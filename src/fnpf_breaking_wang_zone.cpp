/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Jacques Amiot
--------------------------------------------------------------------*/

// Breaking zone of Wang & Ducrozet (2025), see fnpf_breaking_wang.h.

#include "fdm_fnpf.h"
#include "fnpf_breaking_wang.h"
#include "ghostcell.h"
#include "lexer.h"
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

fnpf_breaking_wang::fnpf_breaking_wang(lexer *p, fdm_fnpf *c, ghostcell *pgc) : T_breaking(p), t_break_end(p), R_break(p), nu_break(p), eta_t_wang(p), G_wang(p), seed_old(p)
{
    SLICELOOP4
    {
        T_breaking(i, j) = 0.0;
        t_break_end(i, j) = -1.0e20;
        R_break(i, j) = 0.0;
        nu_break(i, j) = 0.0;
        eta_t_wang(i, j) = 0.0;
        G_wang(i, j) = 0.0;
        seed_old(i, j) = 0;
    }

    // onset log, one file per rank
    if(p->P312 > 0)
    {
        mkdir("./REEF3D_FNPF_Breaking_Onset", 0777);

        char name[500];
        sprintf(name, "./REEF3D_FNPF_Breaking_Onset/REEF3D_FNPF-Breaking-Onset-%06i.r3d", p->mpirank + 1);
        onset_out.open(name);

        onset_out << "time x y eta B k c omega theta Lbr Tbr nu_eddy" << endl;
    }
}

fnpf_breaking_wang::~fnpf_breaking_wang()
{
    if(onset_out.is_open())
        onset_out.close();
}

// Jacobsen ramp: 0 at the outer edge, 1 beyond L
static inline double wang_ramp(double d, double L)
{
    if(L <= 1.0e-12)
        return (d > 0.0) ? 1.0 : 0.0;
    const double xn = MAX(1.0 - d / L, 0.0);
    return 1.0 - (exp(pow(xn, 3.5)) - 1.0) / (EE - 1.0);
}

// Arbitration of the triples (W = nu_eddy*R, T = Tbr, R): largest W, then T, then R.
// A total order, so the result does not depend on the number of ranks.
static inline bool wang_wins(const double *a, const double *b)
{
    if(a[0] != b[0])
        return a[0] > b[0];
    if(a[1] != b[1])
        return a[1] > b[1];
    return a[2] > b[2];
}

// length, duration and rectangle of every onset cell
void fnpf_breaking_wang::breakers(lexer *p, ghostcell *pgc, const bart_config &cfg,
                                  const std::vector<bart_seed> &seeds, int theta_refresh, wang_call &wc)
{
    wc.n_onset_glob = pgc->globalisum((int)seeds.size());

    // new onsets for the log: cells that did not seed at the previous timestep
    std::vector<int> onset(seeds.size(), 0);
    if(p->P312 > 0 && theta_refresh)
    {
        for(size_t q = 0; q < seeds.size(); ++q)
            onset[q] = (seed_old(seeds[q].i, seeds[q].j) == 0) ? 1 : 0;

        SLICELOOP4
        seed_old(i, j) = 0;

        for(size_t q = 0; q < seeds.size(); ++q)
            seed_old(seeds[q].i, seeds[q].j) = 1;
    }

    if(wc.n_onset_glob <= 0)
        return;

    wc.brk.reserve(seeds.size());
    for(size_t q = 0; q < seeds.size(); ++q)
    {
        const bart_seed &s = seeds[q];
        wang_breaker b;
        b.seed = s;
        b.onset = onset[q];

        const double k_p = MAX(s.k, 1.0e-10);
        const double w_p = MAX(s.omega, 1.0e-10);

        b.Lbr = 2.0 * PI / k_p;
        b.Lbr = MIN(b.Lbr, (double)cfg.Nx * cfg.dx0);

        b.Tbr = p->A383 * 2.0 * PI / w_p;

        const double Tref = MAX(p->wTp, p->wT);
        if(Tref > 0.0)
            b.Tbr = MIN(b.Tbr, WANG_TBR_MAX_TP * Tref);

        b.geom = zone_geometry(p, s.gi, s.gj, s.theta, b.Lbr, cfg);
        b.nu = 1.0;
        b.I = 0.0;

        wc.brk.push_back(b);
    }

    // exchange 1: own zones to the subdomains they meet; nu is set by intensity
    std::vector<double> out, in;
    for(size_t q = 0; q < wc.brk.size(); ++q)
        wang_pack((int)q, wc.brk[q].geom, wc.brk[q].Tbr, 1.0, out);

    pgc->deliver_quads(p, out, WANG_REC, in);

    for(size_t k = 0; k < in.size(); k += WANG_REC)
        wc.foreign.push_back(wang_unpack(&in[k]));
}

// rectangle of one breaker, its corners in global cell indices
wang_zone_geom fnpf_breaking_wang::zone_geometry(lexer *p, int gi0, int gj0, double theta,
                                                 double Lbr, const bart_config &cfg)
{
    wang_zone_geom g;
    g.gi0 = gi0;
    g.gj0 = gj0;
    g.ct = cos(theta);
    g.st = sin(theta);
    g.dx0 = cfg.dx0;
    g.dy0 = cfg.dy0;

    g.Lrp_t = WANG_ZONE_RAMP_T * Lbr;
    g.Lrp_n = cfg.Lrp_n;

    g.s0 = -g.Lrp_t;
    g.s1 = Lbr + g.Lrp_t;
    g.t1 = 0.5 * cfg.width + g.Lrp_n;

    // corners counterclockwise; in 2D the zone has no width and lies on the row gj0
    const double cs[4] = {g.s0, g.s1, g.s1, g.s0};
    const double ctt[4] = {-g.t1, -g.t1, g.t1, g.t1};
    for(int k = 0; k < 4; ++k)
    {
        const double t = (p->j_dir == 1) ? ctt[k] : 0.0;
        g.ci[k] = gi0 + (cs[k] * g.ct - t * g.st) / g.dx0;
        g.cj[k] = (p->j_dir == 1) ? gj0 + (cs[k] * g.st + t * g.ct) / g.dy0 : (double)gj0;
    }

    return g;
}

// ramp weight R = R(r^t)*R(r^n) at the centre of a global cell, 0 outside the zone
double fnpf_breaking_wang::zone_ramp_weight(lexer *p, const wang_zone_geom &g, int gi, int gj)
{
    const double Xc = (double)(gi - g.gi0) * g.dx0;
    const double Yc = (p->j_dir == 1) ? (double)(gj - g.gj0) * g.dy0 : 0.0;

    const double sc = Xc * g.ct + Yc * g.st;
    const double tc = -Xc * g.st + Yc * g.ct;

    const double d_t = MIN(sc - g.s0, g.s1 - sc);
    if(d_t <= 0.0)
        return 0.0;

    double R = wang_ramp(d_t, g.Lrp_t);

    if(p->j_dir == 1)
    {
        const double d_n = g.t1 - fabs(tc);
        if(d_n <= 0.0)
            return 0.0;
        R *= wang_ramp(d_n, g.Lrp_n);
    }

    return R;
}

// cells of the zone inside this subdomain with R > 0, in the order gi, then gj. The loop bounds
// come from the four edges in cell indices, moved one cell outwards; R decides.
void fnpf_breaking_wang::zone_cells(lexer *p, const wang_zone_geom &g, std::vector<wang_cell> &cells)
{
    cells.clear();

    const double *ci = g.ci, *cj = g.cj;
    const double i_lo = p->origin_i, i_hi = p->origin_i + p->knox - 1;
    const double j_lo = p->origin_j, j_hi = p->origin_j + p->knoy - 1;

    // columns crossed by the rectangle
    const int i0 = (int)ceil(MAX(*std::min_element(ci, ci + 4) - 1.0, i_lo));
    const int i1 = (int)floor(MIN(*std::max_element(ci, ci + 4) + 1.0, i_hi));

    for(int gi = i0; gi <= i1; ++gi)
    {
        // rows of the column on the inner side of every edge, within one cell
        double lo = j_lo, hi = j_hi;
        for(int k = 0; k < 4; ++k)
        {
            const int m = (k + 1) % 4;
            const double ei = ci[m] - ci[k], ej = cj[m] - cj[k];
            const double L = sqrt(ei * ei + ej * ej);
            if(L == 0.0 || ei == 0.0)
                continue;

            // ei*(gj - cj[k]) - ej*(gi - ci[k]) >= -L
            const double b = cj[k] + (ej * (gi - ci[k]) - L) / ei;
            if(ei > 0.0)
                lo = MAX(lo, b);
            else
                hi = MIN(hi, b);
        }
        if(lo > hi)
            continue;

        for(int gj = (int)ceil(lo); gj <= (int)floor(hi); ++gj)
        {
            const double R = zone_ramp_weight(p, g, gi, gj);
            if(R > 0.0)
                cells.push_back({gi - p->origin_i, gj - p->origin_j, R});
        }
    }
}

// the triple (nu*R, Tbr, R) of one zone on the cells of this subdomain, the winning triple as a whole
void fnpf_breaking_wang::zone_stamp(lexer *p, const wang_zone_geom &g, double nu, double Tbr,
                                    double *stamp)
{
    std::vector<wang_cell> cells;
    zone_cells(p, g, cells);

    for(const wang_cell &e : cells)
    {
        const size_t q = 3 * ((size_t)e.i * p->knoy + e.j);
        const double cand[3] = {nu * e.R, Tbr, e.R};
        if(wang_wins(cand, stamp + q))
        {
            stamp[q] = cand[0];
            stamp[q + 1] = cand[1];
            stamp[q + 2] = cand[2];
        }
    }
}

// stamping of the own and foreign zones on this subdomain, breaking duration and decision
void fnpf_breaking_wang::zone(lexer *p, fdm_fnpf *c, ghostcell *pgc, const bart_config &cfg,
                              wang_call &wc)
{
    if(wc.n_onset_glob > 0)
    {
        wc.stamp.assign(3 * (size_t)p->knox * p->knoy, 0.0);

        for(const wang_breaker &b : wc.brk)
            zone_stamp(p, b.geom, b.nu, b.Tbr, wc.stamp.data());

        for(const wang_foreign &z : wc.foreign)
            zone_stamp(p, z.geom, z.nu, z.Tbr, wc.stamp.data());
    }

    SLICELOOP4
    {
        double Tg = 0.0, Rg = 0.0;
        if(!wc.stamp.empty())
        {
            const size_t q = 3 * ((size_t)i * p->knoy + j);
            Tg = wc.stamp[q + 1];
            Rg = wc.stamp[q + 2];
        }

        // breaking duration: latest end time of all zones stamped on this cell
        if(Rg > 0.0)
            t_break_end(i, j) = MAX(t_break_end(i, j), p->simtime + Tg);

        T_breaking(i, j) = MAX(t_break_end(i, j) - p->simtime, 0.0);

        // ramp R kept for the whole event
        if(Rg > 0.0)
            R_break(i, j) = MAX(R_break(i, j), Rg);

        if(!(T_breaking(i, j) > 0.0))
            R_break(i, j) = 0.0;

        c->breaking(i, j) = (T_breaking(i, j) > 0.0) ? 1 : 0;
    }

    pgc->gcsl_start4int(p, c->breaking, 50);
    pgc->gcsl_start4(p, T_breaking, 1);
    pgc->gcsl_start4(p, R_break, 1);
}
