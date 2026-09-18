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

// Breaking dissipation of Wang & Ducrozet (2025, Eq. 1-4), see fnpf_breaking_wang.h.

#include "fdm_fnpf.h"
#include "fnpf_breaking_wang.h"
#include "ghostcell.h"
#include "lexer.h"
#include <cmath>
#include <vector>

// d(eta)/dt at the first RK substage of the timestep, kept for the whole timestep
void fnpf_breaking_wang::freeze_eta_t(lexer *p, slice &eta, slice &eta_n, double alpha,
                                      int theta_refresh)
{
    if(theta_refresh)
        SLICELOOP4 eta_t_wang(i, j) = (eta(i, j) - eta_n(i, j)) / MAX(alpha * p->dt, 1.0e-10);
}

// integral I and nu_eddy of every breaker, onset log
void fnpf_breaking_wang::intensity(lexer *p, fdm_fnpf *c, ghostcell *pgc, const bart_config &cfg, slice &Fifsf,
                                   int theta_refresh, wang_call &wc)
{
    if(wc.n_onset_glob <= 0)
        return;

    if(cfg.diss_on)
    {
        // integrand lap_h(phi_s)*d(eta)/dt on this subdomain
        SLICELOOP4
        {
            // nothing is dissipated on a dry cell
            if(!bart_wet(p, c, i, j))
            {
                G_wang(i, j) = 0.0;
                continue;
            }

            // the laplacian inverted by fnpf_fsfbc::damping(): a dry neighbour takes the value of
            // the cell, no flux across the shoreline
            const double f = Fifsf(i, j);
            const double f_ip = bart_wet(p, c, i + 1, j) ? Fifsf(i + 1, j) : f;
            const double f_im = bart_wet(p, c, i - 1, j) ? Fifsf(i - 1, j) : f;
            const double f_jp = bart_wet(p, c, i, j + 1) ? Fifsf(i, j + 1) : f;
            const double f_jm = bart_wet(p, c, i, j - 1) ? Fifsf(i, j - 1) : f;

            const double lap =
                ((f_ip - f) / (p->DXP[IP] * p->DXN[IP]) + (f_im - f) / (p->DXP[IM1] * p->DXN[IP])) * p->x_dir + ((f_jp - f) / (p->DYP[JP] * p->DYN[JP]) + (f_jm - f) / (p->DYP[JM1] * p->DYN[JP])) * p->y_dir;

            G_wang(i, j) = lap * eta_t_wang(i, j);
        }

        // part of I on this subdomain
        for(wang_breaker &b : wc.brk)
            b.I = zone_sum(p, b.geom);

        // exchange 2: part of I of every foreign zone to its owner, addressed to the onset cell
        std::vector<double> out, in;
        for(const wang_foreign &z : wc.foreign)
        {
            const double gi0 = z.geom.gi0, gj0 = z.geom.gj0;
            const double r[WANG_CONTRIB] = {gi0, gi0, gi0, gi0, gj0, gj0, gj0, gj0,
                                            (double)z.id, zone_sum(p, z.geom)};
            out.insert(out.end(), r, r + WANG_CONTRIB);
        }

        pgc->deliver_quads(p, out, WANG_CONTRIB, in);

        for(size_t k = 0; k < in.size(); k += WANG_CONTRIB)
            wc.brk[(size_t)in[k + 8]].I += in[k + 9];
    }

    for(wang_breaker &b : wc.brk)
    {
        const bart_seed &s = b.seed;

        b.nu = eddy_viscosity(p, b.I, s.C, cfg);

        if(p->P312 > 0 && theta_refresh && b.onset)
        {
            onset_out << p->simtime << " " << p->Xout(p->XP[s.i + marge], p->YP[s.j + marge]) << " "
                      << p->Yout(p->XP[s.i + marge], p->YP[s.j + marge]) << " "
                      << s.eta << " " << s.B << " " << s.k << " " << s.C << " " << s.omega << " " << s.theta << " "
                      << b.Lbr << " " << b.Tbr << " " << (cfg.diss_on ? b.nu : 0.0) << endl;
        }
    }

    // exchange 3: own zones with nu_eddy to the subdomains they meet
    if(cfg.diss_on)
    {
        std::vector<double> out, in;
        for(size_t q = 0; q < wc.brk.size(); ++q)
            wang_pack((int)q, wc.brk[q].geom, wc.brk[q].Tbr, wc.brk[q].nu, out);

        pgc->deliver_quads(p, out, WANG_REC, in);

        wc.foreign.clear();
        for(size_t k = 0; k < in.size(); k += WANG_REC)
            wc.foreign.push_back(wang_unpack(&in[k]));
    }
}

// sum of R*G*dA over the cells of the zone inside this subdomain, in the order gi, then gj
double fnpf_breaking_wang::zone_sum(lexer *p, const wang_zone_geom &g)
{
    // per unit width in 2D
    const double dA = (p->j_dir == 1) ? g.dx0 * g.dy0 : g.dx0;

    std::vector<wang_cell> cells;
    zone_cells(p, g, cells);

    double I = 0.0;
    for(const wang_cell &e : cells)
        I += e.R * G_wang(e.i, e.j) * dA;

    return I;
}

// nu_eddy = -b*C^5*l_cl / (2*g*I); 1 in detection-only mode
double fnpf_breaking_wang::eddy_viscosity(lexer *p, double I, double C, const bart_config &cfg)
{
    if(!cfg.diss_on)
        return 1.0;

    // crest length of one onset cell: the normal plateau
    const double l_cl = (p->j_dir == 1) ? cfg.width : 1.0;

    // I >= 0 would give a negative viscosity
    double nu = 0.0;
    if(I < -1.0e-30)
    {
        nu = -(cfg.b_str * pow(C, 5.0) * l_cl) / (2.0 * 9.81 * I);
        nu = MIN(nu, 0.5 * WANG_VB_MAX);
        nu = MAX(nu, 0.0);
    }

    return nu;
}

// c->vb = 2*nu_eddy*R, applied for the whole event
void fnpf_breaking_wang::apply(lexer *p, fdm_fnpf *c, ghostcell *pgc, const bart_config &cfg,
                               const wang_call &wc)
{
    SLICELOOP4
    {
        double Wg = 0.0;
        if(!wc.stamp.empty())
            Wg = wc.stamp[3 * ((size_t)i * p->knoy + j)];

        if(Wg > 0.0)
            nu_break(i, j) = MAX(nu_break(i, j), Wg);

        if(!(T_breaking(i, j) > 0.0))
            nu_break(i, j) = 0.0;

        c->vb(i, j) = cfg.diss_on ? MIN(2.0 * nu_break(i, j), WANG_VB_MAX) : 0.0;
    }

    pgc->gcsl_start4(p, nu_break, 1);
    pgc->gcsl_start4(p, c->vb, 1);
}
