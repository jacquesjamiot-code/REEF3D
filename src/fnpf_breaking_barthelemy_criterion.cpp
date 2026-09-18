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

// Angular frequency, crest celerity c = omega/k, criterion B = |u|/c and onset cells.
//
// omega, A384:
//   1  linear dispersion relation, omega = sqrt(g*k*tanh(kh)), h the local water depth
//   2  -dphi/dt by the four-node Lagrange derivative of Wang & Ducrozet (2025, Appendix B)
//        omega = -sum_i phi(t_i) dl_i/dt,  dl_i/dt = sum_{j!=i} 1/(t_i-t_j) prod_{k!=i,j} (t-t_k)/(t_i-t_k)
// The estimate is bounded to [0.5, 1.3] times the dispersion relation.

#include "fdm_fnpf.h"
#include "fnpf_breaking_barthelemy.h"
#include "ghostcell.h"
#include "lexer.h"

// Weights of the Lagrange derivative on the phase increments d1, d2, d3, oldest first.
// d1, d2, d3 are the time intervals; the nodes are t = 0, d1, d1+d2, d1+d2+d3.
static void bart_lagrange_weights(double d1, double d2, double d3, int eval,
                                  double *W, int &ok)
{
    ok = 0;
    W[0] = W[1] = W[2] = 0.0;

    // refuse a collapsed interval
    const double dmin = MIN(d1, MIN(d2, d3));
    const double dmax = MAX(d1, MAX(d2, d3));
    if(dmax <= 0.0 || dmin <= BART_LAG_DTMIN * dmax)
        return;

    const double t[4] = {0.0, d1, d1 + d2, d1 + d2 + d3};
    const double te = t[eval];

    double w[4];
    for(int i4 = 0; i4 < 4; ++i4)
    {
        double sum = 0.0;
        for(int j4 = 0; j4 < 4; ++j4)
        {
            if(j4 == i4)
                continue;

            double prod = 1.0 / (t[i4] - t[j4]);
            for(int k4 = 0; k4 < 4; ++k4)
            {
                if(k4 == i4 || k4 == j4)
                    continue;
                prod *= (te - t[k4]) / (t[i4] - t[k4]);
            }
            sum += prod;
        }
        w[i4] = sum;
    }

    // the weights sum to zero
    double wmax = 0.0;
    for(int i4 = 0; i4 < 4; ++i4)
        wmax = MAX(wmax, fabs(w[i4]));
    if(wmax <= 0.0)
        return;
    if(fabs(w[0] + w[1] + w[2] + w[3]) > 1.0e-6 * wmax)
        return;

    // on phase differences
    W[0] = w[1] + w[2] + w[3];
    W[1] = w[2] + w[3];
    W[2] = w[3];
    ok = 1;
}

// shift the time intervals and update the weights, once per timestep
void fnpf_breaking_barthelemy::lag_advance(lexer *p)
{
    if(p->A384 != 2 || theta_refresh == 0)
        return;

    const double dtn = p->simtime - lag_t_prev;
    lag_t_prev = p->simtime;

    if(theta_ini == 0)
        return;

    lag_dt1 = lag_dt2;
    lag_dt2 = lag_dt3;
    lag_dt3 = dtn;

    if(lag_hist_n < 3)
        lag_hist_n++;

    lag_ok = 0;
    if(lag_hist_n == 3)
        bart_lagrange_weights(lag_dt1, lag_dt2, lag_dt3, BART_LAG_EVAL, lag_w, lag_ok);
}

void fnpf_breaking_barthelemy::omega_celerity(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    SLICELOOP4
    {
        const double k = ch.k(i, j);
        const double h = MAX(c->WL(i, j), 0.01);

        double omega;
        if(p->A384 == 2 && lag_ok)
        {
            const double dphidt = ch.d1(i, j) * lag_w[0] + ch.d2(i, j) * lag_w[1] + ch.d3(i, j) * lag_w[2];
            const double omega_est = -dphidt;
            omega = (omega_est > BART_OMEGA_EPS) ? omega_est : omega_dispersion(k, h);
        }
        else
            omega = omega_dispersion(k, h);

        const double w_ref = omega_dispersion(k, h);
        omega = MAX(omega, BART_OMEGA_LO * w_ref);
        omega = MIN(omega, BART_OMEGA_HI * w_ref);

        ch.omega(i, j) = omega;

        ch.c(i, j) = omega / MAX(k, 1.0e-10);
        ch.c(i, j) = MAX(ch.c(i, j), BART_C_MIN);
        ch.c(i, j) = MIN(ch.c(i, j), BART_C_MAX);
    }
    pgc->gcsl_start4(p, ch.c, 1);
    pgc->gcsl_start4(p, ch.omega, 1);
}

// phases of reference and wet state for the increment of the next timestep
void fnpf_breaking_barthelemy::update_phase_refs(lexer *p, fdm_fnpf *c)
{
    if(theta_refresh)
    {
        SLICELOOP4
        {
            ch.phase_sig_old(i, j) = (*ch.phase_sig)(i, j);
            ch.wet_old(i, j) = bart_wet(p, c, i, j) ? 1 : 0;
        }

        if(pyr_num > 0)
            for(int l = pyr_lev_min; l <= pyr_lev_max + 1; ++l)
            {
                SLICELOOP4
                pyr_level[l]->phase_sig_old(i, j) = pyr_level[l]->phase_sig(i, j);
            }

        theta_ini = 1;
    }
}

// B = |u_surface| / c
void fnpf_breaking_barthelemy::criterion(lexer *p, fdm_fnpf *c)
{
    SLICELOOP4
    {
        const double u_x = c->Fx(i, j);
        const double u_y = c->Fy(i, j);
        const double u_mag = sqrt(u_x * u_x + u_y * u_y);

        const double cp = ch.c(i, j);
        ch.B(i, j) = (cp > BART_OMEGA_EPS) ? u_mag / cp : 0.0;
        ch.B(i, j) = MIN(ch.B(i, j), BART_B_MAX);
    }
}

// cells with B >= B_on, outside the seeding margins A386 and the wet-edge band A387; with
// wetting-drying the cell and its neighbours are wet, the centered derivatives of the cell do
// not read a dry value
void fnpf_breaking_barthelemy::seeds(lexer *p, fdm_fnpf *c, const bart_config &cfg,
                                     std::vector<bart_seed> &sd)
{
    SLICELOOP4
    if(ch.B(i, j) >= cfg.B_on)
    {
        if(p->A343 > 0)
        {
            if(!bart_wet(p, c, i, j) || !bart_wet(p, c, i - 1, j) || !bart_wet(p, c, i + 1, j))
                continue;
            if(p->j_dir == 1 && (!bart_wet(p, c, i, j - 1) || !bart_wet(p, c, i, j + 1)))
                continue;
        }

        // band along the wet edge, A387
        if(wet_edge_dist(i, j) < wet_edge_band)
            continue;

        const int gi = i + p->origin_i;
        const int gj = j + p->origin_j;

        const double d_xs = (double)gi * cfg.dx0;
        const double d_xe = (double)(cfg.Nx - 1 - gi) * cfg.dx0;
        if(cfg.m_xs > 0.0 && d_xs < cfg.m_xs)
            continue;
        if(cfg.m_xe > 0.0 && d_xe < cfg.m_xe)
            continue;

        if(p->j_dir == 1)
        {
            const double d_ys = (double)gj * cfg.dy0;
            const double d_ye = (double)(cfg.Ny - 1 - gj) * cfg.dy0;
            if(cfg.m_ys > 0.0 && d_ys < cfg.m_ys)
                continue;
            if(cfg.m_ye > 0.0 && d_ye < cfg.m_ye)
                continue;
        }

        bart_seed s;
        s.i = i;
        s.j = j;
        s.gi = gi;
        s.gj = gj;
        s.k = ch.k(i, j);
        s.omega = ch.omega(i, j);
        s.theta = ch.theta(i, j);
        s.C = ch.c(i, j);
        s.B = ch.B(i, j);
        s.eta = c->eta(i, j);
        sd.push_back(s);
    }
}
