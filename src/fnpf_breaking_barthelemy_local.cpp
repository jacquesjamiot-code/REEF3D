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

// Local amplitude, phase, wavenumber and propagation direction.
// The direction is the axis of the wavenumber vector, turned by pi when the filtered phase
// increment is positive: for a progressive wave the phase decreases in time.

#include "fdm_fnpf.h"
#include "fnpf_breaking_barthelemy.h"
#include "ghostcell.h"
#include "lexer.h"

// Hilbert method, Wang & Ducrozet (Ocean Eng. 331, 2025): the wavenumber vector is the gradient
// of the phase phase12 built from the two half-signals (eta -/+ Hxy, Hx +/- Hy).
void fnpf_breaking_barthelemy::local_hilbert(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    slice4 &HX = *Hx, &HY = *Hy, &HXY = *Hxy;

    SLICELOOP4
    {
        const double Hx_c = HX(i, j);
        const double Hy_c = HY(i, j);
        const double Hxy_c = HXY(i, j);

        const double dHxdx = centered_dx(HX, i, j, p, marge);
        const double dHydx = centered_dx(HY, i, j, p, marge);
        const double dHxydx = centered_dx(HXY, i, j, p, marge);

        const double eta_c = c->eta(i, j);
        const double detadx = centered_dx(c->eta, i, j, p, marge);

        double dHxdy = 0.0, dHydy = 0.0, dHxydy = 0.0, detady = 0.0;
        if(p->j_dir == 1)
        {
            dHxdy = centered_dy(HX, i, j, p, marge);
            dHydy = centered_dy(HY, i, j, p, marge);
            dHxydy = centered_dy(HXY, i, j, p, marge);
            detady = centered_dy(c->eta, i, j, p, marge);
        }

        // amplitudes and phases of the two half-signals
        double A1 = sqrt((eta_c - Hxy_c) * (eta_c - Hxy_c) + (Hx_c + Hy_c) * (Hx_c + Hy_c));
        double A2 = sqrt((eta_c + Hxy_c) * (eta_c + Hxy_c) + (Hx_c - Hy_c) * (Hx_c - Hy_c));
        double A1s = MAX(A1, 1.0e-10);
        double A2s = MAX(A2, 1.0e-10);

        double cosphase1 = (eta_c - Hxy_c) / A1s;
        double cosphase2 = (eta_c + Hxy_c) / A2s;
        double sinphase1 = (Hx_c + Hy_c) / A1s;
        double sinphase2 = (Hx_c - Hy_c) / A2s;

        double dphase1dx = (cosphase1 * (dHxdx + dHydx) - sinphase1 * (detadx - dHxydx)) / A1s;
        double dphase1dy = (cosphase1 * (dHxdy + dHydy) - sinphase1 * (detady - dHxydy)) / A1s;

        double dphase2dx = (cosphase2 * (dHxdx - dHydx) - sinphase2 * (detadx + dHxydx)) / A2s;
        double dphase2dy = (cosphase2 * (dHxdy - dHydy) - sinphase2 * (detady + dHxydy)) / A2s;

        double dA1dx = cosphase1 * (detadx - dHxydx) + sinphase1 * (dHxdx + dHydx);
        double dA1dy = cosphase1 * (detady - dHxydy) + sinphase1 * (dHxdy + dHydy);

        double dA2dx = cosphase2 * (detadx + dHxydx) + sinphase2 * (dHxdx - dHydx);
        double dA2dy = cosphase2 * (detady + dHxydy) + sinphase2 * (dHxdy - dHydy);

        // phase12 = atan2(N, D) and its gradient
        double N = A1 * sinphase1 + A2 * sinphase2;
        double D = A1 * cosphase1 + A2 * cosphase2;

        double dNdx = dA1dx * sinphase1 + A1 * cosphase1 * dphase1dx + dA2dx * sinphase2 + A2 * cosphase2 * dphase2dx;
        double dNdy = dA1dy * sinphase1 + A1 * cosphase1 * dphase1dy + dA2dy * sinphase2 + A2 * cosphase2 * dphase2dy;

        double dDdx = dA1dx * cosphase1 - A1 * sinphase1 * dphase1dx + dA2dx * cosphase2 - A2 * sinphase2 * dphase2dx;
        double dDdy = dA1dy * cosphase1 - A1 * sinphase1 * dphase1dy + dA2dy * cosphase2 - A2 * sinphase2 * dphase2dy;

        double phase12 = atan2(N, D);
        double ND2 = MAX(D * D + N * N, 1.0e-20);
        double dphase12dx = (dNdx * D - N * dDdx) / ND2;
        double dphase12dy = (dNdy * D - N * dDdy) / ND2;

        ch.phase(i, j) = phase12;
        ch.A(i, j) = 0.5 * sqrt(N * N + D * D);

        double k = sqrt(dphase12dx * dphase12dx + dphase12dy * dphase12dy);
        ch.k(i, j) = k;

        // axis of the wavenumber vector, in ]-pi/2, pi/2]
        double th = (fabs(dphase12dx) > 1.0e-12)
                        ? atan(dphase12dy / dphase12dx)
                        : ((dphase12dy >= 0.0) ? 0.5 * PI : -0.5 * PI);

        // phase increment, once per timestep
        if(theta_refresh && theta_ini)
        {
            double dphi = ch.phase(i, j) - ch.phase_sig_old(i, j);
            while(dphi > PI)
                dphi -= 2.0 * PI;
            while(dphi < -PI)
                dphi += 2.0 * PI;
            ch.dphi_filt(i, j) = (1.0 - BART_DIR_FILTER) * ch.dphi_filt(i, j) + BART_DIR_FILTER * dphi;

            if(p->A384 == 2)
            {
                ch.d1(i, j) = ch.d2(i, j);
                ch.d2(i, j) = ch.d3(i, j);
                ch.d3(i, j) = dphi;
            }
        }

        if(ch.dphi_filt(i, j) > 0.0)
            th = (th > 0.0) ? (th - PI) : (th + PI);

        ch.theta(i, j) = th;
    }

    pgc->gcsl_start4(p, ch.phase, 1);
    pgc->gcsl_start4(p, ch.A, 1);
    pgc->gcsl_start4(p, ch.k, 1);
    pgc->gcsl_start4(p, ch.theta, 1);
}

// Riesz method, monogenic signal (Felsberg & Sommer 2001):
//   A = sqrt(eta^2 + |fo|^2),  phase = atan2(|fo|, eta),  k = |grad phase|
// with (fo1, fo2) the Riesz pair from the FFT or from the pyramid.
void fnpf_breaking_barthelemy::local_riesz(lexer *p, fdm_fnpf *c, ghostcell *pgc,
                                           slice4 &fo1_src, slice4 &fo2_src)
{
    SLICELOOP4
    {
        const double fo1IJ = fo1_src(i, j);
        const double fo2IJ = fo2_src(i, j);
        const double feIJ = c->eta(i, j);

        const double dfo1_dx = centered_dx(fo1_src, i, j, p, marge);
        const double dfo2_dx = centered_dx(fo2_src, i, j, p, marge);
        const double dfe_dx = centered_dx(c->eta, i, j, p, marge);

        double dfo1_dy = 0.0, dfo2_dy = 0.0, dfe_dy = 0.0;
        if(p->j_dir == 1)
        {
            dfo1_dy = centered_dy(fo1_src, i, j, p, marge);
            dfo2_dy = centered_dy(fo2_src, i, j, p, marge);
            dfe_dy = centered_dy(c->eta, i, j, p, marge);
        }

        // |fo| and its gradient
        const double foIJ = sqrt(fo1IJ * fo1IJ + fo2IJ * fo2IJ);
        const double foSafe = MAX(foIJ, 1.0e-10);
        const double dfo_dx = (fo1IJ * dfo1_dx + fo2IJ * dfo2_dx) / foSafe;
        const double dfo_dy = (fo1IJ * dfo1_dy + fo2IJ * dfo2_dy) / foSafe;

        const double A_R = sqrt(feIJ * feIJ + foIJ * foIJ);
        const double local_phase = atan2(foIJ, feIJ);

        const double A2safe = MAX(A_R * A_R, 1.0e-20);
        const double dphase_dx = (feIJ * dfo_dx - foIJ * dfe_dx) / A2safe;
        const double dphase_dy = (feIJ * dfo_dy - foIJ * dfe_dy) / A2safe;

        ch.A(i, j) = A_R;
        ch.phase(i, j) = local_phase;
        ch.k(i, j) = sqrt(dphase_dx * dphase_dx + dphase_dy * dphase_dy);

        // axis from fo2/fo1 = ky/kx
        double th = (fabs(fo1IJ) > 1.0e-12)
                        ? atan(fo2IJ / fo1IJ)
                        : ((fo2IJ >= 0.0) ? 0.5 * PI : -0.5 * PI);

        // signed phase: the Riesz vector projected on the axis
        const double q = fo1IJ * cos(th) + fo2IJ * sin(th);
        (*ch.phase_sig)(i, j) = atan2(q, feIJ);

        // phase increment, once per timestep
        if(theta_refresh && theta_ini)
        {
            double dphi = (*ch.phase_sig)(i, j) - ch.phase_sig_old(i, j);
            while(dphi > PI)
                dphi -= 2.0 * PI;
            while(dphi < -PI)
                dphi += 2.0 * PI;
            ch.dphi_filt(i, j) = (1.0 - BART_DIR_FILTER) * ch.dphi_filt(i, j) + BART_DIR_FILTER * dphi;

            if(p->A384 == 2)
            {
                ch.d1(i, j) = ch.d2(i, j);
                ch.d2(i, j) = ch.d3(i, j);
                ch.d3(i, j) = dphi;
            }
        }

        if(ch.dphi_filt(i, j) > 0.0)
            th = (th > 0.0) ? (th - PI) : (th + PI);

        ch.theta(i, j) = th;
    }

    pgc->gcsl_start4(p, ch.A, 1);
    pgc->gcsl_start4(p, ch.phase, 1);
    pgc->gcsl_start4(p, ch.k, 1);
    pgc->gcsl_start4(p, ch.theta, 1);
}

// direction of every pyramid level, P313, same construction on the level's own fields
void fnpf_breaking_barthelemy::pyramid_level_theta(lexer *p)
{
    for(int l = pyr_lev_min; l <= pyr_lev_max + 1; ++l)
    {
        bart_pyr_level &LV = *pyr_level[l];

        SLICELOOP4
        {
            const double f1 = LV.fo1(i, j), f2 = LV.fo2(i, j);

            double th = (fabs(f1) > 1.0e-12)
                            ? atan(f2 / f1)
                            : ((f2 >= 0.0) ? 0.5 * PI : -0.5 * PI);

            const double q = f1 * cos(th) + f2 * sin(th);
            LV.phase_sig(i, j) = atan2(q, LV.band(i, j));

            if(theta_refresh && theta_ini)
            {
                double dphi = LV.phase_sig(i, j) - LV.phase_sig_old(i, j);
                while(dphi > PI)
                    dphi -= 2.0 * PI;
                while(dphi < -PI)
                    dphi += 2.0 * PI;
                LV.dphi_filt(i, j) = (1.0 - BART_DIR_FILTER) * LV.dphi_filt(i, j) + BART_DIR_FILTER * dphi;
            }

            if(LV.dphi_filt(i, j) > 0.0)
                th = (th > 0.0) ? (th - PI) : (th + PI);

            LV.theta(i, j) = th;
        }
    }
}
