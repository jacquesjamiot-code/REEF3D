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

// Breaking onset criterion, A350 4.
// X. Barthelemy, M.L. Banner, W.L. Peirson, F. Fedele, M. Allis, F. Dias, "On a unified breaking
// onset threshold for gravity waves in deep and intermediate depth water", J. Fluid Mech. 841,
// 463-488, 2018.
//
//   B = |u_s| / c >= B_on (0.855)
//
// u_s is the horizontal velocity at the free surface and c the local crest celerity omega/k.
// The local wavenumber vector k and the direction theta come from a spatial transform of eta,
// A380:
//   1  Hilbert transform by FFT. Y. Wang, G. Ducrozet, "Computing local crest phase speed in
//      multi-directional sea states through the spatial Hilbert transform", Ocean Eng. 331,
//      121357, 2025. Valid for slowly varying wave trains.
//   2  Riesz transform by FFT, monogenic signal. M. Felsberg, G. Sommer, "The monogenic signal",
//      IEEE Trans. Signal Process. 49(12), 2001.
//   3  Riesz transform rebuilt from the Riesz pyramid of Wadhwa et al. (2014), local
//      communication only.
// The cells reaching B_on seed the breaking zones and the dissipation of fnpf_breaking_wang.

#include "fnpf_breaking_barthelemy.h"
#include "fdm_fnpf.h"
#include "fnpf_breaking.h"
#include "fnpf_breaking_wang.h"
#include "ghostcell.h"
#include "lexer.h"
#include <iostream>

void fnpf_breaking::breaking_barthelemy(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    pbart->step_begin(p);

    // Fx, Fy are not available at the first timestep
    if(p->count <= 1)
        return;

    bart_config cfg = pbart->config(p);
    std::vector<bart_seed> seeds;
    wang_call wc;

    pbart->diagnostic(p, c, pgc, cfg, seeds);

    // intensity before zone: overlapping zones keep the largest nu_eddy*R
    pwang->breakers(p, pgc, cfg, seeds, pbart->theta_refresh, wc);
    pwang->freeze_eta_t(p, eta, eta_n, alpha, pbart->theta_refresh);
    pwang->intensity(p, pgc, cfg, Fifsf, pbart->theta_refresh, wc);
    pwang->zone(p, c, pgc, cfg, wc);
    pwang->apply(p, c, pgc, cfg, wc);

    SLICELOOP4
    c->breaklog(i, j) = 0;

    // breaklog
    int count = 0;
    SLICELOOP4
    if(c->breaking(i, j) > 0)
    {
        c->breaklog(i, j) = 1;
        ++count;
    }

    count = pgc->globalisum(count);

    if(p->mpirank == 0 && (p->count % p->P12 == 0))
        cout << "breaking: " << count << endl;
}

void bart_check_parameters(lexer *p, ghostcell *pgc)
{
    const char *name = nullptr, *detail = nullptr;

    // the first failed check is reported
    auto check = [&](bool fail, const char *par, const char *msg)
    {
        if(fail && !name)
        {
            name = par;
            detail = msg;
        }
    };

    if(p->A350 == 4)
    {
        check(p->A380 < 1 || p->A380 > 3, "A 380", "(transform) expected 1, 2 or 3");
#if !USE_FFTW
        check(p->A380 == 1 || p->A380 == 2, "A 380",
              "1 and 2 (FFT transforms) require REEF3D built with USE_FFTW=1, or use A 380 3");
#endif
        check(p->A381 <= 0.0, "A 381", "(onset threshold) must be positive");
        check(p->A382 < 0.0, "A 382", "(dissipation strength) cannot be negative");
        check(p->A383 <= 0.0, "A 383", "(factor on breaking duration) must be positive");
        check(p->A384 < 1 || p->A384 > 2, "A 384", "(omega estimate) expected 1 or 2");
        check(p->A385_nmin < 0 || p->A385_nmin > p->A385_nmax || p->A385_nmax > RIESZ_PYR_MAX - 2,
              "A 385", "(pyramid levels) requires 0 <= N_min <= N_max <= 5");
        check(p->A386_xs < 0.0 || p->A386_xe < 0.0 || p->A386_ys < 0.0 || p->A386_ye < 0.0,
              "A 386", "(seeding margins) cannot be negative");
    }

    check(p->P312 < 0 || p->P312 > 1, "P 312", "expected 0 or 1");
    check(p->P313 < 0 || p->P313 > 1, "P 313", "expected 0 or 1");
    check(p->P314 < 0 || p->P314 > 1, "P 314", "expected 0 or 1");

    check(p->P313 > 0 && (p->A350 != 4 || p->A380 != 3), "P 313",
          "(pyramid levels) requires A 350 4 and A 380 3");
    check(p->P314 > 0 && (p->A350 != 4 || p->P51 <= 0), "P 314",
          "(breaking criterion at gauges) requires A 350 4 and P 51 gauges");

    if(name)
    {
        if(p->mpirank == 0)
            cout << "\n"
                 << "!!! wrong input error for " << name << " !!!\n\n"
                 << name << " " << detail << "\n\n"
                 << "!!! please check the REEF3D User Guide !!!\n\n\n"
                 << endl;

        pgc->final(true);
    }
}

bart_chain::bart_chain(lexer *p) : A(p), phase(p), k(p), theta(p),
                                   phase_sig(&phase_sig_own), phase_sig_own(p), phase_sig_old(p), dphi_filt(p),
                                   d1(p), d2(p), d3(p), omega(p), c(p), B(p)
{
}

bart_pyr_level::bart_pyr_level(lexer *p) : k(p), fo1(p), fo2(p), A(p), phase(p), band(p),
                                           theta(p), phase_sig(p), phase_sig_old(p), dphi_filt(p)
{
}

fnpf_breaking_barthelemy::fnpf_breaking_barthelemy(lexer *p, fdm_fnpf *c, ghostcell *pgc) : ch(p)
{
    Hx = Hy = Hxy = fo1 = fo2 = nullptr;

    if(p->A380 == 1)
    {
        Hx = new slice4(p);
        Hy = new slice4(p);
        Hxy = new slice4(p);
        ch.phase_sig = &ch.phase;
    }
    else
    {
        fo1 = new slice4(p);
        fo2 = new slice4(p);
    }

    SLICELOOP4
    {
        ch.A(i, j) = 0.0;
        ch.phase(i, j) = 0.0;
        ch.k(i, j) = 0.0;
        ch.theta(i, j) = 0.0;
        ch.phase_sig_own(i, j) = 0.0;
        ch.phase_sig_old(i, j) = 0.0;
        ch.dphi_filt(i, j) = 0.0;
        ch.d1(i, j) = 0.0;
        ch.d2(i, j) = 0.0;
        ch.d3(i, j) = 0.0;
        ch.omega(i, j) = 0.0;
        ch.c(i, j) = 0.0;
        ch.B(i, j) = 0.0;
    }

    if(p->A380 == 1)
    {
        SLICELOOP4
        {
            (*Hx)(i, j) = 0.0;
            (*Hy)(i, j) = 0.0;
            (*Hxy)(i, j) = 0.0;
        }
    }
    else
    {
        SLICELOOP4
        {
            (*fo1)(i, j) = 0.0;
            (*fo2)(i, j) = 0.0;
        }
    }

    window = nullptr;
    window_ini = 0;

    lag_dt1 = lag_dt2 = lag_dt3 = 0.0;
    lag_t_prev = 0.0;
    lag_hist_n = 0;
    lag_w[0] = lag_w[1] = lag_w[2] = 0.0;
    lag_ok = 0;

    theta_count = -1;
    theta_refresh = 0;
    theta_ini = 0;

    // pyramid levels, and the levels exported to the vtp with P313
    pyr_lev_min = p->A385_nmin;
    pyr_lev_max = p->A385_nmax;
    pyr_num = 0;
    pyr_gather_level = -1;
    pyr_setup_done = 0;

    for(int l = 0; l < RIESZ_PYR_MAX; ++l)
        pyr_level[l] = nullptr;

    exp_A = exp_phase = exp_k = exp_theta = exp_band = exp_fo1 = exp_fo2 = nullptr;

    if(p->A380 == 3 && p->P313 > 0)
    {
        pyr_num = pyr_lev_max - pyr_lev_min + 2;

        exp_A = new slice4 *[pyr_num];
        exp_phase = new slice4 *[pyr_num];
        exp_k = new slice4 *[pyr_num];
        exp_theta = new slice4 *[pyr_num];
        exp_band = new slice4 *[pyr_num];
        exp_fo1 = new slice4 *[pyr_num];
        exp_fo2 = new slice4 *[pyr_num];

        for(int l = pyr_lev_min; l <= pyr_lev_max + 1; ++l)
        {
            bart_pyr_level *LV = new bart_pyr_level(p);
            pyr_level[l] = LV;

            SLICELOOP4
            {
                LV->k(i, j) = 0.0;
                LV->fo1(i, j) = 0.0;
                LV->fo2(i, j) = 0.0;
                LV->A(i, j) = 0.0;
                LV->phase(i, j) = 0.0;
                LV->band(i, j) = 0.0;
                LV->theta(i, j) = 0.0;
                LV->phase_sig(i, j) = 0.0;
                LV->phase_sig_old(i, j) = 0.0;
                LV->dphi_filt(i, j) = 0.0;
            }

            const int n = l - pyr_lev_min;
            exp_A[n] = &LV->A;
            exp_phase[n] = &LV->phase;
            exp_k[n] = &LV->k;
            exp_theta[n] = &LV->theta;
            exp_band[n] = &LV->band;
            exp_fo1[n] = &LV->fo1;
            exp_fo2[n] = &LV->fo2;
        }
    }

    // fields read by the printer
    c->brk_B = &ch.B;
    c->brk_k = &ch.k;
    c->brk_c = &ch.c;
    c->brk_omega = &ch.omega;

    c->brk_pyr_num = pyr_num;
    c->brk_pyr_first = pyr_lev_min;
    c->brk_pyr_A = exp_A;
    c->brk_pyr_phase = exp_phase;
    c->brk_pyr_k = exp_k;
    c->brk_pyr_theta = exp_theta;
    c->brk_pyr_band = exp_band;
    c->brk_pyr_fo1 = exp_fo1;
    c->brk_pyr_fo2 = exp_fo2;
}

fnpf_breaking_barthelemy::~fnpf_breaking_barthelemy()
{
    delete Hx;
    delete Hy;
    delete Hxy;
    delete fo1;
    delete fo2;
    delete[] window;

    for(int l = 0; l < RIESZ_PYR_MAX; ++l)
        delete pyr_level[l];

    delete[] exp_A;
    delete[] exp_phase;
    delete[] exp_k;
    delete[] exp_theta;
    delete[] exp_band;
    delete[] exp_fo1;
    delete[] exp_fo2;
}

// once per call; the phase increment is measured at the first RK substage of a timestep
void fnpf_breaking_barthelemy::step_begin(lexer *p)
{
    theta_refresh = (p->count != theta_count) ? 1 : 0;
    if(theta_refresh)
        theta_count = p->count;

    // before the phases of this call
    lag_advance(p);
}

bart_config fnpf_breaking_barthelemy::config(lexer *p)
{
    bart_config cfg;

    cfg.B_on = p->A381;
    cfg.diss_on = (p->A382 > 0.0);
    cfg.b_str = cfg.diss_on ? p->A382 : 0.0;

    // p->dx is not set on the FNPF path
    cfg.dx0 = p->DXP[marge + p->knox / 2];
    if(cfg.dx0 <= 1.0e-20)
        cfg.dx0 = 1.0;
    cfg.dy0 = (p->j_dir == 1) ? p->DYP[marge + p->knoy / 2] : cfg.dx0;
    if(cfg.dy0 <= 1.0e-20)
        cfg.dy0 = cfg.dx0;

    const double hmin = MIN(cfg.dx0, cfg.dy0);
    cfg.width = WANG_ZONE_WIDTH * hmin;
    cfg.Lrp_n = WANG_ZONE_RAMP_N * hmin;

    cfg.Nx = p->gknox;
    cfg.Ny = p->gknoy;
    cfg.N = cfg.Nx * cfg.Ny;

    cfg.m_xs = p->A386_xs;
    cfg.m_xe = p->A386_xe;
    cfg.m_ys = p->A386_ys;
    cfg.m_ye = p->A386_ye;

    return cfg;
}

// spatial transform, local variables, criterion, onset cells
void fnpf_breaking_barthelemy::diagnostic(lexer *p, fdm_fnpf *c, ghostcell *pgc, const bart_config &cfg,
                                          std::vector<bart_seed> &sd)
{
    if(p->A380 == 1 || p->A380 == 2)
        spatial_transforms(p, c, pgc);

    if(p->A380 == 3)
        riesz_pyramid(p, c, pgc, cfg);

    if(p->A380 == 1)
        local_hilbert(p, c, pgc);
    else
        local_riesz(p, c, pgc, *fo1, *fo2);

    if(pyr_num > 0)
        pyramid_level_theta(p);

    // omega reads the previous phases, update_phase_refs overwrites them
    omega_celerity(p, c, pgc);
    update_phase_refs(p);

    criterion(p, c);
    seeds(p, c, cfg, sd);
}
