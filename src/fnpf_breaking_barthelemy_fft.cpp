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

// Spatial Hilbert (A380 1) or Riesz (A380 2) transform of eta by FFT on the global grid.
//   Hilbert: Hx = -i*sign(kx)*eta_hat, Hy = -i*sign(ky)*eta_hat, Hxy = -sign(kx)*sign(ky)*eta_hat
//   Riesz:   fo1 = -i*(kx/|k|)*eta_hat, fo2 = -i*(ky/|k|)*eta_hat
// eta is zero on dry cells and tapered by a Jacobsen window near non-fluid cells before the transform.
// A380 1 and 2 need REEF3D built with USE_FFTW=1.

#include "fdm_fnpf.h"
#include "fnpf_breaking_barthelemy.h"
#include "ghostcell.h"
#include "lexer.h"

#if USE_FFTW
#include <fftw3.h>

void fnpf_breaking_barthelemy::spatial_transforms(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const bool do_H = (p->A380 == 1);
    const bool do_R = (p->A380 == 2);

    const int Nx = p->gknox;
    const int Ny = p->gknoy;
    const int N = Nx * Ny;
    const int Nyh = Ny / 2 + 1;
    const double inv_N = 1.0 / (double)N;

    // global eta on every rank, zero on dry cells
    double *global_eta = new double[N]();
    SLICELOOP4
    global_eta[(i + p->origin_i) * Ny + (j + p->origin_j)] = bart_wet(p, c, i, j) ? c->eta(i, j) : 0.0;

    pgc->globalsumV(global_eta, N);

    // window
    if(!window_ini)
    {
        precompute_window(p, pgc);
        window_ini = 1;
    }
    for(int q = 0; q < N; ++q)
        global_eta[q] *= window[q];

    // forward FFT
    fftw_complex *eta_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
    fftw_plan plan_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, global_eta, eta_hat, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);
    fftw_destroy_plan(plan_fwd);
    delete[] global_eta;

    fftw_complex *Hx_hat = nullptr, *Hy_hat = nullptr, *Hxy_hat = nullptr;
    fftw_complex *fo1_hat = nullptr, *fo2_hat = nullptr;
    double *Hx_real = nullptr, *Hy_real = nullptr, *Hxy_real = nullptr;
    double *fo1_real = nullptr, *fo2_real = nullptr;

    if(do_H)
    {
        Hx_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
        Hy_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
        Hxy_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
        Hx_real = (double *)fftw_malloc(sizeof(double) * N);
        Hy_real = (double *)fftw_malloc(sizeof(double) * N);
        Hxy_real = (double *)fftw_malloc(sizeof(double) * N);
    }
    if(do_R)
    {
        fo1_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
        fo2_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
        fo1_real = (double *)fftw_malloc(sizeof(double) * N);
        fo2_real = (double *)fftw_malloc(sizeof(double) * N);
    }

    // spectral filters; r2c layout: kx = gi or gi-Nx, ky = gj >= 0
    for(int gi = 0; gi < Nx; gi++)
    {
        const int kx = (gi <= Nx / 2) ? gi : gi - Nx;
        const int sx = (kx > 0) ? 1 : (kx < 0) ? -1
                                               : 0;

        for(int gj = 0; gj < Nyh; gj++)
        {
            const int idx = gi * Nyh + gj;
            const int ky = gj;
            const int sy = (ky > 0) ? 1 : 0;

            const double re = eta_hat[idx][0];
            const double im = eta_hat[idx][1];

            if(do_H)
            {
                Hx_hat[idx][0] = (double)sx * im;
                Hx_hat[idx][1] = -(double)sx * re;
                Hy_hat[idx][0] = (double)sy * im;
                Hy_hat[idx][1] = -(double)sy * re;
                Hxy_hat[idx][0] = -(double)(sx * sy) * re;
                Hxy_hat[idx][1] = -(double)(sx * sy) * im;
            }
            if(do_R)
            {
                const double knorm = sqrt((double)(kx * kx + ky * ky));
                double fx, fy;
                if(knorm < 1.0e-10)
                {
                    fx = 0.0;
                    fy = 0.0;
                }
                else
                {
                    fx = (double)kx / knorm;
                    fy = (double)ky / knorm;
                }

                fo1_hat[idx][0] = fx * im;
                fo1_hat[idx][1] = -fx * re;
                fo2_hat[idx][0] = fy * im;
                fo2_hat[idx][1] = -fy * re;
            }
        }
    }

    fftw_free(eta_hat);

    // inverse FFTs, back to the local slices
    if(do_H)
    {
        fftw_plan plan_Hx = fftw_plan_dft_c2r_2d(Nx, Ny, Hx_hat, Hx_real, FFTW_ESTIMATE);
        fftw_plan plan_Hy = fftw_plan_dft_c2r_2d(Nx, Ny, Hy_hat, Hy_real, FFTW_ESTIMATE);
        fftw_plan plan_Hxy = fftw_plan_dft_c2r_2d(Nx, Ny, Hxy_hat, Hxy_real, FFTW_ESTIMATE);
        fftw_execute(plan_Hx);
        fftw_execute(plan_Hy);
        fftw_execute(plan_Hxy);
        fftw_destroy_plan(plan_Hx);
        fftw_destroy_plan(plan_Hy);
        fftw_destroy_plan(plan_Hxy);

        slice4 &HX = *Hx, &HY = *Hy, &HXY = *Hxy;
        SLICELOOP4
        {
            const int gi = i + p->origin_i;
            const int gj = j + p->origin_j;
            HX(i, j) = Hx_real[gi * Ny + gj] * inv_N;
            HY(i, j) = Hy_real[gi * Ny + gj] * inv_N;
            HXY(i, j) = Hxy_real[gi * Ny + gj] * inv_N;
        }

        pgc->gcsl_start4(p, HX, 1);
        pgc->gcsl_start4(p, HY, 1);
        pgc->gcsl_start4(p, HXY, 1);

        fftw_free(Hx_hat);
        fftw_free(Hy_hat);
        fftw_free(Hxy_hat);
        fftw_free(Hx_real);
        fftw_free(Hy_real);
        fftw_free(Hxy_real);
    }

    if(do_R)
    {
        fftw_plan plan_fo1 = fftw_plan_dft_c2r_2d(Nx, Ny, fo1_hat, fo1_real, FFTW_ESTIMATE);
        fftw_plan plan_fo2 = fftw_plan_dft_c2r_2d(Nx, Ny, fo2_hat, fo2_real, FFTW_ESTIMATE);
        fftw_execute(plan_fo1);
        fftw_execute(plan_fo2);
        fftw_destroy_plan(plan_fo1);
        fftw_destroy_plan(plan_fo2);

        slice4 &F1 = *fo1, &F2 = *fo2;
        SLICELOOP4
        {
            const int gi = i + p->origin_i;
            const int gj = j + p->origin_j;
            F1(i, j) = fo1_real[gi * Ny + gj] * inv_N;
            F2(i, j) = fo2_real[gi * Ny + gj] * inv_N;
        }

        pgc->gcsl_start4(p, F1, 1);
        pgc->gcsl_start4(p, F2, 1);

        fftw_free(fo1_hat);
        fftw_free(fo2_hat);
        fftw_free(fo1_real);
        fftw_free(fo2_real);
    }
}
#else
// without FFTW3, bart_check_parameters stops A380 1 and 2 before any call
void fnpf_breaking_barthelemy::spatial_transforms(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}
#endif

// Two-pass chamfer distance [cells] from every cell of the global grid to the nearest non-fluid
// or dry cell, 1e9 if there is none. p->wet is not valid at construction time.
std::vector<double> bart_dry_distance(lexer *p, ghostcell *pgc)
{
    int i, j;  // indices of SLICELOOP4
    const int Nx = p->gknox;
    const int Ny = p->gknoy;
    const int N = Nx * Ny;

    // global fluid mask, 0 or 1: each cell is written by its owner only, so the sum is the mask
    double *mask = new double[N]();
    SLICELOOP4
    mask[(i + p->origin_i) * Ny + (j + p->origin_j)] = (p->wet[IJ] > 0) ? 1.0 : 0.0;
    pgc->globalsumV(mask, N);

    // chamfer distance to the nearest non-fluid cell
    const double LARGE = 1.0e9, d1 = 1.0, dd = sqrt(2.0);
    double *dist = new double[N];
    for(int q = 0; q < N; ++q)
        dist[q] = (mask[q] == 0.0) ? 0.0 : LARGE;

    for(int gi = 0; gi < Nx; ++gi)
        for(int gj = 0; gj < Ny; ++gj)
        {
            const int idx = gi * Ny + gj;
            double v = dist[idx];
            if(gi > 0)
                v = MIN(v, dist[(gi - 1) * Ny + gj] + d1);
            if(gj > 0)
                v = MIN(v, dist[gi * Ny + (gj - 1)] + d1);
            if(gi > 0 && gj > 0)
                v = MIN(v, dist[(gi - 1) * Ny + (gj - 1)] + dd);
            if(gi > 0 && gj < Ny - 1)
                v = MIN(v, dist[(gi - 1) * Ny + (gj + 1)] + dd);
            dist[idx] = v;
        }
    for(int gi = Nx - 1; gi >= 0; --gi)
        for(int gj = Ny - 1; gj >= 0; --gj)
        {
            const int idx = gi * Ny + gj;
            double v = dist[idx];
            if(gi < Nx - 1)
                v = MIN(v, dist[(gi + 1) * Ny + gj] + d1);
            if(gj < Ny - 1)
                v = MIN(v, dist[gi * Ny + (gj + 1)] + d1);
            if(gi < Nx - 1 && gj < Ny - 1)
                v = MIN(v, dist[(gi + 1) * Ny + (gj + 1)] + dd);
            if(gi < Nx - 1 && gj > 0)
                v = MIN(v, dist[(gi + 1) * Ny + (gj - 1)] + dd);
            dist[idx] = v;
        }

    std::vector<double> out(dist, dist + N);

    delete[] dist;
    delete[] mask;

    return out;
}

// Jacobsen relaxation profile over BART_WET_EDGE_FFT cells from the nearest non-fluid or dry cell,
// computed once.
void fnpf_breaking_barthelemy::precompute_window(lexer *p, ghostcell *pgc)
{
    const int N = p->gknox * p->gknoy;
    const double L_relax = BART_WET_EDGE_FFT;
    const std::vector<double> dist = bart_dry_distance(p, pgc);

    // w = 0 at the non-fluid cell, 1 beyond L_relax
    window = new double[N];
    for(int q = 0; q < N; ++q)
    {
        double xn = MAX(1.0 - dist[q] / L_relax, 0.0);
        window[q] = 1.0 - (exp(pow(xn, 3.5)) - 1.0) / (EE - 1.0);
    }
}
