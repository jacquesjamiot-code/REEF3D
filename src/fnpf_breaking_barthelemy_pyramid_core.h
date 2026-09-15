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

// Riesz pyramid, numerical core.
// N. Wadhwa, M. Rubinstein, F. Durand, W.T. Freeman, "Riesz pyramids for fast phase-based
// video magnification", IEEE ICCP 2014, and its supplemental material.
//
// One level splits the field low into a high-pass band H9*low and a low-pass L9*low. The Riesz
// transform of the band is taken with the 3x3 kernels FO1 and FO2 = FO1^T, and the low-pass
// is decimated by 2 to give the next level. The 9x9 kernels are the published 1-D filters
// (Table 1a) through the McClellan transform, so a run with one cell in y applies the 1-D
// filters directly, together with the published 3x1 Riesz kernel (Table 2a).
//
// No lexer, ghostcell or MPI. The operators write the interior of their output only and never
// a halo: the caller fills the halo by exchange, or leaves it at zero (zero padding).

#ifndef FNPF_BREAKING_BARTHELEMY_PYRAMID_CORE_H_
#define FNPF_BREAKING_BARTHELEMY_PYRAMID_CORE_H_

#include <cmath>
#include <cstddef>
#include <vector>

// highest level index + 2: the band stops at RIESZ_PYR_MAX-2, the residual takes the next slot
#define RIESZ_PYR_MAX 7

// reach of the 9x9 kernels, width of the halo a level needs
inline constexpr int PYR_HALO = 4;

// ---------------------------------------------------------------------------------------------
// Kernels, tables row by row (row = i offset)
// ---------------------------------------------------------------------------------------------

inline constexpr double PYR_L9[81] = {
    -0.0001, -0.0007, -0.0023, -0.0046, -0.0057, -0.0046, -0.0023, -0.0007, -0.0001,
    -0.0007, -0.0030, -0.0047, -0.0025, -0.0003, -0.0025, -0.0047, -0.0030, -0.0007,
    -0.0023, -0.0047, 0.0054, 0.0272, 0.0387, 0.0272, 0.0054, -0.0047, -0.0023,
    -0.0046, -0.0025, 0.0272, 0.0706, 0.0910, 0.0706, 0.0272, -0.0025, -0.0046,
    -0.0057, -0.0003, 0.0387, 0.0910, 0.1138, 0.0910, 0.0387, -0.0003, -0.0057,
    -0.0046, -0.0025, 0.0272, 0.0706, 0.0910, 0.0706, 0.0272, -0.0025, -0.0046,
    -0.0023, -0.0047, 0.0054, 0.0272, 0.0387, 0.0272, 0.0054, -0.0047, -0.0023,
    -0.0007, -0.0030, -0.0047, -0.0025, -0.0003, -0.0025, -0.0047, -0.0030, -0.0007,
    -0.0001, -0.0007, -0.0023, -0.0046, -0.0057, -0.0046, -0.0023, -0.0007, -0.0001};
inline constexpr double PYR_H9[81] = {
    0.0000, 0.0003, 0.0011, 0.0022, 0.0027, 0.0022, 0.0011, 0.0003, 0.0000,
    0.0003, 0.0020, 0.0059, 0.0103, 0.0123, 0.0103, 0.0059, 0.0020, 0.0003,
    0.0011, 0.0059, 0.0151, 0.0249, 0.0292, 0.0249, 0.0151, 0.0059, 0.0011,
    0.0022, 0.0103, 0.0249, 0.0402, 0.0469, 0.0402, 0.0249, 0.0103, 0.0022,
    0.0027, 0.0123, 0.0292, 0.0469, -0.9455, 0.0469, 0.0292, 0.0123, 0.0027,
    0.0022, 0.0103, 0.0249, 0.0402, 0.0469, 0.0402, 0.0249, 0.0103, 0.0022,
    0.0011, 0.0059, 0.0151, 0.0249, 0.0292, 0.0249, 0.0151, 0.0059, 0.0011,
    0.0003, 0.0020, 0.0059, 0.0103, 0.0123, 0.0103, 0.0059, 0.0020, 0.0003,
    0.0000, 0.0003, 0.0011, 0.0022, 0.0027, 0.0022, 0.0011, 0.0003, 0.0000};

// Riesz 3x3. Signs follow the FFT convention R1 = -i*kx/|k|.
inline constexpr double PYR_FO1[9] = {
    0.12, 0.34, 0.12,
    0.00, 0.00, 0.00,
    -0.12, -0.34, -0.12};
inline constexpr double PYR_FO2[9] = {
    0.12, 0.00, -0.12,
    0.34, 0.00, -0.34,
    0.12, 0.00, -0.12};

// published 1-D filters, used when the run has one cell in y
inline constexpr double PYR_L9_1D[9] = {-0.0209, -0.0219, 0.0900, 0.2723, 0.3611, 0.2723, 0.0900, -0.0219, -0.0209};
inline constexpr double PYR_H9_1D[9] = {0.0099, 0.0492, 0.1230, 0.2020, -0.7633, 0.2020, 0.1230, 0.0492, 0.0099};
inline constexpr double PYR_FO1_1D[3] = {0.49, 0.00, -0.49};
inline constexpr double PYR_FO2_1D[3] = {0.00, 0.00, 0.00};

// kernel of radius r: (2r+1)x(2r+1) table and the 1-D filter of the runs with one cell in y
struct pyr_kernel
{
    const double *k2;
    const double *k1;
    int r;
};

inline constexpr pyr_kernel PYR_L9K = {PYR_L9, PYR_L9_1D, 4};
inline constexpr pyr_kernel PYR_H9K = {PYR_H9, PYR_H9_1D, 4};
inline constexpr pyr_kernel PYR_FO_K[2] = {{PYR_FO1, PYR_FO1_1D, 1}, {PYR_FO2, PYR_FO2_1D, 1}};

// ---------------------------------------------------------------------------------------------
// Level grids and fields
// ---------------------------------------------------------------------------------------------

inline constexpr int PYR_TILE = 0;    // the level on this rank, halo filled by exchange
inline constexpr int PYR_GLOBAL = 1;  // the level gathered on every rank, zero halo

// Grid of one level, named after the lexer. Cell (i,j) has the global level index
// gi = i + origin_i and sits on the fine global cell s*gi. knoy = 1, origin_j = 0 when !j2d.
struct pyr_grid
{
    int knox = 0, knoy = 0;
    int origin_i = 0, origin_j = 0;
    int s = 1;                  // 2^l
    double dx = 1.0, dy = 1.0;  // spacing of the level
    bool j2d = false;           // more than one cell in y
    bool tile = true;           // false: grid gathered on every rank
};

// Next level on the same rank. The rank owns the coarse cells gi whose fine cell s*gi lies in its
// fine tile [origin, origin+kno): gi from ceil(origin/s) to floor((origin+kno-1)/s).
inline pyr_grid pyr_grid_coarser_tile(const pyr_grid &g, int origin_i, int knox, int origin_j, int knoy)
{
    pyr_grid c = g;
    c.s = 2 * g.s;
    c.dx = 2.0 * g.dx;
    c.dy = 2.0 * g.dy;

    auto owned = [&](int origin, int kno, int &first, int &count)
    {
        first = (origin + c.s - 1) / c.s;
        const int last = (origin + kno - 1) / c.s;
        count = (last >= first) ? (last - first + 1) : 0;
    };

    owned(origin_i, knox, c.origin_i, c.knox);

    c.origin_j = 0;
    c.knoy = 1;
    if(g.j2d)
        owned(origin_j, knoy, c.origin_j, c.knoy);

    return c;
}

// next level of the gathered grid
inline pyr_grid pyr_grid_coarser_global(const pyr_grid &g)
{
    pyr_grid c = g;
    c.s = 2 * g.s;
    c.dx = 2.0 * g.dx;
    c.dy = 2.0 * g.dy;
    c.knox = (g.knox + 1) / 2;
    c.knoy = g.j2d ? (g.knoy + 1) / 2 : 1;
    return c;
}

// Field on a level grid, stored like slice_base: f(i,j), halo at negative indices,
// flat and contiguous so that the buffer can be exchanged by ghostcell::gcsl_thick_halo.
struct pyr_field
{
    pyr_grid g;
    int hw = 0;
    int jmax = 0;  // row stride, knoy + 2*hw
    std::vector<double> V;

    void allocate(const pyr_grid &grid, int halo)
    {
        g = grid;
        hw = halo;
        jmax = grid.knoy + 2 * halo;
        V.assign((size_t)(grid.knox + 2 * halo) * jmax, 0.0);
    }

    inline double &operator()(int i, int j) noexcept { return V[(size_t)((i + hw) * jmax + (j + hw))]; }
    inline double operator()(int i, int j) const noexcept { return V[(size_t)((i + hw) * jmax + (j + hw))]; }

    // row(i)[j] = f(i,j)
    inline const double *row(int i) const noexcept { return V.data() + (i + hw) * jmax + hw; }

    double *data() noexcept { return V.data(); }
};

#define PYRLOOP(F)                      \
    for(int i = 0; i < (F).g.knox; ++i) \
        for(int j = 0; j < (F).g.knoy; ++j)

// fields of one stage of the recursion, allocated once, only those the level uses
struct pyr_stage
{
    pyr_field low;                      // halo 4, low-pass field entering the level
    pyr_field band;                     // halo 1, H9*low, or copy of the residual
    pyr_field band_fo[2];               // halo 4, Riesz pair of the band
    pyr_field up_fo[2];                 // halo 4, result of the coarser stage, zero-upsampled
    pyr_field fo[2];                    // halo 1, Riesz pair rebuilt at this level
    pyr_field out_A, out_phase, out_k;  // halo 1, monogenic variables for P313
};

// ---------------------------------------------------------------------------------------------
// Convolution
// ---------------------------------------------------------------------------------------------

// Kernel K applied to f at cell (i,j). Requires: f halo of width K.r filled.
// even_only, for a zero-upsampled f: only the cells of even global index can be non-zero, the
// others hold 0.0 and are skipped. A skipped term is +-0.0 and the sum starts at +0.0, so it
// never becomes -0.0: the result is the full sum, bit for bit.
inline double pyr_conv(const pyr_field &f, int i, int j, const pyr_kernel &K, bool even_only = false)
{
    const int r = K.r;
    const int w = 2 * r + 1;
    const int gi = i + f.g.origin_i;
    const int gj = j + f.g.origin_j;
    double s = 0.0;

    for(int p = -r; p <= r; ++p)
    {
        if(even_only && ((gi + p) & 1))
            continue;
        const double *row = f.row(i + p);

        if(!f.g.j2d)
        {
            s += K.k1[p + r] * row[j];
            continue;
        }

        for(int q = -r; q <= r; ++q)
        {
            if(even_only && ((gj + q) & 1))
                continue;
            s += K.k2[(p + r) * w + (q + r)] * row[j + q];
        }
    }

    return s;
}

// ---------------------------------------------------------------------------------------------
// Operators of one level
// ---------------------------------------------------------------------------------------------

// dst = src on the interior
void pyr_copy(const pyr_field &src, pyr_field &dst);

// band = H9*low.  Requires: low halo filled.
void pyr_highpass(const pyr_field &low, pyr_field &band);

// fo[c] = FO_c*band.  Requires: band halo filled.
void pyr_riesz(const pyr_field &band, pyr_field fo[2]);

// coarse = L9*low taken at the cell 2*gc of low for every cell gc of coarse.  Requires: low halo filled.
void pyr_down(const pyr_field &low, pyr_field &coarse);

// fine = coarse on the even global cells of the finer grid, 0 elsewhere
void pyr_upsample(const pyr_field &coarse, pyr_field &fine);

// fo = gain*L9*up, plus H9*band_fo when the level has a band; gain = 4 in 3D, 2 in 2D.
// Requires: up and band_fo halos filled.
void pyr_rebuild(const pyr_field &up, const pyr_field &band_fo, bool has_band, pyr_field &fo);

// A = sqrt(band^2 + |fo|^2), phase = atan2(|fo|, band), no smoothing
void pyr_monogenic(const pyr_field &band, const pyr_field fo[2], pyr_field &A, pyr_field &phase);

// k = |grad phase| on the spacing of the level.  Requires: phase halo filled.
void pyr_wavenumber(const pyr_field &phase, pyr_field &k);

#endif
