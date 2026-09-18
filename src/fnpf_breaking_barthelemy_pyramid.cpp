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

// Riesz transform of eta rebuilt from the Riesz pyramid, A380 3.
//
// One recursive stage per level l, on the tile of the rank (PYR_TILE) or on the grid gathered
// on every rank (PYR_GLOBAL):
//   analysis   low -> band H9*low and its Riesz pair (levels N_min..N_max, from A385), low -> L9 decimated
//   recursion  the coarser stage returns its rebuilt pair fo_{l+1}
//   synthesis  fo_l = gain * L9(zero-upsampled fo_{l+1}) + H9(Riesz pair of the band)
// The base case is the low-pass residual N_max+1, whose Riesz pair seeds the synthesis.
// The coarse cell of global index gi at level l sits on the fine global cell 2^l*gi, so every
// rank owns a contiguous block of coarse cells and exchanges a 4-cell halo with its neighbours.
// When the tiles of a level become narrower than the halo, the level is gathered and the
// remaining levels run in serial, then the result is cut back to the tile.
// With P313 the monogenic variables of every level and of the residual are kept for output,
// constant on each coarse cell. Uniform spacing is assumed.
//
// Every halo exchange and collective of the pyramid is called from this file, through ghostcell.

#include "fdm_fnpf.h"
#include "fnpf_breaking_barthelemy.h"
#include "ghostcell.h"
#include "lexer.h"
#include <algorithm>
#include <vector>

// Collective on a tile field: halo exchange. Nothing on a gathered field, whose halo stays zero.
static void pyr_exchange(lexer *p, ghostcell *pgc, pyr_field &f)
{
    if(f.g.tile)
        pgc->gcsl_thick_halo(p, f.data(), f.hw, f.g.knox, f.g.knoy);
}

// Collective: tile field -> gathered field, through the flat buffer buf of the gathered size
static void pyr_gather(ghostcell *pgc, const pyr_field &tile, pyr_field &global, std::vector<double> &buf)
{
    const int gny = global.g.knoy;

    std::fill(buf.begin(), buf.end(), 0.0);
    PYRLOOP(tile)
    buf[(size_t)(i + tile.g.origin_i) * gny + (j + tile.g.origin_j)] = tile(i, j);

    pgc->globalsumV(buf.data(), (int)buf.size());

    PYRLOOP(global)
    global(i, j) = buf[(size_t)i * gny + j];
}

// rebuilt pair of the gathered stage -> tile of the rank
static void pyr_scatter(const pyr_stage &global, pyr_stage &tile)
{
    for(int c = 0; c < 2; ++c)
    {
        pyr_field &fo = tile.fo[c];
        PYRLOOP(fo)
        fo(i, j) = global.fo[c](i + fo.g.origin_i, j + fo.g.origin_j);
    }
}

// fields the stage of level l uses on grid g
static void allocate_stage(pyr_stage &S, const pyr_grid &g, bool has_band, bool residual, bool output)
{
    S.low.allocate(g, PYR_HALO);

    for(int c = 0; c < 2; ++c)
    {
        S.fo[c].allocate(g, 1);
        if(has_band)
            S.band_fo[c].allocate(g, PYR_HALO);
        if(!residual)
            S.up_fo[c].allocate(g, PYR_HALO);
    }

    if(has_band || residual)
        S.band.allocate(g, 1);

    if(output && (has_band || residual))
    {
        S.out_A.allocate(g, 1);
        S.out_phase.allocate(g, 1);
        S.out_k.allocate(g, 1);
    }
}

void fnpf_breaking_barthelemy::riesz_pyramid(lexer *p, fdm_fnpf *c, ghostcell *pgc, const bart_config &cfg)
{
    if(!pyr_setup_done)
        pyr_setup(p, cfg);

    pyr_stage &S0 = pyr_stages[PYR_TILE][0];

    // eta, zero on dry cells
    SLICELOOP4
    S0.low(i, j) = bart_wet(p, c, i, j) ? c->eta(i, j) : 0.0;

    riesz_level(p, pgc, 0, PYR_TILE);

    slice4 &F1 = *fo1, &F2 = *fo2;
    SLICELOOP4
    {
        F1(i, j) = S0.fo[0](i, j);
        F2(i, j) = S0.fo[1](i, j);
    }
    pgc->gcsl_start4(p, F1, 1);
    pgc->gcsl_start4(p, F2, 1);
}

// grids of all levels and their fields, once: the decomposition does not change
void fnpf_breaking_barthelemy::pyr_setup(lexer *p, const bart_config &cfg)
{
    const int l_residual = pyr_lev_max + 1;
    const bool output = (pyr_num > 0);

    auto allocate = [&](int mode, int l, const pyr_grid &g)
    {
        allocate_stage(pyr_stages[mode][l], g, (l >= pyr_lev_min && l <= pyr_lev_max),
                       (l == l_residual), output);
    };

    pyr_grid tile_grid;
    tile_grid.knox = p->knox;
    tile_grid.knoy = p->knoy;
    tile_grid.j2d = (p->j_dir == 1);
    tile_grid.origin_i = p->origin_i;
    tile_grid.origin_j = tile_grid.j2d ? p->origin_j : 0;
    tile_grid.dx = cfg.dx0;
    tile_grid.dy = cfg.dy0;

    pyr_gather_level = -1;

    for(int l = 0; l <= l_residual; ++l)
    {
        allocate(PYR_TILE, l, tile_grid);

        // global test, identical on every rank; the residual is never gathered
        const int s = tile_grid.s;
        const int gnx = (p->gknox + s - 1) / s;
        const int gny = tile_grid.j2d ? (p->gknoy + s - 1) / s : 1;

        if(l < l_residual && (gnx < PYR_HALO * p->mx || (tile_grid.j2d && gny < PYR_HALO * p->my)))
        {
            pyr_gather_level = l;
            pyr_gather_buf.assign((size_t)gnx * gny, 0.0);

            pyr_grid global_grid = tile_grid;
            global_grid.knox = gnx;
            global_grid.knoy = gny;
            global_grid.origin_i = 0;
            global_grid.origin_j = 0;
            global_grid.tile = false;

            for(int l2 = l; l2 <= l_residual; ++l2)
            {
                allocate(PYR_GLOBAL, l2, global_grid);
                global_grid = pyr_grid_coarser_global(global_grid);
            }
            break;
        }

        tile_grid = pyr_grid_coarser_tile(tile_grid, p->origin_i, p->knox, p->origin_j, p->knoy);
    }

    pyr_setup_done = 1;
}

// One stage of the pyramid on level l of mode PYR_TILE or PYR_GLOBAL; fills its rebuilt pair fo.
// Requires: the interior of pyr_stages[mode][l].low.
void fnpf_breaking_barthelemy::riesz_level(lexer *p, ghostcell *pgc, int l, int mode)
{
    pyr_stage &S = pyr_stages[mode][l];

    // base case: Riesz pair of the low-pass residual
    if(l == pyr_lev_max + 1)
    {
        pyr_copy(S.low, S.band);
        pyr_exchange(p, pgc, S.band);
        pyr_riesz(S.band, S.fo);

        if(pyr_num > 0)
            pyr_export(p, pgc, l, S, S.fo);
        return;
    }

    // tiles narrower than the halo: the remaining levels on the gathered grid
    if(mode == PYR_TILE && l == pyr_gather_level)
    {
        pyr_stage &G = pyr_stages[PYR_GLOBAL][l];

        pyr_gather(pgc, S.low, G.low, pyr_gather_buf);
        riesz_level(p, pgc, l, PYR_GLOBAL);
        pyr_scatter(G, S);
        return;
    }

    const bool has_band = (l >= pyr_lev_min);
    pyr_stage &Sc = pyr_stages[mode][l + 1];

    // analysis
    pyr_exchange(p, pgc, S.low);

    if(has_band)
    {
        pyr_highpass(S.low, S.band);
        pyr_exchange(p, pgc, S.band);
        pyr_riesz(S.band, S.band_fo);

        if(pyr_num > 0)
            pyr_export(p, pgc, l, S, S.band_fo);
    }

    pyr_down(S.low, Sc.low);

    // coarser levels
    riesz_level(p, pgc, l + 1, mode);

    // synthesis
    for(int c = 0; c < 2; ++c)
    {
        pyr_upsample(Sc.fo[c], S.up_fo[c]);
        pyr_exchange(p, pgc, S.up_fo[c]);
    }

    if(has_band)
        for(int c = 0; c < 2; ++c)
            pyr_exchange(p, pgc, S.band_fo[c]);

    for(int c = 0; c < 2; ++c)
        pyr_rebuild(S.up_fo[c], S.band_fo[c], has_band, S.fo[c]);
}

// Monogenic variables of level l on its own grid, replicated on the fine grid of the rank.
// fo is the Riesz pair of the band, or the rebuilt pair of the residual.
// Requires: S.band halo filled.
void fnpf_breaking_barthelemy::pyr_export(lexer *p, ghostcell *pgc, int l, pyr_stage &S, pyr_field fo[2])
{
    const pyr_grid &g = S.band.g;

    pyr_monogenic(S.band, fo, S.out_A, S.out_phase);
    pyr_exchange(p, pgc, S.out_phase);
    pyr_wavenumber(S.out_phase, S.out_k);

    // a fine cell may belong to a coarse cell of the neighbour: halos of the replicated fields
    pyr_exchange(p, pgc, S.out_k);
    pyr_exchange(p, pgc, S.out_A);
    pyr_exchange(p, pgc, fo[0]);
    pyr_exchange(p, pgc, fo[1]);

    if(g.knox * g.knoy == 0)
        return;

    bart_pyr_level &LV = *pyr_level[l];
    SLICELOOP4
    {
        const int gi = i + p->origin_i, gj = j + p->origin_j;
        int ic = gi / g.s - g.origin_i;
        if(ic < -1)
            ic = -1;
        if(ic > g.knox)
            ic = g.knox;
        int jc = 0;
        if(g.j2d)
        {
            jc = gj / g.s - g.origin_j;
            if(jc < -1)
                jc = -1;
            if(jc > g.knoy)
                jc = g.knoy;
        }
        LV.k(i, j) = S.out_k(ic, jc);
        LV.A(i, j) = S.out_A(ic, jc);
        LV.phase(i, j) = S.out_phase(ic, jc);
        LV.fo1(i, j) = fo[0](ic, jc);
        LV.fo2(i, j) = fo[1](ic, jc);
        LV.band(i, j) = S.band(ic, jc);
    }
}
