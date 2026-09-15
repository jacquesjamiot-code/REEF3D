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

#include "fnpf_breaking_barthelemy_pyramid_core.h"

void pyr_copy(const pyr_field &src, pyr_field &dst)
{
    PYRLOOP(dst)
    dst(i, j) = src(i, j);
}

void pyr_highpass(const pyr_field &low, pyr_field &band)
{
    PYRLOOP(band)
    band(i, j) = pyr_conv(low, i, j, PYR_H9K);
}

void pyr_riesz(const pyr_field &band, pyr_field fo[2])
{
    PYRLOOP(fo[0])
    for(int c = 0; c < 2; ++c)
        fo[c](i, j) = pyr_conv(band, i, j, PYR_FO_K[c]);
}

void pyr_down(const pyr_field &low, pyr_field &coarse)
{
    const pyr_grid &g = low.g, &gc = coarse.g;

    PYRLOOP(coarse)
    {
        int ii = 2 * (i + gc.origin_i) - g.origin_i;
        if(ii < 0)
            ii = 0;
        if(ii > g.knox - 1)
            ii = g.knox - 1;
        int jj = 0;
        if(g.j2d)
        {
            jj = 2 * (j + gc.origin_j) - g.origin_j;
            if(jj < 0)
                jj = 0;
            if(jj > g.knoy - 1)
                jj = g.knoy - 1;
        }
        coarse(i, j) = pyr_conv(low, ii, jj, PYR_L9K);
    }
}

void pyr_upsample(const pyr_field &coarse, pyr_field &fine)
{
    const pyr_grid &g = fine.g, &gc = coarse.g;

    // the parent of an even cell gi is gi/2, if this rank owns it
    PYRLOOP(fine)
    {
        const int gi = i + g.origin_i;
        const int gj = g.j2d ? (j + g.origin_j) : 0;
        double v = 0.0;
        if((gi % 2 == 0) && (!g.j2d || gj % 2 == 0))
        {
            const int ic = gi / 2 - gc.origin_i;
            const int jc = g.j2d ? (gj / 2 - gc.origin_j) : 0;
            if(ic >= 0 && ic < gc.knox && jc >= 0 && jc < gc.knoy)
                v = coarse(ic, jc);
        }
        fine(i, j) = v;
    }
}

void pyr_rebuild(const pyr_field &up, const pyr_field &band_fo, bool has_band, pyr_field &fo)
{
    const double gain = up.g.j2d ? 4.0 : 2.0;

    PYRLOOP(fo)
    {
        fo(i, j) = gain * pyr_conv(up, i, j, PYR_L9K, true);
        if(has_band)
            fo(i, j) += pyr_conv(band_fo, i, j, PYR_H9K);
    }
}

void pyr_monogenic(const pyr_field &band, const pyr_field fo[2], pyr_field &A, pyr_field &phase)
{
    PYRLOOP(A)
    {
        const double h = band(i, j);
        const double rmag = sqrt(fo[0](i, j) * fo[0](i, j) + fo[1](i, j) * fo[1](i, j));
        A(i, j) = sqrt(h * h + rmag * rmag);
        phase(i, j) = atan2(rmag, h);
    }
}

void pyr_wavenumber(const pyr_field &phase, pyr_field &k)
{
    const pyr_grid &g = phase.g;

    PYRLOOP(k)
    {
        const double dphx = (phase(i + 1, j) - phase(i - 1, j)) / (2.0 * g.dx);
        double dphy = 0.0;
        if(g.j2d)
            dphy = (phase(i, j + 1) - phase(i, j - 1)) / (2.0 * g.dy);
        k(i, j) = sqrt(dphx * dphx + dphy * dphy);
    }
}
