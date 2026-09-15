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

// Halo exchange of width HW for a flat row-major buffer (nx+2HW)x(ny+2HW), used by the Riesz
// pyramid. buf[a*(ny+2HW)+b], interior cell (i,j) at (i+HW, j+HW). The halo is zero-filled by
// the caller and stays zero at the domain boundary. The corners are filled by exchanging in x
// over the interior, then in y over the extended x range. Face neighbours only.
//
// The buffer is horizontal: the z directions carry nothing. Sendrecv_3D takes the six directions
// and only uses the neighbours of the decomposition, whatever the number of decomposed directions.

#include "ghostcell.h"
#include "lexer.h"
#include <vector>

// collective MPI communication
void ghostcell::gcsl_thick_halo(lexer *p, double *buf, int HW, int nx, int ny)
{
    if(!do_comms)
        return;

    starttime = timer();

    const int W2 = ny + 2 * HW;
    double dummy = 0.0;

    // x direction, interior j: directions 0 (x-) and 1 (x+)
    {
        const int cnt = HW * ny;
        std::vector<double> sL(cnt, 0.0), sR(cnt, 0.0), rL(cnt, 0.0), rR(cnt, 0.0);

        for(int s = 0; s < HW; ++s)
            for(int j = 0; j < ny; ++j)
            {
                sL[s * ny + j] = buf[(HW + s) * W2 + (HW + j)];
                sR[s * ny + j] = buf[(nx + s) * W2 + (HW + j)];
            }

        const void *sp[6] = {cnt ? sL.data() : &dummy, cnt ? sR.data() : &dummy, &dummy, &dummy, &dummy, &dummy};
        void *rp[6] = {cnt ? rL.data() : &dummy, cnt ? rR.data() : &dummy, &dummy, &dummy, &dummy, &dummy};
        int sc[6] = {cnt, cnt, 0, 0, 0, 0};
        int rc[6] = {cnt, cnt, 0, 0, 0, 0};
        Sendrecv_3D(sp, sc, rp, rc, MPI_DOUBLE);

        for(int s = 0; s < HW; ++s)
            for(int j = 0; j < ny; ++j)
                buf[(nx + HW + s) * W2 + (HW + j)] = rR[s * ny + j];

        for(int s = 0; s < HW; ++s)
            for(int j = 0; j < ny; ++j)
                buf[s * W2 + (HW + j)] = rL[s * ny + j];
    }

    // y direction, extended i including the x halos: directions 2 (y-) and 3 (y+)
    {
        const int na = nx + 2 * HW;
        const int cnt = na * HW;
        std::vector<double> sB(cnt, 0.0), sT(cnt, 0.0), rB(cnt, 0.0), rT(cnt, 0.0);

        for(int a = 0; a < na; ++a)
            for(int s = 0; s < HW; ++s)
            {
                sB[a * HW + s] = buf[a * W2 + (HW + s)];
                sT[a * HW + s] = buf[a * W2 + (ny + s)];
            }

        const void *sp[6] = {&dummy, &dummy, cnt ? sB.data() : &dummy, cnt ? sT.data() : &dummy, &dummy, &dummy};
        void *rp[6] = {&dummy, &dummy, cnt ? rB.data() : &dummy, cnt ? rT.data() : &dummy, &dummy, &dummy};
        int sc[6] = {0, 0, cnt, cnt, 0, 0};
        int rc[6] = {0, 0, cnt, cnt, 0, 0};
        Sendrecv_3D(sp, sc, rp, rc, MPI_DOUBLE);

        for(int a = 0; a < na; ++a)
            for(int s = 0; s < HW; ++s)
                buf[a * W2 + (ny + HW + s)] = rT[a * HW + s];

        for(int a = 0; a < na; ++a)
            for(int s = 0; s < HW; ++s)
                buf[a * W2 + s] = rB[a * HW + s];
    }

    endtime = timer();
    p->xtime += endtime - starttime;
}
