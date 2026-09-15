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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include<algorithm>
#include<vector>

void ghostcell::gather_int(int *sendbuf, int sendcount, int *recvbuf, int recvcount)
{
    MPI_Gather(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, 0, mpi_comm);
}

void ghostcell::gather_double(double *sendbuf, int sendcount, double *recvbuf, int recvcount)
{
    MPI_Gather(sendbuf,sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, 0, mpi_comm);
}

void ghostcell::gatherv_int(int *sendbuf, int sendcount, int *recvbuf, int *recvcount, int *recvdispl)
{
    MPI_Gatherv(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, recvdispl, MPI_INT, 0, mpi_comm);
}

void ghostcell::gatherv_double(double *sendbuf, int sendcount, double *recvbuf, int *recvcount, int *recvdispl)
{
    MPI_Gatherv(sendbuf,sendcount, MPI_DOUBLE, recvbuf, recvcount, recvdispl, MPI_DOUBLE, 0, mpi_comm);
}

void ghostcell::allgather_int(int *sendbuf, int sendcount, int *recvbuf, int recvcount)
{
    MPI_Allgather(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, MPI_INT, mpi_comm);
}

void ghostcell::allgatherv_int(int *sendbuf, int sendcount, int *recvbuf, int *recvcount, int *recvdispl)
{
    MPI_Allgatherv(sendbuf,sendcount, MPI_INT, recvbuf, recvcount, recvdispl, MPI_INT, mpi_comm);
}

void ghostcell::bcast_int(int *sendbuf, int sendcount)
{
    MPI_Bcast(sendbuf,sendcount, MPI_INT, 0, mpi_comm);
}

void ghostcell::bcast_double(double *sendbuf, int sendcount, int root)
{
    MPI_Bcast(sendbuf,sendcount, MPI_DOUBLE, root, mpi_comm);
}

// the convex quadrilateral q (i of the 4 vertices, then j, in contour order) meets the area of the
// cells of the subdomain t (i_lo, i_hi, j_lo, j_hi): separating axis test on the grid axes and on
// the normals of the edges of q. A point (4 equal vertices) is tested against the area only.
static bool quad_meets(const double *q, const int *t)
{
    const double *qi = q, *qj = q + 4;
    const double ti[2] = {t[0] - 0.5, t[1] + 0.5};
    const double tj[2] = {t[2] - 0.5, t[3] + 0.5};

    // grid axes
    if(*std::max_element(qi, qi + 4) < ti[0] || *std::min_element(qi, qi + 4) > ti[1])
        return false;
    if(*std::max_element(qj, qj + 4) < tj[0] || *std::min_element(qj, qj + 4) > tj[1])
        return false;

    // normals of the edges of q
    for(int k = 0; k < 4; ++k)
    {
        const int m = (k + 1) % 4;
        const double ni = -(qj[m] - qj[k]), nj = qi[m] - qi[k];
        if(ni == 0.0 && nj == 0.0)
            continue;

        double qmin = ni * qi[0] + nj * qj[0], qmax = qmin;
        for(int v = 1; v < 4; ++v)
        {
            const double d = ni * qi[v] + nj * qj[v];
            qmin = std::min(qmin, d);
            qmax = std::max(qmax, d);
        }

        double tmin = ni * ti[0] + nj * tj[0], tmax = tmin;
        for(int a = 0; a < 2; ++a)
            for(int b = 0; b < 2; ++b)
            {
                const double d = ni * ti[a] + nj * tj[b];
                tmin = std::min(tmin, d);
                tmax = std::max(tmax, d);
            }

        if(qmax < tmin || qmin > tmax)
            return false;
    }

    return true;
}

// collective MPI communication: each record, len doubles starting with its address (a convex
// quadrilateral in global cell indices), is sent directly to every other rank whose subdomain meets
// the address. recv: records received from the other ranks, grouped by rank, in the order of the ranks.
void ghostcell::deliver_quads(lexer *p, const std::vector<double> &mine, int len, std::vector<double> &recv)
{
    const int P = p->mpi_size;

    // bounds of every subdomain, once; every rank reaches its first call at the same time
    if(tiles.empty())
    {
        int b[4] = {p->origin_i, p->origin_i + p->knox - 1, p->origin_j, p->origin_j + p->knoy - 1};
        tiles.resize(4 * P);
        MPI_Allgather(b, 4, MPI_INT, tiles.data(), 4, MPI_INT, mpi_comm);
    }

    // records for each rank, in their order
    std::vector<std::vector<double>> out(P);

    for(size_t k = 0; k < mine.size(); k += len)
        for(int r = 0; r < P; ++r)
            if(r != p->mpirank && quad_meets(&mine[k], &tiles[4 * r]))
                out[r].insert(out[r].end(), mine.begin() + k, mine.begin() + k + len);

    // one flat send buffer, with the count and the start of the part of every rank
    std::vector<int> scount(P, 0), sdispl(P, 0), rcount(P, 0), rdispl(P, 0);
    std::vector<double> sbuf;
    for(int r = 0; r < P; ++r)
    {
        sdispl[r] = (int)sbuf.size();
        scount[r] = (int)out[r].size();
        sbuf.insert(sbuf.end(), out[r].begin(), out[r].end());
    }

    // the counts first, so that every rank can size its receive buffer
    MPI_Alltoall(scount.data(), 1, MPI_INT, rcount.data(), 1, MPI_INT, mpi_comm);

    int total = 0;
    for(int r = 0; r < P; ++r)
    {
        rdispl[r] = total;
        total += rcount[r];
    }
    recv.assign(total, 0.0);

    MPI_Alltoallv(sbuf.data(), scount.data(), sdispl.data(), MPI_DOUBLE,
                  recv.data(), rcount.data(), rdispl.data(), MPI_DOUBLE, mpi_comm);
}
