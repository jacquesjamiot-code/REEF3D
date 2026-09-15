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

#include "fnpf_print_breaking_gauge.h"
#include "fdm_fnpf.h"
#include "ghostcell.h"
#include "lexer.h"
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

fnpf_print_breaking_gauge::fnpf_print_breaking_gauge(lexer *p, fdm_fnpf *c) : fileFlushMaxCount(100)
{
    gauge_num = p->P51;
    x = p->P51_x;
    y = p->P51_y;

    // Create Folder
    if(p->mpirank == 0)
        mkdir("./REEF3D_FNPF_WSF", 0777);

    if(p->mpirank == 0 && gauge_num > 0)
    {
        open_file(p, out_B, "REEF3D-FNPF-Breaking-B.dat");
        open_file(p, out_k, "REEF3D-FNPF-Breaking-k.dat");
        open_file(p, out_c, "REEF3D-FNPF-Breaking-c.dat");
        open_file(p, out_omega, "REEF3D-FNPF-Breaking-omega.dat");
    }

    p->Iarray(iloc, gauge_num);
    p->Iarray(jloc, gauge_num);
    p->Iarray(flag, gauge_num);

    ini_location(p);
}

fnpf_print_breaking_gauge::~fnpf_print_breaking_gauge()
{
    out_B.close();
    out_k.close();
    out_c.close();
    out_omega.close();
}

void fnpf_print_breaking_gauge::open_file(lexer *p, ofstream &out, const char *name)
{
    std::string path = std::string("./REEF3D_FNPF_WSF/") + name;
    out.open(path.c_str());

    out << "number of gauges:  " << gauge_num << endl
        << endl;
    out << "x_coord     y_coord" << endl;
    for(n = 0; n < gauge_num; ++n)
        out << n + 1 << "\t " << p->Xout(x[n], y[n]) << "\t " << p->Yout(x[n], y[n]) << endl;

    out << endl
        << endl;

    out << "time";
    for(n = 0; n < gauge_num; ++n)
        out << "\t P" << n + 1;

    out << endl;
}

void fnpf_print_breaking_gauge::start(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    write_slice(p, pgc, out_B, c->brk_B);
    write_slice(p, pgc, out_k, c->brk_k);
    write_slice(p, pgc, out_c, c->brk_c);
    write_slice(p, pgc, out_omega, c->brk_omega);
}

void fnpf_print_breaking_gauge::write_slice(lexer *p, ghostcell *pgc, ofstream &out, slice4 *f)
{
    double val;

    if(p->mpirank == 0)
        out << setprecision(9) << p->simtime << "\t";

    for(n = 0; n < gauge_num; ++n)
    {
        val = -1.0e20;

        if(flag[n] > 0 && f != nullptr)
            val = (*f)(iloc[n], jloc[n]);

        val = pgc->globalmax(val);

        if(p->mpirank == 0)
        {
            out << setprecision(9) << val << "\t";
            // flush print to disc limited to prevent data loss for many gauges
            if(n % fileFlushMaxCount == 0 && n != 0)
                out << std::flush;
        }
    }

    if(p->mpirank == 0)
        out << endl;
}

void fnpf_print_breaking_gauge::ini_location(lexer *p)
{
    for(n = 0; n < gauge_num; ++n)
    {
        iloc[n] = p->posc_i(x[n]);

        if(p->j_dir == 0)
            jloc[n] = 0;

        if(p->j_dir == 1)
            jloc[n] = p->posc_j(y[n]);

        if(iloc[n] >= 0 && iloc[n] < p->knox)
            if(jloc[n] >= 0 && jloc[n] < p->knoy)
                flag[n] = 1;
    }
}
