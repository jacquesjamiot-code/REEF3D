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

#ifndef FNPF_PRINT_BREAKING_GAUGE_H_
#define FNPF_PRINT_BREAKING_GAUGE_H_

#include "increment.h"
#include <fstream>
#include <iostream>

class lexer;
class fdm_fnpf;
class ghostcell;
class slice4;

using namespace std;

// breaking criterion B, wavenumber k, celerity c and omega at the P51 gauges, P314
class fnpf_print_breaking_gauge : public increment
{
public:
    fnpf_print_breaking_gauge(lexer *, fdm_fnpf *);
    virtual ~fnpf_print_breaking_gauge();

    void start(lexer *, fdm_fnpf *, ghostcell *);

private:
    void ini_location(lexer *);
    void open_file(lexer *, ofstream &, const char *);
    void write_slice(lexer *, ghostcell *, ofstream &, slice4 *);

    double *x, *y;
    int gauge_num;

    int *iloc, *jloc, *flag;
    int n;
    ofstream out_B, out_k, out_c, out_omega;
    const int fileFlushMaxCount;
};

#endif
