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

#include"fnpf_vtp_fsf.h"
#include"lexer.h"

void fnpf_vtp_fsf::pvtp(lexer *p, int num)
{
    sprintf(name,"./REEF3D_FNPF_VTP_FSF/REEF3D-FNPF-FSF-%08i.pvtp",num);

    ofstream result;
    result.open(name);

    vtp3D::beginningParallel(p,result);

    vtp3D::pointsParallel(result);

    result<<"<PPointData>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"Fifsf\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"eta\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"depth\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"breaking\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"coastline\"/>\n";
    if(p->P23==1)
        result<<"<PDataArray type=\"Float32\" Name=\"test\"/>\n";
    if(p->P110==1)
        result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>\n";

    // riesz pyramid levels
    if(p->P313 > 0)
    {
        const int pyr_num = p->A385_nmax - p->A385_nmin + 2;
        const char *pyr_name[7] = {"A_pyramid", "phase_pyramid", "k_pyramid", "theta_pyramid",
                                   "band_pyramid", "fo1_pyramid", "fo2_pyramid"};

        for(int q = 0; q < 7; ++q)
        {
            result << "<PDataArray type=\"Float32\" Name=\"" << pyr_name[q] << "\" NumberOfComponents=\"" << pyr_num << "\"";
            for(int m = 0; m < pyr_num; ++m)
            {
                if(m < pyr_num - 1)
                    result << " ComponentName" << m << "=\"level" << p->A385_nmin + m << "\"";
                else
                    result << " ComponentName" << m << "=\"residual\"";
            }
            result << "/>\n";
        }
    }
    result<<"</PPointData>\n";

    char pname[200];
    for(n=0; n<p->M10; ++n)
    {
        sprintf(pname,"REEF3D-FNPF-FSF-%08i-%06i.vtp",num,n+1);
        result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    vtp3D::endingParallel(result);

    result.close();
}
