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

// Breaking zone and dissipation.
// Y. Wang, G. Ducrozet, "Modeling of wave breaking in short-crested seas",
// 40th International Workshop on Water Waves and Floating Bodies, 2025, hal-05533115.
//
// Zone: a node with B >= B_on is the origin of a rectangle oriented along the propagation
// direction theta, of length L^t_br = 2*pi/k and plateau width L^n_br = min(dx,dy), with a ramp
// on every side. Every covered node breaks for Tbr = 2*pi/omega; where zones overlap the
// breaker with the largest nu_eddy*R wins. The loops over a zone are bounded by its four edges
// in cell indices; the ramp weight R > 0 decides.
//
// Parallel: every rank stamps the zones on its own cells. Three exchanges through
// ghostcell::deliver_quads, addressed by the corners of a zone or by its onset cell:
//     1  breakers   own zones to the ranks they touch, nu_eddy = 1
//     2  intensity  part of I on every touched rank, back to the owner
//     3  intensity  own zones with nu_eddy to the ranks they touch
//
// Dissipation (Eq. 1-4): a viscous term 2*nu_eddy*R*lap_h is added to both free-surface
// conditions, solved by fnpf_fsfbc::damping() with
//     c->vb   = 2*nu_eddy*R
//     nu_eddy = -b*C^5*l_cl / (2*g*I),   I = int_S R*lap_h(phi_s)*d(eta)/dt dS
// with b = A382 and l_cl = L^n_br. With b = 0 the model detects only and vb = 0.
//
// Order of the calls: breakers, freeze_eta_t, intensity, zone, apply. The zone needs nu_eddy.

#ifndef FNPF_BREAKING_WANG_H_
#define FNPF_BREAKING_WANG_H_

#include "fnpf_breaking_barthelemy.h"
#include "increment.h"
#include "slice4.h"
#include "sliceint4.h"
#include <fstream>
#include <vector>

class lexer;
class fdm_fnpf;
class ghostcell;
class slice;

using namespace std;

static const double WANG_ZONE_WIDTH = 1.0;    // normal plateau, in units of min(dx,dy)
static const double WANG_ZONE_RAMP_N = 0.25;  // normal ramp on each side, same units
static const double WANG_ZONE_RAMP_T = 0.25;  // tangential ramp, fraction of Lbr
static const double WANG_TBR_MAX_TP = 3.0;    // maximum breaking duration, in peak periods
static const double WANG_VB_MAX = 10.0;       // ceiling on vb

// rectangle of one breaker: s along theta, t across, edges s0, s1, -t1, +t1
struct wang_zone_geom
{
    int gi0, gj0;   // origin, global cell
    double ct, st;  // cos and sin of theta
    double s0, s1, t1;
    double Lrp_t, Lrp_n;  // tangential and normal ramps
    double dx0, dy0;      // spacings of the owner, used by every rank for this zone
    double ci[4], cj[4];  // corners (s0,-t1) (s1,-t1) (s1,+t1) (s0,+t1) in global cell indices, counterclockwise
};

// cell of a zone inside this subdomain
struct wang_cell
{
    int i, j;  // local indices
    double R;  // ramp weight, > 0
};

// zone of another rank that touches this subdomain
struct wang_foreign
{
    int id;  // index of the breaker in the list of its owner
    wang_zone_geom geom;
    double Tbr, nu;
};

// zone record exchanged between ranks, WANG_REC doubles:
// ci[0..3], cj[0..3] (address), id, gi0, gj0, ct, st, s0, s1, t1, Lrp_t, Lrp_n, dx0, dy0, Tbr, nu
static const int WANG_REC = 22;

// contribution record, WANG_CONTRIB doubles: gi0 x4, gj0 x4 (address: the onset cell), id, part of I
static const int WANG_CONTRIB = 10;

inline void wang_pack(int id, const wang_zone_geom &g, double Tbr, double nu, std::vector<double> &buf)
{
    const double r[WANG_REC] = {g.ci[0], g.ci[1], g.ci[2], g.ci[3], g.cj[0], g.cj[1], g.cj[2], g.cj[3],
                                (double)id, (double)g.gi0, (double)g.gj0,
                                g.ct, g.st, g.s0, g.s1, g.t1, g.Lrp_t, g.Lrp_n, g.dx0, g.dy0, Tbr, nu};
    buf.insert(buf.end(), r, r + WANG_REC);
}

inline wang_foreign wang_unpack(const double *r)
{
    wang_foreign z;
    for(int k = 0; k < 4; ++k)
    {
        z.geom.ci[k] = r[k];
        z.geom.cj[k] = r[4 + k];
    }
    z.id = (int)r[8];
    z.geom.gi0 = (int)r[9];
    z.geom.gj0 = (int)r[10];
    z.geom.ct = r[11];
    z.geom.st = r[12];
    z.geom.s0 = r[13];
    z.geom.s1 = r[14];
    z.geom.t1 = r[15];
    z.geom.Lrp_t = r[16];
    z.geom.Lrp_n = r[17];
    z.geom.dx0 = r[18];
    z.geom.dy0 = r[19];
    z.Tbr = r[20];
    z.nu = r[21];
    return z;
}

struct wang_breaker
{
    bart_seed seed;
    double Lbr, Tbr;
    wang_zone_geom geom;
    double nu, I;  // nu_eddy and its integral
    int onset;     // 1 if the cell did not seed at the previous timestep
};

// data of one call of the breaking routine
struct wang_call
{
    int n_onset_glob = 0;               // breakers on all ranks
    std::vector<wang_breaker> brk;      // breakers of this rank
    std::vector<wang_foreign> foreign;  // zones of other ranks touching this subdomain
    std::vector<double> stamp;          // (nu_eddy*R, Tbr, R) per cell of this subdomain, empty if none
};

class fnpf_breaking_wang : public increment
{
public:
    fnpf_breaking_wang(lexer *, fdm_fnpf *, ghostcell *);
    virtual ~fnpf_breaking_wang();

    // collective MPI communication
    void breakers(lexer *, ghostcell *, const bart_config &, const std::vector<bart_seed> &,
                  int theta_refresh, wang_call &);

    void freeze_eta_t(lexer *, slice &eta, slice &eta_n, double alpha, int theta_refresh);

    // collective MPI communication
    void intensity(lexer *, fdm_fnpf *, ghostcell *, const bart_config &, slice &Fifsf, int theta_refresh,
                   wang_call &);

    // collective MPI communication
    void zone(lexer *, fdm_fnpf *, ghostcell *, const bart_config &, wang_call &);

    void apply(lexer *, fdm_fnpf *, ghostcell *, const bart_config &, const wang_call &);

private:
    wang_zone_geom zone_geometry(lexer *, int gi0, int gj0, double theta, double Lbr,
                                 const bart_config &);
    double zone_ramp_weight(lexer *, const wang_zone_geom &, int gi, int gj);
    void zone_cells(lexer *, const wang_zone_geom &, std::vector<wang_cell> &);
    void zone_stamp(lexer *, const wang_zone_geom &, double nu, double Tbr, double *stamp);
    double zone_sum(lexer *, const wang_zone_geom &);
    double eddy_viscosity(lexer *, double I, double C, const bart_config &);

    slice4 T_breaking;   // remaining breaking time
    slice4 t_break_end;  // end of breaking, absolute time
    slice4 R_break;      // ramp weight
    slice4 nu_break;     // nu_eddy*R
    slice4 eta_t_wang;   // d(eta)/dt at the first substage of the timestep
    slice4 G_wang;       // lap_h(phi_s)*d(eta)/dt, integrand of I, on this subdomain
    sliceint4 seed_old;  // 1 where a breaker was seeded at the previous timestep

    ofstream onset_out;  // onset log, P312
};

#endif
