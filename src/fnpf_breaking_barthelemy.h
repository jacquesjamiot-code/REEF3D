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

#ifndef FNPF_BREAKING_BARTHELEMY_H_
#define FNPF_BREAKING_BARTHELEMY_H_

#include "fdm_fnpf.h"
#include "fnpf_breaking_barthelemy_pyramid_core.h"
#include "increment.h"
#include "lexer.h"
#include "slice4.h"
#include "sliceint4.h"
#include <cmath>
#include <vector>

class ghostcell;

using namespace std;

static const double BART_KH_DEEP = 10.0;      // kh above which tanh(kh) is taken as 1
static const double BART_OMEGA_EPS = 1.0e-6;  // smallest omega estimate taken as a wave
static const double BART_C_MIN = 0.1;         // celerity floor [m/s]
static const double BART_C_MAX = 50.0;        // celerity ceiling [m/s]
static const double BART_B_MAX = 10.0;        // ceiling on B
static const double BART_OMEGA_LO = 0.5;      // omega floor, fraction of the dispersion relation
static const double BART_OMEGA_HI = 1.3;      // omega ceiling, fraction of the dispersion relation
static const double BART_DIR_FILTER = 0.2;    // weight of the filter on the phase increment
static const int BART_LAG_EVAL = 3;           // Lagrange derivative taken at the newest node
static const double BART_LAG_DTMIN = 1.0e-3;  // smallest admissible min/max ratio of the time intervals
static const double BART_HMIN_WD = 1.5;       // minimum water depth of the model, in wetting criteria A344

// centered derivatives on a non-uniform grid
static inline double centered_dx(slice4 &f, int i, int j, lexer *p, int mg)
{
    return ((f(i + 1, j) - f(i, j)) / p->DXP[i + mg] + (f(i, j) - f(i - 1, j)) / p->DXP[i - 1 + mg]) * 0.5;
}

static inline double centered_dy(slice4 &f, int i, int j, lexer *p, int mg)
{
    return ((f(i, j + 1) - f(i, j)) / p->DYP[j + mg] + (f(i, j) - f(i, j - 1)) / p->DYP[j - 1 + mg]) * 0.5;
}

// cell seen by the breaking model: under wetting-drying, wet and deeper than BART_HMIN_WD times the
// wetting criterion, which leaves out the film of A343 1; every cell without wetting-drying
static inline bool bart_wet(lexer *p, fdm_fnpf *c, int i, int j)
{
    return p->A343 == 0 || (p->wet[IJ] == 1 && c->WL(i, j) > BART_HMIN_WD * c->wd_criterion);
}

// linear dispersion relation
static inline double omega_dispersion(double k, double h)
{
    const double kh = k * h;
    return (kh > BART_KH_DEEP) ? sqrt(9.81 * k) : sqrt(9.81 * k * tanh(kh));
}

// configuration of one call
struct bart_config
{
    double B_on, b_str;             // A381, A382
    bool diss_on;                   // A382 > 0
    double dx0, dy0;                // representative spacings
    double width, Lrp_n;            // zone plateau and normal ramp
    int Nx, Ny, N;                  // global grid
    double m_xs, m_xe, m_ys, m_ye;  // seeding margins at the domain edges, A386
};

// breaking onset detected at one cell
struct bart_seed
{
    int i, j, gi, gj;           // local and global cell
    double k, omega, theta, C;  // wavenumber, angular frequency, direction, celerity
    double B, eta;              // criterion and elevation
};

// local variables and criterion of the selected transform
struct bart_chain
{
    bart_chain(lexer *);

    slice4 A, phase, k, theta;
    slice4 *phase_sig;  // signed phase: phase for Hilbert, phase_sig_own for Riesz
    slice4 phase_sig_own;
    slice4 phase_sig_old;  // signed phase at the previous timestep
    slice4 dphi_filt;      // filtered phase increment, its sign gives the direction
    slice4 d1, d2, d3;     // last three phase increments
    slice4 omega, c, B;
    sliceint4 wet_old;     // bart_wet at the previous timestep
};

// one level of the Riesz pyramid, P313
struct bart_pyr_level
{
    bart_pyr_level(lexer *);

    slice4 k, fo1, fo2, A, phase, band, theta;
    slice4 phase_sig, phase_sig_old, dphi_filt;
};

void bart_check_parameters(lexer *, ghostcell *);

class fnpf_breaking_barthelemy : public increment
{
public:
    fnpf_breaking_barthelemy(lexer *, fdm_fnpf *, ghostcell *);
    virtual ~fnpf_breaking_barthelemy();

    void step_begin(lexer *);
    bart_config config(lexer *);

    // collective MPI communication
    void diagnostic(lexer *, fdm_fnpf *, ghostcell *, const bart_config &, std::vector<bart_seed> &);

    int theta_refresh;  // 1 at the first call of a timestep

private:
    // collective MPI communication
    void spatial_transforms(lexer *, fdm_fnpf *, ghostcell *);
    // collective MPI communication
    void precompute_window(lexer *, ghostcell *);
    // collective MPI communication
    void riesz_pyramid(lexer *, fdm_fnpf *, ghostcell *, const bart_config &);
    void pyr_setup(lexer *, const bart_config &);
    // collective MPI communication
    void riesz_level(lexer *, ghostcell *, int l, int mode);
    // collective MPI communication
    void pyr_export(lexer *, ghostcell *, int l, pyr_stage &, pyr_field fo[2]);

    void local_hilbert(lexer *, fdm_fnpf *, ghostcell *);
    void local_riesz(lexer *, fdm_fnpf *, ghostcell *, slice4 &fo1_src, slice4 &fo2_src);
    void phase_increment(lexer *, fdm_fnpf *, int i, int j);
    void pyramid_level_theta(lexer *);

    void lag_advance(lexer *);
    void omega_celerity(lexer *, fdm_fnpf *, ghostcell *);
    void update_phase_refs(lexer *, fdm_fnpf *);
    void criterion(lexer *, fdm_fnpf *);
    void seeds(lexer *, fdm_fnpf *, const bart_config &, std::vector<bart_seed> &);

    bart_chain ch;

    slice4 *Hx, *Hy, *Hxy;  // Hilbert transforms of eta, A380 1
    slice4 *fo1, *fo2;      // Riesz transforms of eta, A380 2 and 3

    double *window;  // Jacobsen window on the global grid
    int window_ini;

    int pyr_lev_min, pyr_lev_max;  // active pyramid levels, A385
    int pyr_num;                   // exported components, active levels and residual
    bart_pyr_level *pyr_level[RIESZ_PYR_MAX];
    slice4 **exp_A, **exp_phase, **exp_k, **exp_theta, **exp_band, **exp_fo1, **exp_fo2;

    pyr_stage pyr_stages[2][RIESZ_PYR_MAX];  // fields of every level, [PYR_TILE or PYR_GLOBAL][l]
    std::vector<double> pyr_gather_buf;      // flat gathered level, for the Allreduce
    int pyr_gather_level;                    // first gathered level, -1 if none
    int pyr_setup_done;                      // 0 until pyr_setup has run

    // Lagrange time derivative, A384 2
    double lag_dt1, lag_dt2, lag_dt3;  // last three time intervals
    double lag_t_prev;                 // simtime of the previous refresh
    int lag_hist_n;                    // valid intervals, up to 3
    double lag_w[3];                   // weights of the three increments
    int lag_ok;                        // 0: window not usable, dispersion relation instead

    int theta_count;  // last timestep refreshed
    int theta_ini;    // 0 until the previous phases carry a value
};

#endif
