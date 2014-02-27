#ifndef CONST_H
#define CONST_H

#include <complex>
#include <stdint.h>
// System Globals;
const int n_sites = 8;
const double J[3] = {1.0, 1.0, 1.0};
const bool fixed_ends = true;
extern int n_states;
const double de = 1e-16;
const double PI = 3.1415926535897932;
const int dos_its = 100000;
const int n_moments = 16;
const int max_time = 1000;
const double time_step = 0.01;
const double eta_g = 0.0005;
const int N_its = 1000;
typedef std::complex<double> cplxd;
const cplxd I_c(0.0,1.0);
extern int n_bits;
extern int depth;
extern int init_basis;
const cplxd initial_coeff = 1.0;

#endif /* CONST_H */
