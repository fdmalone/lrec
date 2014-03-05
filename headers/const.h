#ifndef CONST_H
#define CONST_H

#include <complex>
#include <stdint.h>

// Some typedefs.
typedef std::complex<double> cplxd;

// Actual constants.
const double de = 1e-16;
const double PI = 3.1415926535897932;
const cplxd I_c(0.0,1.0);

// System Globals;

// Number of sites.
extern int n_sites;
// Heisenberg coupling constant.
extern double J[3];
// Number of states, for exact diagonalisation.
extern int n_states;
// Number of iterations when calculating density of states.
extern int dos_its;
extern double dos_step;
// Starting point for density of states calc.
extern double omega;
// Number of moments you want to calculate.
extern int n_moments;
// Correlation functions: get until max_time*time_step;
extern int max_time;
extern double time_step;
// Positive complex number for avoiding poles.
extern double eta_g;
// Number of iterations for random noise calculations, monte carlo steps effectively.
extern int N_its;
// Number of bits to store operator string.
extern int n_bits;
// Recursion depth.
extern int depth;
// Initial basis operator.
extern int init_basis;
extern cplxd initial_coeff;
// How many overlap matrix elements we want to use.
extern int overlap_depth;
extern double noise_factor;

// Calculation options.

// Do recursion..
extern bool recursion;
extern bool exact_diag;
extern bool inf;
extern bool T0;
// Fixed end (true) or periodic (false) boundary conditions.
extern bool fixed_ends;
extern bool find_overlap;
extern bool dos;
extern bool moments;
extern bool corr_func;
extern bool keep_files;
//extern bool time;

#endif /* CONST_H */
