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
extern double max_time;
extern double time_step;
// Positive complex number for avoiding poles.
extern double eta_g;
// Number of bits to store operator string.
extern int n_bits;
// Recursion depth.
extern int depth;
// Initial basis operator.
extern int init_basis;
extern cplxd initial_coeff;
// How many overlap matrix elements we want to use.
extern int overlap_depth;
// Randomisation strength.
extern double noise_factor;
// Number of iterations for random noise calculations, monte carlo steps effectively.
extern int noise_its;


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
extern bool random_sim;
extern bool keep_files;
extern bool poly_spec;
extern bool convert_moments;
//extern bool time;

// Neighbours.
enum nearest {
    Left,
    Right
};

// Bit masks.
const uint16_t bit_mask = 0XF, on_site_mask = 3, nn_mask = 0XC;

const uint16_t bit_cycle[2] = {1, 2};

const int n_neigh[2] = {Left, Right};
// xor_array:
// arranged: {I, sx, sy, sz}
// sx and sy flip bits.
const uint16_t xor_array[4] = {0,1,1,0};
// spin_coeff:
// arranged as above for rows, columns = {down, up} to conincide with
// definition of basis in exact diagonalisation i.e. 0 = down.
// These are the coefficintes which result from acting on an up or
// down with one of the four matrices.
const cplxd spin_coeff[4][2] = {{1.0,1.0},{1.0,1.0},{-I_c,I_c},{-1.0,1.0}};

#endif /* CONST_H */
