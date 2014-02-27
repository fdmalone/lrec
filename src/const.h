#include <complex>
#include <stdint.h>
// System Globals;
const int n_sites = 8;
const double J[3] = {1.0, 1.0, 1.0};
const bool fixed_ends = true;
const int n_states = (int)pow(2.0, n_sites);
const double de = 1e-16;
const double PI = 3.1415926535897932;
const int dos_its = 100000;
const int n_moments = 16;
const int max_time = 1000;
const double time_step = 0.01;
typedef std::complex<double> cplxd;
extern int n_bits;
