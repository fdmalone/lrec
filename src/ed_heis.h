#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

int look_up(uint16_t *I, uint16_t);

double h_zz(uint16_t *I, uint16_t state);

void hamiltonian(mat &H, uint16_t *I);

int count_bits(uint16_t inp);

void states(uint16_t *I);

void transition_freq(vector<double> input, double diff[]);

void non_zero_overlap(mat input, uint16_t *I, double diff[]);

double vec_noise(vector<double> &input1, vector<double> &input2, double factor);

double energy_noise(vector<double> input1, vector<double> input2, int it, double run, double e_run, double n_run);

void diag_heis(vector<double> &eigen, uint16_t *I, bool boundary_type);
