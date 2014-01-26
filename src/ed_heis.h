#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

int look_up(uint16_t *I, uint16_t);

double h_zz(uint16_t *I, uint16_t state);

void hamiltonian(mat &H, uint16_t *I);

int count_bits(uint16_t inp);

void states(uint16_t *I);

double vec_noise(vector<double> &input, double factor);

void energy_noise(mat H, vec input);

void diag_heis(vector<double> &eigen, uint16_t *I);
