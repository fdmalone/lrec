#ifndef ED_HEIS_H
#define ED_HEIS_H

#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

int look_up(uint16_t *I, uint16_t);

double h_zz(uint16_t *I, uint16_t state);

void hamiltonian(mat &H, uint16_t *I);

int calc_ms(vector<double> input, uint16_t *I);

int count_bits(uint16_t inp);

void states(uint16_t *I);

void transition_freq(vector<double> input, double diff[]);

void non_zero_overlap(mat input, uint16_t *I, double diff[], double mag[]);

double vec_noise(vector<double> &input1, vector<double> &input2, double factor);

void correlation(vector<double> input1, vector<double> input2, vector<double> &av1, vector<double> &av2, mat &corr1, mat &corr2);

double energy_noise(vector<double> input1, vector<double> input2, int it, double run, double e_run, double n_run);

void exact_moments(double trans[], double mag[], int n);

void calc_moments(vector<double> &mom_vec);

void correlation_function_exact(double trans[], double mag[]);

void correlation_function_calc(vector<double> trans, vector<double> mag);

int count_non_zero(vec input);

void non_zero_all(mat input, uint16_t *I, vec evec, double mag[]);

void diag_heis(vector<double> &eigen, uint16_t *I);

void correlation_function_moments(vector<double> mom_vec);

#endif /* ED_HEIS_H */
