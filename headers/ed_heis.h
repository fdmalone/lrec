#ifndef ED_HEIS_H
#define ED_HEIS_H

#include <armadillo>
#include <vector>
#include "const.h"

using namespace std;
using namespace arma;

int look_up(bitint *I, bitint);

double h_zz(bitint *I, bitint state);

void hamiltonian(mat &H, bitint *I);

int calc_ms(vector<double> input, bitint *I);

int count_bits(bitint inp);

void states(bitint *I);

void transition_freq(vector<double> input, double diff[]);

void non_zero_overlap(mat input, bitint *I, double diff[], double mag[]);

double vec_noise(vector<double> &input1, vector<double> &input2, double factor);

void correlation(vector<double> input1, vector<double> input2, vector<double> &av1, vector<double> &av2, mat &corr1, mat &corr2);

double energy_noise(vector<double> input1, vector<double> input2, int it, double run, double e_run, double n_run);

void exact_moments(double trans[], double mag[], int n);

void conv_moments(vector<double> &a_c, vector<double> &b_c);

void calc_moments(vector<double> &mom_vec);

void correlation_function_exact(double trans[], double mag[]);

void correlation_function_calc(vector<double> trans, vector<double> mag);

int count_non_zero(vec input);

void calc_moments_poly();

void non_zero_all(mat input, bitint *I, vec evec, double mag[]);

void diag_heis(vector<double> &eigen, bitint *I);

void correlation_function_moments(vector<double> mom_vec);

#endif /* ED_HEIS_H */
