#ifndef RECURSION_H
#define RECURSION_H

#include <stdint.h>
#include <vector>
#include "const.h"

using namespace std;

// Commutation/recursion.
bitint merge_bits(bitint mod_bits, bitint input_bit_str, int pos, int nn);
int mask_out_left(bitint inp_bit_str, int pos);
bitint comm_bits(bitint onsite_bit_str, bitint nn_bit_str, cplxd &curr_coeff, int iter, int pos, int nn);
void commute_wrapper(vector<double> ground_state, bitint *configs, vector<double> &cf_a, vector<double> &cf_b);
void commute_wrapper_inf(vector<double> &cf_a, vector<double> &cf_b);
int boundary(int pos, int nn);
void add_new_bit_str(bitint bits[], cplxd coeffs[], bitint rank[], int length, vector<bitint> &bit_str_mod, vector<cplxd> &coeff_mod, int &max);
void print(vector<bitint> input);
void print_c(vector<cplxd> input);
void remove_zeros(vector<bitint> &input, vector<cplxd> &coeffs);
void insert_element_ulong(vector<bitint> &a, int pos, int max, bitint val);
void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val);
double inf_trace(vector<bitint> bit_str_a, vector<bitint> bit_str_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b);
void merge_lists(vector<bitint> &bit_str_new, vector<bitint> bit_str_old, vector<cplxd> &coeff_new, vector<cplxd> coeff_old, double mult_a);
cplxd gs_trace(vector<bitint> input_a, vector<bitint> input_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b, bitint *ground_state, vector<double> gs_coeff);
void divide(vector<double> &input, double divisor);
void divide_c(vector<cplxd> &input, double divisor);

#endif /* RECURSION_H */
