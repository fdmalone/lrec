#ifndef RECURSION_H
#define RECURSION_H

#include <stdint.h>
#include <vector>

using namespace std;

// Commutation/recursion.
uint16_t merge_bits(uint16_t mod_bits, uint16_t input_bit_str, int pos, int nn);
int mask_out_left(uint16_t inp_bit_str, int pos);
uint16_t comm_bits(uint16_t onsite_bit_str, uint16_t nn_bit_str, cplxd &curr_coeff, int iter, int pos, int nn);
void commute_wrapper(vector<double> ground_state, uint16_t *configs, double *cf_a, double *cf_b);
void commute_wrapper_inf(double *cf_a, double *cf_b);
int boundary(int pos, int nn);
void add_new_bit_str(uint16_t bits[], cplxd coeffs[], uint16_t rank[], int length, vector<uint16_t> &bit_str_mod, vector<cplxd> &coeff_mod, int &max);
void print(vector<uint16_t> input);
void print_c(vector<cplxd> input);
void remove_zeros(vector<uint16_t> &input, vector<cplxd> &coeffs);
void insert_element(vector<uint16_t> &a, int pos, int res, int max, uint16_t val);
void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val);
double inf_trace(vector<uint16_t> bit_str_a, vector<uint16_t> bit_str_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b);
void merge_lists(vector<uint16_t> &bit_str_new, vector<uint16_t> bit_str_old, vector<cplxd> &coeff_new, vector<cplxd> coeff_old, double mult_a);
cplxd gs_trace(vector<uint16_t> input_a, vector<uint16_t> input_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b, uint16_t *ground_state, vector<double> gs_coeff);
void divide(vector<double> &input, double divisor);
void divide_c(vector<cplxd> &input, double divisor);

#endif /* RECURSION_H */
