#ifndef READ_VEC_H
#define READ_VEC_H

#include <iostream>
#include <vector>
#include "const.h"

using namespace std;

void read_output(vector<int> &lengths, vector<bitint> &bas_el_d, vector<cplxd> &bas_coeff_d);

void find_overlap_product(vector<bitint> &bit_str, vector<cplxd> &coeffs);

void sort_operator_lists(vector<int> lengths, vector<bitint> bit_str, vector<cplxd> coeffs);

void write_operator_file();

#endif /* READ_VEC_H */
