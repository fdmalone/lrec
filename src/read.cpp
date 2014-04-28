#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <vector>
#include <sstream>
#include <complex>
#include <fstream>
#include "const.h"

using namespace std;

void read_output(vector<int> &lengths, vector<bitint> &bas_el_d, vector<cplxd> &bas_coeff_d) {

    ifstream input_basis, input_coeffs, input_lengths;

    input_basis.open("basis_elements.dat");
    input_coeffs.open("basis_coeffs.dat");
    input_lengths.open("basis_lengths.dat");

    int inp, num_el = 0;
    bitint bas;
    complex<double> coeff;

    while (input_lengths >> inp) {
        lengths.push_back(inp);
    }
    for (int i = 0; i < lengths.size(); i++) {
        num_el += lengths[i];
    }
    cout << num_el << endl;
    for (int i = 0; i < num_el; i++) {
        input_basis >> bas;
        input_coeffs >> coeff;
        bas = (bitint)bas;
        coeff = (complex<double>)coeff;
        bas_el_d.push_back(bas);
        bas_coeff_d.push_back(coeff);
    }

    input_lengths.close();
    input_basis.close();
    input_coeffs.close();

}

void write_operator_file() {

    vector<int> lengths;
    vector<bitint> bit_str;
    vector<cplxd> coeffs;
    read_output(lengths, bit_str, coeffs);

    std::stringstream composite;
    string dot = ".";

    composite << "ops_file.";
    composite << init_basis;
    composite << dot;
    composite << n_sites;
    composite << dot;
    composite << J[0];
    composite << dot;
    composite << J[1];
    composite << dot;
    composite << J[2];
    composite << dot;
    composite << depth;

    string file_string = composite.str();
    const char* file_name = file_string.c_str();

    ofstream file;
    file.open(file_name);

    file << "Lengths: ";
    for (int i = 0; i < lengths.size(); i++) {
        file << lengths[i] << "  ";
    }
    file << endl;

    for (int i = 0; i < bit_str.size(); i++) {
        file << bit_str[i] << "   " << coeffs[i] << endl;
    }
    file.close();

}
