#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <vector>
#include <complex>
#include <fstream>

using namespace std;

void read_output(vector<int> &lengths, vector<uint16_t> &bas_el_d, vector<complex <double> > &bas_coeff_d) {

    ifstream input_basis, input_coeffs, input_lengths;

    input_basis.open("basis_elements.dat");
    input_coeffs.open("basis_coeffs.dat");
    input_lengths.open("basis_lengths.dat");

    int inp, num_el = 0;
    uint16_t bas;
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
        bas = (uint16_t)bas;
        coeff = (complex<double>)coeff;
        bas_el_d.push_back(bas);
        bas_coeff_d.push_back(coeff);
    }

    input_lengths.close();
    input_basis.close();
    input_coeffs.close();

}
