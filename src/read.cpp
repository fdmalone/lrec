#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <vector>
#include <complex>
#include <fstream>
#include "const.h"
#include "sorting.h"

using namespace std;

string dump_file;

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

void find_overlap_product(vector<bitint> &bit_str, vector<cplxd> &coeffs) {

    // Work out u_0^{\dagger}u_n for the overlap elements.

    // In/Out:
    //     bit_str: array of bit strings corresponding to u_n.
    //     coeffs: array of coefficients of u_n.

    bitint onsite_bits, new_bits;
    double sgn;

    for (int i = 0; i < bit_str.size(); i++) {
        // Get onsite Pauli matrix.
        onsite_bits = bit_str[i] & on_site_mask;
        // Find u_0 * \sigma_0^{\alpha}
        new_bits = permute_norm(init_basis, onsite_bits, sgn);
        // Mask back in.
        bit_str[i] &= ~on_site_mask;
        bit_str[i] |= new_bits;
        // \sigma_a \sigma_b = i \varepsilon_{abc} \sigma_c + \delta{ab} I.
        if (new_bits != 0) {
            coeffs[i] *= I_c*sgn;
        }
    }

}

void find_num_spin_flips(vector<int> &n_flips, vector<bitint> bit_str) {

    // Find the number of bit flipping matrices in each basis element.

    // In:
    //     bit_str: array of basis "elements" or strings of bits.
    // Out:
    //     n_flips: array containing number of flipping matrices for each element.

    bitint onsite_bits;
    int bit_flips;

    for (int i = 0; i < bit_str.size(); i++) {
        bit_flips = 0;
        for (int j = 0; j < n_sites; j++) {
            onsite_bits = (bit_str[i] >> 2*j) & on_site_mask;
            if (onsite_bits == 1 || onsite_bits == 2) {
                // Current pauli matrix is either \sigma_x or \sigma_y.
                bit_flips++;
            }
        }
        n_flips.push_back(bit_flips);
    }

}

struct CmpBits {

    // Use vec values as comparison function for c++ sort.
    // This is get a sorted list of indices. Taken from stackexchange.

    CmpBits(vector<int>& vec) : values(vec){}
    bool operator() (const int& a, const int& b) const
    {
        return values[a] < values[b];
    }
    vector<int>& values;
};

void sort_operator_list(vector<int> lengths, vector<bitint> bit_str, vector<int> &index, vector<int> &n_flips) {

    // Sort the operators according to the number of bit flipping matrices there are.

    // In:
    //     lenghts: number of elements in each basis function.
    //     bit_str: elements of basis functions.
    // Out:
    //     index: array containing index corresponding to sorted list.
    //     n_flips: array containing the number of bit flipping matrices.

    int step = 0;
    n_flips.reserve(bit_str.size());

    for (int i = 0; i < bit_str.size(); i++) index.push_back(i);

    find_num_spin_flips(n_flips, bit_str);

    for (int i = 0; i < lengths.size(); i++) {
        sort(index.begin()+step, index.begin()+step+lengths[i]+1, CmpBits(n_flips));
        step += lengths[i];
    }

}

void construct_file_name() {

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

    dump_file = composite.str();

}

void write_operator_file() {

    vector<int> lengths, n_flips, index;
    vector<bitint> bit_str;
    vector<cplxd> coeffs;
    read_output(lengths, bit_str, coeffs);
    construct_file_name();
    const char* file_name = dump_file.c_str();
    ofstream file;
    file.open(file_name);

    file << "Lengths: ";
    for (int i = 0; i < lengths.size(); i++) {
        file << lengths[i] << "   ";
    }
    file << endl;

    find_overlap_product(bit_str, coeffs);

    sort_operator_list(lengths, bit_str, index, n_flips);

    for (int i = 0; i < bit_str.size(); i++) {
        file << bit_str[index[i]] << "   " << coeffs[index[i]] << "   " << n_flips[index[i]]<< endl;
    }
    file.close();

}
