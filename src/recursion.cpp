#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <stdint.h>
#include <fstream>
#include <bitset>
#include <armadillo>
#include <time.h>
#include "const.h"
#include "bit_utils.h"
#include "sorting.h"

using namespace std;
using namespace arma;


// reduntant.
void insert_element(vector<bitint> &a, int pos, int res, int max, bitint val) {

    int i, k;

    a.insert(a.begin() + pos, val);
}

void insert_element_ulong(vector<bitint> &a, int pos, int max, bitint val) {

    int i, k;

    for (i = max; i >= pos; i--) {
        k = i + 1;
        a[k] = a[i];
    }
    a[pos] = val;

}

// Might need these when moving to different arrays.
void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val) {

    int i, k;

    for (i = max; i >= pos; i--) {
        k = i + 1;
        a[k] = a[i];
    }
    a[pos] = val;

}
void divide(vector<double> &input, double divisor) {

    for (int i = 0; i < input.size(); i++ ) {
        input[i] = input[i]/divisor;
    }

}

void divide_c(vector<cplxd> &input, double divisor) {

    for (int i = 0; i < input.size(); i++ ) {
        input[i] = input[i]/divisor;
    }

}

double inf_trace(vector<bitint> bit_str_a, vector<bitint> bit_str_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b) {

    // Infinite temperature inner product.
    // For Heisenberg model Tr(sigma_i^{\alpha}sigma_j^{\beta}) = \delta_{i,j}\delta^{\alpha,\beta}.
    // So we only need do (c_i*)(c_i)

    // In:
    //   bit_strings: lists of bit strings we multiply together.
    //   coeff_: coefficients of bit strings.
    // Out:
    //   trace: infinite temperature trace.

    int i, j;
    cplxd trace = 0.0;

    i = 0;
    j = 0;

    do {
        if (bit_str_a[i] == bit_str_b[j]) {
            trace += conj(coeff_a[i])*coeff_b[j];
            i++;
            j++;
        }
        else if (bit_str_a[i] > bit_str_b[j]) {
            i++;
        }
        else {
            j++;
        }
    } while (i < bit_str_a.size() && j < bit_str_b.size());

    return(trace.real());
}

cplxd gs_trace(vector<bitint> input_a, vector<bitint> input_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b, bitint *ground_state, vector<double> gs_coeff) {

    int i, j, k;
    int tmp_el, bit;
    bitint onsite_sigma_a, onsite_sigma_b, onsite_sigma, basis_element;
    cplxd trace = 0;
    double sgn;
    cplxd reduced_coeff, basis_coeff;
    // Loop over Pauli matrices in product state.
    for (int iter = 0; iter < n_states; iter++) {
    for (i = 0; i < input_a.size(); i++) {
    for (j = 0; j < input_b.size(); j++) {
            basis_element = ground_state[iter];
            basis_coeff = gs_coeff[iter];
            if (abs(basis_coeff) > de) {
            for (k = 0; k < n_sites; k++) {
                // Find if spin up or down at site k.
                bit = (basis_element >> k) & 1;
                // Operate on basis element on site k with Pauli matrix.
                onsite_sigma_a = (input_a[i] >> 2*k) & on_site_mask;
                onsite_sigma_b = (input_b[j] >> 2*k) & on_site_mask;
                onsite_sigma = permute_norm(onsite_sigma_b, onsite_sigma_a, sgn);
                basis_element ^= (xor_array[onsite_sigma] << k);
                // Work out coefficient = a_i*
                if (onsite_sigma_a != 0 && onsite_sigma_b != 0 && onsite_sigma != 0) {
                    basis_coeff *= I_c*sgn*spin_coeff[onsite_sigma][bit];
                }
                else {
                    basis_coeff *= sgn*spin_coeff[onsite_sigma][bit];
                }
            }
            trace += conj(coeff_b[j])*coeff_a[i]*basis_coeff*gs_coeff[look_up_table(basis_element, ground_state)];
        }
        }
    }
    }
    if (abs(trace) < de) {
        return(0.0);
    }
    else {
        return(trace);
    }
}

void merge_lists(vector<bitint> &bit_str_new, vector<bitint> bit_str_old, vector<cplxd> &coeff_new, vector<cplxd> coeff_old, double mult_a) {

    // Merge two sorted lists. Makes use of binary searching, very similar to annihlation in HANDE.

    // In/Out:
    //   bit_str_new: List we merge into.
    //   bit_str_old: list we merge from.
    //   coeff_new: coefficients of basis elements in new list.
    //   coeff_old: coefficients of old elements.
    // In:
    //   mult_a: prefactor for adding two lists. See recursion.

    int new_min, res, pos;

    new_min = 0;

    for (int i = 0; i < bit_str_old.size(); i++) {
        res = binary_search(bit_str_new, new_min, bit_str_new.size()-1, bit_str_old[i], pos);
        // Elements is in the list.
        if (res == 1) {
            coeff_new[pos] += mult_a*coeff_old[i];
        }
        else {
            bit_str_new.insert(bit_str_new.begin() + pos, bit_str_old[i]);
            coeff_new.insert(coeff_new.begin() + pos, mult_a*coeff_old[i]);
        }
        // Since both lists are sorted we know the next element to insert is greater than the previous one.
        // So we can do binary search on a shrinking list.
        new_min = pos;
    }
}

void add_new_bit_str(bitint bits[], cplxd coeffs[], bitint rank[], int length, vector<bitint> &bit_str_mod, vector<cplxd> &coeff_mod, int &max) {

    // If bit string is already present in list add coefficients else need to insert new bit string in appropriate position in list.

    int i;
    int res, pos, min;

    min = 0;

    for (i = 0; i < length; i++) {
        // There is probably a way around this.
        // I think there is an issue with sorting of lists of size 0 so need to
        // explicitly insert the first few elements.
        if (max < 1) {
            bit_str_mod.push_back(bits[rank[i]]);
            coeff_mod.push_back(coeffs[rank[i]]);
        }
        else {
            res = binary_search(bit_str_mod, min, bit_str_mod.size()-1, bits[rank[i]], pos);
            // Basis element is already in list.
            if (res == 1) {
                coeff_mod[pos] += coeffs[rank[i]];
            }
            // List is sorted so this should go in pos found from binary search.
            // Insert takes care of all the moving about.
            else {
                bit_str_mod.insert(bit_str_mod.begin() + pos, bits[rank[i]]);
                coeff_mod.insert(coeff_mod.begin() + pos, coeffs[rank[i]]);
            }
        }
        min = pos + 1;
        max = max + 1;
    }
}

bitint merge_bits(bitint mod_bits, bitint inp_bit_str, int pos, int nn) {

    // Merge bits which are the results of [H, local_bits] back into original string.

    // In/out:
    //      mod_bits: bits we want to merge back into the original string.
    // In:
    //      inp_bit_str: original bit string.
    //      pos: the position in bit string of the active operator i.e. the one we commuted H with.
    //      nn: nearest neighbour is either on the left or the right of the active operator.

    if (nn == Left) {
        // Shift bits back to their position, they might need moving around boundaries.
        mod_bits = rotl(mod_bits, pos);
        // Merge back into input bit string.
        mod_bits |= (inp_bit_str & ~rotl(bit_mask, pos));
    }
    else {
        // The right nearest neighbour is slightly more complicated for various reasons.
        // I can't remember why this is the correct shift of the bits and amn't thinking
        // about it now. I remember this taking a bit of time.
        mod_bits = rotr(swap_bits(mod_bits), (n_bits-pos+2)%n_bits);
        mod_bits |= (inp_bit_str & ~rotr(bit_mask, (n_bits-pos+2)%n_bits));
    }
    return(mod_bits);

}

bitint comm_bits(bitint onsite_bit_str, bitint nn_bit_str, cplxd &curr_coeff, int iter, int pos, int nn) {

    // Commute active bits with the Hamiltonian. This is slightly subtle for a number of reasons.

    // In:
    //   onsite_bit_str: active bit rep of operator we want to commute with Hamiltonian.
    //   nn_bit_str: nearest neighbour of onsite_bit_str.
    //   iter: Commute with one of the pauli matrices sort of defined through this index.
    //   pos: position of something
    //   nn: Either Left of Right nearest neighbour.
    // In/Out:
    //   curr_coeff: current prefactor in front of our basis element.

    int i;
    double sgn;
    bitint sigma_nn = 0, onsite_tmp, tmp_nn;
    cplxd I(0.0,1.0), onsite_coeff;

    // Perform commutation.
    // This isn't necessarily deterministic.
    // The bit_cycle mask has two entries. This NAND operation cycles through
    // the possible outcomes of [H, \sigma_{\alpha}. We then later find out what
    // the particular sigma matrix was by knowing the input and output matrix.
    onsite_tmp = ~(onsite_bit_str & bit_cycle[iter]);
    onsite_tmp &= on_site_mask;
    // Work out resulting coefficient.
    // = J_nn/4 * 2 * I * epsilon_{nn,onsite,res}
    sigma_nn = permute(onsite_bit_str, onsite_tmp, sgn);
    curr_coeff *= 0.25*2*I*sgn*J[sigma_nn-1];
    // Deal with nearest neighbour pauli matrix reduction.
    // Three possibilites:
    // 1. Identity on nearest neighbour.
    // 2. Identical pauli matrix.
    // 3. Different Pauli Matrix.
    // sigma^r_{pos-1} [H, sigma^a_pos] ~ sigma^r_{pos-1}sigma^b_pos sigma^c_{pos-1}
    // sigma^r sigma^c = i epsilon^{rcs} sigma^s
    // Here find sigma^r sigma^c and as tmp_nn while finding contribution to
    // overall coefficient at the same time.
    // It matters if the neighbour is on our left or right.
    if (nn_bit_str == 0) {
       tmp_nn = sigma_nn;
    }
    else if (nn_bit_str == sigma_nn) {
        // sigma^2 = I
        tmp_nn = 0;
    }
    else {
        // Reduce the product of nearest neighbour pauli matrices.
        tmp_nn = permute_norm(nn_bit_str, sigma_nn, sgn);
        // This is annoying. While we have periodic boundary conditions on our lattice
        // the same doesn't hold for our string of pauli operators. So when moving the operator
        // back so it can be reduced we don't loop circularly around the end. This is very important.
        // Deserves a better explanation.
        if ((pos == 0 && nn == 1) || (pos == n_bits - 2 && nn == 0)) {
            curr_coeff *= -1.0*I*sgn*pow(-1.0, nn);
        }
        else {
            curr_coeff *= 1.0*I*sgn*pow(-1.0, nn);
        }
    }
    // If not using periodic boundary conditions we zero the coefficients which arise
    // from going over the boundary.
    if (fixed_ends) {
        if ((pos == 0 && nn == 1) || (pos == n_bits - 2 && nn == 0)) {
            curr_coeff = 0.0;
        }
    }
    // Shift back to correct place.
    tmp_nn <<= 2;
    // Merge with nearest neighbour.
    onsite_tmp |= tmp_nn;

    return(onsite_tmp);
}

int boundary(int pos, int nn) {
    if (nn == Left) {
        return((pos+2)%n_bits);
    }
    else {
        return((pos-2+n_bits)%n_bits);
    }
}

void remove_zeros(vector<bitint> &input, vector<cplxd> &coeffs) {

    // Removed basis elements from list which have zero coefficient.

    int start = 0;

    do {
        if (abs(coeffs[start]) < 1e-14) {
            coeffs.erase(coeffs.begin() + start);
            input.erase(input.begin() + start);
            start = start;
        }
        else {
            start++;
        }
    } while(start < coeffs.size());

}

void commute_wrapper(vector<double> ground_state, bitint *configs, vector<double> &cf_a, vector<double> &cf_b) {

    ofstream file1, file2, file3;
    file1.open("basis_lengths.dat");
    file2.open("basis_elements.dat");
    file3.open("basis_coeffs.dat");

    int bits, pos, nn, sig, i, num_bit_str, disp_start, disp_end, dep, disp, shift, max;
    bitint rank[4];
    bitint onsite_bits, nn_bits;
    bitint new_bits, bits_sig[4];
    double lanc_a[1000], lanc_b[1000], J = 1.0, delta;
    cplxd new_coeff, tmp_coeff, coeff_sig[4], check = 0;
    vector<bitint> bit_str_0, bit_str_i, bit_str_old;
    vector<cplxd> coeff_array_0, coeff_array_i, coeff_array_old;
    bit_str_0.push_back(init_basis);
    coeff_array_0.push_back(initial_coeff);
    i = 0;
    delta = 1;
    lanc_b[0] = gs_trace(bit_str_0, bit_str_0, coeff_array_0, coeff_array_0, configs, ground_state).real();
    cf_b[0] = lanc_b[0];
    divide_c(coeff_array_0,sqrt(lanc_b[0]));
    // Max size of space ~ (dep+1)*2Z*N_s ~ (number of matrices)*(2*connectivity)*(number of bit strings at last iteration)
    // Hopefully should reduce on reallocation of array, although probably too large at the same time.
    bit_str_i.reserve(1e5);
    coeff_array_i.reserve(1e5);

    for (dep = 0; dep < depth; dep++) {
        max = -1;
        for (bits = 0; bits < bit_str_0.size(); bits++) {
            for (pos = 0; pos < n_bits; pos = pos + 2) {
                onsite_bits = (bit_str_0[bits] >> pos) & on_site_mask;
                // [H, I] = 0.
                if (onsite_bits != 0) {
                    i = 0;
                    // Loop over neighbours.
                    tmp_coeff = coeff_array_0[bits];
                    for (nn = 0; nn < 2; nn++) {
                        nn_bits = ((bit_str_0[bits] >> boundary(pos, nn)) & on_site_mask);
                        // Perform commutation of input sigma matrix with the two other types.
                        for (sig = 0; sig < 2; sig++) {
                            // Find result of [H,sigma].
                            new_bits = comm_bits(onsite_bits, nn_bits, tmp_coeff, sig, pos, nn);
                            new_bits = merge_bits(new_bits, bit_str_0[bits], pos, nn);
                            bits_sig[i] = new_bits;
                            coeff_sig[i] = tmp_coeff;
                            tmp_coeff = coeff_array_0[bits];
                            i = i + 1;
                        }
                    }
                    // Rank new bits.
                    insertion_rank(bits_sig, rank, 4);
                    // Add new bits and coeffecients to list.
                    add_new_bit_str(bits_sig, coeff_sig, rank, 4, bit_str_i, coeff_array_i, max);
                }
            }
        }
        // a_i = Tr(Lu, u)
        remove_zeros(bit_str_i, coeff_array_i);
        //lanc_a[dep] = inf_trace(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0);
        lanc_a[dep] = gs_trace(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0, configs, ground_state).real();
        // Calculate Lu - a_i u.
        merge_lists(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0, -1.0*lanc_a[dep]);
        // Caluculate V = Lu_i - a_i u_i - b[i] u_i-1
        merge_lists(bit_str_i, bit_str_old, coeff_array_i, coeff_array_old, -1.0*lanc_b[dep]);
        // b_{i+1} = Tr(V_{i+1}, V_{i+1})
        //check = inf_trace(bit_str_i, bit_str_i, coeff_array_i, coeff_array_i);
        check = gs_trace(bit_str_i, bit_str_i, coeff_array_i, coeff_array_i, configs, ground_state).real();
        cout << dep << "   " << lanc_a[dep] << "   " << lanc_b[dep]<<"  " << lanc_b[dep]*lanc_b[dep] << endl;
        cf_a[dep] = lanc_a[dep];
        if (check.real() < 0) {
            cout << "imag" << endl;
            break;
        }
        lanc_b[dep+1] = sqrt(check.real());
        cf_b[dep+1] = lanc_b[dep+1];
        file1 << bit_str_0.size() << endl;
        for (int iter = 0; iter < bit_str_0.size(); iter++) {
            file2 << bit_str_0[iter] << endl;
            file3 << setprecision(16) << coeff_array_0[iter] << endl;
        }
        if (abs(lanc_b[dep+1]) < de) {
            cout << "terminated" << endl;
            break;
        }
        divide_c(coeff_array_i, lanc_b[dep+1]);
        remove_zeros(bit_str_i, coeff_array_i);
        bit_str_old = bit_str_0;
        bit_str_0 = bit_str_i;
        coeff_array_old = coeff_array_0;
        coeff_array_0 = coeff_array_i;
        bit_str_i.resize(0);
        coeff_array_i.resize(0);
    }

    file1.close();
    file2.close();
    file3.close();
}

void commute_wrapper_inf(vector<double> &cf_a, vector<double> &cf_b) {

    ofstream file1, file2, file3;
    file1.open("basis_lengths.dat");
    file2.open("basis_elements.dat");
    file3.open("basis_coeffs.dat");

    int bits, pos, nn, sig, i, num_bit_str, disp_start, disp_end, dep, disp, shift, max;
    bitint rank[4];
    bitint onsite_bits, nn_bits;
    bitint new_bits, bits_sig[4];
    double lanc_a[1000], lanc_b[1000], delta;
    cplxd new_coeff, tmp_coeff, coeff_sig[4], check = 0;
    vector<bitint> bit_str_0, bit_str_i, bit_str_old;
    vector<cplxd> coeff_array_0, coeff_array_i, coeff_array_old;

    bit_str_0.push_back(init_basis);
    coeff_array_0.push_back(initial_coeff);
    i = 0;
    delta = 1;
    lanc_b[0] = inf_trace(bit_str_0, bit_str_0, coeff_array_0, coeff_array_0);
    cf_b[0] = lanc_b[0];
    divide_c(coeff_array_0,sqrt(lanc_b[0]));

    clock_t t;
    float total_t = 0;
    // Max size of space ~ (dep+1)*2Z*N_s ~ (number of matrices)*(2*connectivity)*(number of bit strings at last iteration)
    // Hopefully should reduce on reallocation of array, although probably too large at the same time.
    bit_str_i.reserve(1e7);
    coeff_array_i.reserve(1e7);

    for (dep = 0; dep < depth; dep++) {
        max = -1;
         for (bits = 0; bits < bit_str_0.size(); bits++) {
            for (pos = 0; pos < n_bits; pos = pos + 2) {
                onsite_bits = (bit_str_0[bits] >> pos) & on_site_mask;
                // [H, I] = 0.
                if (onsite_bits != 0) {
                    i = 0;
                    // Loop over neighbours.
                    tmp_coeff = coeff_array_0[bits];
                    for (nn = 0; nn < 2; nn++) {
                        nn_bits = ((bit_str_0[bits] >> boundary(pos, nn)) & on_site_mask);
                        // Perform commutation of input sigma matrix with the two other types.
                        for (sig = 0; sig < 2; sig++) {
                            // Find result of [H,sigma].
                            new_bits = comm_bits(onsite_bits, nn_bits, tmp_coeff, sig, pos, nn);
                            new_bits = merge_bits(new_bits, bit_str_0[bits], pos, nn);
                            bits_sig[i] = new_bits;
                            coeff_sig[i] = tmp_coeff;
                            tmp_coeff = coeff_array_0[bits];
                            i = i + 1;
                        }
                    }
                    // Rank new bits.
                    insertion_rank(bits_sig, rank, 4);
                    // Add new bits and coeffecients to list.
                    //t = clock();
                    add_new_bit_str(bits_sig, coeff_sig, rank, 4, bit_str_i, coeff_array_i, max);
                    //t = clock() - t;
                    //total_t += (float)t/CLOCKS_PER_SEC;
                }
            }
        }
        //cout << dep << "  "  <<total_t << endl;
        //total_t = 0;
        //t = clock();
        // a_i = Tr(Lu, u)
        lanc_a[dep] = inf_trace(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0);
        // Calculate Lu - a_i u.
        merge_lists(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0, -1.0*lanc_a[dep]);
        // Caluculate V = Lu_i - a_i u_i - b[i] u_i-1
        merge_lists(bit_str_i, bit_str_old, coeff_array_i, coeff_array_old, -1.0*lanc_b[dep]);
        // b_{i+1} = Tr(V_{i+1}, V_{i+1})
        check = inf_trace(bit_str_i, bit_str_i, coeff_array_i, coeff_array_i);
        cout << dep << "    " << lanc_a[dep] << "    " << lanc_b[dep]<< "    " << lanc_b[dep]*lanc_b[dep] << endl;
        cf_a[dep] = lanc_a[dep];
        if (check.real() < 0) {
            break;
        }
        file1 << bit_str_0.size() << endl;
        //for (int iter = 0; iter < bit_str_0.size(); iter++) {
            //file2 << bit_str_0[iter] << endl;
            //file3 << setprecision(16) << coeff_array_0[iter] << endl;
        //}
        //t = clock() - t;
        //total_t += (float)t/CLOCKS_PER_SEC;
        if (abs(check.real()) < de || dep + 1 == depth) {
            cout << "terminated"<< "  " << dep <<"  " << depth << endl;
            break;
        }
        lanc_b[dep+1] = sqrt(check.real());
        cf_b[dep+1] = lanc_b[dep+1];
        divide_c(coeff_array_i, lanc_b[dep+1]);
        remove_zeros(bit_str_i, coeff_array_i);
        bit_str_old = bit_str_0;
        bit_str_0 = bit_str_i;
        coeff_array_old = coeff_array_0;
        coeff_array_0 = coeff_array_i;
        bit_str_i.resize(0);
        coeff_array_i.resize(0);
    }
    cout << "TOTAL: " << total_t << endl;
    file1.close();
    file2.close();
    file3.close();
}
