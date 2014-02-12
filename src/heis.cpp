#include <iostream>
#include <bitset>
#include <cmath>
#include <vector>
#include <complex>
#include <stdint.h>
#include <fstream>
#include <armadillo>
#include "read_vec.h"
#include "ed_heis.h"

using namespace std;
using namespace arma;

typedef complex<double> cplxd;

// Bit Manipulation.
uint16_t rotl(uint16_t x, int n);
uint16_t rotr(uint16_t x, int n);
// Bit twiddling hacks
uint16_t swap_bits(uint16_t x);
void find_set_bits(uint16_t spins, vector<int>& positions);

// Sorting;
int permute(int a, int b, double &sign);
int permute_norm(int a, int b, double &sign);
int insertion_sort(int array[], int lenght);
int insertion_rank(uint16_t array[], uint16_t rank[], int length);
int binary_search(vector<uint16_t> &a, int min, int max, uint16_t val, int &pos);
void insert_element(vector<uint16_t> &a, int pos, int res, int max, uint16_t val);
void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val);
double inf_trace(vector<uint16_t> bit_str_a, vector<uint16_t> bit_str_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b);
void merge_lists(vector<uint16_t> &bit_str_new, vector<uint16_t> bit_str_old, vector<cplxd> &coeff_new, vector<cplxd> coeff_old, double mult_a);
void divide(vector<double> &input, double divisor);
void divide_c(vector<cplxd> &input, double divisor);
int look_up_table(uint16_t input, uint16_t arr[]);

// Commutation/recursion.
uint16_t merge_bits(uint16_t mod_bits, uint16_t input_bit_str, int pos, int nn);
int mask_out_left(uint16_t inp_bit_str, int pos);
uint16_t comm_bits(uint16_t onsite_bit_str, uint16_t nn_bit_str, cplxd &curr_coeff, int iter, int pos, int nn);
void commute_wrapper(uint16_t initial_bit_str, cplxd initial_coeff);
int boundary(int pos, int nn);
void add_new_bit_str(uint16_t bits[], cplxd coeffs[], uint16_t rank[], int length, vector<uint16_t> &bit_str_mod, vector<cplxd> &coeff_mod, int &max);
void print(vector<uint16_t> input);
void print_c(vector<cplxd> input);
void remove_zeros(vector<uint16_t> &input, vector<cplxd> &coeffs);

// Recursion.
double continued_fraction(double a[], double b[], double num, double omega);
cplxd gs_trace(vector<uint16_t> input_a, vector<uint16_t> input_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b, uint16_t *ground_state, vector<double> gs_coeff);
// Exact Diagonalisation.
void dos_mat(int L, double a_c[], double b_c[], mat input);
void dos_noise(double factor);
void dos_norm(int its, double omega, double step, int depth);
void overlap_matrix(int size, mat &overlap);
void gram_schmidt(mat overlap);

enum nearest {
    Left,
    Right
};

// Bit masks.
uint16_t bit_mask = 0XF, on_site_mask = 3, nn_mask = 0XC;
uint16_t bit_cycle[2] = {1, 2};

// System Globals;
int n_bits = 16;
int n_sites = n_bits/2;
int init_basis = 3;
int n_neigh[2] = {Left, Right};
double J[3] = {1.0, 1.0, 0.0};
int depth = 3;
double a_av[12], b_av[12];
cplxd I_c(0.0,1.0);
double eta_g = 0.01;
vector<double> gs_vec, tmp_vec1, tmp_vec2;
uint16_t *configs;
int num_states = (int)pow(2.0, n_sites);
int N_its = 5000;
bool fixed_ends = false;

double de = 1e-12;
// fin_trace constants

// xor_array:
// arranged: {I, sx, sy, sz}
// sx and sy flip bits.
uint16_t xor_array[4] = {0,1,1,0};
// spin_coeff:
// arranged as above for rows, columns = {down, up} to conincide with
// definition of basis in exact diagonalisation i.e. 0 = down.
// These are the coefficintes which result from acting on an up or
// down with one of the four matrices.
cplxd spin_coeff[4][2] = {{1.0,1.0},{1.0,1.0},{-I_c,I_c},{-1.0,1.0}};

int main() {

    configs = new uint16_t [num_states];
    diag_heis(gs_vec, configs, fixed_ends);
    double div = 0, N0 = 0, erun = 0, noise_factor;
    double dos[1024], d_av[1024], omega_c = 0.0;

    for (int j = 0; j < depth; j++) {
        a_av[j] = 0;
        b_av[j] = 0;
    }
    for (int m = 0; m < 1024; m++) {
        dos[m] = 0.0;
        d_av[m] = 0.0;
    }

    cin >> noise_factor;
    cout << "Peforming recursion method." << endl;
    cout << "Starting vector: " << init_basis << endl;
    cout << "Fixed end boundary conditions: " << fixed_ends << endl;
    cout << "Number of iterations " << N_its << endl;
    cout << "Randomisation factor " << noise_factor << endl;

    tmp_vec1 = gs_vec;
    tmp_vec2 = gs_vec;
    dos_noise(noise_factor);

}

void overlap_matrix(int size, mat &overlap) {

    vector<int> lengths;
    vector<uint16_t> bas_el, tmp_el_a, tmp_el_b;
    vector<cplxd> bas_coeff, tmp_coeff_a, tmp_coeff_b;
    int shift_a, shift_b;

    shift_a = 0;
    shift_b = 0;

    read_output(lengths, bas_el, bas_coeff);

    for (int i = 0; i < size; i++) {
        tmp_el_a.assign(bas_el.begin()+shift_a, bas_el.begin()+lengths[i]+shift_a);
        tmp_coeff_a.assign(bas_coeff.begin()+shift_a, bas_coeff.begin()+lengths[i]+shift_a);
        for (int j = 0; j < size; j++) {
            tmp_el_b.assign(bas_el.begin()+shift_b, bas_el.begin()+lengths[j]+shift_b);
            tmp_coeff_b.assign(bas_coeff.begin()+shift_b, bas_coeff.begin()+lengths[j]+shift_b);
            overlap(i,j) = gs_trace(tmp_el_a, tmp_el_b, tmp_coeff_a, tmp_coeff_b, configs, gs_vec).real();
            shift_b += lengths[j];
        }
        shift_a += lengths[i];
        shift_b = 0;
    }

}

void dos_norm(int its, double omega, double step, int depth) {

    double dos;

    ofstream file;
    file.open("dos_open_act2.dat");

    for (int i = 0; i < its; i++) {
        omega += step;
        dos = continued_fraction(a_av, b_av, depth, omega);
        file << omega << "   " << dos << endl;
    }

    file.close();

}

void dos_noise(double factor) {

    double div, erun, omega_c, N0 = 0, dos[1024];

    ofstream myfile;
    myfile.open("xx_chain_replicas.dat");

    for (int i = 0; i < 1024; i++) {
        dos[i] = 0.0;
    }
    for (int i = 0; i < N_its; i++) {
        for (int k = 0; k < depth; k++) {
            a_av[k] = 0;
            b_av[k] = 0;
        }
        tmp_vec1 = gs_vec;
        tmp_vec2 = gs_vec;
        div = vec_noise(tmp_vec1, tmp_vec2, factor);
        //cout << "div" <<"  "<< div<< endl;
        N0 += div;
        erun = energy_noise(tmp_vec1, tmp_vec2, i, div, erun, N0);
        //cout << erun << endl;
        //cout << N0 << "   " << erun << endl;
        divide(tmp_vec1, sqrt(div));
        divide(tmp_vec2, sqrt(div));
        //cout << "divided" << endl;
        commute_wrapper(3, 1.0);
        omega_c = 0.0;
        for (int l = 0; l < 1024; l++) {
            omega_c += 0.005;
            dos[l] += continued_fraction(a_av, b_av, depth, omega_c);
        }
    }
    omega_c = 0.0;
    for (int j = 0; j < 1024; j++) {
        omega_c += 0.005;
        myfile << omega_c << "   " << dos[j]/N_its << endl;
    }

    myfile.close();
}

void dos_mat(int L, double a_c[], double b_c[], mat input) {

    cx_mat J(L,L), omega(L, L), jinv(L,L);
    double ome = -5.0;
    cplxd ds;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (i == j) {
                J(i,j) =  a_c[i];
            }
            if (j == i+1 ) {
                J(i,j) = b_c[i];
            }
        }
    }

    ofstream out;
    out.open("dos_overlap.dat");
    for (int i = 0; i < 1024; i++) {
        ome += 0.01;
        jinv = inv(ome*omega - eta_g*I_c*omega-J);
        ds = 0.0;
        for (int j = 0; j < L; j++) {
            ds += jinv(j,0)*input(0,j);
        }
        out << ome << "  " << ds.imag() << endl;
    }
    out.close();

}

// Transforming bases is not a good idea.
void gram_schmidt(mat overlap) {

    int L = 7;
    cplxd norm = 0;
    cx_mat S(L,L), J(L,L), omega(L, L), jtmp(L, L), jinv(L,L), gram(L,L), g_tmp(L,L);
    vec v_tmp;
    gram.zeros();
    g_tmp.eye();
    cplxd norm_i[10];

    int step = 0;
    for (int i = 0; i < L; i++) {
        norm = 0.0;
        gram(i,i) = 1.0;
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                if (i > 0) {
                    norm_i[step] += gram(m,i-1)*gram(n,i-1)*overlap(m,n);
                }
                else {
                    norm_i[step] = 1.0;
                }
            }
        }
        norm_i[step] = sqrt(norm_i[step]);
        cout <<"norm: " <<step<<"  "<<norm_i[step] << endl;
        step++;

        for (int j = 0; j < i; j++) {
            for (int k = 0; k < i; k++) {
                for (int l = 0; l < L; l++) {
                    gram(j,i) -= gram(l,k)*overlap(i,l)*gram(j,k)/norm_i[k+1];
                    //cout << gram(j,i) << "  " << l << "   " << k << "   " << i << "   " << j << endl;
                }
            }
        }

    }

}


int boundary(int pos, int nn) {
    if (nn == Left) {
        return((pos+2)%n_bits);
    }
    else {
        return((pos-2+n_bits)%n_bits);
    }
}

cplxd gs_trace(vector<uint16_t> input_a, vector<uint16_t> input_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b, uint16_t *ground_state, vector<double> gs_coeff) {

    int i, j, k;
    int tmp_el, bit;
    uint16_t onsite_sigma_a, onsite_sigma_b, onsite_sigma, basis_element;
    cplxd trace = 0;
    double sgn;
    cplxd reduced_coeff, basis_coeff;

    // Loop over Pauli matrices in product state.
    for (int iter = 0; iter < num_states; iter++) {
    for (i = 0; i < input_a.size(); i++) {
    for (j = 0; j < input_b.size(); j++) {
            basis_element = ground_state[iter];
            //cout <<"initial: " << bitset<16>(basis_element) << endl;
            //cout << iter << endl;
            basis_coeff = gs_coeff[iter];
            //cout << bitset<16>(input_a[i]) << "   " << bitset<16>(input_b[j]) <<"  "<< bitset<4>(basis_element) << endl;
            if (abs(basis_coeff) > de) {
            for (k = 0; k < n_sites; k++) {
                // Find if spin up or down at site k.
                bit = (basis_element >> k)&1;
                //cout << bit << endl;
                // Operate on basis element on site k with Pauli matrix.
                onsite_sigma_a = (input_a[i] >> 2*k) & on_site_mask;
                onsite_sigma_b = (input_b[j] >> 2*k) & on_site_mask;
                onsite_sigma = permute_norm(onsite_sigma_b, onsite_sigma_a, sgn);
                //cout <<"onsite_sigma: "<< (int)onsite_sigma_a <<"   "<<(int)onsite_sigma_b <<"  "<< (int)onsite_sigma <<"  "<<basis_coeff<<"  "<<trace<< endl;
                //cout << bitset<16>(onsite_sigma_a) <<"   "<< bitset<16>(onsite_sigma_b) << "   " << bitset<16>(onsite_sigma) <<"   " <<basis_element << endl;
                basis_element ^= (xor_array[onsite_sigma] << k);
                //cout << bitset<16>(basis_element) <<"   " << coeff_b[j] <<"   "<< reduced_coeff<< "   "<< xor_array[onsite_sigma] << basis_element <<endl;
                // Work out coefficient = a_i*
                //cout << bitset<16>(basis_element) <<"   "<< sgn << "  "<< spin_coeff[onsite_sigma][bit]<< endl;
                if (onsite_sigma_a != 0 && onsite_sigma_b != 0 && onsite_sigma != 0) {
                    basis_coeff *= I_c*sgn*spin_coeff[onsite_sigma][bit];
                }
                else {
                    basis_coeff *= sgn*spin_coeff[onsite_sigma][bit];
                }
                //cout <<"basis_coeff: "<< basis_coeff << endl;
            }
            //cout << i << "  " << j <<"  " <<iter << "   " << trace <<"   "<< coeff_a[i] <<"  "<< coeff_b[j] << endl;
            //cout << basis_coeff << "  " << gs_coeff[look_up_table(basis_element, ground_state)] << endl;
            //cout <<"Basis element: " <<bitset<4>(basis_element) <<endl;
            //cout << basis_coeff <<"   "<< conj(coeff_b[j])*coeff_a[i] << "   " << gs_coeff[look_up_table(basis_element, ground_state)]  <<endl;
            trace += conj(coeff_b[j])*coeff_a[i]*basis_coeff*tmp_vec2[look_up_table(basis_element, ground_state)];
            //cout <<"Trace: "<< coeff_a[i] <<"  " << basis_coeff<<"  "<< look_up_table(basis_element, ground_state)<< endl;
        }
        }
    }
    }
    //cout << trace.real() << "  "  << trace.imag() << endl;
    if (abs(trace) < de) {
        return(0.0);
    }
    else {
        return(trace);
    }
}

int look_up_table(uint16_t input, uint16_t arr[]) {

    for (int i = 0; i < num_states; i++) {
        if (arr[i] == input) {
            return(i);
        }
    }

}

void commute_wrapper(uint16_t initial_bit_str, cplxd initial_coeff) {

    ofstream file1, file2, file3;
    //file.open("../T0_data/rand_coeffs.dat", ios::out | ios::app);
    file1.open("basis_lengths.dat");
    file2.open("basis_elements.dat");
    file3.open("basis_coeffs.dat");

    int bits, pos, nn, sig, i, num_bit_str, disp_start, disp_end, dep, disp, shift, max;
    uint16_t rank[4];
    uint16_t onsite_bits, nn_bits;
    uint16_t new_bits, bits_sig[4];
    double lanc_a[1000], lanc_b[1000], J = 1.0, delta;
    cplxd new_coeff, tmp_coeff, coeff_sig[4], check;
    vector<uint16_t> bit_str_0, bit_str_i, bit_str_old;
    vector<cplxd> coeff_array_0, coeff_array_i, coeff_array_old;

    bit_str_0.push_back(initial_bit_str);
    coeff_array_0.push_back(initial_coeff);
    i = 0;
    delta = 1;
    lanc_b[0] = gs_trace(bit_str_0, bit_str_0, coeff_array_0, coeff_array_0, configs, tmp_vec1).real();
    //lanc_b[0] = inf_trace(bit_str_0, bit_str_0, coeff_array_0, coeff_array_0);
    b_av[0] = lanc_b[0];
    divide_c(coeff_array_0,sqrt(lanc_b[0]));

    for (dep = 0; dep < depth; dep++) {
        max = -1;
        // Max size of space ~ (dep+1)*2Z*N_s ~ (number of matrices)*(2*connectivity)*(number of bit strings at last iteration)
        // Hopefully should reduce on reallocation of array, although probably too large at the same time.
        bit_str_i.reserve((dep+1)*4*bit_str_0.size());
        coeff_array_i.reserve((dep+1)*4*bit_str_0.size());
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
                        //cout << "    " << boundary(pos, nn) <<"  "<< pos << "  " << nn << endl;
                        // Perform commutation of input sigma matrix with the two other types.
                        for (sig = 0; sig < 2; sig++) {
                            // Find result of [H,sigma].
                            new_bits = comm_bits(onsite_bits, nn_bits, tmp_coeff, sig, pos, nn);
                            new_bits = merge_bits(new_bits, bit_str_0[bits], pos, nn);
                            //cout << "finished commutation:" << "   " << bits <<"   "<< pos <<endl;
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
        //print(bit_str_0);
        //cout << endl;
        //print(bit_str_i);
        // a_i = Tr(Lu, u)
        remove_zeros(bit_str_i, coeff_array_i);
        //print_c(coeff_array_i);
        //lanc_a[dep] = inf_trace(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0);
        lanc_a[dep] = gs_trace(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0, configs, tmp_vec1).real();
        // Calculate Lu - a_i u.
        //cout << lanc_a[dep] << endl;
        merge_lists(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0, -1.0*lanc_a[dep]);
        //print_c(coeff_array_i);
        // Caluculate V = Lu_i - a_i u_i - b[i] u_i-1
        merge_lists(bit_str_i, bit_str_old, coeff_array_i, coeff_array_old, -1.0*lanc_b[dep]);
        // b_{i+1} = Tr(V_{i+1}, V_{i+1})
        //lanc_b[dep+1] = sqrt(inf_trace(bit_str_i, bit_str_i, coeff_array_i, coeff_array_i));
        check = gs_trace(bit_str_i, bit_str_i, coeff_array_i, coeff_array_i, configs, tmp_vec1).real();
        //cout << check.real() << endl;
        if (check.real() < 0) {
            break;
        }
        //cout << lanc_a[dep] << "   " << lanc_b[dep]<<"  " <<lanc_b[dep]*lanc_b[dep]<< endl;
        //print_c(coeff_array_i);
        lanc_b[dep+1] = sqrt(check.real());
        a_av[dep] = lanc_a[dep];
        b_av[dep+1] = lanc_b[dep+1];
        //cout <<dep<< "   " <<lanc_a[dep]<< "  " <<lanc_b[dep] <<"  "<<bit_str_0.size() << endl;
        file1 << bit_str_0.size() << endl;
        for (int iter = 0; iter < bit_str_0.size(); iter++) {
            file2 << bit_str_0[iter] << endl;
            file3 << coeff_array_0[iter] << endl;
        }

        if (lanc_b[dep+1] < de) {
            // cout << "terminated" << endl;
            break;
        }
        //recursion(bit_str_old, bit_str_0, bit_str_i, coeff_array_old, coeff_array_0, coeff_array_i);
        divide_c(coeff_array_i, lanc_b[dep+1]);
        remove_zeros(bit_str_i, coeff_array_i);
        //cout <<bit_str_old.size()<<"  " << bit_str_0.size() << "  " << bit_str_i.size() <<"  "<<bit_str_i.capacity() << endl;
        bit_str_old = bit_str_0;
        bit_str_0 = bit_str_i;
        //cout << bit_str_0.size() << "  " << bit_str_i.size() << bit_str_i[0] << "  " <<bit_str_0[0]<< endl;
        coeff_array_old = coeff_array_0;
        coeff_array_0 = coeff_array_i;
        bit_str_i.resize(0);
        coeff_array_i.resize(0);
    }

    file1.close();
    file2.close();
    file3.close();


}

void remove_zeros(vector<uint16_t> &input, vector<cplxd> &coeffs) {

    int start = 0;

    do {
        if (abs(coeffs[start]) < 1e-6) {
            coeffs.erase(coeffs.begin() + start);
            input.erase(input.begin() + start);
            start = start;
        }
        else {
            start++;
        }
    } while(start < coeffs.size());

}

double continued_fraction(double a[], double b[], double num, double omega) {

    int i;
    double eps;
    cplxd eta, f0, c0, d0, delta, tiny_num;

    eta = 0.005*I_c;
    tiny_num = 1e-30;
    eps = 1e-15;
    f0 = tiny_num;
    c0 = f0;
    d0 = 0;

    for (i = 0; i < num; i++) {

        d0 = omega - eta - a[i] - b[i]*b[i]*d0;

        if (abs(d0) < eps) {
            d0 = tiny_num;
        }

        c0 = omega - eta - a[i] - b[i]*b[i]/c0;

        if (abs(c0) < eps) {
            c0 = tiny_num;
        }

        d0 = 1.0/d0;
        delta = c0*d0;
        f0 = f0*delta;

    }

    return(-2.0*f0.imag());

}

void print(vector<uint16_t> input) {
    for (int i = 0; i < input.size(); i++) {
        cout << i << "  " << bitset<16>(input[i]) << endl;
    }
}
void print_c(vector<cplxd> input) {
    for (int i = 0; i < input.size(); i++) {
        cout << i << "  " << input[i] << endl;
    }
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

double inf_trace(vector<uint16_t> bit_str_a, vector<uint16_t> bit_str_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b) {

    int i, j;
    cplxd trace;

    for (i = 0; i < bit_str_a.size(); i++) {
        //cout << coeff_b[i] << endl;
        for (j =0 ; j < bit_str_b.size(); j++) {
            if (bit_str_a[i] == bit_str_b[j]) {
                trace += conj(coeff_a[i])*coeff_b[j];
            }
        }
    }
    return(trace.real());
}

void merge_lists(vector<uint16_t> &bit_str_new, vector<uint16_t> bit_str_old, vector<cplxd> &coeff_new, vector<cplxd> coeff_old, double mult_a) {

    int new_min, res, pos;

    new_min = 0;

    for (int i = 0; i < bit_str_old.size(); i++) {
        res = binary_search(bit_str_new, new_min, bit_str_new.size()-1, bit_str_old[i], pos);
        //cout << mult_a << "  " << res <<"  " <<pos<< coeff_new[pos] << endl;
        if (res == 1) {
            coeff_new[pos] += mult_a*coeff_old[i];
        }
        else {
            bit_str_new.insert(bit_str_new.begin() + pos, bit_str_old[i]);
            coeff_new.insert(coeff_new.begin() + pos, mult_a*coeff_old[i]);
        }
        new_min = pos;
    }
}

void add_new_bit_str(uint16_t bits[], cplxd coeffs[], uint16_t rank[], int length, vector<uint16_t> &bit_str_mod, vector<cplxd> &coeff_mod, int &max) {

    // If bit string is already present in list add coefficients else need to insert new bit string in appropriate position in list.
    int i;
    int res, pos;
    //cout << bit_str_mod.size() << endl;
    for (i = 0; i < length; i++) {
        // There is probably a way around this.
        if (max < 1) {
            bit_str_mod.push_back(bits[rank[i]]);
            coeff_mod.push_back(coeffs[rank[i]]);
        }
        else {
            //cout << bit_str_mod.size() << endl;
            res = binary_search(bit_str_mod, 0, bit_str_mod.size()-1, bits[rank[i]], pos);
            //cout << pos << "    "<< res <<"   "<<bits[rank[i]] << "     " << bit_str_mod[pos] << endl;
            if (res == 1) {
                coeff_mod[pos] += coeffs[rank[i]];
            }
            else {
                bit_str_mod.insert(bit_str_mod.begin() + pos, bits[rank[i]]);
                coeff_mod.insert(coeff_mod.begin() + pos, coeffs[rank[i]]);
                //insert_element(bit_str_mod, pos, res, max, bits[rank[i]]);
                //insert_element_cplxd(coeff_mod, pos, res, max, coeffs[rank[i]]);
                //max = max + 1;
            }
        }
        //for (int j = 0; j < bit_str_mod.size(); j++) {
        //    cout << j << "   " <<bit_str_mod[j] << endl;
        //}
        //cout << endl;
        max = max + 1;
        //cout <<"this: " << bits[rank[i]] << "  " <<bit_str_mod[i] <<"  "<< bit_str_mod.size()<< endl;
    }
}

int binary_search(vector<uint16_t> &a, int min, int max, uint16_t val, int &pos) {

    int mid, lo, hi, safe = 0;

    if (val > a[max]) {
        // C indexing.
        pos = max + 1;
        return 0;

    }
    else {

        lo = min;
        hi = max;

        do {
            pos = lo + ((hi-lo)/2);
            if (a[pos] == val) {
                return 1;
                break;
            }
            else if (a[pos] < val) {
                lo = pos + 1;
            }
            else {
                hi = pos;
            }
        } while (hi != lo);

        if (hi == lo) {
            if (a[hi] == val) {
                pos = hi;
                return 1;
            }
            else if (a[hi] < val) {
                pos = hi + 1;
                return 0;
            }
            else {
                pos = hi;
                return 0;
            }
        }
    }
}

void insert_element(vector<uint16_t> &a, int pos, int res, int max, uint16_t val) {

    int i, k;

    a.insert(a.begin() + pos, val);
}

void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val) {

    int i, k;

    for (i = max; i >= pos; i--) {
        k = i + 1;
        a[k] = a[i];
    }
    a[pos] = val;

}

int insertion_rank(uint16_t array[], uint16_t rank[], int length) {

    int i, j, tmp;

    for (i = 0; i < length; i++) {
        rank[i] = i;
    }
    for (i = 1; i < length; i++) {
        j = i - 1;
        tmp = rank[i];
        do {
            if ((int)(array[rank[j]] - array[tmp]) < 0) break;
            rank[j+1] = rank[j];
            j--;
        } while (j >= 0);
        rank[j+1] = tmp;
    }
}

uint16_t merge_bits(uint16_t mod_bits, uint16_t inp_bit_str, int pos, int nn) {

    if (nn == Left) {
        mod_bits = rotl(mod_bits, pos);
        mod_bits |= (inp_bit_str & ~rotl(bit_mask, pos));
    }
    else {
        mod_bits = rotr(swap_bits(mod_bits), (n_bits-pos+2)%n_bits);
        mod_bits |= (inp_bit_str & ~rotr(bit_mask, (n_bits-pos+2)%n_bits));
    }
    //cout << "mod_bits:  " << bitset<16>(mod_bits) << endl;
    return(mod_bits);

}

uint16_t swap_bits(uint16_t b) {
    unsigned int i = 0, j = 2; // positions of bit sequences to swap
    unsigned int n = 2;    // number of consecutive bits in each sequence
    uint16_t r;    // bit-swapped result goes here

    uint16_t x = ((b >> i) ^ (b >> j)) & ((1U << n) - 1); // XOR temporary
    r = b ^ ((x << i) | (x << j));
    return(r);
}

uint16_t rotl(uint16_t x, int n) {
          return ((x << n) | (x >> (n_bits - n)));
}

uint16_t rotr(uint16_t x, int n) {
          return ((x >> n) | (x << (n_bits - n)));
}

uint16_t comm_bits(uint16_t onsite_bit_str, uint16_t nn_bit_str, cplxd &curr_coeff, int iter, int pos, int nn) {

    int i;
    double sgn;
    uint16_t sigma_nn = 0, onsite_tmp, tmp_nn;
    cplxd I(0.0,1.0), onsite_coeff;

    // Perform commutation.
    onsite_tmp = ~(onsite_bit_str & bit_cycle[iter]);
    onsite_tmp &= on_site_mask;
    //cout << "onsite_tmp: onsite_tmp " << bitset<16> (onsite_bit_str)<< endl;
    // Work out resulting coefficient.
    // = J_nn/4 * 2 * I * epsilon_{nn,onsite,res}
    sigma_nn = permute(onsite_bit_str, onsite_tmp, sgn);
    //cout << sigma_nn <<"   "<<bitset<32>(onsite_bit_str)<<"   "<< J[sigma_nn] << endl;
    curr_coeff *= 0.25*2*I*sgn*J[sigma_nn-1];
    //cout << "onsite_tmp: comm_res " << bitset<16> (sigma_nn)<< endl;
    // Deal with nearest neighbour pauli matrix reduction.
    // Three possibilites:
    // 1. Identity on nearest neighbour.
    // 2. Identical pauli matrix.
    // 3. Different Pauli Matrix.
    // sigma^r_{pos-1} [H, sigma^a_pos] ~ sigma^r_{pos-1}sigma^b_pos sigma^c_{pos-1}
    // sigma^r sigma^c = i epsilon^{rcs} sigma^s
    // Here find sigma^r sigma^c and as tmp_nn while finding contribution to
    // overall coefficient at the same time.
    if (nn_bit_str == 0) {
       tmp_nn = sigma_nn;
    }
    else if (nn_bit_str == sigma_nn) {
        // sigma^2 = I
        tmp_nn = 0;
    }
    else {
        tmp_nn = permute_norm(nn_bit_str, sigma_nn, sgn);
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
    //cout << "onsite_tmp: " <<bitset<16>(onsite_bit_str)<<"   "<<bitset<16>(nn_bit_str)<<"  " <<bitset<16> (onsite_tmp) <<"  " <<curr_coeff<< endl;

    return(onsite_tmp);
}


int permute(int a, int b, double &sign) {

    int epsilon[3], res;
    int i;

    epsilon[1] = a;
    epsilon[2] = b;

    for (i = 1; i < 4; i++) {
        if ((i != a) && (i != b)) {
            epsilon[0] = i;
            res = i;
        }
    }

    sign = pow(-1.0,insertion_sort(epsilon, 3));
    return (res);
}

int permute_norm(int a, int b, double &sign) {

    int epsilon[3], res;

    if (a == 0) {
        sign = 1.0;
        return(b);
    }
    else if (b == 0) {
        sign = 1.0;
        return(a);
    }
    else if (a == b) {
        sign = 1.0;
        return(0);
    }
    else {
        epsilon[0] = a;
        epsilon[1] = b;

        for (int i = 1; i < 4; i++) {
            if((i != a) && (i !=b )) {
                epsilon[2] = i;
                res = i;
            }
        }

        sign = pow(-1.0, insertion_sort(epsilon, 3));
        return(res);
    }
}

int insertion_sort(int array[], int length) {

    int i, j, tmp, counter=0;

    for (i = 1; i<length; i++) {
        j = i;
        while (j > 0 && array[j-1] > array[j]) {
            tmp = array[j];
            array[j] = array[j-1];
            array[j-1] = tmp;
            j--;
            counter++;
        }
    }

    return(counter);
}

void find_set_bits(uint16_t spins, vector<int>& positions) {

    int i;

    bitset<16> SP(spins);

    for (i = 0; i < SP.size(); i++) {
        if (SP.test(i) == 1) {
            positions.push_back(i);
        }
    }

}
