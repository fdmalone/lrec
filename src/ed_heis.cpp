#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <stdint.h>
#include <fstream>
#include <bitset>
#include <armadillo>


// Checked relative to Hande ED. Factor of 4 difference in eigenvalues.
using namespace std;
using namespace arma;

int n_si = 8;
int n_states = (int)pow(2.0, n_si);
double Ji = 1.0;
mat H(n_states, n_states);
vec e_val;
double dee = 1e-7;
int bound;

int look_up(uint16_t *I, uint16_t state) {

    for (int i = 0; i < n_states; i++) {
        if (I[i] == state) {
            return  i;
        }
    }
}

double h_zz(uint16_t *I, uint16_t input) {

    int bit_i, bit_j;
    double sz = 0;

    for (int i = 0; i < bound; i++) {
        bit_i = (input >> i)&1;
        bit_j = (input >> (i+1)%n_si)&1;
        sz +=  Ji/4.0*pow(-1.0, bit_i + bit_j);
    }

    return(0.0*sz);

}

void hamiltonian(mat &H, uint16_t *I) {

    int bit_i, bit_j, mask;
    uint16_t K, L;

    for (int i = 0; i < n_states; i++) {
        // Diagonal.
        H(i,i) = h_zz(I, I[i]);
        for (int j = 0; j < bound; j++) {
            mask = (int)(pow(2.0, j) + pow(2.0, (j+1)%n_si));
            K = I[i] & mask;
            L = K ^ mask;
            L = I[i] - K + L;
            if (K != 0 && K != mask) {
                H(look_up(I,L), i) += Ji/2.0;
            }
        }

    }
}
int count_bits(uint16_t inp) {

    int i;

    for (i = 0; i < inp; i++) {
        inp &= inp - 1;
    }

    return i;

}

void states(uint16_t *I) {

    int r = 0;

    for (int n_up = 0; n_up <= n_si; n_up++) {
        for (int i = 0; i < n_states; i++) {
            if (count_bits(i) == n_up) {
                I[r] = i;
                r++;
            }
        }
    }
 }

double vec_noise(vector<double> &input1, vector<double> &input2, double factor) {

    //ofstream file_r;
    //file_r.open("rand_nums.dat", ios::out | ios::app);

    double r1, r2, norm = 0;

    for (int i = 0; i < n_states; i++) {
        r1 = ((double)rand()/RAND_MAX - 0.5)*factor;
        r2 = ((double)rand()/RAND_MAX - 0.5)*factor;
        // psi = sum (1+\delta_i)*c_i |i>
        input1[i] *= (1.0+r1);
        input2[i] *= (1.0+r2);
        norm += input2[i]*input1[i];
    }
    return(norm);
    //file_r.close();

}

double energy_noise(vector<double> input1, vector<double> input2, int it, double run, double e_run, double n_run) {

    ofstream file;
    file.open("rand_noise_replica.dat", ios::out | ios::app);
    double e_rand, e_exact, e_av = 0;
    vec tmp1 = input1;
    vec tmp2 = input2;

    tmp1 = conv_to< vec >::from(input1);
    tmp2 = conv_to< vec >::from(input2);

    //e_exact = conv_to< double >::from(tmp2.t()*H*tmp1);
    e_rand = conv_to< double >::from(tmp2.t()*H*tmp1);
    e_run += e_rand;
    file << it << "  " << e_val(0)  << "   " << e_rand/run << "   " << e_run/n_run<< endl;
    return (e_run);
    file.close();

}

void transition_freq(vector<double> input, double diff[]) {

    ofstream freq;
    freq.open("frequencies_open.dat");

    for (int i = 0; i < n_states; i++) {
        diff[i] = abs(input[0] - input[i]);
        freq << abs(input[0] - input[i]) << endl;
    }

    freq.close();

}

void non_zero_overlap(mat input, uint16_t *I, double diff[]) {

    uint16_t mask, tmp;
    double factor[2] = {-1.0, 1.0}, res[n_states];
    vec mod_vec, curr_vec, v_j;
    curr_vec = input.col(0);
    mod_vec = curr_vec;
    mask = 1;

    // Work out sigma_0^z |Psi_0>.
    for (int i = 0; i < n_states; i++) {
        tmp = I[i] & mask;
        mod_vec[i] *= factor[tmp];
    }

    ofstream trans;
    trans.open("non_zero_open.dat");

    for (int i = 0; i < n_states; i++) {
        v_j = input.col(i);
        res[i] = dot(v_j, mod_vec);
        trans << i << "  " << res[i] << endl;
    }

    trans.close();

    ofstream nz;
    nz.open("non_zero_freq.dat");

    for (int i = 0; i < n_states; i++) {
        if (abs(res[i]) > dee) {
            nz << diff[i] << endl;
        }
    }

}


void diag_heis(vector<double> &eigen, uint16_t *I, bool boundary_type) {

    states(I);
    H.zeros();
    if (boundary_type) {
        bound = n_si - 1;
    }
    else {
        bound = n_si;
    }
    hamiltonian(H, I);

    vec  gs;
    mat e_vec;
    vector<double> evalues;
    double spec_diff[n_states];

    eig_sym(e_val, e_vec, H);
    //cout << e_val << endl;
    gs = e_vec.col(0);
    eigen = conv_to< vector<double> >::from(gs);
    evalues = conv_to< vector<double> >::from(e_val);
    //transition_freq(evalues, spec_diff);
    //non_zero_overlap(e_vec, I, spec_diff);

}
