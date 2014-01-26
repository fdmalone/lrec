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

    for (int i = 0; i < n_si; i++) {
        bit_i = (input >> i)&1;
        bit_j = (input >> (i+1)%n_si)&1;
        sz +=  Ji/4.0*pow(-1.0, bit_i + bit_j);
    }

    return(0*sz);

}

void hamiltonian(mat &H, uint16_t *I) {

    int bit_i, bit_j, mask;
    uint16_t K, L;

    for (int i = 0; i < n_states; i++) {
        // Diagonal.
        H(i,i) = h_zz(I, I[i]);
        for (int j = 0; j < n_si; j++) {
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

void vec_noise(vec &input, double factor) {

    ofstream file_r;
    file_r.open("rand_nums.dat", ios::out | ios::app);

    double r;

    for (int i = 0; i < n_states; i++) {
        r = ((double)rand()/RAND_MAX - 0.5)*factor;
        file_r << r << endl;
        input[i] += r;
    }

    file_r.close();

}

void energy_noise(mat H, vec input) {

    ofstream file;
    file.open("rand_noise.dat");
    int N_its = 1000;
    double e_rand, e_exact, e_av = 0;
    vec tmp = input;

    for (int i = 0; i < N_its; i++) {
        vec_noise(tmp, 0.1);
        e_exact = conv_to< double >::from(input.t()*H*input);
        e_rand = conv_to< double >::from(tmp.t()*H*tmp);
        file << i << "  " << e_exact << "   " << e_rand << endl;
        tmp = input;
        e_av += e_rand;
    }
    cout << "Average energy: " << e_av/N_its << endl;
    file.close();

}

void diag_heis(vector<double> &eigen, uint16_t *I) {

    states(I);
    mat HAM(n_states, n_states);
    HAM.zeros();
    hamiltonian(HAM, I);
    //cout << HAM << endl;
    vec e_val, gs;
    mat e_vec;
    eig_sym(e_val, e_vec, HAM);
    //cout << e_val << endl;
    gs = e_vec.col(0);
    eigen = conv_to< vector<double> >::from(gs);

}
