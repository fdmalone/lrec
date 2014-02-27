#include <iostream>
#include <iomanip>
#include <bitset>
#include <cmath>
#include <vector>
#include <complex>
#include <stdint.h>
#include <fstream>
#include <armadillo>
#include "const.h"
#include "recursion.h"
#include "read_vec.h"
#include "ed_heis.h"
#include "spect.h"

using namespace std;
using namespace arma;


void print_row(int size, mat input);

void input_output(double &noise);

int depth = 16;
double noise_factor;
int init_basis = 3;
int n_bits = 16;
int n_states = (int)pow(2.0, n_sites);

int main() {

    double cf_a[depth], cf_b[depth];
    vector<double> gs_vec;
    uint16_t *configs;
    configs = new uint16_t [n_states];
    diag_heis(gs_vec, configs);
    input_output(noise_factor);
    commute_wrapper(gs_vec, configs, cf_a, cf_b);
    for (int i = 0; i < depth; i++ ) {
        cout << cf_a[i] << "  " << cf_b[i] << endl;
    }
    dos_norm(dos_its, -5.0, 0.0001, depth, cf_a, cf_b);
    delete[] configs;
}

void input_output(double &noise) {

    cin >> noise;
    cout << "Peforming recursion method." << endl;
    cout << "Starting vector: " << init_basis << endl;
    cout << "Fixed end boundary conditions: " << fixed_ends << endl;
    cout << "Number of sites: " << n_sites << endl;
    cout << "Values of J_x, J_y, J_z: " << J[0] << "   " << J[1] << "   " << J[2] << endl;
    cout << "Number of states: " << n_states << endl;
    cout << "Number of iterations " << N_its << endl;
    cout << "Randomisation factor " << noise << endl;
    cout << "Recursion depth " << depth << endl;
    cout << "Number of moments calculate: " << n_moments << endl;
    cout << endl;

}
