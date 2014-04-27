#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include "const.h"

using namespace std;

// Calculation options.
// Default do nothing.
bool recursion = false;
bool exact_diag = false;
bool inf = false;
bool T0 = false;
bool find_overlap = false;
bool dos = false;
bool moments = false;
bool corr_func = false;
bool keep_files = false;
bool random_sim = false;
bool poly_spec = false;
bool convert_moments = false;

// Global data.
// Some defaults.
int n_sites;
double J[3];
bool fixed_ends = true;
int n_states;
int dos_its = 10000;
double omega;
int n_moments;
double dos_step;
double max_time;
double time_step;
double eta_g = 0.005;
int n_bits = 2*4*sizeof(bit_mask);
int depth;
int init_basis;
cplxd initial_coeff;
int noise_its;
double noise_factor;

void read_input(vector<string> &parsed, vector<double> &val, char *filename) {

    ifstream file;
    file.open(filename);

    string calc;
    string space = " ";
    double str_inp;
    int first;
    vector <string> input;
    string st[4];

    while (getline(file, calc)) {
       input.push_back(calc);
    }

    for (int i = 0; i < input.size(); i++) {
        calc = input[i];
        first = calc.find(space);
        // This is our string.
        parsed.push_back(calc.substr(0, first));
        // Convert numeric part of string to double.
        std::stringstream convert(calc.substr(first+1));
        convert >> str_inp;
        if (convert.fail()) {
            // Won't use this value so doesn't matter.
            val.push_back(0);
        }
        else {
            // Useful number.
            val.push_back(str_inp);
        }
    }
    file.close();
}

bool calc_present(vector<string> parsed, string input) {

    for (int i = 0; i < parsed.size(); i++) {
        if (parsed[i].find(input) != std::string::npos && (parsed[i].compare(input) == 0)) {
            return true;
        }
    }
    return false;
}

double val_present(vector<string> parsed, vector<double> val, string input) {

    for (int i = 0; i < parsed.size(); i++) {
        if (parsed[i].find(input) != std::string::npos) {
            return (val[i]);
        }
    }
    return 0;
}


void set_up_calc(vector<string> parsed) {

    recursion = calc_present(parsed, "recursion");
    exact_diag = calc_present(parsed, "exact_diag");
    inf = calc_present(parsed, "inf");
    T0 = calc_present(parsed, "T0");
    find_overlap = calc_present(parsed, "find_overlap");
    dos = calc_present(parsed, "dos");
    moments = calc_present(parsed, "moments");
    corr_func = calc_present(parsed, "corr_func");
    fixed_ends = calc_present(parsed, "fixed_ends");
    keep_files = calc_present(parsed, "keep_files");
    random_sim = calc_present(parsed, "random_sim");
    poly_spec = calc_present(parsed, "poly_spec");
    convert_moments = calc_present(parsed, "convert_moments");

}

void set_up_globals(vector<string> parsed, vector<double> val) {

    n_sites = (int)val_present(parsed, val, "n_sites");
    J[0] = val_present(parsed, val, "Jx");
    J[1] = val_present(parsed, val, "Jy");
    J[2] = val_present(parsed, val, "Jz");
    dos_its = (int)val_present(parsed, val, "dos_its");
    omega = (int)val_present(parsed, val, "omega");
    n_moments = (int)val_present(parsed, val, "n_moments");
    max_time = (int)val_present(parsed, val, "max_time");
    time_step = val_present(parsed, val, "time_step");
    eta_g = val_present(parsed, val, "eta_g");
    depth = (int)val_present(parsed, val, "depth");
    init_basis = (int)val_present(parsed, val, "init_basis");
    initial_coeff = (cplxd)val_present(parsed, val, "initial_coeff");
    noise_factor = val_present(parsed, val, "noise_factor");
    noise_its = (int)val_present(parsed, val, "noise_its");
    max_time = val_present(parsed, val, "max_time");
    time_step = val_present(parsed, val, "time_step");
    n_states = (int)pow(2.0, n_sites);
    dos_step = abs(2*omega)/dos_its;

}

void print_input() {

    if (recursion) {
        cout << "Peforming recursion method." << endl;
        cout << "Starting vector: " << init_basis << endl;
        if (fixed_ends) {
            cout << "Using fixed end boundary conditions." << endl;
        }
        else {
            cout << "Using periodic boundary conditions." << endl;
        }
        cout << "Number of sites: " << n_sites << endl;
        cout << "Lattice size compiled with: " << n_bits/2 << endl;
        if (n_bits/2 != n_sites) cout << "CONFLICT IN LATTICE SIZE." << endl;
        cout << "Values of J_x, J_y, J_z: " << J[0] << "   " << J[1] << "   " << J[2] << endl;
        cout << "Number of states: " << n_states << endl;
        cout << "Recursion depth: " << depth << endl;
        if (random_sim) {
            cout << "Number of iterations: " << noise_its << endl;
            cout << "Randomisation factor: " << noise_factor << endl;
        }
        if (moments) {
            cout << "Number of moments to calculate: " << n_moments << endl;
        }
        if (dos) {
            cout << "Generating spectral function graph: " << dos << endl;
            cout << "Spectral resolution: " << dos_step << endl;
            cout << "Spectral width: " << abs(2*omega) << endl;
        }
        if (corr_func) {
            cout << "Correlation function: " << corr_func << "  " << max_time/time_step << "   " << max_time << "  " << time_step << endl;
        }
        if (poly_spec) {
            cout << "Polynomials: " << boolalpha << poly_spec << endl;
        }
        cout << endl;
    }

}

int set_up_chain() {

    if (n_sites % 8 != 0 || n_sites > 32) {
        cout << "Chain length is not a multiple of 8 or exceeds max chain length of 32 sites." << endl;
        return 1;
    }
    else {
        return 0;
    }

}

int set_up_system(char *filename) {

    vector<double> values;
    vector<string> calc_type;
    read_input(calc_type, values, filename);
    set_up_calc(calc_type);
    set_up_globals(calc_type, values);
    if (set_up_chain() == 1) {
        return 1;
    }
    else {
        print_input();
    }

}
