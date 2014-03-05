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

// Global data.
// Some defaults, fuck the user if they don't define them.
int n_sites;
double J[3];
bool fixed_ends = true;
int n_states;
int dos_its = 10000;
double omega;
int n_moments;
int max_time;
double time_step;
double eta_g = 0.005;
int N_its;
int n_bits;
int depth;
int init_basis;
cplxd initial_coeff;

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
        if (parsed[i].find(input) != std::string::npos) {
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

}

void set_up_globals(vector<string> parsed, vector<double> val) {

    n_sites = (int)val_present(parsed, val, "n_sites");
    J[0] = val_present(parsed, val, "Jx");
    J[1] = val_present(parsed, val, "Jy");
    J[2] = val_present(parsed, val, "Jz");
    fixed_ends = (bool)val_present(parsed, val, "fixed_ends");
    //n_states = (int)val_present(parsed, val, "n_states");
    dos_its = (int)val_present(parsed, val, "dos_its");
    omega = (int)val_present(parsed, val, "omega");
    n_moments = (int)val_present(parsed, val, "n_moments");
    max_time = (int)val_present(parsed, val, "max_time");
    time_step = val_present(parsed, val, "time_step");
    eta_g = val_present(parsed, val, "eta_g");
    N_its = (int)val_present(parsed, val, "N_its");
    n_bits = (int)val_present(parsed, val, "n_bits");
    depth = (int)val_present(parsed, val, "depth");
    init_basis = (int)val_present(parsed, val, "init_basis");
    initial_coeff = (cplxd)val_present(parsed, val, "initial_coeff");

}

void calc(char *filename) {

    vector<double> values;
    vector<string> calc_type;
    read_input(calc_type, values, filename);
    set_up_calc(calc_type);
    set_up_globals(calc_type, values);

}
void input_output() {

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
