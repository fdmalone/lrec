#include <iostream>
#include <vector>
#include <stdint.h>
#include "const.h"
#include "ed_heis.h"
#include "recursion.h"
#include "spect.h"
#include "get_input.h"
#include "poly.h"
#include "read_vec.h"

using namespace std;

void calc_type() {

    if (recursion) {
        vector<double> cf_a(depth), cf_b(depth);
        if (T0 || find_overlap) {
            vector<double> ground_state;
            vector<double> ovlp(depth);
            bitint *configs;
            configs = new bitint[n_states];
            diag_heis(ground_state, configs);
            if (find_overlap) {
                commute_wrapper_inf(cf_a, cf_b);
                if (ops_file) {
                    write_operator_file();
                    overlap_from_file(ovlp, configs, ground_state);
                }
                else {
                    overlap_matrix(ovlp, configs, ground_state);
                }
                if (random_sim) {
                    random_overlap(ovlp);
                }
                if (poly_spec) {
                    poly_dos_ovlp(cf_a, cf_b, ovlp);
                }
                if (dos) {
                    dos_mat(cf_a, cf_b, ovlp);
                }
            }
            else {
                commute_wrapper(ground_state, configs, cf_a, cf_b);
                //overlap_matrix(ovlp, configs, ground_state);
                if (dos) {
                    dos_norm(cf_a, cf_b);
                }
                if (random_sim) {
                    random_overlap(ovlp);
                }
                if (poly_spec) {
                    poly_dos(cf_a, cf_b, depth);
                }
            }
            if (moments) {
                calc_moments_poly();
                //calc_moments(mom_vec);
            }
            if (convert_moments) {
               conv_moments(cf_a, cf_b);
               poly_dos(cf_a, cf_b, depth/2);
               calc_moments_poly();
            }
            delete[] configs;
        }
        else {
            commute_wrapper_inf(cf_a, cf_b);
            if (ops_file) {
                write_operator_file();
            }
            if (poly_spec) {
                poly_dos(cf_a, cf_b, depth);
            }
            if (dos) {
                dos_norm(cf_a, cf_b);
            }
            if (moments) {
                calc_moments_poly();
            }
        }
    }

}

void finish_calc() {

    // Delete intermediate files.
    if (!keep_files) {
        remove("basis_coeffs.dat");
        remove("basis_elements.dat");
        remove("basis_lengths.dat");
        remove("frequencies_open.dat");
        remove("ms_0.dat");
        remove("non_zero_freq_all_prod.dat");
        remove("non_zero_freq.dat");
        remove("non_zero_open.dat");
    }
}

void do_calc(char *filename) {

    if (set_up_system(filename) == 1) {
       cout << "Input error terminating calculation." << endl;
    }
    else {
        calc_type();
        finish_calc();
    }

}
