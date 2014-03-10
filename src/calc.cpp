#include <iostream>
#include <vector>
#include <stdint.h>
#include "const.h"
#include "ed_heis.h"
#include "recursion.h"
#include "spect.h"
#include "get_input.h"
#include "poly.h"

using namespace std;

void calc_type() {

    if (recursion) {
        double *cf_a, *cf_b;
        cf_a = new double[depth];
        cf_b = new double[depth];
        vector<double> mom_vec;
        if (T0 || find_overlap) {
            vector<double> ground_state;
            uint16_t *configs;
            configs = new uint16_t[n_states];
            diag_heis(ground_state, configs);
            if (find_overlap) {
                commute_wrapper_inf(cf_a, cf_b);
            }
            else {
                commute_wrapper(ground_state, configs, cf_a, cf_b);
            }
            if (dos) {
                if (find_overlap) {
                    vector<double> ovlp(depth);
                    overlap_matrix(ovlp, configs, ground_state);
                    if (random_sim) {
                        random_overlap(ovlp);
                    }
                    if (poly_spec) {
                        poly_dos(cf_a, cf_b, ovlp);
                    }
                    //dos_mat(cf_a, cf_b, ovlp);
                }
                else {
                    dos_norm(cf_a, cf_b);
                }
            }
            if (moments) {
                calc_moments(mom_vec);
            }
            delete[] configs;
        }
        else {
            commute_wrapper_inf(cf_a, cf_b);
            if (dos) {
                dos_norm(cf_a, cf_b);
            }
            if (moments) {
                calc_moments(mom_vec);
            }
        }
        delete[] cf_a;
        delete[] cf_b;
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

   set_up_system(filename);

   calc_type();

   finish_calc();

}
