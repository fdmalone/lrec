#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "const.h"
#include <gsl/gsl_poly.h>

using namespace std;

int p_order = depth + 1;

void mult_const_array(vector<double> input, double mult) {

    for (int i = 0; i < input.size(); i++) {
        input[i] *= mult;
    }

}

void add_two_array(vector<double> &arr_1, vector<double> arr_2) {

    // Output is arr_1.

    for (int i = 0; i < arr_1.size(); i++) {
        arr_1[i] += arr_2[i];
    }

}

void assign_array(vector<double> &output, vector<double> input) {

    for (int i = 0; i < input.size(); i++) {
        output[i] = input[i];
    }

}

void shift_array(vector<double> &input) {

    // Shift elements in array to right.
    // Multiply polynomial by x.
    int k;

    for (int i = p_order-2; i >= 0; i--) {
        k = i + 1;
        input[k] = input[i];
    }
    input[0] = 0;

}

void print_array(vector<double> input) {

    for (int i = 0; i < input.size(); i++) {
        cout << i << "  " << input[i] << endl;
    }

}

void poly_rec(vector<double> &p_nmi, vector<double> &p_N, int p_depth, double a_n[], double b_n[]) {

    vector<double> tmp_array(p_order), p_plus(p_order), p_n(p_order), p_mi(p_order);

    mult_const_array(p_plus, 0.0);
    mult_const_array(p_n, 0.0);
    mult_const_array(p_mi, 0.0);
    p_n[0] = 1.0;
    if (p_depth == 0) {
        assign_array(p_nmi, p_mi);
        assign_array(p_N, p_n);
    }

    for (int i = 0; i < p_depth-1; i++) {
        // shift array
        assign_array(tmp_array, p_n);
        shift_array(tmp_array);
        for (int j = 0; j < p_order; j++) {
            p_plus[j] = 1.0/b_n[i+1]*(tmp_array[j] - a_n[i]*p_n[j] - b_n[i]*p_mi[j]);
        }
        assign_array(p_mi, p_n);
        assign_array(p_n, p_plus);
        if (i == p_depth-1) {
            assign_array(p_nmi, p_mi);
            assign_array(p_N, p_plus);
        }

    }

    //print_array(p_plus, len);

}

void poly_deriv(vector<double> poly) {

    for (int i = 0; i < p_order-1; i++) {
        poly[i] = (i + 1)*poly[i+1];
    }
    poly[p_order-1] = 0.0;

}

void store_norm(vector<double> norm, vector<double> p_N, vector<double> p_Nm, vector<double> eig) {

    for (int i = 0; i < depth; i++) {
        norm[i] = gsl_poly_eval(&p_Nm[0], p_order, eig[2*i])*gsl_poly_eval(&p_N[0], p_order, eig[2*i]);
    }

}

void poly_dos(double *a, double *b, vector<double> overlap) {

    vector<double> deriv(p_order), nminus(p_order);
    vector<double> norm(depth), p_n(p_order), g(depth), empty(p_order);
    vector<double> eigv(2*depth);

    mult_const_array(g, 0.0);
    // Find P_N and P_{N-1}.
    poly_rec(nminus, deriv, p_order, a, b);
    print_array(deriv);
    // Find the roots of P_N, which are the eigenvalues of T.
    gsl_poly_complex_workspace*w = gsl_poly_complex_workspace_alloc(3);
    gsl_poly_complex_solve(&deriv[0], p_order, w, &eigv[0]);
    gsl_poly_complex_workspace_free(w);

    // Work out P'_N.
    poly_deriv(deriv);
    // Find P_{N-1}(E_{\alpha}) * P'_N(E_{\aplha}) and store in norm[E_{\alpha}.
    store_norm(norm, deriv, nminus, eigv);

    // G_00 = sum_n G^_0n S_n0.
    // G^_0n(E) = sum_{\alpha} P_0(E_{\alpha}) * (E - E_{\alpha})^{-1} * P_n(E_{\alpha} / norm(E_{\alpha}).
    for (int i = 0; i < 1; i++) {
        poly_rec(empty, p_n, i, a, b);
        for (int j = 0; j < depth; j++) {
            g[j] += gsl_poly_eval(&p_n[0], p_order, eigv[2*j])*overlap[i];
        }
    }
    // Normalise, no point in doing it many times.
    for (int i = 0; i < depth; i++) {
        g[i] = g[i] / norm[i];
    }

    ofstream file;
    file.open("poly_dos.dat");
    for (int i = 0; i < depth; i++) {
        file << i << "  " << eigv[i] << "  " << g[i] << endl;
    }
    file.close();

}
