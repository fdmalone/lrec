#include <iostream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include <cmath>
#include <stdint.h>
#include "const.h"
#include "read_vec.h"
#include "recursion.h"

using namespace std;
using namespace arma;

double continued_fraction(double *a, double *b) {

    // Lenz's algorithm from numerical recipes.

    int i;
    double eps;
    cplxd eta, f0, c0, d0, delta, tiny_num;

    eta = I_c*eta_g;
    tiny_num = 1e-30;
    eps = 1e-15;
    f0 = tiny_num;
    c0 = f0;
    d0 = 0;

    for (i = 0; i < depth; i++) {

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

void overlap_matrix(mat &overlap, uint16_t *configs, vector<double> gs_vec) {

    vector<int> lengths;
    vector<uint16_t> bas_el, tmp_el_a, tmp_el_b;
    vector<cplxd> bas_coeff, tmp_coeff_a, tmp_coeff_b;
    int shift_a, shift_b;

    shift_a = 0;
    shift_b = 0;

    read_output(lengths, bas_el, bas_coeff);

    for (int i = 0; i < 1; i++) {
        tmp_el_a.assign(bas_el.begin()+shift_a, bas_el.begin()+lengths[i]+shift_a);
        tmp_coeff_a.assign(bas_coeff.begin()+shift_a, bas_coeff.begin()+lengths[i]+shift_a);
        for (int j = 0; j < overlap_depth; j++) {
            tmp_el_b.assign(bas_el.begin()+shift_b, bas_el.begin()+lengths[j]+shift_b);
            tmp_coeff_b.assign(bas_coeff.begin()+shift_b, bas_coeff.begin()+lengths[j]+shift_b);
            overlap(i,j) = gs_trace(tmp_el_a, tmp_el_b, tmp_coeff_a, tmp_coeff_b, configs, gs_vec).real();
            shift_b += lengths[j];
        }
        shift_a += lengths[i];
        shift_b = 0;
    }
}

void random_overlap(mat &overlap) {

    mat av(depth, depth);
    double r;
    av = overlap;

    for (int i = 0; i < noise_its; i++) {
        for (int j = 0; j < overlap_depth; j++) {
           r = (1 + ((double)rand()/RAND_MAX - 0.5)*noise_factor);
           av(0,j) += overlap(0,j)*r;
        }
    }
    //for (int i = 0; i < overlap_depth; i++) {
    //    cout << i << "  " << av(0,i) <<"  " <<0.01*av(0,i) << endl;
    //}
    overlap = 1.0/(double)(noise_its+1) * av;

}

void dos_norm(double *a, double *b) {

    double dos;

    ofstream file;
    file.open("dos.dat");

    for (int i = 0; i < dos_its; i++) {
        omega += dos_step;
        dos = continued_fraction(a, b);
        file << setprecision(16) << omega << "   " << dos << endl;
    }

    file.close();

}
/*
void dos_noise(double factor) {

    double div, erun, omega_c, N0 = 0, dos[1024];
    vector<double> av1(n_states), av2(n_states);
    mat corr1(n_states, n_states), corr2(n_states, n_states), corr_prod(n_states, n_states), diff(n_states, n_states);
    corr1.zeros();
    corr2.zeros();
    corr_prod.zeros();
    diff.zeros();
    for (int i = 0; i < n_states; i++) { av1[i] = 0; av2[i] = 0;}
    ofstream myfile;
    myfile.open("xx_chain_replica.dat");

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
        N0 += div;
         divide(tmp_vec1, sqrt(div));
        //commute_wrapper(init_basis, 1.0);
        omega_c = -5.0;
    }
    cout << "Mean Enegy: " << erun/N0 << endl;

    omega_c = 0.0;
    for (int j = 0; j < 1024; j++) {
        omega_c += 0.005;
        myfile << omega_c << "   " << dos[j]/N_its << endl;
    }

    myfile.close();
}
*/
void dos_mat(double *a_c, double *b_c, mat input) {

    cx_mat J(depth,depth), omega(depth, depth), jinv(depth, depth);
    double ome = -4.0;
    cplxd ds;

    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < depth; j++) {
            if (j == i+1) {
                J(i,j) = b_c[i+1];
            }
        }
    }
    J += J.t();
    for (int i = 0; i < depth; i++) {
        J(i,i) = a_c[i];
        omega(i,i) = 1.0;
    }
    ofstream out;
    out.open("dos.dat");
    for (int i = 0; i < dos_its; i++) {
        ome += dos_step;
        jinv = inv(ome*omega - eta_g*I_c*omega-J);
        ds = 0.0;
        for (int j = 0; j < overlap_depth; j++) {
            ds += jinv(j,0)*input(0,j);
        }
        out << setprecision(16) << ome << "  " << ds.imag() << endl;
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
