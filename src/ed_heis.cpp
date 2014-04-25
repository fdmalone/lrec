#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <stdint.h>
#include <fstream>
#include <cmath>
#include <bitset>
#include <armadillo>
#include "const.h"
#include "bit_utils.h"
#include "sorting.h"

// Checked relative to Hande ED. Factor of 4 difference in eigenvalues.
using namespace std;
using namespace arma;

mat H;
vec e_val;
int bound;

int look_up(bitint *I, bitint state) {

    for (int i = 0; i < n_states; i++) {
        if (I[i] == state) {
            return  i;
        }
    }
}

double h_zz(bitint *I, bitint input) {

    int bit_i, bit_j;
    double sz = 0;

    for (int i = 0; i < bound; i++) {
        bit_i = (input >> i) & 1;
        bit_j = (input >> (i + 1) % n_sites) & 1;
        sz +=  pow(-1.0, bit_i + bit_j);
    }

    return (0.25*J[2]*sz);

}

void hamiltonian(mat &H, bitint *I) {

    int bit_i, bit_j, mask;
    bitint K, L;

    for (int i = 0; i < n_states; i++) {
        // Diagonal.
        H(i,i) = h_zz(I, I[i]);
        for (int j = 0; j < bound; j++) {
            mask = (int)(pow(2.0, j) + pow(2.0, (j+1)%n_sites));
            K = I[i] & mask;
            L = K ^ mask;
            L = I[i] - K + L;
            if (K != 0 && K != mask) {
                H(look_up(I,L), i) += J[1]/2.0;
            }
        }

    }
}

int calc_ms(vec input, bitint *I) {

    int count1 = 0, count2 = 0;

    for (int i = 0; i < input.size(); i++) {
        if (abs(input[i]) > de) {
            count1++;
            if (count_bits(I[i]) == n_sites/2) {
                count2++;
            }
        }
    }

    if (count1 == count2) {
        return 1;
    }
    else {
        return 0;
    }

}

void states(bitint *I) {

    int r = 0;

    for (int n_up = 0; n_up <= n_sites; n_up++) {
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
        input1[i] *= (1.0 + r1);
        input2[i] *= (1.0 + r2);
        norm += input2[i]*input1[i];
    }
    return(norm);
    //file_r.close();

}

void correlation(vector<double> input1, vector<double> input2, vector<double> &av1, vector<double> &av2, mat &corr1, mat &corr2) {

    for (int i = 0; i < input1.size(); i++) {
        av1[i] += input1[i];
        av2[i] += input2[i];
        for (int j = 0; j < input1.size(); j++) {
            corr1(i,j) += input1[i]*input1[j];
            corr2(i,j) += input2[i]*input1[j];
        }
    }

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
    file << it << "  " << e_val(0)  << "   " << e_rand/run << "   " << n_run << endl;
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

int count_non_zero(vec input) {

    int count = 0;

    for (int i = 0; i < input.size(); i++) {
        if (abs(input[i]) > de) {
            count++;
        }
    }

    return(count);

}

void non_zero_overlap(mat input, bitint *I, double diff[], double mag[]) {

    bitint mask, tmp, new_conf;
    double factor[2] = {-1.0, 1.0}, res[n_states], norm = 0, res_dub;
    cplxd result;
    cx_vec mod_vec, curr_vec, v_j;
    curr_vec = conv_to<cx_vec>::from(input.col(0));
    mod_vec = 0*curr_vec;
    mask = 1;
    // Work out sigma_0^z |Psi_0>.
    for (int i = 0; i < n_states; i++) {
        tmp = I[i] & mask;
        new_conf = I[i] ^ xor_array[init_basis];
        mod_vec[look_up_table(new_conf, I)] += spin_coeff[init_basis][tmp]*curr_vec[i];
    }

    ofstream trans;
    trans.open("non_zero_open.dat");

    for (int i = 0; i < n_states; i++) {
        v_j = conv_to<cx_vec>::from(input.col(i));
        result = cdot(v_j, mod_vec);
        res_dub = (conj(result)*result).real();
        trans << i << "  " << res[i]*res[i] << endl;
        mag[i] = res_dub;
    }


    trans.close();

    ofstream nz;
    nz.open("non_zero_freq.dat");

    for (int i = 0; i < n_states; i++) {
        if (abs(res[i]) > de) {
            nz <<diff[i] << endl;
        }
    }

}

void non_zero_all(mat input, bitint *I, vec evec, double mag[]) {

    int it = 0;
    bitint mask, tmp;
    double factor[2] = {-1.0, 1.0}, norm = 0;
    double res[n_states*n_states], prod;
    vec mod_vec, curr_vec, v_j;
    mod_vec = curr_vec;
    mask = 1;

    ofstream file;
    file.open("non_zero_freq_all_prod.dat");

    for (int i = 0; i < n_states; i++) {
        mod_vec = input.col(i);
        for (int j = 0; j < n_states; j++) {
            tmp = I[j] & mask;
            mod_vec[j] *= factor[tmp];
        }
        for (int k = 0; k < n_states; k++) {
            prod = dot(input.col(k), mod_vec);
            prod *= prod;
            if (abs(prod) > 10*de) {
                file << evec[i] - evec[k] << "   "  << prod<< endl;
                it++;
            }
        }
    }

    file.close();
}

void correlation_function_moments(vector<double> mom_vec) {

    ofstream file;
    file.open("corr_func_mom.dat");
    double t = 0.0, corr, sign = 1.0, fact;

    for (int j = 0; j < (int)(max_time/time_step); j++) {
        corr = 0.0;
        fact = 1.0;
        for (int i = 0; i < mom_vec.size(); i++) {
            if ((i % 2) == 0) {
                corr += sign*pow(1.0*t, i)*mom_vec[i]/fact;
                sign *= -1.0;
            }
            fact *= (double)(i+1);
        }
        file << t << "   " << corr << endl;
        t += time_step;
    }

}

void correlation_function_exact(double trans[], double mag[]) {

    double t = 0.0;
    cplxd corr;

    ofstream file;
    file.open("correlation_function_exact.dat");

    for (int i = 0; i < (int)(max_time/time_step); i++) {
        corr = 0.0;
        for (int j = 0; j < n_states; j++) {
            corr.real() += mag[j]*cos(trans[j]*t);
            corr.imag() += mag[j]*sin(trans[j]*t);
        }
        file << t << "  " << corr.real() << "   " << corr.imag() << endl;
        t += time_step;
    }

    file.close();

}

void correlation_function_calc(vector<double> trans, vector<double> mag) {

    double t = 0.0;
    cplxd corr;

    ofstream file;
    file.open("correlation_function_calc.dat");

    for (int i = 0; i < (int)(max_time/time_step); i++) {
        corr = 0.0;
        for (int j = 0; j < mag.size(); j++) {
            corr.real() += mag[j]*cos(trans[j]*t);
            corr.imag() += mag[j]*sin(trans[j]*t);
        }
        file << t << "  " << corr.real() << "   " << corr.imag() << endl;
        t += time_step;
    }

    file.close();

}

void exact_moments(double trans[], double mag[], int n) {

    // Calculate \mu^n = int_{-\infty}^{\infty} \omega^n A(\omega).
    // Uses exact values for A(\omega).

    double mu_n;

    ofstream file;
    file.open("exact_moments.dat");

    for (int i = 0; i < n; i++) {
        mu_n = 0.0;
        for (int j = 0; j < n_states; j++) {
            mu_n += pow(trans[j], i)*mag[j];
        }
        file << i << "   " << setprecision(16) <<mu_n << endl;
    }

    file.close();

}

void calc_moments_poly() {

    ifstream inp;
    inp.open("poly_dos.dat");

    double eigen, magn, mu_n;
    vector<double> eig, mag;

    while (inp >> eigen >> magn) {
        eig.push_back(eigen);
        mag.push_back(magn);
    }

    ofstream out;
    out.open("calc_moments_poly.dat");

    for (int j = 0; j < n_moments; j++) {
        mu_n = 0;
        for (int k = 0; k < eig.size(); k++) {
             mu_n += pow(eig[k], j)*mag[k];
        }
        //mom_vec.push_back(mu_n);
        //mav[j] += mu_n;
        out << j << "  " << setprecision(16) << mu_n << endl;
    }

    out.close();

    if (corr_func) {
        correlation_function_calc(eig, mag);
        //correlation_function_moments(mom_vec);
    }

}

void conv_moments(vector<double> &a_c, vector<double> &b_c) {

    ifstream inp;
    inp.open("calc_moments_poly.dat");

    int it;
    double mom;
    vector<double> moments;

    while (inp >> it >> mom) {
        moments.push_back(mom);
    }

    inp.close();

    vector<double> M(depth), L(depth), M_mi, L_mi;

    for (int i = 0; i < depth; i++) {
        M_mi.push_back(pow(-1.0,i)*moments[i]);
        L_mi.push_back(pow(-1.0,i+1)*moments[i+1]);
    }
    L[0] = L_mi[0];
    M[0] = 1.0;

    for (int i = 1; i < depth; i++) {
        for (int j = i; j < depth-i+1; j++) {
            M[j] = L_mi[j] - L_mi[i-1]*M_mi[j]/M_mi[i-1];
        }
        for (int j = i; j < depth-i+1; j++) {
            L[j] = M[j+1]/M[i] - M_mi[j]/M_mi[i-1];
        }
        for (int j = 0; j < depth; j++) {
            M_mi[j] = M[j];
            L_mi[j] = L[j];
        }
    }

    ofstream out;
    out.open("cont_frac.dat");
    for (int i = 0; i < depth/2; i++) {
        a_c[i] = -1.0*L[i];
        b_c[i] = sqrt(M[i]);
        out << setprecision(16) << -1.0*L[i] << "  " << sqrt(M[i]) << endl;
    }
    out.close();

}

void calc_moments(vector<double> &mom_vec) {

    int it = 0;
    vector < vector<double> > dos(dos_its+1, vector<double> (2,0));
    double mu_n = 0, mu_0 = 0;
    vector<double> sub_sec(10), maxima, freq;
    ifstream inp;
    double a, b, slope1, slope2;
    inp.open("dos.dat");

    while (!inp.eof()) {
        inp >> a >> b;
        dos[it][0] = a;
        dos[it][1] = b;
        it++;
    }
    inp.close();

    ofstream out;

    for (int i = 1; i < dos_its; i++) {
        slope1 = dos[i][1] - dos[i-1][1];
        slope2 = dos[i+1][1] - dos[i][1];
        if (slope1 > 0 && slope2 < 0  || (slope1 < 0 && slope2 > 0 )) {
            maxima.push_back(dos[i][1]);
            freq.push_back(dos[i][0]);
        }
    }
    for (int i = 0; i < maxima.size(); i++) {
        mu_0 += maxima[i];
    }
    for (int i = 0; i < maxima.size(); i++) {
        maxima[i] = maxima[i]/mu_0;
    }

    out.open("calc_moments.dat");

    for (int j = 0; j < n_moments; j++) {
        mu_n = 0;
        for (int k = 0; k < maxima.size(); k++) {
             mu_n += pow(freq[k], j)*maxima[k];
        }
        mom_vec.push_back(mu_n);
        //mav[j] += mu_n;
        out << j << "  " << setprecision(16) << mu_n << endl;
    }

    out.close();
    if (corr_func) {
        correlation_function_calc(freq, maxima);
        correlation_function_moments(mom_vec);
    }


}

void diag_heis(vector<double> &eigen, bitint *I) {

    states(I);
    H.zeros(n_states, n_states);

    if (fixed_ends) {
        bound = n_sites - 1;
    }
    else {
        bound = n_sites;
    }
    hamiltonian(H, I);

    vec  gs;
    mat e_vec;
    vector<double> evalues;
    double spec_diff[n_states], spec_amp[n_states];

    eig_sym(e_val, e_vec, H);
    cout <<"Ground state energy: " << e_val(0) << "  " << e_val(1)<< endl;
    gs = e_vec.col(0);
    eigen = conv_to< vector<double> >::from(gs);
    evalues = conv_to< vector<double> >::from(e_val);
    transition_freq(evalues, spec_diff);
    non_zero_overlap(e_vec, I, spec_diff, spec_amp);
    exact_moments(spec_diff, spec_amp, n_moments);
    if (corr_func) {
        correlation_function_exact(spec_diff, spec_amp);
    }

}
