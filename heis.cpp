#include <iostream>
#include <bitset>
#include <cmath>
#include <vector>
#include <complex>
#include <stdint.h>
#include <fstream>

using namespace std;

typedef complex<double> cplxd;

// Bit Manipulation.
uint32_t rotl(uint32_t x, int n);
uint32_t rotr(uint32_t x, int n);
// Bit twiddling hacks
uint32_t swap_bits(uint32_t x);
void find_set_bits(uint32_t spins, vector<int>& positions);
int count_bits(uint32_t inp);

// Sorting;
int permute(int a, int b, double &sign);
int permute_norm(int a, int b, double &sign);
int insertion_sort(int array[], int lenght);
int insertion_rank(uint32_t array[], uint32_t rank[], int length);
int binary_search(vector<uint32_t> &a, int min, int max, uint32_t val, int &pos);
void insert_element(vector<uint32_t> &a, int pos, int res, int max, uint32_t val);
void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val);
double inf_trace(vector<uint32_t> bit_str_a, vector<uint32_t> bit_str_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b);
void merge_lists(vector<uint32_t> &bit_str_new, vector<uint32_t> bit_str_old, vector<cplxd> &coeff_new, vector<cplxd> coeff_old, double mult_a);
void divide(vector<cplxd> &input, double divisor);
int look_up_table(uint32_t input, uint32_t arr[]);



// Commutation/recursion.
uint32_t merge_bits(uint32_t mod_bits, uint32_t input_bit_str, int pos, int nn);
int mask_out_left(uint32_t inp_bit_str, int pos);
uint32_t comm_bits(uint32_t onsite_bit_str, uint32_t nn_bit_str, cplxd &curr_coeff, int iter, int pos, int nn);
void commute_wrapper(uint32_t initial_bit_str, cplxd initial_coeff);
int boundary(int pos, int nn);
void add_new_bit_str(uint32_t bits[], cplxd coeffs[], uint32_t rank[], int length, vector<uint32_t> &bit_str_mod, vector<cplxd> &coeff_mod, int &max);
void print(vector<uint32_t> input);
void print_c(vector<cplxd> input);
void remove_zeros(vector<uint32_t> &input, vector<cplxd> &coeffs);

// Recursion.
double continued_fraction(double a[], double b[], double num, double omega);
double gs_trace(vector<uint32_t> input_a, vector<uint32_t> input_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b, uint32_t ground_state[], double gs_coeff[], int iter);

// Exact Diagonalisation.

enum nearest {
    Left,
    Right
};

// Bit masks.
uint32_t bit_mask = 0XF, on_site_mask = 3, nn_mask = 0XC;
uint32_t bit_cycle[2] = {1, 2};

// System Globals;
int n_bits = 32;
int n_sites = n_bits/2;
int n_neigh[2] = {Left, Right};
double J[3] = {1.0, 1.0, 0.0};
cplxd I_c(0.0,1.0);

// fin_trace constants

// xor_array:
// arranged: {I, sx, sy, sz}
// sx and sy flip bits.
uint32_t xor_array[4] = {0,1,1,0};
// spin_coeff:
// arranged as above for rows, columns = {up, down}.
// These are the coefficintes which result from acting on up or
// down with one of the four matrices.
cplxd spin_coeff[4][2] = {{1.0,1.0},{1.0,1.0},{I_c,I_c},{1.0,-1.0}};

int main(){

    commute_wrapper(1, 1.0);
    /*
    vector<uint32_t> a, b;
    vector<cplxd> c_a, c_b;
    c_a.push_back(1.0);
    c_b.push_back(2.0);
    a.push_back(1);
    b.push_back(3);
    uint32_t gs[4] = {0, 1, 2, 3};
    double gs_c[4] = {0.5, 0.5, 0.5, 0.5};
    cout << gs_trace(a, b, c_a, c_b, gs, gs_c, 0) << endl;
    */


}

int boundary(int pos, int nn) {
    if (nn == Left) {
        return((pos+2)%n_bits);
    }
    else {
        return((pos-2+n_bits)%n_bits);
    }
}

double gs_trace(vector<uint32_t> input_a, vector<uint32_t> input_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b, uint32_t ground_state[], double gs_coeff[], int iter) {

    int i, j, k;
    int tmp_el, bit;
    uint32_t onsite_sigma_a, onsite_sigma_b, onsite_sigma, basis_element;
    cplxd trace = 0;
    double sgn;
    cplxd reduced_coeff, basis_coeff;


    // Loop over Pauli matrices in product state.
    for (i = 0; i < input_a.size(); i++) {
        for (j = 0; j < input_b.size(); j++) {
            basis_element = ground_state[iter];
            basis_coeff = gs_coeff[iter];
            //reduced_coeff = reduce_product(input_a[i], input_b[j], coeff_a[i], coeff_b[j], reduced_bits);
            //cout << bitset<16>(input_a[i]) << "   " << bitset<16>(input_b[j]) << endl;
            for (k = 0; k < 4; k++) {
                // Find if spin up or down at site k.
                bit = (basis_element >> k)&1;
                //cout << bit << endl;
                // Operate on basis element on site k with Pauli matrix.
                onsite_sigma_a = (input_a[i] >> 2*k)&on_site_mask;
                onsite_sigma_b = (input_b[j] >> 2*k)&on_site_mask;
                onsite_sigma = permute_norm(onsite_sigma_b, onsite_sigma_a, sgn);
                //cout << bitset<16>(onsite_sigma_a) <<"   "<< bitset<16>(onsite_sigma_b) << "   " << bitset<16>(onsite_sigma) <<"   " <<basis_element << endl;
                basis_element ^= xor_array[onsite_sigma];
                //cout << bitset<16>(basis_element) <<"   " << coeff_b[j] <<"   "<< reduced_coeff<< "   "<< xor_array[onsite_sigma] << basis_element <<endl;
                // Work out coefficient = a_i*
                basis_coeff *= sgn*spin_coeff[onsite_sigma][bit];
            }
            // cout << basis_coeff <<"   "<< conj(coeff_b[j])*coeff_a[i] << "   " << gs_coeff[look_up_table(basis_element, ground_state)]  <<endl;
            trace += conj(coeff_b[j])*coeff_a[i]*basis_coeff*gs_coeff[look_up_table(basis_element, ground_state)];
        }
    }
    cout << trace.real() << "  "  << trace.imag() << endl;
    return((double)trace.real());
}

int look_up_table(uint32_t input, uint32_t arr[]) {

    for (int i = 0; i < 4; i++) {
        if (arr[i] == input) {
            return(i);
        }
    }

}

void commute_wrapper(uint32_t initial_bit_str, cplxd initial_coeff) {

    int bits, pos, nn, sig, i, num_bit_str, disp_start, disp_end, dep, disp, shift, max;
    uint32_t rank[4];
    uint32_t onsite_bits, nn_bits;
    uint32_t new_bits, bits_sig[4];
    double lanc_a[1000], lanc_b[1000], J = 1.0, delta;
    cplxd new_coeff, tmp_coeff, coeff_sig[4];
    vector<uint32_t> bit_str_0, bit_str_i, bit_str_old;
    vector<cplxd> coeff_array_0, coeff_array_i, coeff_array_old;

    bit_str_0.push_back(initial_bit_str);
    coeff_array_0.push_back(initial_coeff);
    i = 0;
    delta = 1;
    lanc_b[0] = inf_trace(bit_str_0, bit_str_0, coeff_array_0, coeff_array_0);
    cout << lanc_b[0] << endl;

    for (dep = 0; dep < 20; dep++) {
        max = -1;
        // Max size of space ~ (dep+1)*2Z*N_s ~ (number of matrices)*(2*connectivity)*(number of bit strings at last iteration)
        // Hopefully should reduce on reallocation of array, although probably too large at the same time.
        bit_str_i.reserve((dep+1)*4*bit_str_0.size());
        coeff_array_i.reserve((dep+1)*4*bit_str_0.size());
        //cout << bit_str_0.size() << "  " <<bit_str_i.capacity() <<endl;
        for (bits = 0; bits < bit_str_0.size(); bits++) {
            for (pos = 0; pos < n_bits; pos = pos + 2) {
                onsite_bits = (bit_str_0[bits] >> pos) & on_site_mask;
                // [H, I] = 0.
                if (onsite_bits != 0) {
                    i = 0;
                    // Loop over neighbours.
                    tmp_coeff = coeff_array_0[bits];
                    for (nn = 0; nn < 2; nn++) {
                        nn_bits = ((bit_str_0[bits] >> boundary(pos, nn)) & on_site_mask);
                        // Perform commutation of input sigma matrix with the two other types.
                        for (sig = 0; sig < 2; sig++) {
                            // Find result of [H,sigma].
                            new_bits = comm_bits(onsite_bits, nn_bits, tmp_coeff, sig, pos, nn);
                            new_bits = merge_bits(new_bits, bit_str_0[bits], pos, nn);
                            //cout << "finished commutation:" << "   " << bits <<"   "<< pos <<endl;
                            bits_sig[i] = new_bits;
                            coeff_sig[i] = tmp_coeff;
                            tmp_coeff = coeff_array_0[bits];
                            i = i + 1;
                        }
                    }
                    // Rank new bits.
                    insertion_rank(bits_sig, rank, 4);
                    // Add new bits and coeffecients to list.
                    add_new_bit_str(bits_sig, coeff_sig, rank, 4, bit_str_i, coeff_array_i, max);
                }
            }
        }
        // a_i = Tr(Lu, u)
        remove_zeros(bit_str_i, coeff_array_i);
        lanc_a[dep] = inf_trace(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0);
        // Calculate Lu - a_i u.
        merge_lists(bit_str_i, bit_str_0, coeff_array_i, coeff_array_0, -1.0*lanc_a[dep]);
        // Caluculate V = Lu_i - a_i u_i - b[i] u_i-1
        merge_lists(bit_str_i, bit_str_old, coeff_array_i, coeff_array_old, -1.0*lanc_b[dep]);
        // b_{i+1} = Tr(V_{i+1}, V_{i+1})
        lanc_b[dep+1] = sqrt(inf_trace(bit_str_i, bit_str_i, coeff_array_i, coeff_array_i));
        //cout << lanc_a[dep] << "   " << lanc_b[dep]<<"  " <<lanc_b[dep]*lanc_b[dep]<< endl;
        cout <<dep<<"   " << lanc_a[dep]<<"   " <<lanc_b[dep+1]*lanc_b[dep+1] << endl;
        //recursion(bit_str_old, bit_str_0, bit_str_i, coeff_array_old, coeff_array_0, coeff_array_i);
        divide(coeff_array_i, lanc_b[dep+1]);
        remove_zeros(bit_str_i, coeff_array_i);
        cout <<bit_str_old.size()<<"  " << bit_str_0.size() << "  " << bit_str_i.size() <<"  "<<bit_str_i.capacity() << endl;
        bit_str_old = bit_str_0;
        bit_str_0 = bit_str_i;
        //cout << bit_str_0.size() << "  " << bit_str_i.size() << bit_str_i[0] << "  " <<bit_str_0[0]<< endl;
        coeff_array_old = coeff_array_0;
        coeff_array_0 = coeff_array_i;
        bit_str_i.resize(0);
        coeff_array_i.resize(0);
    }

    /*
    ofstream myfile;
    myfile.open("spect.dat");
    double omega = -5, res;
    for (i = 0; i < 50; i++) {
        lanc_a[i] = 0;
        lanc_b[i] = 0.5*(i+1);
        cout << lanc_a[i] << lanc_b[i] << endl;
    }
    for (i = 0; i < 1024; i++) {
        omega = omega + 0.01;
        res = continued_fraction(lanc_a, lanc_b, 50, omega);
        myfile << omega << "   " << res << "   " << 2*sqrt(3.14159)*exp(-omega*omega) << endl;
    }
    myfile.close();
    */

}

void remove_zeros(vector<uint32_t> &input, vector<cplxd> &coeffs) {

    int start = 0;

    do {
        if (abs(coeffs[start]) < 1e-6) {
            coeffs.erase(coeffs.begin() + start);
            input.erase(input.begin() + start);
            start = start;
        }
        else {
            start++;
        }
    } while(start < coeffs.size());

}

double continued_fraction(double a[], double b[], double num, double omega) {

    int i;
    double eps;
    cplxd eta, f0, c0, d0, delta, tiny_num;

    eta = 0.005*I_c;
    tiny_num = 1e-30;
    eps = 1e-15;
    f0 = tiny_num;
    c0 = f0;
    d0 = 0;

    //cout << eta << "   " << tiny_num << "   " << eps << "   " << d0 << endl;

    for (i = 0; i < num; i++) {

        //cout << a[i] << "   " << b[i] << endl;
        d0 = omega - eta - a[i] - b[i]*d0;

        if (abs(d0) < eps) {
            d0 = tiny_num;
        }

        c0 = omega - eta - a[i] - b[i]/c0;

        if (abs(c0) < eps) {
            c0 = tiny_num;
        }

        d0 = 1.0/d0;
        delta = c0*d0;
        f0 = f0*delta;

    }

    return(-2.0*f0.imag());

}

void print(vector<uint32_t> input) {
    for (int i = 0; i < input.size(); i++) {
        cout << i << "  " << bitset<16>(input[i]) << endl;
    }
}
void print_c(vector<cplxd> input) {
    for (int i = 0; i < input.size(); i++) {
        cout << i << "  " << input[i] << endl;
    }
}


void divide(vector<cplxd> &input, double divisor) {

    for (int i = 0; i < input.size(); i++ ) {
        input[i] = input[i]/divisor;
    }
}

double inf_trace(vector<uint32_t> bit_str_a, vector<uint32_t> bit_str_b, vector<cplxd> coeff_a, vector<cplxd> coeff_b) {

    int i, j;
    cplxd trace;

    for (i = 0; i < bit_str_a.size(); i++) {
        //cout << coeff_b[i] << endl;
        for (j =0 ; j < bit_str_b.size(); j++) {
            if (bit_str_a[i] == bit_str_b[j]) {
                trace += conj(coeff_a[i])*coeff_b[j];
            }
        }
    }
    return(trace.real());
}

void merge_lists(vector<uint32_t> &bit_str_new, vector<uint32_t> bit_str_old, vector<cplxd> &coeff_new, vector<cplxd> coeff_old, double mult_a) {

    int new_min, res, pos;

    new_min = 0;

    for (int i = 0; i < bit_str_old.size(); i++) {
        res = binary_search(bit_str_new, new_min, bit_str_new.size()-1, bit_str_old[i], pos);
        //cout << mult_a << "  " << res <<"  " <<pos<< coeff_new[pos] << endl;
        if (res == 1) {
            coeff_new[pos] += mult_a*coeff_old[i];
        }
        else {
            bit_str_new.insert(bit_str_new.begin() + pos, bit_str_old[i]);
            coeff_new.insert(coeff_new.begin() + pos, mult_a*coeff_old[i]);
        }
        new_min = pos;
    }
}

void add_new_bit_str(uint32_t bits[], cplxd coeffs[], uint32_t rank[], int length, vector<uint32_t> &bit_str_mod, vector<cplxd> &coeff_mod, int &max) {

    // If bit string is already present in list add coefficients else need to insert new bit string in appropriate position in list.
    int i;
    int res, pos;
    //cout << bit_str_mod.size() << endl;
    for (i = 0; i < length; i++) {
        // There is probably a way around this.
        if (max < 1) {
            bit_str_mod.push_back(bits[rank[i]]);
            coeff_mod.push_back(coeffs[rank[i]]);
        }
        else {
            //cout << bit_str_mod.size() << endl;
            res = binary_search(bit_str_mod, 0, bit_str_mod.size()-1, bits[rank[i]], pos);
            //cout << pos << "    "<< res <<"   "<<bits[rank[i]] << "     " << bit_str_mod[pos] << endl;
            if (res == 1) {
                coeff_mod[pos] += coeffs[rank[i]];
            }
            else {
                bit_str_mod.insert(bit_str_mod.begin() + pos, bits[rank[i]]);
                coeff_mod.insert(coeff_mod.begin() + pos, coeffs[rank[i]]);
                //insert_element(bit_str_mod, pos, res, max, bits[rank[i]]);
                //insert_element_cplxd(coeff_mod, pos, res, max, coeffs[rank[i]]);
                //max = max + 1;
            }
        }
        //for (int j = 0; j < bit_str_mod.size(); j++) {
        //    cout << j << "   " <<bit_str_mod[j] << endl;
        //}
        //cout << endl;
        max = max + 1;
        //cout <<"this: " << bits[rank[i]] << "  " <<bit_str_mod[i] <<"  "<< bit_str_mod.size()<< endl;
    }
}

int binary_search(vector<uint32_t> &a, int min, int max, uint32_t val, int &pos) {

    int mid, lo, hi, safe = 0;

    if (val > a[max]) {
        // C indexing.
        pos = max + 1;
        return 0;

    }
    else {

        lo = min;
        hi = max;

        do {
            pos = lo + ((hi-lo)/2);
            if (a[pos] == val) {
                return 1;
                break;
            }
            else if (a[pos] < val) {
                lo = pos + 1;
            }
            else {
                hi = pos;
            }
        } while (hi != lo);

        if (hi == lo) {
            if (a[hi] == val) {
                pos = hi;
                return 1;
            }
            else if (a[hi] < val) {
                pos = hi + 1;
                return 0;
            }
            else {
                pos = hi;
                return 0;
            }
        }
    }
}

void insert_element(vector<uint32_t> &a, int pos, int res, int max, uint32_t val) {

    int i, k;

    a.insert(a.begin() + pos, val);
}

void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val) {

    int i, k;

    for (i = max; i >= pos; i--) {
        k = i + 1;
        a[k] = a[i];
    }
    a[pos] = val;

}

int insertion_rank(uint32_t array[], uint32_t rank[], int length) {

    int i, j, tmp;

    for (i = 0; i < length; i++) {
        rank[i] = i;
    }
    for (i = 1; i < length; i++) {
        j = i - 1;
        tmp = rank[i];
        do {
            if ((int)(array[rank[j]] - array[tmp]) < 0) break;
            rank[j+1] = rank[j];
            j--;
        } while (j >= 0);
        rank[j+1] = tmp;
    }
}

uint32_t merge_bits(uint32_t mod_bits, uint32_t inp_bit_str, int pos, int nn) {

    if (nn == Left) {
        mod_bits = rotl(mod_bits, pos);
        mod_bits |= (inp_bit_str & ~rotl(bit_mask, pos));
    }
    else {
        mod_bits = rotr(swap_bits(mod_bits), (n_bits-pos+2)%n_bits);
        mod_bits |= (inp_bit_str & ~rotr(bit_mask, (n_bits-pos+2)%n_bits));
    }
    //cout << "mod_bits:  " << bitset<16>(mod_bits) << endl;
    return(mod_bits);

}

uint32_t swap_bits(uint32_t b) {
    unsigned int i = 0, j = 2; // positions of bit sequences to swap
    unsigned int n = 2;    // number of consecutive bits in each sequence
    uint32_t r;    // bit-swapped result goes here

    uint32_t x = ((b >> i) ^ (b >> j)) & ((1U << n) - 1); // XOR temporary
    r = b ^ ((x << i) | (x << j));
    return(r);
}

uint32_t rotl(uint32_t x, int n) {
          return ((x << n) | (x >> (n_bits - n)));
}

uint32_t rotr(uint32_t x, int n) {
          return ((x >> n) | (x << (n_bits - n)));
}

uint32_t comm_bits(uint32_t onsite_bit_str, uint32_t nn_bit_str, cplxd &curr_coeff, int iter, int pos, int nn) {

    int i;
    double sgn;
    uint32_t sigma_nn = 0, onsite_tmp, tmp_nn;
    cplxd I(0.0,1.0), onsite_coeff;

    // Perform commutation.
    onsite_tmp = ~(onsite_bit_str & bit_cycle[iter]);
    onsite_tmp &= on_site_mask;
    //cout << "onsite_tmp: onsite_tmp " << bitset<16> (onsite_bit_str)<< endl;
    // Work out resulting coefficient.
    // = J_nn/4 * 2 * I * epsilon_{nn,onsite,res}
    sigma_nn = permute(onsite_bit_str, onsite_tmp, sgn);
    //cout << sigma_nn <<"   "<<bitset<32>(onsite_bit_str)<<"   "<< J[sigma_nn] << endl;
    curr_coeff *= 0.25*2*I*sgn*J[sigma_nn-1];
    //cout << "onsite_tmp: comm_res " << bitset<16> (sigma_nn)<< endl;
    // Deal with nearest neighbour pauli matrix reduction.
    // Three possibilites:
    // 1. Identity on nearest neighbour.
    // 2. Identical pauli matrix.
    // 3. Different Pauli Matrix.
    // sigma^r_{pos-1} [H, sigma^a_pos] ~ sigma^r_{pos-1}sigma^b_pos sigma^c_{pos-1}
    // sigma^r sigma^c = i epsilon^{rcs} sigma^s
    // Here find sigma^r sigma^c and as tmp_nn while finding contribution to
    // overall coefficient at the same time.
    if (nn_bit_str == 0) {
       tmp_nn = sigma_nn;
    }
    else if (nn_bit_str == sigma_nn) {
        // sigma^2 = I
        tmp_nn = 0;
    }
    else {
        tmp_nn = permute_norm(nn_bit_str, sigma_nn, sgn);
        if ((pos == 0 && nn == 1) || (pos == n_bits - 2 && nn == 0)) {
            curr_coeff *= -1.0*I*sgn*pow(-1.0, nn);
        }
        else {
            curr_coeff *= 1.0*I*sgn*pow(-1.0, nn);
        }
    }
    // Shift back to correct place.
    tmp_nn <<= 2;
    // Merge with nearest neighbour.
    onsite_tmp |= tmp_nn;
    //cout << "onsite_tmp: " <<bitset<16>(onsite_bit_str)<<"   "<<bitset<16>(nn_bit_str)<<"  " <<bitset<16> (onsite_tmp) <<"  " <<curr_coeff<< endl;

    return(onsite_tmp);
}

void merge_bits(uint32_t input_bit_str, unsigned long int new_bits[], int pos) {

    int i;
    for (i = 0; i < 2; i++) {
        // Move new operator product back to correct position.
        new_bits[i] <<= 2*pos;
        // Mask back into original bit string.
        new_bits[i] = new_bits[i] & (input_bit_str & ~bit_mask);
    }
}

int permute(int a, int b, double &sign) {

    int epsilon[3], res;
    int i;

    epsilon[1] = a;
    epsilon[2] = b;

    for (i = 1; i < 4; i++) {
        if ((i != a) && (i != b)) {
            epsilon[0] = i;
            res = i;
        }
    }

    sign = pow(-1.0,insertion_sort(epsilon, 3));
    return (res);
}

int permute_norm(int a, int b, double &sign) {

    int epsilon[3], res;

    if (a == 0) {
        sign = 1.0;
        return(b);
    }
    else if (b == 0) {
        sign = 1.0;
        return(a);
    }
    else {
        epsilon[0] = a;
        epsilon[1] = b;

        for (int i = 1; i < 4; i++) {
            if((i != a) && (i !=b )) {
                epsilon[2] = i;
                res = i;
            }
        }

        sign = pow(-1.0, insertion_sort(epsilon, 3));
        return(res);
    }
}

int insertion_sort(int array[], int length) {

    int i, j, tmp, counter=0;

    for (i = 1; i<length; i++) {
        j = i;
        while (j > 0 && array[j-1] > array[j]) {
            tmp = array[j];
            array[j] = array[j-1];
            array[j-1] = tmp;
            j--;
            counter++;
        }
    }

    return(counter);
}

void find_set_bits(uint32_t spins, vector<int>& positions) {

    int i;

    bitset<16> SP(spins);

    for (i = 0; i < SP.size(); i++) {
        if (SP.test(i) == 1) {
            positions.push_back(i);
        }
    }

}

int count_bits(uint32_t inp) {

    int i;

    for (i = 0; i < inp; i++) {
        inp &= inp - 1;
    }

}
