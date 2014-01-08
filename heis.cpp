#include <iostream>
#include <bitset>
#include <cmath>
#include <vector>
#include <complex>
#include <stdint.h>

using namespace std;

typedef complex<double> cplxd;

// Bit Manipulation.
int count_bits(uint16_t inp);
uint16_t rotl(uint16_t x, int n);
uint16_t rotr(uint16_t x, int n);
// Bit twiddling hacks
uint16_t swap_bits(uint16_t x);
void find_set_bits(uint16_t spins, vector<int>& positions);

// Sorting;
int permute(int a, int b, double &sign);
int insertion_sort(int array[], int lenght);
int insertion_rank(uint16_t array[], uint16_t rank[], int length);
int binary_search(vector<uint16_t> &a, int min, int max, uint16_t val, int &pos);
void insert_element(vector<uint16_t> &a, int pos, int res, int max, uint16_t val);
void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val);


// Commutation/recursion.
uint16_t merge_bits(uint16_t mod_bits, uint16_t input_bit_str, int pos, int nn);
int mask_out_left(uint16_t inp_bit_str, int pos);
uint16_t comm_bits(uint16_t onsite_bit_str, uint16_t nn_bit_str, cplxd &curr_coeff, int iter);
void commute_wrapper(uint16_t initial_bit_str, cplxd initial_coeff);
int boundary(int pos, int nn);
void add_new_bit_str(uint16_t bits[], cplxd coeffs[], uint16_t rank[], int length, vector<uint16_t> &bit_str_mod, vector<cplxd> &coeff_mod, int &max);

enum nearest {
    Left,
    Right
};

// Bit masks.
uint16_t bit_mask = 0XF, on_site_mask = 3, nn_mask = 0XC;
uint16_t bit_cycle[2] = {1, 2};

// System Globals;
int n_bits = 16;
int n_sites = n_bits/2;
int n_neigh[2] = {Left, Right};

int main(){

    commute_wrapper(3, 1.0);

}

int boundary(int pos, int nn) {
    if (nn == Left) {
        return((pos+2)%n_bits);
    }
    else {
        return((pos-2+n_bits)%n_bits);
    }
}

void commute_wrapper(uint16_t initial_bit_str, cplxd initial_coeff) {

    int bits, pos, nn, sig, i, num_bit_str, disp_start, disp_end, dep, disp, shift, max;
    uint16_t rank[4];
    uint16_t onsite_bits, nn_bits;
    uint16_t new_bits, bits_sig[4];
    cplxd new_coeff, tmp_coeff, coeff_sig[4];
    vector<uint16_t> bit_str_0, bit_str_i;
    vector<cplxd> coeff_array, coeff_array_i;

    bit_str_0.push_back(initial_bit_str);
    coeff_array.push_back(initial_coeff);
    i = 0;

    for (dep = 0; dep < 6; dep++) {
        max = -1;
        for (bits = 0; bits < bit_str_0.size(); bits++) {
            // cout << bits <<"  " <<bit_str_0[bits] << endl;
            for (pos = 0; pos < n_bits; pos = pos + 2) {
                onsite_bits = (bit_str_0[bits] >> pos) & on_site_mask;
                // [H, I] = 0.
                if (onsite_bits != 0) {
                    // Loop over neighbours.
                    tmp_coeff = coeff_array[bits];
                    for (nn = 0; nn < 2; nn++) {
                        nn_bits = ((bit_str_0[bits] >> boundary(pos, nn)) & on_site_mask);
                        //cout << bitset<16>(onsite_bits) <<"  " <<bitset<16>(nn_bits)<< endl;
                        // Perform commutation of input sigma matrix with the two other types.
                        for (sig = 0; sig < 2; sig++) {
                            i = i + 1;
                            // Find result of [H,sigma].
                            new_bits = comm_bits(onsite_bits, nn_bits, tmp_coeff, sig);
                            new_bits = merge_bits(new_bits, bit_str_0[bits], pos, nn);
                            bits_sig[sig+nn] = new_bits;
                            coeff_sig[sig+nn] = tmp_coeff;
                            tmp_coeff = coeff_array[bits];
                        }
                    }
                    // Rank new bits.
                    insertion_rank(bits_sig, rank, 4);
                    // Add new bits and coeffecients to list.
                    add_new_bit_str(bits_sig, coeff_sig, rank, 4, bit_str_i, coeff_array_i, max);
                    //cout << "got to here" << endl;
                }
            //cout << i << endl;
            }
        }
        //cout << bit_str_0.size() << "  " << bit_str_i.size() << endl;
        bit_str_0 = bit_str_i;
        //bit_str_i.resize(0);
        //cout << i << endl;
    }
    //cout << disp_end << endl;
    //for (i = 0; i < disp_end; i++) {
    //    cout << i << "   " << bitset<16>(bit_str[i]) << endl;
    //}

}

void add_new_bit_str(uint16_t bits[], cplxd coeffs[], uint16_t rank[], int length, vector<uint16_t> &bit_str_mod, vector<cplxd> &coeff_mod, int &max) {

    // If bit string is already present in list add coefficients else need to insert new bit string in appropriate position in list.
    int i;
    int res, pos;

    for (i = 0 ; i < length; i++) {
        // There is probably a way around this.
        if (max < 1) {
            bit_str_mod.push_back(bits[rank[i]]);
        }
        else {
            res = binary_search(bit_str_mod, 0, max, bits[rank[i]], pos);
            if (res == 1) {
                coeff_mod[pos] += coeffs[rank[i]];
            }
            else {
                insert_element(bit_str_mod, pos, res, max, bits[rank[i]]);
                insert_element_cplxd(coeff_mod, pos, res, max, coeffs[rank[i]]);
                max = max + 1;
            }
        }
    }
}

int binary_search(vector<uint16_t> &a, int min, int max, uint16_t val, int &pos) {

    int mid, lo, hi, safe = 0;

    if (val > a[max]) {

        pos = max + 1;
        return 0;

    }
    else {

        lo = min;
        hi = max;

        do {
            pos = lo + ((hi-lo)/2);
            if (a[pos] == val) {
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

void insert_element(vector<uint16_t> &a, int pos, int res, int max, uint16_t val) {

    int i, k;

    for (i = max; i >= pos; i--) {
        k = i + 1;
        a[k] = a[i];
    }
    a[pos] = val;

}

void insert_element_cplxd(vector<cplxd> &a, int pos, int res, int max, cplxd val) {

    int i, k;

    for (i = max; i >= pos; i--) {
        k = i + 1;
        a[k] = a[i];
    }
    a[pos] = val;

}

int insertion_rank(uint16_t array[], uint16_t rank[], int length) {

    int i, j, tmp;

    for (i = 0; i < length; i++) {
        rank[i] = i;
    }
    for (i = 1; i < length; i++) {
        j = i - 1;
        tmp = rank[i];
        do {
            if (array[rank[j]] - array[tmp] < 0) break;
            rank[j+1] = rank[j];
            j--;
        } while (j >= 0);
        rank[j+1] = tmp;
    }
}

uint16_t merge_bits(uint16_t mod_bits, uint16_t inp_bit_str, int pos, int nn) {

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

uint16_t swap_bits(uint16_t b) {
    unsigned int i = 0, j = 2; // positions of bit sequences to swap
    unsigned int n = 2;    // number of consecutive bits in each sequence
    uint16_t r;    // bit-swapped result goes here

    uint16_t x = ((b >> i) ^ (b >> j)) & ((1U << n) - 1); // XOR temporary
    r = b ^ ((x << i) | (x << j));
    return(r);
}

uint16_t rotl(uint16_t x, int n) {
          return ((x << n) | (x >> (n_bits - n)));
}

uint16_t rotr(uint16_t x, int n) {
          return ((x >> n) | (x << (n_bits - n)));
}

uint16_t comm_bits(uint16_t onsite_bit_str, uint16_t nn_bit_str, cplxd &curr_coeff, int iter) {

    int i;
    double sgn;
    uint16_t sigma_nn = 0, onsite_tmp, tmp_nn;
    cplxd I(0.0,1.0), onsite_coeff;

    // Perform commutation.
    onsite_tmp = ~(onsite_bit_str & bit_cycle[iter]);
    onsite_tmp &= on_site_mask;
    //cout << "onsite_tmp: onsite_tmp " << bitset<16> (onsite_bit_str)<< endl;
    // Work out resulting coefficient.
    sigma_nn = permute(onsite_bit_str, onsite_tmp, sgn);
    curr_coeff *= 2.0*I*sgn;
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
    //   cout << "onsite_tmp: identity " << bitset<16> (tmp_nn)<< endl;
    }
    else if (nn_bit_str == sigma_nn) {
        tmp_nn = 0;
    //    cout << "onsite_tmp: same " << bitset<16> (tmp_nn)<< endl;
    }
    else {
        tmp_nn = permute(nn_bit_str, sigma_nn, sgn);
        //cout << "onsite_tmp: diff " << bitset<16> (tmp_nn)<< endl;
        curr_coeff *= I*sgn;
    }
    // Shift back to correct place.
    tmp_nn <<= 2;
    // Merge with nearest neighbour.
    onsite_tmp |= tmp_nn;
    //cout << "onsite_tmp: " << bitset<16> (onsite_tmp) << endl;

    return(onsite_tmp);
}

void merge_bits(uint16_t input_bit_str, unsigned long int new_bits[], int pos) {

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

void find_set_bits(uint16_t spins, vector<int>& positions) {

    int i;

    bitset<16> SP(spins);

    for (i = 0; i < SP.size(); i++) {
        if (SP.test(i) == 1) {
            positions.push_back(i);
        }
    }

}

int count_bits(unsigned long long int inp) {

    int i;

    for (i = 0; i < inp; i++) {
        inp &= inp - 1;
    }

}
