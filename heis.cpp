#include <iostream>
#include <bitset>
#include <cmath>
#include <vector>
#include <complex>
#include <stdint.h>

using namespace std;

typedef complex<double> cplxd;

enum nearest {
    Left,
    Right
};

void find_set_bits(uint16_t spins, vector<int>& positions);
int count_bits(uint16_t inp);
int permute(int a, int b, double &sign);
int insertionSort(int array[], int lenght);
uint16_t merge_bits(uint16_t mod_bits, uint16_t input_bit_str, int pos, int nn);
int mask_out_left(uint16_t inp_bit_str, int pos);
uint16_t comm_bits(uint16_t onsite_bit_str, uint16_t nn_bit_str, cplxd &curr_coeff, int iter);
void commute_wrapper(uint16_t initial_bit_str, cplxd initial_coeff);
uint16_t rotl(uint16_t x, int n);
uint16_t rotr(uint16_t x, int n);
// Bit twiddling hacks
uint16_t swap_bits(uint16_t x);

uint16_t bit_mask = 0XF, on_site_mask = 3, nn_mask = 0XC;
uint16_t bit_cycle[2] = {1, 2};
int n_bits = 16;
int n_sites = n_bits/2;
int n_neigh[2] = {Left, Right};

int main(){

    int i , j, pos;
    uint16_t new_bits[2];
    complex <double> new_coeff[2], coeff_array[10];

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

    int bits, pos, nn, sig, i, num_bit_str, disp_start, disp_end, dep, disp, shift;
    uint16_t onsite_bits, nn_bits;
    uint16_t new_bits;
    cplxd new_coeff, tmp_coeff;
    vector<uint16_t>:: iterator ins;
    vector<uint16_t> bit_str_0, bit_str_i;
    vector<cplxd> coeff_array, coeff_array_i;

    bit_str_0[0] = initial_bit_str;
    coeff_array[0] = initial_coeff;
    disp_start = 0;
    disp_end = 1;
    disp = 0;
    shift = 0;
    for (dep = 0; dep < 2; dep++) {
    for (bits = 0; bits < bit_str_0.size(); bits++) {
        i = 0;
        //cout << bits <<"  " <<disp << endl;
        for (pos = 0; pos < n_bits; pos = pos + 2) {
            onsite_bits = (bit_str_0[bits] >> pos) & on_site_mask;
            if (onsite_bits != 0) {
                // Loop over neighbours.
                tmp_coeff = coeff_array[bits];
                for (nn = 0; nn < 2; nn++) {
                    nn_bits = ((bit_str_0[bits] >> boundary(pos, nn)) & on_site_mask);
                    // Perform commutation of input sigma matrix with the two other types.
                    for (sig = 0; sig < 2; sig++) {
                        i = i + 1;
                        new_bits = comm_bits(onsite_bits, nn_bits, tmp_coeff, sig);
                        new_bits = merge_bits(new_bits, bit_str_0[bits], pos, nn);
                        ins = lower_bound(bit_str_i.begin(), bit_str_i.end(), new_bits);
                        if (bit_str_i[ins-bit_str_i.begin()] == new_bits) {
                            coeff_array_i[ins-bit_str_i.begin()] = tmp_coeff;
                        }
                        else {
                            bit_str_i.insert(ins, new_bits);
                            coeff_array_i.insert(ins, tmp_coeff);
                        }
                        tmp_coeff = coeff_array[bits];
                    }
                }
            }
        }
    }
    bit_str_0 = bit_str_i;
    fill(bit_str_i.begin(), bit_str_i.end(), 0);
    }
    cout << disp_end << endl;
    //for (i = 0; i < disp_end; i++) {
    //    cout << i << "   " << bitset<16>(bit_str[i]) << endl;
    //}

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

    sign = pow(-1.0,insertionSort(epsilon, 3));
    return (res);
}

int insertionSort(int array[], int length) {

    int i, j, tmp, counter=0;

    for(i = 1; i<length; i++) {
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
