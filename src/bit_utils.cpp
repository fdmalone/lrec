#include <iostream>
#include <stdint.h>
#include "const.h"

uint16_t swap_bits(uint16_t b) {

    // Bit twidling hacks: xor swap.
    // Due to how the commutation is organised we sometime need to swap bits in the mod_bits string.

    unsigned int i = 0, j = 2; // positions of bit sequences to swap
    unsigned int n = 2;    // number of consecutive bits in each sequence
    uint16_t r;    // bit-swapped result goes here

    uint16_t x = ((b >> i) ^ (b >> j)) & ((1U << n) - 1); // XOR temporary
    r = b ^ ((x << i) | (x << j));
    return(r);
}

uint16_t rotl(uint16_t x, int n) {

    // Left circular shift. Source: Stack exchange probably.

    // In:
    //    x: integer we want to shift.
    //    n: how much we want to shift.

          return ((x << n) | (x >> (n_bits - n)));
}

uint16_t rotr(uint16_t x, int n) {

    // Right circular shift.

    // In:
    //    x: integer we want to shift.
    //    n: how much we want to shift.

          return ((x >> n) | (x << (n_bits - n)));
}

int count_bits(uint16_t inp) {

    int i;

    for (i = 0; i < inp; i++) {
        inp &= inp - 1;
    }

    return i;

}
