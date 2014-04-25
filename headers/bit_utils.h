#ifndef BIT_UTILS_H
#define BIT_UTILS_H

#include <iostream>
#include <const.h>

// Bit Manipulation.
bitint rotl(bitint x, int n);
bitint rotr(bitint x, int n);
// Bit twiddling hacks
bitint swap_bits(bitint x);
int count_bits(bitint x);

#endif /* BIT_UTILS_H */
