#ifndef BIT_UTILS_H
#define BIT_UTILS_H

#include <iostream>
#include <stdint.h>
// Bit Manipulation.
uint16_t rotl(uint16_t x, int n);
uint16_t rotr(uint16_t x, int n);
// Bit twiddling hacks
uint16_t swap_bits(uint16_t x);
int count_bits(uint16_t x);

#endif /* BIT_UTILS_H */
