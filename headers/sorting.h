#ifndef SORTING_H
#define SORTING_H

#include <stdint.h>
#include <vector>
#include "const.h"

using namespace std;

// Sorting;
int permute(int a, int b, double &sign);
int permute_norm(int a, int b, double &sign);
int insertion_sort(int array[], int lenght);
int insertion_rank(bitint array[], bitint rank[], int length);
int binary_search(vector<bitint> &a, int min, int max, bitint val, int &pos);
int look_up_table(bitint input, bitint arr[]);

#endif /* SORTING_H */
