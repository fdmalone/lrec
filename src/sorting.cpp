#include <iostream>
#include <stdint.h>
#include <cmath>
#include <vector>
#include "const.h"

using namespace std;

int insertion_sort(int array[], int length) {

    // Sort array of integers. Used for working out permutation.

    // In:
    //   array: array of integers we want to sort.
    //   length: length of array.
    // Out:
    //   counter: even/odd for given permutations.

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

    sign = pow(-1.0, insertion_sort(epsilon, 3));
    return (res);
}

int permute_norm(int a, int b, double &sign) {

    // Work out what matrix we commuted with and also the resulting sign.
    // Technically also reduce nearest neighbour matrices. Not good practice two outputs.
    // Should reduce.

    // In:
    //   a: Matrix we input into the commutator.
    //   b: matrix we get out.
    // In/out:
    //   sign: The sign of from permutation.

    int epsilon[3], res;

    if (a == 0) {
        sign = 1.0;
        return(b);
    }
    else if (b == 0) {
        sign = 1.0;
        return(a);
    }
    else if (a == b) {
        sign = 1.0;
        return(0);
    }
    else {
        // Find out what the other matrix is, probably a clever way of doing this.
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

// Sorting / utils
// should be void.
int insertion_rank(bitint array[], bitint rank[], int length) {

    // Rank an array in increasing order result is rank which contains indices or ranked array elements.

    // In:
    //   array: array of unsigned ints to be sorted i.e. the basis vectors.
    //   length: length of array, slightly redundant potentially.
    // Out:
    //    rank: array containing indices of ranked array.

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
int binary_search(vector<bitint> &a, int min, int max, bitint val, int &pos) {

    // Search for position in sorted list for elements.
    // If it's not in the list then we find the best place to put it. Source: HANDE.

    // In/Out:
    //    a: array containing our list.
    // In:
    //    min: where we sort from.
    //    max: where we sort to.
    //    val: value we want to insert into the list.
    //    pos: position in array where we want to insort val.
    // Out:
    //    return value: 0 if element is no in the list.
    //                  1 if element is in the list.

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
int look_up_table(bitint input, bitint arr[]) {

    for (int i = 0; i < n_states; i++) {
        if (arr[i] == input) {
            return(i);
        }
    }

}
