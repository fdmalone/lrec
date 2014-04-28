#ifndef SPECT_H
#define SPECT_H

#include <armadillo>
#include "const.h"

using namespace arma;

// Recursion.
double continued_fraction(vector<double> *a, vector<double> *b);
// Exact Diagonalisation.
void dos_mat(vector<double> a_c, vector<double> b_c, vector<double> input);
//void dos_noise(vector<double> factor);
void dos_norm(vector<double> a, vector<double> b);
void random_overlap(vector<double> &overlap);
void overlap_matrix(vector<double> &overlap, bitint *configs, vector<double> gs_vec);
void overlap_from_file(vector<double> &overlap, bitint *configs, vector<double> gs_vec);
void gram_schmidt(mat overlap);

#endif /* SPECT_H */
