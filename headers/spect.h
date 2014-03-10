#ifndef SPECT_H
#define SPECT_H

#include <armadillo>

using namespace arma;

// Recursion.
double continued_fraction(double *a, double *b);
// Exact Diagonalisation.
void dos_mat(double *a_c, double *b_c, mat input);
//void dos_noise(double factor);
void dos_norm(double *a, double *b);
void random_overlap(mat &overlap);
void overlap_matrix(mat &overlap, uint16_t *configs, vector<double> gs_vec);
void gram_schmidt(mat overlap);

#endif /* SPECT_H */
