#ifndef SPECT_H
#define SPECT_H

#include <armadillo>

using namespace arma;

// Recursion.
double continued_fraction(double a[], double b[], double num, double omega);
// Exact Diagonalisation.
void dos_mat(int L, double a_c[], double b_c[], mat input, double it, int step);
//void dos_noise(double factor);
void dos_norm(int its, double omega, double step, int depth, double a[], double b[]);
void overlap_matrix(int size, mat &overlap, uint16_t *configs, vector<double> gs_vec);
void gram_schmidt(mat overlap);

#endif /* SPECT_H */
