#ifndef POLY_H
#define POLY_H

void mult_const_array(vector<double> &input, double mult);

void add_two_array(vector<double> &arr_1, vector<double> arr_2);

void assign_array(vector<double> &output, vector<double> input);

void shift_array(vector<double> &input);

void print_array(vector<double> input);

void poly_deriv(vector<double> &poly);

void store_norm(vector<double> &norm, vector<double> p_N, vector<double> p_Nm, vector<double> eig);

void poly_rec(vector<double> &p_nmi, vector<double> &p_N, int p_depth, double a_n[], double b_n[]);

void poly_dos(double *a, double *b, vector<double> overlap);

#endif /* POLY_H */
