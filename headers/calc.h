#ifndef CALC_H
#define CALC_H

#include <vector>
#include <iostream>

using namespace std;

void set_up_globals(vector<string> parsed, vector<double> val);
void set_up_calc(vector<string> parsed);
bool calc_present(vector<string> parsed, string input);
double val_present(vector<string> parsed, vector<double> val, string input);
void calc(char *filename);
void read_input(vector<string> &parsed, vector<double> &val, char *filename);
void input_output();

#endif /* CALC_H */
