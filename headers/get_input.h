#ifndef GET_INPUT_H
#define GET_INPUT_H

#include <vector>
#include <iostream>

using namespace std;

void set_up_globals(vector<string> parsed, vector<double> val);
void set_up_calc(vector<string> parsed);
void set_up_system(char *filename);
bool calc_present(vector<string> parsed, string input);
double val_present(vector<string> parsed, vector<double> val, string input);
void read_input(vector<string> &parsed, vector<double> &val, char *filename);
void input_output();

#endif /* GET_INPUT_H */
