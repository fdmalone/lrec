#include <stdint.h>
#include <vector>
// Sorting;
int permute(int a, int b, double &sign);
int permute_norm(int a, int b, double &sign);
int insertion_sort(int array[], int lenght);
int insertion_rank(uint16_t array[], uint16_t rank[], int length);
int binary_search(vector<uint16_t> &a, int min, int max, uint16_t val, int &pos);
int look_up_table(uint16_t input, uint16_t arr[]);
