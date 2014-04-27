#include <iostream>
#include <stdint.h>
#include "calc.h"
#include <time.h>

using namespace std;

int main(int argc, char* argv[]) {

    clock_t t;

    t = clock();

    do_calc(argv[1]);

    t = clock() - t;

    cout << "Time taken: " << (float)t/CLOCKS_PER_SEC << " seconds" << endl;

}
