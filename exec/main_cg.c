#include <stdio.h>
#include <stdlib.h>

#include "../include/cg.h"

int main(int argc, char *argv[])
{
    int NA1K = 150;
    int NONZER = 15;
    int NITER = 75;
    double SHIFT = 110.0; //from [10.0:1500.0]

    if (argc == 5) {
        NA1K = (int)  strtol(argv[1], NULL, 10);
        NONZER = (int)  strtol(argv[2], NULL, 10);
        NITER = (int) strtol(argv[3], NULL, 10);
        SHIFT = (double)  strtol(argv[4], NULL, 10);
        printf("using NA1K=%d NONZER=%d NITER=%d SHIFT=%e\n", NA1K, NONZER, NITER, SHIFT);
    }

    cg_parameters_t parameters = buildCGParameters(NA1K, NONZER, NITER, SHIFT);

    conjugateGradientPacked(parameters);

    freeCGParameters(parameters);


    return 0;
}
