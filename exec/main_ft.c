#include <stdio.h>
#include <stdlib.h>

#include "../include/ft.h"

int main(int argc, char *argv[])
{
    int NX = 64;
    int NY = 64;
    int NZ = 64;
    int NITER = 6;

    if (argc == 5) {
        NX = (int) strtol(argv[1], NULL, 10);
        NY = (int) strtol(argv[2], NULL, 10);
        NZ = (int) strtol(argv[3], NULL, 10);
        NITER = (int) strtol(argv[4], NULL, 10);
        printf("using NX=%d NY=%d NZ=%d NITER=%d", NX, NY, NZ, NITER);
    }

    ft_parameters_t parameters = buildFTParameters(NX, NY, NZ, NITER);

    fourierTransformationPacked(parameters);

    freeFTParameters(parameters);
    

    return 0;
}
