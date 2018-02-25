#include <stdio.h>
#include <stdlib.h>

#include "../include/mg.h"

int main(int argc, char *argv[])
{
    //this segfault comes from earlier! not recent changes!


    int PROBLEM_SIZE = 512;
    int NIT = 20;

    if (argc == 3) {
        PROBLEM_SIZE = (int) strtol(argv[1], NULL, 10);
        NIT = (int) strtol(argv[2], NULL, 10);
        printf("using PROBLEM_SIZE=%d NIT=%d\n", PROBLEM_SIZE, NIT);
    }

    mg_parameters_t parameters = buildMGParameters(PROBLEM_SIZE, NIT);

    multiGridPacked(parameters);

    freeMGParameters(parameters);

    return 0;
}
