#include <stdio.h>
#include <stdlib.h>

#include "../include/ep.h"

int main(int argc, char *argv[])
{
    //[24:40]
    int M = 17;

    if (argc == 2) {
        M = (int) strtol(argv[1], NULL, 10);
        printf("using M=%d\n", M);
    }


    ep_parameter_t parameters = buildEPParameters(M);

    embarassinglyParallelPacked(parameters);

    freeEPParameters(parameters);

    return 0;
}
