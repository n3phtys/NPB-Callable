#include <stdio.h>
#include <stdlib.h>

#include "../include/is.h"

int main(int argc, char *argv[])
{
    long TOTAL_KEYS_LOG_2 = 25; //[16;31]
    long MAX_KEY_LOG_2 = 21; //[11;27]
    long NUM_BUCKETS_LOG_2 = 10; //[9;10]

    if (argc == 4) {
        TOTAL_KEYS_LOG_2 = (size_t) strtol(argv[1], NULL, 10);
        MAX_KEY_LOG_2 = (size_t) strtol(argv[2], NULL, 10);
        NUM_BUCKETS_LOG_2 = (size_t) strtol(argv[3], NULL, 10);
        printf("using TOTAL_KEYS_LOG_2=%ld MAX_KEY_LOG_2=%ld NUM_BUCKETS_LOG_2=%ld\n", TOTAL_KEYS_LOG_2, MAX_KEY_LOG_2, NUM_BUCKETS_LOG_2);
    }


    is_parameters_t parameters = buildISParameters(TOTAL_KEYS_LOG_2, MAX_KEY_LOG_2, NUM_BUCKETS_LOG_2);

    integerSortPacked(parameters);

    freeISParameters(parameters);

    return 0;
}
