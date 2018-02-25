
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "randdp.h"

double MAX(double X, double Y) {
    return (((X) > (Y)) ? (X) : (Y));
}


typedef struct ep_parameter {
    int M;
    double *x;
    double *q;
} ep_parameter_t;


ep_parameter_t buildEPParameters(int M /* = [17:40] */) {
    int x_size = 2 * (1 << 16);
    int q_size = 10;

    ep_parameter_t obj = { .M = M,
            .x = (double*) malloc((x_size) * sizeof(double)),
            .q = (double*) malloc((q_size) * sizeof(double))};

    return obj;
}

void freeEPParameters(ep_parameter_t parameter) {
    free(parameter.q);
    free(parameter.x);
}


/**
 *
 * @param M size parameter (scales exponentially), from [24:40]
 * @param x 2*NK long, where NK = 1 << MK and MK = 16
 * @param q NQ long where NQ = 10
 * @return 0
 */
int embarassinglyParallel(
        const unsigned int M,
        double *x,
        double *q
) {
    if (M < 16) {
        printf("M is below 16, this will break NPB-EP! Try again with M >= 16\n");
    }


    const int MK = 16;
    const long MM = (M - MK);
    const long NN = (1 << MM);
    const long NK = (1 << MK);
    const int NQ = 10;
    const double A = 1220703125.0;
    const double S = 271828183.0;


    double Mops, t1, t2, t3, t4, x1, x2;
    double sx, sy, an, tt, gc;
    int np;
    int i, ik, kk, l, k, nit;
    int k_offset, j;

    double dum[3] = {1.0, 1.0, 1.0};
    char size[16];


    //--------------------------------------------------------------------
    //  Because the size of the problem is too large to store in a 32-bit
    //  integer for some classes, we put it into a string (for printing).
    //  Have to strip off the decimal point put in there by the floating
    //  point print statement (internal file)
    //--------------------------------------------------------------------

    sprintf(size, "%15.0lf", pow(2.0, M + 1));
    j = 14;
    if (size[j] == '.') j--;
    size[j + 1] = '\0';


    //--------------------------------------------------------------------
    //  Compute the number of "batches" of random number pairs generated
    //  per processor. Adjust if the number of processors does not evenly
    //  divide the total number
    //--------------------------------------------------------------------

    np = NN;

    //--------------------------------------------------------------------
    //  Call the random number generator functions and initialize
    //  the x-array to reduce the effects of paging on the timings.
    //  Also, call all mathematical functions that are used. Make
    //  sure these initializations cannot be eliminated as dead code.
    //--------------------------------------------------------------------

    vranlc(0, &dum[0], dum[1], &dum[2]);
    dum[0] = randlc(&dum[1], dum[2]);
    for (i = 0; i < 2 * NK; i++) {
        x[i] = -1.0e99;
    }
    Mops = log(sqrt(fabs(MAX(1.0, 1.0))));


    t1 = A;
    vranlc(0, &t1, A, x);

    //--------------------------------------------------------------------
    //  Compute AN = A ^ (2 * NK) (mod 2^46).
    //--------------------------------------------------------------------

    t1 = A;

    for (i = 0; i < MK + 1; i++) {
        t2 = randlc(&t1, t1);
    }

    an = t1;
    tt = S;
    gc = 0.0;
    sx = 0.0;
    sy = 0.0;

    for (i = 0; i < NQ; i++) {
        q[i] = 0.0;
    }

    //--------------------------------------------------------------------
    //  Each instance of this loop may be performed independently. We compute
    //  the k offsets separately to take into account the fact that some nodes
    //  have more numbers to generate than others
    //--------------------------------------------------------------------

    k_offset = -1;

    for (k = 1; k <= np; k++) {
        kk = k_offset + k;
        t1 = S;
        t2 = an;

        // Find starting seed t1 for this kk.

        for (i = 1; i <= 100; i++) {
            ik = kk / 2;
            if ((2 * ik) != kk) t3 = randlc(&t1, t2);
            if (ik == 0) break;
            t3 = randlc(&t2, t2);
            kk = ik;
        }

        //--------------------------------------------------------------------
        //  Compute uniform pseudorandom numbers.
        //--------------------------------------------------------------------

        vranlc(2 * NK, &t1, A, x);

        //--------------------------------------------------------------------
        //  Compute Gaussian deviates by acceptance-rejection method and
        //  tally counts in concentri//square annuli.  This loop is not
        //  vectorizable.
        //--------------------------------------------------------------------


        for (i = 0; i < NK; i++) {
            x1 = 2.0 * x[2 * i] - 1.0;
            x2 = 2.0 * x[2 * i + 1] - 1.0;
            t1 = x1 * x1 + x2 * x2;
            if (t1 <= 1.0) {
                t2 = sqrt(-2.0 * log(t1) / t1);
                t3 = (x1 * t2);
                t4 = (x2 * t2);
                l = MAX(fabs(t3), fabs(t4));
                q[l] = q[l] + 1.0;
                sx = sx + t3;
                sy = sy + t4;
            }
        }

    }

    for (i = 0; i < NQ; i++) {
        gc = gc + q[i];
    }


    nit = 0;

    return 0;
}


int embarassinglyParallelPacked(ep_parameter_t parameter) {
    return embarassinglyParallel(parameter.M, parameter.x, parameter.q);
}
