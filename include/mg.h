#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "randdp.h"


const int ONE = 1;
const int lb = 1;


typedef struct mg_parameters {

    int PROBLEM_SIZE;
    int NIT;

    int NX_DEFAULT;
    int NY_DEFAULT;
    int NZ_DEFAULT;
    int LM;
    int LT_DEFAULT;
    int NDIM1;
    int NDIM2;
    int NDIM3;
    int lt;


// actual dimension including ghost cells for communications
    long NM;

    long M;

// size of rhs array
    long NV;

// size of residual array
    long NR;

// maximum number of levels
    long MAXLEVEL;

    int NIT_DEFAULT;


    int *nx; //MAXLEVEL+1 long
    int *ny; //MAXLEVEL+1 long
    int *nz; //MAXLEVEL+1 long
    int *m1; //MAXLEVEL+1 long
    int *m2; //MAXLEVEL+1 long
    int *m3; //MAXLEVEL+1 long
    int *ir; //MAXLEVEL+1 long

    double *u; //NR long
    double *v; //NR long
    double *r; //NR long

    int is1, is2, is3, ie1, ie2, ie3;
} mg_parameters_t;


/* integer log base two. Return error is argument isn't
 * a power of two or is less than or equal to zero
 */

int ilog2(int i) {
    int log2;
    int exp2 = 1;
    if (i <= 0) return (-1);

    for (log2 = 0; log2 < 30; log2++) {
        if (exp2 == i) return (log2);
        if (exp2 > i) break;
        exp2 *= 2;
    }
    return (-1);
}


mg_parameters_t buildMGParameters(int problem_size /* = [32:2048] */, int nit /* = [4:50] */ ) {

    const int ndim = ilog2(problem_size);

    if (1 << ndim != problem_size) {
        printf("problem_size is illegal, rest: %d\n", problem_size - (1 << ndim));
    }

    const int n_default = problem_size;
    const int NX_DEFAULT = n_default;
    const int NY_DEFAULT = n_default;
    const int NZ_DEFAULT = n_default;
    const int LM = ndim;
    const int LT_DEFAULT = ndim;
    const int NDIM1 = ndim;
    const int NDIM2 = ndim;
    const int NDIM3 = ndim;
    const int lt = ndim;

    const int NIT_DEFAULT = nit;


// actual dimension including ghost cells for communications
    const long NM = (2 + (1 << LM));

    const long M = NM + 1;

// size of rhs array
    const long NV = (ONE * (2 + (1 << NDIM1)) * (2 + (1 << NDIM2)) * (2 + (1 << NDIM3)));

// size of residual array
    const long NR = (((NV + NM * NM + 5 * NM + 7 * LM + 6) / 7) * 8);

// maximum number of levels
    const long MAXLEVEL = (LT_DEFAULT + 1);


    mg_parameters_t obj = {
            .PROBLEM_SIZE = problem_size,
            .NIT = nit,
            .NX_DEFAULT = NX_DEFAULT,
            .NY_DEFAULT = NY_DEFAULT,
            .NZ_DEFAULT = NZ_DEFAULT,
            .LM = LM,
            .LT_DEFAULT = LT_DEFAULT,
            .NDIM1 = NDIM1,
            .NDIM2 = NDIM2,
            .NDIM3 = NDIM3,
            .lt = lt,
            .NM = NM,
            .M = M,
            .NV = NV,
            .NR = NR,
            .MAXLEVEL = MAXLEVEL,
            .NIT_DEFAULT = NIT_DEFAULT,
            .nx = (int *) malloc((MAXLEVEL + 1) * sizeof(int)),
            .ny = (int *) malloc((MAXLEVEL + 1) * sizeof(int)),
            .nz = (int *) malloc((MAXLEVEL + 1) * sizeof(int)),
            .m1 = (int *) malloc((MAXLEVEL + 1) * sizeof(int)),
            .m2 = (int *) malloc((MAXLEVEL + 1) * sizeof(int)),
            .m3 = (int *) malloc((MAXLEVEL + 1) * sizeof(int)),
            .ir = (int *) malloc((MAXLEVEL + 1) * sizeof(int)),
            .u = (double *) malloc((NR) * sizeof(double)),
            .v = (double *) malloc((NR) * sizeof(double)),
            .r = (double *) malloc((NR) * sizeof(double)),
            .is1 = 0,
            .is2 = 0,
            .is3 = 0,
            .ie1 = 0,
            .ie2 = 0,
            .ie3 = 0
    };


    return obj;

}

void freeMGParameters(mg_parameters_t parameter) {
    free(parameter.nx);
    free(parameter.ny);
    free(parameter.nz);
    free(parameter.m1);
    free(parameter.m2);
    free(parameter.m3);
    free(parameter.ir);
    free(parameter.u);
    free(parameter.v);
    free(parameter.r);
}


//if min = 0 and 0 <= v < max, this will result in v1. Will wrap number between min and max if not.
static long wrap(long v1, long min, long max) {
    const long v2 = v1 < 0 ? (-1 * v1) : v1;
    const long v3 = v2 - min;
    const long v4 = v3 % (max - min);
    const long v5 = v4 + min;
    return v5;
}



//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB MG code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

//---------------------------------------------------------------------
//  program mg
//---------------------------------------------------------------------




void setup(int *n1, int *n2, int *n3, const int lt, //
           const int MAXLEVEL,
           int *nx,//
           int *ny,//
           int *nz, //
           int *m1,//
           int *m2,//
           int *m3,//
           int *ir, //
           int *is1, int *is2, int *is3, int *ie1, int *ie2, int *ie3);

void mg3P(double u[], double v[], double r[],
          double a[4], double c[4], int n1, int n2, int n3, const int lt,
          const int M,
          int *m1,
          int *m2,
          int *m3,
          int *ir);

void psinv(double* r, double* u, int n1, int n2, int n3,
           double c[4], int k, const int M);

void resid(double *u, double *v, double *r, int n1, int n2, int n3,
           double a[4], int k, const int M);

void rprj3(double* r, int m1k, int m2k, int m3k,
           double* s, int m1j, int m2j, int m3j, int k, const int M);

void interp(double *z, int mm1, int mm2, int mm3,
            double *u, int n1, int n2, int n3, int k, const int M);

void norm2u3(double* r, int n1, int n2, int n3,
             double *rnm2, double *rnmu,
             int nx, int ny, int nz);

void comm3(double *u, int n1, int n2, int n3, int kk);

void
zran3(double *z, int n1, int n2, int n3, int nx, int ny, const int k, int *is1, int *is2, int *is3, int *ie1, int *ie2,
      int *ie3);

double power(double a, int n);

void bubble(double ten[][2], int j1[][2], int j2[][2], int j3[][2],
            int m, int ind);

void zero3(double *z, int n1, int n2, int n3);


//-------------------------------------------------------------------------c
// These arrays are in common because they are quite large
// and probably shouldn't be allocated on the stack. They
// are always passed as subroutine args. 
//-------------------------------------------------------------------------c


void multiGrid(
        int problem_size,
        int nit,
        int *nx, //
        int *ny, //
        int *nz, //
        int *m1, //
        int *m2, //
        int *m3, //
        int *ir, //
        double *u, //
        double *v, //
        double *r, //
        int *is1, int *is2, int *is3, int *ie1, int *ie2, int *ie3 //
) {

    const int ndim = ilog2(problem_size);
    const int n_default = problem_size;
    const int NX_DEFAULT = n_default;
    const int NY_DEFAULT = n_default;
    const int NZ_DEFAULT = n_default;
    const int LM = ndim;
    const int LT_DEFAULT = ndim;
    const int NDIM1 = ndim;
    const int NDIM2 = ndim;
    const int NDIM3 = ndim;
    const int lt = ndim;

    const long NM = (2 + (1 << LM));
    const long M = NM + 1;
    const long NV = (ONE * (2 + (1 << NDIM1)) * (2 + (1 << NDIM2)) * (2 + (1 << NDIM3)));
    const long NR = (((NV + NM * NM + 5 * NM + 7 * LM + 6) / 7) * 8);
    const long MAXLEVEL = (LT_DEFAULT + 1);






    //-------------------------------------------------------------------------c
    // k is the current level. It is passed down through subroutine args
    // and is NOT global. it is the current iteration
    //-------------------------------------------------------------------------c
    int it;


    double a[4], c[4];

    double rnm2, rnmu, old2, oldu, epsilon;
    int n1, n2, n3;
    double nn;





    nx[lt] = NX_DEFAULT;
    ny[lt] = NY_DEFAULT;
    nz[lt] = NZ_DEFAULT;


    a[0] = -8.0 / 3.0;
    a[1] = 0.0;
    a[2] = 1.0 / 6.0;
    a[3] = 1.0 / 12.0;

    //---------------------------------------------------------------------
    // Coefficients for the S(a) smoother
    //---------------------------------------------------------------------
    c[0] = -3.0 / 8.0;
    c[1] = +1.0 / 32.0;
    c[2] = -1.0 / 64.0;
    c[3] = 0.0;

    const int k = lt;

    setup(&n1, &n2, &n3, lt, MAXLEVEL, nx, ny, nz, m1, m2, m3, ir, is1, is2, is3, ie1, ie2, ie3);

    printf("n1 = %d n2 = %d n3 = %d , NR = %ld , product = %ld\n", n1, n2, n3, NR, 1L * n1 * n2 * n3);


    zero3(u, n1, n2, n3);
    zran3(v, n1, n2, n3, nx[lt], ny[lt], k, is1, is2, is3, ie1, ie2, ie3);

    norm2u3(v, n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt],
            nz[lt]);


    resid(u, v, r, n1, n2, n3, a, k, M);
    norm2u3(r, n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt],
            nz[lt]);
    old2 = rnm2;
    oldu = rnmu;

    //---------------------------------------------------------------------
    // One iteration for startup
    //---------------------------------------------------------------------
    mg3P(u, v, r, a, c, n1, n2, n3, lt, M, m1, m2, m3, ir);
    resid(u, v, r, n1, n2, n3, a, k, M);
    setup(&n1, &n2, &n3, lt, MAXLEVEL, nx, ny, nz, m1, m2, m3, ir, is1, is2, is3, ie1, ie2, ie3);
    zero3(u, n1, n2, n3);
    zran3(v, n1, n2, n3, nx[lt], ny[lt], k, is1, is2, is3, ie1, ie2, ie3);


    resid(u, v, r, n1, n2, n3, a, k, M);

    norm2u3(r, n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt], nz[lt]);
    old2 = rnm2;
    oldu = rnmu;

    for (it = 1; it <= nit; it++) {
        mg3P(u, v, r, a, c, n1, n2, n3, lt, M, m1, m2, m3, ir);
        resid(u, v, r, n1, n2, n3, a, k, M);
    }

    norm2u3(r, n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt], nz[lt]);


    epsilon = 1.0e-8;


    nn = 1.0 * nx[lt] * ny[lt] * nz[lt];


}


void multiGridPacked(mg_parameters_t parameters) {
    multiGrid(parameters.PROBLEM_SIZE,
                     parameters.NIT,
                     parameters.nx,
                     parameters.ny,
                     parameters.nz,
                     parameters.m1,
                     parameters.m2,
                     parameters.m3,
                     parameters.ir,
                     parameters.u,
                     parameters.v,
                     parameters.r,
                     &parameters.is1, &parameters.is2, &parameters.is3, &parameters.ie1, &parameters.ie2,
                     &parameters.ie3
    );
}

void setup(int *n1, int *n2, int *n3,
           const int lt, //
           const int MAXLEVEL,
           int *nx,//
           int *ny,//
           int *nz, //
           int *m1,//
           int *m2,//
           int *m3,//
           int *ir, //
           int *is1, int *is2, int *is3, int *ie1, int *ie2, int *ie3) {
    int k, j;

    int ax, mi[MAXLEVEL + 1][3];
    int ng[MAXLEVEL + 1][3];

    ng[lt][0] = nx[lt];
    ng[lt][1] = ny[lt];
    ng[lt][2] = nz[lt];
    for (k = lt - 1; k >= 1; k--) {
        for (ax = 0; ax < 3; ax++) {
            ng[k][ax] = ng[k + 1][ax] / 2;
        }
    }
    for (k = lt; k >= 1; k--) {
        nx[k] = ng[k][0];
        ny[k] = ng[k][1];
        nz[k] = ng[k][2];
    }

    for (k = lt; k >= 1; k--) {
        for (ax = 0; ax < 3; ax++) {
            mi[k][ax] = 2 + ng[k][ax];
        }

        m1[k] = mi[k][0];
        m2[k] = mi[k][1];
        m3[k] = mi[k][2];
    }

    k = lt;
    *is1 = 2 + ng[k][0] - ng[lt][0];
    *ie1 = 1 + ng[k][0];
    *n1 = 3 + *ie1 - *is1;
    *is2 = 2 + ng[k][1] - ng[lt][1];
    *ie2 = 1 + ng[k][1];
    *n2 = 3 + *ie2 - *is2;
    *is3 = 2 + ng[k][2] - ng[lt][2];
    *ie3 = 1 + ng[k][2];
    *n3 = 3 + *ie3 - *is3;

    ir[lt] = 0;
    for (j = lt - 1; j >= 1; j--) {
        ir[j] =
                ir[j + 1] + ONE * m1[j + 1] * m2[j + 1] * m3[j + 1];
    }

}


//---------------------------------------------------------------------
// multigrid V-cycle routine
//---------------------------------------------------------------------
void mg3P(double u[], double v[], double r[],
          double a[4], double c[4], int n1, int n2, int n3,
          const int lt,
          const int M,
          int *m1,
          int *m2,
          int *m3,
          int *ir
) {
    int j, k;


    //---------------------------------------------------------------------
    // down cycle.
    // restrict the residual from the find grid to the coarse
    //---------------------------------------------------------------------
    for (k = lt; k >= lb + 1; k--) {
        j = k - 1;
        rprj3(&r[ir[k]], m1[k], m2[k], m3[k],
              &r[ir[j]], m1[j], m2[j], m3[j], k, M);
    }

    k = lb;
    //---------------------------------------------------------------------
    // compute an approximate solution on the coarsest grid
    //---------------------------------------------------------------------
    zero3(&u[ir[k]], m1[k], m2[k], m3[k]);
    psinv(&r[ir[k]], &u[ir[k]], m1[k], m2[k], m3[k], c, k,
          M);

    for (k = lb + 1; k <= lt - 1; k++) {
        j = k - 1;

        //---------------------------------------------------------------------
        // prolongate from level k-1  to k
        //---------------------------------------------------------------------
        zero3(&u[ir[k]], m1[k], m2[k], m3[k]);
        interp(&u[ir[j]], m1[j], m2[j], m3[j], &u[ir[k]],
               m1[k], m2[k], m3[k], k, M);

        //---------------------------------------------------------------------
        // compute residual for level k
        //---------------------------------------------------------------------
        resid(&u[ir[k]], &r[ir[k]], &r[ir[k]], m1[k], m2[k],
              m3[k], a, k, M);

        //---------------------------------------------------------------------
        // apply smoother
        //---------------------------------------------------------------------
        psinv(&r[ir[k]], &u[ir[k]], m1[k], m2[k], m3[k], c, k,
              M);
    }

    j = lt - 1;
    k = lt;
    interp(&u[ir[j]], m1[j], m2[j], m3[j], u, n1, n2, n3, k, M);
    resid(u, v, r, n1, n2, n3, a, k, M);
    psinv(r, u, n1, n2, n3, c, k, M);
}


//---------------------------------------------------------------------
// psinv applies an approximate inverse as smoother:  u = u + Cr
//
// This  implementation costs  15A + 4M per result, where
// A and M denote the costs of Addition and Multiplication.  
// Presuming coefficient c(3) is zero (the NPB assumes this,
// but it is thus not a general case), 2A + 1M may be eliminated,
// resulting in 13A + 3M.
// Note that this vectorizes, and is also fine for cache 
// based machines.  
//---------------------------------------------------------------------
void psinv(double* r, double* u, int n1, int n2, int n3,
           double c[4], int k, const int M) {

    int i3, i2, i1;

    double r1[M], r2[M];


    for (i3 = 1; i3 < n3 - 1; i3++) {
        for (i2 = 1; i2 < n2 - 1; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
                r1[i1] = r[(i3*n2*n1)+((i2 - 1)*n1)+i1] + r[(i3*n2*n1)+((i2+1)*n1)+i1]
                         + r[((i3 - 1)*n2*n1)+(i2*n1)+i1] + r[((i3 + 1)*n2*n1)+(i2*n1)+i1];
                r2[i1] = r[((i3 - 1)*n2*n1)+((i2-1)*n1)+i1] + r[((i3 - 1)*n2*n1)+((i2+1)*n1)+i1]
                         + r[((i3 + 1)*n2*n1)+((i2-1)*n1)+i1] + r[((i3 + 1)*n2*n1)+((i2+1)*n1)+i1];
            }
            for (i1 = 1; i1 < n1 - 1; i1++) {
                u[(i3*n2*n1)+(i2*n1)+i1] = u[(i3*n2*n1)+(i2*n1)+i1]
                                + c[0] * r[(i3*n2*n1)+(i2*n1)+i1]
                                + c[1] * (r[(i3*n2*n1)+(i2*n1)+(i1-1)] + r[(i3*n2*n1)+(i2*n1)+(i1+1)]
                                          + r1[i1])
                                + c[2] * (r2[i1] + r1[i1 - 1] + r1[i1 + 1]);
                //--------------------------------------------------------------------
                // Assume c[3] = 0    (Enable line below if c[3] not= 0)
                //--------------------------------------------------------------------
                //            + c[3] * ( r2[i1-1] + r2[i1+1] )
                //--------------------------------------------------------------------
            }
        }
    }

    //---------------------------------------------------------------------
    // exchange boundary points
    //---------------------------------------------------------------------
    comm3(u, n1, n2, n3, k);

}


//---------------------------------------------------------------------
// resid computes the residual:  r = v - Au
//
// This  implementation costs  15A + 4M per result, where
// A and M denote the costs of Addition (or Subtraction) and 
// Multiplication, respectively. 
// Presuming coefficient a(1) is zero (the NPB assumes this,
// but it is thus not a general case), 3A + 1M may be eliminated,
// resulting in 12A + 3M.
// Note that this vectorizes, and is also fine for cache 
// based machines.  
//---------------------------------------------------------------------
void resid(double *u, double *v, double* r, int n1, int n2, int n3,
           double a[4], int k, const int M) {

    int i3, i2, i1;
    double u1[M], u2[M];


    for (i3 = 1; i3 < n3 - 1; i3++) {
        for (i2 = 1; i2 < n2 - 1; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
                u1[i1] = u[(i3*n2*n1)+((i2-1)*n1)+i1] + u[(i3*n2*n1)+((i2+1)*n1)+i1]
                         + u[((i3 - 1)*n2*n1)+(i2*n1)+i1] + u[((i3 + 1)*n2*n1)+(i2*n1)+i1];
                u2[i1] = u[((i3 - 1)*n2*n1)+((i2-1)*n1)+i1] + u[((i3 - 1)*n2*n1)+((i2+1)*n1)+i1]
                         + u[((i3 + 1)*n2*n1)+((i2-1)*n1)+i1] + u[((i3 + 1)*n2*n1)+((i2+1)*n1)+i1];
            }
            for (i1 = 1; i1 < n1 - 1; i1++) {
                r[(i3*n2*n1)+(i2*n1)+i1] = v[(i3*n2*n1)+(i2*n1)+i1]
                                - a[0] * u[(i3*n2*n1)+(i2*n1)+i1]
                                //-------------------------------------------------------------------
                                //  Assume a[1] = 0      (Enable 2 lines below if a[1] not= 0)
                                //-------------------------------------------------------------------
                                //            - a[1] * ( u[i3][i2][i1-1] + u[i3][i2][i1+1]
                                //                     + u1[i1] )
                                //-------------------------------------------------------------------
                                - a[2] * (u2[i1] + u1[i1 - 1] + u1[i1 + 1])
                                - a[3] * (u2[i1 - 1] + u2[i1 + 1]);
            }
        }
    }

    //---------------------------------------------------------------------
    // exchange boundary data
    //---------------------------------------------------------------------
    comm3(r, n1, n2, n3, k);

}


//---------------------------------------------------------------------
// rprj3 projects onto the next coarser grid, 
// using a trilinear Finite Element projection:  s = r' = P r
//     
// This  implementation costs  20A + 4M per result, where
// A and M denote the costs of Addition and Multiplication.  
// Note that this vectorizes, and is also fine for cache 
// based machines.  
//---------------------------------------------------------------------
void rprj3(double* r, int m1k, int m2k, int m3k,
           double* s, int m1j, int m2j, int m3j, int k, const int M) {

    int j3, j2, j1, i3, i2, i1, d1, d2, d3, j;

    double x1[M], y1[M];
    double x2, y2;


    if (m1k == 3) {
        d1 = 2;
    } else {
        d1 = 1;
    }

    if (m2k == 3) {
        d2 = 2;
    } else {
        d2 = 1;
    }

    if (m3k == 3) {
        d3 = 2;
    } else {
        d3 = 1;
    }

    int n2 = m2j;
    int n1 = m1j;

    for (j3 = 1; j3 < m3j - 1; j3++) {
        i3 = 2 * j3 - d3;
        for (j2 = 1; j2 < m2j - 1; j2++) {
            i2 = 2 * j2 - d2;

            for (j1 = 1; j1 < m1j; j1++) {
                i1 = 2 * j1 - d1;
                x1[i1] = r[((i3 + 1)*n2*n1)+(i2*n1)+i1] + r[((i3 + 1)*n2*n1)+((i2+2)*n1)+i1]
                         + r[(i3*n2*n1)+((i2+1)*n1)+i1] + r[((i3 + 2)*n2*n1)+((i2+1)*n1)+i1];
                y1[i1] = r[(i3*n2*n1)+(i2*n1)+i1] + r[((i3 + 2)*n2*n1)+(i2*n1)+i1]
                         + r[(i3*n2*n1)+((i2+2)*n1)+i1] + r[((i3 + 2)*n2*n1)+((i2+2)*n1)+i1];
            }

            for (j1 = 1; j1 < m1j - 1; j1++) {
                i1 = 2 * j1 - d1;
                y2 = r[(i3*n2*n1)+(i2*n1)+(i1+1)] + r[((i3 + 2)*n2*n1)+(i2*n1)+(i1+1)]
                     + r[(i3*n2*n1)+((i2+2)*n1)+(i1+1)] + r[((i3 + 2)*n2*n1)+((i2+2)*n1)+(i1+1)];
                x2 = r[((i3 + 1)*n2*n1)+(i2*n1)+(i1+1)] + r[((i3 + 1)*n2*n1)+((i2+2)*n1)+(i1+1)]
                     + r[(i3*n2*n1)+((i2+1)*n1)+(i1+1)] + r[((i3 + 2)*n2*n1)+((i2+1)*n1)+(i1+1)];
                s[(j3*m2j*m1j)+(j2*m1j)+j1] =
                        0.5 * r[((i3 + 1)*n2*n1)+((i2+1)*n1)+(i1+1)]
                        + 0.25 * (r[((i3 + 1)*n2*n1)+((i2+1)*n1)+i1] + r[((i3 + 1)*n2*n1)+((i2+1)*n1)+(i1+2)] + x2)
                        + 0.125 * (x1[i1] + x1[i1 + 2] + y2)
                        + 0.0625 * (y1[i1] + y1[i1 + 2]);
            }
        }
    }


    j = k - 1;
    comm3(s, m1j, m2j, m3j, j);

}


//---------------------------------------------------------------------
// interp adds the trilinear interpolation of the correction
// from the coarser grid to the current approximation:  u = u + Qu'
//     
// Observe that this  implementation costs  16A + 4M, where
// A and M denote the costs of Addition and Multiplication.  
// Note that this vectorizes, and is also fine for cache 
// based machines.  Vector machines may get slightly better 
// performance however, with 8 separate "do i1" loops, rather than 4.
//---------------------------------------------------------------------
void interp(double *z, int mm1, int mm2, int mm3,
            double *u, int n1, int n2, int n3, int k, const int M) {


    int i3, i2, i1, d1, d2, d3, t1, t2, t3;

    // note that m = 1037 in globals.h but for this only need to be
    // 535 to handle up to 1024^3
    //      integer m
    //      parameter( m=535 )
    double z1[M], z2[M], z3[M];


    if (n1 != 3 && n2 != 3 && n3 != 3) {
        for (i3 = 0; i3 < mm3 - 1; i3++) {
            for (i2 = 0; i2 < mm2 - 1; i2++) {
                for (i1 = 0; i1 < mm1; i1++) {
                    z1[i1] = z[(i3*mm2*mm1)+((i2+1)*mm1)+i1] + z[(i3*mm2*mm1)+(i2*mm1)+i1];
                    z2[i1] = z[((i3 + 1)*mm2*mm1)+(i2*mm1)+i1] + z[(i3*mm2*mm1)+(i2*mm1)+i1];
                    z3[i1] = z[((i3 + 1)*mm2*mm1)+((i2+1)*mm1)+i1] + z[((i3 + 1)*mm2*mm1)+(i2*mm1)+i1] + z1[i1];
                }

                for (i1 = 0; i1 < mm1 - 1; i1++) {
                    u[((2 * i3)*n2*n1)+((2*i2)*n1)+(2*i1)] = u[((2 * i3)*n2*n1)+((2*i2)*n1)+(2*i1)]
                                                + z[(i3*mm2*mm1)+(i2*mm1)+i1];
                    u[((2 * i3)*n2*n1)+((2*i2)*n1)+(2 * i1 + 1)] = u[((2 * i3)*n2*n1)+((2*i2)*n1)+(2*i1 +1)]
                                                    + 0.5 * (z[(i3*mm2*mm1)+(i2*mm1)+(i1+1)] + z[(i3*mm2*mm1)+(i2*mm1)+i1]);
                }
                for (i1 = 0; i1 < mm1 - 1; i1++) {
                    u[((2 * i3)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1)] = u[((2 * i3)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1)]
                                                    + 0.5 * z1[i1];
                    u[((2 * i3)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1 + 1)] = u[((2 * i3)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1 + 1)]
                                                        + 0.25 * (z1[i1] + z1[i1 + 1]);
                }
                for (i1 = 0; i1 < mm1 - 1; i1++) {
                    u[((2 * i3 + 1)*n2*n1)+((2 * i2)*n1)+(2 * i1)] = u[((2 * i3 + 1)*n2*n1)+((2 * i2)*n1)+(2 * i1)]
                                                    + 0.5 * z2[i1];
                    u[((2 * i3 + 1)*n2*n1)+((2 * i2)*n1)+(2 * i1 + 1)] = u[((2 * i3 + 1)*n2*n1)+((2 * i2)*n1)+(2 * i1 + 1)]
                                                        + 0.25 * (z2[i1] + z2[i1 + 1]);
                }
                for (i1 = 0; i1 < mm1 - 1; i1++) {
                    u[((2 * i3 + 1)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1)] = u[((2 * i3 + 1)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1)]
                                                        + 0.25 * z3[i1];
                    u[((2 * i3 + 1)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1 + 1)] = u[((2 * i3 + 1)*n2*n1)+((2 * i2 + 1)*n1)+(2 * i1 + 1)]
                                                            + 0.125 * (z3[i1] + z3[i1 + 1]);
                }
            }
        }
    } else {
        if (n1 == 3) {
            d1 = 2;
            t1 = 1;
        } else {
            d1 = 1;
            t1 = 0;
        }

        if (n2 == 3) {
            d2 = 2;
            t2 = 1;
        } else {
            d2 = 1;
            t2 = 0;
        }

        if (n3 == 3) {
            d3 = 2;
            t3 = 1;
        } else {
            d3 = 1;
            t3 = 0;
        }

        for (i3 = d3; i3 <= mm3 - 1; i3++) {
            for (i2 = d2; i2 <= mm2 - 1; i2++) {
                for (i1 = d1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - d1 - 1)] =
                            u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - d1 - 1)]
                            + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)];
                }
                for (i1 = 1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - t1 - 1)] =
                            u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - t1 - 1)]
                            + 0.5 * (z[((i3 - 1)*mm2*mm1)+((i2 - 1)+mm1)+i1] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1) +(i1 - 1)]);
                }
            }
            for (i2 = 1; i2 <= mm2 - 1; i2++) {
                for (i1 = d1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1) + (2 * i1 - d1 - 1)] =
                            u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1)+( 2 * i1 - d1 - 1)]
                            + 0.5 * (z[((i3 - 1)*mm2*mm1)+(i2*mm1)+(i1 - 1)] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1) + (i1 - 1)]);
                }
                for (i1 = 1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1)+(2 * i1 - t1 - 1)] =
                            u[((2 * i3 - d3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1)+(2 * i1 - t1 - 1)]
                            + 0.25 * (z[((i3 - 1)*mm2*mm1)+(i2*mm1)+i1] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+i1]
                                      + z[((i3 - 1)*mm2*mm1)+(i2*mm1)+(i1 - 1)] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]);
                }
            }
        }

        for (i3 = 1; i3 <= mm3 - 1; i3++) {
            for (i2 = d2; i2 <= mm2 - 1; i2++) {
                for (i1 = d1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - d1 - 1)] =
                            u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - d1 - 1)]
                            + 0.5 * (z[(i3*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]);
                }
                for (i1 = 1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - t1 - 1)] =
                            u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - d2 - 1)*n1)+(2 * i1 - t1 - 1)]
                            + 0.25 * (z[(i3*mm2*mm1)+((i2 - 1)*mm1)+i1] + z[(i3*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]
                                      + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+i1] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]);
                }
            }
            for (i2 = 1; i2 <= mm2 - 1; i2++) {
                for (i1 = d1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1)+(2 * i1 - d1 - 1)] =
                            u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1)+(2 * i1 - d1 - 1)]
                            + 0.25 * (z[(i3*mm2*mm1)+(i2*mm1)+(i1 - 1)] + z[(i3*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]
                                      + z[((i3 - 1)*mm2*mm1)+(i2*mm1)+(i1 - 1)] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]);
                }
                for (i1 = 1; i1 <= mm1 - 1; i1++) {
                    u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1)+(2 * i1 - t1 - 1)] =
                            u[((2 * i3 - t3 - 1)*n2*n1)+((2 * i2 - t2 - 1)*n1)+(2 * i1 - t1 - 1)]
                            + 0.125 * (z[(i3*mm2*mm1)+(i2*mm1)+i1] + z[(i3*mm2*mm1)+((i2 - 1)*mm1) +i1]
                                       + z[(i3*mm2*mm1)+(i2*mm1)+(i1 - 1)] + z[(i3*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]
                                       + z[((i3 - 1)*mm2*mm1)+(i2*mm1)+i1] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+ i1]
                                       + z[((i3 - 1)*mm2*mm1)+(i2*mm1)+(i1 - 1)] + z[((i3 - 1)*mm2*mm1)+((i2 - 1)*mm1)+(i1 - 1)]);
                }
            }
        }

    }


}


//---------------------------------------------------------------------
// norm2u3 evaluates approximations to the L2 norm and the
// uniform (or L-infinity or Chebyshev) norm, under the
// assumption that the boundaries are periodic or zero.  Add the
// boundaries in with half weight (quarter weight on the edges
// and eighth weight at the corners) for inhomogeneous boundaries.
//---------------------------------------------------------------------
void norm2u3(double* r, int n1, int n2, int n3,
             double *rnm2, double *rnmu,
             int nx, int ny, int nz) {

    double s, a;
    int i3, i2, i1;

    double dn;


    dn = 1.0 * nx * ny * nz;

    s = 0.0;
    *rnmu = 0.0;
    for (i3 = 1; i3 < n3 - 1; i3++) {
        for (i2 = 1; i2 < n2 - 1; i2++) {
            for (i1 = 1; i1 < n1 - 1; i1++) {
                s = s + pow(r[(i3*n2*n1)+(i2*n1)+i1], 2.0);
                a = fabs(r[(i3*n2*n1)+(i2*n1)+i1]);
                if (a > *rnmu) *rnmu = a;
            }
        }
    }

    *rnm2 = sqrt(s / dn);
}


//---------------------------------------------------------------------
// report on norm
//---------------------------------------------------------------------
void rep_nrm(double* u, int n1, int n2, int n3, char *title, int kk, int *nx, int *ny, int *nz) {
    double rnm2, rnmu;

    norm2u3(u, n1, n2, n3, &rnm2, &rnmu, nx[kk], ny[kk], nz[kk]);
}


//---------------------------------------------------------------------
// comm3 organizes the communication on all borders 
//---------------------------------------------------------------------
void comm3(double *u, int n1, int n2, int n3, int kk) {

    for (int i3 = 1; i3 < n3 - 1; i3++) {
        for (int i2 = 1; i2 < n2 - 1; i2++) {

            u[(i3 * n1 * n2) + (i2 * n1) + 0] = u[i3 * n1 * n2 + i2 * n1 + n1 - 2];
            u[i3 * n1 * n2 + i2 * n1 + n1 - 1] = u[(i3 * n1 * n2) + (i2 * n1) + 1];
        }
    }

    for (int i3 = 1; i3 < n3 - 1; i3++) {
        for (int i1 = 0; i1 < n1; i1++) {
            u[(i3 * n1 * n2) + (0 * n1) + i1] = u[(i3 * n2 * n1) + ((n2 -2) * n1) + i1];
            u[(i3 * n1 * n2) + ((n2-1) * n1) + i1] = u[(i3 * n2 * n1) + (1 * n1) + i1];
        }
    }

    for (int i2 = 0; i2 < n2; i2++) {
        for (int i1 = 0; i1 < n1; i1++) {
            u[(0 * n2 * n1) + (i2 * n1) + i1] = u[((n3 - 2)*n2*n1)+(i2*n1)+i1];
            u[((n3 - 1)*n2*n1)+(i2*n1)+i1] = u[(1*n2*n1)+(i2*n1)+i1];
        }
    }

}

#define  MMCONST 10
//---------------------------------------------------------------------
// zran3  loads +1 at ten randomly chosen points,
// loads -1 at a different ten random points,
// and zero elsewhere.
//---------------------------------------------------------------------
void
zran3(double *z, int n1, int n2, int n3, int nx, int ny, const int k, int *is1, int *is2, int *is3, int *ie1, int *ie2,
      int *ie3) {

    int i0;

    int i1, i2, i3, e1, e2, e3;
    double xx, x0, x1, a1, a2, ai;

    int i;

    const int mm = MMCONST;
    const double a = pow(5.0, 13.0);
    const double x = 314159265.0;
    double ten[MMCONST][2];
    double best;
    int j1[MMCONST][2];
    int j2[MMCONST][2];
    int j3[MMCONST][2];
    int jg[4][MMCONST][2];

    double rdummy;

    a1 = power(a, nx);
    a2 = power(a, nx * ny);

    zero3(z, n1, n2, n3);

    i = *is1 - 2 + nx * (*is2 - 2 + ny * (*is3 - 2));

    ai = power(a, i);
    const int d1 = *ie1 - *is1 + 1;
    e1 = *ie1 - *is1 + 2;
    e2 = *ie2 - *is2 + 2;
    e3 = *ie3 - *is3 + 2;
    x0 = x;
    rdummy = randlc(&x0, ai);

    for (i3 = 1; i3 < e3; i3++) {
        x1 = x0;
        for (i2 = 1; i2 < e2; i2++) {
            xx = x1;
            const int idx = (n2 * n1 * i3) + (n1 * i2) + (1);
            vranlc(d1, &xx, a, z + idx);
            rdummy = randlc(&x1, a1);
        }
        rdummy = randlc(&x0, a2);
    }

    //---------------------------------------------------------------------
    // comm3(z,n1,n2,n3);
    // showall(z,n1,n2,n3);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // each processor looks for twenty candidates
    //---------------------------------------------------------------------
    for (int i = 0; i < mm; i++) {
        ten[i][1] = 0.0;
        j1[i][1] = 0;
        j2[i][1] = 0;
        j3[i][1] = 0;
        ten[i][0] = 1.0;
        j1[i][0] = 0;
        j2[i][0] = 0;
        j3[i][0] = 0;
    }

    for (i3 = 1; i3 < n3 - 1; i3++) {
        for (i2 = 1; i2 < n2 - 1; i2++) {
            for (i1 = 1; i1 < n1 - 1; i1++) {
                const int idx = (n2 * n1 * i3) + (n1 * i2) + (i1);
                if (z[idx] > ten[0][1]) {
                    ten[0][1] = z[idx];
                    j1[0][1] = i1;
                    j2[0][1] = i2;
                    j3[0][1] = i3;
                    bubble(ten, j1, j2, j3, mm, 1);
                }
                if (z[idx] < ten[0][0]) {
                    ten[0][0] = z[idx];
                    j1[0][0] = i1;
                    j2[0][0] = i2;
                    j3[0][0] = i3;
                    bubble(ten, j1, j2, j3, mm, 0);
                }
            }
        }
    }


    //---------------------------------------------------------------------
    // Now which of these are globally best?
    //---------------------------------------------------------------------
    i1 = mm - 1;
    i0 = mm - 1;
    for (i = mm - 1; i >= 0; i--) {
        best = 0.0;
        if (best < ten[i1][1]) {
            jg[0][i][1] = 0;
            jg[1][i][1] = *is1 - 2 + j1[i1][1];
            jg[2][i][1] = *is2 - 2 + j2[i1][1];
            jg[3][i][1] = *is3 - 2 + j3[i1][1];
            i1 = i1 - 1;
        } else {
            jg[0][i][1] = 0;
            jg[1][i][1] = 0;
            jg[2][i][1] = 0;
            jg[3][i][1] = 0;
        }

        best = 1.0;
        if (best > ten[i0][0]) {
            jg[0][i][0] = 0;
            jg[1][i][0] = *is1 - 2 + j1[i0][0];
            jg[2][i][0] = *is2 - 2 + j2[i0][0];
            jg[3][i][0] = *is3 - 2 + j3[i0][0];
            i0 = i0 - 1;
        } else {
            jg[0][i][0] = 0;
            jg[1][i][0] = 0;
            jg[2][i][0] = 0;
            jg[3][i][0] = 0;
        }

    }


    for (i3 = 0; i3 < n3; i3++) {
        for (i2 = 0; i2 < n2; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
                const int idx = (n2 * n1 * i3) + (n1 * i2) + (i1);
                z[idx] = 0.0;
            }
        }
    }
    for (i = mm - 1; i >= 0; i--) {
        const long i1 = jg[3][i][0];
        const long i2 = jg[2][i][0];
        const long i3 = jg[1][i][0];

        const long idx = wrap((n2 * n1 * i1) + (n1 * i2) + (i3), 0,
                              n1 * n2 *
                              n3); //wrapping for lower problem sizes where this would lead to negative numbers
        /*if (idx < 0) {
            printf("idx 1 below zero: %ld\n", idx);
        } else {
            printf("idx 1 above zero: %ld\n", idx);
        }*/
        z[idx] = -1.0;
    }
    for (i = mm - 1; i >= 0; i--) {
        const long i1 = jg[3][i][1];
        const long i2 = jg[2][i][1];
        const long i3 = jg[1][i][1];

        const long idx = wrap((n2 * n1 * i1) + (n1 * i2) + (i3), 0,
                              n1 * n2 *
                              n3); //wrapping for lower problem sizes where this would lead to negative numbers
        //TODO: this index is always zero
        /*if (idx < 0) {
            printf("idx 2 below zero: %ld\n", idx);
        } else {
            printf("idx 2 above zero: %ld\n", idx);
        }*/
        z[idx] = 1.0;
    }
    comm3(z, n1, n2, n3, k);

    //---------------------------------------------------------------------
    // showall(z,n1,n2,n3);
    //---------------------------------------------------------------------
}

#undef MMCONST

//---------------------------------------------------------------------
// power  raises an integer, disguised as a double
// precision real, to an integer power
//---------------------------------------------------------------------
double power(double a, int n) {
    double aj;
    int nj;
    double rdummy;
    double power;

    power = 1.0;
    nj = n;
    aj = a;

    while (nj != 0) {
        if ((nj % 2) == 1) rdummy = randlc(&power, aj);
        rdummy = randlc(&aj, aj);
        nj = nj / 2;
    }

    return power;
}


//---------------------------------------------------------------------
// bubble        does a bubble sort in direction dir
//---------------------------------------------------------------------
void bubble(double ten[][2], int j1[][2], int j2[][2], int j3[][2],
            int m, int ind) {
    //all array parameters are [10][2] in dimension

    double temp;
    int j_temp;

    if (ind == 1) {
        for (int i = 0; i < m - 1; i++) {
            if (ten[i][ind] > ten[i + 1][ind]) {
                temp = ten[i + 1][ind];
                ten[i + 1][ind] = ten[i][ind];
                ten[i][ind] = temp;

                j_temp = j1[i + 1][ind];
                j1[i + 1][ind] = j1[i][ind];
                j1[i][ind] = j_temp;

                j_temp = j2[i + 1][ind];
                j2[i + 1][ind] = j2[i][ind];
                j2[i][ind] = j_temp;

                j_temp = j3[i + 1][ind];
                j3[i + 1][ind] = j3[i][ind];
                j3[i][ind] = j_temp;
            } else {
                return;
            }
        }
    } else {
        for (int i = 0; i < m - 1; i++) {
            if (ten[i][ind] < ten[i + 1][ind]) {

                temp = ten[i + 1][ind];
                ten[i + 1][ind] = ten[i][ind];
                ten[i][ind] = temp;

                j_temp = j1[i + 1][ind];
                j1[i + 1][ind] = j1[i][ind];
                j1[i][ind] = j_temp;

                j_temp = j2[i + 1][ind];
                j2[i + 1][ind] = j2[i][ind];
                j2[i][ind] = j_temp;

                j_temp = j3[i + 1][ind];
                j3[i + 1][ind] = j3[i][ind];
                j3[i][ind] = j_temp;
            } else {
                return;
            }
        }
    }
}


void zero3(double *z, int n1, int n2, int n3) {
    int i1, i2, i3;

    for (i3 = 0; i3 < n3; i3++) {
        for (i2 = 0; i2 < n2; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
                const int idx = (n2 * n1 * i3) + (n1 * i2) + (i1);
                z[idx] = 0.0;
            }
        }
    }
}

