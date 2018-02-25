
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "randdp.h"
#include "type.h"


const int FFTBLOCK_DEFAULT    =  16;
const int FFTBLOCKPAD_DEFAULT  = 18;
const int CACHESIZE           =  8192;
const int BLOCKMAX            =  32;

const double PI = 3.141592653589793238;
const double ALPHA = 1.0e-6;


// Cache blocking params. These values are good for most
// RISC processors.
// FFT parameters:
//  fftblock controls how many ffts are done at a time.
//  The default is appropriate for most cache-based machines
//  On vector machines, the FFT can be vectorized with vector
//  length equal to the block size, so the block size should
//  be as large as possible. This is the size of the smallest
//  dimension of the problem: 128 for class A, 256 for class B and
//  512 for class C.

typedef struct ft_parameters {
    const int NX;
    const int NY;
    const int NZ;
    const int NITER_DEFAULT;
    const int NXP;
    const int NYP;
    const unsigned long long NTOTAL;
    const unsigned long long NTOTALP;


    dcomplex *sums; //NITER_DEFAULT+1
    double *twiddle; //NZ|NY|NX+1
    dcomplex *xnt; //NZ|NY|NX+1
    dcomplex *y; //NZ|NY|NX+1
    dcomplex *plane; //(BLOCKMAX+1)*MAXDIM

    dcomplex *scr; //MAXDIM|BLOCKMAX+1

} ft_parameters_t;


//TODO: fix segfault at > 64 on nx or ny or nz
ft_parameters_t buildFTParameters(int nx /* = [64;4096] */, int ny /* = [64;2048] */, int nz /* = [64;2048] */,
                                  int niter /* = [6;25] */) {

    const int maxdim = ny > nx ? (ny > nz ? ny : nz) : (nx > nz ? nx : nz);
    const int nxp = nx + 1;
    const int nyp = ny;
    const unsigned long long ntotal = nx * ny * nz;
    const unsigned long long ntotalp = (nx + 1) * ny * nz;

    ft_parameters_t obj = {.NX = nx, .NY = ny, .NZ = nz, .NITER_DEFAULT = niter, .NXP = nxp, .NYP = nyp, .NTOTAL = ntotal, .NTOTALP = ntotalp,
            .sums = (dcomplex *) malloc((niter+1) * sizeof(dcomplex)),
            .twiddle = (double *) malloc((nz * ny * (nx + 1)) * sizeof(double)),
            .xnt = (dcomplex *) malloc((nz * ny * (nx + 1)) * sizeof(dcomplex)),
            .y = (dcomplex *) malloc((nz * ny * (nx + 1)) * sizeof(dcomplex)),
            .plane =(dcomplex *) malloc(((BLOCKMAX + 1) * maxdim) * sizeof(dcomplex)),
            .scr = (dcomplex *) malloc(((BLOCKMAX + 1) * maxdim) * sizeof(dcomplex))
    };

    return obj;
}

void freeFTParameters(ft_parameters_t parameter) {
    free(parameter.sums);
    free(parameter.twiddle);
    free(parameter.xnt);
    free(parameter.y);
    free(parameter.plane);
    free(parameter.scr);
}


#define dcmplx(r, i)       (dcomplex){r, i}
#define dcmplx_add(a, b)   (dcomplex){(a).real+(b).real, (a).imag+(b).imag}
#define dcmplx_sub(a, b)   (dcomplex){(a).real-(b).real, (a).imag-(b).imag}
#define dcmplx_mul(a, b)   (dcomplex){((a).real*(b).real)-((a).imag*(b).imag),\
                                     ((a).real*(b).imag)+((a).imag*(b).real)}
#define dcmplx_mul2(a, b)  (dcomplex){(a).real*(b), (a).imag*(b)}

static inline dcomplex dcmplx_div(dcomplex z1, dcomplex z2) {
    double a = z1.real;
    double b = z1.imag;
    double c = z2.real;
    double d = z2.imag;

    double divisor = c * c + d * d;
    double real = (a * c + b * d) / divisor;
    double imag = (b * c - a * d) / divisor;
    dcomplex result = (dcomplex) {real, imag};
    return result;
}

#define dcmplx_div2(a, b)  (dcomplex){(a).real/(b), (a).imag/(b)}
#define dcmplx_abs(x)     sqrt(((x).real*(x).real) + ((x).imag*(x).imag))

#define dconjg(x)         (dcomplex){(x).real, -1.0*(x).imag}


int i_log_2(int n) {
    int nn, lg;
    if (n == 1) return 0;

    lg = 1;
    nn = 2;
    while (nn < n) {
        nn = nn * 2;
        lg = lg + 1;
    }
    return lg;
}


//---------------------------------------------------------------------
// compute the roots-of-unity array that will be used for subsequent FFTs.
//---------------------------------------------------------------------
void CompExp(int exponent_array_length, dcomplex *exponent) {
    int m, nu, ku, i, j, ln;
    double t, ti;
    const double pi = 3.141592653589793238;

    const int n = exponent_array_length;

    //printf("CompExp length = %d, pointer = %p\n", exponent_array_length, (void *)exponent);

    nu = n;
    m = i_log_2(n);
    const double md = (double) m;
    exponent[0] = (dcomplex){md, 0.0};
    ku = 2;
    ln = 1;
    for (j = 1; j <= m; j++) {
        t = pi / ln;
        for (i = 0; i <= ln - 1; i++) {
            ti = i * t;
            exponent[i + ku - 1] = dcmplx(cos(ti), sin(ti));
        }
        ku = ku + ln;
        ln = 2 * ln;
    }
}


/* common /blockinfo/ */
//static int fftblock;


//---------------------------------------------------------------------
// Computes NY N-point complex-to-complex FFTs of X using an algorithm due
// to Swarztrauber.  X is both the input and the output array, while Y is a
// scratch array.  It is assumed that N = 2^M.  Before calling
// Swarztrauber to
// perform FFTs
//---------------------------------------------------------------------
static void Swarztrauber(int is, int m, int vlen, int exponent_array_length, int xd1,
                         dcomplex *x, dcomplex *exponent,
                         dcomplex* scr
) {
    const int n = exponent_array_length;

    dcomplex u1, x11, x21;
    int n1, li, lj, lk, ku, i11, i12, i21, i22;

    //---------------------------------------------------------------------
    // Perform one variant of the Stockham FFT.
    //---------------------------------------------------------------------
    n1 = n / 2;
    lj = 1;
    li = 1 << m;
    for (int l = 1; l <= m; l += 2) {
        lk = lj;
        lj = 2 * lk;
        li = li / 2;
        ku = li;

        for (int i = 0; i <= li - 1; i++) {
            i11 = i * lk;
            i12 = i11 + n1;
            i21 = i * lj;
            i22 = i21 + lk;

            if (is >= 1) {
                u1 = exponent[ku + i];
            } else {
                u1 = dconjg(exponent[ku + i]);
            }
            for (int k = 0; k <= lk - 1; k++) {
                for (int j = 0; j < vlen; j++) {
                    x11 = x[(i11 + k)*vlen + j];
                    x21 = x[(i12 + k)*vlen + j];
                    scr[(i21 + k) * (BLOCKMAX+1) + j] = dcmplx_add(x11, x21);
                    scr[(i22 + k) * (BLOCKMAX+1) + j] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
                }
            }
        }

        if (l == m) {
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < vlen; j++) {
                    x[k*vlen + j] = scr[k * vlen + j];
                }
            }
        } else {
            lk = lj;
            lj = 2 * lk;
            li = li / 2;
            ku = li;

            for (int i = 0; i <= li - 1; i++) {
                i11 = i * lk;
                i12 = i11 + n1;
                i21 = i * lj;
                i22 = i21 + lk;

                if (is >= 1) {
                    u1 = exponent[ku + i];
                } else {
                    u1 = dconjg(exponent[ku + i]);
                }
                for (int k = 0; k <= lk - 1; k++) {
                    for (int j = 0; j < vlen; j++) {
                        x11 = scr[(i11 + k) * vlen + j];
                        x21 = scr[(i12 + k) * vlen + j];
                        x[(i21 + k)*vlen + j] = dcmplx_add(x11, x21);
                        x[(i22 + k)*vlen + j] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
                    }
                }
            }
        }
    }
}


void fftXYZ(int sign, int exp1_length, int exp2_length, int exp3_length,
            dcomplex *x, dcomplex *xout,
            dcomplex *exp1, dcomplex *exp2, dcomplex *exp3, int NX,
            int NY, dcomplex* plane, dcomplex* scr) {
    const int n1 = exp1_length;
    const int n2 = exp2_length;
    const int n3 = exp3_length;


    int i, j, k, log;
    int bls, ble;
    int len;
    int blkp;


    int fftblock = CACHESIZE / n1;
    if (fftblock >= BLOCKMAX) fftblock = BLOCKMAX;
    blkp = fftblock + 1;
    log = i_log_2(n1);
    for (k = 0; k < n3; k++) {
        for (bls = 0; bls < n2; bls += fftblock) {
            ble = bls + fftblock - 1;
            if (ble > n2) ble = n2 - 1;
            len = ble - bls + 1;
            for (j = bls; j <= ble; j++) {
                for (i = 0; i < n1; i++) {
                    const int idx = (k * n1 * (ble + 1)) + (j * n1) + i;
                    plane[j - bls + blkp * i] = x[idx];
                }
            }
            Swarztrauber(sign, log, len, n1, blkp, plane, exp1, scr);
            for (j = bls; j <= ble; j++) {
                for (i = 0; i < n1; i++) {
                    const int idx = (k * n1 * (ble + 1)) + (j * n1) + i;
                    x[idx] = plane[j - bls + blkp * i];
                }
            }
        }
    }

    fftblock = CACHESIZE / n2;
    if (fftblock >= BLOCKMAX) fftblock = BLOCKMAX;
    blkp = fftblock + 1;
    log = i_log_2(n2);
    for (k = 0; k < n3; k++) {
        for (bls = 0; bls < n1; bls += fftblock) {
            ble = bls + fftblock - 1;
            if (ble > n1) ble = n1 - 1;
            len = ble - bls + 1;
            Swarztrauber(sign, log, len, n2, n1 + 1,
                         &x[k * NY * (NX + 1) + 0 * (NX + 1) + bls], exp2, scr);
        }
    }

    fftblock = CACHESIZE / n3;
    if (fftblock >= BLOCKMAX) fftblock = BLOCKMAX;
    blkp = fftblock + 1;
    log = i_log_2(n3);
    for (k = 0; k < n2; k++) {
        for (bls = 0; bls < n1; bls += fftblock) {
            ble = bls + fftblock - 1;
            if (ble > n1) ble = n1 - 1;
            len = ble - bls + 1;
            for (i = 0; i < n3; i++) {
                for (j = bls; j <= ble; j++) {
                    const int idx = (i * n2 * (ble + 1)) + (k * (ble + 1)) + j;
                    plane[j - bls + blkp * i] = x[idx];
                }
            }
            Swarztrauber(sign, log, len, n3, blkp, plane, exp3, scr);
            for (i = 0; i <= n3 - 1; i++) {
                for (j = bls; j <= ble; j++) {
                    xout[j + (n1 + 1) * (k + n2 * i)] = plane[j - bls + blkp * i];
                }
            }
        }
    }
}


//---------------------------------------------------------------------
// compute a^exponent mod 2^46
//---------------------------------------------------------------------
static double ipow46(double a, int exponent) {
    double result, dummy, q, r;
    int n, n2;

    //---------------------------------------------------------------------
    // Use
    //   a^n = a^(n/2)*a^(n/2) if n even else
    //   a^n = a*a^(n-1)       if n odd
    //---------------------------------------------------------------------
    result = 1;
    if (exponent == 0) return result;
    q = a;
    r = 1;
    n = exponent;

    while (n > 1) {
        n2 = n / 2;
        if (n2 * 2 == n) {
            dummy = randlc(&q, q);
            n = n2;
        } else {
            dummy = randlc(&r, q);
            n = n - 1;
        }
    }
    dummy = randlc(&r, q);
    result = r;
    return result;
}


void CalculateChecksum(dcomplex *csum, int iterN, int d1, int d2, int d3,
                       dcomplex *u, int NX, int NY) {
    int i, i1, ii, ji, ki;
    dcomplex csum_temp = dcmplx(0.0, 0.0);
    for (i = 1; i <= 1024; i++) {
        i1 = i;
        ii = i1 % d1;
        ji = 3 * i1 % d2;
        ki = 5 * i1 % d3;
        const int index = (ki * NY * (NX + 1)) + (ji * (NX + 1)) + (ii);
        csum_temp = dcmplx_add(csum_temp, u[index]);
    }
    csum_temp = dcmplx_div2(csum_temp, (double) (d1 * d2 * d3));
    *csum = csum_temp;
}


void compute_initial_conditions(int d1, int d2, int d3,
                                dcomplex *u0, int MAXDIM, int NX, int NY
) {
    dcomplex tmp[MAXDIM];
    double x0, start, an, dummy;
    double RanStarts[MAXDIM];

    int i, j, k;
    const double seed = 314159265.0;
    const double a = 1220703125.0;

    start = seed;
    //---------------------------------------------------------------------
    // Jump to the starting element for our first plane.
    //---------------------------------------------------------------------
    an = ipow46(a, 0);
    dummy = randlc(&start, an);
    an = ipow46(a, 2 * d1 * d2);
    //---------------------------------------------------------------------
    // Go through by z planes filling in one square at a time.
    //---------------------------------------------------------------------
    RanStarts[0] = start;
    for (k = 1; k < d3; k++) {
        dummy = randlc(&start, an);
        RanStarts[k] = start;
    }

    for (k = 0; k < d3; k++) {
        x0 = RanStarts[k];
        for (j = 0; j < d2; j++) {
            vranlc(2 * d1, &x0, a, (double *) tmp);
            for (i = 0; i < d1; i++) {
                const int idx = (k * NY * (NX + 1)) + (j * (NX + 1)) + i;
                u0[idx] = tmp[i];
            }
        }
    }
}


void evolve(int nx, int ny, int nz,
            dcomplex *x, dcomplex *y,
            double *twiddle) {
    int i, j, k;
    for (i = 0; i < nz; i++) {
        for (k = 0; k < ny; k++) {
            for (j = 0; j < nx; j++) {
                const int index = (i * (nx * ny)) + (k * nx) + j;
                y[index] = dcmplx_mul2(y[index], twiddle[index]);
                x[(i * ny * nx) + (k * nx) + j] = y[index];
            }
        }
    }
}


//static dcomplex pad1[128], pad2[128];


void fourierTransformation(int niter, int NX,
                           int NY,
                           int NZ,
                           dcomplex* plane,
                           dcomplex* scr,
                           dcomplex* y,
                           dcomplex* xnt,
                           double* twiddle,
                           dcomplex* sums) {
    int i, j, k, kt, n12, n22, n32, ii, jj, kk, ii2, ik2;
    double ap;


    const int MAXDIM = NY > NX ? (NY > NZ ? NY : NZ) : (NX > NZ ? NX : NZ);

    dcomplex* exp1 = (dcomplex*) malloc(NX * sizeof(dcomplex));
    dcomplex* exp2 = (dcomplex*) malloc(NY * sizeof(dcomplex));
    dcomplex* exp3 = (dcomplex*) malloc(NZ * sizeof(dcomplex));


    compute_initial_conditions(NX, NY, NZ, xnt, MAXDIM, NX, NY);

    CompExp(NX, exp1);
    CompExp(NY, exp2);
    CompExp(NZ, exp3);
    fftXYZ(1, NX, NY, NZ, xnt, (dcomplex *) y, exp1, exp2, exp3,
    NX, NY, plane, scr );

    n12 = NX / 2;
    n22 = NY / 2;
    n32 = NZ / 2;
    ap = -4.0 * ALPHA * (PI * PI);
    for (i = 0; i < NZ; i++) {
        ii = i - (i / n32) * NZ;
        ii2 = ii * ii;
        for (k = 0; k < NY; k++) {
            kk = k - (k / n22) * NY;
            ik2 = ii2 + kk * kk;
            for (j = 0; j < NX; j++) {
                jj = j - (j / n12) * NX;
                const int index = i * ((NX + 1) * NY) + k * (NX + 1) + j;
                twiddle[index] = exp(ap * (double) (jj * jj + ik2));
            }
        }
    }
    compute_initial_conditions(NX, NY, NZ, xnt, MAXDIM, NX, NY);
    fftXYZ(1, NX, NY, NZ, xnt, (dcomplex *) y, exp1, exp2, exp3,
           NX, NY, plane, scr );

    for (kt = 1; kt <= niter; kt++) {
        evolve(NX, NY, NZ, xnt, y, twiddle);
        fftXYZ(-1, NX, NY, NZ, xnt, (dcomplex *) xnt, exp1, exp2,
               exp3, NX, NY, plane, scr );
        CalculateChecksum(sums + kt, kt, NX, NY, NZ, xnt,
                          NX, NY);
    }

    free(exp1);
    free(exp2);
    free(exp3);
}


void fourierTransformationPacked(ft_parameters_t parameters) {
    fourierTransformation(parameters.NITER_DEFAULT, parameters.NX, parameters.NY, parameters.NZ, parameters.plane, parameters.scr, parameters.y, parameters.xnt, parameters.twiddle, parameters.sums);
}
