
#include <stdlib.h>
#include <stdio.h>


typedef struct is_parameters {
    const long MAX_KEY_LOG_2;
    const long NUM_BUCKETS_LOG_2;
    const long NUM_KEYS_LOG_2; //=TOTAL_KEYS
    const long SIZE_OF_BUFFERS; //=NUM_KEYS
    const int MAX_ITERATIONS; //=10
    long* key_array; //SIZE_OF_BUFFERS long (TOTAL_KEYS), min 22
    long* key_buff1; //MAX_KEY long
    long* key_buff2; //SIZE_OF_BUFFERS long (TOTAL_KEYS)
    long* bucket_size; //NUM_BUCKETS long
    long* bucket_ptrs; //NUM_BUCKETS long
} is_parameters_t;





double	randlc( double *X, double *A )
{
      static int        KS=0;
      static double	R23, R46, T23, T46;
      double		T1, T2, T3, T4;
      double		A1;
      double		A2;
      double		X1;
      double		X2;
      double		Z;
      int     		i, j;

      if (KS == 0) 
      {
        R23 = 1.0;
        R46 = 1.0;
        T23 = 1.0;
        T46 = 1.0;
    
        for (i=1; i<=23; i++)
        {
          R23 = 0.50 * R23;
          T23 = 2.0 * T23;
        }
        for (i=1; i<=46; i++)
        {
          R46 = 0.50 * R46;
          T46 = 2.0 * T46;
        }
        KS = 1;
      }

/*  Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.  */

      T1 = R23 * *A;
      j  = T1;
      A1 = j;
      A2 = *A - T23 * A1;

/*  Break X into two parts such that X = 2^23 * X1 + X2, compute
    Z = A1 * X2 + A2 * X1  (mod 2^23), and then
    X = 2^23 * Z + A2 * X2  (mod 2^46).                            */

      T1 = R23 * *X;
      j  = T1;
      X1 = j;
      X2 = *X - T23 * X1;
      T1 = A1 * X2 + A2 * X1;
      
      j  = R23 * T1;
      T2 = j;
      Z = T1 - T23 * T2;
      T3 = T23 * Z + A2 * X2;
      j  = R46 * T3;
      T4 = j;
      *X = T3 - T46 * T4;
      return(R46 * *X);
} 




/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/

void	create_seq( double seed, double a, const long MAX_KEY, const long NUM_KEYS, long* key_array)
{
	double x;

    const long k = MAX_KEY/4;

	for (long i=0; i < NUM_KEYS; i++)
	{
	    x = randlc(&seed, &a);
	    x += randlc(&seed, &a);
    	    x += randlc(&seed, &a);
	    x += randlc(&seed, &a);  

            key_array[i] = k*x;
	}
}






/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/


void rank( int iteration ,
const long MAX_KEY_LOG_2, //
const long MAX_KEY, //
const long NUM_BUCKETS_LOG_2, //
const long NUM_BUCKETS, //
const long NUM_KEYS, //
const int MAX_ITERATIONS, //
long* key_array, //
long* key_buff1,//
long* key_buff2,//
long* bucket_size,//
long* bucket_ptrs //
)
{

    const int shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;


    key_array[iteration] = iteration;
    key_array[iteration+MAX_ITERATIONS] = MAX_KEY - iteration;





/*  Initialize */
    for(long i=0; i<NUM_BUCKETS; i++ )
        bucket_size[i] = 0;

/*  Determine the number of keys in each bucket */
    for(long i=0; i<NUM_KEYS; i++ ) {
        const int kai = key_array[i];
        const int idx = kai >> shift;
        const long old = bucket_size[idx];
        if (idx < 0 || idx >= NUM_BUCKETS) {
            printf("Error: while determing number of keys, idx was %d and NUM_BUCKETS was %ld", idx, NUM_BUCKETS);
        }
        bucket_size[idx] = old + 1;
    }


/*  Accumulative bucket sizes are the bucket pointers */
    bucket_ptrs[0] = 0;
    for(long i=1; i< NUM_BUCKETS; i++ )
        bucket_ptrs[i] = bucket_ptrs[i-1] + bucket_size[i-1];


/*  Sort into appropriate bucket */
    for(long i=0; i<NUM_KEYS; i++ )
    {
        const long key = key_array[i];
        key_buff2[bucket_ptrs[key >> shift]++] = key;
    }



/*  Clear the work array */
    for(long i=0; i<MAX_KEY; i++ )
        key_buff1[i] = 0;



/*  In this section, the keys themselves are used as their 
    own indexes to determine how many of each there are: their
    individual population                                       */

    for(long i=0; i<NUM_KEYS; i++ )
        key_buff1[key_buff2[i]]++;  /* Now they have individual key   */
                                       /* population                     */

/*  To obtain ranks of each key, successively add the individual key
    population                                                  */


    for(long i=0; i<MAX_KEY-1; i++ )
        key_buff1[i+1] += key_buff1[i];



}


/*****************************************************************/
/*************             M  A  I  N             ****************/
/*****************************************************************/

int integerSort( long MAX_KEY_LOG_2, //
                 long NUM_BUCKETS_LOG_2, //
                 long NUM_KEYS_LOG_2, //
                 int MAX_ITERATIONS, //
                 long* key_array, //
                 long* key_buff1,//
                 long* key_buff2,//
                 long* bucket_size,//
                 long* bucket_ptrs //
 ) {


    long NUM_KEYS = 1L << NUM_KEYS_LOG_2;
    long MAX_KEY = 1L << MAX_KEY_LOG_2;
    long NUM_BUCKETS = 1L << NUM_BUCKETS_LOG_2;


/*  Generate random number sequence and subsequent keys on all procs */
    create_seq(314159265.00,                    /* Random number gen seed */
               1220703125.00, MAX_KEY, NUM_KEYS, key_array
    );                 /* Random number gen mult */



/*  This is the main iteration */
    for (int iteration = 1; iteration <= MAX_ITERATIONS; iteration++) {
        rank(iteration, MAX_KEY_LOG_2, MAX_KEY, NUM_BUCKETS_LOG_2, NUM_BUCKETS, NUM_KEYS, MAX_ITERATIONS, key_array, key_buff1, key_buff2, bucket_size, bucket_ptrs);
    }


    return 0;
}



int integerSortPacked( is_parameters_t parameters ) {
    return integerSort(parameters.MAX_KEY_LOG_2, parameters.NUM_BUCKETS_LOG_2, parameters.NUM_KEYS_LOG_2, parameters.MAX_ITERATIONS, parameters.key_array, parameters.key_buff1, parameters.key_buff2, parameters.bucket_size, parameters.bucket_ptrs);
}



is_parameters_t buildISParameters(long TK_LOG_2 /* = [16:31] */, long MK_LOG_2 /* = [11:27] */, long NB_LOG_2 /* = [9:10] */) {

    long TK = 1L << TK_LOG_2;
    long MK = 1L << MK_LOG_2;
    long NB = 1L << NB_LOG_2;


    if (TK <= 21) {
        printf("ERROR: TK %ld is not > 21 . This will lead to stack corruption!", TK);
    }


        printf("TK = %ld MK = %ld NB = %ld\n", TK, MK, NB);


    is_parameters_t obj = { .MAX_KEY_LOG_2 = MK_LOG_2, .NUM_BUCKETS_LOG_2 = NB_LOG_2, .NUM_KEYS_LOG_2 = TK_LOG_2, .SIZE_OF_BUFFERS = TK, .MAX_ITERATIONS = 10, .key_array = (long*) malloc((TK)*sizeof(long)), .key_buff1 = (long*) malloc((MK)*sizeof(long)), .key_buff2 = (long*) malloc((TK)*sizeof(long)), .bucket_size = (long*) malloc((NB)*sizeof(long)), .bucket_ptrs = (long*) malloc((NB)*sizeof(long)) };
    return obj;
}



void freeISParameters(is_parameters_t parameters) {
    free(parameters.key_array);
    free(parameters.key_buff1);
    free(parameters.key_buff2);
    free(parameters.bucket_size);
    free(parameters.bucket_ptrs);
}


