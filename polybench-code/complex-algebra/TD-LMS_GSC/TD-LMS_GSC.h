#ifndef __LMS_GSC_H
#define __LMS_GSC_H

/* Constant value */
//K corresponds to the number of sensor
#define K 1000
//P corresponds to the number+1 of delay elements
//associated with each elementary array input
#define P 1000

//Constants
#define LAMBDA 0.5
#define RHO 0.5
//MU is a positive constant that controls the convergence speed of the algo.
#define MU 1.5




/* Default data type */
# if !defined(DATA_TYPE_IS_INT) && !defined(DATA_TYPE_IS_FLOAT) && !defined(DATA_TYPE_IS_DOUBLE)
#  define DATA_TYPE_IS_DOUBLE
# endif

#ifdef DATA_TYPE_IS_INT
#  define DATA_TYPE int
#  define DATA_PRINTF_MODIFIER "%d "
#endif 

#ifdef DATA_TYPE_IS_FLOAT
#  define DATA_TYPE float
#  define DATA_PRINTF_MODIFIER "%0.2f "
#  define SCALAR_VAL(x) x##f
#  define SQRT_FUN(x) sqrtf(x)
#  define EXP_FUN(x) expf(x)
#  define POW_FUN(x,y) powf(x,y)
# endif

#ifdef DATA_TYPE_IS_DOUBLE
#  define DATA_TYPE double
#  define DATA_PRINTF_MODIFIER "%0.2lf "
#  define SCALAR_VAL(x) x
#  define SQRT_FUN(x) sqrt(x)
#  define EXP_FUN(x) exp(x)
#  define POW_FUN(x,y) pow(x,y)
# endif


typedef struct compex_number {
  DATA_TYPE r;
  DATA_TYPE i;
} complex_number_t, * complex_number_p;



/* Default to LARGE_DATASET. */
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif

# if !defined(N)
/* Define sample dataset sizes. */
#  ifdef MINI_DATASET
#   define N 10
#  endif 

#  ifdef SMALL_DATASET
#   define N 50
#  endif 

#  ifdef MEDIUM_DATASET
#   define N 400
#  endif 

#  ifdef LARGE_DATASET
#   define N 1000
#  endif 

#  ifdef EXTRALARGE_DATASET
#   define N 10000
#  endif 

#endif /* !(N) */


#endif
