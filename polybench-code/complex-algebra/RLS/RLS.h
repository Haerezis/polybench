#ifndef __RLS_H
#define __RLS_H

/* Constant value */
#define NB_RUN 300

//ALPHA is a positive constant that controls the convergence speed of the algo.
#define ALPHA 4.2

#define MU 1.0

#define MAX_VALUE 10.0



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

//N corresponds to the number of sensor
//M corresponds to the number of linear constraint in the system
# if !defined(N) || !defined(M)
/* Define sample dataset sizes. */
#  ifdef MINI_DATASET
#   define N 25
#   define M 5
#  endif 

#  ifdef SMALL_DATASET
#   define N 80
#   define M 20
#  endif 

#  ifdef MEDIUM_DATASET
#   define N 230
#   define M 50
#  endif 

#  ifdef LARGE_DATASET
#   define N 900
#   define M 200
#  endif 

#  ifdef EXTRALARGE_DATASET
#   define N 2600
#   define M 550
#  endif 

#endif /* !(N) || !(M) */


#endif
