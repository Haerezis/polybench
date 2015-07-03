#ifndef __LMS_H
#define __LMS_H

/* Constant value */
//N corresponds to the number of sensor
#define N 1000
//M corresponds to the number of linear constraint in the system
#define M 100

//ALPHA is a positive constant that controls the convergence speed of the algo.
#define ALPHA 4.2

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



/* Default to LARGE_DATASET. */
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif

# if !defined(NB_RUN)
/* Define sample dataset sizes. */
#  ifdef MINI_DATASET
#   define NB_RUN 10
#  endif 

#  ifdef SMALL_DATASET
#   define NB_RUN 50
#  endif 

#  ifdef MEDIUM_DATASET
#   define NB_RUN 400
#  endif 

#  ifdef LARGE_DATASET
#   define NB_RUN 1000
#  endif 

#  ifdef EXTRALARGE_DATASET
#   define NB_RUN 10000
#  endif 

#endif /* !(NB_RUN) */


#endif
