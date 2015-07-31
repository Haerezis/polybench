#ifndef _MATRIX_INVERSION_H
#define _MATRIX_INVERSION_H

#define MAX_VALUE 42.0
#define EPSILON 1.0e-7

/* Default to LARGE_DATASET. */
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif

#define MAT_SIZE 3

# if !defined(MAT_SIZE) && !defined(MAT_SIZE)
/* Define sample dataset sizes. */
#  ifdef MINI_DATASET
#   define MAT_SIZE 32
#  endif 

#  ifdef SMALL_DATASET
#   define MAT_SIZE 100
#  endif 

#  ifdef MEDIUM_DATASET
#   define MAT_SIZE 260
#  endif 

#  ifdef LARGE_DATASET
#   define MAT_SIZE 1400
#  endif 

#  ifdef EXTRALARGE_DATASET
#   define MAT_SIZE 3000
#  endif 


#endif /* !(MAT_SIZE) */

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
#  define DATA_PRINTF_MODIFIER "%0.5lf "
#  define SCALAR_VAL(x) x
#  define SQRT_FUN(x) sqrt(x)
#  define EXP_FUN(x) exp(x)
#  define POW_FUN(x,y) pow(x,y)
# endif

typedef struct compex_number {
  DATA_TYPE r;
  DATA_TYPE i;
} complex_number_t, * complex_number_p;

#define float_equals(a, b) (fabs((a)-(b)) < EPSILON)

#endif /* !_COVARIANCE_H */

