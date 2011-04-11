#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "instrument.h"

/* Default problem size. */
#ifndef NX
# define NX 8000
#endif
#ifndef NY
# define NY 8000
#endif

/* Default data type is double. */
#ifndef DATA_TYPE
# define DATA_TYPE double
#endif

/* Array declaration. Enable malloc if POLYBENCH_TEST_MALLOC. */
#ifndef POLYBENCH_TEST_MALLOC
DATA_TYPE A[NX][NY];
DATA_TYPE x[NY];
DATA_TYPE y[NY];
DATA_TYPE tmp[NX];
#else
DATA_TYPE** A = (DATA_TYPE**)malloc(NX * sizeof(DATA_TYPE*));
DATA_TYPE* x = (DATA_TYPE*)malloc(NY * sizeof(DATA_TYPE));
DATA_TYPE* y = (DATA_TYPE*)malloc(NY * sizeof(DATA_TYPE));
DATA_TYPE* tmp = (DATA_TYPE*)malloc(NX * sizeof(DATA_TYPE));
{
  int i;
  for (i = 0; i < NX; ++i)
    A[i] = (DATA_TYPE*)malloc(NY * sizeof(DATA_TYPE));
}
#endif

static inline
void init_array()
{
  int i, j;

  for (i = 0; i < NX; i++)
    {
      x[i] = i * M_PI;
      for (j = 0; j < NY; j++)
	A[i][j] = ((DATA_TYPE) i*j) / NX;
    }
}

/* Define the live-out variables. Code is not executed unless
   POLYBENCH_DUMP_ARRAYS is defined. */
static inline
void print_array(int argc, char** argv)
{
  int i, j;
#ifndef POLYBENCH_DUMP_ARRAYS
  if (argc > 42 && ! strcmp(argv[0], ""))
#endif
    {
      for (i = 0; i < NX; i++) {
	fprintf(stderr, "%0.2lf ", y[i]);
	if (i%80 == 20) fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
    }
}


int main(int argc, char** argv)
{
  int i, j;
  int nx = NX;
  int ny = NY;

  /* Initialize array. */
  init_array();

  /* Start timer. */
  polybench_start_instruments;

#pragma scop
#pragma live-out y
  for (i= 0; i < nx; i++)
    y[i] = 0;
  for (i = 0; i < ny; i++)
    {
      tmp[i] = 0;
      for (j = 0; j < ny; j++)
	tmp[i] = tmp[i] + A[i][j] * x[j];
      for (j = 0; j < ny; j++)
	y[j] = y[j] + A[i][j] * tmp[i];
    }
#pragma endscop

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  print_array(argc, argv);

  return 0;
}
