/**
 * gramschmidt.c: This file is part of the PolyBench 3.0 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 512. */
#include "gramschmidt.h"


/* Array initialization. */
static
void init_array(int ni, int nj,
		DATA_TYPE POLYBENCH_2D(A,NI,NJ),
		DATA_TYPE POLYBENCH_2D(R,NJ,NJ),
		DATA_TYPE POLYBENCH_2D(Q,NI,NJ))
{
  int i, j;

  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
      A[i][j] = ((DATA_TYPE) i*j) / ni;
      Q[i][j] = ((DATA_TYPE) i*(j+1)) / nj;
    }
  for (i = 0; i < nj; i++)
    for (j = 0; j < nj; j++)
      R[i][j] = ((DATA_TYPE) i*(j+2)) / nj;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nj,
		 DATA_TYPE POLYBENCH_2D(A,NI,NJ),
		 DATA_TYPE POLYBENCH_2D(R,NJ,NJ),
		 DATA_TYPE POLYBENCH_2D(Q,NI,NJ))
{
  int i, j;

  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
	fprintf (stderr, DATA_PRINTF_MODIFIER, A[i][j]);
	if (i % 20 == 0) fprintf (stderr, "\n");
    }
  fprintf (stderr, "\n");
  for (i = 0; i < nj; i++)
    for (j = 0; j < nj; j++) {
	fprintf (stderr, DATA_PRINTF_MODIFIER, R[i][j]);
	if (i % 20 == 0) fprintf (stderr, "\n");
    }
  fprintf (stderr, "\n");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
	fprintf (stderr, DATA_PRINTF_MODIFIER, Q[i][j]);
	if (i % 20 == 0) fprintf (stderr, "\n");
    }
  fprintf (stderr, "\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gramschmidt(int ni, int nj,
			DATA_TYPE POLYBENCH_2D(A,NI,NJ),
			DATA_TYPE POLYBENCH_2D(R,NJ,NJ),
			DATA_TYPE POLYBENCH_2D(Q,NI,NJ))
{
  int i, j, k;

  DATA_TYPE nrm;

#pragma scop
  for (k = 0; k < nj; k++)
    {
      nrm = 0;
      for (i = 0; i < ni; i++)
        nrm += A[i][k] * A[i][k];
      R[k][k] = sqrt(nrm);
      for (i = 0; i < ni; i++)
        Q[i][k] = A[i][k] / R[k][k];
      for (j = k + 1; j < nj; j++)
	{
	  R[k][j] = 0;
	  for (i = 0; i < ni; i++)
	    R[k][j] += Q[i][k] * A[i][j];
	  for (i = 0; i < ni; i++)
	    A[i][j] = A[i][j] - Q[i][k] * R[k][j];
	}
    }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int ni = NI;
  int nj = NJ;

  /* Variable declaration/allocation. */
#ifdef POLYBENCH_HEAP_ARRAYS
  /* Heap arrays use variable 'n' for the size. */
  DATA_TYPE POLYBENCH_2D_ARRAY_DECL(A, ni, ni);
  DATA_TYPE POLYBENCH_2D_ARRAY_DECL(R, nj, nj);
  DATA_TYPE POLYBENCH_2D_ARRAY_DECL(Q, ni, nj);
  A = POLYBENCH_ALLOC_2D_ARRAY(ni, nj, DATA_TYPE);
  R = POLYBENCH_ALLOC_2D_ARRAY(nj, nj, DATA_TYPE);
  Q = POLYBENCH_ALLOC_2D_ARRAY(ni, nj, DATA_TYPE);
#else
  /* Stack arrays use the numerical value 'N' for the size. */
  DATA_TYPE POLYBENCH_2D_ARRAY_DECL(A,NI,NJ);
  DATA_TYPE POLYBENCH_2D_ARRAY_DECL(R,NJ,NJ);
  DATA_TYPE POLYBENCH_2D_ARRAY_DECL(Q,NI,NJ);
#endif

  /* Initialize array(s). */
  init_array (ni, nj,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(R),
	      POLYBENCH_ARRAY(Q));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gramschmidt (ni, nj,
		      POLYBENCH_ARRAY(A),
		      POLYBENCH_ARRAY(R),
		      POLYBENCH_ARRAY(Q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nj, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Q);

  return 0;
}
