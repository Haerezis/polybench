/* covariance.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "matrix_covariance.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE *float_n,
		 complex_number_t POLYBENCH_2D(data,M,N,m,n))
{
  unsigned int l = 0, c = 0;

  *float_n = (DATA_TYPE)n;

  for (l = 0; l < M; l++)
    for (c = 0; c < N; c++)
    {
      data[l][c].r = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;
      data[l][c].i = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 complex_number_t POLYBENCH_2D(cov,M,M,m,m))

{
  unsigned int l = 0, c = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("cov");
  for (l = 0; l < m; l++)
  {
    fprintf (POLYBENCH_DUMP_TARGET, "\n");
    for (c = 0; c < m; c++)
    {
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, cov[l][c].r);
      fprintf (POLYBENCH_DUMP_TARGET, "%s", (cov[l][c].i > 0.0) ? "+" : "");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, cov[l][c].i);
    }
  }
  POLYBENCH_DUMP_END("cov");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_covariance(int m, int n,
		       DATA_TYPE float_n,
		       complex_number_t POLYBENCH_2D(data,M,N,m,n),
		       complex_number_t POLYBENCH_2D(cov,M,M,m,m),
		       complex_number_t POLYBENCH_1D(mean,M,m))
{
  unsigned int l = 0, c = 0, k = 0;

#pragma scop
  /*
  for (l = 0; l < _PB_M; l++)
  {
    mean[l].r = SCALAR_VAL(0.0);
    mean[l].i = SCALAR_VAL(0.0);
    for (c = 0; c < _PB_N; c++)
    {
      mean[l].r += data[l][c].r;
      mean[l].i += data[l][c].i;
    }
    mean[l].r /= float_n;
    mean[l].i /= float_n;
  }

  for (l = 0; l < _PB_M; l++)
  {
    for (c = 0; c < _PB_N; c++)
    {
      data[l][c].r -= mean[l].r;
      data[l][c].i -= mean[l].i;
    }
  }
  */

  for (l = 0; l < _PB_M; l++)
  {
    for (c = l; c < _PB_M; c++)
    {
      cov[l][c].r = SCALAR_VAL(0.0);
      cov[l][c].i = SCALAR_VAL(0.0);
      for (k = 0; k < _PB_N; k++)
      {
        cov[l][c].r += data[l][k].r * data[c][k].r - (data[l][k].i * (-data[c][k].i));
        cov[l][c].i += data[l][k].r * (-data[c][k].i) + (data[l][k].r * data[c][k].r);
      }
      cov[l][c].r /= (float_n - SCALAR_VAL(1.0));
      cov[l][c].i /= (float_n - SCALAR_VAL(1.0));
      cov[c][l] = cov[l][c];
    }
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE float_n;
  POLYBENCH_2D_ARRAY_DECL(data, complex_number_t, M, N, m, n);
  POLYBENCH_2D_ARRAY_DECL(cov, complex_number_t, M, M, m, m);
  POLYBENCH_1D_ARRAY_DECL(mean, complex_number_t, M, m);


  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_covariance (m, n, float_n,
		     POLYBENCH_ARRAY(data),
		     POLYBENCH_ARRAY(cov),
		     POLYBENCH_ARRAY(mean));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(cov)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(cov);
  POLYBENCH_FREE_ARRAY(mean);

  return 0;
}
