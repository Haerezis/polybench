#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/*#define MINI_DATASET*/
/* Include benchmark-specific header. */
#include "LMS_GSC.h"


/* Array initialization. */
static
void init_array(int k, int p,
    DATA_TYPE POLYBENCH_2D(V_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(V_i, K, P-1, k, p-1),
    DATA_TYPE d_r,
    DATA_TYPE d_i,
    DATA_TYPE POLYBENCH_2D(X_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(X_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_1D(y_r, K, k),
    DATA_TYPE POLYBENCH_1D(y_i, K, k),
    DATA_TYPE * e_r,
    DATA_TYPE * e_i)
{
  unsigned int l = 0, c = 0;

  d_r = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
  d_i = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
  *e_r = 0.0;
  *e_i = 0.0;

  for(l = 0 ; l<k ; l++)
  {
    y_r[l] = 0.0;
    y_i[l] = 0.0;

    for(c = 0 ; c<(p-1) ; c++)
    {
      //V represents the initial signal received.
      //It is declared and intialize here, but is useless.
      //We kept it if in the future we want to implement
      //completely the algorithm
      V_r[l][c] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      V_i[l][c] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;

      X_r[l][c] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      X_i[l][c] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      
      W_r[l][c] = 0.0;
      W_i[l][c] = 0.0;
    }
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int p,
    DATA_TYPE e_r,
    DATA_TYPE e_i)
{
  fprintf(stderr, "e : ");
  fprintf (stderr, DATA_PRINTF_MODIFIER, e_r);
  if (e_i >= 0.0) fprintf (stderr, "+");
  fprintf (stderr, DATA_PRINTF_MODIFIER, e_r);
  fprintf (stderr, "\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_LMS_GSC(int k, int mu, int p,
    DATA_TYPE d_r,
    DATA_TYPE d_i,
    DATA_TYPE POLYBENCH_2D(X_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(X_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_1D(y_r, K, k),
    DATA_TYPE POLYBENCH_1D(y_i, K, k),
    DATA_TYPE * e_r,
    DATA_TYPE * e_i)
{
  unsigned int i = 0, j = 0;

  unsigned yy_r = 0;
  unsigned yy_i = 0;

#pragma scop
  for (i = 0 ; i < k ; i++)
  {
    y_r[i] = 0;
    y_i[i] = 0;
    for (j = 0 ; j < (p-1) ; j++)
    {
      y_r[i] += W_r[i][j]*X_r[i][j] - (-W_i[i][j])*X_i[i][j];
      y_i[i] += (-W_i[i][j])*X_r[i][j] + W_r[i][j]*X_i[i][j];
    }
    yy_r += y_r[i];
    yy_i += y_i[i];
  }
  *e_r = d_r - yy_r;
  *e_i = d_i - yy_i;

  for (i = 0 ; i < k ; i++)
  {
    for (j = 0 ; j <(p-1) ; j++)
    {
      W_r[i][j] = W_r[i][j] + mu * (X_r[i][j] * (*e_r) - X_i[i][j] * (-*e_i));
      W_i[i][j] = W_i[i][j] + mu * (X_r[i][j] * (-*e_i) - X_i[i][j] * (*e_r));
    }
  }


#pragma endscop

}


int main(int argc, char** argv)
{
  unsigned int i = 0;

  /* Retrieve problem size. */
  int p = P;
  int k = K;
  int mu = mu;
  
  
  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(V_r, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(V_i, DATA_TYPE, K, P-1, k, p-1);

  DATA_TYPE d_r;
  DATA_TYPE d_i;
  
  POLYBENCH_2D_ARRAY_DECL(X_r, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(X_i, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(W_r, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(W_i, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_1D_ARRAY_DECL(y_r, DATA_TYPE, K, k);
  POLYBENCH_1D_ARRAY_DECL(y_i, DATA_TYPE, K, k);
  
  DATA_TYPE e_r;
  DATA_TYPE e_i;

  /* Initialize array(s). */
  init_array (k,
      p,
      POLYBENCH_ARRAY(V_r),
      POLYBENCH_ARRAY(V_i),
      d_r,
      d_i,
      POLYBENCH_ARRAY(X_r),
      POLYBENCH_ARRAY(X_i),
      POLYBENCH_ARRAY(W_r),
      POLYBENCH_ARRAY(W_i),
      POLYBENCH_ARRAY(y_r),
      POLYBENCH_ARRAY(y_i),
      &e_r,
      &e_i);

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  for (i = 0 ; i < N ; i++)
  {
    kernel_LMS_GSC (k,
        mu,
        p,
        d_r,
        d_i,
        POLYBENCH_ARRAY(X_r),
        POLYBENCH_ARRAY(X_i),
        POLYBENCH_ARRAY(W_r),
        POLYBENCH_ARRAY(W_i),
        POLYBENCH_ARRAY(y_r),
        POLYBENCH_ARRAY(y_i),
        &e_r,
        &e_i);
  }

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(p, e_r, e_i));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(V_r);
  POLYBENCH_FREE_ARRAY(V_i);
  POLYBENCH_FREE_ARRAY(X_r);
  POLYBENCH_FREE_ARRAY(X_i);
  POLYBENCH_FREE_ARRAY(W_r);
  POLYBENCH_FREE_ARRAY(W_i);
  POLYBENCH_FREE_ARRAY(y_r);
  POLYBENCH_FREE_ARRAY(y_i);

  return 0;
}
