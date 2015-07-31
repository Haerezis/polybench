#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/*#define MINI_DATASET*/
/* Include benchmark-specific header. */
#include "TD-LMS_GSC.h"


/* Array initialization. */
static
void init_array(int k, int p,
    complex_number_t POLYBENCH_2D(V, K, P-1, k, p-1),
    complex_number_t d,
    complex_number_t POLYBENCH_2D(X, K, P-1, k, p-1),
    complex_number_t POLYBENCH_1D(u, K, k),
    DATA_TYPE POLYBENCH_1D(cst1, P-1, p-1),
    complex_number_t POLYBENCH_2D(F1, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(P_, K, P-1, k, p-1),
    complex_number_t POLYBENCH_2D(F2, K, P-1, k, p-1),
    complex_number_t POLYBENCH_2D(W, K, P-1, k, p-1),
    complex_number_t POLYBENCH_1D(y, K, k),
    complex_number_t * e)
{
  unsigned int l = 0, m = 0;
  double j = 4.2;

  d.r = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
  d.i = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
  e->r = 0.0;
  e->i = 0.0;


  for(m = 0 ; m<(P-1) ; m++)
  {
    cst1[m] = RHO * EXP_FUN(-j * ((2.0 * M_PI * (double)m) / P));
  }

  for(l = 0 ; l<K ; l++)
  {
    y[l].r = 0.0;
    y[l].i = 0.0;
    u[l].r = 0.0;
    u[l].i = 0.0;

    for(m = 0 ; m<(P-1) ; m++)
    {
      //V represents the initial signal remeived.
      //It is demlared and intialize here, but is useless.
      //We kept it if in the future we want to implement
      //mompletely the algorithm
      V[l][m].r = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      V[l][m].i = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;

      X[l][m].r = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      X[l][m].i = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      
      F1[l][m].r = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      F1[l][m].i = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;

      P_[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      
      F2[l][m].r = 0.0;
      F2[l][m].i = 0.0;
      
      W[l][m].r = 0.0;
      W[l][m].i = 0.0;
    }
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(complex_number_t e)
{
  fprintf(stderr, "e : ");
  fprintf (stderr, DATA_PRINTF_MODIFIER, e.r);
  if (e.i >= 0.0) fprintf (stderr, "+");
  fprintf (stderr, DATA_PRINTF_MODIFIER, e.r);
  fprintf (stderr, "\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_TD_LMS_GSC(int k, int mu, int p,
    complex_number_t d,
    complex_number_t POLYBENCH_2D(X, K, P-1, k, p-1),
    complex_number_t POLYBENCH_1D(u, K, k),
    DATA_TYPE POLYBENCH_1D(cst1, P-1, p-1),
    complex_number_t POLYBENCH_2D(F1, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(P_, K, P-1, k, p-1),
    complex_number_t POLYBENCH_2D(F2, K, P-1, k, p-1),
    complex_number_t POLYBENCH_2D(W, K, P-1, k, p-1),
    complex_number_t POLYBENCH_1D(y, K, k),
    complex_number_t * e)
{
  unsigned int i = 0, m = 0;
  unsigned int n = P-1;

  complex_number_t yy = {0.0, 0.0};

#pragma scop

  for (i = 0 ; i < K ; i++)
  {
    u[i].r = X[i][n-1].r - POW_FUN(RHO, P) * X[i][n - P + 1].r;
    u[i].i = X[i][n-1].i - POW_FUN(RHO, P) * X[i][n - P + 1].i;

    y[i].r = 0.0;
    y[i].i = 0.0;

    for (m = 0 ; m < (P-1) ; m++)
    {
      F1[i][m].r = cst1[m] * F1[i][m].r + u[i].r;
      F1[i][m].i = cst1[m] * F1[i][m].i + u[i].i;

      P_[i][m] = LAMBDA * P_[i][m] + (1 - LAMBDA) * POW_FUN(F1[i][m].r, 2) * POW_FUN(F1[i][m].i, 2);

      F2[i][m].r = MU * F1[i][m].r / P_[i][m];
      F2[i][m].i = MU * F1[i][m].i / P_[i][m];

      //Access to W.i are negative because we access to the conjugate
      y[i].r += W[i][m].r*F1[i][m].r - (-W[i][m].i)*F1[i][m].i;
      y[i].i += (-W[i][m].i)*F1[i][m].r + W[i][m].r*F1[i][m].i;
    }
    yy.r += y[i].r;
    yy.i += y[i].i;
  }
  e->r = d.r - yy.r;
  e->i = d.i - yy.i;

  for (i = 0 ; i < k ; i++)
  {
    for (m = 0 ; m <(p-1) ; m++)
    {
      W[i][m].r = W[i][m].r + mu * (F2[i][m].r * (e->r) - F2[i][m].i * (-e->i));
      W[i][m].i = W[i][m].i + mu * (F2[i][m].r * (-e->i) - F2[i][m].i * (e->r));
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
  POLYBENCH_2D_ARRAY_DECL(V, complex_number_t, K, P-1, k, p-1);

  complex_number_t d;
  POLYBENCH_2D_ARRAY_DECL(X, complex_number_t, K, P-1, k, p-1);
  POLYBENCH_1D_ARRAY_DECL(u, complex_number_t, K, k);
  POLYBENCH_1D_ARRAY_DECL(cst1, DATA_TYPE, P-1, p-1);
  POLYBENCH_2D_ARRAY_DECL(F1, complex_number_t, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(P_, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(F2, complex_number_t, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(W, complex_number_t, K, P-1, k, p-1);
  POLYBENCH_1D_ARRAY_DECL(y, complex_number_t, K, k);
  complex_number_t e;


  /* Initialize array(s). */
  init_array (k,
      p,
      POLYBENCH_ARRAY(V),
      d,
      POLYBENCH_ARRAY(X),
      POLYBENCH_ARRAY(u),
      POLYBENCH_ARRAY(cst1),
      POLYBENCH_ARRAY(F1),
      POLYBENCH_ARRAY(P_),
      POLYBENCH_ARRAY(F2),
      POLYBENCH_ARRAY(W),
      POLYBENCH_ARRAY(y),
      &e);

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  for (i = 0 ; i < N ; i++)
  {
    kernel_TD_LMS_GSC (k,
        mu,
        p,
        d,
        POLYBENCH_ARRAY(X),
        POLYBENCH_ARRAY(u),
        POLYBENCH_ARRAY(cst1),
        POLYBENCH_ARRAY(F1),
        POLYBENCH_ARRAY(P_),
        POLYBENCH_ARRAY(F2),
        POLYBENCH_ARRAY(W),
        POLYBENCH_ARRAY(y),
        &e);
  }

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(e));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(V);
  POLYBENCH_FREE_ARRAY(X);
  POLYBENCH_FREE_ARRAY(u);
  POLYBENCH_FREE_ARRAY(cst1);
  POLYBENCH_FREE_ARRAY(F1);
  POLYBENCH_FREE_ARRAY(P_);
  POLYBENCH_FREE_ARRAY(F2);
  POLYBENCH_FREE_ARRAY(W);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
