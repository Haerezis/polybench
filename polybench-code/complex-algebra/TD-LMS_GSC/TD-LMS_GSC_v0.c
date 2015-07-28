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
    DATA_TYPE POLYBENCH_2D(V_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(V_i, K, P-1, k, p-1),
    DATA_TYPE d_r,
    DATA_TYPE d_i,
    DATA_TYPE POLYBENCH_2D(X_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(X_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_1D(u_r, K, k),
    DATA_TYPE POLYBENCH_1D(u_i, K, k),
    DATA_TYPE POLYBENCH_1D(cst1, P-1, p-1),
    DATA_TYPE POLYBENCH_2D(F1_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(F1_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(P_, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(F2_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(F2_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_1D(y_r, K, k),
    DATA_TYPE POLYBENCH_1D(y_i, K, k),
    DATA_TYPE * e_r,
    DATA_TYPE * e_i)
{
  unsigned int l = 0, m = 0;
  double j = 4.2;

  d_r = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
  d_i = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
  *e_r = 0.0;
  *e_i = 0.0;


  for(m = 0 ; m<(P-1) ; m++)
  {
    cst1[m] = RHO * EXP_FUN(-j * ((2.0 * M_PI * (double)m) / P));
  }

  for(l = 0 ; l<K ; l++)
  {
    y_r[l] = 0.0;
    y_i[l] = 0.0;
    u_r[l] = 0.0;
    u_i[l] = 0.0;

    for(m = 0 ; m<(P-1) ; m++)
    {
      //V represents the initial signal remeived.
      //It is demlared and intialize here, but is useless.
      //We kept it if in the future we want to implement
      //mompletely the algorithm
      V_r[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      V_i[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;

      X_r[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      X_i[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      
      F1_r[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      F1_i[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;

      P_[l][m] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      
      F2_r[l][m] = 0.0;
      F2_i[l][m] = 0.0;
      
      W_r[l][m] = 0.0;
      W_i[l][m] = 0.0;
    }
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(DATA_TYPE e_r, DATA_TYPE e_i)
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
    DATA_TYPE POLYBENCH_1D(u_r, K, k),
    DATA_TYPE POLYBENCH_1D(u_i, K, k),
    DATA_TYPE POLYBENCH_1D(cst1, P-1, p-1),
    DATA_TYPE POLYBENCH_2D(F1_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(F1_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(P_, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(F2_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(F2_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_r, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_2D(W_i, K, P-1, k, p-1),
    DATA_TYPE POLYBENCH_1D(y_r, K, k),
    DATA_TYPE POLYBENCH_1D(y_i, K, k),
    DATA_TYPE * e_r,
    DATA_TYPE * e_i)
{
  unsigned int i = 0, m = 0;
  unsigned int n = P-1;

  DATA_TYPE yy_r = 0;
  DATA_TYPE yy_i = 0;

#pragma scop

  for (i = 0 ; i < K ; i++)
  {
    u_r[i] = X_r[i][n-1] - POW_FUN(RHO, P) * X_r[i][n - P + 1];
    u_i[i] = X_i[i][n-1] - POW_FUN(RHO, P) * X_i[i][n - P + 1];

    y_r[i] = 0.0;
    y_i[i] = 0.0;

    for (m = 0 ; m < (P-1) ; m++)
    {
      F1_r[i][m] = cst1[m] * F1_r[i][m] + u_r[i];
      F1_i[i][m] = cst1[m] * F1_i[i][m] + u_i[i];

      P_[i][m] = LAMBDA * P_[i][m] + (1 - LAMBDA) * POW_FUN(F1_r[i][m], 2) * POW_FUN(F1_i[i][m], 2);

      F2_r[i][m] = MU * F1_r[i][m] / P_[i][m];
      F2_i[i][m] = MU * F1_i[i][m] / P_[i][m];

      //Access to W_i are negative because we access to the conjugate
      y_r[i] += W_r[i][m]*F1_r[i][m] - (-W_i[i][m])*F1_i[i][m];
      y_i[i] += (-W_i[i][m])*F1_r[i][m] + W_r[i][m]*F1_i[i][m];
    }
    yy_r += y_r[i];
    yy_i += y_i[i];
  }
  *e_r = d_r - yy_r;
  *e_i = d_i - yy_i;

  for (i = 0 ; i < k ; i++)
  {
    for (m = 0 ; m <(p-1) ; m++)
    {
      W_r[i][m] = W_r[i][m] + mu * (F2_r[i][m] * (*e_r) - F2_i[i][m] * (-*e_i));
      W_i[i][m] = W_i[i][m] + mu * (F2_r[i][m] * (-*e_i) - F2_i[i][m] * (*e_r));
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
  POLYBENCH_1D_ARRAY_DECL(u_r, DATA_TYPE, K, k);
  POLYBENCH_1D_ARRAY_DECL(u_i, DATA_TYPE, K, k);
  POLYBENCH_1D_ARRAY_DECL(cst1, DATA_TYPE, P-1, p-1);
  POLYBENCH_2D_ARRAY_DECL(F1_r, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(F1_i, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(P_, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(F2_r, DATA_TYPE, K, P-1, k, p-1);
  POLYBENCH_2D_ARRAY_DECL(F2_i, DATA_TYPE, K, P-1, k, p-1);
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
      POLYBENCH_ARRAY(u_r),
      POLYBENCH_ARRAY(u_i),
      POLYBENCH_ARRAY(cst1),
      POLYBENCH_ARRAY(F1_r),
      POLYBENCH_ARRAY(F1_i),
      POLYBENCH_ARRAY(P_),
      POLYBENCH_ARRAY(F2_r),
      POLYBENCH_ARRAY(F2_i),
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
        POLYBENCH_ARRAY(u_r),
        POLYBENCH_ARRAY(u_i),
        POLYBENCH_ARRAY(cst1),
        POLYBENCH_ARRAY(F1_r),
        POLYBENCH_ARRAY(F1_i),
        POLYBENCH_ARRAY(P_),
        POLYBENCH_ARRAY(F2_r),
        POLYBENCH_ARRAY(F2_i),
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
  polybench_prevent_dce(print_array(e_r, e_i));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(V_r);
  POLYBENCH_FREE_ARRAY(V_i);
  POLYBENCH_FREE_ARRAY(X_r);
  POLYBENCH_FREE_ARRAY(X_i);
  POLYBENCH_FREE_ARRAY(u_r);
  POLYBENCH_FREE_ARRAY(u_i);
  POLYBENCH_FREE_ARRAY(cst1);
  POLYBENCH_FREE_ARRAY(F1_r);
  POLYBENCH_FREE_ARRAY(F1_i);
  POLYBENCH_FREE_ARRAY(P_);
  POLYBENCH_FREE_ARRAY(F2_r);
  POLYBENCH_FREE_ARRAY(F2_i);
  POLYBENCH_FREE_ARRAY(W_r);
  POLYBENCH_FREE_ARRAY(W_i);
  POLYBENCH_FREE_ARRAY(y_r);
  POLYBENCH_FREE_ARRAY(y_i);

  return 0;
}
