#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/*#define MINI_DATASET*/
/* Include benchmark-specific header. */
#include "LCMP.h"


/* Array update. */
static
void update_array(int n,
    complex_number_t POLYBENCH_1D(x, N, n))
{
  unsigned int i = 0;

  for(i = 0 ; i<n ; i++)
  {
    x[i].r = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;
    x[i].i = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;
  }
}

/* Array initialization. */
static
void init_array(int n, int m,
    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_2D(R, N, N, n, n),
    complex_number_t POLYBENCH_2D(Rinv, N, N, n, n),
    complex_number_t POLYBENCH_2D(C, N, M, n, m),
    complex_number_t POLYBENCH_2D(R_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(Rinv_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(C_R_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(C_R_C_inv, M, N, m, n),
    complex_number_t POLYBENCH_1D(C_R_C_inv_f, M, m),
    complex_number_t POLYBENCH_1D(f, M, m),
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0;

  for(i = 0 ; i<n ; i++)
  {
    x[i] = (complex_number_t) {0.0, 0.0};
    w[i] = (complex_number_t) {0.0, 0.0};
    for (j = 0 ; j<n ; j++)
    {
      R[i][j] = (complex_number_t) {0.0, 0.0};
      Rinv[i][j] = (complex_number_t) {0.0, 0.0};
    }
    for (j = 0 ; j<m ; j++)
    {
      C[i][j] = (complex_number_t) {0.0, 0.0};
      R_C[i][j] = (complex_number_t) {0.0, 0.0};
      Rinv_C[i][j] = (complex_number_t) {0.0, 0.0};
      C_R_C[i][j] = (complex_number_t) {0.0, 0.0};
      C_R_C_inv[j][i] = (complex_number_t) {0.0, 0.0};
    }
  }

  for(i = 0 ; i<m ; i++)
  {
    f[i] = (complex_number_t) {0.0, 0.0};
    C_R_C_inv_f[i] = (complex_number_t) {0.0, 0.0};
  }
  f[0].r = 1.0;
  
  update_array(n,x);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0;

  fprintf(stderr, "w : (");
  for (i = 0 ; i<n ; i++)
  {
    fprintf (stderr, DATA_PRINTF_MODIFIER, w[i].r);
    if (w[i].i >= 0.0) fprintf (stderr, "+");
    fprintf (stderr, DATA_PRINTF_MODIFIER, w[i].r);
    if (i<(n-1)) fprintf (stderr, ", ");
  }
  fprintf (stderr, ")\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_MVDR(int n, int m,
    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_2D(R, N, N, n, n),
    complex_number_t POLYBENCH_2D(Rinv, N, N, n, n),
    complex_number_t POLYBENCH_2D(C, N, M, n, m),
    complex_number_t POLYBENCH_2D(R_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(Rinv_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(C_R_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(C_R_C_inv, M, N, m, n),
    complex_number_t POLYBENCH_1D(C_R_C_inv_f, M, m),
    complex_number_t POLYBENCH_1D(f, M, m),
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0, k = 0;
#pragma scop
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      R[i][j].r = (x[i].r * x[j].r - x[i].i * (-x[j].i)) / N;
      R[i][j].i = (x[i].i * x[j].r + x[i].r * (-x[j].i)) / N;
    }
  }

  //TODO INVERSE R INTO Rinv

  for (i=0 ; i<N ; i++)
  {
    for (j=0 ; j<M ; j++)
    {
      Rinv_C[i][j].r = 0.0;
      Rinv_C[i][j].i = 0.0;
      for (k=0 ; k<N ; k++)
      {
        Rinv_C[i][j].r += (Rinv[i][k].r * C[k][j].r) - (Rinv[i][k].i * C[k][j].i);
        Rinv_C[i][j].i += (Rinv[i][k].i * C[k][j].r) + (Rinv[i][k].r * C[k][j].i);
      }
    }
  }

  for (i=0 ; i<N ; i++)
  {
    for (j=0 ; j<M ; j++)
    {
      R_C[i][j].r = 0.0;
      R_C[i][j].i = 0.0;
      for (k=0 ; k<N ; k++)
      {
        R_C[i][j].r += (R[i][k].r * C[k][j].r) - (R[i][k].i * C[k][j].i);
        R_C[i][j].i += (R[i][k].i * C[k][j].r) + (R[i][k].r * C[k][j].i);
      }
    }
  }

  for (i=0 ; i<N ; i++)
  {
    for (j=0 ; j<M ; j++)
    {
      C_R_C[i][j].r = 0.0;
      C_R_C[i][j].i = 0.0;
    }
  }
  for (k=0 ; k<N ; k++)
  {
    for (i=0 ; i<N ; i++)
    {
      for (j=0 ; j<M ; j++)
      {
        C_R_C[i][j].r += (C[k][i].r * R_C[k][j].r) - ((-C[k][i].i) * R_C[k][j].i);
        C_R_C[i][j].i += ((-C[k][i].i) * R_C[k][j].r) + (C[k][i].r * R_C[k][j].i);
      }
    }
  }

  //TODO INVERSE C_R_C INTO C_R_C_inv

  for (i=0 ; i<M ; i++)
  {
    C_R_C_inv_f[i].r = 0.0;
    C_R_C_inv_f[i].i = 0.0;

    for (j=0 ; j<N ; j++)
    {
      C_R_C_inv_f[i].r += (C_R_C_inv[i][j].r * f[j].r) - (C_R_C_inv[i][j].i * f[j].i);
      C_R_C_inv_f[i].i += (C_R_C_inv[i][j].r * f[j].i) + (C_R_C_inv[i][j].i * f[j].r);
    }
  }


  for (i=0 ; i<M ; i++)
  {
    w[i].r = 0.0;
    w[i].i = 0.0;

    for (j=0 ; j<N ; j++)
    {
      w[i].r += (Rinv_C[i][j].r * C_R_C_inv_f[j].r) - (Rinv_C[i][j].i * C_R_C_inv_f[j].i);
      w[i].i += (Rinv_C[i][j].r * C_R_C_inv_f[j].i) + (Rinv_C[i][j].i * C_R_C_inv_f[j].r);
    }
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  unsigned int i = 0;

  /* Retrieve problem size. */
  int n = N;
  int m = M;
  
  
  /* Variable declaration/allocation. */
  POLYBENCH_1D_ARRAY_DECL(x, complex_number_t, N, n);

  POLYBENCH_2D_ARRAY_DECL(R, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Rinv, complex_number_t, N, N, n, n);
  
  POLYBENCH_2D_ARRAY_DECL(C, complex_number_t, N, M, n, m);
  
  POLYBENCH_2D_ARRAY_DECL(R_C, complex_number_t, N, M, n, m);
  POLYBENCH_2D_ARRAY_DECL(Rinv_C, complex_number_t, N, M, n, m);
  POLYBENCH_2D_ARRAY_DECL(C_R_C, complex_number_t, N, M, n, m);
  POLYBENCH_2D_ARRAY_DECL(C_R_C_inv, complex_number_t, M, N, m, N);
  
  POLYBENCH_1D_ARRAY_DECL(f, complex_number_t, M, m);
  POLYBENCH_1D_ARRAY_DECL(C_R_C_inv_f, complex_number_t, M, m);
  POLYBENCH_1D_ARRAY_DECL(w, complex_number_t, N, n);

  /* Initialize array(s). */
  init_array (n, m,
      POLYBENCH_ARRAY(x),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(Rinv),
      POLYBENCH_ARRAY(C),
      POLYBENCH_ARRAY(R_C),
      POLYBENCH_ARRAY(Rinv_C),
      POLYBENCH_ARRAY(C_R_C),
      POLYBENCH_ARRAY(C_R_C_inv),
      POLYBENCH_ARRAY(C_R_C_inv_f),
      POLYBENCH_ARRAY(f),
      POLYBENCH_ARRAY(w));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  for (i = 0 ; i < NB_RUN ; i++)
  {
    //Just comment this function call if you don't want to 
    //update/modify the input data (x vector) of the radio.
    update_array(n,
        POLYBENCH_ARRAY(x));

    kernel_MVDR (n, m,
      POLYBENCH_ARRAY(x),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(Rinv),
      POLYBENCH_ARRAY(C),
      POLYBENCH_ARRAY(R_C),
      POLYBENCH_ARRAY(Rinv_C),
      POLYBENCH_ARRAY(C_R_C),
      POLYBENCH_ARRAY(C_R_C_inv),
      POLYBENCH_ARRAY(C_R_C_inv_f),
      POLYBENCH_ARRAY(f),
      POLYBENCH_ARRAY(w));
  }

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Rinv);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(R_C);
  POLYBENCH_FREE_ARRAY(Rinv_C);
  POLYBENCH_FREE_ARRAY(C_R_C);
  POLYBENCH_FREE_ARRAY(C_R_C_inv);
  POLYBENCH_FREE_ARRAY(C_R_C_inv_f);
  POLYBENCH_FREE_ARRAY(f);
  POLYBENCH_FREE_ARRAY(w);

  return 0;
}
