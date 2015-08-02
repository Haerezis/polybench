#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* 
 * The LCMP algorithm was taken from the following paper :
 * "A Recursive Least Squares Implementation for
 * LCMP Beamforming Under Quadratic Constraint
 */

#define MINI_DATASET
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
		complex_number_t POLYBENCH_2D(qrN_Q, N, N, n, n),
		complex_number_t POLYBENCH_2D(qrN_R, N, N, n, n),
		complex_number_t POLYBENCH_2D(qrN_R_inv, N, N, n, n),
		complex_number_t POLYBENCH_2D(qrM_Q, M, M, m, n),
		complex_number_t POLYBENCH_2D(qrM_R, M, M, m, n),
		complex_number_t POLYBENCH_2D(qrM_R_inv, M, M, m, n),

    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_2D(R, N, N, n, n),
    complex_number_t POLYBENCH_2D(Rinv, N, N, n, n),
    complex_number_t POLYBENCH_2D(C, N, M, n, m),
    complex_number_t POLYBENCH_2D(R_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(Rinv_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(C_R_C, M, M, m, m),
    complex_number_t POLYBENCH_2D(C_R_C_inv, M, M, m, m),
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
    }
  }

  for(i = 0 ; i<m ; i++)
  {
    f[i] = (complex_number_t) {0.0, 0.0};
    C_R_C_inv_f[i] = (complex_number_t) {0.0, 0.0};
    for (j = 0 ; j<m ; j++)
    {
      C_R_C[i][j] = (complex_number_t) {0.0, 0.0};
      C_R_C_inv[j][i] = (complex_number_t) {0.0, 0.0};
    }
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
		complex_number_t POLYBENCH_2D(qrN_Q, N, N, n, n),
		complex_number_t POLYBENCH_2D(qrN_R, N, N, n, n),
		complex_number_t POLYBENCH_2D(qrN_R_inv, N, N, n, n),
		complex_number_t POLYBENCH_2D(qrM_Q, M, M, m, n),
		complex_number_t POLYBENCH_2D(qrM_R, M, M, m, n),
		complex_number_t POLYBENCH_2D(qrM_R_inv, M, M, m, n),

    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_2D(R, N, N, n, n),
    complex_number_t POLYBENCH_2D(Rinv, N, N, n, n),
    complex_number_t POLYBENCH_2D(C, N, M, n, m),
    complex_number_t POLYBENCH_2D(R_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(Rinv_C, N, M, n, m),
    complex_number_t POLYBENCH_2D(C_R_C, M, M, m, m),
    complex_number_t POLYBENCH_2D(C_R_C_inv, M, M, m, m),
    complex_number_t POLYBENCH_1D(C_R_C_inv_f, M, m),
    complex_number_t POLYBENCH_1D(f, M, m),
    complex_number_t POLYBENCH_1D(w, N, n))
{
  int i = 0, j = 0, k = 0;
  DATA_TYPE r = 0.0;
  complex_number_t qr_c = {0.0, 0.0}, qr_s = {0.0, 0.0};
  complex_number_t Rik = {0.0, 0.0}, Rjk = {0.0, 0.0};
  complex_number_t Qik = {0.0, 0.0}, Qjk = {0.0, 0.0};
  complex_number_t Qki = {0.0, 0.0}, Qkj = {0.0, 0.0};
  complex_number_t tmp = {0.0, 0.0};

#pragma scop
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      R[i][j].r = (x[i].r * x[j].r - x[i].i * (-x[j].i)) / N;
      R[i][j].i = (x[i].i * x[j].r + x[i].r * (-x[j].i)) / N;
    }
  }

  ////////////////////////////////////////////////////////
  //INVERSE R INTO Rinv
  
  //qrN_R = R, qrN_Q = Id
  for (i = 0 ; i < N ; i++)
  {
    for (j = 0 ; j < N ; j++)
    {
      qrN_R[i][j] = R[i][j];
      qrN_Q[i][j].r = 0.0;
      qrN_Q[i][j].i = 0.0;
    }
    qrN_Q[i][i].r = 1.0;
  }

  for (j = 0; j < N-1; j++)
  {
    for (i = j+1; i < N; i++)
    {
      r = copysign(1.0, qrN_R[i][j].r ) * SQRT_FUN(
          (qrN_R[j][j].r*qrN_R[j][j].r + qrN_R[j][j].i*qrN_R[j][j].i) +
          (qrN_R[i][j].r*qrN_R[i][j].r + qrN_R[i][j].i*qrN_R[i][j].i));

      qr_c.r = qrN_R[j][j].r / r;
      qr_c.i = qrN_R[j][j].i / r;
      qr_s.r = -qrN_R[i][j].r / r;
      qr_s.i = -qrN_R[i][j].i / r;

      // Apply the givens rotation:
      for (k = j; k < N ; k++) {
          Rjk = R[j][k];
          Rik = R[i][k];

          qrN_R[j][k].r = (qr_c.r * Rjk.r - qr_c.i * Rjk.i) + ((-qr_s.r) * Rik.r - (-qr_s.i) * Rik.i);
          qrN_R[j][k].i = (qr_c.i * Rjk.r + qr_c.r * Rjk.i) + ((-qr_s.i) * Rik.r + (-qr_s.r) * Rik.i);

          qrN_R[i][k].r = (qr_c.r * Rjk.r - qr_c.i * Rjk.i) + (qr_s.r * Rik.r - qr_s.i * Rik.i);
          qrN_R[i][k].i = (qr_c.i * Rjk.r + qr_c.r * Rjk.i) + (qr_s.i * Rik.r + qr_s.r * Rik.i);
      }

      //Save G in qrN_Q by applying G^t to qrN_Q
      for (k = 0; k < N ; k++) {
          Qki = qrN_Q[k][i];
          Qkj = qrN_Q[k][j];
          /*qrN_Q[j][k].r = Qkj * qr_s + Qki * qr_c;*/
          qrN_Q[j][k].r = (qr_c.r * Qkj.r - qr_c.i * Qkj.i) + ((-qr_s.r) * Qki.r - (-qr_s.i) * Qki.i);
          qrN_Q[j][k].i = (qr_c.i * Qkj.r + qr_c.r * Qkj.i) + ((-qr_s.i) * Qki.r + (-qr_s.r) * Qki.i);
          
          /*qrN_Q[k][i].r = Qkj * qr_s + Qki * qr_c;*/
          qrN_Q[k][i].r = (qr_c.r * Qki.r - qr_c.i * Qki.i) + (qr_s.r * Qkj.r - qr_s.i * Qkj.i);
          qrN_Q[k][i].i = (qr_c.i * Qki.r + qr_c.r * Qki.i) + (qr_s.i * Qkj.r + qr_s.r * Qkj.i);
      }
    }
  }

  //Inverse qrN_R into qrN_R_inv
  for(i=N-1 ; i >= 0  ; i--)
  {
      qrN_R_inv[i][i].r = qrN_R[i][i].r / (qrN_R[i][i].r * qrN_R[i][i].r + qrN_R[i][i].i * qrN_R[i][i].i);
      qrN_R_inv[i][i].i = -qrN_R[i][i].i / (qrN_R[i][i].r * qrN_R[i][i].r + qrN_R[i][i].i * qrN_R[i][i].i);
      for(j=i+1; j<N; j++)
      {
          for(k=i+1; k<N; k++)
          {
              qrN_R_inv[i][j].r += qrN_R[i][k].r * qrN_R_inv[k][j].r ;
              qrN_R_inv[i][j].i += qrN_R[i][k].i * qrN_R_inv[k][j].i ;
          }
          tmp = qrN_R_inv[i][j];
          qrN_R_inv[i][j].r = (tmp.r * (-qrN_R_inv[i][i].r)) - (tmp.i * (-qrN_R_inv[i][i].i));
          qrN_R_inv[i][j].i = (tmp.r * (-qrN_R_inv[i][i].i)) + (tmp.i * (-qrN_R_inv[i][i].r));
      }
  }

  //R^{-1} = qrN_R^{-1} * qrN_Q^H
  for (i=0; i<N ;i++)
  {
    for (j=0; j<N ;j++)
    {
      Rinv[i][j].r = 0.0;
      Rinv[i][j].i = 0.0;
      for (k=0; k<N ;k++)
      {
        Rinv[i][j].r += qrN_R_inv[i][k].r * qrN_Q[j][k].r - (qrN_R_inv[i][k].i * (-qrN_Q[j][k].i));
        Rinv[i][j].i += qrN_R_inv[i][k].i * qrN_Q[j][k].r + qrN_R_inv[i][k].r * qrN_Q[j][k].i;
      }
    }
  }

  //END INVERSE
  ////////////////////////////////////////////////////////

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

  for (i=0 ; i<M ; i++)
  {
    for (j=0 ; j<M ; j++)
    {
      C_R_C[i][j].r = 0.0;
      C_R_C[i][j].i = 0.0;
    }
  }
  for (k=0 ; k<N ; k++)
  {
    for (i=0 ; i<M ; i++)
    {
      for (j=0 ; j<M ; j++)
      {
        C_R_C[i][j].r += (C[k][i].r * R_C[k][j].r) - ((-C[k][i].i) * R_C[k][j].i);
        C_R_C[i][j].i += ((-C[k][i].i) * R_C[k][j].r) + (C[k][i].r * R_C[k][j].i);
      }
    }
  }

  //INVERSE C_R_C INTO C_R_C_inv
  /////////////////////////////////////////////////////////////////////////////

  //qrM_R = R, qrM_Q = Id
  for (i = 0 ; i < M ; i++)
  {
    for (j = 0 ; j < M ; j++)
    {
      qrM_R[i][j] = C_R_C[i][j];
      qrM_Q[i][j].r = 0.0;
      qrM_Q[i][j].i = 0.0;
    }
    qrM_Q[i][i].r = 1.0;
  }

  for (j = 0; j < M-1; j++)
  {
    for (i = j+1; i < M; i++)
    {
      r = copysign(1.0, qrM_R[i][j].r ) * SQRT_FUN(
          (qrM_R[j][j].r*qrM_R[j][j].r + qrM_R[j][j].i*qrM_R[j][j].i) +
          (qrM_R[i][j].r*qrM_R[i][j].r + qrM_R[i][j].i*qrM_R[i][j].i));

      qr_c.r = qrM_R[j][j].r / r;
      qr_c.i = qrM_R[j][j].i / r;
      qr_s.r = -qrM_R[i][j].r / r;
      qr_s.i = -qrM_R[i][j].i / r;

      // Apply the givens rotation:
      for (k = j; k < M ; k++) {
          Rjk = C_R_C[j][k];
          Rik = C_R_C[i][k];

          qrM_R[j][k].r = (qr_c.r * Rjk.r - qr_c.i * Rjk.i) + ((-qr_s.r) * Rik.r - (-qr_s.i) * Rik.i);
          qrM_R[j][k].i = (qr_c.i * Rjk.r + qr_c.r * Rjk.i) + ((-qr_s.i) * Rik.r + (-qr_s.r) * Rik.i);

          qrM_R[i][k].r = (qr_c.r * Rjk.r - qr_c.i * Rjk.i) + (qr_s.r * Rik.r - qr_s.i * Rik.i);
          qrM_R[i][k].i = (qr_c.i * Rjk.r + qr_c.r * Rjk.i) + (qr_s.i * Rik.r + qr_s.r * Rik.i);
      }

      //Save G in qrM_Q by applying G^t to qrM_Q
      for (k = 0; k < M ; k++) {
          Qki = qrM_Q[k][i];
          Qkj = qrM_Q[k][j];
          /*qrM_Q[j][k].r = Qkj * qr_s + Qki * qr_c;*/
          qrM_Q[j][k].r = (qr_c.r * Qkj.r - qr_c.i * Qkj.i) + ((-qr_s.r) * Qki.r - (-qr_s.i) * Qki.i);
          qrM_Q[j][k].i = (qr_c.i * Qkj.r + qr_c.r * Qkj.i) + ((-qr_s.i) * Qki.r + (-qr_s.r) * Qki.i);
          
          /*qrM_Q[k][i].r = Qkj * qr_s + Qki * qr_c;*/
          qrM_Q[k][i].r = (qr_c.r * Qki.r - qr_c.i * Qki.i) + (qr_s.r * Qkj.r - qr_s.i * Qkj.i);
          qrM_Q[k][i].i = (qr_c.i * Qki.r + qr_c.r * Qki.i) + (qr_s.i * Qkj.r + qr_s.r * Qkj.i);
      }
    }
  }

  //Inverse qrM_R into qrM_R_inv
  for(i=M-1 ; i >= 0  ; i--)
  {
      qrM_R_inv[i][i].r = qrM_R[i][i].r / (qrM_R[i][i].r * qrM_R[i][i].r + qrM_R[i][i].i * qrM_R[i][i].i);
      qrM_R_inv[i][i].i = -qrM_R[i][i].i / (qrM_R[i][i].r * qrM_R[i][i].r + qrM_R[i][i].i * qrM_R[i][i].i);
      for(j=i+1; j<M; j++)
      {
          for(k=i+1; k<M; k++)
          {
              qrM_R_inv[i][j].r += qrM_R[i][k].r * qrM_R_inv[k][j].r ;
              qrM_R_inv[i][j].i += qrM_R[i][k].i * qrM_R_inv[k][j].i ;
          }
          tmp = qrM_R_inv[i][j];
          qrM_R_inv[i][j].r = (tmp.r * (-qrM_R_inv[i][i].r)) - (tmp.i * (-qrM_R_inv[i][i].i));
          qrM_R_inv[i][j].i = (tmp.r * (-qrM_R_inv[i][i].i)) + (tmp.i * (-qrM_R_inv[i][i].r));
      }
  }

  //R^{-1} = qrM_R^{-1} * qrM_Q^H
  for (i=0; i<M ;i++)
  {
    for (j=0; j<M ;j++)
    {
      C_R_C_inv[i][j].r = 0.0;
      C_R_C_inv[i][j].i = 0.0;
      for (k=0; k<M ;k++)
      {
        C_R_C_inv[i][j].r += qrM_R_inv[i][k].r * qrM_Q[j][k].r - (qrM_R_inv[i][k].i * (-qrM_Q[j][k].i));
        C_R_C_inv[i][j].i += qrM_R_inv[i][k].i * qrM_Q[j][k].r + qrM_R_inv[i][k].r * qrM_Q[j][k].i;
      }
    }
  }

  //END INVERSE
  /////////////////////////////////////////////////////////////////////////////

  for (i=0 ; i<M ; i++)
  {
    C_R_C_inv_f[i].r = 0.0;
    C_R_C_inv_f[i].i = 0.0;

    for (j=0 ; j<M ; j++)
    {
      C_R_C_inv_f[i].r += (C_R_C_inv[i][j].r * f[j].r) - (C_R_C_inv[i][j].i * f[j].i);
      C_R_C_inv_f[i].i += (C_R_C_inv[i][j].r * f[j].i) + (C_R_C_inv[i][j].i * f[j].r);
    }
  }


  for (i=0 ; i<N ; i++)
  {
    w[i].r = 0.0;
    w[i].i = 0.0;

    for (j=0 ; j<M ; j++)
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
  POLYBENCH_2D_ARRAY_DECL(qrN_Q, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(qrN_R, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(qrN_R_inv, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(qrM_Q, complex_number_t, M, M, m, m);
  POLYBENCH_2D_ARRAY_DECL(qrM_R, complex_number_t, M, M, m, m);
  POLYBENCH_2D_ARRAY_DECL(qrM_R_inv, complex_number_t, M, M, m, m);

  POLYBENCH_1D_ARRAY_DECL(x, complex_number_t, N, n);

  POLYBENCH_2D_ARRAY_DECL(R, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Rinv, complex_number_t, N, N, n, n);
  
  POLYBENCH_2D_ARRAY_DECL(C, complex_number_t, N, M, n, m);
  
  POLYBENCH_2D_ARRAY_DECL(R_C, complex_number_t, N, M, n, m);
  POLYBENCH_2D_ARRAY_DECL(Rinv_C, complex_number_t, N, M, n, m);
  POLYBENCH_2D_ARRAY_DECL(C_R_Ca, complex_number_t, M, M, m, m);
  POLYBENCH_2D_ARRAY_DECL(C_R_C_inv, complex_number_t, M, M, m, m);
  
  POLYBENCH_1D_ARRAY_DECL(f, complex_number_t, M, m);
  POLYBENCH_1D_ARRAY_DECL(C_R_C_inv_f, complex_number_t, M, m);
  POLYBENCH_1D_ARRAY_DECL(w, complex_number_t, N, n);

  /* Initialize array(s). */
  init_array (n, m,
      POLYBENCH_ARRAY(qrN_Q),
      POLYBENCH_ARRAY(qrN_R),
      POLYBENCH_ARRAY(qrN_R_inv),
      POLYBENCH_ARRAY(qrM_Q),
      POLYBENCH_ARRAY(qrM_R),
      POLYBENCH_ARRAY(qrM_R_inv),
      POLYBENCH_ARRAY(x),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(Rinv),
      POLYBENCH_ARRAY(C),
      POLYBENCH_ARRAY(R_C),
      POLYBENCH_ARRAY(Rinv_C),
      POLYBENCH_ARRAY(C_R_Ca),
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
      POLYBENCH_ARRAY(qrN_Q),
      POLYBENCH_ARRAY(qrN_R),
      POLYBENCH_ARRAY(qrN_R_inv),
      POLYBENCH_ARRAY(qrM_Q),
      POLYBENCH_ARRAY(qrM_R),
      POLYBENCH_ARRAY(qrM_R_inv),
      POLYBENCH_ARRAY(x),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(Rinv),
      POLYBENCH_ARRAY(C),
      POLYBENCH_ARRAY(R_C),
      POLYBENCH_ARRAY(Rinv_C),
      POLYBENCH_ARRAY(C_R_Ca),
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
  POLYBENCH_FREE_ARRAY(qrN_Q);
  POLYBENCH_FREE_ARRAY(qrN_R);
  POLYBENCH_FREE_ARRAY(qrN_R_inv);
  POLYBENCH_FREE_ARRAY(qrM_Q);
  POLYBENCH_FREE_ARRAY(qrM_R);
  POLYBENCH_FREE_ARRAY(qrM_R_inv);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Rinv);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(R_C);
  POLYBENCH_FREE_ARRAY(Rinv_C);
  POLYBENCH_FREE_ARRAY(C_R_Ca);
  POLYBENCH_FREE_ARRAY(C_R_C_inv);
  POLYBENCH_FREE_ARRAY(C_R_C_inv_f);
  POLYBENCH_FREE_ARRAY(f);
  POLYBENCH_FREE_ARRAY(w);

  return 0;
}
