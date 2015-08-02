#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/*
 * The MVDR algorithm was taken from the following paper :
 * "Application of Adaptive Beamforming Techniques to HF Radar"
 * written by Peter Vouras and Brian Freburger.
 */

/* Include polybench common header. */
#include <polybench.h>

/*#define MINI_DATASET*/
/* Include benchmark-specific header. */
#include "MVDR.h"


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
void init_array(int n,
		complex_number_t POLYBENCH_2D(qr_Q, N, N, n, n),
		complex_number_t POLYBENCH_2D(qr_R, N, N, n, n),
		complex_number_t POLYBENCH_2D(qr_R_inv, N, N, n, n),
    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_2D(R, N, N, n, n),
    complex_number_t POLYBENCH_2D(Rinv, N, N, n, n),
    complex_number_t POLYBENCH_1D(s, N, n),
    complex_number_t POLYBENCH_1D(Rinv_s, N, n),
    complex_number_t POLYBENCH_1D(R_s, N, n),
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0;

  for(i = 0 ; i<n ; i++)
  {
    x[i] = (complex_number_t) {0.0, 0.0};
    s[i] = (complex_number_t) {0.0, 0.0};
    Rinv_s[i] = (complex_number_t) {0.0, 0.0};
    R_s[i] = (complex_number_t) {0.0, 0.0};
    w[i] = (complex_number_t) {0.0, 0.0};
    for (j = 0 ; j<n ; j++)
    {
      R[i][j] = (complex_number_t) {0.0, 0.0};
      Rinv[i][j] = (complex_number_t) {0.0, 0.0};
      
      qr_Q[i][j].r = 0.0;
      qr_Q[i][j].i = 0.0;

      qr_R[i][j].r = 0.0;
      qr_R[i][j].i = 0.0;

      qr_R_inv[i][j].r = 0.0;
      qr_R_inv[i][j].i = 0.0;

    }
  }
  
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
void kernel_MVDR(int n,
		complex_number_t POLYBENCH_2D(qr_Q, N, N, n, n),
		complex_number_t POLYBENCH_2D(qr_R, N, N, n, n),
		complex_number_t POLYBENCH_2D(qr_R_inv, N, N, n, n),
    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_2D(R, N, N, n, n),
    complex_number_t POLYBENCH_2D(Rinv, N, N, n, n),
    complex_number_t POLYBENCH_1D(s, N, n),
    complex_number_t POLYBENCH_1D(Rinv_s, N, n),
    complex_number_t POLYBENCH_1D(R_s, N, n),
    complex_number_t POLYBENCH_1D(w, N, n))
{
  int i = 0, j = 0, k = 0;
  complex_number_t sH_R_s = {0.0, 0.0};
  double norm = 0.0;

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

  //////////////////////////////////////
  //INVERSE R into Rinv

  //qr_R = R, qr_Q = Id
  for (i = 0 ; i < N ; i++)
  {
    for (j = 0 ; j < N ; j++)
    {
      qr_R[i][j] = R[i][j];
      qr_Q[i][j].r = 0.0;
      qr_Q[i][j].i = 0.0;
    }
    qr_Q[i][i].r = 1.0;
  }

  for (j = 0; j < N-1; j++)
  {
    for (i = j+1; i < N; i++)
    {
      r = copysign(1.0, qr_R[i][j].r ) * SQRT_FUN(
          (qr_R[j][j].r*qr_R[j][j].r + qr_R[j][j].i*qr_R[j][j].i) +
          (qr_R[i][j].r*qr_R[i][j].r + qr_R[i][j].i*qr_R[i][j].i));

      qr_c.r = qr_R[j][j].r / r;
      qr_c.i = qr_R[j][j].i / r;
      qr_s.r = -qr_R[i][j].r / r;
      qr_s.i = -qr_R[i][j].i / r;

      // Apply the givens rotation:
      for (k = j; k < N ; k++) {
          Rjk = R[j][k];
          Rik = R[i][k];

          qr_R[j][k].r = (qr_c.r * Rjk.r - qr_c.i * Rjk.i) + ((-qr_s.r) * Rik.r - (-qr_s.i) * Rik.i);
          qr_R[j][k].i = (qr_c.i * Rjk.r + qr_c.r * Rjk.i) + ((-qr_s.i) * Rik.r + (-qr_s.r) * Rik.i);

          qr_R[i][k].r = (qr_c.r * Rjk.r - qr_c.i * Rjk.i) + (qr_s.r * Rik.r - qr_s.i * Rik.i);
          qr_R[i][k].i = (qr_c.i * Rjk.r + qr_c.r * Rjk.i) + (qr_s.i * Rik.r + qr_s.r * Rik.i);
      }

      //Save G in qr_Q by applying G^t to qr_Q
      for (k = 0; k < N ; k++) {
          Qki = qr_Q[k][i];
          Qkj = qr_Q[k][j];
          /*qr_Q[j][k].r = Qkj * qr_s + Qki * qr_c;*/
          qr_Q[j][k].r = (qr_c.r * Qkj.r - qr_c.i * Qkj.i) + ((-qr_s.r) * Qki.r - (-qr_s.i) * Qki.i);
          qr_Q[j][k].i = (qr_c.i * Qkj.r + qr_c.r * Qkj.i) + ((-qr_s.i) * Qki.r + (-qr_s.r) * Qki.i);
          
          /*qr_Q[k][i].r = Qkj * qr_s + Qki * qr_c;*/
          qr_Q[k][i].r = (qr_c.r * Qki.r - qr_c.i * Qki.i) + (qr_s.r * Qkj.r - qr_s.i * Qkj.i);
          qr_Q[k][i].i = (qr_c.i * Qki.r + qr_c.r * Qki.i) + (qr_s.i * Qkj.r + qr_s.r * Qkj.i);
      }
    }
  }

  //Inverse qr_R into qr_R_inv
  for(i=N-1 ; i >= 0  ; i--)
  {
      qr_R_inv[i][i].r = qr_R[i][i].r / (qr_R[i][i].r * qr_R[i][i].r + qr_R[i][i].i * qr_R[i][i].i);
      qr_R_inv[i][i].i = -qr_R[i][i].i / (qr_R[i][i].r * qr_R[i][i].r + qr_R[i][i].i * qr_R[i][i].i);
      for(j=i+1; j<N; j++)
      {
          for(k=i+1; k<N; k++)
          {
              qr_R_inv[i][j].r += qr_R[i][k].r * qr_R_inv[k][j].r ;
              qr_R_inv[i][j].i += qr_R[i][k].i * qr_R_inv[k][j].i ;
          }
          tmp = qr_R_inv[i][j];
          qr_R_inv[i][j].r = (tmp.r * (-qr_R_inv[i][i].r)) - (tmp.i * (-qr_R_inv[i][i].i));
          qr_R_inv[i][j].i = (tmp.r * (-qr_R_inv[i][i].i)) + (tmp.i * (-qr_R_inv[i][i].r));
      }
  }

  //R^{-1} = qr_R^{-1} * qr_Q^H
  for (i=0; i<N ;i++)
  {
    for (j=0; j<N ;j++)
    {
      Rinv[i][j].r = 0.0;
      Rinv[i][j].i = 0.0;
      for (k=0; k<N ;k++)
      {
        Rinv[i][j].r += qr_R_inv[i][k].r * qr_Q[j][k].r - (qr_R_inv[i][k].i * (-qr_Q[j][k].i));
        Rinv[i][j].i += qr_R_inv[i][k].i * qr_Q[j][k].r + qr_R_inv[i][k].r * qr_Q[j][k].i;
      }
    }
  }

  //END INVERSE
  ///////////////////////////////////////

  //Rinv_s = Rinv * s
  for (i = 0; i < N ; i++)
  {
    Rinv_s[i].r = 0.0;
    Rinv_s[i].i = 0.0;
    for (j= 0; j < N; j++)
    {
      Rinv_s[i].r += Rinv[i][j].r * s[j].r - Rinv[i][j].i * s[j].i;
      Rinv_s[i].i += Rinv[i][j].i * s[j].r + Rinv[i][j].r * s[j].i;
    }
  }

  //R_s = R * s
  for (i = 0; i < N ; i++)
  {
    R_s[i].r = 0.0;
    R_s[i].i = 0.0;
    for (j= 0; j < N; j++)
    {
      R_s[i].r += R[i][j].r * s[j].r - R[i][j].i * s[j].i;
      R_s[i].i += R[i][j].i * s[j].r + R[i][j].r * s[j].i;
    }
  }

  //sH_R_s = s^{H} * R * s
  for (i = 0; i < N ; i++)
  {
      sH_R_s.r += s[i].r * R_s[i].r - s[i].i * R_s[i].i;
      sH_R_s.i += s[i].i * R_s[i].r + s[i].r * R_s[i].i;
  }

  norm = sH_R_s.r * sH_R_s.r + sH_R_s.i * sH_R_s.i;
  for (i = 0 ; i < N ; i++)
  {
    w[i].r = ((R_s[i].r * sH_R_s.r) + (R_s[i].i * sH_R_s.i)) / norm;
    w[i].i = ((R_s[i].i * sH_R_s.r) - (R_s[i].r * sH_R_s.i)) / norm;
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  
  
  /* Variable declaration/allocation. */

  POLYBENCH_2D_ARRAY_DECL(qr_Q, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(qr_R, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(qr_R_inv, complex_number_t, N, N, n, n);

  POLYBENCH_1D_ARRAY_DECL(x, complex_number_t, N, n);
  POLYBENCH_2D_ARRAY_DECL(R, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Rinv, complex_number_t, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(s, complex_number_t, N, n);
  POLYBENCH_1D_ARRAY_DECL(Rinv_s, complex_number_t, N, n);
  POLYBENCH_1D_ARRAY_DECL(R_s, complex_number_t, N, n);
  POLYBENCH_1D_ARRAY_DECL(w, complex_number_t, N, n);

  /* Initialize array(s). */
  init_array (n,
      POLYBENCH_ARRAY(qr_Q),
      POLYBENCH_ARRAY(qr_R),
      POLYBENCH_ARRAY(qr_R_inv),
      POLYBENCH_ARRAY(x),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(Rinv),
      POLYBENCH_ARRAY(s),
      POLYBENCH_ARRAY(Rinv_s),
      POLYBENCH_ARRAY(R_s),
      POLYBENCH_ARRAY(w));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  //Just comment this function call if you don't want to 
  //update/modify the input data (x vector) of the radio.
  update_array(n,
      POLYBENCH_ARRAY(x));

  kernel_MVDR (n,
      POLYBENCH_ARRAY(qr_Q),
      POLYBENCH_ARRAY(qr_R),
      POLYBENCH_ARRAY(qr_R_inv),
      POLYBENCH_ARRAY(x),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(Rinv),
      POLYBENCH_ARRAY(s),
      POLYBENCH_ARRAY(Rinv_s),
      POLYBENCH_ARRAY(R_s),
      POLYBENCH_ARRAY(w));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(qr_Q);
  POLYBENCH_FREE_ARRAY(qr_R);
  POLYBENCH_FREE_ARRAY(qr_R_inv);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Rinv);
  POLYBENCH_FREE_ARRAY(s);
  POLYBENCH_FREE_ARRAY(Rinv_s);
  POLYBENCH_FREE_ARRAY(R_s);
  POLYBENCH_FREE_ARRAY(w);

  return 0;
}
