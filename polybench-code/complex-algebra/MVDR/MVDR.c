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
    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_2D(R, N, N, n, n),
    complex_number_t POLYBENCH_2D(Rinv, N, N, n, n),
    complex_number_t POLYBENCH_1D(s, N, n),
    complex_number_t POLYBENCH_1D(Rinv_s, N, n),
    complex_number_t POLYBENCH_1D(R_s, N, n),
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0;
  complex_number_t sH_R_s = {0.0, 0.0};
  double tmp = 0.0;

#pragma scop
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      R[i][j].r = (x[i].r * x[j].r - x[i].i * (-x[j].i)) / N;
      R[i][j].i = (x[i].i * x[j].r + x[i].r * (-x[j].i)) / N;
    }
  }

  //INVERSE R into Rinv //TODO

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

  tmp = sH_R_s.r * sH_R_s.r + sH_R_s.i * sH_R_s.i;
  for (i = 0 ; i < N ; i++)
  {
    w[i].r = ((R_s[i].r * sH_R_s.r) + (R_s[i].i * sH_R_s.i)) / tmp;
    w[i].i = ((R_s[i].i * sH_R_s.r) - (R_s[i].r * sH_R_s.i)) / tmp;
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  unsigned int i = 0;

  /* Retrieve problem size. */
  int n = N;
  
  
  /* Variable declaration/allocation. */
  POLYBENCH_1D_ARRAY_DECL(x, complex_number_t, N, n);
  POLYBENCH_2D_ARRAY_DECL(R, complex_number_t, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(Rinv, complex_number_t, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(s, complex_number_t, N, n);
  POLYBENCH_1D_ARRAY_DECL(Rinv_s, complex_number_t, N, n);
  POLYBENCH_1D_ARRAY_DECL(R_s, complex_number_t, N, n);
  POLYBENCH_1D_ARRAY_DECL(w, complex_number_t, N, n);

  /* Initialize array(s). */
  init_array (n,
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
  for (i = 0 ; i < NB_RUN ; i++)
  {
    //Just comment this function call if you don't want to 
    //update/modify the input data (x vector) of the radio.
    update_array(n,
        POLYBENCH_ARRAY(x));

    kernel_MVDR (n,
        POLYBENCH_ARRAY(x),
        POLYBENCH_ARRAY(R),
        POLYBENCH_ARRAY(Rinv),
        POLYBENCH_ARRAY(s),
        POLYBENCH_ARRAY(Rinv_s),
        POLYBENCH_ARRAY(R_s),
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
  POLYBENCH_FREE_ARRAY(s);
  POLYBENCH_FREE_ARRAY(Rinv_s);
  POLYBENCH_FREE_ARRAY(R_s);
  POLYBENCH_FREE_ARRAY(w);

  return 0;
}
