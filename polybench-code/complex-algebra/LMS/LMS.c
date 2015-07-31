#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/*#define MINI_DATASET*/
/* Include benchmark-specific header. */
#include "LMS.h"


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
    complex_number_t POLYBENCH_1D(wq, N, n),
    complex_number_t POLYBENCH_2D(B, N, N-M, n, n-m),
    complex_number_t * yc,
    complex_number_t POLYBENCH_1D(z, N-M, n-m),
    complex_number_t POLYBENCH_1D(wa, N-M, n-m),
    complex_number_t * yp,
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0;

  *yc = (complex_number_t) {0.0, 0.0};
  *yp = (complex_number_t) {0.0, 0.0};

  for(i = 0 ; i<n ; i++)
  {
    wq[i].r = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;
    wq[i].i = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;

    w[i].r = 0.0;
    w[i].i = 0.0;

    for (j = 0 ; j<(n-m) ; j++)
    {
      B[i][j].r = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * 10.0;
      B[i][j].i = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * 10.0;
    }
  }
  
  update_array(n,x);

  
  //1. w_a(0) = 0
  for (i = 0 ; i<(n-m) ; i++)
  {
    wa[i].r = 0.0;
    wa[i].i = 0.0;
  }
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
void kernel_LMS(int n, int m,
    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_1D(wq, N, n),
    complex_number_t POLYBENCH_2D(B, N, N-m, n, n-m),
    complex_number_t * yc,
    complex_number_t POLYBENCH_1D(z, N-M, n-m),
    complex_number_t POLYBENCH_1D(wa, N-M, n-m),
    complex_number_t * yp,
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0;

  complex_number_t waz = {0,0};

  complex_number_t Bwa = {0,0};

#pragma scop
  //2. y_c(k) = w_q^{H} * x(k)
  for( i = 0 ; i<n ; i++)
  {
    yc->r += (wq[i].r * x[i].r) - ((-wq[i].i) * x[i].i);
    yc->i += ((-wq[i].i) * x[i].r) + (wq[i].r * x[i].i);
    
  }

  //3. z(k) = B^{H} * x(k)
  for (j = 0 ; j<(n-m) ; j++)
  {
    z[j].r = 0.0;
    z[j].i = 0.0;

    for(i = 0 ; i<n ; i++)
    {
      z[j].r += (B[i][j].r * x[i].r) - ((-B[i][j].i) * x[i].i);
      z[j].i += ((-B[i][j].i) * x[i].r) + (B[i][j].r * x[i].i);
    }
  }

  //4. y_p(k) = y_c(k) - w_a^H(k-1) * z(k)
  for (i = 0 ; i < (n-m) ; i++)
  {
    waz.r += (wa[i].r * z[i].r) - ((-wa[i].i) * z[i].i);
    waz.i += (-wa[i].i)*z[i].r + wa[i].r*z[i].i;
  }
  yp[0].r = yc[0].r - waz.r;
  yp[0].i = yc[0].i - waz.i;
  
  //5. w_a(k) = w_a (k-1) + ALPHA * z(k) * y_p^{*}(k)
  for (i = 0 ; i < (n-m) ; i++)
  {
    wa[i].r = wa[i].r + ALPHA * ((z[i].r * yp[0].r) - (z[i].i * (-yp[0].i)));
    wa[i].i = wa[i].i + ALPHA * ((z[i].i * (-yp[0].i)) - (z[i].r * yp[0].r));
  }

  //6. w(k) = w_q - B * w_a(k)
  for (i = 0 ; i < n ; i++)
  {
    Bwa.r = 0.0;
    Bwa.i = 0.0;
    for (j = 0 ; j < (n-m) ; j++)
    {
      Bwa.r += B[i][j].r * wa[j].r - B[i][j].i * wa[j].i;
      Bwa.i += B[i][j].r * wa[j].i - B[i][j].i * wa[j].r;
    }
    w[i].r = wq[i].r - Bwa.r;
    w[i].i = wq[i].i - Bwa.i;
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
  POLYBENCH_1D_ARRAY_DECL(wq, complex_number_t, N, n);
  POLYBENCH_2D_ARRAY_DECL(B, complex_number_t, N, N-m, n, n-m);
  complex_number_t yc;
  POLYBENCH_1D_ARRAY_DECL(z, complex_number_t, N-M, n-m);
  POLYBENCH_1D_ARRAY_DECL(wa, complex_number_t, N-M, n-m);
  complex_number_t yp;
  POLYBENCH_1D_ARRAY_DECL(w, complex_number_t, N, n);

  /* Initialize array(s). */
  init_array (n,
      m,
      POLYBENCH_ARRAY(x),
      POLYBENCH_ARRAY(wq),
      POLYBENCH_ARRAY(B),
      &yc,
      POLYBENCH_ARRAY(z),
      POLYBENCH_ARRAY(wa),
      &yp,
      POLYBENCH_ARRAY(w))

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  for (i = 0 ; i < NB_RUN ; i++)
  {
    //Just comment this function call if you don't want to 
    //update/modify the input data (x vector) of the radio.
    update_array(n,
        POLYBENCH_ARRAY(x));

    kernel_LMS (n, m,
        POLYBENCH_ARRAY(x),
        POLYBENCH_ARRAY(wq),
        POLYBENCH_ARRAY(B),
        &yc,
        POLYBENCH_ARRAY(z),
        POLYBENCH_ARRAY(wa),
        &yp,
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
  POLYBENCH_FREE_ARRAY(wq);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(z);
  POLYBENCH_FREE_ARRAY(wa);
  POLYBENCH_FREE_ARRAY(w);

  return 0;
}
