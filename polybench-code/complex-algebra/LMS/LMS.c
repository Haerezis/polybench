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
    DATA_TYPE POLYBENCH_1D(x_r, N, n),
    DATA_TYPE POLYBENCH_1D(x_i, N, n))
{
  unsigned int i = 0;

  for(i = 0 ; i<n ; i++)
  {
    x_r[i] = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;
    x_i[i] = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;
  }
}

/* Array initialization. */
static
void init_array(int n, int m,
    DATA_TYPE POLYBENCH_1D(x_r, N, n),
    DATA_TYPE POLYBENCH_1D(x_i, N, n),
    DATA_TYPE POLYBENCH_1D(wq_r, N, n),
    DATA_TYPE POLYBENCH_1D(wq_i, N, n),
    DATA_TYPE POLYBENCH_2D(B_r, N, N-M, n, n-m),
    DATA_TYPE POLYBENCH_2D(B_i, N, N-M, n, n-m),
    DATA_TYPE * yc_r,
    DATA_TYPE * yc_i,
    DATA_TYPE POLYBENCH_1D(z_r, N-M, n-m),
    DATA_TYPE POLYBENCH_1D(z_i, N-M, n-m),
    DATA_TYPE POLYBENCH_1D(wa_r, N-M, n-m),
    DATA_TYPE POLYBENCH_1D(wa_i, N-M, n-m),
    DATA_TYPE * yp_r,
    DATA_TYPE * yp_i,
    DATA_TYPE POLYBENCH_1D(w_r, N, n),
    DATA_TYPE POLYBENCH_1D(w_i, N, n))
{
  unsigned int i = 0, j = 0;

  *yc_r = 0.0;
  *yc_i = 0.0;
  *yp_r = 0.0;
  *yp_i = 0.0;

  for(i = 0 ; i<n ; i++)
  {
    wq_r[i] = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;
    wq_i[i] = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;

    w_r[i] = 0.0;
    w_i[i] = 0.0;

    for (j = 0 ; j<(n-m) ; j++)
    {
      B_r[i][j] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      B_i[i][j] = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
    }
  }
  
  update_array(n,
      x_r,
      x_i);

  
  //1. w_a(0) = 0
  for (j = 0 ; j<(n-m) ; j++)
  {
    wa_r[i] = 0.0;
    wa_i[i] = 0.0;
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
    DATA_TYPE POLYBENCH_1D(w_r, N, n),
    DATA_TYPE POLYBENCH_1D(w_i, N, n))
{
  unsigned int i = 0;

  fprintf(stderr, "w : (");
  for (i = 0 ; i<n ; i++)
  {
    fprintf (stderr, DATA_PRINTF_MODIFIER, w_r[i]);
    if (w_i[i] >= 0.0) fprintf (stderr, "+");
    fprintf (stderr, DATA_PRINTF_MODIFIER, w_r[i]);
    if (i<(n-1)) fprintf (stderr, ", ");
  }
  fprintf (stderr, ")\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_LMS_GSC(int n, int m,
    DATA_TYPE POLYBENCH_1D(x_r, N, n),
    DATA_TYPE POLYBENCH_1D(x_i, N, n),
    DATA_TYPE POLYBENCH_1D(wq_r, N, n),
    DATA_TYPE POLYBENCH_1D(wq_i, N, n),
    DATA_TYPE POLYBENCH_2D(B_r, N, N-m, n, n-m),
    DATA_TYPE POLYBENCH_2D(B_i, N, N-m, n, n-m),
    DATA_TYPE * yc_r,
    DATA_TYPE * yc_i,
    DATA_TYPE POLYBENCH_1D(z_r, N-M, n-m),
    DATA_TYPE POLYBENCH_1D(z_i, N-M, n-m),
    DATA_TYPE POLYBENCH_1D(wa_r, N-M, n-m),
    DATA_TYPE POLYBENCH_1D(wa_i, N-M, n-m),
    DATA_TYPE * yp_r,
    DATA_TYPE * yp_i,
    DATA_TYPE POLYBENCH_1D(w_r, N, n),
    DATA_TYPE POLYBENCH_1D(w_i, N, n))
{
  unsigned int i = 0, j = 0;

  DATA_TYPE waz_r = 0;
  DATA_TYPE waz_i = 0;

  DATA_TYPE Bwa_r = 0;
  DATA_TYPE Bwa_i = 0;

#pragma scop
  //2. y_c(k) = w_q^{H} * x(k)
  for( i = 0 ; i<n ; i++)
  {
    *yc_r += (wq_r[i] * x_r[i]) - ((-wq_i[i]) * x_i[i]);
    *yc_i += ((-wq_i[i]) * x_r[i]) + (wq_r[i] * x_i[i]);
    
  }

  //3. z(k) = B^{H} * x(k)
  for (j = 0 ; j<(n-m) ; j++)
  {
    z_r[j] = 0.0;
    z_i[j] = 0.0;

    for(i = 0 ; i<n ; i++)
    {
      z_r[j] += (B_r[i][j] * x_r[i]) - ((-B_i[i][j]) * x_i[i]);
      z_i[j] += ((-B_i[i][j]) * x_r[i]) + (B_r[i][j] * x_i[i]);
    }
  }

  //4. y_p(k) = y_c(k) - w_a^H(k-1) * z(k)
  for (i = 0 ; i < (n-m) ; i++)
  {
    waz_r += (wa_r[i] * z_r[i]) - ((-wa_i[i]) * z_i[i]);
    waz_i += (-wa_i[i])*z_r[i] + wa_r[i]*z_i[i];
  }
  yp_r[0] = yc_r[0] - waz_r;
  yp_i[0] = yc_i[0] - waz_i;
  
  //5. w_a(k) = w_a (k-1) + ALPHA * z(k) * y_p^{*}(k)
  for (i = 0 ; i < (n-m) ; i++)
  {
    wa_r[i] = wa_r[i] + ALPHA * ((z_r[i] * yp_r[0]) - (z_i[i] * (-yp_i[0])));
    wa_i[i] = wa_i[i] + ALPHA * ((z_i[i] * (-yp_i[0])) - (z_r[i] * yp_r[0]));
  }

  //6. w(k) = w_q - B * w_a(k)
  for (i = 0 ; i < n ; i++)
  {
    Bwa_r = 0.0;
    Bwa_i = 0.0;
    for (j = 0 ; j < (n-m) ; j++)
    {
      Bwa_r += B_r[i][j] * wa_r[j] - B_i[i][j] * wa_i[j];
      Bwa_i += B_r[i][j] * wa_i[j] - B_i[i][j] * wa_r[j];
    }
    w_r[i] = wq_r[i] - Bwa_r;
    w_i[i] = wq_i[i] - Bwa_i;
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
  POLYBENCH_1D_ARRAY_DECL(x_r, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x_i, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(wq_r, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(wq_i, DATA_TYPE, N, n);
  POLYBENCH_2D_ARRAY_DECL(B_r, DATA_TYPE, N, N-m, n, n-m);
  POLYBENCH_2D_ARRAY_DECL(B_i, DATA_TYPE, N, N-m, n, n-m);
  DATA_TYPE yc_r;
  DATA_TYPE yc_i;
  POLYBENCH_1D_ARRAY_DECL(z_r, DATA_TYPE, N-M, n-m);
  POLYBENCH_1D_ARRAY_DECL(z_i, DATA_TYPE, N-M, n-m);
  POLYBENCH_1D_ARRAY_DECL(wa_r, DATA_TYPE, N-M, n-m);
  POLYBENCH_1D_ARRAY_DECL(wa_i, DATA_TYPE, N-M, n-m);
  DATA_TYPE yp_r;
  DATA_TYPE yp_i;
  POLYBENCH_1D_ARRAY_DECL(w_r, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(w_i, DATA_TYPE, N, n);

  /* Initialize array(s). */
  init_array (n,
      m,
      POLYBENCH_ARRAY(x_r),
      POLYBENCH_ARRAY(x_i),
      POLYBENCH_ARRAY(wq_r),
      POLYBENCH_ARRAY(wq_i),
      POLYBENCH_ARRAY(B_r),
      POLYBENCH_ARRAY(B_i),
      &yc_r,
      &yc_i,
      POLYBENCH_ARRAY(z_r),
      POLYBENCH_ARRAY(z_i),
      POLYBENCH_ARRAY(wa_r),
      POLYBENCH_ARRAY(wa_i),
      &yp_r,
      &yp_i,
      POLYBENCH_ARRAY(w_r),
      POLYBENCH_ARRAY(w_i));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  for (i = 0 ; i < NB_RUN ; i++)
  {
    //Just comment this function call if you don't want to 
    //update/modify the input data (x vector) of the radio.
    update_array(n,
        POLYBENCH_ARRAY(x_r),
        POLYBENCH_ARRAY(x_i));

    kernel_LMS_GSC (n, m,
        POLYBENCH_ARRAY(x_r),
        POLYBENCH_ARRAY(x_i),
        POLYBENCH_ARRAY(wq_r),
        POLYBENCH_ARRAY(wq_i),
        POLYBENCH_ARRAY(B_r),
        POLYBENCH_ARRAY(B_i),
        &yc_r,
        &yc_i,
        POLYBENCH_ARRAY(z_r),
        POLYBENCH_ARRAY(z_i),
        POLYBENCH_ARRAY(wa_r),
        POLYBENCH_ARRAY(wa_i),
        &yp_r,
        &yp_i,
        POLYBENCH_ARRAY(w_r),
        POLYBENCH_ARRAY(w_i));
  }

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w_r), POLYBENCH_ARRAY(w_i)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(x_r);
  POLYBENCH_FREE_ARRAY(x_i);
  POLYBENCH_FREE_ARRAY(wq_r);
  POLYBENCH_FREE_ARRAY(wq_i);
  POLYBENCH_FREE_ARRAY(B_r);
  POLYBENCH_FREE_ARRAY(B_i);
  POLYBENCH_FREE_ARRAY(z_r);
  POLYBENCH_FREE_ARRAY(z_i);
  POLYBENCH_FREE_ARRAY(wa_r);
  POLYBENCH_FREE_ARRAY(wa_i);
  POLYBENCH_FREE_ARRAY(w_r);
  POLYBENCH_FREE_ARRAY(w_i);

  return 0;
}
