#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/*#define MINI_DATASET*/
/* Include benchmark-specific header. */
#include "RLS.h"


/* Array update. */
static
void update_array(int n,
    complex_number_t POLYBENCH_1D(x, N, n))
{
  unsigned int i = 0;

  for(i = 0 ; i<n ; i++)
  {
    x[i].r = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;
    x[i].i = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;
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
    complex_number_t POLYBENCH_1D(zP, N-M, n-m),
    complex_number_t POLYBENCH_1D(Pz, N-M, n-m),
    complex_number_t POLYBENCH_1D(g, N-M, n-m),
    complex_number_t POLYBENCH_2D(gz, N-M, N-M, n-m, n-m),
    complex_number_t POLYBENCH_2D(P, N-M, N-M, n-m, n-m),
    complex_number_t POLYBENCH_1D(wa, N-M, n-m),
    complex_number_t * yp,
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0;

  *yc = (complex_number_t) {0.0, 0.0};

  for(i = 0 ; i<n ; i++)
  {
    wq[i].r = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;
    wq[i].i = ((float)(rand()+1)/(float)(RAND_MAX)) * MAX_VALUE;

    w[i].r = 0.0;
    w[i].i = 0.0;

    for (j = 0 ; j<(n-m) ; j++)
    {
      B[i][j].r = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
      B[i][j].i = ((float)(rand()+1)/(float)(RAND_MAX)) * 10.0;
    }
  }
  
  update_array(n, x);

  for (i = 0 ; i<(n-m) ; i++)
  {
    g[i].r = 0.0;
    g[i].i = 0.0;
  }

  
  //1. w_a(0) = 0, P(0) = I_{N-m}
  for (i = 0 ; i<(n-m) ; i++)
  {
    wa[i].r = 0.0;
    wa[i].i = 0.0;

    z[i].r = 0.0;
    z[i].i = 0.0;
    
    Pz[i].r = 0.0;
    Pz[i].i = 0.0;

    zP[i].r = 0.0;
    zP[i].i = 0.0;

    for (j = 0 ; j<(n-m) ; j++)
    {
      P[i][j].r = 0.0;
      P[i][j].i = 0.0;

      gz[i][j].r = 0.0;
      gz[i][j].i = 0.0;
    }
    P[i][i].r = 1.0;
    P[i][i].i = 1.0;
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
void kernel_LMS_GSC(int n, int m,
    complex_number_t POLYBENCH_1D(x, N, n),
    complex_number_t POLYBENCH_1D(wq, N, n),
    complex_number_t POLYBENCH_2D(B, N, N-m, n, n-m),
    complex_number_t * yc,
    complex_number_t POLYBENCH_1D(z, N-M, n-m),
    complex_number_t POLYBENCH_1D(zP, N-M, n-m),
    complex_number_t POLYBENCH_1D(Pz, N-M, n-m),
    complex_number_t POLYBENCH_1D(g, N-M, n-m),
    complex_number_t POLYBENCH_2D(gz, N-M, N-M, n-m, n-m),
    complex_number_t POLYBENCH_2D(P, N-M, N-M, n-m, n-m),
    complex_number_t POLYBENCH_1D(wa, N-M, n-m),
    complex_number_t * yp,
    complex_number_t POLYBENCH_1D(w, N, n))
{
  unsigned int i = 0, j = 0, k = 0;

  complex_number_t tmp = (complex_number_t) {0.0, 0.0};
  DATA_TYPE tmp_norm_pow2 = 0;

  complex_number_t gPz = (complex_number_t) {0.0, 0.0};

  complex_number_t waz = (complex_number_t) {0.0, 0.0};

  complex_number_t Bwa = (complex_number_t) {0.0, 0.0};

#pragma scop
  //2. y_c(k) = w_q^{H} * x(k)
  yc->r = 0.0;
  yc->i = 0.0;
  for( i = 0 ; i<n ; i++)
  {
    yc->r += (wq[i].r * x[i].r) - ((-wq[i].i) * x[i].i);
    yc->i += ((-wq[i].i) * x[i].r) + (wq[i].r * x[i].i);
    
  }

  //3. z(k) = B^{H} * x(k)
  for (i = 0 ; i<(n-m) ; i++)
  {
    z[i].r = 0.0;
    z[i].i = 0.0;
  }
  for(j = 0 ; j<n ; j++)
  {
    for (i = 0 ; i<(n-m) ; i++)
    {
      z[i].r += (B[j][i].r * x[j].r) - ((-B[j][i].i) * x[j].i);
      z[i].i += ((-B[j][i].i) * x[j].r) + (B[j][i].r * x[j].i);
    }
  }

  //4. y_p(k) = y_c(k) - w_a^H(k-1) * z(k)
  waz.r = 0.0;
  waz.i = 0.0;
  for (i = 0 ; i < (n-m) ; i++)
  {
    waz.r += (wa[i].r * z[i].r) - ((-wa[i].i) * z[i].i);
    waz.i += (-wa[i].i)*z[i].r + wa[i].r*z[i].i;
  }
  yp[0].r = yc[0].r - waz.r;
  yp[0].i = yc[0].i - waz.i;

  //5. g(k) = P(k-1) * z(k) / (μ + z^H(k) * P(k-1) * z(k))
  //where μ was set equal to 1
  for (i = 0 ; i < (n-m) ; i++)
  {
    Pz[i].r = 0.0;
    Pz[i].i = 0.0;
  }
  for (j = 0 ; j < (n-m) ; j++)
  {
    for (i = 0 ; i < (n-m) ; i++)
    {
      Pz[i].r += (P[j][i].r * z[j].r) - (P[j][i].i * z[j].i);
      Pz[i].i += (P[j][i].r * z[j].i) + (P[j][i].i * z[j].r);
    }
  }
  for (i = 0 ; i < (n-m) ; i++)
  {
    tmp.r += (z[i].r * Pz[i].r) - ((-z[i].i) * Pz[i].i);
    tmp.i += (-z[i].i)*Pz[i].r + z[i].r*Pz[i].i;
  }
  tmp.r += MU;
  tmp.i += MU;
  tmp_norm_pow2 = tmp.r * tmp.r + tmp.i * tmp.i;

  for (i = 0 ; i < (n-m) ; i++)
  {
    g[i].r = ((Pz[i].r * tmp.r) - (Pz[i].i * (-tmp.i))) / tmp_norm_pow2;
    g[i].i = ((Pz[i].r * (-tmp.i)) + (Pz[i].i * tmp.r)) / tmp_norm_pow2;
  }


  //6. P(k) = μ^{-1} * (P(k-1) - g(k) * z^H(k) * P(k-1))
  for (i = 0 ; i < (n-m) ; i++)
  {
    zP[i].r = 0.0;
    zP[i].i = 0.0;
  }
  for (i = 0 ; i < (n-m) ; i++)
  {
    for (j = 0 ; j < (n-m) ; j++)
    {
      zP[j].r += (z[j].r * P[i][j].r) - ((-z[j].i) * P[i][j].i);
      zP[j].i += ((-z[j].i) * P[i][j].r) + (z[j].r * P[i][j].i);
    }
  }

  for (i = 0 ; i < (n-m) ; i++)
  {
    for (j = 0 ; j < (n-m) ; j++)
    {
      gPz.r = (g[i].r * zP[j].r) - (g[i].i * zP[j].i);
      gPz.i = (g[i].r * zP[j].i) - (g[i].i * zP[j].r);
      P[j][i].r = (P[j][i].r - gPz.r) / MU;
      P[j][i].i = (P[j][i].i - gPz.i) / MU;
    }
  }

  
  //7. w_a(k) = w_a(k-1) + g(k) * y_p^{*}(k)
  for (i = 0 ; i < (n-m) ; i++)
  {
    wa[i].r = wa[i].r + ((g[i].r * yp[0].r) - (g[i].i * (-yp[0].i)));
    wa[i].i = wa[i].i + ((g[i].i * (-yp[0].i)) - (g[i].r * yp[0].r));
  }

  //8. w(k) = w_q - B * w_a(k)
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
  POLYBENCH_1D_ARRAY_DECL(zP, complex_number_t, N-M, n-m);
  POLYBENCH_1D_ARRAY_DECL(Pz, complex_number_t, N-M, n-m);
  POLYBENCH_1D_ARRAY_DECL(g, complex_number_t, N-M, n-m);
  POLYBENCH_2D_ARRAY_DECL(gz, complex_number_t, N, N-M, n-m, n-m);
  POLYBENCH_2D_ARRAY_DECL(P, complex_number_t, N-M, N-M, n-m, n-m);
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
      POLYBENCH_ARRAY(zP),
      POLYBENCH_ARRAY(Pz),
      POLYBENCH_ARRAY(g),
      POLYBENCH_ARRAY(gz),
      POLYBENCH_ARRAY(P),
      POLYBENCH_ARRAY(wa),
      &yp,
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

    kernel_LMS_GSC (n, m,
        POLYBENCH_ARRAY(x),
        POLYBENCH_ARRAY(wq),
        POLYBENCH_ARRAY(B),
        &yc,
        POLYBENCH_ARRAY(z),
        POLYBENCH_ARRAY(zP),
        POLYBENCH_ARRAY(Pz),
        POLYBENCH_ARRAY(g),
        POLYBENCH_ARRAY(gz),
        POLYBENCH_ARRAY(P),
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
  POLYBENCH_FREE_ARRAY(zP);
  POLYBENCH_FREE_ARRAY(Pz);
  POLYBENCH_FREE_ARRAY(g);
  POLYBENCH_FREE_ARRAY(gz);
  POLYBENCH_FREE_ARRAY(P);
  POLYBENCH_FREE_ARRAY(wa);
  POLYBENCH_FREE_ARRAY(w);

  return 0;
}
