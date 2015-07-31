/* M_inversion.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "matrix_inversion.h"


/* Array initialization. */
static
void init_array (unsigned int mat_size,
		 complex_number_t POLYBENCH_2D(M, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
		 complex_number_t POLYBENCH_2D(M_inv, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
		 complex_number_t POLYBENCH_2D(Q, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
		 complex_number_t POLYBENCH_2D(R, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
		 complex_number_t POLYBENCH_2D(R_inv, MAT_SIZE, MAT_SIZE, mat_size, mat_size))
{
  unsigned int l = 0, c = 0;
  unsigned int inc = 0;

  for (l = 0; l < MAT_SIZE ; l++)
    for (c = 0; c < MAT_SIZE ; c++)
    {
      /*M[l][c].r = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;*/
      /*M[l][c].i = ((DATA_TYPE)(rand()+1)/(DATA_TYPE)(RAND_MAX)) * MAX_VALUE;*/
      M[l][c].r = inc++;

      M_inv[c][l].r = 0.0;
      M_inv[c][l].i = 0.0;

      Q[l][c].r = 0.0;
      Q[l][c].i = 0.0;

      R[l][c].r = 0.0;
      R[l][c].i = 0.0;

      R_inv[c][l].r = 0.0;
      R_inv[c][l].i = 0.0;
    }
  M[MAT_SIZE-1][MAT_SIZE-1].r = 0.0;
}


/* DCE code. Must scan the entire live-out M.
   Can be used also to check the correctness of the output. */
static
void print_array(unsigned int mat_size,
    complex_number_t POLYBENCH_2D(M_inv, MAT_SIZE, MAT_SIZE, mat_size, mat_size))

{
  unsigned int l = 0, c = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("M_inv");
  for (l = 0; l < MAT_SIZE; l++)
  {
    fprintf (POLYBENCH_DUMP_TARGET, "\n");
    for (c = 0; c < MAT_SIZE; c++)
    {
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, M_inv[l][c].r);
      /*fprintf (POLYBENCH_DUMP_TARGET, "%s", (M_inv[l][c].i > 0.0) ? "+" : "");*/
      /*fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, M_inv[l][c].i);*/
      fprintf (POLYBENCH_DUMP_TARGET, " ");
    }
  }
  POLYBENCH_DUMP_END("M_inv");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_M_inversion(unsigned int mat_size,
    complex_number_t POLYBENCH_2D(M, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
    complex_number_t POLYBENCH_2D(M_inv, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
		complex_number_t POLYBENCH_2D(Q, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
		complex_number_t POLYBENCH_2D(R, MAT_SIZE, MAT_SIZE, mat_size, mat_size),
		complex_number_t POLYBENCH_2D(R_inv, MAT_SIZE, MAT_SIZE, mat_size, mat_size))
{
  int i = 0, j = 0, k = 0, l = 0;
  DATA_TYPE r = 0.0;
  complex_number_t c = {0.0, 0.0}, s = {0.0, 0.0};
  complex_number_t Rik = {0.0, 0.0}, Rjk = {0.0, 0.0};
  complex_number_t Qik = {0.0, 0.0}, Qjk = {0.0, 0.0};
  complex_number_t Qki = {0.0, 0.0}, Qkj = {0.0, 0.0};
  complex_number_t tmp = {0.0, 0.0};

#pragma scop
  //R = M, Q = Id
  for (i = 0 ; i < MAT_SIZE ; i++)
  {
    for (j = 0 ; j < MAT_SIZE ; j++)
    {
      R[i][j] = M[i][j];
      Q[i][j].r = 0.0;
      Q[i][j].i = 0.0;
    }
    Q[i][i].r = 1.0;
  }

  for (j = 0; j < MAT_SIZE-1; j++)
  {
    for (i = j+1; i < MAT_SIZE; i++)
    {
      r = copysign(1.0, R[i][j].r ) * SQRT_FUN(
          (R[j][j].r*R[j][j].r + R[j][j].i*R[j][j].i) +
          (R[i][j].r*R[i][j].r + R[i][j].i*R[i][j].i));

      c.r = R[j][j].r / r;
      c.i = R[j][j].i / r;
      s.r = -R[i][j].r / r;
      s.i = -R[i][j].i / r;

      // Apply the givens rotation:
      for (k = j; k < MAT_SIZE ; k++) {
          Rjk = R[j][k];
          Rik = R[i][k];

          R[j][k].r = (c.r * Rjk.r - c.i * Rjk.i) + ((-s.r) * Rik.r - (-s.i) * Rik.i);
          R[j][k].i = (c.i * Rjk.r + c.r * Rjk.i) + ((-s.i) * Rik.r + (-s.r) * Rik.i);

          R[i][k].r = (c.r * Rjk.r - c.i * Rjk.i) + (s.r * Rik.r - s.i * Rik.i);
          R[i][k].i = (c.i * Rjk.r + c.r * Rjk.i) + (s.i * Rik.r + s.r * Rik.i);
      }

      //G in Q by applying G^t to Q
      for (k = 0; k < MAT_SIZE ; k++) {
          Qki = Q[k][i];
          Qkj = Q[k][j];
          /*Q[j][k].r = Qkj * s + Qki * c;*/
          Q[j][k].r = (c.r * Qkj.r - c.i * Qkj.i) + ((-s.r) * Qki.r - (-s.i) * Qki.i);
          Q[j][k].i = (c.i * Qkj.r + c.r * Qkj.i) + ((-s.i) * Qki.r + (-s.r) * Qki.i);
          
          /*Q[k][i].r = Qkj * s + Qki * c;*/
          Q[k][i].r = (c.r * Qki.r - c.i * Qki.i) + (s.r * Qkj.r - s.i * Qkj.i);
          Q[k][i].i = (c.i * Qki.r + c.r * Qki.i) + (s.i * Qkj.r + s.r * Qkj.i);
      }
    }
  }

  //Inverse R into R_inv
  for(i=MAT_SIZE-1 ; i >= 0  ; i--)
  {
      R_inv[i][i].r = R[i][i].r / (R[i][i].r * R[i][i].r + R[i][i].i * R[i][i].i);
      R_inv[i][i].i = -R[i][i].i / (R[i][i].r * R[i][i].r + R[i][i].i * R[i][i].i);
      for(j=i+1; j<MAT_SIZE; j++)
      {
          for(k=i+1; k<MAT_SIZE; k++)
          {
              R_inv[i][j].r += R[i][k].r * R_inv[k][j].r ;
              R_inv[i][j].i += R[i][k].i * R_inv[k][j].i ;
          }
          tmp = R_inv[i][j];
          R_inv[i][j].r = (tmp.r * (-R_inv[i][i].r)) - (tmp.i * (-R_inv[i][i].i));
          R_inv[i][j].i = (tmp.r * (-R_inv[i][i].i)) + (tmp.i * (-R_inv[i][i].r));
      }
  }

  //M^{-1} = R^{-1} * Q^H
  for (i=0; i<MAT_SIZE ;i++)
  {
    for (j=0; j<MAT_SIZE ;j++)
    {
      M_inv[i][j].r = 0.0;
      M_inv[i][j].i = 0.0;
      for (k=0; k<MAT_SIZE ;k++)
      {
        M_inv[i][j].r += R_inv[i][k].r * Q[j][k].r - (R_inv[i][k].i * (-Q[j][k].i));
        M_inv[i][j].i += R_inv[i][k].i * Q[j][k].r + R_inv[i][k].r * Q[j][k].i;
      }
    }
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  unsigned int mat_size = MAT_SIZE;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(M, complex_number_t, MAT_SIZE, MAT_SIZE, mat_size, mat_size);
  POLYBENCH_2D_ARRAY_DECL(M_inv, complex_number_t, MAT_SIZE, MAT_SIZE, mat_size, mat_size);
  POLYBENCH_2D_ARRAY_DECL(Q, complex_number_t, MAT_SIZE, MAT_SIZE, mat_size, mat_size);
  POLYBENCH_2D_ARRAY_DECL(R, complex_number_t, MAT_SIZE, MAT_SIZE, mat_size, mat_size);
  POLYBENCH_2D_ARRAY_DECL(R_inv, complex_number_t, MAT_SIZE, MAT_SIZE, mat_size, mat_size);


  /* Initialize array(s). */
  init_array (mat_size,
      POLYBENCH_ARRAY(M),
      POLYBENCH_ARRAY(M_inv),
      POLYBENCH_ARRAY(Q),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(R_inv));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_M_inversion (mat_size,
      POLYBENCH_ARRAY(M),
      POLYBENCH_ARRAY(M_inv),
      POLYBENCH_ARRAY(Q),
      POLYBENCH_ARRAY(R),
      POLYBENCH_ARRAY(R_inv));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out M must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(mat_size, POLYBENCH_ARRAY(M_inv)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(M);
  POLYBENCH_FREE_ARRAY(M_inv);

  return 0;
}
