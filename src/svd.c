/* -*- mode: C; mode: fold -*- */
/* Test with: gcc -DTEST_SVD -O2 $CFLAGS svd.c -o svd_test -lm */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "isismath.h"
#include "svd.h"

/* This implementation of SVD, from NODElib by Gary William Flake,
   is a modification of the SVD routine from Nash (1990).

   Copyright (c) 1992-2000 by Gary William Flake. 

   NODElib is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2 of the License, or (at your option) any later version. 

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 
*/

/* Original comments from Bryant Marks:
 
   This SVD routine is based on pgs 30-48 of "Compact Numerical Methods
   for Computers" by J.C. Nash (1990), used to compute the pseudoinverse.
   Modifications include:
        Translation from Pascal to ANSI C.
        Array indexing from 0 rather than 1.
        Float replaced by double everywhere.
        Support for the Matrix structure.
        I changed the array indexing so that the matricies (float [][])
           could be replaced be a single list (double *) for more
           efficient communication with Mathematica.

  From: bryant@sioux.stanford.edu (Bryant Marks)
  >
  > A couple of things to note: A needs to have twice as much room
  > allocated for it (2*n + 2*m) since the W in the svd function requires
  > this (part of a rotation algorithm).  After the routine has run W
  > contains two maticies of the decomposition A = USV'.  The first nRow
  > rows contain the product US and the next nCol rows contain V (not V').
  > Z is equal to the vector of the sqares of the diagonal elements of S. */

/* Comments from GWF: The note above is not strictly correct.  To compute
   the SVD of an (m x n) matrix, W must be ((m + n) x n) in size.

   12-26-96: I slightly rewrote things to include V so that the interface
   could at least be close to a lapack-like routine. 

   01-02-97: Rewrote it so that U,S,V was actually returned.  V is now optional. */

#define TOLERANCE 1.0e-22

/* Why does glibc use y0 for a Bessel function???? */
#define y0 y0_silly

void svd (double *U, double *S, double *V, int nRow, int nCol) /*{{{*/
{
  int i, j, k, EstColRank, RotCount, SweepCount, slimit;
  double eps, e2, tol, vt, p, x0, y0, q, r, c0, s0, d1, d2;

  eps = TOLERANCE;
  slimit = nCol / 4;
  if (slimit < 6.0)
    slimit = 6;
  SweepCount = 0;
  e2 = 10.0 * nRow * eps * eps;
  tol = eps * .1;
  EstColRank = nCol;
  if(V)
    for (i = 0; i < nCol; i++)
      for (j = 0; j < nCol; j++) {
	V[nCol * i + j] = 0.0;
	V[nCol * i + i] = 1.0;
      }
  RotCount = EstColRank * (EstColRank - 1) / 2;
  while (RotCount != 0 && SweepCount <= slimit) {
    RotCount = EstColRank * (EstColRank - 1) / 2;
    SweepCount++;
    for (j = 0; j < EstColRank - 1; j++) {
      for (k = j + 1; k < EstColRank; k++) {
	p = q = r = 0.0;
	for (i = 0; i < nRow; i++) {
	  x0 = U[nCol * i + j];
	  y0 = U[nCol * i + k];
	  p += x0 * y0;
	  q += x0 * x0;
	  r += y0 * y0;
	}
	S[j] = q;
	S[k] = r;
	if (q >= r) {
	  if (q <= e2 * S[0] || fabs(p) <= tol * q)
	    RotCount--;
	  else {
	    p /= q;
	    r = 1 - r / q;
	    vt = sqrt(4 * p * p + r * r);
	    c0 = sqrt(fabs(.5 * (1 + r / vt)));
	    s0 = p / (vt * c0);
	    for (i = 0; i < nRow; i++) {
	      d1 = U[nCol * i + j];
	      d2 = U[nCol * i + k];
	      U[nCol * i + j] = d1 * c0 + d2 * s0;
	      U[nCol * i + k] = -d1 * s0 + d2 * c0;
	    }
	    if(V)
	      for (i = 0; i < nCol; i++) {
		d1 = V[nCol * i + j];
		d2 = V[nCol * i + k];
		V[nCol * i + j] = d1 * c0 + d2 * s0;
		V[nCol * i + k] = -d1 * s0 + d2 * c0;
	      }
	  }
	}
	else {
	  p /= r;
	  q = q / r - 1;
	  vt = sqrt(4 * p * p + q * q);
	  s0 = sqrt(fabs(.5 * (1 - q / vt)));
	  if (p < 0)
	    s0 = -s0;
	  c0 = p / (vt * s0);
	  for (i = 0; i < nRow; i++) {
	    d1 = U[nCol * i + j];
	    d2 = U[nCol * i + k];
	    U[nCol * i + j] = d1 * c0 + d2 * s0;
	    U[nCol * i + k] = -d1 * s0 + d2 * c0;
	  }
	  if(V)
	    for (i = 0; i < nCol; i++) {
	      d1 = V[nCol * i + j];
	      d2 = V[nCol * i + k];
	      V[nCol * i + j] = d1 * c0 + d2 * s0;
	      V[nCol * i + k] = -d1 * s0 + d2 * c0;
	    }
	}
      }
    }
    while (EstColRank >= 3 && S[(EstColRank - 1)] <= S[0] * tol + tol * tol)
      EstColRank--;
  }
  for(i = 0; i < nCol; i++)
    S[i] = sqrt(S[i]);
  for(i = 0; i < nCol; i++)
    for(j = 0; j < nRow; j++)
      U[nCol * j + i] = U[nCol * j + i] / S[i];
}

/*}}}*/

static int zero_small_singular_values (double *s, int n) /*{{{*/
{
   double cut, max;
   int j;
   
   /* Find the largest singular value */
   max = fabs(s[0]);
   for (j = 1; j < n; j++)
     {
	double a = fabs(s[j]);
	if (a > max) max = a;
     }
   
   /* LAPACK dgelss uses cut = rcond * max(s)
    * where rcond = 1 / (norm(A) * norm(inv(A))) 
    *               where A = original matrix
    *               or DBL_EPSILON if rcond < 0 on input
    *       max(s) = largest singular value */
   cut = 1.e4 * DBL_EPSILON * max;
   
   /* Now zero out the "small" singular values */
   for (j = 0; j < n; j++)
     {
	if (fabs(s[j]) < cut) 
	  s[j] = 0.0;
     }
 
   return 0;
}

/*}}}*/

static int svd_bksb (int m, int n, double *b, double *u, double *s, double *v, double *x) /*{{{*/
{
   double *t;
   double sum;
   int j;
   
   t = ISIS_MALLOC (n * sizeof(*t));
   if (t == NULL)
     return -1;
   
   for (j = 0; j < n; j++)
     {
	int i;
	sum = 0.0;
	if (fabs(s[j]) > 0.0)
	  {
	     for (i = 0; i < m; i++)
	       sum += u [i*m + j] * b[i];
	     sum /= s[j];
	  }
	t[j] = sum;
     }
   
   for (j = 0; j < n; j++)
     {
	int k;
	sum = 0.0;
	for (k = 0; k < n; k++)
	  sum += v[j*n + k] * t[k];
	x[j] = sum;
     }
   
   ISIS_FREE (t);
   
   return 0;
}

/*}}}*/

int svd_solve (double *a, unsigned int m, unsigned int n, double *b, double *x) /*{{{*/
{
   double *u=NULL, *v=NULL, *s=NULL;

   if ((a==NULL) || (b == NULL) || (x == NULL) || (m > n))
     return -1;

   u = ISIS_MALLOC (m * n * sizeof(*u));
   v = ISIS_MALLOC (n * n * sizeof(*v));
   s = ISIS_MALLOC (    n * sizeof(*s));
   if (!(u && v && s))
     {
        ISIS_FREE(u); ISIS_FREE(v); ISIS_FREE(s);
        return -1;
     }

   memcpy ((char *)u, (char *)a, m * n * sizeof(double));

   svd (u, s, v, m, n);
   zero_small_singular_values (s, n);
   svd_bksb (m, n, b, u, s, v, x);

   ISIS_FREE (u); ISIS_FREE (v); ISIS_FREE (s);
   return 0;
}

/*}}}*/

#ifdef TEST_SVD

static double Tol = 500.0 * DBL_EPSILON;

static int init_random_matrix (double *a, int m, int n) /*{{{*/
{
   int i, j;

   /* m rows, n columns  a[i][j] => a[i*n + j] is row i, column j
    * because the last index changes fastest in C */
   for (i = 0; i < m; i++)
     {
	for (j = 0; j < n; j++)
	  {
	     /* a(i,j) row i, column j*/
	     a[i * n + j] = (10.0 * rand ()) / (RAND_MAX + 1.0);
	  }
     }
   
   return 0;
}

/*}}}*/

static int init_random_vector (double *b, int n) /*{{{*/
{
   int j;

   for (j = 0; j < n; j++)
     {
	b[j] = (10.0 * rand ()) / (RAND_MAX + 1.0);
     }
   
   return 0;
}

/*}}}*/

static int verify_svd (double *a, int m, int n, double *u, double *s, double *v) /*{{{*/
{
   int i, j, k;
   double tmp, max;
   
   tmp = max = 0.0;
      
   for (i = 0; i < n; i++)
     {
	for (j = 0; j < m; j++)
	  {
	     tmp = 0.0;
	     for (k = 0; k < n; k++)
	       {
		  /* U*s*transpose(V) = u(j,k) * v(i,k) * s(k) */
		  tmp += u[j*n + k] * v[i*n + k] * s[k];
	       }
	     /* compute infinity norm = max|U*s*tr(V) - a(j,i)| */
	     tmp = fabs (tmp - a[j*n + i]);
	     if (tmp > max) max = tmp;
	  }
     }

   if (fabs(max) > Tol)
     fprintf (stdout, "failed:  *** Size mxn = %dx%d;   Norm of (A - U*s*transpose(V)) = %15.8e\n", m, n, max);
   else
     fputs ("ok:  SVD\n", stdout);
   
   return 0;
}

/*}}}*/

static int verify_solution (double *a, int m, int n, double *b, double *x) /*{{{*/
{
   int i, j, flag;
   
   flag = 0;
   for (j = 0; j < m; j++)
     {
	double sum = 0.0;
	double diff;
	for (i = 0; i < n; i++)
	  {
	     /* sum_i a(j,i) * x[i] */
	     sum += a[j*n + i] * x[i];
	  }
	diff = sum - b[j];
	if (fabs(diff) > Tol)
	  {
	     fprintf (stdout, "*** [%3d]:  A * x - b =  %15.8e\n", j, sum - b[j]);
	     flag = 1;
	  }
     }
   
   if (flag == 0)
     fputs ("ok:  SVD solution\n", stdout);
   else 
     fputs ("failed:  SVD solution\n", stdout);

   return 0;
}

/*}}}*/

int main (void) /*{{{*/
{
   double *a, *b, *u, *s, *v, *x;
   int m, n;
   
   m = n = 16;
   
   a = ISIS_MALLOC (m * n * sizeof(*a));
   u = ISIS_MALLOC (m * n * sizeof(*u));
   v = ISIS_MALLOC (n * n * sizeof(*v));
   s = ISIS_MALLOC (    n * sizeof(*s));
   b = ISIS_MALLOC (    m * sizeof(*b));   
   x = ISIS_MALLOC (    n * sizeof(*b));
   if (!(a && u && v && s && b && x))
     goto finish;
   
   fprintf (stdout, "Tolerance = %15.8e\n", Tol);
   
   (void) init_random_matrix (a, m, n);
   memcpy ((char *)u, (char *)a, m * n * sizeof(double));   
   svd (u, s, v, m, n);
   (void) verify_svd (a, m, n, u, s, v);
   
   (void) init_random_matrix (a, m, n);
   (void) init_random_vector (b, m);
   (void) svd_solve (a, m, n, b, x);
   (void) verify_solution (a, m, n, b, x);   
   
   finish:
   
   ISIS_FREE (a); 
   ISIS_FREE (b); 
   ISIS_FREE (u);    
   ISIS_FREE (v); 
   ISIS_FREE (s);
   
   return 0;
}

/*}}}*/

#endif
