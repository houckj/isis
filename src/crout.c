/* Crout LU factorization */

/* This file is a subset of routines from the JDMath library by
 * John E. Davis, modified by John C. Houck (1/28/99)
 * to work apart from JDMath
 *
 * Copyright (c) 1994, 1998 John E. Davis 
 *
 * You may distribute this file under the terms the GNU General Public
 * License.  See the file COPYING for more information.
 */

/* $Id: crout.c,v 1.2 2003/02/10 00:46:32 houck Exp $ */

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "isis.h"
#include "isismath.h"

#define TOLERANCE 1.e-23

/*
 * Crout's algorithm with partial pivoting.
 * 
 * This algorithm performs a Lower-Upper factorization of a matrix A.  The
 * method is actually quite simple.  Write:
 * 
 *     L_ij = L_ij Q(i>=j)
 *     U_ij = U_ij Q(i<=j)
 * 
 * where Q(x) is 1,0 if x is TRUE,FALSE.  Then, 
 * 
 *     A_ij = \sum_k L_ik U_kj Q(i>=k) Q(k<=j)
 * or
 * 
 *     A_ij = \sum_{k<=min(i,j)} L_ik U_kj
 *
 * or
 * 
 *     A_ij = L_ii U_ij + \sum_{k<i} L_ik U_kj     (i <= j)
 *          = L_ij U_jj + \sum_{k<j} L_ik U_kj     (j <= i)
 * 
 * We are free to choose L_ii = 1 for all i (since there are n^2+n unknown
 * matrix elements.  Then for i <= j, we obtain:
 * 
 *     U_ij = A_ij - \sum{k<i} L_ik U_kj    (i <= j)
 * 
 * and for j > i:
 * 
 *     L_ij = (1/U_jj) (A_ij - \sum{k<j} L_ik U_kj)  (j >= i)
 * 
 * These equations allow a straightforward recursive solution for U_ij and L_ij
 * since the sums on the RHS contain L_ik and U_kj that are known.  In fact,
 * by not storing the diagonal elements of L, which have been set to 1, we can
 * do the factorization in place since each A_ij in the original matrix is
 * explicitly one time.
 * 
 * The only problem is that because of the presence of (1/U_jj), this algorithm
 * is unstable without pivoting.  Pivoting can be easily accomplished by noting
 * that the computation of U_ij involves the element A_ij, the values of
 * L_ij from the same row, and the values of U_kj from earlier rows.  The key 
 * point is that L_ij involves values of L from the same row.  Thus, we may 
 * perform simple row swapping.
 * 
 * Note that the algorithm equilibrates the matrix in order to find the 
 * appropriate pivot.
 */

int JDM_lu_decomp (double **a, unsigned int n, unsigned int *pivot, double epsilon, int *parityp)
{
   unsigned int j;
   int det = 1;
   double *eq;

   if ((n == 0) || (a == NULL) || (pivot == NULL))
     return -1;

   if (NULL == (eq = _JDM_equilibrate (a, n, NULL)))
     return -1;

   for (j = 0; j < n; j++)
     {
	unsigned int i, ipiv;
	double big;

	/* Do i < j case */
	for (i = 0; i < j; i++)
	  {
	     double *ai = a[i];
	     ai[j] -= _JDM_innerprod_col (ai, a, j, i);
	  }
	
	/* Now i >= j.
	 * We also pivot in this step.  Since we do not know what value
	 * goes on the diagonal, compute all possible diagonal values, then
	 * pivot, then divide.
	 */
	
	big = 0;
	ipiv = j;

	for (i = j; i < n; i++)
	  {
	     double *ai = a[i];
	     ai[j] -= _JDM_innerprod_col (ai, a, j, j);
	     
	     if (eq[j] * fabs (ai[j]) > big)
	       {
		  big = eq[j] * fabs (ai[j]);
		  ipiv = i;
	       }
	  }
	if (ipiv != j)
	  {
	     double eq_j;
	     _JDMswap_dvector (a[ipiv], a[j], n);
	     eq_j = eq[j]; eq[j] = eq[ipiv]; eq[ipiv] = eq_j;
	     det = -det;
	  }
	pivot[j] = ipiv;

	big = a[j][j];
	if (fabs (big) < epsilon)
	  {
	     ISIS_FREE (eq);
	     eq = NULL;
	     return -1;
	  }


	big = 1.0/big;
	for (i = j + 1; i < n; i++)
	  a[i][j] *= big;
     }

   if (parityp != NULL)
     *parityp = det;

   ISIS_FREE (eq);
   return 0;
}

static void backsubst (double **a, unsigned int n, double *b)
{
   unsigned int i;

   i = n;
   do
     {
	double *ai;
	unsigned int i1;
	
	i1 = i;
	i--;
	ai = a[i];
	b[i] = (b[i] - _JDM_innerprod (ai + i1, b + i1, n - i1))/ai[i];
     }
   while (i != 0);
}

static void forwsubst (double **a, unsigned int n, double *b, unsigned int *pivot)
{
   unsigned int i;
   /* Do forward substitution and unscramble b in the process */
   for (i = 0; i < n; i++)
     {
	double x_i;
	unsigned int j;
	
	j = pivot[i];
	x_i = b[j] - _JDM_innerprod (a[i], b, i);
	if (i != j)
	  b[j] = b[i];
	b[i] = x_i;
     }
}


int JDM_lu_backsubst (double **a, unsigned int n, unsigned int *pivot, double *b)
{
   if ((a == NULL) || (b == NULL) || (pivot == NULL) || (n == 0))
     return -1;

   forwsubst (a, n, b, pivot);
   backsubst (a, n, b);

   return 0;
}

static void forwsubst_unit_vector (double **a, unsigned int n,
				   double *b, unsigned int *pivot,
				   unsigned int col)
{
   unsigned int i;
   
   for (i = 0; i < n; i++)
     {
	if (i > col)
	  {
	     double *bcol = b + col;
	     while (i < n)
	       {
		  b[i] = -_JDM_innerprod (a[i] + col, bcol, i - col);
		  i++;
	       }
	     break;
	  }
	
	if (col == i)
	  col = pivot[i];
	else if (col == pivot[i])
	  col = i;
	
	b[i] = (i == col) ? 1.0 : 0.0;
     }
}

int JDM_ludecomp_inverse (double **a, unsigned int n)
{
   unsigned int *pivot;
   double *col;
   double **inv_a;
   unsigned int i, j;
   int ret;
   
   if (n == 0)
     return 0;
   
   if (a == NULL)
     return -1;
   
   pivot = NULL;
   col = NULL;
   inv_a = NULL;
   ret = -1;
   
   if ((NULL == (pivot = (unsigned int *) ISIS_MALLOC (n * sizeof(*pivot))))
       || (NULL == (col = JDMdouble_vector (n)))
       || (NULL == (inv_a = JDMdouble_matrix (n, n))))
     goto return_error;
   
   if (-1 == JDM_lu_decomp (a, n, pivot, TOLERANCE, NULL))
     goto return_error;
   
   for (j = 0; j < n; j++)
     {
	forwsubst_unit_vector (a, n, col, pivot, j);
	backsubst (a, n, col);
	
	for (i = 0; i < n; i++)
	  inv_a [i][j] = col[i];
     }
   
   for (i = 0; i < n; i++)
     {
	double *ai, *inv_ai;
	
	ai = a[i];
	inv_ai = inv_a[i];
	for (j = 0; j < n; j++)
	  ai[j] = inv_ai[j];
     }
   
   ret = 0;
   /* drop */
   
   return_error:
   ISIS_FREE (pivot);
   ISIS_FREE (col);
   JDMfree_double_matrix (inv_a, n);
   
   return ret;
}


#ifdef TEST_LUDECOMP
#include <stdio.h>

#define N 100

int main (int argc, char **argv)
{
   unsigned int n = N;
   unsigned int i, j;
   double **a, **a1;
   unsigned int piv[N];
   double b[N], b1[N];
   unsigned int count;
   
   count = 5;
   
   a = JDMdouble_matrix (n, n);
   a1 = JDMdouble_matrix (n, n);
   
   while (count)
     {
	double max;

	count--;

	for (i = 0; i < n; i++)
	  {
	     b[i] = b1[i] = (double) rand();
	     for (j = 0; j < n; j++)
	       a[i][j] = a1[i][j] = (double) rand();
	  }

#ifdef SOLVE 
	if (-1 == JDM_lu_decomp (a, n, piv, TOLERANCE, NULL)
	    || -1 == JDM_lu_backsubst (a, n, piv, b))
	  continue;
	
        for (i = 0; i < n; i++)
	  {
	     double sum = 0;
	     for (j = 0; j < n; j++)
	       sum += a1[i][j] * b[j];
	     
	     fprintf (stdout, "%e\n", sum - b1[i]);
	  }
#else
	if (-1 == JDM_ludecomp_inverse (a, n))
	  continue;
	
	max = 0;

	for (i = 0; i < n; i++)
	  {
	     unsigned int k;
	     for (k = 0; k < n; k++)
	       {
		  double sum = 0;
		  for (j = 0; j < n; j++)
		    sum += a[i][j] * a1[j][k];
		  
		  if (i == k)
		    sum -= 1.0;

		  if (fabs (sum) > max) max = fabs(sum);
	       }
	  }
	fprintf (stdout, "%e\n", max);
#endif	
     }
   
   JDMfree_double_matrix (a, n);
   JDMfree_double_matrix (a1,n);

   return 0;
}

#endif
	
