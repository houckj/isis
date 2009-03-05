
/* This file is a subset of routines from the JDMath library by
 * John E. Davis, modified by John C. Houck (1/28/99)
 * to work apart from JDMath
 *
 * Copyright (c) 1994, 1998 John E. Davis 
 *
 * You may distribute this file under the terms the GNU General Public
 * License.  See the file COPYING for more information.
 */

/* $Id: dblas.c,v 1.1 2003/01/20 16:36:44 houck Exp $ */

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "isis.h"
#include "isismath.h"

void _JDMswap_dvector (double *a, double *b, unsigned int n)
{
   double tmp;
   unsigned int i;
   
   for (i = 0; i < n; i++)
     {
	tmp = a[i];
	a[i] = b[i];
	b[i] = tmp;
     }
}


double _JDM_innerprod (double *a, double *b, unsigned int n)
{
   unsigned int n8;
   double sum;

   n8 = n / 8;
   n = n % 8;
   sum = 0.0;
   while (n8--)
     {
	sum += a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
	  + a[4]*b[4] + a[5]*b[5] + a[6]*b[6] + a[7]*b[7];
	a += 8;
	b += 8;
     }

   while (n--)
     {
	sum += a[0] * b[0];
	a++;
	b++;
     }
   return sum;
}

double _JDM_innerprod_col (double *a, double **b, unsigned int col, unsigned int n)
{
   unsigned int n8;
   double sum;

   n8 = n / 8;
   n = n % 8;
   sum = 0.0;
   while (n8--)
     {
	sum += a[0]*b[0][col] + a[1]*b[1][col] + a[2]*b[2][col] 
	  + a[3]*b[3][col] + a[4]*b[4][col] + a[5]*b[5][col] 
	  + a[6]*b[6][col] + a[7]*b[7][col];
	a += 8;
	b += 8;
     }

   while (n--)
     {
	sum += a[0] * b[0][col];
	a++;
	b++;
     }

   return sum;
}


static double _JDM_nvector_max_abs (double *a, unsigned int n)
{
   unsigned int i;
   double max;

   if (n == 0)
     return -1.0;
   
   max = fabs (a[0]);
   for (i = 1; i < n; i++)
     {
	if (fabs(a[i]) > max)
	  max = fabs (a[i]);
     }
   return max;
}

double *_JDM_equilibrate (double **a, unsigned int n, double *eq)
{
   unsigned int i;
   int is_malloced;
   
   is_malloced = 0;
   if (eq == NULL)
     {
	eq = JDMdouble_vector (n+1);
	if (eq == NULL)
	  return eq;
	is_malloced = 1;
     }

   for (i = 0; i < n; i++)
     {
	double max = _JDM_nvector_max_abs (a[i], n);
	if (max == 0)
	  {
	     if (is_malloced) ISIS_FREE (eq);
	     return NULL;
	  }
	eq[i] = 1.0 / max;
     }

   return eq;
}

	
