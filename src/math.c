/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2010  Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Author:  John C. Houck  <houck@space.mit.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/* $Id: math.c,v 1.14 2004/09/10 02:31:09 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>

#include "isis.h"
#include "util.h"
#include "_isis.h"
#include "isismath.h"
#include "errors.h"
/*}}}*/

/* quick-select median finder - order of input array is *not* preserved */
static int find_median (double *x, int n, double *xmed) /*{{{*/
{
   int med, lo, hi;

   if (x == NULL || n == 0 || xmed == NULL)
     return -1;

#define DSWAP(i,j) do {double t_=x[(i)]; x[(i)]=x[(j)]; x[(j)]=t_;} while (0)

   lo = 0;
   hi = n-1;
   med = (lo + hi)/2;

   for (;;)
     {
        int l, h, mid;

        /* one item left? */
        if (hi <= lo)
          {
             *xmed = x[med];
             return 0;
          }

        /* two items left? */
        if (hi == lo+1)
          {
             if (x[lo] > x[hi]) DSWAP(lo, hi);
             *xmed = x[med];
             return 0;
          }

        /* order the current (lo,hi) values and their midpoint */
        mid = (lo + hi) / 2;
        if (x[mid] > x[hi]) DSWAP(mid, hi);
        if (x[lo]  > x[hi]) DSWAP(lo, hi);
        if (x[mid] > x[lo]) DSWAP(mid, lo);

        /* swap middle value to a safe place (lo+1)
         * and process the items between [lo+2:hi]
         */
        DSWAP(mid, lo+1);

        /* scan the items between [lo+2, hi-1],
         * inward from each end,
         * swapping when necessary.
         */

        l = lo + 1;
        h = hi;
        for (;;)
          {
             while (x[++l] < x[lo])
               ;
             while (x[lo] < x[--h])
               ;

             if (h < l)
               break;

             DSWAP(l, h);
          }

        /* swap the current low value to the spot
         * where the scan pointers met
         */
        DSWAP(lo, h);

        /* update the endpoints, making sure that the median
         * is bracketed.
         */
        if (h <= med) lo = l;
        if (h >= med) hi = h - 1;
     }
#undef DSWAP
}

/*}}}*/

static void median (void) /*{{{*/
{
   SLang_Array_Type *sx = NULL;
   double med = DBL_MAX;
   double *x;
   int n;

   if ((-1 == SLang_pop_array_of_type (&sx, SLANG_DOUBLE_TYPE))
       || (sx == NULL)
       || (sx->num_elements < 1))
     {
        isis_throw_exception (Isis_Error);
        SLang_free_array (sx);
        return;
     }

   x = (double *)sx->data;
   n = sx->num_elements;

   (void) find_median (x, n, &med);

   SLang_push_double (med);
   SLang_free_array (sx);
}

/*}}}*/

typedef struct
{
   int num;
   double ave;
   double var;
   double sdev;
   double sdom;
   double min;
   double max;
}
Moment_Type;

static SLang_CStruct_Field_Type Moment_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Moment_Type, num, "num", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Moment_Type, ave,  "ave", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Moment_Type, var,  "var", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Moment_Type, sdev, "sdev", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Moment_Type, sdom, "sdom", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Moment_Type, min,  "min", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Moment_Type, max,  "max", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

/* subroutine borrowed from John Davis */
static int mean_stddev_doubles (double *x, unsigned int num, double *s) /*{{{*/
{
   unsigned int i;
   double mean_i, variance_i;

   mean_i = variance_i = 0.0;
   i = 0;
   while (i < num)
     {
        double diff, x_i;

        x_i = x[i];
        diff = x_i - mean_i;
        i++;
        mean_i += diff / i;
        variance_i += diff * (x_i - mean_i);
     }

   s[0] = mean_i;
   if (num > 1)
     s[1] = sqrt (variance_i / (num - 1));
   else
     s[1] = 0;

   return 0;
}

/*}}}*/

static void moment (void) /*{{{*/
{
   SLang_Array_Type *sx = NULL;
   Moment_Type m;
   double s[2];  /* mean, stddev */
   double xmax, xmin;
   double *x;
   int i, n;

   if ((-1 == SLang_pop_array_of_type (&sx, SLANG_DOUBLE_TYPE))
       || (sx == NULL))
     {
        isis_throw_exception (Isis_Error);
        goto error_return;
     }

   if (sx->num_elements < 1)
     {
        m.num = 0;
        m.ave = m.var = m.sdev = m.sdom = m.min = m.max = isis_nan();
        SLang_push_cstruct ((VOID_STAR)&m, Moment_Layout);
        return;
     }

   x = (double *)sx->data;
   n = sx->num_elements;

   m.num = n;
   xmin = xmax = x[0];

   for (i = 0; i < n; i++)
     {
        double x_i = x[i];
        if (x_i > xmax)
          xmax = x_i;
        else if (x_i < xmin)
          xmin = x_i;
     }
   m.min = xmin;
   m.max = xmax;

   if (-1 == mean_stddev_doubles (x, n, s))
     {
        isis_throw_exception (Isis_Error);
        goto error_return;
     }

   m.ave = s[0];
   m.sdev = s[1];
   m.var = m.sdev * m.sdev;
   m.sdom = m.sdev / sqrt(n);

   error_return:

   SLang_free_array (sx);
   SLang_push_cstruct ((VOID_STAR)&m, Moment_Layout);
}

/*}}}*/

static int pop_two_darrays (SLang_Array_Type **a, SLang_Array_Type **b) /*{{{*/
{
   *a = NULL;
   *b = NULL;

   if (-1 == SLang_pop_array_of_type (b, SLANG_DOUBLE_TYPE)
       || -1 == SLang_pop_array_of_type (a, SLANG_DOUBLE_TYPE)
       || (*a == NULL || *b == NULL))
     {
        isis_throw_exception (Isis_Error);
        SLang_free_array (*a);
        SLang_free_array (*b);
        *a = NULL;
        *b = NULL;
        return -1;
     }

   return 0;
}

/*}}}*/

/* reverse index converter from John Davis */
static SLang_Array_Type *convert_reverse_indices (SLindex_Type *r, SLindex_Type num_r, SLindex_Type num_h)
{
   SLang_Array_Type *new_r;
   SLang_Array_Type **new_r_data;
   SLindex_Type i, *lens;

   if (NULL == (new_r = SLang_create_array (SLANG_ARRAY_TYPE, 0, NULL, &num_h, 1)))
     return NULL;

   if (NULL == (lens = (SLindex_Type *)SLmalloc (num_h * sizeof (SLindex_Type))))
     {
        SLang_free_array (new_r);
        return NULL;
     }
   memset ((char *)lens, 0, num_h*sizeof(SLindex_Type));

   for (i = 0; i < num_r; i++)
     {
        SLindex_Type r_i = r[i];

        if (r_i >= 0)
          lens[r_i]++;
     }

   new_r_data = (SLang_Array_Type **) new_r->data;
   for (i = 0; i < num_h; i++)
     {
        if (NULL == (new_r_data[i] = SLang_create_array (SLANG_ARRAY_INDEX_TYPE, 0, NULL, &lens[i], 1)))
          goto return_error;

        lens[i] = 0;
     }

   for (i = 0; i < num_r; i++)
     {
        SLang_Array_Type *at;
        SLindex_Type r_i = r[i];

        if (r_i < 0)
          continue;

        at = new_r_data[r_i];

        ((SLindex_Type *)at->data)[lens[r_i]] = i;
        lens[r_i]++;
     }

   SLfree ((char *)lens);
   return new_r;

   return_error:
   SLfree ((char *) lens);
   SLang_free_array (new_r);
   return NULL;
}

static void make_1d_histogram (int *reverse) /*{{{*/
{
   SLang_Array_Type *v, *lo, *hi, *b, *rev;
   double *xlo, *xhi, *bv;
   unsigned int *num;
   SLindex_Type i, n, nbins;
   SLindex_Type *r = NULL;

   v = lo = hi = b = rev = NULL;

   if ((-1 == pop_two_darrays (&lo, &hi))
       || -1 == SLang_pop_array_of_type (&v, SLANG_DOUBLE_TYPE)
       || (v == NULL))
     goto push_result;

   if (lo->num_elements != hi->num_elements)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto push_result;
     }

   n = v->num_elements;
   nbins = lo->num_elements;

   if (n < 1 || nbins < 1)
     goto push_result;

   if (*reverse == 0)
     r = NULL;
   else
     {
        if (NULL == (r = (SLindex_Type *) ISIS_MALLOC (n * sizeof(SLindex_Type))))
          {
             isis_throw_exception (Isis_Error);
             goto push_result;
          }
        for (i = 0; i < n; i++)
          r[i] = -1;
     }

   if (NULL == (b = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &nbins, 1)))
     {
        isis_throw_exception (Isis_Error);
        goto push_result;
     }

   num = (unsigned int *)b->data;
   memset ((char *)num, 0, nbins * sizeof(unsigned int));

   bv = (double *)v->data;
   xlo = (double *)lo->data;
   xhi = (double *)hi->data;

   /* If the (lo,hi) grid has holes, this algorithm will
    * give the wrong answer because every item will go
    * into a bin.  But what if the grid has holes by
    * accident because it was poorly constructed?
    * Perhaps that is a strong reason to deprecate this
    * interface.
    */

   for (i = 0; i < n; i++)
     {
        double t = bv[i];
        int k = find_bin (t, xlo, xhi, (int) nbins);
        if (k >= 0)
          {
             num[k] += 1;
             if (r != NULL) r[i] = k;
          }
     }

   if ((r != NULL)
       && (NULL == (rev = convert_reverse_indices (r, n, nbins))))
     goto push_result;

   push_result:

   SLang_free_array (v);
   SLang_free_array (hi);
   SLang_free_array (lo);
   ISIS_FREE(r);

   SLang_push_array (b, 1);
   SLang_push_array (rev, 1);
}

/*}}}*/

static void make_2d_histogram (int *reverse) /*{{{*/
{
   SLang_Array_Type *grid_x, *grid_y, *sl_x, *sl_y, *b;
   SLang_Array_Type *rev;
   double *x, *y, *bx, *by;
   double xmax, ymax;
   SLindex_Type *num;
   SLindex_Type dims[2];
   SLindex_Type i, n, nx, ny, nbins;
   SLindex_Type *r = NULL;

   grid_x = grid_y = sl_x = sl_y = b = rev = NULL;

   if (-1 == pop_two_darrays (&grid_x, &grid_y))
     goto push_result;

   /* need at least 1 point */
   if ((-1 == pop_two_darrays (&sl_x, &sl_y))
       || (sl_x->num_elements != sl_y->num_elements)
       || (sl_x->num_elements < 1))
     goto push_result;

   n = sl_x->num_elements;
   nx = grid_x->num_elements;
   ny = grid_y->num_elements;

   if (*reverse == 0)
     r = NULL;
   else
     {
        if (NULL == (r = (SLindex_Type *) ISIS_MALLOC (n * sizeof (SLindex_Type))))
          {
             isis_throw_exception (Isis_Error);
             goto push_result;
          }
        for (i = 0; i < n; i++)
          {
             r[i] = -1;
          }
     }

   dims[0] = nx;
   dims[1] = ny;
   nbins = dims[0] * dims[1];
   if (NULL == (b = SLang_create_array (SLANG_INT_TYPE, 0, NULL, dims, 2)))
     {
        isis_throw_exception (Isis_Error);
        goto push_result;
     }

   num = (SLindex_Type *)b->data;
   memset ((char *)num, 0, nbins * sizeof(SLindex_Type));

   bx = (double *)sl_x->data;
   by = (double *)sl_y->data;
   x = (double *)grid_x->data;
   y = (double *)grid_y->data;

   xmax = x[nx-1];
   ymax = y[ny-1];

   for (i = 0; i < n; i++)
     {
        double b_x = bx[i];
        double b_y = by[i];
        SLindex_Type ix, iy, k;

        if (b_x >= xmax)
          ix = nx-1;
        else if ((ix = find_bin (b_x, x, x+1, nx-1)) < 0)
          continue;

        if (b_y >= ymax)
          iy = ny-1;
        else if ((iy = find_bin (b_y, y, y+1, ny-1)) < 0)
          continue;

        k = iy + ny * ix;

        num[k] += 1;
        if (r != NULL) r[i] = k;
     }

   if ((r != NULL)
       && (NULL == (rev = convert_reverse_indices (r, n, nx*ny))))
     goto push_result;

   push_result:

   SLang_free_array (sl_x);
   SLang_free_array (sl_y);
   SLang_free_array (grid_x);
   SLang_free_array (grid_y);

   ISIS_FREE(r);

   SLang_push_array (b, 1);
   SLang_push_array (rev, 1);
}

/*}}}*/

static void _fft1d (int *isign, double *scaling) /*{{{*/
{
   SLang_Array_Type *re, *im;
   int dims;
   int ndim = 1;

   re = im = NULL;

   if (-1 == SLang_pop_array_of_type (&im, SLANG_DOUBLE_TYPE)
       || im == NULL
       || -1 == SLang_pop_array_of_type (&re, SLANG_DOUBLE_TYPE)
       || re == NULL
       || re->num_elements != im->num_elements
       || abs(*isign) != 1)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid input to FFT");
        goto push_values;
     }

   dims = (int) re->num_elements;

   if (-1 == JDMfftn (ndim, &dims, (double *)re->data,
                      (double *)im->data, *isign, *scaling))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "computing FFT");

   JDMfft_free ();

   push_values:
   (void) SLang_push_array (re, 1);
   (void) SLang_push_array (im, 1);
}

/*}}}*/

static int dsort (const void *v1, const void *v2) /*{{{*/
{
   const double *a = (const double *) v1;
   const double *b = (const double *) v2;

   if (*a < *b) return -1;
   else if (*a > *b) return 1;
   return 0;
}

/*}}}*/

static double ks_difference (void) /*{{{*/
{
   SLang_Array_Type *sl1, *sl2;
   double *a1, *a2;
   int n1, n2, i1, i2;
   double f01, f02, diff;

   sl1 = sl2 = NULL;

   if (-1 == SLang_pop_array_of_type (&sl1, SLANG_DOUBLE_TYPE)
       || sl1 == NULL
       || -1 == SLang_pop_array_of_type (&sl2, SLANG_DOUBLE_TYPE)
       || sl2 == NULL)
     return -1.0;

   a1 = (double *)sl1->data;
   a2 = (double *)sl2->data;

   n1 = sl1->num_elements;
   n2 = sl2->num_elements;

   qsort (a1, n1, sizeof(double), dsort);
   qsort (a2, n2, sizeof(double), dsort);

   i1 = i2 = 0;
   f01 = f02 = 0.0;
   diff = 0.0;

   while ((i1 < n1) && (i2 < n2))
     {
        double fn1, fn2, dt;

        if (a1[i1] < a2[i2])
          {
             fn1 = (double) i1 / n1;
             dt = MAX (fabs(fn1-f02), fabs(f01-f02));
             if (dt > diff) diff = dt;
             f01 = fn1;
             i1++;
          }
        else
          {
             fn2 = (double) i2 / n2;
             dt = MAX (fabs(fn2-f01), fabs(f02-f01));
             if (dt > diff) diff = dt;
             f02 = fn2;
             i2++;
          }
     }

   diff *= sqrt (n1*n2/((double)(n1+n2)));

   SLang_free_array (sl1);
   SLang_free_array (sl2);

   return diff;
}

/*}}}*/

static double ks_probability (double *lambda) /*{{{*/
{
   double eps1, eps2;
   double a2, fac, termbf, p;
   int i;

   eps1 = 1.e-3;
   eps2 = 1.e-8;

   a2 = -2 * (*lambda)*(*lambda);
   fac = 2.0;
   termbf = 0.0;

   p = 0.0;

   for (i = 1; i <= 100; i++)
     {
        double term = fac * exp (a2*i*i);
        p += term;
        if ((fabs(term) < eps1*termbf) || (fabs(term) < eps2*p))
          return p;
        fac *= -1;
        termbf = fabs(term);
     }

   return 1.0;
}

/*}}}*/

static void seed_random (unsigned long *seed) /*{{{*/
{
   random_seed (*seed);
}

/*}}}*/

static void rand_array (SLindex_Type num, double (*rand_fun)(void)) /*{{{*/
{
   SLang_Array_Type *at = NULL;
   double *ad;
   SLindex_Type i;

   if (num <= 0)
     return;
   else if (num == 1)
     {
        SLang_push_double ((*rand_fun)());
        return;
     }

   if (NULL == (at = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &num, 1)))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "creating array of random values");
        return;
     }

   ad = (double *) at->data;
   for (i = 0; i < num; i++)
     ad[i] = (*rand_fun) ();

   SLang_push_array (at, 1);
}

/*}}}*/

static void urand_array (SLindex_Type *num) /*{{{*/
{
   /* Uniform random numbers */
   rand_array (*num, urand);
}

/*}}}*/

static void grand_array (SLindex_Type *num) /*{{{*/
{
   /* Gaussian random numbers */
   rand_array (*num, grand);
}

/*}}}*/

static void prand_vec (void) /*{{{*/
{
   SLang_Array_Type *rate = NULL;
   double *r;
   SLindex_Type i, n;

   if (-1 == SLang_pop_array_of_type (&rate, SLANG_DOUBLE_TYPE)
       || rate == NULL)
     {
        SLang_free_array (rate); rate = NULL;
        SLang_push_array (rate, 1);
        return;
     }

   n = rate->num_elements;
   r = (double *)rate->data;

   for (i = 0; i < n; i++)
     {
        r[i] = prand (r[i]);
     }

   SLang_push_array (rate, 1);
}

/*}}}*/

static void prand_array (double *rate, SLindex_Type *num) /*{{{*/
{
   SLang_Array_Type *at = NULL;
   int *ai;
   SLindex_Type i, n;

   n = *num;

   if (n == 0)
     return;
   else if (n == 1)
     {
        SLang_push_integer ((int) prand (*rate));
        return;
     }

   if (NULL == (at = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &n, 1)))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "creating array of random values");
        return;
     }

   ai = (int *) at->data;
   for (i = 0; i < n; i++)
     ai[i] = prand (*rate);

   SLang_push_array (at, 1);
}

/*}}}*/

/* solving linear systems */

typedef struct
{
   double **a;
   double *b;
   SLindex_Type n;
}
Linear_System_Type;

static int pop_linear_system (Linear_System_Type *t) /*{{{*/
{
   SLang_Array_Type *sl_a=NULL, *sl_b=NULL;
   double **a=NULL, *b=NULL;
   SLindex_Type i, j, n, dims[2];
   int status = -1;

   t->a = NULL;
   t->b = NULL;
   t->n = 0;

   if ((-1 == SLang_pop_array_of_type (&sl_b, SLANG_DOUBLE_TYPE))
       || (-1 == SLang_pop_array_of_type (&sl_a, SLANG_DOUBLE_TYPE)))
     goto return_error;

   n = sl_b->num_elements;

   if (sl_a->num_elements != (unsigned int) n*n)
     goto return_error;

   if ((NULL == (a = JDMdouble_matrix (n, n)))
       || (NULL == (b = JDMdouble_vector (n))))
     goto return_error;

   memcpy ((char *)b, (char *)sl_b->data, n * sizeof(double));
   for (i = 0; i < n; i++)
     {
        double *ai = a[i];
        if (-1 == SLang_get_array_element (sl_b, &i, &b[i]))
          goto return_error;
        dims[0] = i;
        for (j = 0; j < n; j++)
          {
             double aij;
             dims[1] = j;
             if (-1 == SLang_get_array_element (sl_a, dims, &aij))
               goto return_error;
             ai[j] = aij;
          }
     }

   t->a = a;
   t->b = b;
   t->n = n;

   status = 0;
return_error:
   SLang_free_array (sl_a);
   SLang_free_array (sl_b);
   return status;
}

/*}}}*/

static void free_linear_system (Linear_System_Type *t)
{
   if (t == NULL)
     return;
   JDMfree_double_matrix (t->a, t->n);
   t->a = NULL;
   ISIS_FREE(t->b);
   t->n = 0;
}

static void lu_solve_intrin (void)
{
   Linear_System_Type t;
   SLang_Array_Type *sl_b = NULL;
   unsigned int *piv = NULL;

   if ((-1 == pop_linear_system (&t))
       || (NULL == (piv = (unsigned int *) ISIS_MALLOC (t.n * sizeof(unsigned int)))))
     {
        isis_throw_exception (Isis_Error);
        goto the_return;
     }

   if (-1 == isis_lu_solve (t.a, t.n, piv, t.b))
     goto the_return;

   sl_b = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &t.n, 1);
   if (sl_b != NULL)
     {
        memcpy ((char *)sl_b->data, (char *)t.b, t.n * sizeof (double));
     }

the_return:
   SLang_push_array (sl_b, 1);
   free_linear_system (&t);
   ISIS_FREE(piv);
}

static void svd_solve_intrin (void)
{
   Linear_System_Type t;
   SLang_Array_Type *sl_b = NULL;

   if (-1 == pop_linear_system (&t))
     {
        isis_throw_exception (Isis_Error);
        goto the_return;
     }

   if (-1 == isis_svd_solve (t.a, t.n, t.b))
     goto the_return;

   sl_b = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &t.n, 1);
   if (sl_b != NULL)
     {
        memcpy ((char *)sl_b->data, (char *)t.b, t.n * sizeof (double));
     }

the_return:
   SLang_push_array (sl_b, 1);
   free_linear_system (&t);
}

/*{{{ intrinsics */

#define V SLANG_VOID_TYPE
#define AI SLANG_ARRAY_INDEX_TYPE
#define I SLANG_INT_TYPE
#define U SLANG_UINT_TYPE
#define F SLANG_FLOAT_TYPE
#define D SLANG_DOUBLE_TYPE
#define S SLANG_STRING_TYPE

static SLang_Intrin_Fun_Type Math_Intrinsics [] =
{
   MAKE_INTRINSIC_1("_make_1d_histogram", make_1d_histogram, V, I),
   MAKE_INTRINSIC_1("_make_2d_histogram", make_2d_histogram, V, I),
   MAKE_INTRINSIC_2("_fft1d", _fft1d, V, I, D),
   MAKE_INTRINSIC("_moment", moment, V, 0),
   MAKE_INTRINSIC("_median", median, V, 0),
   MAKE_INTRINSIC("_ks_difference", ks_difference, D, 0),
   MAKE_INTRINSIC_1("_ks_probability", ks_probability, D, D),
   MAKE_INTRINSIC_1("_seed_random", seed_random, V, SLANG_ULONG_TYPE),
   MAKE_INTRINSIC_1("_urand_array", urand_array, V, AI),
   MAKE_INTRINSIC_1("_grand_array", grand_array, V, AI),
   MAKE_INTRINSIC_2("_prand_array", prand_array, V, D, AI),
   MAKE_INTRINSIC("_prand_vec", prand_vec, V, 0),
   MAKE_INTRINSIC("lu_solve_intrin", lu_solve_intrin, V, 0),
   MAKE_INTRINSIC("svd_solve_intrin", svd_solve_intrin, V, 0),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef V
#undef AI
#undef I
#undef F
#undef D
#undef S
#undef U

SLANG_MODULE(math);
int init_math_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if (-1 == SLns_add_intrin_fun_table (ns, Math_Intrinsics, NULL))
     return isis_trace_return(-1);

   return 0;
}

void deinit_math_module (void)
{
}

/*}}}*/
