/* -*- mode: C; mode: fold -*- */

/* This routine is based on marq.sl, an S-Lang implementation of
 * Levenberg-Marquardt by John E. Davis.
 * 1/27/99, John C. Houck did the C translation, with modifications
 *
 * Copyright (c) 1999, 1998 John E. Davis
 *
 * You may distribute this file under the terms the GNU General Public
 * License.  See the file COPYING for more information.
 */

/* $Id: marq.c,v 1.15 2003/09/10 22:48:46 houck Exp $ */

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdarg.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#define ISIS_FIT_ENGINE_PRIVATE_DATA \
   double tol;             /* chi-square acceptance tolerance */       \
   double delta;           /* initial dx numerical derivative calc */  \
   double inc;             /* lambda adjustment */                     \
   double dec;             /* lambda adjustment */                     \
   double add;             /* additive damping factor */               \
   double max_step_frac;   /* max relative param step */               \
   unsigned int max_loops;

#include "isis.h"
#include "isismath.h"
#include "util.h"

/* JCH: some features based on Marquardt implementation in
 *  Nash, J.C., 1990, "Compact Numerical Methods for Computers,
 *        Linear Algebra and Function Minimization", 2nd ed.
 *        (Adam Hilger; New York).
 */ 

static int marquardt_range_hook (void *clientdata, double *amin, double *amax, double *a, unsigned int n) /*{{{*/
{
   unsigned int i;

   (void) clientdata;

   for (i=0; i < n; i++)
     {
        if ((a[i] < amin[i]) || (amax[i] < a[i]))
          {
             a[i] = amin[i] + (amax[i] - amin[i]) * urand ();
          }
     }

   return 0;
}

/*}}}*/

static void marquardt_warn_hook (void *clientdata, const char * fmt, ...) /*{{{*/
{
   char buf[256];
   va_list ap;

   (void) clientdata;

   if (fmt == NULL)
     return;

   va_start (ap, fmt);
   if (-1 == isis_vsnprintf (buf, sizeof(buf), fmt, ap))
     fputs ("**** String buffer overflow in marquardt_warn_hook\n", stderr);
   va_end (ap);

   fputs (buf, stdout);
}

/*}}}*/

/*}}}*/

static void marquardt_verbose_hook (void *clientdata, double statistic, /*{{{*/
                                    double *par, unsigned int n)
{
   unsigned int i;
   (void) clientdata;
   fprintf (stdout, "statistic: %e\n", statistic);
   for (i = 0; i < n; i++)
     fprintf (stdout, "\tp[%u]=%e", i, par[i]);
   (void) fputs ("\n", stdout);
}

/*}}}*/

/* decide how long to iterate to get the derivative */
static int howmany_iterations (Isis_Fit_Engine_Type *e) /*{{{*/
{
   double max_delta = 0.5;
   int num_doublings;

   if (e->delta <= 0.0)
     return -1;

   (void) frexp (max_delta/e->delta, &num_doublings);

   return (num_doublings < 5) ? 5 : num_doublings;
}

/*}}}*/

static int _marquardt_compute_alpha_beta (Isis_Fit_Type *ft, void *clientdata, /*{{{*/
                                          double **alpha, double *beta,
                                          double *x, double *y_dat,
                                          double *weight, unsigned int ny,
                                          double *y_a, double *a,
                                          unsigned int nparms)
{
   Isis_Fit_Engine_Type *e = ft->engine;
   Isis_Fit_Statistic_Type *fs = ft->stat;
   double *y_1 = NULL;
   double **dyda = NULL;
   int i, count_max, num_pars = nparms;
   unsigned int k;
   int ret = -1;
   
   if ((count_max = howmany_iterations (e)) < 0)
     return -1;
   
   if (NULL == (dyda = JDMdouble_matrix (nparms,ny))
       || NULL == (y_1 = JDMdouble_vector (ny)))
     goto finish;

   for (i=0; i < num_pars; i++)
     {
        double da, sumsq_diff;
        double ai = a[i];
        double a_min = e->par_min[i];
        double a_max = e->par_max[i];
        int count = 0;

        /* recommended by Nash */
        da = (fabs(ai) + sqrt(e->delta)) * sqrt(e->delta);

        do
          {
             double a_test = ai + da;

             /* parameter value must stay in bounds */
             if (a_test < a_min)
               a_test = a_min;
             else if (a_max < a_test)
               a_test = a_max;

             a[i] = a_test;

             if (-1 == ft->compute_model (fs->opt_data, x, ny, a, nparms, y_1))
               {
                  e->warn_hook (clientdata, "function evaluation failed\n");
                  goto finish;
               }
             
             sumsq_diff = 0.0;
             for (k=0; k < ny; k++)
               {
                  double diff = y_1[k] - y_a[k];
                  dyda[i][k] = diff / da;
                  sumsq_diff += diff * diff;
               }

             a[i] = ai;
             da *= 2;
             count++;
          }
        while (count < count_max && sumsq_diff == 0.0);
     }

   for (i=0; i < num_pars; i++)
     {
        double sum;
        int j;

        sum = 0.0;
        for (k=0; k < ny; k++)
          sum += weight[k] * (y_dat[k] - y_1[k]) * dyda[i][k];

        beta[i] = sum;

        for (j=0; j <= i; j++)
          {
             sum = 0.0;
             for (k=0; k < ny; k++)
               sum += weight[k] * dyda[i][k] * dyda[j][k];

             alpha[i][j] = alpha[j][i] = sum;
          }
     }

   ret = 0;

   finish:

   ISIS_FREE (y_1);
   JDMfree_double_matrix (dyda,nparms);
      
   return ret;
}

/*}}}*/

static int marquardt (Isis_Fit_Type *ft, void *clientdata, /*{{{*/
                      double *x, double *y, double *weight, unsigned int ny,
                      double *a, unsigned int nparms)
{
   Isis_Fit_Engine_Type *e = ft->engine;
   Isis_Fit_Statistic_Type *fs = ft->stat;
   double **alpha = NULL;
   double **cov = NULL;
   double *beta, *a1, *y_1, *alpha_diag, *vec;
   unsigned int *piv = NULL;
   double chisqr, chisqr0, chisqr1, lambda, dchisqr = 0.0;
   unsigned int i, j, k;
   int ret = -1;

   beta = a1 = y_1 = alpha_diag = vec = NULL;

   if (ny == 0 || nparms == 0)
     return -1;

   if (NULL == y || NULL == weight || NULL == a)
     {
        e->warn_hook (clientdata, "marquardt: internal error (null ptr input)\n");
        return -1;
     }

   if (NULL == (alpha_diag = JDMdouble_vector (nparms))
       || NULL == (alpha = JDMdouble_matrix (nparms,nparms))
       || NULL == (cov = JDMdouble_matrix (nparms,nparms))
       || NULL == (piv = (unsigned int *) ISIS_MALLOC (nparms * sizeof (unsigned int)))
       || NULL == (beta = JDMdouble_vector (nparms))
       || NULL == (a1 = JDMdouble_vector (nparms))
       || NULL == (y_1 = JDMdouble_vector (ny))
       || NULL == (vec = JDMdouble_vector (ny)))
     {
        e->warn_hook (clientdata, "marquardt: allocation failed\n");
        goto finish;
     }

   if (-1 == ft->compute_model (fs->opt_data, x, ny, a, nparms, y_1))
     {
        e->warn_hook (clientdata, "function evaluation failed\n");
        goto finish;
     }

   (void) fs->compute_statistic (fs, y, y_1, weight, ny, vec, &chisqr);
   chisqr1 = chisqr;
   chisqr0 = chisqr;

   lambda = 0.1;

   for (i=0; i < e->max_loops; i++)
     {
        if (-1 == _marquardt_compute_alpha_beta (ft, clientdata,
                                                 alpha, beta, x, y, weight, ny,
                                                 y_1, a, nparms))
          {
             e->warn_hook (clientdata, "marquardt:  Failed computing matrix of derivatives d(chisqr)/d(param)\n");
             ft->statistic = chisqr;
             goto finish;
          }

        for (j=0; j < nparms; j++)
          alpha_diag[j] = alpha[j][j];

        for (;;)
          {
             if (e->add == 0.0)
               {
                  /* multiplicative damping */
                  for (j=0; j < nparms; j++)
                    alpha[j][j] = alpha_diag[j] * (1.0 + lambda);
               }
             else
               {
                  /* additive damping */
                  for (j=0; j < nparms; j++)
                    alpha[j][j] = alpha_diag[j] + lambda * e->add;
               }

             /* First try solving with LU decomposition.
              * if that fails, use SVD */
             for (j=0; j < nparms; j++)
               for (k=0; k < nparms; k++)
                 cov[j][k] = alpha[j][k];
             memcpy ((char *)a1, (char *)beta, nparms * sizeof(double));

             if (-1 == isis_lu_solve (cov, nparms, piv, a1))
               {
                  for (j=0; j < nparms; j++)
                    for (k=0; k < nparms; k++)
                      cov[j][k] = alpha[j][k];
                  memcpy ((char *)a1, (char *)beta, nparms * sizeof(double));

                  if (-1 == isis_svd_solve (cov, nparms, a1))
                    {
                       if (e->verbose > 0)
                         e->warn_hook (clientdata, "marquardt: matrix solution failed\n");
                       ft->statistic = chisqr;
                       goto finish;
                    }
               }

             /* apply parameter step, within limits */
             if (e->max_step_frac == 0.0)
               {
                  for (j=0; j < nparms; j++)
                    a1[j] += a[j];
               }
             else
               {
                  for (j=0; j < nparms; j++)
                    {
                       double da_max = (e->max_step_frac
                                        * ((a[j] != 0.0) ? fabs(a[j]) : 1.0));
                       double da = a1[j];
                       if (fabs(da) > da_max)
                         {
                            da = da_max * ((da < 0.0) ? -1.0 : 1.0);
                         }
                       a1[j] = a[j] + da;
                    }
               }

             (void) e->range_hook (e->range_hook_client_data, e->par_min, e->par_max, a1, nparms);

             if (-1 == ft->compute_model (fs->opt_data, x, ny, a1, nparms, y_1))
               {
                  e->warn_hook (clientdata, "function evaluation failed\n");
                  goto finish;
               }

             (void) fs->compute_statistic (fs, y, y_1, weight, ny, vec, &chisqr1);

             if (e->verbose > 0)
               e->verbose_hook (clientdata, chisqr1, a1, nparms);

             if (chisqr1 < chisqr)
               break;

             lambda *= e->inc;
             if (lambda > 1.e10)
               {
                  if (e->verbose > 0)
                    e->warn_hook (clientdata, "marquardt:  Too many iterations:  no decrease in fit-statistic\n");
                  ft->statistic = chisqr;
                  if (chisqr < chisqr0) ret = 0;
                  goto finish;
               }
          }

        lambda *= e->dec;
        if (lambda < DBL_EPSILON)
          lambda = DBL_EPSILON;

        memcpy ((char *)a, (char *)a1, nparms * sizeof(double));

        if (e->verbose > 0)
          e->warn_hook (clientdata, "delta-chisqr = %11.4e   accept = %11.4e\n",
                         chisqr - chisqr1, e->tol * chisqr);

        if (chisqr - chisqr1 <= e->tol * chisqr)
          {
             ft->statistic = chisqr1;
             ret = 0;
             goto finish;
          }

        dchisqr = chisqr - chisqr1;          /* for output only */
        chisqr = chisqr1;
     }

   if (e->verbose > 0)
     e->warn_hook (clientdata, "marquardt:  Failed to converge:  last delta-chisqr = %11.4e\n", dchisqr);

   ft->statistic = chisqr1;

   finish:

#if 0
   if (!ret)
     {
        for (j=0; j < nparms; j++)
          for (k=0; k < nparms; k++)
            cov[j][k] = alpha[j][k];

        if (-1 == JDM_ludecomp_inverse (cov, nparms))
          vmessage ("singular covariance matrix\n");
     }
#endif

   JDMfree_double_matrix (cov, nparms);
   JDMfree_double_matrix (alpha, nparms);
   ISIS_FREE (alpha_diag);
   ISIS_FREE (piv);
   ISIS_FREE (beta);
   ISIS_FREE (a1);
   ISIS_FREE (y_1);
   ISIS_FREE (vec);
   return ret;
}

/*}}}*/

static int handle_double_option (char *subsystem, char *optname, char *value, double *d) /*{{{*/
{
   if (1 != sscanf (value, "%lf", d))
     {
        fprintf (stderr, "%s;%s option requires a double\n", subsystem, optname);
        return -1;
     }
   return 0;
}

/*}}}*/

static int handle_tol_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double tol;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &tol))
     return -1;
   e->tol = tol;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static int handle_delta_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double delta;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &delta))
     return -1;
   e->delta = fabs(delta);
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static int handle_inc_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double inc;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &inc))
     return -1;
   e->inc = inc;
   return isis_update_option_string (&e->option_string, optname, value);   
}

/*}}}*/

static int handle_add_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double add;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &add))
     return -1;
   e->add = add;
   return isis_update_option_string (&e->option_string, optname, value);   
}

/*}}}*/

static int handle_dec_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double dec;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &dec))
     return -1;
   e->dec = dec;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static int handle_max_step_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double max_step;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &max_step))
     return -1;
   e->max_step_frac = max_step;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static int handle_loops_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   unsigned int max_loops;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (1 != sscanf (value, "%u", &max_loops))
     {
        fprintf (stderr, "%s;%s option requires an unsigned int\n", subsystem, optname);
        return -1;
     }
   e->max_loops = max_loops;
   return isis_update_option_string (&e->option_string, optname, value);   
}

/*}}}*/

static Isis_Option_Table_Type Option_Table [] =
{
     {"tol", handle_tol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-4", "Fractional chisqr acceptance tolerance"},
     {"delta", handle_delta_option, ISIS_OPT_REQUIRES_VALUE, "1.e-6", "Initial numerical derivative step size"},
     {"dec", handle_dec_option, ISIS_OPT_REQUIRES_VALUE, "0.1", "Lambda decrease adjustment factor"},   
     {"inc", handle_inc_option, ISIS_OPT_REQUIRES_VALUE, "10.0", "Lambda increase adjustment factor"},
     {"add", handle_add_option, ISIS_OPT_REQUIRES_VALUE, "0.0", "Additive damping factor"},
     {"max_step_frac", handle_max_step_option, ISIS_OPT_REQUIRES_VALUE, "0.0", "Max relative param step"},
     {"max_loops", handle_loops_option, ISIS_OPT_REQUIRES_VALUE, "50", "Max number of iterations"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_options (Isis_Fit_Engine_Type *e, Isis_Option_Type *opts) /*{{{*/
{
   return isis_process_options (opts, Option_Table, (void *)e, 1);
}

/*}}}*/

static int set_range_hook (Isis_Fit_Engine_Type *e, Isis_Fit_Range_Hook_Type *r) /*{{{*/
{
   if (e == NULL)
     return -1;

   e->range_hook = r ? r : &marquardt_range_hook;

   return 0;
}

/*}}}*/

static void deallocate (Isis_Fit_Engine_Type *e)
{
   ISIS_FREE (e->engine_name);
   ISIS_FREE (e->default_statistic_name);
   ISIS_FREE (e->option_string);
}

ISIS_FIT_ENGINE_METHOD(marquardt,name,sname)
{
   Isis_Fit_Engine_Type *e;

   if (NULL == (e = (Isis_Fit_Engine_Type *) ISIS_MALLOC (sizeof(Isis_Fit_Engine_Type))))
     return NULL;
   memset ((char *)e, 0, sizeof (*e));

   if ((NULL == (e->engine_name = isis_make_string (name)))
       || (NULL == (e->default_statistic_name = isis_make_string (sname))))
     {
        deallocate (e);
        ISIS_FREE (e);
        return NULL;
     }

   e->method = marquardt;
   e->set_options = set_options;
   e->deallocate = deallocate;
   e->set_range_hook = set_range_hook;
   e->range_hook = marquardt_range_hook;
   e->verbose_hook = marquardt_verbose_hook;
   e->warn_hook = marquardt_warn_hook;

   /* Private */
   e->tol = 1.e-4;
   e->delta = 1.e-6;
   e->dec = 0.1;
   e->inc = 10.0;
   e->add = 0.0;
   e->max_step_frac = 0.0;
   e->max_loops = 50;

   e->option_string = isis_make_default_option_string ("marquardt", Option_Table);
   if (e->option_string == NULL)
     {
        deallocate (e);
        return NULL;
     }

   return e;
}

