/* -*- mode: C; mode: fold -*- */
/*
    This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2009 Massachusetts Institute of Technology

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

/* $Id: subplex.c,v 1.6 2004/02/09 11:14:25 houck Exp $ */

/*{{{ Includes  */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdarg.h>
#include <string.h>

#include <setjmp.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

/*}}}*/

#define ISIS_FIT_ENGINE_PRIVATE_DATA \
   double *scale;       \
   double scale_factor; \
   double tol;          \
   int mode;            \
   int maxnfe;          \
   int nfe;             \
   int iflag;

#define PARAM_OUTSIDE_BOUND_CONSTRAINTS 1.e30

#include "isis.h"
#include "isismath.h"

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
}
#endif
#define SUBPLX_FC FC_FUNC(subplx,SUBPLX)
typedef double subplx_fun_type (int *num, double *pars);
extern void SUBPLX_FC(subplx_fun_type *f, int *n, double *tol, int *maxnfe,
                      int *mode, double *scale, double *x, double *fx,
                      int *nfe, double *work, int *iwork, int *iflag);
#if 0
{
#endif
#ifdef __cplusplus
}
#endif

typedef struct
{
   double *x;
   double *y;
   double *weights;
   double *fx;
   unsigned int npts;
}
Fun_Info_Type;

static Fun_Info_Type Fun_Info;
static Isis_Fit_Type *Isis_Fit_Info;
static double *Stat_Vec;
static void *Client_Data;

static jmp_buf Jumpbuf;

static double fun (int *num, double *pars) /*{{{*/
{
   Isis_Fit_Type *ift = Isis_Fit_Info;
   Isis_Fit_Engine_Type *e = ift->engine;
   Isis_Fit_Statistic_Type *fs = ift->stat;
   Fun_Info_Type *fi = &Fun_Info;
   void *clientdata = Client_Data;
   double statistic;
   unsigned int npars = (unsigned int) *num;

   if (e->range_hook (clientdata, e->par_min, e->par_max, pars, npars))
     return PARAM_OUTSIDE_BOUND_CONSTRAINTS;

   if (-1 == ift->fun (clientdata, fi->x, fi->npts, pars, npars, fi->fx))
     {
        /* subplex doesn't test for error returns */
        longjmp (Jumpbuf, 1);
        /* return -1.0; */
     }

   (void)fs->fun (fs, fi->y, fi->fx, fi->weights, fi->npts, Stat_Vec, &statistic);

   if (e->verbose > 0)
     e->verbose_hook (clientdata, statistic, pars, npars);

   return statistic;
}

/*}}}*/

static int alloc_work (unsigned int npars, double **work, int **iwork) /*{{{*/
{
   unsigned int nsmin, nsmax, worksize, iworksize, i;

   nsmin = (npars < 2) ? npars : 2;
   nsmax = (npars < 5) ? npars : 5;

   worksize = 2*npars + nsmax * (nsmax + 4) + 1;
   *work = (double *) malloc (worksize * sizeof(double));

   iworksize = npars + npars / nsmin;
   *iwork = (int *) malloc (iworksize * sizeof(int));

   if ((*work == NULL) || (*iwork == NULL))
     {
        free (*work);
        free (*iwork);
        return -1;
     }

   /* Unfortunately, subplex seems to assume these arrays
    * have been initialized to zero -- iwork anyway.
    */
   for (i = 0; i < worksize; i++)
     (*work)[i] = 0.0;
   for (i = 0; i < iworksize; i++)
     (*iwork)[i] = 0;

   return 0;
}

/*}}}*/

static void print_status (int status, int nfe) /*{{{*/
{
   char *msg;

   switch (status)
     {
      case -2: msg = "invalid input";
        break;
      case -1: msg = "too many function evaluations";
        break;
      case 0:  msg = "succeeded";
        break;
      case 1: msg = "result limited by machine precision";
        break;
      case 2: msg = "fstop reached";
        break;
      default:
        msg = "internal error; unknown exit status code from subplex";
        break;
     }

   fprintf (stdout, "subplex[%d]:  %s\n", nfe, msg);
}

/*}}}*/

static void set_initial_stepsizes (Isis_Fit_Engine_Type *e, /*{{{*/
                                  double *pars, unsigned int npars)
{
   double *step = e->par_step;
   unsigned int i;

   for (i = 0; i < npars; i++)
     {
        double scale;

        if (step[i] > 0)
          {
             scale = step[i];
          }
        else
          {
             double a = fabs (pars[i]);
             double s = (a > 100.0 * DBL_EPSILON) ? a : 1.0;
             scale = s * e->scale_factor;
          }

        e->scale[i] = scale;
     }
}

/*}}}*/

static int call_subplx (Isis_Fit_Engine_Type *e, double *pars, int num_pars, double *work, int *iwork, double *st) /*{{{*/
{
   double result;

   /* subplex doesnt check function return values */
   if (setjmp (Jumpbuf))
     return -1;

   SUBPLX_FC(fun,&num_pars,&e->tol,&e->maxnfe,&e->mode,e->scale,pars,&result,&e->nfe,work,iwork,&e->iflag);

   *st = result;

   return 0;
}

/*}}}*/

static int subplex (Isis_Fit_Type *ift, void *clientdata, /*{{{*/
                    double *x, double *y, double *weights, unsigned int npts,
                    double *pars, unsigned int npars)
{
   Isis_Fit_Engine_Type *e;
   Fun_Info_Type *fi = &Fun_Info;
   double *work = NULL;
   int *iwork = NULL;
   int num_pars = (int) npars;
   volatile int ret = -1;

   Client_Data = clientdata;

   if (ift == NULL)
     return -1;

   if (npars == 0)
     return -1;

   e = ift->engine;
   Isis_Fit_Info = ift;

   fi->x = x;
   fi->y = y;
   fi->weights = weights;
   fi->npts = npts;

   if ((NULL == (fi->fx = (double *) malloc (npts * sizeof(double))))
       || (NULL == (Stat_Vec = (double *) malloc (npts * sizeof(double)))))
     goto finish;
   memset ((char *)fi->fx, 0, npts * sizeof(double));

   if (NULL == (e->scale = (double *) malloc (npars * sizeof(double))))
     goto finish;

   set_initial_stepsizes (e, pars, npars);

   if (-1 == alloc_work (npars, &work, &iwork))
     goto finish;

   e->iflag = 0;
   if (e->maxnfe <= 0)
     e->maxnfe = 1024 * npars;

   if (-1 == call_subplx (e, pars, num_pars, work, iwork, &ift->statistic))
     goto finish;

   if (e->verbose > 0)
     print_status (e->iflag, e->nfe);

   ret = 0;
   finish:

   free (fi->fx);
   free (e->scale);
   free (work);
   free (iwork);
   free (Stat_Vec);

   if ((ret != 0) || (e->iflag < 0))
     return -1;

   return 0;
}

/*}}}*/

static void warn_hook (void *clientdata, const char * fmt, ...) /*{{{*/
{
   char buf[256];
   va_list ap;

   (void) clientdata;

   if (fmt == NULL)
     return;

   va_start (ap, fmt);
   if (-1 == isis_vsnprintf (buf, sizeof(buf), fmt, ap))
     fputs ("**** String buffer overflow in warn_hook\n", stderr);
   va_end (ap);

   fputs (buf, stdout);
}

/*}}}*/

static void verbose_hook (void *clientdata, double statistic, /*{{{*/
                          double *par, unsigned int n)
{
   unsigned int i;
   (void) clientdata;
   fprintf (stdout, "statistic: %e", statistic);
   for (i = 0; i < n; i++)
     fprintf (stdout, "\tp[%u]=%e", i, par[i]);
   (void) fputs ("\n", stdout);
}

/*}}}*/

static int range_hook (void *clientdata, double *amin, double *amax, double *a, unsigned int n) /*{{{*/
{
   unsigned int i;

   (void) clientdata;

   for (i = 0; i < n; i++)
     {
        if ((a[i] < amin[i]) || (amax[i] < a[i]))
          return 1;
     }

   return 0;
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

static int handle_scale_factor_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double scale_factor;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &scale_factor))
     return -1;
   e->scale_factor = scale_factor;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static int handle_maxnfe_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   int maxnfe;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (1 != sscanf (value, "%d", &maxnfe))
     {
        fprintf (stderr, "%s;%s option requires an unsigned int\n", subsystem, optname);
        return -1;
     }
   e->maxnfe = maxnfe;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static Isis_Option_Table_Type Option_Table [] =
{
     {"tol", handle_tol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-4", "Relative error tolerance"},
     {"scale_factor", handle_scale_factor_option, ISIS_OPT_REQUIRES_VALUE, "1.0", "Parameter scale factor"},
     {"maxnfe", handle_maxnfe_option, ISIS_OPT_REQUIRES_VALUE, "0", "Max number of function evaluations"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_options (Isis_Fit_Engine_Type *e, Isis_Option_Type *opts)
{
   return isis_process_options (opts, Option_Table, (void *)e, 1);
}

static int set_range_hook (Isis_Fit_Engine_Type *e, Isis_Fit_Range_Hook_Type *r)
{
   if (e == NULL)
     return -1;

   e->range_hook = r ? r : &range_hook;

   return 0;
}

static void deallocate (Isis_Fit_Engine_Type *e)
{
   ISIS_FREE (e->engine_name);
   ISIS_FREE (e->default_statistic_name);
   ISIS_FREE (e->option_string);
}

ISIS_FIT_ENGINE_METHOD(subplex,name,sname)
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

   e->method = subplex;
   e->set_options = set_options;
   e->deallocate = deallocate;
   e->set_range_hook = set_range_hook;
   e->range_hook = range_hook;
   e->verbose_hook = verbose_hook;
   e->warn_hook = warn_hook;

   e->scale = NULL;
   e->scale_factor = 1.0;
   e->tol = 1.e-4;
   e->maxnfe = 0;
   e->mode = 0;
   e->nfe = 0;
   e->iflag = 0;

   e->option_string = isis_make_default_option_string ("subplex", Option_Table);
   if (e->option_string == NULL)
     {
        deallocate (e);
        return NULL;
     }

   return e;
}
