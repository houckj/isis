/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008 Massachusetts Institute of Technology

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

/*{{{ Includes  */

#include "config.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>

/*}}}*/

#define ISIS_FIT_ENGINE_PRIVATE_DATA \
   double tol; \
   double epsfcn; \
   int maxfev;

#include "isis.h"

#define LMDIF_FC FC_FUNC(lmdif1,LMDIF1)
typedef void lmdif_fun_type
       (int *npts, int *npars, double *pars, double *fvec, int *iflag);
extern void LMDIF_FC(lmdif_fun_type *fun, int *npts, int *npars, double *pars,
                     double *y, double *stat, double *tol, double *epsfcn, int *maxfev,
                     int *info, int *iwa, double *wa, int *lwa);

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
static double *Tmp_Params;
static void *Client_Data;
static int Error_Occurred;

static void fun (int *npts, int *npars, double *pars, double *fvec, int *iflag) /*{{{*/
{
   Isis_Fit_Type *ift = Isis_Fit_Info;
   Isis_Fit_Engine_Type *e = ift->engine;
   Isis_Fit_Statistic_Type *fs = ift->stat;
   Fun_Info_Type *fi = &Fun_Info;
   int i, params_in_range = 1;

   for (i = 0; i < *npars; i++)
     {
        double p = pars[i];

        if (p < e->par_min[i])
          p = e->par_min[i];
        else if (e->par_max[i] < p)
          p = e->par_max[i];

        Tmp_Params[i] = p;
        params_in_range = params_in_range && (p == pars[i]);
     }

   if (params_in_range == 0)
     {
        /* If one or more params is out of range, return a model
         * value guaranteed to be distant from the data */
        double *fx = fi->fx;
        double *y = fi->y;
        double huge = 0.5 * sqrt(DBL_MAX / *npts);
        for (i = 0; i < *npts; i++)
          fx[i] = y[i] + huge;
     }
   else if (-1 == ift->fun (Client_Data, fi->x, fi->npts, Tmp_Params, *npars, fi->fx))
     {
        *iflag = -1;
        Error_Occurred = 1;
        return;
     }

   (void)fs->fun (fs, fi->y, fi->fx, fi->weights, fi->npts, fvec, &ift->statistic);

   if (e->verbose > 0)
     e->verbose_hook (Client_Data, ift->statistic, pars, *npars);
}

/*}}}*/

static void print_status (unsigned int info) /*{{{*/
{
   char *messages[] = {
    "improper input parameters."
  , "algorithm estimates that the relative error in the sum of squares is at most tol."
  , "algorithm estimates that the relative error between x and the solution is at most tol."
  , "conditions for info = 1 and info = 2 both hold."
  , "fvec is orthogonal to the columns of the jacobian to machine precision."
  , "number of calls to fcn has reached or exceeded maxfev."
  , "tol is too small. no further reduction in the sum of squares is possible."
  , "tol is too small. no further improvement in the approximate solution x is possible."
   };
   char *msg;

   if (info == INT_MAX)
     return;

   if (info < sizeof(messages)/sizeof(*messages) )
     msg = messages[info];
   else
     {
        msg = "Unknown info value (internal error)";
        (void) fprintf (stderr, "lmdif1: info = %d\n", info);
     }

   (void) fprintf (stderr, "lmdif1: %s\n", msg);
}

/*}}}*/

static int lmdif (Isis_Fit_Type *ift, void *clientdata, /*{{{*/
                  double *x, double *y, double *weights, unsigned int npts,
                  double *pars, unsigned int npars)
{
   Isis_Fit_Engine_Type *e;
   Fun_Info_Type *fi = &Fun_Info;
   double *fvec = NULL;
   double *wa = NULL;
   int *iwa = NULL;
   int npars_int=npars, npts_int=npts;
   int info = INT_MAX, lwa, status=-1;

   Client_Data = clientdata;

   if (ift == NULL)
     return -1;

   e = ift->engine;
   Isis_Fit_Info = ift;

   fi->x = x;
   fi->y = y;
   fi->weights = weights;
   fi->npts = npts;

   lwa = (npts + 5) * npars + npts;

   fi->fx = NULL;
   Tmp_Params = NULL;

   if ((NULL == (iwa = ISIS_MALLOC (npars * sizeof(int))))
       || (NULL == (Tmp_Params = ISIS_MALLOC (npars * sizeof(double))))
       || (NULL == (wa = ISIS_MALLOC (lwa * sizeof(double))))
       || (NULL == (fi->fx = ISIS_MALLOC (npts * sizeof(double))))
       || (NULL == (fvec = ISIS_MALLOC (npts * sizeof(double)))))
     goto finish;

   Error_Occurred = 0;
   LMDIF_FC(fun,&npts_int,&npars_int,pars,fvec,&ift->statistic,&e->tol,&e->epsfcn,&e->maxfev,&info,iwa,wa,&lwa);

   if (e->verbose > 0)
     print_status (info);

   if ((info >= 1) && (info <= 4) && (Error_Occurred == 0))
     status = 0;
   else status = -1;

   finish:

   ISIS_FREE (wa);
   ISIS_FREE (iwa);
   ISIS_FREE (fvec);
   ISIS_FREE (fi->fx);
   ISIS_FREE (Tmp_Params);

   return status;
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

static int handle_epsfcn_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double epsfcn;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &epsfcn))
     return -1;
   e->epsfcn = epsfcn;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static int handle_maxfev_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   int maxfev;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (1 != sscanf (value, "%d", &maxfev))
     {
        fprintf (stderr, "%s;%s option requires an unsigned int\n", subsystem, optname);
        return -1;
     }
   if (maxfev < 1)
     {
        fprintf (stderr, "%s;%s option must be at least 1\n", subsystem, optname);
        return -1;
     }
   e->maxfev = maxfev;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

static Isis_Option_Table_Type Option_Table [] =
{
     {"tol", handle_tol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-6", "Relative error tolerance"},
     {"epsfcn", handle_epsfcn_option, ISIS_OPT_REQUIRES_VALUE, "1.e-6", "Relative error in function values"},
     {"maxfev", handle_maxfev_option, ISIS_OPT_REQUIRES_VALUE, "500", "Max number of function evaluations"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_options (Isis_Fit_Engine_Type *e, Isis_Option_Type *opts) /*{{{*/
{
   return isis_process_options (opts, Option_Table, (void *)e, 1);
}

/*}}}*/

static int set_range_hook (Isis_Fit_Engine_Type *e, Isis_Fit_Range_Hook_Type *r) /*{{{*/
{
   (void) e; (void) r;
   return 0;
}

/*}}}*/

static void deallocate (Isis_Fit_Engine_Type *e)
{
   ISIS_FREE (e->engine_name);
   ISIS_FREE (e->default_statistic_name);
   ISIS_FREE (e->option_string);
}

ISIS_FIT_ENGINE_METHOD(lmdif,name,sname)
{
   Isis_Fit_Engine_Type *e;

   if (NULL == (e = ISIS_MALLOC (sizeof(Isis_Fit_Engine_Type))))
     return NULL;
   memset ((char *)e, 0, sizeof (*e));

   if ((NULL == (e->engine_name = isis_make_string (name)))
       || (NULL == (e->default_statistic_name = isis_make_string (sname))))
     {
        deallocate (e);
        ISIS_FREE (e);
        return NULL;
     }

   e->method = lmdif;
   e->set_options = set_options;
   e->deallocate = deallocate;
   e->set_range_hook = set_range_hook;
   e->range_hook = NULL;
   e->verbose_hook = verbose_hook;
   e->warn_hook = warn_hook;

   e->tol = 1.e-6;
   e->epsfcn = 1.e-6;
   e->maxfev = 500;

   e->option_string = isis_make_default_option_string ("lmdif", Option_Table);
   if (e->option_string == NULL)
     {
        deallocate (e);
        return NULL;
     }

   return e;
}
