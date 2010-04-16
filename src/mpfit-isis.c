/* -*- mode: C; mode: fold -*- */

/*
    Copyright (C) 1998-2010 Massachusetts Institute of Technology

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

/* #include "config.h" */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>

/*}}}*/

#include "mpfit.h"

#define ISIS_FIT_ENGINE_PRIVATE_DATA \
   struct mp_config_struct mpfit_config;

#include "isis.h"

static void *Isis_Client_Data;

typedef struct
{
   double *x;
   double *y;
   double *weights;
   double *fx;
   unsigned int npts;
   unsigned int npars;
}
Fun_Info_Type;
static Fun_Info_Type Fun_Info;

static int mpfit_objective (int num_fvec_values, int num_pars, double *pars, /*{{{*/
                            double *fvec, double **deriv_vec, void *client_data)
{
   Isis_Fit_Type *ift = (Isis_Fit_Type *)client_data;
   Isis_Fit_Engine_Type *e = ift->engine;
   Isis_Fit_Statistic_Type *fs = ift->stat;
   Fun_Info_Type *fi = &Fun_Info;

   (void) num_fvec_values;  (void) deriv_vec;

   if (-1 == ift->fun (Isis_Client_Data, fi->x, fi->npts, pars, num_pars, fi->fx))
     return -1;

   (void)fs->fun (fs, fi->y, fi->fx, fi->weights, fi->npts, fvec, &ift->statistic);

   if (e->verbose > 0)
     e->verbose_hook (Isis_Client_Data, ift->statistic, pars, fi->npars);

   return 0;
}

/*}}}*/

static int mpfit_config_set_defaults (struct mp_config_struct *s) /*{{{*/
{
   if (s == NULL)
     return -1;

   s->ftol = 1.e-10;
   /* Relative chi-square convergence criterion */

   s->xtol = 1.e-10;
   /* Relative parameter convergence criterion */

   s->gtol = 1.e-10;
   /* Orthogonality convergence criterion */

   s->stepfactor = 100.0;
   /* Initial step bound */

   s->epsfcn = MP_MACHEP0;
   /* Finite derivative step size */

   s->covtol = 1.e-14;
   /* Range tolerance for covariance calculation */

   s->maxiter = 200;
   /* Maximum number of iterations.  If maxiter == 0,
    then basic error checking is done, and parameter
    errors/covariances are estimated based on input
    parameter values, but no fitting iterations are done.
    */

   s->maxfev = 0;
   /* Maximum number of function evaluations */

   s->nprint = 0;

   s->douserscale = 0;
   /* Scale variables by user values?
    1 = yes, user scale values in diag;
    0 = no, variables scaled internally */

   s->nofinitecheck = 0;
   /* Disable check for infinite quantities from user?
    0 = do not perform check
    1 = perform check
    */

   s->iterproc = 0;
   /* Placeholder pointer - must set to 0 */

   return 0;
}

/*}}}*/

static int mpfit_method (Isis_Fit_Type *ift, void *clientdata, /*{{{*/
                          double *x, double *y, double *weights, unsigned int npts,
                          double *pars, unsigned int npars)
{
   Isis_Fit_Engine_Type *e;
   struct mp_result_struct mpfit_result;
   struct mp_par_struct *mpfit_pars;
   Fun_Info_Type *fi = &Fun_Info;
   unsigned int i;

   Isis_Client_Data = clientdata;

   if (ift == NULL)
     return -1;

   e = ift->engine;

   if (NULL == (mpfit_pars = (struct mp_par_struct *) malloc (npars * sizeof *mpfit_pars)))
     return -1;

   for (i = 0; i < npars; i++)
     {
        struct mp_par_struct *ps = &mpfit_pars[i];
        ps->fixed = 0;
        ps->limited[0] = e->par_min[i] > -DBL_MAX;
        ps->limited[1] = e->par_max[i] <  DBL_MAX;
        ps->limits[0] = e->par_min[i];
        ps->limits[1] = e->par_max[i];
        ps->parname = 0;
        ps->step = e->par_step[i];
        ps->relstep = 0.01;  /* FIXME !! */
        ps->side = 0;        /* FIXME !! */
        /* Note that mpfit's two-sided numerical derivative option
         * does not fully support parameter bounds.  The method used
         * to support parameter bounds with one-sided numerical
         * derivatives could easily be extended to the two-sided
         * derivative case, but as of 4/2010 that hasn't been done.
         */
        ps->deriv_debug = 0;
        ps->deriv_reltol = 0.0; /* rel. tolerance for deriv debug printout */
        ps->deriv_abstol = 0.0; /* abs. tolerance for deriv debug printout */
     }

   fi->x = x;
   fi->y = y;
   fi->weights = weights;
   fi->npts = npts;
   fi->npars = npars;
   fi->fx = NULL;

   memset ((void *)&mpfit_result, 0, sizeof (struct mp_result_struct));

   if (NULL == (fi->fx = malloc (npts * sizeof(double))))
     goto finish;

   (void) mpfit (mpfit_objective, npts, npars, pars,
                 mpfit_pars, &e->mpfit_config, ift, &mpfit_result);
   ift->statistic = mpfit_result.bestnorm;

   finish:

   free(fi->fx);
   free(mpfit_pars);

   return (mpfit_result.status > 0) ? 0 : -1;
}

/*}}}*/

static void warn_hook (void *clientdata, const char * fmt, ...) /*{{{*/
{
   char buf[1024];
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

static int handle_int_option (char *subsystem, char *optname, char *value, int *d) /*{{{*/
{
   if (1 != sscanf (value, "%d", d))
     {
        fprintf (stderr, "%s;%s option requires a int\n", subsystem, optname);
        return -1;
     }
   return 0;
}

/*}}}*/

#define OPT_HANDLER(name,type) \
  static int handle_##name##_option (char *subsystem, char *optname, char *value, void *clientdata) \
  { \
     Isis_Fit_Engine_Type *e; \
     e = (Isis_Fit_Engine_Type *) clientdata; \
     if ((-1 == handle_##type##_option (subsystem, optname, value, &e->mpfit_config.name)) \
          || (isis_update_option_string (&e->option_string, optname, value))) \
       return -1; \
     return 0; \
  }

OPT_HANDLER(ftol,double)
OPT_HANDLER(xtol,double)
OPT_HANDLER(gtol,double)
OPT_HANDLER(stepfactor,double)
OPT_HANDLER(epsfcn,double)
OPT_HANDLER(covtol,double)
OPT_HANDLER(maxiter,int)
OPT_HANDLER(maxfev,int)

static Isis_Option_Table_Type Option_Table [] =
{
   {"ftol", handle_ftol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-10", "Relative chi-square convergence criterion"},
   {"xtol", handle_xtol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-10", "Relative parameter convergence criterion"},
   {"gtol", handle_gtol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-10", "Orthogonality convergence criterion"},
   {"stepfactor", handle_stepfactor_option, ISIS_OPT_REQUIRES_VALUE, "100.0", "Initial step bound"},
   {"epsfcn", handle_epsfcn_option, ISIS_OPT_REQUIRES_VALUE, "2.2204460e-16", "Finite derivative step size"},
   {"covtol", handle_covtol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-14", "Range tolerance for covariance calculation"},
   {"maxiter", handle_maxiter_option, ISIS_OPT_REQUIRES_VALUE, "200", "Maximum number of iterations"},
   {"maxfev", handle_maxfev_option, ISIS_OPT_REQUIRES_VALUE, "0", "Maximum number of function evaluations"},
   ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_options (Isis_Fit_Engine_Type *e, Isis_Option_Type *opts) /*{{{*/
{
   return isis_process_options (opts, Option_Table, (void *)e, 1);
}

/*}}}*/

static int set_range_hook (Isis_Fit_Engine_Type *e, Isis_Fit_Range_Hook_Type r) /*{{{*/
{
   (void) e; (void) r;
   return 0;
}

/*}}}*/

static void deallocate (Isis_Fit_Engine_Type *e) /*{{{*/
{
   ISIS_FREE (e->engine_name);
   ISIS_FREE (e->default_statistic_name);
   ISIS_FREE (e->option_string);
}

/*}}}*/

ISIS_FIT_ENGINE_METHOD(mpfit,name,sname)
{
   Isis_Fit_Engine_Type *e;

   if (NULL == (e = ISIS_MALLOC (sizeof(Isis_Fit_Engine_Type))))
     return NULL;
   memset ((char *)e, 0, sizeof (*e));

   if ((NULL == (e->engine_name = isis_make_string (name)))
       || (NULL == (e->default_statistic_name = isis_make_string (sname))))
     {
        deallocate(e);
        ISIS_FREE(e);
        return NULL;
     }

   e->method = mpfit_method;
   e->set_options = set_options;
   e->deallocate = deallocate;
   e->set_range_hook = set_range_hook;
   e->range_hook = NULL;
   e->verbose_hook = verbose_hook;
   e->warn_hook = warn_hook;

   mpfit_config_set_defaults (&e->mpfit_config);

   e->option_string = isis_make_default_option_string ("mpfit", Option_Table);
   if (e->option_string == NULL)
     {
        deallocate (e);
        ISIS_FREE(e);
        return NULL;
     }

   return e;
}

static int handle_tol_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   double tol;
   e = (Isis_Fit_Engine_Type *) clientdata;
   if (-1 == handle_double_option (subsystem, optname, value, &tol))
     return -1;
   e->mpfit_config.ftol = tol;
   e->mpfit_config.xtol = tol;
   return isis_update_option_string (&e->option_string, optname, value);
}

/*}}}*/

/* An lmdif interface is provided for back-compatibility */ /*{{{*/

static Isis_Option_Table_Type Lmdif_Option_Table [] =
{
   {"tol", handle_tol_option, ISIS_OPT_REQUIRES_VALUE, "1.e-6", "Relative error tolerance"},
   {"epsfcn", handle_epsfcn_option, ISIS_OPT_REQUIRES_VALUE, "1.e-6", "Relative error in function values"},
   {"maxfev", handle_maxfev_option, ISIS_OPT_REQUIRES_VALUE, "500", "Max number of function evaluations"},
   ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_lmdif_options (Isis_Fit_Engine_Type *e, Isis_Option_Type *opts) /*{{{*/
{
   return isis_process_options (opts, Lmdif_Option_Table, (void *)e, 1);
}

/*}}}*/

ISIS_FIT_ENGINE_METHOD(lmdif,name,sname)
{
   Isis_Fit_Engine_Type *e;

   if (NULL == (e = ISIS_MALLOC (sizeof(Isis_Fit_Engine_Type))))
     return NULL;
   memset ((char *)e, 0, sizeof (*e));

   if ((NULL == (e->engine_name = isis_make_string (name)))
       || (NULL == (e->default_statistic_name = isis_make_string (sname))))
     {
        deallocate(e);
        ISIS_FREE(e);
        return NULL;
     }

   e->method = mpfit_method;
   e->set_options = set_lmdif_options;
   e->deallocate = deallocate;
   e->set_range_hook = set_range_hook;
   e->range_hook = NULL;
   e->verbose_hook = verbose_hook;
   e->warn_hook = warn_hook;

   mpfit_config_set_defaults (&e->mpfit_config);

   e->mpfit_config.ftol = 1.e-6;
   e->mpfit_config.xtol = 1.e-6;
   e->mpfit_config.gtol = 0.0;
   e->mpfit_config.epsfcn = 1.e-6;
   e->mpfit_config.maxfev = 500;

   e->option_string = isis_make_default_option_string ("lmdif", Lmdif_Option_Table);
   if (e->option_string == NULL)
     {
        deallocate (e);
        ISIS_FREE(e);
        return NULL;
     }

   return e;
}

/*}}}*/
