/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2005 Massachusetts Institute of Technology

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

/* $Id$ */

/*{{{ Includes  */

#include "config.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#include <setjmp.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <string.h>

/*}}}*/

#define ISIS_FIT_ENGINE_PRIVATE_DATA \
   double initial_step, eps, t, rt; \
   int ns, nt, neps, maxevl, iprint; \
   int iseed1, iseed2;

#include "isis.h"

#define SIMANN_FC FC_FUNC(simann,SIMANN)
typedef void simann_fun_type (int npars, double *pars, double *val);
extern void SIMANN_FC (simann_fun_type *fcn,
                       int *n, double *x, int *max, double *rt, double *eps,
                       int *ns, int *nt, int *neps, int *maxevl,
                       double *lb, double *ub, double *c,
                       int *iprint, int *iseed1, int *iseed2,
                       double *t, double *vm, double *xopt, double *fopt,
                       int *nacc, int *nfcnev, int *nobds, int *ier,
                       double *fstar, double *xp, int *nacp);

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
static Isis_Fit_Type *Isis_Fit_Info;
static void *Client_Data;
static int Error_Occurred;
static double *Fvec;

static jmp_buf Jumpbuf;

static void fcn (int npars, double *pars, double *val) /*{{{*/
{
   Isis_Fit_Type *ift = Isis_Fit_Info;
   Isis_Fit_Engine_Type *e = ift->engine;
   Isis_Fit_Statistic_Type *fs = ift->stat;
   Fun_Info_Type *fi = &Fun_Info;
   double statistic;

   (void) npars;

   *val = DBL_MAX;

   if (-1 == ift->fun (Client_Data, fi->x, fi->npts, pars, fi->npars, fi->fx))
     {
        Error_Occurred = 1;
        /* simann doesn't test for error returns */
        longjmp (Jumpbuf, 1);
        return;
     }

   (void)fs->fun (fs, fi->y, fi->fx, fi->weights, fi->npts, Fvec, &statistic);
   *val = statistic;

   if (e->verbose > 0)
     e->verbose_hook (Client_Data, statistic, pars, fi->npars);
}

/*}}}*/

static void print_status (unsigned int info) /*{{{*/
{
   static const char *messages[] = {
      "normal return"
     ,"exceeded maximum allowed number of function evaluations"
     ,"starting parameter values lie outside specified bounds"
     ,"initial temperature is not positive"
   };
   const char *msg;

   if (info == INT_MAX)
     return;

   if (info < sizeof(messages)/sizeof(*messages) )
     msg = messages[info];
   else
     {
        msg = "*** simann internal error";
        (void) fprintf (stderr, "simann: ier = %d\n", info);
     }

   (void) fprintf (stderr, "simann: %s\n", msg);
}

/*}}}*/

static int call_simann (Isis_Fit_Engine_Type *e, double *pars, int npars,
                        int nt, double *c, double *vm, double *xopt, double *fopt,
                        int *nacc, int *nfcnev, int *nobds, int *ier, double *fstar,
                        double *xp, int *nacp)
{
   double x_fopt;
   int x_nacc, x_nfcnev, x_nobds, x_ier, zero=0;

   /* simann doesnt check function return values */
   if (setjmp (Jumpbuf))
     return -1;

   SIMANN_FC(fcn,&npars,pars,&zero,&e->rt,&e->eps,&e->ns,&nt,&e->neps,&e->maxevl,
             e->par_min,e->par_max,c,&e->iprint,&e->iseed1,&e->iseed2,
             &e->t,vm,xopt,&x_fopt,&x_nacc,&x_nfcnev,&x_nobds,&x_ier,fstar,xp,nacp);

   *fopt = x_fopt;
   *nacc = x_nacc;
   *nfcnev = x_nfcnev;
   *nobds = x_nobds;
   *ier = x_ier;

   return 0;
}

static int simann (Isis_Fit_Type *ift, void *clientdata, /*{{{*/
                  double *x, double *y, double *weights, unsigned int npts,
                  double *pars, unsigned int npars)
{
   Isis_Fit_Engine_Type *e;
   Fun_Info_Type *fi = &Fun_Info;
   double *fstar=NULL, *xp=NULL, *vm=NULL, *xopt=NULL;
   double *c=NULL;
   int *nacp=NULL;
   double fopt;
   int nt, nacc, nfcnev, nobds, ier;
   unsigned int i;
   int ret=-1;

   Client_Data = clientdata;

   if (ift == NULL)
     return -1;

   e = ift->engine;
   Isis_Fit_Info = ift;

   fi->x = x;
   fi->y = y;
   fi->weights = weights;
   fi->npts = npts;
   fi->npars = npars;
   fi->fx = NULL;

   if ((NULL == (Fvec = malloc (npts * sizeof(double))))
       || (NULL == (fi->fx = malloc (npts * sizeof(double))))
       || (NULL == (vm = malloc (npars * sizeof(double))))
       || (NULL == (xopt = malloc (npars * sizeof(double))))
       || (NULL == (c = malloc (npars * sizeof(double))))
       || (NULL == (xp = malloc (npars * sizeof(double))))
       || (NULL == (fstar = malloc (e->neps * sizeof(double))))
       || (NULL == (nacp = malloc (npars * sizeof(int)))))
     goto finish;

   if (e->nt > 0)
     nt = e->nt;
   else nt = (5 * npars < 100) ? 100 : 5 * npars;

   for (i = 0; i < npars; i++)
     {
        c[i] = 2.0;
     }

   if (e->initial_step > 0)
     {
        for (i = 0; i < npars; i++)
          {
             vm[i] = (fabs(pars[i]) + sqrt(e->initial_step))
               * sqrt(e->initial_step);
          }
     }
   else
     {
        for (i = 0; i < npars; i++)
          {
             double lo, hi;
             lo = pars[i] - e->par_min[i];
             hi = e->par_max[i] - pars[i];
             vm[i] = (lo < hi) ? lo : hi;
          }
     }

   Error_Occurred = 0;
   if (-1 == call_simann (e, pars, npars, nt, c, vm, xopt, &fopt, &nacc, &nfcnev,
                          &nobds, &ier, fstar, xp, nacp))
     goto finish;

   if (e->verbose > 0)
     print_status (ier);

   ret = ((ier == 0) && (Error_Occurred == 0)) ? 0 : -1;

   for (i = 0; i < npars; i++)
     {
        pars[i] = xopt[i];
     }

   finish:

   free(xp);
   free(vm);
   free(xopt);
   free(c);
   free(fstar);
   free(nacp);
   free(Fvec);
   free(fi->fx);

   return ret;
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

#define HANDLE_IOPTION(name) \
static int handle_##name##_option (char *subsystem, char *optname, char *value, void *clientdata) \
{\
   Isis_Fit_Engine_Type *e; \
   e = (Isis_Fit_Engine_Type *) clientdata; \
   return handle_int_option (subsystem, optname, value, &e->name); \
}

#define HANDLE_DOPTION(name) \
static int handle_##name##_option (char *subsystem, char *optname, char *value, void *clientdata) \
{\
   Isis_Fit_Engine_Type *e; \
   e = (Isis_Fit_Engine_Type *) clientdata; \
   return handle_double_option (subsystem, optname, value, &e->name); \
}

HANDLE_DOPTION(initial_step)
HANDLE_DOPTION(eps)
HANDLE_DOPTION(t)
HANDLE_DOPTION(rt)
HANDLE_IOPTION(ns)
HANDLE_IOPTION(nt)
HANDLE_IOPTION(neps)
HANDLE_IOPTION(maxevl)
HANDLE_IOPTION(iprint)
HANDLE_IOPTION(iseed1)
HANDLE_IOPTION(iseed2)

static Isis_Option_Table_Type Option_Table [] =
{
     {"initial_step", handle_initial_step_option, ISIS_OPT_REQUIRES_VALUE, "0.05", "Initial step (<= 0 means scale using specified min/max)"},
     {"eps", handle_eps_option, ISIS_OPT_REQUIRES_VALUE, "1.e-6", "Error tolerance for termination"},
     {"t", handle_t_option, ISIS_OPT_REQUIRES_VALUE, "1.0", "Initial Temperature"},
     {"rt", handle_rt_option, ISIS_OPT_REQUIRES_VALUE, "0.85", "Temperature reduction factor"},
     {"ns", handle_ns_option, ISIS_OPT_REQUIRES_VALUE, "20", "Number of cycles (suggest 20)"},
     {"nt", handle_nt_option, ISIS_OPT_REQUIRES_VALUE, "-1", "Number of iterations before temperature reduction (suggest max(100,5*n))"},
     {"neps", handle_neps_option, ISIS_OPT_REQUIRES_VALUE, "4", "Number of final function values to decide termination (suggest 4)"},
     {"maxevl", handle_maxevl_option, ISIS_OPT_REQUIRES_VALUE, "1024", "Maximum number of function evaluations"},
     {"iprint", handle_iprint_option, ISIS_OPT_REQUIRES_VALUE, "0", "print control [0|1|2|3]"},
     {"iseed1", handle_iseed1_option, ISIS_OPT_REQUIRES_VALUE, "1", "first RNG seed (0 <= iseed1 <= 31328)"},
     {"iseed2", handle_iseed2_option, ISIS_OPT_REQUIRES_VALUE, "2", "second RNG seed (0 <= iseed2 <= 30081)"},
     {NULL, NULL, 0, NULL, NULL}
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

static void deallocate (Isis_Fit_Engine_Type *e)
{
   ISIS_FREE (e->engine_name);
   ISIS_FREE (e->default_statistic_name);
   ISIS_FREE (e->option_string);
}

ISIS_FIT_ENGINE_METHOD(simann,name,sname)
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

   e->method = simann;
   e->set_options = set_options;
   e->deallocate = deallocate;
   e->set_range_hook = set_range_hook;
   e->range_hook = NULL;
   e->verbose_hook = verbose_hook;
   e->warn_hook = warn_hook;

   e->initial_step = 0.05;
   e->t = 1.0;
   e->rt = 0.85;
   e->eps = 1.e-6;
   e->ns = 20;
   e->nt = -1;
   e->neps = 4;
   e->maxevl = 1024;
   e->iprint = 0;
   e->iseed1 = 1;
   e->iseed2 = 2;

   e->option_string = isis_make_default_option_string ("simann", Option_Table);
   if (e->option_string == NULL)
     {
        deallocate (e);
        return NULL;
     }

   return e;
}
