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

/* $Id: fit-chisqrconf.c,v 1.28 2004/09/10 02:31:09 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <limits.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>

#include "isis.h"
#include "util.h"
#include "fit.h"
#include "errors.h"
#include "_isis.h"

/*}}}*/

enum
{
   EVAL_ERROR    = -3,
   EVAL_FAILED   = -2,
   EVAL_INVALID  = -1,
   EVAL_OK       =  0,
   EVAL_IMPROVED =  1
};

struct Search_Info_Type
{
   double min_chisqr;
   double delt;
   double tolerance;
   int find_best;
};

struct Sample_Type
{
   double v;
   double chisqr;
   double *par;
};

static int run_fit_improved_hook (void) /*{{{*/
{
   static char hook_name[] = "isis_fit_improved_hook";
   int status = -1;

   if (2 != SLang_is_defined (hook_name))
     return 0;

   SLang_execute_function (hook_name);

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "getting return value from isis_fit_improved_hook");
     }

   return status;
}

/*}}}*/

static int examine_fit_statistic (Param_t *pt, unsigned int idx, double p, double *chisqr,  /*{{{*/
                                  Search_Info_Type *sinfo, Fit_Object_Type *fo)
{
   Fit_Info_Type *info = fo->info;
   /* should have
    *  examine_fit_statistic:improve_tol <= find_limit:accept_tol
    */
   double improve_tol = 0.5 * sinfo->tolerance * sinfo->delt;
   double dchisqr;

   if (info->verbose)
     verbose_warn_hook (NULL, "Minimizing for par[%d]= %15.8e\n", idx, p);

   /* frozen at this value */
   Fit_set_param_value (pt, idx, p);
   if (run_fit_improved_hook ())
     return EVAL_IMPROVED;

   if ((-1 == fit_object_config (fo, pt, 1))
       || (-1 == fit_statistic (fo, sinfo->find_best, chisqr, NULL)))
     return EVAL_FAILED;
   dchisqr = *chisqr - sinfo->min_chisqr;

   if (info->verbose)
     verbose_warn_hook (NULL, "par[%d]= %15.8e  stat= %11.4e  dstat= %11.4e\n",
                        idx, p, *chisqr, dchisqr);

   if ((dchisqr < 0.0) && (-dchisqr > improve_tol))
     {
        verbose_warn_hook (NULL, "Found improved fit, stat= %0.6g for param[%d] = %0.6g\n", *chisqr, idx, p);
        return EVAL_IMPROVED;
     }

   return EVAL_OK;
}

/*}}}*/

enum
{
   USE_QUADRATIC = 0,
   USE_LINEAR    = 1,
   USE_BISECT    = 2,
   MAX_BRACKET_TRIES = 30,
   MAX_CONVERGE_TRIES = 30
};

static int interp_par (int mode, Sample_Type *a, Sample_Type *b, /*{{{*/
                       double min_chisqr, double target_chisqr,
                       double *ptest, Fit_Param_t *par)
{
   double r, factor;
   int i;

   switch (mode)
     {
      case USE_BISECT:
        *ptest = 0.5 * (a->v + b->v);
        for (i = 0; i < par->npars; i++)
          par->par[i] = 0.5 * (a->par[i] + b->par[i]);
        mode = USE_LINEAR;
        break;

      case USE_LINEAR:
        r = (target_chisqr - b->chisqr) / (a->chisqr - b->chisqr);
        *ptest = r * a->v + (1.0 - r) * b->v;
        factor = 1-r;
        for (i = 0; i < par->npars; i++)
          par->par[i] = r * a->par[i] + (1.0 - r) * b->par[i];
        break;

      case USE_QUADRATIC:
        factor = sqrt ((target_chisqr - min_chisqr)/(a->chisqr - min_chisqr));
        *ptest = b->v + (a->v - b->v)*factor;
        for (i = 0; i < par->npars; i++)
          par->par[i] = b->par[i] + (a->par[i] - b->par[i]) * factor;
        mode = USE_LINEAR;
        break;
     }

   return mode;
}

/*}}}*/

static int bracket_target (Param_t *pt, int idx, double *conf_limit, /*{{{*/
                           double pbest, double ptest, double prange,
                           Sample_Type *a, Sample_Type *b,
                           Search_Info_Type *sinfo, Fit_Object_Type *fo)
{
   Fit_Info_Type *info = fo->info;
   Fit_Param_t *par = info->par;
   double chisqr, target_chisqr;
   int i, k, status;

   /* 'b' starts at the minimum chisqr.
    * Move 'b' away from the minimum until the target chisqr is exceeded.
    * After the loop exits, 'a' and 'b' should more-or-less narrowly bracket
    * the target chisqr with 'b' closest to the minimum and 'a' just outside.
    */

   b->v = pbest;
   b->chisqr = sinfo->min_chisqr;
   for (i = 0; i < par->npars; i++)
     b->par[i] = par->par[i];

   status = examine_fit_statistic (pt, idx, ptest, &chisqr, sinfo, fo);
   if (status != EVAL_OK)
     {
        *conf_limit = ptest;
        return status;
     }

   target_chisqr = sinfo->min_chisqr + sinfo->delt;

   k = 0;
   while (chisqr < target_chisqr)
     {
        k++;
        if ((k > MAX_BRACKET_TRIES) || (ptest == prange))
          {
             verbose_warn_hook (NULL, "**** Parameter range endpoint %0.6g is inside the confidence limit\n", prange);
             *conf_limit = prange;
             return EVAL_INVALID;
          }

        if (SLang_get_error())
          return EVAL_INVALID;

        for (i = 0; i < par->npars; i++)
          b->par[i] = par->par[i];

        /* Try to avoid the endpoint if possible */
        ptest = 0.5 * (prange + ptest);
        status = examine_fit_statistic (pt, idx, ptest, &chisqr, sinfo, fo);
        if (status != EVAL_OK)
          {
             *conf_limit = ptest;
             return status;
          }

        if (chisqr >= target_chisqr)
          break;

        b->v = ptest;
        b->chisqr = chisqr;
     }

   a->v = ptest;
   a->chisqr = chisqr;
   for (i = 0; i < par->npars; i++)
     a->par[i] = par->par[i];

   return 0;
}

/*}}}*/

static int find_limit (double *conf_limit, unsigned int idx, /*{{{*/
                       double ptest, double prange, double pbest,
                       Search_Info_Type *sinfo, Fit_Object_Type *fo,
                       Param_t *pt, char *which)
{
   Fit_Info_Type *info = fo->info;
   double test, chisqr, target_chisqr, accept_tol;
   int i, k, count, status, verbose;
   Sample_Type a, b;
   Fit_Param_t *par;
   int mode = USE_QUADRATIC;

   /* on input:
    *  pbest = best fit parameter value (corresponding to min_chisqr),
    *  ptest = parameter value to test first
    *  prange = extreme allowed parameter value (biggest or smallest;
    *           pbest provides the other limit)
    *  idx = index of (frozen) param
    *  par -> best fit values of all variable parameters
    */

   *conf_limit = prange;

   if (info->par == NULL)
     {
        status = EVAL_ERROR;
        goto finish;
     }

   par = info->par;
   verbose = info->verbose;

   sinfo->find_best = (par->npars > 0) ? 1 : 0;

   a.par = ISIS_MALLOC (par->npars * sizeof(double));
   b.par = ISIS_MALLOC (par->npars * sizeof(double));
   if ((a.par == NULL) || (b.par == NULL))
     {
        status = EVAL_ERROR;
        goto finish;
     }

   status = bracket_target (pt, idx, conf_limit, pbest, ptest, prange, &a, &b, sinfo, fo);
   if (status != EVAL_OK)
     goto finish;

   /* At this point, 'b' has the last value tested < target_chisqr
    * and 'a' has the first value tested > target_chisqr
    * With target_chisqr bracketed, interpolation should work well:
    *
    * Note that should have
    *  examine_fit_statistic:improve_tol <= find_limit:accept_tol
    *
    * Also, while testing convergence of the fit-parameter would
    * be nice, its difficult to do since we don't know anything
    * about the scale size of the parameter.  Parameters for which
    * zero is a valid value are particularly troublesome.
    */

   target_chisqr = sinfo->min_chisqr + sinfo->delt;
   accept_tol = sinfo->tolerance * sinfo->delt;

   count = 0;
   k = 0;
   do
     {
        Sample_Type *update = &a;

        if (fabs(a.chisqr - b.chisqr) < accept_tol)
          {
             verbose_warn_hook (NULL, "*** fit statistic isnt changing...\n");
             *conf_limit = ptest;
             status = EVAL_INVALID;
             goto finish;
          }

        if (SLang_get_error())
          {
             *conf_limit = ptest;
             status = EVAL_INVALID;
             goto finish;
          }

        mode = interp_par (mode, &a, &b, sinfo->min_chisqr, target_chisqr, &ptest, par);
        status = examine_fit_statistic (pt, idx, ptest, &chisqr, sinfo, fo);
        if (status != EVAL_OK)
          {
             *conf_limit = ptest;
             goto finish;
          }

        test = chisqr - target_chisqr;

        if (test > 0.0)
          {
             count++;
             update = &a;
          }
        else if (test < 0.0)
          {
             count--;
             update = &b;
          }

        update->v = ptest;
        update->chisqr = chisqr;
        for (i = 0; i < par->npars; i++)
          update->par[i] = par->par[i];

        if ((abs (count) > 2) || (fabs(test) > 3))
          {
             count = 0;
             mode = USE_BISECT;
          }

     } while ((fabs(test) > accept_tol) && (k++ < MAX_CONVERGE_TRIES));

   *conf_limit = ptest;

   if (fabs(test) > accept_tol)
     {
        status = EVAL_INVALID;
        goto finish;
     }

   if (verbose > 0)
     verbose_warn_hook (NULL,"limit found: %15.8e  dstat= %11.4e\n",
                        ptest, chisqr - sinfo->min_chisqr);

   status = EVAL_OK;

   finish:

   ISIS_FREE (a.par);
   ISIS_FREE (b.par);

   switch (status)
     {
      case EVAL_OK:
        break;
      case EVAL_IMPROVED:
        verbose_warn_hook (NULL, "**** Found improved fit\n");
        Fit_unpack_variable_params (pt, info->par->par);
        Fit_set_param_value (pt, idx, *conf_limit);
        break;
      case EVAL_FAILED:
        verbose_warn_hook (NULL, "**** %s confidence limit: fit failed\n", which);
        break;
      case EVAL_INVALID:
        verbose_warn_hook (NULL, "**** %s confidence limit didn't converge[%d]:  allow wider parameter ranges?\n", which, idx);
        break;
      case EVAL_ERROR:
        verbose_warn_hook (NULL, "**** internal error\n");
        break;
      default:
        isis_vmesg (FAIL, I_INTERNAL, __FILE__, __LINE__, "find_limit: invalid status for parameter %d", idx);
        break;
     }

   return status;
}

/*}}}*/

static Param_Info_t *get_valid_param_info (Param_t *pt, int idx) /*{{{*/
{
   Param_Info_t *p;

   if (NULL == (p = Fit_param_info (pt, idx)))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "invalid param index = %d", idx);
        return NULL;
     }

   if (p->tie_param_name)
     {
        Param_Info_t *tie_info;
        tie_info = Fit_find_param_info_by_full_name (pt, p->tie_param_name);
        if (tie_info != NULL)
          isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "parameter %d is tied to parameter %d",
                      idx, tie_info->idx);
        else
          isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__,
                      "parameter %d is tied to a nonexistent parameter(!)", idx);
        return NULL;
     }

   if ((p->min == -DBL_MAX) || (p->max == DBL_MAX))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "param %d:  please set min/max values for search.", idx);
        return NULL;
     }

   return p;
}

/*}}}*/

static int init_conf_limit_search (Param_t *pt, int idx, /*{{{*/
                                   Param_Info_t **par_info, Fit_Param_t **pars)
{
   unsigned int num_all, num_vary;

   if (NULL == (*par_info = get_valid_param_info (pt, idx)))
     return -1;

   if (-1 == Fit_count_params (pt, &num_all, &num_vary))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (NULL == (*pars = new_fit_param_type (num_all)))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return 0;
}

/*}}}*/

int get_confidence_limits (Fit_Object_Type *fo, Param_t *pt, Isis_Fit_CLC_Type *ctrl,  /*{{{*/
                           int idx, double *pconf_min, double *pconf_max)
{
   Fit_Info_Type *info = NULL;
   Search_Info_Type sinfo;
   Param_Info_t *par_info = NULL;
   Fit_Param_t *initial_pars = NULL;
   double pstart, pbest, conf_min, conf_max;
   unsigned int initial_state;
   int ret = 0;

   if ((fo == NULL) || (ctrl == NULL) || (pconf_min == NULL) || (pconf_max == NULL))
     return -1;

   *pconf_min = conf_min = 0.0;
   *pconf_max = conf_max = 0.0;
   initial_state = 0;

   sinfo.delt = ctrl->delta_stat;
   sinfo.tolerance = ctrl->tol;

   init_verbose_hook ();

   if (-1 == fit_object_config (fo, pt, 0))
     goto free_and_return;
   info = fo->info;

   if (-1 == init_conf_limit_search (pt, idx, &par_info, &initial_pars))
     return -1;

   /* Unless the confidence limit search turns up a better fit,
    * we'll restore the initial state of the parameter table
    * when the search ends.
    * Its best to use the given initial values as the reference point,
    * rather than those returned by a search for the best fit.
    * Compute the current chi^2 = the "best" fit (by assumption)
    * then freeze param for confidence limit search
    */

   /* Save the initial state */
   initial_state = par_info->freeze;
   if (-1 == Fit_pack_all_params (pt, initial_pars))
     goto free_and_return;

   /* First freeze param of interest, then pack remaining variable params */
   Fit_set_freeze (pt, idx, 1);
   Fit_get_param_value (pt, idx, &pbest);

   if (-1 == fit_object_config (fo, pt, 1))
     goto restore_old_params;

   /* Evaluate the current fit-statistic */
   if (-1 == fit_statistic (fo, 0, &sinfo.min_chisqr, NULL))
     goto restore_old_params;

   if (ctrl->verbose > 0)
     verbose_warn_hook (NULL, "Best fit par[%d]= %15.8e: stat=%15.8e\n",
                        idx, pbest, sinfo.min_chisqr);

   /* Search for the lower limit */
   pstart = 0.5 * (pbest + par_info->min);

   ret = find_limit (&conf_min, idx, pstart, par_info->min, pbest, &sinfo, fo, pt, "Lower");
   if (ret == EVAL_ERROR || ret == EVAL_FAILED)
     goto restore_old_params;
   else if (ret == EVAL_IMPROVED)
     {
        conf_max = conf_min;
        goto immediate_return;
     }

   /* Got the lower one -- reset to the initial params */
   Fit_unpack_all_params (pt, initial_pars->par);

   if (-1 == fit_object_config (fo, pt, 1))
     goto restore_old_params;

   /* Search for the upper limit */
   pstart = pbest + 2.0*(pbest - conf_min);
   if (pstart >= par_info->max)
     pstart = 0.5 * (pbest + par_info->max);

   ret = find_limit (&conf_max, idx, pstart, par_info->max, pbest, &sinfo, fo, pt, "Upper");
   if (ret == EVAL_ERROR || ret == EVAL_FAILED)
     goto restore_old_params;
   else if (ret == EVAL_IMPROVED)
     {
        conf_min = conf_max;
        goto immediate_return;
     }

   /* For better or worse, we're done except for cleanup */

   restore_old_params:
   Fit_unpack_all_params (pt, initial_pars->par);

   immediate_return:
   if (par_info) Fit_set_freeze (pt, idx, initial_state);

   /* Evaluate the model once more to leave it consistent
    * with the restored parameter values */
   if (0 == fit_object_config (fo, pt, 1))
     {
        double not_used;
        (void) fit_statistic (fo, 0, &not_used, NULL);
     }

   free_and_return:
   deinit_verbose_hook ();
   free_fit_param_type (initial_pars);

   *pconf_min = conf_min;
   *pconf_max = conf_max;

   return 0;
}

/*}}}*/

/*}}}*/

