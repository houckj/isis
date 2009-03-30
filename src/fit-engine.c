/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Author: John E. Davis <davis@space.mit.edu>
            John C. Houck  <houck@space.mit.edu>

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

/* $Id: fit-engine.c,v 1.51 2004/09/07 14:50:03 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <slang.h>

#include "_isis.h"
#include "isis.h"
#include "util.h"
#include "fit.h"
#include "errors.h"

/*}}}*/

extern void list_statistics_and_engines (void);
extern void _add_slang_statistic (char *);

char *Fit_Default_Fit_Method = "lmdif";
char *Fit_Default_Fit_Statistic = "chisqr";
int Isis_Fit_In_Progress = 0;

static Isis_Fit_Engine_Type *Engine_List = NULL;

typedef struct _Statistic_List_Type
{
   struct _Statistic_List_Type *next;
   char *sname;
   Isis_Fit_Statistic_Type *stat;
}
Statistic_List_Type;
static Statistic_List_Type *Statistic_List;

/* Fit engine */

static Isis_Fit_Engine_Type *find_fit_engine (char *name) /*{{{*/
{
   Isis_Fit_Engine_Type *e;

   e = Engine_List;
   while (e != NULL)
     {
        if (0 == strcmp (name, e->engine_name))
          return e;

        e = e->next;
     }
   return NULL;
}

/*}}}*/

Isis_Fit_Engine_Type *isis_find_fit_engine (char *name) /*{{{*/
{
   Isis_Option_Type *opts;
   Isis_Fit_Engine_Type *e;

   if (name == NULL)
     name = Fit_Default_Fit_Method;

   if (NULL == (opts = isis_parse_option_string (name)))
     return NULL;

   name = opts->subsystem;

   if (NULL == (e = find_fit_engine (name)))
     {
        fprintf (stderr, "fit method %s does not exist\n", name);
        isis_free_options (opts);
        return NULL;
     }

   if (e->set_options != NULL)
     {
        if (-1 == e->set_options (e, opts))
          {
             isis_free_options (opts);
             return NULL;
          }
     }

   isis_free_options (opts);

   return e;
}

/*}}}*/

void isis_fit_free_fit_engine (Isis_Fit_Engine_Type *e) /*{{{*/
{
   if (e == NULL)
     return;

   if (e->deallocate != NULL)
     e->deallocate (e);
   ISIS_FREE (e);
}

/*}}}*/

int isis_fit_add_engine (char *name, char *sname, Isis_Fit_Engine_Init_Type *init) /*{{{*/
{
   Isis_Fit_Engine_Type *e;

   if (sname == NULL)
     sname = Fit_Default_Fit_Statistic;

   e = (*init)(name, sname);
   if (e == NULL)
     return -1;

   e->next = Engine_List;
   Engine_List = e;

   return 0;
}

/*}}}*/

int isis_fit_load_fit_routine (char *file, char *name, char *sname) /*{{{*/
{
   Isis_Fit_Engine_Init_Type *f;

   f = (Isis_Fit_Engine_Init_Type *)isis_load_function (file, name, "feng");
   if (f == NULL)
     return -1;

   return isis_fit_add_engine (name, sname, f);
}

/*}}}*/

/* Fit statistic */

static Isis_Fit_Statistic_Fun_Type penalty_statistic;

static int penalty_statistic (Isis_Fit_Statistic_Type *st, /*{{{*/
                              double *y, double *fx, double *w,
                              unsigned int npts, double *vec, double *stat)
{
   SLang_Array_Type *sl_pars = NULL;
   double *pars = NULL;
   double penalty;
   unsigned int num;

   if (-1 == (*st->assigned_fun)(st, y, fx, w, npts, vec, stat))
     return -1;

   if (st->constraint_fun == NULL)
     return 0;

   if (-1 == Fit_copy_fun_params ("constraint", 1, &pars, &num))
     return -1;

   if (num > 0)
     {
        int n = num;
        sl_pars = SLang_create_array (SLANG_DOUBLE_TYPE, 0, pars, &n, 1);
        if (NULL == sl_pars)
          return -1;
     }

   SLang_start_arg_list ();
   /* NULL sl_pars is ok */
   if ((-1 == SLang_push_double (*stat))
       || (-1 == SLang_push_array (sl_pars, 0)))
     {
        SLang_end_arg_list ();
        SLang_free_array (sl_pars);
        return -1;
     }
   SLang_end_arg_list ();

   if ((-1 == SLexecute_function ((SLang_Name_Type *)st->constraint_fun))
       || -1 == SLang_pop_double (&penalty))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "evaluating fit-constraint function");
        SLang_free_array (sl_pars);
        return -1;
     }

   SLang_free_array (sl_pars);

   *stat = penalty;

   return 0;
}

/*}}}*/

static Isis_Fit_Statistic_Type *find_statistic (char *name) /*{{{*/
{
   Statistic_List_Type *s;

   s = Statistic_List;
   while (s != NULL)
     {
        if (0 == strcmp (name, s->sname))
          return s->stat;

        s = s->next;
     }
   return NULL;
}

/*}}}*/

Isis_Fit_Statistic_Type *isis_find_fit_statistic (char *name) /*{{{*/
{
   Isis_Option_Type *opts;
   Isis_Fit_Statistic_Type *s;

   if (name == NULL)
     name = Fit_Default_Fit_Method;

   if (NULL == (opts = isis_parse_option_string (name)))
     return NULL;

   name = opts->subsystem;

   if (NULL == (s = find_statistic (name)))
     {
        fprintf (stderr, "fit statistic %s does not exist\n", name);
        isis_free_options (opts);
        return NULL;
     }

   s->message_string = NULL;

   if (s->set_options != NULL)
     {
        if (-1 == s->set_options (s, opts))
          {
             isis_free_options (opts);
             return NULL;
          }
     }

   isis_free_options (opts);

   return s;
}

/*}}}*/

void list_statistics_and_engines (void) /*{{{*/
{
   Statistic_List_Type *s;
   Isis_Fit_Engine_Type *e;

   fprintf (stdout, "Fit Statistics:\n");
   s = Statistic_List;
   while (s != NULL)
     {
        fprintf (stdout, "\t%s\n", s->sname);
        s = s->next;
     }

   fprintf (stdout, "Minimization Algorithms:\n");
   e = Engine_List;
   while (e != NULL)
     {
        fprintf (stdout, "\t%s\n", e->engine_name);
        e = e->next;
     }

   fflush (stdout);
}

/*}}}*/

void isis_fit_free_fit_statistic (Isis_Fit_Statistic_Type *s) /*{{{*/
{
   if (s == NULL)
     return;

   if (s->deallocate != NULL)
     s->deallocate (s);
   if (s->sl_fun != NULL)
     SLang_free_function ((SLang_Name_Type *)s->sl_fun);
   if (s->sl_report != NULL)
     SLang_free_function ((SLang_Name_Type *)s->sl_report);

   ISIS_FREE (s->symbol);   
   ISIS_FREE (s);
}

/*}}}*/

static void free_statistic_type (Statistic_List_Type *s) /*{{{*/
{
   if (s == NULL)
     return;

   isis_fit_free_fit_statistic (s->stat);
   ISIS_FREE (s->sname);
   ISIS_FREE (s);
}

/*}}}*/

int isis_fit_add_statistic (char *name, Isis_Fit_Statistic_Init_Type *init) /*{{{*/
{
   Statistic_List_Type *s;

   if ((name == NULL)
       || (init == NULL))
     return -1;

   if (NULL == (s = ISIS_MALLOC (sizeof *s)))
     return -1;
   memset ((char *)s, 0, sizeof (*s));

   if (NULL == (s->sname = isis_make_string (name)))
     {
        free_statistic_type (s);
        return -1;
     }

   if (NULL == (s->stat = (*init)()))
     {
        free_statistic_type (s);
        return -1;
     }

   /* Keep a read-only backup copy of the s->fun pointer
    * so we can temporarily over-write it and then restore
    * the original (e.g. to introduce Lagrange multiplier
    * constraints).
    */

   s->stat->assigned_fun = s->stat->fun;

   s->next = Statistic_List;
   Statistic_List = s;

   return 0;
}

/*}}}*/

int isis_fit_load_statistic (char *file, char *sname) /*{{{*/
{
   Isis_Fit_Statistic_Init_Type *x;

   x = (Isis_Fit_Statistic_Init_Type *)isis_load_function (file, sname, "stat");
   if (x == NULL)
     return -1;

   return isis_fit_add_statistic (sname, x);
}

/*}}}*/

static int sl_statistic_function (Isis_Fit_Statistic_Type *s,/*{{{*/
                                  double *y, double *fx, double *w, unsigned int npts,
                                  double *vec, double *stat)
{
   SLang_Array_Type *a_fx, *a_w, *a_y, *sl_vec;
   double st;
   int n, ret = -1;

   *stat = -1.0;

   if (s == NULL || s->sl_fun == NULL)
     return -1;

   sl_vec = NULL;
   a_fx = a_w = a_y = NULL;
   n = (int) npts;
   st = -1.0;

   if ((NULL == (a_y = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
       || (NULL == (a_fx = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
       || ((NULL == (a_w = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))))
     goto free_and_return;

   memcpy (a_y->data, y, npts*sizeof(double));
   memcpy (a_fx->data, fx, npts*sizeof(double));
   memcpy (a_w->data, w, npts*sizeof(double));

   /* (vec, stat) = slang_statistic (y, fx, w) */

   SLang_start_arg_list ();
   if ((-1 == SLang_push_array (a_y, 0))
       || (-1 == SLang_push_array (a_fx, 0))
       || (-1 == SLang_push_array (a_w, 0)))
     goto free_and_return;
   SLang_end_arg_list ();

   if (-1 != SLexecute_function ((SLang_Name_Type *)s->sl_fun))
     {
        (void) SLang_pop_double (&st);

        if ((-1 ==  SLang_pop_array_of_type (&sl_vec, SLANG_DOUBLE_TYPE))
            || (sl_vec == NULL)
            || (sl_vec->num_elements != npts))
          {
             isis_throw_exception (Isis_Error);
          }
        else
          {
             ret = 0;
             memcpy ((char *)vec, (char *)sl_vec->data, npts * sizeof(double));
          }
     }
   /* drop */

   free_and_return:
   SLang_free_array (a_y);
   SLang_free_array (a_fx);
   SLang_free_array (a_w);
   SLang_free_array (sl_vec);

   *stat = st;

   return ret;
}

/*}}}*/

static int sl_report_function (Isis_Fit_Statistic_Type *s, void *pfp, double stat, unsigned int npts, unsigned int nvpars) /*{{{*/
{
   FILE *fp = (FILE *)pfp;
   char *str;

   if (s == NULL || s->sl_report == NULL)
     return -1;

   SLang_start_arg_list ();
   if ((-1 == SLang_push_double (stat))
       || (-1 == SLang_push_integer ((int) npts))
       || (-1 == SLang_push_integer ((int) nvpars)))
     return -1;
   SLang_end_arg_list ();

   if (-1 == SLexecute_function ((SLang_Name_Type *)s->sl_report))
     return -1;

   if (-1 == SLang_pop_slstring (&str))
     return -1;

   if (EOF == fputs (str, fp))
     {
        SLang_free_slstring (str);
        return -1;
     }
   SLang_free_slstring (str);
   return 0;
}

/*}}}*/

static Isis_Fit_Statistic_Type *init_sl_statistic (void) /*{{{*/
{
   SLang_Name_Type *statistic_fun, *report_fun;
   Isis_Fit_Statistic_Type *s;

   if (NULL == (s = ISIS_MALLOC (sizeof(Isis_Fit_Statistic_Type))))
     return NULL;
   memset ((char *)s, 0, sizeof (*s));

   if (NULL == (report_fun = SLang_pop_function ()))
     {
        ISIS_FREE (s);
        return NULL;
     }

   if (NULL == (statistic_fun = SLang_pop_function ()))
     {
        ISIS_FREE (s);
        SLang_free_function (report_fun);
        return NULL;
     }

   if (NULL == (s->symbol = isis_make_string (statistic_fun->name)))
     {
        ISIS_FREE(s);
        SLang_free_function (report_fun);
        SLang_free_function (statistic_fun);
        return NULL;
     }

   s->fun = sl_statistic_function;
   s->report = sl_report_function;
   s->sl_fun = (isis_fptr_type) statistic_fun;
   s->sl_report = (isis_fptr_type) report_fun;

   return s;
}

/*}}}*/

void _add_slang_statistic (char *name) /*{{{*/
{
   if (-1 == isis_fit_add_statistic (name, init_sl_statistic))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "adding fit statistic %s", name);
}

/*}}}*/

/* open and close */

int init_fit_engine (void) /*{{{*/
{
   if (-1 == isis_fit_add_engine ("subplex", "chisqr", Isis_subplex_feng))
     return -1;

   if (-1 == isis_fit_add_engine ("simann", "chisqr", Isis_simann_feng))
     return -1;

   if (-1 == isis_fit_add_engine ("marquardt", "chisqr", Isis_marquardt_feng))
     return -1;

   if (-1 == isis_fit_add_engine ("lmdif", "chisqr", Isis_lmdif_feng))
     return -1;

   if (-1 == isis_fit_add_statistic ("cash", Isis_cash_stat))
     return -1;

   if (-1 == isis_fit_add_statistic ("chisqr", Isis_chisqr_stat))
     return -1;

   return 0;
}

/*}}}*/

void deinit_fit_engine (void) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   Statistic_List_Type *s;

   e = Engine_List;
   while (e != NULL)
     {
        Engine_List = e->next;
        isis_fit_free_fit_engine (e);
        e = Engine_List;
     }

   s = Statistic_List;
   while (s != NULL)
     {
        Statistic_List = s->next;
        free_statistic_type (s);
        s = Statistic_List;
     }
}

/*}}}*/

Isis_Fit_Type *isis_fit_open_fit (char *name, char *sname, Isis_Fit_Fun_Type *fun, /*{{{*/
                                  SLang_Name_Type *constraint_fun)
{
   Isis_Fit_Engine_Type *e;
   Isis_Fit_Statistic_Type *s;
   Isis_Fit_Type *f;

   if (fun == NULL)
     {
        fprintf (stderr, "fit_open_fit: invalid parameters\n");
        return NULL;
     }

   if (NULL == (e = isis_find_fit_engine (name)))
     return NULL;

   if ((sname == NULL)
       || (0 == strcmp (sname, "default")))
     sname = e->default_statistic_name;

   if (NULL == (s = isis_find_fit_statistic (sname)))
     {
        fprintf (stderr, "fit_open_fit: statistic %s does not exist\n", sname);
        return NULL;
     }

   if (NULL == (f = ISIS_MALLOC (sizeof(Isis_Fit_Type))))
     return NULL;

   /* When a constraint function is provided, we call it
    * by introducing an extra indirection.  When the
    * constraint is removed, we restore the original
    * function pointer.
    */

   if (constraint_fun)
     {
        s->constraint_fun = constraint_fun;
        s->fun = &penalty_statistic;
     }
   else
     {
        s->constraint_fun = NULL;
        s->fun = s->assigned_fun;
     }

   f->fun = fun;
   f->engine = e;
   f->stat = s;
   f->statistic = DBL_MAX;

   Isis_Fit_In_Progress = 1;

   return f;
}

/*}}}*/

void isis_fit_close_fit (Isis_Fit_Type *f) /*{{{*/
{
   Isis_Fit_Statistic_Type *s;
   
   Isis_Fit_In_Progress = 0;
   if (f == NULL)
     return;

   s = f->stat;
   if (s->message_string)
     {
        fprintf (stdout, "%s\n", s->message_string);
        s->message_string = NULL;
     }

   ISIS_FREE (f);
}

/*}}}*/

/* perform fits, evaluate statistics */

int isis_invalid_params (Isis_Fit_Engine_Type *e, double *pars, unsigned int npars) /*{{{*/
{
   unsigned int i;

   /* validate returned parameter values */
   for (i = 0; i < npars; i++)
     {
        if ((0 == isfinite (pars[i]))
            || (pars[i] < e->par_min[i])
            || (e->par_max[i] < pars[i]))
          {
             return 1;
          }
     }

   return 0;
}

/*}}}*/

int isis_fit_perform_fit (Isis_Fit_Type *f, void *clientdata, /*{{{*/
                          double *x, double *y, double *weights, unsigned int npts,
                          double *pars, unsigned int npars, double *statistic)
{
   Isis_Fit_Engine_Type *e;
   double *save_pars = NULL;
   unsigned int i;
   int status;

   if (f == NULL)
     return -1;

   e = f->engine;

   if (NULL == (save_pars = ISIS_MALLOC (npars * sizeof(double))))
     return -1;
   memcpy ((char *)save_pars, (char *)pars, npars * sizeof(double));

   status = e->method (f, clientdata, x, y, weights, npts, pars, npars);
   *statistic = f->statistic;

   /* validate returned parameter values */
   if (isis_invalid_params (e, pars, npars))
     {
        memcpy ((char *)pars, (char *)save_pars, npars * sizeof(double));
        *statistic = DBL_MAX;
        status = -1;
     }

   ISIS_FREE(save_pars);

   return status;
}

/*}}}*/

int isis_fit_report_statistic (Isis_Fit_Type *f, FILE *fp, double stat, unsigned int npts, unsigned int nvpars) /*{{{*/
{
   if ((f == NULL) || (f->stat->report == NULL))
     return -1;

   return f->stat->report (f->stat, fp, stat, npts, nvpars);
}

/*}}}*/

int isis_delta_stat_is_chisqr_distrib (Isis_Fit_Type *f) /*{{{*/
{
   if (f == NULL)
     return -1;

   return f->stat->delta_is_chisqr_distributed;
}

/*}}}*/

/* Set fit properties */

int isis_fit_set_ranges (Isis_Fit_Type *f, double *par_min, double *par_max) /*{{{*/
{
   if ((f == NULL) || (par_min == NULL) || (par_max == NULL))
     return -1;
   f->engine->par_min = par_min;
   f->engine->par_max = par_max;
   return 0;
}

/*}}}*/

int isis_fit_set_warn_hook (Isis_Fit_Type *f, Isis_Fit_Warning_Hook_Type *v) /*{{{*/
{
   if (f == NULL)
     return -1;
   f->engine->warn_hook = v;
   return 0;
}

/*}}}*/

int isis_fit_set_verbose_hook (Isis_Fit_Type *f, Isis_Fit_Verbose_Hook_Type *v) /*{{{*/
{
   if (f == NULL)
     return -1;
   f->engine->verbose_hook = v;
   return 0;
}

/*}}}*/

int isis_fit_set_range_hook (Isis_Fit_Type *f, Isis_Fit_Range_Hook_Type *r) /*{{{*/
{
   if ((f == NULL) || (f->engine->set_range_hook == NULL))
     return -1;
   return f->engine->set_range_hook (f->engine, r);
}

/*}}}*/

int isis_fit_set_verbose_level (Isis_Fit_Type *f, int verbose) /*{{{*/
{
   if (f == NULL)
     return -1;
   f->engine->verbose = verbose;
   return 0;
}

/*}}}*/

