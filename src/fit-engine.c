/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2018 Massachusetts Institute of Technology

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

#include "isis.h"
#include "_isis.h"
#include "util.h"
#include "fit.h"
#include "errors.h"

/*}}}*/

extern void list_statistics_and_engines (void);
extern void _add_slang_statistic (char *);

char *Fit_Default_Fit_Method = "mpfit";
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
   double old_stat, new_stat, penalty, b, vec_penalty;
   SLindex_Type num;
   unsigned int i, n;

   if (-1 == (*st->assigned_fun)(st, y, fx, w, npts, vec, stat))
     return -1;

   if (st->constraint_fun == NULL)
     return 0;

   if (-1 == Fit_copy_fun_params ("constraint", 1, &pars, &n))
     return -1;

   num = n;

   if (num > 0)
     {
        sl_pars = SLang_create_array (SLANG_DOUBLE_TYPE, 0, pars, &num, 1);
        if (NULL == sl_pars)
          return -1;
     }

   old_stat = *stat;

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
       || -1 == SLang_pop_double (&new_stat))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "evaluating fit-constraint function");
        SLang_free_array (sl_pars);
        return -1;
     }

   SLang_free_array (sl_pars);

   *stat = new_stat;

   /* The penalty must also affect the vector statistic somehow.
    * Try spreading it uniformly over all bins, assuming that
    * the base statistic is the Euclidean norm = sum (vec^2);
    * We want to define:
    *   new_vec = old_vec + vec_penalty
    * so that
    *   new_stat = sum(new_vec^2) = old_stat + penalty.
    *
    * FIXME?  This seems ok for chi-square, but perhaps something
    * different would be better for max. likelihood or cash statistic?
    * Maybe the statistic object should have a vec_penalty() method?
    */
   b = 2 * isis_kahan_sum (vec, npts) / npts;
   penalty = new_stat - old_stat;
   vec_penalty = -0.5*b + sqrt (0.25*b*b + penalty / npts);
   for (i = 0; i < npts; i++)
     {
        vec[i] += vec_penalty;
     }

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

static void free_statistic_opt_data (Isis_Fit_Statistic_Optional_Data_Type *opt_data) /*{{{*/
{
   if (opt_data == NULL)
     return;
   ISIS_FREE(opt_data->bkg);
   ISIS_FREE(opt_data->bkg_at);
   ISIS_FREE(opt_data->src_at);
   ISIS_FREE(opt_data);
}

/*}}}*/

static Isis_Fit_Statistic_Optional_Data_Type *allocate_statistic_opt_data (int n) /*{{{*/
{
   Isis_Fit_Statistic_Optional_Data_Type *opt_data;

   if (NULL == (opt_data = (Isis_Fit_Statistic_Optional_Data_Type *)ISIS_MALLOC (sizeof *opt_data)))
     return NULL;
   memset ((char *)opt_data, 0, sizeof *opt_data);

   opt_data->num = n;
   opt_data->malloced = 1;

   if ((NULL == (opt_data->bkg = (double *)ISIS_MALLOC (n * sizeof(double))))
       ||(NULL == (opt_data->bkg_at = (double *)ISIS_MALLOC (n * sizeof(double))))
       ||(NULL == (opt_data->src_at = (double *)ISIS_MALLOC (n * sizeof(double))))
       )
     {
        free_statistic_opt_data (opt_data);
        return NULL;
     }
   memset ((char *)opt_data->bkg, 0, n*sizeof(double));
   memset ((char *)opt_data->bkg_at, 0, n*sizeof(double));
   memset ((char *)opt_data->src_at, 0, n*sizeof(double));

   return opt_data;
}

/*}}}*/

void isis_fit_free_fit_statistic (Isis_Fit_Statistic_Type *s) /*{{{*/
{
   if (s == NULL)
     return;

   if (s->deallocate != NULL)
     s->deallocate (s);
   if (s->sl_fun != NULL)
     SLang_free_function (s->sl_fun);
   if (s->sl_report != NULL)
     SLang_free_function (s->sl_report);

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

   if (NULL == (s = (Statistic_List_Type *) ISIS_MALLOC (sizeof *s)))
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

   s->stat->assigned_fun = s->stat->compute_statistic;

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

typedef struct
{
   SLang_Array_Type *bkg;
   SLang_Array_Type *bkg_at;
   SLang_Array_Type *src_at;
}
Optional_Data_Type;

static SLang_CStruct_Field_Type Optional_Data_Type_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Optional_Data_Type, bkg, "bkg", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Optional_Data_Type, bkg_at, "bkg_at", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Optional_Data_Type, src_at, "src_at", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int push_opt_data (Isis_Fit_Statistic_Optional_Data_Type *opt_data) /*{{{*/
{
   Optional_Data_Type odt;
   int n, status=-1;

   if (opt_data == NULL)
     {
        SLang_push_null ();
        return 0;
     }

   memset ((char *)&odt, 0, sizeof odt);

   n = opt_data->num;

   if ((NULL == (odt.bkg = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
       || (NULL == (odt.bkg_at = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
       || ((NULL == (odt.src_at = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))))
     goto free_and_return;

   memcpy ((char *)odt.bkg->data, (char *)opt_data->bkg, n*sizeof(double));
   memcpy ((char *)odt.bkg_at->data, (char *)opt_data->bkg_at, n*sizeof(double));
   memcpy ((char *)odt.src_at->data, (char *)opt_data->src_at, n*sizeof(double));

   if (-1 == SLang_push_cstruct ((VOID_STAR)&odt, Optional_Data_Type_Layout))
     goto free_and_return;

   status = 0;
free_and_return:
   SLang_free_array (odt.bkg);
   SLang_free_array (odt.bkg_at);
   SLang_free_array (odt.src_at);

   return status;
}

/*}}}*/

static int sl_statistic_function (Isis_Fit_Statistic_Type *s,/*{{{*/
                                  double *y, double *fx, double *w, unsigned int npts,
                                  double *vec, double *stat)
{
   SLang_Array_Type *a_fx, *a_w, *a_y, *sl_vec;
   SLindex_Type n;
   double st;
   int ret = -1;

   *stat = -1.0;

   if (s == NULL || s->sl_fun == NULL)
     return -1;

   sl_vec = NULL;
   a_fx = a_w = a_y = NULL;
   n = npts;
   st = -1.0;

   if ((NULL == (a_y = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
       || (NULL == (a_fx = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
       || ((NULL == (a_w = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))))
     goto free_and_return;

   memcpy (a_y->data, y, npts*sizeof(double));
   memcpy (a_fx->data, fx, npts*sizeof(double));
   memcpy (a_w->data, w, npts*sizeof(double));

   /* (vec, stat) = slang_statistic (y, fx, w)
    *   OR, if opt_data is used:
    * (vec, stat) = slang_statistic (y, fx, w, opt_data)
    */

   SLang_start_arg_list ();
   if ((-1 == SLang_push_array (a_y, 0))
       || (-1 == SLang_push_array (a_fx, 0))
       || (-1 == SLang_push_array (a_w, 0))
       || ((s->uses_opt_data != 0) && (-1 == push_opt_data (s->opt_data))))
     goto free_and_return;
   SLang_end_arg_list ();

   if (-1 != SLexecute_function (s->sl_fun))
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

static void sl_deallocate_function (Isis_Fit_Statistic_Type *s)
{
   ISIS_FREE (s->symbol);
   ISIS_FREE (s->option_string);
}

static Isis_Fit_Statistic_Type *init_sl_statistic (void) /*{{{*/
{
   SLang_Name_Type *statistic_fun, *report_fun;
   Isis_Fit_Statistic_Type *s;

   if (NULL == (s = (Isis_Fit_Statistic_Type *) ISIS_MALLOC (sizeof(Isis_Fit_Statistic_Type))))
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

   s->compute_statistic = sl_statistic_function;
   s->deallocate = sl_deallocate_function;
   s->report = sl_report_function;
   s->sl_fun = statistic_fun;
   s->sl_report = report_fun;

   return s;
}

/*}}}*/

static int fixup_sl_statistic_name (char *name)
{
   Isis_Fit_Statistic_Type *s;

   if (NULL == (s = find_statistic (name)))
     return -1;

   if (NULL == (s->option_string = isis_make_string (name)))
     return -1;

   ISIS_FREE(s->symbol);
   if (NULL == (s->symbol = isis_make_string (name)))
     return -1;

   return 0;
}

void _add_slang_statistic (char *name) /*{{{*/
{
   if ((-1 == isis_fit_add_statistic (name, init_sl_statistic))
       || (-1 == fixup_sl_statistic_name (name)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "adding fit statistic %s", name);
     }
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

   if (-1 == isis_fit_add_engine ("mpfit", "chisqr", Isis_mpfit_feng))
     return -1;

   if (-1 == isis_fit_add_engine ("lmdif", "chisqr", Isis_lmdif_feng))
     return -1;

   if (-1 == isis_fit_add_statistic ("ml", Isis_ml_stat))
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
                                  SLang_Name_Type *constraint_fun, int n)
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

   if (NULL == (f = (Isis_Fit_Type *) ISIS_MALLOC (sizeof(Isis_Fit_Type))))
     return NULL;

   if (s->uses_opt_data)
     {
        if (NULL == (s->opt_data = allocate_statistic_opt_data (n)))
          {
             ISIS_FREE(f);
             return NULL;
          }
     }

   /* When a constraint function is provided, we call it
    * by introducing an extra indirection.  When the
    * constraint is removed, we restore the original
    * function pointer.
    */

   if (constraint_fun)
     {
        s->constraint_fun = constraint_fun;
        s->compute_statistic = &penalty_statistic;
     }
   else
     {
        s->constraint_fun = NULL;
        s->compute_statistic = s->assigned_fun;
     }

   f->compute_model = fun;
   f->engine = e;
   f->stat = s;
   f->statistic = DBL_MAX;
   f->covariance_matrix = NULL;

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

   if (s->opt_data)
     {
        free_statistic_opt_data (s->opt_data);
        s->opt_data = NULL;
     }

   ISIS_FREE (f->covariance_matrix);
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
   int status;

   if (f == NULL)
     return -1;

   e = f->engine;

   if (NULL == (save_pars = (double *) ISIS_MALLOC (npars * sizeof(double))))
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

   if (Isis_Verbose >= 0)
     {
        (void) SLang_run_hooks ("_isis->check_optimizer_context", 0);
     }

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

int isis_fit_set_param_step (Isis_Fit_Type *f, double *step, double *relstep) /*{{{*/
{
   if ((f == NULL) || (step == NULL) || (relstep == NULL))
     return -1;
   f->engine->par_step = step;
   f->engine->par_relstep = relstep;
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

int isis_fit_set_range_hook (Isis_Fit_Type *f, Isis_Fit_Range_Hook_Type *r, void *client_data) /*{{{*/
{
   if ((f == NULL) || (f->engine->set_range_hook == NULL))
     return -1;
   f->engine->range_hook_client_data = client_data;
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

