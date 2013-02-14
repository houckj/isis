/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2012 Massachusetts Institute of Technology

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

/* $Id: fit-params.c,v 1.42 2004/09/07 14:50:04 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>
#include "isis.h"
#include "util.h"
#include "fit.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

/* parameter table definition */

enum
{
   VARIABLE_PARS = 0,
   ALL_PARS = 1
};

struct Param_t
{
   Param_t *next_fun;
   Param_Info_t *info;
   unsigned int num_params;
   unsigned int fun_id;
   unsigned int fun_type;
   unsigned int fun_version;
};

int Isis_Evaluating_Derived_Param;

static int update_derived_params (Param_t *, int out_of_range_severity);

static void free_param_function (Param_Info_t *p)
{
   if (p == NULL)
     return;
   SLang_free_function (p->fun_ptr);
   ISIS_FREE (p->fun_str);
   p->fun_str = NULL;
}

static void free_param_info (Param_Info_t *i, unsigned int num_params) /*{{{*/
{
   Param_Info_t *pi;

   if (i == NULL)
     return;

   pi = i;
   while (num_params-- > 0)
     {
        free_param_function (pi);
        ISIS_FREE (pi->tie_param_name);
        ISIS_FREE (pi->param_name);
        pi++;
     }

   ISIS_FREE (i);
}

/*}}}*/

static Param_Info_t *new_param_info (unsigned int num_params) /*{{{*/
{
   Param_Info_t *i = (Param_Info_t *) ISIS_MALLOC (num_params * sizeof(Param_Info_t));

   if (i != NULL)
     memset ((char *)i, 0, num_params * sizeof (*i));

   return i;
}

/*}}}*/

void Fit_free_param_table (Param_t *pt) /*{{{*/
{
   while (pt)
     {
        Param_t *tmp = pt->next_fun;
        free_param_info (pt->info, pt->num_params);
        ISIS_FREE (pt);
        pt = tmp;
     }
}

/*}}}*/

Param_t *Fit_new_param_table (unsigned int num_params,  /*{{{*/
                              unsigned int fun_type, unsigned int fun_id, unsigned int fun_version)
{
   Param_t *p;

   if (NULL == (p = (Param_t *) ISIS_MALLOC (sizeof(Param_t))))
     return NULL;
   memset ((char *)p, 0, sizeof (*p));

   p->next_fun = NULL;
   p->num_params = num_params;
   p->fun_type = fun_type;
   p->fun_id = fun_id;
   p->fun_version = fun_version;

   if (num_params == 0)
     p->info = NULL;
   else
     {
        p->info = new_param_info (num_params);
        if (p->info == NULL)
          {
             free_param_info (p->info, p->num_params);
             Fit_free_param_table (p);
             p = NULL;
          }
     }

   return p;
}

/*}}}*/

/* build and search parameter table */

static Param_t *locate_fun_params (Param_t *pt, unsigned int fun_type, unsigned int fun_id) /*{{{*/
{
   while (pt)
     {
        if ((pt->fun_type == fun_type) && (pt->fun_id == fun_id))
          return pt;
        pt = pt->next_fun;
     }

   return NULL;
}

/*}}}*/

static int map_table (Param_t *pt, int which, int (*fun)(Param_Info_t *, void *), void *cl) /*{{{*/
{
   while (pt)
     {
        unsigned int i;
        for (i = 0; i < pt->num_params; i++)
          {
             Param_Info_t *p = &pt->info[i];
             if ((which == ALL_PARS)
                 || ((p->tie_param_name == NULL) && (p->freeze == 0)))
               {
                  if ((*fun)(p, cl))
                    return -1;
               }
          }
        pt = pt->next_fun;
     }

   return 0;
}

/*}}}*/

static Param_Info_t *find_param_by_index (Param_t *pt, unsigned int idx) /*{{{*/
{
   while (pt)
     {
        unsigned int i;
        for (i = 0; i < pt->num_params; i++)
          {
             Param_Info_t *p = &pt->info[i];
             if (p->in_use && (p->idx == idx))
               return p;
          }
        pt = pt->next_fun;
     }

   return NULL;
}

/*}}}*/

static Param_Info_t *find_param_by_vary_index (Param_t *pt, unsigned int vary_idx) /*{{{*/
{
   while (pt)
     {
        unsigned int i;
        for (i = 0; i < pt->num_params; i++)
          {
             Param_Info_t *p = &pt->info[i];
             if (p->in_use
                 && (p->vary_idx == vary_idx)
                 && (p->freeze == 0)
                 && (p->tie_param_name == NULL))
               return p;
          }
        pt = pt->next_fun;
     }

   return NULL;
}

/*}}}*/

Param_Info_t *Fit_param_info2 (Param_t *pt, unsigned int fun_type, unsigned int fun_id, unsigned int fun_par) /*{{{*/
{
   while (pt)
     {
        if ((pt->fun_type == fun_type) && (pt->fun_id == fun_id))
          {
             unsigned int i;
             for (i = 0; i < pt->num_params; i++)
               {
                  Param_Info_t *p = &pt->info[i];
                  if (p->in_use && (p->fun_par == fun_par))
                    return p;
               }
          }
        pt = pt->next_fun;
     }

   return NULL;
}

/*}}}*/

Param_Info_t *Fit_find_param_info_by_full_name (Param_t *pt, char *name) /*{{{*/
{
   Fit_Fun_t *ff;
   char fun_name[MAX_NAME_SIZE], par_name[MAX_NAME_SIZE];
   int n, fun_type, fun_id, fun_par;

   if (name == NULL)
     return NULL;

   n = sscanf (name, "%[^(<]%*1[(<]%d%*1[)>].%s", fun_name, &fun_id, par_name);
   if (n != 3)
     return NULL;

   if (-1 == (fun_type = Fit_get_fun_type (fun_name)))
     return NULL;

   if (NULL == (ff = Fit_get_fit_fun (fun_type)))
     return NULL;

   if (-1 == (fun_par = Fit_get_fun_par (ff, par_name)))
     return NULL;

   return Fit_param_info2 (pt, fun_type, fun_id, fun_par);
}

/*}}}*/

int Fit_get_fun_params (Param_t *pt, unsigned int fun_type, unsigned int fun_id, /*{{{*/
                        double *par)
{
   unsigned int i;

   /* update derived params here to support dataset dependent
    * parameter values. */
   if (-1 == update_derived_params (pt, INTR))
     return -1;

   if (NULL == (pt = locate_fun_params (pt, fun_type, fun_id)))
     return -1;

   for (i = 0; i < pt->num_params; i++)
     {
        par[i] = pt->info[i].value;
     }

   return 0;
}

/*}}}*/

static int append_fun (Param_t *pt, Param_t *fun) /*{{{*/
{
   if (pt == NULL) return -1;

   while (pt->next_fun != NULL)
     pt = pt->next_fun;

   pt->next_fun = fun;
   return 0;
}

/*}}}*/

static int set_param_minmax (Param_Info_t *p, double p_min, double p_max) /*{{{*/
{
   if (p == NULL)
     return -1;

   if ((p_min < p->hard_min) || (p_min > p->hard_max))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "parameter %d minimum=%g violates hard limits (%g, %g)",
                    p->idx, p_min, p->hard_min, p->hard_max);
        return -1;
     }
   if ((p_max < p->hard_min) || (p_max > p->hard_max))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "parameter %d maximum=%g violates hard limits (%g, %g)",
                    p->idx, p_max, p->hard_min, p->hard_max);
        return -1;
     }

   p->min = p_min;
   p->max = p_max;

   if ((p->min == 0.0) && (p->max == 0.0))
     {
        p->min = (p->hard_min > -DBL_MAX) ? p->hard_min : -DBL_MAX;
        p->max = (p->hard_max <  DBL_MAX) ? p->hard_max :  DBL_MAX;
     }

   return 0;
}

/*}}}*/

static int validate_param_minmax (Param_Info_t *p) /*{{{*/
{
   if (p == NULL)
     return -1;

   if ((p->value < p->min) || (p->max < p->value))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid range for parameter %d: value=%g is outside [%g,%g)",
                    p->idx, p->value, p->min, p->max);
        return -1;
     }

   return 0;
}

/*}}}*/

static int init_parameter_of_type (Fit_Fun_t *ff, Param_Info_t *p) /*{{{*/
{
   char name[OUT_PARAM_NAME_SIZE];
   char ch, *s;

   if (p->default_applied)
     return 0;

   s = ff->name[p->fun_par+1];
   ch = (*s != 0) ? '.' : ' ';
   sprintf (name, "%s(%d)%c%s", ff->name[0], p->fun_id, ch, s);

   if (NULL == (p->param_name = isis_make_string (name)))
     return -1;

   if (-1 == (*ff->set_param_default) (ff, p))
     return -1;

   if ((-1 == set_param_minmax (p, p->min, p->max))
       ||(-1 == validate_param_minmax (p)))
     return -1;

   p->default_applied = 1;

   return 0;
}
/*}}}*/

static int param_is_a_norm (Param_Info_t *p, Isis_User_Source_t *s) /*{{{*/
{
   unsigned int k;

   if ((s->norm_indexes == NULL) || (s->num_norms == 0))
     return 0;

   for (k = 0; k < s->num_norms; k++)
     {
        if (p->fun_par == s->norm_indexes[k])
          return 1;
     }

   return 0;
}

/*}}}*/

static int mark_unused (Param_Info_t *p, void *cl) /*{{{*/
{
   (void) cl;
   p->in_use = 0;
   return 0;
}

/*}}}*/

int Fit_mark_params_unused (Param_t *pt) /*{{{*/
{
   return map_table (pt, ALL_PARS, mark_unused, NULL);
}

/*}}}*/

static int delete_fun_params (Param_t *pt, unsigned int fun_type, unsigned int fun_id) /*{{{*/
{
   while (pt)
     {
        Param_t *tmp = pt->next_fun;
        if (tmp == NULL) return 0;
        if ((tmp->fun_type == fun_type) && (tmp->fun_id == fun_id))
          {
             pt->next_fun = tmp->next_fun;
             free_param_info (tmp->info, tmp->num_params);
             ISIS_FREE (tmp);
             return 0;
          }
        pt = pt->next_fun;
     }

   return -1;
}

/*}}}*/

int Fit_register_fun (Param_t *pt, Fit_Fun_t *ff, unsigned int fun_id,  /*{{{*/
                      unsigned int *idx)
{
   Param_t *fun = NULL;
   Param_Info_t *p = NULL;
   unsigned int i = 0;

   /* no duplicates */
   fun = locate_fun_params (pt, ff->fun_type, fun_id);
   if (fun)
     {
        if (fun->fun_version == ff->fun_version)
          {
             for (i = 0; i < fun->num_params; i++)
               {
                  p = &fun->info[i];
                  p->in_use++;
                  if (p->in_use == 1)
                    p->idx = ++(*idx);
               }

             return 0;
          }
        delete_fun_params (pt, ff->fun_type, fun_id);
     }

   if (NULL == (fun = Fit_new_param_table (ff->nparams, ff->fun_type, fun_id, ff->fun_version)))
     return -1;

   for (i = 0; i < ff->nparams; i++)
     {
        p = &fun->info[i];
        p->fun_type = ff->fun_type;
        p->fun_version = ff->fun_version;
        p->fun_id = fun_id;
        p->fun_par = i;
        p->idx = ++(*idx);
        p->in_use = 1;
        p->tie_param_name = NULL;
        p->freeze = 0;
        p->is_a_norm = param_is_a_norm (p, &ff->s);

        p->hard_min = -Isis_Inf;
        p->hard_max =  Isis_Inf;
        p->relstep = Isis_Default_Relstep;

        if (-1 == init_parameter_of_type (ff, p))
          {
             Fit_free_param_table (fun);
             return -1;
          }
     }

   if (-1 == append_fun (pt, fun))
     {
        Fit_free_param_table (fun);
        return -1;
     }

   return 0;
}
/*}}}*/

/* parameter table access */

static int update_tie (Param_Info_t *p, void *cl) /*{{{*/
{
   Param_t *pt = (Param_t *)cl;

   if (p->in_use && p->tie_param_name)
     {
        Param_Info_t *tie_info;

        tie_info = Fit_find_param_info_by_full_name (pt, p->tie_param_name);
        if ((NULL == tie_info)
            || (tie_info->in_use == 0))
          {
             isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "parameter %d tie is invalid",
                         p->idx);
             return -1;
          }

        p->value = tie_info->value;
        return set_param_minmax (p, tie_info->min, tie_info->max);
     }

   return 0;
}

/*}}}*/

static int update_tied_params (Param_t *pt) /*{{{*/
{
   return map_table (pt, ALL_PARS, update_tie, pt);
}

/*}}}*/

static int update_derived (Param_Info_t *p, void *cl) /*{{{*/
{
   int out_of_range_severity = cl ? *(int *)cl : FAIL;
   double value;
   int err = 0;

   if (Isis_Evaluating_Derived_Param)
     return 0;

   if (p->in_use && p->fun_str != NULL)
     {
        Isis_Evaluating_Derived_Param = 1;

        if ((-1 == SLexecute_function (p->fun_ptr))
            ||(-1 == SLang_pop_double (&value)))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "computing derived parameter: %s",
                         p->fun_str);
             err = -1;
          }

        Isis_Evaluating_Derived_Param = 0;

        /* Out of range values are unacceptable. */
        if ((value < p->min || p->max < value)
           || (0 != isnan(value)))
          {
             isis_vmesg (out_of_range_severity, I_ERROR, __FILE__, __LINE__,
                         "parameter %s:  derived value=%g ignored (violates limits min=%g max=%g)",
                         p->param_name, value, p->min, p->max);
             if (out_of_range_severity == INTR)
               err = -1;
          }
        else p->value = value;
     }

   return err;
}

/*}}}*/

static int update_derived_params (Param_t *pt, int out_of_range_severity) /*{{{*/
{
   return map_table (pt, ALL_PARS, update_derived, (void *)&out_of_range_severity);
}

/*}}}*/

int Fit_set_fun_params (Param_t *pt, unsigned int fun_type, unsigned int fun_id, /*{{{*/
                        double *par, double *par_min, double *par_max)
{
   Param_t *fun = locate_fun_params (pt, fun_type, fun_id);
   unsigned int i;

   if (fun == NULL)
     return -1;

   for (i = 0; i < fun->num_params; i++)
     {
        Param_Info_t *p = &fun->info[i];
        p->value = par[i];
        if (-1 == set_param_minmax (p, par_min[i], par_max[i]))
          return -1;
     }

   return update_derived_params (pt, INFO);
}
/*}}}*/

Param_Info_t *Fit_param_info (Param_t *pt, unsigned int idx) /*{{{*/
{
   return find_param_by_index (pt, idx);
}

/*}}}*/

Param_Info_t *Fit_variable_param_info (Param_t *pt, unsigned int vary_idx) /*{{{*/
{
   return find_param_by_vary_index (pt, vary_idx);
}

/*}}}*/

int Fit_set_param_hard_limits1 (Param_t *pt, int idx,
                                int hard_limits_only, Param_Info_t *pnew)
{
   Param_Info_t *p;

   if (NULL == (p = find_param_by_index (pt, idx)))
     return -1;

   if (hard_limits_only)
     {
        if ((pnew->hard_min > p->min) || (p->max > pnew->hard_max))
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                         "parameter %d: new hard limits do not contain current min/max range", idx);
             return -1;
          }

        p->hard_min = pnew->hard_min;
        p->hard_max = pnew->hard_max;

        return 0;
     }

   if ((pnew->hard_min > pnew->min) || (pnew->min > pnew->value)
       || (pnew->value > pnew->max) || (pnew->max > pnew->hard_max))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "parameter %d: assigning inconsistent parameter limits", idx);
        return -1;
     }

   p->hard_min = pnew->hard_min;
   p->hard_max = pnew->hard_max;
   p->min = pnew->min;
   p->max = pnew->max;
   p->value = pnew->value;

   return update_derived_params (pt, FAIL);
}

int Fit_set_param_control (Param_t *pt, unsigned int idx, int update_minmax,  /*{{{*/
                           double min, double max, int freeze, char *tie,
                           double step, double relstep)
{
   Param_Info_t *p = find_param_by_index (pt, idx);
   if (p == NULL)
     return -1;

   if (update_minmax)
     {
        if (-1 == set_param_minmax (p, min, max))
          return -1;
        p->set_minmax = 1;
     }

   if (p->set_minmax == 0)
     {
        if (-1 == set_param_minmax (p, 0.0, 0.0))
          return -1;
        p->set_minmax = 1;
     }

   /* step/relstep values <= 0 may be interpreted as signals
    * to the optimizer, so we have to accept them.
    */
   if (0 == isnan(step))
     p->step = step;
   if (0 == isnan(relstep))
     p->relstep = relstep;

   ISIS_FREE(p->tie_param_name);
   if ((tie != NULL) && (*tie != 0))
     {
        if (NULL == (p->tie_param_name = isis_make_string (tie)))
          return -1;
     }

   if (Fit_Loading_Parameters_From_File)
     {
        free_param_function (p);
     }

   if (p->fun_str != NULL)
     {
        if (freeze == 0)
          isis_vmesg (WARN, I_ERROR, __FILE__, __LINE__, "param %d is derived; freeze=1 is implied.", idx);
        return 0;
     }

   if (freeze >= 0)
     p->freeze = (freeze == 0) ? 0 : 1;

   return 0;
}

/*}}}*/

int Fit_set_freeze (Param_t *pt, unsigned int idx, unsigned int freeze) /*{{{*/
{
   Param_Info_t *p = find_param_by_index (pt, idx);
   if (p == NULL)
     return -1;

   if (p->fun_str != NULL)
     {
        if (freeze == 0)
          isis_vmesg (WARN, I_ERROR, __FILE__, __LINE__, "param %d is derived; freeze=1 is implied.", idx);
        return 0;
     }

   p->freeze = (freeze == 0) ? 0 : 1;
   return 0;
}

/*}}}*/

int Fit_tie (Param_t *pt, unsigned int idx_a, unsigned int idx_b) /*{{{*/
{
   Param_Info_t *p, *x;
   if (NULL == (p = find_param_by_index (pt, idx_a)))
     return -1;
   /* tie (0, idx) is the same as untie(idx) */
   x = find_param_by_index (pt, idx_b);
   ISIS_FREE(p->tie_param_name);
   if (x != NULL)
     {
        if (NULL == (p->tie_param_name = isis_make_string (x->param_name)))
          return -1;
     }
   else if (idx_b > 0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "cannot tie parameter %d to nonexistent parameter %d",
                    idx_a, idx_b);
        return -1;
     }

   return 0;
}

/*}}}*/

int Fit_untie (Param_t *pt, unsigned int idx) /*{{{*/
{
   Param_Info_t *p;
   if (NULL == (p = find_param_by_index (pt, idx)))
     return -1;
   ISIS_FREE(p->tie_param_name);
   return 0;
}

/*}}}*/

int Fit_set_param_function (Param_t *pt, unsigned int idx, char *str) /*{{{*/
{
   static char *fmt = "define _pf_%u() {return %s;}";
   Param_Info_t *p;
   char fun_name[32];
   char *buf;
   int len;

   if (NULL == (p = find_param_by_index (pt, idx)))
     return -1;

   if (str == NULL)
     {
        if (p->fun_str != NULL) p->freeze = 0;
        free_param_function (p);
        return 0;
     }

   len = strlen (fmt) + strlen (str) + 16;
   if (NULL == (buf = (char *) ISIS_MALLOC (len * sizeof(char))))
     return -1;
   sprintf (buf, fmt, idx, str);
   sprintf (fun_name, "_pf_%u", idx);
   if (-1 == SLang_load_string (buf))
     {
        ISIS_FREE (buf);
        return -1;
     }
   ISIS_FREE (buf);

   SLang_free_function (p->fun_ptr);
   if (NULL == (p->fun_ptr = SLang_get_function (fun_name)))
     return -1;

   ISIS_FREE (p->fun_str);
   if (NULL == (p->fun_str = isis_make_string (str)))
     {
        SLang_free_function (p->fun_ptr);
        p->fun_ptr = NULL;
        return -1;
     }

   p->freeze = 1;

   return update_derived_params (pt, INFO);
}

/*}}}*/

int Fit_set_param_value (Param_t *pt, unsigned int idx, double value) /*{{{*/
{
   Param_Info_t *p;
   if (NULL == (p = find_param_by_index (pt, idx)))
     return -1;
   p->value = value;
   return update_derived_params (pt, INFO);
}

/*}}}*/

int Fit_get_param_value (Param_t *pt, unsigned int idx, double *value) /*{{{*/
{
   Param_Info_t *p;
   if (NULL == (p = find_param_by_index (pt, idx)))
     return -1;
   *value = p->value;
   return 0;
}

/*}}}*/

double *Fit_get_kernel_params (Param_t *pt, int id, Isis_Kernel_Def_t *def) /*{{{*/
{
   double *kp;
   Param_t *fun;
   unsigned int i;

   /* std kernel has kernel_id = 0 */
   if ((def == NULL)
       || (def->kernel_id == 0)
       || (def->num_kernel_parms == 0))
     return NULL;

   if (NULL == (fun = locate_fun_params (pt, def->fun_type, id)))
     return NULL;

   if (NULL == (kp = (double *) ISIS_MALLOC (def->num_kernel_parms * sizeof(double))))
     return NULL;

   for (i = 0; i < def->num_kernel_parms; i++)
     {
        kp[i] = fun->info[i].value;
     }

   return kp;
}

/*}}}*/

/*}}}*/

typedef struct
{
   int all, vary;
}
Count_Type;

static int count (Param_Info_t *p, void *cl) /*{{{*/
{
   Count_Type *num = (Count_Type *)cl;

   if (p->in_use)
     {
        num->all++;
        if ((p->tie_param_name == NULL) && (p->freeze == 0))
          num->vary++;
     }

   return 0;
}

/*}}}*/

int Fit_count_params (Param_t *pt, int *num_all, int *num_vary) /*{{{*/
{
   Count_Type num;

   num.all = 0;
   num.vary = 0;

   if (-1 == map_table (pt, ALL_PARS, count, &num))
     return -1;

   *num_all = num.all;
   *num_vary = num.vary;

   return 0;
}
/*}}}*/

/* pack/unpack parameters for minimizer */

static int pack (Param_Info_t *p, void *cl) /*{{{*/
{
   Fit_Param_t *par = (Fit_Param_t *)cl;
   unsigned int n;

   if (p->in_use)
     {
        n = par->npars;

        p->vary_idx = n;
        par->par[n] = p->value;
        par->par_min[n] = p->min;
        par->par_max[n] = p->max;
        par->step[n] = p->step;
        par->relstep[n] = p->relstep;
        par->idx[n] = p->idx;

        par->npars++;
     }

   return 0;
}

/*}}}*/

int Fit_pack_variable_params (Param_t *pt, Fit_Param_t *p) /*{{{*/
{
   p->npars = 0;

   if (-1 == update_derived_params (pt, INTR))
     return -1;

   if (-1 == map_table (pt, VARIABLE_PARS, pack, p))
     {
        isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "no variable parameters");
        return -1;
     }

   return 0;
}

/*}}}*/

int Fit_pack_all_params (Param_t *pt, Fit_Param_t *p) /*{{{*/
{
   p->npars = 0;

   if (-1 == update_derived_params (pt, INTR))
     return -1;

   if (-1 == map_table (pt, ALL_PARS, pack, p))
     {
        isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "no parameters");
        return -1;
     }

   return 0;
}

/*}}}*/

typedef struct
{
   double *par;
   unsigned int n;
}
Pack_Type;

static int unpack (Param_Info_t *p, void *cl) /*{{{*/
{
   Pack_Type *x = (Pack_Type *) cl;

   if (p->in_use)
     {
        p->value = x->par[x->n];
        x->n++;
     }

   return 0;
}
/*}}}*/

int Fit_unpack_variable_params (Param_t *pt, double *par) /*{{{*/
{
   Pack_Type x;

   x.par = par;
   x.n = 0;

   if (-1 == map_table (pt, VARIABLE_PARS, unpack, &x))
     return -1;

   if (-1 == update_derived_params (pt, INTR))
     return -1;

   return update_tied_params (pt);
}

/*}}}*/

int Fit_unpack_all_params (Param_t *pt, double *par) /*{{{*/
{
   Pack_Type x;

   x.par = par;
   x.n = 0;

   if (-1 == map_table (pt, ALL_PARS, unpack, &x))
     return -1;

   if (-1 == update_derived_params (pt, INTR))
     return -1;

   return update_tied_params (pt);
}

/*}}}*/

int Fit_sync_tied_params (Param_t *pt)
{
   return update_tied_params (pt);
}

int Fit_sync_derived_params (Param_t *pt)
{
   return update_derived_params (pt, INTR);
}
