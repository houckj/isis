/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
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

/* $Id: fit-cmds.c,v 1.266 2004/09/10 02:36:58 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <limits.h>
#include <errno.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>

#include "isis.h"
#include "util.h"
#include "histogram.h"
#include "fit.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

/*{{{ imported function prototypes */

extern Hist_t *get_histogram_list_head (void);

/*}}}*/

static int (*Unpack)(Param_t *, double *);
int Fit_Loading_Parameters_From_File = 0;

int Isis_Voigt_Is_Normalized = 1;

/*{{{ internal globals */

static int Isis_Active_Function_Id;
static int Num_Statistic_Evaluations;

static enum Mode_type
{
   INVALID_MODE   = -1,
   BIN_EVAL_MODE  = 0,
   INIT_MODE      = 1,
   DIFF_EVAL_MODE = 2,
   NAME_EVAL_MODE = 3
}
Mode = BIN_EVAL_MODE;

typedef struct
{
   SLang_Name_Type *fun_ptr;    /* function pointer */
   char *fun_string;            /* definition string */
}
User_Function_Type;

typedef struct Cached_Grid_Type Cached_Grid_Type;
struct Cached_Grid_Type
{
   Cached_Grid_Type *next;
   Isis_Hist_t grid;
   unsigned int id;
   unsigned int type;
   unsigned int updated_cached_model_values;  /* boolean */
};

struct Fit_Data_t
{
   Hist_t **datasets;
   int *offsets;
   int *offset_for_marked;
   Cached_Grid_Type *cache;
   double *data;
   double *weight;
   double *tmp;
   int nbins;
   int num_datasets;
   int nbins_after_datasets_combined;
};

static int Response_Type = USE_ASSIGNED_ARF_RMF;

static char *Fit_Method;
static char *Fit_Statistic;
static SLang_Name_Type *Fit_Constraint;
static SLang_Name_Type *Fit_Range_Hook;
static SLang_Name_Type *Define_Model_Hook;
static Fit_Data_t *Current_Fit_Data_Info;

static Kernel_Table_t _Kernel_Table;
static Kernel_Table_t *Kernel_Table = &_Kernel_Table;

static Isis_Hist_t Eval_Grid;
static Isis_User_Grid_t Differential_Grid;

static int Fit_Verbose;
static int Fit_Store_Model;
static int Use_Interactive_Param_Init;
static int Looking_For_Confidence_Limits;
static int Computing_Statistic_Only;
static unsigned int Fit_Data_Type;

/* total number of parameters */
static unsigned int Num_Params;
static Param_t *Param;
static User_Function_Type User_Function;

static SLang_Name_Type *Array_Fit_Fun;

static int push_kernel (Kernel_Table_t *t, unsigned int hist_index,
                        Isis_Kernel_Def_t *def, char *params);
static int pop_kernel (Kernel_Table_t *t, unsigned int hist_index, Kernel_Info_t *ki);
static Cached_Grid_Type *find_cached_grid (Cached_Grid_Type *t, unsigned int type, unsigned int id);
static void cached_models_need_updating (Fit_Data_t *d);

/*}}}*/

static Isis_Hist_t *get_evaluation_grid (void) /*{{{*/
{
   return &Eval_Grid;
}

/*}}}*/

static int map_datasets (int (*fun)(Hist_t *, void *), void *cl) /*{{{*/
{
   enum {check_exclude = 1};
   Hist_t *head;

   if (NULL == (head = get_histogram_list_head ()))
     return 0;

   return Hist_map (head, fun, cl, check_exclude);
}

/*}}}*/

/*{{{ kernel handling */

Kernel_Table_t *get_kernel_table (void) /*{{{*/
{
   return Kernel_Table;
}

/*}}}*/

static int pre_fit (Hist_t *h, void *cl) /*{{{*/
{
   (void) cl;

   if (is_flux (Fit_Data_Type))
     {
        Kernel_Table_t *t = Kernel_Table;
        Isis_Kernel_Def_t *def;
        if (NULL == (def = Fit_find_kernel_by_name (t->kernel_defs, "std")))
          return -1;
        if (-1 == push_kernel (t, Hist_get_index(h), def, NULL))
          return -1;
     }

   return 0;
}

/*}}}*/

static int post_fit (Hist_t *h, void *cl) /*{{{*/
{
   (void) cl;

   if (is_flux (Fit_Data_Type))
     {
        Kernel_Info_t p;
        return pop_kernel (Kernel_Table, Hist_get_index(h), &p);
     }

   return 0;
}

/*}}}*/

static Isis_Kernel_Def_t *get_histogram_kernel (Kernel_Table_t *t, Hist_t *h) /*{{{*/
{
   Kernel_Info_t *p;

   if (NULL == (p = find_kernel_info (t, Hist_get_index (h))))
     return NULL;

   return Fit_find_kernel (t->kernel_defs, p->kernel_id);
}

/*}}}*/

static char *get_kernel_params (Kernel_Table_t *t, Hist_t *h) /*{{{*/
{
   Kernel_Info_t *p;

   if (NULL == (p = find_kernel_info (t, Hist_get_index (h))))
     return NULL;

   return p->kernel_params;
}

/*}}}*/

int get_kernel_params_for_hist (Hist_t *h, double **params, unsigned int *num) /*{{{*/
{
   Isis_Kernel_t *k;
   Isis_Kernel_Def_t *def;

   if (h == NULL)
     return -1;

   k = Hist_get_kernel (h);
   if ((k == NULL) || (k->kernel_def == NULL))
     return -1;

   def = k->kernel_def;
   *num = def->num_kernel_parms;
   *params = Fit_get_kernel_params (Param, Hist_get_index (h), def);

   return 0;
}

/*}}}*/

static int assign_kernel_to_histogram (Kernel_Table_t *t, unsigned int hist_index, /*{{{*/
                                       Isis_Kernel_Def_t *def, char *params)
{
   Kernel_Info_t *p;

   if (NULL == (p = find_kernel_info (t, hist_index)))
     return -1;

   p->kernel_id = def->kernel_id;
   ISIS_FREE (p->kernel_params);
   if (NULL == (p->kernel_params = isis_make_string (params)))
     return -1;

   return 0;
}

/*}}}*/

static int pop_kernel (Kernel_Table_t *t, unsigned int hist_index, Kernel_Info_t *ki)
{
   Kernel_Info_t *p, *saved;

   if (NULL == (p = find_kernel_info (t, hist_index)))
     return -1;

   if (p->saved == NULL)
     return -1;

   ki->kernel_id = p->kernel_id;
   ki->kernel_params = p->kernel_params;

   saved = p->saved;
   p->kernel_id = saved->kernel_id;
   p->kernel_params = saved->kernel_params;

   ISIS_FREE(p->saved);

   return 0;
}

static int push_kernel (Kernel_Table_t *t, unsigned int hist_index,
                        Isis_Kernel_Def_t *def, char *params)
{
   Kernel_Info_t *p, *saved;
   char *s = NULL;

   if (NULL == (p = find_kernel_info (t, hist_index)))
     return -1;

   if ((params != NULL) && (NULL == (s = isis_make_string (params))))
     return -1;

   if (NULL == (p->saved = (Kernel_Info_t *) ISIS_MALLOC(sizeof(*saved))))
     return -1;

   saved = p->saved;
   saved->kernel_id = p->kernel_id;
   saved->kernel_params = p->kernel_params;

   p->kernel_id = def->kernel_id;
   p->kernel_params = s;

   return 0;
}

static void _print_kernel (int * hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   (void) Hist_print_kernel (h);
}

/*}}}*/

static int _load_kernel (char *libfile, char *init_name) /*{{{*/
{
   Kernel_Table_t *t = Kernel_Table;
   Isis_Kernel_Def_t *new_def = NULL;
   char *init_args = NULL;

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     SLdo_pop ();
   else if (-1 == SLpop_string (&init_args))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if ((NULL == libfile)
       || (NULL == init_name))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if ((NULL == (new_def = Fit_new_kernel ()))
       || (-1 == Hist_load_aux_kernel (libfile, init_name, init_args, new_def)))
     {
        Fit_free_kernel (new_def);
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (-1 == Fit_append_kernel (t->kernel_defs, new_def))
     {
        Fit_free_kernel (new_def);
        new_def = NULL;
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return 0;
}
/*}}}*/

static void _list_kernels (void) /*{{{*/
{
   Kernel_Table_t *t = Kernel_Table;
   Fit_push_kernel_names (t->kernel_defs);
}

/*}}}*/

/*}}}*/

static User_Function_Type *get_user_function (void) /*{{{*/
{
   return &User_Function;
}

/*}}}*/

static void reset_user_function (void) /*{{{*/
{
   Num_Params = 0;
}

/*}}}*/

static void free_user_function (void) /*{{{*/
{
   ISIS_FREE (User_Function.fun_string);
}

/*}}}*/

static void delete_user_model (void) /*{{{*/
{
   ISIS_FREE(User_Function.fun_string);
   User_Function.fun_ptr = NULL;
}

/*}}}*/

static int reset_param_lookup_table (Param_t **pt) /*{{{*/
{
   if (pt == NULL)
     return -1;

   if (*pt == NULL)
     {
        *pt = Fit_new_param_table (0, UINT_MAX, UINT_MAX, UINT_MAX);
        if (NULL == *pt)
          return -1;
     }

   reset_user_function ();
   return Fit_mark_params_unused (*pt);
}

/*}}}*/

/*{{{ set/freeze/thaw/tie/untie/list/save/load/edit params */

static void push_num_params (void) /*{{{*/
{
   SLang_push_integer (Num_Params);
}

/*}}}*/

static void _get_par (unsigned int *idx) /*{{{*/
{
   double value;

   if (-1 == Fit_get_param_value (Param, *idx, &value))
     {
        isis_vmesg (INTR, I_INVALID, __FILE__, __LINE__, "parameter index=%d", *idx);
        return;
     }

   SLang_push_double (value);
}
/*}}}*/

static int get_param_name (unsigned int idx, char *name, int len) /*{{{*/
{
   Param_Info_t *p;
   if (NULL == (p = Fit_param_info (Param, idx)))
     return -1;
   isis_strcpy (name, p->param_name, len);
   return 0;
}

/*}}}*/

static int get_variable_param_name (int vary_idx, char *name, int len) /*{{{*/
{
   Param_Info_t *p;
   if (NULL == (p = Fit_variable_param_info (Param, vary_idx)))
     return -1;
   isis_strcpy (name, p->param_name, len);
   return 0;
}

/*}}}*/

static int get_index_for_param_name (char *name) /*{{{*/
{
   Param_Info_t *p = Fit_find_param_info_by_full_name (Param, name);
   if (p == NULL)
     return -1;
   return p->idx;
}

/*}}}*/

typedef struct
{
   char *name;
   char *fun_str;
   char *tie;
   char *units;
   int idx, freeze, is_a_norm;
   double value, step;
   double min, max, hard_min, hard_max;
   int malloced_fun_str;
}
_Param_Info_Type;

static SLang_CStruct_Field_Type _Param_Info_Type_Layout [] =
{
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, idx, "index", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, value, "value", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, min, "min", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, max, "max", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, hard_min, "hard_min", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, hard_max, "hard_max", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, step, "step", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, freeze, "freeze", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, tie, "tie", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, units, "units", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, is_a_norm, "is_a_norm", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (_Param_Info_Type, fun_str, "fun", SLANG_STRING_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static void _free_param_info (_Param_Info_Type *pi) /*{{{*/
{
   if (pi == NULL)
     return;

   ISIS_FREE(pi->name);
   if (pi->malloced_fun_str)
     ISIS_FREE(pi->fun_str);
   if (pi->tie)
     ISIS_FREE(pi->tie);
   if (pi->units)
     ISIS_FREE(pi->units);
}

/*}}}*/

static int _copy_param_info (_Param_Info_Type *pi, unsigned int idx) /*{{{*/
{
   unsigned int name_size = OUT_PARAM_NAME_SIZE;
   Param_Info_t *p;
   Fit_Fun_t *ff;

   if (pi == NULL)
     return -1;

   memset ((char *)pi, 0, sizeof(*pi));

   if (NULL == (p = Fit_param_info (Param, idx)))
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "parameter index=%d", idx);
        return -1;
     }

   if (p->fun_str)
     {
        if (NULL == (pi->fun_str = isis_make_string (p->fun_str)))
          return -1;
        pi->malloced_fun_str = 1;
     }
   else pi->fun_str = "";

   if (NULL == (pi->name = (char *) ISIS_MALLOC(name_size)))
     return -1;

   if (-1 == get_param_name (p->idx, pi->name, name_size))
     return -1;

   pi->idx = p->idx;
   pi->value = p->value;
   pi->min = p->min;
   pi->max = p->max;
   pi->hard_min = p->hard_min;
   pi->hard_max = p->hard_max;
   pi->freeze = p->freeze;
   pi->is_a_norm = p->is_a_norm;
   pi->step = p->step;

   if (p->tie_param_name)
     {
        if (NULL == (pi->tie = isis_make_string (p->tie_param_name)))
          return -1;
     }

   if (NULL == (ff = Fit_get_fit_fun (p->fun_type)))
     return -1;

   if (NULL == (pi->units = isis_make_string (ff->unit[p->fun_par])))
     return -1;

   return 0;
}

/*}}}*/

static void _get_param_info (int *idx) /*{{{*/
{
   _Param_Info_Type pi;

   if (-1 == _copy_param_info (&pi, *idx))
     {
        isis_throw_exception (Isis_Error);
        _free_param_info (&pi);
        return;
     }

   (void) SLang_push_cstruct ((VOID_STAR)&pi, _Param_Info_Type_Layout);

   _free_param_info (&pi);
}

/*}}}*/

static int do_set_par (_Param_Info_Type *pi, int update_minmax) /*{{{*/
{
   Param_Info_t *p;

   if (-1 == Fit_set_param_control (Param, pi->idx, update_minmax, pi->min, pi->max, pi->freeze, pi->tie))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting parameter %d", pi->idx);
        return -1;
     }

   if (NULL == (p = Fit_param_info (Param, pi->idx)))
     return -1;

   if ((p->min <= pi->value) && (pi->value <= p->max))
     return Fit_set_param_value (Param, pi->idx, pi->value);

   return 1;
}

/*}}}*/

static int run_set_par_hook (char *hook_name, _Param_Info_Type *pi) /*{{{*/
{
   _Param_Info_Type pit;
   int status;

   if (pi->tie == NULL)
     {
        if (NULL == (pi->tie = (char *) ISIS_MALLOC (1)))
          return -1;
        pi->tie[0] = 0;
     }

   SLang_start_arg_list ();
   (void) SLang_push_cstruct ((VOID_STAR)pi, _Param_Info_Type_Layout);
   SLang_end_arg_list ();

   SLang_execute_function (hook_name);

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&pit, _Param_Info_Type_Layout))
     return -1;

   status = do_set_par (&pit, 1);
   SLang_free_cstruct ((VOID_STAR)&pit, _Param_Info_Type_Layout);

   return status;
}

/*}}}*/

static int set_par (unsigned int idx, int p_tie, int p_freeze, /*{{{*/
                    double p_value, double p_min, double p_max,
                    int update_minmax, double p_step)
{
   static char hook_name[] = "isis_set_par_hook";
   char tie_name[OUT_PARAM_NAME_SIZE];
   int name_size = OUT_PARAM_NAME_SIZE;
   _Param_Info_Type pi;
   Param_Info_t *p;
   int status;

   if (-1 == _copy_param_info (&pi, idx))
     {
        _free_param_info (&pi);
        isis_throw_exception (Isis_Error);
        return -1;
     }

   pi.idx = idx;
   pi.value = p_value;
   pi.freeze = p_freeze;

   if (p_step >= 0.0)
     pi.step = p_step;

   if (update_minmax)
     {
        pi.min = p_min;
        pi.max = p_max;
     }

   if (p_tie <= 0)
     ISIS_FREE(pi.tie);
   else if ((-1 == get_param_name (p_tie, tie_name, name_size))
            || (NULL == (pi.tie = isis_make_string (tie_name))))
     {
        _free_param_info (&pi);
        isis_throw_exception (Isis_Error);
        return -1;
     }

   status = do_set_par (&pi, update_minmax);

   if (status > 0)
     {
        if (2 == SLang_is_defined (hook_name))
          status = run_set_par_hook (hook_name, &pi);
        else
          {
             p = Fit_param_info (Param, idx);
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "param %d not set:  value %0.5g lies outside [min, max] interval [%0.5g, %0.5g]",
                         idx, p_value, p->min, p->max);
          }
     }

   _free_param_info (&pi);
   if (status) isis_throw_exception (Isis_Error);
   return status;
}
/*}}}*/

static void set_params (void) /*{{{*/
{
   int idx, p_tie, p_freeze, update_minmax;
   double p_value, p_min, p_max, p_step;

   /* (idx, tie, freeze, value, min, max, update_minmax, step) */

   SLang_pop_double (&p_step);
   SLang_pop_integer (&update_minmax);

   SLang_pop_double (&p_max);
   SLang_pop_double (&p_min);
   SLang_pop_double (&p_value);

   SLang_pop_integer (&p_freeze);
   SLang_pop_integer (&p_tie);
   SLang_pop_integer (&idx);

   if (set_par (idx, p_tie, p_freeze, p_value, p_min, p_max, update_minmax, p_step))
     isis_throw_exception (Isis_Error);
}

/*}}}*/

/*}}}*/

static void _get_fit_fun (void) /*{{{*/
{
   User_Function_Type *f;
   char *s = NULL;

   f = get_user_function ();
   if (f == NULL)
     isis_throw_exception (Isis_Error);
   else
     {
        s = f->fun_string;
        if (s != NULL)
          {
             char *p = strchr (s, '\n');
             if (p) *p = 0;
          }
     }

   if (-1 == SLang_push_string (s))
     isis_throw_exception (Isis_Error);
}

/*}}}*/

static void _get_par_fun (int *idx) /*{{{*/
{
   Param_Info_t *p = Fit_param_info (Param, *idx);

   if ((p == NULL)
       || (-1 == SLang_push_string (p->fun_str)))
     {
        isis_throw_exception (Isis_Error);
     }
}

/*}}}*/

static void _set_par_fun (int *idx, char *fun_str) /*{{{*/
{
   unsigned int uidx = (unsigned int) *idx;

   if (*fun_str == 0)
     fun_str = NULL;

   if (Fit_set_param_function (Param, uidx, fun_str))
     isis_throw_exception (Isis_Error);
}

/*}}}*/

static int _define_user_model (char * fun_body);

static int parse_comment_line (char *line, int line_num, char *fname, unsigned int idx) /*{{{*/
{
   char *p;

   if (0 != strncmp (line, "#=>", 3))
     return 0;

   p = strchr (line, '\n');
   if (p) *p = 0;

   if (-1 == Fit_set_param_function (Param, idx, &line[3]))
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "parsing line %d, file %s",
                    line_num, fname);
        return -1;
     }

   return 0;
}

/*}}}*/

int Fit_copy_fun_params (char *fun_name, unsigned int fun_id, double **par, unsigned int *num_pars) /*{{{*/
{
   Fit_Fun_t *ff;
   double *p;
   int fun_type;

   *num_pars = 0;
   *par = NULL;

   if (-1 == (fun_type = Fit_get_fun_type (fun_name)))
     return -1;

   if (NULL == (ff = Fit_get_fit_fun (fun_type)))
     return -1;

   if (ff->nparams == 0)
     return 0;

   if (NULL == (p = (double *) ISIS_MALLOC (ff->nparams * sizeof(double))))
     return -1;

   if (-1 == Fit_get_fun_params (Param, fun_type, fun_id, p))
     {
        ISIS_FREE (p);
        return -1;
     }

   *par = p;
   *num_pars = ff->nparams;

   return 0;
}

/*}}}*/

static int parse_param_info (char *line, int line_num, char *fname, unsigned int *idx) /*{{{*/
{
   Param_Info_t *p;
   Fit_Fun_t *ff;
   double p_value, p_min, p_max;
   int n, fun_type, fun_id, fun_par, p_tie, p_freeze;
   char fun_name[MAX_NAME_SIZE], par_name[MAX_NAME_SIZE];

   /* Although identifiers of the form foo(n).bar are expected,
    * identifiers of the form pileup<n>.foo are recognized for
    * back-compatibility
    */
   n = sscanf (line, "%*u %[^(<]%*1[(<]%d%*1[)>].%s %d %d %le %le %le",
               fun_name, &fun_id, par_name,
               &p_tie, &p_freeze, &p_value, &p_min, &p_max);

   if (n == 0)
     return 1;
   else if (n != 8)
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "parsing line %d, file %s",
                    line_num, fname);
        return -1;
     }

   if (-1 == (fun_type = Fit_get_fun_type (fun_name)))
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "unrecognized function name `%s'",
                    fun_name);
        return -1;
     }

   if (NULL == (ff = Fit_get_fit_fun (fun_type)))
     return -1;

   if (-1 == (fun_par = Fit_get_fun_par (ff, par_name)))
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "unrecognized parameter name `%s'",
                    par_name);
        return -1;
     }

   if (NULL == (p = Fit_param_info2 (Param, fun_type, fun_id, fun_par)))
     return -1;

   *idx = p->idx;

   /* Currently no support for setting param step from parameter file. */
   return set_par (*idx, p_tie, p_freeze, p_value, p_min, p_max, 1, -1.0);
}

/*}}}*/

static int load_params (char *fname) /*{{{*/
{
   FILE *fp;
   char line[BUFSIZE];
   unsigned int i, num_params;
   unsigned int idx = 0;
   int line_num = 0;
   int ret = 0;

   if ((0 == is_regular_file (fname))
       || (NULL == (fp = fopen (fname, "r"))))
     {
        isis_vmesg (INTR, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", fname);
        return -1;
     }

   /* skip comment lines */
   for (;;)
     {
        if (NULL == fgets (line, sizeof (line), fp))
          {
             fclose (fp);
             isis_vmesg (INTR, I_READ_FAILED, __FILE__, __LINE__, "%s", fname);
             return -1;
          }
        if (line[0] != COMMENT_CHAR)
          break;
        line_num++;
     }

   if (NULL == strchr (line, '\n'))
     {
        fclose (fp);
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "model string is too long");
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__,
                    "*** You may wish to consider expressing the model as a S-Lang function");
        return -1;
     }

   if (-1 == _define_user_model (line))
     {
        fclose (fp);
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "%s:  invalid model definition", fname);
        return -1;
     }

   num_params = Num_Params;
   Fit_Loading_Parameters_From_File = 1;

   i = 0;
   do
     {
        line_num++;

        if (NULL == fgets (line, sizeof (line), fp))
          break;

        if (*line == COMMENT_CHAR)
          ret = parse_comment_line (line, line_num, fname, idx);
        else if (i < num_params)
          {
             ret = parse_param_info (line, line_num, fname, &idx);
             if (ret == 0) i++;
          }
        else break;

        /* don't stop on (i < num_params) because trailing
         * comment line may hold parameter linkage information
         */
     } while (ret >= 0);

   Fit_Loading_Parameters_From_File = 0;
   fclose(fp);

   if (ret == -1)
     isis_throw_exception (Isis_Error);

   return ret;
}
/*}}}*/

static int save_params (char *fname) /*{{{*/
{
   return SLang_run_hooks ("save_par", 1, fname);
}
/*}}}*/

static void edit_params (char *file) /*{{{*/
{
   if (-1 == edit_temp_file (save_params, load_params,
                             (*file == 0) ? NULL : file))
     {
        isis_throw_exception (Isis_Error);
     }
}

/*}}}*/

static void _freeze (int *idx) /*{{{*/
{
   if (-1 == Fit_set_freeze (Param, *idx, 1))
     {
        isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "param %d not frozen", *idx);
     }
}

/*}}}*/

static void _thaw (int *idx) /*{{{*/
{
   if (-1 == Fit_set_freeze (Param, *idx, 0))
     {
        isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "param %d not thawed", *idx);
     }
}

/*}}}*/

static void _tie (int *idx_a, int *idx_b) /*{{{*/
{
   /* tie (x,x) is a no-op */
   if (*idx_a == *idx_b)
     return;

   /* tie (idx, 0) is the same as untie(idx) */
   if (-1 == Fit_tie (Param, *idx_a, *idx_b))
     {
        isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "params (%d, %d) not tied",
                    *idx_a, *idx_b);
     }
}

/*}}}*/

static void _untie (int *idx) /*{{{*/
{
   if (-1 == Fit_untie (Param, *idx))
     {
        isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "failed to untie param %d", *idx);
     }
}

/*}}}*/

/*}}}*/

/*{{{ mode switching, etc. */

static void bin_eval_mode (void) /*{{{*/
{
   Mode = BIN_EVAL_MODE;
}

/*}}}*/

static int do_init_mode_eval (SLang_Name_Type *fun_ptr) /*{{{*/
{
   enum Mode_type tmp = Mode;
   double d1;

   if (fun_ptr == NULL)
     goto error_return;

   Mode = INIT_MODE;
   if (-1 == SLexecute_function (fun_ptr)
       || -1 == SLang_pop_double (&d1))
     goto error_return;
   Mode = tmp;
   return 0;

   error_return:
   bin_eval_mode ();
   reset_param_lookup_table (&Param);
   isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "initializing user model");
   return -1;
}

/*}}}*/

static int do_mode_eval (SLang_Name_Type *fun_ptr, enum Mode_type mode) /*{{{*/
{
   enum Mode_type tmp = Mode;
   int ret = 0;

   if (fun_ptr == NULL)
     return -1;

   Mode = mode;

   if (-1 == SLexecute_function (fun_ptr))
     ret = -1;

   Mode = tmp;

   return ret;
}

/*}}}*/

static int init_fun_params (Fit_Fun_t *ff, unsigned int fun_id) /*{{{*/
{
   double *par = NULL;
   double *par_min = NULL;
   double *par_max = NULL;
   int ret = -1;

   if (NULL == ff->s.init_params_from_screen)
     return -1;

   if (NULL == (par = (double *) ISIS_MALLOC (3 * ff->nparams * sizeof(double))))
     goto finish;

   par_min = par + ff->nparams;
   par_max = par + 2*ff->nparams;

   if ((-1 == (*ff->s.init_params_from_screen) (fun_id, par, par_min, par_max, ff->nparams))
       || (-1 == Fit_set_fun_params (Param, ff->fun_type, fun_id, par, par_min, par_max)))
     goto finish;

   ret = 0;

   finish:
   ISIS_FREE (par);

   return ret;
}
/*}}}*/

static int push_parameter_value_struct (Fit_Fun_t *ff, unsigned int fun_id) /*{{{*/
{
   SLtype *types = NULL;
   char **names = NULL;
   double *par = NULL;
   double **pvalues = NULL;
   unsigned int i, n;
   int ret = -1;

   if (ff == NULL)
     return -1;

   n = ff->nparams;

   if ((NULL == (par = (double *) ISIS_MALLOC (n * sizeof(double))))
       || (NULL == (names = (char **) ISIS_MALLOC (n * sizeof(char*))))
       || (NULL == (types = (SLtype *) ISIS_MALLOC (n * sizeof(SLtype))))
       || (NULL == (pvalues = (double **) ISIS_MALLOC (n * sizeof(double *)))))
     goto finish;

   memset ((char *)par, 0, n * sizeof(double));
   ret = Fit_get_fun_params (Param, ff->fun_type, fun_id, par);

   for (i = 0; i < n; i++)
     {
        names[i] = ff->name[i+1];
        types[i] = SLANG_DOUBLE_TYPE;
        pvalues[i] = &par[i];
     }

   finish:

   (void) SLstruct_create_struct (n, names, types, (void **) pvalues);

   ISIS_FREE (par);
   ISIS_FREE (names);
   ISIS_FREE (types);
   ISIS_FREE (pvalues);

   return ret;
}
/*}}}*/

static int init_eval (Fit_Fun_t *ff, unsigned int fun_id, unsigned int num_extra_args) /*{{{*/
{
   if (NULL == ff)
     return -1;

   if (-1 == Fit_register_fun (Param, ff, fun_id, &Num_Params))
     return -1;

   if (num_extra_args > 0)
     SLdo_pop_n (num_extra_args);

   SLang_push_double (1.0);

   if (!Use_Interactive_Param_Init)
     return 0;

   if (current_plot_shows_bin_density () <= 0)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "Interactive param init requires bin-density plot");
        return -1;
     }

   return init_fun_params (ff, fun_id);
}

/*}}}*/

static int name_eval (Fit_Fun_t *ff, unsigned int fun_id, unsigned int num_extra_args) /*{{{*/
{
   char buf[2*MAX_NAME_SIZE];
   int ret;

   if (ff == NULL)
     return -1;

   sprintf (buf, "%s(%d)", ff->name[0], fun_id);

   ret = SLang_run_hooks ("_isis->name_mode_eval_hook", 1, buf);

   if (num_extra_args > 0)
     SLdo_pop_n (num_extra_args);

   SLang_push_double (1.0);

   return ret;
}
/*}}}*/

static int diff_eval (Fit_Fun_t *ff, unsigned int fun_id, unsigned int num_extra_args) /*{{{*/
{
   double *par = NULL;

   (void) num_extra_args;

   if (ff == NULL)
     return -1;

   if (NULL == (par = (double *) ISIS_MALLOC (ff->nparams * sizeof(double))))
     return -1;

   if (-1 == Fit_get_fun_params (Param, ff->fun_type, fun_id, par))
     {
        ISIS_FREE (par);
        return -1;
     }

   if (-1 == (*ff->diff_eval_method)(ff, &Differential_Grid, par))
     {
        ISIS_FREE (par);
        SLang_push_double (1.0);
        return -1;
     }

   ISIS_FREE (par);

   return 0;
}
/*}}}*/

static int bin_eval (Fit_Fun_t *ff, unsigned int fun_id, unsigned int num_extra_args) /*{{{*/
{
   double *par = NULL;
   Isis_Hist_t *g;
   int ret;

   (void) num_extra_args;

   if ((NULL == ff) || (NULL == ff->fun.c))
     return -1;

   if (NULL == (g = get_evaluation_grid ()))
     return -1;

   if ((g->n_notice > 0)
       && ((g->notice_list == NULL)
           || (g->bin_lo == NULL)
           || (g->bin_hi == NULL)))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Invalid model evaluation grid");
        return -1;
     }

   if (ff->nparams)
     {
        par = (double *) ISIS_MALLOC (ff->nparams * sizeof(double));
        if (NULL == par)
          return -1;

        if (-1 == Fit_get_fun_params (Param, ff->fun_type, fun_id, par))
          {
             ISIS_FREE (par);
             return -1;
          }
     }
   else par = NULL;

   ret = (*ff->bin_eval_method)(ff, g, par);
   if (ret) SLang_push_double (1.0);
   ISIS_FREE (par);
   return ret;
}
/*}}}*/

typedef struct
{
   int (*task)(Fit_Fun_t *, unsigned int, unsigned int);
   enum Mode_type mode;
   char *err_msg;
}
Mode_Task_Type;

/* Each fit intrinsic (e.g. Gaussian) supplies a function pointer
 * that gets called when that intrinsic is evaluated in the
 * corresponding mode.  Context determines the mode.
 *
 * INIT_MODE
 * allocates and initializes space in the parameter table for that
 * instance of the intrinsic and optionally obtains initial parameter
 * values by reading cursor coordinates off an interactive plot.
 * e.g. fit_fun ("gauss(1)"); implies (INIT_MODE, fun_type, fun_id)
 *
 * BIN_EVAL_MODE
 * evaluates the wavelength-bin-integrated form of the fit-intrinsic
 * and pushes the result onto the S-Lang stack.
 * e.g.  () = fit_counts (); implies BIN_EVAL_MODE
 *
 * DIFF_EVAL_MODE
 * evaluates the differential form of the fit-intrinsic and pushes the
 * result onto the S-Lang stack.
 * e.g. y = get_cfun (x);  implies DIFF_EVAL_MODE
 *
 * NAME_EVAL_MODE
 * compiles a list of function component names as a S-Lang array
 * e.g.  "gauss(1)*gauss(2)"   =>   ["gauss(1)",  "gauss(2)"]
*/

static Mode_Task_Type Mode_Table[] =
{
     {init_eval,  INIT_MODE,      "function definition failed"},
     {bin_eval,   BIN_EVAL_MODE,  "function evaluation failed"},
     {diff_eval,  DIFF_EVAL_MODE, "function evaluation failed"},
     {name_eval,  NAME_EVAL_MODE, "failed retrieving function name"},
     {NULL, INVALID_MODE, "internal error:  invalid mode in mode_switch"}
};

/* I'm told that this code is hard to understand. Is it? */
static void mode_switch (unsigned int *fun_id, unsigned int *fun_type, unsigned int *num_extra_args) /*{{{*/
{
   Fit_Fun_t *ff = Fit_get_fit_fun (*fun_type);
   Mode_Task_Type *t = Mode_Table;

   if (Isis_Evaluating_Derived_Param)
     {
        push_parameter_value_struct (ff, *fun_id);
        return;
     }

   while ((t->task != NULL) && (t->mode != Mode))
     t++;

   Isis_Active_Function_Id = *fun_id;
   if ((t->task == NULL) || (-1 == (*t->task)(ff, *fun_id, *num_extra_args)))
     {
        int severity = Looking_For_Confidence_Limits ? FAIL : INTR;
        isis_vmesg (severity, I_ERROR, __FILE__, __LINE__, "%s", t->err_msg);
     }
   Isis_Active_Function_Id = 0;
}

/*}}}*/

/*}}}*/

/*{{{ model and kernel definition */

static int show_kernel_help (Isis_Option_Type *o) /*{{{*/
{
   unsigned int k = 0;

   if ((o == NULL) || (o->num_options == 0))
     return 0;

   while (k < o->num_options)
     {
        if (0 == isis_strcasecmp ("help", o->option_names[k++]))
          return 1;
     }

   return 0;
}

/*}}}*/

static void _set_kernel (unsigned int *hist_index) /*{{{*/
{
   Kernel_Table_t *t = Kernel_Table;
   Isis_Kernel_Def_t *def = NULL;
   Isis_Option_Type *o = NULL;
   char *options = NULL;

   if ((*hist_index == 0)
       || (-1 == SLpop_string (&options)))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "data set %d, kernel '%s', kernel not set",
                    *hist_index,
                    options ? options : "<null>");
        return;
     }

   if (NULL == (o = isis_parse_option_string (options)))
     {
        SLfree (options);
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "parsing kernel specifier for data set %d, kernel not set",
                    *hist_index);
        return;
     }

   if (NULL == (def = Fit_find_kernel_by_name (t->kernel_defs, o->subsystem)))
     {
        isis_vmesg (INTR, I_INVALID, __FILE__, __LINE__, "kernel %s", o->subsystem);
        SLfree (options);
        isis_free_options (o);
        return;
     }

   if (show_kernel_help (o))
     {
        Hist_t *h = find_hist (*hist_index);
        if (h != NULL)
          (void) Hist_reallocate_kernel (h, options, def);
     }
   else if (-1 == assign_kernel_to_histogram (Kernel_Table, *hist_index, def, options))
     {
        isis_throw_exception (Isis_Error);
        SLfree (options);
        isis_free_options (o);
        return;
     }

   SLfree (options);
   isis_free_options (o);
   update_user_model ();
}

/*}}}*/

static SLang_Name_Type *create_registered_model (char *fun_name, char *fun_body) /*{{{*/
{
   static char fmt[] = "define %s() {return %s;}";
   SLang_Name_Type *fun_ptr = NULL;
   char *model = NULL;
   int model_size, status;

   model_size = strlen (fun_body) + strlen (fun_name) + sizeof(fmt);
   if (NULL == (model = (char *) ISIS_MALLOC (model_size * sizeof(char))))
     return NULL;
   sprintf (model, fmt, fun_name, fun_body);

   status = SLang_load_string (model);
   ISIS_FREE (model);
   if (status == -1) return NULL;

   if ((NULL == (fun_ptr = SLang_get_function (fun_name)))
       || (-1 == do_init_mode_eval (fun_ptr)))
     return NULL;

   return fun_ptr;
}

/*}}}*/

static int register_kernel_pars (Hist_t *h, void *cl) /*{{{*/
{
   char name[MAX_NAME_SIZE+16], s[MAX_NAME_SIZE+16];
   Isis_Kernel_Def_t *def = NULL;
   int idx;

   (void) cl;

   if (Hist_num_data_noticed (h) < 1)
     return 0;

   if (NULL == (def = get_histogram_kernel (Kernel_Table, h)))
     return -1;

   /* nothing to do for the standard kernel */
   if (def->kernel_id == 0)
     return 0;

   idx = Hist_get_index (h);
   sprintf (name, "_kern_%d", idx);
   sprintf (s, "%s(%d)", def->kernel_name, idx);

   /* We won't actually call this function, we just
    * want the parameters registered */
   if (NULL == create_registered_model (name, s))
     return -1;

   return 0;
}

/*}}}*/

static int register_iback_pars (Hist_t *h, void *cl) /*{{{*/
{
   SLang_Name_Type *fun_ptr;
   char name[MAX_NAME_SIZE];
   char *s;

   (void) cl;

   if (Hist_num_data_noticed (h) < 1)
     return 0;

   /* its ok if there's no instrumental background function */
   if (NULL == (s = Hist_get_instrumental_background_hook_name (h)))
     return 0;

   sprintf (name, "_iback_%d", Hist_get_index (h));
   fun_ptr = create_registered_model (name, s);

   return Hist_set_instrumental_background_hook (h, fun_ptr);
}

/*}}}*/

typedef struct
{
   char *fun_name;
   char *fun_body;
   User_Function_Type *f;
}
Fitfun_Register_Type;

static int register_fitfun_pars (Hist_t *h, void *cl) /*{{{*/
{
   Fitfun_Register_Type *fr = (Fitfun_Register_Type *)cl;
   User_Function_Type *f;
   char *name;
   char *body;

   if (h != NULL)
     {
        SLang_Name_Type *fun_ptr;
        if (NULL != (fun_ptr = Hist_assigned_model (h)))
          {
             Isis_Arg_Type *args;
             if (NULL != (args = Hist_assigned_model_args (h)))
               {
                  (void) isis_push_args (args);
               }
             return do_init_mode_eval (fun_ptr);
          }
     }

   if (fr == NULL)
     return -1;

   name = fr->fun_name;
   body = fr->fun_body;
   f = fr->f;

   if (NULL == (f->fun_ptr = create_registered_model (name, body)))
     return -1;

   return 0;
}

/*}}}*/

static int no_noticed_data (void) /*{{{*/
{
  return (0 == Hist_num_noticed_histograms (get_histogram_list_head()));
}

/*}}}*/

static int register_constraint_pars (void) /*{{{*/
{
   static char name[] = "__isis_fit_constraint__";
   static char body[] = "constraint(1)";

   if (Fit_Constraint == NULL)
     return 0;

   /* We won't actually call this function, we just
    * want the parameters registered */
   if (NULL == create_registered_model (name, body))
     return -1;

   return 0;
}

/*}}}*/

static int do_define_user_model (char *fun_body) /*{{{*/
{
   User_Function_Type *f = get_user_function ();
   Fitfun_Register_Type fr;

   if (-1 == reset_param_lookup_table (&Param))
     return -1;

   if ((fun_body == NULL) || (NULL == f))
     return -1;
   ISIS_FREE (f->fun_string);
   if (NULL == (f->fun_string = isis_make_string (fun_body)))
     return -1;

   fr.f = f;
   fr.fun_name = "__isis_tmp_ffname__";
   fr.fun_body = fun_body;

   if (no_noticed_data())
     {
        if (-1 == register_fitfun_pars (NULL, &fr))
          return -1;
     }
   else if (-1 == map_datasets (register_fitfun_pars, &fr))
     return -1;

   if (-1 == map_datasets (register_kernel_pars, NULL)
       || -1 == map_datasets (register_iback_pars, NULL))
     return -1;

   if (-1 == register_constraint_pars ())
     return -1;

   return 0;
}

/*}}}*/

static int Defining_User_Model;
static int _define_user_model (char *fun_body) /*{{{*/
{
   int status;

   if (Defining_User_Model)
     {
        SLang_verror (Isis_Error, "invalid model definition");
        return -1;
     }

   Defining_User_Model = 1;
   status = do_define_user_model (fun_body);
   Defining_User_Model = 0;

   if ((status == 0) && (Define_Model_Hook != NULL))
     status = (SLexecute_function (Define_Model_Hook) == -1) ? -1 : 0;

   return status;
}

/*}}}*/

static void define_user_model (int *interactive, char *fun_body) /*{{{*/
{
   if (interactive)
     Use_Interactive_Param_Init = *interactive ? 1 : 0;

   if (*fun_body != 0)
     _define_user_model (fun_body);
   else
     delete_user_model ();

   Use_Interactive_Param_Init = 0;
}

/*}}}*/

int update_user_model (void) /*{{{*/
{
   User_Function_Type *f;
   int ret;
   char *copy;
   char *s;

   if (NULL == (f = get_user_function ()))
     return -1;

   /* We can't use f->fun_string directly because that
    * pointer will be freed before being re-defined.
    * The string literal is a silly hack to enable flux_corr
    * with user-defined kernels (like the grating pileup kernel).
    */

   s = f->fun_string ? f->fun_string : (char *) "null";
   if (NULL == (copy = isis_make_string (s)))
     return -1;
   ret = _define_user_model (copy);
   ISIS_FREE (copy);

   if (ret)
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "updating fit-function definition");

   return ret;
}

/*}}}*/

/*}}}*/

static void set_fit_type (int *response_type, int *fit_to_flux_corrected) /*{{{*/
{
   if (*fit_to_flux_corrected)
     bit_set (Fit_Data_Type, H_FLUX);
   else
     bit_clear (Fit_Data_Type, H_FLUX);

   switch (*response_type)
     {
      case USE_IDEAL_ARF:
      case USE_IDEAL_RMF:
      case USE_IDEAL_ARF_RMF:
      case USE_ASSIGNED_ARF_RMF:
        Response_Type = *response_type;
        break;

      default:
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "response type %d:  response type not set",
                    *response_type);
        return;
     }
}

/*}}}*/

/*{{{ evaluate model function and apply response */

static int pop_model_result (Isis_Hist_t *x) /*{{{*/
{
   SLang_Array_Type *value = NULL;
   int ret = 0;

   if ((SLANG_ARRAY_TYPE == SLang_peek_at_stack ())
       && (0 == SLang_pop_array_of_type (&value, SLANG_DOUBLE_TYPE))
       && (value != NULL)
       && (x->n_notice == (int) value->num_elements))
     {
        memcpy ((char *)x->val, (char *) value->data,
                x->n_notice * sizeof(double));
     }
   else
     {
        SLdo_pop ();
        memset ((char *)x->val, 0, x->n_notice * sizeof(double));
        verbose_warn_hook (NULL, "Failed evaluating user fit function\n");
        if (!Looking_For_Confidence_Limits)
          isis_throw_exception (Isis_Error);
        ret = -1;
     }

   SLang_free_array (value);

   return ret;
}

/*}}}*/

static int eval_model_using_global_grid (Hist_t *h, Isis_Hist_t *g) /*{{{*/
{
   User_Function_Type *f;
   SLang_Name_Type *fun_ptr;

   if (NULL != (fun_ptr = Hist_assigned_model (h)))
     {
        Isis_Arg_Type *args;
        if (NULL != (args = Hist_assigned_model_args (h)))
          {
             (void) isis_push_args (args);
          }
     }
   else
     {
        f = get_user_function ();
        if ((f == NULL) || (f->fun_ptr == NULL))
          {
             isis_vmesg (INTR, I_INFO, __FILE__, __LINE__, "fit function not defined");
             return -1;
          }
        fun_ptr = f->fun_ptr;
     }

   SLexecute_function (fun_ptr);

   return pop_model_result (g);
}

/*}}}*/

static double *temp_workspace_for_eval (Isis_Hist_t *g) /*{{{*/
{
   if (g == NULL)
     return NULL;
   /* Access temporary work space without a malloc. */
   return g->val + g->n_notice;
}

/*}}}*/

static int eval_model_using_cached_grid (Hist_t *h, Isis_Hist_t *g) /*{{{*/
{
   Fit_Data_t *d = Current_Fit_Data_Info;
   Hist_Eval_Grid_Method_Type *egm;
   Cached_Grid_Type *m;
   Isis_Hist_t *x;

   if (NULL == (egm = Hist_eval_grid_method (h)))
     return -1;

   if (NULL == (m = find_cached_grid (d->cache, egm->type, egm->id)))
     return -1;

   if (m->updated_cached_model_values == 0)
     {
        Isis_Hist_t tmp;
        int ret;

        /* FIXME!!! This is ugly. */
        tmp = *g;
        *g = m->grid;
        ret = eval_model_using_global_grid (h, g);
        *g = tmp;
        if (ret == -1) return -1;

        m->updated_cached_model_values = 1;
     }

   /* x is an arbitrary grid with notice flags such that the noticed
    * bin values of x are packed into x->val[0:x->n_notice-1].
    * g is a different arbitrary grid with different notice flags.
    * We must map from x => g such that the noticed bin values
    * of g are packed into g->val[0:g->n_notice-1].
    */

   x = &m->grid;
   return Isis_Hist_rebin_noticed (x, temp_workspace_for_eval(x),
                                   g, temp_workspace_for_eval(g));
}

/*}}}*/

static int evaluate_model (Isis_Hist_t *g) /*{{{*/
{/* Don't change the signature of this function without
  * changing the prototype of kernel->compute_kernel()
  * in isis.h
  */
   Hist_t *h = Hist_Current;
   Hist_Eval_Grid_Method_Type *m;
   int status;

   /* global evaluation grid pointer should now point to
    * the source model grid (usually the arf grid)
    */
   if (NULL == g)
     return -1;

   if (NULL == (m = Hist_eval_grid_method (h)))
     return -1;

   if (m->eval_model == NULL)
     m->eval_model = &eval_model_using_global_grid;

   status = (*m->eval_model)(h, g);

   if (Fit_Store_Model)
     {
        if (-1 == Hist_set_model (h, H_FLUX, g->val))
          return -1;
     }

   return status;
}

/*}}}*/

static int allocate_notice_arrays (Isis_Hist_t *g) /*{{{*/
{
   int i, n = g->nbins;

   g->n_notice = n;

   if ((NULL == (g->notice = (int *) ISIS_MALLOC (n * sizeof(int))))
       || (NULL == (g->notice_list = (int *) ISIS_MALLOC (n * sizeof(int)))))
     return -1;

   for (i = 0; i < n; i++)
     {
        g->notice[i] = 1;
        g->notice_list[i] = i;
     }

   return 0;
}

/*}}}*/

static void free_notice_arrays (Isis_Hist_t *g) /*{{{*/
{
   if (g == NULL)
     return;
   ISIS_FREE (g->notice);
   ISIS_FREE (g->notice_list);
   g->n_notice = 0;
}

/*}}}*/

static int pop_stored_background (SLang_Array_Type **sl_bgd, Hist_t *h) /*{{{*/
{
   double *bgd = NULL;
   int orig_nbins;

   *sl_bgd = NULL;

   if (-1 == Hist_copy_scaled_background (h, &bgd))
     return -1;

   /* its ok if no background data is available. */
   if (bgd == NULL)
     return 0;

   orig_nbins = Hist_orig_hist_size (h);

   *sl_bgd = SLang_create_array (SLANG_DOUBLE_TYPE, 0, bgd, &orig_nbins, 1);
   if (*sl_bgd == NULL)
     return -1;

   return 0;
}

/*}}}*/

int pop_instrumental_background (SLang_Array_Type **bgd, Hist_t *h) /*{{{*/
{
   SLang_Name_Type *hook = Hist_get_instrumental_background_hook (h);
   Isis_Hist_t *g;
   int ret = -1;

   *bgd = NULL;

   /* If available, use background data by default.
    * Probably this isn't rigorously correct, but it
    * its the way most people treat the background.
    */
   if (hook == NULL)
     return pop_stored_background (bgd, h);

   /* g is a pointer to a global structure */
   if (NULL == (g = get_evaluation_grid ()))
     return -1;
   g->val = NULL;

   /* point evaluation grid 'g' at the DATA/channel grid */
   if ((-1 == _Hist_get_orig_hist_grid (h, g))
       || (-1 == allocate_notice_arrays (g)))
     return -1;

   SLexecute_function (hook);

   if ((SLANG_ARRAY_TYPE == SLang_peek_at_stack ())
       && (0 == SLang_pop_array_of_type (bgd, SLANG_DOUBLE_TYPE))
       && (*bgd != NULL)
       && (g->n_notice == (int) (*bgd)->num_elements))
     {
        ret = 0;
     }
   else
     {
        SLdo_pop ();
        verbose_warn_hook (NULL, "Failed evaluating instrument background function\n");
        if (!Looking_For_Confidence_Limits)
          isis_throw_exception (Isis_Error);
        SLang_free_array (*bgd);
        *bgd = NULL;
        ret = -1;
     }

   free_notice_arrays (g);

   return ret;
}

/*}}}*/

static int add_instrumental_background (double *cts, Hist_t *h) /*{{{*/
{
   SLang_Array_Type *bgd = NULL;
   unsigned int i, n;
   double *bd;

   if (0 == is_flux (Fit_Data_Type))
     {
        if (-1 == pop_instrumental_background (&bgd, h))
          return -1;
     }

   if (NULL != Hist_post_model_hook (h))
     {
        int status = Hist_run_post_model_hook (h, cts, bgd);
        SLang_free_array (bgd);
        return status;
     }

   if (bgd == NULL)
     return 0;

   n = bgd->num_elements;
   bd = (double *)bgd->data;

   for (i = 0; i < n; i++)
     cts[i] += bd[i];

   SLang_free_array (bgd);

   return 0;
}

/*}}}*/

static int apply_response (double *result, Isis_Hist_t *g, Hist_t *h) /*{{{*/
{
   Isis_Kernel_Def_t *def;
   Isis_Kernel_t *k;
   double *kp = NULL;
   int hist_index, ret;

   if (NULL == h || NULL == g || NULL == result)
     return -1;

   k = Hist_get_kernel (h);
   if ((NULL == k) || (NULL == k->compute_kernel))
     return -1;

   def = k->kernel_def;
   hist_index = Hist_get_index (h);

   kp = Fit_get_kernel_params (Param, hist_index, def);
   if ((NULL == kp) && (def->num_kernel_parms > 0))
     return -1;

   ret = k->compute_kernel (k, result, g, kp, def->num_kernel_parms,
                            evaluate_model);
   ISIS_FREE (kp);

   return ret;
}

/*}}}*/

static int compute_hist_model (Hist_t *h, double *bincts) /*{{{*/
{
   Isis_Hist_t *g = NULL;
   double *cts = NULL;
   int orig_nbins;
   int ret = -1;

   /* allocate space for the full-resolution result */
   if ((-1 == (orig_nbins = Hist_orig_hist_size (h)))
       || (NULL == (cts = (double *) ISIS_MALLOC (orig_nbins * sizeof(double)))))
     return -1;
   memset ((char *)cts, 0, orig_nbins * sizeof(double));

   /* g is a pointer to a global structure */
   if (NULL == (g = get_evaluation_grid ()))
     return -1;
   g->val = NULL;

   /* point evaluation grid 'g' at the appropriate model grid */
   if (-1 == Hist_get_model_grid (g, h))
     goto finish;

   /* allocate space for the noticed bins in the source model,
    * plus some temporary work space which may or may not be used.
    * (use temp_workspace_for_eval() to get a pointer to this space)
    */

   if (NULL == (g->val = (double *) ISIS_MALLOC ((g->n_notice + g->nbins) * sizeof(double))))
     goto finish;

   if (-1 == apply_response (cts, g, h))
     goto finish;

   ISIS_FREE (g->val);
   memset ((char *)g, 0, sizeof (*g));

   if (-1 == add_instrumental_background (cts, h))
     goto finish;

   if (-1 == Hist_apply_rebin_and_notice_list (bincts, cts, h))
     goto finish;

   if (Fit_Store_Model)
     {
        unsigned int model_type = 0;
        if (is_flux(Fit_Data_Type))
          bit_set (model_type, H_FLUX | H_CONVOLVED);
        if (-1 == Hist_set_model (h, model_type, bincts))
          goto finish;
     }

   ret = 0;
   finish:

   ISIS_FREE (g->val);
   ISIS_FREE (cts);

   return ret;
}

/*}}}*/

static int compute_hist_model_hook (Hist_t *h, void *cl) /*{{{*/
{
   double **model = (double **)cl;
   int ret;

   if (Hist_num_data_noticed (h) < 1)
     return 0;

   ret = compute_hist_model (h, *model);
   *model += Hist_num_data_noticed (h);

   return ret;
}

/*}}}*/

static int copy_hist_model_hook (Hist_t *h, void *cl) /*{{{*/
{
   double **model = (double **)cl;
   unsigned int model_type = 0;
   int status;

   if (Hist_num_data_noticed (h) < 1)
     return 0;

   if (is_flux(Fit_Data_Type))
     {
        bit_set (model_type, H_FLUX | H_CONVOLVED);
     }

   /* Assume the model has already been computed */
   status = Hist_copy_packed_model (h, model_type, *model);
   *model += Hist_num_data_noticed (h);

   return status;
}

/*}}}*/

static int compute_model (double *model, double *par_list, int npars_vary) /*{{{*/
{
   static char hook_name[] = "isis_start_eval_hook";
   double *model_end;
   Fit_Data_t *d = Current_Fit_Data_Info;
   int severity;

   if ((NULL == model) || (d == NULL)
       || (NULL == par_list) || (npars_vary < 0))
     return -1;

   if (Unpack == NULL)
     {
        isis_vmesg (FAIL, I_INTERNAL, __FILE__, __LINE__, "compute_model: param_unpack_method = NULL");
        return -1;
     }

   if (-1 == (*Unpack)(Param, par_list))
     return -1;

   /* Check for user-defined hook function
    * Now that the current trial parameters are unpacked,
    * the user can get the values via get_params()
    */
   if (2 == SLang_is_defined (hook_name))
     (void) SLang_run_hooks (hook_name, 0);

   cached_models_need_updating(d);
   model_end = model + d->nbins;

   if (Computing_Statistic_Only == 0)
     {
        (void) map_datasets (&compute_hist_model_hook, &model);
     }
   else
     {
        (void) map_datasets (&copy_hist_model_hook, &model);
     }

   Num_Statistic_Evaluations++;
   if (model == model_end)
     return 0;

   severity = Looking_For_Confidence_Limits ? FAIL : INTR;

   isis_vmesg (severity, I_FAILED, __FILE__, __LINE__, "evaluating model");
   return -1;
}

/*}}}*/

static int combine_marked_datasets (Fit_Data_t *d, int apply_weights, double *y, double *yc) /*{{{*/
{
   double w = 1.0;
   int i;

   if (d == NULL || y == NULL || yc == NULL)
     return -1;

   memset ((char *)yc, 0, d->nbins_after_datasets_combined * sizeof(double));

   for (i = 0; i < d->num_datasets; i++)
     {
        Hist_t *h = d->datasets[i];
        double *yt, *yf;
        int k, n;

        n = Hist_num_data_noticed (h);

        if (apply_weights)
          {
             (void) Hist_combination_weight (h, &w);
          }

        yf = y + d->offsets[i];
        yt = yc + d->offset_for_marked[i];

        for (k = 0; k < n; k++)
          {
             yt[k] += w * yf[k];
          }
     }

   return 0;
}

/*}}}*/

static int _fitfun (void *cl, double *x,  unsigned int nbins, /*{{{*/
                    double *par, unsigned int npars, double *model)
{
   Fit_Data_t *d = Current_Fit_Data_Info;
   int no_combined_datasets, num_pars = npars;
   (void) x; (void) cl; (void) nbins;

   if (d == NULL)
     {
        isis_vmesg (FAIL, I_INTERNAL, __FILE__, __LINE__, "failed getting fit data info");
        return -1;
     }

   no_combined_datasets = (d->nbins == d->nbins_after_datasets_combined);

   if (no_combined_datasets)
     return compute_model (model, par, num_pars);

   if (-1 == compute_model (d->tmp, par, num_pars))
     return -1;

   return combine_marked_datasets (d, 0, d->tmp, model);
}

/*}}}*/

/*}}}*/

/*{{{ assemble data to fit */

static int prepare_kernel_for_fit (Hist_t *h) /*{{{*/
{
   Isis_Kernel_Def_t *def;
   char *params;

   /* Prevent re-initializing the fit-kernel
    * while the fit is in progress -- e.g. if the
    * user-defined fit function calls flux_corr()
    */
   if (Isis_Fit_In_Progress
       && (NULL != Hist_get_kernel (h)))
     return 0;

   /* Can't use Hist_get_kernel() here because we're
    * defining the kernel type. */
   if (NULL == (def = get_histogram_kernel (Kernel_Table, h)))
     return -1;

   /* ok if params == NULL (e.g. std kernel) */
   params = get_kernel_params (Kernel_Table, h);

   if (-1 == Hist_reallocate_kernel (h, params, def))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "allocating kernel for data set %d",
                    Hist_get_index(h));
        return -1;
     }

   return 0;
}

/*}}}*/

static int prepare_for_fit (Hist_t *h, void *cl) /*{{{*/
{
   (void) cl;

   if (Hist_num_data_noticed (h) < 1)
     return 0;

   if (-1 == Hist_set_fit_responses (h, Fit_Data_Type, Response_Type))
     return -1;

   return prepare_kernel_for_fit (h);
}

/*}}}*/

int sync_model_with_data (void) /*{{{*/
{
   if (-1 == map_datasets (prepare_for_fit, NULL))
     return -1;

   /* By updating the model, we initialize the global
    * parameter table 'Param' */
   if (-1 == update_user_model ())
     return -1;

   /* Sync tied params here to make sure tieing behaves
    * properly with flux_corr() */
   return Fit_sync_tied_params (Param);
}

/*}}}*/

typedef struct
{
   Fit_Data_t *data;
   int nbins;
   int index;
   unsigned int version;
}
Load_Data_Ptr_Type;

static int compute_dataset_combination_offsets (Fit_Data_t *d) /*{{{*/
{
   int i, j, n, gid, npts, offset;

   /* First, assume each dataset goes at the end,
    * then check to see if it might go elsewhere (brute force).
    * We already know that the grids and ignore/notice match
    */

   npts = 0;

   for (i = 0; i < d->num_datasets; i++)
     {
        Hist_t *h = d->datasets[i];

        n = Hist_num_data_noticed (h);
        gid = Hist_combination_id (h);

        offset = npts;

        if (gid != 0)
          {
             for (j = 0; j < i; j++)
               {
                  if (gid == Hist_combination_id (d->datasets[j]))
                    offset = d->offset_for_marked[j];
               }
          }

        d->offset_for_marked[i] = offset;

        if (offset == npts)
          npts += n;
     }

   d->nbins_after_datasets_combined = npts;

   return 0;
}

/*}}}*/

static int load_data_hook (Hist_t *h, void *cl) /*{{{*/
{
   Load_Data_Ptr_Type *ld = (Load_Data_Ptr_Type *)cl;
   Fit_Data_t *d = ld->data;
   double *pdata, *pweight;
   int n;

   n = Hist_num_data_noticed (h);

   if (n < 1)
     return 0;

   d->datasets[ld->index] = h;
   d->offsets[ld->index] = ld->nbins;

   pdata = d->data + ld->nbins;
   pweight = d->weight + ld->nbins;

   if (-1 == Hist_copy_noticed_data (h, ld->version, pdata, pweight))
     return -1;

   ld->index += 1;
   ld->nbins += n;

   if (ld->nbins > d->nbins)
     return -1;

   return 0;
}

/*}}}*/

static int Fit_load_data (Fit_Data_t *d) /*{{{*/
{
   Load_Data_Ptr_Type ld;
   double *wt;
   int i;

   if (NULL == d)
     return -1;

   ld.data = d;
   ld.nbins = 0;
   ld.index = 0;
   ld.version = 0;
   bit_set (ld.version, H_DATA | Fit_Data_Type);

   if (-1 == map_datasets (&load_data_hook, &ld))
     return -1;

   if (ld.nbins != d->nbins)
     return -1;

   if (-1 == compute_dataset_combination_offsets (d))
     return -1;

   if (NULL == (d->tmp = (double *) ISIS_MALLOC (d->nbins * sizeof(double))))
     return -1;

   /* wt != 0 guaranteed elsewhere */
   wt = d->weight;
   for (i=0; i < d->nbins; i++)
     wt[i] = 1.0/(wt[i] * wt[i]);

   return 0;
}

/*}}}*/

static int get_num_noticed_bins (void) /*{{{*/
{
   int nbins;
   Hist_t *head = get_histogram_list_head ();

   if (NULL == head)
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "no data has been loaded");
        return -1;
     }

   nbins = Hist_all_data_noticed (head);
   if (nbins <= 0)
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "all data bins are being ignored");
        return -1;
     }

   return nbins;
}

/*}}}*/

void free_fit_data (Fit_Data_t *d) /*{{{*/
{
   Cached_Grid_Type *t;

   if (d == NULL)
     return;

   ISIS_FREE (d->data);
   ISIS_FREE (d->weight);
   ISIS_FREE (d->datasets);
   ISIS_FREE (d->offsets);
   ISIS_FREE (d->offset_for_marked);
   ISIS_FREE (d->tmp);

   t = d->cache;
   while (t)
     {
        Cached_Grid_Type *next = t->next;
        Isis_Hist_free (&t->grid);
        ISIS_FREE(t);
        t = next;
     }

   ISIS_FREE (d);
}

/*}}}*/

static Fit_Data_t *new_fit_data (int nbins) /*{{{*/
{
   Fit_Data_t *d;

   if (NULL == (d = (Fit_Data_t *) ISIS_MALLOC (sizeof(Fit_Data_t))))
     return NULL;
   memset ((char *)d, 0, sizeof (*d));

   d->nbins = nbins;
   d->num_datasets = Hist_num_noticed_histograms (get_histogram_list_head());
   d->tmp = NULL;

   if (NULL == (d->data = (double *) ISIS_MALLOC (nbins * sizeof(double)))
       || NULL == (d->weight = (double *) ISIS_MALLOC (nbins * sizeof(double)))
       || NULL == (d->datasets = (Hist_t **) ISIS_MALLOC (d->num_datasets * sizeof(Hist_t *)))
       || NULL == (d->offsets = (int *) ISIS_MALLOC (d->num_datasets * sizeof(int)))
       || NULL == (d->offset_for_marked = (int *) ISIS_MALLOC (d->num_datasets * sizeof(int)))
       )
     {
        free_fit_data (d);
        return NULL;
     }

   memset ((char *)d->offset_for_marked, 0, d->num_datasets * sizeof(int));

   return d;
}

/*}}}*/

static int grids_match_in_dataset_combinations (Fit_Data_t *d) /*{{{*/
{
   Hist_t *head = get_histogram_list_head ();
   int i, j, k, n, gid;
   int *ok;

   /* quick return when no datasets are combined */
   if (0 == Hist_any_noticed_dataset_combinations (head))
     return 1;

   if (NULL == (ok = (int *) ISIS_MALLOC (d->num_datasets * sizeof(int))))
     return 0;
   memset ((char *)ok, 0, d->num_datasets * sizeof(int));
   n = 0;

   for (i = 0; i < d->num_datasets; i++)
     {
        Hist_t *hi = d->datasets[i];
        gid = Hist_combination_id (hi);

        /* is this dataset a combination member? */
        if (gid == 0)
          continue;

        for (k = 0; k < n; k++)
          {
             if (gid == ok[k])
               break;
          }

        /* has this combination been checked already? */
        if ((n > 0) && (k < n))
          continue;

        for (j = 0; j < d->num_datasets; j++)
          {
             Hist_t *hj = d->datasets[j];
             if ((gid == Hist_combination_id (hj))
                 && (0 == Hist_grids_are_identical (hi, hj)))
               {
                  ISIS_FREE (ok);
                  return 0;
               }
          }

        ok[n++] = gid;
     }

   ISIS_FREE (ok);
   return 1;
}

/*}}}*/

/* Cached model-eval grids */

static Cached_Grid_Type *find_cached_grid (Cached_Grid_Type *t, unsigned int type, unsigned int id) /*{{{*/
{
   for ( ; t != NULL; t = t->next)
     {
        if ((t->type == type) && (t->id == id))
          break;
     }
   return t;
}

/*}}}*/

static void cached_models_need_updating (Fit_Data_t *d) /*{{{*/
{
   Cached_Grid_Type *t;

   for (t = d->cache; t != NULL; t = t->next)
     {
        t->updated_cached_model_values = 0;
     }
}

/*}}}*/

static int make_merged_eval_grid (Hist_t *h, void *pv) /*{{{*/
{
   Cached_Grid_Type *m = (Cached_Grid_Type *)pv;
   Hist_t *head = get_histogram_list_head ();
   (void) h;
   if (m == NULL)
     return -1;
   return Hist_make_merged_eval_grid (head, m->id, &m->grid);
}

/*}}}*/

static int make_user_eval_grid (Hist_t *h, void *pv) /*{{{*/
{
   Cached_Grid_Type *m = (Cached_Grid_Type *)pv;
   if (m == NULL)
     return -1;
   return Hist_make_user_eval_grid (h, &m->grid);
}

/*}}}*/

static int make_cached_eval_grid (Hist_t *h, void *cl) /*{{{*/
{
   Fit_Data_t *d = (Fit_Data_t *)cl;
   Hist_Eval_Grid_Method_Type *egm;
   Cached_Grid_Type *m;

   if (Hist_num_data_noticed (h) < 1)
     return 0;

   if (NULL == (egm = Hist_eval_grid_method (h)))
     return -1;

   if (egm->make_grid == NULL)
     return 0;

   if (NULL != (m = find_cached_grid (d->cache, egm->type, egm->id)))
     return 0;

   if (NULL == (m = (Cached_Grid_Type *) ISIS_MALLOC(sizeof *m)))
     return -1;

   m->type = egm->type;
   m->id = egm->id;

   if (-1 == (*egm->make_grid)(h, m))
     {
        ISIS_FREE(m);
        return -1;
     }

   m->next = d->cache;
   d->cache = m;

   return 0;
}

/*}}}*/

static int init_cached_eval_grids (Fit_Data_t *d) /*{{{*/
{
   return map_datasets (&make_cached_eval_grid, (void *)d);
}

/*}}}*/

static int init_model_structs (Hist_t *h, void *cl) /*{{{*/
{
   (void) cl;

   if (Hist_num_data_noticed (h) < 1)
     return 0;

   return Hist_init_model_structs (h);
}

/*}}}*/

Fit_Data_t *get_fit_data (void) /*{{{*/
{
   Fit_Data_t *d;
   int nbins;

   nbins = get_num_noticed_bins ();
   if (nbins <= 0)
     return NULL;

   /* FIXME?  sync model and data _before_ loading data,
    * so that scaled background gets the right ARF exposure
    * (it would be really nice to find a way to avoid this
    * tricky/subtle side-effect....)
    */
   if (-1 == sync_model_with_data ())
     return NULL;

   if (-1 == map_datasets (init_model_structs, NULL))
     return NULL;

   if (NULL == (d = new_fit_data (nbins)))
     return NULL;

   if (-1 == Fit_load_data (d))
     {
        free_fit_data (d);
        return NULL;
     }

   if (0 == grids_match_in_dataset_combinations (d))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "combined datasets have mismatched grids");
        free_fit_data (d);
        return NULL;
     }

   if (-1 == init_cached_eval_grids (d))
     {
        free_fit_data (d);
        return NULL;
     }

   return d;
}

/*}}}*/

/*}}}*/

/*{{{ fit verbose hook */

static int Open_Fit_Verbose_Hook;

void init_verbose_hook (void) /*{{{*/
{
   if (!Fit_Verbose) return;
   if (Open_Fit_Verbose_Hook) return;

   SLsig_block_signals ();
   SLang_run_hooks ("open_fit_verbose_hook", 0);
   SLsig_unblock_signals ();

   Open_Fit_Verbose_Hook = 1;
}

/*}}}*/

void deinit_verbose_hook (void) /*{{{*/
{
   if (!Open_Fit_Verbose_Hook) return;

   SLsig_block_signals ();
   SLang_run_hooks ("close_fit_verbose_hook", 0);
   SLsig_unblock_signals ();

   Open_Fit_Verbose_Hook = 0;
}

/*}}}*/

void verbose_warn_hook (void *cl, const char * fmt, ...) /*{{{*/
{
   char buf[BUFSIZE];
   va_list ap;

   (void) cl;

   SLsig_block_signals ();

   if (fmt == NULL)
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   va_start (ap, fmt);
   if (-1 == isis_vsnprintf (buf, sizeof(buf), fmt, ap))
     fprintf (stderr, "**** String buffer overflow in verbose_warn_hook\n");
   va_end (ap);

   finish:
   SLang_run_hooks ("fit_verbose_warn_hook", 1, buf);
   SLsig_unblock_signals ();
}

/*}}}*/

void verbose_info_hook (void *cl, double statistic, /*{{{*/
                        double *par, unsigned int npars_vary)
{
   SLang_Array_Type *sl_param = NULL;
   SLang_Array_Type *sl_param_name = NULL;
   int i, npar = npars_vary;

   (void) cl;

   SLsig_block_signals ();

   sl_param = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &npar, 1);
   if (sl_param == NULL)
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   memcpy ((char *)sl_param->data, (char *)par, npar * sizeof(double));

   sl_param_name = SLang_create_array (SLANG_STRING_TYPE, 1, NULL, &npar, 1);
   if (sl_param == NULL)
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   for (i = 0; i < npar; i++)
     {
        char name[OUT_PARAM_NAME_SIZE];
        char *namep = name;

        if (Array_Fit_Fun)
          sprintf (name, "p[%d]", i);
        else if (-1 == get_variable_param_name (i, name, OUT_PARAM_NAME_SIZE))
          break;

        if (SLang_set_array_element (sl_param_name, &i, &namep))
          break;
     }

   finish:

   SLang_start_arg_list ();
   SLang_push_double (statistic);
   SLang_push_array (sl_param, 1);
   SLang_push_array (sl_param_name, 1);
   SLang_end_arg_list ();

   SLang_execute_function ("fit_verbose_info_hook");

   SLsig_unblock_signals ();
}

/*}}}*/

/*}}}*/

static SLang_Array_Type *make_sl_darray_copy (double *x, int n) /*{{{*/
{
   SLang_Array_Type *s;

   if (x == NULL)
     return NULL;

   if (NULL == (s = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
     return NULL;

   memcpy ((char *)s->data, (char *)x, n * sizeof(double));

   return s;
}

/*}}}*/

void set_hook_from_stack (SLang_Name_Type **hook) /*{{{*/
{
   SLang_free_function (*hook);

   if (SLANG_REF_TYPE == SLang_peek_at_stack ())
     {
        *hook = SLang_pop_function ();
     }
   else if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        (void) SLang_pop_null();
        *hook = NULL;
     }
   else
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "expected a function reference or NULL");
        isis_throw_exception (Isis_Error);
     }
}

/*}}}*/

/*{{{ alternate fit methods and statistics */

static int load_fit_method (char *file, char *name) /*{{{*/
{
   return isis_fit_load_fit_routine (file, name, NULL);
}

/*}}}*/

static int load_fit_statistic (char *file, char *sname) /*{{{*/
{
   return isis_fit_load_statistic (file, sname);
}

/*}}}*/

static void set_fit_constraint_fun (void) /*{{{*/
{
   set_hook_from_stack (&Fit_Constraint);
   (void) update_user_model ();
}

/*}}}*/

static void set_fit_method_name (char *name) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   char *str;

   if (NULL == (e = isis_find_fit_engine (name)))
     return;

   if (e->option_string)
     name = e->option_string;

   if (NULL == (str = isis_make_string (name)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   ISIS_FREE (Fit_Method);
   Fit_Method = str;
}

/*}}}*/

static void set_fit_statistic_name (char *name) /*{{{*/
{
   Isis_Fit_Statistic_Type *s;
   char *str;

   if (NULL == (s = isis_find_fit_statistic (name)))
     return;

   if (s->option_string)
     name = s->option_string;

   if (NULL == (str = isis_make_string (name)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   ISIS_FREE (Fit_Statistic);
   Fit_Statistic = str;
}

/*}}}*/

Isis_Fit_Statistic_Type *current_fit_statistic (void) /*{{{*/
{
   Isis_Fit_Statistic_Type *s;
   char *sname = Fit_Statistic;

   if (NULL == (s = isis_find_fit_statistic (sname)))
     {
        fprintf (stderr, "statistic %s does not exist\n",
                 sname ? sname : "<null>");
        return NULL;
     }

   return s;
}

/*}}}*/

Isis_Fit_Engine_Type *current_fit_method (void) /*{{{*/
{
   Isis_Fit_Engine_Type *e;
   char *name = Fit_Method;

   if (NULL == (e = isis_find_fit_engine (name)))
     {
        fprintf (stderr, "method %s does not exist\n",
                 name ? name : "<null>");
        return NULL;
     }

   return e;
}

/*}}}*/

static void get_fit_method_name (void) /*{{{*/
{
   Isis_Fit_Engine_Type *e;

   if (NULL == (e = current_fit_method ()))
     return;

   SLang_push_string (e->option_string);
}

/*}}}*/

static void get_fit_statistic_name (void) /*{{{*/
{
   Isis_Fit_Statistic_Type *s;

   if (NULL == (s = current_fit_statistic()))
     return;

   SLang_push_string (s->option_string);
}

/*}}}*/

/*}}}*/

static int sl_fit_range_hook (void *cl, double *par_min, double *par_max, double *par, unsigned int npars) /*{{{*/
{
   SLang_Array_Type *sl_pmin, *sl_pmax, *sl_par, *sl_idx;
   int *idx = (int *)cl;
   int num_pars = npars;

   sl_pmin = sl_pmax = sl_par = sl_idx = NULL;

   if (idx != NULL)
     {
        sl_idx = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &num_pars, 1);
        if (sl_idx == NULL)
          return -1;
        memcpy ((char *)sl_idx->data, (char *)idx, num_pars * sizeof(int));
     }

   sl_pmin = make_sl_darray_copy (par_min, num_pars);
   sl_pmax = make_sl_darray_copy (par_max, num_pars);
   sl_par = make_sl_darray_copy (par, num_pars);

   if ((sl_pmin == NULL) || (sl_pmax == NULL) || (sl_par == NULL))
     {
        SLang_free_array (sl_pmin);
        SLang_free_array (sl_pmax);
        SLang_free_array (sl_par);
        SLang_free_array (sl_idx);
        return -1;
     }

   SLang_start_arg_list ();
   if (-1 == SLang_push_array (sl_par, 1)
       ||-1 == SLang_push_array (sl_pmin, 1)
       ||-1 == SLang_push_array (sl_pmax, 1)
       ||-1 == SLang_push_array (sl_idx, 1))
     {
        SLang_free_array (sl_pmin);
        SLang_free_array (sl_pmax);
        SLang_free_array (sl_par);
        SLang_free_array (sl_idx);
        return -1;
     }
   SLang_end_arg_list ();

   /* => fit_range_hook (par, min, max, idx) */
   SLexecute_function (Fit_Range_Hook);

   /* (par, min, max) <= */
   if (-1 == Isis_pop_double_array (par_max, num_pars)
       ||-1 == Isis_pop_double_array (par_min, num_pars)
       || -1 == Isis_pop_double_array (par, num_pars))
     return -1;

   return 0;
}

/*}}}*/

static void set_fit_range_hook (void) /*{{{*/
{
   set_hook_from_stack (&Fit_Range_Hook);
}

/*}}}*/

static void set_define_model_hook_intrin (void)
{
   set_hook_from_stack (&Define_Model_Hook);
}

static int set_fit_method_hooks (Isis_Fit_Type *f) /*{{{*/
{
   if (Fit_Range_Hook)
     isis_fit_set_range_hook (f, sl_fit_range_hook);
   else
     isis_fit_set_range_hook (f, NULL);

   /* FIXME -- this is ugly */
   if (0 == isis_strcasecmp (Fit_Method, Fit_Default_Fit_Method))
     {
        (void) isis_fit_set_warn_hook (f, verbose_warn_hook);
        (void) isis_fit_set_verbose_hook (f, verbose_info_hook);
     }

   return 0;
}

/*}}}*/

static void set_slangfun_param_default_hook (char *fun_name) /*{{{*/
{
   Fit_Fun_t *ff;
   int fun_type;

   if ((-1 == (fun_type = Fit_get_fun_type (fun_name)))
       ||(NULL == (ff = Fit_get_fit_fun (fun_type))))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "setting parameter default hook: %s",
                    fun_name ? fun_name : "<null>");
        isis_throw_exception (Isis_Error);
        return;
     }

   set_hook_from_stack (&ff->slangfun_param_default);

   isis_free_args (ff->slangfun_param_default_args);
   ff->slangfun_param_default_args = isis_pop_args (SLang_Num_Function_Args-2);
}

/*}}}*/

static int combine_datasets (Fit_Data_t *d, double *y, double *w, int n, /*{{{*/
                             double *yc, double *wg)
{
   double *tmp = d->tmp;
   int i, ng;

   ng = d->nbins_after_datasets_combined;

   if (-1 == combine_marked_datasets (d, 1, y, yc))
     return -1;

   /* hack for adding errors in quadrature */

   for (i = 0; i < n; i++)
     tmp[i] = 1.0/w[i];

   if (-1 == combine_marked_datasets (d, 1, tmp, wg))
     return -1;

   for (i = 0; i < ng; i++)
     wg[i] = 1.0/wg[i];

   return 0;
}

/*}}}*/

static void free_fit_object_data (Fit_Object_Data_Type *dt) /*{{{*/
{
   if (dt == NULL)
     return;

   if (dt->malloced)
     {
        ISIS_FREE(dt->data);
        ISIS_FREE(dt->weight);
     }
}

/*}}}*/

static Fit_Object_Data_Type *setup_fit_object_data (Fit_Data_t *d) /*{{{*/
{
   Fit_Object_Data_Type *dt = NULL;

   if (d == NULL)
     return NULL;

   if (NULL == (dt = (Fit_Object_Data_Type *) ISIS_MALLOC (sizeof *dt)))
     return NULL;

   if (0 == Hist_any_noticed_dataset_combinations (get_histogram_list_head ()))
     {
        dt->malloced = 0;
        dt->num = d->nbins;
        dt->data = d->data;
        dt->weight = d->weight;
     }
   else
     {
        dt->malloced = 1;
        dt->num = d->nbins_after_datasets_combined;

        if ((NULL == (dt->data = (double *) ISIS_MALLOC (dt->num * sizeof(double))))
            || (NULL == (dt->weight = (double *) ISIS_MALLOC (dt->num * sizeof(double))))
            || (-1 == combine_datasets (d, d->data, d->weight, d->nbins,
                                        dt->data, dt->weight)))
          {
             free_fit_object_data (dt);
             ISIS_FREE(dt);
             return NULL;
          }
     }

   return dt;
}

/*}}}*/

static void fit_object_close (Fit_Object_Type *fo) /*{{{*/
{
   if (fo == NULL)
     return;

   deinit_verbose_hook ();
   (void) map_datasets (post_fit, NULL);

   Current_Fit_Data_Info = NULL;
   isis_fit_close_fit (fo->ft);
   free_fit_data (fo->d);
   free_fit_object_data (fo->dt);
   if (fo->info) free_fit_param_type (fo->info->par);
   ISIS_FREE (fo->dt);
   ISIS_FREE (fo->info);
   ISIS_FREE (fo);
}

/*}}}*/

int fit_object_config (Fit_Object_Type *fo, Param_t *pt, int unpack_variable) /*{{{*/
{
   int (*pack)(Param_t *, Fit_Param_t *);
   Fit_Info_Type *info;

   if ((fo == NULL) || (fo->info == NULL))
     return -1;

   info = fo->info;

   if (unpack_variable)
     {
        pack = Fit_pack_variable_params;
        Unpack = Fit_unpack_variable_params;
     }
   else
     {
        pack = Fit_pack_all_params;
        Unpack = Fit_unpack_all_params;
     }

   fo->unpack = Unpack;

   if (-1 == (*pack)(pt, info->par))
     return -1;

   isis_fit_set_verbose_level (fo->ft, info->verbose);
   isis_fit_set_ranges (fo->ft, info->par->par_min, info->par->par_max);
   isis_fit_set_param_step (fo->ft, info->par->step);

   return 0;
}

/*}}}*/

static int fit_object_keep_params (Fit_Object_Type *fo, Param_t *pt) /*{{{*/
{
   if (fo == NULL)
     return -1;

   if (fo->unpack == NULL)
     return -1;

   if (fo->info == NULL || fo->info->par == NULL)
     return -1;

   return fo->unpack (pt, fo->info->par->par);
}

/*}}}*/

static Fit_Object_Type *fit_object_open (void) /*{{{*/
{
   Fit_Object_Type *fo = NULL;
   Fit_Info_Type *info;

   (void) map_datasets (pre_fit, NULL);
   init_verbose_hook ();

   SLang_run_hooks ("isis_prefit_hook", 0);

   if (NULL == (fo = (Fit_Object_Type *) ISIS_MALLOC (sizeof *fo)))
     return NULL;
   memset ((char *)fo, 0, sizeof *fo);

   if (NULL == (fo->info = (Fit_Info_Type *) ISIS_MALLOC (sizeof (Fit_Info_Type))))
     goto return_error;
   memset ((char *)fo->info, 0, sizeof (Fit_Info_Type));
   info = fo->info;

   if (-1 == Fit_count_params (Param, &info->num_pars, &info->num_vary))
     goto return_error;

   if (NULL == (info->par = new_fit_param_type (info->num_pars)))
     goto return_error;

   info->verbose = Fit_Verbose;
   info->cl_verbose = NULL;

   if (NULL == (fo->d = get_fit_data ()))
     goto return_error;
   Current_Fit_Data_Info = fo->d;

   if (NULL == (fo->dt = setup_fit_object_data (fo->d)))
     goto return_error;

   if (NULL == (fo->ft = isis_fit_open_fit (Fit_Method, Fit_Statistic, _fitfun, Fit_Constraint)))
     goto return_error;

   set_fit_method_hooks (fo->ft);

   return fo;
return_error:
   fit_object_close (fo);
   return NULL;
}

/*}}}*/

static int eval_fit_stat2 (Isis_Fit_Statistic_Type *s,  int enable_copying, /*{{{*/
                           double *y, double *w, int n, double *par, int npars,
                           double *stat, double **vec, int which)
{
   double *fx=NULL, *vstat=NULL;
   int ret;

   *stat = 0;

   if (NULL == (fx = (double *) ISIS_MALLOC (2*n * sizeof(double))))
     return -1;
   vstat = fx + n;

   Fit_Store_Model = enable_copying;

   if (-1 == _fitfun (NULL, NULL, n, par, npars, fx))
     {
        verbose_warn_hook (NULL, "Failed evaluating fit-function\n");
        ISIS_FREE (fx);
        return -1;
     }

   ret = s->compute_statistic (s, y, fx, w, n, vstat, stat);

   if (vec)
     {
        double *cpy;
        int i;

        if (NULL == (cpy = (double *) ISIS_MALLOC (n * sizeof(double))))
          {
             ISIS_FREE(fx);
             return -1;
          }

        if (which == 0)
          {
             /* vector statistic */
             memcpy ((char *)cpy, (char *)vstat, n * sizeof(double));
          }
        else if (which == 1)
          {
             /* fit residuals */
             for (i = 0; i < n; i++)
               cpy[i] = y[i] - fx[i];
          }
        else
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                         "unrecognized option which=%d", which);
          }
        *vec = cpy;
     }

   ISIS_FREE (fx);

   return ret;
}

/*}}}*/

static int eval_fit_stat (Isis_Fit_Statistic_Type *s,  int enable_copying, /*{{{*/
                          double *y, double *w, int n,
                          double *par, int npars, double *stat)
{
   return eval_fit_stat2 (s, enable_copying, y, w, n, par, npars, stat, NULL, 0);
}

/*}}}*/

static double *eval_fit_residual (Isis_Fit_Statistic_Type *s, int enable_copying, /*{{{*/
                                  double *y, double *w, int n,
                                  double *par, int npars)
{
   double *resid = NULL;
   double stat;

   if (-1 == eval_fit_stat2 (s, enable_copying, y, w, n, par, npars, &stat, &resid, 1))
     {
        ISIS_FREE(resid);
     }

   return resid;
}

/*}}}*/

static SLang_MMT_Type *create_fit_object_mmt_type (Fit_Object_Type *fo);

int fit_statistic (Fit_Object_Type *fo, int optimize, double *stat, int *num_bins) /*{{{*/
{
   Isis_Fit_Type *ft = NULL;
   Fit_Info_Type *info = NULL;
   Fit_Object_Data_Type *dt = NULL;
   Fit_Param_t *par = NULL;
   int say_results;
   int status = -1;
   int fit_ret = 0;

   if ((NULL == fo) || (stat == NULL))
     return -1;

   *stat = 0.0;
   if (num_bins)
     *num_bins = 0;

   info = fo->info;
   dt = fo->dt;
   ft = fo->ft;

   par = info->par;

   if (optimize)
     {
        /* FIXME!!! this is an ugly hack... */
        int slang_optimizer = is_slang_optimizer (ft);
        if (slang_optimizer)
          set_slopt_fit_object (create_fit_object_mmt_type (fo));

        /* disable model copying during the fit */
        Fit_Store_Model = 0;
        fit_ret = isis_fit_perform_fit (ft, par->idx, NULL, dt->data, dt->weight, dt->num,
                                        par->par, par->npars, stat);

        if (slang_optimizer)
          set_slopt_fit_object (NULL);

        if ((fit_ret != 0)
            && (Looking_For_Confidence_Limits == 0))
          {
             verbose_warn_hook (NULL, "warning: minimization failed\n");
          }
     }

   if (-1 == eval_fit_stat (ft->stat, 1, dt->data, dt->weight, dt->num, par->par, par->npars, stat))
     goto return_error;

   if (num_bins)
     *num_bins = dt->num;

   say_results = ((Fit_Verbose >= 0) && (Isis_Verbose >= WARN)
                  && (Looking_For_Confidence_Limits == 0));

   if (say_results)
     {
        fprintf (stdout, " Parameters[Variable] = %d[%d]\n", Num_Params, info->num_vary);
        fprintf (stdout, "            Data bins = %d\n", dt->num);
        isis_fit_report_statistic (ft, stdout, *stat, dt->num, info->num_vary);
     }

   status = 0;
   return_error:

   if (status < 0)
     return -1;

   return (fit_ret ? 1 : 0);
}

/*}}}*/

/* fit data supplied directly as slang arrays */

static SLang_Array_Type *X_Array;
static int eval_array_fit_fun (void *cl, double *x, unsigned int nbins,  /*{{{*/
                               double *par, unsigned int npars, double *model)
{
   SLang_Array_Type *sl_pars = NULL;
   SLang_Array_Type *sl_model = NULL;

   (void) cl;

   if (x == NULL || par == NULL || model == NULL
       || Array_Fit_Fun == NULL || X_Array == NULL)
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   sl_pars = make_sl_darray_copy (par, npars);
   if (sl_pars == NULL)
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_array (X_Array, 0);
   SLang_push_array (sl_pars, 1);
   SLang_end_arg_list ();

   SLexecute_function (Array_Fit_Fun);

   if (-1 == SLang_pop_array_of_type (&sl_model, SLANG_DOUBLE_TYPE)
       || (sl_model == NULL)
       || (sl_model->num_elements != (unsigned int) nbins))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "evaluating fit function");
        SLang_free_array (sl_model);
        return -1;
     }

   memcpy ((char *)model, (char *)sl_model->data, nbins * sizeof(double));
   SLang_free_array (sl_model);

   return 0;
}

/*}}}*/

static int pop3_double_arrays (SLang_Array_Type **sa, SLang_Array_Type **sb, /*{{{*/
                               SLang_Array_Type **sc)
{
   if (sa == NULL || sb == NULL || sc == NULL)
     return -1;

   *sa = NULL;
   *sb = NULL;
   *sc = NULL;

   if (-1 == SLang_pop_array_of_type (sc, SLANG_DOUBLE_TYPE)
       || *sc == NULL
       || -1 == SLang_pop_array_of_type (sb, SLANG_DOUBLE_TYPE)
       || *sb == NULL
       || -1 == SLang_pop_array_of_type (sa, SLANG_DOUBLE_TYPE)
       || *sa == NULL)
     {
        SLang_free_array (*sa);
        SLang_free_array (*sb);
        SLang_free_array (*sc);
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return 0;
}

/*}}}*/

static void array_fit (void) /*{{{*/
{
#define DDATA(x) ((double *)x->data)

   Isis_Fit_Type *f;
   SLang_Array_Type *par, *par_min, *par_max;
   SLang_Array_Type *x, *y, *wt;
   double *par_step = NULL;
   double statistic = 0.0;
   int i, nx, npars, ret = -1;

   f = NULL;
   x = y = wt = par = par_min = par_max = NULL;

   Array_Fit_Fun = SLang_pop_function();
   if (Array_Fit_Fun == NULL)
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (-1 == pop3_double_arrays (&par, &par_min, &par_max)
       || -1 == pop3_double_arrays (&x, &y, &wt))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   X_Array = x;
   nx = x->num_elements;
   npars = par->num_elements;

   if (npars == 0 || nx == 0 || nx < npars)
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (NULL == (f = isis_fit_open_fit (Fit_Method, Fit_Statistic, eval_array_fit_fun,
                                       Fit_Constraint)))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   /* FIXME !! */
   if (is_slang_optimizer (f))
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__,
                    "array_fit does not support using optimizers implemented in S-Lang");
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (NULL == (par_step = ISIS_MALLOC(npars * sizeof(double))))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }
   for (i = 0; i < npars; i++)
     {
        double value;
        (void) SLang_get_array_element (par, &i, &value);
        par_step[i] = 0.01 * ((value == 0) ? 1.0 : fabs(value));
     }

   isis_fit_set_verbose_level (f, Fit_Verbose);
   isis_fit_set_ranges (f, DDATA(par_min), DDATA(par_max));
   isis_fit_set_param_step (f, par_step);

   ret = isis_fit_perform_fit (f, NULL, DDATA(x), DDATA(y), DDATA(wt), nx,
                               DDATA(par), npars, &statistic);
   if (ret)
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "fitting array data");

   finish:

   SLang_free_function (Array_Fit_Fun);
   Array_Fit_Fun = NULL;

   SLang_push_array (par, 1);
   SLang_push_double (statistic);

   isis_fit_close_fit (f);
   X_Array = NULL;
   SLang_free_array (x);
   SLang_free_array (y);
   SLang_free_array (wt);
   SLang_free_array (par_min);
   SLang_free_array (par_max);
   ISIS_FREE(par_step);
#undef DDATA
}

/*}}}*/

/*{{{ high-level eval-model or fit-model */

void free_fit_param_type (Fit_Param_t *p) /*{{{*/
{
   if (p == NULL) return;

   ISIS_FREE (p->par);
   ISIS_FREE (p->idx);
   ISIS_FREE (p);
}

/*}}}*/

Fit_Param_t *new_fit_param_type (unsigned int num) /*{{{*/
{
   Fit_Param_t *p = NULL;

   if (NULL == (p = (Fit_Param_t *) ISIS_MALLOC (sizeof(Fit_Param_t))))
     return NULL;
   memset ((char *)p, 0, sizeof (*p));

   p->par = (double *) ISIS_MALLOC (num*4 * sizeof(double));
   p->idx = (int *) ISIS_MALLOC (num * sizeof(int));
   if ((p->par == NULL) || (p->idx == NULL))
     {
        free_fit_param_type (p);
        p = NULL;
     }
   memset ((char *) p->par, 0, 4*num * sizeof(double));
   memset ((char *) p->idx, 0, num * sizeof(int));

   p->par_min = p->par + num;
   p->par_max = p->par + num*2;
   p->step   = p->par + num*3;
   p->npars = 0;

   return p;
}

/*}}}*/

typedef struct
{
   Fit_Object_Type *fo;
}
Fit_Object_MMT_Type;
static int Fit_Object_MMT_Type_Id = -1;

static void destroy_fit_object_mmt_type (SLtype type, VOID_STAR f) /*{{{*/
{
   Fit_Object_MMT_Type *mt = (Fit_Object_MMT_Type *) f;
   (void) type;
   if (mt != NULL)
     fit_object_close (mt->fo);
   SLfree ((char *)f);
}

/*}}}*/

static SLang_MMT_Type *create_fit_object_mmt_type (Fit_Object_Type *fo) /*{{{*/
{
   SLang_MMT_Type *mmt;
   Fit_Object_MMT_Type *mt;

   if (fo == NULL)
     return NULL;

   if (NULL == (mt = (Fit_Object_MMT_Type *)SLmalloc (sizeof *mt)))
     return NULL;

   mt->fo = fo;

   mmt = SLang_create_mmt (Fit_Object_MMT_Type_Id, (VOID_STAR) mt);
   if (NULL == mmt)
     {
        SLfree ((char *)mt);
        return NULL;
     }

   return mmt;
}

/*}}}*/

static void open_fit_object_mmt_intrin (void) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   Fit_Object_Type *fo = NULL;

   if (NULL == (fo = fit_object_open ()))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == fit_object_config (fo, Param, 1))
     {
        isis_throw_exception (Isis_Error);
        fit_object_close (fo);
        return;
     }

   if (NULL == (mmt = create_fit_object_mmt_type (fo)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == SLang_push_mmt (mmt))
     SLang_free_mmt (mmt);
}

/*}}}*/

static int pop_params (Isis_Fit_Type *ft, Fit_Param_t *par) /*{{{*/
{
   SLang_Array_Type *sl_pars = NULL;

   if ((-1 == SLang_pop_array_of_type (&sl_pars, SLANG_DOUBLE_TYPE))
       || (sl_pars == NULL))
     return -1;

   if (sl_pars->num_elements != (unsigned int) par->npars)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "expecting %d parameter values, got %d",
                    par->npars, sl_pars->num_elements);
        SLang_free_array (sl_pars);
        return -1;
     }

   if (isis_invalid_params (ft->engine, (double *)sl_pars->data, par->npars))
     {
        SLang_free_array (sl_pars);
        return -1;
     }

   memcpy ((char *)par->par, (char *)sl_pars->data, par->npars * sizeof(double));

   SLang_free_array (sl_pars);

   return 0;
}

/*}}}*/

static void fobj_eval_statistic (Fit_Object_MMT_Type *mmt, int *enable_copying) /*{{{*/
{
   Fit_Object_Type *fo = mmt->fo;
   Fit_Info_Type *info = fo->info;
   Isis_Fit_Type *ft = fo->ft;
   Fit_Object_Data_Type *dt = fo->dt;
   Fit_Param_t *par = info->par;
   double stat = DBL_MAX;
   int status = -1;

   if (-1 == pop_params (ft, par))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   status = eval_fit_stat (ft->stat, *enable_copying, dt->data, dt->weight, dt->num, par->par, par->npars, &stat);
   ft->statistic = stat;

   SLang_push_integer (status ? -1 : 0);
   SLang_push_double (stat);
   SLang_push_integer ((int) info->num_vary);
   SLang_push_integer (dt->num);
}

/*}}}*/

static void fobj_eval_residuals (Fit_Object_MMT_Type *mmt, int *enable_copying) /*{{{*/
{
   Fit_Object_Type *fo = mmt->fo;
   Fit_Info_Type *info = fo->info;
   Isis_Fit_Type *ft = fo->ft;
   Fit_Object_Data_Type *dt = fo->dt;
   Fit_Param_t *par = info->par;
   SLang_Array_Type *sl_res = NULL;
   double *res = NULL;

   if (-1 == pop_params (ft, par))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   res = eval_fit_residual (ft->stat, *enable_copying, dt->data, dt->weight, dt->num, par->par, par->npars);

   if ((res == NULL)
       || (NULL == (sl_res = SLang_create_array (SLANG_DOUBLE_TYPE, 0, res, &dt->num, 1))))
     {
        ISIS_FREE(res);
     }

   SLang_push_array (sl_res, 1);
}

/*}}}*/

static void fobj_get_data_weights (Fit_Object_MMT_Type *mmt) /*{{{*/
{
   Fit_Object_Type *fo = mmt->fo;
   Fit_Object_Data_Type *dt = fo->dt;
   SLang_Array_Type *sl_wt = NULL;

   if (NULL == (sl_wt = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &dt->num, 1)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   memcpy ((char *)sl_wt->data, (char *)dt->weight, dt->num * sizeof(double));

   SLang_push_array (sl_wt, 1);
}

/*}}}*/

typedef struct
{
   double stat;
   SLang_Array_Type *covar;
   unsigned int num_vary;
   int num_bins;
}
Fit_Result_Type;

static SLang_CStruct_Field_Type Fit_Result_Type_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Fit_Result_Type, stat, "statistic", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Fit_Result_Type, covar, "covariance_matrix", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Fit_Result_Type, num_vary, "num_variable_params", SLANG_UINT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Fit_Result_Type, num_bins, "num_bins", SLANG_INT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int copy_covariance_matrix (Fit_Result_Type *fr, Fit_Object_Type *fo) /*{{{*/
{
   Isis_Fit_Type *ift = fo->ft;
   int dims[2] = {0, 0};

   if (ift->covariance_matrix)
     {
        dims[0] = fr->num_vary;
        dims[1] = fr->num_vary;
     }

   if (NULL == (fr->covar = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, dims, 2)))
     return -1;

   if (ift->covariance_matrix)
     {
        int i, j, k = 0, pos[2];
        for (i = 0; i < dims[0]; i++)
          {
             pos[0] = i;
             for (j = 0; j < dims[1]; j++)
               {
                  double m = ift->covariance_matrix[k];
                  pos[1] = j;
                  SLang_set_array_element (fr->covar, pos, &m);
                  k++;
               }
          }
     }

   return 0;
}

/*}}}*/

static int iterate_fit_fun (int optimize, SLang_Ref_Type *ref) /*{{{*/
{
   Fit_Result_Type fr;
   Fit_Object_Type *fo = NULL;
   unsigned int num_pars = 0;
   int status = -1;

   memset ((char *)&fr, 0, sizeof(fr));
   fr.stat = isis_nan();

   Num_Statistic_Evaluations = 0;

   if ((NULL != (fo = fit_object_open ()))
       && (0 == fit_object_config (fo, Param, optimize)))
     {
        status = fit_statistic (fo, optimize, &fr.stat, &fr.num_bins);
        fit_object_keep_params (fo, Param);
     }

   if ((fo != NULL) && (fo->info != NULL))
     {
        fr.num_vary = fo->info->num_vary;
     }
   else Fit_count_params (Param, &num_pars, &fr.num_vary);

   copy_covariance_matrix (&fr, fo);
   fit_object_close (fo);

   if (ref != NULL)
     (void) SLang_assign_cstruct_to_ref (ref, &fr, Fit_Result_Type_Layout);

   SLang_free_cstruct (&fr, Fit_Result_Type_Layout);

   return status ? -1 : 0;
}

/*}}}*/

static int find_best_fit (SLang_Ref_Type *ref) /*{{{*/
{
   Computing_Statistic_Only = 0;
   return iterate_fit_fun (1, ref);
}
/*}}}*/

static int eval_model (SLang_Ref_Type *ref) /*{{{*/
{
   Computing_Statistic_Only = 0;
   return iterate_fit_fun (0, ref);
}

/*}}}*/

static int eval_statistic_only (SLang_Ref_Type *ref) /*{{{*/
{
   int status;
   Computing_Statistic_Only = 1;
   status = iterate_fit_fun (0, ref);
   Computing_Statistic_Only = 0;
   return status;
}

/*}}}*/

static int delta_stat_is_chisqr_distributed (void) /*{{{*/
{
   Isis_Fit_Type *f;
   int flag;

   if (NULL == (f = isis_fit_open_fit (Fit_Method, Fit_Statistic, _fitfun,
                                       Fit_Constraint)))
     return 0;

   flag = isis_delta_stat_is_chisqr_distrib (f);
   isis_fit_close_fit (f);

   return flag;
}

/*}}}*/

static void _set_conf_limit_search (int *flag) /*{{{*/
{
   Looking_For_Confidence_Limits = *flag;
}

/*}}}*/

static void confidence_limits (int *idx, double *delta_chisqr, int *verbose, double *tolerance) /*{{{*/
{
   Fit_Object_Type *fo = NULL;
   Param_Info_t *p;
   Isis_Fit_CLC_Type c;
   double conf_min, conf_max;

   if (0 == delta_stat_is_chisqr_distributed())
     {
        isis_vmesg (INTR, I_INFO, __FILE__, __LINE__, "statistic = '%s': confidence limits not supported",
                    Fit_Statistic);
        return;
     }

   c.delta_stat = *delta_chisqr;

   if (NULL == (p = Fit_param_info (Param, *idx)))
     {
        isis_vmesg (INTR, I_INVALID, __FILE__, __LINE__, "parameter %d", *idx);
        return;
     }
   if (p->freeze != 0)
     {
        isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "parameter %d is frozen", *idx);
        return;
     }

   c.tol = *tolerance;
   c.verbose = *verbose;

   if (NULL == (fo = fit_object_open ()))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "initializing fit engine");
        return;
     }

   Looking_For_Confidence_Limits = 1;
   (void) get_confidence_limits (fo, Param, &c, *idx, &conf_min, &conf_max);
   Looking_For_Confidence_Limits = 0;

   fit_object_close (fo);

   SLang_push_double (conf_min);
   SLang_push_double (conf_max);
}

/*}}}*/

/*}}}*/

/* etc */

static void get_instrumental_background (int *hist_index) /*{{{*/
{
   SLang_Array_Type *bgd = NULL;
   Hist_t *h;
   double *rebin_bgd;
   int nbins;

   if (NULL == (h = find_hist (*hist_index)))
     {
        SLang_push_array (bgd, 1);
        return;
     }

   if (-1 == pop_instrumental_background (&bgd, h))
     {
        isis_vmesg (INTR, I_INFO, __FILE__, __LINE__, "no background for dataset %d", *hist_index);
        SLang_push_array (bgd, 1);
        return;
     }

   if (bgd == NULL)
     {
        SLang_push_array (bgd, 1);
        return;
     }

   if (-1 == Hist_apply_rebin ((double *)bgd->data, h, &rebin_bgd, &nbins))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "rebinning background for dataset %d", *hist_index);
        SLang_push_array (bgd, 1);
        return;
     }

   SLang_free_array (bgd);
   bgd = SLang_create_array (SLANG_DOUBLE_TYPE, 0, rebin_bgd, &nbins, 1);

   SLang_push_array (bgd, 1);
}

/*}}}*/

/* Fit_Fun_t MMT */

typedef struct
{
   Fit_Fun_t *ff;
}
Fit_Fun_MMT_Type;
static int Fit_Fun_Type_Id = -1;

static SLang_MMT_Type *create_mmt_fitfun_type (Fit_Fun_t *ff) /*{{{*/
{
   SLang_MMT_Type *mmt;
   Fit_Fun_MMT_Type *mt;

   if (ff == NULL)
     return NULL;

   if (NULL == (mt = (Fit_Fun_MMT_Type *)SLmalloc (sizeof *mt)))
     return NULL;

   mt->ff = ff;

   mmt = SLang_create_mmt (Fit_Fun_Type_Id, (VOID_STAR) mt);
   if (NULL == mmt)
     {
        SLfree ((char *)mt);
        return NULL;
     }

   return mmt;
}

/*}}}*/

static void push_mmt_fitfun_type_intrin (char *fname) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   Fit_Fun_t *ff;
   int fun_type;

   if ((-1 == (fun_type = Fit_get_fun_type (fname)))
       || (NULL == (ff = Fit_get_fit_fun (fun_type))))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (NULL == (mmt = create_mmt_fitfun_type (ff)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == SLang_push_mmt (mmt))
     SLang_free_mmt (mmt);
}

/*}}}*/

static void eval_fitfun_using_handle_intrin (Fit_Fun_MMT_Type *mmt) /*{{{*/
{
   SLang_Array_Type *sl_par=NULL;
   Isis_Hist_t g;
   Fit_Fun_t *ff;

   ff = mmt->ff;

   memset ((char *)&g, 0, sizeof g);
   if (-1 == Isis_Hist_pop_valid_grid (&g))
     {
        Isis_Hist_free (&g);
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == SLang_pop_array_of_type (&sl_par, SLANG_DOUBLE_TYPE))
     {
        Isis_Hist_free (&g);
        isis_throw_exception (Isis_Error);
        return;
     }

   if (sl_par->num_elements != ff->nparams)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__,
                    "%s requires %d parameters but %d were provided",
                    ff->name[0], ff->nparams, sl_par->num_elements);
        SLang_free_array (sl_par);
        Isis_Hist_free (&g);
        return;
     }

   if (-1 == (*ff->bin_eval_method)(ff, &g, (double *)sl_par->data))
     isis_throw_exception (Isis_Error);

   SLang_free_array (sl_par);
   Isis_Hist_free (&g);
}

/*}}}*/

static void eval_diff_fitfun_using_handle_intrin (Fit_Fun_MMT_Type *mmt) /*{{{*/
{
   SLang_Array_Type *sl_par=NULL, *sl_x=NULL;
   Isis_User_Grid_t g;
   Fit_Fun_t *ff;

   ff = mmt->ff;

   if ((-1 == SLang_pop_array_of_type (&sl_x, SLANG_DOUBLE_TYPE)
       || (sl_x == NULL)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }
   g.x = (double *)sl_x->data;
   g.npts = (int) sl_x->num_elements;

   if ((-1 == SLang_pop_array_of_type (&sl_par, SLANG_DOUBLE_TYPE)
       || (sl_par == NULL)))
     {
        SLang_free_array (sl_x);
        isis_throw_exception (Isis_Error);
        return;
     }

   if (sl_par->num_elements != ff->nparams)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__,
                    "%s requires %d parameters but %d were provided",
                    ff->name[0], ff->nparams, sl_par->num_elements);
        SLang_free_array (sl_par);
        SLang_free_array (sl_x);
        return;
     }

   if (-1 == (*ff->diff_eval_method)(ff, &g, (double *)sl_par->data))
     isis_throw_exception (Isis_Error);

   SLang_free_array (sl_par);
   SLang_free_array (sl_x);
}

/*}}}*/

/* This is intended for use by kernels needing to compute the
 * spectral model on a specific grid, different from the default
 * ARF grid */
int isis_eval_model_on_alt_grid (Isis_Hist_t *x) /*{{{*/
{
   User_Function_Type *f;
   Isis_Hist_t save_g;
   Isis_Hist_t *g;
   int ret = -1;

   if (SLang_get_error ())
     return -1;

   /* Validate x only to check for NULL pointers */
   if (!Isis_Hist_has_grid(x))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid grid passed to isis_eval_model_on_alt_grid");
        return -1;
     }

   f = get_user_function ();
   if ((NULL == f) || (NULL == f->fun_ptr))
     return -1;

   /* Save the current global evaluation grid,
    * then, replace the global evaluation grid with 'x'
    */

   if (NULL == (g = get_evaluation_grid ()))
     return -1;
   /* copy structs */
   save_g = *g;

   /* copy structs */
   *g = *x;

   Mode = BIN_EVAL_MODE;
   if (-1 == SLexecute_function (f->fun_ptr))
     isis_throw_exception (Isis_Error);
   /* BIN_EVAL_MODE is where we want to stay */

   ret = pop_model_result (x);

   /* Restore the global eval grid to its initial state */
   *g = save_g;

   if (Fit_Store_Model)
     {
        double *wrk, *x_wrk, *g_wrk;
        int wrk_size = g->nbins + x->nbins;
        if (NULL == (wrk = (double *) ISIS_MALLOC (wrk_size * sizeof(double))))
          return -1;
        x_wrk = wrk;
        g_wrk = x_wrk + x->nbins;
        if ((-1 == Isis_Hist_rebin_noticed (x, x_wrk, g, g_wrk))
            || (-1 == Hist_set_model (Hist_Current, H_FLUX, g->val)))
          {
             ISIS_FREE(wrk);
             return -1;
          }
        ISIS_FREE(wrk);
     }

   return ret;
}

/*}}}*/

static void get_model_on_user_grid (void) /*{{{*/
{
   Isis_Hist_t *g;
   User_Function_Type *f;

   if (NULL == (g = get_evaluation_grid ()))
     return;

   f = get_user_function ();
   if ((NULL == f) || (NULL == f->fun_ptr))
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "fit function not defined");
        return;
     }

   /* Any dataset-dependence will be resolved using
    * the user-supplied value of Isis_Active_Dataset.
    * If the user forgot to set it, the current value
    * will be used, for better or worse...
    */
   if ((-1 == Fit_sync_derived_params (Param))
       || (-1 == Fit_sync_tied_params (Param)))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "evaluating derived parameters");
        goto finish;
     }

   if (-1 == Isis_Hist_pop_valid_grid (g))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "failed initializing grid");
        goto finish;
     }

   /* Evaluate the model and leave the result on the stack */

   Mode = BIN_EVAL_MODE;
   if (-1 == SLexecute_function (f->fun_ptr))
     isis_throw_exception (Isis_Error);
   /* BIN_EVAL_MODE is the default, so dont change it back */

   finish:

   Isis_Hist_free (g);
}

/*}}}*/

static void get_differential_model (void) /*{{{*/
{
   SLang_Array_Type *slx = NULL;
   User_Function_Type *f;

   f = get_user_function ();
   if ((NULL == f) || (NULL == f->fun_ptr))
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "fit function not defined");
        return;
     }

   if ((-1 == Fit_sync_derived_params (Param))
       || (-1 == Fit_sync_tied_params (Param)))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "evaluating derived parameters");
        goto finish;
     }

   if (-1 == SLang_pop_array_of_type (&slx, SLANG_DOUBLE_TYPE)
       || slx == NULL)
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   Differential_Grid.x = (double *) slx->data;
   Differential_Grid.npts = (int) slx->num_elements;

   do_mode_eval (f->fun_ptr, DIFF_EVAL_MODE);

   Differential_Grid.x = NULL;
   Differential_Grid.npts = 0;

   finish:
   SLang_free_array (slx);
}

/*}}}*/

static int name_eval_mode_dataset_hook (Hist_t *h, void *cl) /*{{{*/
{
   SLang_Name_Type *fun_ptr;
   double d;

   (void) cl;

   if (NULL == (fun_ptr = Hist_get_instrumental_background_hook (h)))
     return 0;

   do_mode_eval (fun_ptr, NAME_EVAL_MODE);
   (void) SLang_pop_double (&d);
   return 0;
}

/*}}}*/

static void _do_name_eval_mode (void) /*{{{*/
{
   User_Function_Type *f;
   double d;

   f = get_user_function ();
   if ((NULL == f) || (NULL == f->fun_ptr))
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "fit function not defined");
        return;
     }

   do_mode_eval (f->fun_ptr, NAME_EVAL_MODE);
   (void) SLang_pop_double (&d);

   (void) map_datasets (&name_eval_mode_dataset_hook, NULL);
}

/*}}}*/

/* combining datasets */

static int init_eval_grid_method (Hist_Eval_Grid_Method_Type *m, int type, /*{{{*/
                                  void *options, void (*destroy_options)(void *))
{
   if (m == NULL)
     return -1;

   switch (type)
     {
      default:
        /* drop */
      case ISIS_EVAL_GRID_SEPARATE:
        m->type = ISIS_EVAL_GRID_SEPARATE;
        m->cache_model_values = 0;
        m->make_grid = NULL;
        m->eval_model = &eval_model_using_global_grid;
        break;

      case ISIS_EVAL_GRID_MERGED:
        m->type = ISIS_EVAL_GRID_MERGED;
        m->cache_model_values = 1;
        m->make_grid = &make_merged_eval_grid;
        m->eval_model = &eval_model_using_cached_grid;
        break;

      case ISIS_EVAL_GRID_USER:
        m->type = ISIS_EVAL_GRID_USER;
        /* m->cache_model_values is defined in the calling routine */
        m->make_grid = &make_user_eval_grid;
        m->eval_model = &eval_model_using_cached_grid;
        if (options == NULL)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "user method requires a function reference");
             return -1;
          }

        break;
     }

   m->options = options;
   m->destroy_options = destroy_options;

   return 0;
}

/*}}}*/

static int mark_dataset_combination (void) /*{{{*/
{
   Hist_Eval_Grid_Method_Type egm;
   SLang_Array_Type *sl_members = NULL;
   SLang_Array_Type *sl_weights = NULL;
   Hist_t *head;
   double *weights;
   unsigned int *members;
   unsigned int num_members;
   int gid = -1;

   if (-1 == SLang_pop_array_of_type (&sl_weights, SLANG_DOUBLE_TYPE)
       || sl_weights == NULL
       || -1 == SLang_pop_array_of_type (&sl_members, SLANG_UINT_TYPE)
       || sl_members == NULL)
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (NULL == (head = get_histogram_list_head ()))
     goto return_error;

   weights = (double *)sl_weights->data;
   members = (unsigned int *)sl_members->data;
   num_members = sl_members->num_elements;

   gid = Hist_set_combination (head, members, weights, num_members);

   (void) init_eval_grid_method (&egm, ISIS_EVAL_GRID_MERGED, NULL, NULL);
   (void) Hist_set_eval_grid_info (head, members, num_members, &egm);

   return_error:

   SLang_free_array (sl_members);
   SLang_free_array (sl_weights);

   return (gid ? gid : -1);
}

/*}}}*/

static void break_dataset_combination (unsigned int *gid) /*{{{*/
{
   Hist_t *head = get_histogram_list_head();
   if (-1 == Hist_break_combination (head, *gid))
     isis_throw_exception (Isis_Error);
}

/*}}}*/

static int hist_is_combined (int *hist_index, unsigned int *gid) /*{{{*/
{
   Hist_t *h;

   if (NULL == (h = find_hist (*hist_index)))
     return 0;

   return Hist_dataset_is_combined (h, *gid);
}

/*}}}*/

/* merging eval grids */

static void set_eval_grid_method (int *method, int *cache_model_values) /*{{{*/
{
   Hist_Eval_Grid_Method_Type egm;
   Hist_t *head = get_histogram_list_head();
   SLang_Array_Type *sl_members = NULL;
   unsigned int *members;
   unsigned int num_members;
   SLang_Name_Type *hook = NULL;

   if ((-1 == SLang_pop_array_of_type (&sl_members, SLANG_UINT_TYPE))
       || (sl_members == NULL))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   set_hook_from_stack (&hook);

   members = (unsigned int *)sl_members->data;
   num_members = sl_members->num_elements;

   /* This value is used only for user-defined methods.
    * For internal methods, it gets over-written.
    */
   egm.cache_model_values = *cache_model_values;

   if (-1 == init_eval_grid_method (&egm, *method, (void *)hook,
                                    (void (*)(void *))&SLang_free_function))
     {
        isis_throw_exception (Isis_Error);
        SLang_free_array (sl_members);
        SLang_free_function (hook);
        return;
     }

   (void) Hist_set_eval_grid_info (head, members, num_members, &egm);

   SLang_free_array (sl_members);
}

/*}}}*/

/*{{{ SLang Intrinsics */

/* DUMMY_FITFUN_MMT_TYPE is a temporary hack that will be modified to the true
  * id once the interpreter provides it when the class is registered.  See below
  * for details.  The reason for this is simple: for a module, the type-id
  * must be assigned dynamically.
  */
#define DUMMY_FITFUN_MMT_TYPE   255
#define MT DUMMY_FITFUN_MMT_TYPE
#define DUMMY_FITOBJ_MMT_TYPE   254
#define MTO DUMMY_FITOBJ_MMT_TYPE
#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define UI SLANG_UINT_TYPE
#define D SLANG_DOUBLE_TYPE
#define S SLANG_STRING_TYPE
#define R SLANG_REF_TYPE

static SLang_Intrin_Var_Type Fit_Intrin_Vars [] =
{
   MAKE_VARIABLE("Fit_Verbose", &Fit_Verbose, I, 0),
   MAKE_VARIABLE("Fit_Statistic", &Fit_Statistic, S, 0),
   MAKE_VARIABLE("Fit_Method", &Fit_Method, S, 0),
   MAKE_VARIABLE("Isis_Fit_In_Progress", &Isis_Fit_In_Progress, I, 1),
   MAKE_VARIABLE("Isis_Active_Dataset", &Isis_Active_Dataset, I, 0),
   /* Note that Isis_Active_Dataset cannot be read-only,
    * because, in some cases, it must be set by the user before
    * calling get_model_on_user_grid()=eval_fun() or
    * get_differential_model()=get_cfun().
    * This is necessary only when more than one dataset is loaded
    * and individual datasets are assigned different functions.
    */
   MAKE_VARIABLE("Isis_Active_Function_Id", &Isis_Active_Function_Id, I, 1),
   MAKE_VARIABLE("Isis_Voigt_Is_Normalized", &Isis_Voigt_Is_Normalized, I, 0),
   /* FIXME:  In isis-2 the Voigt profile will be unit-normalized.
    *         This is also the new default.  This global switch is
    *         a temporary hack to provide back-compatibility in isis-1
    */
   MAKE_VARIABLE("_num_statistic_evaluations", &Num_Statistic_Evaluations, I, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_DConstant_Type Fit_Intrin_Dconst [] =
{
   MAKE_DCONSTANT("DBL_MAX", DBL_MAX),
   MAKE_DCONSTANT("DBL_MIN", DBL_MIN),
   SLANG_END_DCONST_TABLE
};

static SLang_IConstant_Type Fit_Public_Intrin_Iconst [] =
{
   MAKE_ICONSTANT("Assigned_ARFRMF", USE_ASSIGNED_ARF_RMF),
   MAKE_ICONSTANT("Assigned_RMFARF", USE_ASSIGNED_ARF_RMF),
   MAKE_ICONSTANT("Assigned_RMF", USE_ASSIGNED_ARF_RMF),
   MAKE_ICONSTANT("Ideal_ARFRMF", USE_IDEAL_ARF_RMF),
   MAKE_ICONSTANT("Ideal_RMFARF", USE_IDEAL_ARF_RMF),
   MAKE_ICONSTANT("Ideal_ARF", USE_IDEAL_ARF),
   MAKE_ICONSTANT("Ideal_RMF", USE_IDEAL_RMF),
   MAKE_ICONSTANT("ISIS_FUN_ADDMUL", ISIS_FUN_ADDMUL),
   MAKE_ICONSTANT("ISIS_FUN_OPERATOR", ISIS_FUN_OPERATOR),
   SLANG_END_ICONST_TABLE
};

static SLang_Intrin_Fun_Type Fit_Intrinsics [] =
{
   MAKE_INTRINSIC_3("_mode_switch", mode_switch, V, UI, UI, UI),
   MAKE_INTRINSIC("_sync_model_with_data", sync_model_with_data, I, 0),
   MAKE_INTRINSIC_2("_set_fit_fun", define_user_model, V, I, S),
   MAKE_INTRINSIC_1("_load_par", load_params, I, S),
   MAKE_INTRINSIC_1("_edit_par", edit_params, V, S),
   MAKE_INTRINSIC("_set_par", set_params, V, 0),
   MAKE_INTRINSIC_I("_get_par", _get_par, V),
   MAKE_INTRINSIC_2("_set_par_fun", _set_par_fun, V, I, S),
   MAKE_INTRINSIC_I("_get_par_fun", _get_par_fun, V),
   MAKE_INTRINSIC_I("_get_param_info", _get_param_info, V),
   MAKE_INTRINSIC_1("_get_index_for_param_name", get_index_for_param_name, I, S),
   MAKE_INTRINSIC("_get_num_params", push_num_params, V, 0),
   MAKE_INTRINSIC("_do_name_eval_mode", _do_name_eval_mode, V, 0),
   MAKE_INTRINSIC_I("_freeze", _freeze, V),
   MAKE_INTRINSIC_I("_thaw", _thaw, V),
   MAKE_INTRINSIC_II("_tie", _tie, V),
   MAKE_INTRINSIC_I("_untie", _untie, V),
   MAKE_INTRINSIC_1("_fit", find_best_fit, I, R),
   MAKE_INTRINSIC_1("_eval_model", eval_model, I, R),
   MAKE_INTRINSIC_1("_eval_statistic_only", eval_statistic_only, I, R),
   MAKE_INTRINSIC_4("_confidlev", confidence_limits, V, I, D, I, D),
   MAKE_INTRINSIC_I("_set_conf_limit_search", _set_conf_limit_search, V),
   MAKE_INTRINSIC("_array_fit", array_fit, V, 0),
   MAKE_INTRINSIC("_get_differential_model", get_differential_model, V, 0),
   MAKE_INTRINSIC("_get_model_on_user_grid", get_model_on_user_grid, V, 0),
   MAKE_INTRINSIC_I("_get_instrumental_background", get_instrumental_background, V),
   MAKE_INTRINSIC_2("_set_fit_type", set_fit_type, V, I, I),
   MAKE_INTRINSIC_2("_load_fit_method", load_fit_method, I, S, S),
   MAKE_INTRINSIC_2("_load_fit_statistic", load_fit_statistic, I, S, S),
   MAKE_INTRINSIC("_add_slangfun_intrin", add_slangfun_intrin, V, 0),
   MAKE_INTRINSIC("_add_cfun_intrin", add_cfun_intrin, V, 0),
   MAKE_INTRINSIC_1("_del_function", del_function, V, S),
   MAKE_INTRINSIC_2("_set_function_category", set_function_category, V, S, UI),
   MAKE_INTRINSIC("_function_list", function_list, V, 0),
   MAKE_INTRINSIC("_get_fit_fun", _get_fit_fun, V, 0),
   MAKE_INTRINSIC("_get_fit_method_name", get_fit_method_name, V, 0),
   MAKE_INTRINSIC("_get_fit_statistic_name", get_fit_statistic_name, V, 0),
   MAKE_INTRINSIC_1("_set_fit_method_name", set_fit_method_name, V, S),
   MAKE_INTRINSIC_1("_set_fit_statistic_name", set_fit_statistic_name, V, S),
   MAKE_INTRINSIC("_set_fit_range_hook", set_fit_range_hook, V, 0),
   MAKE_INTRINSIC("_set_fit_constraint_fun", set_fit_constraint_fun, V, 0),
   MAKE_INTRINSIC("set_define_model_hook_intrin", set_define_model_hook_intrin, V, 0),
   MAKE_INTRINSIC_1("_set_slangfun_param_default_hook", set_slangfun_param_default_hook, V, S),
   MAKE_INTRINSIC_1("_set_kernel", _set_kernel, V, UI),
   MAKE_INTRINSIC_2("_load_kernel", _load_kernel, I, S, S),
   MAKE_INTRINSIC_I("_print_kernel", _print_kernel, V),
   MAKE_INTRINSIC("_list_kernels", _list_kernels, V, 0),
   MAKE_INTRINSIC_1("_add_slang_statistic", _add_slang_statistic, V, S),
   MAKE_INTRINSIC_2("_add_slang_optimizer", add_slang_fit_engine_intrin, V, S, S),
   MAKE_INTRINSIC("list_statistics_and_engines", list_statistics_and_engines, V, 0),
   MAKE_INTRINSIC("_define_hist_combination", mark_dataset_combination, I, 0),
   MAKE_INTRINSIC_1("_undefine_hist_combination", break_dataset_combination, V, UI),
   MAKE_INTRINSIC_2("_hist_is_combined", hist_is_combined, I, I, UI),
   MAKE_INTRINSIC_2("_set_eval_grid_method", set_eval_grid_method, V, I, I),
   MAKE_INTRINSIC_1("eval_fitfun_using_handle_intrin", eval_fitfun_using_handle_intrin, V, MT),
   MAKE_INTRINSIC_1("eval_diff_fitfun_using_handle_intrin", eval_diff_fitfun_using_handle_intrin, V, MT),
   MAKE_INTRINSIC_1("get_fitfun_handle_intrin", push_mmt_fitfun_type_intrin, V, S),
   MAKE_INTRINSIC_1("get_fitfun_info", Fit_get_fun_info, V, S),
   MAKE_INTRINSIC_4("set_hard_limits", Fit_set_hard_limits, I, S, S, D, D),
   MAKE_INTRINSIC("open_fit_object_mmt_intrin", open_fit_object_mmt_intrin, V, 0),
   MAKE_INTRINSIC_2("fobj_eval_statistic", fobj_eval_statistic, V, MTO, I),
   MAKE_INTRINSIC_2("fobj_eval_residuals", fobj_eval_residuals, V, MTO, I),
   MAKE_INTRINSIC_1("fobj_get_data_weights", fobj_get_data_weights, V, MTO),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef V
#undef I
#undef UI
#undef D
#undef S
#undef R
#undef MT
#undef MTO

/*}}}*/

#ifdef __cplusplus
extern "C"
#endif
void deinit_fit_module (void);

void deinit_fit_module (void) /*{{{*/
{
   free_user_function ();
   free_kernel_table (Kernel_Table);
   Fit_free_param_table (Param);  Param = NULL;

   ISIS_FREE (Fit_Method);  Fit_Method = NULL;
   ISIS_FREE (Fit_Statistic); Fit_Statistic = NULL;
   SLang_free_function (Fit_Range_Hook); Fit_Range_Hook = NULL;
   SLang_free_function (Define_Model_Hook); Define_Model_Hook = NULL;
   SLang_free_function (Fit_Constraint);  Fit_Constraint = NULL;

   deinit_fit_functions ();
   deinit_fit_engine ();

   bin_eval_mode ();
   Fit_Verbose = 0;
   Fit_Data_Type = 0;
   Use_Interactive_Param_Init = 0;
   Response_Type = USE_ASSIGNED_ARF_RMF;
}

/*}}}*/

int init_fit_module_internals (void) /*{{{*/
{
   char *minimizer = Fit_Default_Fit_Method;
   char *statistic = Fit_Default_Fit_Statistic;

   ISIS_FREE (Fit_Method);
   if (NULL == (Fit_Method = isis_make_string (minimizer)))
     return -1;

   ISIS_FREE (Fit_Statistic);
   if (NULL == (Fit_Statistic = isis_make_string (statistic)))
     {
        ISIS_FREE (Fit_Method);
        return -1;
     }

   if ((-1 == init_fit_engine ())
       || (-1 == init_fit_functions ())
       || (-1 == init_kernel_table (Kernel_Table)))
     return isis_trace_return(-1);

   return 0;
}

/*}}}*/

static void patchup_intrinsic_table (unsigned int dummy_id, unsigned int assigned_id) /*{{{*/
{
   SLang_Intrin_Fun_Type *f;

   f = Fit_Intrinsics;
   while (f->name != NULL)
     {
        unsigned int i, nargs;
        SLtype *args;

        nargs = f->num_args;
        args = f->arg_types;
        for (i = 0; i < nargs; i++)
          {
             if (args[i] == dummy_id)
               args[i] = assigned_id;
          }

        /* For completeness */
        if (f->return_type == dummy_id)
          f->return_type = assigned_id;

        f++;
     }
}

/*}}}*/

static void destroy_fitfun_mmt_type (SLtype type, VOID_STAR f) /*{{{*/
{
   Fit_Fun_MMT_Type *mt = (Fit_Fun_MMT_Type *) f;
   (void) type;
   SLfree ((char *) mt);
}

/*}}}*/

SLANG_MODULE(fit);
int init_fit_module_ns (char *ns_name) /*{{{*/
{
   SLang_Class_Type *cl;
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if (Fit_Fun_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("Fitfun_Type")))
          return isis_trace_return(-1);
        (void) SLclass_set_destroy_function (cl, destroy_fitfun_mmt_type);

        /* By registering as SLANG_VOID_TYPE,
         * slang will dynamically allocate a type
         */
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (Fit_Fun_MMT_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return isis_trace_return(-1);

        Fit_Fun_Type_Id = SLclass_get_class_id (cl);
        patchup_intrinsic_table (DUMMY_FITFUN_MMT_TYPE, Fit_Fun_Type_Id);
     }

   if (Fit_Object_MMT_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("Fit_Object_Type")))
          return isis_trace_return(-1);
        (void) SLclass_set_destroy_function (cl, destroy_fit_object_mmt_type);

        /* By registering as SLANG_VOID_TYPE,
         * slang will dynamically allocate a type
         */
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (Fit_Object_MMT_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return isis_trace_return(-1);

        Fit_Object_MMT_Type_Id = SLclass_get_class_id (cl);
        patchup_intrinsic_table (DUMMY_FITOBJ_MMT_TYPE, Fit_Object_MMT_Type_Id);
     }

   if ((-1 == SLns_add_intrin_fun_table (ns, Fit_Intrinsics, NULL))
       || (-1 == SLns_add_dconstant_table (ns, Fit_Intrin_Dconst, NULL))
       || (-1 == SLns_add_intrin_var_table (pub_ns, Fit_Intrin_Vars, NULL))
       || (-1 == SLns_add_iconstant_table (pub_ns, Fit_Public_Intrin_Iconst, NULL)))
     return isis_trace_return(-1);

   return init_fit_module_internals ();
}

/*}}}*/

