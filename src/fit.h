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

/* $Id: fit.h,v 1.43 2004/03/01 13:18:12 houck Exp $ */

#ifndef ISIS_FIT_H
#define ISIS_FIT_H

#define MAX_NAME_SIZE  32

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

enum
{
   OUT_PARAM_NAME_SIZE = (2*MAX_NAME_SIZE + 5)
};

extern char *Fit_Default_Fit_Method;
extern char *Fit_Default_Fit_Statistic;

extern int Isis_Evaluating_Derived_Param;
extern int Isis_Fit_In_Progress;
extern int Fit_Loading_Parameters_From_File;

extern double Isis_Inf;
extern double Isis_Nan;

typedef struct _Fit_Param_t Fit_Param_t;
typedef struct _Param_Info_t Param_Info_t;
typedef struct Param_t Param_t;
typedef struct _Fit_Fun_t Fit_Fun_t;
typedef char Fit_Fun_Name_t[MAX_NAME_SIZE];

extern Fit_Fun_t *Fit_get_fit_fun (int fun_type);
extern int Fit_get_fun_type (char *name);
extern int Fit_get_fun_par (Fit_Fun_t *ff, char *par_name);
extern void Fit_get_fun_info (char *name);
extern int Fit_is_valid_fit_fun (Fit_Fun_t *ff_test);
extern int Fit_set_fun_post_hook (char *fun_name);
extern int Fit_set_fun_trace_hook (char *fun_name);

struct _Fit_Param_t
{
   double *par;
   double *par_min;
   double *par_max;
   double *step;
   double *relstep;
   int *idx;
   int npars;
};

struct _Fit_Fun_t
{
   Fit_Fun_t *next;
   unsigned int fun_type;      /* e.g. Gaussian, Lorentzian, mekal, ... */
   unsigned int fun_version;
   unsigned int nparams;

   union
     {
        Isis_Binned_Function_t *c;  /* ptr to C function */
        SLang_Name_Type *sl;        /* ptr to S-Lang function */
     }
   fun;
   int (*bin_eval_method)(Fit_Fun_t *, Isis_Hist_t *, double *, SLang_Struct_Type *);

   int (*diff_eval_method)(Fit_Fun_t *, Isis_User_Grid_t *, double *, SLang_Struct_Type *);

   Isis_User_Source_t s;       /* optional function ptrs    */
   Fit_Fun_Name_t *name;       /* = { fcn_name, par_name0, par_name1, ... } */
   Fit_Fun_Name_t *unit;       /* = {           par_unit0, par_unit1, ... } */

   int (*set_param_default)(Fit_Fun_t *, Param_Info_t *);
   SLang_Name_Type *slangfun_param_default;
   Isis_Arg_Type *slangfun_param_default_args;
   SLang_Name_Type *slangfun_diff_eval;

   SLang_Name_Type *trace_hook;
   SLang_Name_Type *post_hook;

   int (*set_param_hard_limits)(Fit_Fun_t *, unsigned int, int, double, double, double, double, double);

   void (*destroy_fun)(Fit_Fun_t *);
   void *client_data;
};

struct _Param_Info_t
{
   double value;
   double min, max;             /* max allowed range (user-defined) */
   double hard_min, hard_max;   /* max allowed range (model validity) */
   double step;                 /* suggested absolute parameter step */
   double relstep;              /* suggested relative parameter step */

   char *fun_str;
   SLang_Name_Type *fun_ptr;     /* param may be a function of other params */
   char *tie_param_name;         /* fname(i).parname */

   char *param_name;             /* fname(i).parname */
   unsigned int fun_type;
   unsigned int fun_version;
   unsigned int fun_id;
   unsigned int fun_par;
   unsigned int idx;
   unsigned int vary_idx;
   unsigned int freeze;         /* boolean */
   unsigned int default_applied;
   unsigned int set_minmax;
   unsigned int is_a_norm;      /* boolean */
   unsigned int in_use;         /* boolean */
};

extern Param_t *Fit_new_param_table (unsigned int num_params,
                                     unsigned int fun_type, unsigned int fun_id, unsigned int fun_version);
extern void Fit_free_param_table (Param_t *pt);
extern void Fit_reset_param_lookup_table (Param_t *pt);
extern int Fit_mark_params_unused (Param_t *pt);
extern Param_Info_t *Fit_param_info (Param_t *pt, unsigned int idx);
extern Param_Info_t *Fit_param_info2 (Param_t *pt, unsigned int fun_type, unsigned int fun_id, unsigned int fun_par);
extern Param_Info_t *Fit_variable_param_info (Param_t *pt, unsigned int vary_idx);
extern Param_Info_t *Fit_find_param_info_by_full_name (Param_t *pt, char *name);
extern int Fit_pack_all_params (Param_t *pt, Fit_Param_t *p);
extern int Fit_unpack_all_params (Param_t *pt, double *par);
extern int Fit_sync_tied_params (Param_t *pt);
extern int Fit_sync_derived_params (Param_t *pt);
extern int Fit_count_params (Param_t *pt, int *num_all, int *num_vary);
extern int Fit_pack_variable_params (Param_t *pt, Fit_Param_t *p);
extern int Fit_unpack_variable_params (Param_t *pt, double * par);

extern int Fit_register_fun (Param_t *pt, Fit_Fun_t *ff, unsigned int fun_id,
                             unsigned int *addr);
extern int Fit_set_fun_params (Param_t *pt, unsigned int fun_type, unsigned int fun_id,
                               double *par, double *par_min, double *par_max);
extern int Fit_get_fun_params (Param_t *pt, unsigned int fun_type, unsigned int fun_id,
                               double *par);
extern int Fit_copy_fun_params (char *fun_name, unsigned int fun_id, double **par, unsigned int *num_pars);

extern int Fit_get_param_value (Param_t *pt, unsigned int idx, double *value);
extern int Fit_set_param_value (Param_t *pt, unsigned int idx, double value);
extern int Fit_set_param_function (Param_t *pt, unsigned int idx, char *fun_str);
extern int Fit_tie (Param_t *pt, unsigned int idx_a, int unsigned idx_b);
extern int Fit_untie (Param_t *pt, unsigned int idx);
extern int Fit_set_freeze (Param_t *pt, unsigned int idx, unsigned int freeze_value);
extern int Fit_set_param_control (Param_t *pt, unsigned int idx,
                                  int update_minmax, double min, double max,
                                  int freeze, char *tie, double step, double relstep);

extern int Fit_set_hard_limits (char *fun_name, char *par_name, int *idx, int *hard_limits_only);
extern int Fit_set_param_hard_limits (int idx, int hard_limits_only, Param_Info_t *p);
extern int Fit_set_param_hard_limits1 (Param_t *pt, int idx, int hard_limits_only, 
                                       Param_Info_t *pnew);

extern int Fit_Append_Builtin_Functions (void);

/* kernel */

extern Isis_Kernel_Def_t * Fit_new_kernel (void);
extern int Fit_add_kernel_function (Isis_Kernel_Def_t *def);
extern void Fit_free_kernel (Isis_Kernel_Def_t *def);
extern void Fit_free_aux_kernels (Isis_Kernel_Def_t *t);
extern void Fit_push_kernel_names (Isis_Kernel_Def_t *t);
extern int Fit_append_kernel (Isis_Kernel_Def_t *head, Isis_Kernel_Def_t *def);
extern Isis_Kernel_Def_t * Fit_find_kernel (Isis_Kernel_Def_t *t, unsigned int kernel_id);
extern Isis_Kernel_Def_t * Fit_find_kernel_by_name (Isis_Kernel_Def_t *t, char *kernel_name);
extern double *Fit_get_kernel_params (Param_t *pt, int hist_index, Isis_Kernel_Def_t *def);
extern int Fit_set_kernel_param_default (Isis_Kernel_Def_t *def,
                                         int fun_par, Param_Info_t *p);

typedef struct Kernel_Info_t Kernel_Info_t;
typedef struct Kernel_Table_t Kernel_Table_t;

struct Kernel_Info_t
{
   Kernel_Info_t *next;
   Kernel_Info_t *saved;
   char *kernel_params;
   unsigned int kernel_id;
   int hist_index;
};

struct Kernel_Table_t
{
   Isis_Kernel_Def_t *kernel_defs;
   Kernel_Info_t kernel_info;
};

extern int init_kernel_table (Kernel_Table_t *t);
extern void free_kernel_table (Kernel_Table_t *t);
extern Kernel_Table_t *get_kernel_table (void);
extern Kernel_Info_t *find_kernel_info (Kernel_Table_t *t, int hist_index);

extern int init_fit_engine (void);
extern void deinit_fit_engine (void);

/* conf */

typedef struct Isis_Fit_CLC_Type Isis_Fit_CLC_Type;
typedef struct Fit_Data_t Fit_Data_t;
typedef struct Search_Info_Type Search_Info_Type;
typedef struct Sample_Type Sample_Type;
typedef struct Fit_Info_Type Fit_Info_Type;

struct Fit_Info_Type
{
   Fit_Param_t *par;
   int num_pars;
   int num_vary;
   int verbose;
   void *cl_verbose;
};

/* CLC=confidence limit control type. */
struct Isis_Fit_CLC_Type
{
   double tol;
   double delta_stat;
   int verbose;
};

typedef struct
{
   double *data;
   double *weight;
   int num;
   int malloced;
}
Fit_Object_Data_Type;

typedef struct
{
   Isis_Fit_Type *ft;
   /* Isis_Fit_Type:  optimizer + statistic + objective function */

   Fit_Info_Type *info;
   /* Fit_Info_Type:  parameter info */

   int (*unpack)(Param_t *, double *);
   /* param unpack method */

   Fit_Data_t *d;
   /* Fit_Data_t: raw data being fit */

   Fit_Object_Data_Type *dt;
   /* Fit_Object_Data_Type:  pointers to support combining datasets
    * If no datasets are being combined, Fit_Object_Data_Type contains
    * copies of the Fit_Data_t pointers.
    * If datasets are being combined,  Fit_Object_Data_Type contains
    * pointers to malloced workspace.
    */
}
Fit_Object_Type;

extern Fit_Param_t *new_fit_param_type (unsigned int num);
extern void free_fit_param_type (Fit_Param_t *p);
extern void free_fit_data (Fit_Data_t *d);

extern int fit_object_config (Fit_Object_Type *fo, Param_t *pt, int unpack_variable);
extern int fit_statistic (Fit_Object_Type *fo, int optimize, double *statistic, int *num_data_bins);

extern void verbose_warn_hook (void *clientdata, const char * fmt, ...);
extern void init_verbose_hook (void);
extern void deinit_verbose_hook (void);
extern void verbose_info_hook (void *clientdata, double chisqr,
                               double *par, unsigned int npars_vary);

extern int get_confidence_limits (Fit_Object_Type *fo, Param_t *pt, Isis_Fit_CLC_Type *ctrl,
                                  int idx, double *pconf_min, double *pconf_max);

/* engine */

extern Isis_Fit_Engine_Type *isis_find_fit_engine (char *name);
extern Isis_Fit_Statistic_Type *isis_find_fit_statistic (char *name);
extern void isis_fit_close_fit (Isis_Fit_Type *);
extern Isis_Fit_Type *isis_fit_open_fit
  (char *method_name, char *statistic_name, Isis_Fit_Fun_Type *f,
   SLang_Name_Type *constraint_fun);
extern int isis_fit_perform_fit
  (Isis_Fit_Type *f, void *clientdata,
   double *x, double *y, double *weights, unsigned int npts,
   double *pars, unsigned int npars, double *statistic);
extern int isis_fit_report_statistic
  (Isis_Fit_Type *f, FILE *fp, double statistic, unsigned int npts, unsigned int nvpars);
extern int isis_invalid_params (Isis_Fit_Engine_Type *e,
                                double *pars, unsigned int npars);

extern void set_slopt_fit_object (SLang_MMT_Type *mmt);
extern int is_slang_optimizer (Isis_Fit_Type *f);

extern int isis_fit_set_ranges (Isis_Fit_Type *f, double *par_min, double *par_max);
extern int isis_fit_set_param_step (Isis_Fit_Type *f, double *step, double *relstep);
extern int isis_fit_set_warn_hook (Isis_Fit_Type *f, Isis_Fit_Warning_Hook_Type *v);
extern int isis_fit_set_verbose_hook (Isis_Fit_Type *f, Isis_Fit_Verbose_Hook_Type *v);
extern int isis_fit_set_range_hook (Isis_Fit_Type *f, Isis_Fit_Range_Hook_Type *r);
extern int isis_fit_set_verbose_level (Isis_Fit_Type *f, int verbose);
extern int isis_delta_stat_is_chisqr_distrib (Isis_Fit_Type *f);

extern int isis_fit_load_statistic (char *file, char *name);
extern int isis_fit_load_fit_routine (char *file, char *name, char *sname);

extern void add_slang_fit_engine_intrin (char *eng_name, char *stat_name);
extern int isis_fit_add_engine (char *name, char *sname, Isis_Fit_Engine_Init_Type *init);
extern int isis_fit_add_statistic (char *name, Isis_Fit_Statistic_Init_Type *x);

/* Specific fitting methods */
extern Isis_Fit_Engine_Type *Isis_lmdif_feng (char *name, char *sname);
extern Isis_Fit_Engine_Type *Isis_marquardt_feng (char *name, char *sname);
extern Isis_Fit_Engine_Type *Isis_mpfit_feng (char *name, char *sname);
extern Isis_Fit_Engine_Type *Isis_subplex_feng (char *name, char *sname);
extern Isis_Fit_Engine_Type *Isis_simann_feng (char *name, char *sname);
extern Isis_Fit_Statistic_Type *Isis_chisqr_stat (void);
extern Isis_Fit_Statistic_Type *Isis_cash_stat (void);
extern Isis_Fit_Statistic_Type *Isis_ml_stat (void);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
