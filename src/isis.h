#ifndef ISIS_H /* -*- mode: C; mode: fold -*- */
#define ISIS_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2011  Massachusetts Institute of Technology

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

/* $Id: isis.h,v 1.30 2004/06/06 02:42:22 houck Exp $ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <slang.h>

#define ISIS_VERSION          10601
#define ISIS_VERSION_STRING  "1.6.1-46"
#define ISIS_VERSION_PREFIX   1.6.1

#define ISIS_API_VERSION 5

enum
{
   ISIS_MAX_PROTON_NUMBER = 36
};

/* globals */

extern int Isis_List_Filenames;

extern unsigned int Isis_Version;
extern int Isis_Error;
extern int Isis_Trace;

extern void (*Isis_Errno_Hook)(int);
extern int isis_user_break (void);

extern char *Isis_Srcdir;
extern char *Isis_Pager;
extern char *Isis_Public_Namespace_Name;

/*{{{ memory handling */

extern void *isis_realloc(void *ptr, size_t size);
extern void *isis_malloc(size_t size);

#define ISIS_MALLOC  isis_malloc
#define ISIS_REALLOC  isis_realloc
#define ISIS_FREE(p) do {if (p) free ((void *)p); p = NULL;} while (0)

/*}}}*/

/*{{{ dynamic linking */

#define _x_0(a,b)  a##_##b
#define _x_(a,b)   _x_0(a,b)

typedef void (*isis_fptr_type)(void);
extern isis_fptr_type isis_load_function (const char *path, const char *name, const char *type);

/*}}}*/

extern double isis_gpf (double x);

extern int isis_strcasecmp (const char *pa, const char *pb);
extern char *isis_strcpy (char *dest, const char *src, int size);
extern int isis_strcat (char *dest, int size, ...);
extern char *isis_mkstrcat (const char *arg, ...);
extern char *isis_make_string (const char *str);

/*{{{ RMFs */

typedef struct
{
   double *bin_lo;
   double *bin_hi;
   unsigned int nbins;
   int units;
}
Isis_Rmf_Grid_Type;

extern void Isis_free_rmf_grid (Isis_Rmf_Grid_Type *eb);
extern Isis_Rmf_Grid_Type *Isis_new_rmf_grid (unsigned int nbins, double *lo, double *hi);

typedef struct Isis_Rmf_t Isis_Rmf_t;

struct Isis_Rmf_t
{
   Isis_Rmf_t *next;
   int index;
   int ref_count;

   int method;
   int includes_effective_area;

   int (*set_arf_grid)(Isis_Rmf_t *, double *, double *, unsigned int);
   int (*set_data_grid)(Isis_Rmf_t *, double *, double *, unsigned int);
   /* grids are always in angstrom units, in increasing order;
    * data (int *chan) array is assumed to have offset zero. */

   int (*init) (Isis_Rmf_t *);

   int (*get_arf_grid)(Isis_Rmf_t *, double **, double **, unsigned int *);
   int (*get_data_grid)(Isis_Rmf_t *, double **, double **, unsigned int *, int *);
   int (*redistribute)(Isis_Rmf_t *, unsigned int, double, double *, unsigned int);
   int (*set_noticed_model_bins)(Isis_Rmf_t *, int, int *, int, int *);
   void (*delete_client_data)(Isis_Rmf_t *);

   int (*rebin_rmf) (Isis_Rmf_t *, double *, double *, unsigned int);
   int (*factor_rsp) (Isis_Rmf_t *, double *);

#define ISIS_RMF_BUFSIZE  72
   int order;
   char    grating[ISIS_RMF_BUFSIZE];
   char instrument[ISIS_RMF_BUFSIZE];
   char *arg_string;		       /* may be NULL */
   void *client_data;
};

typedef int Isis_Rmf_Load_Method_t (Isis_Rmf_t *, void *);

#ifdef __cplusplus
# define ISIS_RMF_METHOD(name,a,b) \
   extern "C" int Isis_##name##_api_version; \
   int Isis_##name##_api_version = ISIS_API_VERSION; \
   extern "C" int Isis_##name##_rmf(Isis_Rmf_t *, char *); \
   int Isis_##name##_rmf(Isis_Rmf_t *a, char *b)
#else
# define ISIS_RMF_METHOD(name,a,b) \
   int Isis_##name##_api_version = ISIS_API_VERSION; \
   extern int Isis_##name##_rmf(Isis_Rmf_t *, char *); \
   int Isis_##name##_rmf(Isis_Rmf_t *a, char *b)
#endif

extern int Isis_Rmf_OGIP_Compliance;

/*}}}*/

/*{{{ ARFs */

typedef struct Isis_Arf_t Isis_Arf_t;

struct Isis_Arf_t
{
   Isis_Arf_t *next;

   double exposure;              /* exposure time [sec] */

   double *bin_lo;               /* bin left edge */
   double *bin_hi;               /* bin right edge */
   double *arf;                  /* Effective area [count cm^2 / photon] */
   double *arf_err;              /* uncertainty in eff. area [count cm^2 / photon] */
   union
     {
        double s;
        double *v;
     }
   fracexpo;
   int fracexpo_is_vector;       /* boolean */

   int nbins;                    /* number of bins */

   int index;                    /* id number */
   int ref_count;                /* count how many Hist_t's point to here */
   int is_identity;

   int order;                    /* TG_M  = spectral order */
   int part;                     /* TG_PART = code indicating heg, meg, leg, etc. */
   int srcid;                    /* TG_SRCID = code indicating which source in the fov */

#define ISIS_ARF_VALUE_SIZE      72
   char object[ISIS_ARF_VALUE_SIZE];      /* target name */
   char grating[ISIS_ARF_VALUE_SIZE];     /* HETG or LETG */
   char instrument[ISIS_ARF_VALUE_SIZE];  /* ACIS-S or HRC-S */

   char *file;     /* name of input ARF file */
};

/*}}}*/

typedef struct Isis_Rsp_t Isis_Rsp_t;

struct Isis_Rsp_t
{
   Isis_Rsp_t *next;
   Isis_Rmf_t *rmf;
   Isis_Arf_t *arf;
   unsigned int is_flux_data;
};

typedef struct
{
   double *val, *val_err;
   double *bin_lo, *bin_hi;
   double *sys_err_frac;
   int nbins;
   int *notice, *notice_list;
   int n_notice;
}
Isis_Hist_t;
#define ISIS_HIST_INIT {NULL,NULL,NULL,NULL,NULL,0,NULL,NULL,0}

extern int Isis_Hist_allocate (unsigned int n, Isis_Hist_t *x);
extern void Isis_Hist_free (Isis_Hist_t *x);
extern int Isis_Hist_push_noticed_grid (Isis_Hist_t *g);
extern int Isis_Hist_rebin_noticed (Isis_Hist_t *x, double *x_tmp,
                                    Isis_Hist_t *g, double *g_tmp);
extern int Isis_Hist_pop_valid_grid (Isis_Hist_t *g);
extern int Isis_Hist_has_grid (Isis_Hist_t *x);

/*{{{ Fit-kernel */

typedef struct Isis_Kernel_t Isis_Kernel_t;
typedef struct Isis_Kernel_Def_t Isis_Kernel_Def_t;

typedef struct
{
   double exposure_time;
   double frame_time;
   double dtcor;
   Isis_Rsp_t rsp;
   int (*apply_rmf)(Isis_Rmf_t *, double *, int, double *, int *, int);
   int num_orig_data;
   int tg_part;
   int tg_m;
}
Isis_Obs_t;

struct Isis_Kernel_Def_t
{
   char *kernel_name;
   Isis_Kernel_t *(*allocate_kernel) (Isis_Obs_t *, char *);
   int (*allows_ignoring_model_intervals)(Isis_Kernel_t *);
   unsigned int num_kernel_parms;
   char **kernel_parm_names;
   char **kernel_parm_units;
   double *default_min;
   double *default_max;
   double *default_hard_min;
   double *default_hard_max;
   double *default_value;
   double *default_step;
   double *default_relstep;
   unsigned int *default_freeze;

#ifdef ISIS_SRC
   Isis_Kernel_Def_t *next;
   unsigned int fun_type;
   unsigned int kernel_id;
   unsigned int malloced_kernel;
#endif
};

#define ISIS_NULL_KERNEL_DEF {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}

struct Isis_Kernel_t
{
   int (*compute_kernel)(Isis_Kernel_t *, double *, Isis_Hist_t *, double *, unsigned int,
                         int (*)(Isis_Hist_t *));
   int (*compute_flux)(Isis_Kernel_t *, double *, unsigned int, Isis_Hist_t *,
                       double *, double *, double *, double **, char *);
   void (*delete_kernel)(Isis_Kernel_t *);
   int (*print_kernel)(Isis_Kernel_t *);
   Isis_Kernel_Def_t *kernel_def;

   double frame_time;
   double exposure_time;
   Isis_Rsp_t rsp;
   int (*apply_rmf)(Isis_Rmf_t *, double *, int, double *, int *, int);
   unsigned int num_orig_data;

   char *params;

   int tg_part;
   int tg_m;

#ifdef ISIS_KERNEL_PRIVATE_DATA
   ISIS_KERNEL_PRIVATE_DATA
#endif
};

#ifdef __cplusplus
# define ISIS_USER_KERNEL_MODULE(name,b,c) \
   extern "C" int Isis_##name##_api_version; \
   int Isis_##name##_api_version = ISIS_API_VERSION; \
   extern "C" int Isis_##name##_kernel(Isis_Kernel_Def_t *, char *); \
   int Isis_##name##_kernel(Isis_Kernel_Def_t *b, char *c)
#else
# define ISIS_USER_KERNEL_MODULE(name,b,c) \
   int Isis_##name##_api_version = ISIS_API_VERSION; \
   extern int Isis_##name##_kernel(Isis_Kernel_Def_t *, char *); \
   int Isis_##name##_kernel(Isis_Kernel_Def_t *b, char *c)
#endif

typedef int Isis_Kernel_Init_t (Isis_Kernel_Def_t *, char *);
extern Isis_Kernel_t *isis_init_kernel (Isis_Kernel_t *, unsigned int, Isis_Obs_t *);
extern double *isis_unit_source_response (Isis_Kernel_t *k);

/*}}}*/

/*{{{ option strings */

typedef struct
{
   char *subsystem;
   unsigned int num_options;
   char **option_names;
   char **option_values;
}
Isis_Option_Type;
extern Isis_Option_Type *isis_parse_option_string (char *str);
extern void isis_free_options (Isis_Option_Type *);

typedef struct
{
   char *optname;
   int (*fun)(char *, char *, char *, void *);
   unsigned int value_flags;
#define ISIS_OPT_REQUIRES_VALUE  1
#define ISIS_OPT_NO_VALUE        2
   char *default_value_string;
   char *help_string;
}
Isis_Option_Table_Type;
#define ISIS_OPTION_TABLE_TYPE_NULL  {NULL, NULL, 0, NULL, NULL}

extern int isis_process_options (Isis_Option_Type *opt,
                                 Isis_Option_Table_Type *table,
                                 void *client_data,
                                 int err_on_unsupported);
extern int isis_is_option_present (char *options, char *option_name, char **value);
/* *value is a malloced pointer */

extern char *isis_add_option (char *subsys, char *option);

extern int isis_get_option_d (Isis_Option_Type *opt, char *option, double *d);
extern int isis_get_option_i (Isis_Option_Type *opt, char *option, int *d);
extern int isis_get_option_u (Isis_Option_Type *opt, char *option, unsigned int *d);

extern int isis_update_option_string (char **optstring, char *optname, char *optvalue);
extern char *isis_make_default_option_string (const char *subsystem, Isis_Option_Table_Type *table);

/*}}}*/

/*{{{ User-defined fit-functions */

typedef struct
{
   double *x;
   int npts;
}
Isis_User_Grid_t;

typedef int Isis_Binned_Function_t (double *, Isis_Hist_t *, double *, unsigned int);
typedef int Isis_Unbinned_Function_t (double *, Isis_User_Grid_t *, double *, unsigned int);

enum
{
   ISIS_FUN_ADDMUL = 0,
   ISIS_FUN_OPERATOR = 1
};

typedef struct
{
   Isis_Binned_Function_t *binned;
   Isis_Unbinned_Function_t *unbinned;
   void (*function_exit) (void);
   int (*init_params_from_screen) (unsigned int, double *, double *, double *, unsigned int);
   unsigned int *norm_indexes;
   unsigned int num_norms;
   char **parameter_names;
   char **parameter_units;
   double *default_min;
   double *default_max;
   double *default_hard_min;
   double *default_hard_max;
   double *default_value;
   double *default_step;
   double *default_relstep;
   unsigned int *default_freeze;
   unsigned int num_parameters;
   unsigned int category;          /* addmul, operator, etc. */
}
Isis_User_Source_t;

#ifdef __cplusplus
# define ISIS_USER_SOURCE_MODULE(name,p,q) \
    extern "C" int Isis_##name##_api_version; \
    int Isis_##name##_api_version = ISIS_API_VERSION; \
    extern "C" int Isis_##name##_function(Isis_User_Source_t *, char *); \
    int Isis_##name##_function(Isis_User_Source_t *p, char *q)
#else
# define ISIS_USER_SOURCE_MODULE(name,p,q) \
    int Isis_##name##_api_version = ISIS_API_VERSION; \
    extern int Isis_##name##_function(Isis_User_Source_t *, char *); \
    int Isis_##name##_function(Isis_User_Source_t *p, char *q)
#endif

typedef int Isis_User_Source_Init_Fun_t (Isis_User_Source_t *, char *);

extern int isis_eval_model_on_alt_grid (Isis_Hist_t *x);
extern int Isis_Add_Static_Fun (Isis_User_Source_Init_Fun_t *us_init, char *us_name);

extern double Isis_Default_Relstep;

/*}}}*/

typedef struct Isis_Fit_Type Isis_Fit_Type;
typedef struct Isis_Fit_Engine_Type Isis_Fit_Engine_Type;
typedef struct Isis_Fit_Statistic_Type Isis_Fit_Statistic_Type;

/*{{{ User-defined fit-statistic */

typedef int Isis_Fit_Statistic_Fun_Type
  (Isis_Fit_Statistic_Type *st, double *y, double *fx, double *w, unsigned int npts,
   double *vec, double *val);
typedef int Isis_Fit_Statistic_Report_Type
  (Isis_Fit_Statistic_Type *st, void *fp, double statistic,  unsigned int npts, unsigned int nvpars);
typedef Isis_Fit_Statistic_Type *Isis_Fit_Statistic_Init_Type (void);

struct Isis_Fit_Statistic_Type
{
   Isis_Fit_Statistic_Fun_Type *compute_statistic;
   SLang_Name_Type *sl_fun;

   Isis_Fit_Statistic_Report_Type *report;
   SLang_Name_Type *sl_report;

   int (*set_options)(Isis_Fit_Statistic_Type *, Isis_Option_Type *);
   void (*deallocate)(Isis_Fit_Statistic_Type *);

   int delta_is_chisqr_distributed;
   char *symbol;

   /* These two pointers are used to impose Lagrange multiplier
    * fit constraints */
   void *constraint_fun;
   Isis_Fit_Statistic_Fun_Type *assigned_fun;

   char *option_string;
   char *message_string;

#ifdef ISIS_FIT_STATISTIC_PRIVATE_DATA
   ISIS_FIT_STATISTIC_PRIVATE_DATA
#endif
};

#ifdef __cplusplus
# define ISIS_FIT_STATISTIC_METHOD(method) \
   extern "C" int Isis_##method##_api_version; \
   int Isis_##method##_api_version = ISIS_API_VERSION; \
   extern "C" Isis_Fit_Statistic_Type *Isis_##method##_stat(void); \
   Isis_Fit_Statistic_Type *Isis_##method##_stat(void)
#else
# define ISIS_FIT_STATISTIC_METHOD(method) \
   int Isis_##method##_api_version = ISIS_API_VERSION; \
   extern Isis_Fit_Statistic_Type *Isis_##method##_stat(void); \
   Isis_Fit_Statistic_Type *Isis_##method##_stat(void)
#endif

/*}}}*/

/*{{{ Fit-method */

typedef int Isis_Fit_Fun_Type(void *clientdata,
                              double *x, unsigned int nbins,
                              double *par, unsigned int npars,
                              double *fx);
struct Isis_Fit_Type
{
   Isis_Fit_Fun_Type *compute_model;
   Isis_Fit_Engine_Type *engine;
   Isis_Fit_Statistic_Type *stat;
   double statistic;
   double *covariance_matrix;   /* npars x npars */
};

extern void isis_fit_free_fit_engine (Isis_Fit_Engine_Type *e);
extern void isis_fit_free_fit_statistic (Isis_Fit_Statistic_Type *s);

/*}}}*/

/*{{{ User-defined fit-engine */

typedef Isis_Fit_Engine_Type *Isis_Fit_Engine_Init_Type (char *, char *);

typedef int Isis_Fit_Routine_Type(Isis_Fit_Type *, void *clientdata,
                                  double *x, double *y, double *weights, unsigned int npts,
                                  double *pars, unsigned int npars);

typedef int Isis_Fit_Range_Hook_Type(void *clientdata, double *par_min, double *par_max, double *par, unsigned int npars);
typedef void Isis_Fit_Verbose_Hook_Type(void *clientdata, double statistic, double *par, unsigned int n);
typedef void Isis_Fit_Warning_Hook_Type(void *clientdata, const char *fmt, ...);

struct Isis_Fit_Engine_Type
{
   Isis_Fit_Engine_Type *next;

   Isis_Fit_Routine_Type *method;
   void (*deallocate)(Isis_Fit_Engine_Type *);
   char *engine_name;
   char *default_statistic_name;
   double *par_min;
   double *par_max;
   double *par_step;
   double *par_relstep;

   int (*set_options) (Isis_Fit_Engine_Type *, Isis_Option_Type *);
   int (*set_range_hook) (Isis_Fit_Engine_Type *, Isis_Fit_Range_Hook_Type *);

   Isis_Fit_Range_Hook_Type *range_hook;
   Isis_Fit_Verbose_Hook_Type *verbose_hook;
   Isis_Fit_Warning_Hook_Type *warn_hook;
   int verbose;

   char *option_string;

#ifdef ISIS_FIT_ENGINE_PRIVATE_DATA
   ISIS_FIT_ENGINE_PRIVATE_DATA
#endif
};

#ifdef __cplusplus
# define ISIS_FIT_ENGINE_METHOD(method,name,sname) \
    extern "C" int Isis_##method##_api_version; \
    int Isis_##method##_api_version = ISIS_API_VERSION; \
    extern "C" Isis_Fit_Engine_Type *Isis_##method##_feng(char *, char *); \
    Isis_Fit_Engine_Type *Isis_##method##_feng(char *name, char *sname)
#else
# define ISIS_FIT_ENGINE_METHOD(method,name,sname) \
    int Isis_##method##_api_version = ISIS_API_VERSION; \
    extern Isis_Fit_Engine_Type *Isis_##method##_feng(char *, char *); \
    Isis_Fit_Engine_Type *Isis_##method##_feng(char *name, char *sname)
#endif

/*}}}*/

enum
{
 ISIS_EVAL_GRID_SEPARATE = 0,
 ISIS_EVAL_GRID_MERGED   = 1,
 ISIS_EVAL_GRID_USER     = 2
};

enum
{
   ISIS_STAT_RESID = 1,
   ISIS_DIFF_RESID,
   ISIS_RATIO_RESID
};
extern int Isis_Residual_Plot_Type;

extern int Isis_Active_Dataset;
extern int Isis_Batch_Mode;
extern int Isis_Silent_Errors;
extern int Isis_Verbose;
extern int Isis_Load_File_Verbose_Mask;

typedef int Isis_Line_Profile_Type (Isis_Hist_t *, double, double, double, int,
                                    double *, int, void *);
#ifdef __cplusplus
# define ISIS_LINE_PROFILE_MODULE(name,grid,flux,wl,atwt,mid,pars,npars,options) \
    extern "C" int Isis_##name##_api_version; \
    int Isis_##name##_api_version = ISIS_API_VERSION; \
    extern "C" int Isis_##name##_line_profile(Isis_Hist_t *, double, double, double, int, double *, int, void *); \
    int Isis_##name##_line_profile (Isis_Hist_t *grid, double flux, double wl, double atwt, int mid, double *pars, int npars, void *options)
#else
# define ISIS_LINE_PROFILE_MODULE(name,grid,flux,wl,atwt,mid,pars,npars,options) \
    int Isis_##name##_api_version; \
    int Isis_##name##_api_version = ISIS_API_VERSION; \
    int Isis_##name##_line_profile(Isis_Hist_t *, double, double, double, int, double *, int, void *); \
    int Isis_##name##_line_profile (Isis_Hist_t *grid, double flux, double wl, double atwt, int mid, double *pars, int npars, void *options)
#endif

extern int isis_system (char *cmd);

#include <stdarg.h>
extern int isis_vsnprintf (char *buf, unsigned int bufsize,
                           const char *fmt, va_list ap);

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#if !(defined(HAVE_ISNAN) || defined(isnan))
# define isnan(x) ((x) != (x))
#endif

#ifndef isfinite
# ifdef HAVE_FINITE
#  define isfinite(x)  finite(x)
# else
#  define isfinite(x)  ((0 == isnan(x)) && (fabs(x) <= DBL_MAX))
# endif
#endif

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif

