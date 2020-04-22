/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2020 Massachusetts Institute of Technology

    Author:  John C. Houck <houck@space.mit.edu>

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

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

/*{{{ Includes */

#include "config.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>
#include <signal.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <slang.h>

#include "isis.h"

/*}}}*/

typedef struct _Xspec_Param_t
{
   union {double *d; float *f;} ear;
   union {double *d; float *f;} param;
   union {double *d; float *f;} photar;
   union {double *d; float *f;} photer;
   int ne;
   int ifl;
   char *filename;
}
Xspec_Param_t;

typedef void Xspec_Fun_t (Xspec_Param_t *p);

/*{{{ Eval */

/*
 * NIST 1998 CODATA recommended values of physical constants:
 */

#define PLANCK      ((double) 6.62606876e-27)       /* Planck's constant (erg s) */
#define CLIGHT      ((double) 2.99792458e10)        /* speed of light (cm/s) */
#define ERG_PER_EV  ((double) 1.602176462e-12)
#define KEV_ANGSTROM (((PLANCK * CLIGHT) / (ERG_PER_EV * 1.e3)) * 1.e8)

#define TOL  (10 * FLT_MIN)

static char *Table_Model_Filename = NULL;
static int Table_Model_Number_Of_Parameters;
static char *Table_Model_Type = NULL;
static char *Xspec_Model_Names_File = NULL;
static int Xspec_Version;

static volatile sig_atomic_t Signal_In_Progress;
static void sig_segv (int signo) /*{{{*/
{
   static char msg[] =
"\n**** XSPEC signal: segmentation fault (SIGSEGV) while in an XSPEC function.\n";

   (void) signo;
   if (Signal_In_Progress)
     return;
   Signal_In_Progress = 1;
   write (STDERR_FILENO, msg, sizeof(msg));
   /* so more signals won't interfere with exit() */
   SLsignal (SIGSEGV, SIG_DFL);
   exit (EXIT_FAILURE);
}

/*}}}*/

static void sig_abrt (int signo) /*{{{*/
{
   static char msg[] =
"\n**** XSPEC signal: abort signal (SIGABRT) generated while in an XSPEC function.\n";

   (void) signo;
   if (Signal_In_Progress)
     return;
   Signal_In_Progress = 1;
   write (STDERR_FILENO, msg, sizeof(msg));
   /* so more signals won't interfere with exit() */
   SLsignal (SIGABRT, SIG_DFL);
   exit (EXIT_FAILURE);
}

/*}}}*/

typedef struct
{
   SLSig_Fun_Type *sig_func;
   SLSig_Fun_Type *sig_func_prev;
   const char *sig_name;
   int sig_type;
}
Signal_Type;
#define SIGNAL_TABLE_END {NULL,NULL,NULL,0}
#define SIGNAL_ENTRY(fun,typ,name) {fun,NULL,name,typ}

static Signal_Type Signal_Table[] =
{
   SIGNAL_ENTRY(sig_segv, SIGSEGV, "segmentation violation"),
   SIGNAL_ENTRY(sig_abrt, SIGABRT, "abort"),
   SIGNAL_TABLE_END
};

static void set_signal_handlers (void)
{
   Signal_Type *st;
   for (st = Signal_Table; st->sig_func != NULL; st++)
     {
        if (SIG_ERR == (st->sig_func_prev = SLsignal (st->sig_type, st->sig_func)))
          {
             fprintf (stderr, "warning: failed initializing signal handler for signal=%d\n",
                      st->sig_type);
          }
     }
}

static void unset_signal_handlers (void)
{
   Signal_Type *st;
   for (st = Signal_Table; st->sig_func != NULL; st++)
     {
        if (SLsignal (st->sig_type, st->sig_func_prev) == SIG_ERR)
          {
             fprintf (stderr, "warning: failed to re-set signal handler for signal=%d\n",
                      st->sig_type);
          }
     }
}

static void call_xspec_fun (Xspec_Fun_t *fun, Xspec_Param_t *p) /*{{{*/
{
   set_signal_handlers ();
   (*fun)(p);
   unset_signal_handlers ();
}

/*}}}*/

typedef struct
{
   union {double *d; float *f;} ebins;
   union {double *d; float *f;} photar;
   int *keep;
   int nbins;
}
Xspec_Info_Type;

#define FREE_XIT(s) \
static void free_##s##_xspec_info_type (Xspec_Info_Type *x) \
{ \
   if (x == NULL) \
     return; \
 \
   ISIS_FREE (x->ebins.s); \
   ISIS_FREE (x->photar.s); \
   ISIS_FREE (x->keep); \
   ISIS_FREE (x); \
}
FREE_XIT(f)
FREE_XIT(d)
#if 0
}
#endif

#define NEW_XIT(s,type) \
static Xspec_Info_Type *new_##s##_xspec_info_type (int nbins) \
{ \
   Xspec_Info_Type *x; \
 \
   if (nbins <= 0) \
     return NULL; \
 \
   if (NULL == (x = (Xspec_Info_Type *) ISIS_MALLOC (sizeof *x))) \
     return NULL; \
 \
   if (NULL == (x->ebins.s = (type *) ISIS_MALLOC ((nbins+1) * sizeof(type))) \
       || NULL == (x->photar.s = (type *) ISIS_MALLOC (nbins * sizeof(type))) \
       || NULL == (x->keep = (int *) ISIS_MALLOC (nbins * sizeof(int)))) \
     { \
        free_##s##_xspec_info_type (x); \
        return NULL; \
     } \
 \
   memset ((char *) x->photar.s, 0, nbins * sizeof(type)); \
   x->nbins = nbins; \
 \
   return x; \
}
NEW_XIT(f,float)
NEW_XIT(d,double)
#if 0
}
#endif

/*
 *   The ISIS data grid might have holes in it, but the XSPEC
 *   grid cannot.  To work around this, I'll generate a hole-free
 *   grid than spans the full range, then pick out the relevant
 *   bin values later.
 *
 *   Also, the input data grid is in Angstrom, but XSPEC wants keV.
 */

#define MAKE_XG(s,type) \
static Xspec_Info_Type *make_##s##_xspec_grid (Isis_Hist_t *g) \
{ \
   Xspec_Info_Type *x; \
   int i, n, k, nbins, n_notice; \
   int *notice_list, *keep; \
   type *ebins; \
   double *bin_lo, *bin_hi; \
 \
   if ((NULL == g) || (g->notice_list == NULL)) \
     { \
        fprintf (stderr, "*** internal error:  got NULL ptr in make_xspec_grid\n"); \
        return NULL; \
     } \
 \
   if (g->n_notice < 1) \
     { \
        fprintf (stderr, "*** no noticed bins\n"); \
        return NULL; \
     } \
 \
   bin_lo = g->bin_lo; \
   bin_hi = g->bin_hi; \
   n_notice = g->n_notice; \
   notice_list = g->notice_list; \
   nbins = 1; \
 \
   for (i=1; i < n_notice; i++) \
     { \
        int n1 = notice_list[i  ]; \
        int n0 = notice_list[i-1]; \
        double diff = fabs (bin_lo[n1] - bin_hi[n0]); \
        double avg  = 0.5 * fabs (bin_lo[n1] + bin_hi[n0]); \
        nbins += ((diff < TOL * avg) ? 1 : 2); \
     } \
 \
   if (NULL == (x = new_##s##_xspec_info_type (nbins))) \
     return NULL; \
 \
   n = g->notice_list[0]; \
 \
   k = nbins; \
   x->ebins.s[k] = (type) KEV_ANGSTROM / g->bin_lo[n]; \
 \
   ebins = x->ebins.s; \
   keep = x->keep; \
 \
   for (i=1; i < n_notice; i++) \
     { \
        int n1 = notice_list[i  ]; \
        int n0 = notice_list[i-1]; \
        double diff = fabs (bin_lo[n1] - bin_hi[n0]); \
        double avg  = 0.5 * fabs (bin_lo[n1] + bin_hi[n0]); \
 \
        k--; \
        ebins[k] = (type) KEV_ANGSTROM / bin_hi[n0]; \
        keep[k]  = 1; \
 \
        if (diff > TOL * avg) \
          { \
             k--; \
             ebins[k] = (type) KEV_ANGSTROM / bin_lo[n1]; \
             keep[k]  = 0; \
          } \
     } \
 \
   k--;             /* low edge of first ENERGY bin */ \
 \
   if (k != 0) \
     { \
        fprintf (stderr, "Invalid xspec grid\n"); \
        free_##s##_xspec_info_type (x); \
        return NULL; \
     } \
 \
   n = g->notice_list[g->n_notice-1]; \
   x->ebins.s[k] = (type) KEV_ANGSTROM / g->bin_hi[n]; \
   x->keep[k]  = 1; \
 \
   return x; \
}
MAKE_XG(f,float)
MAKE_XG(d,double)
#if 0
}
#endif

/*    to unpack the xspec result (on energy grid),
 *    reverse array order consistent with the input wavelength grid
 *
 *    And, apply the normalization.
 *    If norm isnt relevant for this function, the caller should
 *    just set it to 1.0
 */
#define EVAL_XF(s,type) \
static int eval_##s##_xspec_fun (Xspec_Fun_t *fun, double *val, Isis_Hist_t *g, \
                                 type *param, type norm, int category) \
{ \
   Xspec_Param_t p; \
   Xspec_Info_Type *x; \
   int i, k; \
   int ret = -1; \
 \
   if (NULL == (x = make_##s##_xspec_grid (g))) \
     goto finish; \
 \
   p.ear.s = x->ebins.s; \
   p.ne = x->nbins; \
   p.param.s = param; \
   p.ifl = 0; \
   p.photar.s = x->photar.s; \
   if (NULL == (p.photer.s = (type *) ISIS_MALLOC (x->nbins * sizeof(type)))) \
     goto finish;  \
   memset ((char *)p.photer.f, 0, x->nbins * sizeof(type)); \
 \
   p.filename = Table_Model_Filename; \
 \
   if (category == ISIS_FUN_OPERATOR) \
     { \
        k = g->n_notice; \
        for (i=0; i < x->nbins; i++) \
          { \
             if (x->keep[i]) \
               x->photar.s[i] = (type) val[--k]; \
          } \
     } \
 \
   call_xspec_fun (fun, &p); \
 \
   ISIS_FREE (p.photer.s); \
   p.photer.s = NULL; \
 \
   k = g->n_notice; \
   for (i=0; i < x->nbins; i++) \
     { \
        if (x->keep[i]) \
          val[--k] = (double) (norm * x->photar.s[i]); \
     } \
 \
   if (k == 0) \
     ret = 0; \
   else \
     fprintf (stderr, "Inconsistent grid while evaluating XSPEC function\n"); \
 \
   finish: \
 \
   free_##s##_xspec_info_type (x); \
   return ret; \
}
EVAL_XF(f,float)
EVAL_XF(d,double)
#if 0
}
#endif
/*}}}*/

typedef void fptr_type (void);
static fptr_type *Generic_Fptr;
static char *Model_Init_String;

typedef int Hook_Type (double *, Isis_Hist_t *, double *, unsigned int);

typedef struct
{
   char *name;
   fptr_type *symbol;
   char *hook_name;
   char *init_string;
   int malloced;
}
Xspec_Type;
static int Xspec_Type_Id = -1;

typedef void Fcn_f_Type (float *, int *, float *, int *, float *, float *);
typedef void Fcn_fn_Type (float *, int *, float *, int *, float *);
typedef void Fcn_F_Type (double *, int *, double *, int *, double *, double *);
typedef void Fcn_C_Type (double *, int, double *, int, double *, double *, const char *);
typedef void Fcn_c_Type (double *, int, double *, int, double *, double *, const char *);

static void f_sub (Xspec_Param_t *p) /*{{{*/
{
   int ne, ifl;
   if (p == NULL || Generic_Fptr == NULL)
     return;
   ne = p->ne;
   ifl = p->ifl;
   (*((Fcn_f_Type *)Generic_Fptr))(p->ear.f, &ne, p->param.f, &ifl, p->photar.f, p->photer.f);
}

/*}}}*/

static void fn_sub (Xspec_Param_t *p) /*{{{*/
{
   int ne, ifl;
   if (p == NULL || Generic_Fptr == NULL)
     return;
   ne = p->ne;
   ifl = p->ifl;
   (*((Fcn_fn_Type *)Generic_Fptr))(p->ear.f, &ne, p->param.f, &ifl, p->photar.f);
}

/*}}}*/

static void F_sub (Xspec_Param_t *p) /*{{{*/
{
   int ne, ifl;
   if (p == NULL || Generic_Fptr == NULL)
     return;
   ne = p->ne;
   ifl = p->ifl;
   (*((Fcn_F_Type *)Generic_Fptr))(p->ear.d, &ne, p->param.d, &ifl, p->photar.d, p->photer.d);
}

/*}}}*/

static void C_sub (Xspec_Param_t *p) /*{{{*/
{
   int ne, ifl;
   if (p == NULL || Generic_Fptr == NULL)
     return;
   ne = p->ne;
   ifl = p->ifl;
   (*((Fcn_C_Type *)Generic_Fptr))(p->ear.d, ne, p->param.d, ifl, p->photar.d, p->photer.d, Model_Init_String);
}

/*}}}*/

static int mul_f (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   float *param = NULL;
   int ret;

   if (npar > 0)
     {
        unsigned int i;

        if (NULL == (param = (float *) ISIS_MALLOC (npar * sizeof(float))))
          return -1;

        for (i = 0; i < npar; i++)
          param[i] = (float) par[i];
     }

   ret = eval_f_xspec_fun (f_sub, val, g, param, 1.0, ISIS_FUN_ADDMUL);
   ISIS_FREE (param);

   return ret;
}

/*}}}*/

static int con_f (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   float *param = NULL;
   int ret;

   if (npar > 0)
     {
        unsigned int i;

        if (NULL == (param = (float *) ISIS_MALLOC (npar * sizeof(float))))
          return -1;

        for (i = 0; i < npar; i++)
          param[i] = (float) par[i];
     }

   ret = eval_f_xspec_fun (f_sub, val, g, param, 1.0, ISIS_FUN_OPERATOR);
   ISIS_FREE (param);

   return ret;
}

/*}}}*/

static int add_f (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   float tmp[2] = {0.0};
   float *param = tmp;
   unsigned int i;
   int ret, malloced = 0;

   if (npar > 2)
     {
        if (NULL == (param = (float *) ISIS_MALLOC (npar * sizeof(float))))
          return -1;
        malloced = 1;
     }

   for (i = 0; i < npar; i++)
     param[i] = (float) par[i];

   ret = eval_f_xspec_fun (f_sub, val, g, &param[1], param[0], ISIS_FUN_ADDMUL);

   if (malloced)
     ISIS_FREE(param);

   return ret;
}

/*}}}*/

static int mul_fn (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   float *param = NULL;
   int ret;

   if (npar > 0)
     {
        unsigned int i;

        if (NULL == (param = (float *) ISIS_MALLOC (npar * sizeof(float))))
          return -1;

        for (i = 0; i < npar; i++)
          param[i] = (float) par[i];
     }

   ret = eval_f_xspec_fun (fn_sub, val, g, param, 1.0, ISIS_FUN_ADDMUL);
   ISIS_FREE (param);

   return ret;
}

/*}}}*/

static int con_fn (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   float *param = NULL;
   int ret;

   if (npar > 0)
     {
        unsigned int i;

        if (NULL == (param = (float *) ISIS_MALLOC (npar * sizeof(float))))
          return -1;

        for (i = 0; i < npar; i++)
          param[i] = (float) par[i];
     }

   ret = eval_f_xspec_fun (fn_sub, val, g, param, 1.0, ISIS_FUN_OPERATOR);
   ISIS_FREE (param);

   return ret;
}

/*}}}*/

static int add_fn (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   float tmp[2] = {0.0};
   float *param = tmp;
   unsigned int i;
   int ret, malloced = 0;

   if (npar > 2)
     {
        if (NULL == (param = (float *) ISIS_MALLOC (npar * sizeof(float))))
          return -1;
        malloced = 1;
     }

   for (i = 0; i < npar; i++)
     param[i] = (float) par[i];

   ret = eval_f_xspec_fun (fn_sub, val, g, &param[1], param[0], ISIS_FUN_ADDMUL);

   if (malloced)
     ISIS_FREE(param);

   return ret;
}

/*}}}*/

static int mul_F (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   int ret;
   (void) npar;
   ret = eval_d_xspec_fun (F_sub, val, g, par, 1.0, ISIS_FUN_ADDMUL);
   return ret;
}

/*}}}*/

static int con_F (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   int ret;
   (void) npar;
   ret = eval_d_xspec_fun (F_sub, val, g, par, 1.0, ISIS_FUN_OPERATOR);
   return ret;
}

/*}}}*/

static int add_F (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double tmp[2] = {0.0};
   double *param = tmp;
   int ret, malloced = 0;
   unsigned int i;

   if (npar > 2)
     {
        if (NULL == (param = (double *) ISIS_MALLOC (npar * sizeof(double))))
          return -1;
        malloced = 1;
     }

   for (i = 0; i < npar; i++)
     param[i] = par[i];

   ret = eval_d_xspec_fun (F_sub, val, g, &param[1], param[0], ISIS_FUN_ADDMUL);

   if (malloced)
     ISIS_FREE(param);

   return ret;
}

/*}}}*/

static int mul_C (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   int ret;
   (void) npar;
   ret = eval_d_xspec_fun (C_sub, val, g, par, 1.0, ISIS_FUN_ADDMUL);
   return ret;
}

/*}}}*/

static int con_C (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   int ret;
   (void) npar;
   ret = eval_d_xspec_fun (C_sub, val, g, par, 1.0, ISIS_FUN_OPERATOR);
   return ret;
}

/*}}}*/

static int add_C (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double tmp[2] = {0.0};
   double *param = tmp;
   int ret, malloced = 0;
   unsigned int i;

   if (npar > 2)
     {
        if (NULL == (param = (double *) ISIS_MALLOC (npar * sizeof(double))))
          return -1;
        malloced = 1;
     }

   for (i = 0; i < npar; i++)
     param[i] = par[i];

   ret = eval_d_xspec_fun (C_sub, val, g, &param[1], param[0], ISIS_FUN_ADDMUL);

   if (malloced)
     ISIS_FREE(param);

   return ret;
}

/*}}}*/

static int pop_2_matched_arrays (int type, SLang_Array_Type **ap, SLang_Array_Type **bp) /*{{{*/
{
   SLang_Array_Type *a, *b;

   *ap = *bp = NULL;

   if (-1 == SLang_pop_array_of_type (&b, type))
     return -1;

   if (-1 == SLang_pop_array_of_type (&a, type))
     {
        SLang_free_array (b);
        return -1;
     }

   if (a->num_elements == b->num_elements)
     {
        *ap = a;
        *bp = b;
        return 0;
     }

   fprintf (stderr, "*** inconsistent array sizes\n");
   SLang_set_error (Isis_Error);
   SLang_free_array (a);
   SLang_free_array (b);

   return -1;
}

/*}}}*/

static void _xspec_hook (Xspec_Type *xt, Hook_Type *hook) /*{{{*/
{
   SLang_Array_Type *sl_lo, *sl_hi, *sl_par, *sl_val, *sl_arg;
   double *val, *par;
   Isis_Hist_t g;
   int *notice_list = NULL;
   int i, nbins, npar;
   int ret = -1;

   sl_lo = sl_hi = sl_par = sl_val = sl_arg = NULL;

   /* stack should contain  lo, hi, par [, arg] */

   if (hook == con_f || hook == con_fn || hook == con_F
       || hook == con_C)
     {
        if ((-1 == SLang_pop_array_of_type (&sl_arg, SLANG_DOUBLE_TYPE))
             || sl_arg == NULL)
          goto finish;
     }

   if (-1 == SLang_pop_array_of_type (&sl_par, SLANG_DOUBLE_TYPE))
     goto finish;

   if (-1 == pop_2_matched_arrays (SLANG_DOUBLE_TYPE, &sl_lo, &sl_hi))
     goto finish;

   nbins = sl_lo->num_elements;

   if (NULL == (notice_list = (int *) ISIS_MALLOC (nbins * sizeof(int))))
     goto finish;
   for (i = 0; i < nbins; i++)
     notice_list[i] = i;

   sl_val = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nbins, 1);
   if (sl_val == NULL)
     goto finish;

   if (sl_arg != NULL)
     {
        if (sl_arg->num_elements != (unsigned int) nbins)
          {
             fprintf (stderr, "*** inconsistent array size for operator arg\n");
             goto finish;
          }
        memcpy (sl_val->data, sl_arg->data, nbins * sizeof(double));
        SLang_free_array (sl_arg);
     }

   g.bin_lo = (double *)sl_lo->data;
   g.bin_hi = (double *)sl_hi->data;
   g.nbins = nbins;
   g.n_notice = nbins;
   g.notice_list = notice_list;

   val = (double *)sl_val->data;
   par = (double *)sl_par->data;
   npar = sl_par->num_elements;

   /* set global function pointer */
   Generic_Fptr = xt->symbol;
   Model_Init_String = xt->init_string;
   ret = (*hook) (val, &g, par, npar);

   finish:
   if (ret) SLang_set_error (Isis_Error);

   ISIS_FREE(notice_list);
   SLang_free_array (sl_lo);
   SLang_free_array (sl_hi);
   SLang_free_array (sl_par);

   SLang_push_array (sl_val, 1);
}

/*}}}*/

#define XS_HOOK(type) \
static void _xspec_##type##_hook (Xspec_Type *xt) \
{ \
   _xspec_hook (xt, type); \
}
XS_HOOK(add_f)
XS_HOOK(mul_f)
XS_HOOK(con_f)
XS_HOOK(add_fn)
XS_HOOK(mul_fn)
XS_HOOK(con_fn)
XS_HOOK(add_F)
XS_HOOK(mul_F)
XS_HOOK(con_F)
XS_HOOK(add_C)
XS_HOOK(mul_C)
XS_HOOK(con_C)
#if 0
}
#endif

static void _xspec_model_init_string (Xspec_Type *xt, char *init)
{
   if (xt == NULL)
     return;
   ISIS_FREE(xt->init_string);
   if ((init != NULL)
       && (NULL == (xt->init_string = isis_make_string (init))))
     {
        SLang_set_error (Isis_Error);
     }
}

static void handle_link_error (char *path, char *name) /*{{{*/
{
   const char *error = dlerror ();

   if (error == NULL)
     error = "unknown";

   (void) SLang_run_hooks ("_isis->xspec_register_link_error", 3,
                           name, path, error);
}
/*}}}*/

static fptr_type *load_function (char *path, char *name) /*{{{*/
{
   fptr_type *f = NULL;
   void *handle = NULL;

   if (path == NULL || name == NULL)
     return NULL;

#ifndef RTLD_GLOBAL
# define RTLD_GLOBAL 0
#endif
#ifdef RTLD_NOW
# define DLOPEN_FLAG  (RTLD_NOW | RTLD_GLOBAL)
#else
# define DLOPEN_FLAG  (RTLD_LAZY | RTLD_GLOBAL)
#endif

   handle = dlopen (path, DLOPEN_FLAG);
   if (handle == NULL)
     {
        handle_link_error (path, name);
        return NULL;
     }

   f = (fptr_type *) dlsym (handle, name);
   if (f == NULL)
     {
        handle_link_error (path, name);
        /* dlclose (handle); */
        return NULL;
     }

   return f;
}
/*}}}*/

static void load_xspec_fun (char *file, char *fun_name) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   Xspec_Type *xt = NULL;
   fptr_type *fptr;

   if ((NULL == (fptr = (fptr_type *) load_function (file, fun_name)))
       || (NULL == (xt = (Xspec_Type *) ISIS_MALLOC (sizeof(Xspec_Type)))))
     goto push_null;

   memset ((char *)xt, 0, sizeof (*xt));

   xt->symbol = (fptr_type *)fptr;
   xt->init_string = NULL;
   xt->malloced = 1;

   if ((NULL == (mmt = SLang_create_mmt (Xspec_Type_Id, (void *) xt)))
       || (-1 == SLang_push_mmt (mmt)))
     goto push_null;

   return;

push_null:

   ISIS_FREE (xt);
   SLang_free_mmt (mmt);
   SLang_push_null();
}

/*}}}*/

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
}
#endif

#if defined(HAVE_XSPEC_11)
# include "_model_externs_xspec11.inc"
# define XSPEC_MODEL_NAMES_FILE "_names_xspec11.dat"
# define XSPEC_VERSION 11
#elif defined(HAVE_XSPEC_12)
# include "_model_externs_xspec12.inc"
# define XSPEC_MODEL_NAMES_FILE "_names_xspec12.dat"
# define XSPEC_VERSION 12
#endif

static Xspec_Type Static_Fun_Table[] =
{
#if defined(HAVE_XSPEC_11)
#  include "_model_table_xspec11.inc"
#elif defined(HAVE_XSPEC_12)
#  include "_model_table_xspec12.inc"
#endif
     {NULL, NULL, NULL, NULL, 0}
};

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

static void find_xspec_fun (char *fun_name) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   Xspec_Type *xt = NULL;

   /* silence compiler warnings about unused symbols */
   (void) mul_C;  (void) con_C;  (void) add_C;
   (void) mul_F;  (void) con_F;  (void) add_F;

   for (xt = Static_Fun_Table; xt->name != NULL; xt++)
     {
        if (0 == strcmp (fun_name, xt->name))
          break;
     }

   if (xt->name == NULL)
     {
        SLang_push_null ();
        return;
     }

   if (NULL == (mmt = SLang_create_mmt (Xspec_Type_Id, (void *) xt)))
     {
        SLang_push_null ();
        return;
     }

   if (-1 == SLang_push_mmt (mmt))
     {
        SLang_free_mmt (mmt);
        return;
     }
}

/*}}}*/

/* Intrinsics */

/* DUMMY_MODEL_TYPE is a temporary hack that will be modified to the true
  * id once the interpreter provides it when the class is registered.  See below
  * for details.  The reason for this is simple: for a module, the type-id
  * must be assigned dynamically.
  */
#define DUMMY_MODEL_TYPE   255
#define MT DUMMY_MODEL_TYPE
#define I SLANG_INT_TYPE
#define R SLANG_REF_TYPE
#define S SLANG_STRING_TYPE
#define V SLANG_VOID_TYPE

static SLang_Intrin_Fun_Type Intrinsics [] =
{
   MAKE_INTRINSIC_2("load_xspec_fun", load_xspec_fun, V, S, S),
   MAKE_INTRINSIC_1("find_xspec_fun", find_xspec_fun, V, S),
   MAKE_INTRINSIC_1("_xspec_mul_fn_hook", _xspec_mul_fn_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_add_fn_hook", _xspec_add_fn_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_con_fn_hook", _xspec_con_fn_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_mul_f_hook", _xspec_mul_f_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_add_f_hook", _xspec_add_f_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_con_f_hook", _xspec_con_f_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_mul_F_hook", _xspec_mul_F_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_add_F_hook", _xspec_add_F_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_con_F_hook", _xspec_con_F_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_mul_C_hook", _xspec_mul_C_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_add_C_hook", _xspec_add_C_hook, V, MT),
   MAKE_INTRINSIC_1("_xspec_con_C_hook", _xspec_con_C_hook, V, MT),
   MAKE_INTRINSIC_2("_xspec_model_init_string", _xspec_model_init_string, V, MT, S),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef I
#undef R
#undef S
#undef V
#undef MT

static void free_xspec_fun_type (SLtype type, void *f) /*{{{*/
{
   Xspec_Type *xt = (Xspec_Type *)f;
   (void) type;
   (void) xt;
   if (xt == NULL)
     return;
   ISIS_FREE(xt->init_string);
   if (xt->malloced)
     ISIS_FREE(xt);
}

/*}}}*/

#include "xsFortran.h"

/* In the xspec source code,
 * fpdatd is defined in xspec/src/functions/fpunc.f
 * fgdatd is defined in xspec/src/functions/fgunc.f
 */

static int xs_set_datadir (char *name) /*{{{*/
{
   if (name == NULL)
     return -1;
   FPDATD(name);
   return 0;
}

/*}}}*/

static void xs_get_datadir (void) /*{{{*/
{
   SLang_push_string (FGDATD());
}

/*}}}*/

/* In the xspec source code,
 * fpsolr is defined in xspec/src/functions/fpunc.f
 * fgsolr is defined in xspec/src/functions/fgunc.f
 */

static int xs_set_abundance_table (char *name) /*{{{*/
{
   int ierr = 0;
   if (name == NULL)
     return -1;
   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     SLdo_pop();
   else
     {
        SLang_Array_Type *sl_abund = NULL;
        int len = 30;
        if ((-1 == SLang_pop_array_of_type (&sl_abund, SLANG_FLOAT_TYPE))
            || (sl_abund == NULL))
          return -1;
        if (sl_abund->num_elements != (unsigned int) len)
          {
             fprintf (stderr, "*** Error:  abundance array length=%d (should be %d)\n",
                      sl_abund->num_elements, len);
             SLang_free_array (sl_abund);
             return -1;
          }
        FPSLFL((float *)sl_abund->data,len,&ierr);
        SLang_free_array (sl_abund);
        if (ierr)
          {
             fprintf (stderr, "*** Error:  failed initializing new abundance table\n");
             return -1;
          }
        name = "file";
     }
   /* FPSOLR modifies ierr on return */
   FPSOLR(name, &ierr);
   return ierr ? -1 : 0;
}

/*}}}*/

static void xs_get_abundance_table (void) /*{{{*/
{
   SLang_push_string (FGSOLR());
}

/*}}}*/

/* In the xspec source code,
 * fpxsct is defined in xspec/src/functions/fpunc.f
 * fgxsct is defined in xspec/src/functions/fgunc.f
 */

static int xs_set_xsection_table (char *name) /*{{{*/
{
   int ierr = -1;
   if (name == NULL)
     return -1;
   /* FPXSCT modifies ierr on return */
   FPXSCT(name, &ierr);
   return ierr ? -1 : 0;
}

/*}}}*/

static void xs_get_xsection_table (void) /*{{{*/
{
   SLang_push_string (FGXSCT());
}

/*}}}*/

static void xs_fpmstr (char *p, char *v) /*{{{*/
{
   if (p == NULL || v == NULL)
     return;
   FPMSTR(p, v);
}

/*}}}*/

#define XSPEC_DEFAULT_H0 70.0
#define XSPEC_DEFAULT_Q0 0.0
#define XSPEC_DEFAULT_L0 0.73

static void xs_set_cosmo_hubble (float *h0)
{
   float h = h0 ? *h0 : XSPEC_DEFAULT_H0;
   csmph0(h);
}
static void xs_set_cosmo_decel (float *q0)
{
   float q = q0 ? *q0 : XSPEC_DEFAULT_Q0;
   csmpq0(q);
}
static void xs_set_cosmo_lambda (float *l0)
{
   float l = l0 ? *l0 : XSPEC_DEFAULT_L0;
   csmpl0(l);
}

static double xs_get_cosmo_hubble (void)
{
   return (double) csmgh0();
}
static double xs_get_cosmo_decel (void)
{
   return (double) csmgq0();
}
static double xs_get_cosmo_lambda (void)
{
   return (double) csmgl0();
}

static double xs_get_element_solar_abundance (char *element) /*{{{*/
{
   double x = FGABND(element);
   return x;
}

/*}}}*/

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
}
#endif

#define PHOTO_FC FC_FUNC(photo,PHOTO)
extern float PHOTO_FC(float *kev1, float *kev2, int *Z, int *versn, long *status);

#define GPHOTO_FC FC_FUNC(gphoto,GPHOTO)
extern float GPHOTO_FC(float *kev1, float *kev2, int *Z, long *status);

#define PHFIT2_FC FC_FUNC(phfit2,PHFIT2)
extern float PHFIT2_FC(int *nz, int *ne, int *is, float *e, float *s);

/* As of heasoft-6.16 (July 2014), these symbols are no longer exported.
 * These interfaces can be updated if someone's using them, otherwise,
 * they can eventually be removed.
 */
#undef OBSOLETE_NEI_SYMBOLS
#ifdef OBSOLETE_NEI_SYMBOLS
#define INITNEI_FC FC_FUNC(initnei,INITNEI)
extern void INITNEI_FC(int *ni, int *nz);

#define IONSNEQR_FC FC_FUNC(ionsneqr,IONSNEQR)
extern void IONSNEQR_FC(float *tmp, float *tau, int *n, int*nzmax, int *nionp,
                        float *fout, int *ionel, int *ionstage);

#define NONEQ_FC FC_FUNC(noneq,NONEQ)
extern void NONEQ_FC(float *tempr, int *ntp, float *tau, int *n, float *weight,
                     int *nzmax, int *nzpmax, int *nionp, double *vr, double *vl,
                     double *eig, double *feqb, double *feqs, double *fs,
                     double *fspec, double *work, int *lt, float *fout,
                     int *ionel, int *ionstage);
#endif

/* To pass a C string to Fortran:  for each string, append to the Fortran
 * function's parameter list a 'long' containing the length of the string.
 * This works with gcc/gfortran, but may not be totally portable.
 */
#ifdef XSPEC_OLDTABLE

#define XSPEC11_TABLE_FUN(name,xsname,XSNAME)                              \
   extern void FC_FUNC(xsname,XSNAME)(float *,int *,float *,char *,int *,float *,float *,long); \
   static void name (Xspec_Param_t *p)                                         \
   {                                                                           \
      FC_FUNC(xsname,XSNAME)(p->ear.f,&p->ne,p->param.f,p->filename,&p->ifl,p->photar.f,p->photer.f,(long)strlen(p->filename));   \
   }

XSPEC11_TABLE_FUN(xs_atbl,xsatbl,XSATBL)
XSPEC11_TABLE_FUN(xs_mtbl,xsmtbl,XSMTBL)

#else

/* Since the tabint and tabint_ functions are not exported
 * by xspec from the libXSFunction library the compile can not
 * report type missmatches.
 * The current version of heasoft v6.27.1 does change void tabint
 * to extern "C" void tabint so better use this directly than that string + long
 * hack for fortran.
 */
#define XSPEC12_TABLE_FUN(name,xsname,XSNAME)                              \
   extern void xsname (float *,int, float *, int, char *,int, char *, float *,float *); \
   static void name (Xspec_Param_t *p)                                         \
   {                                                                           \
     int npar = Table_Model_Number_Of_Parameters;\
     xsname (p->ear.f, p->ne,p->param.f, npar, p->filename, p->ifl, Table_Model_Type, p->photar.f,p->photer.f);   \
    }

XSPEC12_TABLE_FUN(xs_atbl,tabint,TABINT)
XSPEC12_TABLE_FUN(xs_mtbl,tabint,TABINT)

#endif

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

static double xs_photo (float *kev1, float *kev2, int *Z, int *versn) /*{{{*/
{
   long status = 0;
   double xsect = (double) PHOTO_FC(kev1,kev2,Z,versn,&status);
   if (status)
     {
        SLang_set_error(Isis_Error);
        xsect = 0.0;
     }
   return xsect;
}

/*}}}*/

static double xs_gphoto (float *kev1, float *kev2, int *Z) /*{{{*/
{
   long status = 0;
   double xsect = (double) GPHOTO_FC(kev1,kev2,Z,&status);
   if (status)
     {
        SLang_set_error(Isis_Error);
        xsect = 0.0;
     }
   return xsect;
}

/*}}}*/

static double xs_phfit2 (int *nz, int *ne, int *is, float *e) /*{{{*/
{
   float s;
   PHFIT2_FC(nz, ne, is, e, &s);
   return (double)s;
}

/*}}}*/

#ifdef OBSOLETE_NEI_SYMBOLS
static void xs_initnei (int *nionp,int *nzmax) /*{{{*/
{
   static int nei_is_initialized = 0;
   static int ni, nz;
   if (nei_is_initialized == 0)
     {
        INITNEI_FC(&ni,&nz);
        nei_is_initialized = 1;
     }
   *nionp = ni;
   *nzmax = nz;
}

/*}}}*/

static void xs_ionsneqr(void) /*{{{*/
{
   SLang_Array_Type *sl_tmp=NULL, *sl_tau=NULL, *sl_fout=NULL;
   SLang_Array_Type *sl_ionel=NULL, *sl_ionstage=NULL;
   float *tmp, *tau, *fout;
   int *ionel, *ionstage;
   int n, nzmax, nionp;

   xs_initnei(&nionp, &nzmax);

   if (-1 == pop_2_matched_arrays (SLANG_FLOAT_TYPE, &sl_tmp, &sl_tau))
     return;

   n = sl_tmp->num_elements;

   sl_fout = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &nionp, 1);
   sl_ionel= SLang_create_array (SLANG_INT_TYPE, 0, NULL, &nionp, 1);
   sl_ionstage= SLang_create_array (SLANG_INT_TYPE, 0, NULL, &nionp, 1);
   if ((sl_fout == NULL) || (sl_ionel == NULL) || (sl_ionstage == NULL))
     {
        SLang_set_error (Isis_Error);
        goto push_results;
     }

   tmp = (float *)sl_tmp->data;
   tau = (float *)sl_tau->data;
   fout = (float *)sl_fout->data;
   ionel = (int *)sl_ionel->data;
   ionstage = (int *)sl_ionstage->data;

   IONSNEQR_FC(tmp, tau, &n, &nzmax, &nionp, fout, ionel, ionstage);

   push_results:
   SLang_free_array (sl_tmp);
   SLang_free_array (sl_tau);
   SLang_push_array (sl_fout,1);
   SLang_push_array (sl_ionel,1);
   SLang_push_array (sl_ionstage,1);
}

/*}}}*/

static void xs_noneq (void) /*{{{*/
{
   SLang_Array_Type *sl_tempr=NULL, *sl_tau=NULL, *sl_weight=NULL;
   SLang_Array_Type *sl_fout=NULL, *sl_ionel=NULL, *sl_ionstage=NULL;
   float *tempr=NULL, *tau=NULL, *weight=NULL, *fout=NULL;
   double *vr=NULL, *vl=NULL, *eig=NULL, *feqb=NULL, *feqs=NULL, *fs=NULL,
     *fspec=NULL, *work=NULL;
   double *feqs_p1=NULL, *fs_p1=NULL;
   int *ionel=NULL, *ionstage=NULL, *lt=NULL;
   int nionp, nzmax, nzpmax, ntp, n;

   if (-1 == pop_2_matched_arrays (SLANG_FLOAT_TYPE, &sl_tau, &sl_weight))
     goto push_results;

   if (-1 == SLang_pop_array_of_type (&sl_tempr, SLANG_FLOAT_TYPE))
     goto push_results;

   ntp = sl_tempr->num_elements;
   n = sl_tau->num_elements;
   xs_initnei (&nionp, &nzmax);
   nzpmax = nzmax + 1;

   sl_fout = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &nionp, 1);
   sl_ionel= SLang_create_array (SLANG_INT_TYPE, 0, NULL, &nionp, 1);
   sl_ionstage= SLang_create_array (SLANG_INT_TYPE, 0, NULL, &nionp, 1);
   if ((sl_fout == NULL) || (sl_ionel == NULL) || (sl_ionstage == NULL))
     {
        SLang_set_error (Isis_Error);
        goto push_results;
     }

   tempr = (float *)sl_tempr->data;
   tau = (float *)sl_tau->data;
   weight = (float *)sl_weight->data;
   fout = (float *)sl_fout->data;
   ionel = (int *)sl_ionel->data;
   ionstage = (int *)sl_ionstage->data;

   memset ((char *)fout, 0, nionp*sizeof(float));

   if (   (NULL == (vr    = (double *) ISIS_MALLOC (sizeof(double) * nzmax*nzmax*ntp)))
       || (NULL == (vl    = (double *) ISIS_MALLOC (sizeof(double) * nzmax*nzmax*ntp)))
       || (NULL == (eig   = (double *) ISIS_MALLOC (sizeof(double) * nzmax*ntp)))
       || (NULL == (feqb  = (double *) ISIS_MALLOC (sizeof(double) * nzpmax*ntp)))
       || (NULL == (feqs  = (double *) ISIS_MALLOC (sizeof(double) * (nzpmax + 1))))
       || (NULL == (fs    = (double *) ISIS_MALLOC (sizeof(double) * (nzpmax + 1))))
       || (NULL == (fspec = (double *) ISIS_MALLOC (sizeof(double) * nzmax)))
       || (NULL == (work  = (double *) ISIS_MALLOC (sizeof(double) * nzmax)))
       || (NULL == (lt    = (int *) ISIS_MALLOC (sizeof(int) * ntp)))
      )
     {
        SLang_set_error (Isis_Error);
        goto push_results;
     }

   /* The fortran code will access feqs(0) and fs(0) */
   feqs_p1 = feqs + 1;
   fs_p1 = fs + 1;

   NONEQ_FC(tempr,&ntp,tau,&n,weight,&nzmax,&nzpmax,&nionp,vr,vl,
            eig,feqb,feqs_p1,fs_p1,fspec,work,lt,fout,ionel,ionstage);

push_results:
   ISIS_FREE(vr);
   ISIS_FREE(vl);
   ISIS_FREE(eig);
   ISIS_FREE(feqb);
   ISIS_FREE(feqs);
   ISIS_FREE(fs);
   ISIS_FREE(fspec);
   ISIS_FREE(work);
   ISIS_FREE(lt);

   SLang_free_array (sl_tempr);
   SLang_free_array (sl_tau);
   SLang_free_array (sl_weight);

   SLang_push_array (sl_fout,1);
   SLang_push_array (sl_ionel,1);
   SLang_push_array (sl_ionstage,1);
}

/*}}}*/

#endif

static int xs_gchat (void)
{
#ifdef HAVE_XSPEC_12
   return FGCHAT();
#else
   return 0;
#endif
}

static void xs_pchat (int *lev)
{
#ifdef HAVE_XSPEC_12
   FPCHAT(*lev);
#endif
}

static int xs_init (void)
{
   float h0=XSPEC_DEFAULT_H0;
   float q0=XSPEC_DEFAULT_Q0;
   float l0=XSPEC_DEFAULT_L0;
   FNINIT();
   csmph0(h0);
   csmpq0(q0);
   csmpl0(l0);
   return 0;
}

/* XSPEC table models */

static void set_table_model_filename (char *filename) /*{{{*/
{
   char *t;

   if (filename == NULL)
     {
        fputs ("*** error: filename not set", stderr);
        return;
     }

   if (NULL == (t = (char *) ISIS_MALLOC (1 + strlen(filename))))
     {
        fputs ("*** error: malloc failed", stderr);
        return;
     }

   strcpy (t, filename);
   ISIS_FREE (Table_Model_Filename);
   Table_Model_Filename = t;
}

/*}}}*/

static void set_table_model_number_of_parameters (int *npar) /*{{{*/
{
  if (npar == NULL)
  {
    fputs ("*** error: number of parameters not set", stderr);
    return;
  }

  Table_Model_Number_Of_Parameters = *npar;
}
/* }}} */

static void set_table_model_type (char *tabtype) /*{{{*/
{
   char *t;

   if (tabtype == NULL)
     {
        fputs ("*** error: table type not set", stderr);
        return;
     }

   if (NULL == (t = (char *) ISIS_MALLOC (1 + strlen(tabtype))))
     {
        fputs ("*** error: malloc failed", stderr);
        return;
     }

   strcpy (t, tabtype);
   ISIS_FREE (Table_Model_Type);
   Table_Model_Type = t;
}
/* }}} */

static int evaluate_table_model (Xspec_Fun_t *fun) /*{{{*/
{
   Isis_Hist_t g;
   SLang_Array_Type *sl_lo, *sl_hi, *sl_val, *sl_par;
   double *val = NULL;
   float *param = NULL;
   int *notice_list = NULL;
   int *notice = NULL;
   int i, nbins, ret = -1;

   sl_lo = sl_hi = sl_val = sl_par = NULL;

   if (Table_Model_Filename == NULL)
     {
        fprintf (stderr, "Internal error in xspec module - table model filename not set\n");
        return -1;
     }

   if (-1 == SLang_pop_array_of_type (&sl_par, SLANG_FLOAT_TYPE)
       ||-1 == SLang_pop_array_of_type (&sl_hi, SLANG_DOUBLE_TYPE)
       ||-1 == SLang_pop_array_of_type (&sl_lo, SLANG_DOUBLE_TYPE)
       || (sl_par == NULL) || (sl_hi == NULL) || (sl_lo == NULL))
     goto finish;

   nbins = sl_lo->num_elements;
   if (nbins != (int) sl_hi->num_elements)
     goto finish;

   if ((NULL == (notice_list = (int *) ISIS_MALLOC (nbins * sizeof (int))))
       || (NULL == (notice = (int *) ISIS_MALLOC (nbins * sizeof (int))))
       || (NULL == (val = (double *) ISIS_MALLOC (nbins * sizeof (double)))))
     goto finish;

   for (i = 0; i < nbins; i++)
     {
        notice_list[i] = i;
        notice[i] = 1;
     }

   g.val = val;
   g.bin_lo = (double *)sl_lo->data;
   g.bin_hi = (double *)sl_hi->data;
   g.nbins = nbins;
   g.n_notice = nbins;
   g.notice = notice;
   g.notice_list = notice_list;

   param = (float *)sl_par->data;

   ret = eval_f_xspec_fun (fun, val, &g, param, 1.0, ISIS_FUN_ADDMUL);

   finish:
   SLang_free_array (sl_par);
   SLang_free_array (sl_hi);
   SLang_free_array (sl_lo);
   ISIS_FREE (notice_list);
   ISIS_FREE (notice);

   sl_val = SLang_create_array (SLANG_DOUBLE_TYPE, 0, val, &nbins, 1);
   SLang_push_array (sl_val, 1);

   return ret;
}

/*}}}*/

/*}}}*/

static int atbl (void) /*{{{*/
{
   return evaluate_table_model (xs_atbl);
}

/*}}}*/

static int mtbl (void) /*{{{*/
{
   return evaluate_table_model (xs_mtbl);
}

/*}}}*/

static SLang_Intrin_Fun_Type Table_Model_Intrinsics [] =
{
   MAKE_INTRINSIC_S("_set_table_model_filename", set_table_model_filename, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_S("_set_table_model_type", set_table_model_type, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_I("_set_table_model_number_of_parameters", set_table_model_number_of_parameters, SLANG_VOID_TYPE),
   MAKE_INTRINSIC("_atbl", atbl, SLANG_VOID_TYPE, 0),
   MAKE_INTRINSIC("_mtbl", mtbl, SLANG_VOID_TYPE, 0),
   SLANG_END_INTRIN_FUN_TABLE
};

#ifndef HEADAS
  #define HEADAS "xxx"
#endif
static char *Compiled_Headas_Path = HEADAS ;

#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define F SLANG_FLOAT_TYPE
#define D SLANG_DOUBLE_TYPE

static SLang_Intrin_Fun_Type Private_Intrinsics [] =
{
   MAKE_INTRINSIC_S("_xs_set_datadir", xs_set_datadir, I),
   MAKE_INTRINSIC_0("_xs_get_datadir", xs_get_datadir, V),
   MAKE_INTRINSIC_S("_xs_set_abundances", xs_set_abundance_table, I),
   MAKE_INTRINSIC_0("_xs_get_abundances", xs_get_abundance_table, V),
   MAKE_INTRINSIC_S("_xs_set_xsections", xs_set_xsection_table, I),
   MAKE_INTRINSIC_0("_xs_get_xsections", xs_get_xsection_table, V),
   MAKE_INTRINSIC_4("_xs_photo", xs_photo, D, F,F,I,I),
   MAKE_INTRINSIC_3("_xs_gphoto", xs_gphoto, D, F,F,I),
   MAKE_INTRINSIC_4("_xs_phfit2", xs_phfit2, D, I,I,I,F),
#ifdef OBSOLETE_NEI_SYMBOLS
   MAKE_INTRINSIC_0("_xs_ionsneqr", xs_ionsneqr, V),
   MAKE_INTRINSIC_0("_xs_noneq", xs_noneq, V),
#endif
   MAKE_INTRINSIC_S("_xs_get_element_solar_abundance", xs_get_element_solar_abundance, D),
   MAKE_INTRINSIC_2("_xs_fpmstr", xs_fpmstr, V, S, S),
   MAKE_INTRINSIC_1("_xs_set_cosmo_hubble", xs_set_cosmo_hubble, V, F),
   MAKE_INTRINSIC_1("_xs_set_cosmo_decel", xs_set_cosmo_decel, V, F),
   MAKE_INTRINSIC_1("_xs_set_cosmo_lambda", xs_set_cosmo_lambda, V, F),
   MAKE_INTRINSIC_0("_xs_get_cosmo_hubble", xs_get_cosmo_hubble, D),
   MAKE_INTRINSIC_0("_xs_get_cosmo_decel", xs_get_cosmo_decel, D),
   MAKE_INTRINSIC_0("_xs_get_cosmo_lambda", xs_get_cosmo_lambda, D),
   MAKE_INTRINSIC_1("_xs_pchat", xs_pchat, V, I),
   MAKE_INTRINSIC_0("_xs_gchat", xs_gchat, I),
   SLANG_END_INTRIN_FUN_TABLE
};

static char *Fc_Mangle_Suffix = FC_MANGLE_SUFFIX;

static SLang_Intrin_Var_Type Private_Vars [] =
{
   MAKE_VARIABLE("Xspec_Compiled_Headas_Path", &Compiled_Headas_Path, S, 1),
   MAKE_VARIABLE("Xspec_Model_Names_File", &Xspec_Model_Names_File, S, 1),
   MAKE_VARIABLE("Xspec_Version", &Xspec_Version, I, 1),
   MAKE_VARIABLE("__FC_MANGLE_SUFFIX", &Fc_Mangle_Suffix, S, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_IConstant_Type Private_IConst [] =
{
   MAKE_ICONSTANT("__FC_MANGLE_UPCASE", FC_MANGLE_UPCASE),
   SLANG_END_ICONST_TABLE
};

#undef V
#undef I
#undef S
#undef F
#undef D

/*}}}*/

/*}}}*/

/* init */

static char *Xanadu_Setenv = NULL;
static char *Headas_Setenv = NULL;

static void free_env (void)
{
   ISIS_FREE(Xanadu_Setenv);
   ISIS_FREE(Headas_Setenv);
}

static char *copy_and_set_env (char *env_name, char *env_builtin_value) /*{{{*/
{
   struct stat st;
   char *env, *env_set;

   env = getenv (env_name);

   if (env != NULL)
     {
        if (-1 == stat (env, &st))
          {
             fprintf (stderr, "*** %s environment variable provides an invalid path.\n", env_name);
             fprintf (stderr, "    Falling back to compiled-in path %s=%s\n",
                     env_name, env_builtin_value);
             env = env_builtin_value;
          }
     }
   else env = env_builtin_value;

   if ((env == env_builtin_value)
       && (-1 == stat (env, &st)))
     {
        fprintf (stderr, "*** Invalid path: %s=%s\n", env_name, env);
        return NULL;
     }

   if (NULL == (env_set = isis_mkstrcat (env_name, "=", env, NULL)))
     return NULL;

   if (-1 == putenv (env_set))
     {
        fprintf (stderr, "Failed setting %s environment variable: %s\n", env_name, env_set);
        ISIS_FREE(env_set);
        return NULL;
     }

   return env_set;
}

/*}}}*/

#ifdef __cplusplus
extern "C"
#endif
void deinit_xspec_module (void);
void deinit_xspec_module (void) /*{{{*/
{
   ISIS_FREE (Table_Model_Filename);
   free_env();
}

/*}}}*/

SLANG_MODULE(xspec);
int init_xspec_module_ns (char *ns_name) /*{{{*/
{
   SLang_Class_Type *cl;
   SLang_NameSpace_Type *ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return -1;

   if (Xspec_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("Xspec_Type")))
          return -1;

        (void) SLclass_set_destroy_function (cl, free_xspec_fun_type);

        /* By registering as SLANG_VOID_TYPE, slang will dynamically allocate a
         * type.
         */
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (Xspec_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return -1;

        Xspec_Type_Id = SLclass_get_class_id (cl);
        SLclass_patch_intrin_fun_table1 (Intrinsics, DUMMY_MODEL_TYPE, Xspec_Type_Id);
     }

   if (-1 == SLns_add_intrin_fun_table (NULL, Intrinsics, NULL))
     return -1;

#ifdef HAVE_XSPEC_11
   if (NULL == (Xanadu_Setenv = copy_and_set_env ("XANADU", HEADAS "/..")))
     goto return_error;
#endif

   if (NULL == (Headas_Setenv = copy_and_set_env ("HEADAS", HEADAS)))
     goto return_error;

   if (-1 == SLns_add_intrin_fun_table (NULL, Table_Model_Intrinsics, "__HAVE_XSPEC_TABLE_MODELS__"))
     {
        fprintf (stderr, "Failed initializing XSPEC table-model intrinsics\n");
        goto return_error;
     }

   Xspec_Model_Names_File = XSPEC_MODEL_NAMES_FILE;
   Xspec_Version = XSPEC_VERSION;

   if ((-1 == SLns_add_intrin_fun_table (ns, Private_Intrinsics, NULL))
      || (-1 == SLns_add_intrin_var_table (ns, Private_Vars, NULL))
      || (-1 == SLns_add_iconstant_table (ns, Private_IConst, NULL)))
     {
        fprintf (stderr, "Failed initializing XSPEC intrinsics\n");
        goto return_error;
     }

   if (-1 == xs_init())
     {
        fprintf (stderr, "Failed initializing XSPEC module\n");
        goto return_error;
     }

   (void) SLdefine_for_ifdef ("__XSPEC__");

#ifdef WITH_XSPEC_STATIC_LINKED
   (void) SLdefine_for_ifdef ("__XSPEC_STATIC_LINKED__");
#endif

#ifdef HAVE_XSPEC_12
   (void) SLdefine_for_ifdef ("__HAVE_XSPEC_12__");
#endif

   return 0;

   return_error:
   deinit_xspec_module();
   return -1;
}

/*}}}*/

