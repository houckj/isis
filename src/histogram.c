/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2012  Massachusetts Institute of Technology

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

/* $Id: histogram.c,v 1.52 2004/05/20 00:27:44 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <slang.h>

#include "isis.h"
#include "_isis.h"
#include "cfits.h"
#include "util.h"
#include "plot.h"
#include "histogram.h"
#include "rmf.h"
#include "arf.h"
#include "keyword.h"
#include "errors.h"

/*}}}*/

#undef TEST_MODEL_NOTICE

/*{{{ types */

Hist_t *Hist_Current;
int Isis_Active_Dataset;
int Isis_Residual_Plot_Type = ISIS_STAT_RESID;
int Hist_Ignore_PHA_Response_Keywords;
int Hist_Ignore_PHA_Backfile_Keyword;
int Hist_Allow_Multiple_Arf_Factors;
int Hist_Warn_Invalid_Uncertainties;

/* array must be NULL terminated */
static const char *Spectrum_Hdu_Names[] = {"SPECTRUM", NULL};
static const char *Spectrum_Hdu_Names_Hook = "_nonstandard_spectrum_hdu_names";

int Isis_List_Filenames = 1;
double Hist_Min_Stat_Err = -1.0;        /* used only if > 0 */
double Hist_Min_Model_Spacing = 0.0;
double Hist_Rmf_Grid_Match_Tol = 1.0e-4;
#define RMF_GRID_TOL  fabs(Hist_Rmf_Grid_Match_Tol)

static double get_assigned_exposure (Hist_t *h);
static int get_exposure_time (Hist_t *h, double *t);
static Isis_Rmf_t *matching_identity_rmf (Hist_t *h, Isis_Arf_t *a);
static int rebin_flux_using_weights (Hist_t *h, double *, double *, double *, double *);
static int run_stat_error_hook (Hist_t *h);
static int copy_d_array (double **to, double *from, int nbins);
static int load_background_from_file (Hist_t *h, char *file);

typedef struct
{
   union
     {
        double *v;
        double s;
     } value;
   int is_vector;
   unsigned int n;
}
Area_Type;

static int rebin_backscale (Area_Type *at, Hist_t *h, double *lo, double *hi, int nbins);

static void area_free (Area_Type *a);

enum
{
   PHA_TYPE_I = 1,
   PHA_TYPE_II = 2,
   MIN_WIDTH = 11
};

/* combo_id=0 means 'not combined with any other dataset' */
static unsigned int Next_Combo_Id = 1;
static unsigned int Next_Eval_Grid_Id = 0;

struct _Hist_t
{
   Hist_t *next;
   int index;                    /* id number */

   unsigned int combo_id;        /* for combining datasets to improve statistics */
   double combo_weight;
   SLang_Name_Type *pre_combine;

   int exclude;                  /* != 0 means exclude from fit */
   Isis_Rsp_t f_rsp;             /* fit-response (used in fit) */
   Isis_Rsp_t a_rsp;             /* assigned response */
   Isis_Kernel_t *kernel;        /* fit kernel */

   SLang_Name_Type *post_model_hook;
   void (*post_model_hook_delete)(SLang_Name_Type *);

   SLang_Name_Type *instrumental_background_hook;          /* instrumental background */
   char *instrumental_background_hook_name;

   SLang_Name_Type *assigned_model;
   Isis_Arg_Type *assigned_model_args;

   SLang_Name_Type *stat_error_hook;
   void (*stat_error_hook_delete)(SLang_Name_Type *);

   SLang_Any_Type *user_meta;    /* user-defined metadata */

   char *bgd_file;               /* name of background file */
   char *file;                   /* name of data file */

   double min_stat_err;

   /* [R] marks things directly or indirectly affected by re-binning */

   double *bin_lo;               /* [R] bin left edge */
   double *bin_hi;               /* [R] bin right edge */
   double *counts;               /* [R] number of counts in bin */
   double *stat_err;             /* [R] statistical error in number of counts */
   double *sys_err_frac;         /* [R] systematic err is (sys_err_frac * data) */
   double *flux;                 /* [R] flux corrected data (photons/s/cm^2/bin) */
   double *flux_err;             /* [R] statistical error in flux corrected data */
   double *model_counts;         /* [R] number of model counts in bin */
   double *convolved_model_flux; /* [R] \int dE R(h,E)S(E) */
   int nbins;                    /* [R] number of data bins */

   Isis_Hist_t model_flux;       /* [R'] \int dE S(E) = model flux in bin (photons/s/cm^2) */
   double *scaled_bgd;           /* area-, exposure-scaled background */

   /* as-input spectrum (not re-binned) */

   double *orig_bin_lo;          /* bin left edge */
   double *orig_bin_hi;          /* bin right edge */
   double *orig_counts;          /* number of counts in bin */
   double *orig_stat_err;        /* statistical error in number of counts */
   double *orig_sys_err_frac;    /* systematic error is (sys_err_frac * data) */
   double *orig_bgd;
   double *orig_flux;            /* full-resolution flux-corrected data */
   double *orig_flux_err;        /* full-resolution fcd uncertainty */
   double *flux_weights;         /* response-weights for full-res fcd */
   int orig_nbins;               /* original number of bins */
   int *rebin;                   /* index array for rebinning scheme */
   int *orig_notice;             /* ONLY for recording ignore/notice on unbinned data */
   int *quality;                 /* quality flags, for ignore_bad */

   int  n_notice;                /* [R] number of noticed bins */
   int *notice;                  /* [R] boolean */
   int *notice_list;             /* [R] index list of noticed bins */
   int *color;                   /* plot color - if not NULL, overrides color in Plot_t fmt */

   /* Header keywords: */

   char     object[CFLEN_VALUE];      /* target name */
   char    grating[CFLEN_VALUE];      /* HETG or LETG */
   char instrument[CFLEN_VALUE];      /* ACIS-S or HRC-S */
   char  *ancrfile;                   /* ARF file for this spectrum */
   char  *respfile;                   /* RMF file for this spectrum */
   double exposure;                   /* sum of good time intervals [sec] */
   double totcts;                     /* total counts */
   Area_Type area;                    /* area of extraction region [usually pixels^2] */
   Area_Type bgd_area;                /* area of bgd extraction region */
   double bgd_exposure;               /* bgd exposure time [sec] */

   double frame_time;                 /* CCD frame time (EXPTIME, sec) */
   double timedel, exptime;           /* TIMEDEL = EXPTIME + 0.04104 s */
   double tstart;

   Hist_Eval_Grid_Method_Type eval_grid_method;

   int spec_num;
   int order;                    /* TG_M  = spectral order */
   int part;                     /* TG_PART = code indicating heg, meg, leg, etc. */
   int srcid;                    /* TG_SRCID = code indicating which source in the fov */
   int has_grid;
   int swapped_bin_order;
   int is_fake_data;
};

/*}}}*/

static int update_notice_list (Hist_t *h);
static int init_rebin (Hist_t *h);
static void free_kernel (Hist_t *h);
static int initialize_rmf (Isis_Rmf_t *rmf, Hist_t *h);
static int attach_rsp (Hist_t *h, Isis_Rsp_t *rsp);
static int assign_matching_rsp (Hist_t *h, Isis_Rmf_t **rmf, Isis_Arf_t **arf);
static int set_model_grid (Hist_t *h);

/* new / free */

static void incr_rsp_refcount (Isis_Rsp_t *rsp) /*{{{*/
{
   if (rsp == NULL)
     return;

   rsp->arf->ref_count++;
   rsp->rmf->ref_count++;
}

/*}}}*/

static void release_arf (Isis_Arf_t *a) /*{{{*/
{
   if (a == NULL)
     return;

   a->ref_count--;
   if (Arf_is_identity (a))
     Arf_free_arf (a);
}

/*}}}*/

static void release_rmf (Isis_Rmf_t *r) /*{{{*/
{
   if (r == NULL)
     return;

   r->ref_count--;
   if (Rmf_is_identity (r))
     Rmf_free_rmf (r);
}

/*}}}*/

static void release_rsp (Isis_Rsp_t *rsp) /*{{{*/
{
   if (rsp == NULL)
     return;

   release_arf (rsp->arf);
   release_rmf (rsp->rmf);

   rsp->arf = NULL;
   rsp->rmf = NULL;
}

/*}}}*/

static void map_rsp_list (Isis_Rsp_t *rsp, void (*fcn)(Isis_Rsp_t *)) /*{{{*/
{
   Isis_Rsp_t *next;

   while (rsp != NULL)
     {
        next = rsp->next;
        (*fcn) (rsp);
        rsp = next;
     }
}

/*}}}*/

static void free_rsp (Isis_Rsp_t *rsp) /*{{{*/
{
   ISIS_FREE (rsp);
}

/*}}}*/

static void free_eval_grid_method (Hist_Eval_Grid_Method_Type *m) /*{{{*/
{
   if (m == NULL)
     return;
   if (m->destroy_options != NULL)
     (*m->destroy_options)(m->options);
}

/*}}}*/

static int init_eval_grid_method (Hist_Eval_Grid_Method_Type *m) /*{{{*/
{
   if (m == NULL)
     return -1;

   m->make_grid = NULL;
   m->eval_model = NULL;
   m->options = NULL;
   m->type = 0;
   m->id = 0;
   m->cache_model_values = 0;

   return 0;
}

/*}}}*/

static void free_hist (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return;

   ISIS_FREE (h->bin_lo);
   ISIS_FREE (h->bin_hi);
   ISIS_FREE (h->notice);
   ISIS_FREE (h->notice_list);
   ISIS_FREE (h->counts);
   ISIS_FREE (h->stat_err);
   ISIS_FREE (h->sys_err_frac);
   ISIS_FREE (h->flux);
   ISIS_FREE (h->flux_err);
   ISIS_FREE (h->orig_flux);
   ISIS_FREE (h->orig_flux_err);
   ISIS_FREE (h->flux_weights);
   ISIS_FREE (h->model_counts);
   ISIS_FREE (h->convolved_model_flux);
   ISIS_FREE (h->color);
   ISIS_FREE (h->rebin);
   ISIS_FREE (h->orig_bin_lo);
   ISIS_FREE (h->orig_bin_hi);
   ISIS_FREE (h->orig_counts);
   ISIS_FREE (h->orig_stat_err);
   ISIS_FREE (h->orig_sys_err_frac);
   ISIS_FREE (h->orig_bgd);
   ISIS_FREE (h->bgd_file);
   ISIS_FREE (h->ancrfile);
   ISIS_FREE (h->respfile);
   ISIS_FREE (h->file);
   ISIS_FREE (h->orig_notice);
   ISIS_FREE (h->quality);

   SLang_free_function (h->pre_combine);

   ISIS_FREE (h->instrumental_background_hook_name);
   SLang_free_function (h->instrumental_background_hook);

   SLang_free_function (h->assigned_model);
   isis_free_args (h->assigned_model_args);

   area_free (&h->area);
   area_free (&h->bgd_area);

   Isis_Hist_free (&h->model_flux);
   map_rsp_list (&h->a_rsp, &release_rsp);
   map_rsp_list (h->a_rsp.next, &free_rsp);
   if (h->f_rsp.arf && Arf_is_identity (h->f_rsp.arf))
     {
        /* FIXME? -
         * This is an ugly hack but it does prevent a memory leak
         * caused by assign_rsp (0,rmf,pha); when
         * using slang RMFs or user-defined RMFs (e.g. RMF_USER)
         * I guess this means my reference counting scheme needs work.
         */
        h->f_rsp.arf->ref_count--;
     }
   release_rsp (&h->f_rsp);
   free_kernel (h);

   free_eval_grid_method (&h->eval_grid_method);

   if (h->post_model_hook && h->post_model_hook_delete)
     (*h->post_model_hook_delete)(h->post_model_hook);

   if (h->stat_error_hook && h->stat_error_hook_delete)
     (*h->stat_error_hook_delete)(h->stat_error_hook);

   SLang_free_anytype (h->user_meta);

   ISIS_FREE (h);
}

/*}}}*/

static void area_init (Area_Type *a) /*{{{*/
{
   if (a == NULL)
     return;
   a->value.s = 1.0;
   a->is_vector = 0;
   a->n = 1;
}

/*}}}*/

static void area_free (Area_Type *a) /*{{{*/
{
   if (a == NULL)
     return;
   if (a->is_vector)
     ISIS_FREE(a->value.v);
   area_init (a);
}

/*}}}*/

static int area_set_vector_copy (Area_Type *at, double *area, unsigned int nbins) /*{{{*/
{
   double *area_copy;
   unsigned int i;

   if (NULL == (area_copy = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     return -1;

   for (i = 0; i < nbins; i++)
     {
        area_copy[i] = area[i];
     }

   area_free (at);
   at->value.v = area_copy;
   at->n = nbins;
   at->is_vector = 1;

   return 0;
}

/*}}}*/

static int area_copy (Area_Type *to, Area_Type *from) /*{{{*/
{
   if (to == NULL || from == NULL)
     return -1;

   if (from->is_vector)
     {
        Area_Type x;
        area_init(&x);
        if (-1 == area_set_vector_copy (&x, from->value.v, from->n))
          return -1;
        area_free (to);
        *to = x;
     }
   else
     {
        area_free (to);
        /* struct copy */
        *to = *from;
     }

   return 0;
}

/*}}}*/

static int area_set (Area_Type *a, double *area, int num) /*{{{*/
{
   if (a == NULL || num <= 0 || area == NULL)
     return -1;

   if (num == 1)
     {
        /* negative areas are right out. */
        if (*area < 0)
          return -1;
        /* Background area==0 is taken to mean that the background
         * is exact and that it should not be scaled.
         * Data area==0 is taken to mean that the background doesn't
         * matter (although this is kind of a silly case...).
         */
        area_free (a);
        a->is_vector = 0;
        a->value.s = *area;
     }
   else if (-1 == area_set_vector_copy (a, area, num))
     {
        return -1;
     }

   return 0;
}

/*}}}*/

static int area_get (Area_Type *at, double **area, int *nbins) /*{{{*/
{
   double *a;

   if (at == NULL)
     return -1;

   if (NULL == (a = (double *) ISIS_MALLOC (at->n * sizeof(double))))
     return -1;

   if (at->is_vector)
     {
        memcpy ((char *)a, (char *)at->value.v, at->n * sizeof(double));
     }
   else a[0] = at->value.s;

   *area = a;
   *nbins = at->n;

   return 0;
}

/*}}}*/

static int area_force_vector (Area_Type *a, int nbins) /*{{{*/
{
   double *av;
   double sa;
   int i;

   if (a->is_vector)
     return 0;

   if (NULL == (av = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     return -1;

   sa = a->value.s;

   for (i = 0; i < nbins; i++)
     {
        av[i] = sa;
     }

   a->value.v = av;
   a->is_vector = 1;
   a->n = nbins;

   return 0;
}

/*}}}*/

static Hist_t *Hist_new_hist (int nbins) /*{{{*/
{
   Hist_t *h = NULL;
   int i;

   if (NULL == (h = (Hist_t *) ISIS_MALLOC (sizeof(Hist_t))))
     return NULL;
   memset ((char *)h, 0, sizeof(*h));

   h->min_stat_err = Hist_Min_Stat_Err;

   area_init (&h->area);
   h->exposure = 1.0;

   area_init (&h->bgd_area);
   h->bgd_exposure = 1.0;

   h->timedel = -1.0;
   h->exptime = -1.0;
   h->frame_time = -1.0;
   h->tstart = -1.0;

   h->respfile = NULL;
   h->ancrfile = NULL;

   h->a_rsp.next = NULL;
   h->f_rsp.next = NULL;
   h->next = NULL;
   h->kernel = NULL;
   h->orig_bgd = NULL;

   h->flux_weights = NULL;

   h->pre_combine = NULL;

   h->post_model_hook = NULL;
   h->post_model_hook_delete = NULL;

   h->assigned_model = NULL;
   h->assigned_model_args = NULL;

   h->user_meta = NULL;

   h->stat_error_hook = NULL;
   h->stat_error_hook_delete = NULL;

   h->bgd_file = NULL;
   h->file = NULL;
   h->color = NULL;

   h->has_grid = 0;
   h->swapped_bin_order = 0;
   h->is_fake_data = 0;

   init_eval_grid_method (&h->eval_grid_method);

   /* 0 means 'not combined with any other dataset' */
   h->combo_id = 0;
   h->combo_weight = 1.0;

   /* allowed for null head node */
   if (nbins == 0)
     return h;

   if ((NULL == (h->bin_lo = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->bin_hi = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->notice = (int *) ISIS_MALLOC (nbins * sizeof(int))))
       || (NULL == (h->notice_list = (int *) ISIS_MALLOC (nbins * sizeof(int))))
       || (NULL == (h->counts = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->model_counts = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->convolved_model_flux = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->model_flux.val = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->stat_err = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->flux = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->flux_err = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->rebin = (int *) ISIS_MALLOC (nbins * sizeof(int))))
       || (NULL == (h->orig_bin_lo = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->orig_bin_hi = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->orig_counts = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->orig_stat_err = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->orig_notice = (int *) ISIS_MALLOC (nbins * sizeof(int))))
       || (NULL == (h->quality = (int *) ISIS_MALLOC (nbins * sizeof(int))))
       || (NULL == (h->orig_flux = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (h->orig_flux_err = (double *) ISIS_MALLOC (nbins * sizeof(double)))))
     {
        free_hist (h);
        return NULL;
     }

   memset ((char *)h->quality, 0, nbins*sizeof(int));
   memset ((char *)h->notice, 0, nbins*sizeof(int));
   memset ((char *)h->model_flux.val, 0, nbins*sizeof(double));
   memset ((char *)h->model_counts, 0, nbins*sizeof(double));
   memset ((char *)h->convolved_model_flux, 0, nbins*sizeof(double));
   memset ((char *)h->counts, 0, nbins*sizeof(double));
   memset ((char *)h->flux, 0, nbins*sizeof(double));
   memset ((char *)h->orig_flux, 0, nbins*sizeof(double));

   for (i = 0; i < nbins; i++)
     {
        h->flux_err[i] = -1.0;
        h->stat_err[i] = -1.0;
        h->orig_flux_err[i] = -1.0;
     }

   h->n_notice = nbins;
   h->nbins = nbins;
   h->orig_nbins = nbins;

   return h;
}

/*}}}*/

/* append, delete, find */

static int init_rebin (Hist_t *h) /*{{{*/
{
   int k, sign;
   size_t d_size;

   if (NULL == h)
     return -1;

   h->orig_nbins = h->nbins;

   sign = 1;
   for (k = 0; k < h->orig_nbins; k++)
     {
        sign *= -1;
        h->rebin[k] = sign;
        h->orig_notice[k] = 1;
     }

   if (h->flux_weights != NULL)
     {
        for (k = 0; k < h->orig_nbins; k++)
          {
             h->flux_weights[k] = 1.0;
          }
     }

   d_size = h->orig_nbins * sizeof(double);

   if (h->sys_err_frac)
     {
        memcpy ((char *)h->orig_sys_err_frac, (char *)h->sys_err_frac, d_size);
     }

   memcpy ((char *)h->orig_bin_lo, (char *)h->bin_lo, d_size);
   memcpy ((char *)h->orig_bin_hi, (char *)h->bin_hi, d_size);
   memcpy ((char *)h->orig_counts, (char *)h->counts, d_size);
   memcpy ((char *)h->orig_stat_err, (char *)h->stat_err, d_size);
   memcpy ((char *)h->orig_flux, (char *)h->flux, d_size);
   memcpy ((char *)h->orig_flux_err, (char *)h->flux_err, d_size);

   return 0;
}

/*}}}*/

static int finish_hist_init (Hist_t *h) /*{{{*/
{
   Isis_Rsp_t *a_rsp;
   int i;

   /* The call to init_rebin() will force h->nbins = h->orig_nbins */

   if (h->has_grid == 0)
     {
        for (i = 0; i < h->nbins; i++)
          {
             h->bin_lo[i] = i + 1.0;
             h->bin_hi[i] = h->bin_lo[i] + 1.0;
          }
        h->has_grid = -1;
     }

   if (-1 == init_rebin (h))
     return -1;

   h->totcts = 0.0;
   for (i=0; i < h->nbins; i++)
     {
        h->totcts += h->counts[i];
        h->notice[i] = 1;
     }

   if (-1 == update_notice_list (h))
     return -1;

   a_rsp = &h->a_rsp;

   if (a_rsp->arf == NULL)
     {
        Isis_Rmf_t *rmf = a_rsp->rmf;
        double *lo, *hi;
        unsigned int n;
        int malloced;

        /* If we have an RMF, the ARF should match it */
        if ((rmf == NULL)
            || (rmf->get_arf_grid == NULL)
            || (-1 == rmf->get_arf_grid (rmf, &lo, &hi, &n)))
          {
             malloced = 0;
             lo = h->bin_lo;
             hi = h->bin_hi;
             n = h->nbins;
          }
        else malloced = 1;

        if (NULL == (a_rsp->arf = Arf_make_identity_arf (lo, hi, n)))
          return -1;
        a_rsp->arf->ref_count++;

        if (malloced)
          {
             ISIS_FREE (lo);
             ISIS_FREE (hi);
          }
     }
   if (a_rsp->rmf == NULL)
     {
        a_rsp->rmf = matching_identity_rmf (h, a_rsp->arf);
        if (a_rsp->rmf == NULL)
          return -1;
        a_rsp->rmf->ref_count++;
     }

   /* struct copy */
   h->f_rsp = h->a_rsp;
   incr_rsp_refcount (&h->f_rsp);

   return set_model_grid (h);
}

/*}}}*/

static int __histogram_list_append (Hist_t *head, Hist_t *h) /*{{{*/
{
   Hist_t *prev;

   if (-1 == finish_hist_init (h))
     return -1;

   for (prev=head; prev->next != NULL; prev=prev->next)
     ;

   prev->next = h;
   h->next = NULL;
   h->index = prev->index + 1;

   return h->index;
}

/*}}}*/

static int histogram_list_append (Hist_t *head, Hist_t *h)
{
   int indx = __histogram_list_append (head, h);
   /* Changing the number of datasets may change the fit-function */
   update_user_model();
   return indx;
}

int Hist_map (Hist_t *head, int (*fun)(Hist_t *, void *), void *cl, int check_exclude) /*{{{*/
{
   Hist_t *h;

   if (head == NULL)
     return -1;

   Isis_Active_Dataset = 0;

   for (h = head->next; h != NULL; h = h->next)
     {
        if (check_exclude && h->exclude)
          continue;

        Isis_Active_Dataset = h->index;
        Hist_Current = h;

        if ((*fun)(h, cl))
          return -1;
     }

   return 0;
}

/*}}}*/

static int delete_after_hist (Hist_t *h) /*{{{*/
{
   Hist_t *dead;

   if (h == NULL)
     return -1;

   dead = h->next;
   h->next = dead->next;
   free_hist (dead);

   /* Changing the number of datasets may change the fit-function */
   update_user_model();

   return 0;
}

/*}}}*/

Hist_t *_Hist_find_hist_index (Hist_t * head, int hist_index) /*{{{*/
{
   Hist_t *h;

   if (head == NULL)
     return NULL;

   for (h=head->next; h != NULL; h = h->next)
     {
        if (h->index == hist_index)
          return h;
     }

   return NULL;
}

/*}}}*/

Hist_t *Hist_find_hist_index (Hist_t *head, int hist_index) /*{{{*/
{
   Hist_t *h = _Hist_find_hist_index (head, hist_index);

   if (h == NULL)
     isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "data set %d not found", hist_index);

   return h;
}

/*}}}*/

static int get_hist_version (Hist_t *h, unsigned int version, Isis_Hist_t *x) /*{{{*/
{
   if (h == NULL)
     return -1;

   x->bin_lo = h->bin_lo;
   x->bin_hi = h->bin_hi;
   x->nbins = h->nbins;
   x->notice = h->notice;
   x->notice_list = h->notice_list;
   x->n_notice = h->n_notice;
   x->sys_err_frac = h->sys_err_frac;

   if (is_data(version))
     {
        if (is_flux(version))
          {
             x->val = h->flux;
             x->val_err = h->flux_err;
          }
        else
          {
             x->val = h->counts;
             x->val_err = h->stat_err;
          }
     }
   else
     {
        x->val_err = NULL;

        if (is_convolved(version))
          x->val = h->convolved_model_flux;
        else if (is_flux(version))
          *x = h->model_flux;   /* struct copy */
        else
          x->val = h->model_counts;
     }

   return 0;
}

/*}}}*/

/* free, delete */

void Hist_free_list (Hist_t *head) /*{{{*/
{
   Hist_t *new_head;

   while (head != NULL)
     {
        new_head = head->next;
        free_hist (head);
        head = new_head;
     }
}

/*}}}*/

static int _delete_hist (Hist_t *head, int hist_index) /*{{{*/
{
   Hist_t *h;

   if (head == NULL)
     return -1;

   for (h=head; h->next != NULL; h = h->next)
     {
        if (h->next->index == hist_index)
          return delete_after_hist (h);
     }

   isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "data set %d not found", hist_index);
   return -1;
}

/*}}}*/

int Hist_delete_hist (Hist_t *head, int hist_index) /*{{{*/
{
   return _delete_hist (head, hist_index);
}

/*}}}*/

/* init, read, list */

Hist_t *Hist_init_list (void) /*{{{*/
{
   return Hist_new_hist (0);
}

/*}}}*/

static double counts_uncertainty (double counts, double min_stat_err) /*{{{*/
{
   if (counts >= 1.0)
     {
        return sqrt (counts);
     }
   else if (min_stat_err > 0.0)
     {
        return min_stat_err;
     }

   return 1.0;
}

/*}}}*/

static int fixup_stat_err (double *counts, double *stat_err, int n, double min_stat_err) /*{{{*/
{
   int i, reset = 0;

   if (min_stat_err > 0.0)
     {
        for (i = 0; i < n; i++)
          {
             if (stat_err[i] >= min_stat_err)
               continue;
             reset = 1;
             stat_err[i] = min_stat_err;
          }
     }
   else
     {
        for (i=0; i < n; i++)
          {
             if (stat_err[i] > 0.9)
               continue;
             reset = 1;
             stat_err[i] = counts_uncertainty (counts[i], min_stat_err);
          }
     }

   return reset;
}

/*}}}*/

static int validate_stat_err (Hist_t *h) /*{{{*/
{
   if (NULL == h)
     return -1;

   return fixup_stat_err (h->counts, h->stat_err, h->nbins, h->min_stat_err);
}
/*}}}*/

static int fixup_flux_err (double *flux, double *flux_err, int n) /*{{{*/
{
   int i, reset = 0;

   (void) flux;

   /* I'm assuming that INDEFS get set to DBL_MIN
    * on input and that any real flux_err value
    * would be larger than DBL_MIN.
    * Is there a cleaner approach?
    */

   for (i=0; i < n; i++)
     {
        if (flux_err[i] > DBL_MIN)
          continue;

        reset = 1;
        flux_err[i] = 1.0;
     }

   return reset;
}

/*}}}*/

static int validate_flux_err (Hist_t *h) /*{{{*/
{
   if (NULL == h)
     return -1;

   return fixup_flux_err (h->flux, h->flux_err, h->nbins);
}

/*}}}*/

static int get_canonical_hist_coordinates (Hist_t *h, int input_units) /*{{{*/
{
   Area_Type *a, *b;
   unsigned int h_nbins;
   int status;

   status = get_canonical_coordinates (h->bin_lo, h->bin_hi, &h->nbins, input_units);

   /* in case nbins got truncated */
   if (h->n_notice != h->nbins)
     {
        int i;
        h->n_notice = 0;
        ISIS_FREE (h->notice_list);
        h->notice_list = (int *) ISIS_MALLOC (h->nbins * sizeof(int));
        if (h->notice_list == NULL)
          return -1;
        h->n_notice = h->nbins;
        for (i = 0; i < h->nbins; i++)
          {
             h->notice_list[i] = i;
             h->notice[i] = 1;
          }
     }

   switch (status)
     {
      case BINS_REVERSED:
        h_nbins = h->nbins;
        /* DONT reverse notice_list */
        if (-1 == reverse_d (h->counts, h_nbins)
            || -1 == reverse_d (h->stat_err, h_nbins)
            || -1 == reverse_d (h->flux, h_nbins)
            || -1 == reverse_d (h->flux_err, h_nbins)
            || -1 == reverse_i (h->notice, h_nbins)
            || -1 == reverse_i (h->quality, h_nbins))
          return -1;
        a = &h->area;
        if (a->is_vector)
          {
             if (-1 == reverse_d (a->value.v, h_nbins))
               return -1;
          }
        b = &h->bgd_area;
        if (b->is_vector)
          {
             if (-1 == reverse_d (b->value.v, h_nbins))
               return -1;
          }
        if (h->orig_bgd)
          {
             if (-1 == reverse_d (h->orig_bgd, h_nbins))
               return -1;
          }
        if (h->sys_err_frac)
          {
             if (-1 == reverse_d (h->sys_err_frac, h_nbins))
               return -1;
          }
        /* drop */
      case BINS_OK:
        break;

      case BINS_INVALID:
      default:
        return -1;
     }

   h->has_grid = 1;

   return 0;
}

/*}}}*/

static void invalid_uncertainties_replaced (void) /*{{{*/
{
   if (Hist_Warn_Invalid_Uncertainties)
     isis_vmesg (WARN, I_INVALID_UNCERT, __FILE__, __LINE__, NULL);
}

/*}}}*/

static int do_instrument_specific_hacks (Hist_t *h) /*{{{*/
{
   if ((0 == strncmp (h->instrument, "ACIS", 4))
       || (0 == strncmp (h->instrument, "acis", 4)))
     {
        if (h->exptime > 0.0)
          h->frame_time = h->exptime;
        else if (h->timedel > 0.0)
          h->frame_time = h->timedel - 0.04104;
        else /* use standard EXPTIME for ACIS timed mode */
          h->frame_time = 3.2;

        if (h->frame_time <= 0.0)
          {
             /* time to clock out 3 pixels */
             h->frame_time = 0.00285*3;
             isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "Assuming frame time = %g", h->frame_time);
             isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "If this is incorrect, use set_frame_time() to provide a value");
          }
     }

   return 0;
}

/*}}}*/

/*{{{ keyword table */

/* When updating this table,
 * remember to update copy_histogram_keywords ()
 */

#define OFFSETOF(x)  (offsetof(struct _Hist_t,x))

static Keyword_t Hist_Keyword_Table[] =
{
     { "OBJECT",     OPTIONAL_KEY,   STRING_KEY,   OFFSETOF(object)     },
     { "INSTRUME",   OPTIONAL_KEY,   STRING_KEY,   OFFSETOF(instrument) },
     { "GRATING",    OPTIONAL_KEY,   STRING_KEY,   OFFSETOF(grating)    },
     { "EXPOSURE",   OPTIONAL_KEY,   DOUBLE_KEY,   OFFSETOF(exposure)   },
     { "EXPTIME",    OPTIONAL_KEY,   DOUBLE_KEY,   OFFSETOF(exptime)    },
     { "TSTART",     OPTIONAL_KEY,   DOUBLE_KEY,   OFFSETOF(tstart)     },
     { "TIMEDEL",    OPTIONAL_KEY,   DOUBLE_KEY,   OFFSETOF(timedel)    },
     { "RFLORDER",   OPTIONAL_KEY,   INT_KEY,      OFFSETOF(order)      },
     { "TG_M",       OPTIONAL_KEY,   INT_KEY,      OFFSETOF(order)      },
     { "TG_PART",    OPTIONAL_KEY,   INT_KEY,      OFFSETOF(part)       },
     { "TG_SRCID",   OPTIONAL_KEY,   INT_KEY,      OFFSETOF(srcid)      },
     { NULL,         OPTIONAL_KEY,   NULL_KEY,     0 }
};

/*}}}*/

static Hist_t *read_ascii (char *filename, double min_stat_err) /*{{{*/
{
   FILE *fp = NULL;
   Hist_t *h = NULL;
   Keyword_t *keytable = Hist_Keyword_Table;
   unsigned int bin_type = 0;
   double *bin_val, *bin_err;
   double backscal;
   int (*validate_fcn)(Hist_t *) = NULL;   /* ptr to validate fcn */
   int i, nbins, valid;
   int input_units = U_ANGSTROM;
   unsigned int num_ascii_keys;

   /* data type defaults to counts */
   bit_clear(bin_type, H_FLUX);

   if (filename == NULL)
     return NULL;

   /* Scan the file to count the number of valid data lines
    * and the number of keywords.  (yuk)
    */

   if ((0 == is_regular_file (filename)) /* paranoia */
       || (NULL == (fp = fopen (filename, "r"))))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", filename);
        return NULL;
     }

   num_ascii_keys = 0;
   nbins = 0;

   while (!feof (fp))
     {
        double lo, hi, val, unc;
        char buf[BUFSIZE];

        if (NULL == fgets (buf, BUFSIZE, fp))
          break;

        if (!feof (fp) && (NULL == strchr (buf, '\n')))
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s, lines too long", filename);
             return NULL;
          }

        if (buf[0] == COMMENT_CHAR)
          continue;
        else if (buf[0] == KEYWORD_CHAR)
          num_ascii_keys++;
        else if (4 == sscanf (buf, "%le %le %le %le", &lo, &hi, &val, &unc))
          nbins++;
        else if (strlen (buf) != strspn (buf, " \t\n"))
          {
             nbins = 0;
             break;
          }
     }

   fclose(fp);

   if (nbins == 0)
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "bad format, %s", filename);
        return NULL;
     }

   /* Now, allocate space and read the file */

   if (NULL == (h = Hist_new_hist (nbins)))
     return NULL;

   h->min_stat_err = min_stat_err;

   if (NULL == (fp = fopen (filename, "r")))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", filename);
        free_hist (h);
        return NULL;
     }

   if (num_ascii_keys == 0)
     {
        isis_vmesg (WARN, I_FILE_NOKEYS, __FILE__, __LINE__, "%s", filename);
        isis_strcpy (h->object, filename, CFLEN_VALUE);
     }
   else
     {
        char value[CFLEN_VALUE];
        Ascii_Keybuf_t * keybuf = NULL;

        if ((NULL == (keybuf = Key_ascii_allocate_keybuf (num_ascii_keys)))
            || (-1 == Key_ascii_load_keybuf (fp, keybuf)))
          {
             isis_vmesg (WARN, I_READ_KEY_FAILED, __FILE__, __LINE__, NULL);
             fclose (fp);
             Key_ascii_free_keybuf (keybuf);
             free_hist (h);
             return NULL;
          }

        (void) Key_read_header_keywords (keybuf, (char *)h, keytable, ASCII_FILE);

        if (-1 == Key_ascii_read_string_keyword (value, "xunit", keybuf))
          input_units = U_ANGSTROM;
        else if (-1 == (input_units = unit_id (value)))
          {
             fclose (fp);
             Key_ascii_free_keybuf (keybuf);
             free_hist (h);
             return NULL;
          }

        if (0 == Key_ascii_read_string_keyword (value, "bintype", keybuf)
            && 0 == isis_strcasecmp (value, "FLUX"))
          bit_set(bin_type,H_FLUX);

        if (0 == Key_ascii_read_double_keyword (&backscal, "backscal", keybuf))
          {
             if (-1 == area_set (&h->area, &backscal, 1))
               {
                  fclose (fp);
                  Key_ascii_free_keybuf (keybuf);
                  free_hist (h);
                  return NULL;
               }
          }

        Key_ascii_free_keybuf (keybuf);
     }

   if (is_flux(bin_type))
     {
        memset ((char *)h->counts, 0, nbins * sizeof(double));
        for (i = 0; i < nbins; i++)
          h->stat_err[i] = 1.0;
        bin_val = h->flux;
        bin_err = h->flux_err;
        validate_fcn = validate_flux_err;
     }
   else
     {
        memset ((char *)h->flux, 0, nbins * sizeof(double));
        for (i = 0; i < nbins; i++)
          h->flux_err[i] = 1.0;
        bin_val = h->counts;
        bin_err = h->stat_err;
        validate_fcn = validate_stat_err;
     }

   /* I'm not using a for() loop because we might need to skip
    * lines without incrementing i (comments or format errors)
    */

   i=0;

   while (!feof(fp) && i < nbins)
     {
        char buf[BUFSIZE];

        if (NULL == fgets (buf, BUFSIZE, fp))
          break;

        if (buf[0] == COMMENT_CHAR)
          continue;

        if (4 != sscanf(buf, "%le %le %le %le", &h->bin_lo[i], &h->bin_hi[i],
                        &bin_val[i], &bin_err[i]))
          {
             isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "skipped: %s", buf);
             continue;
          }

        i++;
     }

   fclose(fp);

   if (i != nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "truncated file");
        if (i==0)
          {
             free_hist (h);
             return NULL;
          }
        h->nbins = nbins = i;
     }

   if (-1 == get_canonical_hist_coordinates (h, input_units))
     {
        free_hist (h);
        return NULL;
     }

   valid = (*validate_fcn) (h);

   if (-1 == valid)
     {
        free_hist (h);
        return NULL;
     }
   else if (1 == valid)
     invalid_uncertainties_replaced ();

   (void) do_instrument_specific_hacks (h);

   isis_vmesg (INFO, I_READ_OK, __FILE__, __LINE__, "%d bins from %s",
               nbins, filename);
   return h;
}

/*}}}*/

int Hist_read_ascii (Hist_t *head, char *filename, double min_stat_err) /*{{{*/
{
   Hist_t *h;

   if (NULL == (h = read_ascii (filename, min_stat_err)))
     return -1;

   if (NULL == (h->file = isis_make_string (filename)))
     {
        free_hist (h);
        return -1;
     }

   return histogram_list_append (head, h);
}

/*}}}*/

static Hist_t *create_hist_from_grid (Hist_t *head, Isis_Hist_t *x, unsigned int x_unit) /*{{{*/
{
   Hist_t *h = NULL;
   unsigned int i, n, size;

   if (x == NULL)
     return NULL;

   n = x->nbins;

   if (NULL == (h = Hist_new_hist (n)))
     return NULL;

   h->is_fake_data = 1;

   size = n * sizeof(double);

   memcpy ((char *)h->bin_lo, (char *)x->bin_lo, size);
   memcpy ((char *)h->bin_hi, (char *)x->bin_hi, size);
   if (-1 == get_canonical_hist_coordinates (h, x_unit))
     {
        free_hist (h);
        return NULL;
     }

   memset ((char *)h->counts, 0, size);
   memset ((char *)h->flux, 0, size);

   for (i = 0; i < n; i++)
     {
        h->stat_err[i] = 1.0;
        h->flux_err[i] = 1.0;
     }

   if (-1 == histogram_list_append (head, h))
     {
        free_hist (h);
        return NULL;
     }

   return h;
}

/*}}}*/

static int update_hist_sort_order (Hist_t *head, Hist_t *h) /*{{{*/
{
   Hist_t *q, *h_prev;

   for (q=head; q->next != NULL; q=q->next)
     {
        if (q->next->index > h->index)
          {
             break;
          }
        else if ((q->next->index == h->index) && (q->next != h))
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "%d already exists", h->index);
             return -1;
          }
     }

   if (q->next == NULL)
     return 0;

   for (h_prev=q; h_prev->next != NULL; h_prev=h_prev->next)
     {
        if (h_prev->next == h)
          break;
     }

   h_prev->next = h->next;
   h->next = q->next;
   q->next = h;

   return 0;
}

/*}}}*/

static Hist_t *create_hist_from_arf (Hist_t *head, int hist_index, Isis_Arf_t *a) /*{{{*/
{
   Hist_t *h;
   Isis_Hist_t x;

   if (hist_index <= 0)
     return NULL;

   x.bin_lo = a->bin_lo;
   x.bin_hi = a->bin_hi;
   x.nbins = a->nbins;

   if (NULL == (h = create_hist_from_grid (head, &x, U_ANGSTROM)))
     return NULL;

   h->index = hist_index;
   if (-1 == update_hist_sort_order (head, h))
     {
        free_hist (h);
        return NULL;
     }

   return h;
}

/*}}}*/

static Hist_t *create_hist_from_rmf (Hist_t *head, int hist_index, Isis_Rmf_t *r) /*{{{*/
{
   Hist_t *h;
   Isis_Hist_t x;
   int energy_ordered_ebounds;
   unsigned int n;

   if (hist_index <= 0)
     return NULL;

   if (-1 == r->get_data_grid (r, &x.bin_lo, &x.bin_hi, &n, &energy_ordered_ebounds))
     return NULL;

   x.nbins = n;
   if (NULL == (h = create_hist_from_grid (head, &x, U_ANGSTROM)))
     return NULL;

   ISIS_FREE (x.bin_lo);
   ISIS_FREE (x.bin_hi);

   h->order = r->order;

   h->index = hist_index;
   if (-1 == update_hist_sort_order (head, h))
     {
        free_hist (h);
        return NULL;
     }

   return h;
}

/*}}}*/

static int redefine_fake_data_grid (Hist_t *head, Isis_Rmf_t *r, Hist_t **h) /*{{{*/
{
   int indx = (*h)->index;
   Isis_Rsp_t rsp;
   Isis_Arf_t *arf;
   SLang_Name_Type *instrumental_background_hook;
   char *instrumental_background_hook_name;

   /* preserve some of the existing configuration */
   arf = (*h)->a_rsp.arf;
   if (Arf_is_identity (arf))
     arf = NULL;
   instrumental_background_hook = (*h)->instrumental_background_hook;

   instrumental_background_hook_name = isis_make_string ((*h)->instrumental_background_hook_name); /* NULL ok */

   if (-1 == _delete_hist (head, indx))
     {
        ISIS_FREE (instrumental_background_hook_name);
        return -1;
     }

   *h = create_hist_from_rmf (head, indx, r);
   if (*h == NULL)
     {
        ISIS_FREE (instrumental_background_hook_name);
        return -1;
     }

   rsp.next = NULL;
   rsp.rmf = r;
   /* restore pre-existing configuration */
   if (arf != NULL)
     rsp.arf = arf;
   else
     rsp.arf = (*h)->a_rsp.arf;
   (*h)->instrumental_background_hook_name = instrumental_background_hook_name;
   (*h)->instrumental_background_hook = instrumental_background_hook;

   return attach_rsp (*h, &rsp);
}

/*}}}*/

int Hist_define_data (Hist_t *head, Isis_Hist_t *x, unsigned int bin_type, unsigned int has_grid, /*{{{*/
                      double min_stat_err)
{
   int (*validate_fcn)(Hist_t *) = NULL;
   Hist_t *h = NULL;
   unsigned int size;
   int valid;
   int id = -1;

   if (x == NULL || head == NULL)
     return -1;

   if (NULL == (h = Hist_new_hist (x->nbins)))
     return -1;

   size = x->nbins * sizeof(double);

   if (is_flux(bin_type))
     {
        validate_fcn = validate_flux_err;
        memcpy ((char *)h->flux, (char *)x->val, size);
        if (x->val_err)
          memcpy ((char *)h->flux_err, (char *)x->val_err, size);
        else
          memset ((char *)h->flux_err, 0, size);
     }
   else
     {
        h->min_stat_err = min_stat_err;
        validate_fcn = validate_stat_err;
        memcpy ((char *)h->counts, (char *)x->val, size);
        if (x->val_err)
          memcpy ((char *)h->stat_err, (char *)x->val_err, size);
        else
          memset ((char *)h->stat_err, 0, size);
     }

   if (has_grid)
     {
        memcpy ((char *)h->bin_lo, (char *)x->bin_lo, size);
        memcpy ((char *)h->bin_hi, (char *)x->bin_hi, size);
        if (-1 == get_canonical_hist_coordinates (h, U_ANGSTROM))
          {
             isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "bad grid");
             goto finish;
          }
     }

   valid = (*validate_fcn) (h);
   if (-1 == valid)
     goto finish;
   else if (1 == valid)
     invalid_uncertainties_replaced ();

   id = histogram_list_append (head, h);

   finish:
   if (id < 0) ISIS_FREE (h);

   return id;
}

/*}}}*/

static void delete_bgd (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return;

   h->bgd_exposure = 1.0;
   area_free (&h->bgd_area);
   ISIS_FREE (h->orig_bgd);
}

/*}}}*/

static double *background_scale_factor (Hist_t *h, int do_rebin) /*{{{*/
{
   Area_Type da, ba;
   double *da_v, *ba_v, *s;
   double src_exposure, bgd_exposure, f;
   unsigned int i, n;

   if (h == NULL)
     return NULL;

   if (-1 == get_exposure_time (h, &src_exposure))
     return NULL;

   /* Background exposure may not be up to date */

   if (h->bgd_exposure <= 0.0)
     h->bgd_exposure = src_exposure;
   bgd_exposure = h->bgd_exposure;

   f = src_exposure / bgd_exposure;

   n = do_rebin ? h->nbins : h->orig_nbins;

   if (NULL == (s = (double *) ISIS_MALLOC (n * sizeof(double))))
     return NULL;

   /* quick return for exact background */
   if ((h->bgd_area.is_vector == 0) && (h->bgd_area.value.s == 0.0))
     {
        /* background is exact -- don't scale it */
        for (i = 0; i < n; i++)
          {
             s[i] = 1.0;
          }
        return s;
     }

   area_init (&da);
   area_init (&ba);
   area_copy (&da, &h->area);
   area_copy (&ba, &h->bgd_area);

   /* quick return if both areas are scalars */
   if ((da.is_vector == 0) && (ba.is_vector == 0))
     {
        double r = da.value.s / ba.value.s;
        for (i = 0; i < n; i++)
          {
             s[i] = f * r;
          }

        area_free (&da);
        area_free (&ba);
        return s;
     }

   /* quick return if rebinning is not needed */
   if (do_rebin == 0)
     {
        if (ba.is_vector && da.is_vector)
          {
             da_v = da.value.v;
             ba_v = ba.value.v;
             for (i = 0; i < n; i++)
               {
                  if (ba_v[i] > 0.0)
                    s[i] = f * da_v[i] / ba_v[i];
                  else s[i] = 0.0;
               }
          }
        else if (ba.is_vector)
          {
             f *= da.value.s;
             ba_v = ba.value.v;
             for (i = 0; i < n; i++)
               {
                  if (ba_v[i] > 0.0)
                    s[i] = f / ba_v[i];
                  else s[i] = 0.0;
               }
          }
        else /* da.is_vector */
          {
             f /= ba.value.s;
             da_v = da.value.v;
             for (i = 0; i < n; i++)
               {
                  s[i] = f * da_v[i];
               }
          }

        area_free (&da);
        area_free (&ba);
        return s;
     }

   /* ok, we're rebinning -- simplify processing by forcing the
    * areas to be vectors
    */

   if ((-1 == rebin_backscale (&da, h, h->bin_lo, h->bin_hi, h->nbins))
       || (-1 == rebin_backscale (&ba, h, h->bin_lo, h->bin_hi, h->nbins)))
     {
        area_free (&da);
        area_free (&ba);
        ISIS_FREE(s);
        return NULL;
     }

   da_v = da.value.v;
   ba_v = ba.value.v;

   for (i = 0; i < n; i++)
     {
        if (ba_v[i] > 0)
          s[i] = f * da_v[i] / ba_v[i];
        else s[i] = 0.0;
     }

   area_free (&da);
   area_free (&ba);
   return s;
}

/*}}}*/

static int make_scaling_vector (Hist_t *h, int do_rebin, int pack_noticed, /*{{{*/
                                Area_Type *a, double t, double **pvat, int *num_vat)
{
   double *vat = NULL, *av = NULL;
   int i, k, n;

   n = do_rebin ? h->nbins : h->orig_nbins;
   *num_vat = n;

   if (NULL == (vat = (double *) ISIS_MALLOC (n * sizeof(double))))
     return -1;
   *pvat = vat;

#if 0
   (void) a;
   for (i = 0; i < n; i++)
     {
        vat[i] = t;
     }
#else
   if (do_rebin == 0)
     {
        double at;
        if (a->is_vector == 0)
          {
             at = t;
             if (a->value.s > 0) at *= a->value.s;

             for (i = 0; i < n; i++)
               {
                  vat[i] = at;
               }
          }
        else
          {
             av = a->value.v;
             for (i = 0; i < n; i++)
               {
                  vat[i] = t * av[i];
               }
          }
     }
   else
     {
        Area_Type a_cpy;

        area_init (&a_cpy);
        area_copy (&a_cpy, a);
        if (-1 == rebin_backscale (&a_cpy, h, h->bin_lo, h->bin_hi, h->nbins))
          {
             area_free (&a_cpy);
             ISIS_FREE(vat);
             return -1;
          }

        av = a_cpy.value.v;
        for (i = 0; i < n; i++)
          {
             vat[i] = t * av[i];
          }
        area_free (&a_cpy);
     }
#endif

   if (pack_noticed == 0)
     return 0;

   if (NULL == (av = (double *)ISIS_MALLOC (h->n_notice * sizeof(double))))
     return -1;
   for (k = 0; k < h->n_notice; k++)
     {
        i = h->notice_list[k];
        av[k] = vat[i];
     }

   *pvat = av;
   *num_vat = h->n_notice;

   ISIS_FREE(vat);

   return 0;
}

/*}}}*/

int Hist_scaling_vectors (Hist_t *h, int do_rebin, int pack_noticed, /*{{{*/
                          double **psrc_at, double **pbkg_at, int *pnum)
{
   double src_exposure, bgd_exposure;
   double *src_at=NULL, *bkg_at=NULL;
   int n;

   if (-1 == get_exposure_time (h, &src_exposure))
     return -1;

   /* Background exposure may not be up to date */

   if (h->bgd_exposure <= 0.0)
     h->bgd_exposure = src_exposure;
   bgd_exposure = h->bgd_exposure;

   if ((-1 == make_scaling_vector (h, do_rebin, pack_noticed, &h->bgd_area, bgd_exposure, &bkg_at, &n))
       || (-1 == make_scaling_vector (h, do_rebin, pack_noticed, &h->area, src_exposure, &src_at, &n)))
     {
        ISIS_FREE(bkg_at);
        ISIS_FREE(src_at);
        return -1;
     }

   *pnum = n;
   *pbkg_at = bkg_at;
   *psrc_at = src_at;

   return 0;
}

/*}}}*/

static int copy_input_background (Hist_t *h, int do_rebin, double **bc, int *nbc) /*{{{*/
{
   double *b = NULL;
   int n;

   if (h == NULL || h->orig_bgd == NULL
       || bc == NULL || nbc == NULL)
     return -1;

   *bc = NULL;
   *nbc = 0;

   n = do_rebin ? h->nbins : h->orig_nbins;

   if (NULL == (b = (double *) ISIS_MALLOC (n * sizeof(double))))
     return -1;

   if (do_rebin)
     {
        if (-1 == apply_rebin (h->orig_bgd, h->orig_nbins, h->rebin, h->nbins, b))
          {
             ISIS_FREE(b);
             return -1;
          }
     }
   else
     {
        memcpy ((char *)b, (char *)h->orig_bgd, n * sizeof(double));
     }

   *bc = b;
   *nbc = n;

   return 0;
}

/*}}}*/

static int scale_background (Hist_t *h, int do_rebin, /*{{{*/
                             double **scaled_bgd, double **scaled_bgd_err)
{
   double *scale=NULL, *b=NULL, *berr=NULL;
   int i, n, status = -1;

   if (h == NULL)
     return -1;

   if (h->orig_bgd == NULL)
     return 0;

   if (-1 == copy_input_background (h, do_rebin, &b, &n))
     return -1;

   if (NULL == (scale = background_scale_factor (h, do_rebin)))
     goto error_return;

   if (NULL == (berr = (double *) ISIS_MALLOC (n * sizeof(double))))
     goto error_return;

   for (i = 0; i < n; i++)
     {
        /* note - order is important here */
        berr[i] = scale[i] * counts_uncertainty (b[i], h->min_stat_err);
        b[i] *= scale[i];
     }

   status = 0;
   error_return:

   ISIS_FREE(scale);

   if (scaled_bgd)
     {
        *scaled_bgd = b;
     }
   else ISIS_FREE(b);

   if (scaled_bgd_err)
     {
        *scaled_bgd_err = berr;
     }
   else ISIS_FREE(berr);

   return status;
}

/*}}}*/

int Hist_define_background (Hist_t *h, double bgd_exposure, /*{{{*/
                            double *bgd_area, int area_is_vector,
                            double *bgd, unsigned int nbins)
{
   Area_Type at;
   double *bgd_copy = NULL;

   if (h == NULL)
     return -1;

   /* bgd==NULL is ok */
   if (bgd == NULL)
     {
        delete_bgd (h);
        return 0;
     }

   if (h->nbins != h->orig_nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Can't define background for binned data");
        return -1;
     }

   if (nbins != (unsigned int) h->orig_nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "size mismatch: data has %d bins",
                    h->orig_nbins);
        return -1;
     }

   area_init (&at);

   if (area_is_vector)
     {
        if (-1 == area_set_vector_copy (&at, bgd_area, nbins))
          return -1;
     }
   else if (*bgd_area > 0)
     {
        at.is_vector= 0;
        at.value.s = *bgd_area;
     }
   else if (-1 == area_copy (&at, &h->area))
     {
        return -1;
     }

   if (NULL == (bgd_copy = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     {
        area_free (&at);
        return -1;
     }

   /* I'm assuming the bins are in the right order.. */
   memcpy ((char *)bgd_copy, (char *)bgd, nbins * sizeof(double));

   delete_bgd (h);
   h->orig_bgd = bgd_copy;

   /* struct copy */
   h->bgd_area = at;

   if (bgd_exposure > 0.0)
     h->bgd_exposure = bgd_exposure;
   else
     {
        /* ARF may not be assigned yet */
        h->bgd_exposure = -1.0;
     }

   return 0;
}

/*}}}*/

static int alloc_sys_err (Hist_t *h, int num) /*{{{*/
{
   ISIS_FREE(h->sys_err_frac);
   ISIS_FREE(h->orig_sys_err_frac);

   if ((NULL == (h->sys_err_frac = (double *) ISIS_MALLOC (num * sizeof(double))))
       || (NULL == (h->orig_sys_err_frac = (double *) ISIS_MALLOC (num * sizeof(double)))))
     return -1;

   return 0;
}

/*}}}*/

static int assign_sys_err (Hist_t *h, double *sys_err, int num) /*{{{*/
{
   int i;

   if (num == 1)
     {
        for (i = 0; i < h->nbins; i++)
          {
             h->sys_err_frac[i] = *sys_err;
          }
     }
   else if (num == h->nbins)
     {
        for (i = 0; i < h->nbins; i++)
          {
             h->sys_err_frac[i] = sys_err[i];
          }
     }
   else
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "sys_err_frac array size mismatch, len = %d", num);
        return -1;
     }

   return 0;
}

/*}}}*/

int Hist_define_sys_err_frac (Hist_t *h, double *sys_err_frac, int nbins) /*{{{*/
{
   if (h == NULL)
     return -1;

   if (sys_err_frac == NULL || nbins == 0)
     {
        ISIS_FREE(h->sys_err_frac);
        ISIS_FREE(h->orig_sys_err_frac);
        return 0;
     }

   if (h->nbins != h->orig_nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Can't define systematic error for binned data");
        return -1;
     }

   if ((nbins != h->orig_nbins) && (nbins != 1))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "size mismatch: data has %d bins",
                    h->orig_nbins);
        return -1;
     }

   if (h->sys_err_frac == NULL)
     {
        if (-1 == alloc_sys_err (h, h->orig_nbins))
          return -1;
     }

   /* nbins = 1 or h->orig_nbins */
   if (-1 == assign_sys_err (h, sys_err_frac, nbins))
     return -1;

   memcpy ((char *)h->orig_sys_err_frac, (char *)h->sys_err_frac,
           h->orig_nbins * sizeof(double));

   return 0;
}

/*}}}*/

int Hist_copy_sys_err_frac (Hist_t *h, double **sys_err_frac, int *nbins) /*{{{*/
{
   if (h == NULL)
     return -1;

   if (h->sys_err_frac == NULL)
     {
        *sys_err_frac = NULL;
        *nbins = 0;
        return 0;
     }

   if (-1 == copy_d_array (sys_err_frac, h->sys_err_frac, h->nbins))
     return -1;
   *nbins = h->nbins;

   return 0;
}

/*}}}*/

int Hist_copy_histogram_keywords (Hist_t *dst, Hist_t *src) /*{{{*/
{
   if (NULL == dst || NULL == src)
     return -1;

   isis_strcpy (dst->object, src->object, CFLEN_VALUE);
   isis_strcpy (dst->instrument, src->instrument, CFLEN_VALUE);
   isis_strcpy (dst->grating, src->grating, CFLEN_VALUE);
   dst->exposure = src->exposure;
   dst->exptime = src->exptime;
   dst->tstart = src->tstart;
   dst->timedel = src->timedel;
   dst->frame_time = src->frame_time;
   dst->order = src->order;
   dst->part = src->part;
   dst->srcid = src->srcid;

   if (-1 == area_copy (&dst->area, &src->area))
     return -1;

   return 0;
}

/*}}}*/

int Hist_run_rmf_post_fit_method (Hist_t *h)
{
   return Rmf_run_post_fit_method (h->a_rsp.rmf);
}

static int set_hist_grid_using_rmf (Isis_Rmf_t *rmf, Hist_t *h) /*{{{*/
{
   double *lo = NULL;
   double *hi = NULL;
   unsigned int i, n;
   int energy_ordered_ebounds;

   if ((NULL == rmf) || (NULL == h))
     return -1;

   /* Its too much trouble to first define the 'orig' grid
    * with the RMF and then reconstruct the binned grid.
    * Much undo the grouping, assign the RMF and restore
    * the grouping (all at a higher level)
    */
   if (h->nbins != h->orig_nbins)
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "operation not supported for rebinned data");
        return -1;
     }

   if (-1 == rmf->get_data_grid (rmf, &lo, &hi, &n, &energy_ordered_ebounds))
     return -1;

   if (n != (unsigned int) h->orig_nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "grid mismatch, data=%d RMF=%d channels",
                    h->orig_nbins, n);
        ISIS_FREE (lo);
        ISIS_FREE (hi);
        return -1;
     }

   for (i = 0; i < n; i++)
     {
        h->bin_lo[i] = lo[i];
        h->bin_hi[i] = hi[i];
     }

   ISIS_FREE (lo);
   ISIS_FREE (hi);

   /* I'm assuming the data arrays and the ebounds grid
    * are in the same order.  This might break.. FIXME?
    */

   if (energy_ordered_ebounds
       && (h->has_grid <= 0)
       && (h->swapped_bin_order == 0))
     {
        h->swapped_bin_order = 1;
        reverse_d (h->counts, n);
        reverse_d (h->stat_err, n);
        reverse_d (h->flux, n);
        reverse_d (h->flux_err, n);
        reverse_i (h->quality, n);
        if (h->sys_err_frac) reverse_d (h->sys_err_frac, n);
     }

   h->order = rmf->order;

   if (-1 == get_canonical_hist_coordinates (h, U_ANGSTROM))
     return -1;

   return init_rebin (h);
}

/*}}}*/

static int get_grid_from_rmf (Hist_t *h, Isis_Rmf_t *rmf_head) /*{{{*/
{
   int rmf_index = Rmf_load_rmf (rmf_head, RMF_FILE, h->respfile);
   if (-1 == rmf_index)
     return -1;

   h->a_rsp.rmf = Rmf_find_rmf_index (rmf_head, rmf_index);
   if ((h->a_rsp.rmf == NULL)
       || (-1 == initialize_rmf (h->a_rsp.rmf, h)))
     return -1;

   return set_hist_grid_using_rmf (h->a_rsp.rmf, h);
}

/*}}}*/

static int read_fits_backscal_keywords (cfitsfile *fp, Hist_t *h) /*{{{*/
{
   double backscal;

   if (0 == cfits_read_double_keyword (&backscal, "BACKSCAL", fp))
     {
        if (-1 == area_set (&h->area, &backscal, 1))
          return -1;
     }

   if (0 == cfits_read_double_keyword (&backscal, "BG_AREA", fp))
     {
        if (-1 == area_set (&h->bgd_area, &backscal, 1))
          return -1;
     }

   return 0;
}

/*}}}*/

static int read_fits_backscal_column (cfitsfile *fp, int k, Hist_t *h) /*{{{*/
{
   double *b = NULL;
   int num;

   if (-1 == cfits_get_repeat_count (&num, "BACKSCAL", fp))
     return -1;

   if ((num != 1) && (num != h->nbins))
     return -1;

   if (NULL == (b = (double *) ISIS_MALLOC (h->nbins * sizeof(double))))
     return -1;

   if (0 == cfits_read_optional_double_col (b, h->nbins, k, "BACKSCAL", fp))
     {
        if (-1 == area_set (&h->area, b, h->nbins))
          {
             ISIS_FREE(b);
             return -1;
          }
     }

   ISIS_FREE(b);
   return 0;
}

/*}}}*/

static int apply_areascale (cfitsfile *fp, int k, Hist_t *h, double areascal_keyword_value, int have_areascal_keyword) /*{{{*/
{
   double *a = NULL;
   int i;

   /* XMM RGS uses AREASCAL.  Why? */

   if (have_areascal_keyword == 0)
     {
        int num;

        if (-1 == cfits_get_repeat_count (&num, "AREASCAL", fp))
          return -1;

        if ((num != 1) && (num != h->nbins))
          return -1;
     }

   if (NULL == (a = (double *) ISIS_MALLOC (h->nbins * sizeof(double))))
     return -1;

   if (have_areascal_keyword)
     {
        for (i = 0; i < h->nbins; i++)
          {
             a[i] = areascal_keyword_value;
          }
     }
   else
     {
        if (-1 == cfits_read_optional_double_col (a, h->nbins, k, "AREASCAL", fp))
          {
             ISIS_FREE(a);
             return -1;
          }
     }

   for (i = 0; i < h->nbins; i++)
     {
        if (isfinite(a[i]) && a[i] != 0)
          {
             h->counts[i] /= a[i];
             h->stat_err[i] /= a[i];
          }
     }

   ISIS_FREE(a);
   return 0;
}

/*}}}*/

static int read_fits_bg_area_column (cfitsfile *fp, int k, Hist_t *h) /*{{{*/
{
   double *b = NULL;

   if (NULL == (b = (double *) ISIS_MALLOC (h->nbins * sizeof(double))))
     return -1;

   if (0 == cfits_read_optional_double_col (b, h->nbins, k, "BG_AREA", fp))
     {
        if (-1 == area_set (&h->bgd_area, b, h->nbins))
          {
             ISIS_FREE(b);
             return -1;
          }
     }

   ISIS_FREE(b);
   return 0;
}

/*}}}*/

static int read_fits_bg_counts_column (cfitsfile *fp, int k, Hist_t *h) /*{{{*/
{
   Area_Type *a;
   double *bg_counts, *area;

   if (NULL == (bg_counts = (double *) ISIS_MALLOC (h->nbins * sizeof(double))))
     return -1;

   if (-1 == cfits_read_double_col (bg_counts, h->nbins, k, "BG_COUNTS", fp))
     {
        ISIS_FREE(bg_counts);
        return -1;
     }

   a = &h->bgd_area;
   area = a->is_vector ? a->value.v : &a->value.s;
   if (-1 == Hist_define_background (h, h->exposure, area, a->is_vector,
                                     bg_counts, h->nbins))
     {
        ISIS_FREE(bg_counts);
        return -1;
     }

   ISIS_FREE(bg_counts);
   return 0;
}

/*}}}*/

static int read_fits_background_updown_columns (cfitsfile *fp, int k, Hist_t *h) /*{{{*/
{
   Area_Type *a;
   double *up, *down, *area;
   int i, n = h->nbins;

   if (NULL == (up = (double *) ISIS_MALLOC (2*h->nbins * sizeof(double))))
     return -1;
   down = up + h->nbins;

   memset ((char *)up, 0, 2*h->nbins*sizeof(double));

   if (cfits_col_exist ("BACKGROUND_UP", fp))
     {
        double backscup;
        if (-1 == cfits_read_double_col (up, h->nbins, k, "BACKGROUND_UP", fp))
          {
             ISIS_FREE(up);
             return -1;
          }
        if (-1 == cfits_read_double_keyword (&backscup, "BACKSCUP", fp))
          {
             ISIS_FREE(up);
             return -1;
          }
        if (backscup != 0)
          {
             for (i = 0; i < h->nbins; i++)
               {
                  up[i] /= backscup;
               }
          }
     }
   if (cfits_col_exist ("BACKGROUND_DOWN", fp))
     {
        double backscdn;
        if (-1 == cfits_read_double_col (down, h->nbins, k, "BACKGROUND_DOWN", fp))
          {
             ISIS_FREE(up);
             return -1;
          }
        if (-1 == cfits_read_double_keyword (&backscdn, "BACKSCDN", fp))
          {
             ISIS_FREE(up);
             return -1;
          }
        if (backscdn != 0)
          {
             for (i = 0; i < h->nbins; i++)
               {
                  down[i] /= backscdn;
               }
          }
     }

   for (i = 0; i < n; i++)
     {
        up[i] += down[i];
     }

   a = &h->bgd_area;
   area = a->is_vector ? a->value.v : &a->value.s;
   if (-1 == Hist_define_background (h, h->exposure, area, a->is_vector,
                                     up, h->nbins))
     {
        ISIS_FREE(up);
        return -1;
     }

   ISIS_FREE(up);
   return 0;
}

/*}}}*/

static int is_whitespace (const char *s) /*{{{*/
{
   if (s == NULL)
     return 0;

   for ( ; *s != 0; s++)
     {
        if (0 == isspace(*s))
          return 0;
     }

   return 1;
}

/*}}}*/

static char *read_file_keyword (cfitsfile *fp, const char *keyname, const char *filename) /*{{{*/
{
   char buf[CFLEN_FILENAME];
   char *path=NULL, *dir=NULL, *slash=NULL;

   *buf = 0;

   if (0 != cfits_read_string_keyword (buf, keyname, fp)
       || (0 != is_whitespace (buf))
       || (0 == isis_strcasecmp (buf, "NONE")))
     return NULL;

#ifndef HAVE_UNISTD_H
   return isis_make_string (buf);
#else
   if ((0 == access (buf, F_OK))
       || (NULL == (slash = strrchr ((char *)filename, '/')))
       || (NULL == (dir = isis_make_string (filename))))
     {
        return isis_make_string (buf);
     }

   *(dir + (slash - filename)) = 0;
   path = isis_mkstrcat (dir, "/", buf, NULL);
   ISIS_FREE(dir);
   return path;
#endif
}

/*}}}*/

static Hist_t *_read_typeI_pha (char *pha_filename, double min_stat_err, int use_bkg_updown) /*{{{*/
{
   Keyword_t *keytable = Hist_Keyword_Table;
   Hist_t *h = NULL;
   cfitsfile *fp = NULL;
   double sys_err_keyword, areascal_keyword;
   int nbins, val_stat, val_flux, have_sys_err_keyword, have_areascal_keyword;
   int i, ret = -1;

   if (pha_filename == NULL)
     return NULL;

   if (NULL == (fp = cfits_open_file_readonly (pha_filename)))
     return NULL;

   if (-1 == cfits_move_to_matching_hdu (fp, Spectrum_Hdu_Names, Spectrum_Hdu_Names_Hook, NULL))
     {
        isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__,
                    "No recognized spectrum HDU found in %s", pha_filename);
        goto finish;
     }

   if (-1 == cfits_read_int_keyword (&nbins, "NAXIS2", fp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "NAXIS2 => %s", pha_filename);
        goto finish;
     }

   have_sys_err_keyword =
     (0 == cfits_read_double_keyword (&sys_err_keyword, "SYS_ERR", fp));
   have_areascal_keyword =
     (0 == cfits_read_double_keyword (&areascal_keyword, "AREASCAL", fp));

   if (NULL == (h = Hist_new_hist (nbins)))
     goto finish;

   h->min_stat_err = min_stat_err;

   (void) Key_read_header_keywords (fp, (char *)h, keytable, FITS_FILE);

   if (-1 == read_fits_backscal_keywords (fp, h))
     goto finish;

   if (Hist_Ignore_PHA_Response_Keywords == 0)
     {
        h->respfile = read_file_keyword (fp, "RESPFILE", pha_filename);
        h->ancrfile = read_file_keyword (fp, "ANCRFILE", pha_filename);
     }

   if (Hist_Ignore_PHA_Backfile_Keyword == 0)
     {
        h->bgd_file = read_file_keyword (fp, "BACKFILE", pha_filename);
     }

   h->spec_num = 1;

   if ((-1 == cfits_read_optional_double_col (h->stat_err, nbins, 1, "STAT_ERR", fp)
        || (-1 == cfits_read_optional_int_col (h->quality, nbins, 1, "QUALITY", fp))))
     goto finish;

   if (cfits_col_exist ("SYS_ERR", fp))
     {
        if (-1 == alloc_sys_err (h, nbins)
            || -1 == cfits_read_double_col (h->sys_err_frac, h->nbins, 1, "SYS_ERR", fp))
          goto finish;
     }
   else if (have_sys_err_keyword)
     {
        if (-1 == alloc_sys_err (h, nbins)
            || -1 == assign_sys_err (h, &sys_err_keyword, 1))
          goto finish;
     }

   /* supercede backscal keywords */
   if (cfits_col_exist ("BACKSCAL", fp))
     {
        if (-1 == read_fits_backscal_column (fp, 1, h))
          goto finish;
     }

   if (cfits_col_exist ("BG_AREA", fp))
     {
        if (-1 == read_fits_bg_area_column (fp, 1, h))
          goto finish;
     }

   if (cfits_col_exist ("BG_COUNTS", fp))
     {
        if (-1 == read_fits_bg_counts_column (fp, 1, h))
          goto finish;
     }
   else if (use_bkg_updown
            && (cfits_col_exist ("BACKGROUND_UP", fp)
                || cfits_col_exist ("BACKGROUND_DOWN", fp)))
     {
        if (-1 == read_fits_background_updown_columns (fp, 1, h))
          goto finish;
     }

   if (cfits_col_exist ("RATE", fp))
     {
        if (-1 == cfits_read_double_col (h->counts, nbins, 1, "RATE", fp))
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "RATE => %s", pha_filename);
             goto finish;
          }
        for (i = 0; i < nbins; i++)
          {
             h->counts[i] *= h->exposure;
             h->stat_err[i] *= h->exposure;
          }
     }
   else if (cfits_col_exist ("COUNTS", fp))
     {
        if (-1 == cfits_read_double_col (h->counts, nbins, 1, "COUNTS", fp))
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "COUNTS => %s", pha_filename);
             goto finish;
          }
     }
   else
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "no COUNTS or RATE column in %s",
                    pha_filename);
        goto finish;
     }

   if (have_areascal_keyword || cfits_col_exist ("AREASCAL", fp))
     {
        if (-1 == apply_areascale (fp, 1, h, areascal_keyword, have_areascal_keyword))
          goto finish;
     }

   /* Load the data grid if available */
   if (cfits_col_exist ("BIN_LO", fp)
       && cfits_col_exist ("BIN_HI", fp))
     {
        char bin_units[CFLEN_VALUE];
        int input_units;
        if (-1 == cfits_get_colunits (bin_units, "BIN_LO", fp))
          input_units = U_ANGSTROM;
        else if (-1 == (input_units = unit_id (bin_units)))
          goto finish;

        if (-1 == cfits_read_double_col (h->bin_lo, nbins, 1, "BIN_LO", fp)
            || -1 == cfits_read_double_col (h->bin_hi, nbins, 1, "BIN_HI", fp))
          goto finish;

        if (-1 == get_canonical_hist_coordinates (h, input_units))
          goto finish;
     }

   memset ((char *)h->flux, 0, nbins*sizeof(double));
   memset ((char *)h->flux_err, 0, nbins*sizeof(double));

   val_stat = validate_stat_err (h);
   val_flux = validate_flux_err (h);
   if (-1 == val_stat || -1 == val_flux)
     goto finish;
   else if (1 == val_stat || 1 == val_flux)
     invalid_uncertainties_replaced ();

   (void) do_instrument_specific_hacks (h);

   ret = 0;
   finish:

   cfits_close_file (fp);
   if (ret)
     {
        free_hist (h);
        return NULL;
     }

   return h;
}

/*}}}*/

static int read_typeI_pha (Hist_t *head, Isis_Arf_t *arf_head, /*{{{*/
                           Isis_Rmf_t *rmf_head, char *file,
                           int *hist_index, int strict,
                           double min_stat_err, int use_bkg_updown)
{
   Hist_t *h;
   int id;

   (void) strict;

   if (NULL == (h = _read_typeI_pha (file, min_stat_err, use_bkg_updown)))
     return -1;

   if (h->respfile)
     {
        if (-1 == get_grid_from_rmf (h, rmf_head))
          {
             free_hist (h);
             return -1;
          }
     }

   if (h->respfile && h->bgd_file)
     {
        /* try loading the background only if the RMF has defined the grid */
        if (-1 == load_background_from_file (h, h->bgd_file))
          {
             isis_vmesg (WARN, I_ERROR, __FILE__, __LINE__,
                         "failed loading BACKFILE=%s",
                         h->bgd_file ? h->bgd_file : "<null");
             free_hist (h);
             return -1;
          }
     }

   if (NULL == (h->file = isis_make_string (file)))
     {
        free_hist (h);
        return -1;
     }

   /* finish initializing the histogram structure */
   if ((id = histogram_list_append (head, h)) < 0)
     return -1;

   /* now assign the ARF if there is one. */
   if (h->ancrfile)
     {
        Isis_Arf_t *a;
        int arf_index;

        arf_index = Arf_read_arf (arf_head, h->ancrfile);
        if ((arf_index > 0)
            && (NULL != (a = Arf_find_arf_index (arf_head, arf_index))))
          {
             /* Because 'h' has been fully initialized,
              * we know we've got a valid RMF */
             if (-1 == assign_matching_rsp (h, &h->a_rsp.rmf, &a))
               {
                  isis_vmesg (WARN, I_ERROR, __FILE__, __LINE__,
                              "ARF not compatible with RMF: %s", h->ancrfile);
               }
          }
     }

   *hist_index = id;

   return 0;
}

/*}}}*/

static int load_background_from_file (Hist_t *h, char *file) /*{{{*/
{
   Hist_t *b;
   cfitsfile *fp;
   double *bgd_area;
   int area_is_vector;

   if (h == NULL)
     return -1;

   /* recognize special file-names */
   if (file[0] == COMMENT_CHAR)
     return 0;

   /* First try opening as FITS, then as ASCII */
   if (NULL == (fp = cfits_open_file_readonly_silent (file)))
     {
        if (NULL == (b = read_ascii (file, h->min_stat_err)))
          return -1;
     }
   else
     {
        cfits_close_file (fp);
        if (NULL == (b = _read_typeI_pha (file, h->min_stat_err, 0)))
          return -1;
        if (-1 == set_hist_grid_using_rmf (h->a_rsp.rmf, b))
          {
             free_hist (b);
             return -1;
          }
     }

   if (b->area.is_vector)
     {
        bgd_area = b->area.value.v;
        area_is_vector = 1;
     }
   else
     {
        bgd_area = &b->area.value.s;
        area_is_vector = 0;
     }

   /* Define the background using the counts vector */
   if (-1 == Hist_define_background (h, b->exposure, bgd_area, area_is_vector,
                                     b->counts, b->orig_nbins))
     {
        free_hist (b);
        return -1;
     }

   free_hist (b);
   return 0;
}

/*}}}*/

int Hist_set_background_name (Hist_t *h, char *name) /*{{{*/
{
   char *s;

   if (name == NULL)
     {
        ISIS_FREE(h->bgd_file);
        return 0;
     }

   /* use string copy in case file == h->bgd_file */
   if (NULL == (s = isis_make_string (name)))
     return -1;
   ISIS_FREE (h->bgd_file);
   h->bgd_file = s;
   return 0;
}

/*}}}*/

int Hist_set_background_from_file (Hist_t *h, char *file) /*{{{*/
{
   if (h == NULL)
     return -1;

   if (-1 == load_background_from_file (h, file))
     return -1;

   return Hist_set_background_name (h, file);
}

/*}}}*/

static int read_typeII_pha (Hist_t *head, char * filename, int **indices, int *num_spectra, int just_one, /*{{{*/
                            double min_stat_err, int use_bkg_updown)
{
   Keyword_t *keytable = Hist_Keyword_Table;
   Hist_t *h1 = NULL;
   Hist_t *h = NULL;
   cfitsfile *cfp = NULL;
   char bin_units[CFLEN_VALUE];
   int k, num, nbins, reset = 0;
   int have_backscal_col, have_bg_area_col, have_bg_counts_col;
   int have_bkg_up_col, have_bkg_down_col;
   int have_areascal_col, have_rate, have_bin_lohi, have_exposure_col;
   int have_sys_err_col, have_sys_err_keyword, have_areascal_keyword;
   int input_units;
   double sys_err_keyword, areascal_keyword;
   char *s;
   int ret = -1;

   if (filename == NULL)
     return -1;

   if (NULL == (cfp = cfits_open_file_readonly_silent (filename)))
     return NOT_FITS_FORMAT;

   if (-1 == cfits_move_to_matching_hdu (cfp, Spectrum_Hdu_Names, Spectrum_Hdu_Names_Hook, NULL))
     {
        isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__,
                    "No recognized spectrum HDU found in %s", filename);
        goto finish;
     }

   if (-1 == cfits_get_colunits (bin_units, "BIN_LO", cfp))
     input_units = U_ANGSTROM;
   else if (-1 == (input_units = unit_id (bin_units)))
     goto finish;

   have_sys_err_keyword =
     (0 == cfits_read_double_keyword (&sys_err_keyword, "SYS_ERR", cfp));
   have_areascal_keyword =
     (0 == cfits_read_double_keyword (&areascal_keyword, "AREASCAL", cfp));

   if (just_one)
     *num_spectra = 1;
   else if (-1 == cfits_read_int_keyword (num_spectra, "NAXIS2", cfp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "NAXIS2 => %s", filename);
        goto finish;
     }

   if (NULL == (*indices = (int *) ISIS_MALLOC (*num_spectra * sizeof(int))))
     goto finish;
   memset ((char *) *indices, 0, (*num_spectra) * sizeof(int));

   /* Get number of elements from first column, and assume all data */
   /* columns have the same number of elements */

   have_exposure_col = cfits_col_exist ("EXPOSURE", cfp);
   have_bin_lohi = cfits_col_exist ("BIN_LO", cfp);
   have_rate = cfits_col_exist ("RATE", cfp);
   have_sys_err_col = cfits_col_exist ("SYS_ERR", cfp);

   have_backscal_col = cfits_col_exist ("BACKSCAL", cfp);
   have_bg_counts_col = cfits_col_exist ("BG_COUNTS", cfp);
   have_bg_area_col = cfits_col_exist ("BG_AREA", cfp);

   have_bkg_up_col = cfits_col_exist ("BACKGROUND_UP", cfp);
   have_bkg_down_col = cfits_col_exist ("BACKGROUND_DOWN", cfp);

   have_areascal_col = cfits_col_exist ("AREASCAL", cfp);

   s = have_rate ? (char *) "RATE" : (char *) "COUNTS";
   if (-1 == cfits_get_repeat_count(&nbins, s, cfp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "%s repeat count => %s",
                    s, filename);
        goto finish;
     }

   if ((!just_one) && (Isis_Verbose >= WARN))
     fputs ("Reading: ", stderr);

   k = just_one ? just_one : 1;

   for (num = 0; num < *num_spectra; num++)
     {
        int val_stat, val_flux;

        if (NULL == (h = Hist_new_hist (nbins)))
          goto finish;

        h->min_stat_err = min_stat_err;

        /* read keywords first pass only */
        if (num > 0)
          Hist_copy_histogram_keywords (h, h1);
        else
          {
             h1 = h;
             (void) Key_read_header_keywords (cfp, (char *)h, keytable, FITS_FILE);
             if (have_backscal_col == 0)
               {
                  if (-1 == read_fits_backscal_keywords (cfp, h))
                    goto finish;
               }
          }

        if (-1 == cfits_read_int_col (&h->spec_num, 1, k, "SPEC_NUM", cfp)
            || -1 == cfits_read_optional_int_col (&h->order, 1, k, "TG_M", cfp)
            || -1 == cfits_read_optional_int_col (&h->part, 1, k, "TG_PART", cfp)
            || -1 == cfits_read_optional_int_col (&h->srcid, 1, k, "TG_SRCID", cfp)
            || -1 == cfits_read_optional_int_col (h->quality, nbins, k, "QUALITY", cfp)
            || -1 == cfits_read_optional_double_col (h->stat_err, nbins, k, "STAT_ERR", cfp)
            || -1 == cfits_read_optional_double_col (h->flux, nbins, k, "FLUX", cfp)
            || -1 == cfits_read_optional_double_col (h->flux_err, nbins, k, "FLUX_ERR", cfp)
            )
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "%s", filename);
             goto finish;
          }

        if (have_sys_err_col)
          {
             if (-1 == alloc_sys_err (h, nbins)
                 || -1 == cfits_read_double_col (h->sys_err_frac, h->nbins, k, "SYS_ERR", cfp))
               goto finish;
          }
        else if (have_sys_err_keyword)
          {
             if (-1 == alloc_sys_err (h, nbins)
                 || -1 == assign_sys_err (h, &sys_err_keyword, 1))
               goto finish;
          }

        if (have_exposure_col)
          {
             if (-1 == cfits_read_double_col (&h->exposure, 1, k, "EXPOSURE", cfp))
               {
                  isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "%s", filename);
                  goto finish;
               }
          }

        if (have_bin_lohi)
          {
             if ((-1 == cfits_read_double_col (h->bin_lo, nbins, k, "BIN_LO", cfp))
                 || (-1 == cfits_read_double_col (h->bin_hi, nbins, k, "BIN_HI", cfp)))
               {
                  isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "%s", filename);
                  goto finish;
               }
          }

        if (have_rate)
          {
             int i;
             if (-1 == cfits_read_double_col (h->counts, nbins, 1, "RATE", cfp))
               {
                  isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "RATE => %s", filename);
                  goto finish;
               }
             for (i = 0; i < nbins; i++)
               {
                  h->counts[i] *= h->exposure;
                  h->stat_err[i] *= h->exposure;
               }
          }
        else if (-1 == cfits_read_double_col (h->counts, nbins, k, "COUNTS", cfp))
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "COUNTS => %s", filename);
             goto finish;
          }

        if (have_backscal_col)
          {
             if (-1 == read_fits_backscal_column (cfp, k, h))
               goto finish;
          }

        if (have_bg_area_col)
          {
             if (-1 == read_fits_bg_area_column (cfp, k, h))
               goto finish;
          }

        if (have_bg_counts_col)
          {
             if (-1 == read_fits_bg_counts_column (cfp, k, h))
               goto finish;
          }
        else if (use_bkg_updown && (have_bkg_up_col || have_bkg_down_col))
          {
             if (-1 == read_fits_background_updown_columns (cfp, k, h))
               goto finish;
          }

        if (have_areascal_col || have_areascal_keyword)
          {
             if (-1 == apply_areascale (cfp, k, h, areascal_keyword, have_areascal_keyword))
               goto finish;
          }

        if (have_bin_lohi)
          {
             if (-1 == get_canonical_hist_coordinates (h, input_units))
               goto finish;
          }

        val_stat = validate_stat_err (h);
        val_flux = validate_flux_err (h);
        if (-1 == val_stat || -1 == val_flux)
          goto finish;
        else if (1 == val_stat || 1 == val_flux)
          reset = 1;

        (void) do_instrument_specific_hacks (h);

        if ((!just_one) && (Isis_Verbose >= WARN))
          fputc ('.', stderr);

        if (NULL == (h->file = isis_make_string (filename)))
          goto finish;

        if (-1 == ((*indices)[num] = histogram_list_append (head, h)))
          goto finish;

        k++;
     }

   ret = 0;

   finish:

   if ((!just_one) && (Isis_Verbose >= WARN))
     fputc ('\n',stderr);

   if (reset)
     invalid_uncertainties_replaced ();

   if (ret)
     {
        free_hist (h);
        ISIS_FREE (*indices);
     }

   (void) cfits_close_file (cfp);
   return ret;
}

/*}}}*/

static int get_pha_type (cfitsfile *fp) /*{{{*/
{
   int repeat_count;

   if (fp == NULL)
     return -1;

   if (-1 == cfits_move_to_matching_hdu (fp, Spectrum_Hdu_Names, Spectrum_Hdu_Names_Hook, NULL))
     {
        isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__,
                    "File contains no recognized spectrum HDU");
        return -1;
     }

   if (cfits_col_exist ("COUNTS", fp))
     {
        if (-1 == cfits_get_repeat_count(&repeat_count, "COUNTS", fp))
          repeat_count = 1;
     }
   else if (cfits_col_exist ("RATE", fp))
     {
        if (-1 == cfits_get_repeat_count(&repeat_count, "RATE", fp))
          repeat_count = 1;
     }

   if (repeat_count > 1)
     return PHA_TYPE_II;

   if (repeat_count == 1)
     return PHA_TYPE_I;

   return -1;
}

/*}}}*/

int Hist_read_fits (Hist_t *head, Isis_Arf_t *arf_head, Isis_Rmf_t *rmf_head, char * pha_filename, /*{{{*/
                    int **indices, int *num_spectra, int strict, int just_one,
                    double min_stat_err, int use_bkg_updown)
{
   cfitsfile *cfp = NULL;
   int pha_type;

   if (pha_filename == NULL)
     return -1;

   if (NULL == (cfp = cfits_open_file_readonly (pha_filename)))
     return NOT_FITS_FORMAT;

   pha_type = get_pha_type (cfp);

   (void) cfits_close_file (cfp);

   switch (pha_type)
     {
      case PHA_TYPE_I:
        *num_spectra = 1;
        if (NULL == (*indices = (int *) ISIS_MALLOC (sizeof(int))))
          return -1;
        return read_typeI_pha (head, arf_head, rmf_head, pha_filename, *indices, strict,
                               min_stat_err, use_bkg_updown);

      case PHA_TYPE_II:
        return read_typeII_pha (head, pha_filename, indices, num_spectra, just_one,
                                min_stat_err, use_bkg_updown);

      default:
        break;
     }

   /* should never happen */
   return -1;
}

/*}}}*/

/* get / copy */

static Isis_Arf_t *identity_arf_on_arf_grid (Hist_t *h) /*{{{*/
{
   Isis_Arf_t *arf = h->a_rsp.arf;
   return Arf_make_identity_arf (arf->bin_lo, arf->bin_hi, arf->nbins);
}

/*}}}*/

static Isis_Rmf_t *matching_identity_rmf (Hist_t *h, Isis_Arf_t *a) /*{{{*/
{
   Isis_Rmf_Grid_Type arf, ebounds;

   if ((h == NULL) || (a == NULL))
     return NULL;

   ebounds.nbins = h->orig_nbins;
   ebounds.bin_lo = h->orig_bin_lo;
   ebounds.bin_hi = h->orig_bin_hi;
   ebounds.units = U_ANGSTROM;

   arf.nbins = a->nbins;
   arf.bin_hi = a->bin_hi;
   arf.bin_lo = a->bin_lo;
   arf.units = U_ANGSTROM;

   return Rmf_create_delta_rmf (&arf, &ebounds);
}

/*}}}*/

int Hist_get_data_region_area (Hist_t *h, double **area, int *num) /*{{{*/
{
   if (h == NULL || area == NULL || num == NULL)
     return -1;

   return area_get (&h->area, area,  num);
}

/*}}}*/

int Hist_get_back_region_area (Hist_t *h, double **area, int *num) /*{{{*/
{
   if (h == NULL || area == NULL || num == NULL)
     return -1;

   return area_get (&h->bgd_area, area,  num);
}

/*}}}*/

int Hist_set_data_region_area (Hist_t *h, double *area, int num) /*{{{*/
{
   if (h == NULL)
     return -1;

   if ((num > 1) && (num != h->orig_nbins))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Background vector area has size %d, should be %d",
                    num, h->orig_nbins);
        return -1;
     }

   return area_set (&h->area, area, num);
}

/*}}}*/

int Hist_set_back_region_area (Hist_t *h, double *area, int num) /*{{{*/
{
   if (h == NULL)
     return -1;

   if ((num > 1) && (num != h->orig_nbins))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Background vector area has size %d, should be %d",
                    num, h->orig_nbins);
        return -1;
     }

   return area_set (&h->bgd_area, area, num);
}

/*}}}*/

int Hist_set_frame_time (Hist_t *h, double frame_time) /*{{{*/
{
   if (NULL == h)
     return -1;

   h->frame_time = frame_time;

   return 0;
}

/*}}}*/

int Hist_get_frame_time (Hist_t *h, double *frame_time) /*{{{*/
{
   if (NULL == h)
     return -1;

   *frame_time = h->frame_time;

   return 0;
}

/*}}}*/

/* instrumental background hook */

char *Hist_get_instrumental_background_hook_name (Hist_t *h) /*{{{*/
{
   if (h == NULL) return NULL;
   return h->instrumental_background_hook_name;
}

/*}}}*/

int Hist_set_instrumental_background_hook_name (Hist_t *h, char *name) /*{{{*/
{
   char *s;

   if (NULL == h)
     return -1;

   if ((name == NULL) || (*name == 0))
     s = NULL;
   else
     {
        if (NULL == (s = isis_make_string (name)))
          return -1;
     }

   ISIS_FREE (h->instrumental_background_hook_name);
   h->instrumental_background_hook_name = s;

   /* old function pointer is no longer valid */
   h->instrumental_background_hook = NULL;

   return 0;
}

/*}}}*/

SLang_Name_Type *Hist_get_instrumental_background_hook (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return NULL;
   return h->instrumental_background_hook;
}

/*}}}*/

int Hist_set_instrumental_background_hook (Hist_t *h, SLang_Name_Type *hook) /*{{{*/
{
   if (h == NULL)
     return -1;

   h->instrumental_background_hook = hook;

   return 0;
}

/*}}}*/

/* user-defined metadata */

SLang_Any_Type *Hist_get_metadata (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return NULL;
   return h->user_meta;
}

/*}}}*/

int Hist_set_metadata (Hist_t *h, SLang_Any_Type *meta) /*{{{*/
{
   if (h == NULL)
     return -1;
   if (h->user_meta != NULL)
     SLang_free_anytype (h->user_meta);
   h->user_meta = meta;
   return 0;
}

/*}}}*/

/* assigned model */
SLang_Name_Type *Hist_assigned_model (Hist_t *h)
{
   if (h == NULL)
     return NULL;
   return h->assigned_model;
}

Isis_Arg_Type *Hist_assigned_model_args (Hist_t *h)
{
   if (h == NULL)
     return NULL;
   return h->assigned_model_args;
}

int Hist_assign_model (Hist_t *h, SLang_Name_Type *fun_ptr, Isis_Arg_Type *args)
{
   if (h == NULL)
     return -1;
   SLang_free_function (h->assigned_model);
   isis_free_args (h->assigned_model_args);
   h->assigned_model = fun_ptr;
   h->assigned_model_args = args;
   return 0;
}

/* post model-evaluation  hook */

SLang_Name_Type *Hist_post_model_hook (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return NULL;
   return h->post_model_hook;
}

/*}}}*/

int Hist_set_post_model_hook (Hist_t *h, SLang_Name_Type *hook, void (*delete_hook)(SLang_Name_Type *)) /*{{{*/
{
   if (h == NULL)
     return -1;
   if (h->post_model_hook_delete)
     (*h->post_model_hook_delete) (h->post_model_hook);
   h->post_model_hook = hook;
   h->post_model_hook_delete = delete_hook;
   return 0;
}

/*}}}*/

int Hist_run_post_model_hook (Hist_t *h, double *cts, SLang_Array_Type *sl_bgd) /*{{{*/
{
   SLindex_Type orig_nbins;

   SLang_Array_Type *sl_cts = NULL;
   SLang_Array_Type *sl_lo = NULL;
   SLang_Array_Type *sl_hi = NULL;

   if ((h == NULL) || (h->post_model_hook == NULL))
     return -1;

   orig_nbins = h->orig_nbins;

   if ((NULL == (sl_cts = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &orig_nbins, 1)))
       || (NULL == (sl_lo = SLang_create_array (SLANG_DOUBLE_TYPE, 1, h->orig_bin_lo, &orig_nbins, 1)))
       || (NULL == (sl_hi = SLang_create_array (SLANG_DOUBLE_TYPE, 1, h->orig_bin_hi, &orig_nbins, 1))))
     return -1;

   memcpy ((char *)sl_cts->data, (char *)cts, h->orig_nbins * sizeof(double));

   SLang_start_arg_list ();
   SLang_push_array (sl_lo, 0);
   SLang_push_array (sl_hi, 0);
   SLang_push_array (sl_cts, 1);
   if (sl_bgd != NULL)
     {
        SLang_push_array (sl_bgd, 0);
     }
   else
     {
        SLang_push_double (0.0);
     }
   SLang_end_arg_list ();

   SLexecute_function (h->post_model_hook);

   sl_lo->data = NULL; SLang_free_array (sl_lo);
   sl_hi->data = NULL; SLang_free_array (sl_hi);

   return Isis_pop_double_array (cts, h->orig_nbins);
}

/*}}}*/

/* access */

int Hist_get_rsp_list (Hist_t *h, int **arfs, int **rmfs, int *num) /*{{{*/
{
   Isis_Rsp_t *rsp;
   int n, *a, *r;

   *arfs = NULL;
   *rmfs = NULL;
   *num = 0;

   if (h == NULL)
     return -1;

   n = 0;
   for (rsp = &h->a_rsp; rsp != NULL; rsp = rsp->next)
     n++;

   if (n < 1)
     return -1;

   a = (int *) ISIS_MALLOC (n * sizeof(int));
   r = (int *) ISIS_MALLOC (n * sizeof(int));
   if ((a == NULL) || (r == NULL))
     {
        ISIS_FREE (a);
        ISIS_FREE (r);
        return -1;
     }

   n = 0;
   for (rsp = &h->a_rsp; rsp != NULL; rsp = rsp->next)
     {
        a[n] = rsp->arf->index;
        r[n] = rsp->rmf->index;
        n++;
     }

   *num = n;
   *arfs = a;
   *rmfs = r;

   return 0;
}

/*}}}*/

int Hist_get_info (Hist_t *h, Hist_Info_Type *info) /*{{{*/
{
   if ((NULL == h) || (info == NULL))
     return -1;

   info->tstart = h->tstart;
   info->min_stat_err = h->min_stat_err;
   info->spec_num = h->spec_num;
   info->order = h->order;
   info->part = h->part;
   info->srcid = h->srcid;
   info->exclude = h->exclude;
   info->combo_id = h->combo_id;
   info->combo_weight = h->combo_weight;
   isis_strcpy (info->object, h->object, sizeof(info->object));
   isis_strcpy (info->instrument, h->instrument, sizeof(info->instrument));
   isis_strcpy (info->grating, h->grating, sizeof(info->grating));
   info->file = h->file ? h->file : (char *) "";
   info->bgd_file = h->bgd_file ? h->bgd_file : (char *) "";

   return 0;
}

/*}}}*/

int Hist_set_info (Hist_t *h, Hist_Info_Type *info) /*{{{*/
{
   char *s;

   if ((NULL == h) || (info == NULL))
     return -1;

   h->tstart = info->tstart;
   h->min_stat_err = info->min_stat_err;
   h->spec_num = info->spec_num;
   h->order = info->order;
   h->part = info->part;
   h->srcid = info->srcid;
   h->exclude = info->exclude;
   h->combo_id = info->combo_id;
   h->combo_weight = info->combo_weight;

   if (info->file)
     {
        if (NULL == (s = isis_make_string (info->file)))
          return -1;
        ISIS_FREE(h->file);
        h->file = s;
     }

   if (info->bgd_file)
     {
        if (-1 == Hist_set_background_name (h, info->bgd_file))
          return -1;
     }

   return 0;
}

/*}}}*/

int Hist_set_object_name (Hist_t *h, char *object) /*{{{*/
{
   if (h == NULL || object == NULL)
     return -1;

   isis_strcpy (h->object, object, sizeof(h->object));
   return 0;
}

/*}}}*/

int Hist_set_fake (Hist_t *h, int value) /*{{{*/
{
   if (NULL == h)
     return -1;

   h->is_fake_data = value ? 1 : 0;

   return 0;
}

/*}}}*/

int Hist_is_fake (Hist_t *h) /*{{{*/
{
   if (NULL == h)
     return -1;

   return h->is_fake_data;
}

/*}}}*/

int Hist_get_back_exposure (Hist_t *h, double *back_exposure) /*{{{*/
{
   if (NULL == h)
     return -1;

   *back_exposure = h->bgd_exposure;

   return 0;
}
/*}}}*/

int Hist_set_back_exposure (Hist_t *h, double back_exposure) /*{{{*/
{
   if (NULL == h)
     return -1;

   h->bgd_exposure = back_exposure;

   return 0;
}
/*}}}*/

int Hist_get_exposure (Hist_t *h, double * exposure) /*{{{*/
{
   if (NULL == h)
     return -1;

   *exposure = h->exposure;

   return 0;
}
/*}}}*/

int Hist_set_exposure (Hist_t *h, double exposure) /*{{{*/
{
   if (NULL == h)
     return -1;

   h->exposure = exposure;

   return 0;
}

/*}}}*/

int Hist_set_min_stat_err (Hist_t *h, double min_stat_err) /*{{{*/
{
   if (NULL == h)
     return -1;

   h->min_stat_err = min_stat_err;

   return 0;
}

/*}}}*/

int Hist_get_min_stat_err (Hist_t *h, double *min_stat_err) /*{{{*/
{
   if (h == NULL || min_stat_err == NULL)
     return -1;

   *min_stat_err = h->min_stat_err;

   return 0;
}

/*}}}*/

int Hist_hist_size (Hist_t *h, unsigned int version) /*{{{*/
{
   Isis_Hist_t hist;

   if (NULL == h)
     return -1;

   if (-1 == get_hist_version (h, version, &hist))
     return -1;

   return hist.nbins;
}

/*}}}*/

int Hist_orig_hist_size (Hist_t *h) /*{{{*/
{
   if (NULL == h)
     return -1;

   return h->orig_nbins;
}

/*}}}*/

int Hist_get_index (Hist_t *h) /*{{{*/
{
   return ((h == NULL) ? -1 : h->index);
}

/*}}}*/

int _Hist_get_orig_hist_grid (Hist_t *h, Isis_Hist_t *g) /*{{{*/
{
   if ((h == NULL) || (g == NULL))
     return -1;

   g->nbins = h->orig_nbins;
   g->bin_lo = h->orig_bin_lo;
   g->bin_hi = h->orig_bin_hi;
   g->notice_list = NULL;
   g->notice = NULL;
   g->n_notice = 0;

   return 0;
}

/*}}}*/

int Hist_get_hist_grid (Hist_t *h, unsigned int version, Isis_Hist_t *g) /*{{{*/
{
   Isis_Hist_t hist;

   if (-1 == get_hist_version (h, version, &hist))
     return -1;

   g->nbins = hist.nbins;
   g->bin_lo = hist.bin_lo;
   g->bin_hi = hist.bin_hi;
   g->notice_list = hist.notice_list;
   g->notice = hist.notice;
   g->n_notice = hist.n_notice;

   return 0;
}

/*}}}*/

int Hist_get_hist_rebin_info (Hist_t *h, int **rebin, int *orig_nbins) /*{{{*/
{
   *rebin = NULL;
   *orig_nbins = 0;

   if (h == NULL)
     return -1;

   *rebin = h->rebin;
   *orig_nbins = h->orig_nbins;
   return 0;
}

/*}}}*/

int Hist_copy_input_background (Hist_t *h, int do_rebin, double **bgd, int *nbins) /*{{{*/
{
   if (h == NULL)
     return -1;
   *bgd = NULL;
   *nbins = 0;
   if (h->orig_bgd == NULL)
     return 0;
   return copy_input_background (h, do_rebin, bgd, nbins);
}

/*}}}*/

int Hist_background_scale_factor (Hist_t *h, int do_rebin, double **scale_factor, int *nbins) /*{{{*/
{
   if (h == NULL)
     return -1;
   *nbins = do_rebin ? h->nbins : h->orig_nbins;
   *scale_factor = background_scale_factor (h, do_rebin);
   return (scale_factor != NULL) ? 0 : -1;
}

/*}}}*/

int Hist_copy_scaled_background (Hist_t *h, double **bgd) /*{{{*/
{
   if ((h == NULL) || (bgd == NULL))
     return -1;
   *bgd = NULL;
   if (h->orig_bgd == NULL)
     return 0;
   return scale_background (h, 0, bgd, NULL);
}

/*}}}*/

int Hist_get_model_grid (Isis_Hist_t *g, Hist_t *h) /*{{{*/
{
   double *val;

   if (h == NULL || g == NULL)
     return -1;

   val = g->val;
   *g = h->model_flux;  /* struct copy */
   g->val = val;

   return 0;
}

/*}}}*/

int Hist_copy_hist_data (Hist_t *h, unsigned int version, double * val, double * val_err) /*{{{*/
{
   Isis_Hist_t hist;

   if (NULL == h)
     return -1;

   if (-1 == get_hist_version (h, version, &hist))
     return -1;

   if (val != NULL && hist.val != NULL)
     memcpy ((char *)val, (char *)hist.val, hist.nbins * sizeof(double));

   if (val_err != NULL && hist.val_err != NULL)
     memcpy ((char *)val_err, (char *)hist.val_err, hist.nbins * sizeof(double));

   return 0;
}

/*}}}*/

int Hist_copy_packed_model (Hist_t *h, unsigned int version, double *val) /*{{{*/
{
   Isis_Hist_t hist;
   int i;

   if ((NULL == h) || (val == NULL))
     return -1;

   if (-1 == get_hist_version (h, version, &hist))
     return -1;

   for (i = 0; i < hist.n_notice; i++)
     {
        int k = hist.notice_list[i];
        val[i] = hist.val[k];
     }

   return 0;
}

/*}}}*/

int Hist_copy_noticed_data (Hist_t *h, unsigned int version, double *val, double *val_err) /*{{{*/
{
   Isis_Hist_t hist;
   double *scaled_bgd = NULL;
   double *scaled_bgd_err = NULL;
   int i;

   if (is_model(version))
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__, NULL);
        return -1;
     }

   /* FIXME!? - background error propagation.
    *    h->stat_error_hook should manage background error propagation,
    *    so the interface should provide the necessary information.
    *    If h->instrumental_background_hook is defined, the background
    *    might change during the fit.  In that case, it might be
    *    better to re-do the background error propagation at each
    *    fit step.
    */

   if (h->stat_error_hook != NULL)
     {
        if (-1 == run_stat_error_hook (h))
          return -1;
     }

   if (val == NULL || val_err == NULL || h == NULL
       || -1 == get_hist_version (h, version, &hist))
     return -1;

   if ((h->stat_error_hook == NULL)
       && ((h->orig_bgd != NULL) && is_counts(version)))
     {
        /* Scaling the background here ensures it is up-to-date
         * for fitting in case the exposure times or extraction
         * area values have changed since it was initially defined.
         */
        if (-1 == scale_background (h, 1, &scaled_bgd, &scaled_bgd_err))
          {
             return -1;
          }
     }

   /* Background-subtracted, FLUX-CORRECTED data computed by ISIS
    * already includes the background subtraction uncertainty in flux_err,
    * so we don't need to include it now.  If the data was
    * flux-corrected in some other way, I assume the provided
    * flux_err value is correct and that any background provided
    * is exact.
    */

   if (scaled_bgd == NULL)
     {
        for (i=0; i < hist.n_notice; i++)
          {
             int k = hist.notice_list[i];
             val[i] = hist.val[k];
             val_err[i] = hist.val_err[k];
          }
     }
   else
     {
        for (i=0; i < hist.n_notice; i++)
          {
             int k = hist.notice_list[i];
             val[i] = hist.val[k];
             val_err[i] = isis_hypot (hist.val_err[k], scaled_bgd_err[k]);
          }

        ISIS_FREE(scaled_bgd);
        ISIS_FREE(scaled_bgd_err);
     }

   /* Include xspec-style systematic error */
   if (hist.sys_err_frac != NULL)
     {
        for (i=0; i < hist.n_notice; i++)
          {
             int k = hist.notice_list[i];
             val_err[i] = isis_hypot (val_err[i], hist.sys_err_frac[k] * hist.val[k]);
          }
     }

   return 0;
}

/*}}}*/

/* replace */

static int handle_grid_change (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return -1;

   *h->object = 0;
   *h->grating = 0;
   *h->instrument = 0;
   h->ancrfile = NULL;
   h->respfile = NULL;
   h->exposure = 1.0;

   h->frame_time = -1.0;
   h->timedel = -1.0;
   h->exptime = -1.0;

   h->spec_num = 0;
   h->order = 0;
   h->part = 0;
   h->srcid = 0;

   release_rsp (&h->a_rsp);
   release_rsp (&h->f_rsp);

   if (-1 == finish_hist_init (h))
     return -1;

   isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "reset data grid");

   return 0;
}

/*}}}*/

int Hist_replace_hist_grid (Hist_t *h, double * bin_lo, double * bin_hi, int nbins) /*{{{*/
{
   int i, size;

   if (NULL == h)
     return -1;

   /* should have nbins <= h->nbins <= h->orig_nbins */
   if (nbins > h->orig_nbins)
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__, "new size must be <= old size");
        return -1;
     }

   if (-1 == validate_wavelength_grid (nbins, bin_lo, bin_hi))
     {
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "Kept existing grid");
        return -1;
     }

   for (i = 0; i < nbins; i++)
     {
        if ( (fabs(1.0 - h->bin_lo[i]/bin_lo[i]) > GRID_TOL)
             || (fabs(1.0 - h->bin_hi[i]/bin_hi[i]) > GRID_TOL) )
          break;
     }

   size = nbins * sizeof(double);
   memcpy ((char *)h->bin_lo, (char *)bin_lo, size);
   memcpy ((char *)h->bin_hi, (char *)bin_hi, size);

   if (i == h->nbins)
     return 0;

   h->nbins = nbins;
   return handle_grid_change (h);
}

/*}}}*/

int Hist_replace_hist_data (Hist_t *h, unsigned int version, double *val, double *val_err, int nbins) /*{{{*/
{
   Isis_Hist_t hist;
   int valid, size, i;

#if 0
   if (is_model(version))
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__, NULL);
        return -1;
     }
#endif

   if (NULL == h)
     return -1;

   /* If the grid just changed, we should have
    *   nbins = h->nbins = h->orig_nbins
    * otherwise, it may be that
    *   nbins = h->nbins <= h->orig_nbins
    */
   if (nbins != h->nbins)
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__,
                    "new size must match current size");
        return -1;
     }

   if (-1 == get_hist_version (h, version, &hist))
     return -1;

   size = h->nbins * sizeof(double);

   if (val == NULL || hist.val == NULL)
     return -1;
   memcpy ((char *)hist.val, (char *)val, size);

   if (val_err == NULL || hist.val_err == NULL)
     return 0;

   memcpy ((char *)hist.val_err, (char *)val_err, size);

   if (is_flux(version))
     valid = validate_flux_err (h);
   else
     valid = validate_stat_err (h);

   if (valid == -1)
     return -1;
   else if (valid == 1)
     invalid_uncertainties_replaced ();

   if (is_counts(version))
     {
        double tot, *oc = h->orig_counts;
        memcpy ((char *)h->orig_counts, (char *)h->counts, size);
        memcpy ((char *)h->orig_stat_err, (char *)h->stat_err, size);

        tot = 0.0;
        for (i = 0; i < h->nbins; i++)
          {
             tot += oc[i];
          }
        h->totcts = tot;
     }
   else if (is_flux(version))
     {
        memcpy ((char *)h->orig_flux, (char *)h->flux, size);
        memcpy ((char *)h->orig_flux_err, (char *)h->flux_err, size);
        if (h->flux_weights != NULL)
          {
             for (i = 0; i < h->orig_nbins; i++)
               {
                  h->flux_weights[i] = 1.0;
               }
          }
     }

   return 0;
}

/*}}}*/

static int copy_d_array (double **to, double *from, int nbins) /*{{{*/
{
   double *cpy;

   if (NULL == (cpy = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     return -1;

   memcpy ((char *)cpy, (char *)from, nbins * sizeof(double));

   ISIS_FREE (*to);
   *to = cpy;

   return 0;
}

/*}}}*/

static int rebin_rmf (Hist_t *h, double *lo, double *hi, int nbins) /*{{{*/
{
   Isis_Rmf_t *r;
   int ret;

   /* this will get re-set by finish_hist_init() */
   release_rmf (h->f_rsp.rmf);
   release_arf (h->f_rsp.arf);

   /* rmf rebinning is allowed only if the RMF is not in use */
   r = h->a_rsp.rmf;
   r->ref_count--;
   ret = Rmf_rebin_rmf (r, lo, hi, nbins);
   r->ref_count++;

   return ret;
}

/*}}}*/

static int rebin_backscale (Area_Type *at, Hist_t *h, double *lo, double *hi, int nbins) /*{{{*/
{
   double *v;

   if (-1 == area_force_vector (at, h->nbins))
     return -1;

   if (NULL == (v = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     return -1;

   if (-1 == rebin_histogram (at->value.v, h->bin_lo, h->bin_hi, h->nbins,
                              v, lo, hi, nbins))
     {
        ISIS_FREE(v);
        return -1;
     }

   ISIS_FREE(at->value.v);
   at->value.v = v;
   at->n = nbins;

   return 0;
}

/*}}}*/

static double *rebin_sys_err_frac (Hist_t *h, double *lo, double *hi, int nbins) /*{{{*/
{
   double *sys_err_frac_dlam=NULL, *sys_err_frac=NULL;
   int i;

   if ((NULL == (sys_err_frac_dlam = (double *) ISIS_MALLOC (h->orig_nbins * sizeof(double))))
       || (NULL == (sys_err_frac = (double *) ISIS_MALLOC (nbins * sizeof(double)))))
     {
        ISIS_FREE(sys_err_frac_dlam);
        return NULL;
     }

   for (i = 0; i < h->orig_nbins; i++)
     {
        sys_err_frac_dlam[i] = (h->orig_sys_err_frac[i]
                                * (h->orig_bin_hi[i] - h->orig_bin_lo[i]));
     }

   if (-1 == rebin_histogram (sys_err_frac_dlam, h->orig_bin_lo, h->orig_bin_hi, h->orig_nbins,
                              sys_err_frac, lo, hi, nbins))
     {
        ISIS_FREE(sys_err_frac_dlam);
        ISIS_FREE(sys_err_frac);
        return NULL;
     }

   ISIS_FREE(sys_err_frac_dlam);

   /* use the mean sys_err within each new bin */
   for (i = 0; i < nbins; i++)
     {
        sys_err_frac[i] /= hi[i] - lo[i];
     }

   return sys_err_frac;
}

/*}}}*/

int Hist_rebin (Hist_t *h, double *lo, double *hi, int nbins) /*{{{*/
{
   double *bgd_cts = NULL;
   double *cts = NULL;
   double *sys_err_frac = NULL;
   union {double *d; int *i;} tmp;
   int i;

   if ((h == NULL) || (lo == NULL) || (hi == NULL) || (nbins < 1))
     return -1;

   if (h->nbins != h->orig_nbins)
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__,
                    "please revert to the original input data grid");
        return -1;
     }

   if (h->a_rsp.rmf->ref_count > 2)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "rebin failed: RMF in use by multiple datasets");
        return -1;
     }

   if (-1 == validate_wavelength_grid (nbins, lo, hi))
     return -1;

   if (h->sys_err_frac != NULL)
     {
        /* rebin this before we modify the original grid, counts, etc. */
        if (NULL == (sys_err_frac = rebin_sys_err_frac (h, lo, hi, nbins)))
          goto error_return;
     }

#define D_REALLOC_OR_FAIL(x,n) \
   if (NULL == (tmp.d = (double *) ISIS_REALLOC((x), (n) * sizeof(double)))) \
       goto error_return; \
   else (x) = tmp.d;

#define I_REALLOC_OR_FAIL(x,n) \
   if (NULL == (tmp.i = (int *) ISIS_REALLOC((x), (n) * sizeof(int)))) \
       goto error_return; \
   else (x) = tmp.i;

   I_REALLOC_OR_FAIL(h->notice, nbins);
   I_REALLOC_OR_FAIL(h->notice_list, nbins);
   D_REALLOC_OR_FAIL(h->model_counts, nbins);
   D_REALLOC_OR_FAIL(h->convolved_model_flux, nbins);
   D_REALLOC_OR_FAIL(h->model_flux.val, nbins);
   D_REALLOC_OR_FAIL(h->stat_err, nbins);
   D_REALLOC_OR_FAIL(h->flux, nbins);
   D_REALLOC_OR_FAIL(h->flux_err, nbins);
   D_REALLOC_OR_FAIL(h->orig_flux, nbins);
   D_REALLOC_OR_FAIL(h->orig_flux_err, nbins);
   if (h->flux_weights != NULL)
     {
        D_REALLOC_OR_FAIL(h->flux_weights, nbins);
     }
   if (h->sys_err_frac != NULL)
     {
        D_REALLOC_OR_FAIL(h->orig_sys_err_frac, nbins);
     }
   I_REALLOC_OR_FAIL(h->rebin, nbins);
   D_REALLOC_OR_FAIL(h->orig_bin_lo, nbins);
   D_REALLOC_OR_FAIL(h->orig_bin_hi, nbins);
   D_REALLOC_OR_FAIL(h->orig_counts, nbins);
   D_REALLOC_OR_FAIL(h->orig_stat_err, nbins);
   I_REALLOC_OR_FAIL(h->orig_notice, nbins);
   I_REALLOC_OR_FAIL(h->quality, nbins);

   memset ((char *)h->stat_err, 0, nbins*sizeof(double));
   memset ((char *)h->quality, 0, nbins*sizeof(int));
   memset ((char *)h->notice, 0, nbins*sizeof(int));
   memset ((char *)h->model_flux.val, 0, nbins*sizeof(double));
   memset ((char *)h->model_counts, 0, nbins*sizeof(double));
   memset ((char *)h->convolved_model_flux, 0, nbins*sizeof(double));

   for (i = 0; i < nbins; i++)
     {
        h->flux[i] = 0.0;
        h->flux_err[i] = -1.0;
     }

   if (NULL == (cts = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     goto error_return;

   if (-1 == rebin_histogram (h->counts, h->bin_lo, h->bin_hi, h->nbins,
                              cts, lo, hi, nbins))
     goto error_return;

   if (h->orig_bgd != NULL)
     {
        if (NULL == (bgd_cts = (double *) ISIS_MALLOC (nbins * sizeof(double))))
          goto error_return;

        if (-1 == rebin_histogram (h->orig_bgd, h->bin_lo, h->bin_hi, h->nbins,
                                   bgd_cts, lo, hi, nbins))
          goto error_return;
     }

   if ((-1 == rebin_backscale (&h->area, h, lo, hi, nbins))
       || (-1 == rebin_backscale (&h->bgd_area, h, lo, hi, nbins)))
     goto error_return;

   if ((-1 == copy_d_array (&h->bin_lo, lo, nbins))
       || (-1 == copy_d_array (&h->bin_hi, hi, nbins)))
     goto error_return;

   ISIS_FREE (h->counts);
   h->counts = cts;

   ISIS_FREE (h->orig_bgd);
   h->orig_bgd = bgd_cts;

   ISIS_FREE (h->sys_err_frac);
   h->sys_err_frac = sys_err_frac;

   h->nbins = nbins;
   h->orig_nbins = nbins;
   h->n_notice = nbins;

   if (validate_stat_err (h) < 0)
     goto error_return;

   if (0 == Rmf_is_identity (h->a_rsp.rmf))
     {
        if (-1 == rebin_rmf (h, lo, hi, nbins))
          goto error_return;
     }
   else
     {
    /* finish_hist_init() will assign a new delta-RMF */
        release_rmf (h->a_rsp.rmf);
        h->a_rsp.rmf = NULL;
    /* this will get re-set by finish_hist_init() */
        release_rmf (h->f_rsp.rmf);
        release_arf (h->f_rsp.arf);
        h->f_rsp.rmf = NULL;
        h->f_rsp.arf = NULL;
     }

   return finish_hist_init (h);

error_return:
   if (h->counts == cts || h->orig_bgd == bgd_cts)
     {
        isis_vmesg (FAIL, I_INTERNAL, __FILE__, __LINE__,
                    "*** rebinned dataset was left in an inconsistent state");
     }
   else
     {
        ISIS_FREE(cts);
        ISIS_FREE(bgd_cts);
     }
   return -1;
}

/*}}}*/

int Hist_set_model (Hist_t *h, unsigned int version, double *model_value) /*{{{*/
{
   Isis_Hist_t hist;

   if (is_data(version))
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__, NULL);
        return -1;
     }

   if (h == NULL)
     return -1;

   if (-1 == get_hist_version (h, version, &hist))
     return -1;

   memset ((char *)hist.val, 0, hist.nbins * sizeof(double));

   return unpack_noticed (model_value, hist.notice_list,
                          hist.n_notice, hist.nbins, hist.val);
}

/*}}}*/

/* Feature measurements */

int Hist_print_stats (FILE * fp, Hist_Stat_t *s) /*{{{*/
{
   static const char *header =
     "#   w_min           w_max       nbins     centroid     ctr_err      sum       sum_err      net       net_err      contin    con_slope      eq_width   eqw_err";

   if (NULL == fp)
     return -1;

   if (s == NULL)
     {
        fprintf (fp, "%s\n", header);
        fflush(fp);
        return 0;
     }

   fprintf (fp, "%15.8e %15.8e %5d %15.8e %9.2e %12.5e %9.2e %12.5e %9.2e %12.5e %12.5e %12.5e %9.2e\n",
            s->min, s->max,
            s->nbins,
            s->centroid, s->centroid_err,
            s->sum, s->sum_err,
            s->net, s->net_err,
            s->contin, s->slope,
            s->eqwidth, s->eqwidth_err
            );
   fflush(fp);
   return 0;
}

/*}}}*/

int Hist_region_stats (Hist_t *h, unsigned int version, Hist_Stat_t *s, double x[/*2*/], double y[/*2*/]) /*{{{*/
{
   Isis_Hist_t hist;
   int i, imin, imax;
   double a, b, lo, hi, val, err2, contin, sxy, sxy_err;
   int eqwidth = 1;

   /* x and y should contain the coordinates of 2 points which mark
    * the region boundaries at the intended continuum DENSITY level --- the
    * continuum is assumed to be a straight line connecting the two points.
    */
   if (s == NULL || h == NULL || x[0] >= x[1])
     return -1;

   memset ((char *) s, 0, sizeof(*s));

   if (-1 == get_hist_version (h, version, &hist))
     return -1;

   imin = find_bin (x[0], hist.bin_lo, hist.bin_hi, hist.nbins);
   if (imin < 0)
     imin = 0;

   imax = find_bin (x[1], hist.bin_lo, hist.bin_hi, hist.nbins);
   if (imax < 0)
     imax = hist.nbins - 1;

   /* linear continuum DENSITY:  y = a * x + b */

   a = (y[1] - y[0]) / (x[1] - x[0]);
   b = y[0] - a * x[0];

   s->min = x[0];
   s->max = x[1];

   for (i=imin; i <= imax; i++)
     {
        double width, mid;

        if (hist.notice[i])
          {
             lo = hist.bin_lo[i];
             hi = hist.bin_hi[i];
             width = hi - lo;
             mid   = 0.5 * (hi + lo);

             /* stored histogram values are always some bin-integrated
              * quantity vs. wavelength in angstroms */

             val  = hist.val[i];
             if (hist.val_err != NULL)
               err2 = hist.val_err[i] * hist.val_err[i];
             else
               err2 = 0.0;

             /* bin integrated continuum */
             contin = width * (a * mid + b);

             s->sum += val;
             s->net += val - contin;
             s->centroid += mid * fabs(val - contin);

             /* rms error */
             s->sum_err += err2;
             s->net_err += err2;
             s->centroid_err += err2 * mid * mid;

             /* negative of usual equivalent width definition,
              * with rms error */
             if (contin > 0.0)
               {
                  double contin_dens = contin / width;
                  s->eqwidth += width * fabs(1.0 - val / contin);
                  s->eqwidth_err += err2 / contin_dens / contin_dens;
               }
             else
               eqwidth = 0;

             s->nbins++;
          }
     }

   if (s->nbins <= 0 || s->net == 0.0)
     return 0;

   s->sum_err = sqrt (s->sum_err / s->nbins);
   s->net_err = sqrt (s->net_err / s->nbins);

   sxy = s->centroid;
   s->centroid = sxy / fabs(s->net);

   sxy_err = sqrt (s->centroid_err / s->nbins);
   if (sxy != 0.0)
     s->centroid_err = s->centroid * fabs (sxy_err / sxy + s->net_err / s->net);
   else
     s->centroid_err = s->centroid * (s->net_err / s->net);

   if (eqwidth)
     s->eqwidth_err = sqrt (s->eqwidth_err / s->nbins);
   else
     {
        s->eqwidth = 0.0;
        s->eqwidth_err = 0.0;
     }

   /* continuum density and slope */
   s->contin = a * s->centroid + b;
   s->slope  = a;

   return 0;
}

/*}}}*/

/* set notice flags, get/put using notice flags */

unsigned int Hist_num_noticed_histograms (Hist_t *head) /*{{{*/
{
   Hist_t *h;
   unsigned int num;

   if (head == NULL)
     return 0;

   num = 0;
   for (h=head->next; h != NULL; h = h->next)
     {
        if ((h->exclude == 0) && (h->n_notice > 0))
          num++;
     }

   return num;
}

/*}}}*/

int Hist_id_list (Hist_t *head, int noticed, unsigned int **ids, unsigned int *num) /*{{{*/
{
   Hist_t *h;
   unsigned int n;

   *ids = NULL;
   *num = 0;

   if (head == NULL)
     return 0;

   for (h = head->next; h != NULL; h = h->next)
     {
        if ((noticed != 0) && ((h->n_notice == 0) || (h->exclude != 0)))
          continue;
        *num += 1;
     }

   if (*num == 0)
     return 0;

   if (NULL == (*ids = (unsigned int *) ISIS_MALLOC (*num * sizeof(unsigned int))))
     return -1;

   n = 0;
   for (h = head->next; h != NULL; h = h->next)
     {
        if ((noticed != 0) && ((h->n_notice == 0) || (h->exclude != 0)))
          continue;
        (*ids)[n++] = h->index;
     }

   return 0;
}

/*}}}*/

int Hist_num_data_noticed (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return -1;

   if (h->exclude != 0)
     return 0;

   return h->n_notice;
}

/*}}}*/

int Hist_all_data_noticed (Hist_t *head) /*{{{*/
{
   Hist_t *h;
   int n = 0;

   if (head == NULL)
     return -1;

   for (h=head->next; h != NULL; h = h->next)
     {
        if (h->exclude == 0)
          n += h->n_notice;
     }

   return n;
}

/*}}}*/

int Hist_all_model_noticed (Hist_t *head) /*{{{*/
{
   Hist_t *h;
   int num = 0;

   if (head == NULL)
     return -1;

   for (h = head->next; h != NULL; h = h->next)
     {
        if (h->exclude == 0)
          num += h->model_flux.n_notice;
     }

   return num;
}

/*}}}*/

static int update_notice_list (Hist_t *h) /*{{{*/
{
   int n_notice1;

   if (h == NULL)
     return -1;

   n_notice1 = h->n_notice;

   if (-1 == _update_notice_list (h->notice, &h->notice_list, &h->n_notice, h->nbins))
     return -1;

   /* Changing the number of noticed datasets might change
    * the fit-function (e.g. when components of the fit-function
    * apply only to specific datasets and those specific datasets
    * go from being ignored to being noticed, or vice-versa.)
    */
   if (((h->n_notice > 0) && (n_notice1 == 0))
       || ((h->n_notice == 0) && (n_notice1 > 0)))
     {
        update_user_model();
     }

   /* nbins == orig_nbins implies that the data
    * are currently not re-binned, so we should record this
    * ignore/notice state -- NOTE that h->orig_notice ONLY gets
    * used when rebinning the data.
    */

   if (h->nbins == h->orig_nbins)
     memcpy ((char *)h->orig_notice, (char *)h->notice, h->orig_nbins * sizeof(int));

   return 0;
}

/*}}}*/

static int apply_notice_value (int value, double bin_lo, double bin_hi, /*{{{*/
                               int *notice, double *grid_lo, double *grid_hi, int nbins)
{
   int i;

   if (NULL == notice || NULL == grid_lo || NULL == grid_hi)
     return -1;

   value = (value == 0) ? 0 : 1;

   for (i = 0; i < nbins; i++)
     {
        if (INTERVALS_OVERLAP(bin_lo, bin_hi, grid_lo[i], grid_hi[i]))
          notice[i] = value;
     }

   return 0;
}

/*}}}*/

int Hist_ignore_bad (Hist_t *h) /*{{{*/
{
   int i, n;

   if (NULL == h)
     return -1;

   n = h->nbins;

   if (n != h->orig_nbins)
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__, "data must be unbinned");
        return -1;
     }

   for (i = 0; i < n; i++)
     {
        if (h->quality[i] != 0)
          h->notice[i] = 0;
     }

   return update_notice_list (h);
}

/*}}}*/

int Hist_set_notice_using_mask (Hist_t *h, int *mask, int nbins) /*{{{*/
{
   int i;

   if ((h == NULL) || (mask == NULL) || (nbins != h->nbins))
     return -1;

   for (i = 0; i < nbins; i++)
     {
        h->notice[i] = mask[i];
     }

   return update_notice_list (h);
}

/*}}}*/

int Hist_set_notice_using_list (Hist_t *h, int value, unsigned int *list, unsigned int n) /*{{{*/
{
   unsigned int i;

   if ((h == NULL) || (list == NULL) || ((int) n > h->nbins))
     return -1;

   for (i = 0; i < n; i++)
     {
        h->notice[list[i]] = value;
     }

   return update_notice_list (h);
}

/*}}}*/

int Hist_set_notice (Hist_t *h, int value, double bin_lo, double bin_hi) /*{{{*/
{
   if (NULL == h)
     return -1;

   if ((h->has_grid == -1)
       && ((h->bin_lo[0] < bin_lo) || (bin_hi < h->bin_hi[h->nbins-1])))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Data set %d has no grid", h->index);
        return -1;
     }

   if (bin_lo > bin_hi)
     {
        double tmp = bin_lo;
        bin_lo = bin_hi;
        bin_hi = tmp;
     }

   if (-1 == apply_notice_value (value, bin_lo, bin_hi, h->notice, h->bin_lo, h->bin_hi, h->nbins))
     return -1;

   return update_notice_list (h);
}

/*}}}*/

#ifdef TEST_MODEL_NOTICE
static void print_notice_match (Hist_t *h, int *notice) /*{{{*/
{
   FILE *fp = NULL;
   static char file [] = "/tmp/isis_tmp.out";
   Isis_Hist_t *m = &h->model_flux;
   int i;

   fp = fopen (file, "w");
   if (fp == NULL)
     {
        fprintf (stderr, "*** Failed opening test output file %s\n", file);
        return;
     }

   fprintf (fp, "###############  Data set %d\n", h->index);
   for (i = 0; i < h->orig_nbins; i++)
     {
        if (notice[i] != 0)
          {
             int j, jstart, jend;

             fprintf (fp, "data bin %d:  %g - %g Angstrom\n",
                      i, h->orig_bin_lo[i], h->orig_bin_hi[i]);

             jstart = jend = -1;
             for (j = 0; j < m->nbins; j++)
               {
                  if (m->notice[j] != 0)
                    {
                       if (jstart < 0) jstart = j;
                       jend = j;
                    }
               }

             if (jstart > 0 && jend > 0)
               {
                  fprintf (fp, "model bins %d-%d:  %g - %g Angstroms\n",
                           jstart, jend,
                           m->bin_lo[jstart], m->bin_hi[jend]);
               }
             else fprintf (fp, "--> no model bins noticed !!");

          }
     }

   fclose (fp);
}
/*}}}*/
#endif

static int set_model_notice (Hist_t *h) /*{{{*/
{
   Isis_Hist_t *m = &h->model_flux;
   Isis_Rsp_t *r;
   int have_multi_rsp;
   int *notice = NULL;
   int ret = -1;

   /* Apply current notice flags to original bins,
    * then map onto model grid using the RMF
    */
   if (NULL == (notice = (int *) ISIS_MALLOC (h->orig_nbins * sizeof(int))))
     return -1;

   memset ((char *)notice, 0, h->orig_nbins * sizeof(int));
   memset ((char *)m->notice, 0, m->nbins * sizeof(int));

   if (-1 == transfer_notice (h->bin_lo, h->bin_hi, h->notice_list, h->n_notice,
                              h->orig_bin_lo, h->orig_bin_hi, h->orig_nbins, notice))
     goto finish;

   have_multi_rsp = (h->f_rsp.next != NULL) ? 1 : 0;

   for (r = &h->f_rsp; r != NULL; r = r->next)
     {
        Isis_Rmf_t *rmf = r->rmf;
        Isis_Hist_t x, *px = NULL;
        int malloced_px = 0;

        /* easy case:  The model is on the ARF grid */
        if (have_multi_rsp == 0)
          px = m;
        else
          {
             /* harder case:
              *   The model grid is different from the ARF grid,
              *   so we'll have to take an extra step. */
             Isis_Arf_t *arf = r->arf;
             px = &x;
             px->bin_lo = arf->bin_lo;
             px->bin_hi = arf->bin_hi;
             px->nbins = arf->nbins;
             px->notice_list = NULL;
             px->n_notice = 0;
             if (NULL == (px->notice = (int *) ISIS_MALLOC (arf->nbins * sizeof(int))))
               goto finish;
             malloced_px = 1;
          }

        if (-1 == rmf->set_noticed_model_bins (rmf, h->orig_nbins, notice, px->nbins, px->notice))
          {
             if (malloced_px)
               ISIS_FREE (px->notice);
             goto finish;
          }

        if (have_multi_rsp != 0)
          {
             /* Use the ARF-grid notice flags to update the model-grid notice flags */
             if ((-1 == _update_notice_list (px->notice, &px->notice_list, &px->n_notice, px->nbins))
                 || (-1 == transfer_notice (px->bin_lo, px->bin_hi, px->notice_list, px->n_notice,
                                            m->bin_lo, m->bin_hi, m->nbins, m->notice)))
               {
                  ISIS_FREE (px->notice);
                  goto finish;
               }
             ISIS_FREE (px->notice_list);
             ISIS_FREE (px->notice);
          }
     }

   ret = 0;
   finish:

   ISIS_FREE (notice);

   return ret;
}

/*}}}*/

static int update_model_notice_list (Hist_t *h) /*{{{*/
{
   Isis_Kernel_Def_t *def = h->kernel->kernel_def;
   Isis_Hist_t *m = &h->model_flux;
   int allow_ignore;

   if (def->allows_ignoring_model_intervals)
     allow_ignore = (*def->allows_ignoring_model_intervals)(h->kernel);
   else allow_ignore = 0;

   if (allow_ignore)
     {
        if (-1 == set_model_notice (h))
          return -1;
     }
   else
     {
        int i;
        for (i = 0; i < m->nbins; i++)
          m->notice[i] = 1;
     }

   if (-1 == _update_notice_list (m->notice, &m->notice_list, &m->n_notice, m->nbins))
     return -1;

   return update_notice_list (h);
}

/*}}}*/

int Hist_init_model_structs (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return -1;

   if (-1 == set_model_grid (h))
     return -1;

   return update_model_notice_list (h);
}

/*}}}*/

int Hist_set_exclude_flag (Hist_t *h, int exclude) /*{{{*/
{
   if (h == NULL)
     return -1;

   h->exclude = exclude ? 1 : 0;

   /* Changing the number of noticed datasets might change
    * the fit-function (e.g. when components of the fit-function
    * apply only to specific datasets and those specific datasets
    * go from being ignored to being noticed, or vice-versa.)
    */
   update_user_model ();

   return 0;
}

/*}}}*/

/* combining datasets */

int Hist_any_noticed_dataset_combinations (Hist_t *head) /*{{{*/
{
   Hist_t *h;

   if (head == NULL)
     return 0;

   /* more subtle definitions are possible */

   for (h = head->next; h != NULL; h = h->next)
     {
        if ((h->exclude == 0)
            && (h->n_notice > 0)
            && (h->combo_id != 0))
          return 1;
     }

   return 0;
}

/*}}}*/

int Hist_break_combination (Hist_t *head, unsigned int gid) /*{{{*/
{
   Hist_t *h;

   if (head == NULL)
     return 0;

   /* gid = 0 means break up all dataset combinations */

   for (h = head->next; h != NULL; h = h->next)
     {
        if ((gid == 0) || (h->combo_id == gid))
          {
             h->combo_id = 0;
             init_eval_grid_method (&h->eval_grid_method);
          }
     }

   if (gid == 0)
     Next_Combo_Id = 1;

   return 0;
}

/*}}}*/

int Hist_combination_id (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return 0;

   return h->combo_id;
}

/*}}}*/

int Hist_combination_weight (Hist_t *h, double *weight) /*{{{*/
{
   if (h == NULL || weight == NULL)
     return -1;

   *weight = h->combo_weight;
   return 0;
}

/*}}}*/

int Hist_dataset_is_combined (Hist_t *h, unsigned int gid) /*{{{*/
{
   if (h == NULL)
     return 0;

   return (h->combo_id == gid);
}

/*}}}*/

int Hist_grids_are_identical (Hist_t *a, Hist_t *b) /*{{{*/
{
   int i;

   if (a == NULL || b == NULL)
     return 0;

   if (a == b)
     return 1;

   /* check the grid sizes and bin edges */

   if (a->nbins != b->nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grid size mismatch, datasets %d and %d",
                    a->index, b->index);
        return 0;
     }

   for (i = 0; i < a->nbins; i++)
     {
        if ((a->bin_lo[i] != b->bin_lo[i])
            || (a->bin_hi[i] != b->bin_hi[i]))
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grid mismatch, datasets %d and %d",
                         a->index, b->index);
             return 0;
          }
     }

   /* check that the same bins are noticed */

   if (a->n_notice != b->n_notice)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "different ranges noticed, datasets %d and %d",
                    a->index, b->index);
        return 0;
     }

   for (i = 0; i < a->n_notice; i++)
     {
        if (a->notice_list[i] != b->notice_list[i])
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "different ranges noticed, datasets %d and %d",
                         a->index, b->index);
             return 0;
          }
     }

   return 1;
}

/*}}}*/

static int all_member_grids_are_identical (Hist_t *head, unsigned int *members, /*{{{*/
                                           unsigned int num_members)
{
   Hist_t *first;
   unsigned int k;

   if ((members == NULL) || (num_members < 2))
     return 0;

   if (NULL == (first = Hist_find_hist_index (head, members[0])))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "dataset %d does not exist",
                   members[0]);
        return 0;
     }

   for (k = 1; k < num_members; k++)
     {
        Hist_t *h;

        if (NULL == (h = Hist_find_hist_index (head, members[k])))
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "dataset %d does not exist",
                         members[k]);
             return 0;
          }

        if (0 == Hist_grids_are_identical (first, h))
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grid mismatch; run match_dataset_grids() first.");
             return 0;
          }
     }

   return 1;
}

/*}}}*/

int Hist_set_combination (Hist_t *head, unsigned int *members, /*{{{*/
                          double *weights, unsigned int num_members)
{
   Hist_t *h;
   unsigned int k;
   int gid;

   if ((head == NULL) || (weights == NULL) || (num_members < 1))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "dataset combination not defined");
        return 0;
     }

   if (num_members > 1)
     {
        if (0 == all_member_grids_are_identical (head, members, num_members))
          return -1;
     }

   gid = Next_Combo_Id;

   for (k = 0; k < num_members; k++)
     {
        if (NULL == (h = Hist_find_hist_index (head, members[k])))
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Dataset %d not found -- dataset combination not defined",
                        members[k]);
             return -1;
          }
        h->combo_id = gid;
        h->combo_weight = weights[k];
     }

   Next_Combo_Id++;

   return gid;
}

/*}}}*/

int Hist_set_pre_combine_hook (Hist_t *h, SLang_Name_Type *hook) /*{{{*/
{
   if (h == NULL)
     return -1;
   SLang_free_function (h->pre_combine);
   h->pre_combine = hook;
   return 0;
}

/*}}}*/

int Hist_call_pre_combine_hook (Hist_t *h, double *y, SLindex_Type n, double **y_new) /*{{{*/
{
   SLang_Array_Type *sl_y = NULL, *sl_new_y = NULL;
   double *new_y = NULL;
   int status = -1;

   *y_new = NULL;

   if (h->pre_combine == NULL)
     {
        *y_new = y;
        return 0;
     }

   /* When a pre-combine hook is used, the eval_grid_method should
    * probably be ISIS_EVAL_GRID_SEPARATE, but I don't want to rule
    * out the possibility that ISIS_EVAL_GRID_USER might work too.
    * ISIS_EVAL_GRID_MERGED definitely makes no sense, so I'll
    * warn about that.
    */
   if (h->eval_grid_method.type == ISIS_EVAL_GRID_MERGED)
     {
        isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__,
                    "dataset %d: a pre-combine hook is assigned, but eval_grid_method=MERGED_GRID",
                    h->index);
        return -1;
     }

   if (NULL == (sl_y = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1)))
     return -1;

   memcpy ((char *)sl_y->data, (char *)y, n * sizeof(double));

   SLang_start_arg_list ();
   SLang_push_integer (h->index);
   SLang_push_array (sl_y, 1);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (h->pre_combine))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "executing pre_combine hook (dataset id = %d)", h->index);
        goto return_error;
     }

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLdo_pop ();
        goto return_error;
     }

   if ((-1 == SLang_pop_array_of_type (&sl_new_y, SLANG_DOUBLE_TYPE))
       || (sl_new_y == NULL))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "in pre_combine hook (dataset id = %d)", h->index);
        goto return_error;
     }

   if (sl_new_y->num_elements != (unsigned int) n)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "pre_combine hook (dataset id = %d) returned wrong length array (was %d should be %d)",
                    h->index, sl_new_y->num_elements, n);
        goto return_error;
     }

   if (NULL == (new_y = (double *) ISIS_MALLOC (n * sizeof(*new_y))))
     goto return_error;

   memcpy ((char *)new_y, (char *)sl_new_y->data, n * sizeof(double));
   *y_new = new_y;
   status = 0;

return_error:
   SLang_free_array (sl_new_y);
   if (status) ISIS_FREE (new_y);
   return status;
}

/*}}}*/

/* make histogram plots */

int Hist_unset_color (Hist_t *h) /*{{{*/
{
   if (NULL == h)
     return -1;

   ISIS_FREE (h->color);

   return 0;
}

/*}}}*/

int Hist_set_color (Hist_t *h, int color) /*{{{*/
{
   if (NULL == h)
     return -1;

   /* I know allocating a single int looks weird, but
    * I'm using a NULL ptr to indicate that no color has been chosen -
    * that way any valid integer can map to a color.
    * too extreme?
    */

   if (h->color == NULL)
     {
        if (NULL == (h->color = (int *) ISIS_MALLOC (sizeof(int))))
          return -1;
     }

   *h->color = color;

   return 0;
}

/*}}}*/

static int set_default_labels (unsigned int version, Hist_Plot_Tune_Type *info, Plot_t *fmt) /*{{{*/
{
   char xlabel[PLOT_LABEL_STRING_SIZE], ylabel[PLOT_LABEL_STRING_SIZE];
   char xunit[PLOT_LABEL_STRING_SIZE], yunit[PLOT_LABEL_STRING_SIZE];
   const char *rtype, *ytype, *xlog, *ylog;
   int x_unit_id = Plot_x_unit (fmt);
   int x_unit_is_wavelength = unit_is_wavelength (x_unit_id);

   xlabel[0] = 0;
   ylabel[0] = 0;
   yunit[0] = 0;

   if (-1 == unit_name (xunit, x_unit_id))
     return -1;

   /* see plot.c cpgenv/cpgbox usage */
#ifndef PLOT_NUMBERS_ARE_LOGARITHMS
# define PLOT_NUMBERS_ARE_LOGARITHMS 0
#endif
#if PLOT_NUMBERS_ARE_LOGARITHMS
   if (Plot_get_logx (fmt))
     xlog = "log ";
   else
     xlog = "";

   if (Plot_get_logy (fmt))
     ylog = "log ";
   else
     ylog = "";
#else
   xlog = "";
   ylog = "";
#endif

   if (x_unit_is_wavelength)
     {
        if (-1 == isis_strcat (xlabel, sizeof(xlabel), xlog, "Wavelength [", xunit, "]", NULL))
          return -1;
     }
   else
     {
        if (-1 == isis_strcat (xlabel, sizeof(xlabel), xlog, "Energy [", xunit, "]", NULL))
          return -1;
     }

   if (is_flux(version))
     ytype = "Photons/s/cm\\u2\\d";
   else
     ytype = "Counts";

   if ((info->rate != 0) && is_counts(version))
     rtype = "/s";
   else
     rtype = "";

   if (Plot_bin_density(fmt))
     {
        if (-1 == isis_strcat (yunit, sizeof(yunit), ytype, rtype, "/", xunit, NULL))
          return -1;
     }
   else
     {
        if (-1 == isis_strcat (yunit, sizeof(yunit), ytype, rtype, "/bin", NULL))
          return -1;
     }

   if (is_flux(version))
     {
        if (-1 == isis_strcat (ylabel, sizeof(ylabel), ylog, "Flux [", yunit, "]", NULL))
          return -1;
     }
   else
     {
        if (-1 == isis_strcat (ylabel, sizeof(ylabel), ylog, yunit, NULL))
          return -1;
     }

   Plot_set_label_type (fmt, PLOT_DEFAULT_LABEL);
   return Plot_set_labels (fmt, xlabel, ylabel, NULL);
}

/*}}}*/

static void free_plot_copy (Isis_Hist_t *cpy) /*{{{*/
{
   ISIS_FREE (cpy->val);
   ISIS_FREE (cpy->val_err);
   ISIS_FREE (cpy->notice);
   ISIS_FREE (cpy->bin_lo);
   ISIS_FREE (cpy->bin_hi);
   ISIS_FREE (cpy->sys_err_frac);
}

/*}}}*/

static int make_plot_copy (Isis_Hist_t *src, Isis_Hist_t *dest) /*{{{*/
{
   unsigned int d_size, i_size;

   if ((NULL == src) || (NULL == dest))
     return -1;

   if (NULL == (dest->val = (double *) ISIS_MALLOC (src->nbins * sizeof(double)))
       || NULL == (dest->bin_lo = (double *) ISIS_MALLOC (src->nbins * sizeof(double)))
       || NULL == (dest->bin_hi = (double *) ISIS_MALLOC (src->nbins * sizeof(double)))
       || NULL == (dest->val_err = (double *) ISIS_MALLOC (src->nbins * sizeof(double)))
       || ((src->sys_err_frac != NULL) && (NULL == (dest->sys_err_frac = (double *) ISIS_MALLOC (src->nbins * sizeof(double)))))
       || ((src->notice != NULL) && (NULL == (dest->notice = (int *) ISIS_MALLOC (src->nbins * sizeof(int)))))
      )
     {
        free_plot_copy (dest);
        return -1;
     }

   d_size = src->nbins * sizeof(double);
   i_size = src->nbins * sizeof(int);

   memcpy ((char *) dest->val,    (char *) src->val,    d_size);
   memcpy ((char *) dest->bin_lo, (char *) src->bin_lo, d_size);
   memcpy ((char *) dest->bin_hi, (char *) src->bin_hi, d_size);

   if (src->val_err != NULL)
     memcpy ((char *) dest->val_err, (char *) src->val_err, d_size);
   else
     memset ((char *) dest->val_err, 0, d_size);

   if (src->sys_err_frac != NULL)
     memcpy ((char *) dest->sys_err_frac, (char *) src->sys_err_frac, d_size);

   if (src->notice != NULL)
     memcpy ((char *) dest->notice, (char *) src->notice, i_size);

   dest->nbins = src->nbins;

   return 0;
}
/*}}}*/

static int convert_from_angstrom (int x_unit, Isis_Hist_t *g) /*{{{*/
{
   if (g == NULL)
     return -1;

   if (0 == unit_is_wavelength (x_unit))
     {
        double *tmp = g->bin_lo;
        g->bin_lo = g->bin_hi;
        g->bin_hi = tmp;

        if (-1 == reverse_d (g->val, g->nbins)
            || -1 == reverse_d (g->bin_lo, g->nbins)
            || -1 == reverse_d (g->bin_hi, g->nbins)
            || ((NULL != g->val_err) && (-1 == reverse_d (g->val_err, g->nbins)))
            || ((NULL != g->notice) && (-1 == reverse_i (g->notice, g->nbins))))
          {
             return -1;
          }
     }

   if (-1 == unit_convert_x (g->bin_lo, g->nbins, x_unit, U_ANGSTROM)
       || -1 == unit_convert_x (g->bin_hi, g->nbins, x_unit, U_ANGSTROM))
     {
        return -1;
     }

   return 0;
}

/*}}}*/

static double get_assigned_exposure (Hist_t *h) /*{{{*/
{
   double t = 1.0;

   if ((0 == Arf_is_identity (h->a_rsp.arf))
       && (h->a_rsp.arf->exposure > 0.0))
     t = h->a_rsp.arf->exposure;
   else if (h->exposure > 0.0)
     t = h->exposure;

   return t;
}

/*}}}*/

static int compute_rate (double t, Isis_Hist_t *cpy) /*{{{*/
{
   double *val;
   int i, n;

   val = cpy->val;
   n  = cpy->nbins;

   for (i = 0; i < n; i++)
     val[i] /= t;

   if (cpy->val_err != NULL)
     {
        double *val_err = cpy->val_err;
        for (i = 0; i < n; i++)
          val_err[i] /= t;
     }

   return 0;
}

/*}}}*/

static void free_opt_data (Isis_Fit_Statistic_Optional_Data_Type *opt_data) /*{{{*/
{
   if (opt_data == NULL)
     return;
   ISIS_FREE(opt_data->bkg);
   ISIS_FREE(opt_data->bkg_at);
   ISIS_FREE(opt_data->src_at);
}

/*}}}*/

static int init_opt_data (Hist_t *h, unsigned int version, /*{{{*/
                          Isis_Fit_Statistic_Optional_Data_Type *opt_data)
{
   int nb, n = h->nbins;

   opt_data->num = h->nbins;
   opt_data->bkg = NULL;

   if (is_counts(version))
     {
        if (h->instrumental_background_hook != NULL)
          {
             /* FIXME? */ fprintf (stderr, "*** init_opt_data:  sorry, this "
                                   "functionality has not been implemented.\n");
             return -1;
          }
        if (-1 == copy_input_background (h, 0, &opt_data->bkg, &nb))
          return -1;
        /* opt_data->bkg may be NULL here */
     }

   if (opt_data->bkg == NULL)
     {
        if (NULL == (opt_data->bkg = (double *)ISIS_MALLOC (n*sizeof(double))))
          return -1;
        memset ((char *)opt_data->bkg, 0, n*sizeof(double));
     }

   if (-1 == Hist_scaling_vectors (h, 0, 0, &opt_data->src_at, &opt_data->bkg_at, &nb))
     {
        free_opt_data (opt_data);
        return -1;
     }

   return 0;
}

/*}}}*/

static int compute_plot_residuals (Hist_t *h, unsigned int version, /*{{{*/
                                   Hist_Plot_Tune_Type *info,
                                   Isis_Fit_Statistic_Type *s,
                                   Isis_Hist_t *cpy, char **label)
{
   Isis_Hist_t *m, *d;
   Isis_Hist_t other;
   unsigned int i, n;
   double t, statistic;
   double *wt, *total_err;
   int status = -1;

   if (is_model(version))
     {
        bit_set(version, H_DATA);
        if (-1 == get_hist_version (h, version, &other))
          return -1;

        m = cpy;
        d = &other;
     }
   else
     {
        bit_clear(version, H_DATA);
        if (-1 == get_hist_version (h, version, &other))
          return -1;

        m = &other;
        d = cpy;
     }

   n = cpy->nbins;

   if (d->sys_err_frac == NULL)
     total_err = d->val_err;
   else
     {
        if (NULL == (total_err = (double *)ISIS_MALLOC (n * sizeof(double))))
          return -1;
        for (i = 0; i < n; i++)
          {
             total_err[i] = isis_hypot (d->val_err[i],
                                        d->sys_err_frac[i] * d->val[i]);
          }
     }

   switch (info->residuals)
     {
      default:
        if (s == NULL)
          goto return_error;
        *label = s->symbol;
        /* use cpy->val_err as temporary work space for weights */
        wt = cpy->val_err;
        for (i = 0; i < n; i++)
          {
             /* val_err != 0 guaranteed elsewhere */
             double dv = total_err[i];
             wt[i] = 1.0/(dv*dv);
          }
          {
             Isis_Fit_Statistic_Optional_Data_Type opt_data;
             int stat_status;

             if (s->uses_opt_data == 0)
               s->opt_data = NULL;
             else
               {
                  memset ((char *)&opt_data, 0, sizeof opt_data);
                  if (-1 == init_opt_data (h, version, &opt_data))
                    goto return_error;
                  s->opt_data = &opt_data;
               }

             stat_status = s->compute_statistic (s, d->val, m->val, wt, n, cpy->val, &statistic);
             s->opt_data = NULL;

             if (s->uses_opt_data)
               free_opt_data (&opt_data);

             if (stat_status)
               goto return_error;
          }
        for (i = 0; i < n; i++)
          cpy->val_err[i] = 1.0;
        break;

      case ISIS_DIFF_RESID:
        *label = "D-M";

        if (info->rate)
          t = get_assigned_exposure (h);
        else t = 1.0;

        for (i = 0; i < n; i++)
          {
             cpy->val[i] = (d->val[i] - m->val[i]) / t;
             cpy->val_err[i] = total_err[i] / t;
          }
        break;

      case ISIS_RATIO_RESID:
        *label = "D/M";
        for (i = 0; i < n; i++)
          {
             if (m->val[i] == 0.0)
               {
                  cpy->val[i] = FLT_MAX;
                  cpy->val_err[i] = FLT_MAX;
               }
             else
               {
                  cpy->val[i] = d->val[i] / m->val[i];
                  cpy->val_err[i] = total_err[i] / m->val[i];
               }
          }
        break;
     }

   status = 0;
return_error:

   if (total_err != d->val_err)
     ISIS_FREE(total_err);

   return status;
}

/*}}}*/

static void plot_reference_line (double ref) /*{{{*/
{
   float x[2], y[2];
   float xmin, xmax, ymin, ymax;

   Plot_query_plot_limits (&xmin, &xmax, &ymin, &ymax);
   x[0] = xmin;
   x[1] = xmax;
   y[0] = ref;
   y[1] = ref;

   _Plot_set_color (PLOT_LINE_COLOR);
   Plot_line (2, x, y);
}

/*}}}*/

int Hist_plot (Hist_t *h, unsigned int version, Hist_Plot_Tune_Type *info, /*{{{*/
               Plot_t *fmt, Isis_Fit_Statistic_Type *s)
{
   int x_unit, eb_save, ylog_save, bd_save;
   Isis_Hist_t hist, cpy;
   int ret = -1;
   char *label = NULL;

   if (NULL == fmt)
     return -1;

   if ((NULL == h)
       || (-1 == get_hist_version (h, version, &hist)))
     return -1;

   if (hist.nbins == 0)
     return -1;

   /* If the user has not explicitly specified the plot units,
    * choose the default units based on the dispersion order.
    * Grating data defaults to wavelength units,
    * CCD data defaults to keV
    */
   if (Plot_x_unit_is_default (fmt))
     Plot_set_x_unit (fmt, h->order ? U_ANGSTROM : U_KEV);

   x_unit = Plot_x_unit (fmt);

   memset ((char *)&cpy, 0, sizeof(cpy));
   if (-1 == make_plot_copy (&hist, &cpy))
     return -1;

   /* The order of these transformations matters! */
   if (info->residuals != 0)
     {
        if (-1 == compute_plot_residuals (h, version, info, s, &cpy, &label))
          goto finish;
     }
   else if (info->rate != 0)
     {
        if (!is_flux(version))
          {
             double t = get_assigned_exposure (h);
             if (-1 == Plot_save_exposure_time (fmt, t)) /* ugly hack */
               goto finish;
             if (-1 == compute_rate (t, &cpy))
               goto finish;
          }
     }

   if ((cpy.sys_err_frac != NULL) && (info->residuals == 0))
     {
        int i;
        for (i = 0; i < cpy.nbins; i++)
          {
             cpy.val_err[i] = isis_hypot (cpy.val_err[i],
                                          cpy.sys_err_frac[i] * cpy.val[i]);
          }
     }

   if (x_unit != U_ANGSTROM)
     {
        if (-1 == convert_from_angstrom (x_unit, &cpy))
          goto finish;
     }

   if (info->default_labels)
     {
        if ((-1 == set_default_labels (version, info, fmt))
            || ((info->residuals != 0)
                && (label != NULL)
                && (-1 == Plot_set_labels (fmt, NULL, label, NULL))))
          {
             goto finish;
          }
     }

   /* Don't plot errorbars on models */
   eb_save = Plot_use_errorbars (fmt);
   if (is_model(version) && (info->residuals == 0))
     Plot_set_errorbars (fmt, 0, -1.0);

   /* usually don't plot residuals as logs or as bin density.. */
   ylog_save = Plot_get_logy (fmt);
   bd_save = Plot_bin_density (fmt);
   if (info->residuals != 0)
     {
        if (info->residuals != ISIS_RATIO_RESID)
          Plot_set_logy (fmt, 0);
        if (info->residuals != ISIS_DIFF_RESID)
          Plot_set_bin_density (fmt, 0);
     }

   if (info->use_style == NULL)
     info->use_style = h->color;

   ret = Plot_hist (fmt, info->use_style, info->overlay, &cpy);

   /* plot a horizontal reference line for residuals */
   switch (info->residuals)
     {
      case ISIS_STAT_RESID:
        /* drop */
      case ISIS_DIFF_RESID:
        plot_reference_line (0.0);
        break;
      case ISIS_RATIO_RESID:
        plot_reference_line (1.0);
        break;
      default:
        /* do nothing (should never happen) */
        break;
     }

   /* restore initial format state */
   Plot_set_errorbars (fmt, eb_save, -1.0);
   Plot_set_logy (fmt, ylog_save);
   Plot_set_bin_density (fmt, bd_save);

   finish:

   if (info->default_labels)
     (void) Plot_set_labels (fmt, "", "", NULL);

   free_plot_copy (&cpy);

   return ret;
}

/*}}}*/

/* Rebinning */

static int iterate_rebin_mask (Hist_t *h, int *rebin, int num_new, /*{{{*/
                               int (*start_bin)(Hist_t *, int, int),
                               int (*incr_bin)(Hist_t *, int, int))
{
   int k, n;

   if (rebin == NULL || h == NULL)
     return -1;

   k = 0;
   while (rebin[k] == 0) k++;

   for (n = 0; k < h->orig_nbins && n < num_new; n++)
     {
        int sign = rebin[k];

        if (start_bin)
          {
             if (-1 == (*start_bin)(h, k, n))
               return -1;
          }

        do
          {
             if (rebin[k] == 0)
               continue;
             if (rebin[k] != sign)
               break;

             if (incr_bin)
               {
                  if (-1 == (*incr_bin)(h, k, n))
                    {
                       return -1;
                    }
               }
          }
        while (++k < h->orig_nbins);
     }

   return n;
}

/*}}}*/

int Hist_rebin_index (Hist_t *h, void *s) /*{{{*/
{
   int *rebin = (int *)s;
   int n;

   /* s MUST point to an integer array of the same size
    * as the original histogram!
    */
   if (rebin == NULL || h == NULL)
     return -1;

   /* Count the number of bins (n) after rebinning. */
   if ((n = iterate_rebin_mask (h, rebin, h->orig_nbins, NULL, NULL)) < 0)
     return -1;

   memcpy ((char *)h->rebin, (char *)rebin, h->orig_nbins * sizeof(int));

   return n;
}

/*}}}*/

int Hist_rebin_min_counts (Hist_t *h, void * s) /*{{{*/
{
   double min_bin_counts, tot, remaining;
   int n, k, sign;

   if (s == NULL || h == NULL)
     return -1;

   min_bin_counts = *((double *)s);

   if (min_bin_counts < 0.0)
     return -1;

   /* revert to original binning */
   if (min_bin_counts == 0.0)
     {
        sign = 1;
        for (k = 0; k < h->orig_nbins; k++)
          {
             sign *= -1;
             h->rebin[k] = sign;
          }

        return h->orig_nbins;
     }

   remaining = 0.0;

   for (k = 0; k < h->orig_nbins; k++)
     {
        if (h->orig_notice[k])
          remaining += h->orig_counts[k];
     }

   if (remaining == 0.0)
     isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "zero total counts in dataset %d", h->index);

   n = 0;
   k = 0;
   sign = 1;

   for (n = 0; (n < h->orig_nbins) && (min_bin_counts < remaining); n++)
     {
        tot = 0.0;
        sign *= -1;

        for ( ; (k < h->orig_nbins) && (tot < min_bin_counts); k++)
          {
             if (h->orig_notice[k])
               {
                  h->rebin[k] = sign;
                  tot += h->orig_counts[k];
                  remaining -= h->orig_counts[k];
               }
             else h->rebin[k] = 0;
          }
     }

   if (remaining > 0.0)
     {
#if 0
        if (n < h->orig_nbins)
          n++;

        sign *= -1;
#else
        if (remaining >= min_bin_counts)
          {
             n++;
             sign *= -1;
          }
#endif
        for (; k < h->orig_nbins; k++)
          h->rebin[k] = sign;
     }
   else
     {
        /* include any remaining noticed zero content bins in the
         * last good bin, marking ignored bins as usual */
        for (; k < h->orig_nbins; k++)
          {
             h->rebin[k] = h->orig_notice[k] ? sign : 0;
          }
     }

   return n;
}

/*}}}*/

SLang_Name_Type *Hist_get_stat_error_hook (Hist_t *h) /*{{{*/
{
   return h ? h->stat_error_hook : NULL;
}

/*}}}*/

int Hist_set_stat_error_hook (Hist_t *h, SLang_Name_Type *hook, void (*delete_hook)(SLang_Name_Type *)) /*{{{*/
{
   if (h == NULL)
     return -1;
   if (h->stat_error_hook_delete)
     (*h->stat_error_hook_delete) (h->stat_error_hook);
   h->stat_error_hook = hook;
   h->stat_error_hook_delete = delete_hook;
   return 0;
}

/*}}}*/

static int slang_stat_error_hook (SLang_Name_Type *hook, /*{{{*/
                                  double *orig_counts, double *orig_stat_err,
                                  int *rebin, SLindex_Type orig_nbins,
                                  double *stat_err, int nbins)
{
   SLang_Array_Type *oc, *ose, *reb, *err;
   int ret = -1;

   oc = ose = reb = err = NULL;

   oc = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &orig_nbins, 1);
   ose = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &orig_nbins, 1);
   reb = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &orig_nbins, 1);
   if ((oc == NULL) || (ose == NULL) || (reb == NULL))
     goto finish;

   memcpy ((char *)oc->data, (char *)orig_counts, orig_nbins * sizeof(double));
   memcpy ((char *)ose->data, (char *)orig_stat_err, orig_nbins * sizeof(double));
   memcpy ((char *)reb->data, (char *)rebin, orig_nbins * sizeof(int));

   SLang_start_arg_list ();
   SLang_push_array (oc, 1);
   SLang_push_array (ose, 1);
   SLang_push_array (reb, 1);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (hook))
     goto finish;

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLdo_pop ();
        goto finish;
     }

   if (-1 == SLang_pop_array_of_type (&err, SLANG_DOUBLE_TYPE)
       || err == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "in rebin error hook");
        goto finish;
     }

   if ((int) err->num_elements != nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "wrong size array returned from rebin error hook", hook);
        goto finish;
     }

   memcpy ((char *)stat_err, (char *)err->data, nbins * sizeof(double));

   ret = 0;

   finish:
   if (ret)
     isis_throw_exception (Isis_Error);
   SLang_free_array (err);

   return ret;
}

/*}}}*/

static int run_stat_error_hook (Hist_t *h) /*{{{*/
{
   double *oc, *ose;
   int on, val_stat;

   if (h == NULL)
     return -1;

   on = h->orig_nbins;
   oc = h->orig_counts;
   ose = h->orig_stat_err;

   Isis_Active_Dataset = h->index;
   if (-1 == slang_stat_error_hook (h->stat_error_hook,
                                    oc, ose, h->rebin, on, h->stat_err, h->nbins))
     return -1;

   val_stat = validate_stat_err (h);
   if (-1 == val_stat)
     return -1;
   else if (1 == val_stat)
     invalid_uncertainties_replaced ();

   return 0;
}

/*}}}*/

static int update_stat_err (Hist_t *h) /*{{{*/
{
   int k, val_stat;

   if (h->stat_error_hook != NULL)
     return run_stat_error_hook (h);

   /* Assume Poisson */
   for (k = 0; k < h->nbins; k++)
     {
        h->stat_err[k] = sqrt (fabs(h->counts[k]));
     }

   val_stat = validate_stat_err (h);
   if (-1 == val_stat)
     return -1;
#if 0
   else if (1 == val_stat)
     invalid_uncertainties_replaced ();
#endif

   return 0;
}

/*}}}*/

static int bin_start (Hist_t *h, int old_index, int new_index) /*{{{*/
{
   h->counts[new_index] = 0.0;
   h->bin_lo[new_index] = h->orig_bin_lo[old_index];
   return 0;
}

/*}}}*/

static int bin_incr (Hist_t *h, int old_index, int new_index) /*{{{*/
{
   h->bin_hi[new_index] = h->orig_bin_hi[old_index];
   h->counts[new_index] += h->orig_counts[old_index];
   return 0;
}

/*}}}*/

int Hist_do_rebin (Hist_t *h, int (*rebin_fcn)(Hist_t *, void *), void *s)  /*{{{*/
{
   int n, k;

   if ((NULL == rebin_fcn) || (NULL == h))
     return -1;

   n = (*rebin_fcn)(h, s);
   if (n <= 0)
     return -1;

   h->n_notice = h->nbins = n;

   /* If the data is at the original binning, restore
    * the original notice flags;  otherwise, notice
    * everything by default.
    */

   if (h->nbins == h->orig_nbins)
     memcpy ((char *)h->notice, (char *)h->orig_notice, h->orig_nbins * sizeof(int));
   else
     {
        for (k = 0; k < h->nbins; k++)
          h->notice[k] = 1;
     }

   if (-1 == update_notice_list (h))
     return -1;

   if (-1 == iterate_rebin_mask (h, h->rebin, h->nbins, &bin_start, &bin_incr))
     return -1;

   if (-1 == update_stat_err (h))
     return -1;

   if (h->sys_err_frac)
     {
        double *sys_err_frac;
        if (NULL == (sys_err_frac = rebin_sys_err_frac (h, h->bin_lo, h->bin_hi, h->nbins)))
          return -1;
        ISIS_FREE(h->sys_err_frac);
        h->sys_err_frac = sys_err_frac;
     }

   if (h->flux_weights)
     {
        int valid;
        if (-1 == rebin_flux_using_weights (h, h->orig_flux, h->orig_flux_err, h->flux, h->flux_err))
          return -1;
        valid = validate_flux_err (h);
        if (valid == -1)
          return -1;
#if 0
        else if (valid == 1)
          invalid_uncertainties_replaced ();
#endif
     }

   return 0;
}

/*}}}*/

/* ARF/RMF check, assign and apply */

static int cmp_doubles (const void *va, const void *vb) /*{{{*/
{
   const double *a = (const double *) va;
   const double *b = (const double *) vb;

   if (*a < *b) return -1;
   else if (*a > *b) return 1;
   else return 0;
}
/*}}}*/

static int merge_grids (double *y, unsigned int num, double min_dy, Isis_Hist_t *x) /*{{{*/
{
   unsigned int i, j, num_uniq;
   double last_y;

   qsort (y, num, sizeof(double), &cmp_doubles);

   /* Provide user-adjustable override. */
   if (Hist_Min_Model_Spacing > 0.0)
     min_dy = Hist_Min_Model_Spacing;

   num_uniq = 1;
   last_y = y[0];
   for (j = 1; j < num; j++)
     {
        if (y[j] - last_y < min_dy)
          continue;

        last_y = y[j];
        y[num_uniq] = last_y;
        num_uniq++;
     }

   num_uniq--;

   if (-1 == Isis_Hist_allocate(num_uniq, x))
     return -1;

   for (i = 0; i < num_uniq; i++)
     {
        x->bin_lo[i] = y[i];
        x->bin_hi[i] = y[i+1];
     }

   return 0;
}

/*}}}*/

static int set_model_grid (Hist_t *h) /*{{{*/
{
   Isis_Hist_t x = ISIS_HIST_INIT;
   Isis_Rsp_t *rsp = &h->f_rsp;
   Isis_Rsp_t *r;
   unsigned int i, num;
   double *y, min_dy;

   if (rsp->next == NULL)
     {
        Isis_Hist_t *m = &h->model_flux;
        Isis_Arf_t *a = rsp->arf;

        /* allocate new before freeing the old */
        if (-1 == Isis_Hist_allocate (a->nbins, &x))
          return -1;
        Isis_Hist_free (m);

        /* struct copy */
        *m = x;

        memcpy ((char *)m->bin_lo, (char *)a->bin_lo, a->nbins*sizeof(double));
        memcpy ((char *)m->bin_hi, (char *)a->bin_hi, a->nbins*sizeof(double));

        return 0;
     }

   /* Generate a grid by combining all the ARFs.
    * Do this by sorting and keeping the unique points
    */
   num = 0;
   for (r = rsp; r != NULL; r = r->next)
     {
        num += r->arf->nbins + 1;
     }

   if (NULL == (y = (double *) ISIS_MALLOC (num * sizeof(double))))
     return -1;

   min_dy = DBL_MAX;
   i = 0;

   for (r = rsp; r != NULL; r = r->next)
     {
        double *bin_lo = r->arf->bin_lo;
        double *bin_hi = r->arf->bin_hi;
        unsigned int j, ny = r->arf->nbins;

        for (j = 0; j < ny; j++)
          {
             double dy = bin_hi[j] - bin_lo[j];
             if (dy < min_dy)
               min_dy = dy;
             y[i++] = bin_lo[j];
          }
        y[i++] = bin_hi[ny-1];
     }

   if (-1 == merge_grids (y, num, min_dy, &x))
     {
        ISIS_FREE(y);
        return -1;
     }
   ISIS_FREE (y);

   Isis_Hist_free (&h->model_flux);
   h->model_flux = x; /* struct copy */

   return 0;
}

/*}}}*/

static int initialize_rmf (Isis_Rmf_t *rmf, Hist_t *h) /*{{{*/
{
   Isis_Rmf_Grid_Type arf, ebounds;

   if (h == NULL)
     {
        memset ((char *)&ebounds, 0, sizeof(Isis_Rmf_Grid_Type));
        memset ((char *)&arf, 0, sizeof(Isis_Rmf_Grid_Type));
     }
   else
     {
        Isis_Arf_t *a = h->a_rsp.arf;

        ebounds.nbins = h->orig_nbins;
        ebounds.bin_lo = h->orig_bin_lo;
        ebounds.bin_hi = h->orig_bin_hi;

        if (a == NULL)
          {
             arf.nbins = 0;
             arf.bin_hi = NULL;
             arf.bin_lo = NULL;
          }
        else
          {
             arf.nbins = a->nbins;
             arf.bin_hi = a->bin_hi;
             arf.bin_lo = a->bin_lo;
          }
     }

   ebounds.units = U_ANGSTROM;
   arf.units = U_ANGSTROM;

   return Rmf_init_rmf (rmf, &arf, &ebounds);
}

/*}}}*/

static int grid_mismatch (double *bin_lo, double *bin_hi, unsigned int n, int order, /*{{{*/
                          double *rlo, double *rhi, unsigned int rn, int rorder,
                          double *max_lo, double *max_hi)
{
   unsigned int i;
   double of = 1.0;

   *max_lo = *max_hi = 0.0;

   if (n != rn)
     return 1;

   /* RMF EBOUNDS grid may need re-scaling */
   if (order != rorder)
     {
        double o1, o2;
        o1 = (double) ((rorder == 0) ? 1 : abs(rorder));
        o2 = (double) ((order == 0) ? 1 : abs(order));
        of = o1 / o2;
     }

   for (i = 0; i < rn; i++)
     {
        double lo, hi;

        hi = fabs(1.0 - of * rhi[i] / bin_hi[i]);
        lo = fabs(1.0 - of * rlo[i] / bin_lo[i]);
        *max_lo = MAX(lo, *max_lo);
        *max_hi = MAX(hi, *max_hi);
     }

   if ((fabs(*max_lo) > RMF_GRID_TOL)
       || (fabs(*max_hi) > RMF_GRID_TOL))
     return 1;

   return 0;
}

/*}}}*/

static int try_rmf_data_match (Hist_t *h, Isis_Rmf_t *rmf) /*{{{*/
{
   double err_lo, err_hi;
   double *lo, *hi;
   unsigned int n;
   int mismatch, energy_ordered_ebounds;

   if (h == NULL || rmf == NULL)
     return -1;

   lo = hi = NULL;

   if (-1 == rmf->get_data_grid (rmf, &lo, &hi, &n, &energy_ordered_ebounds))
     return -1;

   mismatch = grid_mismatch (h->orig_bin_lo, h->orig_bin_hi, h->orig_nbins, h->order,
                             lo, hi, n, rmf->order, &err_lo, &err_hi);
   ISIS_FREE (lo);
   ISIS_FREE (hi);

   if ((mismatch == 1) && (0 == Rmf_is_identity (rmf)))
     {
        mismatch = -1;
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "%d RMF bins vs %d data bins",
                    n, h->orig_nbins);
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__,
                    "Largest fractional errors:  %g (low edge)  %g (high edge)",
                    err_lo, err_hi);
     }

   return mismatch;
}

/*}}}*/

static int try_rmf_arf_match (Isis_Rmf_t *rmf, Isis_Arf_t *arf) /*{{{*/
{
   double err_lo, err_hi;
   double *lo, *hi;
   unsigned int n;
   int mismatch;

   if (rmf == NULL || arf == NULL)
     return -1;

   lo = hi = NULL;

   if (-1 == rmf->get_arf_grid (rmf, &lo, &hi, &n))
     return -1;

   mismatch = grid_mismatch (arf->bin_lo, arf->bin_hi, arf->nbins, 1,
                             lo, hi, n, 1, &err_lo, &err_hi);
   ISIS_FREE (lo);
   ISIS_FREE (hi);

   if ((mismatch == 1)
       && ((0 == Arf_is_identity(arf)) && (0 == Rmf_is_identity (rmf))))
     {
        mismatch = -1;
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "%d RMF bins vs %d ARF bins",
                    n, arf->nbins);
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__,
                    "Largest fractional errors:  %g (low edge)  %g (high edge)",
                    err_lo, err_hi);
     }

   return mismatch;
}

/*}}}*/

static int adjust_identity_rsp (Hist_t *h, Isis_Rmf_t **rmf, Isis_Arf_t **arf) /*{{{*/
{
   /* An attempt has been made to assign the responses in rsp
    * to histogram 'h', but the full set (h, rmf, arf) didnt match.
    * If either 'rmf' or 'arf' is an identity, it can be
    * re-created to guarantee a matching set.
    *
    * It is an error to call this routine with
    * both rmf != identity and arf != identity
    */

   if (Rmf_is_identity (*rmf))
     {
        release_rmf (*rmf);
        *rmf = matching_identity_rmf (h, *arf);
        if (NULL == *rmf)
          return -1;
        (*rmf)->ref_count++;
     }
   else if (Arf_is_identity (*arf))
     {
        double *lo, *hi;
        unsigned int n;

        if (-1 == (*rmf)->get_arf_grid (*rmf, &lo, &hi, &n))
          return -1;

        release_arf (*arf);
        *arf = Arf_make_identity_arf (lo, hi, n);
        ISIS_FREE (lo); ISIS_FREE (hi);
        if (*arf == NULL)
          return -1;
        (*arf)->ref_count++;
     }
   else
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Internal error: [adjust_identity_rsp]");
        return -1;
     }

   return 0;
}

/*}}}*/

static int set_fit_response (Hist_t *h, Isis_Rsp_t *rsp) /*{{{*/
{
   if ((h == NULL) || (rsp == NULL)
       || (rsp->arf == NULL) || (rsp->rmf == NULL))
     return -1;

   release_rsp (&h->f_rsp);
   h->f_rsp = *rsp;  /* struct copy */
   incr_rsp_refcount (&h->f_rsp);

   return 0;
}

/*}}}*/

static void warn_single_rsp (Isis_Rsp_t *rsp) /*{{{*/
{
   if (rsp->next != NULL)
     fprintf (stdout, "** Fit will use a single ARF/RMF pair\n");
   rsp->next = NULL;
}

/*}}}*/

int Hist_set_fit_responses (Hist_t *h, unsigned int fit_data_type, /*{{{*/
                            unsigned int response_mask)
{
   Isis_Rsp_t rsp;
   unsigned int is_flux_data = is_flux(fit_data_type);

   if (h == NULL)
     return -1;

   rsp = h->a_rsp;  /* struct copy */

   rsp.is_flux_data = is_flux_data;

   /* The assigned response might be a list of ARF/RMF pairs.
    * However, if the assigned response is not used, the fit
    * will use only a single ARF/RMF pair.
    */
   if ((/* is_flux_data || */ want_ideal_arf(response_mask))
        && (0 == Arf_is_identity (rsp.arf)))
     {
        rsp.arf = identity_arf_on_arf_grid (h);
        warn_single_rsp (&rsp);
     }

   if (want_ideal_rmf(response_mask) && (0 == Rmf_is_identity (rsp.rmf)))
     {
        rsp.rmf = matching_identity_rmf (h, rsp.arf);
        warn_single_rsp (&rsp);
     }

#if 0
   if (is_flux_data
       && Rmf_includes_effective_area (rsp.rmf))
     {
        isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "instrument response inconsistent with fit data type");
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, " ==> Use factor_rsp() to factor out the effective area.");
        return -1;
     }
#endif

   return set_fit_response (h, &rsp);
}

/*}}}*/

static int attach_rsp (Hist_t *h, Isis_Rsp_t *rsp) /*{{{*/
{
   /* Keep the head of the list, free the rest */
   if (h->a_rsp.next != NULL)
     {
        map_rsp_list (h->a_rsp.next, &release_rsp);
        map_rsp_list (h->a_rsp.next, &free_rsp);
        h->a_rsp.next = NULL;
     }

   /* Handle the head of the list */
   if (rsp->rmf != h->a_rsp.rmf)
     {
        release_rmf (h->a_rsp.rmf);
        h->a_rsp.rmf = rsp->rmf;
        rsp->rmf->ref_count++;
     }

   if (rsp->arf != h->a_rsp.arf)
     {
        release_arf (h->a_rsp.arf);
        h->a_rsp.arf = rsp->arf;
        rsp->arf->ref_count++;
     }

   /* Bump the reference count on the rest of the list */
   if (rsp->next != NULL)
     {
        h->a_rsp.next = rsp->next;
        map_rsp_list (h->a_rsp.next, &incr_rsp_refcount);
     }

   if (-1 == set_fit_response (h, rsp))
     return -1;

   return set_model_grid (h);
}

/*}}}*/

static int match_rsp (Hist_t *h, Isis_Rsp_t *rsp) /*{{{*/
{
   int ma, md;

   if ((ma = try_rmf_arf_match (rsp->rmf, rsp->arf)) < 0)
     return -1;

   if ((md = try_rmf_data_match (h, rsp->rmf)) < 0)
     return -1;

   return (ma || md);
}

/*}}}*/

static int assign_matching_rsp (Hist_t *h, Isis_Rmf_t **rmf, Isis_Arf_t **arf) /*{{{*/
{
   Isis_Rsp_t rsp;
   int adjust;

   rsp.next = NULL;
   rsp.arf = *arf;
   rsp.rmf = *rmf;

   if ((0 != Rmf_includes_effective_area (*rmf))
       && (0 == Arf_is_identity (*arf)))
     {
        if (Hist_Allow_Multiple_Arf_Factors < 0)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "RMF already includes effective area; ARF not assigned");
             isis_vmesg (WARN, I_INFO, __FILE__, __LINE__,
                         "To force the ARF assignment, set Allow_Multiple_Arf_Factors=0 or 1");
             return -1;
          }
        else if (Hist_Allow_Multiple_Arf_Factors == 0)
          {
             isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "RMF already includes effective area; assigning ARF anyway");
          }
     }

   if ((adjust = match_rsp (h, &rsp)) < 0)
     return -1;

   if (adjust)
     {
        if (-1 == adjust_identity_rsp (h, rmf, arf))
          return -1;
     }

   /* arf or rmf might have changed. */
   rsp.next = NULL;
   rsp.arf = *arf;
   rsp.rmf = *rmf;

   return attach_rsp (h, &rsp);
}

/*}}}*/

static Isis_Arf_t *get_arf (Hist_t *h, Isis_Arf_t *ahead, int arf_index, int *malloced) /*{{{*/
{
   Isis_Arf_t *a;

   if (malloced != NULL)
     *malloced = 0;

   if (arf_index == 0)
     {
        a = Arf_make_identity_arf (h->orig_bin_lo, h->orig_bin_hi, h->orig_nbins);
        if (a == NULL)
          return NULL;
        if (malloced != NULL)
          *malloced = 1;
        return a;
     }

   return Arf_find_arf_index (ahead, arf_index);
}

/*}}}*/

int Hist_assign_arf_to_hist (Hist_t *h_head, int hist_index, Isis_Arf_t *ahead, int arf_index) /*{{{*/
{
   Hist_t *h;
   Isis_Arf_t *a;
   int fake_data = 0;
   int malloced = 0;

   if (NULL == (h = _Hist_find_hist_index (h_head, hist_index)))
     {
        if (arf_index <= 0)
          return -1;
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__,
                    "generating empty data set to match ARF %d",
                    arf_index);
        fake_data = 1;
     }

   if (NULL == (a = get_arf (h, ahead, arf_index, &malloced)))
     return -1;

   if (fake_data != 0)
     {
        if (NULL == (h = create_hist_from_arf (h_head, hist_index, a)))
          goto fail;
     }

   if (0 == assign_matching_rsp (h, &h->a_rsp.rmf, &a))
     return 0;

   fail:
   if (malloced)
     Arf_free_arf (a);
   return -1;
}

/*}}}*/

static Isis_Rmf_t *get_rmf (Hist_t *h, Isis_Rmf_t *rhead, int rmf_index, int *malloced) /*{{{*/
{
   Isis_Rmf_t *r;

   if (malloced != NULL)
     *malloced = 0;

   if (rmf_index == 0)
     {
        r = matching_identity_rmf (h, h->a_rsp.arf);
        if (malloced != NULL)
          *malloced = 1;
        return r;
     }

   r = Rmf_find_rmf_index (rhead, rmf_index);
   if (r == NULL)
     return NULL;
   if (-1 == initialize_rmf (r, h))
     return NULL;

   return r;
}

/*}}}*/

static Hist_t *fake_hist_with_rmf (Hist_t *h, Hist_t *h_head, int hist_index, Isis_Rmf_t *rmf) /*{{{*/
{
   if (h == NULL)
     return create_hist_from_rmf (h_head, hist_index, rmf);

   if (-1 == redefine_fake_data_grid (h_head, rmf, &h))
     return NULL;

   return h;
}

/*}}}*/

static int rmf_chan_grid_size (Isis_Rmf_t *r) /*{{{*/
{
   double *lo, *hi;
   unsigned int n;
   int eoe;

   if (r == NULL)
     return -1;

   lo = hi = NULL;

   if (-1 == r->get_data_grid (r, &lo, &hi, &n, &eoe))
     return -1;

   ISIS_FREE (lo);
   ISIS_FREE (hi);

   return n;
}

/*}}}*/

int Hist_assign_rmf_to_hist (Hist_t *h_head, int hist_index, Isis_Rmf_t *rhead, int rmf_index) /*{{{*/
{
   Hist_t *h;
   Isis_Rmf_t *r;
   int malloced = 0;
   int fake_data = 0;

   if (NULL == (h = _Hist_find_hist_index (h_head, hist_index)))
     {
        if (rmf_index <= 0)
          return -1;
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__,
                    "generating empty data set to match RMF %d",
                    rmf_index);
        fake_data = 1;
     }

   /* If this is a no-op, return early for two reasons:
    *  1) save time
    *  2) avoid possibility that the data grid order might
    *     get erroneously swapped
    */
   if ((rmf_index == 0) && (0 != Rmf_is_identity (h->a_rsp.rmf)))
     return 0;

   if (NULL == (r = get_rmf (h, rhead, rmf_index, &malloced)))
     return -1;

   /* If fake data got initialized with a wrong-sized ARF grid and then
    * gets assigned an RMF grid, we'll just force it to fit.
    */
   if ((fake_data != 0)
       || ((h->is_fake_data != 0)
           && (h->orig_nbins != rmf_chan_grid_size (r))))
     {
        if (NULL == (h = fake_hist_with_rmf (h, h_head, hist_index, r)))
          goto fail;
     }
   else
     {
        /* FIXME(not?)  Using the RMF to auto-set the data grid
         * may lead to some subtle problems.  I'm not yet sure if I
         * want to leave it this way.
         */
        if (/* (h->has_grid == -1) && */
            (-1 == set_hist_grid_using_rmf (r, h)))
          goto fail;
     }

   if (0 != assign_matching_rsp (h, &r, &h->a_rsp.arf))
     goto fail;

   if ((h->bgd_file == NULL) || (*h->bgd_file == 0))
     return 0;

   /* Re-load the background file to make sure the grids still match */
   if (-1 == Hist_set_background_from_file (h, h->bgd_file))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "re-assigning background spectrum");
        goto fail;
     }

   return 0;

   fail:
   if (malloced)
     Rmf_free_rmf (r);
   return -1;
}

/*}}}*/

int Hist_assign_rsp_list (Hist_t *hhead, int hist_index, /*{{{*/
                          Isis_Arf_t *ahead, int *arfs,
                          Isis_Rmf_t *rhead, int *rmfs, int num)
{
   Hist_t *h;
   Isis_Rsp_t *x = NULL;
   Isis_Rsp_t **last = &x;
   Isis_Rsp_t *rsp = NULL;
   int i, fake_data = 0;
   int first = 1;

   if (num == 0)
     return 0;

   if ((hhead == NULL)
       || ((ahead == NULL) || (arfs == NULL))
       || ((rhead == NULL) || (rmfs == NULL))
       || (num < 0))
     return -1;

   h = _Hist_find_hist_index (hhead, hist_index);
   if (NULL == h)
     {
        if (rmfs[0] <= 0)
          return -1;
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__,
                    "generating empty data set to match assigned RMFs");
        fake_data = 1;
     }

   for (i = 0; i < num; i++)
     {
        if (NULL == (rsp = (Isis_Rsp_t *) ISIS_MALLOC (sizeof(Isis_Rsp_t))))
          goto fail;

        rsp->next = NULL;
        rsp->arf = get_arf (h, ahead, arfs[i], NULL);
        rsp->rmf = get_rmf (h, rhead, rmfs[i], NULL);
        if ((rsp->arf == NULL) || (rsp->rmf == NULL))
          goto fail;

        if (first == 0)
          {
             if (0 != match_rsp (h, rsp))
               goto fail;
          }
        else
          {
             first = 0;
             if (fake_data || h->is_fake_data)
               {
                  h = fake_hist_with_rmf (h, hhead, hist_index, rsp->rmf);
                  if (h == NULL)
                    goto fail;
               }
             else
               {
                  if (-1 == set_hist_grid_using_rmf (rsp->rmf, h))
                    goto fail;
               }
          }

        *last = rsp;
        last = &rsp->next;
     }

   if (-1 == attach_rsp (h, x))
     goto fail;

   free_rsp (x);

   return 0;

   fail:
   release_rsp (rsp);
   free_rsp (rsp);
   map_rsp_list (x, &release_rsp);
   map_rsp_list (x, &free_rsp);
   return -1;
}

/*}}}*/

/* Support model caching during fit */

Hist_Eval_Grid_Method_Type *Hist_eval_grid_method (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return NULL;
   return &h->eval_grid_method;
}

/*}}}*/

int Hist_set_eval_grid_info (Hist_t *head, unsigned int *indices, unsigned int num_indices, /*{{{*/
                             Hist_Eval_Grid_Method_Type *m)
{
   unsigned int i, id;
   Hist_t *h;

   if (head == NULL || indices == NULL)
     return -1;

   if (m->cache_model_values && ISIS_EVAL_GRID_SEPARATE)
     {
        if (0 == all_member_grids_are_identical (head, indices, num_indices))
          {
             isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "To cache model values, all grids must be identical");
             return -1;
          }
     }

   id = Next_Eval_Grid_Id;

   for (i = 0; i < num_indices; i++)
     {
        Hist_Eval_Grid_Method_Type *x;
        if (NULL == (h = (_Hist_find_hist_index (head, indices[i]))))
          continue;
        x = &h->eval_grid_method;
        free_eval_grid_method (x);
        *x = *m;                        /* struct copy */
        x->id = x->cache_model_values ? id : id++;
     }

   Next_Eval_Grid_Id = id+1;

   return 0;
}

/*}}}*/

static int merge_model_grids (Hist_t *head, unsigned int *indices, unsigned int num_indices, Isis_Hist_t *m) /*{{{*/
{
   Hist_t *h;
   double *y;
   double min_dy;
   unsigned int i, j, k, num;

   if (head == NULL || indices == NULL || m == NULL)
     return -1;

   num = 0;
   for (i = 0; i < num_indices; i++)
     {
        if (NULL == (h = (_Hist_find_hist_index (head, indices[i]))))
          continue;
        num += h->model_flux.nbins + 1;
     }

   if (NULL == (y = (double *) ISIS_MALLOC (num * sizeof(double))))
     return -1;

   min_dy = DBL_MAX;
   k = 0;

   for (i = 0; i < num_indices; i++)
     {
        Isis_Hist_t *f;
        double *bin_lo, *bin_hi;
        unsigned int ny;

        if (NULL == (h = (_Hist_find_hist_index (head, indices[i]))))
          continue;

        f = &h->model_flux;
        bin_lo = f->bin_lo;
        bin_hi = f->bin_hi;
        ny = f->nbins;

        for (j = 0; j < ny; j++)
          {
             double dy = bin_hi[j] - bin_lo[j];
             if (dy < min_dy)
               min_dy = dy;
             y[k++] = bin_lo[j];
          }
        y[k++] = bin_hi[ny-1];
     }

   if (-1 == merge_grids (y, num, min_dy, m))
     {
        ISIS_FREE(y);
        return -1;
     }
   ISIS_FREE(y);

   for (i = 0; i < num_indices; i++)
     {
        Isis_Hist_t *x;
        if (NULL == (h = (_Hist_find_hist_index (head, indices[i]))))
          continue;

        x = &h->model_flux;
        if (-1 == transfer_notice (x->bin_lo, x->bin_hi, x->notice_list, x->n_notice,
                                   m->bin_lo, m->bin_hi, m->nbins, m->notice))
          return -1;
     }

   return _update_notice_list (m->notice, &m->notice_list, &m->n_notice, m->nbins);
}

/*}}}*/

static int allocate_temporary_workspace (Isis_Hist_t *g) /*{{{*/
{
   double *tmp;
   if (g == NULL)
     return -1;
   if (NULL == (tmp = (double *) ISIS_REALLOC (g->val, 2*g->nbins*sizeof(double))))
     return -1;
   g->val = tmp;
   memset ((char *)g->val, 0, 2*g->nbins*sizeof(double));
   return 0;
}

/*}}}*/

int Hist_make_merged_eval_grid (Hist_t *head, int eval_grid_id, Isis_Hist_t *grid) /*{{{*/
{
   Hist_t *h;
   unsigned int *id;
   unsigned int num;

   if ((head == NULL) || (grid == NULL))
     return -1;

   if (0 == (num = Hist_num_noticed_histograms (head)))
     return -1;

   if (NULL == (id = (unsigned int *) ISIS_MALLOC (num * sizeof(unsigned int))))
     return -1;

   num = 0;
   for (h = head->next; h != NULL; h = h->next)
     {
        Hist_Eval_Grid_Method_Type *m = &h->eval_grid_method;

        if ((h->exclude == 0) && (h->n_notice > 0)
            && (m->type == ISIS_EVAL_GRID_MERGED)
            && (m->id == eval_grid_id))
          {
             id[num++] = h->index;
          }
     }

   if (-1 == merge_model_grids (head, id, num, grid))
     {
        ISIS_FREE(id);
        return -1;
     }

   /* need extra temp space during the fit. */
   if (-1 == allocate_temporary_workspace (grid))
     {
        ISIS_FREE(id);
        return -1;
     }

   ISIS_FREE(id);
   return 0;
}

/*}}}*/

int Hist_finalize_user_eval_grid (Hist_t *h, Isis_Hist_t *g) /*{{{*/
{
   Hist_Eval_Grid_Method_Type *m;
   Isis_Hist_t *x;
   double xmin, xmax;
   int i;

   if (h == NULL || g == NULL)
     return -1;

   if ((h->exclude != 0) || (h->n_notice == 0))
     return 0;

   m = &h->eval_grid_method;

   if (m->type != ISIS_EVAL_GRID_USER)
     return -1;

   x = &h->model_flux;

   if (-1 == transfer_notice (x->bin_lo, x->bin_hi, x->notice_list, x->n_notice,
                              g->bin_lo, g->bin_hi, g->nbins, g->notice))
     return -1;

   /* make sure bins extending beyond the original grid
    * are also noticed.
    */

   xmin = x->bin_lo[0];
   for (i = 0; g->bin_lo[i] < xmin; i++)
     g->notice[i] = 1;

   xmax = x->bin_hi[x->nbins-1];
   for (i = g->nbins-1; xmax < g->bin_hi[i]; i--)
     g->notice[i] = 1;

   if (-1 == _update_notice_list (g->notice, &g->notice_list, &g->n_notice, g->nbins))
     return -1;

   /* need extra temp space during the fit. */
   return allocate_temporary_workspace (g);
}

/*}}}*/

typedef struct
{
   SLang_Array_Type *bin_lo;
   SLang_Array_Type *bin_hi;
}
Model_Eval_Grid_Type;

static SLang_CStruct_Field_Type Model_Eval_Grid_Type_Layout[] =
{
   MAKE_CSTRUCT_FIELD (Model_Eval_Grid_Type, bin_lo, "bin_lo", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Model_Eval_Grid_Type, bin_hi, "bin_hi", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

int Hist_make_user_eval_grid (Hist_t *h, Isis_Hist_t *m_grid) /*{{{*/
{
   Hist_Eval_Grid_Method_Type *egm;
   Model_Eval_Grid_Type grid;
   double *lo, *hi;
   int n;
   int status = -1;

   if ((h == NULL) || (m_grid == NULL))
     return -1;

   if (NULL == (egm = Hist_eval_grid_method (h)))
     return -1;

   if (NULL == egm->options)
     return -1;

   memset ((char *)&grid, 0, sizeof(grid));

   SLang_start_arg_list ();
   (void) SLang_push_integer (h->index);
   (void) SLang_push_cstruct ((VOID_STAR)&grid, Model_Eval_Grid_Type_Layout);
   SLang_end_arg_list ();

   (void) SLexecute_function ((SLang_Name_Type *)egm->options);

   if ((-1 == SLang_pop_cstruct ((VOID_STAR)&grid, Model_Eval_Grid_Type_Layout))
       || (grid.bin_lo == NULL)
       || (grid.bin_hi == NULL))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "handling Struct_Type returned from user eval grid function");
        goto error_return;
     }

   if (grid.bin_lo->num_elements != grid.bin_hi->num_elements)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "mismatched array sizes");
        goto error_return;
     }

   n = (int) grid.bin_lo->num_elements;
   lo = (double *)grid.bin_lo->data;
   hi = (double *)grid.bin_hi->data;

   if (-1 == validate_wavelength_grid (n, lo, hi))
     goto error_return;

   if (-1 == Isis_Hist_allocate (n, m_grid))
     goto error_return;

   memcpy ((char *)m_grid->bin_lo, (char *)lo, n*sizeof(double));
   memcpy ((char *)m_grid->bin_hi, (char *)hi, n*sizeof(double));

   if (-1 == Hist_finalize_user_eval_grid (h, m_grid))
     goto error_return;

   status = 0;

   error_return:
   SLang_free_array (grid.bin_lo);
   SLang_free_array (grid.bin_hi);

   return status;
}

/*}}}*/

/* apply rebin and notice */

int Hist_apply_rebin (double *x, Hist_t *h, double **rebinned, int *nbins) /*{{{*/
{
   double *xn = NULL;

   if (h == NULL)
     return -1;

   if (NULL == (xn = (double *) ISIS_MALLOC (h->nbins * sizeof(double))))
     return -1;

   if (-1 == apply_rebin (x, h->orig_nbins, h->rebin, h->nbins, xn))
     {
        ISIS_FREE (xn);
        return -1;
     }

   *rebinned = xn;
   if (nbins != NULL) *nbins = h->nbins;

   return 0;
}

/*}}}*/

int Hist_apply_rebin_and_notice_list (double *bin_and_notice_result, double *x, Hist_t *h)/*{{{*/
{
   double *xn;
   int allocated = 0;
   int k;

   /* optionally rebin the result */
   if (h->nbins == h->orig_nbins)
     xn = x;
   else
     {
        allocated = 1;
        if (-1 == Hist_apply_rebin (x, h, &xn, NULL))
          return -1;
     }

   /* apply notice list */
   for (k = 0; k < h->n_notice; k++)
     {
        int j = h->notice_list[k];
        bin_and_notice_result[k] = xn[j];
     }

   if (allocated)
     ISIS_FREE (xn);

   return 0;
}

/*}}}*/

/* response kernel */

static void free_kernel (Hist_t *h) /*{{{*/
{
   Isis_Kernel_t *k;

   if (h == NULL || h->kernel == NULL)
     return;

   k = h->kernel;
   ISIS_FREE(k->params);

   if (k->delete_kernel != NULL)
     k->delete_kernel (h->kernel);

   h->kernel = NULL;
}

/*}}}*/

Isis_Kernel_t *isis_init_kernel (Isis_Kernel_t *k, unsigned int sizeof_kernel, /*{{{*/
                                 Isis_Obs_t *o)
{
   if (o == NULL)
     return NULL;

   if (k == NULL)
     {
        if (NULL == (k = (Isis_Kernel_t *) ISIS_MALLOC (sizeof_kernel)))
          return NULL;
        memset ((char *) k, 0, sizeof_kernel);
     }

   /* struct copy */
   k->rsp = o->rsp;

   k->exposure_time = o->exposure_time;
   k->apply_rmf = o->apply_rmf;
   k->num_orig_data = o->num_orig_data;
   k->params = NULL;

   k->frame_time = o->frame_time;
   k->tg_part = o->tg_part;
   k->tg_m = o->tg_m;

   return k;
}
/*}}}*/

static int get_exposure_time (Hist_t *h, double *t) /*{{{*/
{
   *t = h->f_rsp.arf->exposure;

   if (Rmf_includes_effective_area (h->f_rsp.rmf)
       || Arf_is_identity (h->f_rsp.arf))
     {
        if (h->f_rsp.is_flux_data)
          *t = 1.0;
        else *t = h->exposure;
     }

   if (*t <= 0.0)
     {
        isis_vmesg (INFO, I_WARNING, __FILE__, __LINE__, "dataset %d, ARF exposure time not set", h->index);
        if (h->exposure > 0.0)
          {
             isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, " => defaulting to data exposure time = %0.4g", h->exposure);
             *t = h->exposure;
          }
        else
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "dataset %d, invalid exposure time = %0.4g sec",
                         h->index, h->exposure);
             return -1;
          }
     }

   return 0;
}

/*}}}*/

static int set_observation_info (Isis_Obs_t *o, Hist_t *h) /*{{{*/
{
   /* struct copy (fit-responses) */
   o->rsp = h->f_rsp;

   o->apply_rmf = &Rmf_apply_rmf;
   o->num_orig_data = h->orig_nbins;
   o->frame_time = h->frame_time;
   o->tg_part = h->part;
   o->tg_m = h->order;

   if (-1 == get_exposure_time (h, &o->exposure_time))
     return -1;

   return 0;
}

/*}}}*/

int Hist_reallocate_kernel (Hist_t *h, char *params, Isis_Kernel_Def_t *def) /*{{{*/
{
   Isis_Kernel_t *k = NULL;
   Isis_Obs_t o;

   if ((h == NULL) || (def == NULL) || (def->allocate_kernel == NULL))
     return -1;

   if (-1 == set_observation_info (&o, h))
     return -1;

   if (NULL == (k = def->allocate_kernel (&o, params)))
     return -1;

   if (params == NULL)
     k->params = NULL;
   else if (NULL == (k->params = isis_make_string (params)))
     {
        k->delete_kernel (k);
        return -1;
     }

   k->kernel_def = def;

   if (h->kernel != NULL)
     free_kernel (h);

   h->kernel = k;

   return 0;
}

/*}}}*/

int Hist_print_kernel (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return -1;

   if (h->kernel == NULL)
     return -1;

   if (h->kernel->print_kernel == NULL)
     return 0;

   return (*h->kernel->print_kernel) (h->kernel);
}

/*}}}*/

int Hist_load_aux_kernel (char *libfile, char *init_name, char *init_args,  /*{{{*/
                          Isis_Kernel_Def_t *def)
{
   Isis_Kernel_Init_t *init_kernel = NULL;

   if (libfile == NULL || init_name == NULL || def == NULL)
     return -1;

   def->kernel_name = "";
   def->num_kernel_parms = 0;
   def->allocate_kernel = NULL;

   init_kernel = (Isis_Kernel_Init_t *) isis_load_function (libfile, init_name, "kernel");
   if (init_kernel == NULL)
     return -1;

   if (0 != (*init_kernel)(def, init_args))
     return -1;

   return 0;
}

/*}}}*/

Isis_Kernel_t *Hist_get_kernel (Hist_t *h) /*{{{*/
{
   if (h == NULL)
     return NULL;

   return h->kernel;
}

/*}}}*/

char *Hist_get_kernel_params (Hist_t *h) /*{{{*/
{
   if (h == NULL || h->kernel == NULL)
     return NULL;

   return h->kernel->params;
}

/*}}}*/

/* flux-correct */

int Hist_get_flux_corr_weights (Hist_t *h, double **weights, int *num_weights) /*{{{*/
{
   int len;

   if ((h == NULL) || (weights == NULL) || (num_weights == NULL))
     return -1;

   if (h->flux_weights == NULL)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "flux-correction weights have not yet been computed for dataset %d", h->index);
        return -1;
     }

   len = h->orig_nbins * sizeof(double);
   if (NULL == (*weights = (double *) ISIS_MALLOC (len)))
     return -1;
   *num_weights = h->orig_nbins;

   memcpy ((char *)*weights, (char *)h->flux_weights, len);

   return 0;
}

/*}}}*/

static int setup_flux_correct (Hist_t *h, int using_model, Isis_Hist_t *c, /*{{{*/
                               double **flux, double **flux_err)
{
   double *f=NULL, *f_err=NULL;
   int i, n;

   if (h == NULL || c == NULL)
     return -1;

   n = h->orig_nbins;
   c->nbins = h->orig_nbins;
   c->bin_lo = h->orig_bin_lo;
   c->bin_hi = h->orig_bin_hi;

   if (using_model == 0)
     {
        /* numerator is observed counts */

        c->val = h->orig_counts;
        c->val_err = h->orig_stat_err;
        *flux = h->orig_flux;
        *flux_err = h->orig_flux_err;

        return 0;
     }

   /* numerator is model counts */

   if (h->nbins != h->orig_nbins)
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "*** operation not implemented for binned data");
        return -1;
     }

   if ((NULL == (f = (double *) ISIS_MALLOC (n * sizeof(double))))
       || (NULL == (f_err = (double *) ISIS_MALLOC (n * sizeof(double))))
       || (NULL == (c->val_err = (double *) ISIS_MALLOC (n* sizeof(double)))))
     {
        ISIS_FREE(f);
        ISIS_FREE(f_err);
        ISIS_FREE(c->val_err);
        return -1;
     }

   *flux = f;
   *flux_err = f_err;
   c->val = h->model_counts;

   for (i = 0; i < n; i++)
     {
        c->val_err[i] = sqrt(fabs(c->val[i]));
     }

   return 0;
}

/*}}}*/

static int set_flux_floor (Hist_t *h, double frac, double *flux) /*{{{*/
{
   double *net=NULL, *err=NULL;
   int i, malloced = 0;;

   if (h->orig_bgd == NULL)
     {
        net = h->counts;
        err = h->stat_err;
     }
   else
     {
        double *b=NULL, *berr=NULL;

        if (-1 == scale_background (h, 1, &b, &berr))
          return -1;

        if (NULL == (net = (double *) ISIS_MALLOC (2 * h->nbins * sizeof(double))))
          {
             ISIS_FREE(b);
             ISIS_FREE(berr);
             return -1;
          }
        malloced = 1;
        err = net + h->nbins;

        for (i = 0; i < h->nbins; i++)
          {
             net[i] = h->counts[i] - b[i];
             err[i] = isis_hypot (h->stat_err[i], berr[i]);
          }

        if (h->sys_err_frac)
          {
             for (i = 0; i < h->nbins; i++)
               {
                  err[i] = isis_hypot (err[i], h->sys_err_frac[i] * h->counts[i]);
               }
          }

        ISIS_FREE(b);
        ISIS_FREE(berr);
     }

   for (i = 0; i < h->nbins; i++)
     {
        if (net[i] < frac * err[i])
          flux[i] = 0.0;
     }

   if (malloced)
     ISIS_FREE(net);

   return 0;
}

/*}}}*/

double *isis_unit_source_response (Isis_Kernel_t *k) /*{{{*/
{
   double *arf_rmf_dlam = NULL;
   double *arf_dlam = NULL;
   double *arf = NULL;
   int i, n, arf_nbins;
   Isis_Arf_t *arft;

   if (k == NULL)
     return NULL;

   if (k->rsp.next != NULL)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "flux-correction for multiple responses is not supported");
        return NULL;
     }

   arft = k->rsp.arf;
   arf = arft->arf;
   arf_nbins = arft->nbins;

   if (NULL == (arf_dlam = (double *) ISIS_MALLOC (arf_nbins * sizeof(double))))
     return NULL;

   for (i = 0; i < arf_nbins; i++)
     {
        arf_dlam[i] = arf[i] * (arft->bin_hi[i] - arft->bin_lo[i]);
     }

   n = k->num_orig_data;
   if (NULL == (arf_rmf_dlam = (double *) ISIS_MALLOC (n * sizeof(double))))
     {
        ISIS_FREE (arf_dlam);
        return NULL;
     }

   memset ((char *) arf_rmf_dlam, 0, n * sizeof (double));

   if (-1 == k->apply_rmf (k->rsp.rmf, arf_rmf_dlam, n, arf_dlam, NULL, arf_nbins))
     {
        ISIS_FREE (arf_rmf_dlam);
        ISIS_FREE (arf_dlam);
        return NULL;
     }

   ISIS_FREE (arf_dlam);

   for (i = 0; i < n; i++)
     {
        arf_rmf_dlam[i] *= k->exposure_time;
     }

   /* t \int dy R(h,y) A(y) */
   return arf_rmf_dlam;
}

/*}}}*/

static int rebin_flux_using_weights (Hist_t *h, /*{{{*/
                                     double *old_flux, double *old_flux_err,
                                     double *new_flux, double *new_flux_err)
{
   double *lo, *hi;
   int i, n;

   if (h->flux_weights == NULL)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "cannot rebin flux-corrected data for unknown response");
        return -1;
     }

   /* f(h) = (D(h)-B(h)) / t \int dy R(h,y) A(y)   [photons/s/cm^2/A]
    * The weights are defined so that the rebinned result yields
    * F(H) [photons/s/cm^2/A] on the big bins.
    */
   if (-1 == apply_rebin_weights (old_flux, h->flux_weights, 1,
                                  h->orig_nbins, h->rebin, h->nbins, new_flux))
     return -1;

   n = h->nbins;
   lo = h->bin_lo;
   hi = h->bin_hi;

   /* multiply by the bin width, dy, to yield the rebinned flux
    * F(H) in bin-integrated units, photons/s/cm^2/bin
    */
   for (i = 0; i < n; i++)
     {
        new_flux[i] *= hi[i] - lo[i];
     }

   if (new_flux_err == NULL)
     return 0;

   if (-1 == apply_rebin_weights (old_flux_err, h->flux_weights, 2,
                                  h->orig_nbins, h->rebin, h->nbins, new_flux_err))
     return -1;

   /* new_flux_err is quadrature sum so it can't be < 0 */
   for (i = 0; i < n; i++)
     {
        new_flux_err[i] = (hi[i] - lo[i]) * sqrt(new_flux_err[i]);
     }

   return 0;
}

/*}}}*/

static int verify_flux_corr_kernel (Isis_Kernel_t *k) /*{{{*/
{
   if (NULL == k->compute_flux)
     {
        isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__, "operation not supported by assigned kernel");
        return -1;
     }

   if ((k->rsp.arf != NULL)
       && (0 != Arf_is_identity (k->rsp.arf))
       && (0 == Rmf_includes_effective_area (k->rsp.rmf)))
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "ARF not assigned");
     }

   return 0;
}

/*}}}*/

/* Define W(h) = t \int dy R(h,y)A(y)
 *
 * In the linear regime, the flux-corrected data is defined to be
 *
 *        f(h) = (C(h) - B(h)) / W(h)
 *     df^2(h) = (C(h) + B(h)) / W^2(h)
 *
 * To rebin, summing over several bins h, to make a big bin H,
 * the correct result is:
 *
 *             \sum_h C(h) - \sum_h B(h)
 *      f(H) = -------------------------
 *                  \sum_h W(h)
 *
 *             \sum_h C(h) + \sum_h B(h)
 *   df^2(H) = -------------------------
 *                 (\sum_h W(h))^2
 *
 * If we define weights
 *
 *               W(h)
 *    w(h) = -------------
 *           \sum_h' W(h')
 *
 * we can rewrite f(H) and df(H) as:
 *
 *      f(H) = \sum_h w(h) f(h)
 *   df^2(H) = \sum_h w^2(h) df^2(h)
 *
 * This approach allows one to perform arbitrary rebinning
 * after computing the f(h), df(h) and W(h) values only once.
 */
int Hist_flux_correct (Hist_t *h, double frac, double *bkg, int using_model, /*{{{*/
                       double *kernel_params, unsigned int num_kernel_params,
                       char *options)
{
   Isis_Hist_t cts;
   Isis_Kernel_t *k;
   double *wt = NULL;
   double *f, *f_err, *rf, *rf_err;
   int status = -1;

   if ((NULL == h) || (NULL == h->kernel))
     return -1;

   k = h->kernel;

   if (-1 == verify_flux_corr_kernel (k))
     return -1;

   if (-1 == setup_flux_correct (h, using_model, &cts, &f, &f_err))
     return -1;

   /* f(h) = (D(h)-B(h)) / t \int dy R(h,y) A(y)   [photons/s/cm^2/A]  */
   if (-1 == k->compute_flux (k, kernel_params, num_kernel_params,
                              &cts, bkg, f, f_err, &wt, options))
     goto return_status;

   if (wt == NULL)
     {
        if (NULL == (wt = isis_unit_source_response (k)))
          goto return_status;
     }

   ISIS_FREE (h->flux_weights);
   h->flux_weights = wt;

   if (using_model)
     {
        rf = h->convolved_model_flux;
        rf_err = NULL;
     }
   else
     {
        rf = h->flux;
        rf_err = h->flux_err;
     }

   if (-1 == rebin_flux_using_weights (h, f, f_err, rf, rf_err))
     goto return_status;

   if (-1 == set_flux_floor (h, frac, rf))
     goto return_status;

   if (using_model == 0)
     {
        int valid = validate_flux_err (h);
        if (valid == -1)
          goto return_status;
        else if (valid == 1)
          invalid_uncertainties_replaced ();
     }

   status = 0;
   return_status:

   if (status)
     isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "flux correct failed");

   if (using_model)
     {
        ISIS_FREE (f);
        ISIS_FREE (f_err);
        ISIS_FREE (cts.val_err);
     }

   return status;
}

/*}}}*/

int Hist_flux_corr_model (Hist_t *h, double *bincts) /*{{{*/
{
   double *wt_sum = NULL, *lo, *hi;
   int i, n_notice, *notice_list;

   if (h->flux_weights == NULL)
     {
        if (NULL == (h->flux_weights = isis_unit_source_response (h->kernel)))
          return -1;
     }

   if (NULL == (wt_sum = (double *) ISIS_MALLOC (h->nbins * sizeof(double))))
     return -1;
   memset ((char *)wt_sum, 0, h->nbins * sizeof(double));

   if (-1 == Hist_apply_rebin_and_notice_list (wt_sum, h->flux_weights, h))
     {
        ISIS_FREE(wt_sum);
        return -1;
     }

   lo = h->bin_lo;
   hi = h->bin_hi;
   n_notice = h->n_notice;
   notice_list = h->notice_list;

   /* [bincts] = counts
    * [wt_sum] = sec Angstrom (cm^2 counts/photon)
    * [bincts/wt_sum] = photons/sec/cm^2/A
    * so multiply by bin width to get bin-integrated photon flux,
    *   photons/sec/cm^2/bin
    *
    * We're operating on the data in place for efficiency,
    * so even though the array name is no longer mnemonic, the
    * final result has flux units.
    */

   for (i = 0; i < n_notice; i++)
     {
        int k = notice_list[i];
        if (wt_sum[i] != 0.0)
          {
             bincts[i] *= (hi[k] - lo[k]) / wt_sum[i];
          }
     }

   ISIS_FREE(wt_sum);

   return 0;
}

/*}}}*/
