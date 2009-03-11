#ifndef ISIS_HISTOGRAM_H
#define ISIS_HISTOGRAM_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008  Massachusetts Institute of Technology

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

/* $Id: histogram.h,v 1.21 2004/05/07 17:22:42 houck Exp $ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include "cfits.h"

#define NOT_FITS_FORMAT (-2)

/* Histogram version bit-masks:
 *
 *    DATA    bit indicates "data" or "model"
 *    FLUX    bit indicates "flux" or "counts"
 */

#define H_CONVOLVED ((unsigned int) 1<<0)
#define H_DATA      ((unsigned int) 1<<1)
#define H_FLUX      ((unsigned int) 1<<2)

#define is_data(v)       ((v) & H_DATA)
#define is_model(v)    (!((v) & H_DATA))
#define is_flux(v)       ((v) & H_FLUX)
#define is_counts(v)   (!((v) & H_FLUX))
#define is_convolved(v)  ((v) & H_CONVOLVED)

#define USE_ASSIGNED_ARF_RMF 0x0
#define USE_IDEAL_ARF        0x1
#define USE_IDEAL_RMF        0x2
#define USE_IDEAL_ARF_RMF    (USE_IDEAL_ARF | USE_IDEAL_RMF)

#define want_ideal_arf(v) ((v) & USE_IDEAL_ARF)
#define want_ideal_rmf(v) ((v) & USE_IDEAL_RMF)

extern double Hist_Rmf_Grid_Match_Tol;
extern double Hist_Min_Stat_Err;
extern double Hist_Gehrels_Threshold;
extern double Hist_Min_Model_Spacing;
extern int Hist_Ignore_PHA_Response_Keywords;
extern int Hist_Ignore_PHA_Backfile_Keyword;

typedef int Hist_Stat_Error_Hook_Type (void *, double *, double *, int *, int,
                                       double *, int);
extern Hist_Stat_Error_Hook_Type *Hist_Stat_Error_Hook;

typedef struct _Hist_t Hist_t;
typedef struct _Hist_Stat_t Hist_Stat_t;

struct _Hist_Stat_t
{
   double min, max;                  /* region boundaries */
   double sum, sum_err;
   double net, net_err;              /* continuum subtracted */
   double centroid, centroid_err;    /* continuum subtracted */
   double eqwidth, eqwidth_err;      /* equivalent width */
   double contin, slope;             /* continuum density, slope at centroid position */
   int nbins;
};

typedef struct
{
   double tstart;
   double combo_weight;
   char object[CFLEN_VALUE];
   char instrument[CFLEN_VALUE];
   char grating[CFLEN_VALUE];
   char *file;
   char *bgd_file;
   int spec_num;
   int order;
   int part;
   int srcid;
   int exclude;
   int combo_id;
}
Hist_Info_Type;

/* mostly boolean flags */
typedef struct
{
   int *use_style;       /* override current line style? */
   int overlay;          /* this plot is an overlay? */
   int residuals;        /* plot residuals? */
   int default_labels;   /* use default labels? */
   int rate;             /* divide by the exposure time? */
}
Hist_Plot_Tune_Type;
extern Hist_t *Hist_Current;
extern Hist_t *Hist_init_list (void);

/* Things that use the list head pointer */

/* List handling */
extern void Hist_free_list (Hist_t *head);
extern unsigned int Hist_num_noticed_histograms (Hist_t *head);
extern int Hist_all_data_noticed (Hist_t *head);
extern int Hist_all_model_noticed (Hist_t *head);
extern int Hist_delete_hist (Hist_t *head, int hist_index);
extern Hist_t *_Hist_find_hist_index (Hist_t * head, int hist_index);
extern Hist_t *Hist_find_hist_index (Hist_t *head, int hist_index);
extern int Hist_id_list (Hist_t *head, int noticed, unsigned int **ids, unsigned int *num);
extern int Hist_map (Hist_t *head, int (*fun)(Hist_t *, void *), void *cl, int check_exclude);

/* data I/O */
extern int Hist_read_ascii (Hist_t *head, char *filename);
extern int Hist_read_fits (Hist_t *head, Isis_Arf_t *arf_head, Isis_Rmf_t *rmf_head, char *pha_filename,
                           int **indices, int *num_spectra, int strict, int just_one);
extern int Hist_define_data (Hist_t *head, Isis_Hist_t *x, unsigned int bin_type, unsigned int has_grid);
extern int Hist_assign_arf_to_hist (Hist_t *h_head, int hist_index, Isis_Arf_t *r_head, int arf_index);
extern int Hist_assign_rmf_to_hist (Hist_t *h_head, int hist_index, Isis_Rmf_t *r_head, int rmf_index);
extern int Hist_get_rsp_list (Hist_t *h, int **arfs, int **rmfs, int *num);
extern int Hist_assign_rsp_list (Hist_t *hhead, int hist_index,
                                 Isis_Arf_t *ahead, int *arfs,
                                 Isis_Rmf_t *rhead, int *rmfs, int num);

/* dataset combinations */
extern int Hist_any_noticed_dataset_combinations (Hist_t *head);
extern int Hist_grids_are_identical (Hist_t *a, Hist_t *b);
extern int Hist_break_combination (Hist_t *head, unsigned int gid);
extern int Hist_dataset_is_combined (Hist_t *h, unsigned int gid);
extern int Hist_combination_id (Hist_t *h);
extern int Hist_combination_weight (Hist_t *h, double *weight);
extern int Hist_set_combination (Hist_t *head, unsigned int *members,
                                 double *weights, unsigned int num_members);

/* model evaluation grid */
typedef struct
{
   int (*make_grid)(Hist_t *, void *);   
   int (*eval_model)(Hist_t *, Isis_Hist_t *);
   void *options;
   void (*destroy_options)(void *);
   int type;
   int id;
   int cache_model_values;
   /* boolean:  non-zero means one grid shared between datasets
    *           with the same type/id */   
}
Hist_Eval_Grid_Method_Type;

extern Hist_Eval_Grid_Method_Type *Hist_eval_grid_method (Hist_t *h);
extern int Hist_set_eval_grid_info (Hist_t *head, unsigned int *indices, unsigned int num_indices,
                                    Hist_Eval_Grid_Method_Type *m);
extern int Hist_make_merged_eval_grid (Hist_t *head, int eval_grid_id, Isis_Hist_t *grid);
extern int Hist_finalize_user_eval_grid (Hist_t *h, Isis_Hist_t *g);
extern int Hist_make_user_eval_grid (Hist_t *, Isis_Hist_t *m);

/* Things that use 'version' */
extern int Hist_hist_size (Hist_t *h, unsigned int version);
extern int Hist_get_hist_grid (Hist_t *h, unsigned int version, Isis_Hist_t *g);
extern int Hist_copy_hist_data (Hist_t *h, unsigned int version, double *val, double *val_err);
extern int Hist_replace_hist_data (Hist_t *h, unsigned int version, double *val, double *val_err, int nbins);
extern int Hist_copy_noticed_data (Hist_t *h, unsigned int version, double *val, double *val_err);
extern int Hist_region_stats (Hist_t *h, unsigned int version,
                              Hist_Stat_t *s, double x[/*2*/], double y[/*2*/]);
extern int Hist_copy_packed_model (Hist_t *h, unsigned int version, double *val);
extern int Hist_set_model (Hist_t *h, unsigned int version, double *model_value);
#ifdef ISIS_PLOT_H
extern int Hist_plot (Hist_t *h, unsigned int version, Hist_Plot_Tune_Type *info, /*{{{*/
                      Plot_t *fmt, Isis_Fit_Statistic_Type *s);
#endif

/* Things that operate on a specific histogram (without version) */

/* background */
extern int Hist_set_background_from_file (Hist_t *h, char *file);
extern int Hist_set_background_name (Hist_t *h, char *name);
extern int Hist_copy_scaled_background (Hist_t *h, double **bgd);
extern int Hist_set_instrumental_background_hook_name (Hist_t *h, char *hook_name);
extern int Hist_set_instrumental_background_hook (Hist_t *h, void *hook);
extern char *Hist_get_instrumental_background_hook_name (Hist_t *h);
extern void *Hist_get_instrumental_background_hook(Hist_t *h);
extern int Hist_define_background (Hist_t *h, double bgd_exposure, 
                                   double *bgd_area, int area_is_vector,
                                   double *bgd, unsigned int nbins);

/* systematic errors */
extern int Hist_define_sys_err_frac (Hist_t *h, double *sys_err_frac, int nbins);
extern int Hist_copy_sys_err_frac (Hist_t *h, double **sys_err_frac, int *nbins);

/* generic access */
extern int Hist_get_info (Hist_t *h, Hist_Info_Type *info);
extern int Hist_set_info (Hist_t *h, Hist_Info_Type *info);
extern int Hist_set_frame_time (Hist_t *h, double frame_time);
extern int Hist_get_frame_time (Hist_t *h, double *frame_time);
extern int Hist_get_exposure (Hist_t *h, double *exposure);
extern int Hist_set_exposure (Hist_t *h, double exposure);
extern int Hist_get_data_region_area (Hist_t *h, double **region_area, int *num);
extern int Hist_get_back_region_area (Hist_t *h, double **region_area, int *num);
extern int Hist_set_data_region_area (Hist_t *h, double *region_area, int num);
extern int Hist_set_back_region_area (Hist_t *h, double *region_area, int num);
extern int Hist_get_back_exposure (Hist_t *h, double *back_exposure);
extern int Hist_set_color (Hist_t *h, int color);
extern int Hist_unset_color (Hist_t *h);
extern int Hist_set_object_name (Hist_t *h, char *object);
extern int Hist_get_flux_corr_weights (Hist_t *h, double **weights, int *num_weights);

/* misc.. */
extern int Hist_is_fake (Hist_t *h);
extern int Hist_set_fake (Hist_t *h, int value);
extern int Hist_print_stats (FILE *fp, Hist_Stat_t *s);
extern int Hist_get_index (Hist_t *h);
extern int Hist_orig_hist_size (Hist_t *h);
extern int Hist_init_model_structs (Hist_t *h);
extern int Hist_get_model_grid (Isis_Hist_t *g, Hist_t *h);
extern int _Hist_get_orig_hist_grid (Hist_t *h, Isis_Hist_t *g);
extern int Hist_replace_hist_grid (Hist_t *h, double *bin_lo, double *bin_hi, int nbins);
extern int Hist_copy_histogram_keywords (Hist_t *dst, Hist_t *src);

/* responses */
extern int Hist_set_fit_responses
  (Hist_t *h, unsigned int fit_data_type, unsigned int response_mask);

/* ignore/notice */
extern int Hist_num_data_noticed (Hist_t *h);
extern int Hist_ignore_bad (Hist_t *h);
extern int Hist_set_notice (Hist_t *h, int value, double bin_lo, double bin_hi);
extern int Hist_set_notice_using_mask (Hist_t *h, int *mask, int nbins);
extern int Hist_set_notice_using_list (Hist_t *h, int value, unsigned int *list, unsigned int n);
extern int Hist_set_exclude_flag (Hist_t *h, int exclude);

/* rebin */
extern int Hist_rebin_min_counts (Hist_t *h, void *s);
extern int Hist_rebin_index (Hist_t *h, void *s);
extern int Hist_get_hist_rebin_info (Hist_t *h, int **rebin, int *orig_nbins);
extern int Hist_do_rebin (Hist_t *h, int (*rebin_fcn)(Hist_t *, void *), void *s);
extern int Hist_apply_rebin_and_notice_list (double *bin_and_notice_result, double *x, Hist_t *h);
extern int Hist_apply_rebin (double *x, Hist_t *h, double **rebinned, int *nbins);
extern int Hist_rebin (Hist_t *h, double *lo, double *hi, int nbins);
extern int Hist_set_stat_error_hook (Hist_t *h, void *hook, void (*delete_hook)(void *));

/* flux correct */
extern int Hist_flux_correct (Hist_t *h, double threshold, double *bgd, int method,
                              double *kernel_params, unsigned int num_kernel_params,
                              char *options);
/* Kernel handling */
extern int Hist_load_aux_kernel (char *libfile, char *init_name, char *init_args, Isis_Kernel_Def_t *def);
extern int Hist_reallocate_kernel (Hist_t *h, char *params, Isis_Kernel_Def_t *def);
extern int Hist_print_kernel (Hist_t *h);
extern Isis_Kernel_t *Hist_get_kernel (Hist_t *h);
extern char *Hist_get_kernel_params (Hist_t *h);
extern Isis_Kernel_t *Hist_index_get_kernel (Hist_t *h);
extern int Hist_compute_flux (Isis_Kernel_t *k, Isis_Hist_t *counts, double threshold,
                              Isis_Hist_t *flux);

/* post model hook */
extern void *Hist_post_model_hook (Hist_t *h);
extern int Hist_set_post_model_hook (Hist_t *h, void *hook, void (*delete_hook)(void *));
extern int Hist_run_post_model_hook (Hist_t *h, double *cts, SLang_Array_Type *sl_bgd);

/* user-defined metadata */
extern void *Hist_get_metadata (Hist_t *h);
extern int Hist_set_metadata (Hist_t *h, void *meta);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
