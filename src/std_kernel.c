/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2010  Massachusetts Institute of Technology

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

/* $Id: std_kernel.c,v 1.3 2004/02/09 11:14:25 houck Exp $ */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#define ISIS_KERNEL_PRIVATE_DATA \
   int allows_ignoring_model_intervals;

#include "isis.h"
#include "util.h"
#include "errors.h"

static void delete_kernel (Isis_Kernel_t *k) /*{{{*/
{
   ISIS_FREE (k);
}

/*}}}*/

static int print_kernel (Isis_Kernel_t *k) /*{{{*/
{
   (void) k;
   fprintf (stderr, "function not implemented\n");
   return 0;
}

/*}}}*/

static int match_arf_grid (Isis_Kernel_t *k, Isis_Rsp_t *rsp, Isis_Hist_t *g, Isis_Hist_t *m) /*{{{*/
{
   Isis_Arf_t *a;
   double *g_val, *m_val;
   int i;

   g_val = m_val = NULL;

   /* When there's only one ARF, the model is computed on that grid.
    * If we've got multiple responses, we interpolate the model
    * onto the current ARF grid
    */

   if (k->rsp.next == NULL)
     {  /* easy case:  The model is on the ARF grid */
        m->val = g->val;
        m->notice_list = g->notice_list;
        m->n_notice = g->n_notice;
        return 0;
     }

   /* harder case: map the model onto the ARF grid */
   a = rsp->arf;
   m->nbins = a->nbins;
   m->bin_lo = a->bin_lo;
   m->bin_hi = a->bin_hi;
   m->n_notice = 0;
   m->val = NULL;
   m->notice = NULL;
   m->notice_list = NULL;
   if ((NULL == (m_val = (double *) ISIS_MALLOC (m->nbins * sizeof(double))))
       || (NULL == (m->notice = (int *) ISIS_MALLOC (m->nbins * sizeof(int))))
       || (NULL == (g_val = (double *) ISIS_MALLOC (g->nbins * sizeof(double)))))
     {
        goto fail;
     }
   memset ((char *)m->notice, 0, m->nbins * sizeof(int));
   memset ((char *)g_val, 0, g->nbins * sizeof(double));

   if (-1 == unpack_noticed (g->val, g->notice_list, g->n_notice, g->nbins, g_val))
     goto fail;

   if ((-1 == rebin_histogram (g_val, g->bin_lo, g->bin_hi, g->nbins,
                                m_val, m->bin_lo, m->bin_hi, m->nbins))
       || (-1 == transfer_notice (g->bin_lo, g->bin_hi, g->notice_list, g->n_notice,
                                  m->bin_lo, m->bin_hi, m->nbins, m->notice))
       || (-1 == _update_notice_list (m->notice, &m->notice_list, &m->n_notice, m->nbins)))
     {
        goto fail;
     }

   /* Finish by packing according to the notice list */
   if (NULL == (m->val = (double *) ISIS_MALLOC (m->n_notice * sizeof(double))))
     goto fail;

   for (i = 0; i < m->n_notice; i++)
     {
        int n = m->notice_list[i];
        m->val[i] = m_val[n];
     }

   ISIS_FREE (g_val);
   ISIS_FREE (m_val);
   ISIS_FREE (m->notice);
   return 1;

   fail:
   ISIS_FREE (g_val);
   ISIS_FREE (m_val);
   ISIS_FREE (m->notice);
   ISIS_FREE (m->val);
   ISIS_FREE (m->notice_list);
   m->n_notice = 0;
   return -1;
}

/*}}}*/

static int compute_kernel (Isis_Kernel_t *k, double *result, Isis_Hist_t *g, double *par, unsigned int num, /*{{{*/
                           int (*fun)(Isis_Hist_t *))
{
   Isis_Rsp_t *rsp;
   int ret = -1;

   (void) par; (void) num;

   if ((k == NULL) || (g == NULL) || (NULL == fun))
     return -1;

   /* Evaluate the model once */
   if (-1 == (*fun)(g))
     return -1;

   /* Fold the model through each of the responses,
    * incrementing 'result' for each such contribution
    */
   for (rsp = &k->rsp; rsp != NULL; rsp = rsp->next)
     {
        double *arf = rsp->arf->arf;
        Isis_Hist_t m;
        int i, malloced;

        if ((malloced = match_arf_grid (k, rsp, g, &m)) < 0)
          return -1;

        for (i = 0; i < m.n_notice; i++)
          {
             int n = m.notice_list[i];
             m.val[i] *= (arf[n] * k->exposure_time);
          }

        ret = k->apply_rmf (rsp->rmf, result, k->num_orig_data,
                            m.val, m.notice_list, m.n_notice);
        if (malloced)
          {
             ISIS_FREE (m.notice_list);
             ISIS_FREE (m.val);
          }

        if (ret == -1)
          return ret;
     }

   return ret;
}

/*}}}*/

static int compute_flux (Isis_Kernel_t *k, double *kernel_params, unsigned int num_kernel_params, /*{{{*/
                         Isis_Hist_t *counts, double *bgd,
                         double *f, double *df, double **weights, char *options)
{
   double *cts, *ar;
   int i, n;

   (void) kernel_params; (void) num_kernel_params; (void) options;

   if (k == NULL)
     return -1;

   if (k->rsp.next != NULL)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "flux-correction for multiple responses is not supported");
        return -1;
     }

   n = k->num_orig_data;
   cts = counts->val;

   /* Do background subtraction */

   if (bgd == NULL)
     {
        for (i = 0; i < n; i++)
          {
             f[i] = cts[i];
             /* zero is ok */
             df[i] = sqrt(fabs(cts[i]));
          }
     }
   else
     {
        for (i = 0; i < n; i++)
          {
             f[i] = cts[i] - bgd[i];
             /* zero is ok */
             df[i] = sqrt(fabs(cts[i] + bgd[i]));
          }
     }

   /* flux-correct */

   if (NULL == (ar = isis_unit_source_response (k)))
     return -1;

   for (i = 0; i < n ; i++)
     {
        if (ar[i] != 0.0)
          {
             f[i] /= ar[i];
             df[i] /= ar[i];
          }
     }

   *weights = ar;

   return 0;
}

/*}}}*/

static int allows_ignoring_model_intervals (Isis_Kernel_t *k) /*{{{*/
{
   if (k == NULL)
     return 0;
   return k->allows_ignoring_model_intervals;
}

/*}}}*/

static int eval_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Kernel_t *k = (Isis_Kernel_t *)clientdata;

   (void) subsystem;
   (void) optname;

   if (k == NULL)
     return -1;

   if (0 == isis_strcasecmp (value, "all"))
     k->allows_ignoring_model_intervals = 0;
   else if (0 == isis_strcasecmp (value, "noticed"))
     k->allows_ignoring_model_intervals = 1;
   else
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "unrecognized kernel option '%s'",
                    value ? value : "<null>");
        return -1;
     }

   return 0;
}

/*}}}*/

static Isis_Option_Table_Type Std_Option_Table [] = /*{{{*/
{
     {"eval", eval_option, ISIS_OPT_REQUIRES_VALUE, "noticed", "specify energies to evaluate model: (all | noticed)"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

/*}}}*/

static int process_options (Isis_Kernel_t *k, char *options) /*{{{*/
{
   Isis_Option_Type *o;
   int status;

   if (options == NULL)
     return 0;

   if (NULL == (o = isis_parse_option_string (options)))
     return -1;

   status = isis_process_options (o, Std_Option_Table, (void *)k, 1);
   isis_free_options (o);

   return status;
}

/*}}}*/

static Isis_Kernel_t *allocate_std_kernel (Isis_Obs_t *o, char *options) /*{{{*/
{
   Isis_Kernel_t *k = NULL;

   (void) options;

   if (NULL == (k = isis_init_kernel (NULL, sizeof(*k), o)))
     return NULL;

   k->allows_ignoring_model_intervals = 1;

   if (-1 == process_options (k, options))
     {
        delete_kernel (k);
        return NULL;
     }

   k->delete_kernel = delete_kernel;
   k->compute_kernel = compute_kernel;
   k->compute_flux = compute_flux;
   k->print_kernel = print_kernel;

   return k;
}

/*}}}*/

ISIS_USER_KERNEL_MODULE(std,def,options)
{
   (void) options;

   def->kernel_name = "std";
   def->allocate_kernel = allocate_std_kernel;
   def->allows_ignoring_model_intervals = allows_ignoring_model_intervals;
   def->num_kernel_parms = 0;
   def->kernel_parm_names = NULL;
   def->default_min = NULL;
   def->default_max = NULL;
   def->default_value = NULL;
   def->default_freeze = NULL;

   return 0;
}

static int compute_yshift_kernel (Isis_Kernel_t *k, double *result, Isis_Hist_t *g, double *par, unsigned int num, /*{{{*/
                                 int (*fun)(Isis_Hist_t *))
{
   Isis_Rmf_t *rmf = k->rsp.rmf;
   double *ylo=NULL, *yhi=NULL;
   double dy = par[0];
   unsigned int n;
   int status = -1;

   if (-1 == compute_kernel (k, result, g, par, num, fun))
     return -1;

   if (dy == 0.0)
     return 0;

   if (-1 == rmf->get_data_grid (rmf, &ylo, &yhi, &n, NULL))
     return -1;

   if (ylo[0] + dy > 0.0)
     {
        double *tmp, *shift_lo, *shift_hi;
        unsigned int i, len;

        len = 3 * n * sizeof(double);
        if (NULL == (tmp = (double *) ISIS_MALLOC (len)))
          goto return_error;
        shift_lo  = tmp + n;
        shift_hi  = tmp + 2*n;

        /* dy > 0 moves features to longer wavelengths */

        shift_lo[0] = ylo[0] - dy;
        for (i = 1; i < n; i++)
          {
             shift_lo[i] = ylo[i] - dy;
             shift_hi[i-1] = shift_lo[i];
          }
        shift_hi[n-1] = yhi[n-1] - dy;

        if (-1 == rebin_histogram (result, ylo, yhi, n,
                                   tmp, shift_lo, shift_hi, n))
          {
             ISIS_FREE(tmp);
             isis_vmesg(FAIL, I_ERROR, __FILE__, __LINE__,
                        "shift kernel failed while rebinning histogram");
             goto return_error;
          }

        memcpy ((char *)result, (char *)tmp, n * sizeof(double));
        ISIS_FREE(tmp);
     }
   else
     {
        isis_vmesg(FAIL, I_ERROR, __FILE__, __LINE__, "offset=%g yields invalid grid", dy);
        goto return_error;
     }

   status = 0;
   return_error:
   ISIS_FREE(ylo);
   ISIS_FREE(yhi);

   return status;
}

/*}}}*/

static Isis_Kernel_t *allocate_yshift_kernel (Isis_Obs_t *o, char *options) /*{{{*/
{
   Isis_Kernel_t *k = NULL;

   (void) options;

   if (NULL == (k = isis_init_kernel (NULL, sizeof(*k), o)))
     return NULL;

   k->allows_ignoring_model_intervals = 0;

   k->delete_kernel = delete_kernel;
   k->compute_kernel = compute_yshift_kernel;
   k->compute_flux = NULL;
   k->print_kernel = print_kernel;

   return k;
}

/*}}}*/

ISIS_USER_KERNEL_MODULE(yshift,def,options)
{
   static char *parm_names[] = {"offset", NULL};
   static char *parm_units[] = {"A", NULL};
   static double default_min [] = {0.0};
   static double default_max [] = {0.1};
   static double default_value [] = {0.0};
   static unsigned int default_freeze [] = {1};
   (void) options;

   def->kernel_name = "yshift";
   def->allocate_kernel = allocate_yshift_kernel;
   def->allows_ignoring_model_intervals = NULL;
   def->num_kernel_parms = 1;
   def->kernel_parm_names = parm_names;
   def->kernel_parm_units = parm_units;
   def->default_min = default_min;
   def->default_max = default_max;
   def->default_value = default_value;
   def->default_freeze = default_freeze;

   return 0;
}

static int compute_gainshift_kernel (Isis_Kernel_t *k, double *result, Isis_Hist_t *g, double *par, unsigned int num, /*{{{*/
                                     int (*fun)(Isis_Hist_t *))
{
   Isis_Rmf_t *rmf = k->rsp.rmf;
   double *ylo=NULL, *yhi=NULL;
   double *tmp, *shift_lo, *shift_hi;
   double r0 = par[0]/KEV_ANGSTROM, slope = par[1];
   unsigned int i, len, n;
   int status = -1;

   if (-1 == compute_kernel (k, result, g, par, num, fun))
     return -1;

   if ((par[0] < 0) || (par[1] == 0.0))
     {
        isis_vmesg(FAIL, I_ERROR, __FILE__, __LINE__,
                   "gainshift kernel:  parameters (%g, %g) define an invalid grid",
                   par[0], par[1]);
        return -1;
     }

   if (-1 == rmf->get_data_grid (rmf, &ylo, &yhi, &n, NULL))
     return -1;

   len = 3 * n * sizeof(double);
   if (NULL == (tmp = (double *) ISIS_MALLOC (len)))
     goto return_error;
   shift_lo  = tmp + n;
   shift_hi  = tmp + 2*n;

#define NEW_LAMBDA(y)    (1.0/(1.0/y/slope - r0))

   shift_lo[0] = NEW_LAMBDA(ylo[0]);
   for (i = 1; i < n; i++)
     {
        shift_lo[i] = NEW_LAMBDA(ylo[i]);
        shift_hi[i-1] = shift_lo[i];
     }
   shift_hi[n-1] = NEW_LAMBDA(yhi[n-1]);

   if (-1 == rebin_histogram (result, ylo, yhi, n,
                              tmp, shift_lo, shift_hi, n))
     {
        ISIS_FREE(tmp);
        isis_vmesg(FAIL, I_ERROR, __FILE__, __LINE__,
                   "gainshift kernel failed while rebinning histogram");
        goto return_error;
     }

   memcpy ((char *)result, (char *)tmp, n * sizeof(double));
   ISIS_FREE(tmp);

   status = 0;
   return_error:
   ISIS_FREE(ylo);
   ISIS_FREE(yhi);

   return status;
}

/*}}}*/

static Isis_Kernel_t *allocate_gainshift_kernel (Isis_Obs_t *o, char *options) /*{{{*/
{
   Isis_Kernel_t *k = NULL;

   (void) options;

   if (NULL == (k = isis_init_kernel (NULL, sizeof(*k), o)))
     return NULL;

   k->allows_ignoring_model_intervals = 0;

   k->delete_kernel = delete_kernel;
   k->compute_kernel = compute_gainshift_kernel;
   k->compute_flux = NULL;
   k->print_kernel = print_kernel;

   return k;
}

/*}}}*/

ISIS_USER_KERNEL_MODULE(gainshift,def,options)
{
   static char *parm_names[] = {"intercept", "slope", NULL};
   static char *parm_units[] = {"eV", "", NULL};
   static double default_min [] = {-1.0, 0.8};
   static double default_max [] = { 1.0, 1.2};
   static double default_value [] = {0, 1.0};
   static unsigned int default_freeze [] = {1, 0};
   (void) options;

   def->kernel_name = "gainshift";
   def->allocate_kernel = allocate_gainshift_kernel;
   def->allows_ignoring_model_intervals = NULL;
   def->num_kernel_parms = 2;
   def->kernel_parm_names = parm_names;
   def->kernel_parm_units = parm_units;
   def->default_min = default_min;
   def->default_max = default_max;
   def->default_value = default_value;
   def->default_freeze = default_freeze;

   return 0;
}

