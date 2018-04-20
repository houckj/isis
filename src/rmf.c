/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2018  Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Authors:  John C. Houck  <houck@space.mit.edu>

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

/* $Id: rmf.c,v 1.18 2004/06/06 20:56:07 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "util.h"
#include "rmf.h"
#include "errors.h"

/*}}}*/

/* new/free */

void Isis_free_rmf_grid (Isis_Rmf_Grid_Type *eb) /*{{{*/
{
   if (NULL == eb)
     return;

   ISIS_FREE (eb->bin_lo);
   ISIS_FREE (eb->bin_hi);
   ISIS_FREE (eb);
}
/*}}}*/

Isis_Rmf_Grid_Type *Isis_new_rmf_grid (unsigned int nbins, double *lo, double *hi) /*{{{*/
{
   Isis_Rmf_Grid_Type *g;

   if (NULL == (g = (Isis_Rmf_Grid_Type *) ISIS_MALLOC (sizeof(Isis_Rmf_Grid_Type))))
     return NULL;
   memset ((char *)g, 0, sizeof (*g));

   g->nbins = nbins;
   g->units = -1;

   if (NULL == (g->bin_lo = (double *) ISIS_MALLOC (nbins * sizeof(double)))
       || NULL == (g->bin_hi = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     {
        Isis_free_rmf_grid (g);
        g = NULL;
        return g;
     }
   if (lo == NULL)
     memset ((char *)g->bin_lo, 0, nbins*sizeof(double));
   else
     memcpy ((char *)g->bin_lo, (char *)lo, nbins * sizeof(double));

   if (hi == NULL)
     memset ((char *)g->bin_hi, 0, nbins*sizeof(double));
   else
     memcpy ((char *)g->bin_hi, (char *)hi, nbins * sizeof(double));

   return g;
}
/*}}}*/

void Rmf_free_rmf (Isis_Rmf_t *rmf) /*{{{*/
{
   if (rmf == NULL)
     return;

   /* Free memory only if nobody is pointing to it */

   if (rmf->ref_count > 0)
     return;

   if (rmf->delete_client_data)
     rmf->delete_client_data (rmf);

   ISIS_FREE (rmf->arg_string);
   ISIS_FREE (rmf);
}

/*}}}*/

static int default_set_noticed_model_bins (Isis_Rmf_t *rmf, int num_chan, int *chan_notice, /*{{{*/
                                           int num_model, int *model_notice)
{
   int k;

   (void) num_chan; (void) chan_notice;

   if (rmf == NULL || model_notice == NULL)
     return -1;

   for (k = 0; k < num_model; k++)
     {
        model_notice[k] = 1;
     }

   return 0;
}

/*}}}*/

static int default_rebin_rmf (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int num) /*{{{*/
{
   (void) rmf; (void) lo; (void) hi; (void) num;

   isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__,
               "rebinning is not supported for this RMF type");
   return -1;
}

/*}}}*/

static int default_factor_rsp (Isis_Rmf_t *rmf, double *arf) /*{{{*/
{
   (void) rmf; (void) arf;

   isis_vmesg (FAIL, I_NOT_IMPLEMENTED, __FILE__, __LINE__,
               "factorization is not supported for this RMF type");
   return -1;
}

/*}}}*/

static Isis_Rmf_t *new_rmf (void) /*{{{*/
{
   Isis_Rmf_t *rmf;

   if (NULL == (rmf = (Isis_Rmf_t *) ISIS_MALLOC (sizeof(Isis_Rmf_t))))
     return NULL;
   memset ((char *)rmf, 0, sizeof(*rmf));

   rmf->next = NULL;
   rmf->index = 0;
   rmf->ref_count = 0;

   rmf->order = 0;
   rmf->grating[0] = 0;
   rmf->instrument[0] = 0;
   rmf->arg_string = NULL;
   rmf->method = RMF_DELTA;
   rmf->includes_effective_area = 0;

   rmf->set_arf_grid = NULL;
   rmf->get_arf_grid = NULL;
   rmf->set_data_grid = NULL;
   rmf->get_data_grid = NULL;
   rmf->redistribute = NULL;
   rmf->delete_client_data = NULL;

   rmf->set_noticed_model_bins = default_set_noticed_model_bins;
   rmf->rebin_rmf = default_rebin_rmf;
   rmf->factor_rsp = default_factor_rsp;

   return rmf;
}

/*}}}*/

/* linked list maintenance */

static int rmf_list_append (Isis_Rmf_t *head, Isis_Rmf_t *rmf) /*{{{*/
{
  Isis_Rmf_t *a;

   if (rmf == NULL)
     return -1;

  for (a = head; a->next != NULL; a = a->next)
     ;

  a->next = rmf;
  rmf->next = NULL;
  rmf->index = a->index + 1;

  return rmf->index;
}

/*}}}*/

static int delete_after_rmf (Isis_Rmf_t *a) /*{{{*/
{
   Isis_Rmf_t *dead;
   Isis_Rmf_t *next;

   if (a == NULL)
     return -1;

   dead = a->next;

   if (dead->ref_count > 0)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "RMF in use");
        return 0;
     }

   isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "Deleting RMF %d", dead->index);

   next = dead->next;
   Rmf_free_rmf (dead);
   a->next = next;

  return 0;
}

/*}}}*/

Isis_Rmf_t *Rmf_find_rmf_index (Isis_Rmf_t *head, int rmf_index) /*{{{*/
{
   Isis_Rmf_t *a;

   if (head == NULL)
     return NULL;

   for (a = head->next; a != NULL; a = a->next)
     {
        if (a->index == rmf_index)
          return a;
     }

   isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "RMF %d not found", rmf_index);
   return NULL;
}

/*}}}*/

void Rmf_free_rmf_list (Isis_Rmf_t *head) /*{{{*/
{
   while (head != NULL)
     {
        Isis_Rmf_t *n = head->next;
        Rmf_free_rmf (head);
        head = n;
     }
}

/*}}}*/

int Rmf_delete_rmf (Isis_Rmf_t *head, int rmf_index) /*{{{*/
{
   Isis_Rmf_t *a;

   if (head == NULL)
     return -1;

   for (a = head; a->next != NULL; a = a->next)
     {
        if (a->next->index == rmf_index)
          return delete_after_rmf (a);
     }

   isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "RMF %d not found", rmf_index);
   return -1;
}

/*}}}*/

Isis_Rmf_t *Rmf_init_rmf_list (void) /*{{{*/
{
   return new_rmf ();
}

/*}}}*/

void Rmf_free_info (Rmf_Info_Type *info) /*{{{*/
{
   if (info == NULL)
     return;

   ISIS_FREE(info->grating);
   ISIS_FREE(info->instrument);
   ISIS_FREE(info->arg_string);
}

/*}}}*/

int Rmf_get_info (Isis_Rmf_t *rmf, Rmf_Info_Type *info) /*{{{*/
{
   if ((rmf == NULL) || (info == NULL))
     return -1;

   info->index = rmf->index;
   info->method = rmf->method;
   info->order = rmf->order;

   if (NULL == (info->grating = isis_make_string (rmf->grating)))
     return -1;
   if (NULL == (info->instrument = isis_make_string (rmf->instrument)))
     return -1;
   if (NULL == (info->arg_string = isis_make_string (rmf->arg_string)))
     return -1;

   return 0;
}

/*}}}*/

int Rmf_set_info (Isis_Rmf_t *rmf, Rmf_Info_Type *info) /*{{{*/
{
   if ((rmf == NULL) || (info == NULL))
     return -1;

   /* don't overwrite arg-string, index or method fields */

   rmf->order = info->order;

   isis_strcpy(rmf->grating,info->grating, sizeof(rmf->grating));
   isis_strcpy(rmf->instrument,info->instrument, sizeof(rmf->instrument));

   return 0;
}

/*}}}*/

/*{{{ RMF access */

int Rmf_includes_effective_area (Isis_Rmf_t *rmf) /*{{{*/
{
   if (rmf == NULL)
     return -1;

   return rmf->includes_effective_area;
}

/*}}}*/

int Rmf_is_identity (Isis_Rmf_t *rmf) /*{{{*/
{
   if (rmf == NULL)
     return -1;

   return (rmf->method == RMF_DELTA);
}

/*}}}*/

char *Rmf_name (Isis_Rmf_t *rmf) /*{{{*/
{
   if (rmf == NULL)
     return NULL;
   return rmf->arg_string;
}

/*}}}*/

int Rmf_index (Isis_Rmf_t *rmf) /*{{{*/
{
   if (rmf == NULL)
     return -1;
   return rmf->index;
}

/*}}}*/

int Rmf_id_list (Isis_Rmf_t *head, unsigned int **ids, unsigned int *num) /*{{{*/
{
   Isis_Rmf_t *r;
   unsigned int n;

   *ids = NULL;
   *num = 0;

   if (head == NULL)
     return 0;

   for (r = head->next; r != NULL; r = r->next)
     {
        *num += 1;
     }

   if (*num == 0)
     return 0;

   if (NULL == (*ids = (unsigned int *) ISIS_MALLOC (*num * sizeof (unsigned int))))
     return -1;

   n = 0;
   for (r = head->next; r != NULL; r = r->next)
     {
        (*ids)[n++] = r->index;
     }

   return 0;
}

/*}}}*/

int Rmf_run_post_fit_method (Isis_Rmf_t *rmf)
{
   if (rmf->post_fit == NULL)
     return 0;
   return (*rmf->post_fit)(rmf);
}

/*}}}*/

int Rmf_apply_rmf (Isis_Rmf_t *rmf, double *x, int num_orig_data, /*{{{*/
                   double *arf_src, int *arf_notice_list,
                   int num_arf_noticed)
{
   int k;

   if (NULL == rmf || NULL == arf_src )
     return -1;

   if (rmf->pre_apply != NULL)
     {
        if (-1 == (*rmf->pre_apply)(rmf))
          return -1;
     }

   for (k = 0; k < num_arf_noticed; k++)
     {
        int nk;
        if (arf_src[k] == 0.0)
          continue;
        nk = (arf_notice_list != NULL) ? arf_notice_list[k] : k;
        if (-1 == rmf->redistribute (rmf, nk, arf_src[k], x, num_orig_data))
          return -1;
     }

   if (rmf->post_apply != NULL)
     {
        if (-1 == (*rmf->post_apply)(rmf))
          return -1;
     }

   return 0;
}

/*}}}*/

int Rmf_find_peaks (Isis_Rmf_t *rmf, double **h_P, int *num) /*{{{*/
{
   double *arf_lo, *arf_hi, *ebounds_lo, *ebounds_hi, *profile;
   unsigned int arf_n, ebounds_n, k;
   unsigned int *indices;
   int status = -1;

   if (rmf == NULL || h_P == NULL || num == NULL)
     return -1;

   profile = arf_lo = arf_hi = ebounds_lo = ebounds_hi = NULL;
   indices = NULL;
   *h_P = NULL;

   /* wavelength grids */
   if ((-1 == rmf->get_arf_grid (rmf, &arf_lo, &arf_hi, &arf_n))
       || (-1 == rmf->get_data_grid (rmf, &ebounds_lo, &ebounds_hi, &ebounds_n, NULL)))
     goto return_error;

   *num = arf_n;

   if ((NULL == (*h_P = (double *) ISIS_MALLOC (arf_n * sizeof(double))))
       || (NULL == (indices = (unsigned int *) ISIS_MALLOC (ebounds_n * sizeof(unsigned int))))
       || (NULL == (profile = (double *) ISIS_MALLOC (ebounds_n * sizeof(double)))))
     goto return_error;

   /* For each arf wavelength, fold a delta-function source
    * through the rmf, find the peak in the output, and record
    * the bin-center wavelength of that peak
    */

   for (k = 0; k < arf_n; k++)
     {
        double peak_value;
        unsigned int i, num_indices;

        for (i = 0; i < ebounds_n; i++)
          profile[i] = 0.0;

        /* set default value */
        (*h_P)[k] = (double) k;

        if (-1 == rmf->redistribute (rmf, k, 1.0, profile, ebounds_n))
          goto return_error;

        peak_value = 0.0;
        for (i = 0; i < ebounds_n; i++)
          {
             if (profile[i] > peak_value)
               peak_value = profile[i];
          }

        if (peak_value == 0.0)
          continue;

        num_indices = 0;
        for (i = 0; i < ebounds_n; i++)
          {
             if (profile[i] == peak_value)
               indices[num_indices++] = i;
          }

        if (num_indices != 0)
          {
#if 0
             (*h_P)[k] = indices[num_indices/2];
#else
             unsigned int hw = 3;
             i = indices[num_indices/2];
             if ((hw <= i) && (i + hw < ebounds_n))
               {
                  double s, xp;
                  unsigned int j;
                  s = xp = 0.0;
                  for (j=i-hw; j<=i+hw; j++)
                    {
                       xp += j * profile[j];
                       s += profile[j];
                    }
                  (*h_P)[k] = xp / s;
               }
#endif
          }
     }

   status = 0;

   return_error:
   ISIS_FREE (indices);
   ISIS_FREE (profile);
   ISIS_FREE (arf_lo);
   ISIS_FREE (arf_hi);
   ISIS_FREE (ebounds_lo);
   ISIS_FREE (ebounds_hi);

   if ((status != 0) && (h_P != NULL))
     ISIS_FREE (*h_P);

   return status;
}

/*}}}*/

static Isis_Rmf_Load_Method_t *get_user_rmf_load_method (char *options) /*{{{*/
{
   static const char *delim = ":;";
   Isis_Rmf_Load_Method_t *load_method;
   char *s, *file, *init_name;

   if (options == NULL)
     return NULL;

   /* options string has form
    *   "libname.so:init_name[;opt1;opt2;opt3;...]"
    */

   if (NULL == (s = isis_make_string (options)))
     return NULL;

   if ((NULL == (file = strtok (s, delim)))
       || (NULL == (init_name = strtok (NULL, delim))))
     {
        ISIS_FREE(s);
        return NULL;
     }

   load_method = (Isis_Rmf_Load_Method_t *)isis_load_function (file, init_name, "rmf");

   ISIS_FREE(s);

   return load_method;
}

/*}}}*/

static int rmf_load_user (Isis_Rmf_t *rmf, void *opt)
{
   Isis_Rmf_Load_Method_t *load_method;
   char *options;

   options = (char *)opt;
   if ((NULL == (load_method = get_user_rmf_load_method (options)))
       || (-1 == (*load_method)(rmf, options)))
     return -1;

   rmf->arg_string = isis_make_string (options);   /* NULL ok */
   return 0;
}

static Isis_Rmf_t *open_rmf (int method, void *options) /*{{{*/
{
   Isis_Rmf_Load_Method_t *load_method = NULL;
   Isis_Rmf_t *rmf;

   switch (method)
     {
      case RMF_FILE:
        load_method = Rmf_load_file;
        break;

      case RMF_DELTA:
        load_method = Rmf_load_delta;
        break;

      case RMF_USER:
        load_method = rmf_load_user;
        break;

      case RMF_SLANG:
        load_method = Rmf_load_slang;
        break;

      default:
        break;
     }

   if (load_method == NULL)
     return NULL;

   if (NULL == (rmf = new_rmf ()))
     return NULL;

   rmf->method = method;
   if (-1 == (*load_method)(rmf, options))
     {
        Rmf_free_rmf (rmf);
        return NULL;
     }
   return rmf;
}
/*}}}*/

int Rmf_load_rmf (Isis_Rmf_t *head, int method, void *options) /*{{{*/
{
   Isis_Rmf_t *rmf;

   if (head == NULL)
     return -1;

   rmf = open_rmf (method, options);
   if (rmf == NULL)
     return -1;

   return rmf_list_append (head, rmf);
}

/*}}}*/

int Rmf_init_rmf (Isis_Rmf_t *rmf, Isis_Rmf_Grid_Type *arf, Isis_Rmf_Grid_Type *ebounds) /*{{{*/
{
   if ((ebounds == NULL) || (arf == NULL) || (rmf == NULL))
     return -1;

   if ((-1 == rmf->set_arf_grid (rmf, arf->bin_lo, arf->bin_hi, arf->nbins))
       || (-1 == rmf->set_data_grid (rmf, ebounds->bin_lo, ebounds->bin_hi, ebounds->nbins)))
     return -1;

   return rmf->init (rmf);
}

/*}}}*/

Isis_Rmf_t *Rmf_create_delta_rmf (Isis_Rmf_Grid_Type *arf, Isis_Rmf_Grid_Type *ebounds) /*{{{*/
{
   Isis_Rmf_t *rmf = open_rmf (RMF_DELTA, NULL);
   if (rmf == NULL)
     return NULL;

   if (-1 == Rmf_init_rmf (rmf, arf, ebounds))
     {
        Rmf_free_rmf (rmf);
        return NULL;
     }

   return rmf;
}

/*}}}*/

int Rmf_rebin_rmf (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int num) /*{{{*/
{
   int swapped;

   if ((rmf == NULL) || (lo == NULL) || (hi == NULL))
     return -1;

   if (num == 0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "0 bins");
        return -1;
     }

   if (rmf->ref_count > 0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "RMF in use, rebin failed");
        return -1;
     }

#if 1
   if (-1 == _isis_fixup_lo_hi_grids (lo, hi, num, &swapped, NULL))
     return -1;
#else
   if ((num > 1)
       && ((lo[0] >= lo[1]) || (hi[0] >= hi[1]) || (hi[0] <= lo[0])))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grids must be in ascending order");
        return -1;
     }
#endif

   if (-1 == rmf->rebin_rmf (rmf, lo, hi, num))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "rebin failed");
        return -1;
     }

   return 0;
}

/*}}}*/

int Rmf_get_data_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   if ((rmf == NULL)
       || (lo == NULL) || (hi == NULL) || (n == NULL))
     return -1;

   return rmf->get_data_grid (rmf, lo, hi, n, NULL);
}

/*}}}*/

int Rmf_get_arf_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   if ((rmf == NULL)
       || (lo == NULL) || (hi == NULL) || (n == NULL))
     return -1;

   return rmf->get_arf_grid (rmf, lo, hi, n);
}

/*}}}*/

