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

/* $Id: hist-cmds.c,v 1.132 2004/09/09 11:31:56 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <slang.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "rmf.h"
#include "arf.h"
#include "plot.h"
#include "util.h"
#include "histogram.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

static Hist_t *Data_List_Head = NULL;
static Isis_Arf_t *Arf_List_Head = NULL;
static Isis_Rmf_t *Rmf_List_Head = NULL;

extern Hist_t *get_histogram_list_head (void);

static unsigned int Label_By_Default = 1;

Hist_t *get_histogram_list_head (void)
{
   return Data_List_Head;
}

/*{{{ start/load/end */

static int init_data_list (void) /*{{{*/
{
   if (NULL == Data_List_Head
       && NULL == (Data_List_Head = Hist_init_list()))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return 0;
}

/*}}}*/

static int init_lists (void) /*{{{*/
{
   if (-1 == init_data_list ())
     return -1;

   if (NULL == Rmf_List_Head
       && NULL == (Rmf_List_Head = Rmf_init_rmf_list()))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (NULL == Arf_List_Head
       && NULL == (Arf_List_Head = Arf_init_arf_list()))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return 0;
}

/*}}}*/

static void _load_data (char * pha_filename, int *just_one) /*{{{*/
{
   int *indices = NULL;
   int num_spectra = 0;
   int ret = -1;

   if (-1 == init_lists ())
     goto error_return;

   ret = Hist_read_fits (Data_List_Head, Arf_List_Head, Rmf_List_Head, pha_filename,
                         &indices, &num_spectra, 1, *just_one);
   if (ret == NOT_FITS_FORMAT)
     ret = Hist_read_ascii (Data_List_Head, pha_filename);

   if ((ret == -1) || (indices == NULL))
     goto error_return;

   if (num_spectra == 1)
     SLang_push_integer (indices[0]);
   else if (num_spectra > 1)
     {
        SLang_Array_Type *at = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &num_spectra, 1);
        if (NULL == at)
          {
             ISIS_FREE (indices);
             goto error_return;
          }
        memcpy ((char *)at->data, (char *)indices, num_spectra * sizeof(int));
        SLang_push_array (at, 1);
     }
   else
     {
        isis_vmesg (WARN, I_READ_FAILED, __FILE__, __LINE__, "%s", pha_filename);
        goto error_return;
     }

   ret = 0;
   error_return:

   if (ret) SLang_push_integer (ret);
   ISIS_FREE (indices);
}

/*}}}*/

static int load_arf (char * filename) /*{{{*/
{
   int id;

   if (NULL == Arf_List_Head
       && NULL == (Arf_List_Head = Arf_init_arf_list()))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (-1 == (id = Arf_read_arf (Arf_List_Head, filename)))
     return -1;

   return id;
}

/*}}}*/

static int load_rmf_internal (int method, void *opt) /*{{{*/
{
   int id;

   if (NULL == Rmf_List_Head
       && NULL == (Rmf_List_Head = Rmf_init_rmf_list()))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (-1 == (id = Rmf_load_rmf (Rmf_List_Head, method, opt)))
     return -1;

   return id;
}

/*}}}*/

static int load_user_rmf (char *options) /*{{{*/
{
   return load_rmf_internal (RMF_USER, options);
}

/*}}}*/

static int load_file_rmf (char *options) /*{{{*/
{
   return load_rmf_internal (RMF_FILE, options);
}

/*}}}*/

static SLang_CStruct_Field_Type Rmf_SLang_Info_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Rmf_SLang_Info_Type, func_ref, "func", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_SLang_Info_Type, arf_bin_lo, "arf_bin_lo", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_SLang_Info_Type, arf_bin_hi, "arf_bin_hi", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_SLang_Info_Type, data_bin_lo, "data_bin_lo", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_SLang_Info_Type, data_bin_hi, "data_bin_hi", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_SLang_Info_Type, threshold, "threshold", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_SLang_Info_Type, client_data, "parms", SLANG_ANY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int load_slang_rmf (void) /*{{{*/
{
   Rmf_SLang_Info_Type info;
   int ret = -1;

   memset ((char *) &info, 0, sizeof(Rmf_SLang_Info_Type));

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&info, Rmf_SLang_Info_Layout))
     return -1;

   if ((info.func_ref == NULL)
       || (NULL == (info.func = SLang_get_fun_from_ref (info.func_ref)))
       || (info.arf_bin_lo == NULL)
       || (info.arf_bin_hi == NULL)
       || (info.data_bin_lo == NULL)
       || (info.data_bin_hi == NULL)
       || (info.threshold < 0))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "One or more load_slang_rmf parmeters are invalid");
        goto free_and_return;
     }

   if ((info.arf_bin_lo->num_elements != info.arf_bin_hi->num_elements)
       || (info.data_bin_lo->num_elements != info.data_bin_hi->num_elements))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "load_slang_rmf: bin_lo/hi grids are incompatible");
        goto free_and_return;
     }

   ret = load_rmf_internal (RMF_SLANG, &info);

free_and_return:

   SLang_free_function (info.func);
   SLang_free_cstruct ((VOID_STAR)&info, Rmf_SLang_Info_Layout);
   return ret;
}

/*}}}*/

/*}}}*/

/*{{{ access */

Hist_t *find_hist (int hist_index) /*{{{*/
{
   if (-1 == init_data_list ())
     return NULL;

   return Hist_find_hist_index (Data_List_Head, hist_index);
}

/*}}}*/

static Isis_Arf_t *find_arf (int arf_index) /*{{{*/
{
   return Arf_find_arf_index (Arf_List_Head, arf_index);
}

/*}}}*/

static Isis_Rmf_t *find_rmf (int rmf_index) /*{{{*/
{
   return Rmf_find_rmf_index (Rmf_List_Head, rmf_index);
}

/*}}}*/

static void _assign_rmf_to_hist (int *rmf_index, int *hist_index) /*{{{*/
{
   /* Might be using RMF to generate fake data */
   if (-1 == init_data_list ())
     return;

   if (-1 == Hist_assign_rmf_to_hist (Data_List_Head, *hist_index,
                                      Rmf_List_Head, *rmf_index))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "assigning RMF %d to data set %d",
                    *rmf_index, *hist_index);
     }
}

/*}}}*/

static void _assign_arf_to_hist (int *arf_index, int *hist_index) /*{{{*/
{
   /* Might be using ARF to generate fake data */
   if (-1 == init_data_list ())
     return;

   if (-1 == Hist_assign_arf_to_hist (Data_List_Head, *hist_index,
                                      Arf_List_Head, *arf_index))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "assigning ARF %d to data set %d",
                    *arf_index, *hist_index);
     }
}

/*}}}*/

static void _assign_rsp_list (int *hist_index) /*{{{*/
{
   SLang_Array_Type *sl_arfs, *sl_rmfs;
   int *arfs, *rmfs;
   int num;

   sl_arfs = sl_rmfs = NULL;

   if (-1 == init_lists ())
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == SLang_pop_array_of_type (&sl_rmfs, SLANG_INT_TYPE)
       || sl_rmfs == NULL
       || -1 == SLang_pop_array_of_type (&sl_arfs, SLANG_INT_TYPE)
       || sl_arfs == NULL
       || sl_arfs->num_elements != sl_rmfs->num_elements)
     {
        SLang_free_array (sl_arfs);
        SLang_free_array (sl_rmfs);
        isis_throw_exception (Isis_Error);
        return;
     }

   num = sl_arfs->num_elements;
   arfs = (int *)sl_arfs->data;
   rmfs = (int *)sl_rmfs->data;

   if (-1 == Hist_assign_rsp_list (Data_List_Head, *hist_index,
                                   Arf_List_Head, arfs, Rmf_List_Head, rmfs, num))
     {
        SLang_free_array (sl_arfs);
        SLang_free_array (sl_rmfs);
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "assigning responses to data set %d", *hist_index);
        return;
     }

   SLang_free_array (sl_arfs);
   SLang_free_array (sl_rmfs);
   return;
}

/*}}}*/

static void _delete_arf (int * arf_index) /*{{{*/
{
   if (-1 == Arf_delete_arf (Arf_List_Head, *arf_index))
     isis_vmesg (INFO, I_FAILED, __FILE__, __LINE__, "deleting ARF %d", *arf_index);
}

/*}}}*/

static void _delete_rmf (int * rmf_index) /*{{{*/
{
   if (-1 == Rmf_delete_rmf (Rmf_List_Head, *rmf_index))
     isis_vmesg (INFO, I_FAILED, __FILE__, __LINE__, "deleting RMF %d", *rmf_index);
}

/*}}}*/

static void _delete_hist (int * hist_index) /*{{{*/
{
   if (-1 == Hist_delete_hist (Data_List_Head, *hist_index))
     isis_vmesg (INFO, I_FAILED, __FILE__, __LINE__, "deleting data set %d", *hist_index);

   (void) update_user_model ();
}

/*}}}*/

static int define_bgd_file (int *hist_index, char *file) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);

   if (-1 == Hist_set_background_from_file (h, file))
     {
        isis_vmesg (INFO, I_FAILED, __FILE__, __LINE__, "loading background from %s", file);
        return -1;
     }

   return 0;
}

/*}}}*/

static int define_bgd (int *hist_index, double *exposure) /*{{{*/
{
   SLang_Array_Type *sl_bgd = NULL;
   SLang_Array_Type *sl_area = NULL;
   unsigned int nbins = 0;
   unsigned int nbins_area = 0;
   double *bgd = NULL;
   double *area = NULL;
   Hist_t *h = NULL;
   char name[32];
   int area_is_vector = 0;

   if (NULL == (h = find_hist (*hist_index)))
     return -1;

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLdo_pop_n (SLang_Num_Function_Args-2);
     }
   else
     {
        if (-1 == SLang_pop_array_of_type (&sl_bgd, SLANG_DOUBLE_TYPE)
            || sl_bgd == NULL
            || sl_bgd->num_elements == 0)
          {
             SLang_free_array (sl_bgd);
             return -1;
          }
        bgd = (double *)sl_bgd->data;
        nbins = sl_bgd->num_elements;
        if (-1 == SLang_pop_array_of_type (&sl_area, SLANG_DOUBLE_TYPE)
            || sl_area == NULL
            || sl_area->num_elements == 0)
          {
             SLang_free_array (sl_area);
             SLang_free_array (sl_bgd);
             return -1;
          }
        area = (double *)sl_area->data;
        nbins_area = sl_area->num_elements;

        if (nbins_area == 1)
          area_is_vector = 0;
        else if (nbins_area == nbins)
          area_is_vector = 1;
        else
          {
             SLang_free_array (sl_area);
             SLang_free_array (sl_bgd);
             return -1;
          }
     }

   if (-1 == Hist_define_background (h, *exposure, area, area_is_vector, bgd, nbins))
     {
        SLang_free_array (sl_bgd);
        SLang_free_array (sl_area);
        return -1;
     }

   sprintf (name, "%c_define_bgd()", COMMENT_CHAR);
   if (-1 == Hist_set_background_name (h, name))
     {
        SLang_free_array (sl_bgd);
        SLang_free_array (sl_area);
        return -1;
     }

   SLang_free_array (sl_bgd);
   SLang_free_array (sl_area);

   return 0;
}

/*}}}*/

static void set_sys_err_frac_intrin (int *hist_index) /*{{{*/
{
   Hist_t *h;
   SLang_Array_Type *sl_sef = NULL;
   double *sef;
   unsigned int n;

   if (NULL == (h = find_hist (*hist_index)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     SLdo_pop ();
   else if (-1 == SLang_pop_array_of_type (&sl_sef, SLANG_DOUBLE_TYPE))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (sl_sef != NULL)
     {
        sef = (double *)sl_sef->data;
        n = sl_sef->num_elements;
     }
   else
     {
        sef = NULL;
        n = 0;
     }

   if (-1 == Hist_define_sys_err_frac (h, sef, n))
     isis_throw_exception (Isis_Error);

   SLang_free_array (sl_sef);
}

/*}}}*/

static void get_sys_err_frac_intrin (int *hist_index) /*{{{*/
{
   SLang_Array_Type *sl_sef = NULL;
   Hist_t *h;
   double *sef = NULL;
   int nbins;

   if (NULL == (h = find_hist (*hist_index)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == Hist_copy_sys_err_frac (h, &sef, &nbins))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (NULL == (sl_sef = SLang_create_array (SLANG_DOUBLE_TYPE, 0, sef, &nbins, 1)))
     isis_throw_exception (Isis_Error);

   SLang_push_array (sl_sef, 1);
}

/*}}}*/

static int _define_hist (void) /*{{{*/
{
   SLang_Array_Type *binlo, *binhi, *value, *uncert;
   Isis_Hist_t x;
   int hist_index;
   unsigned int bin_type, has_grid, nbins;

   binlo = binhi = value = uncert = NULL;
   memset ((char *)&x, 0, sizeof(x));

   if (-1 == init_data_list ())
     return -1;

   if (-1 == SLang_pop_array_of_type (&uncert, SLANG_DOUBLE_TYPE)
       || uncert == NULL
       || -1 == SLang_pop_array_of_type (&value, SLANG_DOUBLE_TYPE)
       || value == NULL
       || -1 == SLang_pop_array_of_type (&binhi, SLANG_DOUBLE_TYPE)
       || binhi == NULL
       || -1 == SLang_pop_array_of_type (&binlo, SLANG_DOUBLE_TYPE)
       || binlo == NULL
       || -1 == SLang_pop_uinteger (&bin_type)
       || value->num_elements == 0)
     {
        /* use slang wrapper to print usage message */
        isis_throw_exception (Isis_Error);
        hist_index = -1;
        goto free_and_return;
     }

   x.nbins = value->num_elements;
   nbins = x.nbins;

   x.val = (double *)value->data;

   if (uncert->num_elements == nbins)
     x.val_err = (double *)uncert->data;

   has_grid = ((binlo->num_elements == nbins)
               && (binhi->num_elements == nbins));
   if (has_grid)
     {
        x.bin_lo = (double *)binlo->data;
        x.bin_hi = (double *)binhi->data;
     }
   else
     {
        /* might define the data set with counts only, and
         * later define the grid from elsewhere */
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "valid data grid not supplied, using default");
     }

   hist_index = Hist_define_data (Data_List_Head, &x, bin_type, has_grid);

   free_and_return:

   SLang_free_array (binlo);
   SLang_free_array (binhi);
   SLang_free_array (value);
   SLang_free_array (uncert);

   return hist_index;
}

/*}}}*/

static void _copy_hist_keywords (int *dst, int *src) /*{{{*/
{
   Hist_t *d, *s;

   d = find_hist (*dst);
   s = find_hist (*src);

   if ((d == NULL) || (s == NULL)
       || (-1 == Hist_copy_histogram_keywords (d, s)))
     isis_throw_exception (Isis_Error);
}

/*}}}*/

static void put_hist (void) /*{{{*/
{
   Hist_t *h;
   SLang_Array_Type *binlo, *binhi, *value, *uncert;
   int hist_index;
   unsigned int n, version;

   binlo = binhi = value = uncert = NULL;

   if (-1 == SLang_pop_array_of_type (&uncert, SLANG_DOUBLE_TYPE)
       || uncert == NULL
       || -1 == SLang_pop_array_of_type (&value, SLANG_DOUBLE_TYPE)
       || value == NULL
       || -1 == SLang_pop_array_of_type (&binhi, SLANG_DOUBLE_TYPE)
       || binhi == NULL
       || -1 == SLang_pop_array_of_type (&binlo, SLANG_DOUBLE_TYPE)
       || binlo == NULL
       || -1 == SLang_pop_uinteger (&version)
       || -1 == SLang_pop_integer (&hist_index))
     {
        /* use slang wrapper to print usage message */
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   if (NULL == (h = find_hist (hist_index)))
     {
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   n = value->num_elements;

   if ((n == 0)
       || (uncert->num_elements != n)
       || (binlo->num_elements != n)
       || (binhi->num_elements != n))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto free_and_return;
     }

   /* FIXME this should be one subroutine call! */
   if (-1 == Hist_replace_hist_grid (h, (double *) binlo->data,
                                     (double *) binhi->data, n)
       || -1 == Hist_replace_hist_data (h, version, (double *) value->data,
                                        (double *) uncert->data, n))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "replacing data set %d", hist_index);
     }

   free_and_return:
   SLang_free_array (binlo);
   SLang_free_array (binhi);
   SLang_free_array (value);
   SLang_free_array (uncert);
}

/*}}}*/

static void put_model_intrin (void) /*{{{*/
{
   SLang_Array_Type *value;
   Hist_t *h;
   int hist_index;
   unsigned int version;

   value = NULL;

   if (-1 == SLang_pop_array_of_type (&value, SLANG_DOUBLE_TYPE)
       || value == NULL
       || -1 == SLang_pop_uinteger (&version)
       || -1 == SLang_pop_integer (&hist_index))
     {
        /* use slang wrapper to print usage message */
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   if (NULL == (h = find_hist (hist_index)))
     {
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   if (value->num_elements == 0)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array size");
        goto free_and_return;
     }

   if (-1 == Hist_replace_hist_data (h, version, (double *) value->data,
                                     NULL, value->num_elements))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "replacing model for data set %d", hist_index);
     }

   free_and_return:
   SLang_free_array (value);
}

/*}}}*/

static int is_fake_data (int *hist_index) /*{{{*/
{
   return Hist_is_fake (find_hist (*hist_index));
}

/*}}}*/

static int set_fake (int *hist_index, int *value) /*{{{*/
{
   return Hist_set_fake (find_hist (*hist_index), *value);
}

/*}}}*/

static void get_hist (int *hist_index, unsigned int *version) /*{{{*/
{
   Hist_t *h;
   Isis_Hist_t g;
   SLang_Array_Type *bin_lo = NULL;
   SLang_Array_Type *bin_hi = NULL;
   SLang_Array_Type *value = NULL;
   SLang_Array_Type *uncert = NULL;
   int nbins;

   if (NULL == (h = find_hist (*hist_index)))
     goto push_values;

   nbins = Hist_hist_size (h, *version);
   if (-1 == nbins)
     goto push_values;

   if (NULL == (bin_lo = SLang_create_array (SLANG_DOUBLE_TYPE, 0,
                                             NULL, &nbins, 1))
       || NULL == (bin_hi = SLang_create_array (SLANG_DOUBLE_TYPE, 0,
                                                NULL, &nbins, 1))
       || NULL == (value = SLang_create_array (SLANG_DOUBLE_TYPE, 0,
                                                NULL, &nbins, 1))
       || NULL == (uncert = SLang_create_array (SLANG_DOUBLE_TYPE, 0,
                                                  NULL, &nbins, 1)))
     goto push_values;

   if (-1 == Hist_get_hist_grid (h, *version, &g)
       || -1 == Hist_copy_hist_data (h, *version, (double *) value->data,
                                     (double *) uncert->data))
     goto push_values;

   if (g.nbins != nbins)
     {
        isis_vmesg (FAIL, I_INTERNAL, __FILE__, __LINE__, "inconsistent grid");
        goto push_values;
     }

   memcpy ((char *)bin_lo->data, (char *)g.bin_lo, nbins * sizeof(double));
   memcpy ((char *)bin_hi->data, (char *)g.bin_hi, nbins * sizeof(double));

   /* always push values onto the run-time stack (NULLs are ok),
    * otherwise we might generate an interpreter error (stack underflow)
    */

   push_values:

   SLang_push_array (bin_lo, 1);
   SLang_push_array (bin_hi, 1);
   SLang_push_array (value, 1);

   if (is_data(*version))
     SLang_push_array (uncert, 1);
   else
     SLang_free_array (uncert);
}

/*}}}*/

static int have_data (int *hist_index) /*{{{*/
{
   if (NULL != _Hist_find_hist_index (Data_List_Head, *hist_index))
     return 1;

   return 0;
}

/*}}}*/

static void get_hist_notice_info (int * hist_index, unsigned int *version) /*{{{*/
{
   SLang_Array_Type * notice_list = NULL;
   SLang_Array_Type * notice = NULL;
   SLang_Array_Type * rebin = NULL;
   int *prebin = NULL;
   int orig_nbins;
   Isis_Hist_t g;
   Hist_t *h;

   if ((NULL == (h = find_hist (*hist_index)))
       || -1 == Hist_get_hist_grid (h, *version, &g)
       || -1 == Hist_get_hist_rebin_info (h, &prebin, &orig_nbins))
     goto push_values;

   if ((NULL == (notice_list = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &g.n_notice, 1)))
       || (NULL == (notice = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &g.nbins, 1)))
       || (NULL == (rebin = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &orig_nbins, 1))))
     goto push_values;

   memcpy ((char *)notice_list->data, (char *)g.notice_list, g.n_notice * sizeof(int));
   memcpy ((char *)notice->data, (char *)g.notice, g.nbins * sizeof(int));
   memcpy ((char *)rebin->data, (char *)prebin, orig_nbins * sizeof(int));

   push_values:
   SLang_push_array (notice, 1);
   SLang_push_array (notice_list, 1);
   SLang_push_array (rebin, 1);
}

/*}}}*/

static void put_arf (void) /*{{{*/
{
   SLang_Array_Type *binlo, *binhi, *arf, *arf_err;
   int nbins, arf_index;
   Isis_Arf_t *a;

   binlo = binhi = arf = arf_err = NULL;

   if (-1 == SLang_pop_array_of_type (&arf_err, SLANG_DOUBLE_TYPE)
       || arf_err == NULL
       || -1 == SLang_pop_array_of_type (&arf, SLANG_DOUBLE_TYPE)
       || arf == NULL
       || -1 == SLang_pop_array_of_type (&binhi, SLANG_DOUBLE_TYPE)
       || binhi == NULL
       || -1 == SLang_pop_array_of_type (&binlo, SLANG_DOUBLE_TYPE)
       || binlo == NULL
       || -1 == SLang_pop_integer (&arf_index))
     {
        /* use slang wrapper to print usage message */
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   if (NULL == (a = find_arf (arf_index)))
     goto free_and_return;

   if (-1 == (nbins = Arf_arf_size (a)))
     {
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   if (arf->num_elements != (unsigned) nbins
       || arf_err->num_elements != (unsigned) nbins
       || binlo->num_elements != (unsigned) nbins
       || binhi->num_elements != (unsigned) nbins)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto free_and_return;
     }

   if (-1 == Arf_put_arf (a, (double *) arf->data, (double *) arf_err->data,
                          (double *) binlo->data, (double *) binhi->data))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "Couldn't reset ARF data");
     }

   free_and_return:
   SLang_free_array (binlo);
   SLang_free_array (binhi);
   SLang_free_array (arf);
   SLang_free_array (arf_err);
}

/*}}}*/

static void get_arf (int * arf_index) /*{{{*/
{
   SLang_Array_Type * bin_lo = NULL;
   SLang_Array_Type * bin_hi = NULL;
   SLang_Array_Type * arf = NULL;
   SLang_Array_Type * arf_err = NULL;
   Isis_Arf_t *a = find_arf (*arf_index);
   int nbins;

   if ((a == NULL)
       || -1 == (nbins = Arf_arf_size (a)))
     goto push_values;

   bin_lo  = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nbins, 1);
   bin_hi  = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nbins, 1);
   arf     = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nbins, 1);
   arf_err = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nbins, 1);
   if (NULL == bin_lo || NULL == bin_hi || NULL == arf || NULL == arf_err)
     goto push_values;

   if (-1 == Arf_get_arf (a, (double *) arf->data, (double *) arf_err->data,
                           (double *) bin_lo->data, (double *) bin_hi->data))
     goto push_values;

   /* always push values onto the run-time stack (NULLs are ok),
    * otherwise we might generate an interpreter error (stack underflow)
    */

   push_values:
   SLang_push_array (bin_lo, 1);
   SLang_push_array (bin_hi, 1);
   SLang_push_array (arf, 1);
   SLang_push_array (arf_err, 1);
}

/*}}}*/

static double _get_arf_exposure_time (int *arf_index) /*{{{*/
{
   Isis_Arf_t *a = find_arf (*arf_index);
   double exposure = 1.0;

   if ((a == NULL)
       || -1 == Arf_get_arf_exposure (a, &exposure))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't get exposure from ARF %d",
                    *arf_index);
     }

   return exposure;
}

/*}}}*/

static void _set_arf_exposure_time (int *arf_index, double *exposure) /*{{{*/
{
   Isis_Arf_t *a = find_arf (*arf_index);

   if ((a == NULL)
       || -1 == Arf_set_arf_exposure (a, *exposure))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set exposure for ARF %d",
                    *arf_index);
     }
}
/*}}}*/

static int _define_arf (void) /*{{{*/
{
   SLang_Array_Type *sl_binlo, *sl_binhi, *sl_arf, *sl_arf_err;
   double *binlo, *binhi, *arf, *arf_err;
   int arf_index;
   unsigned int nbins;

   binlo = binhi = arf = arf_err = NULL;
   sl_binlo = sl_binhi = sl_arf = sl_arf_err = NULL;

   if (NULL == Arf_List_Head
       && NULL == (Arf_List_Head = Arf_init_arf_list()))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (-1 == SLang_pop_array_of_type (&sl_arf_err, SLANG_DOUBLE_TYPE)
       || sl_arf_err == NULL
       || -1 == SLang_pop_array_of_type (&sl_arf, SLANG_DOUBLE_TYPE)
       || sl_arf == NULL
       || -1 == SLang_pop_array_of_type (&sl_binhi, SLANG_DOUBLE_TYPE)
       || sl_binhi == NULL
       || -1 == SLang_pop_array_of_type (&sl_binlo, SLANG_DOUBLE_TYPE)
       || sl_binlo == NULL)

     {
        /* use slang wrapper to print usage message */
        isis_throw_exception (Isis_Error);
        arf_index = -1;
        goto free_and_return;
     }

   nbins = sl_arf->num_elements;

   if ((nbins == 0)
       || (sl_binlo->num_elements != nbins)
       || (sl_binhi->num_elements != nbins)
       || (sl_arf_err->num_elements != nbins))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        arf_index = -1;
        goto free_and_return;
     }

   binlo = (double *)sl_binlo->data;
   binhi = (double *)sl_binhi->data;
   arf = (double *)sl_arf->data;
   arf_err = (double *)sl_arf_err->data;

   arf_index = Arf_define_arf (Arf_List_Head, nbins, binlo, binhi, arf, arf_err);

   free_and_return:

   SLang_free_array (sl_binlo);
   SLang_free_array (sl_binhi);
   SLang_free_array (sl_arf);
   SLang_free_array (sl_arf_err);

   return arf_index;
}

/*}}}*/

static SLang_CStruct_Field_Type Arf_Info_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, exposure, "exposure", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, order, "order", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, part, "part", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, srcid, "srcid", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, object, "object", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, grating, "grating", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, instrument, "instrument", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Arf_Info_Type, file, "file", SLANG_STRING_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int pop_arf_info (Arf_Info_Type *info) /*{{{*/
{
   SLang_Array_Type *f;

   if (info == NULL)
     return -1;

   if (-1 == SLang_pop_cstruct ((VOID_STAR)info, Arf_Info_Layout))
     return -1;

   if ((-1 == SLang_pop_array_of_type (&f, SLANG_DOUBLE_TYPE))
       || f == NULL)
     return -1;

   if (f->num_elements > 1)
     {
        int size;
        info->fracexpo_is_vector = 1;
        info->nbins = f->num_elements;
        size = info->nbins * sizeof(double);
        if (NULL == (info->fracexpo.v = (double *) ISIS_MALLOC(size)))
          {
             SLang_free_array (f);
             return -1;
          }
        memcpy ((char *)info->fracexpo.v, (char *)f->data, size);
     }
   else
     {
        info->fracexpo_is_vector = 0;
        info->nbins = 1;
        if (f->num_elements == 1)
          info->fracexpo.s = *((double *)f->data);
        else
          info->fracexpo.s = 1.0;
     }

   SLang_free_array (f);

   return 0;
}

/*}}}*/

static int push_arf_info (Arf_Info_Type *info) /*{{{*/
{
   SLang_Array_Type *f;

   if (info == NULL)
     return -1;

   if (-1 == SLang_push_cstruct ((VOID_STAR)info, Arf_Info_Layout))
     return -1;

   f = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &info->nbins, 1);
   if (f == NULL)
     {
        SLang_push_array (f, 1);
        return -1;
     }

   if (info->nbins > 1)
     {
        int size;
        size = info->nbins * sizeof(double);
        memcpy ((char *)f->data, (char *)info->fracexpo.v, size);
     }
   else
     {
        int i = 0;
        SLang_set_array_element (f, &i, &info->fracexpo.s);
     }

   SLang_push_array (f, 1);

   return 0;
}

/*}}}*/

static void _set_arf_info (int *arf_index) /*{{{*/
{
   Arf_Info_Type info ;
   Isis_Arf_t *a;

   memset ((char *)&info, 0, sizeof(info));

   if (-1 == pop_arf_info (&info))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't pop arf info",
                    *arf_index);
        goto free_and_return;
     }

   if (NULL == (a = find_arf (*arf_index)))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't find arf %d",
                    *arf_index);
        goto free_and_return;
     }

   if (-1 == Arf_set_arf_info (a, &info))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set info for arf %d",
                    *arf_index);
        goto free_and_return;
     }

   free_and_return:
   if (info.fracexpo_is_vector)
     {
        ISIS_FREE (info.fracexpo.v);
     }
   SLang_free_cstruct ((VOID_STAR)&info, Arf_Info_Layout);
}

/*}}}*/

static int _get_arf_info (int *arf_index) /*{{{*/
{
   Arf_Info_Type info;
   Isis_Arf_t *a;
   int ret = 0;

   memset ((char *)&info, 0, sizeof(info));

   if ((NULL == (a = find_arf (*arf_index)))
       || (-1 == Arf_get_arf_info (a, &info)))
     ret = -1;

   if (ret == 0)
     ret = push_arf_info (&info);

   Arf_free_arf_info (&info);

   return ret;
}

/*}}}*/

static SLang_CStruct_Field_Type Rmf_Info_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Rmf_Info_Type, index, "index", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_Info_Type, method, "method", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_Info_Type, order, "order", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_Info_Type, grating, "grating", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_Info_Type, instrument, "instrument", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Rmf_Info_Type, arg_string, "arg_string", SLANG_STRING_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static void _set_rmf_info (int *rmf_index) /*{{{*/
{
   Rmf_Info_Type info ;
   Isis_Rmf_t *rmf;

   memset ((char *)&info, 0, sizeof(info));

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&info, Rmf_Info_Layout))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't pop rmf info",
                    *rmf_index);
        goto free_and_return;
     }

   if (NULL == (rmf = find_rmf (*rmf_index)))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't find rmf %d",
                    *rmf_index);
        goto free_and_return;
     }

   if (-1 == Rmf_set_info (rmf, &info))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set info for rmf %d",
                    *rmf_index);
        goto free_and_return;
     }

   free_and_return:
   SLang_free_cstruct ((VOID_STAR)&info, Rmf_Info_Layout);
}

/*}}}*/

static void _get_rmf_info (int *rmf_index) /*{{{*/
{
   Rmf_Info_Type info;
   Isis_Rmf_t *rmf;
   int status;

   memset ((char *)&info, 0, sizeof(info));

   if ((NULL == (rmf = find_rmf (*rmf_index)))
       || (-1 == Rmf_get_info (rmf, &info)))
     {
        Rmf_free_info (&info);
        SLang_push_null ();
        return;
     }

   status = SLang_push_cstruct ((VOID_STAR)&info, Rmf_Info_Layout);
   Rmf_free_info (&info);
}

/*}}}*/

static void _set_hist_frame_time (int *hist_index, double *frame_time) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_set_frame_time (h, *frame_time))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set frame-time for data set %d",
                    *hist_index);
     }
}

/*}}}*/

static double _get_hist_frame_time (int *hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   double ft = 1.0;

   if (-1 == Hist_get_frame_time (h, &ft))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't get frame-time for data set %d",
                    *hist_index);
     }

   return ft;
}

/*}}}*/

static void _get_kernel (int *hist_index) /*{{{*/
{
   Hist_t *h;
   Isis_Kernel_t *k;
   char *s;

   if (NULL == (h = find_hist (*hist_index)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (NULL == (k = Hist_get_kernel (h)))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__,
                    "dataset %d:  kernel uninitialized (until model is evaluated)",
                    *hist_index);
        return;
     }

   s = Hist_get_kernel_params (h);
   if (s == NULL) s = k->kernel_def->kernel_name;

   SLang_push_string (s);
}

/*}}}*/

static void get_flux_corr_weights_intrin (int *hist_index) /*{{{*/
{
   SLang_Array_Type *sl_wt = NULL;
   double *wt = NULL;
   Hist_t *h;
   int n;

   if (NULL == (h = find_hist (*hist_index)))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "data set %d not found", *hist_index);
        goto return_result;
     }

   if (-1 == Hist_get_flux_corr_weights (h, &wt, &n))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "retrieving flux_corr weights for data set %d", *hist_index);
        goto return_result;
     }

   if (NULL == (sl_wt = SLang_create_array (SLANG_DOUBLE_TYPE, 0, wt, &n, 1)))
     {
        ISIS_FREE(wt);
        goto return_result;
     }

   return_result:

   SLang_push_array (sl_wt, 1);
}

/*}}}*/

static void do_flux_correct (int hist_index, int method, double threshold) /*{{{*/
{
   SLang_Array_Type *bgd = NULL;
   Hist_t *h;
   char *options = NULL;
   double *params = NULL;
   unsigned int num_params;
   double *bd;

   /* To initialize the kernel parameters, etc.
    * slang wrapper must call sync_model_with_data()
    * before calling this routine
    */

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     SLdo_pop ();
   else if (-1 == SLpop_string (&options))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (NULL == (h = find_hist (hist_index)))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "data set %d not found", hist_index);
        return;
     }

   if (0 == Hist_num_data_noticed (h))
     return;

   if (-1 == pop_instrumental_background (&bgd, h))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "flux_corr failed retrieving background for dataset %d",
                    hist_index);
        SLang_free_array (bgd);
        return;
     }

   if (-1 == get_kernel_params_for_hist (h, &params, &num_params))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "flux_corr failed retrieving kernel params for dataset %d",
                    hist_index);
        SLang_free_array (bgd);
        ISIS_FREE (params);
        return;
     }

   if (bgd)
     bd = (double *)bgd->data;
   else
     bd = NULL;

   if (-1 == Hist_flux_correct (h, threshold, bd, method, params, num_params, options))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "flux-correcting dataset %d", hist_index);
        SLang_free_array (bgd);
        ISIS_FREE (params);
        return;
     }

   SLang_free_array (bgd);
   ISIS_FREE (params);
}

/*}}}*/

static void _flux_correct (int *hist_index, double *threshold) /*{{{*/
{
   do_flux_correct (*hist_index, 0, *threshold);
}

/*}}}*/

static void _flux_correct_model_counts (int *hist_index, double *threshold) /*{{{*/
{
   do_flux_correct (*hist_index, 1, *threshold);
}

/*}}}*/

static SLang_CStruct_Field_Type Hist_Info_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, spec_num, "spec_num", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, order, "order", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, part, "part", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, srcid, "srcid", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, exclude, "exclude", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, combo_id, "combo_id", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, combo_weight, "combo_weight", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, tstart, "tstart", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, file, "file", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Info_Type, bgd_file, "bgd_file", SLANG_STRING_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int _get_hist_info (int *hist_index) /*{{{*/
{
   SLang_Array_Type *sl_arfs, *sl_rmfs;
   int *arfs, *rmfs, num;
   Hist_Info_Type info;
   Hist_t *h;
   int ret = 0;

   sl_arfs = sl_rmfs = NULL;
   arfs = rmfs = NULL;

   memset ((char *)&info, 0, sizeof(info));

   if ((NULL == (h = find_hist (*hist_index)))
       || (-1 == Hist_get_info (h, &info))
       || (-1 == Hist_get_rsp_list (h, &arfs, &rmfs, &num)))
     {
        ret = -1;
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "getting info for dataset %d", *hist_index);
     }

   if (num > 0)
     {
        sl_arfs = SLang_create_array (SLANG_INT_TYPE, 0, arfs, &num, 1);
        sl_rmfs = SLang_create_array (SLANG_INT_TYPE, 0, rmfs, &num, 1);
     }

   if (-1 == SLang_push_cstruct ((VOID_STAR)&info, Hist_Info_Layout))
     {
        ret = -1;
        isis_vmesg (FAIL, I_INTERNAL, __FILE__, __LINE__, "creating info struct");
     }

   SLang_push_string (info.object);
   SLang_push_string (info.instrument);
   SLang_push_string (info.grating);
   SLang_push_array (sl_arfs, 1);
   SLang_push_array (sl_rmfs, 1);

   return ret;
}

/*}}}*/

static void _set_hist_object_name (int *hist_index, char *object) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   (void) Hist_set_object_name (h, object);
}

/*}}}*/

static void _set_hist_info (int *hist_index) /*{{{*/
{
   Hist_Info_Type info;
   Hist_t *h;

   memset ((char *)&info, 0, sizeof(info));

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&info, Hist_Info_Layout))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't pop cstruct",
                    *hist_index);
     }

   if (NULL == (h = find_hist (*hist_index)))
     {
        SLang_free_cstruct ((VOID_STAR)&info, Hist_Info_Layout);
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't find dataset %d",
                    *hist_index);
        return;
     }

   if (-1 == Hist_set_info (h, &info))
     {
        SLang_free_cstruct ((VOID_STAR)&info, Hist_Info_Layout);
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set info for dataset %d", *hist_index);
     }
}

/*}}}*/

static int is_grouped_data (int *hist_index) /*{{{*/
{
   Hist_t *h;
   int orig_nbins, nbins;

   if (-1 == init_data_list ())
     return -1;

   h = _Hist_find_hist_index (Data_List_Head, *hist_index);
   if (h == NULL)
     return 0;

   orig_nbins = Hist_orig_hist_size (h);
   nbins = Hist_hist_size (h, H_DATA | ~H_FLUX);

   if (orig_nbins < 0 || nbins < 0)
     {
        isis_throw_exception (Isis_Error);
        return 0;
     }

   return (orig_nbins != nbins);
}

/*}}}*/

static void _get_histogram_exposure_time (int *hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   double exposure = 1.0;

   if (-1 == Hist_get_exposure (h, &exposure))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't get exposure time for dataset %d",
                    *hist_index);
     }

   SLang_push_double (exposure);
}

/*}}}*/

static void _get_back_exposure (int *hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   double back_exposure = 1.0;

   if (-1 == Hist_get_back_exposure (h, &back_exposure))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't get background exposure time for dataset %d",
                    *hist_index);
     }

   SLang_push_double (back_exposure);
}

/*}}}*/

static void _set_back_exposure (int *hist_index, double *t) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_set_back_exposure (h, *t))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set exposure for spectrum %d background", *hist_index);
     }
}

/*}}}*/

static void get_region_area (int hist_index, /*{{{*/
                             int (*getfun)(Hist_t *, double **, int *))
{
   SLang_Array_Type *sl_area = NULL;
   Hist_t *h = find_hist (hist_index);
   double *region_area = NULL;
   int num;

   if (-1 == (*getfun)(h, &region_area, &num))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't get backscale dataset %d", hist_index);
     }

   sl_area = SLang_create_array (SLANG_DOUBLE_TYPE, 0, region_area, &num, 1);
   SLang_push_array (sl_area, 1);
}

/*}}}*/

static void set_region_area (int hist_index, /*{{{*/
                             int (*setfun)(Hist_t *, double *, int))
{
   SLang_Array_Type *sl_area = NULL;
   Hist_t *h = find_hist (hist_index);
   double *area;
   int nbins_area;
   int ret = -1;

   if (-1 == SLang_pop_array_of_type (&sl_area, SLANG_DOUBLE_TYPE)
       || sl_area == NULL
       || sl_area->num_elements == 0)
     goto return_error;

   area = (double *)sl_area->data;
   nbins_area = sl_area->num_elements;

   if (-1 == (*setfun)(h, area, nbins_area))
     goto return_error;

   ret = 0;
   return_error:

   if (ret)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set backscale for dataset %d", hist_index);
     }

   SLang_free_array (sl_area);
}

/*}}}*/

static void _get_data_region_area (int *hist_index) /*{{{*/
{
   get_region_area (*hist_index, &Hist_get_data_region_area);
}

/*}}}*/

static void _get_back_region_area (int *hist_index) /*{{{*/
{
   get_region_area (*hist_index, &Hist_get_back_region_area);
}

/*}}}*/

static void _set_data_region_area (int *hist_index) /*{{{*/
{
   set_region_area (*hist_index, &Hist_set_data_region_area);
}

/*}}}*/

static void _set_back_region_area (int *hist_index) /*{{{*/
{
   set_region_area (*hist_index, &Hist_set_back_region_area);
}

/*}}}*/

static void _set_histogram_exposure_time (int *hist_index, double *t) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_set_exposure (h, *t))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set exposure for spectrum %d", *hist_index);
     }
}

/*}}}*/

static void _set_instrument_background (int *hist_index, char *fun_string) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_set_instrumental_background_hook_name (h, fun_string))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set background function for dataset %d",
                      *hist_index);
     }

   (void) sync_model_with_data ();
}

/*}}}*/

static void _get_instrumental_background_hook_name (int *hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   char *s = Hist_get_instrumental_background_hook_name (h);
   SLang_push_string (s);
}

/*}}}*/

/*}}}*/

static void _all_data (int *noticed) /*{{{*/
{
   SLang_Array_Type *sl_ids = NULL;
   unsigned int *ids = NULL;
   unsigned int num;

   if (0 == Hist_id_list (Data_List_Head, *noticed, &ids, &num))
     {
        int n = num;
        if (ids != NULL)
          sl_ids = SLang_create_array (SLANG_UINT_TYPE, 0, ids, &n, 1);
     }

   SLang_push_array (sl_ids, 1);
}

/*}}}*/

static void _all_arfs (void) /*{{{*/
{
   SLang_Array_Type *sl_ids = NULL;
   unsigned int *ids = NULL;
   unsigned int num;

   if (0 == Arf_id_list (Arf_List_Head, &ids, &num))
     {
        int n = num;
        if (ids != NULL)
          sl_ids = SLang_create_array (SLANG_UINT_TYPE, 0, ids, &n, 1);
     }

   SLang_push_array (sl_ids, 1);
}

/*}}}*/

static void _all_rmfs (void) /*{{{*/
{
   SLang_Array_Type *sl_ids = NULL;
   unsigned int *ids = NULL;
   unsigned int num;

   if (0 == Rmf_id_list (Rmf_List_Head, &ids, &num))
     {
        int n = num;
        if (ids != NULL)
          sl_ids = SLang_create_array (SLANG_UINT_TYPE, 0, ids, &n, 1);
     }

   SLang_push_array (sl_ids, 1);
}

/*}}}*/

/*{{{ notice / ignore */

static void _ignore_bad (int *hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   Hist_ignore_bad (h);
}

/*}}}*/

static void _set_exclude_flag (int *hist_index, int *exclude) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_set_exclude_flag (h, *exclude))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't set exclude state for spectrum %d", *hist_index);
     }
}

/*}}}*/

static void set_notice (int *value, double *lo, double *hi, int *hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);

   if (-1 == Hist_set_notice (h, *value, *lo, *hi))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't notice/ignore bins for dataset %d",
                    *hist_index);
     }
}

/*}}}*/

static void set_notice_using_mask (int *hist_index) /*{{{*/
{
   SLang_Array_Type *sl_mask = NULL;
   Hist_t *h;
   int *mask;
   int nbins;

   if (NULL == (h = find_hist (*hist_index)))
     return;

   if (-1 == SLang_pop_array_of_type (&sl_mask, SLANG_INT_TYPE)
       || sl_mask == NULL)
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   mask = (int *)sl_mask->data;
   nbins = sl_mask->num_elements;

   if (-1 == Hist_set_notice_using_mask (h, mask, nbins))
     {
        SLang_free_array (sl_mask);
        sl_mask = NULL;
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't mask bins for dataset %d",
                    *hist_index);
     }

   SLang_free_array (sl_mask);
}

/*}}}*/

static void set_notice_using_list (int *hist_index, int *value) /*{{{*/
{
   SLang_Array_Type *sl_list = NULL;
   Hist_t *h;
   unsigned int *list;
   unsigned int n;

   if (NULL == (h = find_hist (*hist_index)))
     return;

   if (-1 == SLang_pop_array_of_type (&sl_list, SLANG_UINT_TYPE)
       || sl_list == NULL)
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   list = (unsigned int *)sl_list->data;
   n = sl_list->num_elements;

   if (-1 == Hist_set_notice_using_list (h, *value, list, n))
     {
        SLang_free_array (sl_list);
        sl_list = NULL;
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't reset noticed bins for dataset %d",
                    *hist_index);
     }

   SLang_free_array (sl_list);
}

/*}}}*/

/*}}}*/

static void set_dataset_metadata (int *hist_index) /*{{{*/
{
   SLang_Any_Type *meta = NULL;
   Hist_t *h;

   if (NULL == (h = find_hist (*hist_index)))
     return;

   if ((-1 == SLang_pop_anytype (&meta))
       || (-1 == Hist_set_metadata (h, meta)))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "failed setting metadata for dataset %d", *hist_index);
     }
}

/*}}}*/

static void get_dataset_metadata (int *hist_index) /*{{{*/
{
   Hist_t *h;

   if (NULL == (h = find_hist (*hist_index)))
     return;

   SLang_push_anytype (Hist_get_metadata (h));
}

/*}}}*/

static void set_post_model_hook (void) /*{{{*/
{
   SLang_Name_Type *post_model_hook = NULL;
   int hist_index;
   Hist_t *h;

   set_hook_from_stack (&post_model_hook);

   if (-1 == SLang_pop_integer (&hist_index))
     return;

   if (NULL == (h = find_hist (hist_index)))
     return;

   (void) Hist_set_post_model_hook (h, post_model_hook, SLang_free_function);
}

/*}}}*/

static void set_pre_combine_hook (void) /*{{{*/
{
   SLang_Name_Type *pre_combine_hook = NULL;
   int hist_index;
   Hist_t *h;

   set_hook_from_stack (&pre_combine_hook);

   if (-1 == SLang_pop_integer (&hist_index))
     return;

   if (NULL == (h = find_hist (hist_index)))
     return;

   (void) Hist_set_pre_combine_hook (h, pre_combine_hook);
}

/*}}}*/

static void assign_model_intrin (int *num_extra_args) /*{{{*/
{
   SLang_Name_Type *fun_ptr = NULL;
   Isis_Arg_Type *args = NULL;
   int hist_index;
   Hist_t *h;

   if (*num_extra_args > 0)
     args = isis_pop_args (*num_extra_args);

   set_hook_from_stack (&fun_ptr);

   if (-1 == SLang_pop_integer (&hist_index))
     return;

   if (NULL == (h = find_hist (hist_index)))
     return;

   (void) Hist_assign_model (h, fun_ptr, args);

   (void) sync_model_with_data ();
}

/*}}}*/

/*{{{ plotting */

static void _plot_hist (int * hist_index, unsigned int *version, int *style, /*{{{*/
                        int *overlay, int *residuals)
{
   Hist_Plot_Tune_Type info;
   Hist_t *h = find_hist (*hist_index);
   Plot_t *fmt;
   Isis_Fit_Statistic_Type *s;

   if (h == NULL)
     return;

   if ((NULL == (fmt = current_format ())
        || (force_open_plot_device() < 0)))
     return;

   if (*overlay == 0)
     reset_current_pane (fmt);

   /* overrides auto-incr on overlay
    * but won't override fmt style
    */
   info.use_style = (*style > 0) ? style : NULL;
   info.residuals = (*residuals) ? Isis_Residual_Plot_Type : 0;
   info.overlay = *overlay;
   info.default_labels = Label_By_Default;
   info.rate = Plot_bin_density (fmt);

   s = info.residuals ? current_fit_statistic () : NULL;

   if (-1 == Hist_plot (h, *version, &info, fmt, s))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "plot failed");
        return;
     }
   else
     Plot_auto_incr_line_type (fmt);

}
/*}}}*/

static void _set_hist_color (int *hist_index, int *color) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_set_color (h, *color))
     isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "color not set for dataset %d", *hist_index);
}

/*}}}*/

static void _unset_hist_color (int *hist_index) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_unset_color (h))
     isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed to unset color for dataset %d", *hist_index);
}

/*}}}*/

/*}}}*/

/*{{{ region stats */

static SLang_CStruct_Field_Type Region_Stat_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, min, "min", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, max, "max", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, sum, "sum", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, sum_err, "sum_err", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, net, "net", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, net_err, "net_err", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, centroid, "centroid", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, centroid_err, "centroid_err", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, eqwidth, "eqwidth", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, eqwidth_err, "eqwidth_err", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, contin, "contin", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, slope, "slope", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Hist_Stat_t, nbins, "nbins", SLANG_INT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static void _get_region_sum (int *hist_index, unsigned int *version, /*{{{*/
                             double *xmin, double *xmax,
                             double *ymin, double *ymax)
{
   Hist_t *h = find_hist (*hist_index);
   Hist_Stat_t s;
   double x[2], y[2];

   if (h == NULL)
     return;

   x[0] = *xmin;      y[0] = *ymin;
   x[1] = *xmax;      y[1] = *ymax;

   if (-1 == Hist_region_stats (h, *version, &s, x, y))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "couldn't get region sum for dataset %d", *hist_index);
        return;
     }

   if (-1 == SLang_push_cstruct ((VOID_STAR)&s, Region_Stat_Layout))
     {
        isis_throw_exception (Isis_Error);
        return;
     }
}
/*}}}*/

static void _cursor_region_stats (int * hist_index, int *version, /*{{{*/
                                  int * subtract, char * file)
{
   Hist_t *h = find_hist (*hist_index);
   FILE * fp = NULL;
   Plot_t *fmt;
   char *msg;

   if (h == NULL)
     return;

   if (NULL == (fmt = current_format ()))
     return;

   if (!Plot_bin_density (fmt))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "data not plotted as bin-density");
        return;
     }

   errno = 0;
   if (NULL == (fp = fopen (file, "a")))
     {
        isis_vmesg (FAIL, I_WRITE_OPEN_FAILED, __FILE__, __LINE__, "%s", file);
        return;
     }

   if (*subtract)
     msg = "appending continuum-subtracted results to %s";
   else
     msg = "appending results to %s (NOT continuum subtracted)";

   isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, msg, file);

   /* print header */
   if (-1 == Hist_print_stats (fp, NULL))
     goto finish;

   for (;;)
     {
        Hist_Stat_t s;
        float fx[2], fy[2];
        double x[2], y[2];

        fputc ('\n', stdout); fflush(stdout);

        if (0 != Plot_read_line_endpoints (fx, fy, *subtract)
            || -1 == change_xunits_to_angstrom (fx, fy))
          goto finish;

        x[0] = (double) fx[0];
        x[1] = (double) fx[1];

        if (*subtract)
          {
             y[0] = (double) fy[0];
             y[1] = (double) fy[1];
          }
        else
          {
             y[0] = y[1] = 0.0;
          }

        if (-1 == Hist_region_stats (h, *version, &s, x, y))
          isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "couldn't get region sum for data set %d", *hist_index);
        else if (-1 == Hist_print_stats (fp, &s))
          goto finish;
        else
          {
             isis_vmesg (WARN, I_INFO, __FILE__, __LINE__,
                         "saved region [%12.6e, %12.6e] angstrom",
                         x[0], x[1]);
          }
     }

   finish:
   isis_fclose (fp);
}
/*}}}*/

/*}}}*/

/*{{{ rebinning */

static void get_stat_error_hook (int *indx) /*{{{*/
{
   Hist_t *h;

   if (NULL == (h = find_hist (*indx)))
     return;

   (void) SLang_push_function (Hist_get_stat_error_hook (h));
}

/*}}}*/

static void set_stat_error_hook (int *indx) /*{{{*/
{
   SLang_Name_Type *stat_error_hook = NULL;
   Hist_t *h;

   if (NULL == (h = find_hist (*indx)))
     return;

   set_hook_from_stack (&stat_error_hook);

   if (-1 == Hist_set_stat_error_hook (h, stat_error_hook, SLang_free_function))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "setting dataset %d rebin error hook", *indx);
     }
}

/*}}}*/

static void _rebin_min_counts (int *hist_index, double * min_bin_counts) /*{{{*/
{
   Hist_t *h = find_hist (*hist_index);
   if (-1 == Hist_do_rebin (h, Hist_rebin_min_counts,(void *) min_bin_counts))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "rebinning data set %d", *hist_index);
     }
}

/*}}}*/

static void _rebin_index (void) /*{{{*/
{
   Hist_t *h;
   SLang_Array_Type *sl_rebin_index = NULL;
   int hist_index, orig_nbins;

   if (-1 == SLang_pop_array_of_type (&sl_rebin_index, SLANG_INT_TYPE)
       || sl_rebin_index == NULL
       || -1 == SLang_pop_integer (&hist_index))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (NULL == (h = find_hist (hist_index)))
     goto finish;

   if (-1 == (orig_nbins = Hist_orig_hist_size (h)))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (sl_rebin_index->num_elements != (unsigned) orig_nbins)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto finish;
     }

   if (-1 == Hist_do_rebin (h, Hist_rebin_index, (void *) sl_rebin_index->data))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "rebinning data set %d", hist_index);
     }

   finish:
   SLang_free_array (sl_rebin_index);
}

/*}}}*/

static void _rebin_histogram (void) /*{{{*/
{
   SLang_Array_Type *inlo, *inhi, *inval, *outlo, *outhi, *outval;
   int n;

   inlo = inhi = inval = NULL;
   outlo = outhi = outval = NULL;

   if (-1 == SLang_pop_array_of_type (&inval, SLANG_DOUBLE_TYPE)
       || inval == NULL
       || -1 == SLang_pop_array_of_type (&inhi, SLANG_DOUBLE_TYPE)
       || inhi == NULL
       || -1 == SLang_pop_array_of_type (&inlo, SLANG_DOUBLE_TYPE)
       || inlo == NULL
       || -1 == SLang_pop_array_of_type (&outhi, SLANG_DOUBLE_TYPE)
       || outhi == NULL
       || -1 == SLang_pop_array_of_type (&outlo, SLANG_DOUBLE_TYPE)
       || outlo == NULL)
     {
        /* use slang wrapper to print usage message */
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   if ((inval->num_elements != inhi->num_elements)
       || (inval->num_elements != inlo->num_elements)
       || (outhi->num_elements != outlo->num_elements))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto free_and_return;
     }

   n = (int) outhi->num_elements;

   outval = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1);
   if (NULL == outval)
     {
        isis_throw_exception (Isis_Error);
        goto free_and_return;
     }

   rebin_histogram ((double *)inval->data,
                     (double *)inlo->data, (double *)inhi->data,
                     (int) inval->num_elements,
                     (double *)outval->data,
                     (double *)outlo->data, (double *)outhi->data,
                     (int) outval->num_elements);

   free_and_return:

   SLang_free_array (inlo);
   SLang_free_array (inhi);
   SLang_free_array (inval);
   SLang_free_array (outlo);
   SLang_free_array (outhi);
   SLang_push_array (outval, 1);
}

/*}}}*/

static void _rebin_array (void) /*{{{*/
{
   SLang_Array_Type *sl_rebin = NULL;
   SLang_Array_Type *sl_a = NULL;
   SLang_Array_Type *sl_ra = NULL;
   int *rebin, *pend, *p;
   double *a;
   unsigned int n;
   int sgn, nr;

   /* pop the data array and a matching array of rebin flags */
   if ((-1 == SLang_pop_array_of_type (&sl_rebin, SLANG_INT_TYPE))
       || (sl_rebin == NULL)
       || (-1 == SLang_pop_array_of_type (&sl_a, SLANG_DOUBLE_TYPE))
       || (sl_a == NULL)
       || (sl_a->num_elements != sl_rebin->num_elements))
     goto finish;

   n = sl_a->num_elements;
   a = (double *) sl_a->data;
   rebin = (int *) sl_rebin->data;

   pend = rebin + n;
   /* count the number of bins in the rebinned result */
   for (p = rebin; (p < pend) && (*p == 0); p++)
     ;
   if (p == pend)
     goto finish;

   nr = 1;
   sgn = *p;
   p++;

   for ( ; p < pend; p++)
     {
        if (*p == 0 || *p == sgn)
          continue;
        nr++;
        sgn = *p;
     }

   /* allocate space for the result */
   sl_ra = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nr, 1);
   if (sl_ra == NULL)
     goto finish;

   /* rebin, and we're done */
   if (-1 == apply_rebin (a, n, rebin, nr, (double *)sl_ra->data))
     goto finish;

   finish:
   SLang_free_array (sl_a);
   SLang_free_array (sl_rebin);
   SLang_push_array (sl_ra, 1);
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

   isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
   SLang_free_array (a);
   SLang_free_array (b);

   return -1;
}

/*}}}*/

static void _rebin_dataset (int *hist_index) /*{{{*/
{
   SLang_Array_Type *sl_lo, *sl_hi;
   double *lo, *hi;
   Hist_t *h;
   int nbins;

   if ((NULL == (h = find_hist (*hist_index)))
       || (-1 == pop_2_matched_arrays (SLANG_DOUBLE_TYPE, &sl_lo, &sl_hi)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   lo = (double *) sl_lo->data;
   hi = (double *) sl_hi->data;
   nbins = sl_lo->num_elements;

   if (-1 == Hist_rebin (h, lo, hi, nbins))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "rebinning dataset %d", *hist_index);
        return;
     }

   SLang_free_array (sl_lo);
   SLang_free_array (sl_hi);
}

/*}}}*/

/*}}}*/

/* Three arrays are expected with usage like:
 *    new_y = interpol (new_x, old_x, old_y);
 * <initial implementation borrowed from John Davis's jdl>
 */
static void array_interp (void) /*{{{*/
{
   SLang_Array_Type *x, *y, *new_x, *new_y;
   unsigned int new_npts, npts;

   if (SLang_Num_Function_Args != 3)
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   x = y = new_x = new_y = NULL;

   if (-1 == pop_2_matched_arrays (SLANG_DOUBLE_TYPE, &x, &y))
     return;

   if (-1 == SLang_pop_array_of_type (&new_x, SLANG_DOUBLE_TYPE))
     goto return_error;

   npts = x->num_elements;

   if (npts < 2)
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "Need at least two points for interpolation.");
        goto return_error;
     }

   if (NULL == (new_y = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL,
                                            new_x->dims, new_x->num_dims)))
     goto return_error;

   new_npts = new_x->num_elements;

   (void) interpolate_dvector ((double *)new_x->data, (double *)new_y->data, new_npts,
                               (double *)x->data, (double *)y->data, npts);

   SLang_push_array (new_y, 0);
   /* drop */

   return_error:

   SLang_free_array (x);
   SLang_free_array (y);
   SLang_free_array (new_x);
   SLang_free_array (new_y);
}

/*}}}*/

/* Three arrays are expected with usage like:
 *    new_y = interpol_points (new_x, old_x, old_y);
 */
static void array_interp_points (void) /*{{{*/
{
   SLang_Array_Type *x, *y, *new_x, *new_y;
   double *xx, *yy;
   unsigned int new_npts, npts, i;

   if (SLang_Num_Function_Args != 3)
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   x = y = new_x = new_y = NULL;

   if (-1 == pop_2_matched_arrays (SLANG_DOUBLE_TYPE, &x, &y))
     return;

   if (-1 == SLang_pop_array_of_type (&new_x, SLANG_DOUBLE_TYPE))
     goto return_error;

   npts = x->num_elements;

   if (npts < 2)
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "Need at least two points for interpolation.");
        goto return_error;
     }

   if (NULL == (new_y = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL,
                                            new_x->dims, new_x->num_dims)))
     goto return_error;

   new_npts = new_x->num_elements;
   xx = (double *)new_x->data;
   yy = (double *)new_y->data;

   for (i = 0; i < new_npts; i++)
     {
        yy[i] = interpolate_d (xx[i], (double *)x->data, (double *)y->data, npts);
     }

   SLang_push_array (new_y, 0);
   /* drop */

   return_error:

   SLang_free_array (x);
   SLang_free_array (y);
   SLang_free_array (new_x);
   SLang_free_array (new_y);
}

/*}}}*/

/*{{{ RMFs */

static Isis_Rmf_t *get_rmf_from_index (int i) /*{{{*/
{
   Isis_Rmf_t *r;

   if ((NULL == Rmf_List_Head)
       && (NULL == (Rmf_List_Head = Rmf_init_rmf_list())))
     {
        isis_throw_exception (Isis_Error);
        return NULL;
     }

   if (NULL == (r = Rmf_find_rmf_index (Rmf_List_Head, i)))
     {
        isis_throw_exception (Isis_Error);
        return NULL;
     }

   return r;
}

/*}}}*/

static Isis_Rmf_t *pop_rmf (void) /*{{{*/
{
   int i;

   if (-1 == SLang_pop_integer (&i))
     return NULL;

   return get_rmf_from_index (i);
}

/*}}}*/

static void release_rmf (Isis_Rmf_t *r) /*{{{*/
{
   (void) r;
}

/*}}}*/

static void _rebin_rmf (void) /*{{{*/
{
   SLang_Array_Type *lo, *hi;
   Isis_Rmf_t *rmf = NULL;

   lo = hi = NULL;

   if (-1 == pop_2_matched_arrays (SLANG_DOUBLE_TYPE, &lo, &hi))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (NULL == (rmf = pop_rmf ()))
     {
        SLang_free_array (lo);
        SLang_free_array (hi);
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == Rmf_rebin_rmf (rmf, (double *) lo->data, (double *) hi->data, lo->num_elements))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "rebinning RMF");
        /* drop */
     }

   SLang_free_array (lo);
   SLang_free_array (hi);
   release_rmf (rmf);
}

/*}}}*/

static Isis_Arf_t *_factor_rsp (Isis_Rmf_t *rmf) /*{{{*/
{
   Isis_Arf_t *arf;
   double *lo, *hi;
   unsigned int n;

   lo = hi = NULL;

   if (rmf == NULL)
     return NULL;

   if (rmf->ref_count > 0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "RMF in use, factorization failed");
        return NULL;
     }

   if (-1 == rmf->get_arf_grid (rmf, &lo, &hi, &n))
     return NULL;

   if (NULL == (arf = Arf_new_arf (n)))
     {
        ISIS_FREE(lo);
        ISIS_FREE(hi);
        return NULL;
     }

   memcpy ((char *)arf->bin_lo, (char *)lo, n * sizeof (double));
   memcpy ((char *)arf->bin_hi, (char *)hi, n * sizeof (double));

   ISIS_FREE(lo);
   ISIS_FREE(hi);

   arf->nbins = n;
   arf->order = rmf->order;
   isis_strcpy (arf->grating, rmf->grating, ISIS_ARF_VALUE_SIZE);
   isis_strcpy (arf->instrument, rmf->instrument, ISIS_ARF_VALUE_SIZE);

   if (-1 == rmf->factor_rsp (rmf, arf->arf))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "factoring response matrix");
        Arf_free_arf (arf);
        return NULL;
     }

   memset ((char *)arf->arf_err, 0, n * sizeof(double));

   return arf;
}

/*}}}*/

static int factor_rsp (int *rmf_index) /*{{{*/
{
   Isis_Rmf_t *r;
   Isis_Arf_t *a;

   if (NULL == (r = get_rmf_from_index (*rmf_index)))
     return -1;

   if (NULL == Arf_List_Head
       && NULL == (Arf_List_Head = Arf_init_arf_list()))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (NULL == (a = _factor_rsp (r)))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return Arf_append_arf (Arf_List_Head, a);
}

/*}}}*/

typedef struct
{
   SLang_Array_Type *bin_lo;
   SLang_Array_Type *bin_hi;
}
Grid_Type;

static SLang_CStruct_Field_Type Grid_Type_Layout[] =
{
   MAKE_CSTRUCT_FIELD (Grid_Type, bin_lo, "bin_lo", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Grid_Type, bin_hi, "bin_hi", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int get_rmf_grid (int rmf_index, /*{{{*/
                         int (*get_grid_fun)(Isis_Rmf_t *, double **, double **, unsigned int *))
{
   Grid_Type g;
   Isis_Rmf_t *rmf;
   double *lo, *hi;
   unsigned int n;
   int num, status = -1;

   memset ((char *)&g, 0, sizeof g);

   if (get_grid_fun == NULL)
     goto return_error;

   if (NULL == (rmf = Rmf_find_rmf_index (Rmf_List_Head, rmf_index)))
     goto return_error;

   if (-1 == (*get_grid_fun)(rmf, &lo, &hi, &n))
     goto return_error;

   num = n;

   g.bin_lo = SLang_create_array (SLANG_DOUBLE_TYPE, 0, lo, &num, 1);
   g.bin_hi = SLang_create_array (SLANG_DOUBLE_TYPE, 0, hi, &num, 1);

   status = 0;

   return_error:
   SLang_push_cstruct ((VOID_STAR)&g, Grid_Type_Layout);
   return status;
}

/*}}}*/

static void _get_rmf_data_grid (int *rmf_index) /*{{{*/
{
   (void) get_rmf_grid (*rmf_index, &Rmf_get_data_grid);
}

/*}}}*/

static void _get_rmf_arf_grid (int *rmf_index) /*{{{*/
{
   (void) get_rmf_grid (*rmf_index, &Rmf_get_arf_grid);
}

/*}}}*/

static void _find_rmf_peaks (int *rmf_index) /*{{{*/
{
   SLang_Array_Type *sl_peaks = NULL;
   Isis_Rmf_t *rmf;
   double *h_P = NULL;
   int n;

   if (NULL == (rmf = Rmf_find_rmf_index (Rmf_List_Head, *rmf_index)))
     goto return_error;

   if (-1 == Rmf_find_peaks(rmf, &h_P, &n))
     goto return_error;

   sl_peaks = SLang_create_array (SLANG_DOUBLE_TYPE, 0, h_P, &n, 1);

   return_error:
   SLang_push_array (sl_peaks, 1);
}

/*}}}*/

/*}}}*/

/*{{{ Slang intrinsic tables */

#define V SLANG_VOID_TYPE
#define S SLANG_STRING_TYPE
#define I SLANG_INT_TYPE
#define UI SLANG_UINT_TYPE
#define D SLANG_DOUBLE_TYPE

static SLang_Intrin_Fun_Type Hist_Intrinsics [] =
{
   MAKE_INTRINSIC_0("_rebin_rmf", _rebin_rmf, V),
   MAKE_INTRINSIC_I("_factor_rsp", factor_rsp, I),
   MAKE_INTRINSIC_2("_copy_hist_keywords", _copy_hist_keywords, V, I, I),
   MAKE_INTRINSIC("_rebin_histogram", _rebin_histogram, V, 0),
   MAKE_INTRINSIC_S("_load_arf", load_arf, I),
   MAKE_INTRINSIC_2("_assign_arf_to_hist", _assign_arf_to_hist, V, I, I),
   MAKE_INTRINSIC_I("_delete_arf", _delete_arf, V),
   MAKE_INTRINSIC_I("_get_arf", get_arf, V),
   MAKE_INTRINSIC("_put_arf", put_arf, V, 0),
   MAKE_INTRINSIC_1("_load_user_rmf", load_user_rmf, I, S),
   MAKE_INTRINSIC_1("_load_file_rmf", load_file_rmf, I, S),
   MAKE_INTRINSIC_0("_load_slang_rmf", load_slang_rmf, I),
   MAKE_INTRINSIC_2("_assign_rmf_to_hist", _assign_rmf_to_hist, V, I, I),
   MAKE_INTRINSIC_I("_delete_rmf", _delete_rmf, V),
   MAKE_INTRINSIC_I("_get_rmf_data_grid", _get_rmf_data_grid, V),
   MAKE_INTRINSIC_I("_get_rmf_arf_grid", _get_rmf_arf_grid, V),
   MAKE_INTRINSIC_I("_find_rmf_peaks", _find_rmf_peaks, V),
   MAKE_INTRINSIC_I("_assign_rsp_list", _assign_rsp_list, V),
   MAKE_INTRINSIC_I("get_flux_corr_weights_intrin", get_flux_corr_weights_intrin, V),
   MAKE_INTRINSIC_2("_flux_correct", _flux_correct, V, I, D),
   MAKE_INTRINSIC_2("_flux_correct_model_counts", _flux_correct_model_counts, V, I, D),
   MAKE_INTRINSIC_2("_rebin_min_counts", _rebin_min_counts, V, I, D),
   MAKE_INTRINSIC("_rebin_index", _rebin_index, V, 0),
   MAKE_INTRINSIC("_rebin_array", _rebin_array, V, 0),
   MAKE_INTRINSIC_I("_rebin_dataset", _rebin_dataset, V),
   MAKE_INTRINSIC("_array_interp", array_interp, V, 0),
   MAKE_INTRINSIC("_array_interp_points", array_interp_points, V, 0),
   MAKE_INTRINSIC_1("_set_stat_error_hook", set_stat_error_hook, V, I),
   MAKE_INTRINSIC_1("_get_stat_error_hook", get_stat_error_hook, V, I),
   MAKE_INTRINSIC_2("_load_data", _load_data, V, S, I),
   MAKE_INTRINSIC_I("_delete_hist", _delete_hist, V),
   MAKE_INTRINSIC_2("_get_hist", get_hist, V, I, UI),
   MAKE_INTRINSIC_2("_get_hist_notice_info", get_hist_notice_info, V, I, UI),
   MAKE_INTRINSIC("_put_hist", put_hist, V, 0),
   MAKE_INTRINSIC("put_model_intrin", put_model_intrin, V, 0),
   MAKE_INTRINSIC("_define_hist", _define_hist, I, 0),
   MAKE_INTRINSIC("_define_arf", _define_arf, I, 0),
   MAKE_INTRINSIC_2("_define_bgd", define_bgd, I, I,D),
   MAKE_INTRINSIC_2("_define_bgd_file", define_bgd_file, I, I,S),
   MAKE_INTRINSIC_I("_get_back_exposure", _get_back_exposure, V),
   MAKE_INTRINSIC_2("_set_back_exposure", _set_back_exposure, V, I,D),
   MAKE_INTRINSIC_I("_get_histogram_exposure_time", _get_histogram_exposure_time, V),
   MAKE_INTRINSIC_2("_set_histogram_exposure_time", _set_histogram_exposure_time, V, I, D),
   MAKE_INTRINSIC_I("_get_data_region_area", _get_data_region_area, V),
   MAKE_INTRINSIC_1("_set_data_region_area", _set_data_region_area, V, I),
   MAKE_INTRINSIC_I("_get_back_region_area", _get_back_region_area, V),
   MAKE_INTRINSIC_1("_set_back_region_area", _set_back_region_area, V, I),
   MAKE_INTRINSIC_2("_set_instrument_background", _set_instrument_background, V, I, S),
   MAKE_INTRINSIC_1("_get_instrumental_background_hook_name", _get_instrumental_background_hook_name, V, I),
   MAKE_INTRINSIC_1("_all_data", _all_data, V, I),
   MAKE_INTRINSIC("_all_arfs", _all_arfs, V, 0),
   MAKE_INTRINSIC("_all_rmfs", _all_rmfs, V, 0),
   MAKE_INTRINSIC_1("_is_fake_data", is_fake_data, I, I),
   MAKE_INTRINSIC_2("_set_fake", set_fake, I, I, I),
   MAKE_INTRINSIC_1("_is_grouped_data", is_grouped_data, I, I),
   MAKE_INTRINSIC_1("have_data", have_data, I, I),
   MAKE_INTRINSIC_I("_get_hist_info", _get_hist_info, I),
   MAKE_INTRINSIC_I("_set_hist_info", _set_hist_info, V),
   MAKE_INTRINSIC_I("_get_rmf_info", _get_rmf_info, V),
   MAKE_INTRINSIC_I("_set_rmf_info", _set_rmf_info, V),
   MAKE_INTRINSIC_I("_get_arf_info", _get_arf_info, I),
   MAKE_INTRINSIC_I("_set_arf_info", _set_arf_info, V),
   MAKE_INTRINSIC_1("_get_arf_exposure_time", _get_arf_exposure_time, D, I),
   MAKE_INTRINSIC_2("_set_arf_exposure_time", _set_arf_exposure_time, V, I, D),
   MAKE_INTRINSIC_2("_set_hist_frame_time", _set_hist_frame_time, V, I, D),
   MAKE_INTRINSIC_1("_get_hist_frame_time", _get_hist_frame_time, D, I),
   MAKE_INTRINSIC_1("_get_kernel", _get_kernel, V, I),
   MAKE_INTRINSIC_5("_plot_hist", _plot_hist, V, I, UI, I, I, I),
   MAKE_INTRINSIC_2("_set_exclude_flag", _set_exclude_flag, V, I, I),
   MAKE_INTRINSIC_4("_set_notice", set_notice, V, I, D, D, I),
   MAKE_INTRINSIC_I("_set_notice_using_mask", set_notice_using_mask, V),
   MAKE_INTRINSIC_II("_set_notice_using_list", set_notice_using_list, V),
   MAKE_INTRINSIC_I("_ignore_bad", _ignore_bad, V),
   MAKE_INTRINSIC_6("_get_region_sum", _get_region_sum, V, I, UI, D, D, D, D),
   MAKE_INTRINSIC_4("_cursor_region_stats", _cursor_region_stats, V, I, I, I, S),
   MAKE_INTRINSIC_II("_set_hist_color", _set_hist_color, V),
   MAKE_INTRINSIC_I("_unset_hist_color", _unset_hist_color, V),
   MAKE_INTRINSIC_2("_set_hist_object_name", _set_hist_object_name, V, I, S),
   MAKE_INTRINSIC("_set_post_model_hook", set_post_model_hook, V, 0),
   MAKE_INTRINSIC("_set_pre_combine_hook", set_pre_combine_hook, V, 0),
   MAKE_INTRINSIC_I("_assign_model", assign_model_intrin, V),
   MAKE_INTRINSIC_I("_set_dataset_metadata", set_dataset_metadata, V),
   MAKE_INTRINSIC_I("_get_dataset_metadata", get_dataset_metadata, V),
   MAKE_INTRINSIC_I("_set_sys_err_frac_intrin", set_sys_err_frac_intrin, V),
   MAKE_INTRINSIC_I("_get_sys_err_frac_intrin", get_sys_err_frac_intrin, V),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef V
#undef S
#undef I
#undef UI
#undef D

static SLang_IConstant_Type Hist_Intrin_Const [] =
{
   MAKE_ICONSTANT("_Flux", H_FLUX),
   MAKE_ICONSTANT("_Data", H_DATA),
   MAKE_ICONSTANT("_Convolved", H_CONVOLVED),
   MAKE_ICONSTANT("RMF_DELTA", RMF_DELTA),
   MAKE_ICONSTANT("RMF_FILE", RMF_FILE),
   MAKE_ICONSTANT("RMF_USER", RMF_USER),
   MAKE_ICONSTANT("RMF_SLANG", RMF_SLANG),
   SLANG_END_ICONST_TABLE
};

static SLang_IConstant_Type Hist_Pub_Intrin_Const [] =
{
   MAKE_ICONSTANT("DELCHI", ISIS_STAT_RESID),  /* kept for back-compatibility */
   MAKE_ICONSTANT("STAT", ISIS_STAT_RESID),
   MAKE_ICONSTANT("DIFF", ISIS_DIFF_RESID),
   MAKE_ICONSTANT("RATIO", ISIS_RATIO_RESID),
   MAKE_ICONSTANT("SEPARATE_GRID", ISIS_EVAL_GRID_SEPARATE),
   MAKE_ICONSTANT("MERGED_GRID", ISIS_EVAL_GRID_MERGED),
   MAKE_ICONSTANT("USER_GRID", ISIS_EVAL_GRID_USER),
   SLANG_END_ICONST_TABLE
};

#define UI SLANG_UINT_TYPE
#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE

static SLang_Intrin_Var_Type Hist_Intrin_Vars [] =
{
   MAKE_VARIABLE("Label_By_Default", &Label_By_Default, UI, 0),
   MAKE_VARIABLE("Isis_Residual_Plot_Type", &Isis_Residual_Plot_Type, I, 0),
   MAKE_VARIABLE("Isis_List_Filenames", &Isis_List_Filenames, I, 0),
   MAKE_VARIABLE("Ignore_PHA_Response_Keywords", &Hist_Ignore_PHA_Response_Keywords, I, 0),
   MAKE_VARIABLE("Ignore_PHA_Backfile_Keyword", &Hist_Ignore_PHA_Backfile_Keyword, I, 0),
   MAKE_VARIABLE("Allow_Multiple_Arf_Factors", &Hist_Allow_Multiple_Arf_Factors, I, 0),
   MAKE_VARIABLE("Warn_Invalid_Uncertainties", &Hist_Warn_Invalid_Uncertainties, I, 0),
   MAKE_VARIABLE("Rmf_Grid_Tol", &Hist_Rmf_Grid_Match_Tol, D, 0),
   MAKE_VARIABLE("Rmf_OGIP_Compliance", &Isis_Rmf_OGIP_Compliance, I, 0),
   MAKE_VARIABLE("Minimum_Stat_Err", &Hist_Min_Stat_Err, D, 0),
   MAKE_VARIABLE("Min_Model_Spacing", &Hist_Min_Model_Spacing, D, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

#undef UI
#undef D
#undef S
#undef I

/*}}}*/

static void quit_hist (void) /*{{{*/
{
   Label_By_Default = 1;
   Hist_free_list (Data_List_Head); Data_List_Head=NULL;
   Arf_free_arf_list (Arf_List_Head); Arf_List_Head=NULL;
   Rmf_free_rmf_list (Rmf_List_Head); Rmf_List_Head=NULL;
}

/*}}}*/

SLANG_MODULE(hist);
int init_hist_module_ns (char *ns_name) /*{{{*/
{
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if ((-1 == SLns_add_intrin_fun_table (ns, Hist_Intrinsics, NULL))
       || (-1 == SLns_add_iconstant_table (ns, Hist_Intrin_Const, NULL))
       || (-1 == SLns_add_iconstant_table (pub_ns, Hist_Pub_Intrin_Const, NULL))
       || (-1 == SLns_add_intrin_var_table (pub_ns, Hist_Intrin_Vars, NULL)))
     return isis_trace_return(-1);

   return 0;
}

/*}}}*/

void deinit_hist_module (void) /*{{{*/
{
   quit_hist();
}

/*}}}*/
