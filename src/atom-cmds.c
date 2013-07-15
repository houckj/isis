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

/* $Id: atom-cmds.c,v 1.47 2004/09/10 02:36:58 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <slang.h>
#include <errno.h>

#include "isis.h"
#include "plot.h"
#include "util.h"
#include "db-atomic.h"
#include "db-em.h"
#include "db-display.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

typedef struct
{
   DB_t *ptr;
}
DB_MMT_Type;
static int DB_MMT_Type_Id = -1;

static void destroy_atomic_db_mmt_type (SLtype type, VOID_STAR f) /*{{{*/
{
   DB_MMT_Type *mt = (DB_MMT_Type *)f;
   (void) type;

   if (mt != NULL)
     {
        DB_end (mt->ptr);
     }

   SLfree ((char *)f);
}

/*}}}*/

static SLang_MMT_Type *create_atomic_db_mmt_type (DB_t *ptr) /*{{{*/
{
   SLang_MMT_Type *mmt;
   DB_MMT_Type *mt;

   if (ptr == NULL)
     return NULL;

   if (NULL == (mt = (DB_MMT_Type *)SLmalloc (sizeof *mt)))
     return NULL;

   mt->ptr = ptr;

   mmt = SLang_create_mmt (DB_MMT_Type_Id, (VOID_STAR) mt);
   if (NULL == mmt)
     {
        SLfree ((char *)mt);
        return NULL;
     }

   return mmt;
}

/*}}}*/

static void push_atomic_db_pointer_intrin (DB_t *ptr) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;

   if (NULL == (mmt = create_atomic_db_mmt_type (ptr)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == SLang_push_mmt (mmt))
     SLang_free_mmt (mmt);
}

/*}}}*/

static DB_t *pop_atomic_db_pointer_intrin (void) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   DB_MMT_Type *mt;

   SLang_run_hooks ("_isis->get_atomic_db_pointer", 0, NULL);

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLang_pop_null();
        return NULL;
     }

   if (NULL == (mmt = SLang_pop_mmt (DB_MMT_Type_Id)))
     return NULL;

   mt = (DB_MMT_Type *) SLang_object_from_mmt (mmt);
   SLang_free_mmt (mmt);

   return mt->ptr;
}

/*}}}*/

/*{{{ utilities */

/* null terminated pointer array */

static char **copy_string_array (SLang_Array_Type *s) /*{{{*/
{
   char **a;
   SLindex_Type i;

   if (NULL == (a = (char **) ISIS_MALLOC ((s->num_elements + 1) * sizeof(char *))))
     return NULL;

   for (i = 0; i < (SLindex_Type) s->num_elements; i++)
     {
        if (-1 == SLang_get_array_element (s, &i, (VOID_STAR) &a[i]))
          {
             ISIS_FREE (a);
             return NULL;
          }
     }

   /* ensure NULL terminated */
   a[s->num_elements] = NULL;

   return a;
}

/*}}}*/

static DB_t *try_loading_filemap (char *filemap) /*{{{*/
{
   DB_t *db = NULL;
   SLang_Array_Type *sl_e = NULL;
   SLang_Array_Type *sl_w = NULL;
   char **elev_files = NULL;
   char **wavelen_files = NULL;

   if ((filemap == NULL) || (*filemap == 0))
     goto finish;

   if (1 == SLang_run_hooks ("_isis->atomdb_start_hook", 1, filemap))
     {
        if (-1 == SLang_pop_array_of_type (&sl_e, SLANG_STRING_TYPE)
            || -1 == SLang_pop_array_of_type (&sl_w, SLANG_STRING_TYPE))
          {
             SLang_free_array (sl_e);
             SLang_free_array (sl_w);
             goto finish;
          }

        /* NULL ok */
        elev_files = copy_string_array (sl_e);
        wavelen_files = copy_string_array (sl_w);
     }

   finish:

   /* NULL ok */
   db = DB_start (elev_files, wavelen_files);

   SLang_free_array (sl_e);
   SLang_free_array (sl_w);
   ISIS_FREE (elev_files);
   ISIS_FREE (wavelen_files);

   return db;
}

/*}}}*/
#ifdef __cplusplus
extern "C"
#endif
void deinit_emis_module (void);
static void db_start (char *filemap) /*{{{*/
{
   DB_t *db;

   isis_vmesg (WARN, I_INITIALIZING, __FILE__, __LINE__, "atomic data");

   db = try_loading_filemap (filemap);
   if (NULL == db)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "Atomic database not initialized");
        return;
     }

   if (NULL != _ptr_to_emissivity_db ())
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__,
                    "Replaced internal data structures: please reload emissivity tables");
        (void) deinit_emis_module ();
     }

   push_atomic_db_pointer_intrin (db);
}

/*}}}*/

static void quit_atom (void) /*{{{*/
{
   DB_Ion_Format = FMT_ROMAN;
}

/*}}}*/

/* the silent version */
DB_t *_ptr_to_atomic_db (void) /*{{{*/
{
   return pop_atomic_db_pointer_intrin ();
}

/*}}}*/

DB_t *ptr_to_atomic_db (void) /*{{{*/
{
   DB_t *db = pop_atomic_db_pointer_intrin ();
   if (NULL == db)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "atomic database not loaded");
        isis_throw_exception (Isis_Error);
     }
   return db;
}

/*}}}*/

/*}}}*/

/*{{{ Line group definition */

static void get_unblended (void) /*{{{*/
{
   SLang_Array_Type *at = NULL, *sllist = NULL;
   DB_t *db = ptr_to_atomic_db ();
   int n, *t = NULL;
   unsigned int type;
   float f, wl_f;

   if (NULL == db)
     return;

   if (-1 == SLang_pop_array_of_type (&sllist, SLANG_INT_TYPE)
       || sllist == NULL
       || sllist->num_elements == 0
       || -1 == SLang_pop_uinteger (&type)
       || -1 == SLang_pop_float (&wl_f)
       || -1 == SLang_pop_float (&f))
     goto fail;

   if (NULL == (t = (int *) ISIS_MALLOC (sllist->num_elements * sizeof(int))))
     goto fail;

   if (-1 == DB_get_unblended (t, &n, f, wl_f, type, (int *)sllist->data,
                                 sllist->num_elements, db))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "searching for unblended lines");

   if (NULL == (at = SLang_create_array (SLANG_INT_TYPE, 1, NULL, &n, 1)))
     goto fail;

   memcpy ((char *)at->data, (char *) t, n*sizeof(int));

   if (-1 == SLang_push_array (at, 1))
     goto fail;

   fail:
   ISIS_FREE (t);
   SLang_free_array (sllist);
}

/*}}}*/

static void get_k_brightest_lines (void) /*{{{*/
{
   SLang_Array_Type *at = NULL, *sllist = NULL;
   DB_t *db = ptr_to_atomic_db ();
   int n, *t = NULL;
   int k;

   if (NULL == db)
     return;

   if (-1 == SLang_pop_array_of_type (&sllist, SLANG_INT_TYPE)
       || sllist == NULL
       || -1 == SLang_pop_integer (&k)
       || k <= 0)
     goto fail;

   if (NULL == (t = (int *) ISIS_MALLOC (k * sizeof(int))))
     goto fail;

   if (-1 == DB_get_k_brightest (t, &n, k, (int *)sllist->data,
                                 sllist->num_elements, db))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "searching for %d brightest lines", k);

   if (NULL == (at = SLang_create_array (SLANG_INT_TYPE, 1, NULL, &n, 1)))
     goto fail;

   memcpy ((char *)at->data, (char *) t, n*sizeof(int));

   if (-1 == SLang_push_array (at, 1))
     goto fail;

   fail:
   ISIS_FREE (t);
   SLang_free_array (sllist);
}

/*}}}*/

typedef int (*set_filter_fcn_t)(DB_line_filter_t *, double, double);

static void filter_line_by_double_range (double xmin, double xmax, /*{{{*/
                                         set_filter_fcn_t set_filter_fcn,
                                         DB_line_filter_function_t filter_fcn)
{
   SLang_Array_Type *ct = NULL;
   DB_t *db = ptr_to_atomic_db ();
   DB_line_filter_t *filter = NULL;
   int n;

   if (NULL == db
       || NULL == (filter = DB_new_line_filter()))
     return;

   n = DB_get_nlines (db);
   if (n <= 0)
     {
        isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "no lines in database");
        goto fail;
     }

   /* will push read-only array onto stack */
   if (NULL == (ct = SLang_create_array (SLANG_CHAR_TYPE, 1, NULL, &n, 1)))
     goto fail;

   if (-1 == (*set_filter_fcn) (filter, xmin, xmax)
       || -1 == DB_apply_line_filter ( (char *)ct->data, filter_fcn, filter, db)
       || -1 == SLang_push_array (ct, 1))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "applying line filter");

   fail:
   DB_free_line_filter (filter);
}

/*}}}*/

static void filter_line_by_flux (double *flux_min, double *flux_max) /*{{{*/
{
   filter_line_by_double_range (*flux_min, *flux_max,
                                DBf_set_flux, DBf_flux);
}

/*}}}*/

static void filter_line_by_wavelength (double *wl_min, double *wl_max) /*{{{*/
{
   filter_line_by_double_range (*wl_min, *wl_max,
                                DBf_set_wavelength, DBf_wavelength);
}

/*}}}*/

static int pop_two_int_arrays (SLang_Array_Type **a, SLang_Array_Type **b) /*{{{*/
{
   if (-1 == SLang_pop_array_of_type (b, SLANG_INT_TYPE)
       || -1 == SLang_pop_array_of_type (a, SLANG_INT_TYPE))
     return -1;

   return 0;
}

/*}}}*/

static void filter_line_by_el_ion (void) /*{{{*/
{
   SLang_Array_Type *ct, *slZ, *slq;
   DB_t *db = ptr_to_atomic_db ();
   DB_line_filter_t *filter = NULL;
   int nions, nelem, *proton_number, *ion_charge;
   int n;

   ct = slZ = slq = NULL;

   if (NULL == db)
     return;

   if ((-1 == pop_two_int_arrays (&slZ, &slq))
       || (slZ == NULL) || (slq == NULL))
     goto fail;

   ion_charge = (int *) slq->data;
   nions = slq->num_elements;
   if (ion_charge[0] < 0)
     nions = 0;

   proton_number = (int *) slZ->data;
   nelem = slZ->num_elements;
   if (proton_number[0] <= 0)
     nelem = 0;

   if (NULL == (filter = DB_new_line_filter()))
     goto fail;

   n = DB_get_nlines (db);
   if (n <= 0)
     {
        isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "no lines in database");
        goto fail;
     }

   /* will push read-only array onto stack */
   if (NULL == (ct = SLang_create_array (SLANG_CHAR_TYPE, 1, NULL, &n, 1)))
     goto fail;

   if (-1 == DBf_set_el_ion (filter, nelem, proton_number, nions, ion_charge)
       || -1 == DB_apply_line_filter ((char *)ct->data, DBf_el_ion, filter, db)
       || -1 == SLang_push_array (ct, 1))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "applying line filter");

   fail:

   DB_free_line_filter (filter);
   SLang_free_array (slZ);
   SLang_free_array (slq);

   return;
}

/*}}}*/

static void filter_line_by_trans (int *Z, int *q) /*{{{*/
{
   SLang_Array_Type *ct, *sl_up, *sl_lo;
   DB_t *db = ptr_to_atomic_db ();
   DB_line_filter_t *filter = NULL;
   int nup, nlo, n;
   int *up, *lo;

   ct = sl_up = sl_lo = NULL;

   if (NULL == db)
     return;

   if ((-1 == pop_two_int_arrays (&sl_up, &sl_lo))
       || (sl_lo == NULL) || (sl_lo == NULL))
     goto fail;

   up = (int *) sl_up->data;
   nup = sl_up->num_elements;
   if (up[0] < 0)
     nup = 0;

   lo = (int *) sl_lo->data;
   nlo = sl_lo->num_elements;
   if (lo[0] <= 0)
     nlo = 0;

   if (NULL == (filter = DB_new_line_filter()))
     goto fail;

   n = DB_get_nlines (db);
   if (n <= 0)
     {
        isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "no lines in database");
        goto fail;
     }

   /* will push read-only array onto stack */
   if (NULL == (ct = SLang_create_array (SLANG_CHAR_TYPE, 1, NULL, &n, 1)))
     goto fail;

   if (-1 == DBf_set_trans (filter, *Z, *q, nup, up, nlo, lo)
       || -1 == DB_apply_line_filter ((char *)ct->data, DBf_trans, filter, db)
       || -1 == SLang_push_array (ct, 1))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "applying line filter");

   fail:

   DB_free_line_filter (filter);
   SLang_free_array (sl_up);
   SLang_free_array (sl_lo);

   return;
}

/*}}}*/

static void make_group_from_list (void) /*{{{*/
{
   SLang_Array_Type *at = NULL;
   DB_t *db = ptr_to_atomic_db ();
   int n, *list;
   int group;

   if (NULL == db)
     return;

   if (-1 == SLang_pop_array_of_type (&at, SLANG_INT_TYPE)
       || at == NULL
       || -1 ==  SLang_pop_integer (&group))
     {
        isis_throw_exception (Isis_Error);
        goto fail;
     }

   list = (int *) at->data;
   n = at->num_elements;

   if (-1 == DB_list_to_group (group, list, n, db))
     isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "undefined group (%d)", group);

   fail:
   SLang_free_array (at);
}

/*}}}*/

static void name_group (int * group, char *name) /*{{{*/
{
   DB_line_group_t *cl;
   DB_t *db = ptr_to_atomic_db ();

   if (NULL == db)
     return;

   if (NULL == (cl = DB_find_group (*group, db))
       || -1 == DB_set_group_name (name, cl))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting group name");
}

/*}}}*/

/*}}}*/

/*{{{ Manipulating line groups:  page/delete/list/plot */

static void page_group (int * group) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();
   DB_line_group_t *g;
   FILE * fp;

   if (db == NULL || group == NULL)
     return;

   if (NULL == (g = DB_find_group (*group, db)))
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "undefined group (%d)", *group);
        return;
     }

   if (NULL == (fp = isis_open_pager ()))
     return;

   (void) DB_print_line_group (fp, g, db);

   isis_close_pager (fp);
}

/*}}}*/

static void page_line_list (void) /*{{{*/
{
   SLang_Array_Type *at = NULL;
   DB_t *db = ptr_to_atomic_db ();
   DB_line_group_t *g = NULL;
   FILE * fp;

   if (db == NULL)
     return;

   if (-1 == SLang_pop_array_of_type (&at, SLANG_INT_TYPE))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (NULL == (g = DB_make_group_from_list (0, (int *)at->data,
                                             at->num_elements,
                                             db)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "generating line group from list");
        goto finish;
     }

   if (NULL == (fp = isis_open_pager()))
     goto finish;

   (void) DB_print_line_group (fp, g, db);

   isis_close_pager (fp);

   finish:
   DB_free_line_group (g);
   SLang_free_array (at);
}

/*}}}*/

static void _save_group (int *group, char *fname) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();
   DB_line_group_t *g;
   FILE *fp = NULL;

   if (db == NULL)
     return;

   if (NULL == (g = DB_find_group (*group, db)))
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "undefined group (%d)", *group);
        return;
     }

   if (NULL == fname)
     return;

   errno = 0;
   fp = fopen (fname, "a");
   if (NULL == fp)
     {
        isis_vmesg (FAIL, I_WRITE_OPEN_FAILED, __FILE__, __LINE__, "%s", fname);
        return;
     }

   if (-1 == DB_print_line_group (fp, g, db))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "saving group %d", *group);

   isis_fclose(fp);
}

/*}}}*/

static void _save_line_list (void) /*{{{*/
{
   SLang_Array_Type *at = NULL;
   DB_line_group_t *g = NULL;
   DB_t *db = ptr_to_atomic_db ();
   FILE *fp = NULL;
   char *fname = NULL;

   if (db == NULL)
     return;

   if (-1 == SLpop_string (&fname)
       || -1 == SLang_pop_array_of_type (&at, SLANG_INT_TYPE))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (NULL == (g = DB_make_group_from_list (0, (int *)at->data,
                                             at->num_elements,
                                             db)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "generating line group from list");
        goto finish;
     }

   if (NULL == fname)
     goto finish;

   errno = 0;
   fp = fopen (fname, "a");
   if (NULL == fp)
     {
        isis_vmesg (FAIL, I_WRITE_OPEN_FAILED, __FILE__, __LINE__, "%s", fname);
        goto finish;
     }

   if (-1 == DB_print_line_group (fp, g, db))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "saving line group");

   finish:
   DB_free_line_group (g);
   SLang_free_array (at);
   SLfree (fname);

   isis_fclose(fp);
}

/*}}}*/

static void delete_group (int * group) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();

   if (db != NULL && group != NULL)
     DB_delete_line_group (*group, db);
}

/*}}}*/

static void list_group (void) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();

   if (db != NULL)
     DB_list_line_group_stats (stdout, db);
}

/*}}}*/

static int pop_line_label_style (Line_Label_Style_Type *s) /*{{{*/
{
   if (s == NULL) return -1;

   s->label_length = SHORT_LINE_NAME_SIZE;

   if (-1 == SLang_pop_integer (&s->label_type)
       || -1 == SLang_pop_float (&s->char_height)
       || -1 == SLang_pop_float (&s->offset)
       || -1 == SLang_pop_float (&s->bottom_frac)
       || -1 == SLang_pop_float (&s->top_frac)
       || -1 == SLang_pop_float (&s->justify)
       || -1 == SLang_pop_float (&s->angle))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting label format parameters");
        return -1;
     }

   return 0;
}

/*}}}*/

static char *Linelabel_Hook = "isis_linelabel_hook";

static char *massage_label (char *s) /*{{{*/
{
   char *ms = NULL;

   if (1 == SLang_run_hooks (Linelabel_Hook, 1, s))
     {
        if (-1 == SLpop_string (&ms))
          return NULL;
     }

   return ms;
}

/*}}}*/

static void plot_line_group (void) /*{{{*/
{
   Line_Label_Style_Type s;
   DB_t *db = ptr_to_atomic_db ();
   DB_line_group_t *g;
   Plot_t *fmt = NULL;
   int group, use_color;
   float redshift;
   char *(*label_massage_hook)(char *);

   if (db == NULL)
     return;

   if (-1 == pop_line_label_style (&s))
     return;

   if (-1 == SLang_pop_integer (&use_color)
       || -1 == SLang_pop_float (&redshift)
       || -1 == SLang_pop_integer (&group))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (NULL == (g = DB_find_group (group, db)))
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "undefined group (%d)", group);
        return;
     }

   if (NULL == (fmt = current_format ()))
     return;

   if (2 == SLang_is_defined (Linelabel_Hook))
     label_massage_hook = massage_label;
   else
     label_massage_hook = NULL;

   s.use_color = (use_color > 0) ? &use_color : NULL;

   (void) DB_plot_line_group (g, &s, redshift, label_massage_hook, fmt, db);
   Plot_auto_incr_line_type (fmt);
}

/*}}}*/

static void plot_line_list (void) /*{{{*/
{
   SLang_Array_Type *at = NULL;
   Line_Label_Style_Type s;
   DB_t *db = ptr_to_atomic_db ();
   DB_line_group_t *g = NULL;
   Plot_t *fmt = NULL;
   int use_color;
   float redshift;
   char *(*label_massage_hook)(char *);

   if (db == NULL)
     return;

   if (-1 == pop_line_label_style (&s))
     return;

   if (-1 == SLang_pop_integer (&use_color)
       || -1 == SLang_pop_float (&redshift)
       || -1 == SLang_pop_array_of_type (&at, SLANG_INT_TYPE))
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (NULL == (g = DB_make_group_from_list (0, (int *)at->data,
                                             at->num_elements,
                                             db)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "generating line group from list");
        goto finish;
     }

   if (NULL == (fmt = current_format ()))
     goto finish;

   if (2 == SLang_is_defined (Linelabel_Hook))
     label_massage_hook = massage_label;
   else
     label_massage_hook = NULL;

   s.use_color = (use_color > 0) ? &use_color : NULL;

   (void) DB_plot_line_group (g, &s, redshift, label_massage_hook, fmt, db);
   Plot_auto_incr_line_type (fmt);

   finish:
   DB_free_line_group (g);
   SLang_free_array (at);
}

/*}}}*/

static void plot_line_list2 (void) /*{{{*/
{
   SLang_Array_Type *sl_lambdas = NULL;
   SLang_Array_Type *sl_labels = NULL;
   Line_Label_Style_Type s;
   Plot_t *fmt = NULL;
   float *lambdas = NULL;
   char **labels = NULL;
   char *label_text = NULL;
   SLindex_Type i, nlines;
   unsigned int label_size;
   int use_color;
   float redshift;

   if (-1 == pop_line_label_style (&s))
     return;

   if (-1 == SLang_pop_integer (&use_color)
       || -1 == SLang_pop_float (&redshift)
       || -1 == SLang_pop_array_of_type (&sl_labels, SLANG_STRING_TYPE)
       || sl_labels == NULL
       || -1 == SLang_pop_array_of_type (&sl_lambdas, SLANG_FLOAT_TYPE)
       || sl_lambdas == NULL
       || sl_labels->num_elements != sl_lambdas->num_elements)
     {
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   if (NULL == (fmt = current_format ()))
     goto finish;

   s.use_color = (use_color > 0) ? &use_color : NULL;

   lambdas = (float *)sl_lambdas->data;
   nlines = sl_lambdas->num_elements;

   label_size = s.label_length * sizeof(char);

   if ((NULL == (labels =(char **) ISIS_MALLOC (nlines * sizeof (char *))))
       || (NULL == (label_text = (char *) ISIS_MALLOC (nlines * label_size))))
     goto finish;

   memset ((char *)label_text, 0, nlines * label_size);
   for (i = 0; i < nlines; i++)
     {
        labels[i] = label_text + i * s.label_length;
        if (-1 == SLang_get_array_element (sl_labels, &i, (VOID_STAR) &labels[i]))
          goto finish;
     }

   (void) DB_plot_line_list (fmt, &s, redshift, lambdas, labels, nlines);
   Plot_auto_incr_line_type (fmt);

   finish:

   ISIS_FREE(labels);
   ISIS_FREE(label_text);
   SLang_free_array (sl_labels);
   SLang_free_array (sl_lambdas);
}

/*}}}*/

static void _group_members (int * group) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();
   SLang_Array_Type *at = NULL;
   int *indx = NULL;
   int n = 0;

   if (db == NULL)
     goto push;

   if (-1 == DB_get_group_members (&indx, &n, *group, db))
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "undefined group (%d)", *group);
        goto push;
     }

   if (NULL == (at = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &n, 1)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "creating index array");
        goto push;
     }

   memcpy ((char *)at->data, (char *)indx, n * sizeof(int));

   push:

   (void) SLang_push_array (at, 1);
   ISIS_FREE (indx);
}

/*}}}*/

#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define F SLANG_FLOAT_TYPE
#define D SLANG_DOUBLE_TYPE
#define S SLANG_STRING_TYPE
#define UC SLANG_UCHAR_TYPE

typedef struct
{
   double flux;
   float wavelen, wavelen_err;
   float A, A_err;
   int indx;
   int proton_number, ion_charge;
   int   upper_index, lower_index;
   float upper_g, lower_g;   /* statistical weights */
   float upper_E, lower_E;
   int   upper_L, upper_S;
   int   lower_L, lower_S;
   char *upper_label;
   char *lower_label;
}
Line_Info_Type;
#define LINE_INFO_TYPE_INIT {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0}

static SLang_CStruct_Field_Type Line_Info_Table [] =
{
   MAKE_CSTRUCT_FIELD(Line_Info_Type, flux, "flux", D, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, wavelen, "lambda", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, wavelen_err, "lambda_err", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, A, "A", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, A_err, "A_err", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, indx, "id", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, ion_charge, "ion", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, proton_number, "Z", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, upper_index, "upper", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, lower_index, "lower", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, upper_g, "upper_g", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, lower_g, "lower_g", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, upper_E, "upper_E", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, lower_E, "lower_E", F, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, upper_L, "upper_L", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, upper_S, "upper_S", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, lower_L, "lower_L", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, lower_S, "lower_S", I, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, upper_label, "up_name", S, 0),
   MAKE_CSTRUCT_FIELD(Line_Info_Type, lower_label, "lo_name", S, 0),
   SLANG_END_CSTRUCT_TABLE
};

#undef V
#undef I
#undef F
#undef D
#undef S
#undef UC

static void copy_line_info (Line_Info_Type *li, DB_line_t *p) /*{{{*/
{
   if (li == NULL || p == NULL)
     return;

   li->flux = p->flux;
   li->wavelen = p->wavelen;
   li->wavelen_err = p->wavelen_err;
   li->A = p->A;
   li->A_err = p->A_err;
   li->indx = p->indx;
   li->proton_number = p->proton_number;
   li->ion_charge = p->ion_charge;
   li->upper_index = p->upper_level;
   li->lower_index = p->lower_level;
}

/*}}}*/

static void _get_line_info (int *id) /*{{{*/
{
   Line_Info_Type li = LINE_INFO_TYPE_INIT;
   DB_t *db = ptr_to_atomic_db ();
   DB_line_t *line = NULL;
   DB_ion_t *ion = NULL;
   DB_level_t *e_upper, *e_lower;
   int Z, q, up, lo;
   char upper_label[LEVEL_NAME_SIZE];
   char lower_label[LEVEL_NAME_SIZE];

   if (NULL == db
       || NULL == (line = DB_get_line_from_index (*id, db)))
     goto finish;

   Z = line->proton_number;
   q = line->ion_charge;
   up = line->upper_level;
   lo = line->lower_level;

   copy_line_info (&li, line);

   (void) DB_get_level_label (upper_label, Z, q, up, db);
   (void) DB_get_level_label (lower_label, Z, q, lo, db);

   li.upper_label = upper_label;
   li.lower_label = lower_label;

   ion = DB_get_ion (db, Z, q);
   if (NULL != (e_upper = DB_get_ion_level (up, ion)))
     {
        li.upper_g = e_upper->stat_weight;
        li.upper_E = e_upper->energy;
        li.upper_L = e_upper->L;
        li.upper_S = e_upper->S;
     }

   if (NULL != (e_lower = DB_get_ion_level (lo, ion)))
     {
        li.lower_g = e_lower->stat_weight;
        li.lower_E = e_lower->energy;
        li.lower_L = e_lower->L;
        li.lower_S = e_lower->S;
     }

   finish:

   if (-1 == SLang_push_cstruct ((VOID_STAR)&li, Line_Info_Table))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't push cstruct");
     }
}

/*}}}*/

/*}}}*/

/*{{{ wavelength correction */

static void correct_line_wavelengths (char * filename) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();
   FILE *fp;
   char buf[BUFSIZE];
   int nchanged = 0;

   if (db == NULL)
     return;

   if (filename == NULL)
     return;

   if (NULL == (fp = fopen (filename, "r")))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", filename);
        return;
     }

   while (!feof(fp))
     {
        DB_line_t *line;
        float lambda, lambda_err;
        int proton_number, ion, upper_level, lower_level;

        if (NULL == fgets(buf, BUFSIZE, fp))
             break;

        if (buf[0] == COMMENT_CHAR)
          continue;

        if (6 != sscanf (buf, "%f %f %d %d %d %d\n",
                        &lambda, &lambda_err, &proton_number, &ion,
                        &upper_level, &lower_level))
          continue;

        line = DB_get_line_by_indices (proton_number, ion-1,
                                      upper_level, lower_level, db);
        if (NULL == line)
          {
             isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "line %s not found", buf);
             continue;
          }

        if (-1 == DB_set_line_wavelength (lambda, lambda_err, line))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "changing line %s", buf);
             continue;
          }

        nchanged++;
     }

   fclose(fp);

   if (0 == DB_sort_line_list (db))
     {
        fprintf (stderr,
"Changed %d wavelength values and re-sorted the line list.\n\
Note that these wavelength changes have not been\n\
propagated to any currently defined line groups.\n\
You might want to re-load the emissivity tables now.\n", nchanged);
     }
   else
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "changed %d wavelength values", nchanged);
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "resorting the line new list");
     }
}

/*}}}*/

/*}}}*/

/*{{{ display energy levels */

static void _plot_elev (void) /*{{{*/
{
   SLang_Array_Type *slline = NULL;
   DB_t *db = ptr_to_atomic_db ();
   Plot_t *fmt = NULL;
   int Z, q, subset, overlay;
   int *lines = NULL;
   int nlines;

   if (NULL == db)
     return;

   if (-1 == SLang_pop_integer (&overlay)
       || -1 == SLang_pop_integer (&subset)
       || -1 == SLang_pop_array_of_type (&slline, SLANG_INT_TYPE)
       || slline == NULL
       || -1 == SLang_pop_integer (&q)
       || -1 == SLang_pop_integer (&Z))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        goto fail;
     }

   nlines = slline->num_elements;
   lines = (int *) slline->data;

   if ((nlines == 1) && (lines[0] == -1))
     nlines = 0;

   if ((NULL == (fmt = current_format ())
        || (force_open_plot_device() < 0)))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        goto fail;
     }

   if (overlay == 0)
     Plot_restart_style_cycle (fmt);

   if (-1 == DB_plot_levels (Z, q, lines, nlines, subset, overlay, fmt, db))
     isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");

   fail:
   SLang_free_array (slline);
}

/*}}}*/

static void _oplot_transitions (void) /*{{{*/
{
   SLang_Array_Type *slline = NULL;
   DB_t *db = ptr_to_atomic_db ();
   Plot_t *fmt = NULL;
   int Z, q, style = 0;
   int *lines = NULL;
   int nlines;
   int *pstyle;

   if (NULL == db)
     return;

   if (-1 == SLang_pop_integer (&style)
       || -1 == SLang_pop_array_of_type (&slline, SLANG_INT_TYPE)
       || slline == NULL
       || -1 == SLang_pop_integer (&q)
       || -1 == SLang_pop_integer (&Z))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        goto fail;
     }

   nlines = slline->num_elements;
   lines = (int *) slline->data;

   if (nlines == 1 && lines[0] == -1)
     nlines = 0;

   if (NULL == (fmt = current_format ()))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        goto fail;
     }

   /* overrides auto-incr but won't override fmt style */
   pstyle = (style > 0) ? &style : NULL;

   if (-1 == DB_plot_transitions (Z, q, lines, nlines, pstyle, fmt, db))
     isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
   else
     Plot_auto_incr_line_type (fmt);

   fail:
   SLang_free_array (slline);
}

/*}}}*/

static void _list_elev (int *proton_number, int * ion_charge) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();
   FILE *fp;

   if (db == NULL)
     return;

   fp = isis_open_pager ();
   if ((fp != NULL)
       && (0 != DB_list_ion_levels (fp, *proton_number, *ion_charge, db)))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "listing ion levels");

   isis_close_pager (fp);
}

/*}}}*/

static void _list_branching (int *Z, int *q) /*{{{*/
{
   DB_t *db = ptr_to_atomic_db ();
   FILE *fp;

   if (db == NULL)
     return;

   fp = isis_open_pager ();

   if ((fp != NULL)
       && (0 != DB_print_branching_for_ion (fp, *Z, *q, db)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "listing branching ratios for Z=%d q=%d",
                    *Z, *q);
     }

   isis_close_pager (fp);
}

/*}}}*/

/*}}}*/

/*{{{ Slang intrinsic function table */

#define V SLANG_VOID_TYPE
#define F SLANG_FLOAT_TYPE
#define D SLANG_DOUBLE_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE

static SLang_IConstant_Type Atom_Intrin_Const [] =
{
   MAKE_ICONSTANT("SAME_ION", SAME_ION),
   MAKE_ICONSTANT("SAME_ELEM", SAME_ELEM),
   MAKE_ICONSTANT("FMT_ROMAN", FMT_ROMAN),
   MAKE_ICONSTANT("FMT_CHARGE", FMT_CHARGE),
   MAKE_ICONSTANT("FMT_INT_ROMAN", FMT_INT_ROMAN),
   SLANG_END_ICONST_TABLE
};

static SLang_Intrin_Var_Type Atom_Intrin_Vars [] =
{
   MAKE_VARIABLE("Ion_Format", &DB_Ion_Format, SLANG_INT_TYPE, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

#define DUMMY_DBOBJ_MMT_TYPE 255
#define MT DUMMY_DBOBJ_MMT_TYPE

static SLang_Intrin_Fun_Type Atom_Intrinsics [] =
{
   MAKE_INTRINSIC_S("_db_start", db_start, V),
   MAKE_INTRINSIC("make_group_from_list", make_group_from_list, V, 0),
   MAKE_INTRINSIC("filter_line_by_el_ion", filter_line_by_el_ion, V, 0),
   MAKE_INTRINSIC_2("filter_line_by_trans", filter_line_by_trans, V, I, I),
   MAKE_INTRINSIC_2("filter_line_by_flux", filter_line_by_flux, V, D, D),
   MAKE_INTRINSIC_2("filter_line_by_wavelength", filter_line_by_wavelength, V, D, D),
   MAKE_INTRINSIC("get_k_brightest_lines", get_k_brightest_lines, V, 0),
   MAKE_INTRINSIC("get_unblended", get_unblended, V, 0),
   MAKE_INTRINSIC_IS("_name_group", name_group, V),
   MAKE_INTRINSIC_I("_delete_group", delete_group, V),
   MAKE_INTRINSIC_I("_page_group", page_group, V),
   MAKE_INTRINSIC("_page_line_list", page_line_list, V, 0),
   MAKE_INTRINSIC_2("_save_group", _save_group, V, I, S),
   MAKE_INTRINSIC("_save_line_list", _save_line_list, V, 0),
   MAKE_INTRINSIC("_list_group", list_group, V, 0),
   MAKE_INTRINSIC("plot_line_group", plot_line_group, V, 0),
   MAKE_INTRINSIC("_plot_line_list", plot_line_list, V, 0),
   MAKE_INTRINSIC("_plot_line_list2", plot_line_list2, V, 0),
   MAKE_INTRINSIC_I("_group_members", _group_members, V),
   MAKE_INTRINSIC_I("_get_line_info", _get_line_info, V),
   MAKE_INTRINSIC_II("_list_elev", _list_elev, V),
   MAKE_INTRINSIC("_plot_elev", _plot_elev, V, 0),
   MAKE_INTRINSIC("_oplot_transitions", _oplot_transitions, V, 0),
   MAKE_INTRINSIC_S("_change_wl", correct_line_wavelengths, V),
   MAKE_INTRINSIC_II("_list_branching", _list_branching, V),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef MT
#undef V
#undef F
#undef D
#undef I
#undef S

SLANG_MODULE(atom);
int init_atom_module_ns (char *ns_name)
{
   SLang_Class_Type *cl;
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if((-1 == SLns_add_intrin_fun_table (ns, Atom_Intrinsics, NULL))
      || (-1 == SLns_add_iconstant_table (pub_ns, Atom_Intrin_Const, NULL))
      || (-1 == SLns_add_intrin_var_table (pub_ns, Atom_Intrin_Vars, NULL)))
     return isis_trace_return(-1);

   if (DB_MMT_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("DB_Object_Type")))
          return isis_trace_return(-1);
        (void) SLclass_set_destroy_function (cl, destroy_atomic_db_mmt_type);
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (DB_MMT_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return isis_trace_return(-1);
        DB_MMT_Type_Id = SLclass_get_class_id (cl);
        SLclass_patch_intrin_fun_table1 (Atom_Intrinsics, DUMMY_DBOBJ_MMT_TYPE, DB_MMT_Type_Id);
     }

   return 0;
}

void deinit_atom_module (void)
{
   quit_atom ();
}

/*}}}*/
