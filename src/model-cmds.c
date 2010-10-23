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

/* $Id: model-cmds.c,v 1.41 2004/09/09 11:31:57 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>

#include "isis.h"
#include "util.h"
#include "db-atomic.h"
#include "db-em.h"
#include "db-cie.h"
#include "model.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

/* DUMMY_MODEL_TYPE is a temporary hack that will be modified to the true
  * id once the interpreter provides it when the class is registered.  See below
  * for details.  The reason for this is simple: for a module, the type-id
  * must be assigned dynamically.
  */
#define DUMMY_MODEL_TYPE   255

/* Line profile functions */
typedef struct
{
   Isis_Line_Profile_Type *profile;
}
Line_Profile_Type;
static int Line_Profile_Type_Id = -1;

static SLang_MMT_Type *create_mmt_profile_type (Isis_Line_Profile_Type *x) /*{{{*/
{
   Line_Profile_Type *pt;
   SLang_MMT_Type *mmt;

   if (x == NULL)
     return NULL;

   if (NULL == (pt = (Line_Profile_Type *)SLmalloc (sizeof *pt)))
     return NULL;

   pt->profile = x;

   mmt = SLang_create_mmt (Line_Profile_Type_Id, (VOID_STAR) pt);
   if (NULL == mmt)
     {
        SLfree ((char *)pt);
        return NULL;
     }

   return mmt;
}

/*}}}*/

static void load_line_profile_function_intrin (char *path, char *name) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   isis_fptr_type p;

   if (NULL == (p = isis_load_function (path, name, "line_profile")))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (NULL == (mmt = create_mmt_profile_type ((Isis_Line_Profile_Type *)p)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == SLang_push_mmt (mmt))
     SLang_free_mmt (mmt);
}

/*}}}*/

/* Model definition */

/* This type is for the S-Lang model class */
typedef struct
{
   Model_t *model;
}
Model_Type;
static int Model_Type_Id = -1;

/* This points to the default internal model */
static Model_t *Model_List = NULL;

/*{{{ Stuff for operating on the internal model */

static void quit_model (void) /*{{{*/
{
   Model_end (Model_List);
   Model_List = NULL;
}

/*}}}*/

static int _load_ascii_model (char *filename) /*{{{*/
{
   Model_t *m = Model_load_ascii_file (filename);

   if (NULL == m)
     {
        if (Model_List != NULL)
          isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "kept existing model");
        return -1;
     }

   Model_end (Model_List);
   Model_List = m;

   return 0;
}

/*}}}*/

static int _save_model (Model_t *m, char *fname) /*{{{*/
{
   FILE *fp;
   int ret;

   errno = 0;
   if (NULL == (fp = fopen (fname, "w")))
     {
        isis_vmesg (FAIL, I_WRITE_OPEN_FAILED, __FILE__, __LINE__, "%s", fname);
        return -1;
     }

   ret = Model_print_model (fp, m);

   isis_fclose(fp);

   return ret;
}

/*}}}*/

static int _save_ascii_model (char *fname) /*{{{*/
{
   return _save_model (Model_List, fname);
}

/*}}}*/

static void _list_model (void) /*{{{*/
{
   (void) Model_print_model (stdout, Model_List);
   fflush (stdout);
}

/*}}}*/

static void _edit_model (char *file) /*{{{*/
{
   (void) edit_temp_file (_save_ascii_model, _load_ascii_model,
                          (*file == 0) ? NULL : file);
}

/*}}}*/

static void _set_model_profile (int *idx) /*{{{*/
{
   Model_set_profile_function (*idx);
}

/*}}}*/

static int pop_line_profile_info (Model_Info_Type *info) /*{{{*/
{
   SLang_Array_Type *a;
   SLang_MMT_Type *mmt;
   Line_Profile_Type *lp;

   if (-1 == SLang_pop_array_of_type (&a, SLANG_DOUBLE_TYPE))
     return -1;

   if (NULL == (mmt = SLang_pop_mmt (Line_Profile_Type_Id)))
     {
        SLang_free_array (a);
        return -1;
     }

   lp = (Line_Profile_Type *)SLang_object_from_mmt (mmt);
   if (lp->profile == NULL)
     {
        SLang_free_array (a);
        SLang_free_mmt (mmt);
        return -1;
     }

   info->profile = lp->profile;
   SLang_free_mmt (mmt);

   SLang_free_array (info->profile_params);
   info->profile_params = a;

   return 0;
}

/*}}}*/

static int pop_line_modifier_info (Model_Info_Type *info) /*{{{*/
{
   SLang_Array_Type *a;
   SLang_Name_Type *f;
   int num_extra_args;

   if (-1 == SLang_pop_array_of_type (&a, SLANG_DOUBLE_TYPE))
     return -1;

   if (NULL == (f = SLang_pop_function()))
     {
        SLang_free_array (a);
        return -1;
     }

   if (-1 == SLang_pop_integer (&num_extra_args))
     {
        SLang_free_array (a);
        return -1;
     }

   isis_free_args (info->line_emis_modifier_args);
   info->line_emis_modifier_args = NULL;
   if (num_extra_args)
     {
        info->line_emis_modifier_args = isis_pop_args (num_extra_args);
     }

   SLang_free_array (info->line_emis_modifier_params);
   info->line_emis_modifier_params = a;

   SLang_free_function (info->line_emis_modifier);
   info->line_emis_modifier = f;

   return 0;
}

/*}}}*/

static int pop_ionpop_modifier_info (Model_Info_Type *info) /*{{{*/
{
   SLang_Array_Type *a;
   SLang_Name_Type *f;
   int num_extra_args;

   if (-1 == SLang_pop_array_of_type (&a, SLANG_DOUBLE_TYPE))
     return -1;

   if (NULL == (f = SLang_pop_function()))
     {
        SLang_free_array (a);
        return -1;
     }

   if (-1 == SLang_pop_integer (&num_extra_args))
     {
        SLang_free_array (a);
        return -1;
     }

   isis_free_args (info->ionpop_args);
   info->ionpop_args = NULL;
   if (num_extra_args)
     {
        info->ionpop_args = isis_pop_args (num_extra_args);
     }

   SLang_free_array (info->ionpop_params);
   info->ionpop_params = a;

   SLang_free_function (info->ionpop_modifier);
   info->ionpop_modifier = f;

   return 0;
}

/*}}}*/

static SLang_CStruct_Field_Type Model_Info_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Model_Info_Type, contrib_flag, "contrib_flag", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Model_Info_Type, line_list, "line_list", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int init_model_info (Model_Info_Type *info) /*{{{*/
{
   info->db = ptr_to_atomic_db ();
   info->em = ptr_to_emissivity_db ();
   if ((NULL == info->db) || (NULL == info->em))
     return -1;

   info->contrib_flag = MODEL_LINES_AND_CONTINUUM;
   info->line_list = NULL;

   for (;;)
     {
        int status;
        char *s;

        if ((0 == SLstack_depth())
            || (SLANG_STRING_TYPE != SLang_peek_at_stack ()))
          {
             break;
          }

        if (-1 == SLpop_string (&s))
          return -1;

        if (0 == strcmp (s, "aped_profile"))
          {
             status = pop_line_profile_info (info);
          }
        else if (0 == strcmp (s, "aped_line_modifier"))
          {
             status = pop_line_modifier_info (info);
          }
        else if (0 == strcmp (s, "aped_ionpop_modifier"))
          {
             status = pop_ionpop_modifier_info (info);
          }
        else if (0 == strcmp (s, "aped_hook"))
          {
             status = SLang_pop_cstruct ((VOID_STAR)info, Model_Info_Layout);
             if (info->line_list != NULL)
               {
                  SLang_Array_Type *at = info->line_list;
                  SLtype data_type = at->data_type;

                  if ((data_type != SLANG_NULL_TYPE)
                      && (data_type != SLANG_INT_TYPE)
                      && (at->num_elements > 0))
                    {
                       isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__,
                                   "aped_hook: expected Integer_Type line_list array (data_type=%d), got data_type=%d",
                                   SLANG_INT_TYPE, data_type);
                       status = -1;
                    }
               }
          }
        else
          {
             status = -1;
             SLfree(s);
             break;
          }

        if (status)
          {
             isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__,
                         "processing %s model qualifiers", s ? s : "(null)");
             SLfree(s);
             return status;
          }

        SLfree(s);
     }

   return 0;
}

/*}}}*/

static void free_model_info (Model_Info_Type *info) /*{{{*/
{
   isis_free_args (info->line_emis_modifier_args);
   SLang_free_array (info->line_emis_modifier_params);
   SLang_free_function (info->line_emis_modifier);
   SLang_free_array (info->profile_params);
   SLang_free_cstruct ((VOID_STAR)info, Model_Info_Layout);
   isis_free_args (info->ionpop_args);
   SLang_free_array (info->ionpop_params);
   SLang_free_function (info->ionpop_modifier);
}

/*}}}*/

typedef struct
{
   SLang_Array_Type *lo;
   SLang_Array_Type *hi;
   SLang_Array_Type *val;
}
Model_Result_Type;

static void free_model_grid (Model_Result_Type *m) /*{{{*/
{
   if (m == NULL) return;
   SLang_free_array (m->lo);
   SLang_free_array (m->hi);
}

/*}}}*/

static int get_model_grid (Model_Result_Type *m) /*{{{*/
{
   SLindex_Type n;

   if (-1 == SLang_pop_array_of_type (&m->hi, SLANG_DOUBLE_TYPE)
       || m->hi == NULL
       || -1 == SLang_pop_array_of_type (&m->lo, SLANG_DOUBLE_TYPE)
       || m->lo == NULL
       || m->lo->num_elements != m->hi->num_elements)
     {
        free_model_grid (m);
        return -1;
     }

   n = m->lo->num_elements;

   m->val = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1);
   if (NULL == m->val)
     {
        free_model_grid (m);
        return -1;
     }

   return 0;
}

/*}}}*/

static void do_model_calculation (Model_t *m) /*{{{*/
{
   Model_Result_Type gr;
   Model_Info_Type info;
   double *lo, *hi, *val;
   int n;

   memset ((char *)&gr, 0, sizeof gr);
   memset ((char *)&info, 0, sizeof info);

   if (-1 == get_model_grid (&gr))
     {
        isis_vmesg (INTR, I_WARNING, __FILE__, __LINE__, "wavelength grid not set");
        goto finish;
     }

   if (-1 == init_model_info (&info))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "initializing model");
        free_model_info (&info);
        goto finish;
     }

   n = gr.lo->num_elements;
   lo = (double *)gr.lo->data;
   hi = (double *)gr.hi->data;
   val = (double *)gr.val->data;

   if (-1 == Model_spectrum (m, &info, lo, hi, n, val))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "computing model spectrum");
     }

   finish:

   SLang_push_array (gr.val, 1);
   free_model_grid (&gr);
   free_model_info (&info);
}

/*}}}*/

static int handle_abundance_list (SLang_Array_Type **sl_abun, SLang_Array_Type **sl_elem) /*{{{*/
{
   *sl_abun = *sl_elem = NULL;

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLdo_pop ();
        SLdo_pop ();
        return 1;
     }
   else
     {
        if ((-1 == SLang_pop_array_of_type (sl_elem, SLANG_UINT_TYPE))
            || (-1 == SLang_pop_array_of_type (sl_abun, SLANG_FLOAT_TYPE))
            || (*sl_elem == NULL) || (*sl_abun == NULL)
            || ((*sl_elem)->num_elements != (*sl_abun)->num_elements))
          {
             SLang_free_array (*sl_abun);
             SLang_free_array (*sl_elem);
             *sl_abun = *sl_elem = NULL;
             return -1;
          }
     }

   return 0;
}

/*}}}*/

static SLang_CStruct_Field_Type Model_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Model_t, norm, "norm", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Model_t, metal_abund, "metal_abund", SLANG_FLOAT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Model_t, redshift, "redshift", SLANG_FLOAT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Model_t, vturb, "vturb", SLANG_FLOAT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Model_t, density, "density", SLANG_FLOAT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Model_t, temperature, "temperature", SLANG_FLOAT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int pop_model_component (Model_t *m, SLang_Array_Type **sl_abun, /*{{{*/
                                SLang_Array_Type **sl_elem)
{
   int ret;

   if (m == NULL)
     return -1;

   if (-1 == (ret = handle_abundance_list (sl_abun, sl_elem)))
     return -1;

   if (-1 == SLang_pop_cstruct ((VOID_STAR)m, Model_Layout))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return ret;
}

/*}}}*/

static Model_t *_add_component (Model_t *head) /*{{{*/
{
   SLang_Array_Type *sl_abun, *sl_elem;
   float *elem_abun = NULL;
   unsigned int *elem = NULL;
   unsigned int num_elems = 0;
   Model_t m;
   Model_t *h;
   int ret;

   sl_abun = sl_elem = NULL;

   if ((ret = pop_model_component (&m, &sl_abun, &sl_elem)) < 0)
     {
        SLang_free_array (sl_abun);
        SLang_free_array (sl_elem);
        isis_throw_exception (Isis_Error);
        return NULL;
     }

   if (ret == 0)
     {
        elem_abun = (float *)sl_abun->data;
        elem = (unsigned int *)sl_elem->data;
        num_elems = sl_elem->num_elements;
     }

   h = Model_add_component (head, &m, elem, elem_abun, num_elems);
   if (h == NULL)
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "updating model");
     }
   else if (head == NULL)
     {
        head = h;
     }

   SLang_free_array (sl_abun);
   SLang_free_array (sl_elem);

   return head;
}

/*}}}*/

static void add_component (int *init) /*{{{*/
{
   Model_t *head;

   head = *init ? NULL : Model_List;
   head = _add_component (head);
   if ((head != NULL) && (*init != 0))
     {
        Model_end (Model_List);
        Model_List = head;
     }
}

/*}}}*/

static void calc_model_list (void) /*{{{*/
{
   Model_t *m = Model_List;

   if (NULL == m)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "model not defined");
        return;
     }

   do_model_calculation (m);
}
/*}}}*/

/*}}}*/

/* Stuff for operating on the S-Lang model class */

static SLang_MMT_Type *create_mmt_model_type (Model_t *x) /*{{{*/
{
   Model_Type *mt;
   SLang_MMT_Type *mmt;

   if (x == NULL)
     return NULL;

   if (NULL == (mt = (Model_Type *)SLmalloc (sizeof *mt)))
     return NULL;

   memset ((char *) mt, 0, sizeof *mt);
   mt->model = x;

   mmt = SLang_create_mmt (Model_Type_Id, (VOID_STAR) mt);
   if (NULL == mmt)
     {
        SLfree ((char *)mt);
        Model_end (x);
        return NULL;
     }

   return mmt;
}

/*}}}*/

static void cl_add_component (void) /*{{{*/
{
   SLang_MMT_Type *mmt;
   Model_Type *mt;
   Model_t *x, *m;

   if (Model_Type_Id == SLang_peek_at_stack())
     {
        if ((NULL == (mmt = SLang_pop_mmt (Model_Type_Id)))
            || (NULL == (mt = (Model_Type *) SLang_object_from_mmt (mmt))))
          goto return_error;

        m = mt->model;
     }
   else
     {
        SLdo_pop ();
        mmt = NULL;
        m = NULL;
     }

   if (NULL == (x = _add_component (m)))
     goto return_error;

   if (mmt == NULL)
     {
        if (NULL == (mmt = create_mmt_model_type (x)))
          goto return_error;
        if (-1 == SLang_push_mmt (mmt))
          SLang_free_mmt (mmt);
        return;
     }

   SLang_push_mmt (mmt);
   SLang_free_mmt (mmt);
   return;

   return_error:
   isis_throw_exception (Isis_Error);
   SLang_free_mmt (mmt);               /* NULL ok */
}

/*}}}*/

static void cl_calc_model (Model_Type *mt) /*{{{*/
{
   if (NULL == mt || NULL == mt->model)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "model not defined");
        return;
     }

   do_model_calculation (mt->model);
}

/*}}}*/

static void cl_get_model_details (Model_Type *mt)
{
   Model_t *m;

   if (NULL == mt || NULL == mt->model)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "model not defined");
        return;
     }

   for (m = mt->model ; m != NULL; m = m->next)
     {
        Model_t m_cpy;

        memset ((char *)&m_cpy, 0, sizeof m_cpy);
        m_cpy.norm = m->norm;
        m_cpy.metal_abund = m->metal_abund;
        m_cpy.redshift = m->redshift;
        m_cpy.vturb = m->vturb;
        m_cpy.density = m->density;
        m_cpy.temperature = m->temperature;

        SLang_push_cstruct ((VOID_STAR)&m_cpy, Model_Layout);
        SLang_push_array (m->line_flux, 0);
     }
}

static void cl_load_ascii_model (char *filename) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   Model_t *m;

   m = Model_load_ascii_file (filename);
   if (m == NULL)
     isis_throw_exception (Isis_Error);

   mmt = create_mmt_model_type (m);
   if (mmt == NULL)
     isis_throw_exception (Isis_Error);

   SLang_push_mmt (mmt);
}
/*}}}*/

static int cl_save_ascii_model (Model_Type *mt, char *fname) /*{{{*/
{
   if (mt->model == NULL)
     return -1;
   return _save_model (mt->model, fname);
}

/*}}}*/

static void cl_list_model (Model_Type *mt) /*{{{*/
{
   Model_print_model (stdout, mt->model);
   fflush (stdout);
}

/*}}}*/

/*{{{ Intrinsics */

#define MT DUMMY_MODEL_TYPE
#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define F SLANG_FLOAT_TYPE
#define D SLANG_DOUBLE_TYPE

static SLang_IConstant_Type Model_Intrin_Const [] = /*{{{*/
{
   MAKE_ICONSTANT("MODEL_LINES_AND_CONTINUUM", MODEL_LINES_AND_CONTINUUM),
   MAKE_ICONSTANT("MODEL_LINES", MODEL_LINES),
   MAKE_ICONSTANT("MODEL_CONTIN", MODEL_CONTIN),
   MAKE_ICONSTANT("MODEL_CONTIN_PSEUDO", MODEL_CONTIN_PSEUDO),
   MAKE_ICONSTANT("MODEL_CONTIN_TRUE", MODEL_CONTIN_TRUE),
   MAKE_ICONSTANT("_isis_max_proton_number", ISIS_MAX_PROTON_NUMBER),
   SLANG_END_ICONST_TABLE
};
/*}}}*/

static SLang_Intrin_Fun_Type Model_Intrinsics [] = /*{{{*/
{
   MAKE_INTRINSIC_1("_load_ascii_model", _load_ascii_model, I, S),
   MAKE_INTRINSIC_1("_save_ascii_model", _save_ascii_model, I, S),
   MAKE_INTRINSIC_1("_edit_model", _edit_model, V, S),
   MAKE_INTRINSIC_I("_set_model_profile", _set_model_profile, V),
   MAKE_INTRINSIC("_list_model", _list_model, V, 0),
   MAKE_INTRINSIC("_calc_model_list", calc_model_list, V, 0),
   MAKE_INTRINSIC_1("_add_component", add_component, V, I),
   MAKE_INTRINSIC("_cl_add_component", cl_add_component, V, 0),
   MAKE_INTRINSIC_1("_cl_calc_model", cl_calc_model, V, MT),
   MAKE_INTRINSIC_1("_cl_get_model_details", cl_get_model_details, V, MT),
   MAKE_INTRINSIC_1("_cl_load_ascii_model", cl_load_ascii_model, V, S),
   MAKE_INTRINSIC_2("_cl_save_ascii_model", cl_save_ascii_model, I, MT, S),
   MAKE_INTRINSIC_1("_cl_list_model", cl_list_model, V, MT),
   MAKE_INTRINSIC_2("load_line_profile_function_intrin", load_line_profile_function_intrin, V, S, S),
   SLANG_END_INTRIN_FUN_TABLE
};

/*}}}*/

#undef V
#undef S
#undef I
#undef F
#undef D
#undef MT

/*}}}*/

static void destroy_model_type (SLtype type, VOID_STAR f) /*{{{*/
{
   Model_Type *mt = (Model_Type *) f;
   (void) type;

   Model_end (mt->model);
   SLfree ((char *) mt);
}

/*}}}*/

static void destroy_line_profile_type (SLtype type, VOID_STAR f) /*{{{*/
{
   Line_Profile_Type *pt = (Line_Profile_Type *) f;
   (void) type;
   SLfree ((char *) pt);
}

/*}}}*/

SLANG_MODULE(model);
int init_model_module_ns (char *ns_name) /*{{{*/
{
   SLang_Class_Type *cl;
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   ns = SLns_create_namespace (ns_name);
   if (NULL == ns)
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if (Model_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("Model_Type")))
          return isis_trace_return(-1);
        (void) SLclass_set_destroy_function (cl, destroy_model_type);

        /* By registering as SLANG_VOID_TYPE,
         * slang will dynamically allocate a type
         */
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (Model_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return isis_trace_return(-1);

        Model_Type_Id = SLclass_get_class_id (cl);
        SLclass_patch_intrin_fun_table1 (Model_Intrinsics, DUMMY_MODEL_TYPE, Model_Type_Id);
     }

   if (Line_Profile_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("Line_Profile_Type")))
          return isis_trace_return(-1);
        (void) SLclass_set_destroy_function (cl, destroy_line_profile_type);

        /* By registering as SLANG_VOID_TYPE,
         * slang will dynamically allocate a type
         */
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (Line_Profile_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return isis_trace_return(-1);

        Line_Profile_Type_Id = SLclass_get_class_id (cl);
     }

   if (-1 == SLns_add_iconstant_table (pub_ns, Model_Intrin_Const, NULL))
     return isis_trace_return(-1);

   if (-1 == SLns_add_intrin_fun_table (ns, Model_Intrinsics, NULL))
     return isis_trace_return(-1);

   return 0;
}

/*}}}*/

void deinit_model_module (void) /*{{{*/
{
   quit_model ();
}

/*}}}*/
