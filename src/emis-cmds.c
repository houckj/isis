/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2016  Massachusetts Institute of Technology

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

/* $Id: emis-cmds.c,v 1.31 2004/09/09 11:31:56 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <slang.h>

#include "isis.h"
#include "db-atomic.h"
#include "db-em.h"
#include "db-cie.h"
#include "util.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

typedef struct
{
   EM_t *ptr;
}
EM_MMT_Type;
static int EM_MMT_Type_Id = -1;

static void destroy_emis_db_mmt_type (SLtype type, VOID_STAR f) /*{{{*/
{
   EM_MMT_Type *mt = (EM_MMT_Type *)f;
   (void) type;

   if (mt != NULL)
     {
        EM_end (mt->ptr);
     }

   SLfree ((char *)f);
}

/*}}}*/

static SLang_MMT_Type *create_emis_db_mmt_type (EM_t *ptr) /*{{{*/
{
   SLang_MMT_Type *mmt;
   EM_MMT_Type *mt;

   if (ptr == NULL)
     return NULL;

   if (NULL == (mt = (EM_MMT_Type *)SLmalloc (sizeof *mt)))
     return NULL;

   mt->ptr = ptr;

   mmt = SLang_create_mmt (EM_MMT_Type_Id, (VOID_STAR) mt);
   if (NULL == mmt)
     {
        SLfree ((char *)mt);
        return NULL;
     }

   return mmt;
}

/*}}}*/

static void push_emis_db_pointer_intrin (EM_t *ptr) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;

   if (NULL == (mmt = create_emis_db_mmt_type (ptr)))
     {
        isis_throw_exception (Isis_Error);
        return;
     }

   if (-1 == SLang_push_mmt (mmt))
     SLang_free_mmt (mmt);
}

/*}}}*/

static EM_t *pop_emis_db_pointer_intrin (void) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   EM_MMT_Type *mt;

   SLang_run_hooks ("_isis->get_emis_db_pointer", 0, NULL);

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLang_pop_null();
        return NULL;
     }

   if (NULL == (mmt = SLang_pop_mmt (EM_MMT_Type_Id)))
     return NULL;

   mt = (EM_MMT_Type *) SLang_object_from_mmt (mmt);
   SLang_free_mmt (mmt);

   return mt->ptr;
}

/*}}}*/

/* start / stop / list */

static void quit_emis (void) /*{{{*/
{
   EM_Use_Memory = EM_USE_MEMORY_DEFAULT;
   EM_Maybe_Missing_Lines = 1;
}

/*}}}*/

#define V SLANG_VOID_TYPE
#define S SLANG_STRING_TYPE

static SLang_CStruct_Field_Type EM_File_Table [] =
{
   MAKE_CSTRUCT_FIELD(EM_File_Type, contin_emis, "contin_emis", S, 0),
   MAKE_CSTRUCT_FIELD(EM_File_Type, line_emis, "line_emis", S, 0),
   MAKE_CSTRUCT_FIELD(EM_File_Type, ionization, "ionization", S, 0),
   MAKE_CSTRUCT_FIELD(EM_File_Type, abundance, "abundance", S, 0),
   MAKE_CSTRUCT_FIELD(EM_File_Type, filemap, "filemap", S, 0),
   SLANG_END_CSTRUCT_TABLE
};

#undef V
#undef S

static void em_start (float *tmin, float *tmax, float *dmin, float *dmax) /*{{{*/
{
   EM_File_Type f = NULL_EM_FILE_TYPE;
   DB_t *db = _ptr_to_atomic_db();
   EM_t *em = NULL;
   EM_Range_Type r;

   if (1 == SLang_run_hooks ("_isis->emisdb_start_hook", 0))
     {
        if (-1 == SLang_pop_cstruct ((VOID_STAR)&f, EM_File_Table))
          goto finish;
     }

   if (NULL == db)
     goto finish;

   r.trange[0] = *tmin;   r.trange[1] = *tmax;
   r.drange[0] = *dmin;   r.drange[1] = *dmax;

   if (NULL == (em = EM_start (&f, (void *) &r, db)))
     goto finish;

   push_emis_db_pointer_intrin (em);

   finish:
   SLang_free_cstruct ((VOID_STAR)&f, EM_File_Table);

   if (!em)
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing emissivity database");
}

/*}}}*/

EM_t *_ptr_to_emissivity_db (void) /*{{{*/
{
   return pop_emis_db_pointer_intrin ();
}

/*}}}*/

EM_t *ptr_to_emissivity_db (void) /*{{{*/
{
   EM_t *ptr = pop_emis_db_pointer_intrin ();
   if (NULL == ptr)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "Emissivity database not initialized");
        isis_throw_exception (Isis_Error);
     }
   return ptr;
}

/*}}}*/

static void _get_filemap (char *emis_file) /*{{{*/
{
   EM_t *em = _ptr_to_emissivity_db ();
   SLang_Array_Type *sl_t = NULL;
   SLang_Array_Type *sl_d = NULL;
   double *temp = NULL;
   double *dens = NULL;
   unsigned int num = 0;
   SLindex_Type n;

   if (-1 == EM_get_filemap (em, emis_file, NULL, &num, &temp, &dens))
     goto finish;

   n = num;
   sl_t = SLang_create_array (SLANG_DOUBLE_TYPE, 1, temp, &n, 1);
   sl_d = SLang_create_array (SLANG_DOUBLE_TYPE, 1, dens, &n, 1);
   if (sl_t == NULL || sl_d == NULL)
     {
        ISIS_FREE (temp);
        ISIS_FREE (dens);
        SLang_free_array (sl_t);
        SLang_free_array (sl_d);
        isis_throw_exception (Isis_Error);
        goto finish;
     }

   finish:

   SLang_push_array (sl_t, 1);
   SLang_push_array (sl_d, 1);
}

/*}}}*/

/* get ionization balance at given T, dens */

static int get_ioniz_balance (float * charge, float * frac, /*{{{*/
                              int Z, float etemp, float edens,
                              int ion_table_id)
{
   EM_t *em = ptr_to_emissivity_db ();
   int q;

   if (NULL == em)
     return -1;

   for (q=0; q <= Z; q++)
     {
        charge[q] = (float) q;
        if (-1 == EM_get_ion_fraction (&frac[q], etemp, edens, Z, q, ion_table_id, em))
          return -1;
     }

   return 0;
}

/*}}}*/

static void retrieve_ionization_bal (int *Z, float *temp, /*{{{*/
                                     float *dens, int *ion_table_id)
{
   SLang_Array_Type *slcharge, *slfrac;
   SLindex_Type dims[1];

   slcharge = slfrac = NULL;

   dims[0] = *Z + 1;

   slfrac   = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, dims, 1);
   slcharge = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, dims, 1);

   if (NULL == slfrac
       || NULL == slcharge
       || -1 == get_ioniz_balance ((float *) slcharge->data,
                                   (float *) slfrac->data,
                                   *Z, *temp, *dens, *ion_table_id))
     goto done;

   done:
   (void) SLang_push_array (slcharge, 1);
   (void) SLang_push_array (slfrac, 1);
}

/*}}}*/

/* get ionization fraction vs. T  */

static int get_ioniz_vs_temp (float *frac, float *temp, int ntemp, /*{{{*/
                             int Z, int charge, int ion_table_id)
{
   EM_t *em = ptr_to_emissivity_db ();
   int i;

   if (NULL == em)
     return -1;

   for (i=0; i < ntemp; i++)
     {
        /* FIXME:  density not used */
        if (-1 == EM_get_ion_fraction (&frac[i], temp[i], -1.0, Z, charge, ion_table_id, em))
          return -1;
     }

   return 0;
}

/*}}}*/

static void retrieve_ionization_fcn (void) /*{{{*/
{
   SLang_Array_Type *slfrac, *sltemp;
   int Z, ion_charge, dims[1], ion_table_id;

   sltemp = slfrac = NULL;

   if (-1 == SLang_pop_integer (&ion_table_id)
       || -1 == SLang_pop_array_of_type (&sltemp, SLANG_FLOAT_TYPE)
       || sltemp == NULL
       || -1 == SLang_pop_integer (&ion_charge)
       || -1 == SLang_pop_integer (&Z))
     goto free_and_return;

   dims[0] = sltemp->num_elements;

   slfrac = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, dims, 1);

   if (NULL == slfrac
       || -1 == get_ioniz_vs_temp ((float *) slfrac->data,
                                   (float *) sltemp->data, dims[0],
                                   Z, ion_charge, ion_table_id))
     goto free_and_return;

   free_and_return:
   SLang_free_array (sltemp);
   (void) SLang_push_array (slfrac, 1);
}

/*}}}*/

/* get line emissivity vs. T, dens */

static int get_emissivity_fcn (float **emis, int *nemis, /*{{{*/
                               float *temp, int ntemp, float *dens, int ndens,
                               int *line_index, int nlines)
{
   DB_t *db = ptr_to_atomic_db ();
   EM_t *em = ptr_to_emissivity_db ();
   int t, d;

   if (em == NULL || db == NULL)
     return -1;

   *nemis = ntemp * ndens;

   if (NULL == (*emis = (float *) ISIS_MALLOC (*nemis * sizeof(float))))
     return -1;

   memset ((char *) *emis, 0, *nemis * sizeof(float));

   for (d=0; d < ndens; d++)
     {
        for (t=0; t < ntemp; t++)
          {
             float emis_sum;

             if (-1 == EM_sum_line_emissivity (&emis_sum, temp[t], dens[d],
                                               line_index, nlines, em))
               return -1;

             (*emis)[ t + d*ntemp ] = emis_sum;
          }
     }

   return 0;
}

/*}}}*/

static void retrieve_emissivity_fcn (void) /*{{{*/
{
   float *emis = NULL;
   SLindex_Type dims[2];
   SLang_Array_Type *sltemp, *sldens, *slemis, *slindx;
   int nemis, problem = -1;

   sltemp = sldens = slemis = slindx = NULL;

   if (-1 == SLang_pop_array_of_type (&sldens, SLANG_FLOAT_TYPE)
       || sldens == NULL
       || -1 == SLang_pop_array_of_type (&sltemp, SLANG_FLOAT_TYPE)
       || sltemp == NULL
       || -1 == SLang_pop_array_of_type (&slindx, SLANG_INT_TYPE)
       || slindx == NULL)
     goto free_and_return;

   if (-1 == get_emissivity_fcn (&emis, &nemis,
                                 (float *) sltemp->data, sltemp->num_elements,
                                 (float *) sldens->data, sldens->num_elements,
                                 (int *) slindx->data, (int) slindx->num_elements))
     goto free_and_return;

   if (1 == sltemp->num_elements
       || 1 == sldens->num_elements)
     {
        dims[0] = nemis;

        if (NULL == (slemis = SLang_create_array (SLANG_FLOAT_TYPE, 0,
                                                  NULL, dims, 1)))
          goto free_and_return;
     }
   else if (sltemp->num_elements > 1
            && sldens->num_elements > 1)
     {
        dims[0] = sldens->num_elements;
        dims[1] = sltemp->num_elements;

        if (NULL == (slemis = SLang_create_array (SLANG_FLOAT_TYPE, 0,
                                                  NULL, dims, 2)))
          goto free_and_return;
     }
   else
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto free_and_return;
     }

   memcpy ((char *)slemis->data, (char *)emis, nemis * sizeof(float));
   problem = 0;

   free_and_return:
   if (problem)
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "retrieving emissivity function");

   ISIS_FREE (emis);
   SLang_free_array (sltemp);
   SLang_free_array (sldens);
   SLang_free_array (slindx);
   (void) SLang_push_array (slemis, 1);
}

/*}}}*/

static void _line_emissivities_on_datafile_grid (int *idx) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();
   float *emis, *temps, *densities;
   SLang_Array_Type *sl_emis, *sl_temps, *sl_densities;
   int num_points;
   int line_index = (int) *idx;

   emis = temps = densities = NULL;
   sl_emis = sl_temps = sl_densities = NULL;

   if (NULL == em)
     goto push;

   if (0 == EM_get_line_emissivity_function (&emis, &temps, &densities, &num_points,
                                              line_index, em))
     {
        sl_emis = SLang_create_array (SLANG_FLOAT_TYPE, 1, emis, &num_points, 1);
        sl_temps = SLang_create_array (SLANG_FLOAT_TYPE, 1, temps, &num_points, 1);
        sl_densities = SLang_create_array (SLANG_FLOAT_TYPE, 1, densities, &num_points, 1);
     }

   push:

   /* ok to push NULLs */
   (void) SLang_push_array (sl_emis, 1);
   (void) SLang_push_array (sl_temps, 1);
   (void) SLang_push_array (sl_densities, 1);
}

/*}}}*/

/* get continuum emissivity for given T, dens, Z, q */

static void _get_continuum (void) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();
   EM_cont_type_t *p = NULL;
   EM_cont_select_t s = {0, 0, NULL};
   SLang_Array_Type *sl_lo, *sl_hi, *sl_true, *sl_pseudo;
   float rel_abun[ISIS_MAX_PROTON_NUMBER+1];
   float temp, dens;
   int nbins, q, Z, size;

   if (NULL == em)
     return;

   for (Z=0; Z <= ISIS_MAX_PROTON_NUMBER; Z++)
     rel_abun[Z] = 1.0;

   sl_lo = sl_hi = sl_true = sl_pseudo = NULL;

   if (-1 == SLang_pop_integer (&q)
       || -1 == SLang_pop_integer (&Z)
       || -1 == SLang_pop_float (&dens)
       || -1 == SLang_pop_float (&temp)
       || -1 == SLang_pop_array_of_type (&sl_hi, SLANG_DOUBLE_TYPE)
       || sl_hi == NULL
       || -1 == SLang_pop_array_of_type (&sl_lo, SLANG_DOUBLE_TYPE)
       || sl_lo == NULL)
     goto finish;

   if (sl_hi->num_elements != sl_lo->num_elements)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto finish;
     }

   nbins = (int) sl_hi->num_elements;
   size = nbins * sizeof(double);

   if (NULL == (p = EM_new_continuum (nbins)))
     goto finish;

   memcpy ((char *)p->wlhi, (char *) sl_hi->data, size);
   memcpy ((char *)p->wllo, (char *) sl_lo->data, size);

   s.Z = Z;
   s.q = q;
   s.rel_abun = rel_abun;
   if (-1 == EM_get_continuum (p, temp, dens, NULL /* FIXME - ionpop_new? */, &s, em))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "retrieving continuum data");
        goto finish;
     }

   if ((NULL == (sl_true = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nbins, 1)))
       || (NULL == (sl_pseudo = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nbins, 1)))
       )
     goto finish;

   memcpy ((char *)sl_true->data, (char *) p->true_contin, size);
   memcpy ((char *)sl_pseudo->data, (char *) p->pseudo, size);

   finish:

   EM_free_continuum (p);
   SLang_free_array (sl_hi);
   SLang_free_array (sl_lo);
   (void) SLang_push_array (sl_true, 1);
   (void) SLang_push_array (sl_pseudo, 1);
}

/*}}}*/

/* alt abundance and ionization tables */

typedef struct
{
   SLang_Array_Type *abun;
   SLang_Array_Type *z;
   char *name;
}
Abund_Table_Type;

static SLang_CStruct_Field_Type Abund_Table_Type_Layout [] =
{
   MAKE_CSTRUCT_FIELD(Abund_Table_Type, abun, "abun", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Abund_Table_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Abund_Table_Type, z, "z", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static void _get_abundance_table (int *k) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();
   Abund_Table_Type t;
   char *name = NULL;
   float *a = NULL;
   int *z = NULL;
   int n, i;

   memset ((char *)&t, 0, sizeof t);

   i = (*k < 0) ? EM_get_standard_abundance (em) : *k;

   if (0 == EM_get_abundance_table (em, i, &name, &a, &z, &n))
     {
        t.abun = SLang_create_array (SLANG_FLOAT_TYPE, 0, a, &n, 1);
        t.z = SLang_create_array (SLANG_INT_TYPE, 0, z, &n, 1);
        t.name = name;
     }

   SLang_push_cstruct ((VOID_STAR)&t, Abund_Table_Type_Layout);
   ISIS_FREE(name);
}

/*}}}*/

static void _choose_abs_abundance_table (int *k) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();

   if (-1 == EM_set_chosen_abundance (em, *k))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "abundance table not changed");
}

/*}}}*/

static int _index_for_abundance_table (char *name) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();
   int k;

   if (-1 == (k = EM_get_abundance_table_index_by_name (em, name)))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "abundance table '%s' not found", name);

   return k;
}

/*}}}*/

static void _list_abund_tables (int *verbose) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();

   if (NULL == em)
     return;

   (void) EM_list_abundance_tables (stdout, em, *verbose);
}

/*}}}*/

static int _add_abund_table (void) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();
   Abund_Table_Type t;
   int num_abun, status;
   float *abun;
   int *z;

   if (em == NULL)
     return -1;

   memset ((char *)&t, 0, sizeof (Abund_Table_Type));

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&t, Abund_Table_Type_Layout))
     return -1;

   if (t.z->num_elements != t.abun->num_elements)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "mismatched array sizes in abundance table");
        return -1;
     }

   if ((-1 == isis_coerce_array_to_type (&t.z, SLANG_INT_TYPE))
       ||(-1 == isis_coerce_array_to_type (&t.abun, SLANG_FLOAT_TYPE)))
     {
        SLang_free_cstruct ((VOID_STAR)&t, Abund_Table_Type_Layout);
        return -1;
     }

   z = (int *)t.z->data;
   abun = (float *)t.abun->data;
   num_abun = t.z->num_elements;

   status = EM_add_abundance_table (em, t.name, abun, z, num_abun);
   SLang_free_cstruct ((VOID_STAR)&t, Abund_Table_Type_Layout);

   return status;
}

/*}}}*/

static void _free_alt_ionization_table (void) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();

   if (NULL == em)
     return;

   EM_free_alt_ionization_table (em);
}

/*}}}*/

static void _load_alt_ionization_table (char *file) /*{{{*/
{
   EM_t *em = ptr_to_emissivity_db ();

   if (NULL == em)
     return;

   if (-1 == EM_load_alt_ionization_table (file, em))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "loading ionization table %s", file);
}

/*}}}*/

/* SLang intrinsics */

#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define F SLANG_FLOAT_TYPE

static SLang_Intrin_Var_Type Emis_Intrin_Vars [] = /*{{{*/
{
   MAKE_VARIABLE("Use_Memory", &EM_Use_Memory, SLANG_UINT_TYPE, 0),
   MAKE_VARIABLE("Incomplete_Line_List", &EM_Maybe_Missing_Lines, SLANG_INT_TYPE, 0),
   MAKE_VARIABLE("EM_Hash_Table_Size_Hint", &EM_Hash_Table_Size_Hint, SLANG_UINT_TYPE, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

/*}}}*/

#define DUMMY_EMOBJ_MMT_TYPE 255
#define MT DUMMY_EMOBJ_MMT_TYPE

static SLang_Intrin_Fun_Type Emis_Intrinsics [] = /*{{{*/
{
   MAKE_INTRINSIC_4("_em_start", em_start, V, F, F, F, F),
   MAKE_INTRINSIC_S("_get_filemap", _get_filemap, V),
   MAKE_INTRINSIC_4("_ioniz_bal", retrieve_ionization_bal, V, I, F, F, I),
   MAKE_INTRINSIC("_ioniz_fcn", retrieve_ionization_fcn, V, 0),
   MAKE_INTRINSIC("_line_em", retrieve_emissivity_fcn, V, 0),
   MAKE_INTRINSIC_I("_line_emissivities_on_datafile_grid", _line_emissivities_on_datafile_grid, V),
   MAKE_INTRINSIC("_get_continuum", _get_continuum, V, 0),
   MAKE_INTRINSIC_I("_list_abund_tables", _list_abund_tables, V),
   MAKE_INTRINSIC_I("_choose_abs_abundance_table", _choose_abs_abundance_table, V),
   MAKE_INTRINSIC_S("_index_for_abundance_table", _index_for_abundance_table, I),
   MAKE_INTRINSIC_I("_get_abundance_table", _get_abundance_table, V),
   MAKE_INTRINSIC("_add_abund_table", _add_abund_table, I, 0),
   MAKE_INTRINSIC_S("_load_alt_ionization_table", _load_alt_ionization_table, V),
   MAKE_INTRINSIC("_free_alt_ionization_table", _free_alt_ionization_table, V, 0),
   SLANG_END_INTRIN_FUN_TABLE
};

/*}}}*/

#undef V
#undef I
#undef F

SLANG_MODULE(emis);
int init_emis_module_ns (char *ns_name) /*{{{*/
{
   SLang_Class_Type *cl;
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if ((-1 == SLns_add_intrin_fun_table (ns, Emis_Intrinsics, NULL))
       || (-1 == SLns_add_intrin_var_table (pub_ns, Emis_Intrin_Vars, NULL)))
     return isis_trace_return(-1);

   if (EM_MMT_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("EM_Object_Type")))
          return isis_trace_return(-1);
        (void) SLclass_set_destroy_function (cl, destroy_emis_db_mmt_type);
        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (EM_MMT_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return isis_trace_return(-1);
        EM_MMT_Type_Id = SLclass_get_class_id (cl);
        SLclass_patch_intrin_fun_table1 (Emis_Intrinsics, DUMMY_EMOBJ_MMT_TYPE, EM_MMT_Type_Id);
     }

   return 0;
}

/*}}}*/

#ifdef __cplusplus
extern "C"
#endif
void deinit_emis_module (void);
void deinit_emis_module (void) /*{{{*/
{
   quit_emis ();
}

/*}}}*/

