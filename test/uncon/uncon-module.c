/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2009 Massachusetts Institute of Technology

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

/* Purpose: Provide a S-Lang interface to the 'uncon'
 *          suite of optimizer tests from netlib.org
 * Author:  John C. Houck <houck@space.mit.edu>
 *   Date:  12/2004
 */

#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <math.h>
#include <float.h>

#include "slang.h"

#if (SLANG_VERSION<20000)
void SLang_set_error (int err) /*{{{*/
{
   SLang_Error = err;
}

/*}}}*/
#endif

static int Uncon_Type_Id = -1;

typedef void Fptr_Type (double *,int *,double *,int *,double *,
                        double [][1],int *, double *, int *);

static Fptr_Type *Fortran_Sub_Ptr;

typedef struct
{
   Fptr_Type *fun;
}
Uncon_Type;

typedef struct
{
   SLang_Array_Type *x, *f;
   double ftf;
   int n, m, mode;
}
Uncon_Interface_Type;

static int unsupported_mode (int mode) /*{{{*/
{
   return ((mode > 0) && (mode % 100));
}

/*}}}*/

static void fortran_sub (Uncon_Interface_Type *u) /*{{{*/
{
   double *x, *f;
   int mode;

   if ((Fortran_Sub_Ptr == NULL) || (u == NULL)
       || unsupported_mode (u->mode))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   x = ((double *)u->x->data);
   f = ((double *)u->f->data);
   mode = u->mode;

     {
        double fj[1][1], g[1];
        int lfj;
        /* unused parameters (for unsupported modes) */
        lfj = 1; g[0] = 0.0; fj[0][0] = 0.0;

        (*Fortran_Sub_Ptr)(x, &u->n, f, &u->m, &u->ftf, fj, &lfj, g, &mode);
     }

   if (mode == 999)
     {
        int i, m = u->m;
        double huge_value = 0.5 * sqrt(DBL_MAX/m);
        for (i = 0; i < m; i++)
          {
             f[i] = huge_value;
          }
     }
}

/*}}}*/

static SLang_CStruct_Field_Type Uncon_Interface_Layout [] = /*{{{*/
{
   MAKE_CSTRUCT_FIELD (Uncon_Interface_Type, x, "x", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Uncon_Interface_Type, n, "n", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Uncon_Interface_Type, f, "f", SLANG_ARRAY_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Uncon_Interface_Type, m, "m", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Uncon_Interface_Type, ftf, "ftf", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Uncon_Interface_Type, mode, "mode", SLANG_INT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

/*}}}*/

static void _uncon_hook (Uncon_Type *ut) /*{{{*/
{
   Uncon_Interface_Type u;

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&u, Uncon_Interface_Layout))
     {
        SLang_set_error (SL_INTRINSIC_ERROR);
        return;
     }

   Fortran_Sub_Ptr = ut->fun;
   fortran_sub (&u);

   if (-1 == SLang_push_cstruct ((VOID_STAR)&u, Uncon_Interface_Layout))
        SLang_set_error (SL_INTRINSIC_ERROR);
   
   (void) SLang_free_cstruct ((VOID_STAR)&u, Uncon_Interface_Layout);
}

/*}}}*/

static void handle_link_error (char *path, char *name) /*{{{*/
{
   const char *error = dlerror ();

   if (error != NULL)
     fprintf (stderr, "Link error:  %s\n", error);
   else
     fprintf (stderr, "Link error:  failed loading %s from %s\n",
              name, path);
}
/*}}}*/

static Fptr_Type *load_function (char *path, char *name) /*{{{*/
{
   Fptr_Type *f = NULL;
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

   if (NULL == (handle = dlopen (path, DLOPEN_FLAG)))
     {
        handle_link_error (path, name);
        return NULL;
     }

   if (NULL == (f = (Fptr_Type *) dlsym (handle, name)))
     {
        handle_link_error (path, name);
        dlclose (handle);
        return NULL;
     }

   return f;
}
/*}}}*/

static int _load_uncon_fun (SLang_Ref_Type *ref, char *file, char *fun_name) /*{{{*/
{
   SLang_MMT_Type *mmt = NULL;
   Uncon_Type *ut = NULL;
   Fptr_Type *fptr;

   if (-1 == SLang_assign_to_ref (ref, SLANG_NULL_TYPE, NULL))
     return -1;

   if (NULL == (fptr = load_function (file, fun_name)))
     return -1;

   if (NULL == (ut = malloc(sizeof(*ut))))
     return -1;

   ut->fun = (Fptr_Type *)fptr;

   if (NULL == (mmt = SLang_create_mmt (Uncon_Type_Id, (void *) ut)))
     {
        free (ut);
        return -1;
     }

   if (-1 == SLang_assign_to_ref (ref, Uncon_Type_Id, &mmt))
     {
        free (ut);
        SLang_free_mmt (mmt);
        return -1;
     }

   return 0;
}

/*}}}*/

/*{{{ Intrinsics */

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
   MAKE_INTRINSIC_3("_load_uncon_fun", _load_uncon_fun, I, R, S, S),
   MAKE_INTRINSIC_1("_uncon_hook", _uncon_hook, V, MT),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef I
#undef R
#undef S
#undef V
#undef MT

static SLang_IConstant_Type Intrin_Const [] =
{
   MAKE_ICONSTANT("UNCON_SET_BESTF",   -2),
   MAKE_ICONSTANT("UNCON_SET_SIZE",    -1),
   MAKE_ICONSTANT("UNCON_SET_INITPT",   0),
   MAKE_ICONSTANT("UNCON_COMPUTE",   1100),
   SLANG_END_ICONST_TABLE
};

/*}}}*/

static void patchup_intrinsic_table (void) /*{{{*/
{
   SLang_Intrin_Fun_Type *f;

   f = Intrinsics;
   while (f->name != NULL)
     {
        unsigned int i, nargs;
        SLtype *args;

        nargs = f->num_args;
        args = f->arg_types;
        for (i = 0; i < nargs; i++)
          {
             if (args[i] == DUMMY_MODEL_TYPE)
               args[i] = Uncon_Type_Id;
          }

        /* For completeness */
        if (f->return_type == DUMMY_MODEL_TYPE)
          f->return_type = Uncon_Type_Id;

        f++;
     }
}

/*}}}*/

static void free_uncon_fun_type (SLtype type, void *f) /*{{{*/
{
   Uncon_Type *ut = (Uncon_Type *)f;
   (void) type; (void) ut;
   free (ut);
}

/*}}}*/

SLANG_MODULE(uncon);
int init_uncon_module_ns (char *ns_name) /*{{{*/
{
   SLang_Class_Type *cl;
   SLang_NameSpace_Type *ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return -1;

   if (Uncon_Type_Id == -1)
     {
        if (NULL == (cl = SLclass_allocate_class ("Uncon_Type")))
          return -1;

        (void) SLclass_set_destroy_function (cl, free_uncon_fun_type);

        if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (Uncon_Type),
                                          SLANG_CLASS_TYPE_MMT))
          return -1;

        Uncon_Type_Id = SLclass_get_class_id (cl);
        patchup_intrinsic_table ();
     }

   if ((-1 == SLns_add_intrin_fun_table (ns, Intrinsics, NULL))
       || (-1 == SLns_add_iconstant_table (ns, Intrin_Const, NULL)))
     return -1;

   return 0;
}

/*}}}*/
