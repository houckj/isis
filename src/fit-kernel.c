/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2019 Massachusetts Institute of Technology

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

/* $Id: fit-kernel.c,v 1.27 2004/02/09 11:14:20 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>
#include "util.h"
#include "isis.h"
#include "fit.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

/* Manage list of fit kernels */

Isis_Kernel_Def_t * Fit_new_kernel (void) /*{{{*/
{
   Isis_Kernel_Def_t *def;

   if (NULL == (def = (Isis_Kernel_Def_t *) ISIS_MALLOC (sizeof(Isis_Kernel_Def_t))))
     return NULL;

   memset ((char *)def, 0, sizeof(*def));

   def->fun_type = UINT_MAX;
   def->kernel_id = UINT_MAX;
   def->malloced_kernel = 1;
   def->next = NULL;

   return def;
}

/*}}}*/

void Fit_free_kernel (Isis_Kernel_Def_t *def) /*{{{*/
{
   if ((NULL == def) || (0 == def->malloced_kernel))
     return;

   ISIS_FREE (def);
}

/*}}}*/

void Fit_free_aux_kernels (Isis_Kernel_Def_t *t) /*{{{*/
{
   if (t == NULL) return;

   while (t != NULL)
     {
        Isis_Kernel_Def_t *n = t->next;
        Fit_free_kernel (t);
        t = n;
     }
}

/*}}}*/

void Fit_push_kernel_names (Isis_Kernel_Def_t *t) /*{{{*/
{
   Isis_Kernel_Def_t *def = NULL;

   for (def = t; def != NULL; def = def->next)
     {
        SLang_push_string (def->kernel_name);
     }
}

/*}}}*/

Isis_Kernel_Def_t * Fit_find_kernel (Isis_Kernel_Def_t *t, unsigned int kernel_id) /*{{{*/
{
   Isis_Kernel_Def_t *def = NULL;

   for (def = t; def != NULL; def = def->next)
     {
        if (def->kernel_id == kernel_id)
          break;
     }

   if (NULL == def)
     isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Kernel %d not found", kernel_id);

   return def;
}

/*}}}*/

Isis_Kernel_Def_t * Fit_find_kernel_by_name (Isis_Kernel_Def_t *t, char *kernel_name) /*{{{*/
{
   Isis_Kernel_Def_t *def = NULL;

   if (NULL == kernel_name)
     return NULL;

   for (def = t; def != NULL; def = def->next)
     {
        if (!strcmp(def->kernel_name, kernel_name))
          break;
     }

   return def;
}

/*}}}*/

int Fit_append_kernel (Isis_Kernel_Def_t *head, Isis_Kernel_Def_t *def) /*{{{*/
{
   Isis_Kernel_Def_t *p;
   Isis_Kernel_Def_t *last = head;

   if (NULL == def)
     return -1;

   /* check for duplicates (by name) */

   if (!strcmp (head->kernel_name, def->kernel_name))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "can't replace the standard kernel");
        return -1;
     }

   if (-1 == Fit_add_kernel_function (def))
     return -1;

   for (p = head; p->next != NULL; p = p->next)
     {
        if (!strcmp (def->kernel_name, p->next->kernel_name))
          {
             /* replace if duplicate name, keeping old kernel id */
             def->kernel_id = p->next->kernel_id;
             def->next = p->next->next;
             Fit_free_kernel (p->next);
             p->next = def;
             return 0;
          }
        last = p->next;
     }

   /* not duplicate name, so append */
   last->next = def;
   def->next = NULL;
   def->kernel_id = last->kernel_id + 1;

   return 0;
}

/*}}}*/

/* Kernel table management */

static Kernel_Info_t *new_kernel_info (int hist_index) /*{{{*/
{
   Kernel_Info_t *p;

   p = (Kernel_Info_t *) ISIS_MALLOC (sizeof(Kernel_Info_t));
   if (p)
     {
        p->next = NULL;
        p->saved = NULL;
        p->kernel_params = NULL;
        p->hist_index = hist_index;
        p->kernel_id = 0;
     }

   return p;
}

/*}}}*/

static void free_kernel_info (Kernel_Info_t *p) /*{{{*/
{
   ISIS_FREE (p->kernel_params);
   ISIS_FREE (p);
}

/*}}}*/

Kernel_Info_t *find_kernel_info (Kernel_Table_t *t, int hist_index) /*{{{*/
{
   Kernel_Info_t *p = &t->kernel_info;
   Kernel_Info_t *last = p;

   while (p)
     {
        if (p->hist_index == hist_index)
          return p;

        last = p;
        p = p->next;
     }

   last->next = new_kernel_info (hist_index);
   return last->next;
}

/*}}}*/

void free_kernel_table (Kernel_Table_t *t) /*{{{*/
{
   Isis_Kernel_Def_t *def;
   Kernel_Info_t *p;

   def = t->kernel_defs;
   while (def)
     {
        if (def->malloced_kernel)
          break;
        def = def->next;
     }
   Fit_free_aux_kernels (def);

   p = &t->kernel_info;
   p = p->next;
   while (p)
     {
        Kernel_Info_t *tmp = p->next;
        free_kernel_info (p);
        p = tmp;
     }
}

/*}}}*/

typedef struct Static_Kernel_Table_t Static_Kernel_Table_t;
struct Static_Kernel_Table_t
{
   Isis_Kernel_Init_t *initfun;
   Isis_Kernel_Def_t def;
};

#define ISIS_KERNEL_NAME(s) Isis_##s##_kernel

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
}
#endif
extern int ISIS_KERNEL_NAME(std) (Isis_Kernel_Def_t *, char *);
extern int ISIS_KERNEL_NAME(pileup) (Isis_Kernel_Def_t *, char *);
extern int ISIS_KERNEL_NAME(yshift) (Isis_Kernel_Def_t *, char *);
extern int ISIS_KERNEL_NAME(gainshift) (Isis_Kernel_Def_t *, char *);
#if 0
{
#endif
#ifdef __cplusplus
}
#endif

static Static_Kernel_Table_t Static_Kernels [] =
{
  /* std kernel must be first */
  {ISIS_KERNEL_NAME(std),     ISIS_NULL_KERNEL_DEF},
  {ISIS_KERNEL_NAME(pileup),  ISIS_NULL_KERNEL_DEF},
  {ISIS_KERNEL_NAME(yshift),  ISIS_NULL_KERNEL_DEF},
  {ISIS_KERNEL_NAME(gainshift),  ISIS_NULL_KERNEL_DEF},
  {NULL,            ISIS_NULL_KERNEL_DEF}
};

static int init_static_kernels (Kernel_Table_t *t) /*{{{*/
{
   Static_Kernel_Table_t *sk = Static_Kernels;
   Isis_Kernel_Def_t **kd = &t->kernel_defs;
   unsigned int n = 0;

   while (sk->initfun)
     {
        if (-1 == (*sk->initfun)(&sk->def, NULL))
          return -1;

        *kd = &sk->def;
        (*kd)->kernel_id = n;
        (*kd)->malloced_kernel = 0;
        (*kd)->next = NULL;

        kd = &(*kd)->next;
        n++;
        sk++;
     }

   return 0;
}

/*}}}*/

int init_kernel_table (Kernel_Table_t *t) /*{{{*/
{
   Isis_Kernel_Def_t *def;

   memset ((char *)t, 0, sizeof (*t));

   if (-1 == init_static_kernels (t))
     return -1;

   for (def = t->kernel_defs; def != NULL; def = def->next)
     {
        if (-1 == Fit_add_kernel_function (def))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing fit-kernels");
             return -1;
          }
     }

   return 0;
}

/*}}}*/
