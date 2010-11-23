/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2010 Massachusetts Institute of Technology

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

/* $Id: model.h,v 1.5 2004/02/09 11:14:23 houck Exp $ */

#ifndef ISIS_MODEL_H
#define ISIS_MODEL_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct _Model_t Model_t;

enum
{
   MODEL_LINES_AND_CONTINUUM,
   MODEL_LINES,
   MODEL_CONTIN,
   MODEL_CONTIN_PSEUDO,
   MODEL_CONTIN_TRUE
};

struct _Model_t
{
   Model_t *next;          /* next component in linked list */
   double norm;            /* emission measure /(4*pi*D^2)*/
   float temperature;      /* electron temperature [K] */
   float density;          /* electron density [cm^-3] */
   float metal_abund;      /* abundance relative to cosmic */
   float rel_abund[ISIS_MAX_PROTON_NUMBER+1];
   float vturb;            /* turbulent velocity width [cm/s]*/
   float redshift;
   int id;                 /* id number of this component in the list */

   SLang_Array_Type *line_flux;       /* line fluxes */
   SLang_Array_Type *last_ionpop;
   /* last user-provided ionization balance */
};

typedef struct
{
   DB_t *db;
   EM_t *em;

   int contrib_flag;
   SLang_Array_Type *line_list;

   /* line emissivity modifier */
   SLang_Name_Type *line_emis_modifier;
   SLang_Array_Type *line_emis_modifier_params;
   Isis_Arg_Type *line_emis_modifier_args;
   SLang_Struct_Type *line_emis_modifier_qualifiers;

   /* user-defined line profile */
   Isis_Line_Profile_Type *profile;
   SLang_Array_Type *profile_params;
   char *profile_options;

   /* ionization balance hook */
   SLang_Name_Type *ionpop_modifier;
   SLang_Array_Type *ionpop_params;
   Isis_Arg_Type *ionpop_args;
   SLang_Struct_Type *ionpop_qualifiers;
}
Model_Info_Type;

extern int Model_start (void);
extern void Model_end (Model_t *model);

extern Model_t *Model_load_ascii_file (char *filename);
extern int Model_print_model (FILE *fp, Model_t *h);
extern void Model_set_profile_function (int idx);

extern Model_t *Model_add_component (Model_t *head, Model_t *x,
       unsigned int *elem, float *elem_abund, unsigned int num_elems);

extern int Model_spectrum (Model_t *h, Model_Info_Type *info,
                           double *wllo, double *wlhi, int nbins, double *val);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
