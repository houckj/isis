/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008   Massachusetts Institute of Technology

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

/* $Id: db-em.h,v 1.6 2004/06/06 02:42:21 houck Exp $ */

#ifndef ISIS_DBEM_H
#define ISIS_DBEM_H

/* This header represents the generic interface between ISIS and
 * the emissivity database.  It should be kept independent of the
 * details of any particular physical model.
 */

#define EM_PARAM_NAME_SIZE  32
#define EM_PARAM_UNITS_SIZE 32

#define ABUND_NAME_SIZE     16

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

enum
{
  EM_USE_MEMORY_DEFAULT=1
};

extern unsigned int EM_Use_Memory;
extern int EM_Maybe_Missing_Lines;
extern unsigned int EM_Hash_Table_Size_Hint;

typedef struct _EM_t EM_t;
typedef struct _EM_ioniz_table_t EM_ioniz_table_t;
typedef struct _EM_line_emis_t EM_line_emis_t;
typedef struct _EM_cont_type_t EM_cont_type_t;

typedef struct
{
   char *contin_emis;
   char *line_emis;
   char *ionization;
   char *abundance;
   char *filemap;
}
EM_File_Type;
#define NULL_EM_FILE_TYPE {NULL,NULL,NULL,NULL,NULL}   

typedef struct
{
   int Z, q;
   float *rel_abun;
}
EM_cont_select_t;

struct _EM_cont_type_t
{
   double *wlhi;         /* high edge of bin (angstrom) */
   double *wllo;         /* low edge of bin (angstrom) */
   double *true_contin;  /* true continuum */
   double *pseudo;       /* weak lines */
   int nbins;            /* number of bins */
};

extern EM_t *EM_start (EM_File_Type *f, void *cl, DB_t *db);
extern void EM_end (EM_t *em);

extern int EM_get_filemap (EM_t *em, char *emis_file, void *cl, unsigned int *num_hdus, double **temp, double **dens);

extern int EM_list_abundance_tables (FILE *fp, EM_t *em, int verbose);
extern int EM_set_chosen_abundance (EM_t *em, int k);
extern int EM_set_standard_abundance (EM_t *em, int k);
extern int EM_get_standard_abundance (EM_t *em);
extern int EM_get_abundance_table_index_by_name (EM_t *em, char *name);
extern int EM_get_abundance_table (EM_t *em, int k, char **name, float **a, int **z, int *n);
extern int EM_add_abundance_table (EM_t *em, char *name, float *abun, int *Z, int num_abun);

extern int EM_get_ion_fraction (float * frac, float * par_val,
                                int proton_number, int charge, int table_id,
                                EM_t *em);
extern int EM_load_alt_ionization_table (char *file, EM_t *em);
extern void EM_free_alt_ionization_table (EM_t *em);

extern int EM_get_nlines (EM_line_emis_t *t);
extern void EM_free_line_emis_list (EM_line_emis_t *p);
extern EM_line_emis_t *EM_get_line_spectrum (char *flag, float *par, EM_t *em);

extern EM_cont_type_t *EM_new_continuum (int nbins);
extern void EM_free_continuum (EM_cont_type_t *t);
extern int EM_get_continuum (EM_cont_type_t *cont, float *par, EM_cont_select_t *s, EM_t *em);

extern int EM_sum_line_emissivity (float *emis, float *par,
                                   int *list, int nlines, EM_t *em);

extern int _EM_get_line_emis_wl (DB_line_t **line, float *emis, float *wl, int k,
                                EM_line_emis_t *t);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
