#ifndef ISIS_DBATOMIC_H
#define ISIS_DBATOMIC_H

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

/* $Id: db-atomic.h,v 1.13 2004/02/09 11:14:18 houck Exp $ */

#define SHORT_LINE_NAME       0
#define LONG_LINE_NAME        1

#define LONG_LINE_NAME_SIZE  64
#define SHORT_LINE_NAME_SIZE 48
#define LEVEL_NAME_SIZE      48

/* codes to indicate the blend type */ 

#define  SAME_ION   ((unsigned int)0x01)
#define  SAME_ELEM  ((unsigned int)0x02)

#define WAVELEN_TOL  1.e-5    /* Two wavelengths are "different" if 
                              * fabs(lam/lam0 - 1) >= WAVELEN_TOL */ 

/* I hope this is enough primes... */
#define PRIME_LIST \
     { \
           16381,    32749,    65521,    131071,    262139, \
          524287,  1048573,  2097143,   4194301,   8388593, \
        16777213, 33554393, 67108859, 134217689, 268435399, \
       536870909, 1073741789, UINT_MAX \
     }

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct _DB_t DB_t;

typedef struct _DB_ion_t DB_ion_t;
typedef struct _DB_level_t DB_level_t; 
typedef struct _DB_line_t DB_line_t;

typedef struct _DB_line_filter_t DB_line_filter_t;
typedef struct _DB_line_group_t DB_line_group_t;

extern int DB_Ion_Format;
#define FMT_CHARGE     0
#define FMT_ROMAN      1
#define FMT_INT_ROMAN  2

struct _DB_line_t
{
   double flux;            /* flux in current model */
   float wavelen;          /* Wavelength [angstrom] */
   float wavelen_err;      /* Wavelength uncertainty [angstrom] */
   float A;                /* Einstein A value [s^-1] */
   float A_err;            /* uncertainty in Einstein A value [s^-1] */
   int indx;               /* array index in database list */
   int upper_level;        /* energy level indices */
   int lower_level;
   int source_index;       /* index to literature source */
   int type_index;
   unsigned char ion_charge;
   unsigned char proton_number;
   unsigned char have_emissivity_data;
};

/* IMPORTANT:
 * The ground state is listed as level = 1 in the data files, 
 * but I'm storing it at array index 0 (in the array of 
 * DB_level_t structures).
 */
struct _DB_level_t
{
   float energy;                 /* excitation energy [eV] */
   float stat_weight;            /* statistical weight */
   float S;                      /* quantum numbers */
   int n, L;                     
   int indx;                     /* ground state = level 1, etc. */
   char label[LEVEL_NAME_SIZE];
   int ndown;
   DB_line_t **down;             /* list of downward transitions */
};

typedef struct
{
   int n;
   float *lambda, *lambda_err, *A_value, *A_err;
   int *upper_lev, *lower_lev;
   unsigned char *proton_number, *ion_charge;
}
DB_Merge_Type;

extern DB_t *DB_start (char **energy_level_files, char **wavelength_files);
extern void DB_end (DB_t *db);

/* general db */

extern int DB_sort_line_list (DB_t *db);
extern int DB_get_nlines (DB_t *db);
extern DB_line_t *DB_get_line_from_index (int indx, DB_t *db);
extern DB_line_t *DB_get_line (float wavelen, int proton_number, int ion_charge,
                              int upper_level, int lower_level, DB_t *db);
extern DB_line_t *DB_get_line_by_indices (int proton_number, int ion_charge,
                                        int upper_level, int lower_level, DB_t *db);
extern int _DB_get_element_name (char *name, int proton_number);
extern int _DB_get_element_Z (int *Z, char *name);
extern int DB_get_ion_name (char *name, int size, int Z, int q, int roman);
extern int DB_get_atomic_weight_amu (float *wt, DB_line_t *line);
extern int DB_print_branching_for_ion (FILE *fp, int Z, int q, DB_t *db);
extern int DB_merge_lines (DB_t *db, DB_Merge_Type *m);

/* element */

extern DB_ion_t *DB_get_ion (DB_t *db, int Z, int charge);

/* ion */

extern DB_level_t *DB_get_ion_level (int indx, DB_ion_t *ion);
extern int DB_get_ion_nlevels (int *nlevels, DB_ion_t *ion);

/* level */

extern int DB_get_level_label (char *label, int proton_number,
                               int ion_charge, int indx, DB_t *db);

/* Line filtering and classification */

extern DB_line_filter_t *DB_new_line_filter (void);
extern void DB_free_line_filter (DB_line_filter_t *filter);

typedef int (*DB_line_filter_function_t) (DB_line_t *line, DB_line_filter_t *f);

extern int DBf_trans (DB_line_t *line, DB_line_filter_t *f);
extern int DBf_set_trans (DB_line_filter_t *f, int Z, int q,
                          int nup, int *upper, int nlo, int *lower);

extern int DBf_el_ion (DB_line_t *line, DB_line_filter_t *f);
extern int DBf_set_el_ion (DB_line_filter_t *filter, 
                           int n_elem, int *proton_number,
                           int n_ion, int *ion_charge);

extern int DBf_wavelength (DB_line_t *line, DB_line_filter_t *f);
extern int DBf_set_wavelength (DB_line_filter_t *f,
                               double wavelen_min, double wavelen_max);

extern int DBf_flux (DB_line_t *line, DB_line_filter_t *f);
extern int DBf_set_flux (DB_line_filter_t *f, double flux_min, double flux_max);

extern int DB_apply_line_filter (char *t,
                                 DB_line_filter_function_t line_in_group,
                                 DB_line_filter_t *filter, DB_t *db);

extern int DB_list_to_group (int group, int *list, int n, DB_t *db);

extern int DB_get_group_members (int **index_list, int *nmembers, 
                                 int group, DB_t *db);

extern int DB_get_k_brightest (int *t, int *n, int k, 
                               int *list, int nlines, DB_t *db);

extern int DB_get_unblended (int *t, int *n, float f, float wl_f,
                             unsigned int type,
                             int *list, int nlines, DB_t *db);

/* line group */

extern void DB_free_line_group (DB_line_group_t *t);

extern void DB_delete_line_group (int group, DB_t *db);
extern DB_line_group_t *DB_get_group_table_head (DB_t *db);
extern DB_line_group_t *DB_next_group (DB_line_group_t *cl);
extern DB_line_group_t *DB_find_group (int group, DB_t *db);
extern int DB_edit_group (DB_line_group_t **g, int *list, int n, int add,
                          DB_t *db);

extern int DB_get_group_line_list (int *nlines, DB_line_t ***line, 
                                   DB_line_group_t *cl);
extern int DB_get_group_name (char **name, DB_line_group_t *cl);
extern int DB_get_group_index (int *group, DB_line_group_t *cl);
extern int DB_set_group_name (char *name, DB_line_group_t *cl);

extern char *DB_flag_array_from_list (int *list, int n, DB_t *db);
extern DB_line_group_t *DB_make_group_from_list (int group, int *list, int n,
                                                 DB_t *db);

/* line */

extern int DB_zero_line_flux (DB_t *db);
extern int DB_set_line_wavelength (float wavelen, float wavelen_err,
                                   DB_line_t *line);
extern int DB_get_line_ion (int *proton_number, int *charge, DB_line_t *line);
extern int DB_get_line_namestring (char *s, int size, int type, 
                                   DB_line_t *line, DB_t *db);
extern unsigned int DB_hash (float wl, int Z, int q, int upper, int lower,
                              unsigned int size);
extern unsigned int DB_hash2 (int Z, int q);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
