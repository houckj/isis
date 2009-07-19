/* -*- mode: C; mode: fold -*- */

/*
    This file is part of ISIS, the Interactive Spectral Interpretation System
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

/* $Id: db-cie.c,v 1.18 2004/06/06 20:56:06 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "cfits.h"
#include "errors.h"
#include "db-atomic.h"
#include "db-em.h"
#include "db-cie.h"

#undef MAX
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

#undef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/*
 * NIST 1998 CODATA recommended values of physical constants:
 */

#define PLANCK      ((double) 6.62606876e-27)       /* Planck's constant (erg s) */
#define BOLTZ       ((double) 1.3806503e-16)        /* Boltzmann's constant (erg/K) */
#define CLIGHT      ((double) 2.99792458e10)        /* speed of light (cm/s) */
#define AMU         ((double) 1.66053873e-24)       /* atomic mass unit (g) */
#define ERG_PER_EV  ((double) 1.602176462e-12)

#define ERG_ANGSTROM    (PLANCK * CLIGHT / 1.e-8)
#define KEV_ANGSTROM    (ERG_ANGSTROM / ERG_PER_EV / 1.e3)

/*}}}*/

/*{{{ globals & private data type definitions */

#define DENSITY_TOL   1.e-4
#define TEMP_TOL      1.e-4

enum
{
   EM_TEMP = 0,               /* temperature */
   EM_DENS = 1                /* density */
};

unsigned int EM_Use_Memory = EM_USE_MEMORY_DEFAULT;
/* EM_Use_Memory is a bitmap.  See set_memory_usage_level() for details */

static int EM_Load_Cont_Emis;
static int EM_Load_Line_Emis;

/*EM_Load_Line_Emis != 0 means line emissivities all loaded into memory at once
 *                       (the emissivity database is "memory resident")
 *                   = 0 means line emissivities are loaded from disk as needed
 *                       (the emissivity database is "disk resident")
 */

int EM_Maybe_Missing_Lines = 1;
/*EM_Maybe_Missing_Lines != 0 means that the atomic database might not contain
 *                          every line listed in the emissivity tables.
 *                        = 0 means that every line that appears in the emissivity
 *                          tables is known to appear in the atomic data tables.
 */

typedef struct EM_filemap_t EM_filemap_t;
typedef struct _EM_line_data_t EM_line_data_t;
typedef struct _EM_cont_data_t EM_cont_data_t;
typedef struct _EM_cont_emis_t EM_cont_emis_t;
typedef struct _EM_ionfrac_t EM_ionfrac_t;
typedef struct _EM_abund_t EM_abund_t;

struct _EM_t
{
   DB_t *db;
   EM_ioniz_table_t *ioniz_table[2];
   EM_line_data_t *line_data;
   EM_cont_data_t *cont_data;
   EM_abund_t *abund;
   int chosen_abund_table;       /* user-specified abund table */
   int standard_abund_table;     /* the standard abund table */
};

struct EM_filemap_t
{
   float *densities;
   float *temps;
   int *hdu;
   int num_hdus, num_temps, num_densities;
   char filename[CFLEN_FILENAME];
   char abund_table[CFLEN_KEYWORD];
};

struct _EM_line_data_t
{
   EM_filemap_t *map;
   EM_line_emis_t **emis;   /* storage for memory resident mode */
};

struct _EM_line_emis_t
{
   DB_line_t **line;          /* vector of ptrs to atomic data for each line */
   float *emissivity;         /* vector of line emissivities */
   int *lookup;
   float par[2];              /* e.g (T, density) values for these emissivities */
   int nlines;
};
/* contains all the line emissivities for e.g. a given (T, density) pair */

struct _EM_cont_data_t
{
   EM_filemap_t *map;
   EM_cont_emis_t **emis;  /* storage for memory resident mode */
};

struct _EM_cont_emis_t        /* mirrors structure of FITS file extension */
{                             /* units are photons cm^3 s^-1 Angstrom^-1 */
   double *g_true_contin, *true_contin;             /* true continuum */
   double *g_pseudo, *pseudo;                       /* weak lines */
   int ntrue_contin, npseudo;
   double temp, dens;
   int Z, q;
   EM_cont_emis_t *next;
};
/*  Z==0  q==-1 (rmJ=0) means  node contains sum over elements/ions.
 *  Z >0  q==-1 (rmJ=0) means  node contains sum over ions of a single element.
 *  Z >0  q>= 0 (rmJ>0) means  node contains values for a single ion.
 */

struct _EM_ionfrac_t
{
   float *fraction;    /* packed vector of ionization fractions */
   float par[2];       /* e.g (T, density) values for this ioniz. structure */
};

struct _EM_ioniz_table_t
{
   EM_ionfrac_t **ionfrac;
   int *offset;          /* offset[i] is offset to ith element in each ionfrac vector */
   int num_td_pairs;     /* number of temp/density grid points in ionfrac */
};

struct _EM_abund_t
{
   EM_abund_t *next;
   float abundance[ISIS_MAX_PROTON_NUMBER+1];
   char name[ABUND_NAME_SIZE];
   int abund_id;
};

typedef struct
{
   float temperature;
   float density;
   float *epsilon;
   float *lambda;
   int *Z;
   int *rmJ;
   int *up;
   int *lo;
}
Load_Linefile_Type;

/* abundance[Z] = log abundance of element with proton_number = Z relative
 * to Hydrogen, with Hydrogen = 12.00 */

/*}}}*/

/*{{{ abundances */

static void free_abund_table (EM_abund_t *t) /*{{{*/
{
   ISIS_FREE (t);
}

/*}}}*/

static void free_abund_list (EM_abund_t *head) /*{{{*/
{
   while (head)
     {
        EM_abund_t *t = head->next;
        free_abund_table (head);
        head = t;
     }        
}

/*}}}*/

static EM_abund_t *new_abund_table (void) /*{{{*/
{
   EM_abund_t *t;
   if (NULL != (t = (EM_abund_t *) ISIS_MALLOC(sizeof *t)))
     {
        memset ((char *)t, 0, sizeof *t);
     }   
   return t;
}

/*}}}*/

static EM_abund_t *make_abundance_table (char *name, float *abun, int *Z, int num_abun) /*{{{*/
{
   EM_abund_t *t;
   int i;
   
   if (NULL == name || NULL == abun || NULL == Z)
     return NULL;
   
   if (NULL == (t = new_abund_table ()))
     return NULL;
   
   isis_strcpy (t->name, name, ABUND_NAME_SIZE);
   for (i = 0; i < num_abun; i++)
     {
        if (0 <= Z[i] && Z[i] <= ISIS_MAX_PROTON_NUMBER)
          t->abundance[Z[i]] = abun[i];
     }
   
   return t;
}

/*}}}*/

int EM_add_abundance_table (EM_t *em, char *name, float *abun, int *Z, int num_abun) /*{{{*/
{
   EM_abund_t *t, *a;
   
   if (em == NULL)
     return -1;
   
   if (NULL == (t = make_abundance_table (name, abun, Z, num_abun)))
     return -1;
   
   if (em->abund)
     {
        for (a = em->abund; a != NULL; a = a->next)
          {
             if (a->next == NULL)
               {
                  a->next = t;
                  t->abund_id = a->abund_id + 1;
                  break;
               }        
          }
     } 
   else em->abund = t;
   
   return t->abund_id;
}

/*}}}*/

static void free_string_array (char **p, int n) /*{{{*/
{
   if (p == NULL)
     return;

   while (n-- > 0)
     ISIS_FREE (p[n]);

   ISIS_FREE (p);
}

/*}}}*/

static char **new_string_array (int n, int len) /*{{{*/
{
   char **p;
   int i;

   if (n <= 0 || len <= 0)
     return NULL;

   if (NULL == (p = (char **) ISIS_MALLOC (n * sizeof(char *))))
     return NULL;
   memset ((char *)p, 0, n * sizeof(char *));

   for (i=0; i < n; i++)
     {
        if (NULL == (p[i] = (char *) ISIS_MALLOC (len * sizeof(char))))
          {
             free_string_array (p, n);
             return NULL;
          }
     }

   return p;
}

/*}}}*/

static EM_abund_t *load_abundance_file (char *path) /*{{{*/
{
   char **tmp_name = NULL;
   float *tmp_abun = NULL;
   cfitsfile * fp = NULL;
   EM_abund_t *t = NULL;
   EM_abund_t *head = NULL;
   int i, Z, nabund, ok = 0;

   if (NULL == path)
     return NULL;

   if (NULL == (fp = cfits_open_file_readonly (path)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", path);
        goto close_and_return;
     }

   if (-1 == cfits_movabs_hdu (2, fp))
     {
        isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__, "hdu=2, %s", path);
        goto close_and_return;
     }

   if (-1 == cfits_read_int_keyword (&nabund, "NAXIS2", fp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "NAXIS2 => %s", path);
        goto close_and_return;
     }

   if (NULL == (tmp_abun = (float *) ISIS_MALLOC (nabund * sizeof(float)))
       || NULL == (tmp_name = new_string_array (nabund, ABUND_NAME_SIZE)))
     goto close_and_return;

   for (i = 0; i < nabund; i++)
     {
        if (NULL == (t = new_abund_table ()))
          goto close_and_return;        
        t->abund_id = nabund - i - 1;
        t->next = head;
        head = t;
     }

   if (-1 == cfits_read_string_col (tmp_name, nabund, 1L, "Source", fp))
     {
        isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "Source => %s", path);
        goto close_and_return;
     }

   t = head;
   for (i = 0; i < nabund; i++)
     {
        isis_strcpy (t->name, tmp_name[i], ABUND_NAME_SIZE);
        t = t->next;
     }   

   for (Z = 1; Z <= ISIS_MAX_PROTON_NUMBER; Z++)
     {
        char el_name[3];

        if (-1 == _DB_get_element_name (el_name, Z))
          goto close_and_return;

        if (-1 == cfits_read_float_col (tmp_abun, nabund, 1L, el_name, fp))
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "%s", path);
             goto close_and_return;
          }

        t = head;
        for (i = 0; i < nabund; i++)
          {
             t->abundance[Z] = tmp_abun[i];
             t = t->next;
          }        
     }

   isis_vmesg (INFO, I_READ_OK, __FILE__, __LINE__, "%s", path);
   ok = 1;

   close_and_return:

   ISIS_FREE (tmp_abun);
   free_string_array (tmp_name, nabund);

   (void) cfits_close_file (fp);

   if (!ok)
     {
        free_abund_list (head);
        head = NULL;
     }   

   return head;
}

/*}}}*/

static int have_abundance_tables (EM_t *em) /*{{{*/
{
   if (NULL == em)
     return 0;

   return (NULL != em->abund);
}

/*}}}*/

int EM_list_abundance_tables (FILE *fp, EM_t *em, int verbose) /*{{{*/
{
   EM_abund_t *t;

   if (NULL == em)
     return -1;

   if (0 == have_abundance_tables (em))
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "abundance tables not loaded?");
        return 0;
     }

   fprintf (fp, "Abundance Tables:\n");
   if (verbose)
     fprintf (fp, "Log cosmic abundance by number, relative to Hydrogen = 12.0\n");

   for (t = em->abund; t != NULL; t = t->next)
     {
        int Z;
        if (fprintf (fp, "%d  %s %s\n",
                     t->abund_id,
                     t->name,
                     (t->abund_id == em->chosen_abund_table) ? "<Current Definition>" : "")
            < 0)
          return -1;

        if (verbose == 0)
          continue;

        for (Z = 1; Z <= ISIS_MAX_PROTON_NUMBER; Z++)
          {
             char el_name[3];
             if (-1 == _DB_get_element_name (el_name, Z))
               return -1;
             if (fprintf (fp, " %2s: %7.4f%c",
                          el_name, t->abundance[Z],
                          (Z % 6) ? ' ' : '\n')
                 < 0)
               return -1;
          }
        fputc ('\n', fp);
     }

   return 0;
}

/*}}}*/

static EM_abund_t *find_abund_table (EM_t *em, int k) /*{{{*/
{
   EM_abund_t *t;
   
   if (em == NULL)
     return 0;
   
   for (t = em->abund; t != NULL; t = t->next)
     {
        if (t->abund_id == k)
          return t;
     }
   
   return NULL;
}

/*}}}*/

int EM_set_chosen_abundance (EM_t *em, int k) /*{{{*/
{
   if (NULL == find_abund_table (em, k))
     return -1;

   em->chosen_abund_table = k;
   return 0;
}
/*}}}*/

int EM_set_standard_abundance (EM_t *em, int k) /*{{{*/
{
   if (NULL == find_abund_table (em, k))
     return -1;

   em->standard_abund_table = k;
   return 0;
}
/*}}}*/

int EM_get_standard_abundance (EM_t *em) /*{{{*/
{
   if (!have_abundance_tables (em))
     return -1;

   return em->standard_abund_table;
}
/*}}}*/

int EM_get_abundance_table_index_by_name (EM_t *em, char *name) /*{{{*/
{
   EM_abund_t *t;

   if ((NULL == em) || (name == NULL))
     return -1;

   if (0 == have_abundance_tables (em))
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "abundance tables not loaded?");
        return -1;
     }

   for (t = em->abund; t != NULL; t = t->next)
     {
        if (0 == isis_strcasecmp (t->name, name))
          return t->abund_id;
     }

   return -1;
}

/*}}}*/

int EM_get_abundance_table (EM_t *em, int k, char **name, float **a, int **z, int *n) /*{{{*/
{
   EM_abund_t *t;
   float *abund;
   int zz, num = ISIS_MAX_PROTON_NUMBER;

   if (a == NULL || z == NULL || n == NULL)
     return -1;

   if ((0 == have_abundance_tables(em))
        || (NULL == (t = find_abund_table (em, k))))
     return -1;
   
   *name = NULL;
   *a = NULL;
   *z = NULL;

   if ((NULL == (*name = isis_make_string (t->name)))
       || (NULL == (*a = (float *) ISIS_MALLOC (num * sizeof(float))))
       || (NULL == (*z = (int *) ISIS_MALLOC (num * sizeof(int)))))
     {
        ISIS_FREE (*z);
        ISIS_FREE (*a);
        ISIS_FREE (*name);
        return -1;
     }   

   abund = t->abundance;
   
   for (zz = 1; zz <= num; zz++)
     {
        (*a)[zz-1] = abund[zz];
        (*z)[zz-1] = zz;
     }

   *n = num;

   return 0;
}

/*}}}*/

static int set_standard_abund_table (EM_t *em) /*{{{*/
{
   char *abund_line, *abund_cont, *std;
   int k, status;

   if ((em == NULL) || (0 == have_abundance_tables(em)))
     return 0;

   abund_line = abund_cont = std = NULL;

   if (em->line_data && em->line_data->map)
     abund_line = em->line_data->map->abund_table;

   if (em->cont_data && em->cont_data->map)
     abund_cont = em->cont_data->map->abund_table;

   if (abund_line && abund_cont)
     {
        if (0 == isis_strcasecmp (abund_line, abund_cont))
          std = abund_line;
        else
          {
             isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "abundance table mismatch:  lines->'%s',  continua->'%s'",
                         abund_line, abund_cont);
             return -1;
          }
     }
   else if (abund_line)
     {
        std = abund_line;
     }
   else if (abund_cont)
     {
        std = abund_cont;
     }

   k = EM_get_abundance_table_index_by_name (em, std);

   /* default chosen to standard */
   status = EM_set_standard_abundance (em, k);
   EM_set_chosen_abundance (em, k);

   return status;
}

/*}}}*/

static int use_alt_abund (EM_t *em) /*{{{*/
{
   if (NULL == em
       || (em->standard_abund_table == em->chosen_abund_table))
     return 0;
   else
     return 1;
}
/*}}}*/

static int get_abundance_factor (float *f_abund, EM_t *em) /*{{{*/
{
   EM_abund_t *s, *a;
   int Z;

   if ((NULL == f_abund) || (0 == have_abundance_tables (em)))
     return -1;

   s = find_abund_table (em, em->standard_abund_table);
   a = find_abund_table (em, em->chosen_abund_table);

   if ((s == NULL) || (a == NULL))
     return -1;

   /* abundance[Z] = log abundance of element with proton_number = Z relative
    * to Hydrogen, with Hydrogen = 12.00 */

   f_abund[0] = 0.0;
   for (Z = 1; Z <= ISIS_MAX_PROTON_NUMBER; Z++)
     {
        float xp = a->abundance[Z] - s->abundance[Z];
        f_abund[Z] = pow (10.0, xp);
     }

   return 0;
}

/*}}}*/

static int scale_line_abundance (EM_line_emis_t *t, EM_t *em) /*{{{*/
{
   float f_abund[ISIS_MAX_PROTON_NUMBER+1];
   int k;

   if (!use_alt_abund (em))
     return 0;

   if (-1 == get_abundance_factor (f_abund, em))
     return -1;

   for (k = 0; k < t->nlines; k++)
     {
        int Z, q;
        if (-1 == DB_get_line_ion (&Z, &q, t->line[k]))
          return -1;
        t->emissivity[k] *= f_abund[Z];
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ ionization balance tables */

static void free_ionfrac (EM_ionfrac_t *p) /*{{{*/
{
   if (p == NULL)
     return;

   ISIS_FREE (p->fraction);
   ISIS_FREE (p);
}

/*}}}*/

static void free_ioniz_table (EM_ioniz_table_t *p) /*{{{*/
{
   int i;

   if (p == NULL)
     return;

   if (p->ionfrac)
     {
        for (i=0; i < p->num_td_pairs; i++)
          free_ionfrac (p->ionfrac[i]);
        ISIS_FREE (p->ionfrac);
     }

   ISIS_FREE (p->offset);
   ISIS_FREE (p);
}
/*}}}*/

static EM_ioniz_table_t *new_ioniz_table (int num_ions, int num_td_pairs) /*{{{*/
{
   EM_ioniz_table_t *t = NULL;
   int i;

   if (num_td_pairs <= 0
       || num_ions <= 0)
     return NULL;

   if (NULL == (t = (EM_ioniz_table_t *) ISIS_MALLOC (sizeof(EM_ioniz_table_t))))
     return NULL;
   memset ((char *)t, 0, sizeof (*t));

   t->num_td_pairs = num_td_pairs;

   if (NULL == (t->offset = (int *) ISIS_MALLOC ((ISIS_MAX_PROTON_NUMBER + 1) * sizeof(int)))
       || NULL == (t->ionfrac = (EM_ionfrac_t **) ISIS_MALLOC (num_td_pairs * sizeof(EM_ionfrac_t *))))
     goto free_and_return;
   memset ((char *)t->ionfrac, 0, num_td_pairs * sizeof (EM_ionfrac_t *));

   for (i=0; i < num_td_pairs; i++)
     {
        EM_ionfrac_t *p = NULL;

           if (NULL == (p = (EM_ionfrac_t *) ISIS_MALLOC (sizeof(EM_ionfrac_t))))
          goto free_and_return;
        memset ((char *)p, 0, sizeof (*p));

        if (NULL == (p->fraction = (float *) ISIS_MALLOC (num_ions * sizeof(float))))
          {
             free_ionfrac (p);
             goto free_and_return;
          }

        t->ionfrac[i] = p;
     }

   return t;

   free_and_return:
   free_ioniz_table (t);
   return NULL;
}
/*}}}*/

static int compute_offsets (int *element_offset, int *Z, int num_elements) /*{{{*/
{
   int j;
   int off;

   if (NULL == element_offset || NULL == Z)
     return -1;

   for (j=0; j <= ISIS_MAX_PROTON_NUMBER; j++)
     element_offset[j] = -1;

   off = 0;
   for (j=0; j < num_elements; j++)
     {
        if (Z[j] > ISIS_MAX_PROTON_NUMBER)
          {
             isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__,
                         "truncating ionization balance table\n\t(number of elements in ionization table exceeds\n\t ISIS_MAX_PROTON_NUMBER=%d)",
                         ISIS_MAX_PROTON_NUMBER);
             break;
          }
        element_offset[ Z[j] ] = off;
        off += Z[j] + 1;
     }

   return 0;
}

/*}}}*/

static EM_ioniz_table_t *load_ionization_table (char * filename) /*{{{*/
{
   EM_ioniz_table_t *t = NULL;
   cfitsfile * fp;
   int ret = -1;
   int i, ntemp, ndens, num_td_pairs;
   int num_elements;              /* # of elements for which data exists */
   int num_ions;                  /* # of ions per temp/density pair*/
   int *Z = NULL;                 /* array of proton numbers for which data exists */

   if (NULL == filename)
     return NULL;

   if (NULL == (fp = cfits_open_file_readonly (filename)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", filename);
        goto close_and_return;
     }

   /*   (void) cfits_test_keyword ("ION_BAL", "HDUCLAS1", fp); */

   if (-1 == cfits_movabs_hdu (2, fp))
     {
        isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__, "hdu=2, %s", filename);
        goto close_and_return;
     }

   if (-1 == cfits_read_int_keyword (&ntemp, "T_NUMBER", fp)
       || -1 == cfits_read_int_keyword (&ndens, "N_NUMBER", fp)
       || -1 == cfits_read_int_keyword (&num_elements, "N_ELEMEN", fp)
       || -1 == cfits_read_int_keyword (&num_ions, "N_IONS", fp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "%s", filename);
        goto close_and_return;
     }

   num_td_pairs = ntemp * ndens;

   if (NULL == (t = new_ioniz_table (num_ions, num_td_pairs)))
     goto close_and_return;

   /* Read the atomic numbers vector from the first row only,
      assuming all rows are identical */

   if (NULL == (Z = (int *) ISIS_MALLOC (num_elements * sizeof(int)))
       || -1 == cfits_read_int_col (Z, num_elements, 1, "Z_ELEMENT", fp))
     {
        isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "Z_ELEMENT => %s", filename);
        goto close_and_return;
     }

   if (-1 == compute_offsets (t->offset, Z, num_elements))
     goto close_and_return;

   for (i=0; i < num_td_pairs; i++)
     {
        EM_ionfrac_t *p;
        int start_row = i+1;

        p = t->ionfrac[i];

        if (-1 == cfits_read_float_col(&p->par[EM_TEMP], 1, start_row, "Temperature", fp)
            || -1 == cfits_read_float_col(&p->par[EM_DENS], 1, start_row, "Density", fp)
            || -1 == cfits_read_float_col(p->fraction, num_ions, start_row, "X_IONPOP", fp))
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "%s", filename);
             goto close_and_return;
          }
     }

   ret = 0;
   isis_vmesg (INFO, I_READ_OK, __FILE__, __LINE__, "%s", filename);

   close_and_return:

   if (ret)
     {
        free_ioniz_table (t);
        t = NULL;
     }

   (void) cfits_close_file(fp);
   ISIS_FREE (Z);

   return t;
}
/*}}}*/

static int get_ion_fraction (float *frac, float *par, int Z, int q, EM_ioniz_table_t *t) /*{{{*/
{
   float etemp, edens;
   int i, n;
   int off;

   if (NULL == t)
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "ionization table not loaded?");
        return -1;
     }

   if (Z < 1 || Z > ISIS_MAX_PROTON_NUMBER
       || q < 0 || q > Z)
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "ion Z = %d, q = %d", Z, q);
        return -1;
     }

   if (t->offset[ Z ] < 0)
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "Z = %d ionization", Z);
        return -1;
     }

   off = t->offset[ Z ] + q;

   etemp = par[ EM_TEMP ];
   edens = par[ EM_DENS ];          /* FIXME: ignoring density dependence of ionization */
   (void) edens;

   n = t->num_td_pairs;

   for (i=0; i < n-1; i++)
     {
        EM_ionfrac_t *p  = t->ionfrac[i  ];
        EM_ionfrac_t *pn = t->ionfrac[i+1];
        float t1 = p->par[ EM_TEMP ];
        float t2 = pn->par[ EM_TEMP ];
        float x, f1, f2;

        if (etemp < t1 || t2 <= etemp)
          continue;

        x = (etemp - t1) / (t2 - t1);

        f1 = p->fraction[ off ];
        f2 = pn->fraction[ off ];

        *frac = (1.0 - x) * f1 + x * f2;
        return 0;
     }
                            /* temperature not found */
   *frac = 0.0;
   isis_vmesg (FAIL, I_RANGE_ERROR, __FILE__, __LINE__, 
               "%11.4e K out of range [%11.4e, %11.4e]",
               etemp,
               t->ionfrac[0  ]->par[ EM_TEMP ],
               t->ionfrac[n-1]->par[ EM_TEMP ]);
   return -1;
}
/*}}}*/

int EM_get_ion_fraction (float *frac, float *par, int Z, int q, int k, EM_t *em) /*{{{*/
{
   if (k != 0 && k != 1)
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "table index %d", k);
        return -1;
     }

   return get_ion_fraction (frac, par, Z, q, em->ioniz_table[k]);
}

/*}}}*/

static int use_alt_ioniz (EM_t *em) /*{{{*/
{
   if (NULL == em
       || NULL == em->ioniz_table[0]
       || NULL == em->ioniz_table[1])
     return 0;
   else
     return 1;
}

/*}}}*/

static int get_ioniz_factor (float f_ioniz[ISIS_MAX_PROTON_NUMBER+1][ISIS_MAX_PROTON_NUMBER+1], /*{{{*/
                             float * par, EM_ioniz_table_t *t_new, EM_ioniz_table_t *t_old)
{
   float f_old, f_new;
   int size = (ISIS_MAX_PROTON_NUMBER+1)*(ISIS_MAX_PROTON_NUMBER+1);
   int Z;

   if (NULL == t_old || NULL == t_new || par == NULL)
     return -1;

   memset ((char *) f_ioniz, 0, size * sizeof(float));

   for (Z = 1; Z <= ISIS_MAX_PROTON_NUMBER; Z++)
     {
        int q;
        if (-1 == t_old->offset[Z]
            || -1 == t_new->offset[Z])
          {
             for (q = 0; q <= Z; q++)
               f_ioniz[Z][q] = 1.0;
          }
        else
          {
             for (q = 0; q <= Z; q++)
               {
                  if (0 == get_ion_fraction (&f_old, par, Z, q, t_old)
                      && 0 == get_ion_fraction (&f_new, par, Z, q, t_new)
                      && f_old > 0.0)
                    f_ioniz[Z][q] = f_new / f_old;
                  else
                    f_ioniz[Z][q] = 1.0;
               }
          }
     }

   return 0;
}

/*}}}*/

static int scale_line_ionization (EM_line_emis_t *t, float * par, EM_t *em) /*{{{*/
{
   EM_ioniz_table_t *t_old;
   EM_ioniz_table_t *t_new;
   float f_ioniz[ISIS_MAX_PROTON_NUMBER+1][ISIS_MAX_PROTON_NUMBER+1];
   int k;

   if (NULL == em || NULL == t || par == NULL
       || em->ioniz_table == NULL)
     return -1;

   t_old = em->ioniz_table[0];
   t_new = em->ioniz_table[1];

   if (NULL == t_old || NULL == t_new)
     return 0;

   if (-1 == get_ioniz_factor (f_ioniz, par, t_new, t_old))
     return -1;

   for (k = 0; k < t->nlines; k++)
     {
        int Z, q;
        if (-1 == DB_get_line_ion (&Z, &q, t->line[k]))
          return -1;
        t->emissivity[k] *= f_ioniz[Z][q];
     }

   return 0;
}
/*}}}*/

void EM_free_alt_ionization_table (EM_t *em) /*{{{*/
{
   if (NULL == em || NULL == em->ioniz_table || NULL == em->ioniz_table[1])
     return;

   free_ioniz_table (em->ioniz_table[1]);
   em->ioniz_table[1] = NULL;
}

/*}}}*/

int EM_load_alt_ionization_table (char *file, EM_t *em) /*{{{*/
{
   EM_ioniz_table_t *t;

   if (NULL == em || NULL == file)
     return -1;

   if (NULL == (t = load_ionization_table (file)))
     return -1;

   free_ioniz_table (em->ioniz_table[1]);
   em->ioniz_table[1] = t;

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ Merge emission lines into atomic database */

typedef struct
{
   float lambda;                 /* Wavelength [angstrom] */
   int Z, q;                     /* proton number, ion charge */
   int up, lo;                   /* energy level indices */
   int is_new;
}
Line_t;

typedef struct
{
   Line_t *table;
   unsigned int *entries;
   unsigned int num_entries;
   unsigned int size;
}
Hash_t;

enum
{
   MAX_ALLOWED_MISSES = 32,
   DEFAULT_HASH_TABLE_SIZE_HINT = 50000
};

unsigned int EM_Hash_Table_Size_Hint = DEFAULT_HASH_TABLE_SIZE_HINT;
static Hash_t *Hash_Table;

static void free_hash_table (Hash_t *h) /*{{{*/
{
   if (h == NULL)
     return;
   ISIS_FREE (h->table);
   ISIS_FREE (h->entries);
   ISIS_FREE (h);
}
/*}}}*/

static Hash_t *new_hash_table (unsigned int size) /*{{{*/
{
   Hash_t *h;

   if (NULL == (h = (Hash_t *) ISIS_MALLOC (sizeof(Hash_t))))
     return NULL;
   memset ((char *)h, 0, sizeof (*h));

   h->table = (Line_t *) ISIS_MALLOC (size * sizeof(Line_t));
   h->entries = (unsigned int *) ISIS_MALLOC (size * sizeof(unsigned int));
   if ((h->table == NULL) || (h->entries == NULL))
     {
        free_hash_table (h);
        return NULL;
     }

   memset ((char *)h->table, 0, size * sizeof(Line_t));
   memset ((char *)h->entries, 0, size * sizeof(unsigned int));
   h->size = size;
   h->num_entries = 0;

   return h;
}
/*}}}*/

static int init_hash_table (unsigned int size) /*{{{*/
{
   static unsigned int prime[] = PRIME_LIST;
   int k;

   for (k = 0; prime[k] < size; k++)
     ;
   size = prime[k];

   Hash_Table = new_hash_table (size);
   if (Hash_Table == NULL)
     return -1;

   return 0;
}

/*}}}*/

static int start_hashing (void ***vp, int nlines) /*{{{*/
{
   (void) vp; (void) nlines;
   return 0;
}

/*}}}*/

static int finish_hashing (void **vp, int status) /*{{{*/
{
   (void) vp;
   return status;
}

/*}}}*/

static int do_hashing (void **vp, int nread, int start_row, Load_Linefile_Type *x, DB_t *db) /*{{{*/
{
   int i;

   (void) vp; (void) db; (void) start_row;

   for (i = 0; i < nread; i++)
     {
        unsigned int h, step;
        int q, miss;
        Line_t *p, *ph;

        q = x->rmJ[i] - 1;

        h = DB_hash (x->lambda[i], x->Z[i], q, x->up[i], x->lo[i], Hash_Table->size);
        step = DB_hash2 (x->Z[i], q);

        miss = 0;
        p = Hash_Table->table;

        for (;;)
          {
             int seen;

             ph = &p[h];

             if (ph->Z == 0)
               break;

             seen = ((ph->Z == x->Z[i]) && (ph->q == q)
                     && (ph->up == x->up[i]) && (ph->lo == x->lo[i])
                     && (ph->lambda == x->lambda[i]));

             if (seen)
               goto already_seen;

             h = (h + step) % Hash_Table->size;
             if (miss++ > MAX_ALLOWED_MISSES)
               {
                  /* Try increasing EM_Hash_Table_Size_Hint */
                  isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, 
                              "hash table overflow\n Try setting EM_Hash_Table_Size_Hint [currently = %d]\n to a value at least 25%% larger than the number of lines\n in your database\n",
                              EM_Hash_Table_Size_Hint);
                  return -1;
               }
          }

        ph->Z = x->Z[i];
        ph->q = q;
        ph->lambda = x->lambda[i];
        ph->up = x->up[i];
        ph->lo = x->lo[i];
        ph->is_new = 1;     /* assume not in the database */

        Hash_Table->entries[ Hash_Table->num_entries ] = h;
        Hash_Table->num_entries++;

        already_seen: ;
     }

   return 0;
}

/*}}}*/

static void free_merge_space (DB_Merge_Type *m) /*{{{*/
{
   if (m == NULL)
     return;
   ISIS_FREE (m->lambda); ISIS_FREE (m->lambda_err);
   ISIS_FREE (m->A_value); ISIS_FREE (m->A_err);
   ISIS_FREE (m->upper_lev); ISIS_FREE (m->lower_lev);
   ISIS_FREE (m->proton_number); ISIS_FREE (m->ion_charge);

   memset ((char *)m, 0, sizeof(*m));
}

/*}}}*/

static int allocate_merge_space (DB_Merge_Type *m, unsigned int n) /*{{{*/
{
   if (m == NULL)
     return -1;

   memset ((char *)m, 0, sizeof(*m));

   if (NULL == (m->lambda = (float *) ISIS_MALLOC (n * sizeof(float)))
       || NULL == (m->lambda_err = (float *) ISIS_MALLOC (n * sizeof(float)))
       || NULL == (m->A_value = (float *) ISIS_MALLOC (n * sizeof(float)))
       || NULL == (m->A_err = (float *) ISIS_MALLOC (n * sizeof(float)))
       || NULL == (m->upper_lev = (int *) ISIS_MALLOC (n * sizeof(int)))
       || NULL == (m->lower_lev = (int *) ISIS_MALLOC (n * sizeof(int)))
       || NULL == (m->proton_number = (unsigned char *) ISIS_MALLOC (n * sizeof(unsigned char)))
       || NULL == (m->ion_charge = (unsigned char *) ISIS_MALLOC (n * sizeof(unsigned char))))
     {
        free_merge_space (m);
        return -1;
     }

   m->n = 0;

   return 0;
}
/*}}}*/

static void copy_to_merge_space (DB_Merge_Type *m, Line_t *p) /*{{{*/
{
   int n;
   if (m == NULL || p == NULL)
     return;

   n = m->n;
   m->lambda[n] = p->lambda;
   m->lambda_err[n] = 0.0;
   m->A_value[n] = 0.0;
   m->A_err[n] = 0.0;
   m->proton_number[n] = (unsigned char) p->Z;
   m->ion_charge[n] = (unsigned char) p->q;
   m->upper_lev[n] = p->up;
   m->lower_lev[n] = p->lo;
   m->n++;
}
/*}}}*/

static int merge_hash_table (DB_t *db) /*{{{*/
{
   DB_Merge_Type m;
   unsigned int i, h, num_new_lines;
   Line_t *p = Hash_Table->table;
   Line_t *ph;
   int ret = -1;

   num_new_lines = 0;
   for (i = 0; i < Hash_Table->num_entries; i++)
     {
        h = Hash_Table->entries[i];
        ph = &p[h];
        if (NULL == DB_get_line (ph->lambda, ph->Z, ph->q, ph->up, ph->lo, db))
          num_new_lines++;
        else ph->is_new = 0;
     }

   if (num_new_lines == 0)
     return 0;

   if (-1 == allocate_merge_space (&m, num_new_lines))
     return -1;

   for (i = 0; i < Hash_Table->num_entries; i++)
     {
        h = Hash_Table->entries[i];
        ph = &p[h];
        if (ph->is_new)
          copy_to_merge_space (&m, ph);
     }

   ret = DB_merge_lines (db, &m);
   free_merge_space (&m);

   return ret;
}
/*}}}*/

static void deallocate_hash_table (void) /*{{{*/
{
   free_hash_table (Hash_Table);
}
/*}}}*/

/*}}}*/

/*{{{ file map definition */

static void free_filemap (EM_filemap_t *map) /*{{{*/
{
   if (NULL == map)
     return;

   ISIS_FREE (map->densities);
   ISIS_FREE (map->temps);
   ISIS_FREE (map->hdu);
   ISIS_FREE (map);
}

/*}}}*/

static EM_filemap_t *new_filemap (void) /*{{{*/
{
   EM_filemap_t *map;

   if (NULL == (map = (EM_filemap_t *) ISIS_MALLOC (sizeof(EM_filemap_t))))
     return NULL;

   map->densities = NULL;
   map->temps = NULL;
   map->hdu = NULL;
   map->num_hdus = 0;

   return map;
}

/*}}}*/

static int filter_filemap (EM_filemap_t *map, void *cl) /*{{{*/
{
   EM_Range_Type *r = (EM_Range_Type *)cl;
   int hdu, keep = 0;

   if (r == NULL) return 0;

   for (hdu = 0; hdu < map->num_hdus; hdu++)
     {
        if (((r->trange[0] <= map->temps[hdu])
             && (map->temps[hdu] <= r->trange[1]))
            && ((r->drange[0] <= map->densities[hdu])
                && (map->densities[hdu] <= r->drange[1])))
          {
             if (keep != hdu)
               {
                  map->hdu[keep] = map->hdu[hdu];
                  map->temps[keep] = map->temps[hdu];
                  map->densities[keep] = map->densities[hdu];
               }
             keep++;
          }
     }

   map->num_hdus = keep;
   if (map->num_densities == 1) map->num_temps = keep;
   if (map->num_temps == 1) map->num_densities = keep;

   return 0;
}

/*}}}*/

static EM_filemap_t *get_filemap (char *filename, void *cl) /*{{{*/
{
   EM_filemap_t *map = NULL;
   cfitsfile *fp = NULL;
   int hdu, num_temps, num_densities;
   int ret = -1;

   if (filename == NULL)
     return NULL;

   if (NULL == (fp = cfits_open_file_readonly (filename)))
     return NULL;

   if (NULL == (map = new_filemap ()))
     goto finish;

   isis_strcpy (map->filename, filename, CFLEN_FILENAME);

   if (-1 == cfits_read_int_keyword (&num_temps, "INUM_TEMP", fp)
       || -1 == cfits_read_int_keyword (&num_densities, "INUM_DENSITIES", fp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "%s", filename);
        goto finish;
     }

   *map->abund_table = 0;
   (void) cfits_read_string_keyword (map->abund_table, "SABUND_SOURCE", fp);

   if ((-1 == cfits_movabs_hdu (2, fp))
       || (-1 == cfits_test_keyword ("PARAMETERS", "EXTNAME", fp)))
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", filename);
        goto finish;
     }

   map->num_temps = num_temps;
   map->num_densities = num_densities;
   map->num_hdus = num_temps * num_densities;  /* initial guess */

   if (NULL == (map->temps = (float *) ISIS_MALLOC (map->num_hdus * sizeof(float)))
       || NULL == (map->densities = (float *) ISIS_MALLOC (map->num_hdus * sizeof(float)))
       || NULL == (map->hdu = (int *) ISIS_MALLOC (map->num_hdus * sizeof(int))))
     goto finish;

   if (-1 == cfits_read_float_col (map->temps, map->num_hdus, 1, "kT", fp)
       || -1 == cfits_read_float_col (map->densities, map->num_hdus, 1, "EDensity", fp))
     {
        isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "%s", filename);
        goto finish;
     }

   /* Convert keV -> Kelvin  and update num_hdus if necessary */
   for (hdu = 0; hdu < map->num_hdus; hdu++)
     {
        /* hardwired offset -- FIXME ? */
        enum {MAP_HDU_OFFSET=3};
        map->temps[hdu] *= (1000.0 * ERG_PER_EV / BOLTZ);
        map->hdu[hdu] = hdu + MAP_HDU_OFFSET;

        if ((map->temps[hdu] <= 0.0) || (map->densities[hdu] <= 0.0))
          {
             map->num_hdus = hdu + 1;
             break;
          }
     }

   if (-1 == filter_filemap (map, cl))
     goto finish;

   ret = 0;
   finish:

   (void) cfits_close_file (fp);
   if (ret) 
     {
        free_filemap (map);
        map = NULL;
     }

   return map;
}

/*}}}*/

int EM_get_filemap (EM_t *em, char *emis_file, void *cl, unsigned int *num_hdus, double **temp, double **dens) /*{{{*/
{
   EM_filemap_t *map;
   int i, malloced = 0;

   if ((emis_file != NULL) && (*emis_file != 0))
     {
        if (NULL == (map = get_filemap (emis_file, cl)))
          return -1;
        malloced = 1;
     }
   else
     {
        if (em == NULL || em->line_data == NULL)
          return -1;
        /* assuming line and continuum maps are the same */
        map = em->line_data->map;
     }

   *num_hdus = map->num_hdus;

   *temp = (double *) ISIS_MALLOC (map->num_hdus * sizeof(double));
   *dens = (double *) ISIS_MALLOC (map->num_hdus * sizeof(double));
   if (*temp == NULL || *dens == NULL)
     {
        ISIS_FREE (*temp);
        ISIS_FREE (*dens);
        free_filemap (map);
        return -1;
     }

   for (i = 0; i < map->num_hdus; i++)
     {
        (*temp)[i] = (double) map->temps[i];
        (*dens)[i] = (double) map->densities[i];
     }

   if (malloced) free_filemap (map);

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ line emissivity input */

void EM_free_line_emis_list (EM_line_emis_t *p) /*{{{*/
{
   if (p == NULL)
     return;

   ISIS_FREE (p->line);
   ISIS_FREE (p->lookup);
   ISIS_FREE (p->emissivity);
   ISIS_FREE (p);
}

/*}}}*/

static EM_line_emis_t *new_line_emis_list (int nlines) /*{{{*/
{
   EM_line_emis_t *emis = NULL;

   if (NULL == (emis = (EM_line_emis_t *) ISIS_MALLOC (sizeof(EM_line_emis_t))))
     return NULL;
   memset ((char *)emis, 0, sizeof (*emis));

   if (NULL == (emis->line = (DB_line_t **) ISIS_MALLOC (nlines * sizeof(DB_line_t *)))
       || NULL == (emis->emissivity = (float *) ISIS_MALLOC (nlines * sizeof(float))))
     {
        EM_free_line_emis_list (emis);
        return NULL;
     }

   emis->par[0] = emis->par[1] = 0.0;
   emis->nlines = nlines;
   emis->lookup = NULL;

   return emis;
}
/*}}}*/

static int read_temperature_keyword (float *temp, cfitsfile *fp) /*{{{*/
{
   if (temp == NULL || fp == NULL)
     return -1;

   if (-1 == cfits_read_float_keyword (temp, "TEMPERAT", fp)
       && -1 == cfits_read_float_keyword (temp, "TEMPERATURE", fp))
     return -1;
   else
     return 0;
}

/*}}}*/

typedef struct
{
   int (*init)(void ***, int);
   int (*process)(void **, int, int, Load_Linefile_Type *, DB_t *);
   int (*fini)(void **, int);
}
Linefile_Action_Type;

static int Unidentified_Lines;

static int apply_to_linefile_hdu (cfitsfile *fp, DB_t *db, void **user_data, /*{{{*/
                                  Linefile_Action_Type *a)
{
   Load_Linefile_Type x;
   int nlines, k, opt;
   int ret = -1;

   memset ((char *)&x, 0, sizeof (x));

   if (isis_user_break ())
     return -1;

   if (-1 == cfits_read_int_keyword (&nlines, "NAXIS2", fp))
     return -1;

   if (-1 == (*a->init)(&user_data, nlines))
     return -1;

   if (-1 == read_temperature_keyword (&x.temperature, fp)
       || -1 == cfits_read_float_keyword (&x.density, "DENSITY", fp))
     goto free_and_return;

   if (-1 == (opt = cfits_optimal_numrows (fp)))
     goto free_and_return;

   if (NULL == (x.lambda = (float *) ISIS_MALLOC (opt * sizeof(float)))
       || NULL == (x.epsilon = (float *) ISIS_MALLOC (opt * sizeof(float)))
       || NULL == (x.Z = (int *) ISIS_MALLOC (opt * sizeof(int)))
       || NULL == (x.rmJ = (int *) ISIS_MALLOC (opt * sizeof(int)))
       || NULL == (x.up = (int *) ISIS_MALLOC (opt * sizeof(int)))
       || NULL == (x.lo = (int *) ISIS_MALLOC (opt * sizeof(int))))
     goto free_and_return;

   k = 0;
   while (nlines - k*opt > 0)
     {
        int nread, start_row;

        start_row = k * opt;
        nread = MIN (opt, nlines - start_row);

        /* FIXME:  Look for "Plambda" column only for backward compatibility;
         *         this should go away soon.
         */

        if ((-1 == cfits_read_float_col (x.lambda, nread, 1+start_row, "Lambda", fp))
            || ((-1 == cfits_read_float_col (x.epsilon, nread, 1+start_row, "Epsilon", fp)
                 && -1 == cfits_read_float_col (x.epsilon, nread, 1+start_row, "Plambda", fp)))
            || (-1 == cfits_read_int_col (x.Z, nread, 1+start_row, "Element", fp))
            || (-1 == cfits_read_int_col (x.rmJ, nread, 1+start_row, "Ion", fp))
            || (-1 == cfits_read_int_col (x.up, nread, 1+start_row, "UpperLev", fp))
            || (-1 == cfits_read_int_col (x.lo, nread, 1+start_row, "LowerLev", fp)))
          goto free_and_return;

        if (-1 == (*a->process)(user_data, nread, start_row, &x, db))
          goto free_and_return;

        k++;
     }

   ret = 0;

   free_and_return:

   ISIS_FREE (x.lambda);
   ISIS_FREE (x.epsilon);
   ISIS_FREE (x.Z);
   ISIS_FREE (x.rmJ);
   ISIS_FREE (x.up);
   ISIS_FREE (x.lo);

   return (*a->fini)(user_data, ret);
}

/*}}}*/

static int init_load_linefile_hdu (void ***vp, int nlines) /*{{{*/
{
   EM_line_emis_t **p = *(EM_line_emis_t ***) vp;

   *p = new_line_emis_list (nlines);

   return ((NULL == *p) ? -1 : 0);
}

/*}}}*/

static int clean_load_linefile_hdu (void **vp, int ret) /*{{{*/
{
   if (ret)
     EM_free_line_emis_list (*(EM_line_emis_t **)vp);

   return ret;
}

/*}}}*/

static int do_load_linefile_hdu (void **vp, int nread, int start_row, /*{{{*/
                                 Load_Linefile_Type *x, DB_t *db)
{
   EM_line_emis_t *p = *(EM_line_emis_t **)vp;
   float *emis = p->emissivity;
   DB_line_t **line = p->line;
   int i;

   p->par[EM_TEMP] = x->temperature;
   p->par[EM_DENS] = x->density;

   for (i=0; i < nread; i++)
     {
        int o = i + start_row;

        /* ion charge from roman numeral */
        int q = x->rmJ[i] - 1;

        /* Emissivities already in photons cm^3 s^-1 */
        emis[o] = x->epsilon[i];
        line[o] = DB_get_line (x->lambda[i], x->Z[i], q, x->up[i], x->lo[i], db);

        if (line[o])
          line[o]->have_emissivity_data = 1;
        else
          Unidentified_Lines++;
     }

   return 0;
}

/*}}}*/

static void free_line_data (EM_line_data_t *ld) /*{{{*/
{
   if (NULL == ld)
     return;

   if (ld->emis)
     {
        int i;
        for(i = 0; i < ld->map->num_hdus; i++)
          {
             if (ld->emis[i] != NULL)
               EM_free_line_emis_list (ld->emis[i]);
          }
        ISIS_FREE (ld->emis);
     }

   free_filemap (ld->map);
   ISIS_FREE (ld);
}

/*}}}*/

static EM_line_data_t *new_line_data (EM_filemap_t *map) /*{{{*/
{
   EM_line_data_t *ld = NULL;

   if (map == NULL)
     return NULL;

   if (NULL == (ld = (EM_line_data_t *) ISIS_MALLOC (sizeof(EM_line_data_t))))
     return NULL;
   memset ((char *)ld, 0, sizeof (*ld));

   ld->map = map;

   if (EM_Load_Line_Emis)
     {
        ld->emis = (EM_line_emis_t **) ISIS_MALLOC (map->num_hdus * sizeof(EM_line_emis_t *));
        if (ld->emis == NULL)
          {
             ISIS_FREE (ld);
             return NULL;
          }
        memset ((char *)ld->emis, 0, map->num_hdus * sizeof(EM_line_emis_t *));
     }

   return ld;
}

/*}}}*/

static int update_line_atomic_data (EM_filemap_t *map, DB_t *db) /*{{{*/
{
   Linefile_Action_Type hash =
     {
        &start_hashing, &do_hashing, &finish_hashing
     };

   FILE *progress = stderr;
   cfitsfile *fp = NULL;
   int j, ret = -1;
   static const char *fmt = "hdu:  %d/%d\r";

   if (db == NULL || map == NULL)
     return -1;

   if (NULL == (fp = cfits_open_file_readonly (map->filename)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", map->filename);
        goto finish;
     }

   if (-1 == init_hash_table (EM_Hash_Table_Size_Hint))
     goto finish;

   for (j=0; j < map->num_hdus; j++)
     {
        int hdu = map->hdu[j];

        if (-1 == cfits_movabs_hdu (hdu, fp))
          {
             isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__, "hdu=%d, %s", hdu, map->filename);
             goto finish;
          }

        if (Isis_Verbose >= WARN)
          fprintf (progress, fmt, j+1, map->num_hdus);

        if (-1 == apply_to_linefile_hdu (fp, db, NULL, &hash))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "building line list");
             goto finish;
          }
     }
   if (Isis_Verbose >= WARN)
     fputc ('\n', progress);

   if (-1 == merge_hash_table (db))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "updating line list");
        goto finish;
     }

   ret = 0;
   finish:

   deallocate_hash_table ();
   (void) cfits_close_file (fp);

   return ret;
}

/*}}}*/

static int load_line_spectrum_hdu (cfitsfile *fp, int hdu, EM_t *em, EM_line_emis_t **x) /*{{{*/
{
   Linefile_Action_Type load =
     {
        &init_load_linefile_hdu,
        &do_load_linefile_hdu,
        &clean_load_linefile_hdu
     };
   EM_filemap_t *map;

   if (fp == NULL || em == NULL
       || em->db == NULL || em->line_data == NULL)
     return -1;

   map = em->line_data->map;

   if (-1 == cfits_movabs_hdu (hdu, fp))
     return -1;

   Unidentified_Lines = 0;
   if (-1 == apply_to_linefile_hdu (fp, em->db, (void **) x, &load))
     return -1;

   if (Unidentified_Lines > 0)
     isis_vmesg (WARN, I_ERROR, __FILE__, __LINE__, "%d unidentified lines!",
                 Unidentified_Lines);
   return 0;
}

/*}}}*/

static int load_line_emissivity_data (EM_line_data_t *ld, DB_t *db) /*{{{*/
{
   Linefile_Action_Type load =
     {
        &init_load_linefile_hdu,
        &do_load_linefile_hdu,
        &clean_load_linefile_hdu
     };

   FILE *progress = stderr;
   EM_filemap_t *map = ld->map;
   cfitsfile *fp = NULL;
   int j, hdu, ret = -1;
   static const char *fmt = "hdu:  %d/%d\r";

   if (NULL == (fp = cfits_open_file_readonly (map->filename)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", map->filename);
        goto finish;
     }

   isis_vmesg (WARN, I_LOADING, __FILE__, __LINE__, "line emissivity tables [%d hdu%s]",
               map->num_hdus,
               (map->num_hdus > 1) ? "s" : "");

   for (j = 0; j < map->num_hdus; j++)
     {
        hdu = map->hdu[j];

        if (-1 == cfits_movabs_hdu (hdu, fp))
          {
             isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__, "hdu=%d, %s", hdu, map->filename);
             goto finish;
          }

        if (Isis_Verbose >= WARN)
          fprintf (progress, fmt, j+1, map->num_hdus);

        Unidentified_Lines = 0;
        if (-1 == apply_to_linefile_hdu (fp, db, (void **) &ld->emis[j], &load))
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", map->filename);
             goto finish;
          }
        if (Unidentified_Lines > 0)
          isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, 
                      "%d unidentified lines in extension %d",
                      Unidentified_Lines, hdu);
     }

   if (Isis_Verbose >= WARN)
     fputc ('\n', progress);
   isis_vmesg (INFO, I_READ_OK, __FILE__, __LINE__, "%s", map->filename);

   ret = 0;

   finish:
   (void) cfits_close_file (fp);

   return ret;
}

/*}}}*/

static EM_line_data_t *load_line_data (char *filename, void *cl, DB_t *db) /*{{{*/
{
   EM_line_data_t *ld = NULL;
   EM_filemap_t *map;

   if (NULL == db)
     return NULL;

   map = get_filemap (filename, cl);
   if (NULL == map)
     return NULL;

   ld = new_line_data (map);
   if (ld == NULL)
     return NULL;

   /* Because some lines in the emissivity tables
    * aren't listed in the atomic data tables (!) we have to
    *   1) scan all the emissivity tables to assemble a complete
    *      line list (with no duplicates)
    *   2) update our internal run-time list.
    */

   if (EM_Maybe_Missing_Lines)
     {
        isis_vmesg (WARN, I_SCANNING, __FILE__, __LINE__, 
                    "line emissivity tables [%d hdu%s]",
                    ld->map->num_hdus,
                    (ld->map->num_hdus > 1) ? "s" : "");

        if (-1 == update_line_atomic_data (ld->map, db))
          {
             free_line_data (ld);
             return NULL;
          }
     }

   if (EM_Load_Line_Emis)
     {
        if (-1 == load_line_emissivity_data (ld, db))
          {
             free_line_data (ld);
             return NULL;
          }
     }

   return ld;
}

/*}}}*/

/*}}}*/

/*{{{ continuum emissivity input */

static void free_cont_type (EM_cont_emis_t *p) /*{{{*/
{
   if (p == NULL)
     return;

   ISIS_FREE (p->g_true_contin);
   ISIS_FREE (p->true_contin);
   ISIS_FREE (p->g_pseudo);
   ISIS_FREE (p->pseudo);
   ISIS_FREE (p);
}

/*}}}*/

static EM_cont_emis_t *new_cont_type (int ntrue_contin, int npseudo) /*{{{*/
{
   EM_cont_emis_t *p;

   if (ntrue_contin <= 0 || npseudo <= 0)
     return NULL;

   if (NULL == (p = (EM_cont_emis_t *) ISIS_MALLOC (sizeof(EM_cont_emis_t))))
     return NULL;
   memset ((char *)p, 0, sizeof (*p));

   if (NULL == (p->g_true_contin = (double *) ISIS_MALLOC (ntrue_contin * sizeof(double)))
       || NULL == (p->true_contin = (double *) ISIS_MALLOC (ntrue_contin * sizeof(double)))
       || NULL == (p->g_pseudo = (double *) ISIS_MALLOC (npseudo * sizeof(double)))
       || NULL == (p->pseudo = (double *) ISIS_MALLOC (npseudo * sizeof(double))))
     {
        free_cont_type (p);
        return NULL;
     }

   p->ntrue_contin = ntrue_contin;
   p->npseudo = npseudo;

   p->Z = -1;
   p->q = -1;

   p->next = NULL;

   return p;
}

/*}}}*/

static void free_cont_list (EM_cont_emis_t *head) /*{{{*/
{
   EM_cont_emis_t *new_head;

   if (head == NULL)
     return;

   while (head != NULL)
     {
        new_head = head->next;
        free_cont_type (head);
        head = new_head;
     }
}

/*}}}*/

static EM_cont_emis_t *find_cont_type (EM_cont_emis_t *head, int Z, int q) /*{{{*/
{
   EM_cont_emis_t *h;

   if (head == NULL)
     return NULL;

   for (h=head; h != NULL; h = h->next)
     {
        if (h->Z == Z && h->q == q)
          return h;
     }

   return NULL;
}

/*}}}*/

void EM_free_continuum (EM_cont_type_t *p) /*{{{*/
{
   if (p == NULL)
     return;

   ISIS_FREE (p->wlhi);
   ISIS_FREE (p->wllo);
   ISIS_FREE (p->pseudo);
   ISIS_FREE (p->true_contin);
   ISIS_FREE (p);
}

/*}}}*/

EM_cont_type_t *EM_new_continuum (int nbins) /*{{{*/
{
   EM_cont_type_t *p = NULL;

   if (NULL == (p = (EM_cont_type_t *) ISIS_MALLOC (sizeof(EM_cont_type_t))))
     return NULL;
   memset ((char *)p, 0, sizeof(*p));

   if (NULL == (p->wlhi = (double *) ISIS_MALLOC (nbins * sizeof(double)))
       || NULL == (p->wllo = (double *) ISIS_MALLOC (nbins * sizeof(double)))
       || NULL == (p->pseudo = (double *) ISIS_MALLOC (nbins * sizeof(double)))
       || NULL == (p->true_contin = (double *) ISIS_MALLOC (nbins * sizeof(double))))
     {
        EM_free_continuum (p);
        return NULL;
     }

   p->nbins = nbins;

   return p;
}

/*}}}*/

static int reverse_dbl (double *x, int n) /*{{{*/
{
   int i;

   if (NULL == x || n <= 0)
     return -1;

   for (i=0; i < n/2; i++)
     {
        int j = n - i - 1;
        double tmp = x[j];
        x[j] = x[i];
        x[i] = tmp;
     }

   return 0;
}

/*}}}*/

static int cvt_energy_to_wavelength (double *grid, double *val, int *nbins) /*{{{*/
{
   int i, n;

   /* count non-zero energy values */

   n = 0;
   for (i = 0; i < *nbins; i++)
     {
        if (grid[i] > 0.0)
          n++;
        else
          break;
     }

   *nbins = n;

   /* convert grid to wavelength units
    * and bin value from photons / keV to photons / Angstrom
    * and make sure there are no negative bin values
    * since that's unphysical.
    */

   if (-1 == reverse_dbl (grid, n)
       || -1 == reverse_dbl (val, n))
     return -1;

   while (n-- > 0)
     {
        if (grid[n] <= 0.0)
          return -1;

        grid[n] = KEV_ANGSTROM / grid[n];
        val[n]  = fabs(val[n]) * (KEV_ANGSTROM / grid[n] / grid[n]);
     }

   return 0;
}

/*}}}*/

static int make_canonical_continuum (EM_cont_emis_t *p) /*{{{*/
{
   if (NULL == p)
     return -1;

   if (-1 == cvt_energy_to_wavelength (p->g_true_contin, p->true_contin, &p->ntrue_contin))
     p->ntrue_contin = 0;
   if (-1 == cvt_energy_to_wavelength (p->g_pseudo, p->pseudo, &p->npseudo))
     p->npseudo = 0;

   return 0;
}

/*}}}*/

static void free_cont_data (EM_cont_data_t *cd) /*{{{*/
{
   if (NULL == cd)
     return;

   if (cd->emis)
     {
        int i;
        for(i = 0; i < cd->map->num_hdus; i++)
          {
             if (cd->emis[i] != NULL)
               free_cont_list (cd->emis[i]);
          }
        ISIS_FREE (cd->emis);
     }

   free_filemap (cd->map);
   ISIS_FREE (cd);
}

/*}}}*/

static EM_cont_data_t *new_cont_data (EM_filemap_t *map) /*{{{*/
{
   EM_cont_data_t *cd = NULL;

   if (map == NULL)
     return NULL;

   if (NULL == (cd = (EM_cont_data_t *) ISIS_MALLOC (sizeof(EM_cont_data_t))))
     return NULL;
   memset ((char *)cd, 0, sizeof (*cd));

   cd->map = map;

   if (EM_Load_Cont_Emis)
     {
        cd->emis = (EM_cont_emis_t **) ISIS_MALLOC (map->num_hdus * sizeof(EM_cont_emis_t *));
        if (cd->emis == NULL)
          {
             ISIS_FREE (cd);
             return NULL;
          }
        memset ((char *)cd->emis, 0, map->num_hdus * sizeof(*cd->emis));
     }

   return cd;
}
/*}}}*/

static EM_cont_emis_t *load_cont_emissivity_hdu (cfitsfile *fp, int Z_req, int q_req, int load_all, EM_t *em) /*{{{*/
{
   EM_cont_emis_t *head;
   EM_cont_emis_t *last;
   int nrows;
   float temp, dens;
   int k, foundit, ret = -1;
   int alt_ioniz = use_alt_ioniz (em);
   int alt_abund = use_alt_abund (em);

   if (NULL == fp || NULL == em)
     return NULL;

   head = last = NULL;

   if (-1 == cfits_read_int_keyword (&nrows, "NAXIS2", fp)
       || -1 == read_temperature_keyword (&temp, fp)
       || -1 == cfits_read_float_keyword (&dens, "DENSITY", fp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "continuum emissivity file");
        goto free_and_return;
     }

   foundit = (Z_req > 0) ? 0 : 1;

   for (k = 1; k <= nrows; k++)
     {
        EM_cont_emis_t *p = NULL;
        int Z, rmj, ntrue_contin, npseudo;

        if (-1 == cfits_read_int_col (&Z, 1, k, "Z", fp)
            || -1 == cfits_read_int_col (&rmj, 1, k, "rmJ", fp)
            || -1 == cfits_read_int_col (&ntrue_contin, 1, k, "N_Cont", fp)
            || -1 == cfits_read_int_col (&npseudo, 1, k, "N_Pseudo", fp))
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "continuum emissivity file");
             goto free_and_return;
          }

        if (NULL == (p = new_cont_type (ntrue_contin, npseudo)))
          goto free_and_return;

        /* e.g. for Fe X    rmJ = 10, charge q = 9 */
        p->Z = Z;
        p->q = rmj - 1;
        p->temp = (double) temp;
        p->dens = (double) dens;

        if (-1 == cfits_read_double_col (p->g_true_contin, p->ntrue_contin, k, "E_Cont", fp)
            || -1 == cfits_read_double_col (p->true_contin, p->ntrue_contin, k, "Continuum", fp)
            || -1 == cfits_read_double_col (p->g_pseudo, p->npseudo, k, "E_Pseudo", fp)
            || -1 == cfits_read_double_col (p->pseudo, p->npseudo, k, "Pseudo", fp)
            )
          {
             isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "continuum emissivity file");
             free_cont_type (p);
             goto free_and_return;
          }

        if (-1 == make_canonical_continuum (p))
          {
             free_cont_type (p);
             goto free_and_return;
          }

        if (k == 1)
          head = p;
        else
          last->next = p;         /* append */

        last = p;

        /* try to bail out early */

        if (!(alt_ioniz || alt_abund || load_all)
            && (Z_req == p->Z && q_req == p->q))
          {
             foundit = 1;
             break;
          }
     }

   if (!foundit)
     {
        isis_vmesg (FAIL, I_NOT_FOUND, __FILE__, __LINE__, "Z=%d, rmJ=%d continuum",
                      Z_req, q_req+1);
        ret = -1;
     }
   else ret = 0;

   free_and_return:

   if (ret)
     {
        free_cont_list (head);
        return NULL;
     }

   return head;                   /* NULL if failure return */
}

/*}}}*/

static int load_cont_emissivity_data (EM_cont_data_t *cd, EM_t *em) /*{{{*/
{
   FILE *progress = stderr;
   EM_filemap_t *map = cd->map;
   cfitsfile *fp = NULL;
   int j, hdu, ret = -1;
   static const char *fmt = "hdu:  %d/%d\r";

   if (NULL == (fp = cfits_open_file_readonly (map->filename)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", map->filename);
        goto finish;
     }

   isis_vmesg (WARN, I_LOADING, __FILE__, __LINE__, "continuum emissivity tables [%d hdu%s]",
               map->num_hdus,
               (map->num_hdus > 1) ? "s" : "");

   for (j = 0; j < map->num_hdus; j++)
     {
        hdu = map->hdu[j];

        if (-1 == cfits_movabs_hdu (hdu, fp))
          {
             isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__, "hdu=%d, %s", hdu, map->filename);
             goto finish;
          }

        if (Isis_Verbose >= WARN)
          fprintf (progress, fmt, j+1, map->num_hdus);

        /* Z_req = 0, q_req = 0, load_all = 1 */
        cd->emis[j] = load_cont_emissivity_hdu (fp, 0, 0, 1, em);
        if (cd->emis[j] == NULL)
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "Reading hdu=%d, %s", hdu, map->filename);
             goto finish;
          }
     }

   if (Isis_Verbose >= WARN)
     fputc ('\n', progress);
   isis_vmesg (INFO, I_READ_OK, __FILE__, __LINE__, "%s", map->filename);

   ret = 0;

   finish:
   (void) cfits_close_file (fp);

   return ret;
}

/*}}}*/

static EM_cont_data_t *load_cont_data (char *filename, void *cl, EM_t *em) /*{{{*/
{
   EM_cont_data_t *cd = NULL;
   EM_filemap_t *map;

   if (NULL == filename)
     return NULL;

   if (NULL == (map = get_filemap (filename, cl)))
     return NULL;

   if (NULL == (cd = new_cont_data (map)))
     return NULL;

   if (EM_Load_Cont_Emis)
     {
        if (-1 == load_cont_emissivity_data (cd, em))
          {
             free_cont_data (cd);
             return NULL;
          }
     }

   return cd;
}

/*}}}*/

/*}}}*/

/*{{{ table search and interpolate */

/* this search assumes t[i] < t[i+1]  */

static int find_in_table (float x, float *t, int n) /*{{{*/
{
   int i;

   if (NULL == t || n <= 0)
     return -1;

   if (x < t[0] || t[n-1] < x)
     {
        isis_vmesg (FAIL, I_RANGE_ERROR, __FILE__, __LINE__, 
                    "%g is outside [%11.4e, %11.4e]",
                    x, t[0], t[n-1]);
        return -1;
     }

   for (i=0; i < n-1; i++)                /* brute force for now */
     {
        if (t[i] <= x && x < t[i+1])
          return i;
     }

   return -1;                             /* should never happen */
}

/*}}}*/

static int linear_interp (float f[/*2*/], int ix[/*2*/], float x, float *t, int n) /*{{{*/
{
   int k;

   if (NULL == t || n <= 0 || NULL == f)
     return -1;

   if (-1 == (k = find_in_table (x, t, n)))
     return -1;

   if ( k == n-1
        || (k == 0 && x < t[0]) )
     {
        f[0] =  f[1] = 0.5;
        ix[0] = ix[1] = k;
     }
   else
     {
        float p = (x - t[k]) / (t[k+1] - t[k]);
        f[0] = 1 - p;      f[1] = p;
        ix[0] = k;         ix[1] = k+1;
     }

   return 0;
}

/*}}}*/

static int bilinear_interp (float c[/*4*/], int ip[/*4*/], /*{{{*/
                            float xp, float yp,
                            float *x, float *y, int n)
{
   float dist[4], fx, fy;
   float tol = 0.0; /* 100.0 * FLT_EPSILON; */
   int i;

   if ((NULL == y) || (NULL == x) || (n <= 0)
       || xp <= 0.0 || yp <= 0.0)
     return -1;

   /* brute force  -- use logs to minimize scale problems
    * ip[1] = x-min, y-max     ip[3] = x-max, y-max
    * ip[0] = x-min, y-min     ip[2] = x-max, y-min
    */

   ip[0] = ip[1] = ip[2] = ip[3] = -1;
   dist[0] = dist[1] = dist[2] = dist[3] = FLT_MAX;

   for (i = 0; i < n; i++)
     {
        float dx, dy, r2;
        unsigned int ii;

        dy = log(yp / y[i]);
        dx = log(xp / x[i]);

        ii = 0;
        if (dy >= tol) ii += 1;
        if (dx >= tol) ii += 2;

        r2 = dx*dx + dy*dy;
        if (r2 < dist[ii])
          {
             dist[ii] = r2;
             ip[ii] = i;
          }
     }

   /* if any corner is unassigned, something failed */
   if (ip[0] < 0 || ip[1] < 0 || ip[2] < 0 || ip[3] < 0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "values out of range:  (%g, %g)", xp, yp);
        return -1;
     }

   if (x[ip[0]] != x[ip[2]])
     fx = (x[ip[2]] - xp) / (x[ip[2]] - x[ip[0]]);
   else fx = 0.5;

   if (y[ip[0]] != y[ip[1]])
     fy = (y[ip[1]] - yp) / (y[ip[1]] - y[ip[0]]);
   else fy = 0.5;

   if ((fx < 0.0) || (1.0 < fx) || (fy < 0.0) || (1.0 < fy))
     {
        isis_vmesg (FAIL, I_RANGE_ERROR, __FILE__, __LINE__, "[%g, %g] is out of range",
                    xp, yp);
        return -1;
     }

   c[0] =        fx  *        fy;
   c[1] =        fx  * (1.0 - fy);
   c[2] = (1.0 - fx) *        fy;
   c[3] = (1.0 - fx) * (1.0 - fy);

   return 0;
}

/*}}}*/

static int interp_coeffs (float *coef, int *idx, int *npoints, /*{{{*/
                          float *par, EM_filemap_t *map)
{
   float temp = par [EM_TEMP];
   float density = par [EM_DENS];

   if (map->num_densities == 1)
     {
        *npoints = 2;
        return linear_interp (coef, idx, temp, map->temps, map->num_temps);
     }
   else if (map->num_temps == 1)
     {
        *npoints = 2;
        return linear_interp (coef, idx, density, map->densities, map->num_densities);
     }
   else
     {
        *npoints = 4;
        return bilinear_interp (coef, idx, temp, density,
                                map->temps, map->densities, map->num_hdus);
     }
}

/*}}}*/

/*}}}*/

/*{{{ interpolate line spectrum using filemap */

static int build_lookup_table (EM_line_emis_t *p, DB_t *db) /*{{{*/
{
   DB_line_t **line;
   int i, nlines;
   int *lookup = NULL;

   if (NULL == p || NULL == db)
     return -1;

   if (-1 == (nlines = DB_get_nlines (db))
       || NULL == (lookup = (int *) ISIS_MALLOC (nlines * sizeof(int))))
     return -1;

   for (i=0; i < nlines; i++)        /* init to invalid values */
     lookup[i] = INT_MIN;

   line = p->line;

   for (i=0; i < p->nlines; i++)
     {
        if (line[i])
          lookup[line[i]->indx] = i;
     }

   p->lookup = lookup;

   return 0;
}

/*}}}*/

static int get_line_interp_points (EM_t *em, int npoints, int *idx, EM_line_emis_t **tbl) /*{{{*/
{
   EM_line_data_t *ld;
   EM_filemap_t *map;
   cfitsfile *fp = NULL;
   int j;

   ld  = em->line_data;
   map = ld->map;

   /* DB in memory */
   if (EM_Load_Line_Emis)
     {
        for (j=0; j < npoints; j++)
          {
             tbl[j] = ld->emis[ idx[j] ];
          }
        return 0;
     }

   /* DB on disk */
   if (NULL == (fp = cfits_open_file_readonly (map->filename)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", map->filename);
        return -1;
     }

   for (j = 0; j < npoints; j++)
     {
        int hdu = map->hdu[ idx[j] ];
        if (-1 == load_line_spectrum_hdu (fp, hdu, em, &tbl[j]))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading %s[%d]", map->filename, hdu);
             (void) cfits_close_file (fp);
             return -1;
          }
     }

   (void) cfits_close_file (fp);

   return 0;
}

/*}}}*/

static EM_line_emis_t *interpolate_line_emis (char *flag, EM_line_emis_t **table, /*{{{*/
                                             float *coef, int n, DB_t *db)
{
   EM_line_emis_t *result = NULL;
   int j, k, nl, ntot;
   int  ret = -1;

   struct line_interp_t
     {
        float emis;
        DB_line_t *line;
     }
   *t = NULL;

   if (NULL == table || NULL == coef
       || (n != 2 && n != 4))
     return NULL;

   /* Use a temporary array of structures to
    *  1) determine the union of the N line lists and
    *  2) to assign an interpolated emissivity to each line.
    *
    * The temp array is used as a lookup table in which the array
    * index is equal to the index of the line in the internal list.
    * (e.g. it is big enough to hold the entire available line list).
    *
    * We waste a bit of memory, but this is a simple algorithm and
    * is guaranteed to work (right?) even if the input line lists
    * are unsorted.  and memory is cheap.... isnt it?
    *
    */

   if (-1 == (ntot = DB_get_nlines (db)))
     return NULL;
   if (NULL == (t = (struct line_interp_t *) ISIS_MALLOC (ntot * sizeof(struct line_interp_t))))
     return NULL;
   memset ((char *)t, 0, ntot * sizeof(struct line_interp_t));

   nl = 0;
   for (j=0; j < n; j++)
     {
        EM_line_emis_t *tbl = table[j];
        int nlines = tbl->nlines;

        for (k=0; k < nlines; k++)
          {
             DB_line_t *p = tbl->line[k];
             int idx;

             if (p == NULL)
               goto fail;

             idx = p->indx;

             if (t[idx].line == NULL)
               {
                  t[idx].line = p;
                  nl++;
               }
             else if (t[idx].line != p)       /* safety check */
               goto fail;

             if (NULL == flag || flag[idx] != 0)
               t[idx].emis += coef[j] * tbl->emissivity[k];
          }
     }

   if (nl == 0
       || NULL == (result = new_line_emis_list (nl)))
     goto fail;

   k = 0;
   for (j=0; j < ntot; j++)
     {
        if (t[j].line == NULL)
          continue;

        result->line[ k ] = t[j].line;
        result->emissivity[ k ] = MAX(t[j].emis, 0.0);
        k++;
     }

   if (k != nl)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Internal error in interpolate_line_emis");
        goto fail;
     }

   if (-1 == build_lookup_table (result, db))
     goto fail;

   ret = 0;

   fail:

   if (ret)
     {
        EM_free_line_emis_list (result);
        result = NULL;
     }

   ISIS_FREE (t);
   return result;
}
/*}}}*/

EM_line_emis_t *EM_get_line_spectrum (char *flag, float *par, EM_t *em) /*{{{*/
{
   EM_line_data_t *ld;
   EM_line_emis_t *line = NULL;
   EM_line_emis_t *tbl[4] = {NULL, NULL, NULL, NULL};
   float coef[4] = {0.0, 0.0, 0.0, 0.0};
   int idx[4], npoints;

   if (NULL == em || NULL == em->line_data)
     return NULL;

   ld  = em->line_data;

   if (-1 == interp_coeffs (coef, idx, &npoints, par, ld->map))
     return NULL;

   if (-1 == get_line_interp_points (em, npoints, idx, tbl))
     goto close_and_return;

   line = interpolate_line_emis (flag, tbl, coef, npoints, em->db);
   if (NULL == line)
     goto close_and_return;

   if (-1 == scale_line_abundance (line, em)
       || -1 == scale_line_ionization (line, par, em))
     goto close_and_return;

   line->par[EM_TEMP] = par[EM_TEMP];
   line->par[EM_DENS] = par[EM_DENS];

   close_and_return:

   if (EM_Load_Line_Emis == 0)
     {
        int j;
        for (j=0; j < npoints; j++)
          {
             EM_free_line_emis_list (tbl[j]);
          }
     }

   return line;
}
/*}}}*/

int EM_get_nlines (EM_line_emis_t *t) /*{{{*/
{
   if (NULL == t)
     return -1;
   else
     return t->nlines;
}

/*}}}*/

int _EM_get_line_emis_wl (DB_line_t **line, float *emis, float *wl, int k, /*{{{*/
                         EM_line_emis_t *t)
{
   DB_line_t *pline;

   if ((NULL == t) || (k < 0) || (k > t->nlines))
     return -1;

   pline = t->line[k];
   if (pline == NULL)
     return -1;

   *line = pline;
   *wl = pline->wavelen;
   *emis = t->emissivity[k];

   return 0;
}

/*}}}*/

int EM_sum_line_emissivity (float * emis, float *par, /*{{{*/
                            int *list, int nlines, EM_t *em)
{
   EM_line_emis_t *t = NULL;
   int i, k, max_index;
   int ret = -1;

   if (em == NULL)
     return -1;

   /* also applies abundance and ion-fraction factors */

   if (NULL == (t = EM_get_line_spectrum (NULL, par, em)))
     return -1;

   max_index = DB_get_nlines (em->db);

   *emis = 0.0;

   for (i=0; i < nlines; i++)
     {
        if (list[i] < 0 || list[i] >= max_index)
          {
             *emis = -1.0;
             goto finish;
          }
        k = t->lookup[ list[i] ];
        if (k < 0)
          continue;
        *emis += t->emissivity[k];
     }

   ret = 0;

   finish:
   EM_free_line_emis_list (t);

   return ret;
}
/*}}}*/

int EM_get_line_emissivity_function (float **emis, float **temps, float **densities, /*{{{*/
                                     int *num_points,
                                     int line_index, EM_t *em)
{
   EM_line_data_t *ld;
   EM_filemap_t *map;
   cfitsfile *fp = NULL;
   int i, not_found = 1;;

   if (NULL == em)
     return -1;

   if (NULL == DB_get_line_from_index (line_index, em->db))
     {
        isis_vmesg (FAIL, I_NOT_FOUND, __FILE__, __LINE__, "line %d", line_index);
        return -1;
     }

   ld = em->line_data;
   map = ld->map;

   *num_points = map->num_hdus;

   if ((NULL == (*emis = (float *) ISIS_MALLOC (*num_points * sizeof(float))))
       || (NULL == (*temps = (float *) ISIS_MALLOC (*num_points * sizeof(float))))
       || (NULL == (*densities = (float *) ISIS_MALLOC (*num_points * sizeof(float)))))
     {
        ISIS_FREE (*emis);
        ISIS_FREE (*temps);
        ISIS_FREE (*densities);
        return -1;
     }

   if (EM_Load_Line_Emis == 0)
     {
        if (NULL == (fp = cfits_open_file_readonly (map->filename)))
          {
             isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", map->filename);
             ISIS_FREE (*emis);
             ISIS_FREE (*temps);
             ISIS_FREE (*densities);
             return -1;
          }
     }

   for (i = 0; i < *num_points; i++)
     {
        EM_line_emis_t *p;
        int k, found;
        DB_line_t **line;

        if (EM_Load_Line_Emis)
          p = ld->emis[i];
        else
          {
             int hdu = map->hdu[i];
             if (-1 == load_line_spectrum_hdu (fp, hdu, em, &p))
               {
                  isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading %s[%d]", map->filename, hdu);
                  ISIS_FREE (*emis);
                  ISIS_FREE (*temps);
                  ISIS_FREE (*densities);
                  (void) cfits_close_file (fp);
                  return -1;
               }
          }

        line = p->line;
        found = -1;

        for (k = 0; k < p->nlines ; k++)
          {
             if ((line[k] != NULL)
                 && (line_index == line[k]->indx))
               {
                  found = k;
                  not_found = 0;
                  break;
               }
          }

        (*temps)[i] = p->par[EM_TEMP];
        (*densities)[i] = p->par[EM_DENS];
        (*emis)[i] = (found < 0) ? 0.0 : p->emissivity[found];
     }

   if (fp)
     (void) cfits_close_file (fp);

   if (not_found)
     isis_vmesg (WARN, I_NO_DATA, __FILE__, __LINE__, "line %d emissivity", line_index);

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ interpolate continuum spectrum using filemap */

#define TEST_CONTIN 0

#if TEST_CONTIN /*{{{*/
static int test_contin (EM_cont_type_t *r)
{
   double tot[] = { 0.0, 0.0 };
   double e;
   int k;

   if (NULL == r)
     return -1;

   for (k = 0; k < r->nbins; k++)
     {
        double wl = 0.5 * (r->wllo[k] + r->wlhi[k]);
        e = ERG_ANGSTROM / wl;
        tot[0] += e * r->pseudo[k];
        tot[1] += e * r->true_contin[k];
     }

   fprintf (stderr, "pseudo = %e  total = %e\n",
            sum, tot[0], tot[1]);

   return 0;
}
#endif
/*}}}*/

static int bracket (int * idx, double *x, int n, double *t, int nt) /*{{{*/
{
   int k, j;

   for (j=0; j < n && x[j] < t[0]; j++)
     idx[j] = -1;

   for (k = 1; k < nt; k++)
     {
        while ((j < n) && (t[k-1] <= x[j]) && (x[j] < t[k]))
          idx[j++] = k-1;
     }

   while (j < n)
     idx[j++] = -1;

   return 0;
}

/*}}}*/

/* integrate from: y(x)
 *             to: (lo, hi, value)
 * where 'value' = either the integral or the bin-average
 */
static int bin_xy (double *y, double *x, int npts, int want_avg, /*{{{*/
                   double *value, double *lo, double *hi, int nbins)
{
   double fac, xm;
   int *k_lo, *k_hi;
   int i, ret = -1;

   k_lo = k_hi = NULL;

   if (NULL == (k_lo = (int *) ISIS_MALLOC (nbins * sizeof(int)))
       || NULL == (k_hi = (int *) ISIS_MALLOC (nbins * sizeof(int))))
     goto finish;

   /*    AREA(a,b,i) computes
    *       F(a,b) = Integral [ f(x), a, b ] when  x[i] <= a < b < x[i+1]
    *       and f(x) is given in tabular form f_i = f(x[i])
    */

#define   AREA(a,b,i)    (fac = ((b) - (a)) / (x[(i)+1] - x[(i)]),   \
                          xm = 0.5 * ((a) + (b)),                    \
                          fac * (y[ (i) ] * (x[(i)+1] - xm)          \
                               + y[(i)+1] * (xm - x[(i)])));

   /*  Bracket bin edges:
    *      x[klo]  <=  lo  <  x[klo+1]
    *      x[khi]  <=  hi  <  x[khi+1]
    */

   if (-1 == bracket (k_lo, lo, nbins, x, npts)
       || -1 == bracket (k_hi, hi, nbins, x, npts))
     goto finish;

   for (i=0; i < nbins; i++)
     {
        int klo = k_lo[i];
        int khi = k_hi[i];

        if (klo < 0 || klo >= npts-1
            || khi < 0 || khi >= npts-1)
          {
             value[i] = 0.0;
          }
        else if (khi == klo)
          {
             value[i] = AREA(lo[i], hi[i], klo);
          }
        else
          {
             double sum;
             int j;

             sum  = AREA( lo[i], x[klo+1], klo);
             sum += AREA(x[khi],    hi[i], khi);

             for (j = klo+1; j < khi; j++)
               {
                  sum += AREA(x[j], x[j+1], j);
               }

             value[i] = sum;
          }

        if (want_avg)
          value[i] /= hi[i] - lo[i];      /* assume grid has been validated */
     }

   ret = 0;
   finish:

   ISIS_FREE (k_lo);
   ISIS_FREE (k_hi);

   return ret;
}

/*}}}*/

static int add_cont_contrib (EM_cont_type_t *r, EM_cont_emis_t *t, float weight) /*{{{*/
{
   double * tmp = NULL;
   double wt = (double) weight;
   int k, ret = -1;

   if (NULL == (tmp = (double *) ISIS_MALLOC (r->nbins * sizeof(double))))
     goto finish;

   if (t->ntrue_contin > 0)
     {
        if (-1 == bin_xy (t->true_contin, t->g_true_contin, t->ntrue_contin, 0,
                          tmp, r->wllo, r->wlhi, r->nbins))
          goto finish;
        for (k = 0; k < r->nbins; k++)
          r->true_contin[k] += wt * tmp[k];
     }

   if (t->npseudo > 0)
     {
        if (-1 == bin_xy (t->pseudo, t->g_pseudo, t->npseudo, 0,
                          tmp, r->wllo, r->wlhi, r->nbins))
          goto finish;
        for (k = 0; k < r->nbins; k++)
          r->pseudo[k] += wt * tmp[k];
     }

   ret = 0;

#if TEST_CONTIN
   (void) test_contin (r);
#endif

   finish:
   ISIS_FREE (tmp);

   return ret;
}
/*}}}*/

static int interpolate_cont_emis (EM_t *em, EM_cont_select_t *s,  /*{{{*/
                                  EM_cont_emis_t **table,
                                  float *coef, int n, float *par,
                                  EM_cont_type_t *r)
{
   EM_cont_emis_t *t;
   float f_ioniz[ISIS_MAX_PROTON_NUMBER+1][ISIS_MAX_PROTON_NUMBER+1];
   float f_abund[ISIS_MAX_PROTON_NUMBER+1];
   int found_Z[ISIS_MAX_PROTON_NUMBER+1];
   int alt_ioniz, alt_abund, missing_total_continuum;
   int iz, iz0, iz1;
   int iq, iq0, iq1;
   int i, found_something = 0;
   int vary_rel_abund;

   if (NULL == r || NULL == table || NULL == coef)
     return -1;

   if ((s->Z == 0) && (NULL == find_cont_type (table[0], s->Z, -1)))
     missing_total_continuum = 1;
   else
     missing_total_continuum = 0;

   alt_abund = use_alt_abund (em);

   if (alt_abund)
     {
        if (-1 == get_abundance_factor (f_abund, em))
          return -1;
     }
   else
     {
        for (i=0; i <= ISIS_MAX_PROTON_NUMBER; i++)
          f_abund[i] = 1.0;
     }

   memset ((char *)found_Z, 0, sizeof(found_Z));

   vary_rel_abund = 0;
   for (i = 1; i <= ISIS_MAX_PROTON_NUMBER; i++)
     {
        if (s->rel_abun[i] != 1.0)
          {
             vary_rel_abund = 1;
             break;
          }
     }

   alt_ioniz = use_alt_ioniz (em);

   if (alt_ioniz)
     {
        EM_ioniz_table_t *t_old = em->ioniz_table[0];
        EM_ioniz_table_t *t_new = em->ioniz_table[1];

        if (-1 == get_ioniz_factor (f_ioniz, par, t_new, t_old))
          return -1;
     }
   else
     {
        int j;

        for (i=0; i <= ISIS_MAX_PROTON_NUMBER; i++)
          for (j=0; j <= ISIS_MAX_PROTON_NUMBER; j++)
            f_ioniz[i][j] = 1.0;
     }

   iz0 = iz1 = s->Z;     /* by default, execute iz, iq loop bodies */
   iq0 = iq1 = s->q;     /* once only */

   if ((s->Z == 0 && (alt_abund || alt_ioniz))
       || missing_total_continuum || vary_rel_abund)
     {
        iz0 = 1; iz1 = ISIS_MAX_PROTON_NUMBER;
     }

   /* calling routine zeros r->(stuff) */

   for (i = 0; i < n; i++)               /* loop over interpolation points */
     {
        iz = iz0;
        do
          {
             if (s->q < 0 && alt_ioniz)
               {
                  iq0 = 0;   iq1 = iz;
               }

             iq = iq0;
             do
               {
                  float weight;

                  if (NULL == (t = find_cont_type (table[i], iz, iq)))
                    continue;

                  found_something = 1;
                  found_Z[iz] = 1;

                  weight = coef[i];

                  if (iz > 0)
                    {
                       weight *= f_abund[iz] * s->rel_abun[iz];
                       if (iq >= 0)
                         weight *= f_ioniz[iz][iq];
                    }

                  if (-1 == add_cont_contrib (r, t, weight))
                    return -1;

                  iq++;
               } while (iq < iq1);        /* end loop over ions */

             iz++;
          } while (iz <= iz1);            /* end loop over elements */
     }

   for (i = 1; i <= ISIS_MAX_PROTON_NUMBER; i++)
     {
        if ((s->rel_abun[i] != 1.0)
            && (s->rel_abun[i] > 0.0) && (found_Z[i] == 0))
          {
             char str[32];
             if (0 != _DB_get_element_name (str, i))
               sprintf (str, "Z=%d", i);
             isis_vmesg (INFO, I_NOT_FOUND, __FILE__, __LINE__, "%s continuum", str);
          }
     }

   if (found_something)
     return 0;

   return -1;
}
/*}}}*/

static int get_cont_interp_points (EM_t *em, EM_cont_select_t *s, /*{{{*/
                                   int npoints, int *idx, EM_cont_emis_t **tbl)
{
   EM_filemap_t *map;
   EM_cont_data_t *cd;
   cfitsfile *fp = NULL;
   int i, vary_rel_abund;

   cd  = em->cont_data;
   map = cd->map;

   /* DB in memory */
   if (EM_Load_Cont_Emis)
     {
        for (i=0; i < npoints; i++)
          {
             tbl[i] = cd->emis[ idx[i] ];
          }
        return 0;
     }

   /* DB on disk */
   vary_rel_abund = 0;
   for (i = 1; i <= ISIS_MAX_PROTON_NUMBER; i++)
     {
        if (s->rel_abun[i] != 1.0)
          {
             vary_rel_abund = 1;
             break;
          }
     }

   if (NULL == (fp = cfits_open_file_readonly (map->filename)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", map->filename);
        return -1;
     }

   for (i=0; i < npoints; i++)
     {
        int load_all = vary_rel_abund ? 1 : 0;
        int hdu = map->hdu[ idx[i] ];

        if (-1 == cfits_movabs_hdu (hdu, fp))
          {
             isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__, "hdu=%d, %s", hdu, map->filename);
             (void) cfits_close_file (fp);
             return -1;
          }

        tbl[i] = load_cont_emissivity_hdu (fp, s->Z, s->q, load_all, em);
        if (NULL == tbl[i])
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "hdu=%d, %s", hdu, map->filename);
             (void) cfits_close_file (fp);
             return -1;
          }
     }
   (void) cfits_close_file (fp);

   return 0;
}

/*}}}*/

int EM_get_continuum (EM_cont_type_t *cont, float *par, EM_cont_select_t *s, EM_t *em) /*{{{*/
{
   EM_cont_data_t *cd;
   EM_cont_emis_t *tbl[4] = {NULL, NULL, NULL, NULL};
   float coef[4] = {0.0, 0.0, 0.0, 0.0};
   int idx[4], npoints;
   int ret = -1;

   if (NULL == em || NULL == cont)
     return -1;

   memset ((char *)cont->true_contin, 0, cont->nbins * sizeof(double));
   memset ((char *)cont->pseudo, 0, cont->nbins * sizeof(double));

   cd  = em->cont_data;

   /* if no continuum data, return zeroed cont structure */
   if (NULL == cd)
     return 0;

   if ((s->q > s->Z)
       || (s->Z > ISIS_MAX_PROTON_NUMBER))
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "Z=%d  q=%d continuum", s->Z, s->q);
        return -1;
     }

   /* This is the hack used in the data tables to indicate a sum over
    * ions/elements.
    * In searching the tables, we have to use this q value to get a match
    * with the appropriate continuum data.
    */

   if (s->Z <= 0) s->Z = 0;
   if (s->q < 0) s->q = -1;

   if (-1 == interp_coeffs (coef, idx, &npoints, par, cd->map))
     return -1;

   if (-1 == get_cont_interp_points (em, s, npoints, idx, tbl))
     goto finish;

   if (-1 == interpolate_cont_emis (em, s, tbl, coef, npoints, par, cont))
     goto finish;

   ret = 0;
   finish:

   if (EM_Load_Cont_Emis == 0)
     {
        int j;
        for (j=0; j < npoints; j++)
          {
             free_cont_list (tbl[j]);
          }
     }

   return ret;
}

/*}}}*/

/*}}}*/

/*{{{ start/end, param queries  */

static EM_t *new_em (void) /*{{{*/
{
   EM_t *em = NULL;

   if (NULL == (em = (EM_t *) ISIS_MALLOC (sizeof(EM_t))))
     return NULL;
   memset ((char *)em, 0, sizeof(*em));

   em->db = NULL;
   em->ioniz_table[0] = NULL;
   em->ioniz_table[1] = NULL;
   em->cont_data = NULL;
   em->line_data = NULL;

   /* default to invalid abundance table */
   em->standard_abund_table = -1;
   em->chosen_abund_table = em->standard_abund_table;

   return em;
}
/*}}}*/

static void set_memory_usage_level (void) /*{{{*/
{
   EM_Load_Line_Emis = EM_Use_Memory & 0x1;
   EM_Load_Cont_Emis = EM_Use_Memory & 0x2;
}

/*}}}*/

EM_t *EM_start (EM_File_Type *f, void *cl, DB_t *db)
{
   EM_t *em = NULL;
   int ret = -1;

   set_memory_usage_level ();

   isis_vmesg (WARN, I_INITIALIZING, __FILE__, __LINE__, "emissivity data");

   if (NULL == db
       || NULL == f
       || isis_user_break()
       || NULL == (em = new_em()))
     {
        /* and unset malloc hook */
        goto free_and_return;
     }

   em->db = db;

#define HAVE_STRING(s) (((s) != NULL) && ((*s) != 0))

   if (HAVE_STRING(f->abundance))
     {
        em->abund = load_abundance_file (f->abundance);
        if (em->abund == NULL)
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", f->abundance);
             goto free_and_return;
          }
     }

   if (HAVE_STRING(f->ionization))
     {
        em->ioniz_table[0] = load_ionization_table (f->ionization);
        if (em->ioniz_table[0] == NULL)
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", f->ionization);
             goto free_and_return;
          }
     }

   if (HAVE_STRING(f->line_emis))
     {
        em->line_data = load_line_data (f->line_emis, cl, db);
        if (em->line_data == NULL)
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", f->line_emis);
             goto free_and_return;
          }
     }

   if (HAVE_STRING(f->contin_emis))
     {
        em->cont_data = load_cont_data (f->contin_emis, cl, em);
        if (em->cont_data == NULL)
          isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", f->contin_emis);
        else
          isis_vmesg (INFO, I_READ_OK, __FILE__, __LINE__, "%s", f->contin_emis);
     }

   set_standard_abund_table (em);

   ret = 0;
   free_and_return:

   if (ret || isis_user_break ())
     {
        EM_end (em);
        em = NULL;
     }

   return em;
}

/*}}}*/

void EM_end (EM_t *em) /*{{{*/
{
   if (em == NULL)
     return;

   free_ioniz_table (em->ioniz_table[0]);
   free_ioniz_table (em->ioniz_table[1]);
   free_abund_list (em->abund);
   free_line_data (em->line_data);
   free_cont_data (em->cont_data);
   ISIS_FREE (em);
}
/*}}}*/

int EM_get_index_for_temperature (void)   /* prototype in db-cie.h */ /*{{{*/
{
   return EM_TEMP;
}

/*}}}*/

int EM_get_index_for_density (void)   /* prototype in db-cie.h */ /*{{{*/
{
   return EM_DENS;
}

/*}}}*/

/*}}}*/
