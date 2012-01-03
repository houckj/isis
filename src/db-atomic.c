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

/* $Id: db-atomic.c,v 1.20 2004/02/09 11:14:17 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "cfits.h"
#include "errors.h"
#include "db-atomic.h"

#undef MAX
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

#undef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/*}}}*/

/*{{{ globals */

/* One danger with double-hashing is the case where the table size is
 * a multiple of the increment derived from hashtwo().  This can be
 * avoided by making the table size a prime (this implementation does that).
 *
 * An advantage of double-hashing is that it works well even when the hash
 * table becomes nearly full. The current implementation was tested with a
 * 90% filled table and the access time wasnt noticeably slower than when the
 * table was only half-filled.
 */

enum
{
   USE_HASHTWO=1,             /* double-hash to resolve collisions */
   NHASH_DIGITS=5,
   WARN_HASH_MISSES=128
};

static const char *Element_Name[] =
{
   "",
   "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
   "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca",
   "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
   "Ga", "Ge", "As", "Se", "Br", "Kr",
   NULL
};

int DB_Ion_Format = FMT_ROMAN;

static const char *Roman [] =
{
   "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
   "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX",
   "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII", "XXVIII", "XXIX", "XXX",
   "XXXI", "XXXII", "XXXIII", "XXXIV", "XXXV", "XXXVI", "XXXVII",
   NULL
};

static int build_internals (DB_t *db);

/*}}}*/

/*{{{ private data type definitions */

#define IONSIZE    (ISIS_MAX_PROTON_NUMBER * ISIS_MAX_PROTON_NUMBER)
#define ION(Z,q)   ((Z-1) * ISIS_MAX_PROTON_NUMBER + (q))

struct _DB_t
{
   DB_ion_t *ion[IONSIZE];
   DB_line_t *line;
   int nlines;
   int *sorted_line_index;
   DB_line_t **hash_table;
   unsigned int hash_table_size;
   int max_hash_misses;
   DB_line_group_t *line_group_table_head;
};

struct _DB_line_filter_t
{
   float wavelen_min, wavelen_max;     /* wavelength range */
   double flux_min, flux_max;        /* flux range */
   int *proton_number;               /* ptr to list of elements */
   int *ion_charge;                  /* ptr to list of ions */
   int *upper;                       /* list of energy levels */
   int *lower;
   int n_elements;
   int n_ions;
   int n_upper;
   int n_lower;
};

struct _DB_line_group_t
{
   DB_line_group_t *next;    /* ptr to next group structure in group table */
   DB_line_t **line;         /* array of ptrs into the big line list */
   char *name;               /* user-defined name string */
   int nlines;               /* number of lines in the group */
   int group;                /* group id number */
};

struct _DB_ion_t
{
   /* DB_ion_t *next; */
   DB_level_t *level;           /* ptr to allocated array of size nlevels */
   int nlevels;
   int nlines;
   int charge;
   int proton_number;
};

/*}}}*/

/*{{{ misc internal fcns */

int DB_get_atomic_weight_amu (float *wt, DB_line_t *line) /*{{{*/
{
   /* http://www.physics.curtin.edu.au/iupac/
    * http://www.chem.qmul.ac.uk/iupac/AtWt/
    * Atomic weights from the 2007 table at Pure Appl. Chem., 81, 2131-2156 (2009)
    * with changes to the values for lutetium, molybdenum, nickel, ytterbium and
    * zinc from the 2005 table
    */
   static float atomic_weight[] =
     {
        -1.0,
         1.00794,     4.002602,  6.941,      9.012182,  10.811,    /* H  - B  */
        12.0107,     14.0067,   15.9994,    18.9984032, 20.1797,   /* C  - Ne */
        22.98976928, 24.3050,   26.9815386, 28.0855,    30.973762, /* Na - P  */
        32.065,      35.453,    39.948,     39.0983,    40.078,    /* S  - Ca */
        44.955912,   47.867,    50.9415,    51.9961,    54.938045, /* Sc - Mn */
        55.845,      58.933195, 58.6934,    63.546,     65.38,     /* Fe - Zn */
        69.723,      72.64,     74.92160,   78.96,      79.904,    /* Ga - Br */
        83.798                                                     /* Kr -    */
     };
   int Z;

   if (NULL == line)
     return -1;

   Z = line->proton_number;

   if (Z < 1 || Z > ISIS_MAX_PROTON_NUMBER)
     return -1;

   *wt = atomic_weight[Z];

   return 0;
}

/*}}}*/

int _DB_get_element_name (char *name, int Z) /*{{{*/
{
   if (name == NULL
       || Z < 1 || Z > ISIS_MAX_PROTON_NUMBER)
     return -1;

   isis_strcpy(name, Element_Name[Z], 3);
   return 0;
}

/*}}}*/

int DB_get_ion_name (char *name, int size, int Z, int q, int format) /*{{{*/
{
   (void) size;

   if (name == NULL
       || Z < 1 || Z > ISIS_MAX_PROTON_NUMBER
       || q < 0 || q > Z)
     {
        isis_vmesg (INFO, I_ERROR, __FILE__, __LINE__, "Unknown ion Z=%d q=%d", Z, q);
        return -1;
     }

   memset (name, 0, size);

   switch (format)
     {
      case FMT_CHARGE:
        (void) sprintf (name, "%2s +%-2d", Element_Name[Z], q);
        break;
      case FMT_INT_ROMAN:
        (void) sprintf (name, "%2s %-2d", Element_Name[Z], q+1);
        break;
      case FMT_ROMAN:
        /* drop */
      default:
        (void) sprintf (name, "%2s %s", Element_Name[Z], Roman[q]);
        break;
     }

   return 0;
}

/*}}}*/

int _DB_get_element_Z (int *Z, char *name) /*{{{*/
{
   int z;

   if (name == NULL)
     return -1;

   for (z = 1; Element_Name[z] != NULL; z++)
     {
        if (0 == isis_strcasecmp(name, Element_Name[z]))
          {
             *Z = z;
             return 0;
          }
     }

   return -1;
}

/*}}}*/

/*}}}*/

/*{{{ branching ratios */

static DB_ion_t *find_ion (DB_t *db, int Z, int q)
{
   int k = ION(Z,q);

   if ((db == NULL) || (k < 0) || (k >= IONSIZE))
     return NULL;

   return db->ion[k];
}

static int get_branching_ratios (DB_t *db) /*{{{*/
{
   enum {NDOWN_SIZE_INCREMENT = 10};   /* wasting some memory ... */
   DB_level_t *level = NULL;
   DB_ion_t *ion = NULL;
   DB_line_t **tmp;
   size_t size, new_size;
   int i, k;

   if (NULL == db)
     return -1;

   /* Its ok if no energy level data is loaded (may be unavailable) */

   if (NULL == db->ion)
     return 0;

   for (i = 0; i < IONSIZE; i++)
     {
        if (db->ion[i] == NULL)
          continue;
        ion = db->ion[i];
        for (k=1; k <= ion->nlevels; k++)
          {
             level = &ion->level[k-1];
             level->ndown = 0;
             if (level->down != NULL)
               ISIS_FREE (level->down);
          }
     }

   for (k = 0; k < db->nlines; k++)
     {
        DB_line_t *p = &db->line[k];
        int up, lo, Z, q;

        Z = p->proton_number;
        q = p->ion_charge;
        up = p->upper_level;
        lo = p->lower_level;

        if (Z <= 0 || Z > ISIS_MAX_PROTON_NUMBER)
          return -1;

        if (q < 0 || q > Z)
          return -1;

        ion = find_ion (db, Z, q);
        if (NULL == ion)
          continue;

        if (up < 1 || up > ion->nlevels       /* FIXME should check for DR */
            || lo < 1 || lo > ion->nlevels)
          continue;

        level = &ion->level[up-1];      /* n=1 is ground */

        if (level->ndown > 0
            && 0 == (level->ndown % NDOWN_SIZE_INCREMENT))
          {
             new_size = ((level->ndown + NDOWN_SIZE_INCREMENT)
                         * sizeof (DB_line_t *));
             if (NULL == (tmp = (DB_line_t **) ISIS_REALLOC (level->down, new_size)))
               return -1;
             level->down = tmp;
          }
        else if (level->ndown == 0)
          {
             size = NDOWN_SIZE_INCREMENT * sizeof(DB_line_t *);
             if (NULL == (level->down = (DB_line_t **) ISIS_MALLOC (size)))
               return -1;
          }

        level->down[level->ndown++] = p;
     }

   /*  Cleanup pass to minimize wasted memory:
    */

   for (i = 0; i < IONSIZE; i++)
     {
        if (db->ion[i] == NULL)
          continue;
        ion = db->ion[i];
        for (k=1; k <= ion->nlevels; k++)
          {
             level = &ion->level[k-1];
             if (level->ndown > 0)
               {
                  size = level->ndown * sizeof(DB_line_t *);
                  if (NULL == (tmp = (DB_line_t **) ISIS_REALLOC (level->down, size)))
                    return -1;
                  level->down = tmp;
               }
          }
     }

   return 0;
}

/*}}}*/

int DB_print_branching_for_ion (FILE *fp, int Z, int q, DB_t *db) /*{{{*/
{
   DB_ion_t *ion;
   char ion_name[10];
   double sum;
   int j, k;

   if (NULL == db || NULL == fp)
     return -1;

   if (-1 == DB_get_ion_name (ion_name, sizeof(ion_name), Z, q, DB_Ion_Format))
     return -1;

   ion = find_ion (db, Z, q);
   if (ion == NULL)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "Energy level data not available for %s", ion_name);
        return 0;
     }

   if (fprintf (fp, "%s\n", ion_name) < 0)
     return -1;

   for (k = 1; k <= ion->nlevels; k++)
     {
        DB_level_t *lev;
        DB_line_t *line;

        if (NULL == ion->level)
          continue;
        lev = &ion->level[k-1];
        if (NULL == lev
            || lev->ndown <= 1)
          continue;

        sum = 0.0;
        for (j = 0; j < lev->ndown; j++)
          {
             line = lev->down[j];
             sum += line->A;
          }
        if (sum <= 0.0)
          continue;

        if (fprintf (fp, "upper level = %4d\n", k) < 0)
          return -1;

        if (fprintf (fp, "lower        wavelen       branch       A        index\n") < 0)
          return -1;

        for (j = 0; j < lev->ndown; j++)
          {
             line = lev->down[j];
             if (fprintf (fp, "%5d  %13.6e  %11.4e  %11.4e  %6d\n",
                          line->lower_level,
                          line->wavelen,
                          line->A / sum,
                          line->A,
                          line->indx)
                 < 0)
               return -1;
          }
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ level new/free */

static void free_level_list (DB_level_t *p, int nlevels) /*{{{*/
{
   int i;

   if (p == NULL)
     return;

   for (i=0; i < nlevels; i++)
     ISIS_FREE (p[i].down);

   ISIS_FREE (p);
}

/*}}}*/

static DB_level_t *new_level_list (int nlevels) /*{{{*/
{
   DB_level_t *p = NULL;
   int i;

   if (NULL == (p = (DB_level_t *) ISIS_MALLOC (nlevels * sizeof(DB_level_t))))
     return NULL;

   for (i=0; i < nlevels; i++)
     {
        p[i].down = NULL;
        p[i].ndown = 0;
     }

   return p;
}

/*}}}*/

static void free_label_list (char **p, int n) /*{{{*/
{
   int i;

   for (i=0; i<n; i++)
     ISIS_FREE (p[i]);

   ISIS_FREE (p);
}

/*}}}*/

static char **new_label_list (int size, int n) /*{{{*/
{
   int i;
   char **p;

   if (NULL == (p = (char **) ISIS_MALLOC (n * sizeof(char *))))
     return NULL;
   memset ((char *)p, 0, n * sizeof(char *));

   for (i=0; i < n; i++)
     {
        if (NULL == (p[i] = (char *) ISIS_MALLOC (size * sizeof(char))))
          {
             free_label_list (p, n);
             break;
          }
     }

   return p;
}

/*}}}*/

/*}}}*/

/*{{{ load energy levels */

static int load_level_list (DB_level_t *level_list, int n_levels, cfitsfile *p) /*{{{*/
{
   float *energy, *stat_weight;
   int *n_quan, *L_quan;
   float *S_quan;
   DB_level_t *level = NULL;
   char **label = NULL;
   int start_row = 1;
   int i;
   int ret = -1;

   energy = stat_weight = NULL;
   n_quan = L_quan = NULL;
   S_quan = NULL;

   if (NULL == level_list)
     return -1;

   if (NULL == (energy = (float *) ISIS_MALLOC (n_levels * sizeof(float)))
       || NULL == (stat_weight = (float *) ISIS_MALLOC (n_levels * sizeof(float)))
       || NULL == (n_quan = (int *) ISIS_MALLOC (n_levels * sizeof(int)))
       || NULL == (L_quan = (int *) ISIS_MALLOC (n_levels * sizeof(int)))
       || NULL == (S_quan = (float *) ISIS_MALLOC (n_levels * sizeof(float))))
     goto free_and_return;

   if (NULL == (label = new_label_list (LEVEL_NAME_SIZE+1, n_levels)))
       goto free_and_return;

   if (-1 == cfits_read_float_col (energy, n_levels, start_row, "Energy", p)
        || -1 == cfits_read_float_col (stat_weight, n_levels, start_row, "Lev_Deg", p)
             || -1 == cfits_read_string_col (label, n_levels, start_row, "Elec_Config", p))
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "energy levels.");
        goto free_and_return;
     }

   if (   cfits_col_exist ("N_quan", p)
       && cfits_col_exist ("L_quan", p)
       && cfits_col_exist ("S_quan", p))
     {
        if (-1 == cfits_read_int_col (n_quan, n_levels, start_row, "N_quan", p)
            || -1 == cfits_read_int_col (L_quan, n_levels, start_row, "L_quan", p)
            || -1 == cfits_read_float_col (S_quan, n_levels, start_row, "S_quan", p))
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "energy level quantum numbers.");
             goto free_and_return;
          }
     }
   else
     {
        for (i=0; i < n_levels; i++)
          {
             n_quan[i] = L_quan[i] = -1;
             S_quan[i] = -1.0;
          }
     }

   /* IMPORTANT:  the ground state is listed as level = 1 in the
    * data files, but we're storing it at array index 0 (in the
    * array of DB_level_t structures).
    */

   for(i = 0; i < n_levels; ++i)
     {
        char *s;
        level = &level_list[i];
        level->indx = i+1;
        level->energy = energy[i];
        level->stat_weight = stat_weight[i];
        level->n = n_quan[i];
        level->L = L_quan[i];
        level->S = S_quan[i];
        s = strchr (label[i], ' ');
        if (s != NULL) *s = 0;
        isis_strcpy (level->label, label[i], LEVEL_NAME_SIZE);
     }

   ret = 0;

   free_and_return:

   ISIS_FREE (energy);
   ISIS_FREE (stat_weight);
   ISIS_FREE (n_quan);
   ISIS_FREE (L_quan);
   ISIS_FREE (S_quan);
   free_label_list (label, n_levels);

   return ret;
}

/*}}}*/

static void free_ion (DB_ion_t *ion) /*{{{*/
{
   if (ion == NULL) return;
   free_level_list (ion->level, ion->nlevels);
   ISIS_FREE (ion);
}

/*}}}*/

static DB_ion_t *new_ion (int nlevels) /*{{{*/
{
   DB_ion_t *ion;

   if (NULL == (ion = (DB_ion_t *) ISIS_MALLOC (sizeof(DB_ion_t))))
     return NULL;
   memset ((char *)ion, 0, sizeof (*ion));

   ion->nlines = 0;
   ion->charge = 0;
   ion->proton_number = 0;
   /* ion->next = NULL; */
   ion->nlevels = nlevels;

   if (nlevels <= 0)
     return ion;

   ion->level = new_level_list (nlevels);
   if (ion->level == NULL)
     {
        free_ion (ion);
        ion = NULL;
     }

   return ion;
}

/*}}}*/

static DB_ion_t *load_ion (char *file) /*{{{*/
{
   DB_ion_t *ion = NULL;
   cfitsfile *p = NULL;
   int Z, q, nlev;
   int ok = 0;

   if (NULL == (p = cfits_open_file_readonly (file)))
     return NULL;

   if (-1 == cfits_movabs_hdu (2, p)
       || -1 == cfits_read_int_keyword (&Z, "ELEMENT", p)
       || -1 == cfits_read_int_keyword (&q, "ION_STAT", p)
       || -1 == cfits_read_int_keyword (&nlev, "N_LEVELS", p))
     {
        goto finish;
     }

   if (NULL == (ion = new_ion (nlev)))
     goto finish;

   ion->proton_number = Z;
   ion->charge = q;

   if ((nlev > 0)
       && (-1 == load_level_list (ion->level, ion->nlevels, p)))
     goto finish;

   ok = 1;
   finish:

   if (!ok) free_ion (ion);
   (void) cfits_close_file (p);

   return ion;
}

/*}}}*/

static int load_energy_levels (DB_t *db, char **files) /*{{{*/
{
   unsigned int n;

   if (files == NULL)
     return 0;

   n = 0;
   while (*files && (0 == isis_user_break()))
     {
        DB_ion_t *ion = load_ion (*files);
        int Z, q, o;

        if (ion == NULL)
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", *files);
             return -1;
          }

        Z = ion->proton_number;
        q = ion->charge;
        o = ION(Z,q);

        if (o < 0 || o >= IONSIZE)
          return -1;

        db->ion[o] = ion;
        files++;
        n++;
        if (Isis_Verbose >= WARN)
          fprintf (stderr, "Read %d energy level files\r", n);
     }
   if (Isis_Verbose >= WARN)
     fputc ('\n', stderr);

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ misc. public functions */

int DB_get_nlines (DB_t *db) /*{{{*/
{
   if (db == NULL)
     return -1;
   else
     return db->nlines;
}

/*}}}*/

DB_ion_t *DB_get_ion (DB_t *db, int Z, int q) /*{{{*/
{
   return find_ion (db, Z, q);
}

/*}}}*/

int DB_get_ion_nlevels (int *nlevels, DB_ion_t *ion) /*{{{*/
{
   if (ion == NULL)
     return -1;

   *nlevels = ion->nlevels;

   return 0;
}

/*}}}*/

DB_level_t *DB_get_ion_level (int indx, DB_ion_t *ion) /*{{{*/
{
   if (ion == NULL || ion->level == NULL)
     return NULL;

   if (indx > ion->nlevels)
     return NULL;

   return &ion->level[ indx - 1 ];    /* index = 1 is ground */
}

/*}}}*/

int DB_get_level_label (char *label, int Z, int ion_charge, /*{{{*/
                        int indx, DB_t *db)
{
   DB_ion_t *ion;
   DB_level_t *lev;

   if (NULL == label)
     return -1;

   *label = '\0';

   if (NULL == db || NULL == db->ion
       || Z < 1 || Z > ISIS_MAX_PROTON_NUMBER
       || ion_charge < 0 || ion_charge > Z)
     return -1;

   ion = find_ion (db, Z, ion_charge);
   if (NULL == ion)
     return 0;

   if (indx < 1
       || indx > ion->nlevels
       || NULL == (lev = &ion->level[indx-1]))     /* index = 1 is ground */
     return -1;

   isis_strcpy (label, lev->label, LEVEL_NAME_SIZE);

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ hash functions */

unsigned int DB_hash2 (int Z, int q) /*{{{*/
{
   return (11 * Z + q) % 53 + 1;
}

/*}}}*/

unsigned int DB_hash (float wl, int Z, int q, int upper, int lower, /*{{{*/
                       unsigned int size)
{
   float xp;
   unsigned int iwl, h;

   /* Dangerous to use _all_ the wavelength bits directly?
    * Probably round-off errors would cause problems..
    * Can't skip the wavelength completely because the energy levels
    * may not be available.
    */

   xp = (float) (NHASH_DIGITS - 1 - ((int) log10 (wl)));
   iwl = (unsigned int) (wl * pow (10.0, xp));

   h = ((    (unsigned int)upper) <<  4) % size;
   h = ((h + (unsigned int)q    ) << 12) % size;
   h = ((h + (unsigned int)lower) <<  4) % size;
   h = ((h + (unsigned int)Z    ) << 12) % size;
   h =  (h + iwl) % size;

   return h;
}

/*}}}*/

static unsigned int hash (DB_line_t *p, unsigned int size) /*{{{*/
{
   return DB_hash (p->wavelen, p->proton_number, p->ion_charge,
                   p->upper_level, p->lower_level, size);
}

/*}}}*/

static int build_hash_table (DB_t *db) /*{{{*/
{
   static unsigned int prime[] = PRIME_LIST;
   DB_line_t **t = NULL;
   unsigned int size;
   int k;

   if (NULL == db || db->nlines <= 0)
     return -1;

   /* Access should be faster if the hash table is less than half-full
    * and has a prime number of elements.  However, if space is at a premium,
    * this double-hashing algorithm should allow the table to become
    * nearly full without causing a large penalty in access time.
    */

   size = (unsigned int) (2 * db->nlines + 1);
   for (k=0; prime[k] < size; k++)
     ;
   size = prime[k];

   if (NULL == (t = (DB_line_t **) ISIS_MALLOC (size * sizeof(DB_line_t *))))
     return -1;
   memset ((char *) t, 0, size * sizeof(DB_line_t *));

   db->max_hash_misses = 0;

   for (k=0; k < db->nlines; k++)
     {
        DB_line_t *line = &db->line[k];
        unsigned int h;
        int miss = 0;
        int incr = 1;

        h = hash (line, size);

        if (USE_HASHTWO)
          incr = DB_hash2 (line->proton_number, line->ion_charge);

        while (t[h] != NULL)
          {
             h = (h + incr) % size;
             if (miss++ > db->nlines)
               {
                  isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                              "Database corruption? (problem building hash table, hash_misses = %d)",
                              miss);
                  ISIS_FREE (t);
                  return -1;
               }
          }
        t[h] = line;

        db->max_hash_misses = MAX((int) miss, db->max_hash_misses);
     }

   if (db->hash_table != NULL)
     ISIS_FREE (db->hash_table);

   db->hash_table = t;
   db->hash_table_size = size;

   if (db->max_hash_misses > WARN_HASH_MISSES)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "Database corruption? (problem building hash table, max_hash_misses = %d)",
                    db->max_hash_misses);
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ load line lists  */

int DB_merge_lines (DB_t *db, DB_Merge_Type *m) /*{{{*/
{
   int i, offset, len;

   if (NULL == db || m == NULL || m->n <= 0)
     return -1;

   if (NULL == db->line)
     {
        offset = 0;
        if (NULL == (db->line = (DB_line_t *) ISIS_MALLOC (m->n * sizeof(DB_line_t))))
          return -1;
     }
   else
     {
        DB_line_t *tmp;
        offset = db->nlines;
        len = (db->nlines + m->n) * sizeof(DB_line_t);
        if (NULL == (tmp = (DB_line_t *) ISIS_REALLOC (db->line, len)))
          return -1;
        db->line = tmp;
     }

   for (i=0; i < m->n; i++)
     {
        DB_line_t *p = &db->line[offset + i];
        p->A = m->A_value[i];
        p->A_err = m->A_err[i];
        p->wavelen = m->lambda[i];
        p->wavelen_err = m->lambda_err[i];
        p->upper_level = m->upper_lev[i];
        p->lower_level = m->lower_lev[i];
        p->proton_number = m->proton_number[i];
        p->ion_charge = m->ion_charge[i];
        p->flux = (double) 0.0;
        p->have_emissivity_data = 0;
        p->source_index = -1;
        p->type_index = -1;
        p->indx = i+offset;
     }

   db->nlines += m->n;

   return build_internals (db);      /* re-build sort tables, etc. */
}

/*}}}*/

static int Last_Line_Index;

static int load_line_list_extension (DB_line_t **line, int *n, cfitsfile *fp) /*{{{*/
{
   int i, nlines, k, opt, start_row, crap;
   float *wavelen_obs, *wavelen, *wavelen_err, *A_value, *A_err;
   int *upper_lev, *lower_lev;
   int element, ion_state;
   int len;
   int ret = -1;

   if (NULL == line)
     return -1;

   wavelen = wavelen_obs = wavelen_err = A_value = A_err = NULL;
   upper_lev = lower_lev = NULL;

   if (-1 == cfits_read_int_keyword (&nlines, "NAXIS2", fp)
       || -1 == cfits_read_int_keyword (&element, "ELEMENT", fp)
       || -1 == cfits_read_int_keyword (&ion_state, "ION_STAT", fp))
     return -1;

   if (nlines == 0)
     return 0;

   if (NULL == *line)
     {
        if (NULL == (*line = (DB_line_t *) ISIS_MALLOC (nlines * sizeof(DB_line_t))))
          return -1;
     }
   else
     {
        DB_line_t *tmp;
        len = (*n + nlines) * sizeof(DB_line_t);
        if (NULL == (tmp = (DB_line_t *) ISIS_REALLOC (*line, len)))
          return -1;
        *line = tmp;
     }

/*
 * Now, read these columns:
 *
 * ___Column_Names_________Formats______Dims______Units___
 * Upper_Lev                  1J
 * Lower_Lev                  1J
 * Wavelen                    1E                  A
 * Wave_Err                   1E                  A
 * Einstein_A                 1E                  s**-1
 * Ein_A_Err                  1E                  s**-1
 */

   if (-1 == (opt = cfits_optimal_numrows (fp)))
     goto free_and_return;

  if (NULL == (wavelen = (float *) ISIS_MALLOC (opt * sizeof(float)))
      || NULL == (wavelen_obs = (float *) ISIS_MALLOC (opt * sizeof(float)))
      || NULL == (wavelen_err = (float *) ISIS_MALLOC (opt * sizeof(float)))
      || NULL == (A_value = (float *) ISIS_MALLOC (opt * sizeof(float)))
      || NULL == (A_err = (float *) ISIS_MALLOC (opt * sizeof(float)))
      || NULL == (upper_lev = (int *) ISIS_MALLOC (opt * sizeof(int)))
      || NULL == (lower_lev = (int *) ISIS_MALLOC (opt * sizeof(int))))
     goto free_and_return;

   crap = 0;
   k = 0;
   while (nlines - k*opt > 0)
     {
        int nread, idx;

        start_row = k * opt;
        nread = MIN (opt, nlines - start_row);

        if (-1 == cfits_read_float_col (wavelen, nread, 1+start_row, "WAVELEN", fp)
            || -1 == cfits_read_float_col (wavelen_obs, nread, 1+start_row, "WAVE_OBS", fp)
            || -1 == cfits_read_float_col (wavelen_err, nread, 1+start_row, "Wave_Err", fp)
            || -1 == cfits_read_float_col (A_value, nread, 1+start_row, "Einstein_A", fp)
            || -1 == cfits_read_float_col (A_err, nread, 1+start_row, "Ein_A_err", fp)
            || -1 == cfits_read_int_col (upper_lev, nread, 1+start_row, "Upper_Lev", fp)
            || -1 == cfits_read_int_col (lower_lev, nread, 1+start_row, "Lower_Lev", fp))
          {
             isis_vmesg (FAIL, I_WARNING, __FILE__, __LINE__, "Error reading wavelength file");
             goto free_and_return;
          }

        idx = Last_Line_Index;

        for (i=0; i < nread; i++)
          {
             DB_line_t *p;

             /* skip lines with wavelen <= 0 or INDEF */
             if (wavelen[i] <= FLT_MIN)
               {
                  crap++;
                  continue;
               }

             p = &(*line)[idx];

             /* FLT_MIN flags NULL values (INDEF) */
             if (wavelen_obs[i] > FLT_MIN)
               p->wavelen = wavelen_obs[i];
             else
               p->wavelen = wavelen[i];

             p->A = A_value[i];
             p->A_err = A_err[i];
             p->wavelen_err = wavelen_err[i];
             p->upper_level = upper_lev[i];
             p->lower_level = lower_lev[i];
             p->proton_number = (unsigned char) element;
             p->ion_charge = (unsigned char) ion_state;
             p->flux = 0.0;              /* initial value */
             p->have_emissivity_data = 0;      /* assume not */
             p->source_index = -1;              /* not used */
             p->type_index = -1;                /* not used */
             p->indx = idx;                    /* save array index */

             idx++;       /* if no crap, same as o+i */
          }

        Last_Line_Index = idx;

        k++;
     }

   *n += (nlines - crap);         /* dont bother resizing if crap > 0 */

   ret = 0;

   free_and_return:

   ISIS_FREE (wavelen);
   ISIS_FREE (wavelen_obs);
   ISIS_FREE (wavelen_err);
   ISIS_FREE (A_value);
   ISIS_FREE (A_err);
   ISIS_FREE (upper_lev);
   ISIS_FREE (lower_lev);

   return ret;
}

/*}}}*/

static int load_line_list_file (DB_line_t **line, int *nlines, char *filename) /*{{{*/
{
   cfitsfile *fp;

   if (NULL == (fp = cfits_open_file_readonly (filename))
       || -1 == cfits_movabs_hdu(1, fp))
     return -1;

   for (;;)
     {
        if ((-1 == cfits_movrel_hdu(1, fp))
            || (-1 == load_line_list_extension (line, nlines, fp)))
          break;
     }

   (void) cfits_close_file (fp);

   return 0;
}

/*}}}*/

static int load_line_list (DB_t *db, char **files) /*{{{*/
{
   DB_line_t *line = NULL;
   unsigned int n;
   int nlines = 0;

   if (NULL == files)
     return 0;

   if (NULL == db)
     return -1;

   Last_Line_Index = 0;

   n = 0;
   while (*files && (0 == isis_user_break()))
     {
        if (-1 == load_line_list_file (&line, &nlines, *files))
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s", *files);
             ISIS_FREE(line);
             return -1;
          }
        files++;
        n++;
        if (Isis_Verbose >= WARN)
          fprintf (stderr, "Read %d wavelength files\r", n);
     }
   if (Isis_Verbose >= WARN)
     fputc ('\n', stderr);

   db->line = line;
   db->nlines = nlines;

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ functions to filter lists of lines */

/*{{{ heap-based selection */

#define SWAP(a,b,t)  { t tmp = (a); (a) = (b); (b) = tmp; }

/* re-establish (minimum-oriented) heap property
 * moving upward from position k */

static void heapify_min_upward (int *idx, int k, double *a) /*{{{*/
{
   while (k > 0 && a[ idx[k] ] < a[ idx[k/2] ])
     {
        SWAP (idx[k], idx[k/2], int)
        k /= 2;
     }
}

/*}}}*/

/* re-establish (minimum-oriented) heap property
 * moving downward from position k */

static void heapify_min_downward (int *idx, int k, int n, double *a) /*{{{*/
{
   int j;

   while (2*k < n)
     {
        j = 2*k;
        if (j+1 < n && a[ idx[j+1] ] < a[ idx[j] ])
          j++;
        if (a[ idx[k] ] <= a[ idx[j] ])
          break;
        SWAP (idx[k], idx[j], int)
        k = j;
     }
}

/*}}}*/

/* See Sedgewick "Algorithms" 3rd ed., p. 379 discussion of
 * selection using heaps.  If k = n, this amounts to a heapsort,
 * although not a maximally efficient one.
 */

static int index_select_k_largest (int *idx, int *k_got, int k_want, int n, double *a) /*{{{*/
{
   int k, i=0;
   int *t = NULL;

   if (k_want <= 0 || n <= 0)
     return -1;

   k = *k_got = MIN(n, k_want);        /* best approximation to request */

   if (NULL == a || NULL == idx)
     return -1;

   if (NULL == (t = (int *) ISIS_MALLOC (k * sizeof(int))))
     return -1;

   /* Make a minimum-oriented heap of size k */

   while (i<k)
     {
        t[i] = i;
        heapify_min_upward (t, i++, a);
     }

   /* Scan the remaining elements:
    * If an element is larger than the current minimum heap element,
    * it replaces the current minimum heap element and the heap is
    * re-established.
    */

   for (i=k; i < n; i++)
     {
        if (a[i] <= a[ t[0] ])
          continue;
        t[0] = i;
        heapify_min_downward (t, 0, k, a);
     }

   /* Heap now contains the k largest -- pop them off in descending order */

   while (k-- > 0)
     {
        idx[k] = t[0];
        SWAP (t[0], t[k], int)
        heapify_min_downward (t, 0, k, a);
     }

   ISIS_FREE (t);

   return 0;
}

/*}}}*/

/*}}}*/

int DB_get_k_brightest (int *t, int *n, int k, int *list, int nlines, /*{{{*/
                        DB_t *db)
{
   double *flux = NULL;
   int *indx = NULL;
   int num_select = k;
   int num_got;
   int ret = -1;
   int i;

   if (db == NULL || t == NULL || n == NULL)
     return -1;

   /* t is assumed to be of size k or larger */

   *n = 0;
   memset ((char *)t, 0, k * sizeof(int));

   if (NULL == (indx = (int *) ISIS_MALLOC (num_select * sizeof(int)))
       || NULL == (flux = (double *) ISIS_MALLOC (nlines * sizeof(double))))
     goto finish;

   for (i=0; i < nlines; i++)
     {
        DB_line_t *line = DB_get_line_from_index (list[i], db);
        if (line == NULL)
          goto finish;
        flux[i] = line->flux;
     }

   if (-1 == index_select_k_largest (indx, &num_got, num_select, nlines, flux))
     goto finish;

   for (i=0; i < num_got; i++)
     {
        if ( flux[ indx[i] ] > 0.0 )
          t[(*n)++] = list[ indx[i] ];
     }

   ret = 0;

   if (*n == 0)
     isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "all lines have flux = 0");
   else if (*n < num_select)
     isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "only %d lines had flux > 0", *n);

   finish:
   ISIS_FREE (flux);
   ISIS_FREE (indx);
   return ret;
}

/*}}}*/

static int bsearch_sorted_index (DB_line_t *t, DB_t *db) /*{{{*/
{
   int n0, n1, n2, k;
   DB_line_t *ln;

   n0 = 0;
   n1 = db->nlines;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;

        k = db->sorted_line_index[n2];
        ln = &db->line[k];

        if (t->wavelen <= ln->wavelen)
          {
             if (ln->wavelen == t->wavelen) return n2;      /* I know -- risky */
             n1 = n2;
          }
        else
          {
             n0 = n2;
          }
     }

   return n0;
}

/*}}}*/

int DB_get_unblended (int *t, int *n, float f, float wl_f, /*{{{*/
                      unsigned int mask_allowed,
                      int *list, int nlines, DB_t *db)
{
   int i;

   if (db == NULL || t == NULL || n == NULL || list == NULL)
     return -1;

   /* t is assumed to be of size >= nlines */

   *n = 0;
   memset ((char *)t, 0, nlines * sizeof(int));

   for (i=0; i < nlines; i++)
     {
        DB_line_t *line;
        int j, idx;
        double contamination;

        if (NULL == (line = DB_get_line_from_index (list[i], db))
            || -1 == (idx = bsearch_sorted_index (line, db)))
          return -1;

        contamination = 0.0;

        for (j = idx-1; 0 < j; j--)
          {
             int k  = db->sorted_line_index[j];
             DB_line_t *ln = &db->line[k];
             if (fabs (1.0 - ln->wavelen / line->wavelen) > wl_f)
               break;
             if (mask_allowed == 0
                 || ((mask_allowed & SAME_ELEM)
                  && (line->proton_number != ln->proton_number))
                 || ((mask_allowed & SAME_ION)
                     && (line->ion_charge != ln->ion_charge)))
               contamination += ln->flux;
          }

        for (j = idx+1; j < db->nlines; j++)
          {
             int k = db->sorted_line_index[j];
             DB_line_t *ln = &db->line[k];
             if (fabs (1.0 - ln->wavelen / line->wavelen) > wl_f)
               break;
             if (mask_allowed == 0
                 || ((mask_allowed & SAME_ELEM)
                     && (line->proton_number != ln->proton_number))
                 || ((mask_allowed & SAME_ION)
                     && (line->ion_charge != ln->ion_charge)))
               contamination += ln->flux;
          }

        if (contamination < f * line->flux)
          t[(*n)++] = list[i];
     }

   isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "got %d unblended lines", *n);

   return 0;
}

/*}}}*/

int DB_apply_line_filter (char *t,  DB_line_filter_function_t line_in_set, /*{{{*/
                          DB_line_filter_t *filter, DB_t *db)
{
   int i;

   if (NULL == db || NULL == line_in_set || NULL == filter || NULL == t)
     return -1;

   /* t is assumed to be a character array of size db->nlines */

   for (i=0; i < db->nlines; i++)
     {
       if (line_in_set (&db->line[i], filter))
          t[i] = 0x01;
        else
          t[i] = 0x00;
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ line group manipulation */

static DB_line_group_t *new_line_group (void) /*{{{*/
{
   DB_line_group_t *table;

   if (NULL == (table = (DB_line_group_t *) ISIS_MALLOC (sizeof(DB_line_group_t))))
     return NULL;
   memset ((char *)table, 0, sizeof (*table));

   table->line = NULL;
   table->name = NULL;
   table->group = INT_MAX;
   table->next = NULL;

   return table;
}

/*}}}*/

void DB_free_line_group (DB_line_group_t *t) /*{{{*/
{
   if (t == NULL)
     return;
   ISIS_FREE (t->name);
   ISIS_FREE (t->line);
   ISIS_FREE (t);
}

/*}}}*/

static void delete_line_group (int group, DB_line_group_t *head) /*{{{*/
{
   DB_line_group_t *t;
   DB_line_group_t *next;

   if (NULL == head)
     return;

   for (t=head; t->next != NULL; t=t->next)
     {
        if (group == t->next->group)
          {
             next = t->next->next;
             DB_free_line_group (t->next);
             t->next = next;
             return;
          }
     }
}

/*}}}*/

void DB_delete_line_group (int group, DB_t *db) /*{{{*/
{
   delete_line_group (group, db->line_group_table_head);
}

/*}}}*/

static void free_line_group_table (DB_line_group_t *head) /*{{{*/
{
   DB_line_group_t *t = head;
   DB_line_group_t *dead;

   while (t != NULL)
     {
        dead = t;
        t = t->next;
        DB_free_line_group (dead);
     }

   head = NULL;
}

/*}}}*/

static DB_line_group_t *find_group(int group, DB_line_group_t *head) /*{{{*/
{
   DB_line_group_t *t = NULL;

   if (NULL == head)
     return NULL;

   for (t=head->next; t != NULL; t = t->next)
     {
        if (group == t->group)
          return t;
     }

   return NULL;
}

/*}}}*/

static void append_group(DB_line_group_t *table, DB_line_group_t *head) /*{{{*/
{
   DB_line_group_t *t;

   for (t=head; t->next != NULL; t=t->next)
     {
     }

   t->next = table;
   table->next = NULL;
}

/*}}}*/

DB_line_group_t *DB_find_group(int group, DB_t *db) /*{{{*/
{
   return find_group (group, db->line_group_table_head);
}

/*}}}*/

DB_line_group_t *DB_get_group_table_head (DB_t *db) /*{{{*/
{
   return db->line_group_table_head;
}

/*}}}*/

DB_line_group_t *DB_next_group (DB_line_group_t *cl) /*{{{*/
{
   if (NULL == cl)
     return NULL;
   else
     return cl->next;
}

/*}}}*/

int DB_set_group_name (char *name, DB_line_group_t *cl) /*{{{*/
{
   if (NULL == cl || NULL == name)
     return -1;

   ISIS_FREE (cl->name);

   if (NULL == (cl->name = isis_make_string (name)))
     return -1;
   else
     return 0;
}

/*}}}*/

int DB_get_group_name (char **name, DB_line_group_t *cl) /*{{{*/
{
   if (cl->name == NULL)
     *name = NULL;

   if (NULL == (*name = isis_make_string (cl->name)))
     return -1;
   else
     return 0;
}

/*}}}*/

int DB_get_group_index (int *group, DB_line_group_t *cl) /*{{{*/
{
   if (NULL == cl)
     return -1;

   *group = cl->group;

   return 0;
}

/*}}}*/

int DB_get_group_line_list (int *nlines, DB_line_t ***line, DB_line_group_t *cl) /*{{{*/
{
   if (cl == NULL)
     return -1;

   *nlines = cl->nlines;
   *line = cl->line;

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ line group member-list manipulation */

int DB_edit_group (DB_line_group_t **g, int *list, int n, int add, DB_t *db) /*{{{*/
{
   DB_line_group_t *g_new = NULL;
   DB_line_t **t = NULL;
   char *flag = NULL;
   char flag_value;
   int idx, n_new;
   int j, k;

   if (NULL == g || NULL == *g || NULL == list || NULL == db || n <= 0 )
     return -1;

   /* maximum anticipated size of new group */

   n_new = (*g)->nlines + (add ? n : 0);

   if (NULL == (g_new = new_line_group ())
       || NULL == (t = (DB_line_t **) ISIS_MALLOC (n_new * sizeof(DB_line_t *)))
       || NULL == (flag = (char *) ISIS_MALLOC (db->nlines * sizeof(char))))
     {
        DB_free_line_group (g_new);
        ISIS_FREE (t);
        ISIS_FREE (flag);
        return -1;
     }

   /* initialize flag array:   1 for current group members, 0 otherwise */

   memset (flag, 0, db->nlines * sizeof(char));
   for (j=0; j < (*g)->nlines; j++)
     {
        idx = (*g)->line[j]->indx;
        flag[idx] = 0x01;
     }

   /* Now add/delete by setting flag array elements -- this avoids
    * duplication problems and the time complexity is O(N).
    * While simple to code and good for manipulating large lists
    * efficiently, this is overkill for small edits.
    * Might need an alternate algorithm for that case if this proves
    * too slow.
    */

   flag_value = add ? 0x01 : 0x00;

   for (j=0; j < n; j++)
     {
        idx = list[j];
        if (0 <= idx && idx < db->nlines)
          flag[ idx ] = flag_value;
        else
          isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "line %d ignored (no data)", idx);
     }

   /* assemble edited list */

   k = 0;
   for (j=0; j < db->nlines; j++)
     {
        if (flag[j])
          t[k++] = &db->line[j];
     }

   ISIS_FREE (flag);

   /* resize if necessary */
   if (k <= 0)
     {
        DB_free_line_group (*g);
        DB_free_line_group (g_new);
        ISIS_FREE (t);
        return 0;
     }
   else if (k != n_new)
     {
        DB_line_t **tmp;
        if (NULL == (tmp = (DB_line_t **) ISIS_REALLOC (t, k * sizeof(DB_line_t *))))
          {
             ISIS_FREE (t);
             DB_free_line_group (g_new);
             return -1;
          }
        t = tmp;
     }

   g_new->line = t;
   g_new->nlines = k;

   /* keep old id number */
   g_new->group = (*g)->group;

   /* new group replaces old group */
   DB_free_line_group (*g);
   *g = g_new;

   return 0;
}

/*}}}*/

char *DB_flag_array_from_list (int *list, int n, DB_t *db) /*{{{*/
{
   char *flag = NULL;
   int i;

   if (NULL == list || n <= 0 || NULL == db || db->nlines <= 0)
     return NULL;

   if (NULL == (flag = (char *) ISIS_MALLOC (db->nlines * sizeof(char))))
     return NULL;

   memset (flag, 0, db->nlines * sizeof(char));

   for (i=0; i < n; i++)
     {
        int idx = list[i];
        if (0 <= idx && idx < db->nlines)
          flag[ idx ] = 0x01;
        else
          isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "line %d ignored (no data)", idx);
     }

   return flag;
}

/*}}}*/

static DB_line_t **Group_Sort_Ptr = NULL;

static int compare_wavelen2 (const void *a, const void *b) /*{{{*/
{
   const int *i = (const int *) a;
   const int *k = (const int *) b;
   float wla = Group_Sort_Ptr[ *i ]->wavelen;
   float wlb = Group_Sort_Ptr[ *k ]->wavelen;
   if (wla > wlb) return +1;
   if (wla < wlb) return -1;
   return 0;
}

/*}}}*/

static int sort_group (DB_line_t **g, int nlines) /*{{{*/
{
   int *index_array = NULL;
   DB_line_t **tmp = NULL;
   int k;

   if (NULL == g || nlines <= 0)
     return -1;

   if (NULL == (index_array = (int *) ISIS_MALLOC (nlines * sizeof(int)))
       || NULL == (tmp = (DB_line_t **) ISIS_MALLOC (nlines * sizeof(DB_line_t *))))
     {
        ISIS_FREE(index_array);
        ISIS_FREE(tmp);
        return -1;
     }

   memcpy ((char *)tmp, (char *)g, nlines * sizeof(DB_line_t *));

   for (k = 0; k < nlines; k++)
     index_array[k] = k;

   Group_Sort_Ptr = g;
   qsort (index_array, nlines, sizeof(int), compare_wavelen2);
   Group_Sort_Ptr = NULL;

   for (k = 0; k < nlines; k++)
     g[k] = tmp[index_array[k]];

   ISIS_FREE (tmp);
   ISIS_FREE (index_array);
   return 0;
}

/*}}}*/

DB_line_group_t *DB_make_group_from_list (int group, int *list, int n, DB_t *db) /*{{{*/
{
   DB_line_group_t *g = NULL;
   DB_line_t **t = NULL;
   int i, k;
   char *flag = NULL;

   if (NULL == db || NULL == list || n <= 0)
     return NULL;

   if (NULL == (g = new_line_group ())
       || NULL == (t = (DB_line_t **) ISIS_MALLOC (n * sizeof(DB_line_t *)))
       || NULL == (flag = DB_flag_array_from_list (list, n, db)))
     {
        DB_free_line_group (g);
        ISIS_FREE (t);
        ISIS_FREE (flag);
        return NULL;
     }

   /* Use flag array algorithm to avoid duplications in new group
    * while retaining O(N) time complexity.  This algorithm is
    * overkill for small groups, but great for big groups.
    * If too slow for the small-group case, we could add an optimized
    * algorithm for that case.
    */

   k = 0;
   for (i=0; i < db->nlines; i++)    /* assemble merged list */
     {
        if (flag[i])
          t[k++] = &db->line[i];
     }

   ISIS_FREE (flag);

   if (k == 0)
     {
        ISIS_FREE (t);
        DB_free_line_group (g);
        return NULL;
     }
   else if (k < n)
     {
        DB_line_t **tmp;
        if (NULL == (tmp = (DB_line_t **) ISIS_REALLOC (t, k * sizeof(DB_line_t *))))
          {
             ISIS_FREE (t);
             DB_free_line_group (g);
             return NULL;
          }
        t = tmp;
     }

   if (-1 == sort_group (t, k))
     {
        ISIS_FREE (t);
        DB_free_line_group (g);
        return NULL;
     }

   g->line = t;
   g->nlines = k;
   g->group = group;

   return g;
}

/*}}}*/

int DB_list_to_group (int group, int *list, int n, DB_t *db) /*{{{*/
{
   DB_line_group_t *g;
   DB_line_group_t *head;

   g = head = NULL;

   if (NULL == db || NULL == list || n <= 0)
     return -1;

   if (NULL == (g = DB_make_group_from_list (group, list, n, db)))
     return -1;

   if (NULL == db->line_group_table_head
       && NULL == (db->line_group_table_head = new_line_group ()))
     return -1;

   head = db->line_group_table_head;

   delete_line_group (group, head);
   append_group (g, head);

   return 0;
}

/*}}}*/

int DB_get_group_members (int **index_list, int *nmembers, int group, DB_t *db) /*{{{*/
{
   DB_line_group_t *g;
   int i;

   if (NULL == db->line_group_table_head
       || NULL == (g = find_group (group, db->line_group_table_head)))
     return -1;

   *nmembers = g->nlines;

   if (NULL == (*index_list = (int *) ISIS_MALLOC (*nmembers * sizeof(int))))
     return -1;

   for (i=0; i < g->nlines; i++)
     {
        DB_line_t *line = g->line[i];

        if (line == NULL)
          {
             ISIS_FREE (*index_list);
             *nmembers = 0;
             return -1;
          }

        (*index_list)[i] = line->indx;
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ line list filtering functions */

DB_line_filter_t *DB_new_line_filter (void) /*{{{*/
{
   DB_line_filter_t *f;

   if (NULL == (f = (DB_line_filter_t *) ISIS_MALLOC (sizeof(DB_line_filter_t))))
     return NULL;
   memset ((char *)f, 0, sizeof(*f));

   f->proton_number = NULL;
   f->ion_charge = NULL;
   f->upper = NULL;
   f->lower = NULL;

   return f;
}

/*}}}*/

void DB_free_line_filter (DB_line_filter_t *f) /*{{{*/
{
   if (f == NULL)
     return;
   ISIS_FREE (f->proton_number);
   ISIS_FREE (f->ion_charge);
   ISIS_FREE (f->upper);
   ISIS_FREE (f->lower);
   ISIS_FREE (f);
}

/*}}}*/

int DBf_wavelength (DB_line_t *line, DB_line_filter_t *f) /*{{{*/
{
   if (f->wavelen_min <= line->wavelen && line->wavelen < f->wavelen_max)
     return 1;
   else
     return 0;
}

/*}}}*/

int DBf_set_wavelength (DB_line_filter_t *f, double wavelen_min, double wavelen_max) /*{{{*/
{
   if (f == NULL)
     return -1;

   if (wavelen_max < 0.0 || wavelen_max > FLT_MAX)
     wavelen_max = FLT_MAX;

   f->wavelen_min = (float) wavelen_min;
   f->wavelen_max = (float) wavelen_max;

   return 0;
}

/*}}}*/

int DBf_flux (DB_line_t *line, DB_line_filter_t *f) /*{{{*/
{
   if (f->flux_min <= line->flux && line->flux < f->flux_max)
     return 1;
   else
     return 0;
}

/*}}}*/

int DBf_set_flux (DB_line_filter_t *f, double flux_min, double flux_max) /*{{{*/
{
   if (f == NULL)
     return -1;

   if (flux_max < 0.0)
     flux_max = DBL_MAX;

   f->flux_min = flux_min;
   f->flux_max = flux_max;

   return 0;
}

/*}}}*/

static int is_int_list_member (int value, int *list, int size) /*{{{*/
{
   int i;

   if ((list == NULL) || (size <= 0)
       || ((size == 1) && (list[0] < 0)))
     return 1;
   /* by convention, the empty list matches everything */

   for (i=0; i<size; i++)
     {
        if (value == list[i])
          return 1;
     }
   return 0;
}

/*}}}*/

static int *dup_int_list (int *from, int num) /*{{{*/
{
   int *to;

   if ((from == NULL) || (num <= 0))
     return NULL;

   to = (int *) ISIS_MALLOC (num * sizeof(int));
   if (to == NULL)
     return NULL;

   memcpy ((char *)to, (char *)from, num * sizeof(int));

   return to;
}

/*}}}*/

int DBf_el_ion (DB_line_t *line, DB_line_filter_t *f) /*{{{*/
{
   if (is_int_list_member (line->proton_number, f->proton_number, f->n_elements)
       && is_int_list_member (line->ion_charge, f->ion_charge, f->n_ions))
     return 1;
   else
     return 0;
}

/*}}}*/

int DBf_set_el_ion (DB_line_filter_t *f, int n_elem, int *Z, /*{{{*/
                    int n_ion, int *ion_charge)
{
   if (f == NULL)
     return -1;

   f->n_elements = (n_elem > 0) ? n_elem : 0;
   f->proton_number = dup_int_list (Z, n_elem);
   if ((n_elem > 0) && (f->proton_number == NULL))
     return -1;

   f->n_ions = (n_ion > 0) ? n_ion : 0;
   f->ion_charge = dup_int_list (ion_charge, n_ion);
   if ((n_ion > 0) && (f->ion_charge == NULL))
     return -1;

   return 0;
}

/*}}}*/

int DBf_trans (DB_line_t *line, DB_line_filter_t *f) /*{{{*/
{
   if (0 == DBf_el_ion (line, f)) return 0;
   if (is_int_list_member (line->upper_level, f->upper, f->n_upper)
       && is_int_list_member (line->lower_level, f->lower, f->n_lower))
     return 1;
   else
     return 0;
}

/*}}}*/

int DBf_set_trans (DB_line_filter_t *f,  int Z, int q,/*{{{*/
                  int nup, int *upper, int nlo, int *lower)
{
   if (f == NULL)
     return -1;

   if (-1 == DBf_set_el_ion (f, 1, &Z, 1, &q))
     return -1;

   f->n_upper = (nup > 0) ? nup : 0;
   f->upper = dup_int_list (upper, nup);
   if ((nup > 0) && (f->upper == NULL))
     return -1;

   f->n_lower = (nlo > 0) ? nlo : 0;
   f->lower = dup_int_list (lower, nlo);
   if ((nlo > 0) && (f->lower == NULL))
     return -1;

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ set/retrieve info on single lines */

int DB_get_line_namestring (char *s, int s_size, int type, /*{{{*/
                            DB_line_t *line, DB_t *db)
{
   char ion_name[10];
   int Z, chg, up, lo;

   *s = '\0';

   if (NULL == line)
     return -1;

   Z = line->proton_number;
   chg = line->ion_charge;
   up = line->upper_level;
   lo = line->lower_level;

   if (-1 == DB_get_ion_name (ion_name, sizeof(ion_name), Z, chg, DB_Ion_Format))
     return -1;

   if (type == SHORT_LINE_NAME)
     {
        isis_strcpy (s, ion_name, s_size);
     }
   else if (type == LONG_LINE_NAME)
     {
        char up_name[LEVEL_NAME_SIZE], lo_name[LEVEL_NAME_SIZE];
        unsigned char pn = line->proton_number;
        unsigned char q = line->ion_charge;

        if (-1 == DB_get_level_label (up_name, pn, q, up, db)
            || -1 == DB_get_level_label (lo_name, pn, q, lo, db))
          isis_strcpy (s, ion_name, s_size);
        else
          {
             isis_strcat (s, s_size, ion_name, " ", up_name, "-", lo_name, NULL);
          }
     }
   else
     {
        /* invalid type */
        return -1;
     }

   return 0;
}

/*}}}*/

DB_line_t *DB_get_line_from_index(int indx, DB_t *db) /*{{{*/
{
   if (db == NULL
       || db->line == NULL
       || indx < 0
       || db->nlines <= indx)
     return NULL;

   return &db->line[indx];
}

/*}}}*/

int DB_get_line_ion (int *Z, int *q, DB_line_t *line) /*{{{*/
{
   if (line == NULL)
     return -1;

   *Z = line->proton_number;
   *q = line->ion_charge;
   return 0;
}

/*}}}*/

int DB_zero_line_flux (DB_t *db) /*{{{*/
{
   int i;

   if (NULL == db)
     return -1;

   for (i = 0; i < db->nlines; i++)
     db->line[i].flux = 0.0;

   return 0;
}

/*}}}*/

int DB_set_line_wavelength (float wavelen, float wavelen_err, DB_line_t *line) /*{{{*/
{
   if (line == NULL)
     return -1;

   line->wavelen = wavelen;
   line->wavelen_err = wavelen_err;
   return 0;
}

/*}}}*/

DB_line_t *DB_get_line (float wavelen, int proton_number, int ion_charge, /*{{{*/
                      int upper_level, int lower_level, DB_t *db)
{
   DB_line_t *line;
   DB_line_t p;
   unsigned int h;
   int k, incr = 1;

   if (wavelen <= 0.0
       || proton_number < 1 || proton_number > ISIS_MAX_PROTON_NUMBER
       || ion_charge < 0 || ion_charge > proton_number)
     return NULL;

   if (db->hash_table_size <= 0)     /* maybe no lines were loaded yet */
     return NULL;

   k = db->max_hash_misses + 1;      /* max allowed hashes for any line */

   p.wavelen = wavelen;
   p.proton_number = (unsigned char) proton_number;
   p.ion_charge = (unsigned char) ion_charge;
   p.upper_level = upper_level;
   p.lower_level = lower_level;

   h = hash (&p, db->hash_table_size);

   if (USE_HASHTWO)
     incr = DB_hash2 (p.proton_number, p.ion_charge);

   while (k-- > 0)
     {
        line = db->hash_table [h];

        if (NULL != line
            && line->proton_number == (unsigned char) proton_number
            && line->ion_charge == (unsigned char) ion_charge
            && fabs(line->wavelen/wavelen - 1.0) < WAVELEN_TOL
            && line->upper_level == upper_level
            && line->lower_level == lower_level)
          return line;

        h = (h + incr) % db->hash_table_size;
     }

   return NULL;
}

/*}}}*/

DB_line_t *DB_get_line_by_indices (int proton_number, int ion_charge, /*{{{*/
                                 int upper_level, int lower_level, DB_t *db)
{
   int i;
   unsigned char z = (unsigned char) proton_number;
   unsigned char q = (unsigned char) ion_charge;

   if (NULL == db || NULL == db->line
       || proton_number < 1 || proton_number > ISIS_MAX_PROTON_NUMBER
       || ion_charge < 0 || ion_charge > proton_number)
     return NULL;

   /* use brute force for now -- if too slow, use fact that the
    * lines are ordered by element and ion (?) */

   for (i=0; i<db->nlines; i++)
     {
        DB_line_t *line = &db->line[i];
        if (line->proton_number == z
            && line->ion_charge == q
            && line->upper_level == upper_level
            && line->lower_level == lower_level)
          return line;
     }

   return NULL;
}

/*}}}*/

/*}}}*/

/*{{{ database initialization and shutdown functions */

static DB_t *new_db (void) /*{{{*/
{
   DB_t *db = NULL;

   if (NULL == (db = (DB_t *) ISIS_MALLOC (sizeof(DB_t))))
     return NULL;
   memset ((char *)db, 0, sizeof(*db));

   db->line = NULL;

   db->sorted_line_index = NULL;
   db->line_group_table_head = NULL;

   return db;
}

/*}}}*/

static float *Sort_Ptr = NULL;    /* used by the sort routine */

static int compare_wavelen (const void *a, const void *b) /*{{{*/
{
   const int *i = (const int *) a;
   const int *k = (const int *) b;
   float wla = Sort_Ptr[*i];
   float wlb = Sort_Ptr[*k];
   if (wla > wlb) return +1;
   if (wla < wlb) return -1;
   return 0;
}

/*}}}*/

int DB_sort_line_list (DB_t *db) /*{{{*/
{
   int nlines = db->nlines;
   DB_line_t *lines = db->line;
   int *index_array = NULL;
   float *lambdas = NULL;
   int k;

   if (nlines <= 0)
     return -1;

   if ((NULL == (index_array = (int *) ISIS_MALLOC (nlines * sizeof(int))))
       || (NULL == (lambdas = (float *) ISIS_MALLOC (nlines * sizeof(float)))))
     {
        ISIS_FREE(index_array);
        ISIS_FREE(lambdas);
        return -1;
     }

   for (k=0; k < nlines; k++)
     {
        index_array[k] = k;
        lambdas[k] = lines[k].wavelen;
     }

   Sort_Ptr = lambdas;
   qsort (index_array, nlines, sizeof(int), compare_wavelen);
   Sort_Ptr = NULL;

   ISIS_FREE(lambdas);

   ISIS_FREE (db->sorted_line_index);
   db->sorted_line_index = index_array;

   return 0;
}

/*}}}*/

static int build_internals (DB_t *db) /*{{{*/
{
   if ((NULL == db) || isis_user_break())
     return -1;

   if ((0 == DB_sort_line_list (db))
       && (db->sorted_line_index != NULL && db->line != NULL))
     {
        int imin = db->sorted_line_index[0];
        int imax = db->sorted_line_index[db->nlines-1];
        float lam_min = db->line[imin].wavelen;
        float lam_max = db->line[imax].wavelen;
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__,
                    "Tables have %d lines between %0.4g and %0.4g Angstrom",
                    db->nlines, lam_min, lam_max);
     }

   if (db->nlines > 0)
     {
        if (-1 == build_hash_table (db))
          return -1;
     }

   if (-1 == get_branching_ratios (db))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "computing branching ratios");

   return 0;
}

/*}}}*/

DB_t *DB_start (char **energy_level_files, char **wavelength_files) /*{{{*/
{
   DB_t *db = NULL;

   if (NULL == (db = new_db()))
     return NULL;

   if (-1 == load_energy_levels (db, energy_level_files))
     isis_vmesg (WARN, I_READ_FAILED, __FILE__, __LINE__, "energy levels");

   if (-1 == load_line_list (db, wavelength_files))
     isis_vmesg (WARN, I_READ_FAILED, __FILE__, __LINE__, "wavelength tables");

   (void) build_internals (db);

   if (isis_user_break())
     {
        DB_end (db);
        db = NULL;
     }

   return db;
}

/*}}}*/

void DB_end (DB_t *db) /*{{{*/
{
   int i;

   if (db == NULL)
     return;

   free_line_group_table (db->line_group_table_head);

   ISIS_FREE (db->sorted_line_index);
   ISIS_FREE (db->line);
   ISIS_FREE (db->hash_table);

   for (i = 0; i < IONSIZE; i++)
     {
        free_ion (db->ion[i]);
     }

   ISIS_FREE (db);
}

/*}}}*/

/*}}}*/

