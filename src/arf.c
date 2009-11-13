/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2009  Massachusetts Institute of Technology

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

/* $Id: arf.c,v 1.12 2004/02/23 16:35:45 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "cfits.h"
#include "util.h"
#include "keyword.h"
#include "arf.h"
#include "errors.h"

/*}}}*/

/* Effective area curves (ARFs) */

/* array must be NULL terminated */
static const char *Arf_Hdu_Names[] = {"SPECRESP", NULL};
static const char *Arf_Hdu_Names_Hook = "_nonstandard_arf_hdu_names";

/*{{{ new/free */

void Arf_free_arf (Isis_Arf_t *a) /*{{{*/
{
   if (NULL == a)
     return;

   /* Free memory only if nobody is pointing to it */

   if (a->ref_count > 0)
     return;

   ISIS_FREE (a->bin_lo);
   ISIS_FREE (a->bin_hi);
   ISIS_FREE (a->arf);
   ISIS_FREE (a->arf_err);

   if (a->fracexpo_is_vector)
     ISIS_FREE (a->fracexpo.v);

   ISIS_FREE(a->file);
   ISIS_FREE(a);
}
/*}}}*/

static Isis_Arf_t *new_arf (int nbins) /*{{{*/
{
   Isis_Arf_t *a;

   if (NULL == (a = (Isis_Arf_t *) ISIS_MALLOC (sizeof(Isis_Arf_t))))
     return NULL;
   memset ((char *) a, 0, sizeof(*a));

   a->next = NULL;
   a->index = 0;
   a->ref_count = 0;
   a->is_identity = 0;
   a->fracexpo_is_vector = 0;
   a->fracexpo.s = 1.0;

   a->exposure = -1.0;

   if (nbins == 0)        /* allowed for null head node */
     return a;

   if ((NULL == (a->bin_lo = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (a->bin_hi = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (a->arf = (double *) ISIS_MALLOC (nbins * sizeof(double))))
       || (NULL == (a->arf_err = (double *) ISIS_MALLOC (nbins * sizeof(double)))))
     {
        Arf_free_arf (a);
        return NULL;
     }

   return a;
}

/*}}}*/

Isis_Arf_t *Arf_new_arf (int nbins)
{
   return new_arf (nbins);
}

/*}}}*/

/*{{{ linked list maintenance */

static int arf_list_append (Isis_Arf_t *head, Isis_Arf_t *arf) /*{{{*/
{
  Isis_Arf_t *a;

  for (a=head; a->next != NULL; a=a->next)
     ;

  a->next = arf;
  arf->next = NULL;
  arf->index = a->index + 1;

  return arf->index;
}

/*}}}*/

static int delete_after_arf (Isis_Arf_t *a) /*{{{*/
{
   Isis_Arf_t *dead;
   Isis_Arf_t *next;

   if (a == NULL)
     return -1;

   dead = a->next;

   if (dead->ref_count > 0)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "ARF in use");
        return 0;
     }

   isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "Deleting ARF %d", dead->index);

   next = dead->next;
   Arf_free_arf (dead);
   a->next = next;

  return 0;
}

/*}}}*/

Isis_Arf_t *Arf_find_arf_index (Isis_Arf_t *head, int arf_index) /*{{{*/
{
  Isis_Arf_t *a;

   if (head == NULL)
     return NULL;

   for (a=head->next; a != NULL; a = a->next)
     {
        if (a->index == arf_index)
          return a;
     }

   /* isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "ARF %d not found", arf_index); */
   return NULL;
}

/*}}}*/

void Arf_free_arf_list (Isis_Arf_t *head) /*{{{*/
{
   Isis_Arf_t *new_head;

   while (head != NULL)
     {
        new_head = head->next;
        Arf_free_arf (head);
        head = new_head;
     }
}

/*}}}*/

int Arf_delete_arf (Isis_Arf_t *head, int arf_index) /*{{{*/
{
   Isis_Arf_t *a;

   if (head == NULL)
     return -1;

   for (a=head; a->next != NULL; a = a->next)
     {
        if (a->next->index == arf_index)
          return delete_after_arf (a);
     }

   isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "ARF %d not found", arf_index);
   return -1;
}

/*}}}*/

Isis_Arf_t *Arf_init_arf_list (void) /*{{{*/
{
   return new_arf (0);
}

/*}}}*/

int Arf_id_list (Isis_Arf_t *head, unsigned int **ids, unsigned int *num) /*{{{*/
{
   Isis_Arf_t *a;
   unsigned int n;

   *ids = NULL;
   *num = 0;

   if (head == NULL)
     return 0;

   for (a = head->next; a != NULL; a = a->next)
     {
        *num += 1;
     }

   if (*num == 0)
     return 0;

   if (NULL == (*ids = (unsigned int *)ISIS_MALLOC (*num * sizeof (unsigned int))))
     return -1;

   n = 0;
   for (a = head->next; a != NULL; a = a->next)
     {
        (*ids)[n++] = a->index;
     }

   return 0;
}

/*}}}*/

/*}}}*/

static int get_canonical_arf_coordinates (Isis_Arf_t *a, int input_units) /*{{{*/
{
   int reversed;
   int status;

   if (-1 == _isis_fixup_lo_hi_grids (a->bin_lo, a->bin_hi, a->nbins,
                                      &reversed, NULL))
     return -1;

   status = get_canonical_coordinates (a->bin_lo, a->bin_hi,
                                       &a->nbins, input_units);

   if (status == BINS_REVERSED)
     {
        /* get_canonical_coordinates reversed the grid when converting to A.
         * If _isis_fixup_lo_hi_grids reversed the grids, then there is
         * no need to reverse the other arrays.
         */
        if (reversed == 0)
          reversed = 1;
        else
          reversed = 0;
     }

   else if (status != BINS_OK)
     return -1;

   if (reversed == 0)
     return 0;

   if (-1 == reverse_d (a->arf, a->nbins)
       || -1 == reverse_d (a->arf_err, a->nbins))
     return -1;

   if ((a->fracexpo_is_vector && (a->fracexpo.v != NULL))
       && (-1 == reverse_d (a->fracexpo.v, a->nbins)))
     return -1;

   return 0;
}

/*}}}*/

/* FITS format ARF files:  */

#define OFFSETOF(x)  (offsetof(struct Isis_Arf_t, x))

static Keyword_t Arf_Keyword_Table[] = /*{{{*/
{
     { "OBJECT",     OPTIONAL_KEY,   STRING_KEY,   OFFSETOF(object)     },
     { "INSTRUME",   OPTIONAL_KEY,   STRING_KEY,   OFFSETOF(instrument) },
     { "GRATING",    OPTIONAL_KEY,   STRING_KEY,   OFFSETOF(grating)    },
     { "EXPOSURE",   REQUIRED_KEY,   DOUBLE_KEY,   OFFSETOF(exposure)   },
     { "TG_M",       OPTIONAL_KEY,   INT_KEY,      OFFSETOF(order)      },
     { "TG_PART",    OPTIONAL_KEY,   INT_KEY,      OFFSETOF(part)       },
     { "TG_SRCID",   OPTIONAL_KEY,   INT_KEY,      OFFSETOF(srcid)      },
     { NULL,         OPTIONAL_KEY,   NULL_KEY,     0 }
};

/*}}}*/

static int validate_arf (int n, double *arf, double * arf_err) /*{{{*/
{
   int i, reset, num_positive;

   if ((n <= 0) || (arf_err == NULL) || (arf == NULL))
     return -1;

   reset = 0;
   num_positive = 0;

   for (i = 0; i < n; i++)
     {
        if (arf[i] > 0.0)
          num_positive++;

        if (arf_err[i] <= 0.0)
          {
             arf_err[i] = 1.0;
             reset=1;
          }
     }

   if (num_positive == 0)
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "all ARF bins <= 0");
        return -1;
     }

   return 0;
}
/*}}}*/

static int read_fexp_info (Isis_Arf_t *a, cfitsfile *fp) /*{{{*/
{
   static char fexp_name[] = "FRACEXPO";

   if (!cfits_col_exist (fexp_name, fp))
     {
        a->fracexpo_is_vector = 0;
        if (-1 == cfits_read_double_keyword (&a->fracexpo.s, fexp_name, fp))
          a->fracexpo.s = 1.0;
        return 0;
     }

   a->fracexpo.v = (double *) ISIS_MALLOC (a->nbins * sizeof(double));
   if (NULL == a->fracexpo.v)
     return -1;

   a->fracexpo_is_vector = 1;

   if (-1 == cfits_read_double_col (a->fracexpo.v, a->nbins, 1, fexp_name, fp))
     return -1;

   return 0;
}

/*}}}*/

static int validate_and_append (Isis_Arf_t *head, Isis_Arf_t *a, int input_units) /*{{{*/
{
   if (-1 == get_canonical_arf_coordinates (a, input_units))
     return -1;

   if (-1 == validate_arf (a->nbins, a->arf, a->arf_err))
     return -1;

   return arf_list_append (head, a);
}

/*}}}*/

int Arf_read_arf (Isis_Arf_t *head, char *filename) /*{{{*/
{
   Keyword_t *keytable = Arf_Keyword_Table;
   Isis_Arf_t *a = NULL;
   cfitsfile *fp = NULL;
   char bin_units[ISIS_ARF_VALUE_SIZE];
   int nbins, input_units;
   int ret = -1;
   int id = -1;

   if (head == NULL || filename == NULL)
     return -1;

   if (NULL == (fp = cfits_open_file_readonly (filename)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", filename);
        return -1;
     }

   if (-1 == cfits_move_to_matching_hdu (fp, Arf_Hdu_Names, Arf_Hdu_Names_Hook, NULL))
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__,
                    "No recognized response HDU found in %s", filename);
        goto finish;
     }

   if (-1 == cfits_get_colunits (bin_units, "ENERG_LO", fp))
     input_units = U_KEV;
   else if (-1 == (input_units = unit_id (bin_units)))
     goto finish;

   if (-1 == cfits_read_int_keyword(&nbins, "NAXIS2", fp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "NAXIS2 => %s", filename);
        goto finish;
     }

   a = new_arf (nbins);
   if (NULL == a)
     goto finish;

   a->nbins = nbins;

   if (NULL == (a->file = isis_make_string(filename)))
     goto finish;

   (void) Key_read_header_keywords (fp, (char *)a, keytable, FITS_FILE);

   if (-1 == read_fexp_info (a, fp)
       || -1 == cfits_read_double_col (a->bin_lo, nbins, 1, "ENERG_LO", fp)
       || -1 == cfits_read_double_col (a->bin_hi, nbins, 1, "ENERG_HI", fp)
       || -1 == cfits_read_double_col (a->arf, nbins, 1, "SPECRESP", fp)
       || -1 == cfits_read_optional_double_col (a->arf_err, nbins, 1, "RESP_ERR", fp)
       )
     {
        isis_vmesg (FAIL, I_READ_COL_FAILED, __FILE__, __LINE__, "%s", filename);
        goto finish;
     }

   if (-1 == (id = validate_and_append (head, a, input_units)))
     goto finish;

   ret = 0;

   finish:

   (void) cfits_close_file (fp);

   if (ret)
     {
        Arf_free_arf (a);
        return ret;
     }

   return id;
}

/*}}}*/

int Arf_define_arf (Isis_Arf_t *head, unsigned int nbins, /*{{{*/
                    double *bin_lo, double *bin_hi, double *arf, double *arf_err)
{
   Isis_Arf_t *a;
   int id;

   if (head == NULL)
     return -1;

   if (NULL == (a = new_arf (nbins)))
     return -1;

   a->nbins = nbins;
   a->file = NULL;
   a->object[0] = 0;
   a->instrument[0] = 0;
   a->grating[0] = 0;
   a->order = 0;
   a->part = 0;
   a->srcid = 0;

   memcpy ((char *)a->bin_lo, (char *)bin_lo, nbins * sizeof(double));
   memcpy ((char *)a->bin_hi, (char *)bin_hi, nbins * sizeof(double));
   memcpy ((char *)a->arf, (char *)arf, nbins * sizeof(double));
   if (arf_err != NULL)
     memcpy ((char *)a->arf_err, (char *)arf_err, nbins * sizeof(double));
   else
     memset ((char *)a->arf_err, 0, nbins * sizeof(double));

   if (-1 == (id = validate_and_append (head, a, U_ANGSTROM)))
     {
        Arf_free_arf (a);
        return -1;
     }

   return id;
}

/*}}}*/

int Arf_append_arf (Isis_Arf_t *head, Isis_Arf_t *a) /*{{{*/
{
   return validate_and_append (head, a, U_ANGSTROM);
}

/*}}}*/

/* ARF get/put/set */

int Arf_arf_size (Isis_Arf_t *a) /*{{{*/
{
   if (NULL == a)
     return -1;

   return a->nbins;
}

/*}}}*/

int Arf_get_arf (Isis_Arf_t *a, double *arf, double *arf_err, double *binlo, double *binhi) /*{{{*/
{
   int size;

   if (NULL == a)
     return -1;

   size = a->nbins * sizeof(double);

   memcpy ((char *) binlo, (char *) a->bin_lo, size);
   memcpy ((char *) binhi, (char *) a->bin_hi, size);
   memcpy ((char *) arf, (char *) a->arf, size);
   memcpy ((char *) arf_err, (char *) a->arf_err, size);

   return 0;
}

/*}}}*/

int Arf_put_arf (Isis_Arf_t *a, double *arf, double *arf_err, double *binlo, double *binhi) /*{{{*/
{
   int i, size;

   if (NULL == a)
     return -1;

   if (-1 == validate_wavelength_grid (a->nbins, binlo, binhi))
     return -1;

   if (-1 == validate_arf (a->nbins, arf, arf_err))
     return -1;

   for (i = 0; i < a->nbins; i++)
     {
        if ( (fabs(1.0 - a->bin_lo[i]/binlo[i]) > GRID_TOL)
            || (fabs(1.0 - a->bin_hi[i]/binhi[i]) > GRID_TOL) )
          break;
     }

   if (i != a->nbins)
     isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "ARF grid change may cause ARF/RMF mismatch");

   size = a->nbins * sizeof(double);
   memcpy ((char *) a->bin_lo, (char *)binlo, size);
   memcpy ((char *) a->bin_hi, (char *)binhi, size);
   memcpy ((char *) a->arf, (char *) arf, size);
   memcpy ((char *) a->arf_err, (char *) arf_err, size);

   return 0;
}

/*}}}*/

int Arf_get_arf_exposure (Isis_Arf_t *a, double * exposure) /*{{{*/
{
   if (NULL == a)
     return -1;

   *exposure = a->exposure;

   return 0;
}

/*}}}*/

int Arf_set_arf_exposure (Isis_Arf_t *a, double exposure)  /*{{{*/
{
   if (NULL == a)
     return -1;

   a->exposure = exposure;

   return 0;
}

/*}}}*/

void Arf_free_arf_info (Arf_Info_Type *ai) /*{{{*/
{
   if (ai == NULL)
     return;

   ISIS_FREE (ai->object);
   ISIS_FREE (ai->instrument);
   ISIS_FREE (ai->grating);
   ISIS_FREE (ai->file);

   if (ai->fracexpo_is_vector)
     {
        ISIS_FREE (ai->fracexpo.v);
     }
}

/*}}}*/

int Arf_set_arf_info (Isis_Arf_t *a, Arf_Info_Type *ai) /*{{{*/
{
   if (a == NULL || ai == NULL)
     return -1;

   if (0 == ai->fracexpo_is_vector)
     {
        if (a->fracexpo_is_vector)
          ISIS_FREE(a->fracexpo.v);
        a->fracexpo.s = ai->fracexpo.s;
     }
   else
     {
        int size;

        if (ai->nbins != a->nbins)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "array size mismatch; fracexpo size %d; arf size %d\n",
                         ai->nbins, a->nbins);
             return -1;
          }

        if (a->fracexpo_is_vector)
          ISIS_FREE(a->fracexpo.v);

        size = ai->nbins * sizeof(double);
        if (NULL == (a->fracexpo.v = (double *) ISIS_MALLOC (size)))
          return -1;
        memcpy ((char *)a->fracexpo.v, (char *)ai->fracexpo.v, size);
     }

   a->fracexpo_is_vector = ai->fracexpo_is_vector;

   if (ai->file)
     {
        char *s;
        if (NULL == (s = isis_make_string (ai->file)))
          return -1;
        ISIS_FREE (a->file);
        a->file = s;
     }

   isis_strcpy(a->object,ai->object, sizeof (a->object));
   isis_strcpy(a->instrument,ai->instrument, sizeof(a->instrument));
   isis_strcpy(a->grating,ai->grating, sizeof(a->grating));

   a->exposure = ai->exposure;
   a->order = ai->order;
   a->part = ai->part;
   a->srcid = ai->srcid;

   return 0;
}

/*}}}*/

int Arf_get_arf_info (Isis_Arf_t *a, Arf_Info_Type *ai) /*{{{*/
{
   if (a == NULL || ai == NULL)
     return -1;

   if (a->fracexpo_is_vector)
     {
        int size = a->nbins * sizeof(double);
        if (NULL == (ai->fracexpo.v = (double *) ISIS_MALLOC (size)))
          return -1;
        memcpy ((char *)ai->fracexpo.v, (char *)a->fracexpo.v, size);
        ai->nbins = a->nbins;
     }
   else
     {
        ai->fracexpo.s = a->fracexpo.s;
        ai->nbins = 1;
     }

   ai->fracexpo_is_vector = a->fracexpo_is_vector;

   if (a->file)
     {
        char *s;
        if (NULL == (s = isis_make_string (a->file)))
          return -1;
        ai->file = s;
     }

   if (NULL == (ai->object = isis_make_string (a->object)))
     return -1;
   if (NULL == (ai->instrument = isis_make_string (a->instrument)))
     return -1;
   if (NULL == (ai->grating = isis_make_string (a->grating)))
     return -1;

   ai->exposure = a->exposure;
   ai->order = a->order;
   ai->part = a->part;
   ai->srcid = a->srcid;

   return 0;
}

/*}}}*/

/* Identity ARF */

Isis_Arf_t *Arf_make_identity_arf (double *bin_lo, double * bin_hi, unsigned int nbins) /*{{{*/
{
   Isis_Arf_t *a;
   double *p;

   if (NULL == (a = new_arf (nbins)))
     return NULL;

   a->is_identity = 1;
   a->exposure = 1.0;
   a->nbins = nbins;

   memcpy ((char *)a->bin_lo, (char *)bin_lo, nbins * sizeof(double));
   memcpy ((char *)a->bin_hi, (char *)bin_hi, nbins * sizeof(double));

   p = a->arf;
   while (p < (a->arf + nbins)) *p++ = 1.0;

   p = a->arf_err;
   while (p < (a->arf_err + nbins)) *p++ = 1.e-6;

   return a;
}

int Arf_is_identity (Isis_Arf_t *a)
{
   if (a == NULL)
     return 0;

   return a->is_identity;
}

/*}}}*/

