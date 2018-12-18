/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2018  Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Authors:  John C. Houck  <houck@space.mit.edu>
              John E. Davis <davis@space.mit.edu>

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

/* $Id: rmf_file.c,v 1.18 2004/06/06 02:42:22 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "cfits.h"
#include "util.h"
#include "rmf.h"
#include "errors.h"

/*}}}*/

int Isis_Rmf_OGIP_Compliance = 2;

/*
 * The input RMF is assumed to have
 * 1) the ARF grid in order of increasing energy
 * 2) the response vectors given in order of increasing
 *    channel number.
 *
 * Both RMF grids are stored in energy units [keV] and are
 * sorted in order of increasing energy.
 *
 * Internally, the RMF has the detector channels
 * numbered in order of increasing energy.  If the channels are
 * numbered backward in the input file, theywill be
 * re-labeled on input and the FCHAN values in the
 * matrix will be updated accordingly.
 *
 * rmf->swapped_channels != 0 indicates that the channels
 * got re-labeled.
 */

/*{{{ struct definitions */

typedef struct Rmf_Vector_t Rmf_Vector_t;
typedef struct Rmf_File_t Rmf_File_t;
typedef struct Rmf_Element_t Rmf_Element_t;

struct Rmf_Element_t
{
   unsigned int first_channel;
   unsigned int num_channels;
   float *response;             /* in order of increasing channel number */
};

struct Rmf_Vector_t
{
   unsigned int num_grps;       /* number of detector channel groups */
   Rmf_Element_t *elem;         /* array of detector channel group responses */
};

typedef struct
{
   double threshold;
   union
     {
        char *file;
        SLang_Name_Type *funct;
     }
   f;
   int type;
#define RMF_TYPE_FILE        1
#define RMF_TYPE_SLANG        2
   int swapped_channels;
   int energy_ordered_ebounds;
   int is_initialized;                       /* set but not unused */
   int strict;
   char *matrix_extname;
   char *ebounds_extname;
   Isis_Rmf_Grid_Type *arf;      /* keV, increasing order */
   Isis_Rmf_Grid_Type *ebounds;  /* keV, increasing order */
   Rmf_Vector_t *v;              /* keV, increasing order */
   unsigned int num_ebins;
   int offset;                   /* F_CHAN TLMIN value */
}
Rmf_Client_Data_t;

struct Rmf_File_t
{
   cfitsfile *ft;
   long num_rows;

   int includes_effective_area;
   int offset;

   int energ_lo_col;
   int energ_hi_col;
   int n_grp_col;
   int f_chan_col;
   int n_chan_col;
   int matrix_col;

   /* If the column values are constant, they may be moved to keyword values
    * to save space.  In that case, the _col fields of this structure will be
    * set to -1 and the actual value is specified below.
    */
   long n_grp_val;
   long f_chan_val;
   long n_chan_val;
};

/*}}}*/

static int print_options_help (Rmf_Client_Data_t *cd);
static int read_order_keyword (cfitsfile *ft, int *order);

static Rmf_Client_Data_t *get_client_data (Isis_Rmf_t *rmf) /*{{{*/
{
   return (Rmf_Client_Data_t *)rmf->client_data;
}

/*}}}*/

static void free_rmf_vector (Rmf_Vector_t *v) /*{{{*/
{
   Rmf_Element_t *elem;
   Rmf_Element_t *elem_max;

   if ((v == NULL) || ((elem = v->elem) == NULL))
     return;

   elem_max = elem + v->num_grps;

   while (elem < elem_max)
     {
        if (elem->response != NULL) ISIS_FREE (elem->response);
        elem++;
     }

   ISIS_FREE (v->elem);
   v->num_grps = 0;
}

/*}}}*/

/* RMF input */

static int check_rmf_extension (cfitsfile *ft) /*{{{*/
{
   char value [CFLEN_KEYWORD];

   if ((-1 == cfits_read_string_keyword (value, "HDUCLAS2", ft))
       || strcmp (value, "RSP_MATRIX"))
     return -1;

   if (-1 == cfits_read_string_keyword (value, "HDUCLAS3", ft))
     return -1;

   if (!strcmp (value, "REDIST"))
     return 1;

   if (!strcmp (value, "DETECTOR"))
     return 2;

   if (!strcmp (value, "FULL"))
     return 3;

   return -1;
}

/*}}}*/

static Rmf_File_t *open_rmf_file (char *file, Isis_Rmf_t *rmf) /*{{{*/
{
   /* This routine is complicated by the fact that CAL/GEN/92-002 allows
    * columns to be specified as a keyword value.
    */
   static const char *names[6] =
     {
        "ENERG_LO", "ENERG_HI", "N_GRP", "F_CHAN", "N_CHAN", "MATRIX"
     };
   static const char *ext_list [] =
     {
        "SPECRESP MATRIX", "MATRIX", NULL
     };
   cfitsfile *ft = NULL;
   int columns[6];
   long allow_as_keyword [6];
   Rmf_File_t *rft = NULL;
   Rmf_Client_Data_t *cd = NULL;
   unsigned int i;
   int arf_status;

   if (NULL == (rft = (Rmf_File_t *) ISIS_MALLOC (sizeof (*rft))))
     return NULL;
   memset ((char *) rft, 0, sizeof (*rft));

   if (NULL == (ft = cfits_open_file_readonly (file)))
     goto return_error;

   /* XMM RMFs keep the order keyword *only* in the primary header */
   (void) read_order_keyword (ft, &rmf->order);

   cd = (Rmf_Client_Data_t *) rmf->client_data;

   if (cd->matrix_extname)
     {
        if (-1 == cfits_movnam_hdu (ft, cd->matrix_extname))
          {
             isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__,
                         cd->matrix_extname);
             goto return_error;
          }
     }
   else if (-1 == cfits_move_to_matching_hdu (ft, ext_list, "_nonstandard_rmf_matrix_hdu_names",
                                              (cd->strict ? check_rmf_extension : NULL)))
     {
        isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__, "RMF Matrix");
        if (cd->strict)
          {
             isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "==> You might try loading: \"%s;strict=0\".", file);
             isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "    or set Rmf_OGIP_Compliance=0.", file);
          }
        (void) print_options_help (cd);
        goto return_error;
     }

   arf_status = check_rmf_extension (ft);
   if (arf_status == 3)
     {
        isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "RMF includes the effective area");
        rft->includes_effective_area = 1;
        if (cd->strict != 2)
          goto return_error;
     }
   else if (arf_status < 0)
     {
        if (cd->strict)
          isis_vmesg (WARN, I_INFO, __FILE__, __LINE__,
                      "**** WARNING:  RMF has non-standard format -- important information is missing.");
     }

   if (-1 == cfits_read_long_keyword (&rft->num_rows, "NAXIS2", ft))
     goto return_error;

   if (-1 == cfits_read_int_colkey (&rft->offset, "TLMIN", "F_CHAN", ft))
     rft->offset = -1;

   if (-1 == cfits_get_column_numbers (ft, 6, names, columns))
     goto return_error;

   /* This is ugly !!!!!!!! */
   allow_as_keyword [0] = 0;
   allow_as_keyword [1] = 0;
   allow_as_keyword [5] = 0;
   allow_as_keyword [2] = 1;
   allow_as_keyword [3] = 1;
   allow_as_keyword [4] = 1;

   for (i = 0; i < 6; i++)
     {
        if (columns[i] != -1) continue;
        if ((allow_as_keyword [i] == 0)
            || (-1 == cfits_read_long_keyword (&allow_as_keyword[i], names[i],ft)))
          {
             isis_vmesg (FAIL, I_COL_NOT_FOUND, __FILE__, __LINE__, "%s", names[i]);
             goto return_error;
          }
     }

   rft->ft = ft;

   rft->energ_lo_col = columns[0];
   rft->energ_hi_col = columns[1];
   rft->matrix_col = columns[5];

   if (-1 == (rft->n_grp_col = columns[2]))
     rft->n_grp_val = allow_as_keyword[2];

   if (-1 == (rft->f_chan_col = columns[3]))
     rft->f_chan_val = allow_as_keyword[3];

   if (-1 == (rft->n_chan_col = columns[4]))
     rft->n_chan_val = allow_as_keyword[4];

   return rft;

   return_error:

   ISIS_FREE (rft);
   cfits_close_file (ft);

   return NULL;
}

/*}}}*/

static void close_rmf_file (Rmf_File_t **rft) /*{{{*/
{
   if (rft == NULL) return;
   cfits_close_file ((*rft)->ft);
   ISIS_FREE (*rft);
}

/*}}}*/

static int read_rmf_vector (Rmf_File_t *rft, int row, Rmf_Vector_t *v, /*{{{*/
                            double *elo, double *ehi,
                            int chan_range[2])
{
   cfitsfile *ft;
   unsigned int *nchan=NULL, *fchan=NULL;
   Rmf_Element_t *elem;
   Rmf_Element_t *elem_max;
   float f_elo, f_ehi;
   unsigned int i;
   unsigned int ngrps;
   long offset;
   int matrix_col;

   if (rft == NULL || v == NULL)
     return -1;

   v->num_grps = 0;
   v->elem = NULL;

   ft = rft->ft;

   if (rft->n_grp_col == -1)
     ngrps = rft->n_grp_val;
   else if (-1 == cfits_read_column_uints (ft, rft->n_grp_col, row, 1, &ngrps, 1))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading number of RMF groups");
        return -1;
     }

   if (-1 == cfits_read_column_floats (ft, rft->energ_lo_col, row, 1, &f_elo, 1)
       || -1 == cfits_read_column_floats (ft, rft->energ_hi_col, row, 1, &f_ehi, 1))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading RMF vector energy grid");
        return -1;
     }

   *elo = (double) f_elo;
   *ehi = (double) f_ehi;

   if ((*elo == FLT_MIN) || (*ehi == FLT_MIN))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Corrupted energy grid in RMF");
        return -1;
     }

   if (ngrps == 0)
     return 0;

   if (NULL == (elem = (Rmf_Element_t *) ISIS_MALLOC (ngrps * sizeof(Rmf_Element_t))))
     return -1;
   memset ((char *) elem, 0, ngrps * sizeof(*elem));

   if (NULL == (nchan = (unsigned int *) ISIS_MALLOC (2*ngrps * sizeof(unsigned int))))
     {
        ISIS_FREE(elem);
        return -1;
     }
   fchan = nchan + ngrps;

   v->num_grps = ngrps;
   v->elem = elem;

   /* Hopefully the cfitsio library will cache things. In general things will
    * be in the heap and whether it has a cache or not may be irrelevant.
    */

   if (-1 == rft->f_chan_col)
     {
        unsigned int imax = ngrps;
        unsigned int f_chan_val = rft->f_chan_val;

        for (i = 0; i < imax; i++) fchan [i] = f_chan_val;
     }
   else if (-1 == cfits_read_column_uints (ft, rft->f_chan_col, row, 1, fchan, ngrps))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading RMF fchan column");
        goto return_error;
     }

   if (-1 == rft->n_chan_col)
     {
        unsigned int imax = ngrps;
        unsigned int n_chan_val = rft->n_chan_val;

        for (i = 0; i < imax; i++) nchan [i] = n_chan_val;
     }
   else if (-1 == cfits_read_column_uints (ft, rft->n_chan_col, row, 1, nchan, ngrps))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading RMF nchan column");
        goto return_error;
     }

   elem = v->elem;
   elem_max = elem + ngrps;
   i = 0;

   matrix_col = rft->matrix_col;
   offset = 1;

   while (elem < elem_max)
     {
        elem->first_channel = fchan[i];
        elem->num_channels = nchan[i];

        /* derive min/max channels included in this mapping */
        if ((int) fchan[i] < chan_range[0])
          {
             chan_range[0] = (int)fchan[i];
          }
        if (chan_range[1] < (int)fchan[i] + (int)nchan[i] - 1)
          {
             chan_range[1] = (int)fchan[i] + (int)nchan[i] - 1;
          }

        elem->response = (float *) ISIS_MALLOC (elem->num_channels * sizeof(float));
        if ((elem->response == NULL)
            || (-1 == cfits_read_column_floats (ft, matrix_col, row, offset,
                                                elem->response, elem->num_channels)))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading RMF response column");
             goto return_error;
          }

        offset += elem->num_channels;
        i++;
        elem++;
     }

   ISIS_FREE(nchan);
   return 0;

   return_error:
   ISIS_FREE(nchan);
   free_rmf_vector (v);
   return -1;
}

/*}}}*/

/* returns energy grid [keV] sorted in increasing order */
static int read_rmf_ebounds (cfitsfile *rmf_fp, int chan_range[2], Isis_Rmf_Grid_Type **gp, /*{{{*/
                             int *energy_ordered_ebounds, int *swap_channels,
                             int *min_chan, Rmf_Client_Data_t *cd)
{
   static const char *ext_names[] = {"EBOUNDS", NULL};
   int num_channels, max_chan;
   int *channel = NULL;
   Isis_Rmf_Grid_Type *g = NULL;
   int swapped;
   unsigned int i;

   *gp = NULL;

   if (cd->ebounds_extname)
     {
        if (-1 == cfits_movnam_hdu (rmf_fp, cd->ebounds_extname))
          {
             isis_vmesg (FAIL, I_HDU_NOT_FOUND, __FILE__, __LINE__,
                         cd->ebounds_extname);
             return -1;
          }
     }
   else if (-1 == cfits_move_to_matching_hdu (rmf_fp, ext_names, "_nonstandard_rmf_ebounds_hdu_names", NULL))
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__,
                    "RMF file contains no recognized ebounds HDU");
        (void) print_options_help (cd);
        return -1;
     }

   if (-1 == cfits_read_int_keyword (&num_channels, "NAXIS2", rmf_fp))
     {
        isis_vmesg (FAIL, I_READ_KEY_FAILED, __FILE__, __LINE__, "NAXIS2 from EBOUNDS extension");
        return -1;
     }

   *gp = Isis_new_rmf_grid (num_channels, NULL, NULL);
   if (*gp == NULL) return -1;
   g = *gp;

   if (NULL == (channel = (int *) ISIS_MALLOC (g->nbins * sizeof(int))))
     return -1;

   /* FIXME -- get units from header */
   g->units = U_KEV;

   if ((-1 == cfits_read_int_col (channel, g->nbins, 1, "CHANNEL", rmf_fp))
       || (-1 == cfits_read_double_col (g->bin_lo, g->nbins, 1, "E_MIN", rmf_fp))
       || (-1 == cfits_read_double_col (g->bin_hi, g->nbins, 1, "E_MAX", rmf_fp)))
     {
        isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "RMF EBOUNDS extension");
        ISIS_FREE (channel);
        return -1;
     }

   *min_chan = INT_MAX;
   max_chan = channel[0];
   for (i = 0; i < g->nbins; i++)
     {
        if (channel[i] < *min_chan)
          *min_chan = channel[i];
        else if (channel[i] > max_chan)
          max_chan = channel[i];
     }

   if ((chan_range[0] < *min_chan) || (max_chan < chan_range[1]))
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "RMF matrix maps to detector channels outside EBOUNDS grid");
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "matrix spans channels: %d to %d   ebounds spans channels: %d to %d",
                    chan_range[0], chan_range[1], *min_chan, max_chan);
     }

   /* Try to fix common sloppiness */
   if (g->bin_lo[0] < 0)
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "negative EBOUNDS value E_MIN=%g, set to zero",
                    g->bin_lo[0]);
        g->bin_lo[0] = 0.0;
     }

   if (-1 == _isis_fixup_lo_hi_grids (g->bin_lo, g->bin_hi, g->nbins, &swapped, NULL))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Invalid EBOUNDS grid in RMF");
        ISIS_FREE (channel);
        return -1;
     }

   if (swapped)
     {
        *energy_ordered_ebounds = 0;
        *swap_channels = (channel[0] < channel[1]);
     }
   else
     {
        *energy_ordered_ebounds = 1;
        *swap_channels = (channel[0] > channel[1]);
     }
   ISIS_FREE (channel);
   return 0;
}

/*}}}*/

static int renumber_detector_channels (Isis_Rmf_t *rmf) /*{{{*/
{
   Rmf_Client_Data_t *cd = get_client_data (rmf);
   unsigned int i, num_channels, swapped;

   if ((rmf == NULL) || (cd == NULL))
     return -1;

   num_channels = cd->ebounds->nbins;
   swapped = cd->swapped_channels;

   for (i = 0; i < cd->num_ebins; i++)
     {
        Rmf_Vector_t *v = &cd->v[i];
        unsigned int g;

        for (g = 0; g < v->num_grps; g++)
          {
             Rmf_Element_t *elem = &v->elem[g];
             elem->first_channel -= cd->offset;
             if (swapped)
               {
                  unsigned int last_chan = elem->first_channel + elem->num_channels - 1;
                  elem->first_channel = num_channels - last_chan - 1;
                  reverse_f (elem->response, elem->num_channels);
               }
          }
     }

   return 0;
}

/*}}}*/

static int validate_rmf (Isis_Rmf_t *rmf) /*{{{*/
{
   unsigned int e, g, num_ok, ignored;
   Rmf_Client_Data_t *cd = get_client_data (rmf);

   ignored = 0;
   num_ok = 0;
   for (e = 0; e < cd->num_ebins; e++)
     {
        Rmf_Vector_t *v = &cd->v[e];
        for (g = 0; g < v->num_grps; g++)
          {
             Rmf_Element_t *elem = &v->elem[g];
             int first = elem->first_channel;

             if (first < 0)
               {
                  ignored++;
                  isis_vmesg (INFO, I_WARNING, __FILE__, __LINE__,
                              "invalid RMF chan: ebin=%d group=%d fchan=%d offset=%d",
                              e, g,
                              elem->first_channel + cd->offset,
                              cd->offset);
                  elem->num_channels = 0;
               }

             if (cd->ebounds->nbins < (first + elem->num_channels))
               {
                  ignored++;
                  isis_vmesg (INFO, I_WARNING, __FILE__, __LINE__,
                              "invalid RMF chan: ebin=%d group=%d fchan=%d nchan=%d  num_ebounds=%d",
                              e, g,
                              elem->first_channel,
                              elem->num_channels,
                              cd->ebounds->nbins);
                  elem->num_channels = 0;
               }

             num_ok += elem->num_channels;
          }
     }

   if (ignored)
     isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "Ignored %u invalid RMF channels", ignored);

   if (num_ok == 0)
     isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "RMF has no valid output channels");

   return num_ok ? 0 : -1;
}

/*}}}*/

static int read_order_keyword (cfitsfile *ft, int *order) /*{{{*/
{
   static const char *keys[] =
     {
        "ORDER", "TG_M", "RFLORDER", NULL
     };
   const char **k;
   int ret = 0;

   for (k = keys; *k != NULL; k++)
     {
        if (0 == (ret = cfits_read_int_keyword (order, *k, ft)))
          break;
     }

   return ret;
}

/*}}}*/

static int read_rmf (Isis_Rmf_t *rmf, char *file) /*{{{*/
{
   Rmf_Client_Data_t *cd;
   Rmf_File_t *rft = NULL;
   Isis_Rmf_Grid_Type *g = NULL;
   int row, ret = -1;
   int reversed, min_chan, rmf_order;
   int chan_range[2];

   if ((file == NULL) || (rmf == NULL))
     return -1;

   if (!isis_strcasecmp (file, "NONE"))
     return -1;

   cd = (Rmf_Client_Data_t *) rmf->client_data;

   if (NULL == (rft = open_rmf_file (file, rmf)))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", file);
        return -1;
     }
   rmf->includes_effective_area = rft->includes_effective_area;

   cd->num_ebins = rft->num_rows;

   if (-1 == cfits_read_string_keyword (rmf->grating, "GRATING", rft->ft))
     isis_vmesg (INFO, I_KEY_NOT_FOUND, __FILE__, __LINE__, "GRATING");

   if (-1 == cfits_read_string_keyword (rmf->instrument, "DETNAM", rft->ft))
     isis_vmesg (INFO, I_KEY_NOT_FOUND, __FILE__, __LINE__, "DETNAM");

   if (0 == read_order_keyword (rft->ft, &rmf_order))
     rmf->order = rmf_order;

   if (-1 == cfits_read_double_keyword (&cd->threshold, "LO_THRES", rft->ft))
     cd->threshold = 1.e-6;

   if (NULL == (cd->arf = Isis_new_rmf_grid (cd->num_ebins, NULL, NULL)))
     goto finish;

   /* FIXME - get arf_grid units from header */
   cd->arf->units = U_KEV;

   if (NULL == (cd->v = (Rmf_Vector_t *) ISIS_MALLOC (cd->num_ebins * sizeof(Rmf_Vector_t))))
     goto finish;
   memset ((char *)cd->v, 0, cd->num_ebins * sizeof (Rmf_Vector_t));

   chan_range[0] = INT_MAX;
   chan_range[1] = -INT_MAX;

   g = cd->arf;
   for (row = 0; row < rft->num_rows; row++)
     {
        if (-1 == read_rmf_vector (rft, row+1, &cd->v[row],
                                   &g->bin_lo[row], &g->bin_hi[row], chan_range))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading RMF vector in row %d", row);
             goto finish;
          }
     }

   /* Try to fix common sloppiness */
   if (g->bin_lo[0] < 0)
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "negative ENERG_LO=%g set to zero",
                    g->bin_lo[0]);
        g->bin_lo[0] = 0.0;
     }

   if (-1 == _isis_fixup_lo_hi_grids (g->bin_lo, g->bin_hi, g->nbins, &reversed, 0))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Invalid grid in RMF");
        goto finish;
     }

   if (reversed ||
       ((g->nbins > 1) && (g->bin_lo[0] > g->bin_lo[1])))
     {
        isis_vmesg (FAIL, I_UNSUPPORTED_FORMAT, __FILE__, __LINE__, "RMF energies in decreasing order");
        goto finish;
     }

   if (-1 == read_rmf_ebounds (rft->ft, chan_range, &cd->ebounds, &cd->energy_ordered_ebounds, &cd->swapped_channels,
                               &min_chan, cd))
     goto finish;

   /* channels might have been numbered starting at zero or one */
   if (rft->offset >= 0)
     cd->offset = rft->offset;
   else
     cd->offset = min_chan;

   if (-1 == renumber_detector_channels (rmf))
     goto finish;

   if (-1 == validate_rmf (rmf))
     goto finish;

   cd->is_initialized = 1;
   ret = 0;

   finish:

   close_rmf_file (&rft);

   return ret;
}

/*}}}*/

/* Method Interface */

static void delete_client_data (Isis_Rmf_t *rmf) /*{{{*/
{
   Rmf_Client_Data_t *cd = get_client_data (rmf);

   if (NULL != cd)
     {
        if (cd->v != NULL)
          {
             unsigned int i;
             for (i = 0; i < cd->num_ebins; i++)
               free_rmf_vector (&cd->v[i]);
             ISIS_FREE (cd->v);
          }
        Isis_free_rmf_grid (cd->arf);
        Isis_free_rmf_grid (cd->ebounds);
        if (cd->type == RMF_TYPE_FILE)
          ISIS_FREE (cd->f.file);
        else if (cd->type == RMF_TYPE_SLANG)
          SLang_free_function (cd->f.funct);

        ISIS_FREE (cd->matrix_extname);
        ISIS_FREE (cd->ebounds_extname);
     }

   ISIS_FREE (rmf->client_data);
}
/*}}}*/

static int redistribute (Isis_Rmf_t *rmf, unsigned int in_lam, double flux, /*{{{*/
                         double *det_chan, unsigned int num_ebounds)
{
   Rmf_Client_Data_t *cd = get_client_data (rmf);
   Rmf_Vector_t *v;
   unsigned int g, in_egy;

   in_egy = cd->num_ebins - in_lam - 1;
   v = &cd->v[in_egy];

   for (g = 0; g < v->num_grps; g++)
     {
        Rmf_Element_t *elem = &v->elem[g];
        float *response = elem->response;
        int k, num_channels = elem->num_channels;
        double *d = det_chan + num_ebounds - elem->first_channel - 1;

        for (k = 0; k < num_channels; k++)
          d[-k] += flux * response[k];
     }

   return 0;
}
/*}}}*/

static int set_noticed_model_bins (Isis_Rmf_t *rmf, int num_chan, int *chan_notice, /*{{{*/
                                   int num_model, int *model_notice)
{
   Rmf_Client_Data_t *cd = get_client_data (rmf);
   unsigned int *undetected_model_bin;
   unsigned int e_model;

   if ((NULL == rmf) || (NULL == cd))
     return -1;

   if (num_chan != (int) cd->ebounds->nbins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "RMF EBOUNDS/data grid mismatch");
        return -1;
     }

   if (num_model != (int) cd->num_ebins)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "RMF EBOUNDS/ARF grid mismatch");
        return -1;
     }

   memset ((char *)model_notice, 0, num_model * sizeof(int));

   undetected_model_bin = (unsigned int *) ISIS_MALLOC (num_model * sizeof(unsigned int));
   if (undetected_model_bin == NULL)
     return -1;
   for (e_model = 0; e_model < (unsigned int) num_model; e_model++)
     undetected_model_bin[e_model] = 1;

   for (e_model = 0; e_model < cd->num_ebins; e_model++)
     {
        Rmf_Vector_t *v = &cd->v[e_model];
        unsigned int g;

        for (g = 0; g < v->num_grps; g++)
          {
             Rmf_Element_t *elem = &v->elem[g];
             int num_channels = elem->num_channels;
             int *noticed = (chan_notice + num_chan - 1) - elem->first_channel;
             int e_ch;

             for (e_ch = 0; e_ch < num_channels; e_ch++)
               {
                  if (elem->response[e_ch] > 0.0)
                    {
                       unsigned int k = num_model - e_model - 1;
                       undetected_model_bin[k] = 0;
                       if (noticed[-e_ch])
                         model_notice[k] = 1;
                    }
               }
          }
     }

   for (e_model = 0; e_model < (unsigned int) num_model; e_model++)
     {
        if (undetected_model_bin[e_model])
          model_notice[e_model] = 1;
     }

   ISIS_FREE (undetected_model_bin);

   return 0;
}

/*}}}*/

static int get_wavelength_grid (Isis_Rmf_Grid_Type *g, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   unsigned int i;

   *lo = (double *) ISIS_MALLOC (g->nbins * sizeof(double));
   *hi = (double *) ISIS_MALLOC (g->nbins * sizeof(double));

   if ((*lo == NULL) || (*hi == NULL))
     {
        ISIS_FREE (*lo);
        ISIS_FREE (*hi);
        return -1;
     }

   *n = g->nbins;

   for (i = 0; i < g->nbins; i++)
     {
        double ghi = g->bin_hi[*n - i - 1];
        double glo = g->bin_lo[*n - i - 1];

        if (ghi <= 0.0)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Negative values in RMF energy grid");
             return -1;
          }
        if (glo <= 0.0) glo = 1.e-4 * ghi;

        /* FIXME:  assumed file values are in keV */
        (*lo)[i] = KEV_ANGSTROM / ghi;
        (*hi)[i] = KEV_ANGSTROM / glo;
     }

   return 0;
}
/*}}}*/

static int set_arf_grid (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int n) /*{{{*/
{
   (void) rmf; (void) lo; (void) hi; (void) n;
   return 0;
}
/*}}}*/

static int get_arf_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   Rmf_Client_Data_t *cd = get_client_data (rmf);

   return get_wavelength_grid (cd->arf, lo, hi, n);
}
/*}}}*/

static int set_data_grid (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int n) /*{{{*/
{
   (void) rmf; (void) lo; (void) hi; (void) n;
   return 0;
}
/*}}}*/

static int get_data_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n, int *energy_ordered_ebounds) /*{{{*/
{
   Rmf_Client_Data_t *cd = get_client_data (rmf);

   if (energy_ordered_ebounds)
     *energy_ordered_ebounds = cd->energy_ordered_ebounds;

   return get_wavelength_grid (cd->ebounds, lo, hi, n);
}
/*}}}*/

static int dummy_init (Isis_Rmf_t *rmf) /*{{{*/
{
   (void) rmf;
   return 0;
}

/*}}}*/

static int factor_rsp (Isis_Rmf_t *rmf, double *arf) /*{{{*/
{
   Rmf_Client_Data_t *cd = get_client_data (rmf);
   Rmf_Element_t *elem;
   float *response;
   int k, num_channels;
   unsigned int in_lam, negative_sum, num_ebins;

   if (rmf == NULL || arf == NULL || cd == NULL)
     return -1;

   /* factor ARF out of RSP matrix
    *   RSP => RMF * ARF
    * where RMF is normalized.
    */

   num_ebins = cd->num_ebins;
   negative_sum = 0;

   for (in_lam = 0; in_lam < num_ebins; in_lam++)
     {
        Rmf_Vector_t *v;
        unsigned int g, in_egy, num_grps;
        double sum;

        in_egy = cd->num_ebins - in_lam - 1;

        v = &cd->v[in_egy];
        num_grps = v->num_grps;

        sum = 0.0;
        arf[in_lam] = 0.0;

        /* sum over RMF groups */
        for (g = 0; g < num_grps; g++)
          {
             elem = &v->elem[g];
             response = elem->response;
             num_channels = elem->num_channels;

             for (k = 0; k < num_channels; k++)
               sum += response[k];
          }

        if (sum == 0.0)
          continue;
        else if (sum < 0.0)
          {
             negative_sum++;
             continue;
          }

        /* normalize */
        for (g = 0; g < num_grps; g++)
          {
             elem = &v->elem[g];
             response = elem->response;
             num_channels = elem->num_channels;

             for (k = 0; k < num_channels; k++)
               response[k] /= sum;
          }

        arf[in_lam] = sum;
     }

   rmf->includes_effective_area = 0;

   if (negative_sum)
     isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "corrupted RSP?  %d RSP groups had norm < 0", negative_sum);

   return 0;
}

/*}}}*/

/* Rebinning RMFs */

static int new_hist (unsigned int num, double **lop, double **hip, double **hp) /*{{{*/
{
   double *lo, *hi, *h;

   if (NULL == (lo = (double *) ISIS_MALLOC (num * sizeof(double))))
     return -1;

   if (NULL == (hi = (double *) ISIS_MALLOC (num * sizeof(double))))
     {
        ISIS_FREE (lo);
        return -1;
     }

   if (NULL == (h = (double *) ISIS_MALLOC (num * sizeof(double))))
     {
        ISIS_FREE (lo);
        ISIS_FREE (hi);
        return -1;
     }

   *lop = lo;
   *hip = hi;
   *hp = h;

   return 0;
}

/*}}}*/

static void free_hist (double *lo, double *hi, double *h) /*{{{*/
{
   ISIS_FREE (lo);
   ISIS_FREE (hi);
   ISIS_FREE (h);
}

/*}}}*/

static int break_into_groups (double *rmf, unsigned int num, /*{{{*/
                              unsigned int *n_grp, int *f_chan, int *n_chan, int num_chan,
                              double threshold)
{
   unsigned int i;
   int num_groups;

   num_groups = 0;
   i = 0;
   while (i < num)
     {
        unsigned int first_chan;

        if (rmf[i] <= threshold)
          {
             i++;
             continue;
          }

        first_chan = i;
        while ((i < num)
               && (rmf[i] > threshold))
          i++;

        if (num_groups == num_chan)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "break_into_groups:  too many groups");
             return -1;
          }

        n_chan[num_groups] = (int) (i - first_chan);
        f_chan[num_groups] = first_chan;

        num_groups++;
     }

   *n_grp = num_groups;
   return 0;
}

/*}}}*/

static int store_rmf_histogram (double *h, unsigned int num, unsigned int offset, /*{{{*/
                                double threshold, int *f_chan, int *n_chan, int num_chan,
                                Rmf_Vector_t *v)
{
   unsigned int i, n_grp;
   Rmf_Element_t *elem;

   if (-1 == break_into_groups (h + offset, num, &n_grp, f_chan, n_chan, num_chan, threshold))
     return -1;

   v->num_grps = n_grp;

   if (n_grp == 0)
     {
        v->elem = NULL;
        return 0;
     }

   if (NULL == (elem = (Rmf_Element_t *) ISIS_MALLOC (n_grp * sizeof(Rmf_Element_t))))
     return -1;

   memset ((char *) elem, 0, n_grp * sizeof (Rmf_Element_t));

   v->elem = elem;

   for (i = 0; i < n_grp; i++)
     {
        unsigned int j, num_channels;
        float *response;
        double *h_i;

        elem->first_channel = f_chan[i] + offset;
        elem->num_channels = num_channels = n_chan[i];

        if (NULL == (response = (float *) ISIS_MALLOC (num_channels * sizeof(float))))
          return -1;
        elem->response = response;

        h_i = h + f_chan[i] + offset;
        for (j = 0; j < num_channels; j++)
          response[j] = (float) h_i[j];

        elem++;
     }

   return 0;
}

/*}}}*/

static void free_rmf_vectors (Rmf_Vector_t *v, unsigned int num) /*{{{*/
{
   unsigned int i;

   if (v == NULL)
     return;

   for (i = 0; i < num; i++)
     free_rmf_vector (&v[i]);
   ISIS_FREE (v);
}

/*}}}*/

static int rebin_rmf (Isis_Rmf_t *rmf, double *wv_lo, double *wv_hi, unsigned int new_num) /*{{{*/
{
   double *new_lo, *new_hi, *new_h;
   double *old_lo, *old_hi, *old_h;
   unsigned int old_num;
   Isis_Rmf_Grid_Type *ebounds;
   Rmf_Client_Data_t *cd;
   Rmf_Vector_t *new_v;
   unsigned int i;
   int *f_chan, *n_chan;
   unsigned int num_rows;

   cd = get_client_data(rmf);
   if (cd == NULL)
     return -1;

   f_chan = n_chan = NULL;
   new_lo = new_hi = new_h = NULL;
   old_h = NULL;
   new_v = NULL;

   /* Note: the lo/hi grid values that are passed in here are in wavelength
    * units (ascending order).  This is checked by the calling
    * routine. So first convert them to an energy grid in ascending
    * order.
    */
   if (-1 == new_hist (new_num, &new_lo, &new_hi, &new_h))
     goto return_error;
   memset ((char *)new_h, 0, new_num * sizeof (double));

   for (i = 0; i < new_num; i++)
     {
        new_hi[new_num-i-1] = KEV_ANGSTROM/wv_lo[i];
        new_lo[new_num-i-1] = KEV_ANGSTROM/wv_hi[i];
     }

   ebounds = cd->ebounds;
   old_num = ebounds->nbins;
   old_lo = ebounds->bin_lo;
   old_hi = ebounds->bin_hi;

   if (NULL == (old_h = (double *) ISIS_MALLOC (old_num * sizeof(double))))
     goto return_error;
   memset ((char *)old_h, 0, old_num * sizeof (double));

   if ((NULL == (f_chan = (int *) ISIS_MALLOC (new_num * sizeof(int))))
       || (NULL == (n_chan = (int *) ISIS_MALLOC (new_num * sizeof(int)))))
     goto return_error;

   num_rows = cd->num_ebins;

   if (NULL == (new_v = (Rmf_Vector_t *) ISIS_MALLOC (num_rows * sizeof(Rmf_Vector_t))))
     goto return_error;
   memset ((char *) new_v, 0, num_rows * sizeof (Rmf_Vector_t));

   for (i = 0; i < num_rows; i++)
     {
        Rmf_Vector_t *v;
        unsigned int g;
        double *old_h_start, *old_h_end;
        unsigned int old_h_num, old_h_offset, new_h_num;
        int i_new_start, i_new_end;

        v = cd->v + i;

        if (v->num_grps == 0)
          continue;

        old_h_start = NULL;
        old_h_end = NULL;

        for (g = 0; g < v->num_grps; g++)
          {
             Rmf_Element_t *elem = &v->elem[g];
             unsigned int num_channels = elem->num_channels;
             float *response = elem->response;
             double *h = old_h + elem->first_channel;
             unsigned int k;

             if ((old_h_end == NULL) || ((h + num_channels - 1) > old_h_end))
               old_h_end = h + num_channels - 1;
             if ((old_h_start == NULL) || (h < old_h_start))
               old_h_start = h;

             for (k = 0; k < num_channels; k++)
               h[k] = response[k];
          }

        old_h_offset = old_h_start - old_h;
        old_h_num = old_h_end - old_h_start + 1;

        i_new_start = find_bin (old_lo[old_h_offset],
                                new_lo, new_hi, new_num);
        if (i_new_start < 0)
          i_new_start = 0;

        i_new_end = find_bin (old_hi[old_h_offset+old_h_num-1],
                              new_lo, new_hi, new_num);
        if (i_new_end < 0)
          i_new_end = new_num-1;

        new_h_num = i_new_end - i_new_start + 1;

        (void) rebin_histogram (old_h + old_h_offset,
                                old_lo + old_h_offset,
                                old_hi + old_h_offset,
                                old_h_num,
                                new_h + i_new_start,
                                new_lo + i_new_start,
                                new_hi + i_new_start,
                                new_h_num);

        /* Re-using the file threshold here causes problems.
         * The simplest solution is to use a zero threshold. */
        if (-1 == store_rmf_histogram (new_h, new_h_num, i_new_start,
                                       0.0, f_chan, n_chan, new_num, new_v + i))
          goto return_error;

        memset ((char *)new_h+i_new_start, 0, new_h_num * sizeof (double));
        memset ((char *)old_h+old_h_offset, 0, old_h_num * sizeof (double));
     }

   /* If we made it this far, then it has been a success. So make the
    * appropriate replacements
    */
   free_rmf_vectors (cd->v, num_rows);

   ISIS_FREE (ebounds->bin_lo);
   ISIS_FREE (ebounds->bin_hi);

   cd->v = new_v;
   ebounds->bin_lo = new_lo;
   ebounds->bin_hi = new_hi;
   ebounds->nbins = new_num;

   ISIS_FREE (f_chan);
   ISIS_FREE (n_chan);
   ISIS_FREE (new_h);
   ISIS_FREE (old_h);

   return 0;

   return_error:

   ISIS_FREE (f_chan);
   ISIS_FREE (n_chan);
   ISIS_FREE (new_h);

   free_hist (new_lo, new_hi, new_h);
   free_rmf_vectors (new_v, new_num);

   return -1;
}

/*}}}*/

/* load options */

static int handle_strict_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Rmf_Client_Data_t *cd = (Rmf_Client_Data_t *)clientdata;
   int strict;

   if (1 != sscanf (value, "%d", &strict))
     {
        fprintf (stderr, "Unknown '%s;%s' option value '%s'\n", subsystem, optname, value);
        return -1;
     }

   cd->strict = strict;

   return 0;
}

/*}}}*/

static int handle_matrix_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Rmf_Client_Data_t *cd = (Rmf_Client_Data_t *)clientdata;

   (void) subsystem; (void) optname;

   ISIS_FREE(cd->matrix_extname);
   if (NULL == (cd->matrix_extname = isis_make_string (value)))
     return -1;

   return 0;
}

/*}}}*/

static int handle_ebounds_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Rmf_Client_Data_t *cd = (Rmf_Client_Data_t *)clientdata;

   (void) subsystem; (void) optname;

   ISIS_FREE(cd->ebounds_extname);
   if (NULL == (cd->ebounds_extname = isis_make_string (value)))
     return -1;

   return 0;
}

/*}}}*/

static Isis_Option_Table_Type Option_Table [] =
{
     {"strict", handle_strict_option, ISIS_OPT_REQUIRES_VALUE, "2", "OGIP strictness"},
     {"ebounds", handle_ebounds_option, ISIS_OPT_REQUIRES_VALUE, "EBOUNDS", "EXTNAME of FITS extension containing EBOUNDS grid"},
     {"matrix", handle_matrix_option, ISIS_OPT_REQUIRES_VALUE, "SPECRESP MATRIX", "EXTNAME of FITS extension containing RMF matrix"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_options (Rmf_Client_Data_t *cd, Isis_Option_Type *opts) /*{{{*/
{
   cd->strict = Isis_Rmf_OGIP_Compliance;
   return isis_process_options (opts, Option_Table, (void *)cd, 1);
}

/*}}}*/

static int parse_options (Rmf_Client_Data_t *cd, char *options) /*{{{*/
{
   Isis_Option_Type *opts;

   if (NULL == (opts = isis_parse_option_string (options)))
     return -1;

   cd->type = RMF_TYPE_FILE;
   if (NULL == (cd->f.file = isis_make_string (opts->subsystem)))
     return -1;

   if (-1 == set_options (cd, opts))
     {
        isis_free_options (opts);
        return -1;
     }

   isis_free_options (opts);

   return 0;
}

/*}}}*/

static int print_options_help (Rmf_Client_Data_t *cd) /*{{{*/
{
   Isis_Option_Type *opts;
   int status;
   char *s;

   if (cd->type == RMF_TYPE_SLANG)
     {
        return 0;
     }

   if (NULL == (s = isis_mkstrcat (cd->f.file, ";help", NULL)))
     return -1;
   if (NULL == (opts = isis_parse_option_string (s)))
     {
        ISIS_FREE(s);
        return -1;
     }
   status = isis_process_options (opts, Option_Table, (void *)cd, 1);
   ISIS_FREE(s);
   isis_free_options (opts);
   return status;
}

/*}}}*/

static int init_rmf_t (Isis_Rmf_t *rmf) /*{{{*/
{
   rmf->set_arf_grid = set_arf_grid;
   rmf->set_data_grid = set_data_grid;
   rmf->get_arf_grid = get_arf_grid;
   rmf->get_data_grid = get_data_grid;
   rmf->init = dummy_init;
   rmf->redistribute = redistribute;
   rmf->set_noticed_model_bins = set_noticed_model_bins;
   rmf->rebin_rmf = rebin_rmf;
   rmf->factor_rsp = factor_rsp;
   rmf->delete_client_data = delete_client_data;

   rmf->client_data = (Rmf_Client_Data_t *) ISIS_MALLOC (sizeof(Rmf_Client_Data_t));
   if (rmf->client_data == NULL)
     return -1;
   memset ((char *)rmf->client_data, 0, sizeof (Rmf_Client_Data_t));

   return 0;
}

/*}}}*/

int Rmf_load_file (Isis_Rmf_t *rmf, void *options)  /*{{{*/
{
   Rmf_Client_Data_t *cd;

   if (-1 == init_rmf_t (rmf))
     return -1;

   if (NULL == (rmf->arg_string = isis_make_string ((const char *)options)))
     return -1;

   cd = (Rmf_Client_Data_t *) rmf->client_data;

   if (-1 == parse_options (cd, (char *)options))
     return -1;

   return read_rmf (rmf, cd->f.file);
}

/*}}}*/

int Rmf_load_slang (Isis_Rmf_t *rmf, void *options) /*{{{*/
{
   Rmf_SLang_Info_Type *info;
   Rmf_Client_Data_t *cd;
   unsigned int i, num_data_bins, num_arf_bins;
   double *en_lo, *en_hi;
   int *f_chan = NULL, *n_chan = NULL;

   if (-1 == init_rmf_t (rmf))
     return -1;

   info = (Rmf_SLang_Info_Type *)options;

   if (NULL == (rmf->arg_string = isis_make_string (info->func->name)))
     return -1;

   cd = (Rmf_Client_Data_t *)rmf->client_data;

   rmf->includes_effective_area = 0;

   num_arf_bins = cd->num_ebins = info->arf_bin_lo->num_elements;
   num_data_bins = info->data_bin_lo->num_elements;
   cd->threshold = info->threshold;

   if (NULL == (cd->arf = Isis_new_rmf_grid (cd->num_ebins,
                                             (double *)info->arf_bin_lo->data,
                                             (double *)info->arf_bin_hi->data)))
     {
        goto return_error;
     }
   cd->arf->units = U_KEV;

   if (NULL == (cd->ebounds = Isis_new_rmf_grid (num_data_bins,
                                                 (double *)info->data_bin_lo->data,
                                                 (double *)info->data_bin_hi->data)))
     {
        goto return_error;
     }
   cd->ebounds->units = U_KEV;

   if (NULL == (cd->v = (Rmf_Vector_t *) ISIS_MALLOC (cd->num_ebins * sizeof(Rmf_Vector_t))))
     {
        goto return_error;
     }
   memset ((char *)cd->v, 0, cd->num_ebins * sizeof (Rmf_Vector_t));

   if ((NULL == (n_chan = (int *) ISIS_MALLOC (num_data_bins*sizeof(int))))
       || (NULL == (f_chan = (int *) ISIS_MALLOC (num_data_bins*sizeof(int)))))
     goto return_error;

   en_lo = (double *) info->arf_bin_lo->data;
   en_hi = (double *) info->arf_bin_hi->data;
   for (i = 0; i < num_arf_bins; i++)
     {
        double en = 0.5*(en_lo[i] + en_hi[i]);
        SLang_Array_Type *at_rmf;

        /* Evaluate \int dE' R(E',E). */
        if ((-1 == SLang_start_arg_list ())
            || (-1 == SLang_push_array (info->data_bin_lo, 0))
            || (-1 == SLang_push_array (info->data_bin_hi, 0))
            || (-1 == SLang_push_double (en))
            || ((info->client_data != NULL)
                && (-1 == SLang_push_anytype (info->client_data)))
            || (-1 == SLang_end_arg_list ())
            || (-1 == SLexecute_function (info->func))
            || (-1 == SLang_pop_array_of_type (&at_rmf, SLANG_DOUBLE_TYPE)))
          goto return_error;

        if (at_rmf->num_elements != num_data_bins)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "The computed RMF profile does not contain the expected number of bins");
             SLang_free_array (at_rmf);
             goto return_error;
          }
        if (-1 == store_rmf_histogram ((double *)at_rmf->data, num_data_bins, 0,
                                       info->threshold, f_chan, n_chan, num_data_bins,
                                       cd->v + i))
          {
             SLang_free_array (at_rmf);
             goto return_error;
          }
        SLang_free_array (at_rmf);
     }

   ISIS_FREE (n_chan);
   ISIS_FREE (f_chan);
   return 0;

return_error:

   ISIS_FREE (n_chan);
   ISIS_FREE (f_chan);
   return -1;
}

/*}}}*/
