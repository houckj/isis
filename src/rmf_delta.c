/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2019  Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Authors:  John C. Houck  <houck@space.mit.edu>

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

/* $Id: rmf_delta.c,v 1.5 2004/05/01 23:05:02 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "util.h"
#include "rmf.h"

/*}}}*/

typedef struct
{
   union
     {
        double s;
        double *v;
     }
   resp;
   union
     {
        unsigned int s;
        unsigned int *v;
     }
   chan;

   unsigned int num_chan;
   unsigned int size;
}
Elem_Type;

typedef struct
{
   Isis_Rmf_Grid_Type *arf;
   Isis_Rmf_Grid_Type *ebounds;
   Elem_Type *elem;
   int num_arf;
}
Rmf_Client_Data_t;

static void free_client_data (Rmf_Client_Data_t *cl) /*{{{*/
{
   Elem_Type *end, *e;
   if (cl == NULL) return;

   e = cl->elem;
   end = cl->elem + cl->num_arf;
   while (e < end)
     {
        if (e->num_chan > 1)
          {
             ISIS_FREE (e->resp.v);
             ISIS_FREE (e->chan.v);
          }
        e++;
     }
   ISIS_FREE (cl->elem);

   Isis_free_rmf_grid (cl->arf);
   Isis_free_rmf_grid (cl->ebounds);
}

/*}}}*/

static void delete_client_data (Isis_Rmf_t *rmf) /*{{{*/
{
   free_client_data ((Rmf_Client_Data_t *)rmf->client_data);
   ISIS_FREE (rmf->client_data);
}

/*}}}*/

static Elem_Type *get_response (Isis_Rmf_t *rmf) /*{{{*/
{
   Rmf_Client_Data_t *cl = (Rmf_Client_Data_t *)rmf->client_data;

   return cl->elem;
}

/*}}}*/

static int redistribute (Isis_Rmf_t *rmf, unsigned int in_lam, double flux, /*{{{*/
                         double *det_chan, unsigned int num_ebounds)
{
   Elem_Type *elem = get_response (rmf);
   Elem_Type *e;
   unsigned int n;

   (void) num_ebounds;

   if (elem == NULL)
     return -1;

   e = &elem[in_lam];
   n = e->num_chan;

   if (n == 1)
     {
        det_chan[ e->chan.s ] += flux * e->resp.s;
     }
   else if (n > 1)
     {
        double *resp = e->resp.v;
        unsigned int i, *chv = e->chan.v;

        for (i = 0; i < n; i++)
          det_chan[chv[i]] += flux * resp[i];
     }

   return 0;
}

/*}}}*/

static int set_noticed_model_bins (Isis_Rmf_t *rmf, int num_chan, int *chan_notice, /*{{{*/
                                   int num_model, int *model_notice)
{
   Elem_Type *elem = get_response (rmf);
   int k;

   (void) num_chan;

   memset ((char *)model_notice, 0, num_model * sizeof(int));

   for (k = 0; k < num_model; k++)
     {
        Elem_Type *e = &elem[k];
        unsigned int n = e->num_chan;

        if (n <= 1)
          {
             if (n == 0 || chan_notice[ e->chan.s ])
               model_notice[k] = 1;
          }
        else if (n > 1)
          {
             unsigned int i, *chv = e->chan.v;
             for (i = 0; i < n; i++)
               {
                  if (chan_notice[chv[i]])
                    model_notice[k] = 1;
               }
          }
     }

   return 0;
}

/*}}}*/

static int record_response (Elem_Type *e, double f, unsigned int ch) /*{{{*/
{
   unsigned int i, start_size = 16;

   e->num_chan++;

   if (e->num_chan == 1)
     {
        e->chan.s = ch;
        e->resp.s = f;
        return 0;
     }
   else if (e->num_chan == 2)
     {
        unsigned int ch0 = e->chan.s;
        double f0 = e->resp.s;

        e->resp.v = (double *) ISIS_MALLOC (start_size * sizeof(double));
        e->chan.v = (unsigned int *) ISIS_MALLOC (start_size * sizeof(unsigned int));
        if (e->resp.v == NULL || e->chan.v == NULL)
          return -1;

        e->size = start_size;
        e->resp.v[0] = f0;
        e->chan.v[0] = ch0;
     }

   if (e->num_chan > e->size)
     {
        unsigned int *i_tmp;
        double *d_tmp;
        e->size *= 2;

        if (NULL == (i_tmp = (unsigned int *) ISIS_REALLOC (e->chan.v, e->size * sizeof(unsigned int))))
          return -1;
        e->chan.v = i_tmp;
        
        if (NULL == (d_tmp = (double *) ISIS_REALLOC (e->resp.v, e->size * sizeof(double))))
          return -1;
        e->resp.v = d_tmp;
     }

   i = e->num_chan - 1;
   e->chan.v[i] = ch;
   e->resp.v[i] = f;

   return 0;
}

/*}}}*/

static double overlap (double alo, double ahi, double dlo, double dhi) /*{{{*/
{
   double max_min, min_max, o;

   if (alo > dlo) max_min = alo;
   else max_min = dlo;

   if (ahi < dhi) min_max = ahi;
   else min_max = dhi;

   o = min_max - max_min;

   return (o > 0.0) ? o : 0.0;
}

/*}}}*/

static int init_client_data (Isis_Rmf_t *rmf) /*{{{*/
{
   Rmf_Client_Data_t *cl = (Rmf_Client_Data_t *) rmf->client_data;
   Isis_Rmf_Grid_Type *arf = cl->arf;
   Isis_Rmf_Grid_Type *det = cl->ebounds;
   Elem_Type *elem;
   unsigned int i, ch, na, nd;
   double *alo, *ahi, *dlo, *dhi;

   elem = (Elem_Type *) ISIS_MALLOC (arf->nbins * sizeof(Elem_Type));
   if (elem == NULL)
     return -1;
   memset ((char *)elem, 0, arf->nbins * sizeof (*elem));

   cl->elem = elem;
   cl->num_arf = arf->nbins;

   na = arf->nbins;
   alo = arf->bin_lo;
   ahi = arf->bin_hi;

   nd = det->nbins;
   dlo = det->bin_lo;
   dhi = det->bin_hi;

   for (ch = 0; (ch < nd) && (dhi[ch] < alo[0]); ch++)
     ;

   for ( ; (ch < nd) && (dlo[ch] < ahi[na-1]); ch++)
     {
        int i1 = find_bin (dlo[ch], alo, ahi, na);
        if (i1 < 0) i1 = 0;

        for (i = i1; (i < na) && (alo[i] < dhi[ch]); i++)
          {
             double o = overlap (alo[i], ahi[i], dlo[ch], dhi[ch]);
             if (o > 0.0)
               {
                  double r = o / (ahi[i] - alo[i]);
                  if (-1 == record_response (&elem[i], r, ch))
                    return -1;
               }
          }
     }

   return 0;
}

/*}}}*/

static int set_grid (Isis_Rmf_Grid_Type **g, double *lo, double *hi, unsigned int n) /*{{{*/
{
   Isis_free_rmf_grid (*g);
   *g = Isis_new_rmf_grid (n, lo, hi);
   if ((*g) == NULL)
     return -1;

   return 0;
}

/*}}}*/

static int get_grid (Isis_Rmf_Grid_Type *g, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   *lo = (double *) ISIS_MALLOC (g->nbins * sizeof(double));
   *hi = (double *) ISIS_MALLOC (g->nbins * sizeof(double));

   if ((*lo == NULL) || (*hi == NULL))
     {
        ISIS_FREE (*lo);
        ISIS_FREE (*hi);
        return -1;
     }

   *n = g->nbins;

   memcpy ((char *)*lo, (char *)g->bin_lo, g->nbins * sizeof(double));
   memcpy ((char *)*hi, (char *)g->bin_hi, g->nbins * sizeof(double));

   return 0;
}

/*}}}*/

static int set_arf_grid (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int n) /*{{{*/
{
   Rmf_Client_Data_t *cl = (Rmf_Client_Data_t *)rmf->client_data;

   return set_grid (&cl->arf, lo, hi, n);
}

/*}}}*/

static int get_arf_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   Rmf_Client_Data_t *cl = (Rmf_Client_Data_t *)rmf->client_data;

   return get_grid (cl->arf, lo, hi, n);
}

/*}}}*/

static int set_data_grid (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int n) /*{{{*/
{
   Rmf_Client_Data_t *cl = (Rmf_Client_Data_t *)rmf->client_data;

   return set_grid (&cl->ebounds, lo, hi, n);
}

/*}}}*/

static int get_data_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n, int *eoeb) /*{{{*/
{
   Rmf_Client_Data_t *cl = (Rmf_Client_Data_t *)rmf->client_data;

   if (eoeb) *eoeb = 0;
   return get_grid (cl->ebounds, lo, hi, n);
}

/*}}}*/

int Rmf_load_delta (Isis_Rmf_t *rmf, void *options) /*{{{*/
{
   Rmf_Client_Data_t *cl;

   if (rmf == NULL)
     return -1;

   (void) options;

   rmf->set_arf_grid = set_arf_grid;
   rmf->get_arf_grid = get_arf_grid;
   rmf->set_data_grid = set_data_grid;
   rmf->get_data_grid = get_data_grid;

   rmf->init = init_client_data;
   rmf->redistribute = redistribute;
   rmf->set_noticed_model_bins = set_noticed_model_bins;

   rmf->delete_client_data = delete_client_data;

   if (NULL == (cl = (Rmf_Client_Data_t *) ISIS_MALLOC (sizeof(Rmf_Client_Data_t))))
     return -1;
   memset ((char *)cl, 0, sizeof(*cl));

   rmf->client_data = (void *) cl;

   return 0;
}
/*}}}*/
