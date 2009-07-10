/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008  Massachusetts Institute of Technology

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

/* $Id: rmf_user.c,v 1.14 2004/02/09 11:14:25 houck Exp $ */

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

#define MALLOC(x)  malloc(x)

/*}}}*/

/*{{{ typedefs */

typedef struct Client_Data_t Client_Data_t;
typedef int Shape_Fun_Type (Client_Data_t *, double, double, double, float *);

typedef struct
{
   union
     {
        float s;
        float *v;
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

struct Client_Data_t
{
   Isis_Rmf_Grid_Type *arf;
   Isis_Rmf_Grid_Type *ebounds;
   Elem_Type *elem;
   Shape_Fun_Type *shape_fun;
   double sigma;
   float threshold;
   int num_arf;
};

/*}}}*/

static int find_bin (double x, double *lo, double *hi, int n) /*{{{*/
{
   int n0, n1, n2;

   if (NULL == lo || NULL == hi || n <= 0)
     return -1;

   if (x < lo[0] || hi[n-1] < x)
     return -1;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;

        if (lo[n2] <= x && x < hi[n2])
          return n2;

        if (x < hi[n2])
          n1 = n2;
        else
          n0 = n2;
     }

   return n0;
}

/*}}}*/

static void free_client_data (Client_Data_t *cl) /*{{{*/
{
   Elem_Type *end, *e;
   if (cl == NULL) return;

   e = cl->elem;
   end = cl->elem + cl->num_arf;
   while (e < end)
     {
        if (e->num_chan > 1)
          {
             ISIS_FREE(e->resp.v);
             ISIS_FREE(e->chan.v);
          }
        e++;
     }
   ISIS_FREE(cl->elem);

   Isis_free_rmf_grid (cl->arf);
   Isis_free_rmf_grid (cl->ebounds);
}

/*}}}*/

static void delete_client_data (Isis_Rmf_t *rmf) /*{{{*/
{
   free_client_data ((Client_Data_t *)rmf->client_data);
   ISIS_FREE(rmf->client_data);
}

/*}}}*/

static Elem_Type *get_response (Isis_Rmf_t *rmf) /*{{{*/
{
   Client_Data_t *cl = (Client_Data_t *)rmf->client_data;

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
        float *resp = e->resp.v;
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

static int record_response (Elem_Type *e, float f, unsigned int ch) /*{{{*/
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
        float f0 = e->resp.s;

        e->resp.v = (float *) MALLOC (start_size * sizeof (float));
        e->chan.v = (unsigned int *) MALLOC (start_size * sizeof (unsigned int));
        if (e->resp.v == NULL || e->chan.v == NULL)
          return -1;

        e->size = start_size;
        e->resp.v[0] = f0;
        e->chan.v[0] = ch0;
     }

   if (e->num_chan > e->size)
     {
        unsigned int *i_tmp;
        float *f_tmp;
        e->size *= 2;

        if (NULL == (i_tmp = (unsigned int *) ISIS_REALLOC (e->chan.v, e->size * sizeof(unsigned int))))
          return -1;
        e->chan.v = i_tmp;

        if (NULL == (f_tmp = (float *) ISIS_REALLOC (e->resp.v, e->size * sizeof(float))))
          return -1;
        e->resp.v = f_tmp;
     }

   i = e->num_chan - 1;
   e->chan.v[i] = ch;
   e->resp.v[i] = f;

   return 0;
}

/*}}}*/

static int init_client_data (Isis_Rmf_t *rmf) /*{{{*/
{
   Client_Data_t *cl = (Client_Data_t *) rmf->client_data;
   Isis_Rmf_Grid_Type *arf = cl->arf;
   Isis_Rmf_Grid_Type *det = cl->ebounds;
   Shape_Fun_Type *shape_fun = cl->shape_fun;
   float threshold = cl->threshold;
   Elem_Type *elem;
   unsigned int i, na, nd;
   double *alo, *ahi, *dlo, *dhi;

   if (NULL == (elem = (Elem_Type *) MALLOC (arf->nbins * sizeof (*elem))))
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

   for (i = 0; i < na; i++)
     {
        double lam = 0.5 * (alo[i] + ahi[i]);
        unsigned int ch;
        float f;
        int ch0;

        if ((ch0 = find_bin (lam, dlo, dhi, nd)) < 0)
          continue;

        ch = ch0;
        while (ch-- > 0)
          {
             if (-1 == (*shape_fun) (cl, lam, dlo[ch], dhi[ch], &f))
               return -1;
             if (f < threshold)
               break;
             if (-1 == record_response (&elem[i], f, ch))
               return -1;
          }

        for (ch = ch0; ch < nd; ch++)
          {
             if (-1 == (*shape_fun) (cl, lam, dlo[ch], dhi[ch], &f))
               return -1;
             if (f < threshold)
               break;
             if (-1 == record_response (&elem[i], f, ch))
               return -1;
          }
     }

   return 0;
}

/*}}}*/

static int set_grid (Isis_Rmf_Grid_Type **g, double *lo, double *hi, unsigned int n) /*{{{*/
{
   Isis_free_rmf_grid (*g);
   *g = Isis_new_rmf_grid (n);
   if ((*g) == NULL)
     return -1;

   memcpy ((char *)(*g)->bin_lo, (char *)lo, n * sizeof(double));
   memcpy ((char *)(*g)->bin_hi, (char *)hi, n * sizeof(double));

   return 0;
}

/*}}}*/

static int get_grid (Isis_Rmf_Grid_Type *g, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   *lo = (double *) MALLOC (g->nbins * sizeof(double));
   *hi = (double *) MALLOC (g->nbins * sizeof(double));

   if ((*lo == NULL) || (*hi == NULL))
     {
        ISIS_FREE(*lo);
        ISIS_FREE(*hi);
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
   Client_Data_t *cl = (Client_Data_t *)rmf->client_data;

   return set_grid (&cl->arf, lo, hi, n);
}

/*}}}*/

static int get_arf_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n) /*{{{*/
{
   Client_Data_t *cl = (Client_Data_t *)rmf->client_data;

   return get_grid (cl->arf, lo, hi, n);
}

/*}}}*/

static int set_data_grid (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int n) /*{{{*/
{
   Client_Data_t *cl = (Client_Data_t *)rmf->client_data;
   Isis_Rmf_Grid_Type *g;
   unsigned int i;
   int order;

   Isis_free_rmf_grid (cl->ebounds);
   cl->ebounds = Isis_new_rmf_grid (n);
   if (cl->ebounds == NULL)
     return -1;

   g = cl->ebounds;

   if (abs(rmf->order) > 1)
     order = rmf->order;
   else
     order = 1;

   for (i = 0; i < n; i++)
     {
        g->bin_lo[i] = lo[i] / order;
        g->bin_hi[i] = hi[i] / order;
     }

   return 0;
}

/*}}}*/

static int get_data_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *n, int *eoeb) /*{{{*/
{
   Client_Data_t *cl = (Client_Data_t *)rmf->client_data;

   *eoeb = 0;
   return get_grid (cl->ebounds, lo, hi, n);
}

/*}}}*/

static int shape_fun (Client_Data_t *cl, double lam, double dlo, double dhi, float *f) /*{{{*/
{
   double sigma = cl->sigma;

   if (sigma > 0.0)
     {
        double dxl, dxh;
        dxh = (dhi - lam) / sigma;
        dxl = (dlo - lam) / sigma;
        *f = (float) (isis_gpf (dxh) - isis_gpf (dxl));
     }
   else *f = 0.0;

   return 0;
}

/*}}}*/

static int handle_order_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Rmf_t *rmf = (Isis_Rmf_t *)clientdata;
   int order;

   if (1 != sscanf (value, "%d", &order))
     {
        fprintf (stderr, "Unknown '%s;%s' option value '%s'\n", subsystem, optname, value);
        return -1;
     }

   rmf->order = order;

   return 0;
}

/*}}}*/

static int handle_sigma_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Rmf_t *rmf = (Isis_Rmf_t *)clientdata;
   Client_Data_t *cl = (Client_Data_t *)rmf->client_data;
   double sigma;

   if (1 != sscanf (value, "%le", &sigma))
     {
        fprintf (stderr, "Unknown '%s;%s' option value '%s'\n", subsystem, optname, value);
        return -1;
     }

   cl->sigma = sigma;

   return 0;
}

/*}}}*/

static Isis_Option_Table_Type Option_Table [] =
{
     {"order", handle_order_option, ISIS_OPT_REQUIRES_VALUE, "1", "dispersion order"},
     {"sigma", handle_sigma_option, ISIS_OPT_REQUIRES_VALUE, "0.025", "Gaussian sigma"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_options (Isis_Rmf_t *rmf, Isis_Option_Type *opts) /*{{{*/
{
   return isis_process_options (opts, Option_Table, (void *)rmf, 1);
}

/*}}}*/

static int parse_options (Isis_Rmf_t *rmf, char *options) /*{{{*/
{
   Isis_Option_Type *opts;

   if (NULL == (opts = isis_parse_option_string (options)))
     return -1;

   if (-1 == set_options (rmf, opts))
     {
        isis_free_options (opts);
        return -1;
     }

   isis_free_options (opts);

   return 0;
}

/*}}}*/

ISIS_RMF_METHOD(gauss_rmf, rmf, options)/*{{{*/
{
   Client_Data_t *cl;

   if (rmf == NULL)
     return -1;

   rmf->set_arf_grid = set_arf_grid;
   rmf->get_arf_grid = get_arf_grid;
   rmf->set_data_grid = set_data_grid;
   rmf->get_data_grid = get_data_grid;

   rmf->init = init_client_data;
   rmf->redistribute = redistribute;
   rmf->set_noticed_model_bins = set_noticed_model_bins;

   rmf->delete_client_data = delete_client_data;

   if (NULL == (cl = (Client_Data_t *) MALLOC (sizeof *cl)))
     return -1;
   memset ((char *)cl, 0, sizeof(*cl));
   rmf->client_data = (void *) cl;

   rmf->order = 1;
   cl->sigma = 0.025;
   cl->threshold = 1.e-4;
   cl->shape_fun = &shape_fun;

   return parse_options (rmf, options);
}
/*}}}*/
