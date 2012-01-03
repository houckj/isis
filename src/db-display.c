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

#include <slang.h>

#include "isis.h"
#include "db-atomic.h"
#include "errors.h"
#include "plot.h"
#include "util.h"

#include "db-display.h"

/*}}}*/

/* plot line groups, energy levels with transitions overlaid  */

/* line groups */

void DB_list_line_group_stats (FILE * fp, DB_t *db) /*{{{*/
{
   DB_line_group_t *head = DB_get_group_table_head (db);
   DB_line_group_t *t;

   if (NULL == DB_next_group (head))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "No line groups are currently defined.");
        return;
     }

   for (t = DB_next_group (head); t != NULL; t = DB_next_group (t))
     {
        DB_line_t **line;
        int group;
        int nlines;
        char *name = NULL;

        if (-1 == DB_get_group_index (&group, t)
            || -1 == DB_get_group_line_list (&nlines, &line, t))
          continue;

        fprintf(fp, "Group %3d [%6d lines]  ", group, nlines);

        if (-1 == DB_get_group_name (&name, t)
            || name == NULL)
          fputc ('\n', fp);
        else
          {
             fprintf (fp, "%s\n", name);
             ISIS_FREE (name);
          }
     }
}

/*}}}*/

static int compress_label (char *label, unsigned int size) /*{{{*/
{
   char *p, *t;
   int one_space = 0;

   t = label;

   for (p = label; (p != NULL) && (*p != 0); p++)
     {
        if (isspace((unsigned char)*p))
          {
             if (one_space)
               continue;
             one_space = 1;
          }
        else one_space = 0;

        *t++ = (*p == '\n') ? ' ' : *p;
     }

   if (t < label + size)
     *t = 0;

   return 0;
}

/*}}}*/

int DB_print_line_group (FILE * fp, DB_line_group_t *t, DB_t *db) /*{{{*/
{
   DB_line_t *line;
   DB_line_t **list;
   int i, nlines;
   char *name = NULL;
   char header[] =
"#  index        ion   lambda    F (ph/cm^2/s)  A(s^-1)  upper lower   label";

   if (NULL == t || NULL == fp || NULL == db)
     return -1;

   (void) DB_get_group_name (&name, t);

   if (name != NULL)
     {
        int ret = fprintf(fp,"%s\n", name);
        ISIS_FREE (name);
        if (ret < 0)
          return -1;
     }

   if (-1 == DB_get_group_line_list (&nlines, &list, t))
     return -1;

   if (fprintf(fp, "%s\n", header) < 0)
     return -1;

   for (i=0; i < nlines; i++)
     {
        unsigned char flg;
        int up, lo, Z, q;
        char s_up[LEVEL_NAME_SIZE], s_lo[LEVEL_NAME_SIZE];
        char label[2*LEVEL_NAME_SIZE+5];
        char ion_name[10];

        line = list[i];
        if (line == NULL)
          continue;

        Z = line->proton_number;
        q = line->ion_charge;
        up = line->upper_level;
        lo = line->lower_level;

        label[0] = 0;

        flg = line->have_emissivity_data ? '*' : ' ';

        (void) DB_get_level_label (s_up, Z, q, up, db);
        (void) DB_get_level_label (s_lo, Z, q, lo, db);

        if (-1 == DB_get_ion_name (ion_name, sizeof(ion_name), Z, q,
                                   DB_Ion_Format))
          (void) sprintf (ion_name, "%d  %d", Z, q);

        if (s_up[0] != '\0' && s_lo[0] != '\0')
          {
             (void) isis_strcat (label, sizeof(label), s_up, " - ", s_lo, NULL);
             compress_label (label, sizeof(label));
          }
        else
          label[0] = '\0';

        if (fprintf(fp, "%7d %c %9s %#9.5g    %9.3e   %9.3e  %4d  %4d %-s\n",
                    line->indx, flg, ion_name,
                    line->wavelen, line->flux, line->A, up, lo, label)
            < 0)
          return -1;          /* check fprintf return code in case fp is a pipe */
     }

   return 0;
}

/*}}}*/

int DB_plot_line_list (Plot_t *fmt, Line_Label_Style_Type *s, float redshift, /*{{{*/
                       float *lambda, char **labels, unsigned int nlines)
{
   double one_plus_z = 1.0 + redshift;
   float xlim[2] = {0.0, 0.0};
   float ylim[2] = {0.0, 0.0};
   float wl[2], y[2];
   int logx, x_unit;
   unsigned int i;
   int ret=-1;

   if (fmt == NULL || s == NULL)
     return -1;

   _Plot_set_charsize (s->char_height);
   Plot_set_style (fmt, s->use_color);
   _Plot_set_linestyle (PLOT_SOLID_LINE);
   Plot_set_line_width (fmt);

   logx = Plot_get_logx (fmt);
   x_unit = Plot_x_unit (fmt);

   Plot_query_plot_limits (&xlim[0], &xlim[1], &ylim[0], &ylim[1]);

   y[0] = s->bottom_frac * (ylim[1] - ylim[0]) + ylim[0];
   y[1] = s->top_frac * (ylim[1] - ylim[0]) + ylim[0];

   for (i=0; i < nlines; i++)
     {
        double z_wavelen = one_plus_z * lambda[i];
        float zw;

        if (x_unit != U_ANGSTROM)
          {
             if (-1 == unit_convert_x (&z_wavelen, 1, x_unit, U_ANGSTROM))
               goto return_error;
          }

        if (logx && z_wavelen > 0.0)
          z_wavelen = log10 (z_wavelen);

        zw = (float) z_wavelen;

        if (zw > xlim[1] || zw < xlim[0])
          continue;

        wl[0] = wl[1] = zw;

        Plot_set_lineclip_state (0);
        Plot_line (2, wl, y);
        Plot_set_lineclip_state (1);

        if (s->label_type >= 0)
          {
             float yl = y[1] + s->offset * (y[1] - y[0]);
             Plot_put_text (wl[1], yl, s->angle, s->justify, labels[i]);
          }
     }

   ret = 0;
   return_error:

   Plot_update_display ();
   Plot_set_linestyle (fmt);
   Plot_set_charsize (fmt);

   return ret;
}

/*}}}*/

int DB_plot_line_group (DB_line_group_t *t,  Line_Label_Style_Type *s, float redshift, /*{{{*/
                        char *(*label_massage_hook)(char *),
                        Plot_t *fmt, DB_t *db)
{
   DB_line_t **list;
   float *lambda = NULL;
   char **labels = NULL;
   char *label_text = NULL;
   unsigned int k, label_size;
   int i, nlines, ret = -1;

   if (-1 == DB_get_group_line_list (&nlines, &list, t))
     goto error_return;

   label_size = s->label_length * sizeof(char);

   if ((NULL == (lambda = (float *) ISIS_MALLOC (nlines * sizeof(float))))
       || (NULL == (labels = (char **) ISIS_MALLOC (nlines * sizeof (char *))))
       || (NULL == (label_text = (char *) ISIS_MALLOC (nlines * label_size))))
     goto error_return;

   memset ((char *)label_text, 0, nlines * label_size);
   for (i = 0; i < nlines; i++)
     {
        labels[i] = label_text + i * s->label_length;
     }

   k = 0;
   for (i=0; i < nlines; i++)
     {
        DB_line_t *p = list[i];

        if ((p == NULL)
            || (-1 == DB_get_line_namestring (labels[k], label_size, s->label_type, p, db)))
          continue;

        /* support user-defined text transformation */
        if (label_massage_hook)
          {
             char *x = (*label_massage_hook)(labels[k]);
             if (x)
               {
                  isis_strcpy (labels[k], x, label_size);
                  ISIS_FREE(x);
               }
          }

        lambda[k] = p->wavelen;
        k++;
     }

   ret = DB_plot_line_list (fmt, s, redshift, lambda, labels, k);

   error_return:

   ISIS_FREE (lambda);
   ISIS_FREE (labels);
   ISIS_FREE (label_text);

   return ret;
}

/*}}}*/

/* energy levels */

int DB_list_ion_levels (FILE * fp, int Z, int q, DB_t *db) /*{{{*/
{
   int lev, nlevels;
   DB_ion_t *ion = NULL;
   char ion_name[10];

   if (NULL == fp || NULL == db)
     return -1;

   if (NULL == (ion = DB_get_ion (db, Z, q))
       || -1 == DB_get_ion_nlevels (&nlevels, ion))
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "for Z=%d, q=%d", Z, q);
        return -1;
     }

   if (-1 == DB_get_ion_name (ion_name, sizeof(ion_name), Z, q,
                              DB_Ion_Format))
     (void) sprintf (ion_name, "%d  %d", Z, q);
   fprintf(fp, "  %s has %3d levels:\n", ion_name, nlevels);

   fprintf (fp, " index    energy   weight        A-sum   ndown    n  L   S      label\n");

   for(lev = 1; lev <= nlevels; ++lev)
     {
        DB_level_t *e;
        double sum;
        int j;

        if (NULL == (e = DB_get_ion_level (lev, ion)))
          continue;

        sum = 0.0;
        for (j = 0; j < e->ndown; j++)
          {
             DB_line_t *line = e->down[j];
             sum += line->A;
          }

        if (fprintf(fp, "%6d %9.3f   %5.2f   %11.4g  %4d     %2d %2d  %3.1f %s\n",
                    lev, e->energy, e->stat_weight, sum, e->ndown,
                    e->n, e->L, e->S, e->label)
            < 0)
          return -1;
     }

   return 0;
}

/*}}}*/

static int plot_level (float xcenter, float energy, float width) /*{{{*/
{
   float x[2], y[2];

   /* left and right endpoints */
   x[0] = xcenter - 0.5 * width;
   x[1] = xcenter + 0.5 * width;

   y[0] = y[1] = energy;

   Plot_line (2, x, y);

   return 0;
}

/*}}}*/

static int plot_transition (float xcu, float eu, float xcl, float el) /*{{{*/
{
   float x[2], y[2];

   x[0] = xcu;    y[0] = eu;        /* upper level */
   x[1] = xcl;    y[1] = el;        /* lower level */

   Plot_line (2, x, y);

   return 0;
}

/*}}}*/

static int get_levels_from_lines (int *lev, int *nplotlev, int nlevels, /*{{{*/
                                  int *line_idx, int nlines, int Z, int q,
                                  DB_t *db)
{
   int m, k, nlev;
   int *tlev = NULL;

   if (NULL == (tlev = (int *) ISIS_MALLOC ((nlevels+1) * sizeof(int))))
     return -1;

   for (k=0; k <= nlevels; k++)
     tlev[k] = -1;

   for (m=0; m < nlines; m++)
     {
        int up, lo;
        DB_line_t *line;

        if (NULL == (line = DB_get_line_from_index (line_idx[m], db)))
          continue;

        if ((line->proton_number != Z)
            || (line->ion_charge != q))
          continue;

        up = line->upper_level;
        lo = line->lower_level;

        tlev[up] = up;
        tlev[lo] = lo;
     }

   nlev = 0;
   for (k=1; k <= nlevels; k++)
     {
        if (tlev[k] >= 0)
          lev[nlev++] = tlev[k];
     }

   *nplotlev = nlev;

   ISIS_FREE (tlev);

   return 0;
}

/*}}}*/

/* energy level's X plot coordinate scales with (L, S) quantum numbers */

#define PLOT_X(scale, L, S)    (2*(S) + (L)/(scale))
#define WIDTH(scale)           (0.9 / (scale))

int DB_plot_levels (int Z, int q, int *line_idx, int nlines, int subset, /*{{{*/
                    int overlay, Plot_t *fmt, DB_t *db)
{
   static const char *symbol = "SPDFGHIJKLMNOQRTUVWXY";  /* think that's enough? */
   DB_ion_t *ion = NULL;
   float *energy, *x, *spin;
   int *orb, *lev;
   float de, emin, emax, xmin, xmax, pxmin, pxmax, scale;
   int k, nlevels, nplotlev;
   int ret = -1;
   char name[3];

   if (NULL == db || NULL == fmt)
     return -1;

   energy = x = spin = NULL;
   orb = lev = NULL;

   xmin = emin = 1.e30;
   xmax = emax = scale = 0.0;

   if (-1 == _DB_get_element_name (name, Z)
       || NULL == (ion = DB_get_ion (db, Z, q))
       || -1 == DB_get_ion_nlevels (&nlevels, ion)
       || nlevels <= 0)
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "for Z=%d, q=%d", Z, q);
        return -1;
     }

   if (NULL == (lev = (int *) ISIS_MALLOC ((nlevels+1) * sizeof(int))))
     return -1;

   if (subset)
     {
        if (nlines <= 0
            || -1 == get_levels_from_lines (lev, &nplotlev, nlevels, line_idx, nlines, Z, q, db))
          goto fail;

        if (nplotlev <= 0)
          {
             isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "%s %d+ level data not found", name, q);
             ret = 0;
             goto fail;
          }
     }
   else
     {
        nplotlev = nlevels;
        for (k=0; k < nlevels; k++)
          lev[k] = k+1;                    /* ground = 1 */
     }

   if (NULL == (energy = (float *) ISIS_MALLOC (nplotlev * sizeof(float)))
       || NULL == (x = (float *) ISIS_MALLOC (nplotlev * sizeof(float)))
       || NULL == (spin = (float *) ISIS_MALLOC (nplotlev * sizeof(float)))
       || NULL == (orb = (int *) ISIS_MALLOC (nplotlev * sizeof(int))))
     goto fail;

   for (k=0; k < nplotlev; k++)
     {
        DB_level_t *e;

        if (NULL == (e = DB_get_ion_level (lev[k], ion)))
          goto fail;

        energy[k] = e->energy;
        orb[k] = e->L;
        spin[k] = e->S;

        scale = (float) MAX(scale, orb[k]);

        emin = MIN(emin, energy[k]);
        emax = MAX(emax, energy[k]);
     }

   Plot_set_charsize (fmt);

   if (overlay)
     Plot_get_scale_for_elev_x_axis (&scale, fmt);
   else
     {
        Plot_next_page (Plot_this_window(fmt));
        scale = MAX(DEFAULT_ELEV_SCALE, scale+1.0);
        Plot_set_scale_for_elev_x_axis (fmt, scale);
     }

   for (k=0; k < nplotlev; k++)
     {
        if (orb[k] < 0.0 || spin[k] < 0.0)
          x[k] = -2.0/scale;
        else x[k] = PLOT_X(scale, orb[k], spin[k]);
        xmin = MIN(xmin, x[k]);
        xmax = MAX(xmax, x[k]);
     }

   if (overlay)
     {
        float *_pxmin;
        float *_pxmax;
        if ((-1 == Plot_get_xrange (fmt, &_pxmin, &_pxmax))
            || (_pxmin == NULL) || (_pxmax == NULL))
          goto fail;
        pxmin = *_pxmin;
        pxmax = *_pxmax;
        xmin = pxmin + WIDTH(scale);
     }
   else
     {
        pxmin = xmin - WIDTH(scale);
        pxmax = xmax + WIDTH(scale);

        de = 0.05 * (emax - emin);
        emin -= de;
        emax += de;

        if (-1 == Plot_set_xrange (fmt, &pxmin, &pxmax)
            || -1 == Plot_set_yrange (fmt, &emin, &emax))
          goto fail;

        (void) Plot_xtick_labels_off (fmt);
        Plot_draw_box_opt (fmt, pxmin, pxmax, emin, emax);
        (void) Plot_set_axis_options (fmt, NULL, "");
     }

   Plot_set_linestyle (fmt);

   for (k=0; k < nplotlev; k++)
     {
        float f;
        char text[16];

        Plot_set_color (fmt);
        Plot_set_line_width (fmt);
        (void) plot_level (x[k], energy[k], (float) WIDTH(scale));

        Plot_set_frame_line_width (fmt);
        f = (x[k] - pxmin) / (pxmax - pxmin);
        if (x[k] < 0)
          {
             (void) sprintf (text, "(X)");
          }
        else
          {
             (void) sprintf (text, "\\u%d\\d%c", (int) (2*spin[k] + 1), symbol[ orb[k] ]);
          }
        Plot_LS_label (f, text);
     }

   Plot_set_frame_line_width (fmt);
   Plot_label_box (fmt);
   ret = 0;

   fail:
   ISIS_FREE (lev);
   ISIS_FREE (x);
   ISIS_FREE (energy);
   ISIS_FREE (spin);
   ISIS_FREE (orb);

   return ret;
}

/*}}}*/

int DB_plot_transitions (int Z, int q, int *line_idx, int nlines, /*{{{*/
                         int *use_style, Plot_t *fmt, DB_t *db)
{
   DB_ion_t *ion = NULL;
   int m;
   float scale;
   char name[3];

   if (NULL == db || NULL == fmt)
     return -1;

   if (-1 == _DB_get_element_name (name, Z)
       || NULL == (ion = DB_get_ion (db, Z, q)))
     {
        isis_vmesg (FAIL, I_NO_DATA, __FILE__, __LINE__, "for element Z=%d, q=%d", Z, q);
        return -1;
     }

   Plot_set_linestyle (fmt);
   Plot_set_line_width (fmt);
   Plot_set_style (fmt, use_style);

   Plot_get_scale_for_elev_x_axis (&scale, fmt);

   for (m=0; m < nlines; m++)
     {
        DB_line_t *line;
        DB_level_t *up, *lo;
        float x_up, x_lo;

        if (NULL == (line = DB_get_line_from_index (line_idx[m], db)))
          continue;

        if ((line->proton_number != Z)
            || (line->ion_charge != q))
          continue;

        if (NULL == (up = DB_get_ion_level (line->upper_level, ion))
            || up->L < 0
            || up->S < 0.0)
          continue;

        if (NULL == (lo = DB_get_ion_level (line->lower_level, ion))
            || lo->L < 0
            || lo->S < 0.0)
          continue;

        x_up = PLOT_X(scale, up->L, up->S);
        x_lo = PLOT_X(scale, lo->L, lo->S);

        if (-1 == plot_transition (x_up, up->energy, x_lo, lo->energy))
          return -1;
     }

   return 0;
}

/*}}}*/
