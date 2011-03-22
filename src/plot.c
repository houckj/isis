/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2011  Massachusetts Institute of Technology

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

/* $Id: plot.c,v 1.22 2004/08/28 13:03:52 houck Exp $ */

/*{{{ includes  */

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>

#include "isis.h"
#include "plot.h"
#include "util.h"
#include "errors.h"

/*}}}*/

Plot_Options_Type Plot_User_Defaults;

/*{{{ Private Data Type Definitions  */

enum
{
   PLOT_STD_WINDOW = 0
};

#define PLOT_LIST_HEAD  (-1)

#define INIT_MIN   (FLT_MAX)
#define INIT_MAX   (-FLT_MAX)

#define INVALID_LOG10  (-100.0)
#define SAFE_LOG10(t)  (((t) > 0.0) ? log10(t) : INVALID_LOG10)

static int Max_Colors = MAX_COLORS;
static int Max_Lines = MAX_LINES;

static int Plot_Use_Auto_Increment = PLOT_DEFAULT_AUTO_INCR;

struct Plot_t
{
/* NOTICE:  Update Plot_copy_formats () if additional pointer
 *          arguments are added to Plot_. */

   /* axis limits (if NULL, use data) */
   float *xmin, *xmax, *ymin, *ymax;

   float xtick, ytick;
   int nxsub, nysub;

   float elev_plot_axis_scale;

   float exposure_time;   /* this is an ugly hack.. */

   int logx, logy;
   int linestyle, pointstyle, start_linestyle;
   float char_height, point_size, ebar_term_length;
   int color, start_color;
   int line_width, frame_line_width;
   int style_means_color;
   int rendered;		/* indicates that at least 1 curve was drawn */
   int use_errorbars;
   int connect_points;
   int use_bin_density;
   int x_unit, plotted_x_unit, ok_to_convert_xrange_units;
   int x_unit_set_by_user;
   int axis;

   int label_type;

   char xlabel[PLOT_LABEL_STRING_SIZE];
   char ylabel[PLOT_LABEL_STRING_SIZE];
   char tlabel[PLOT_LABEL_STRING_SIZE];

   /* axis format option strings */
   char *xopt;
   char *yopt;

   Plot_Node_t *this_window;
};

/* normalized device coordinates of panel viewport
 * used only if (Plot_Node_t->window_type != PLOT_STD_WINDOW) */
typedef struct _Plot_Pane_Type Plot_Pane_Type;
struct _Plot_Pane_Type
{
   float vp_xmin;
   float vp_xmax;
   float vp_ymin;
   float vp_ymax;
};
static Plot_Pane_Type Default_Outer_Viewport =
{
   0.075765, 0.924235,  /* xmin, xmax */
   0.1, 0.9             /* ymin, ymax */
};

struct Plot_Node_t
{
   Plot_Node_t *next;
   Plot_t *fmt;
   Plot_Pane_Type *panes;
   Plot_Pane_Type outer_viewport;
   int *ysizes;
   int nxpanes, nypanes;
   int window_type;
   int id;
   int this_pane;
   int force_clear;
};

/*}}}*/

/*{{{ initialize and manipulate format list */

void Plot_free_options (Plot_Options_Type *p)
{
   if (p == NULL)
     return;
   ISIS_FREE(p->xlabel);
   ISIS_FREE(p->ylabel);
   ISIS_FREE(p->tlabel);
   ISIS_FREE(p->xopt);
   ISIS_FREE(p->yopt);
}

void Plot_init_default_format (void) /*{{{*/
{
   Plot_Options_Type *p = &Plot_User_Defaults;
   Plot_Pane_Type *ovp = &Default_Outer_Viewport;

   p->ovp_xmin = ovp->vp_xmin;
   p->ovp_xmax = ovp->vp_xmax;
   p->ovp_ymin = ovp->vp_ymin;
   p->ovp_ymax = ovp->vp_ymax;

   p->logx = 0;
   p->logy = 0;
   p->pointstyle = PLOT_POINT;
   p->point_size = 1.0;
   p->linestyle = PLOT_SOLID_LINE;
   p->start_linestyle = PLOT_SOLID_LINE;
   p->char_height = 1.0;
   p->ebar_term_length = 0.0;
   p->color = PLOT_LINE_COLOR;
   p->start_color = PLOT_LINE_COLOR;
   p->line_width = 1;
   p->frame_line_width = 1;
   p->use_errorbars = 0;
   p->connect_points = 1;
   p->use_bin_density = 0;
   p->x_unit = U_ANGSTROM;
}

/*}}}*/

int Plot_set_options (Plot_t *fmt, Plot_Options_Type *p) /*{{{*/
{
   Plot_Node_t *w;
   Plot_Pane_Type *ovp;
   float *pmin, *pmax;

   if (p == NULL || fmt == NULL)
     return -1;

   Plot_set_logx (fmt, p->logx);
   Plot_set_logy (fmt, p->logy);

   pmin = (p->xmin == -FLT_MAX) ? NULL : &p->xmin;
   pmax = (p->xmax ==  FLT_MAX) ? NULL : &p->xmax;
   if (-1 == Plot_set_xrange (fmt, pmin, pmax))
     return -1;

   pmin = (p->ymin == -FLT_MAX) ? NULL : &p->ymin;
   pmax = (p->ymax ==  FLT_MAX) ? NULL : &p->ymax;
   if (-1 == Plot_set_yrange (fmt, pmin, pmax))
     return -1;

   w = fmt->this_window;
   ovp = w ? &w->outer_viewport : &Default_Outer_Viewport;
   ovp->vp_xmin = p->ovp_xmin;
   ovp->vp_xmax = p->ovp_xmax;
   ovp->vp_ymin = p->ovp_ymin;
   ovp->vp_ymax = p->ovp_ymax;

   fmt->xtick = p->xtick;
   fmt->ytick = p->ytick;

   fmt->color = p->color;
   fmt->start_color = p->start_color;

   fmt->linestyle = p->linestyle;
   fmt->start_linestyle = p->start_linestyle;

   fmt->line_width = p->line_width;
   fmt->frame_line_width = p->frame_line_width;
   fmt->pointstyle = p->pointstyle;
   fmt->connect_points = p->connect_points;
   fmt->char_height = p->char_height;
   fmt->point_size = p->point_size;
   fmt->ebar_term_length = p->ebar_term_length;
   fmt->use_errorbars = p->use_errorbars;
   fmt->use_bin_density = p->use_bin_density;
   fmt->x_unit = p->x_unit;

#define COPY_STRING(field,size) \
   if (p->field != NULL) \
     isis_strcpy (fmt->field, p->field, size);

   COPY_STRING(xlabel,PLOT_LABEL_STRING_SIZE);
   COPY_STRING(ylabel,PLOT_LABEL_STRING_SIZE);
   COPY_STRING(tlabel,PLOT_LABEL_STRING_SIZE);
#undef COPY_STRING

   if (p->xopt != NULL)
     {
        ISIS_FREE(fmt->xopt);
        if (NULL == (fmt->xopt = isis_make_string (p->xopt)))
          return -1;
     }

   if (p->yopt != NULL)
     {
        ISIS_FREE(fmt->yopt);
        if (NULL == (fmt->yopt = isis_make_string (p->yopt)))
          return -1;
     }

   return 0;
}

/*}}}*/

int Plot_get_options (Plot_t *fmt, Plot_Options_Type *p) /*{{{*/
{
   Plot_Node_t *w;
   Plot_Pane_Type *ovp;
   float *pmin=NULL, *pmax=NULL;
   char *s;

   if (p == NULL || fmt == NULL)
     return -1;

   if (-1 == Plot_get_xrange (fmt, &pmin, &pmax))
     return -1;

   if (pmin == NULL)
     p->xmin = -FLT_MAX;
   else if (fmt->logx == 0)
     p->xmin = *pmin;
   else p->xmin = pow(10.0,*pmin);

   if (pmax == NULL)
     p->xmax = FLT_MAX;
   else if (fmt->logx == 0)
     p->xmax = *pmax;
   else p->xmax = pow(10.0,*pmax);

   if (-1 == Plot_get_yrange (fmt, &pmin, &pmax))
     return -1;

   if (pmin == NULL)
     p->ymin = -FLT_MAX;
   else if (fmt->logy == 0)
     p->ymin = *pmin;
   else p->ymin = pow(10.0,*pmin);

   if (pmax == NULL)
     p->ymax = FLT_MAX;
   else if (fmt->logy == 0)
     p->ymax = *pmax;
   else p->ymax = pow(10.0,*pmax);

   w = fmt->this_window;
   ovp = w ? &w->outer_viewport : &Default_Outer_Viewport;
   p->ovp_xmin = ovp->vp_xmin;
   p->ovp_xmax = ovp->vp_xmax;
   p->ovp_ymin = ovp->vp_ymin;
   p->ovp_ymax = ovp->vp_ymax;

   p->xtick = fmt->xtick;
   p->ytick = fmt->ytick;
   p->logx = fmt->logx;
   p->logy = fmt->logy;
   p->color = fmt->color;
   p->start_color = fmt->start_color;
   p->linestyle = fmt->linestyle;
   p->start_linestyle = fmt->start_linestyle;
   p->line_width = fmt->line_width;
   p->frame_line_width = fmt->frame_line_width;
   p->pointstyle = fmt->pointstyle;
   p->connect_points = fmt->connect_points;
   p->char_height = fmt->char_height;
   p->point_size = fmt->point_size;
   p->ebar_term_length = fmt->ebar_term_length;
   p->use_errorbars = fmt->use_errorbars;
   p->use_bin_density = fmt->use_bin_density;
   p->x_unit = fmt->x_unit;

#define DUP_STRING_OR_RETURN_ERROR(field) \
   if (NULL == (s = isis_make_string (fmt->field))) \
     return -1; \
   ISIS_FREE(p->field); \
   p->field = s;

   DUP_STRING_OR_RETURN_ERROR(xlabel);
   DUP_STRING_OR_RETURN_ERROR(ylabel);
   DUP_STRING_OR_RETURN_ERROR(tlabel);
   DUP_STRING_OR_RETURN_ERROR(xopt);
   DUP_STRING_OR_RETURN_ERROR(yopt);

#undef DUP_STRING_OR_RETURN_ERROR

   return 0;
}

/*}}}*/

static int Plot_set_default_fmt (Plot_t *fmt) /*{{{*/
{
   char *s;
   if (-1 == Plot_set_options (fmt, &Plot_User_Defaults))
     return -1;

   fmt->plotted_x_unit = fmt->x_unit;
   fmt->start_linestyle = fmt->linestyle;
   fmt->start_color = fmt->color;

   if (NULL == (s = Plot_get_option_string_default ()))
     return -1;
   ISIS_FREE(fmt->xopt);
   fmt->xopt = s;

   if (NULL == (s = Plot_get_option_string_default ()))
     return -1;
   ISIS_FREE(fmt->yopt);
   fmt->yopt = s;

   ISIS_FREE(fmt->xmin);
   ISIS_FREE(fmt->xmax);
   ISIS_FREE(fmt->ymin);
   ISIS_FREE(fmt->ymax);

   fmt->xtick = fmt->ytick = 0.0;
   fmt->nxsub = fmt->nysub = 0;
   fmt->elev_plot_axis_scale = DEFAULT_ELEV_SCALE;
   fmt->exposure_time = 1.0;
   fmt->style_means_color = 1;
   fmt->ok_to_convert_xrange_units = 0;
   fmt->label_type = 0;
   fmt->xlabel[0] = '\0';
   fmt->ylabel[0] = '\0';
   fmt->tlabel[0] = '\0';
   fmt->axis = 0;

   return 0;
}
/*}}}*/

static Plot_t *new_fmt (void) /*{{{*/
{
   Plot_t *p;

   if (NULL == (p = (Plot_t *) ISIS_MALLOC (sizeof(Plot_t))))
     return NULL;
   memset ((char *)p, 0, sizeof (*p));

   if (-1 == Plot_set_default_fmt (p))
     return NULL;

   return p;
}

/*}}}*/

static void free_fmt (Plot_t *fmt) /*{{{*/
{
   if (fmt == NULL)
     return;

   ISIS_FREE (fmt->xmin);
   ISIS_FREE (fmt->xmax);
   ISIS_FREE (fmt->ymin);
   ISIS_FREE (fmt->ymax);
   ISIS_FREE (fmt->xopt);
   ISIS_FREE (fmt->yopt);
   ISIS_FREE (fmt);
}

/*}}}*/

Plot_t *Plot_alloc_fmt (void) /*{{{*/
{
   return new_fmt ();
}

/*}}}*/

void Plot_free_fmt (Plot_t *fmt) /*{{{*/
{
   free_fmt (fmt);
}

/*}}}*/

static int free_node (Plot_Node_t *t) /*{{{*/
{
   if (t == NULL)
     return -1;

   if ((t->id > 0) && (t->id != PLOT_NO_DEVICE))
     {
        if (-1 == Plot_select_window (t->id))
          return -1;
        if (-1 == Plot_close ())
          return -1;
     }

   free_fmt (t->fmt);
   ISIS_FREE (t->ysizes);
   ISIS_FREE (t->panes);
   ISIS_FREE (t);

   return 0;
}

/*}}}*/

static int define_panel_viewports (Plot_Node_t *t, int nxpanes, int nypanes) /*{{{*/
{
   Plot_Pane_Type *pane;
   Plot_Pane_Type vp;
   int i, sum;
   float dy, ymin;

   if ((NULL == t) || (NULL == t->ysizes))
     return -1;

   if (nxpanes != 1)
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot style requires nxpanes=1");
        return -1;
     }

   if (nypanes == 1)
     {
        /* struct copy */
        *t->panes = t->outer_viewport;
        return 0;
     }

   sum = 0;
   for (i = 0; i < nypanes; i++)
     sum += t->ysizes[i];

   /* struct copy */
   vp = t->outer_viewport;

   dy = (vp.vp_ymax - vp.vp_ymin) / sum;
   ymin = vp.vp_ymin;

   pane = t->panes + nypanes - 1;

   for (i = nypanes-1; i >= 0; i--)
     {
        pane->vp_xmin = vp.vp_xmin;
        pane->vp_xmax = pane->vp_xmin + (vp.vp_xmax - vp.vp_xmin);

        pane->vp_ymin = ymin;
        pane->vp_ymax = pane->vp_ymin + t->ysizes[i] * dy;

        ymin = pane->vp_ymax;

        pane--;
     }

   return 0;
}

/*}}}*/

static Plot_Node_t *alloc_node (void) /*{{{*/
{
   Plot_Node_t *t = NULL;

   if (NULL == (t = (Plot_Node_t *) ISIS_MALLOC (sizeof(Plot_Node_t))))
     return NULL;
   memset ((char *)t, 0, sizeof (*t));

   t->id = PLOT_NO_DEVICE;
   t->nxpanes = t->nypanes = 0;
   t->next = NULL;
   t->fmt = NULL;
   t->ysizes = NULL;
   t->panes = NULL;
   t->this_pane = 0;
   t->force_clear = 0;
   /* struct copy */
   t->outer_viewport = Default_Outer_Viewport;

   return t;
}

/*}}}*/

static int realloc_node_formats (Plot_Node_t *t, int force, int nxpanes, int nypanes, int type, int *ysizes) /*{{{*/
{
   int npanes = nxpanes * nypanes;

   if (t == NULL)
     return -1;

   /* Don't reallocate unless we really have to. */
   if ((force == 0)
       && (t->nxpanes == nxpanes)
       && (t->nypanes == nypanes)
       && (t->ysizes != NULL) && (ysizes != NULL))
     {
        int i, same = 1;
        for (i = 0; i < npanes; i++)
          {
             if (t->ysizes[i] == ysizes[i])
               continue;
             same = 0;
             break;
          }
        if (same == 1) return 0;
     }

   t->force_clear = 1;
   t->this_pane = 0;
   t->nxpanes = nxpanes;
   t->nypanes = nypanes;
   t->window_type = type;

   if (t->panes != NULL) ISIS_FREE (t->panes);
   if (t->ysizes != NULL) ISIS_FREE (t->ysizes);

   if (t->window_type != PLOT_STD_WINDOW)
     {
        if ((NULL == ysizes)
            || (NULL == (t->ysizes = (int *) ISIS_MALLOC (npanes * sizeof(int))))
            || (NULL == (t->panes = (Plot_Pane_Type *) ISIS_MALLOC (npanes * sizeof(Plot_Pane_Type)))))
          return -1;

        memcpy ((char *)t->ysizes, (char *)ysizes, npanes * sizeof(int));

        if (-1 == define_panel_viewports (t, t->nxpanes, t->nypanes))
          return -1;
     }

   /* If we've got a format structure, leave it alone. */
   if (t->fmt == NULL)
     {
        t->fmt = new_fmt ();
        if (NULL == t->fmt)
          return -1;
        t->fmt->this_window = t;
     }

   return 0;
}

/*}}}*/

int Plot_set_outer_viewport (Plot_Node_t *n, float xmin, float xmax, float ymin, float ymax) /*{{{*/
{
   Plot_Pane_Type *vp;
   Plot_Node_t *t;
   int *ysizes;
   int status;

   t = n->fmt->this_window;
   vp = &t->outer_viewport;

   vp->vp_xmin = xmin;
   vp->vp_xmax = xmax;
   vp->vp_ymin = ymin;
   vp->vp_ymax = ymax;

   if (n->window_type == PLOT_STD_WINDOW)
     return 0;

   if (NULL == (ysizes = (int *) ISIS_MALLOC (t->nypanes * sizeof (int))))
     return -1;

   if (t->nypanes == 1)
     {
        ysizes[0] = 1;
     }
   else if (t->ysizes != NULL)
     {
        memcpy ((char *)ysizes, t->ysizes, t->nypanes * sizeof(int));
     }
   else
     {
        ISIS_FREE(ysizes);
        fprintf (stderr, "*** Error:  Plot_set_outer_viewport:  this should never happen\n");
        return -1;
     }

   status = realloc_node_formats (t, 1, t->nxpanes, t->nypanes, t->window_type, ysizes);

   ISIS_FREE(ysizes);
   return status;
}

/*}}}*/

void Plot_get_outer_viewport (Plot_Node_t *n, float *xmin, float *xmax, float *ymin, float *ymax) /*{{{*/
{
   Plot_Pane_Type *vp = &Default_Outer_Viewport;

   if (n != NULL)
     {
        Plot_Node_t *t= n->fmt->this_window;
        vp = &t->outer_viewport;
     }

   *xmin = vp->vp_xmin;
   *xmax = vp->vp_xmax;
   *ymin = vp->vp_ymin;
   *ymax = vp->vp_ymax;
}

/*}}}*/

Plot_Node_t *Plot_start (void) /*{{{*/
{
   Plot_Node_t *n;

   if (NULL != (n = alloc_node ()))
     n->id = PLOT_LIST_HEAD;

   return n;
}

/*}}}*/

void Plot_end (Plot_Node_t **h) /*{{{*/
{
   Plot_Node_t *nh;

   while (*h != NULL)
     {
        nh = (*h)->next;
        free_node (*h);
        *h = nh;
     }
}

/*}}}*/

static int delete_after_node (Plot_Node_t *n) /*{{{*/
{
   Plot_Node_t *dead;

   if (n == NULL)
     return -1;

   dead = n->next;
   n->next = dead->next;
   free_node (dead);

   return 0;
}

/*}}}*/

static int append_node (Plot_Node_t *head, Plot_Node_t *n) /*{{{*/
{
   Plot_Node_t *t;

   for (t=head; t->next != NULL; t=t->next)
     ;

   t->next = n;
   n->next = NULL;

   return 0;
}

/*}}}*/

Plot_Node_t *Plot_new_window (Plot_Node_t *head) /*{{{*/
{
   Plot_Node_t *n;

   if (NULL == (n = alloc_node ()))
     return NULL;

   if (-1 == append_node (head, n))
     {
        free_node (n);
        return NULL;
     }

   return n;
}

/*}}}*/

/*}}}*/

int Plot_reformat_node (Plot_Node_t *n, int nxpanes, int nypanes, int type, int *ysizes) /*{{{*/
{
   if (n == NULL)
     return -1;

   return realloc_node_formats (n, 0, nxpanes, nypanes, type, ysizes);
}

/*}}}*/

int Plot_open_window (Plot_Node_t *n, char *device) /*{{{*/
{
   if (n == NULL)
     return -1;
   n->id = Plot_open (device);
   return (n->id <= 0) ? -1 : n->id;
}

/*}}}*/

static float *copy_float (float *f) /*{{{*/
{
   float *c = NULL;

   if (f != NULL)
     {
        if (NULL != (c = (float *) ISIS_MALLOC (sizeof(float))))
          *c = *f;
     }

   return c;
}

/*}}}*/

int Plot_copy_formats (Plot_t *src, Plot_t *dst) /*{{{*/
{
   Plot_Node_t *src_this;
   Plot_Node_t *dst_this;

   if ((src == NULL) || (dst == NULL))
     return -1;

   /* preserve "this" pointers */
   src_this = src->this_window;
   dst_this = dst->this_window;

   ISIS_FREE (dst->xmin);  ISIS_FREE (dst->xmax);
   ISIS_FREE (dst->ymin);  ISIS_FREE (dst->ymax);
   ISIS_FREE (dst->xopt);  ISIS_FREE (dst->yopt);

   memcpy ((char *) dst, (char *) src, sizeof (*dst));

   dst->xmin = copy_float (src->xmin);
   dst->xmax = copy_float (src->xmax);
   dst->ymin = copy_float (src->ymin);
   dst->ymax = copy_float (src->ymax);

   dst->xopt = isis_make_string (src->xopt);
   dst->yopt = isis_make_string (src->yopt);

   dst->rendered = 0;             /* never inherit rendered state */

   src->this_window = src_this;
   dst->this_window = dst_this;

   return 0;
}

/*}}}*/

int Plot_copy_window (char *new_device, int id, Plot_Node_t *head) /*{{{*/
{
   Plot_Node_t *t, *n;
   int npanes;

   for (t=head; t != NULL; t=t->next)
     {
        if (id == t->id)
          goto found;
     }

   return -1;

   found:

   if (NULL == (n = Plot_new_window (head)))
     return -1;

   if (-1 == realloc_node_formats (n, 0, t->nxpanes, t->nypanes, t->window_type, t->ysizes))
     return -1;

   if (-1 == Plot_open_window (n, new_device))
     return -1;

   if (n->window_type == PLOT_STD_WINDOW)
     {
        if (-1 == Plot_subdivide (n->nxpanes, n->nypanes))
          return -1;
     }

   npanes = t->nxpanes * t->nypanes;

   if (-1 == Plot_copy_formats (t->fmt, n->fmt))
     return -1;

   return n->id;
}

/*}}}*/

int Plot_delete_window (int id, Plot_Node_t *head) /*{{{*/
{
   Plot_Node_t *t;

   if ((head == NULL) || (id == PLOT_LIST_HEAD))
     return -1;

   if (id == PLOT_NO_DEVICE)
     return 0;

   for (t=head; t->next != NULL; t=t->next)
     {
        if (id == t->next->id)
          return delete_after_node (t);
     }

   return -1;
}
/*}}}*/

Plot_Node_t *Plot_find_node (int id, Plot_Node_t *head) /*{{{*/
{
   Plot_Node_t *n;

   if (head == NULL || id == PLOT_LIST_HEAD)
     return NULL;

   for (n = head->next; n != NULL; n = n->next)
     {
        if (n->id == id)
          return n;
     }

   return NULL;
}

/*}}}*/

Plot_t *Plot_fmt (Plot_Node_t *n) /*{{{*/
{
   return n ? n->fmt : NULL;
}

/*}}}*/

Plot_Node_t *Plot_this_window (Plot_t *fmt)
{
   return fmt ? fmt->this_window : NULL;
}

/*{{{ set format fields */

int Plot_xtick_labels_off (Plot_t *fmt) /*{{{*/
{
   char *xopt;
   if (NULL == (xopt = Plot_configure_axis (fmt->xopt, 0, 0)))
     return -1;
   ISIS_FREE(fmt->xopt);
   fmt->xopt = xopt;
   return 0;
}

/*}}}*/

void Plot_set_axis_type (Plot_t *fmt, int axis) /*{{{*/
{
   fmt->axis = axis;
}

/*}}}*/

int Plot_set_axis_options (Plot_t *fmt, const char *xopt, const char *yopt) /*{{{*/
{
   char *s = NULL;

   /* [xy]opt = NULL means use default
    * [xy]opt = empty string means leave as-is */
   if (xopt == NULL)
     {
        ISIS_FREE(fmt->xopt);
        if (NULL == (fmt->xopt = Plot_get_option_string_default()))
          return -1;
     }
   else if (*xopt != 0)
     {
        if (NULL == (s = isis_make_string (xopt)))
          return -1;
        ISIS_FREE(fmt->xopt);
        fmt->xopt = s;
     }

   if (yopt == NULL)
     {
        ISIS_FREE(fmt->yopt);
        if (NULL == (fmt->yopt = Plot_get_option_string_default()))
          return -1;
     }
   else if (*yopt != 0)
     {
        if (NULL == (s = isis_make_string (yopt)))
          return -1;
        ISIS_FREE(fmt->yopt);
        fmt->yopt = s;
     }

   return 0;
}

/*}}}*/

void Plot_set_scale_for_elev_x_axis (Plot_t *fmt, float scale) /*{{{*/
{
   fmt->elev_plot_axis_scale = scale;
}

/*}}}*/

void Plot_get_scale_for_elev_x_axis (float *scale, Plot_t *fmt) /*{{{*/
{
   *scale = fmt->elev_plot_axis_scale;
}

/*}}}*/

int Plot_save_exposure_time (Plot_t *fmt, float t) /*{{{*/
{
   if (fmt == NULL)
     return -1;

   fmt->exposure_time = t;

   return 0;
}

/*}}}*/

float Plot_get_exposure_time (Plot_t *fmt) /*{{{*/
{
   if (fmt == NULL)
     return 1.0;

   return fmt->exposure_time;
}

/*}}}*/

static int _restart_style_cycle (Plot_t *fmt) /*{{{*/
{
   fmt->color = fmt->start_color;
   fmt->linestyle = fmt->start_linestyle;
   return 0;
}

/*}}}*/

int Plot_restart_style_cycle (Plot_t *fmt)
{
   return _restart_style_cycle (fmt);
}

int Plot_set_fmt_color (Plot_t *fmt, int color) /*{{{*/
{
   fmt->color = color;
   fmt->start_color = color;
   return 0;
}

/*}}}*/

int Plot_get_fmt_color (Plot_t *fmt) /*{{{*/
{
   return fmt->color;
}

/*}}}*/

int Plot_set_fmt_linestyle (Plot_t *fmt, int linestyle) /*{{{*/
{
   fmt->linestyle = linestyle;
   fmt->start_linestyle = linestyle;
   return 0;
}

/*}}}*/

int Plot_set_fmt_charsize (Plot_t *fmt, float char_height) /*{{{*/
{
   fmt->char_height = char_height;
   return 0;
}

/*}}}*/

void Plot_set_auto_increment (int flag) /*{{{*/
{
   Plot_Use_Auto_Increment = (flag == 0) ? 0 : 1;
   isis_vmesg (INFO, I_INFO, __FILE__, __LINE__,
               "overplot style auto-select %s",
               Plot_Use_Auto_Increment ? "ON" : "OFF");
}

/*}}}*/

static void Plot_auto_incr_linestyle (Plot_t *fmt) /*{{{*/
{
   if (!Plot_Use_Auto_Increment)
     return;
   fmt->linestyle = 1 + fmt->linestyle % Max_Lines;
}

/*}}}*/

static void Plot_auto_incr_line_color (Plot_t *fmt) /*{{{*/
{
   if (!Plot_Use_Auto_Increment)
     return;
   fmt->color = 1 + fmt->color % Max_Colors;
}

/*}}}*/

void Plot_set_style_meaning (Plot_t *fmt, int style) /*{{{*/
{
   fmt->style_means_color = style ? 1 : 0;
   isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "overplots use different line %s",
               fmt->style_means_color ? "COLORS" : "STYLES");
}

/*}}}*/

void Plot_auto_incr_line_type (Plot_t *fmt) /*{{{*/
{
   if (fmt->style_means_color)
     Plot_auto_incr_line_color (fmt);
   else
     Plot_auto_incr_linestyle (fmt);
}

/*}}}*/

void Plot_set_style (Plot_t *fmt, int *style) /*{{{*/
{
   if (fmt->style_means_color)
     {
        if (style == NULL)
          Plot_set_color (fmt);
        else
          _Plot_set_color (*style);
     }
   else
     {
        if (style == NULL)
          Plot_set_linestyle (fmt);
        else
          _Plot_set_linestyle (*style);
     }
}

/*}}}*/

int Plot_set_fmt_pointstyle (Plot_t *fmt, int pointstyle) /*{{{*/
{
   fmt->pointstyle = pointstyle;
   return 0;
}

/*}}}*/

static int Plot_get_pointstyle (Plot_t *fmt) /*{{{*/
{
   return fmt->pointstyle;
}

/*}}}*/

int Plot_set_fmt_point_size (Plot_t *fmt, float point_size) /*{{{*/
{
   fmt->point_size = point_size;
   return 0;
}

/*}}}*/

static float Plot_get_point_size (Plot_t *fmt) /*{{{*/
{
   return fmt->point_size;
}

/*}}}*/

int Plot_set_fmt_line_width (Plot_t *fmt, int width) /*{{{*/
{
   fmt->line_width = width;
   return 0;
}

/*}}}*/

int Plot_get_fmt_line_width (Plot_t *fmt) /*{{{*/
{
   return fmt->line_width;
}

/*}}}*/

int Plot_set_fmt_frame_line_width (Plot_t *fmt, int width) /*{{{*/
{
   fmt->frame_line_width = width;
   return 0;
}

/*}}}*/

void Plot_set_bin_density (Plot_t *fmt, int type) /*{{{*/
{
   type = (type == 0) ? 0 : 1;
#if 0
   if (fmt->use_bin_density != type)
     (void) _restart_style_cycle (fmt);
#endif
   fmt->use_bin_density = type;
}

/*}}}*/

int Plot_bin_density (Plot_t *fmt) /*{{{*/
{
   return fmt->use_bin_density;
}

/*}}}*/

int Plot_x_unit (Plot_t *fmt) /*{{{*/
{
   return fmt->x_unit;
}

/*}}}*/

static int convert_plot_limit (float *x, int islog, int out_unit, int in_unit) /*{{{*/
{
   double xd;

   if (NULL == x)
     return 0;

   xd = (double) *x;
   if (islog)
     xd = pow (10.0, xd);

   if (-1 == unit_convert_x (&xd, 1, out_unit, in_unit))
     return -1;

   if (islog)
     xd = log10 (xd);
   *x = (float) xd;

   return 0;
}

/*}}}*/

static int convert_plot_units (Plot_t *fmt) /*{{{*/
{
   int need_unit_conversion;
   int to_unit, from_unit;

   need_unit_conversion = ((fmt->plotted_x_unit > 0)
                           && fmt->ok_to_convert_xrange_units
                           && (fmt->x_unit != fmt->plotted_x_unit));

   if (!need_unit_conversion)
     return 0;

   from_unit = fmt->plotted_x_unit;
   to_unit = fmt->x_unit;

   if (-1 == convert_plot_limit (fmt->xmin, fmt->logx, to_unit, from_unit))
     {
        isis_vmesg (FAIL, I_RANGE_ERROR, __FILE__, __LINE__, "converting plot XMIN units");
        return -1;
     }

   if (-1 == convert_plot_limit (fmt->xmax, fmt->logx, to_unit, from_unit))
     {
        isis_vmesg (FAIL, I_RANGE_ERROR, __FILE__, __LINE__, "converting plot XMAX units");
        return -1;
     }

   if ((fmt->xmax != NULL)
       && (fmt->xmin != NULL)
       && (*fmt->xmin > *fmt->xmax))
     {
        float tmp = *fmt->xmin;
        *fmt->xmin = *fmt->xmax;
        *fmt->xmax = tmp;
     }
   else if (unit_is_wavelength (to_unit) != unit_is_wavelength (from_unit))
     {
        float *p = fmt->xmin;
        fmt->xmin = fmt->xmax;
        fmt->xmax = p;
     }

   fmt->x_unit = to_unit;

   return 0;
}

/*}}}*/

int Plot_x_unit_is_default (Plot_t *fmt) /*{{{*/
{
   if (fmt == NULL)
     return 1;

   return ((0 == fmt->x_unit_set_by_user) && (fmt->x_unit == U_ANGSTROM));
}

/*}}}*/

int Plot_set_x_unit (Plot_t *fmt, int unit) /*{{{*/
{
   if (fmt == NULL)
     return -1;

   fmt->x_unit_set_by_user = 1;

   if (fmt->x_unit == unit)
     return 0;

   if (unit_invalid (unit))
     return -1;

   /* (void) _restart_style_cycle (fmt); */

   fmt->x_unit = unit;

   return 0;
}

/*}}}*/

void Plot_set_errorbars (Plot_t *fmt, int value, float term_length) /*{{{*/
{
   fmt->use_errorbars = value;
   if (term_length >= 0.0)
     fmt->ebar_term_length = term_length;
}

/*}}}*/

int Plot_use_errorbars (Plot_t *fmt) /*{{{*/
{
   return fmt->use_errorbars;
}

/*}}}*/

void Plot_set_connect_points (Plot_t *fmt, int value) /*{{{*/
{
   fmt->connect_points = value;
}

/*}}}*/

int Plot_connect_points (Plot_t *fmt) /*{{{*/
{
   return fmt->connect_points;
}

/*}}}*/

int Plot_get_num_xpanes (int *nxpanes, int id, Plot_Node_t *head) /*{{{*/
{
   Plot_Node_t *t;

   *nxpanes = 0;

   if (head == NULL || id == PLOT_LIST_HEAD)
     return -1;

   for (t=head->next; t != NULL && t->id != id; t=t->next)
     ;

   if (t == NULL)
     return -1;

   *nxpanes = t->nxpanes;
   return 0;
}

/*}}}*/

static void dolog_or_free (float **val) /*{{{*/
{
   if (val == NULL || *val == NULL)
     return;

   if (**val > 0.0)
     **val = log10 (**val);
   else
     ISIS_FREE (*val);
}

/*}}}*/

static void dopow10_or_free (float **val) /*{{{*/
{
   if (val == NULL || *val == NULL)
     return;

   if (**val < FLT_MAX_10_EXP)
     **val = pow (10.0, **val);
   else
     ISIS_FREE (*val);
}

/*}}}*/

static void update_axis_type (Plot_t *fmt) /*{{{*/
{
   Plot_Node_t *w = fmt->this_window;
   char *xopt = NULL;
   char *yopt = NULL;

   if (NULL == (yopt = Plot_configure_axis (fmt->yopt, fmt->logy, 1)))
     return;

   if ((w->this_pane < w->nypanes)
       && (w->window_type != PLOT_STD_WINDOW))
     {
        if (NULL == (xopt = Plot_configure_axis (fmt->xopt, fmt->logx, 0)))
          {
             ISIS_FREE(yopt);
             return;
          }
     }
   else
     {
        if (NULL == (xopt = Plot_configure_axis (fmt->xopt, fmt->logx, 1)))
          {
             ISIS_FREE(yopt);
             return;
          }
     }

   fmt->axis = 0;
   if (fmt->logx) fmt->axis += 10;
   if (fmt->logy) fmt->axis += 20;

   (void) Plot_set_axis_options (fmt, xopt, yopt);

   ISIS_FREE(xopt);
   ISIS_FREE(yopt);
}

/*}}}*/

void Plot_set_logx (Plot_t *fmt, int logx) /*{{{*/
{
   logx = (logx == 0) ? 0 : 1;

   if (logx == fmt->logx)
     return;

   fmt->logx = logx;
   update_axis_type (fmt);
   /* (void) _restart_style_cycle (fmt); */

   if (logx)
     {
        dolog_or_free (&fmt->xmin);
        dolog_or_free (&fmt->xmax);
     }
   else
     {
        dopow10_or_free (&fmt->xmin);
        dopow10_or_free (&fmt->xmax);
     }
}

/*}}}*/

void Plot_set_logy (Plot_t *fmt, int logy) /*{{{*/
{
   logy = (logy == 0) ? 0 : 1;

   if (logy == fmt->logy)
     return;

   fmt->logy = logy;
   update_axis_type (fmt);
   /* (void) _restart_style_cycle (fmt); */

   if (logy)
     {
        dolog_or_free (&fmt->ymin);
        dolog_or_free (&fmt->ymax);
     }
   else
     {
        dopow10_or_free (&fmt->ymin);
        dopow10_or_free (&fmt->ymax);
     }
}

/*}}}*/

int Plot_get_logx (Plot_t *fmt) /*{{{*/
{
   return fmt->logx;
}

/*}}}*/

int Plot_get_logy (Plot_t *fmt) /*{{{*/
{
   return fmt->logy;
}

/*}}}*/

int Plot_set_label_type (Plot_t *fmt, int type) /*{{{*/
{
   fmt->label_type = type;
   return 0;
}

/*}}}*/

int Plot_has_default_labels (Plot_t *fmt) /*{{{*/
{
   return (fmt->label_type == PLOT_DEFAULT_LABEL);
}

/*}}}*/

int Plot_set_title (Plot_t *fmt, const char *title) /*{{{*/
{
   isis_strcpy (fmt->tlabel, title, PLOT_LABEL_STRING_SIZE);
   return 0;
}

/*}}}*/

int Plot_set_xlabel (Plot_t *fmt, const char *xlabel) /*{{{*/
{
   isis_strcpy (fmt->xlabel, xlabel, PLOT_LABEL_STRING_SIZE);
   fmt->label_type = PLOT_USER_LABEL;
   return 0;
}

/*}}}*/

int Plot_set_ylabel (Plot_t *fmt, const char *ylabel) /*{{{*/
{
   isis_strcpy (fmt->ylabel, ylabel, PLOT_LABEL_STRING_SIZE);
   fmt->label_type = PLOT_USER_LABEL;
   return 0;
}

/*}}}*/

int Plot_set_labels (Plot_t *fmt, const char *xlabel, /*{{{*/
                     const char *ylabel, const char *tlabel)
{
   if (xlabel != NULL)
     {
        fmt->label_type = PLOT_USER_LABEL;
        isis_strcpy (fmt->xlabel, xlabel, PLOT_LABEL_STRING_SIZE);
     }
   if (ylabel != NULL)
     {
        fmt->label_type = PLOT_USER_LABEL;
        isis_strcpy (fmt->ylabel, ylabel, PLOT_LABEL_STRING_SIZE);
     }
   if (tlabel != NULL)
     isis_strcpy (fmt->tlabel, tlabel, PLOT_LABEL_STRING_SIZE);

   return 0;
}

/*}}}*/

int Plot_get_xrange (Plot_t *fmt, float **xmin, float **xmax) /*{{{*/
{
   if (fmt == NULL)
     return -1;
   if (xmin != NULL) *xmin = fmt->xmin;
   if (xmax != NULL) *xmax = fmt->xmax;
   return 0;
}

/*}}}*/

int Plot_get_yrange (Plot_t *fmt, float **ymin, float **ymax) /*{{{*/
{
   if (fmt == NULL)
     return -1;
   if (ymin != NULL) *ymin = fmt->ymin;
   if (ymax != NULL) *ymax = fmt->ymax;
   return 0;
}

/*}}}*/

static int set_ptr_or_value (float **dest, float *value, int islog) /*{{{*/
{
   if (value == NULL)
     {
        ISIS_FREE (*dest);
        return 0;
     }

   if (islog && (*value <= 0.0))
     {
        isis_vmesg (FAIL, I_RANGE_ERROR, __FILE__, __LINE__, "log scale requires values > 0");
        return -1;
     }

   if (dest == NULL)
     return -1;

   if ((NULL == *dest)
       && (NULL == (*dest = (float *) ISIS_MALLOC (sizeof(float)))))
     return -1;

   if (islog)
     *(*dest) = log10(*value);
   else
     *(*dest) = *value;

   return 0;
}

/*}}}*/

int Plot_set_xrange (Plot_t *fmt, float *xmin, float *xmax) /*{{{*/
{
   if (-1 == set_ptr_or_value (&fmt->xmin, xmin, fmt->logx)
       || -1 == set_ptr_or_value (&fmt->xmax, xmax, fmt->logx))
     return -1;

   fmt->ok_to_convert_xrange_units = 0;

   /* _restart_style_cycle (fmt); */

   return 0;
}

/*}}}*/

int Plot_set_yrange (Plot_t *fmt, float *ymin, float *ymax) /*{{{*/
{
   if (-1 == set_ptr_or_value (&fmt->ymin, ymin, fmt->logy)
       || -1 == set_ptr_or_value (&fmt->ymax, ymax, fmt->logy))
     return -1;

   /* _restart_style_cycle (fmt); */

   return 0;
}

/*}}}*/

void Plot_reset_xyranges (Plot_t *fmt) /*{{{*/
{
   ISIS_FREE (fmt->xmin); ISIS_FREE (fmt->xmax);
   ISIS_FREE (fmt->ymin); ISIS_FREE (fmt->ymax);
   /* _restart_style_cycle (fmt); */
}

/*}}}*/

/*}}}*/

static void fixup_empty_plot_range (float *min, float *max) /*{{{*/
{
   if (*min == 0.0)
     {
        *min = -1.0;
        *max =  1.0;
     }
   else
     {
        *min = *min - 0.5*fabs(*min);
        *max = *max + 0.5*fabs(*max);
     }
}

/*}}}*/

static int x_in_specified_range (float x, float *xmin, float *xmax) /*{{{*/
{
   if (xmin && xmax)
     {
        float lo, hi;
        lo = MIN(*xmin, *xmax);
        hi = MAX(*xmin, *xmax);
        return (lo <= x && x <= hi);
     }

   if ((xmin != NULL && x < *xmin)
       || (xmax != NULL && x > *xmax))
     return 0;

   return 1;
}

/*}}}*/

static int have_data_in_yrange (Plot_t *fmt, float ymin, float ymax) /*{{{*/
{
   if (fmt->ymin && fmt->ymax)
     {
        float lo, hi;
        lo = MIN(*fmt->ymin, *fmt->ymax);
        hi = MAX(*fmt->ymin, *fmt->ymax);

        return INTERVALS_OVERLAP(lo,hi,ymin,ymax);
     }

   if ((fmt->ymin && (ymax < *fmt->ymin))
       || (fmt->ymax && (*fmt->ymax < ymin)))
     return 0;

   return 1;
}

/*}}}*/

static int have_data_in_xrange (Plot_t *fmt, float xmin, float xmax) /*{{{*/
{
   if (fmt->xmin && fmt->xmax)
     {
        float lo, hi;
        lo = MIN(*fmt->xmin, *fmt->xmax);
        hi = MAX(*fmt->xmin, *fmt->xmax);

        return INTERVALS_OVERLAP(lo,hi,xmin,xmax);
     }

   if ((fmt->xmin && (xmax < *fmt->xmin))
       || (fmt->xmax && (*fmt->xmax < xmin)))
     return 0;

   return 1;
}

/*}}}*/

static int find_plot_range (float *x, int n, int logs, /*{{{*/
                            float *xmin, float *xmax)
{
   float mn, mx;
   int i;

   mn = INIT_MIN;
   mx = INIT_MAX;

   for (i = 0; i < n; i++)
     {
        float t = x[i];

        if (logs && (t <= INVALID_LOG10))
          continue;

        if (t < mn) mn = t;
        if (t > mx) mx = t;
     }

   *xmin = mn;
   *xmax = mx;

   return 0;
}

/*}}}*/

static int yrange_for_xrange (Plot_t *fmt, float *y, float *x, int n, /*{{{*/
                              float *ymin, float *ymax)
{
   float mn, mx;
   int i;

   mn = INIT_MIN;
   mx = INIT_MAX;

   for (i = 0; i < n; i++)
     {
        float tx = x[i];
        float ty;

        if (0 == x_in_specified_range (tx, fmt->xmin, fmt->xmax))
          continue;

        ty = y[i];

        if (fmt->logy && (ty <= INVALID_LOG10))
          continue;

        if (ty < mn) mn = ty;
        if (ty > mx) mx = ty;
     }

   *ymin = mn;
   *ymax = mx;

   return 0;
}

/*}}}*/

static int validate_range (float *vxmin, float *vxmax) /*{{{*/
{
   if (*vxmin == INIT_MIN || *vxmax == INIT_MAX)
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "plot limits");
        return -1;
     }

   if (*vxmin == *vxmax)
     {
        fixup_empty_plot_range (vxmin, vxmax);
     }

   return 0;
}

/*}}}*/

static int pad_range (float *vxmin, float *vxmax) /*{{{*/
{
   float dx = 0.02 * (*vxmax - *vxmin);

   *vxmin -= dx;
   *vxmax += dx;

   return 0;
}

/*}}}*/

static int compute_errorbars (float *y, float *y_err, int n, float *y_bot, float *y_top) /*{{{*/
{
   int i;

   for (i = 0; i < n; i++)
     {
        float t  = y[i];
        float dt = y_err[i];
        y_bot[i] = t - dt;
        y_top[i] = t + dt;
     }

   return 0;
}

/*}}}*/

static int compute_log10 (float *x, int n, float *xlog) /*{{{*/
{
   int i;

   for (i = 0; i < n; i++)
     {
        float t = x[i];
        xlog[i] = SAFE_LOG10(t);
     }

   return 0;
}

/*}}}*/

/*{{{ apply format to plot */

int Plot_visit_pane (Plot_Node_t *n, int pane) /*{{{*/
{
   Plot_Pane_Type *w;
   int status = -1;

   if (n == NULL)
     return status;

   if (n->window_type == PLOT_STD_WINDOW)
     {
        Plot_Pane_Type *vp = &n->outer_viewport;
        return Plot_select_viewport (vp->vp_xmin, vp->vp_xmax,
                                     vp->vp_ymin, vp->vp_ymax);
     }

   n->this_pane = pane % n->nypanes;
   w = &n->panes[n->this_pane];

   return Plot_select_viewport (w->vp_xmin, w->vp_xmax,
                                w->vp_ymin, w->vp_ymax);
}

/*}}}*/

int Plot_next_page (Plot_Node_t *n) /*{{{*/
{
   Plot_Pane_Type *w;
   int num_panes;

   if (n->window_type == PLOT_STD_WINDOW)
     {
        Plot_Pane_Type *vp = &n->outer_viewport;
        _Plot_next_page ();
        return Plot_select_viewport (vp->vp_xmin, vp->vp_xmax,
                                     vp->vp_ymin, vp->vp_ymax);
     }

   num_panes = n->nxpanes * n->nypanes;
   if ((n->this_pane == num_panes) || (n->force_clear != 0))
     {
        _Plot_next_page ();
        n->this_pane = 0;
        n->force_clear = 0;
     }

   w = &n->panes[n->this_pane];
   n->this_pane++;

   return Plot_select_viewport (w->vp_xmin, w->vp_xmax,
                                w->vp_ymin, w->vp_ymax);
}

/*}}}*/

static int get_plot_limits (Plot_t *fmt, float *xmin, float *xmax, float *ymin, float *ymax) /*{{{*/
{
   *xmin = fmt->xmin ? *fmt->xmin : *xmin;
   *xmax = fmt->xmax ? *fmt->xmax : *xmax;
   *ymin = fmt->ymin ? *fmt->ymin : *ymin;
   *ymax = fmt->ymax ? *fmt->ymax : *ymax;

   return 0;
}

/*}}}*/

void Plot_draw_box_opt (Plot_t *fmt, float xmin, float xmax, float ymin, float ymax) /*{{{*/
{
   _Plot_set_plot_limits (xmin, xmax, ymin, ymax);
   _Plot_set_color (PLOT_LINE_COLOR);
   _Plot_set_linestyle (PLOT_SOLID_LINE);
   _Plot_set_line_width (fmt->frame_line_width);
   _Plot_draw_box (fmt->xopt, fmt->xtick, fmt->nxsub,
                   fmt->yopt, fmt->ytick, fmt->nysub);
   _Plot_set_line_width (fmt->line_width);
   _Plot_set_color (fmt->color);
   fmt->rendered = 1;
}

/*}}}*/

static void Plot_draw_box (Plot_t *fmt, float xmin, float xmax, float ymin, float ymax) /*{{{*/
{
   update_axis_type (fmt);
   (void) get_plot_limits (fmt, &xmin, &xmax, &ymin, &ymax);
   _Plot_set_color (PLOT_LINE_COLOR);
   _Plot_set_linestyle (PLOT_SOLID_LINE);
   Plot_draw_box_opt (fmt, xmin, xmax, ymin, ymax);
   _Plot_set_color (fmt->color);
   _Plot_set_linestyle (fmt->linestyle);
}

/*}}}*/

void Plot_label_box (Plot_t *fmt) /*{{{*/
{
   Plot_Node_t *w = fmt->this_window;

   _Plot_set_color (PLOT_LINE_COLOR);
   _Plot_set_line_width (fmt->frame_line_width);

   if (w->window_type == PLOT_STD_WINDOW)
     _Plot_label_box (fmt->xlabel, fmt->ylabel, fmt->tlabel);
   else
     {  /* magic numbers chosen to agree with pglab defaults */
        Plot_put_text_offset ("L", 2.2, 0.5, 0.5, fmt->ylabel);
        if (w->this_pane == 1)
          Plot_put_text_offset ("T", 2.0, 0.5, 0.5, fmt->tlabel);
        if (w->this_pane == w->nypanes)
          Plot_put_text_offset ("B", 3.2, 0.5, 0.5, fmt->xlabel);
     }

   _Plot_set_color (fmt->color);
   _Plot_set_line_width (fmt->line_width);
}

/*}}}*/

void Plot_set_color (Plot_t *fmt) /*{{{*/
{
   _Plot_set_color (fmt->color);
}

/*}}}*/

void Plot_set_charsize (Plot_t *fmt) /*{{{*/
{
   _Plot_set_charsize (fmt->char_height);
}

/*}}}*/

void Plot_set_linestyle (Plot_t *fmt) /*{{{*/
{
   _Plot_set_linestyle (fmt->linestyle);
}

/*}}}*/

void Plot_set_line_width (Plot_t *fmt) /*{{{*/
{
   _Plot_set_line_width (fmt->line_width);
}

/*}}}*/

void Plot_set_frame_line_width (Plot_t *fmt) /*{{{*/
{
   _Plot_set_line_width (fmt->frame_line_width);
}

/*}}}*/

/*}}}*/

static void validate_overlay_mode (Plot_t *fmt, int *overlay) /*{{{*/
{
   /* let oplot behave like plot when nothing has been plotted yet */
   if ((fmt->rendered == 0) && (*overlay != 0))
     {
        *overlay = 0;
     }
}

/*}}}*/

int Plot_box (Plot_t *fmt) /*{{{*/
{
   if ((fmt->xmin == NULL) || (fmt->xmax == NULL)
       || (fmt->ymin == NULL) || (fmt->ymax == NULL))
     return -1;

   if (-1 == Plot_next_page (fmt->this_window))
     return -1;
   Plot_set_charsize (fmt);
   Plot_draw_box (fmt, *fmt->xmin, *fmt->xmax, *fmt->ymin, *fmt->ymax);
   Plot_label_box (fmt);
   Plot_update_display ();

   return 0;
}

/*}}}*/

static int Plot_sized_points (Plot_t *fmt, int n, float *x, float *y) /*{{{*/
{
   _Plot_set_charsize (Plot_get_point_size (fmt));
   Plot_points (n, x, y, Plot_get_pointstyle (fmt));
   Plot_set_charsize (fmt);
   return 0;
}

/*}}}*/

int Plot_curve (float *x, float *y, int npts, int *symbol, /*{{{*/
                int *use_style, int overlay, Plot_t *fmt)
{
   float *px, *py;
   float fmt_xmin, fmt_xmax, fmt_ymin, fmt_ymax;
   int ret = -1;

   px = py = NULL;

   if (fmt == NULL)
     return -1;

   validate_overlay_mode(fmt, &overlay);

   if (!overlay)
     {
        if (-1 == Plot_next_page (fmt->this_window))
          return -1;
     }

   if (!overlay)
     _restart_style_cycle (fmt);

   if (NULL == (px = (float *) ISIS_MALLOC (2 * npts * sizeof(float))))
     return -1;
   memset ((char *)px, 0, 2 * npts * sizeof(float));
   py = px + npts;

   /* compute logs if necessary */
   if (fmt->logx == 0)
     memcpy ((char *)px, (char *)x, npts * sizeof(float));
   else if (-1 == compute_log10 (x, npts, px))
     goto return_error;

   if (fmt->logy == 0)
     memcpy ((char *)py, (char *)y, npts * sizeof(float));
   else if (-1 == compute_log10 (y, npts, py))
     goto return_error;

   /* X-axis range */
   if ((-1 == find_plot_range (px, npts, fmt->logx, &fmt_xmin, &fmt_xmax))
       || (-1 == validate_range (&fmt_xmin, &fmt_xmax)))
     goto return_error;

   if (0 == have_data_in_xrange (fmt, fmt_xmin, fmt_xmax))
     isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "No data in specified X limits");

   /* Y-axis range */
   if ((-1 == yrange_for_xrange (fmt, py, px, npts, &fmt_ymin, &fmt_ymax))
       || (-1 == validate_range (&fmt_ymin, &fmt_ymax)))
     goto return_error;

   if (0 == have_data_in_yrange (fmt, fmt_ymin, fmt_ymax))
     isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "No data in specified Y limits");

   /* leave some padding between the curve and the plot axes */
   if (-1 == pad_range (&fmt_xmin, &fmt_xmax)
       || -1 == pad_range (&fmt_ymin, &fmt_ymax))
     goto return_error;

   Plot_set_charsize (fmt);

   if (!overlay)
     {
        Plot_draw_box (fmt, fmt_xmin, fmt_xmax, fmt_ymin, fmt_ymax);
        Plot_label_box (fmt);
     }

   Plot_set_linestyle (fmt);
   Plot_set_style (fmt, use_style);
   _Plot_set_line_width (fmt->line_width);

   if (fmt->connect_points)
     Plot_line (npts, px, py);

   if (fmt->connect_points >= 0)
     {
        if (symbol)
          Plot_symbol_points (npts, px, py, symbol);
        else
          Plot_sized_points (fmt, npts, px, py);
     }

   ret = 0;
   return_error:

   fmt->plotted_x_unit = -1;
   Plot_update_display ();
   ISIS_FREE (px);

   return ret;
}

/*}}}*/

static float double_to_float (double s) /*{{{*/
{
   if (fabs(s) < FLT_MAX)
     return (float) s;
   return (0.99 * FLT_MAX) * ((s < 0) ? -1.0 : 1.0);
}

/*}}}*/

static int plot_y_errorbars (Plot_t *fmt, int nbins, /*{{{*/
                             float *pxmid, float *pytop, float *pybot)
{
   float *x, *xmid, *ytop, *ybot;
   int k;

   if (fmt->use_errorbars == 0)
     return 0;

   x = NULL;

   if (fmt->use_errorbars == 1)
     {
        /* plot _all_ errorbars */
        k = nbins;
        xmid = pxmid;
        ytop = pytop;
        ybot = pybot;
     }
   else
     {
        /* plot errorbars separated by 'step' bins */
        int i, n, step;

        step = abs(fmt->use_errorbars);
        n = nbins / step;
        if (NULL == (x = (float *) ISIS_MALLOC (3 * n * sizeof(float))))
          return -1;

        xmid = x;
        ytop = x + n;
        ybot = x + 2*n;

        k = 0;
        for (i = 0; i < nbins && k < n; i += step)
          {
             xmid[k] = pxmid[i];
             ytop[k] = pytop[i];
             ybot[k] = pybot[i];
             k++;
          }
     }

   Plot_y_errorbar (k, xmid, ytop, ybot, fmt->ebar_term_length);

   ISIS_FREE (x);

   return 0;
}

/*}}}*/

int Plot_hist (Plot_t *fmt, int *use_style, int overlay, Isis_Hist_t *h) /*{{{*/
{
   float *lo, *hi, *val, *val_err;
   float *xmid, *tmp, *ytop, *ybot;
   float fmt_xmin, fmt_xmax, fmt_ymin, fmt_ymax;
   int i, nbins;
   int ret = -1;

   lo = hi = val = val_err = NULL;
   xmid = tmp = ytop = ybot = NULL;

   if ((h == NULL) || (NULL == fmt))
     return -1;

   validate_overlay_mode(fmt, &overlay);

   if (!overlay)
     {
        if (-1 == Plot_next_page (fmt->this_window))
          return -1;
     }

   nbins = h->nbins;

   if (NULL == (lo = (float *) ISIS_MALLOC (nbins * sizeof(float)))
       || NULL == (hi = (float *) ISIS_MALLOC (nbins * sizeof(float)))
       || NULL == (val = (float *) ISIS_MALLOC (nbins * sizeof(float)))
       || NULL == (val_err = (float *) ISIS_MALLOC (nbins * sizeof(float)))
       || NULL == (ytop = (float *) ISIS_MALLOC (nbins * sizeof(float)))
       || NULL == (ybot = (float *) ISIS_MALLOC (nbins * sizeof(float)))
       || NULL == (xmid = (float *) ISIS_MALLOC (nbins * sizeof(float)))
       || NULL == (tmp = (float *) ISIS_MALLOC ((nbins+1) * sizeof(float))))
     goto return_error;

   /*
    * Default plot is bin-integral values.
    * If the user wants to see bin-density values,
    * they get computed assuming the physical units
    * Y and X are consistent -- e.g. that a simple
    * division by the bin-width is ok and no physical
    * unit conversion is needed.
    */

   for (i = 0; i < nbins; i++)
     {
        val_err[i] = 0.0;
        lo[i] = double_to_float (h->bin_lo[i]);
        hi[i] = double_to_float (h->bin_hi[i]);
        if (hi[i]-lo[i] > 0.0)
          continue;
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "bin %d has single precision width = %g",
                    i, hi[i]-lo[i]);
        goto return_error;
     }

   if (0 == Plot_bin_density(fmt))
     {
        for (i = 0; i < nbins; i++) {
           val[i] = double_to_float (h->val[i]);
        }
        if (h->val_err) {
           for (i=0; i < nbins; i++) {
              val_err[i] = double_to_float (h->val_err[i]);
           }
        }
     }
   else
     {
        double dx;
        /* compute bin-density */
        for (i=0; i < nbins; i++) {
           dx = h->bin_hi[i] - h->bin_lo[i];
           val[i] = double_to_float (h->val[i] / dx);
        }
        if (h->val_err) {
           for (i=0; i < nbins; i++) {
              dx = h->bin_hi[i] - h->bin_lo[i];
              val_err[i] = double_to_float (h->val_err[i] / dx);
           }
        }
     }

   if (NULL != h->notice)
     {
        /* zero out ignored bins */
        for (i=0; i < nbins; i++)
          {
             if (h->notice[i] == 0)
               {
                  val[i] = 0.0;
                  val_err[i] = 0.0;
               }
          }
     }

   if (!overlay)
     _restart_style_cycle (fmt);

   if (-1 == convert_plot_units (fmt))
     goto return_error;

   for (i = 0; i < nbins; i++)
     {
        xmid[i] = 0.5 * (lo[i] + hi[i]);
     }

   if (-1 == compute_errorbars (val, val_err, nbins, ybot, ytop))
     goto return_error;

   /* compute logs if necessary */
   if (fmt->logx)
     {
        if ((-1 == compute_log10 (lo, nbins, lo))
            || (-1 == compute_log10 (hi, nbins, hi))
            || (-1 == compute_log10 (xmid, nbins, xmid)))
          goto return_error;
     }

   if (fmt->logy)
     {
        if ((-1 == compute_log10 (val, nbins, val))
            || (-1 == compute_log10 (ybot, nbins, ybot))
            || (-1 == compute_log10 (ytop, nbins, ytop)))
          goto return_error;
     }

   /* X-axis range */
   memcpy ((char *)tmp, (char *)lo, nbins * sizeof(float));
   tmp[nbins] = hi[nbins-1];

   if ((-1 == find_plot_range (tmp, nbins+1, fmt->logx, &fmt_xmin, &fmt_xmax))
       || (-1 == validate_range (&fmt_xmin, &fmt_xmax)))
     goto return_error;

   if (0 == have_data_in_xrange (fmt, fmt_xmin, fmt_xmax))
     isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "No data in specified X limits");

   /* Y-axis range */
   if ((-1 == yrange_for_xrange (fmt, val, xmid, nbins, &fmt_ymin, &fmt_ymax))
       || (-1 == validate_range (&fmt_ymin, &fmt_ymax)))
     goto return_error;

   if (0 == have_data_in_yrange (fmt, fmt_ymin, fmt_ymax))
     isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "No data in specified Y limits");

   /* leave some padding between the curve and the plot axes */
   if (-1 == pad_range (&fmt_xmin, &fmt_xmax)
       || -1 == pad_range (&fmt_ymin, &fmt_ymax))
     goto return_error;

   Plot_set_charsize (fmt);

   if (!overlay)
     {
        Plot_draw_box (fmt, fmt_xmin, fmt_xmax, fmt_ymin, fmt_ymax);
        Plot_label_box (fmt);
     }

   Plot_set_linestyle (fmt);
   Plot_set_style (fmt, use_style);
   _Plot_set_line_width (fmt->line_width);

   if (fmt->connect_points == 0)
     {
        Plot_sized_points (fmt, nbins, xmid, val);
        plot_y_errorbars (fmt, nbins, xmid, ytop, ybot);
     }
   else if (h->notice == NULL)
     {
        Plot_histogram_data (nbins, lo, hi, val);
        plot_y_errorbars (fmt, nbins, xmid, ytop, ybot);
     }
   else
     {
        /* don't plot ignored bins */
        int k = 0;
        while (k < nbins)
          {
             int n_nbins;
             while ((k < nbins) && (h->notice[k] == 0))
               k++;
             if (k == nbins)
               break;
             n_nbins = 1;
             while ((k + n_nbins < nbins) && (h->notice[k + n_nbins] != 0))
               n_nbins++;
             if (k + n_nbins > nbins)
               n_nbins = nbins - k;
             Plot_histogram_data (n_nbins, lo+k, hi+k, val+k);
             plot_y_errorbars (fmt, n_nbins, xmid+k, ytop+k, ybot+k);
             k += n_nbins;
          }
     }

   ret = 0;

   return_error:

   fmt->plotted_x_unit = fmt->x_unit;
   fmt->ok_to_convert_xrange_units = 1;

   Plot_update_display ();

   ISIS_FREE (lo);   ISIS_FREE (hi);
   ISIS_FREE (val);  ISIS_FREE (val_err);
   ISIS_FREE (xmid); ISIS_FREE (tmp);
   ISIS_FREE (ytop); ISIS_FREE (ybot);

   return ret;
}

/*}}}*/

/* read box cursor */

int Plot_read_box_corners_from_cursor (float bx[/*2*/], float by[/*2*/]) /*{{{*/
{
   float x[2],y[2];
   char ch;

   fprintf (stdout, "click box corners  OR  r)eselect  q)uit\n");
   fflush (stdout);

   do
     {
        x[0] = x[1] = y[0] = y[1] = 0.0;

        if (-1 == Plot_read_crosshair_cursor (&x[0], &y[0], &ch))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading the cursor");
             return -1;
          }

        if (ch == 'q' || ch == 'Q')
          return -1;

        if (-1 == Plot_read_rectangle_cursor (x[0], y[0], &x[1], &y[1], &ch))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading the cursor");
             return -1;
          }

        if (ch == 'q' || ch == 'Q')
          return -1;

     } while (ch == 'r' || ch == 'R');

   bx[0] = MIN (x[0], x[1]);     by[0] = MIN (y[0], y[1]);
   bx[1] = MAX (x[0], x[1]);     by[1] = MAX (y[0], y[1]);

   return 0;
}

/*}}}*/

int Plot_read_line_endpoints (float bx[/*2*/], float by[/*2*/], int draw_line) /*{{{*/
{
   float x[2],y[2];
   char ch;

   fprintf (stdout, "click line endpoints  OR   r)eselect    q)uit\n");
   fflush (stdout);

   do
     {
        x[0] = x[1] = y[0] = y[1] = 0.0;

        if (-1 == Plot_read_crosshair_cursor (&x[0], &y[0], &ch))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading the cursor");
             return -1;
          }

        if (ch == 'q' || ch == 'Q')
          return -1;

        if (-1 == Plot_cursor_read_other_endpoint (x[0], y[0], &x[1], &y[1], &ch))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading the cursor");
             return -1;
          }

        if (ch == 'q' || ch == 'Q')
          return -1;

     } while (ch == 'r' || ch == 'R');

   if (draw_line)
     Plot_line (2, x, y);

   bx[0] = MIN (x[0], x[1]);     by[0] = MIN (y[0], y[1]);
   bx[1] = MAX (x[0], x[1]);     by[1] = MAX (y[0], y[1]);

   return 0;
}

/*}}}*/

int Plot_LS_label (float x, char *text) /*{{{*/
{
   int save;

   save = _Plot_get_color ();
   _Plot_set_color (PLOT_LINE_COLOR);
   Plot_put_text_offset ("B",
                         1.5,     /* dist. from axis in char-heights */
                         x,
                         0.5,     /* 0.5 = center justified */
                         text);
   _Plot_set_color (save);
   return 0;
}

/*}}}*/
