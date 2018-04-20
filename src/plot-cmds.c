/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2018  Massachusetts Institute of Technology

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

/* $Id: plot-cmds.c,v 1.74 2004/09/09 17:06:30 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <slang.h>
#include <float.h>

#include "isis.h"
#include "plot.h"
#include "util.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

/*{{{ globals */

/* ptr to format list */
static Plot_Node_t *Plot_List = NULL;
static Plot_Node_t *Window_sans_device = NULL;
static Plot_t *Saved_Format = NULL;
static unsigned int Num_Open_Plots = 0;
static char Default_Plot_Device [BUFSIZE] = "?";

/* current plot focus */
static int Focus = PLOT_NO_DEVICE;

static void _set_window (int *device);
/*}}}*/

/* Keep a history of the most recent N plot windows
 * so we know where to inherit parameters for next
 * plot window opened
 */

enum
{
   HISTORY_SIZE = 128
};

typedef struct
{
   int where[HISTORY_SIZE];
   unsigned int next;
   unsigned int size;
}
History_Type;

static History_Type History;

static void push_device (int dev) /*{{{*/
{
   History_Type *h = &History;

   Focus = dev;

   if (dev == PLOT_NO_DEVICE)
     return;

   h->where[h->next] = dev;
   h->next++;
   h->next = h->next % HISTORY_SIZE;

   if (h->size < HISTORY_SIZE)
     h->size++;
}

/*}}}*/

static int pop_device (void) /*{{{*/
{
   History_Type *h = &History;

   while (h->size > 0)
     {
        int dev;

        if (h->next == 0)
          h->next += HISTORY_SIZE;

        h->next--;
        dev = h->where[h->next];
        h->size--;

        if (NULL != Plot_find_node (dev, Plot_List))
          return dev;
     }

   return PLOT_NO_DEVICE;
}

/*}}}*/

static int peek_device (void) /*{{{*/
{
   int dev = pop_device ();

   if (dev != PLOT_NO_DEVICE)
     push_device (dev);

   return dev;
}

/*}}}*/

static void update_saved_format (void) /*{{{*/
{
   Plot_t *fmt = current_format ();

   if (fmt)
     {
        Plot_free_fmt (Saved_Format);
        Saved_Format = Plot_alloc_fmt ();
        Plot_copy_formats (fmt, Saved_Format);  /* NULL ok */
     }
}

/*}}}*/

static void new_focus (int dev) /*{{{*/
{
   push_device (dev);
   Focus = dev;
}

/*}}}*/

static int init_plot (void) /*{{{*/
{
   if ((NULL == Plot_List)
       && (NULL == (Plot_List = Plot_start ())))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing plot window");
        return -1;
     }

   return 0;
}

/*}}}*/

static void _quit_plot (void) /*{{{*/
{
   Plot_free_fmt (Saved_Format);
   Saved_Format = NULL;

   Plot_end (&Plot_List);
   Focus = PLOT_NO_DEVICE;
   Num_Open_Plots = 0;

   /* Plot_free_library_interface (); */

   memset ((char *)&History, 0, sizeof(History_Type));
}

/*}}}*/

/* open/close/erase/resize */

static int inherit_plot_settings (Plot_Node_t *to_node) /*{{{*/
{
   Plot_Node_t *from_node;
   Plot_t *from_fmt;
   int from_node_id = peek_device();

   if (from_node_id == PLOT_NO_DEVICE)
     {
        /* keep defaults if nothing to inherit from */
        if (Saved_Format == NULL)
          return 0;
        from_fmt = Saved_Format;
     }
   else
     {
        if (NULL == (from_node = Plot_find_node (from_node_id, Plot_List)))
          return -1;
        from_fmt = Plot_fmt (from_node);
     }

   return Plot_copy_formats (from_fmt, Plot_fmt (to_node));
}

/*}}}*/

static int _open_plot (char *device, int *nxpanes, int *nypanes) /*{{{*/
{
   Plot_Node_t *node;
   int dev;

   if (-1 == init_plot ())
     return -1;

   if (Window_sans_device != NULL)
     {
        node = Window_sans_device;
        Window_sans_device = NULL;
     }
   else
     {
        node = Plot_new_window (Plot_List);
        if ((NULL == node)
            || (-1 == Plot_reformat_node (node, *nxpanes, *nypanes, 0, NULL)))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing plot window");
             return -1;
          }
     }

   if (*device == '\0')
     device = Default_Plot_Device;

   if (-1 == (dev = Plot_open_window (node, device)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing plot window");
        return -1;
     }
   Plot_subdivide (*nxpanes, *nypanes);

   Num_Open_Plots++;

   if (-1 == inherit_plot_settings (node))
     return -1;
   new_focus (dev);

   return dev;
}

/*}}}*/

typedef struct Viewport_Type Viewport_Type;
struct Viewport_Type
{
   float xmin;
   float xmax;
   float ymin;
   float ymax;
};

#define F SLANG_FLOAT_TYPE

static SLang_CStruct_Field_Type Viewport_Table [] =
{
   MAKE_CSTRUCT_FIELD(Viewport_Type, xmin, "xmin", F, 0),
   MAKE_CSTRUCT_FIELD(Viewport_Type, xmax, "xmax", F, 0),
   MAKE_CSTRUCT_FIELD(Viewport_Type, ymin, "ymin", F, 0),
   MAKE_CSTRUCT_FIELD(Viewport_Type, ymax, "ymax", F, 0),
   SLANG_END_CSTRUCT_TABLE
};

#undef F

static void get_outer_viewport_intrin (void) /*{{{*/
{
   Plot_Node_t *n;
   Viewport_Type v;

   n = Plot_find_node (Focus, Plot_List);

   Plot_get_outer_viewport (n, &v.xmin, &v.xmax, &v.ymin, &v.ymax);
   if (-1 == SLang_push_cstruct ((VOID_STAR)&v, Viewport_Table))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't push cstruct");
     }
}

/*}}}*/

static void set_outer_viewport_intrin (void) /*{{{*/
{
   Viewport_Type v;
   Plot_Node_t *n;
   Plot_t *fmt;

   if ((NULL == (fmt = current_format ())
        || (force_open_plot_device() < 0)))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        return;
     }

   n = Plot_find_node (Focus, Plot_List);
   if (n == NULL)
     return;

   memset ((char *)&v, 0, sizeof (v));
   if (-1 == SLang_pop_cstruct ((VOID_STAR)&v, Viewport_Table))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't pop cstruct");
        return;
     }

   if (-1 == Plot_set_outer_viewport (n, v.xmin, v.xmax, v.ymin, v.ymax))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "failed setting outer viewport");
        return;
     }

   SLang_free_cstruct ((VOID_STAR)&v, Viewport_Table);
}

/*}}}*/

static void multiplot (void) /*{{{*/
{
   SLang_Array_Type *sl_ysizes = NULL;
   int *ysizes, num_ysizes;
   Plot_Node_t *n;

   if ((-1 == SLang_pop_array_of_type (&sl_ysizes, SLANG_INT_TYPE))
       || (sl_ysizes == NULL))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "bad multiplot specifier");
        return;
     }

   if ((NULL == current_format ())
       || (force_open_plot_device () < 0))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        SLang_free_array (sl_ysizes);
        return;
     }

   if ((NULL == Plot_List)
       || (NULL == (n = Plot_find_node (Focus, Plot_List))))
     {
        isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "plot device %d", Focus);
        SLang_free_array (sl_ysizes);
        return;
     }

   num_ysizes = sl_ysizes->num_elements;
   ysizes = (int *) sl_ysizes->data;

   if (-1 == Plot_reformat_node (n, 1, num_ysizes, (num_ysizes > 1), ysizes))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing plot window");
        SLang_free_array (sl_ysizes);
        return;
     }

   if (num_ysizes == 1)
     {
        (void) Plot_visit_pane (n, 1);
     }

   SLang_free_array (sl_ysizes);
}

/*}}}*/

int force_open_plot_device (void) /*{{{*/
{
   int dev;

   if (Focus != PLOT_NO_DEVICE)
     return 0;

   dev = Plot_open_window (Window_sans_device, Default_Plot_Device);
   if (dev > 0)
     new_focus (dev);
   Window_sans_device = NULL;

   return dev;
}

/*}}}*/

Plot_t *current_format (void) /*{{{*/
{
   Plot_Node_t *n;
   Plot_t *p;

   n = Plot_find_node (Focus, Plot_List);
   p = Plot_fmt (n);

   if (p == NULL)
     {
        if (-1 == init_plot ())
          return NULL;

        Window_sans_device = Plot_new_window (Plot_List);
        if ((NULL == Window_sans_device)
            || (-1 == Plot_reformat_node (Window_sans_device, 1, 1, 0, NULL)))
          {
             isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing plot window");
             return NULL;
          }

        Num_Open_Plots++;

        if (-1 == inherit_plot_settings (Window_sans_device))
          return NULL;
        new_focus (PLOT_NO_DEVICE);

        p = Plot_fmt (Window_sans_device);
     }

   return p;
}

/*}}}*/

static int _copy_plot (char *new_device, int *id) /*{{{*/
{
   int pid = *id;
   int dev;

   if (NULL == Plot_List)
     {
        isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "no plots to copy");
        return -1;
     }

   if (*new_device == '\0')
     new_device = Default_Plot_Device;

   if (pid <= 0)
     pid = Focus;

   if (-1 == (dev = Plot_copy_window (new_device, pid, Plot_List)))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "initializing plot window");
   else
     {
        Num_Open_Plots++;
        new_focus (dev);
     }

   return dev;
}

/*}}}*/

static void _set_plot_device (char *device_string) /*{{{*/
{
   isis_strcpy (Default_Plot_Device, device_string, BUFSIZE);
}

/*}}}*/

static void _set_auto_increment (int *flag) /*{{{*/
{
   Plot_set_auto_increment (*flag);
}

/*}}}*/

void reset_current_pane (Plot_t *fmt) /*{{{*/
{
   (void) fmt;
}
/*}}}*/

static void _close_plot (void) /*{{{*/
{
   int device = Focus;

   if (SLang_Num_Function_Args > 1
       || (SLang_Num_Function_Args == 1
           && -1 == SLang_pop_integer (&device)))
     {
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "invalid input (plot_close)");
        return;
     }

   if (device == Focus)
     update_saved_format();

   if (-1 == Plot_delete_window (device, Plot_List))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "closing plot window");
   else if (Num_Open_Plots)
     Num_Open_Plots--;

   if (device == Focus)
     Focus = pop_device ();

   _set_window (&Focus);
}

/*}}}*/

static void _erase_plot (void) /*{{{*/
{
   int device = Focus;

   if (Focus == PLOT_NO_DEVICE)
     return;

   if (SLang_Num_Function_Args > 1
       || (SLang_Num_Function_Args == 1
           && -1 == SLang_pop_integer (&device)))
     {
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "invalid input (plot_erase)");
        return;
     }

   if (NULL == current_format ())
     return;

   (void) Plot_erase ();
}

/*}}}*/

static void _resize_window (float *width, float *aspect) /*{{{*/
{
   Plot_select_window (Focus);

   if (*width < 0.0 || *aspect <= 0.0)
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "resizing plot window");
   else
     Plot_reset_viewer_size (*width, *aspect);
}

/*}}}*/

/* get/set various parameters  */

static void _set_window (int *device) /*{{{*/
{
   if (*device == PLOT_NO_DEVICE)
     return;

   if (NULL == Plot_find_node (*device, Plot_List))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "invalid plot device: id = %d", *device);
        return;
     }

   new_focus (*device);
   if (Focus != PLOT_NO_DEVICE)
     Plot_select_window (Focus);
}

/*}}}*/

static void _set_pane (int *pane) /*{{{*/
{
   Plot_Node_t *n;

   n = Plot_find_node (Focus, Plot_List);
   if (n == NULL)
     return;

   /* for multiplots only! */
   (void) Plot_visit_pane (n, *pane);
}

/*}}}*/

static void _set_labels (char * xlabel, char * ylabel, char * tlabel) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     return;

    if (-1 == Plot_set_labels (fmt, xlabel, ylabel, tlabel))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting plot labels");
        return;
     }

   Plot_set_label_type (fmt, PLOT_USER_LABEL);
}

/*}}}*/

static void _set_axis_string (char *s, int (*setter)(Plot_t *, const char *)) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     return;

    if (-1 == (*setter)(fmt, s))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting label string");
        return;
     }
}

/*}}}*/

static void _set_title (char * s) /*{{{*/
{
   _set_axis_string (s, &Plot_set_title);
}

/*}}}*/

static void _set_xlabel (char * s) /*{{{*/
{
   _set_axis_string (s, &Plot_set_xlabel);
}

/*}}}*/

static void _set_ylabel (char * s) /*{{{*/
{
   _set_axis_string (s, &Plot_set_ylabel);
}

/*}}}*/

static void _set_charsize (float * char_height) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ())
       || -1 == Plot_set_fmt_charsize (fmt, *char_height))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting character size");
}

/*}}}*/

static void _set_style_meaning (int * style) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting style");
   else
     Plot_set_style_meaning (fmt, *style);
}

/*}}}*/

static void set_pointstyle (int *pointstyle) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ())
       || -1 == Plot_set_fmt_pointstyle (fmt, *pointstyle))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting point style");
}

/*}}}*/

static void set_point_size (float *point_size) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ())
       || -1 == Plot_set_fmt_point_size (fmt, *point_size))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting point size");
}

/*}}}*/

static int _get_line_width (void) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "getting line width");
        return 0;
     }

   return Plot_get_fmt_line_width (fmt);
}
/*}}}*/

static void _set_line_width (int * width) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting line width");
   else
     (void) Plot_set_fmt_line_width (fmt, *width);
}
/*}}}*/

static void _set_frame_line_width (int * width) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting frame line width");
   else
     (void) Plot_set_fmt_frame_line_width (fmt, *width);
}
/*}}}*/

static void _set_color (int * color) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting line color");
   else
     (void) Plot_set_fmt_color (fmt, *color);
}

/*}}}*/

static void _get_color (void) /*{{{*/
{
   Plot_t *fmt;
   int color = PLOT_LINE_COLOR;

   fmt = current_format ();
   if (fmt)
     color = Plot_get_fmt_color (fmt);

   (void) SLang_push_integer (color);
}

/*}}}*/

static void _axis_is_log (char axis) /*{{{*/
{
   Plot_t *fmt = current_format ();
   int islog;

   /* stack underflow on error */
   if (fmt == NULL)
     return;

   switch (axis)
     {
      case 'x':
      case 'X':
        islog = Plot_get_logx (fmt);
        break;
      case 'y':
      case 'Y':
        islog = Plot_get_logy (fmt);
        break;
      default:
        /* stack underflow on error */
        return;
     }

   SLang_push_integer (islog);
}

/*}}}*/

static void _xaxis_is_log (void) /*{{{*/
{
   _axis_is_log ('x');
}

/*}}}*/

static void _yaxis_is_log (void) /*{{{*/
{
   _axis_is_log ('y');
}

/*}}}*/

static void set_linestyle (int * linestyle) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting line style");
   else
     (void) Plot_set_fmt_linestyle (fmt, *linestyle);
}

/*}}}*/

typedef int (*Set_Range_Fcn) (Plot_t *, float *, float *);

static void _set_range (Set_Range_Fcn set_range) /*{{{*/
{
   float min, max;
   float *pmin = NULL;
   float *pmax = NULL;
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting axis range");
        return;
     }

   switch (SLang_Num_Function_Args)
     {
      case 2:
        if (SLANG_NULL_TYPE == SLang_peek_at_stack())
          SLdo_pop();
        else if ((0 == SLang_pop_float (&max)) && (max < FLT_MAX))
          pmax = &max;
        /* drop */
      case 1:
        if (SLANG_NULL_TYPE == SLang_peek_at_stack())
          SLdo_pop();
        else if ((0 == SLang_pop_float (&min)) && (min > -FLT_MAX))
          pmin = &min;
        break;
      default:
        break;
     }

   if (-1 == (*set_range) (fmt, pmin, pmax))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "changing axis limits");
}

/*}}}*/

static void _set_xrange (void) /*{{{*/
{
   _set_range (Plot_set_xrange);
}

/*}}}*/

static void _set_yrange (void) /*{{{*/
{
   _set_range (Plot_set_yrange);
}

/*}}}*/

static void _get_axis_range (int (*get_range)(Plot_t *, float **, float **)) /*{{{*/
{
   float *pmin, *pmax;
   Plot_t *fmt;

   pmin = pmax = NULL;

   if (NULL == (fmt = current_format ()))
     {
        SLang_push_double(-DBL_MAX);
        SLang_push_double( DBL_MAX);
        return;
     }

   (*get_range)(fmt, &pmin, &pmax);

   if (pmin != NULL)
     SLang_push_float(*pmin);
   else
     SLang_push_double(-DBL_MAX);

   if (pmax != NULL)
     SLang_push_float(*pmax);
   else
     SLang_push_double(DBL_MAX);
}

/*}}}*/

static void _get_yrange (void) /*{{{*/
{
   _get_axis_range (&Plot_get_yrange);
}

/*}}}*/

static void _get_xrange (void) /*{{{*/
{
   _get_axis_range (&Plot_get_xrange);
}

/*}}}*/

static void _reset_plot_ranges (void) /*{{{*/
{
   Plot_t *fmt;

   if (NULL != (fmt = current_format ()))
     Plot_reset_xyranges (fmt);

   return;
}

/*}}}*/

static void set_fmt_boolean (void (*fcn)(Plot_t *, int), int value) /*{{{*/
{
   Plot_t *fmt;

   if (NULL != (fmt = current_format ()))
     (*fcn) (fmt, value);
}

/*}}}*/

static void _logx_on (void) /*{{{*/
{
   set_fmt_boolean (Plot_set_logx, 1);
}

/*}}}*/

static void _logx_off (void) /*{{{*/
{
   set_fmt_boolean (Plot_set_logx, 0);
}

/*}}}*/

static void _logy_on (void) /*{{{*/
{
   set_fmt_boolean (Plot_set_logy, 1);
}

/*}}}*/

static void _logy_off (void) /*{{{*/
{
   set_fmt_boolean (Plot_set_logy, 0);
}

/*}}}*/

static void bin_density_on (void) /*{{{*/
{
   set_fmt_boolean (Plot_set_bin_density, 1);
}

/*}}}*/

static void bin_density_off (void) /*{{{*/
{
   set_fmt_boolean (Plot_set_bin_density, 0);
}

/*}}}*/

int current_plot_shows_bin_density (void) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     return -1;

   return Plot_bin_density (fmt);
}

/*}}}*/

static void _set_errorbar (int *value, float *term_length) /*{{{*/
{
   Plot_t *fmt;

   if (NULL != (fmt = current_format ()))
     Plot_set_errorbars (fmt, *value, *term_length);
}

/*}}}*/

static void _set_connect_points (int *value) /*{{{*/
{
   set_fmt_boolean (Plot_set_connect_points, *value);
}

/*}}}*/

static void _connect_points (void) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     (void) SLang_push_integer (-1);
   else
     (void) SLang_push_integer (Plot_connect_points (fmt));
}

/*}}}*/

static void _errorbar_state (void) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     (void) SLang_push_integer (-1);
   else
     (void) SLang_push_integer (Plot_use_errorbars (fmt));
}

/*}}}*/

static void _set_plot_units (char *name) /*{{{*/
{
   Plot_t *fmt;
   int unit;

   if (NULL == (fmt = current_format ()))
     return;

   if (-1 == (unit = unit_id (name)))
     return;

   if (-1 == Plot_set_x_unit (fmt, unit))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "setting X axis units");
}

/*}}}*/

int change_xunits_to_angstrom (float *x, float *y) /*{{{*/
{
   Plot_t *fmt;
   int x_unit;
   double tx, ty;

   if (NULL == (fmt = current_format ()))
     return -1;

   if (Plot_get_logx (fmt))          /* handle logarithmic axes */
     *x = pow (10.0, *x);

   if (Plot_get_logy (fmt))
     *y = pow (10.0, *y);

   *y *= Plot_get_exposure_time (fmt);

   x_unit = Plot_x_unit (fmt);
   if (x_unit == U_ANGSTROM)
     return 0;

   /* handle non-angstrom units */

   tx = (double) *x;
   ty = (double) *y;            /* y is an x-density, e.g. df/dx */

   if (-1 == unit_convert_dfdx (&ty, tx, U_ANGSTROM, x_unit)
       ||-1 == unit_convert_x (&tx, 1, U_ANGSTROM, x_unit))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "converting X-axis units");
        return -1;
     }

   *x = (float) tx;
   *y = (float) ty;

   return 0;
}

/*}}}*/

#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define F SLANG_FLOAT_TYPE
#define S SLANG_STRING_TYPE

static SLang_IStruct_Field_Type Plot_User_Defaults_Table [] =
{
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, logx, "logx", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, logy, "logy", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, color, "color", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, start_color, "start_color", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, x_unit, "x_unit", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, linestyle, "line_style", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, start_linestyle, "start_line_style", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, line_width, "line_width", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, frame_line_width, "frame_line_width", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, pointstyle, "point_style", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, connect_points, "connect_points", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, char_height, "char_height", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ebar_term_length, "ebar_term_length", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, point_size, "point_size", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, use_errorbars, "use_errorbars", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, use_bin_density, "use_bin_density", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_xmin, "ovp_xmin", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_xmax, "ovp_xmax", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_ymin, "ovp_ymin", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_ymax, "ovp_ymax", F, 0),
   SLANG_END_ISTRUCT_TABLE
};

static SLang_CStruct_Field_Type Plot_Option_Table [] =
{
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, xmin, "xmin", F, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, xmax, "xmax", F, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, ymin, "ymin", F, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, ymax, "ymax", F, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, xlabel, "xlabel", S, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, ylabel, "ylabel", S, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, tlabel, "tlabel", S, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, xopt, "xopt", S, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, yopt, "yopt", S, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, logx, "logx", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, logy, "logy", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, color, "color", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, start_color, "start_color", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, x_unit, "x_unit", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, linestyle, "line_style", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, start_linestyle, "start_line_style", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, line_width, "line_width", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, frame_line_width, "frame_line_width", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, pointstyle, "point_style", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, connect_points, "connect_points", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, char_height, "char_height", F, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, point_size, "point_size", F, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, ebar_term_length, "ebar_term_length", F, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, use_errorbars, "use_errorbars", I, 0),
   MAKE_CSTRUCT_FIELD(Plot_Options_Type, use_bin_density, "use_bin_density", I, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_xmin, "ovp_xmin", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_xmax, "ovp_xmax", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_ymin, "ovp_ymin", F, 0),
   MAKE_ISTRUCT_FIELD(Plot_Options_Type, ovp_ymax, "ovp_ymax", F, 0),
   SLANG_END_CSTRUCT_TABLE
};

#undef V
#undef I
#undef F
#undef S

static void set_plot_options_struct (void) /*{{{*/
{
   Plot_Options_Type p;
   Plot_t *fmt;

   memset ((char *) &p, 0, sizeof(p));

   if (-1 == SLang_pop_cstruct ((VOID_STAR)&p, Plot_Option_Table))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't pop cstruct");
        return;
     }

   if (NULL == (fmt = current_format ()))
     return;

   if (-1 == Plot_set_options (fmt, &p))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "setting plot options");
     }

   SLang_free_cstruct ((VOID_STAR)&p, Plot_Option_Table);
}

/*}}}*/

static void get_plot_options_struct (void) /*{{{*/
{
   Plot_Options_Type p;
   Plot_t *fmt;

   memset ((char *)&p, 0, sizeof (p));

   if (NULL != (fmt = current_format ()))
     {
        if (-1 == Plot_get_options (fmt, &p))
          {
             isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "getting plot options");
          }
     }

   if (-1 == SLang_push_cstruct ((VOID_STAR)&p, Plot_Option_Table))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't push cstruct");
     }

   Plot_free_options(&p);
}

/*}}}*/

/* read cursor, place text */

static int cursor_interval (int (*fcn)(float *, float *, char *), float *cmin, float *cmax) /*{{{*/
{
   int ret;
   char ch = '\0';

   *cmin = *cmax = 0.0;

   fprintf (stdout, "Mark the interval with the cursor:\n");
   fprintf (stdout, "Type 'r' to re-select the first point.\n");
   fprintf (stdout, "Type 'q' to quit without marking any points.\n");
   fflush (stdout);

   do
     {
        ret = (*fcn) (cmin, cmax, &ch);
        if (ret == -1)
          isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading cursor");
     } while (ch == 'r' || ch == 'R');

   return ret;
}

/*}}}*/

static void check_plot_scale (int *is_logx, int *is_logy) /*{{{*/
{
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     return;

   if (is_logx != NULL)
     *is_logx = Plot_get_logx(fmt);
   if (is_logy != NULL)
     *is_logy = Plot_get_logy(fmt);
}

/*}}}*/

static void xinterval (void) /*{{{*/
{
   float xmin,xmax;
   int is_logx = 0;
   check_plot_scale (&is_logx, NULL);
   if (-1 == cursor_interval (Plot_read_x_interval_cursor, &xmin, &xmax))
     return;
   if (is_logx)
     {
        xmin = pow(10.0,xmin);
        xmax = pow(10.0,xmax);
     }
   (void) SLang_push_double ((double) MIN(xmin,xmax));
   (void) SLang_push_double ((double) MAX(xmin,xmax));
}

/*}}}*/

static void yinterval (void) /*{{{*/
{
   float ymin,ymax;
   int is_logy = 0;
   check_plot_scale (NULL, &is_logy);
   if (-1 == cursor_interval (Plot_read_y_interval_cursor, &ymin, &ymax))
     return;
   if (is_logy)
     {
        ymin = pow(10.0,ymin);
        ymax = pow(10.0,ymax);
     }
   (void) SLang_push_double ((double) MIN(ymin,ymax));
   (void) SLang_push_double ((double) MAX(ymin,ymax));
}

/*}}}*/

static void cursor (int *loop) /*{{{*/
{
   float x, y;
   char ch = '\0';
   int is_logx=0, is_logy=0;

   check_plot_scale (&is_logx, &is_logy);

   x = y = 0.0;

   for (;;)
     {
        ch = '\0';
        if (-1 == Plot_read_crosshair_cursor (&x, &y, &ch))
          break;
        if (ch == 'q' || ch == 'Q' || *loop == 0)
          goto done;
        if (is_logx)
          x = pow(10.0,x);
        if (is_logy)
          y = pow(10.0,y);
        fprintf (stdout, "x=%14.6e  y=%14.6e  ch=%c\n",x,y,ch);
        fflush(stdout);
     }

   isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading cursor position");

   done:

   (void) SLang_push_double ((double)x);
   (void) SLang_push_double ((double)y);
   (void) SLang_push_integer ((int)ch);
   SLang_flush_input();
   fputc ('\n', stdout);
   fflush(stdout);
}

/*}}}*/

static void _get_cursor_box (void) /*{{{*/
{
   int is_logx, is_logy;
   float bx[2] = { 0.0, 0.0 };
   float by[2] = { 0.0, 0.0 };

   check_plot_scale (&is_logx, &is_logy);

   if (-1 == Plot_read_box_corners_from_cursor (bx,by))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "reading plot cursor");

   if (is_logx)
     {
        bx[0] = pow(10.0,bx[0]);
        bx[1] = pow(10.0,bx[1]);
     }
   if (is_logy)
     {
        by[0] = pow(10.0,by[0]);
        by[1] = pow(10.0,by[1]);
     }

   (void) SLang_push_double ((double) MIN (bx[0], bx[1]));
   (void) SLang_push_double ((double) MAX (bx[0], bx[1]));
   (void) SLang_push_double ((double) MIN (by[0], by[1]));
   (void) SLang_push_double ((double) MAX (by[0], by[1]));
   SLang_flush_input();
   fputc ('\n', stdout);
   fflush (stdout);
}

/*}}}*/

static int fixup_log (float *x, int is_log, float *new_x) /*{{{*/
{
   if (is_log == 0)
     {
        *new_x = *x;
        return 0;
     }

   if (*x <= 0.0)
     return -1;

   *new_x = log10(*x);
   return 0;
}

/*}}}*/

static void _put_text (float *x, float *y, char *text, float *angle, float *justify) /*{{{*/
{
   float px, py;
   Plot_t *fmt;

   if (NULL == (fmt = current_format ()))
     return;

   Plot_set_color (fmt);
   Plot_set_charsize (fmt);
   Plot_set_line_width (fmt);

   if (-1 == fixup_log (x, Plot_get_logx(fmt), &px)
       || -1 == fixup_log (y, Plot_get_logy(fmt), &py))
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "invalid coordinate");
        return;
     }

   Plot_put_text (px, py, *angle, *justify, text);
}

/*}}}*/

static void _plot_box (void) /*{{{*/
{
   Plot_t *fmt;

   if ((NULL == (fmt = current_format ())
        || (force_open_plot_device () < 0)))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        return;
     }

   if (Plot_has_default_labels (fmt))
     Plot_set_labels (fmt, "", "", NULL);

   if (-1 == Plot_box (fmt))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "drawing box");
}

/*}}}*/

/* plot (x,y) curve or histogram  */

static void _plot_xy (void) /*{{{*/
{
   SLang_Array_Type *slx, *sly, *slsym;
   Plot_t *fmt;
   int style, overlay;
   int *psym;
   int *pstyle;

   slx = sly = slsym = NULL;

   if ((NULL == (fmt = current_format ())
        || (force_open_plot_device () < 0)))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        return;
     }

   if (-1 == SLang_pop_array_of_type (&slsym, SLANG_INT_TYPE)
       || slsym == NULL
       || -1 == SLang_pop_integer (&overlay)
       || -1 == SLang_pop_integer (&style)
       || -1 == SLang_pop_array_of_type (&sly, SLANG_FLOAT_TYPE)
       || sly == NULL
       || -1 == SLang_pop_array_of_type (&slx, SLANG_FLOAT_TYPE)
       || slx == NULL)
     goto free_and_return;

   if (slx->num_elements != sly->num_elements)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto free_and_return;
     }

   if (slsym->num_elements == sly->num_elements
       && sly->num_elements > 1)
     psym = (int *) slsym->data;
   else
     psym = NULL;

   if (!overlay)
     reset_current_pane (fmt);

   if (Plot_has_default_labels (fmt))
     Plot_set_labels (fmt, "", "", NULL);

   /* overrides auto-incr on overlay,
    * but won't override fmt style
    */
   pstyle = (style == INT_MAX) ? NULL : &style;

   if (-1 == Plot_curve ((float *) slx->data, (float *) sly->data,
                         slx->num_elements, psym, pstyle, overlay, fmt))
     isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
   else
     Plot_auto_incr_line_type (fmt);

   free_and_return:
   SLang_free_array (slx);
   SLang_free_array (sly);
   SLang_free_array (slsym);
}

/*}}}*/

static void _plot_hist (void) /*{{{*/
{
   SLang_Array_Type *slxlo, *slxhi, *sly;
   Isis_Hist_t h;
   Plot_t *fmt;
   int i, style, overlay;
   int *pstyle;

   slxlo = slxhi = sly = NULL;

   if ((NULL == (fmt = current_format ())
        || (force_open_plot_device() < 0)))
     {
        isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "plot failed");
        return;
     }

   if ( -1 == SLang_pop_integer (&overlay)
       || -1 == SLang_pop_integer (&style)
       || -1 == SLang_pop_array_of_type (&sly, SLANG_DOUBLE_TYPE)
       || sly == NULL
       || -1 == SLang_pop_array_of_type (&slxhi, SLANG_DOUBLE_TYPE)
       || slxhi == NULL
       || -1 == SLang_pop_array_of_type (&slxlo, SLANG_DOUBLE_TYPE)
       || slxlo == NULL)
     goto free_and_return;

   if (slxlo->num_elements != sly->num_elements
       || slxhi->num_elements != sly->num_elements)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "inconsistent array sizes");
        goto free_and_return;
     }

   h.nbins = (int) sly->num_elements;
   h.val = (double *) sly->data;
   h.val_err = NULL;    /* no bin-uncertainties */
   h.notice = NULL;     /* no bin-notice/ignore */
   h.bin_lo = (double *)slxlo->data;
   h.bin_hi = (double *)slxhi->data;

   for (i=0; i < (int) slxlo->num_elements; i++)
     {
        if (h.bin_lo[i] >= h.bin_hi[i])
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "inconsistent grid");
             goto free_and_return;
          }
     }

   if (!overlay)
     reset_current_pane (fmt);

   if (Plot_has_default_labels (fmt))
     Plot_set_labels (fmt, "", "", NULL);

   /* overrides auto-incr on overlay,
    * but won't override fmt style
    */
   pstyle = (style == INT_MAX) ? NULL : &style;

   if (-1 == Plot_hist (fmt, pstyle, overlay, &h))
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plotting histogram");
   else
     Plot_auto_incr_line_type (fmt);

   free_and_return:
   SLang_free_array (slxlo);
   SLang_free_array (slxhi);
   SLang_free_array (sly);
}

/*}}}*/

/* Plot location of other diffraction orders */
static void lambda_mth_order (float * m, float * lambda_n, float * n) /*{{{*/
{
   enum { ORDER_LABEL_SIZE=32 };
   char label[ORDER_LABEL_SIZE];
   float xlim[2] = {0.0, 0.0};
   float ylim[2] = {0.0, 0.0};
   float wl[2], y[2], yl;
   Plot_t *fmt;
   int logx;

   if ( *m <= 0.0 || *m > 99.0
       || *n <= 0.0 || *n > 99.0 )
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid diffraction order");
        return;
     }

   if (NULL == (fmt = current_format ()))
     return;

   Plot_set_color (fmt);
   Plot_set_charsize (fmt);

   wl[0] = wl[1] = (*lambda_n) * (*m) / (*n);

   logx = Plot_get_logx (fmt);

   if (logx && *lambda_n > 0.0)
     {
        wl[0] = log10 (wl[0]);
        wl[1] = log10 (wl[1]);
     }

   Plot_query_plot_limits (&xlim[0], &xlim[1], &ylim[0], &ylim[1]);
   y[0] = 0.6 * (ylim[1] - ylim[0]) + ylim[0];
   y[1] = 0.7 * (ylim[1] - ylim[0]) + ylim[0];

   Plot_line (2, wl, y);

   sprintf (label, "wl=%0.5f m=%0.2f n=%0.2f", *lambda_n, *m, *n);

   yl = y[1] + 0.25* (y[1] - y[0]);

   Plot_put_text (wl[1], yl,
                  (float) 90.0,        /* angle in degrees */
                  (float) 0.0,         /* justification, 0/0.5/1 = left/center/right */
                  label);
}

/*}}}*/

static void set_plot_library_interface (void) /*{{{*/
{
   if (-1 == Plot_set_library_interface ())
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't set library interface");
        return;
     }
}

/*}}}*/

static void push_plot_library_interface (void) /*{{{*/
{
   Plot_push_library_interface ();
}

/*}}}*/

/*{{{ Slang intrinsic function table */

static SLang_IConstant_Type Plot_Units_Const [] =
{
   MAKE_ICONSTANT("U_ANGSTROM",   U_ANGSTROM),
   MAKE_ICONSTANT("U_NANOMETER",  U_NANOMETER),
   MAKE_ICONSTANT("U_MILLIMETER", U_MILLIMETER),
   MAKE_ICONSTANT("U_CENTIMETER", U_CENTIMETER),
   MAKE_ICONSTANT("U_METER",      U_METER),
   MAKE_ICONSTANT("U_EV",         U_EV),
   MAKE_ICONSTANT("U_KEV",        U_KEV),
   MAKE_ICONSTANT("U_MEV",        U_MEV),
   MAKE_ICONSTANT("U_GEV",        U_GEV),
   MAKE_ICONSTANT("U_TEV",        U_TEV),
   MAKE_ICONSTANT("U_HZ",         U_HZ),
   MAKE_ICONSTANT("U_KHZ",        U_KHZ),
   MAKE_ICONSTANT("U_MHZ",        U_MHZ),
   MAKE_ICONSTANT("U_GHZ",        U_GHZ),
   SLANG_END_ICONST_TABLE
};

#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define F SLANG_FLOAT_TYPE
#define S SLANG_STRING_TYPE

static SLang_Intrin_Fun_Type Plot_Intrinsics [] =
{
   MAKE_INTRINSIC("_set_plot_library_interface", set_plot_library_interface, V, 0),
   MAKE_INTRINSIC("_get_plot_library_interface", push_plot_library_interface, V, 0),
   MAKE_INTRINSIC_3("_open_plot", _open_plot, I, S, I, I),
   MAKE_INTRINSIC_2("_copy_plot", _copy_plot, I, S, I),
   MAKE_INTRINSIC("_multiplot", multiplot, V, 0),
   MAKE_INTRINSIC("_close_plot", _close_plot, V, 0),
   MAKE_INTRINSIC("_plot_quit", _quit_plot, V, 0),
   MAKE_INTRINSIC("_erase_plot", _erase_plot, V, 0),
   MAKE_INTRINSIC("_plot_xy", _plot_xy, V, 0),
   MAKE_INTRINSIC("_plot_box", _plot_box, V, 0),
   MAKE_INTRINSIC("_hplot", _plot_hist, V, 0),
   MAKE_INTRINSIC_1("_cursor", cursor, V, I),
   MAKE_INTRINSIC("_xinterval", xinterval, V, 0),
   MAKE_INTRINSIC("_yinterval", yinterval, V, 0),
   MAKE_INTRINSIC("_get_cursor_box", _get_cursor_box, V, 0),
   MAKE_INTRINSIC_SSS("_set_labels", _set_labels, V),
   MAKE_INTRINSIC_S("_set_title", _set_title, V),
   MAKE_INTRINSIC_S("_set_xlabel", _set_xlabel, V),
   MAKE_INTRINSIC_S("_set_ylabel", _set_ylabel, V),
   MAKE_INTRINSIC_5("_put_text", _put_text, V, F, F, S, F, F),
   MAKE_INTRINSIC("_set_xrange", _set_xrange, V, 0),
   MAKE_INTRINSIC("_set_yrange", _set_yrange, V, 0),
   MAKE_INTRINSIC("_get_xrange", _get_xrange, V, 0),
   MAKE_INTRINSIC("_get_yrange", _get_yrange, V, 0),
   MAKE_INTRINSIC("_reset_plot_ranges", _reset_plot_ranges, V, 0),
   MAKE_INTRINSIC("_logx_on", _logx_on, V, 0),
   MAKE_INTRINSIC("_logx_off", _logx_off, V, 0),
   MAKE_INTRINSIC("_logy_on", _logy_on, V, 0),
   MAKE_INTRINSIC("_logy_off", _logy_off, V, 0),
   MAKE_INTRINSIC("_xaxis_is_log", _xaxis_is_log, V, 0),
   MAKE_INTRINSIC("_yaxis_is_log", _yaxis_is_log, V, 0),
   MAKE_INTRINSIC("_plot_bin_density", bin_density_on, V, 0),
   MAKE_INTRINSIC("_plot_bin_integral", bin_density_off, V, 0),
   MAKE_INTRINSIC_I("_set_connect_points", _set_connect_points, V),
   MAKE_INTRINSIC("_connect_points", _connect_points, V, 0),
   MAKE_INTRINSIC_2("_set_errorbar", _set_errorbar, V, I, F),
   MAKE_INTRINSIC("_errorbar_state", _errorbar_state, V, 0),
   MAKE_INTRINSIC_I("_set_window", _set_window, V),
   MAKE_INTRINSIC_I("_set_pane", _set_pane, V),
   MAKE_INTRINSIC_2("_resize_window", _resize_window, V, F, F),
   MAKE_INTRINSIC_1("_set_charsize", _set_charsize, V, F),
   MAKE_INTRINSIC_I("_set_style_meaning", _set_style_meaning, V),
   MAKE_INTRINSIC_I("_set_color", _set_color, V),
   MAKE_INTRINSIC("_get_color", _get_color, V, 0),
   MAKE_INTRINSIC_I("_set_linestyle", set_linestyle, V),
   MAKE_INTRINSIC_I("_set_pointstyle", set_pointstyle, V),
   MAKE_INTRINSIC_1("_set_point_size", set_point_size, V, F),
   MAKE_INTRINSIC_I("_set_frame_line_width", _set_frame_line_width, V),
   MAKE_INTRINSIC_I("_set_line_width", _set_line_width, V),
   MAKE_INTRINSIC("_get_line_width", _get_line_width, I, 0),
   MAKE_INTRINSIC_1("_set_plot_units", _set_plot_units, V, S),
   MAKE_INTRINSIC_1("_set_plot_device", _set_plot_device, V, S),
   MAKE_INTRINSIC_I("_set_auto_increment", _set_auto_increment, V),
   MAKE_INTRINSIC_3("_lambda_mth_order", lambda_mth_order, V, F, F, F),
   MAKE_INTRINSIC("_set_plot_options_struct", set_plot_options_struct, V, 0),
   MAKE_INTRINSIC("_get_plot_options_struct", get_plot_options_struct, V, 0),
   MAKE_INTRINSIC("set_outer_viewport_intrin", set_outer_viewport_intrin, V, 0),
   MAKE_INTRINSIC("get_outer_viewport_intrin", get_outer_viewport_intrin, V, 0),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef V
#undef I
#undef F
#undef S

SLANG_MODULE(plot);
int init_plot_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;
   static Plot_Options_Type *pud = &Plot_User_Defaults;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if ((-1 == SLns_add_istruct_table (pub_ns, Plot_User_Defaults_Table, &pud, "_isis_plot"))
       || (-1 == SLns_add_iconstant_table (pub_ns, Plot_Units_Const, NULL)))
     return isis_trace_return(-1);

   Plot_init_default_format ();

   return SLns_add_intrin_fun_table (ns, Plot_Intrinsics, NULL);
}

void deinit_plot_module (void)
{
   _quit_plot ();
}

/*}}}*/
