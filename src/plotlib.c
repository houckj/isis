/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008  Massachusetts Institute of Technology

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

#include "config.h"
#include <stdio.h>
#include <signal.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <slang.h>

#include "isis.h"
#include "plot.h"
#include "util.h"
#include "errors.h"

#define PLI_FIELD(n) \
    SLang_Name_Type *n; \
    SLang_Ref_Type *n##_ref

typedef struct
{
   PLI_FIELD(open);
   PLI_FIELD(close);
   PLI_FIELD(subdivide);
   PLI_FIELD(select_window);
   PLI_FIELD(select_viewport);
   PLI_FIELD(set_plot_limits);
   PLI_FIELD(query_plot_limits);
   PLI_FIELD(erase);
   PLI_FIELD(update);
   PLI_FIELD(next_page);
   PLI_FIELD(get_color);
   PLI_FIELD(set_color);
   PLI_FIELD(set_line_style);
   PLI_FIELD(set_line_width);
   PLI_FIELD(set_clipping);
   PLI_FIELD(plot_xy);
   PLI_FIELD(plot_points);
   PLI_FIELD(plot_symbol_points);
   PLI_FIELD(plot_histogram);
   PLI_FIELD(plot_y_errorbar);
   PLI_FIELD(set_viewer_size);
   PLI_FIELD(set_char_size);
   PLI_FIELD(draw_box);
   PLI_FIELD(label_axes);
   PLI_FIELD(put_text_xy);
   PLI_FIELD(put_text_offset);
   PLI_FIELD(default_axis);
   PLI_FIELD(configure_axis);
   PLI_FIELD(read_cursor);
}
Plot_Library_Interface_Type;

#undef PLI_FIELD

static SLang_CStruct_Field_Type Plot_Library_Interface_Table [] =
{
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, open_ref, "open", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, close_ref, "close", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, subdivide_ref, "subdivide", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, select_window_ref, "select_window", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, select_viewport_ref, "select_viewport", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, set_plot_limits_ref, "set_plot_limits", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, query_plot_limits_ref, "query_plot_limits", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, erase_ref, "erase", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, update_ref, "update", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, next_page_ref, "next_page", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, get_color_ref, "get_color", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, set_color_ref, "set_color", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, set_line_style_ref, "set_line_style", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, set_line_width_ref, "set_line_width", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, set_clipping_ref, "set_clipping", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, plot_xy_ref, "plot_xy", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, plot_points_ref, "plot_points", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, plot_symbol_points_ref, "plot_symbol_points", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, plot_histogram_ref, "plot_histogram", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, plot_y_errorbar_ref, "plot_y_errorbar", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, set_viewer_size_ref, "set_viewer_size", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, set_char_size_ref, "set_char_size", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, draw_box_ref, "draw_box", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, label_axes_ref, "label_axes", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, put_text_xy_ref, "put_text_xy", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, put_text_offset_ref, "put_text_offset", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, default_axis_ref, "default_axis", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, configure_axis_ref, "configure_axis", SLANG_REF_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plot_Library_Interface_Type, read_cursor_ref, "read_cursor", SLANG_REF_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static Plot_Library_Interface_Type *PLI;

static int pli_undefined (void) /*{{{*/
{
   if (PLI == NULL)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "plot interface is undefined");
        return 1;
     }
   return 0;
}

/*}}}*/

void Plot_free_library_interface (void) /*{{{*/
{
   if (PLI == NULL)
     return;

#define PLI_FREE(field) \
   SLang_free_function(PLI->field); \
   SLang_free_ref(PLI->field##_ref);

   PLI_FREE(open);
   PLI_FREE(close);
   PLI_FREE(subdivide);
   PLI_FREE(select_window);
   PLI_FREE(select_viewport);
   PLI_FREE(set_plot_limits);
   PLI_FREE(query_plot_limits);
   PLI_FREE(erase);
   PLI_FREE(update);
   PLI_FREE(next_page);
   PLI_FREE(get_color);
   PLI_FREE(set_color);
   PLI_FREE(set_line_style);
   PLI_FREE(set_clipping);
   PLI_FREE(set_line_width);
   PLI_FREE(plot_xy);
   PLI_FREE(plot_points);
   PLI_FREE(plot_symbol_points);
   PLI_FREE(plot_histogram);
   PLI_FREE(plot_y_errorbar);
   PLI_FREE(set_viewer_size);
   PLI_FREE(set_char_size);
   PLI_FREE(draw_box);
   PLI_FREE(label_axes);
   PLI_FREE(put_text_xy);
   PLI_FREE(put_text_offset);
   PLI_FREE(default_axis);
   PLI_FREE(configure_axis);
   PLI_FREE(read_cursor);
#undef PLI_FREE

   ISIS_FREE(PLI);
}

/*}}}*/

int Plot_set_library_interface (void) /*{{{*/
{
   Plot_Library_Interface_Type *pli;
   int status = -1;

   if (NULL == (pli = (Plot_Library_Interface_Type *) ISIS_MALLOC (sizeof *pli)))
     return -1;

   if (-1 == SLang_pop_cstruct ((VOID_STAR)pli, Plot_Library_Interface_Table))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't set library interface");
        return -1;
     }

#define PLI_SET(field) \
    do { \
       if (pli->field##_ref == NULL) \
         goto return_status; \
       else { \
          pli->field = SLang_get_fun_from_ref (pli->field##_ref); \
           if (pli->field == NULL) \
              goto return_status; \
        } \
    } while (0)

   PLI_SET(open);
   PLI_SET(close);
   PLI_SET(subdivide);
   PLI_SET(select_window);
   PLI_SET(select_viewport);
   PLI_SET(set_plot_limits);
   PLI_SET(query_plot_limits);
   PLI_SET(erase);
   PLI_SET(update);
   PLI_SET(next_page);
   PLI_SET(get_color);
   PLI_SET(set_color);
   PLI_SET(set_line_style);
   PLI_SET(set_clipping);
   PLI_SET(set_line_width);
   PLI_SET(plot_xy);
   PLI_SET(plot_points);
   PLI_SET(plot_symbol_points);
   PLI_SET(plot_histogram);
   PLI_SET(plot_y_errorbar);
   PLI_SET(set_viewer_size);
   PLI_SET(set_char_size);
   PLI_SET(draw_box);
   PLI_SET(label_axes);
   PLI_SET(put_text_xy);
   PLI_SET(put_text_offset);
   PLI_SET(default_axis);
   PLI_SET(configure_axis);
   PLI_SET(read_cursor);

#undef PLI_SET

   Plot_free_library_interface ();
   PLI = pli;

   status = 0;
   return_status:

   if (status)
     {
        SLang_free_cstruct ((VOID_STAR)pli, Plot_Library_Interface_Table);
        ISIS_FREE(pli);
     }

   return status;
}

/*}}}*/

void Plot_push_library_interface (void) /*{{{*/
{
   Plot_Library_Interface_Type pli;

   if (PLI == NULL)
     memset ((char *)&pli, 0, sizeof pli);
   else
     {
        /* struct copy */
        pli = *PLI;
     }

   if (-1 == SLang_push_cstruct ((VOID_STAR)&pli, Plot_Library_Interface_Table))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "couldn't set library interface");
        return;
     }
}

/*}}}*/

static volatile sig_atomic_t Signal_In_Progress;
static void sig_segv (int signo) /*{{{*/
{
   static char msg[] =
"\n**** Segmentation Fault occurred while closing plot device\n";
   (void) signo;
   if (Signal_In_Progress)
     return;
   Signal_In_Progress = 1;
   write (STDERR_FILENO, msg, sizeof(msg));
   /* so more SEGVs won't interfere with _exit() */
   SLsignal (SIGSEGV, SIG_DFL);
   _exit (EXIT_FAILURE);
}

/*}}}*/

int Plot_close (void) /*{{{*/
{
   int status;
   SLSig_Fun_Type *sig_func;

   if (pli_undefined())
     return -1;

   if (PLI->close == NULL)
     return -1;

   sig_func = SLsignal (SIGSEGV, sig_segv);
   if (SIG_ERR == sig_func)
     fprintf (stderr, "warning:  failed initializing signal handler for SIGSEGV\n");

   if (-1 == SLexecute_function (PLI->close))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed closing plot device");
        return -1;
     }
   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "function closing plot device did not return status integer");
        return -1;
     }

   if (SLsignal (SIGSEGV, sig_func) == SIG_ERR)
     fprintf (stderr, "warning:  failed to re-set signal handler\n");

   return status;
}

/*}}}*/

int Plot_open (char *device) /*{{{*/
{
   int id;

   if (pli_undefined())
     return -1;

   if (PLI->open == NULL)
     return -1;

   SLang_start_arg_list ();
   SLang_push_string (device);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->open))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed opening plot device");
        return -1;
     }

   if ((-1 == SLang_pop_integer (&id))
       || (id <= 0))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed opening plot device");
        return -1;
     }

   return id;
}

/*}}}*/

int Plot_subdivide (int num_x_subpanels, int num_y_subpanels) /*{{{*/
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->subdivide == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: subdivide operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_integer (num_x_subpanels);
   SLang_push_integer (num_y_subpanels);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->subdivide))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed subdividing plot device");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed subdividing plot device");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_select_window (int device) /*{{{*/
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->select_window == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: select_window operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_integer (device);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->select_window))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed selecting plot device");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed selecting plot device");
        return -1;
     }

   return status;
}

/*}}}*/

#define CALL_FFFF_FUN(op,xmin,xmax,ymin,ymax) \
   int status; \
   if (pli_undefined()) \
     return -1;  \
   if (PLI->op == NULL) \
     { \
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, \
                    "plot: " #op " operation is not supported"); \
        return -1; \
     } \
   SLang_start_arg_list (); \
   SLang_push_float (xmin); \
   SLang_push_float (xmax); \
   SLang_push_float (ymin); \
   SLang_push_float (ymax); \
   SLang_end_arg_list (); \
   if (-1 == SLexecute_function (PLI->op)) \
     { \
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: " #op " failed"); \
        return -1; \
     } \
   if (-1 == SLang_pop_integer (&status)) \
     { \
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: " #op " failed"); \
        return -1; \
     } \
   return status

int Plot_select_viewport (float xmin, float xmax, float ymin, float ymax) /*{{{*/
{
   CALL_FFFF_FUN(select_viewport,xmin,xmax,ymin,ymax);
}

/*}}}*/

int _Plot_set_plot_limits (float xmin, float xmax, float ymin, float ymax) /*{{{*/
{
   CALL_FFFF_FUN(set_plot_limits,xmin,xmax,ymin,ymax);
}

/*}}}*/

#define CALL_V_FUN(op) \
   int status; \
   if (pli_undefined()) \
     return -1;  \
   if (PLI->op == NULL) \
     { \
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, \
                    "plot: " #op " operation is not supported"); \
        return -1; \
     } \
   if (-1 == SLexecute_function (PLI->op)) \
     { \
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: " #op " failed"); \
        return -1; \
     } \
   if (-1 == SLang_pop_integer (&status)) \
     { \
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: " #op " failed"); \
        return -1; \
     } \
   return status

int Plot_erase (void) /*{{{*/
{
   CALL_V_FUN(erase);
}

/*}}}*/

int Plot_update_display (void) /*{{{*/
{
   CALL_V_FUN(update);
}

/*}}}*/

int _Plot_next_page (void) /*{{{*/
{
   CALL_V_FUN(next_page);
}

/*}}}*/

int _Plot_get_color (void) /*{{{*/
{
   CALL_V_FUN(get_color);
}

/*}}}*/

#define SET_INT_PROPERTY(op,prop) \
   int status; \
   if (pli_undefined()) \
     return -1; \
   if (PLI->op == NULL) \
     { \
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, \
                    "plot: " #op " operation is not supported"); \
        return -1; \
     } \
   SLang_start_arg_list (); \
   SLang_push_integer (prop); \
   SLang_end_arg_list (); \
   if (-1 == SLexecute_function (PLI->op)) \
     { \
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed setting " #prop); \
        return -1; \
     } \
   if (-1 == SLang_pop_integer (&status)) \
     { \
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed setting " #prop); \
        return -1; \
     } \
   return status

int _Plot_set_color (int color) /*{{{*/
{
   SET_INT_PROPERTY(set_color,color);
}

/*}}}*/

int _Plot_set_linestyle (int line_style) /*{{{*/
{
   SET_INT_PROPERTY(set_line_style,line_style);
}

/*}}}*/

int Plot_set_lineclip_state (int clipping) /*{{{*/
{
   SET_INT_PROPERTY(set_clipping,clipping);
}

/*}}}*/

int _Plot_set_line_width (int line_width) /*{{{*/
{
   SET_INT_PROPERTY(set_line_width,line_width);
}

/*}}}*/

static int push_two_float_arrays (int n, float *a, float *b) /*{{{*/
{
   SLang_Array_Type *sl_a=NULL, *sl_b=NULL;
   int status = -1;

   if (a == NULL || b == NULL)
     return -1;

   if ((NULL == (sl_a = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &n, 1)))
       || NULL == (sl_b = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &n, 1)))
     goto return_status;

   memcpy ((char *)sl_a->data, (char *)a, n * sizeof(float));
   memcpy ((char *)sl_b->data, (char *)b, n * sizeof(float));

   SLang_push_array (sl_a, 1);
   SLang_push_array (sl_b, 1);

   status = 0;
return_status:
   if (status)
     {
        SLang_free_array (sl_a);
        SLang_free_array (sl_b);
     }

   return status;
}

/*}}}*/

static int push_three_float_arrays (int n, float *a, float *b, float *c) /*{{{*/
{
   SLang_Array_Type *sl_a=NULL, *sl_b=NULL, *sl_c=NULL;
   int status = -1;

   if (a == NULL || b == NULL || c == NULL)
     return -1;

   if ((NULL == (sl_a = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &n, 1)))
       || NULL == (sl_b = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &n, 1))
       || NULL == (sl_c = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &n, 1)))
     goto return_status;

   memcpy ((char *)sl_a->data, (char *)a, n * sizeof(float));
   memcpy ((char *)sl_b->data, (char *)b, n * sizeof(float));
   memcpy ((char *)sl_c->data, (char *)c, n * sizeof(float));

   SLang_push_array (sl_a, 1);
   SLang_push_array (sl_b, 1);
   SLang_push_array (sl_c, 1);

   status = 0;
return_status:
   if (status)
     {
        SLang_free_array (sl_a);
        SLang_free_array (sl_b);
        SLang_free_array (sl_c);
     }

   return status;
}

/*}}}*/

int Plot_line (int n, float *x, float *y) /*{{{*/
{
   int status = -1;

   if (pli_undefined())
     return -1;

   if (PLI->plot_xy == NULL)
     return -1;

   SLang_start_arg_list ();
   status = push_two_float_arrays (n, x, y);
   SLang_end_arg_list ();

   if ((status < 0) || (-1 == SLexecute_function (PLI->plot_xy)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed plotting line");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed plotting line");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_points (int n, float *x, float *y, int symbol) /*{{{*/
{
   int status = -1;

   if (pli_undefined())
     return -1;

   if (PLI->plot_points == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: plot_points operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   status = push_two_float_arrays (n, x, y);
   SLang_push_integer (symbol);
   SLang_end_arg_list ();

   if ((status < 0) || (-1 == SLexecute_function (PLI->plot_points)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed plotting points");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed plotting points");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_symbol_points (int n, float *x, float *y, int *symbol) /*{{{*/
{
   SLang_Array_Type *sl_sym=NULL;
   int status = -1;

   if (pli_undefined())
     return -1;

   if (PLI->plot_symbol_points == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: plot_symbol_points operation is not supported");
        return -1;
     }

   if (NULL == (sl_sym = SLang_create_array (SLANG_INT_TYPE, 0, NULL, &n, 1)))
     return -1;
   memcpy ((char *)sl_sym->data, (char *)symbol, n * sizeof(int));

   SLang_start_arg_list ();
   status = push_two_float_arrays (n, x, y);
   SLang_push_array (sl_sym, 1);
   SLang_end_arg_list ();

   if ((status < 0) || (-1 == SLexecute_function (PLI->plot_symbol_points)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed plotting points");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed plotting points");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_query_plot_limits (float *xmin, float *xmax, float *ymin, float *ymax) /*{{{*/
{
   if (pli_undefined())
     return -1;

   if (PLI->query_plot_limits == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: query_plot_limits operation is not supported");
        return -1;
     }

   if (-1 == SLexecute_function (PLI->query_plot_limits))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed querying plot limits");
        return -1;
     }

   if ((-1 == SLang_pop_float (ymax))
       || (-1 == SLang_pop_float (ymin))
       || (-1 == SLang_pop_float (xmax))
       || (-1 == SLang_pop_float (xmin)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed querying plot limits");
        return -1;
     }

   return 0;
}

/*}}}*/

int Plot_histogram_data (int n, float *lo, float *hi, float *val) /*{{{*/
{
   int status = -1;

   if (pli_undefined())
     return -1;

   if (PLI->plot_histogram == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: plot_histogram operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   status = push_three_float_arrays (n, lo, hi, val);
   SLang_end_arg_list ();

   if ((status < 0) || (-1 == SLexecute_function (PLI->plot_histogram)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed plotting histogram");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed plotting histogram");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_y_errorbar (int n, float *x, float *top, float *bot, /*{{{*/
                     float terminal_length)
{
   int status = -1;

   if (pli_undefined())
     return -1;

   if (PLI->plot_y_errorbar == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: plot_y_errorbar operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   status = push_three_float_arrays (n, x, top, bot);
   SLang_push_float (terminal_length);
   SLang_end_arg_list ();

   if ((status < 0) || (-1 == SLexecute_function (PLI->plot_y_errorbar)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "failed plotting Y errorbar");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "failed plotting Y errorbar");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_reset_viewer_size (float width_cm, float aspect) /*{{{*/
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->set_viewer_size == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: set_viewer_size operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_float (width_cm);
   SLang_push_float (aspect);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->set_viewer_size))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: set_viewer_size failed");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: set_viewer_size failed");
        return -1;
     }

   return status;
}

/*}}}*/

int _Plot_set_charsize (float size) /*{{{*/
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->set_char_size == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: set_char_size operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_float (size);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->set_char_size))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: set_char_size failed");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: set_char_size failed");
        return -1;
     }

   return status;
}

/*}}}*/

int _Plot_draw_box (char *xopt, float xtick, int nxsub, /*{{{*/
                     char *yopt, float ytick, int nysub)
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->draw_box == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: draw_box operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_string (xopt);
   SLang_push_float (xtick);
   SLang_push_integer (nxsub);
   SLang_push_string (yopt);
   SLang_push_float (ytick);
   SLang_push_integer (nysub);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->draw_box))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: draw_box failed");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: draw_box failed");
        return -1;
     }

   return status;
}

/*}}}*/

int _Plot_label_box (char *xlabel, char *ylabel, char *tlabel) /*{{{*/
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->label_axes == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: label_axes operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_string (xlabel);
   SLang_push_string (ylabel);
   SLang_push_string (tlabel);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->label_axes))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: label_axes failed");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: label_axes failed");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_put_text (float x, float y, float angle, float justify, char *txt) /*{{{*/
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->put_text_xy == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: put_text_xy operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_float (x);
   SLang_push_float (y);
   SLang_push_float (angle);
   SLang_push_float (justify);
   SLang_push_string (txt);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->put_text_xy))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: put_text_xy failed");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: put_text_xy failed");
        return -1;
     }

   return status;
}

/*}}}*/

int Plot_put_text_offset (char *where, float offset, float ox, float oy, char *text) /*{{{*/
{
   int status;

   if (pli_undefined())
     return -1;

   if (PLI->put_text_offset == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: put_text_offset operation is not supported");
        return -1;
     }

   SLang_start_arg_list ();
   SLang_push_string (where);
   SLang_push_float (offset);
   SLang_push_float (ox);
   SLang_push_float (oy);
   SLang_push_string (text);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->put_text_offset))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: put_text_offset failed");
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: put_text_offset failed");
        return -1;
     }

   return status;
}

/*}}}*/

char *Plot_get_option_string_default (void) /*{{{*/
{
   char *s;

   if (pli_undefined())
     return NULL;

   if (PLI->default_axis == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: default_axis operation is not supported");
        return NULL;
     }

   if (-1 == SLexecute_function (PLI->default_axis))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: default_axis failed");
        return NULL;
     }

   if (-1 == SLpop_string (&s))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: default_axis failed");
        return NULL;
     }

   return s;
}

/*}}}*/

char *Plot_configure_axis (char *opt, int is_log, int has_numbers) /*{{{*/
{
   char *s;

   if (pli_undefined())
     return NULL;

   if (PLI->configure_axis == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: configure_axis operation is not supported");
        return NULL;
     }

   SLang_start_arg_list ();
   SLang_push_string (opt);
   SLang_push_integer (is_log);
   SLang_push_integer (has_numbers);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (PLI->configure_axis))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: configure_axis failed");
        return NULL;
     }

   if (-1 == SLpop_string (&s))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: configure_axis failed");
        return NULL;
     }

   return s;
}

/*}}}*/

typedef struct
{
   float xanchor;
   float yanchor;
   int position;
   int mode;
}
Cursor_Config_Type;

static SLang_CStruct_Field_Type Cursor_Config_Table [] =
{
   MAKE_CSTRUCT_FIELD(Cursor_Config_Type, xanchor, "xanchor", SLANG_FLOAT_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Cursor_Config_Type, yanchor, "yanchor", SLANG_FLOAT_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Cursor_Config_Type, position, "position", SLANG_INT_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Cursor_Config_Type, mode, "mode", SLANG_INT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int read_plot_cursor (Cursor_Config_Type *cc, float *x, float *y, char *ch) /*{{{*/
{
   if (pli_undefined())
     return -1;

   if (PLI->read_cursor == NULL)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "plot: read_cursor operation is not supported");
        return -1;
     }

   if (-1 == SLang_push_cstruct ((VOID_STAR)cc, Cursor_Config_Table))
     return -1;

   if (-1 == SLexecute_function (PLI->read_cursor))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "plot: read_cursor failed");
        return -1;
     }

   if (SLstack_depth() < 3)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "plot: read_cursor failed");
        return -1;
     }

   SLang_pop_char (ch);
   SLang_pop_float (y);
   SLang_pop_float (x);

   return 0;
}

/*}}}*/

int Plot_read_rectangle_cursor (float xanchor, float yanchor, /*{{{*/
                                float *xfinal, float *yfinal,
                                char *ch)
{
   Cursor_Config_Type c;

   c.xanchor = xanchor;
   c.yanchor = yanchor;
   c.mode = 2;
   /* mode = 2 generates a rectangular box */
   c.position = 0;
   /* posn=0 means after reading the cursor, leave the cursor at
    * its current position */

   return read_plot_cursor (&c, xfinal, yfinal, ch);
}

/*}}}*/

int Plot_cursor_read_other_endpoint (float xanchor, float yanchor, /*{{{*/
                                     float *xfinal, float *yfinal,
                                     char *ch)
{
   Cursor_Config_Type c;

   c.xanchor = xanchor;
   c.yanchor = yanchor;
   c.mode = 1;
   /* mode = 1 generates a line from anchor to cursor */
   c.position = 0;
   /* posn=0 means after reading the cursor, leave the cursor at
    * its current position */

   return read_plot_cursor (&c, xfinal, yfinal, ch);
}

/*}}}*/

int Plot_read_crosshair_cursor (float *xfinal, float *yfinal, /*{{{*/
                                char *ch)
{
   Cursor_Config_Type c;

   c.xanchor = 0.0;
   c.yanchor = 0.0;
   c.mode = 7;
   /* mode = 7 generates a crosshair cursor */
   c.position = 0;
   /* posn=0 means after reading the cursor, leave the cursor at
    * its current position */

   return read_plot_cursor (&c, xfinal, yfinal, ch);
}

/*}}}*/

int Plot_read_y_interval_cursor (float *ymin, float *ymax, char *ch) /*{{{*/
{
   Cursor_Config_Type c;
   float x1=0.0, x2=0.0;

   c.xanchor = 0.0;
   c.yanchor = 0.0;
   c.mode = 5;
   c.position = 0;

   if (-1 == read_plot_cursor (&c, &x1, ymin, ch))
     return -1;

   if (*ch == 'q' || *ch == 'Q')
     return 1;

   c.xanchor = x1;
   c.yanchor = *ymin;
   c.mode = 3;
   c.position = 1;

   *ymax = *ymin;
   x2 = x1;

   return read_plot_cursor (&c, &x2, ymax, ch);
}

/*}}}*/

int Plot_read_x_interval_cursor (float *xmin, float *xmax, char *ch) /*{{{*/
{
   Cursor_Config_Type c;
   float ya=0.0, yb=0.0;

   c.xanchor = 0.0;
   c.yanchor = 0.0;
   c.mode = 6;
   c.position = 0;

   if (-1 == read_plot_cursor (&c, xmin, &ya, ch))
     return -1;

   if (*ch == 'q' || *ch == 'Q')
     return 1;

   c.xanchor = *xmin;
   c.yanchor = ya;
   c.mode = 4;
   c.position = 1;

   *xmax = *xmin;
   yb = ya;

   return read_plot_cursor (&c, xmax, &yb, ch);
}

/*}}}*/

