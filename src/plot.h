#ifndef ISIS_PLOT_H
#define ISIS_PLOT_H

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

/* $Id: plot.h,v 1.13 2004/08/28 13:03:52 houck Exp $ */

#define PLOT_NO_DEVICE INT_MAX

#define PLOT_DEFAULT_AUTO_INCR     1       /* on by default */
#define PLOT_LINE_COLOR            1       /* 1 = white */
#define PLOT_POINT                -1       /* symbol index */
#define PLOT_SOLID_LINE            1       /* line style index */
#define PLOT_LABEL_STRING_SIZE   256

#define MAX_COLORS 16       /* max line colors */
#define MAX_LINES   5       /* max line styles */

#define DEFAULT_ELEV_SCALE  5.0

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

enum
{
   PLOT_USER_LABEL,
   PLOT_DEFAULT_LABEL
};

typedef struct Plot_Node_t Plot_Node_t;
typedef struct Plot_t Plot_t;

typedef struct
{
   float xmin, xmax, ymin, ymax;
   float ovp_xmin, ovp_xmax, ovp_ymin, ovp_ymax;
   float xtick, ytick;
   int logx, logy;
   int color, start_color;
   int linestyle, start_linestyle;
   int line_width, frame_line_width;
   int pointstyle, connect_points;
   float char_height, point_size, ebar_term_length;
   int use_errorbars;
   int use_bin_density;
   int x_unit;
   int label_type;
   char *xlabel;
   char *ylabel;
   char *tlabel;
   char *xopt;
   char *yopt;
}
Plot_Options_Type;

extern Plot_Options_Type Plot_User_Defaults;

extern int Plot_set_library_interface (void);
extern void Plot_push_library_interface (void);
extern void Plot_free_library_interface (void);

/* stuff in plot.c */

extern void Plot_init_default_format (void);
extern Plot_Node_t *Plot_start (void);
extern void Plot_end (Plot_Node_t ** head);

extern void Plot_free_fmt (Plot_t *fmt);
extern Plot_t *Plot_alloc_fmt (void);

extern Plot_Node_t *Plot_new_window (Plot_Node_t *head);
extern int Plot_open_window (Plot_Node_t *n, char *device);
extern int Plot_reformat_node (Plot_Node_t *n, int nxpanes, int nypanes, int type, int *ysizes);
extern int Plot_copy_formats (Plot_t *src, Plot_t *dst);

extern Plot_Node_t *Plot_find_node (int id, Plot_Node_t *head);
extern Plot_t *Plot_fmt (Plot_Node_t *n);
extern Plot_Node_t *Plot_this_window (Plot_t *fmt);
extern int Plot_next_page (Plot_Node_t *n);

extern Plot_t *Plot_find_fmt_for_pane (int id, int xpane, int ypane,
                                      Plot_Node_t *head);
extern int Plot_delete_window (int id, Plot_Node_t *head);
extern int Plot_get_num_xpanes (int *nxpanes, int id, Plot_Node_t *head);

extern int Plot_set_options (Plot_t *fmt, Plot_Options_Type *p);
extern int Plot_get_options (Plot_t *fmt, Plot_Options_Type *p);
extern void Plot_free_options (Plot_Options_Type *p);

extern void Plot_set_scale_for_elev_x_axis (Plot_t *fmt, float scale);
extern void Plot_get_scale_for_elev_x_axis (float *scale, Plot_t *fmt);

extern void Plot_set_logx (Plot_t *fmt, int logx);
extern void Plot_set_logy (Plot_t *fmt, int logy);
extern int Plot_get_logx (Plot_t *fmt);
extern int Plot_get_logy (Plot_t *fmt);

extern int Plot_set_xrange (Plot_t *fmt, float *xmin, float *xmax);
extern int Plot_set_yrange (Plot_t *fmt, float *ymin, float *ymax);

extern int Plot_get_xrange (Plot_t *fmt, float **xmin, float **xmax);
extern int Plot_get_yrange (Plot_t *fmt, float **ymin, float **ymax);

extern void Plot_reset_xyranges (Plot_t *fmt);
extern int Plot_set_fmt_charsize (Plot_t *fmt, float char_height);
extern int Plot_set_labels (Plot_t *fmt,
                            char * xlabel, char * ylabel, char * tlabel);
extern int Plot_has_default_labels (Plot_t *fmt);
extern int Plot_set_label_type (Plot_t *fmt, int type);

extern int Plot_set_xlabel (Plot_t *fmt, char *xlabel);
extern int Plot_set_ylabel (Plot_t *fmt, char *ylabel);
extern int Plot_set_title (Plot_t *fmt, char *title);

extern void Plot_set_linestyle (Plot_t *fmt);
extern void Plot_set_line_width (Plot_t *fmt);
extern void Plot_set_frame_line_width (Plot_t *fmt);
extern void Plot_set_color (Plot_t *fmt);
extern void Plot_set_charsize (Plot_t *fmt);
extern void Plot_label_box (Plot_t *fmt);
extern void Plot_draw_box_opt (Plot_t *fmt, float xmin, float xmax, float ymin, float ymax);
extern int Plot_set_axis_options (Plot_t *fmt, char *xopt, char *yopt);
extern void Plot_set_axis_type (Plot_t *fmt, int axis);

extern void Plot_auto_incr_line_type (Plot_t *fmt);
extern void Plot_set_auto_increment (int flag);

extern void Plot_set_style_meaning (Plot_t *fmt, int style);
extern void Plot_set_style (Plot_t *fmt, int *style);
extern int Plot_set_fmt_pointstyle (Plot_t *fmt, int pointstyle);
extern int Plot_set_fmt_linestyle (Plot_t *fmt, int linestyle);
extern int Plot_set_fmt_point_size (Plot_t *fmt, float point_size);
extern int Plot_set_fmt_color (Plot_t *fmt, int color);
extern int Plot_set_fmt_line_width (Plot_t *fmt, int width);
extern int Plot_set_fmt_frame_line_width (Plot_t *fmt, int width);
extern int Plot_get_fmt_line_width (Plot_t *fmt);

extern int Plot_get_fmt_color (Plot_t *fmt);

extern int Plot_bin_density (Plot_t *fmt);
extern void Plot_set_bin_density (Plot_t *fmt, int type);

extern int Plot_x_unit (Plot_t *fmt);
extern int Plot_set_x_unit (Plot_t *fmt, int unit);
extern int Plot_x_unit_is_default (Plot_t *fmt);
extern int Plot_xtick_labels_off (Plot_t *fmt);

extern void Plot_set_errorbars (Plot_t *fmt, int value, float term_length);
extern int Plot_use_errorbars (Plot_t *fmt);

extern void Plot_set_connect_points (Plot_t *fmt, int value);
extern int Plot_connect_points (Plot_t *fmt);

extern int Plot_curve (float *x, float *y, int npts, int *symbol,
                       int *use_style, int overlay, Plot_t *fmt);
extern int Plot_box (Plot_t *fmt);
extern int Plot_hist (Plot_t *fmt, int *use_style, int overlay, Isis_Hist_t *h);

extern int Plot_read_box_corners_from_cursor (float bx[/*2*/], float by[/*2*/]);
extern int Plot_read_line_endpoints (float bx[/*2*/], float by[/*2*/], int draw_line);

extern int Plot_save_exposure_time (Plot_t *fmt, float t);
extern float Plot_get_exposure_time (Plot_t *fmt);

/* low level plot wrappers */

extern int Plot_open (char * device);
extern int Plot_subdivide (int num_x_subpanels, int num_y_subpanels);
extern int Plot_copy_window (char *new_device, int id, Plot_Node_t *head);

extern int Plot_close (void);
extern int Plot_erase (void);

extern int Plot_select_window (int device_id);
extern int Plot_visit_pane (Plot_Node_t *n, int pane);
extern int Plot_select_viewport (float xmin, float xmax, float ymin, float ymax);
extern void Plot_get_outer_viewport (Plot_Node_t *n, float *xmin, float *xmax, float *ymin, float *ymax);
extern int Plot_set_outer_viewport (Plot_Node_t *n, float xmin, float xmax, float ymin, float ymax);

extern int Plot_update_display (void);

extern int Plot_reset_viewer_size (float width_cm, float aspect);

extern char *Plot_get_option_string_default (void);
extern char *Plot_configure_axis (char *opt, int is_log, int has_numbers);
extern int Plot_set_lineclip_state (int state);
extern int Plot_query_plot_limits (float * xmin, float * xmax,
                                    float * ymin, float * ymax);

extern int Plot_put_text (float x, float y, float angle,
                          float justify, char *txt);

extern int Plot_line (int n, float * x, float * y);
extern int Plot_points (int n, float * x, float * y, int symbol);
extern int Plot_symbol_points (int n, float * x, float * y, int * symbol);

extern int Plot_histogram_data (int nbins, float *lo, float *hi, float *val);
extern int Plot_y_errorbar (int nbins, float *x, float *top, float *bot,
                             float terminal_length);

extern int Plot_read_crosshair_cursor (float * xfinal, float * yfinal,
                                       char * ch);
extern int Plot_read_rectangle_cursor (float xanchor, float yanchor,
                                       float * xfinal, float * yfinal,
                                       char * ch);
extern int Plot_cursor_read_other_endpoint (float xanchor, float yanchor,
                                            float * xfinal, float * yfinal,
                                            char * ch);
extern int Plot_read_x_interval_cursor (float *xmin, float *xmax, char *ch);
extern int Plot_read_y_interval_cursor (float *ymin, float *ymax, char *ch);

/* library internal use only, more or less */

extern int _Plot_next_page (void);
extern int _Plot_set_plot_limits (float xmin, float xmax, float ymin, float ymax);
extern int _Plot_label_box (char * xlabel, char * ylabel, char * tlabel);
extern int Plot_put_text_offset (char *where, float offset, float ox, float oy, char *text);

extern int _Plot_draw_box (char * xopt, float xtick, int nxsub,
                            char * yopt, float ytick, int nysub);

extern int _Plot_set_linestyle (int style);
extern int _Plot_set_color (int color);
extern int _Plot_set_charsize (float size);
extern int _Plot_set_line_width (int width);
extern int _Plot_get_color (void);

extern int Plot_LS_label (float x, int twoS_plus_1, char L);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
