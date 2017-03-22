#ifndef ISIS_DBDISPLAY_H
#define ISIS_DBDISPLAY_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2017   Massachusetts Institute of Technology
 
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

/* $Id$ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* printing/plotting */

extern void DB_list_line_group_stats (FILE * fp, DB_t *db);
extern int DB_print_line_group (FILE * fp, DB_line_group_t *g, DB_t *db);
extern int DB_list_ion_levels (FILE * fp, int Z, int q, DB_t *db);

#ifdef ISIS_PLOT_H

typedef struct
{
   int *use_color;
   int label_type;                  /* label_type = < 0 for no labels, 0 for short labels, 1 for long labels */
   int label_length;                /* sizeof(label) */
   float bottom_frac, top_frac;     /* def: (0.6, 0.7)  indicator line endpoints as fraction of plot limits */
   float offset;                    /* separation between indicator line and text label */
                                    /* as fraction of indicator line length */
   float angle;                     /* degrees (90 ccw is default) */
   float justify;                   /* 0=left[default], 0.5=center, 1.0=right */
   float char_height;               /* character size (1.0 is default) */
}
Line_Label_Style_Type;

extern int DB_plot_line_group (DB_line_group_t *t, Line_Label_Style_Type *s, float redshift, 
                               char *(*label_massage_hook)(char *),
                               Plot_t *fmt, DB_t *db);
extern int DB_plot_line_list (Plot_t *fmt, Line_Label_Style_Type *s, float redshift,
                                float *lambda, char **labels, unsigned int nlines);
extern int DB_plot_levels (int Z, int q, int *line_idx, int nlines, int subset,
                           int overlay, Plot_t *fmt, DB_t *db);
extern int DB_plot_transitions (int Z, int q, int *line_idx, int nlines,
                                int *use_color, Plot_t *fmt, DB_t *db);
#endif /* ISIS_PLOT_H */

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
