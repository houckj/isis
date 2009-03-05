#ifndef ISIS_RMF_H
#define ISIS_RMF_H

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

/* $Id: rmf.h,v 1.8 2004/06/06 02:42:22 houck Exp $ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

extern Isis_Rmf_t *Rmf_init_rmf_list (void);
extern int Rmf_load_rmf (Isis_Rmf_t *head, int method, char *file);

extern Isis_Rmf_t *Rmf_find_rmf_index (Isis_Rmf_t *head, int rmf_index);
extern int Rmf_is_identity (Isis_Rmf_t *rmf);
extern int Rmf_includes_effective_area (Isis_Rmf_t *rmf);
extern int Rmf_index (Isis_Rmf_t *rmf);
extern int Rmf_id_list (Isis_Rmf_t *head, unsigned int **ids, unsigned int *num);

extern int Rmf_get_data_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *num);
extern int Rmf_get_arf_grid (Isis_Rmf_t *rmf, double **lo, double **hi, unsigned int *num);

extern int Rmf_adjust_data_grid (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int num);
extern int Rmf_rebin_rmf (Isis_Rmf_t *rmf, double *lo, double *hi, unsigned int num);

extern int Rmf_find_peaks (Isis_Rmf_t *rmf, double **h_P, int *num);

extern int Rmf_delete_rmf (Isis_Rmf_t *head, int rmf_index);
extern void Rmf_free_rmf (Isis_Rmf_t *rmf);
extern void Rmf_free_rmf_list (Isis_Rmf_t *head);

extern int Rmf_init_rmf (Isis_Rmf_t *rmf, Isis_Rmf_Grid_Type *arf, Isis_Rmf_Grid_Type *ebounds);
extern Isis_Rmf_t *Rmf_create_delta_rmf (Isis_Rmf_Grid_Type *arf, Isis_Rmf_Grid_Type *ebounds); 

extern int Rmf_apply_rmf (Isis_Rmf_t *rmf, double *x, int num_orig_data,
                          double *arf_src, int *arf_notice_list,
                          int num_arf_noticed);

extern char *Rmf_name (Isis_Rmf_t *rmf);

typedef struct
{
   int index;
   int method;
   int order;
   char *grating;
   char *instrument;
   char *arg_string;
}
Rmf_Info_Type;

extern void Rmf_free_info (Rmf_Info_Type *info);
extern int Rmf_get_info (Isis_Rmf_t *rmf, Rmf_Info_Type *info);
extern int Rmf_set_info (Isis_Rmf_t *rmf, Rmf_Info_Type *info);

/* RMF methods */
#define RMF_DELTA 1
#define RMF_FILE  2
#define RMF_USER  3

/* Specific methods */
extern int Rmf_load_delta (Isis_Rmf_t *rmf, char *file);
extern int Rmf_load_file (Isis_Rmf_t *rmf, char *file);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
