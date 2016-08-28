#ifndef ISIS_ARF_H
#define ISIS_ARF_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2016  Massachusetts Institute of Technology

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

/* $Id: arf.h,v 1.8 2004/02/09 11:14:16 houck Exp $ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct
{
   double exposure;
   union
     {
        double s;
        double *v;
     }
   fracexpo;
   int fracexpo_is_vector;       /* boolean */
   int nbins;                    /* number of bins */

   int order;
   int part;
   int srcid;
   char *object;       /* target name */
   char *grating;      /* HETG or LETG */
   char *instrument;   /* ACIS-S or HRC-S */
   char *file;         /* name of input ARF file */
}
Arf_Info_Type;

extern Isis_Arf_t *Arf_init_arf_list (void);
extern void Arf_free_arf_list (Isis_Arf_t *head);
extern void Arf_free_arf (Isis_Arf_t *arf);
extern Isis_Arf_t *Arf_new_arf (int nbins);
extern int Arf_delete_arf (Isis_Arf_t *head, int arf_index);
extern Isis_Arf_t *Arf_find_arf_index (Isis_Arf_t *head, int arf_index);
extern int Arf_id_list (Isis_Arf_t *head, unsigned int **ids, unsigned int *num);

extern int Arf_is_identity (Isis_Arf_t *a);
extern Isis_Arf_t *Arf_make_identity_arf (double *bin_lo, double * bin_hi, unsigned int nbins);
extern int Arf_read_arf (Isis_Arf_t *r_head, char * filename);
extern int Arf_arf_size (Isis_Arf_t *a);
extern int Arf_get_arf (Isis_Arf_t *a, double *arf, double *arf_err,
                         double *binlo, double *binhi);
extern int Arf_put_arf (Isis_Arf_t *a, double *arf, double *arf_err,
                         double *binlo, double *binhi);
extern int Arf_get_arf_exposure (Isis_Arf_t *a, double * exposure);
extern int Arf_set_arf_exposure (Isis_Arf_t *a, double exposure);

extern int Arf_set_arf_info (Isis_Arf_t *a, Arf_Info_Type *ai);
extern int Arf_get_arf_info (Isis_Arf_t *a, Arf_Info_Type *ai);
extern void Arf_free_arf_info (Arf_Info_Type *ai);
extern int Arf_define_arf (Isis_Arf_t *head, unsigned int nbins,
                           double *bin_lo, double *bin_hi,
                           double *arf, double *arf_err);

extern int Arf_append_arf (Isis_Arf_t *head, Isis_Arf_t *a);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
