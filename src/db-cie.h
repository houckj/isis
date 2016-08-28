#ifndef ISIS_DBCIE_H
#define ISIS_DBCIE_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2016  Massachusetts Institute of Technology
 
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

/* $Id: db-cie.h,v 1.2 2004/02/09 11:14:18 houck Exp $ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* this is here instead of in db-em.h because I used temp/density explicitly */
extern int EM_get_line_emissivity_function (float **emis, float **temps, float **densities,
                                            int *num_points, int line_index, EM_t *em);

typedef struct
{
   float trange[2];
   float drange[2];
}
EM_Range_Type;

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif 
