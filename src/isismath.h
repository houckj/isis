
/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2009 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.
    
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

/* $Id: isismath.h,v 1.3 2004/02/09 11:14:22 houck Exp $ */

#ifndef ISISMATH_H
#define ISISMATH_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

extern double *JDMdouble_vector (unsigned int n);
extern void JDMfree_double_matrix (double **matrix, unsigned int n);
extern double **JDMdouble_matrix (unsigned int n, unsigned int m);

extern void _JDMswap_dvector (double *a, double *b, unsigned int n);
extern double _JDM_innerprod (double *a, double *b, unsigned int n);
extern double _JDM_innerprod_col (double *a, double **b, unsigned int col, unsigned int n);
extern double *_JDM_equilibrate (double **a, unsigned int n, double *eq);

extern int JDM_lu_decomp (double **a, unsigned int n, unsigned int *pivot, double epsilon, int *parityp);
extern int JDM_lu_backsubst (double **a, unsigned int n, unsigned int *pivot, double *b);
extern int JDM_ludecomp_inverse (double **a, unsigned int n);

/* Fast Fourier Transforms */

extern int JDMfftn (int ndim, int *dims, double * Re, double * Im,
		    int iSign, double scaling);
extern int JDMfftnf (int ndim, int *dims, float *Re, float *Im,
		     int iSign, double scaling);
extern void JDMfft_free (void);

/* random numbers */

extern void random_seed (unsigned long seed);
extern double urand (void);
extern double grand (void);
unsigned int prand (double rate);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
