#ifndef ISIS_UTIL_H
#define ISIS_UTIL_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2010  Massachusetts Institute of Technology

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

/* $Id: util.h,v 1.9 2004/02/09 11:14:25 houck Exp $ */

#ifndef PI
#define PI 3.14159265358979323846
#endif

/*
 * NIST 1998 CODATA recommended values of physical constants:
 */

#define PLANCK      ((double) 6.62606876e-27)       /* Planck's constant (erg s) */
#define BOLTZ       ((double) 1.3806503e-16)        /* Boltzmann's constant (erg/K) */
#define CLIGHT      ((double) 2.99792458e10)        /* speed of light (cm/s) */
#define AMU         ((double) 1.66053873e-24)       /* atomic mass unit (g) */
#define ERG_PER_EV  ((double) 1.602176462e-12)

#define ERG_ANGSTROM    (PLANCK * CLIGHT / 1.e-8)
#define KEV_ANGSTROM    (ERG_ANGSTROM / ERG_PER_EV / 1.e3)

/*
 *    Codes for supported physical units:
 */
enum
{
   U_ANGSTROM = 1, U_NANOMETER, U_MILLIMETER, U_CENTIMETER, U_METER,
   U_EV, U_KEV, U_MEV, U_GEV, U_TEV,
   U_HZ, U_KHZ, U_MHZ, U_GHZ
};

#define BUFSIZE 1024
#define COMMENT_CHAR  '#'

#define GRID_TOL  1.e-5           /* fractional error tolerance for wavelength grids */

#undef MAX
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

#undef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define INTERVALS_OVERLAP(alo,ahi,blo,bhi) (((alo) < (bhi)) && ((blo) < (ahi)))

#define bit_set(uint,mask)   ((uint) |= (mask))
#define bit_clear(uint,mask)   ((uint) &= ~(mask))

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

extern int Isis_Errno;
extern int Isis_Batch_Mode;
extern int Isis_Remove_Spectrum_Gaps;

#include <stdio.h>
extern int isis_fclose (FILE *fp);
extern FILE *isis_open_pager (void);
extern void isis_close_pager (FILE *fp);

extern int is_regular_file (char *file);
extern void isis_set_errno (int err);
extern int (*Isis_User_Break_Hook) (void);

extern double isis_nan (void);

#ifdef SLANG_VERSION
extern int isis_coerce_array_to_type (SLang_Array_Type **at, int type);
#endif

extern int Isis_pop_double_array (double *x, SLindex_Type n);
extern double isis_hypot (double x, double y);
extern double isis_kahan_sum (double *x, unsigned int n);
extern double isis_kahan_sum_squares (double *x, unsigned int n);
extern int isis_svd_solve (double **a, unsigned int n, double *b);
extern int isis_lu_solve (double **a, unsigned int n, unsigned int *piv, double *b);

extern char **new_string_array (int n, int len);
extern void free_string_array (char **p, int n);
extern int edit_temp_file (int (*save_file)(char *), int (*load_file)(char *), char *file);
extern int bsearch_d (double t, double *x, int n);
extern int find_bin (double x, double *lo, double *hi, int n);

extern int unit_id (char *name);
extern int unit_name (char *name, int unit);
extern int unit_invalid (int unit);
extern int unit_is_wavelength (int unit);
extern int unit_convert_x (double *x, int n, int out_units, int in_units);
extern int unit_convert_dfdx (double *dfdx, double x_in, int out_units, int in_units);

extern int reverse_d (double *x, unsigned int n);
extern int reverse_f (float *x, unsigned int n);
extern int reverse_i (int *x, unsigned int n);
extern int reverse_u (unsigned int *x, unsigned int n);

extern int rebin_histogram (double *fy, double *flo, double *fhi, int nf,
                            double *ty, double *tlo, double *thi, int nt);
extern int interpolate_dvector (double *xp, double *yp, unsigned int n,
                                double *oldxp, double *oldyp, unsigned int oldn);
extern double interpolate_d (double x, double *xp, double *yp, unsigned int n);

#define BINS_OK         0
#define BINS_REVERSED   1
#define BINS_INVALID   -1

extern int get_canonical_coordinates (double *bin_lo, double *bin_hi,
                                      int *nbins, int input_units);
extern int validate_wavelength_grid (int n, double *lo, double *hi);

extern int _isis_fixup_lo_hi_grids (double *lo, double *hi, unsigned int num,
                                    int *reversed, int *tweaked);

extern int transfer_notice (double *flo, double *fhi, int *f_notice_list, int f_num_notice, /*{{{*/
                            double *tlo, double *thi, int tn, int *t_notice);
extern int _update_notice_list (int *notice, int **notice_list, int *n_notice, int nbins);
extern int unpack_noticed (double *packed, int *notice_list, int n_notice,
                           int nbins, double *unpacked);

extern int apply_rebin (double *x, unsigned int num_orig_data,
                        int *rebin_flags, unsigned int num_rebin_data,
                        double *x_rebin);

extern int apply_rebin_weights (double *x, double *wt, unsigned int moment,
                                unsigned int num_orig_data,
                                int *rebin_flags, unsigned int num_rebin_data,
                                double *x_rebin);

extern int pop_qualifiers_arg (SLang_Struct_Type **sp);

typedef struct Isis_Arg_Type Isis_Arg_Type;
struct Isis_Arg_Type
{
   Isis_Arg_Type *next;
   void *arg;
   int type;
};

extern Isis_Arg_Type *isis_pop_args (int n);
extern int isis_push_args (Isis_Arg_Type *at);
extern void isis_free_args (Isis_Arg_Type *at);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif

