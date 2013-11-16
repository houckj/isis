/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2012  Massachusetts Institute of Technology

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

/* $Id: cfits.h,v 1.2 2004/02/09 11:14:17 houck Exp $ */

#ifndef CFITSIO_H
#define CFITSIO_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* these should match the values in fitsio.h */
enum
{
 CFLEN_FILENAME = 161,
 CFLEN_KEYWORD = 72,
 CFLEN_COMMENT = 73,
 CFLEN_VALUE = 71
};

extern int Cfits_Verbose;

typedef void *cfitsfile;

extern cfitsfile *cfits_open_file_readonly (const char *filename);
extern cfitsfile *cfits_open_file_readonly_silent (const char *filename);
extern int cfits_close_file (cfitsfile *fptr);

extern int cfits_get_hdu_num (cfitsfile *fptr, int *hdu_num);
extern int cfits_movabs_hdu (int hdu_num, cfitsfile *fptr);
extern int cfits_movrel_hdu (int nmove, cfitsfile *fptr);
extern int cfits_movnam_hdu (cfitsfile *fptr, const char *extname);
extern int cfits_move_to_matching_hdu (cfitsfile *fptr, const char *std_extnames[],
                                       const char *nonstd_names_hook, int (*check_hdu)(cfitsfile *));

extern int cfits_col_exist (const char *col_name, cfitsfile *fptr);

extern int cfits_read_int_keyword (int *keyvalue, const char *keyname,
                                   cfitsfile *fptr);
extern int cfits_read_long_keyword (long *keyvalue, const char *keyname,
                                    cfitsfile *fptr);
extern int cfits_read_float_keyword (float *keyvalue, const char *keyname,
                                     cfitsfile *fptr);
extern int cfits_read_double_keyword (double *keyvalue, const char *keyname,
                                      cfitsfile *fptr);
extern int cfits_read_string_keyword (char *keyvalue, const char *keyname,
                                      cfitsfile *fptr);
extern int cfits_get_colunits (char *col_units, const char *col_name,
                               cfitsfile *fptr);
extern int cfits_read_int_colkey (int *keyvalue, const char *keyroot,
                                  const char *colname, cfitsfile *fptr);

extern int cfits_test_keyword (const char *test_value, const char *keyword,
                                cfitsfile *fptr);

extern int cfits_get_repeat_count (int *nelems, const char *colname, cfitsfile *fptr);

extern int cfits_read_double_col (double *column_values, long num_values,
                                  long firstrow, const char *colname, cfitsfile *fptr);
extern int cfits_read_float_col (float *column_values, long num_values,
                                 long firstrow, const char *colname, cfitsfile *fptr);
extern int cfits_read_int_col (int *column_values, long num_values,
                               long firstrow, const char *colname, cfitsfile *fptr);
extern int cfits_read_string_col (char **column_values, long num_values,
                                  long firstrow, const char *colname, cfitsfile *fptr);

extern long cfits_optimal_numrows (cfitsfile *fptr);

extern int cfits_locate_vextension (cfitsfile *ft, int argc, const char **ext_names,
                                    int (*fun) (cfitsfile *));
extern int cfits_get_column_numbers (cfitsfile *ft, unsigned int num, const char **names, int *cols);

extern int cfits_read_column_uints (cfitsfile *ft, int col, int row, long ofs,
                                     unsigned int *data, int nrows);
extern int cfits_read_column_floats (cfitsfile *ft, int col, long row, long ofs,
                                     float *data, int nrows);

extern int cfits_read_optional_double_col (double *dat, int nbins,
                                           int k, const char *name, cfitsfile *cfp);
extern int cfits_read_optional_int_col (int *dat, int nbins, int k, const char *name,
                                        cfitsfile *cfp);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
