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

/* $Id: cfits.c,v 1.3 2004/02/09 11:14:17 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

/* workaround: fitsio.h needs off_t but doesn't include sys/types.h */
#ifdef HAVE_SYS_TYPES_H    
#include <sys/types.h>
#endif

#include <fitsio.h>

#include "isis.h"
#include "util.h"
#include "cfits.h"
#include "errors.h"

/*}}}*/

/*{{{ error reporting */

int Cfits_Verbose;

static void cfits_report_error (int status)
{
   if (Cfits_Verbose)
     {
        fits_report_error (stderr, status);
        fflush (stderr);
     }
}

/*}}}*/

/*{{{ open, close, change extensions */

static cfitsfile *cfits_open_file (char *filename, int iomode)
{
   int status = 0;
   fitsfile *fptr = NULL;
   
   if (NULL == filename)
     return NULL;
   
   (void) fits_open_file (&fptr, filename, iomode, &status);
   cfits_report_error (status);
   if (status != 0) 
     return NULL;
   
   return (cfitsfile *) fptr;
}

cfitsfile *cfits_open_file_readonly_silent (char *filename)
{
   int status = 0;
   fitsfile *fptr = NULL;
   
   if (NULL == filename)
     return NULL;
   
   (void) fits_open_file (&fptr, filename, READONLY, &status);
   if (status != 0) 
     return NULL;
   
   return (cfitsfile *) fptr;
}

cfitsfile *cfits_open_file_readonly (char *filename)
{
   return cfits_open_file (filename, READONLY);
}

int cfits_close_file (cfitsfile *fptr)
{
  int status = 0;
   
   if (NULL == fptr)
     return 0;

   (void) fits_close_file((fitsfile *) fptr, &status);
   cfits_report_error (status);
   
  if (status != 0) return -1;

  return 0;
}

int cfits_movabs_hdu (int hdunum, cfitsfile *fptr)
{
   int status = 0;
   int *hdutype = NULL;
   
   (void) fits_movabs_hdu((fitsfile *) fptr, hdunum, hdutype, &status);
   cfits_report_error (status);

   if (status != 0) return -1;
   return 0;
}

int cfits_movrel_hdu (int nmove, cfitsfile *fptr)
{
   int status = 0;
   int *hdutype = NULL;
   
   (void) fits_movrel_hdu((fitsfile *) fptr, nmove, hdutype, &status);
   cfits_report_error (status);

   if (status != 0) return -1;
   return 0;
}

int cfits_movnam_hdu (cfitsfile *fptr, char *extname)
{
   int status = 0;

   (void) fits_movnam_hdu((fitsfile *) fptr, ANY_HDU, extname, 0, &status);
   cfits_report_error (status);

   if (status != 0) return -1;
   return 0;
}

/*}}}*/

/*{{{ read and test keywords */

int cfits_read_int_keyword (int *keyvalue, char *keyname, cfitsfile *fptr)
{
   int status = 0;
   int datatype = TINT;
   char comment[CFLEN_COMMENT];
   
   (void) fits_read_key ((fitsfile *) fptr, datatype, keyname, keyvalue,
                        comment, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_long_keyword (long *keyvalue, char *keyname, 
                             cfitsfile *fptr)
{
   int status = 0;
   int datatype = TLONG;
   char comment[CFLEN_COMMENT];
   
   (void) fits_read_key ((fitsfile *) fptr, datatype, keyname, keyvalue,
                        comment, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_double_keyword (double *keyvalue, char *keyname, 
                                cfitsfile *fptr)
{
   int status = 0;
   int datatype = TDOUBLE;
   char comment[CFLEN_COMMENT];
   
   (void) fits_read_key ((fitsfile *) fptr, datatype, keyname, keyvalue,
                        comment, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_float_keyword (float *keyvalue, char *keyname, 
                                cfitsfile *fptr)
{
   int status = 0;
   int datatype = TFLOAT;
   char comment[CFLEN_COMMENT];
   
   (void) fits_read_key ((fitsfile *) fptr, datatype, keyname, keyvalue,
                        comment, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_string_keyword (char *keyvalue, char *keyname, 
                                cfitsfile *fptr)
{
   int status = 0;
   int datatype = TSTRING;
   char comment[CFLEN_COMMENT];
   
   (void) fits_read_key ((fitsfile *) fptr, datatype, keyname, keyvalue,
                        comment, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

/*returns -1 on error or if value in file is different from 
 *           the string "test_value"
 *return   1 if the keyword is missing
 *returns  0 if the keyword value matches "test_value"
 */ 

int cfits_test_keyword (char *test_value, char *keyword,
                        cfitsfile *fptr)
{
   char value[CFLEN_VALUE];
   
   if (NULL == test_value || NULL == keyword || NULL == fptr)
     return -1;
   
   if (-1 == cfits_read_string_keyword(value, keyword, fptr))
     return -1;

   if (NULL == value)
     {
        isis_vmesg (WARN, I_KEY_NOT_FOUND, __FILE__, __LINE__, "%s", keyword);
        return 1;
     }
   
   return isis_strcasecmp (value, test_value);
}

/*}}}*/

/*{{{ get column info */

static int cfits_get_colnum (int *colnum, char *col_name, cfitsfile *fptr)
{
   int status = 0;
   int casesen = CASEINSEN;
   
   (void) fits_get_colnum ((fitsfile *) fptr, casesen, col_name, colnum,
                          &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_get_colunits (char *col_units, char *col_name, cfitsfile *fptr)
{
   int status = 0;
   int colnum;
   
   if (NULL == col_units || NULL == col_name || NULL == fptr)
     return -1;
   
   if (-1 == cfits_get_colnum (&colnum, col_name, fptr))
     return -1;
   
   (void) fits_get_bcolparms ((fitsfile *) fptr, colnum, NULL, col_units,
                              NULL, NULL, NULL, NULL, NULL, NULL, &status);
   
   cfits_report_error (status);
   if (status != 0) return -1;
   
   return 0;
}

int cfits_col_exist (char *col_name, cfitsfile *fptr)
{
   int status = 0;
   int casesen = CASEINSEN;
   int colnum;
   
   if (NULL == col_name || NULL == fptr)
     return -1;
   
   (void) fits_get_colnum ((fitsfile *) fptr, casesen, col_name, &colnum,
                           &status);
   if (status == 0) 
     return 1;                /*yes, column exists */
   else
     {
        fits_clear_errmsg ();
        return 0;                /*no, column doesn't exist */
     }
}

int cfits_get_repeat_count (int *nelems, char *colname, cfitsfile *fptr)
{
   int status = 0, colnum;
   long ne;

   if (-1 == cfits_get_colnum (&colnum, colname, fptr))
     return -1;
   
   (void) fits_get_coltype((fitsfile *) fptr, colnum, NULL, &ne, NULL, &status);
   cfits_report_error (status);
   
   *nelems = ne;
   
   if (status != 0) return -1;
   return 0;
}

/*}}}*/

/*{{{ read bintable columns */

long cfits_optimal_numrows (cfitsfile *fptr)
{
   int status = 0;
   long nrows;
   
   (void) fits_get_rowsize ((fitsfile *) fptr, &nrows, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return nrows;
}

int cfits_read_int_colkey (int *keyvalue, char *keyroot,
                           char *colname, cfitsfile *fptr)
{
   char keyname[CFLEN_KEYWORD];
   char comment[CFLEN_COMMENT];
   int datatype = TINT;

   int status = 0;
   int colnum;

   if (NULL == fptr || NULL == colname || NULL == keyroot
       || NULL == keyvalue)
     return -1;

   if (-1 == cfits_get_colnum (&colnum, colname, fptr))
     return -1;

   if ((strlen(keyroot) + 3) > CFLEN_KEYWORD)
      return -1;

   (void) sprintf (keyname, "%s%i", keyroot, colnum);
   
   (void) fits_read_key ((fitsfile *) fptr, datatype, keyname, keyvalue,
                        comment, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_double_col (double *column_values, long num_values,
                           long firstrow, char *colname, cfitsfile *fptr)
{
   int datatype = TDOUBLE;
   double nulval = DBL_MIN;

   int status = 0;
   long firstelem = 1;
   int anynul, colnum;
   
   if (-1 == cfits_get_colnum (&colnum, colname, fptr))
     {
       column_values = NULL;
       return -1;
     }

   (void) fits_read_col ((fitsfile *) fptr, datatype, colnum, firstrow,
                        firstelem, num_values, (void *) &nulval,
                        (void *) column_values, &anynul, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_float_col (float *column_values, long num_values,
                         long firstrow, char *colname, cfitsfile *fptr)
{
   int datatype = TFLOAT;
   float nulval = FLT_MIN;

   int status = 0;
   long firstelem = 1;
   int anynul, colnum;
   
   if (-1 == cfits_get_colnum (&colnum, colname, fptr))
     {
       column_values = NULL;
       return -1;
     }

   (void) fits_read_col ((fitsfile *) fptr, datatype, colnum, firstrow,
                        firstelem, num_values, (void *) &nulval,
                        (void *) column_values, &anynul, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_int_col (int *column_values, long num_values,
                       long firstrow, char *colname, cfitsfile *fptr)
{
   int datatype = TINT;
   int nulval = 0;
   
   int status = 0;
   long firstelem = 1;
   int anynul, colnum;
   
   if (-1 == cfits_get_colnum (&colnum, colname, fptr))
     {
       column_values = NULL;
       return -1;
     }

   (void) fits_read_col ((fitsfile *) fptr, datatype, colnum, firstrow,
                        firstelem, num_values, (void *) &nulval,
                        (void *) column_values, &anynul, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_string_col (char **array, long num_values,
                           long firstrow, char *colname, cfitsfile *fptr)
{
   char nulval[CFLEN_VALUE];

   int status = 0;
   long firstelem = 1;
   int anynul, colnum;
   
   nulval[0] = '\0';

   if (-1 == cfits_get_colnum (&colnum, colname, fptr))
     {
       array = NULL;
       return -1;
     }

   (void) fits_read_col_str ((fitsfile *) fptr, colnum, firstrow, firstelem,
                             num_values, nulval, array, &anynul, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

/*}}}*/

int cfits_locate_vextension (cfitsfile *ft, int argc, char **ext_names, /*{{{*/
                             int (*fun) (cfitsfile *))
{
   char buf[CFLEN_VALUE];
   
   if (ft == NULL)
     return -1;
   
   if (-1 == cfits_movabs_hdu (2, ft))
     return -1;
   
   do
     {
        if (0 == cfits_read_string_keyword (buf, "EXTNAME", ft))
          {
             int i;
             
             for (i = 0; i < argc; i++)
               {
                  if (!strcmp (buf, ext_names[i]))
                    {
                       if ((fun == NULL) || (-1 != (*fun) (ft)))
                         return 0;
                    }
               }
          }
     }
   while (-1 != cfits_movrel_hdu (1, ft));
   
   return -1;
}

int cfits_read_column_uints (cfitsfile *ft, int col, int row, long ofs,
                             unsigned int *data, int nrows)
{
   int anynul, status = 0;
   
   if (ft == NULL)
     return -1;

   (void) ofs;

   if (nrows == 0) return 0;
   
   (void) ffgcvuk ((fitsfile *)ft, col, row, 1, nrows, 0, data, &anynul, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_read_column_floats (cfitsfile *ft, int col, long row, long ofs, 
                              float *data, int nrows)
{
   int anynul, status = 0;
   float nulval = FLT_MIN;
   
   if (ft == NULL)
     return -1;
   
   if (nrows == 0) return 0;
   
   (void) ffgcve ((fitsfile *)ft, col, row, ofs, nrows, nulval, data, &anynul, &status);
   cfits_report_error (status);
   
   if (status != 0) return -1;
   return 0;
}

int cfits_get_column_numbers (cfitsfile *ft, unsigned int num, char **names, int *cols)
{
   unsigned int i;
   int status = 0;

   if (NULL == ft)
     return -1;
   
   for (i = 0; i < num; i++)
     {
        if (fits_get_colnum ((fitsfile *)ft, 1, names[i], cols + i, &status))
          {
             status = 0;
             cols[i] = -1;
          }
     }
   return 0;
}


int cfits_read_optional_double_col (double *dat, int nbins, int k, char *name,
                                    cfitsfile *cfp)
{
   if (!cfits_col_exist (name, cfp))
     {
        memset ((char *)dat, 0, nbins *sizeof(double));
        return 0;
     }
   
   return cfits_read_double_col (dat, nbins, k, name, cfp);
}

int cfits_read_optional_int_col (int *dat, int nbins, int k, char *name,
                                 cfitsfile *cfp)
{
   if (!cfits_col_exist (name, cfp))
     {
        memset ((char *)dat, 0, nbins *sizeof(int));
        return 0;
     }
   
   return cfits_read_int_col (dat, nbins, k, name, cfp);
}


