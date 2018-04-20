#ifndef ISIS_ERRORS_H
#define ISIS_ERRORS_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2018  Massachusetts Institute of Technology

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

/* $Id: errors.h,v 1.2 2004/02/09 11:14:19 houck Exp $ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <limits.h>

/*
 * errors.sl parses this file to automatically generate
 * a matching array of message strings.
 * 
 * Be careful how you edit this...
 */ 

enum
{
/*!+!  [don't remove this line] */
   I_INFO                  /* "" */
  ,I_WARNING               /* "Warning" */
  ,I_ERROR                 /* "Error" */
  ,I_RANGE_ERROR           /* "Range Error" */
  ,I_INTERNAL              /* "Internal Error" */     
  ,I_FAILED                /* "Failed" */
  ,I_MALLOC_FAILED         /* "Malloc Failed" */    
  ,I_INVALID               /* "Invalid" */
  ,I_NOT_FOUND             /* "Not Found" */
  ,I_FILE_NOT_FOUND        /* "File not found" */
  ,I_NOT_IMPLEMENTED       /* "Not Implemented" */
  ,I_NO_DATA               /* "No Data" */
  ,I_LOADING               /* "Loading" */
  ,I_SCANNING              /* "Scanning" */
  ,I_INITIALIZING          /* "Initializing" */
  ,I_READ_OK               /* "Read" */
  ,I_READ_FAILED           /* "Read failed" */
  ,I_READ_OPEN_FAILED      /* "Failed opening file for reading" */
  ,I_READ_KEY_FAILED       /* "Failed reading keyword(s) from file" */
  ,I_READ_COL_FAILED       /* "Failed reading column(s) from file" */
  ,I_KEY_NOT_FOUND         /* "Keyword not found" */     
  ,I_COL_NOT_FOUND         /* "Column not found" */
  ,I_HDU_NOT_FOUND         /* "Extension not found" */
  ,I_WRITE_OPEN_FAILED     /* "Failed opening file for writing" */
  ,I_WRITE_FAILED          /* "Failed writing" */
  ,I_FILE_NOKEYS           /* "No keywords in file" */
  ,I_KEY_NOTSET            /* "Keyword not set" */
  ,I_UNSUPPORTED_FORMAT    /* "Unsupported format" */
  ,I_INVALID_UNCERT        /* "Invalid uncertainties replaced" */
/*!-! [don't remove this line] */
};

/* Severity levels:
 * relative to default setting Isis_Verbose=0.
 * Print only messages with severity <= Isis_Verbose
 */
enum
{
   INFO  =  1,
   WARN  =  0,
   FAIL  = -1,
   INTR  = -2,
   FATAL = -INT_MAX
};

extern void isis_vmesg (int severity, unsigned int err, const char *file, int line, const char *fmt, ...);
extern void _isis_throw_exception (int err, const char *file, int line);
extern int _isis_trace_return (int status, const char *file, int line);
extern void _isis_trace (const char *file, int line);

#define isis_throw_exception(err)  _isis_throw_exception((err),__FILE__,__LINE__)
#define isis_trace_return(s)       _isis_trace_return((s),__FILE__,__LINE__)
#define isis_trace                 _isis_trace(__FILE__,__LINE__)

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
