/* -*- mode: C; mode: fold -*- */

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

/* $Id: errors.c,v 1.3 2004/02/09 11:14:19 houck Exp $ */

#include "config.h"
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>

#include <slang.h>

#include "isis.h"
#include "util.h"
#include "errors.h"

static const char *Msg[] = {
#include "errors.inc"
};

#define INVALID_ERR(err)   ((err) >= (sizeof (Msg) / sizeof (*Msg)))

int Isis_Error;
int Isis_Trace;

void isis_vmesg (int severity, unsigned int err, const char *file, int line,
                 const char *fmt, ...)
{
   char buf[4096];

   if (Isis_Verbose < severity)
     return;

   if (INVALID_ERR(err))
     {
        char msg[] = "\n### Internal error [isis_vmesg]:  unknown error code = %d ###\n";
        fprintf (stderr, msg, err);
        return;
     }

   buf[0] = 0;

   if (fmt != NULL)
     {
        va_list ap;
        char *b = buf;
        if (*Msg[err] != 0) 
          {
             *b++ = ':';
             *b++ = ' ';
          }        
        va_start (ap, fmt);
        SLvsnprintf (b, sizeof(buf)-1, (char *) fmt, ap);
        va_end (ap);
     }

   /* For severe errors, throw a slang exception */

   if (severity < FAIL)
     {
        unsigned int print_traceback = SLang_Traceback;
#ifndef SL_TB_FULL
# define SL_TB_FULL  0x1
#endif
        if (SLang_Version >= 20100) print_traceback &= SL_TB_FULL;
        if (print_traceback)
          {
             if (file != NULL)
               fprintf (stderr, "Error: %s:%d\n", file, line);
          }
        SLang_verror (Isis_Error, "%s%s", Msg[err], buf);
     }
   else
     {
        fprintf (stderr, "%s%s\n", Msg[err], buf);
     }
}

void _isis_throw_exception (int err, const char *file, int line)
{
   if ((SLang_Traceback != 0) && (file != NULL))
     {
        SLang_verror (Isis_Error, "Error :%s:%d", file, line);
     }
   else SLang_set_error (err);
}

int _isis_trace_return (int status, const char *file, int line)
{
   if ((Isis_Trace != 0) && (status != 0))
     {
        fprintf (stderr, "tracing: return %d from %s:%d\n", status, file ? file : "??", line);
     }
   return status;
}

void _isis_trace (const char *file, int line)
{
   if (Isis_Trace == 0)
     return;
   fprintf (stderr, "tracing: at %s:%d\n", file ? file : "??", line);
}
