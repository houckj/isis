#ifndef ISIS_KEYWORD_H
#define ISIS_KEYWORD_H

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

/* $Id: keyword.h,v 1.4 2004/02/09 11:14:22 houck Exp $ */

#include <stddef.h>

#define KEYWORD_CHAR  ';'        /* in column 1 to mark keyword lines */
#define KEYLINE_SIZE  1024

#define FITS_FILE  0
#define ASCII_FILE 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct _Keyword_t Keyword_t;

struct _Keyword_t
{
   char *name;
   int required;                /* complain if missing */
#define OPTIONAL_KEY  0
#define REQUIRED_KEY  1
   int key_type;
#define   NULL_KEY     0
#define   STRING_KEY   1        /* of size CFLEN_VALUE */
#define   DOUBLE_KEY   2   
#define   FLOAT_KEY    3
#define   INT_KEY      4
   size_t offset;               /* offsetof (struct, field) */
};

typedef struct _Ascii_Keybuf_t Ascii_Keybuf_t;

extern int Key_read_header_keywords (void *fp, char *ptr, Keyword_t *k, int method);

extern Ascii_Keybuf_t * Key_ascii_allocate_keybuf (unsigned int num_keys);
extern void Key_ascii_free_keybuf (Ascii_Keybuf_t *buf);
extern int Key_ascii_load_keybuf (FILE *fp, Ascii_Keybuf_t *buf);
extern int Key_ascii_read_string_keyword (char *key, char *name, Ascii_Keybuf_t *buf);
extern int Key_ascii_read_double_keyword (double *key, char *name, Ascii_Keybuf_t *buf);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif


#endif

