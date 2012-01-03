/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2012  Massachusetts Institute of Technology

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

/* $Id: keyword.c,v 1.5 2004/02/09 11:14:22 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include "isis.h"
#include "util.h"
#include "cfits.h"
#include "keyword.h"
#include "errors.h"

/*}}}*/

/* methods for ASCII file keyword input */

struct _Ascii_Keybuf_t
{
   char **keyword_line;   /* buffer for ascii file keyword lines */
   unsigned int num_keys;
};

void Key_ascii_free_keybuf (Ascii_Keybuf_t *p) /*{{{*/
{
   unsigned int i;

   if (p == NULL)
     return;

   for (i = 0; i < p->num_keys; i++)
     ISIS_FREE (p->keyword_line[i]);

   ISIS_FREE (p->keyword_line);
   ISIS_FREE (p);
}

/*}}}*/

Ascii_Keybuf_t * Key_ascii_allocate_keybuf (unsigned int num_keys) /*{{{*/
{
   Ascii_Keybuf_t *p;
   unsigned int i;

   if (NULL == (p = (Ascii_Keybuf_t *) ISIS_MALLOC (sizeof(Ascii_Keybuf_t)))
       || NULL == (p->keyword_line = (char **) ISIS_MALLOC (num_keys * sizeof(char *))))
     {
        Key_ascii_free_keybuf (p);
        return NULL;
     }

   p->num_keys = num_keys;

   for (i=0; i < num_keys; i++)
     {
        if (NULL == (p->keyword_line[i] = (char *) ISIS_MALLOC (KEYLINE_SIZE * sizeof(char))))
          {
             Key_ascii_free_keybuf (p);
             return NULL;
          }
     }

   return p;
}

/*}}}*/

int Key_ascii_load_keybuf (FILE *fp, Ascii_Keybuf_t *p) /*{{{*/
{
   char buf[KEYLINE_SIZE];
   unsigned int n;

   n = 0;
   while (!feof(fp) && n < p->num_keys)
     {
        char *k;

        if (NULL == fgets (buf, sizeof(buf), fp))
          break;
        if (buf[0] == COMMENT_CHAR)
          continue;
        if (buf[0] != KEYWORD_CHAR)   /* all keywords must be at top of file */
          break;

        /* skip leading whitespace */

        for (k = buf + 1; *k == ' ' || *k == '\t'; k++)
          ;
        if (*k == '\n' || *k == '\0')
          break;                        /* don't allow blank keyword line */

        isis_strcpy (p->keyword_line[n++], k, KEYLINE_SIZE);
     }

   return ((n < p->num_keys) ? -1 : 0);
}

/*}}}*/

static char *ascii_get_keyword_value (const char *keyname, Ascii_Keybuf_t *p) /*{{{*/
{
   unsigned int i;

   for (i=0; i < p->num_keys; i++)
     {
        char *k = p->keyword_line[i];     /* with no leading whitespace */
        const char *n = keyname;

        while (tolower((unsigned char)*k) == tolower((unsigned char)*n))
          {
             k++; n++;
          }

        if (*n != '\0')
          continue;

        while (!(*k == ' ' || *k == '\t' || *k == '\n' || *k == '\0'))
          k++;

        if (*k == '\n' || *k == '\0')
          return NULL;            /* no whitespace trailing keyword */
        else
          return k;               /* ptr to value string */
     }

   return NULL;
}

/*}}}*/

static int ascii_read_int_keyword (int *key, char *name, Ascii_Keybuf_t *buf) /*{{{*/
{
   char *value;

   if (NULL == (value = ascii_get_keyword_value (name, buf)))
     return -1;

   if (1 != sscanf (value, "%d", key))
     return -1;

   return 0;
}

/*}}}*/

static int ascii_read_float_keyword (float *key, char *name, Ascii_Keybuf_t *buf) /*{{{*/
{
   char *value;

   if (NULL == (value = ascii_get_keyword_value (name, buf)))
     return -1;

   if (1 != sscanf (value, "%e", key))
     return -1;

   return 0;
}

/*}}}*/

int Key_ascii_read_double_keyword (double *key, const char *name, Ascii_Keybuf_t *buf) /*{{{*/
{
   char *value;

   if (NULL == (value = ascii_get_keyword_value (name, buf)))
     return -1;

   if (1 != sscanf (value, "%le", key))
     return -1;

   return 0;
}

/*}}}*/

int Key_ascii_read_string_keyword (char *key, const char *name, Ascii_Keybuf_t *buf) /*{{{*/
{
   char *value, *h, *t, *s;
   int len;

   *key = '\0';

   if (NULL == (value = ascii_get_keyword_value (name, buf)))
     return -1;

   /* Unfortunately, we can't use sscanf() because we
    * want to allow embedded spaces in the value string, so...
    */

   /* skip whitespace */

   for (h = value; *h == ' ' || *h == '\t'; h++)
     ;
   if (*h == '\0' || *h == '\n')
     return 0;

   /* find last non-whitespace char */

   t = h;
   for (s = h; *s != '\0' && *s != '\n'; s++)
     {
        if (*s != ' ' && *s != '\t')
          t = s;
     }

   /* copy trimmed string */

   len = MIN(t-h+1, (int) CFLEN_VALUE);
   isis_strcpy (key, h, len);

   return 0;
}

/*}}}*/

/* keyword input functions */

/*{{{ Various read methods:
 * e.g. ascii:          int ascii_read_int_keyword (&int, "name",  fp);
 * e.g. fits:           int cfits_read_int_keyword (&int, "name", cfp);
 *
 * To add a new method,
 *   1) write functions to read each keyword type
 *   2) add a new preprocessor constant for the file type e.g.:
 *               #define MY_FILE 2
 *   3) modify set_keyword_method () so the global method pointers
 *      (Int_Method, etc.) are set correctly.
 */

typedef int method_t(void *, const char *, void *);

typedef struct
{
   method_t *read_double;
   method_t *read_float;
   method_t *read_int;
   method_t *read_string;
   int method;
}
Methods_t;

static Methods_t Read_Methods[] =
{
     {(method_t *) cfits_read_double_keyword,
      (method_t *) cfits_read_float_keyword,
      (method_t *) cfits_read_int_keyword,
      (method_t *) cfits_read_string_keyword,
      FITS_FILE
     },
     {(method_t *) Key_ascii_read_double_keyword,
      (method_t *) ascii_read_float_keyword,
      (method_t *) ascii_read_int_keyword,
      (method_t *) Key_ascii_read_string_keyword,
      ASCII_FILE
     },
     {NULL,NULL,NULL,NULL,-1}
};

/*}}}*/

static int set_keyword_method (Methods_t *m, int method) /*{{{*/
{
   Methods_t *p = Read_Methods;

   while (p->read_double != NULL)
     {
        if (p->method == method)
          {
             *m = *p;  /* struct copy */
             return 0;
          }
        p++;
     }

   isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "keyword input method");
   memset ((char *)m, 0, sizeof(*m));

   return -1;
}

/*}}}*/

int Key_read_header_keywords (void *fp, char *ptr, Keyword_t *k, int method) /*{{{*/
{
   Methods_t m;

   if (NULL == k || NULL == fp || NULL == ptr)
     return -1;

   if (-1 == set_keyword_method (&m, method))
     return -1;

   for ( ; k->name != NULL; k++)
     {
        method_t *read_key;

        switch (k->key_type)
          {
           case INT_KEY:
             read_key = m.read_int;
             break;

           case FLOAT_KEY:
             read_key = m.read_float;
             break;

           case DOUBLE_KEY:
             read_key = m.read_double;
             break;

           case STRING_KEY:
             read_key = m.read_string;
             break;

           default:
             isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "keyword '%s' of type '%s'",
                         k->name, k->key_type);
             return -1;
          }

        if (-1 == (*read_key) (ptr + k->offset, k->name, fp))
          {
             if (k->required == REQUIRED_KEY)
               isis_vmesg (WARN, I_KEY_NOTSET, __FILE__, __LINE__, "%s", k->name);
          }
     }

   return 0;
}

/*}}}*/

