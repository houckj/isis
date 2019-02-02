/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2019  Massachusetts Institute of Technology

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

/* $Id: miscio.c,v 1.5 2004/09/09 11:31:57 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>

#include "isis.h"
#include "util.h"
#include "_isis.h"
#include "cfits.h"
#include "errors.h"

/*}}}*/

#define ISIS_DEFAULT_PAGER  "${PAGER:-more}"

void isis_set_pager (char *pager) /*{{{*/
{
   char *s;

   if ((pager == NULL) || (*pager == 0))
     pager = ISIS_DEFAULT_PAGER;

   if (NULL == (s = isis_make_string (pager)))
     return;

   ISIS_FREE (Isis_Pager);
   Isis_Pager = s;
}

/*}}}*/

static void isis_get_pager (void) /*{{{*/
{
   (void) SLang_push_string (Isis_Pager);
}

/*}}}*/

static void find_file_in_path (char *path, char *name) /*{{{*/
{
   char *pathname;

   pathname = SLpath_find_file_in_path (path, name);
   (void) SLang_push_string (pathname);
   SLfree(pathname);
}

/*}}}*/

static void set_errno (int *errnum) /*{{{*/
{
   errno = *errnum;
   SLerrno_set_errno (*errnum);
}

/*}}}*/

/* ascii readcol */

typedef struct
{
   char *buf;
   unsigned int bufsize;
   unsigned int len;
}
Buffer_Type;

static void free_buf (Buffer_Type *b) /*{{{*/
{
   ISIS_FREE (b->buf);
   b->len = 0;
   b->bufsize = 0;
}

/*}}}*/

static int init_buf (Buffer_Type *b) /*{{{*/
{
   enum {DEFAULT_BUFSIZE = 8192};

   if (b->bufsize)
     return 0;

   b->buf = (char *) ISIS_MALLOC (DEFAULT_BUFSIZE * sizeof(char));
   if (b->buf == NULL)
     return -1;

   b->bufsize = DEFAULT_BUFSIZE;
   b->len = 0;

   return 0;
}

/*}}}*/

static int skip_to_eol (FILE *fp) /*{{{*/
{
   int ch;

   while (EOF != (ch = getc (fp)))
     {
        if (ch == '\n')
          return 0;
     }

   return -1;
}

/*}}}*/

static int readline (FILE *fp, Buffer_Type *b) /*{{{*/
{
   if (-1 == init_buf (b))
     return -1;

   b->len = 0;

   for (;;)
     {
        int ch = getc (fp);

        switch (ch)
          {
           case COMMENT_CHAR:
             b->buf[b->len] = 0;
             if (-1 == skip_to_eol (fp))
               return -1;
             return b->len ? 0 : COMMENT_CHAR;

           case '\n':
             b->buf[b->len] = 0;
             return 0;

           case EOF:
             b->buf[b->len] = 0;
             return -1;

           default:
             b->buf[b->len++] = ch;
             break;
          }

        if (b->len == b->bufsize)
          {
             unsigned int new_size = b->bufsize * 2 * sizeof(char);
             char *tmp;
             if (NULL == (tmp = (char *) ISIS_REALLOC (b->buf, new_size)))
               return -1;
             b->buf = tmp;
             b->bufsize *= 2;
          }
     }
}

/*}}}*/

/* cols[] must be sorted in increasing order and with no values repeated */
static int ascii_readcol (FILE *fp, unsigned int *cols, unsigned int ncols, /*{{{*/
                          double **px, unsigned int *pnx)
{
   enum {START_ROWS = 1024};
   const char *sepchars = " ,\t";
   Buffer_Type b = {NULL, 0, 0};
   unsigned int *last_col = cols + ncols;
   unsigned int k, line, nx = START_ROWS * ncols;
   double *x = NULL;
   int ret = -1;

   *px = NULL;
   *pnx = 0;

   if (ncols == 0)
     return -1;

   if (NULL == (x = (double *) ISIS_MALLOC (nx * sizeof(double))))
     return -1;

   k = 0;
   line = 0;

   do
     {
        unsigned int *want_col = NULL;
        unsigned int c, k_save;
        char *word = NULL;
        int ch;

        k_save = k;

        if (-1 == (ch = readline (fp, &b)))
          break;

        line++;

        if (ch == COMMENT_CHAR)
          continue;

        word = strtok (b.buf, sepchars);
        want_col = cols;
        c = 1;

        while ((word != NULL) && (want_col < last_col))
          {
             if (c == *want_col)
               {
                  double d;

                  if (1 != sscanf (word, "%le", &d))
                    {
                       char *p = strchr (word, '\n');
                       if (p) *p = 0;
                       isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "invalid input: '%s'", word);
                       goto finish;
                    }

                  if (k >= nx)
                    {
                       unsigned int new_size = 2 * nx * sizeof(double);
                       double *tmp;
                       if (NULL == (tmp = (double *) ISIS_REALLOC (x, new_size)))
                         goto finish;
                       x = tmp;
                       nx *= 2;
                    }

                  x[k++] = d;

                  want_col++;
               }

             word = strtok (NULL, sepchars);
             c++;
          }

        /* support skipping blank lines */
        if ((want_col < last_col) && (c > 1))
          {
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "at line %u: couldn't read all requested columns", line);
             k = k_save;
             break;
          }

     } while (!feof(fp));

   ret = 0;
   finish:

   free_buf (&b);

   if (ret)
     ISIS_FREE (x);
   else
     {
        *px = x;
        *pnx = k;
     }

   return ret;
}

/*}}}*/

static int uint_compar (const void *pa, const void *pb) /*{{{*/
{
   const unsigned int *a = (const unsigned int *) pa;
   const unsigned int *b = (const unsigned int *) pb;

   if (*a < *b) return -1;
   else if (*a > *b) return 1;
   return 0;
}

/*}}}*/

static int uniq_cols (unsigned int *cols, unsigned int *ncols) /*{{{*/
{
   unsigned int i, n;

   qsort (cols, *ncols, sizeof(unsigned int), &uint_compar);

   n = 1;

   for (i = 1; i < *ncols; i++)
     {
        if (cols[i] != cols[n-1])
          cols[n++] = cols[i];
     }

   *ncols = n;

   return 0;
}

/*}}}*/

static int push_cols (double *d, unsigned int n, unsigned int ncols) /*{{{*/
{
   SLindex_Type nrows;
   unsigned int c;

   if ((ncols == 0) || (d == NULL))
     return -1;

   nrows = n / ncols;
   for (c = 0; c < ncols; c++)
     {
        SLang_Array_Type *at;
        unsigned int k, i;
        double *x;

        at = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &nrows, 1);
        if (at == NULL)
          return -1;

        x = (double *) at->data;

        i = 0;
        for (k = c; k < n; k += ncols)
          {
             x[i++] = d[k];
          }

        SLang_push_array (at, 1);
     }

   return 0;
}

/*}}}*/

static void _readcol (void) /*{{{*/
{
   FILE *fp = NULL;
   char *file = NULL;
   SLang_Array_Type *sl_cols = NULL;
   unsigned int *cols = NULL;
   unsigned int ncols, nx;
   double *x = NULL;
   int ret = -1;

   if ((-1 == SLpop_string (&file))
       || (-1 == SLang_pop_array_of_type (&sl_cols, SLANG_UINT_TYPE))
       || (sl_cols == NULL))
     goto finish;

   cols = (unsigned int *) sl_cols->data;
   ncols = sl_cols->num_elements;

   fp = fopen (file, "r");
   if (fp == NULL)
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", file);
        goto finish;
     }

   if ((-1 == uniq_cols (cols, &ncols))
       || (-1 == ascii_readcol (fp, cols, ncols, &x, &nx))
       || (-1 == push_cols (x, nx, ncols)))
     goto finish;

   ret = 0;
   finish:

   if (ret) isis_throw_exception (Isis_Error);
   if (fp) fclose (fp);
   ISIS_FREE (x);
   ISIS_FREE (file);
   SLang_free_array (sl_cols);
}

/*}}}*/

/* printing arrays [code adapted from jdl functions by John Davis] */
static int
pop_matrix (SLang_Array_Type **at_ptr, unsigned int *nr, unsigned int *nc) /*{{{*/
{
   SLang_Array_Type *at;

   if (-1 == SLang_pop_array (&at, 0))
     return -1;

   switch (at->num_dims)
     {
      case 0:
        *nr = *nc = 0;
        break;
      case 1:
        *nr = (unsigned int)at->dims[0];
        *nc = 1;
        break;
      case 2:
        *nr = (unsigned int)at->dims[0];
        *nc = (unsigned int)at->dims[1];
        break;

      default:
        SLang_verror (SL_TYPE_MISMATCH, "operation limited to 2-d arrays");
        SLang_free_array (at);
        *at_ptr = NULL;
        return -1;
     }
   *at_ptr = at;
   return 0;
}

/*}}}*/

static void print_array (void) /*{{{*/
{
   enum
     {
        screen_rows=24,
        screen_cols=10
     };

   SLang_Array_Type *at;
   unsigned int i, num_rows, num_cols;
   unsigned char type;
   VOID_STAR v;
   FILE *fp;
   int just_one_line;
   unsigned int num;

   if (-1 == pop_matrix (&at, &num_rows, &num_cols))
     return;

   type = at->data_type;
   switch (type)
     {
      case SLANG_CHAR_TYPE:
      case SLANG_UCHAR_TYPE:
      case SLANG_SHORT_TYPE:
      case SLANG_USHORT_TYPE:
      case SLANG_INT_TYPE:
      case SLANG_UINT_TYPE:
      case SLANG_LONG_TYPE:
      case SLANG_ULONG_TYPE:
      case SLANG_DOUBLE_TYPE:
      case SLANG_FLOAT_TYPE:
      case SLANG_COMPLEX_TYPE:
        break;

      case SLANG_STRING_TYPE:
        break;

      default:
        SLang_verror (SL_TYPE_MISMATCH, "print_array: %s is not supported",
                      SLclass_get_datatype_name (type));
        SLang_free_array (at);
        return;
     }

   fp = NULL;

   if ((num_rows > screen_rows) || (num_cols > screen_cols))
     fp = isis_open_pager ();

   if (fp == NULL) fp = stdout;

   v = at->data;

   just_one_line = 0;
   num = 0;
   for (i = 0; i < num_rows; i++)
     {
        unsigned int j;

        for (j = 0; j < num_cols; j++)
          {
             int ok;

             switch (type)
               {
                case SLANG_CHAR_TYPE:
                  ok = fprintf (fp, "%d\t", (int)((char *)v)[num]);
                  break;

                case SLANG_UCHAR_TYPE:
                  ok = fprintf (fp, "%d\t", (int)((unsigned char *)v)[num]);
                  break;

                case SLANG_SHORT_TYPE:
                  ok = fprintf (fp, "%hd\t", ((short *)v)[num]);
                  break;

                case SLANG_USHORT_TYPE:
                  ok = fprintf (fp, "%hu\t", ((unsigned short *)v)[num]);
                  break;

                case SLANG_INT_TYPE:
                  ok = fprintf (fp, "%d\t", ((int *)v)[num]);
                  break;

                case SLANG_UINT_TYPE:
                  ok = fprintf (fp, "%u\t", ((unsigned int *)v)[num]);
                  break;

                case SLANG_LONG_TYPE:
                  ok = fprintf (fp, "%ld\t", ((long *)v)[num]);
                  break;

                case SLANG_ULONG_TYPE:
                  ok = fprintf (fp, "%lu\t", ((unsigned long *)v)[num]);
                  break;

                case SLANG_FLOAT_TYPE:
                  ok = fprintf (fp, "%e  ", ((float *)v)[num]);
                  break;

                case SLANG_DOUBLE_TYPE:
                  ok = fprintf (fp, "%e  ", ((double *)v)[num]);
                  break;

                case SLANG_STRING_TYPE:
                  ok = fprintf (fp, "\"%s\"  ", ((char **)v)[num]);
                  break;

                case SLANG_COMPLEX_TYPE:
                  ok = fprintf (fp, "(%e, %e)  ",
                                ((double *)v)[num], ((double *)v)[num+1]);
                  num++;
                  break;
                default:
                  ok = -1;
               }

             if (ok <= 0)
               goto done;

             num++;
          }

        if (fputs ("\n", fp) < 0)
          break;

        if ((Isis_Batch_Mode == 0)
            && (fp == stdout)
            && (((num_rows > screen_rows && (0 == (num_rows % screen_rows))) && (num_rows != 0))
                || just_one_line))
          {
             unsigned int key;

             if (just_one_line == 0)
               fprintf (stdout, "Press SPACE to continue");
             fflush (stdout);

             key = isis_getkey ();
             if (key == ' ')
               just_one_line = 0;
             else if (key == '\r')
               just_one_line = 1;
             else break;
          }
     }

   done:

   if (fp != stdout)
     isis_close_pager (fp);

   fputs ("\n", stdout);
   SLang_free_array (at);
}


/*}}}*/

/*{{{ Intrinsics */

static SLang_Intrin_Var_Type Misc_Intrin_Vars [] =
{
   MAKE_VARIABLE("Verbose_Fitsio", &Cfits_Verbose, SLANG_INT_TYPE, 0),
   MAKE_VARIABLE("Remove_Spectrum_Gaps", &Isis_Remove_Spectrum_Gaps, SLANG_INT_TYPE, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define U SLANG_UINT_TYPE
#define F SLANG_FLOAT_TYPE
#define D SLANG_DOUBLE_TYPE
#define S SLANG_STRING_TYPE

static SLang_Intrin_Fun_Type Misc_Intrinsics [] =
{
   MAKE_INTRINSIC("print_array", print_array, V, 0),
   MAKE_INTRINSIC_1("_isis_set_pager", isis_set_pager, V, S),
   MAKE_INTRINSIC("_isis_get_pager", isis_get_pager, V, 0),
   MAKE_INTRINSIC_2("_find_file_in_path", find_file_in_path, V, S, S),
   MAKE_INTRINSIC("_readcol", _readcol, V, 0),
   MAKE_INTRINSIC_I("_isis_set_errno", set_errno, V),
   SLANG_END_INTRIN_FUN_TABLE
};

static char *Pivs = ISIS_VERSION_STRING;
static char *Install_Prefix = INSTALL_PREFIX ;
static char *Install_Prefix_Input = INSTALL_PREFIX_INPUT ;

static SLang_Intrin_Var_Type Global_Intrin_Vars [] =
{
   MAKE_VARIABLE("_isis_version_string", &Pivs, SLANG_STRING_TYPE, 1),
   MAKE_VARIABLE("_isis_version", &Isis_Version, SLANG_UINT_TYPE, 1),
   MAKE_VARIABLE("_isis_srcdir", &Isis_Srcdir, SLANG_STRING_TYPE, 1),
   MAKE_VARIABLE("_isis_install_prefix", &Install_Prefix, SLANG_STRING_TYPE, 1),
   MAKE_VARIABLE("_isis_install_prefix_sans_subdir", &Install_Prefix_Input, SLANG_STRING_TYPE, 1),
   MAKE_VARIABLE("Isis_Batch_Mode", &Isis_Batch_Mode, SLANG_INT_TYPE, 1),
   MAKE_VARIABLE("Isis_Silent_Errors", &Isis_Silent_Errors, SLANG_INT_TYPE, 0),
   MAKE_VARIABLE("Isis_Verbose", &Isis_Verbose, SLANG_INT_TYPE, 0),
   MAKE_VARIABLE("Isis_Load_File_Verbose_Mask", &Isis_Load_File_Verbose_Mask, SLANG_INT_TYPE, 0),
   MAKE_VARIABLE("Isis_Trace", &Isis_Trace, SLANG_INT_TYPE, 0),   
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_IConstant_Type Misc_Intrin_Const [] =
{
   MAKE_ICONSTANT("_INFO", INFO),
   MAKE_ICONSTANT("_WARN", WARN),
   MAKE_ICONSTANT("_FAIL", FAIL),
   MAKE_ICONSTANT("_FATAL", FATAL),
   MAKE_ICONSTANT("_INT_MAX", INT_MAX),
   SLANG_END_ICONST_TABLE
};

#undef V
#undef I
#undef F
#undef D
#undef S
#undef U

/*}}}*/

SLANG_MODULE(miscio);
int init_miscio_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if ((-1 == SLns_add_intrin_fun_table (ns, Misc_Intrinsics, NULL))
       || (-1 == SLns_add_intrin_var_table (pub_ns, Misc_Intrin_Vars, NULL))
       || (-1 == SLadd_intrin_var_table (Global_Intrin_Vars, NULL))
       || (-1 == SLns_add_iconstant_table (ns, Misc_Intrin_Const, NULL)))
     return isis_trace_return(-1);

   return 0;
}

void deinit_miscio_module (void)
{
   ISIS_FREE (Isis_Pager);
}

