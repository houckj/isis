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

/* $Id: util.c,v 1.22 2004/04/21 23:35:42 houck Exp $ */

/*{{{ includes */

#include "config.h"

/* Need this to get NAN and INFINITY definitions from glibc math.h */
#if defined(__GNUC__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef HAVE_DLFCN_H
#  include <dlfcn.h>
#endif

#include <slang.h>

#include "isis.h"
#include "isismath.h"
#include "svd.h"
#include "util.h"
#include "errors.h"

/*}}}*/

/*{{{ globals */
int Isis_Verbose;
int Isis_Errno;
int Isis_Batch_Mode;
int Isis_Remove_Spectrum_Gaps;

char *Isis_Pager;
char *Isis_Srcdir;
void (*Isis_Errno_Hook)(int);
int (*Isis_User_Break_Hook)(void);

/*}}}*/

void *isis_malloc(size_t size) /*{{{*/
{
   void *m;
   if (NULL == (m = malloc (size)))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "allocating %ld bytes", size);
     }
   return m;
}

/*}}}*/

void *isis_realloc(void *ptr, size_t size) /*{{{*/
{
   void *m;
   if (NULL == (m = realloc (ptr, size)))
     {
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__, "re-allocating %ld bytes", size);
     }
   return m;
}

/*}}}*/

/*{{{ Retrieving NaN */

double isis_nan (void)
{
#if defined(NAN)
   return NAN;
#else
   return DBL_MAX;
#endif
}

/*}}}*/

void isis_set_errno (int err) /*{{{*/
{
   /* always call the error hook */
   if (Isis_Errno_Hook != NULL)
     (*Isis_Errno_Hook)(err);

   /* don't stomp the initial error */
   if ((Isis_Errno != 0) && (err != 0))
     return;

   Isis_Errno = err;
}

/*}}}*/

int isis_user_break (void) /*{{{*/
{
   if (Isis_User_Break_Hook == NULL)
     return 0;

   return (*Isis_User_Break_Hook)();
}

/*}}}*/

int isis_fclose (FILE *fp) /*{{{*/
{
   int ret;

   if (0 != (ret = fclose (fp))
       || (errno != 0))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "write failed (%s)",
                    errno ? strerror (errno) : "unknown error");
     }

   return ret;
}

/*}}}*/

/* was SLhypot from John D. */
double isis_hypot (double x, double y) /*{{{*/
{
   double fr, fi, ratio;

   fr = fabs(x);
   fi = fabs(y);

   if (fr > fi)
     {
        ratio = y / x;
        x = fr * sqrt (1.0 + ratio * ratio);
     }
   else if (fi == 0.0)
     {
        x = 0.0;
     }
   else
     {
        ratio = x / y;
        x = fi * sqrt (1.0 + ratio * ratio);
     }

   return x;
}

/*}}}*/

double isis_kahan_sum (double *x, unsigned int n) /*{{{*/
{
   double s, c, y, t;
   unsigned int i;

/* Beware of optimizers that break this */

   s = x[0];
   c = 0.0;

   for (i = 1; i < n; i++)
     {
        y = x[i] - c;
        t = s + y;
        c = (t - s) - y;
        s = t;
     }

   return s;
}

/*}}}*/

double isis_kahan_sum_squares (double *x, unsigned int n) /*{{{*/
{
   double s, c, y, t;
   unsigned int i;

/* Beware of optimizers that break this */

   s = x[0] * x[0];
   c = 0.0;

   for (i = 1; i < n; i++)
     {
        y = x[i] * x[i] - c;
        t = s + y;
        c = (t - s) - y;
        s = t;
     }

   return s;
}

/*}}}*/

int isis_svd_solve (double **matrix, unsigned int n, double *b) /*{{{*/
{
   double *a, *t;
   unsigned int i;
   int ret;

   a = JDMdouble_vector (n * n);
   t = JDMdouble_vector (n);
   if ((a == NULL) || (t == NULL))
     {
        ISIS_FREE (a);
        ISIS_FREE (t);
        return -1;
     }

   memcpy ((char *)t, (char *)b, n * sizeof(double));
   for (i = 0; i < n; i++)
     memcpy ((char *)&a[i*n], (char *)matrix[i], n * sizeof(double));

   ret = svd_solve (a, n, n, t, b);

   ISIS_FREE (a);
   ISIS_FREE (t);

   return ret;
}

/*}}}*/

/* the input matrix 'matrix' will be over-written */
int isis_lu_solve (double **matrix, unsigned int n, unsigned int *piv, double *b) /*{{{*/
{
#define TOLERANCE 1.e-23
   if (-1 == JDM_lu_decomp (matrix, n, piv, TOLERANCE, NULL)
       || -1 == JDM_lu_backsubst (matrix, n, piv, b))
     {
        return -1;
     }

   return 0;
}

/*}}}*/

/*{{{ strings */

char *isis_make_string (const char *str) /*{{{*/
{
   int len;
   char *cpy;

   if (NULL == str)
     return NULL;

   len = strlen(str) + 1;
   if (NULL == (cpy = (char *) ISIS_MALLOC (len * sizeof(char))))
     return NULL;
   isis_strcpy (cpy, str, len);

   return cpy;
}

/*}}}*/

int isis_strcasecmp (const char *pa, const char *pb) /*{{{*/
{
   /* NULL matches nothing */
   if (NULL == pa || NULL == pb)
     return -1;

   while (*pa != '\0')
     {
        if (tolower ((unsigned char)*pa++) != tolower ((unsigned char)*pb++))
          return -1;
     }

   /* strings different */
   if (*pb != '\0') return -1;

   /* strings the same */
   return 0;
}

/*}}}*/

/* this is much like strncpy except
 *  1) it guarantees a null-terminated result,
 *  2) it returns a ptr to the end of the result string and,
 *  3) it has some tolerance for NULL ptrs
 */

char *isis_strcpy (char *dest, const char *src, int size) /*{{{*/
{
   if (NULL == dest || NULL == src || size < 1)
     return NULL;

   while (*src != '\0' && size-- > 0)
     {
        *dest++ = *src++;
     }

   if (size == 0) dest--;

   *dest = '\0';

   /* return ptr to terminating null char */
   return dest;
}

/*}}}*/

/*
 *  1) last arg of isis_strcat () _must_ be a NULL ptr
 *  2) < size characters will be copied to dest
 *  3) dest is guaranteed to be null-terminated.
 *  4) returns -1 if size is too small to hold all the strings.
 *                In that case, dest contains everything that would
 *                fit in a null-terminated string of the stated size.
 *     returns 0 if the concatenation was completely successful.
 */

int isis_strcat (char *dest, int size, ...) /*{{{*/
{
   char *p, *t, *s;
   va_list ap;

   if (dest == NULL || size < 1)
     return -1;

   p = dest;
   t = dest + size - 1;

   while (*p != 0 && p < t)
     p++;
   if (p == t) return -1;

   va_start (ap, size);
   while (NULL != (s = va_arg (ap, char *)))
     {
        while (*s != '\0' && p < t)
          {
             *p++ = *s++;
          }
        if (p == t) break;
     }
   va_end (ap);

   *p = '\0';

   return ((s == NULL) ? 0 : -1);
}

/*}}}*/

char *isis_mkstrcat (const char *arg, ...) /*{{{*/
{
   va_list ap;
   enum {MAX_ARGS = 20};
   const char *s[MAX_ARGS];
   char *p = NULL;
   char *t;
   int i, n, len;

   if (arg == NULL)
     return NULL;

   va_start (ap, arg);

   s[0] = arg;
   n = 1;
   len = strlen(arg);

   for (i = 1; i < MAX_ARGS; i++)
     {
        s[i] = va_arg (ap, char *);
        if (s[i] == NULL)
          break;
        len += strlen(s[i]);
        n++;
     }

   if (NULL != (p = (char *) ISIS_MALLOC ((len + 1) * sizeof(char))))
     {
        t = p;
        for (i = 0; i < n; i++)
          {
             t = isis_strcpy (t, s[i], len+1);
          }
     }

   va_end (ap);
   return p;
}

/*}}}*/

/*}}}*/

/*{{{ messages */

int isis_vsnprintf (char *buf, unsigned int bufsize, /*{{{*/
                    const char *fmt, va_list ap)
{
#ifdef HAVE_VSNPRINTF
   if (-1 == vsnprintf (buf, bufsize, fmt, ap))
     return -1;
#else
   if (strlen(buf) >= bufsize)
     return -1;
   vsprintf (buf, fmt, ap);
#endif
   return 0;
}

/*}}}*/

/*}}}*/

int isis_system (char *cmd) /*{{{*/
{
   int ret;
   if (cmd == NULL)
     return -1;

   ret = SLsystem (cmd);

#if defined(WIFEXITED) && defined(WEXITSTATUS)
   if ((ret != -1) && WIFEXITED(ret))
     ret = WEXITSTATUS(ret);
#endif

   return ret;
}

/*}}}*/

FILE *isis_open_pager (void) /*{{{*/
{
   FILE * fp;

   if (Isis_Pager == NULL)
     {
        isis_vmesg (FAIL, I_INTERNAL, __FILE__, __LINE__, "null pager; using stdout");
        return stdout;
     }

   fp = (FILE *) popen (Isis_Pager, "w");
   if (NULL == fp)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "opening pager %s; using stdout", Isis_Pager);
        return stdout;
     }

   return fp;
}

/*}}}*/

void isis_close_pager (FILE *fp) /*{{{*/
{
   if ((fp == NULL) || (fp == stdout))
     return;

   pclose (fp);
}

/*}}}*/

int is_regular_file (char *file) /*{{{*/
{
#ifdef HAVE_STAT
   struct stat s;

   if (file == NULL)
     return -1;

   if (-1 == stat (file, &s))
     {
        isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "Couldn't stat file %s", file);
        return 0;
     }

   return S_ISREG(s.st_mode);
#else
   (void) file;
   return 1;
#endif
}

/*}}}*/

/*{{{ make temporary file name */

static unsigned int _Fast_Random;
static double _urand (void)
{
   _Fast_Random = _Fast_Random * 69069U + 1013904243U;
   return ((double)_Fast_Random/4294967296.0);
}

static int make_temp_filename (char *buf, unsigned int size) /*{{{*/
{
   char tmp[BUFSIZE];
   struct stat s;

   if (buf == NULL)
     return -1;

   do
     {
        int pid = (int) getpid ();
        int num = (int) (_urand () * 10000.0);
        sprintf (tmp, "/tmp/.isis.%d.%d", pid, num);
     }
   while (0 == stat (tmp, &s));

   isis_strcpy (buf, tmp, size);

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ spawn editor on file */

static int edit_file (char *fname) /*{{{*/
{
   static char editor[] = "${EDITOR:-vi}";     /* environment or default */
   char *cmd = NULL;
   int len, ret = -1;

   if (NULL == fname)
     return -1;

   len = strlen (editor) + strlen (fname) + 3;

   cmd = (char *) ISIS_MALLOC (len * sizeof(char));
   if (NULL == cmd)
     return -1;

   cmd[0] = 0;
   if (-1 == isis_strcat (cmd, len, editor, " ", fname, NULL))
     {
        ISIS_FREE (cmd);
        return -1;
     }

   ret = isis_system(cmd);
   switch (ret)
     {
      case -1:
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__,
                    "check EDITOR environment variable");
        break;
      case 127:
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__,
                    "The editor could not be found.");
        ret = -1;
        break;
      case 126:
        isis_vmesg (INTR, I_FAILED, __FILE__, __LINE__,
                    "The editor was found, but could not be executed.");
        ret = -1;
        break;
      default:
        if (ret)
          {
             isis_vmesg (WARN, I_INFO, __FILE__, __LINE__,
                         "The editor returned a non-zero status (%d)",
                         ret);
             ret = 0;
          }
        break;
     }

   ISIS_FREE (cmd);
   fputs ("\n", stdout);

   return ret;
}

/*}}}*/

int edit_temp_file (int (*save_file)(char *), int (*load_file)(char *), char *file) /*{{{*/
{
   struct stat buf;
   char fname[BUFSIZE];
   int ret = 0;
   int ok_to_save = 1;
   int load_ret = 0;
   char line[4];

   if ((file == NULL)
       && (-1 == make_temp_filename (fname, sizeof(fname))))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "Creating %s", fname);
        return -1;
     }
   else
     {
        if (file)
          isis_strcpy (fname, file, BUFSIZE);

        if (0 == stat (fname, &buf))
          {
             ok_to_save = 0;
             *line = 0;
             fprintf (stderr, "\nFile %s exists -- overwrite it now (y/n)?:  ", fname);
             if ((line == fgets (line, sizeof(line), stdin))
                 && (*line == 'y' || *line == 'Y'))
               ok_to_save = 1;
          }
     }

   if (ok_to_save && (-1 == (*save_file) (fname)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "Creating %s", fname);
        return -1;
     }

   do
     {
        ret = edit_file (fname);

        if (ret == 0)
          load_ret = (*load_file) (fname);

        if (load_ret == 1)
          {
             *line = 0;
             fputs ("\nRe-edit (y/n)?:  ", stderr);
             if ((line == fgets (line, sizeof(line), stdin))
                 && (*line == 'n' || *line == 'N'))
               {
                  load_ret = -1;
               }
             else
               {
                  if (SLang_get_error())
                    SLang_restart (1);
               }
          }

     } while ((ret == 0) && (load_ret == 1));

   if (file == NULL)
     remove (fname);

   if (ret)
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "Editing %s", fname);

   return ret;
}

/*}}}*/

/*}}}*/

/*{{{ search */

int bsearch_d (double t, double *x, int n) /*{{{*/
{
   int n0, n1, n2;
   double xt;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        xt = x[n2];
        if (t <= xt)
          {
             if (xt == t) return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

/*}}}*/

/* histogram bin containing x: */

int find_bin (double x, double *lo, double *hi, int n) /*{{{*/
{
   int n0, n1, n2;

   if (NULL == lo || NULL == hi || n <= 0)
     return -1;

   if (isnan(x))
     return -1;

   if (x < lo[0] || hi[n-1] < x)
     return -1;

   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        if (x < hi[n2])
          {
             if (lo[n2] <= x)
               return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

/*}}}*/

/*}}}*/

/*{{{ array reversal */

int reverse_d (double *x, unsigned int n) /*{{{*/
{
   unsigned int i;

   if (NULL == x || n == 0)
     return -1;

   for (i=0; i < n/2; i++)
     {
        unsigned int j = n - i - 1;
        double tmp = x[j];
        x[j] = x[i];
        x[i] = tmp;
     }

   return 0;
}

/*}}}*/

int reverse_f (float *x, unsigned int n) /*{{{*/
{
   unsigned int i;

   if (NULL == x || n == 0)
     return -1;

   for (i=0; i < n/2; i++)
     {
        unsigned int j = n - i - 1;
        float tmp = x[j];
        x[j] = x[i];
        x[i] = tmp;
     }

   return 0;
}

/*}}}*/

int reverse_i (int *x, unsigned int n) /*{{{*/
{
   unsigned int i;

   if (NULL == x || n == 0)
     return -1;

   for (i=0; i < n/2; i++)
     {
        unsigned int j = n - i - 1;
        int tmp = x[j];
        x[j] = x[i];
        x[i] = tmp;
     }

   return 0;
}

/*}}}*/

int reverse_u (unsigned int *x, unsigned int n) /*{{{*/
{
   unsigned int i;

   if (NULL == x || n == 0)
     return -1;

   for (i=0; i < n/2; i++)
     {
        unsigned int j = n - i - 1;
        unsigned int tmp = x[j];
        x[j] = x[i];
        x[i] = tmp;
     }

   return 0;
}

/*}}}*/

/*}}}*/

/* Gaussian Probability Function */

double isis_gpf (double x) /*{{{*/
{
  /*
   *  P(x) = (1/sqrt(2 pi)) * Integral( -Inf, x, exp(-0.5*t^2) )
   *
   *  Uses fit coefficients from Abramowitz & Stegun, p. 932
   *
   *  Function value is defined for all x and (according to A&S)
   *  should be accurate to within +/- 7.5e-8
   *
   */

   double t, p;

#define XMAX 6.0

   if (x > XMAX)
     return 1.0;
   else if (x < -XMAX)
     return 0.0;

   t = 1.0 / (1.0 + 0.2316419 * fabs(x));

   p = t * (0.319381530 + t *
            (-0.356563782 + t *
             (1.781477937 + t *
              (-1.821255978 + t * 1.330274429))));

   p *= 0.3989422804 * exp(-0.5 * x * x);         /*  0.3989.. = 1/sqrt(2 * PI) */

   if (x < 0.0)
     return p;
   else
     return (1.0 - p);
}

/*}}}*/

/*{{{ Physical unit conversion for energies and wavelengths */

/* CGS conversion factors for supported physical units: */

typedef struct _Unit_t
{
   double cgs_value;    /* e.g. cm/angstrom or erg/eV */
   int unit;            /* integer key */
   int is_wavelength;   /* 1 if yes, 0 if no */
   const char *name;
}
Unit_t;

#define CM_PER_ANGSTROM    1.0e-8
#define CM_PER_NANOMETER   1.0e-7
#define CM_PER_MILLIMETER   0.1
#define CM_PER_CENTIMETER   1.0
#define CM_PER_METER       100.0
#define ERG_PER_HZ         PLANCK
#define ERG_PER_KHZ   (PLANCK * 1.e3)
#define ERG_PER_MHZ   (PLANCK * 1.e6)
#define ERG_PER_GHZ   (PLANCK * 1.e9)
#define ERG_PER_KEV   (ERG_PER_EV * 1.e3)
#define ERG_PER_MEV   (ERG_PER_EV * 1.e6)
#define ERG_PER_GEV   (ERG_PER_EV * 1.e9)
#define ERG_PER_TEV   (ERG_PER_EV * 1.e12)

static Unit_t Unit_Table[] =
{
     {CM_PER_ANGSTROM,   U_ANGSTROM,  1,  "Angstrom"},
     {CM_PER_ANGSTROM,   U_ANGSTROM,  1,  "A"},
     {CM_PER_NANOMETER,  U_NANOMETER, 1,  "nm"},
     {CM_PER_MILLIMETER, U_MILLIMETER,1,  "mm"},
     {CM_PER_CENTIMETER, U_CENTIMETER,1,  "cm"},
     {CM_PER_METER,      U_METER,     1,  "m"},
     {ERG_PER_EV,        U_EV,        0,  "eV"},
     {ERG_PER_KEV,       U_KEV,       0,  "keV"},
     {ERG_PER_MEV,       U_MEV,       0,  "MeV"},
     {ERG_PER_GEV,       U_GEV,       0,  "GeV"},
     {ERG_PER_TEV,       U_TEV,       0,  "TeV"},
     {ERG_PER_HZ,        U_HZ,        0,  "Hz"},
     {ERG_PER_KHZ,       U_KHZ,       0,  "kHz"},
     {ERG_PER_MHZ,       U_MHZ,       0,  "MHz"},
     {ERG_PER_GHZ,       U_GHZ,       0,  "GHz"},
     {-1.0, -1, 1, NULL}
};

static Unit_t * unit_data (int unit) /*{{{*/
{
   Unit_t * p = Unit_Table;

   while (p->cgs_value > 0.0)
     {
        if (p->unit == unit)
          return p;
        p++;
     }

   return NULL;
}

/*}}}*/

int unit_is_wavelength (int unit) /*{{{*/
{
   Unit_t *p = unit_data (unit);

   if (p == NULL)
     return 1;
   else
     return p->is_wavelength;
}

/*}}}*/

int unit_id (char *name) /*{{{*/
{
   Unit_t *p = Unit_Table;

   while (p->cgs_value > 0.0)
     {
        if (0 == isis_strcasecmp (name, p->name))
          return p->unit;
        p++;
     }

   isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "Unrecognized unit %s", name);

   return -1;
}

/*}}}*/

int unit_invalid (int unit) /*{{{*/
{
   if (NULL != unit_data (unit))
     return 0;

   isis_vmesg (FAIL, I_INFO, __FILE__, __LINE__, "Unrecognized unit");
   return 1;
}

/*}}}*/

int unit_name (char * name, int unit) /*{{{*/
{
   Unit_t *p;

   if (NULL == (p = unit_data(unit)))
     return -1;

   strcpy (name, p->name);
   return 0;
}

/*}}}*/

/* unit conversion for derivative df/dx when x-units change */
int unit_convert_dfdx (double *dfdx, double x_in, int out_units, int in_units) /*{{{*/
{
   Unit_t *p_in, *p_out;
   double factor;

   if (dfdx == NULL
       || NULL == (p_in = unit_data (in_units))
       || NULL == (p_out = unit_data (out_units)))
     return -1;

   if (p_in->is_wavelength == p_out->is_wavelength)
     {
        factor = p_in->cgs_value / p_out->cgs_value;
     }
   else
     {
        if (x_in <= 0.0)
          return -1;
        factor = (PLANCK * CLIGHT / p_in->cgs_value) / p_out->cgs_value;
        factor = factor / x_in / x_in;
     }

   *dfdx /= factor;

   return 0;
}

/*}}}*/

int unit_convert_x (double *x, int n, int out_units, int in_units) /*{{{*/
{
   Unit_t *p_in, *p_out;
   double factor;
   int i;

   if (NULL == x || n <= 0
       || NULL == (p_in = unit_data (in_units))
       || NULL == (p_out = unit_data (out_units)))
     return -1;

   if (p_in->is_wavelength == p_out->is_wavelength)
     {
        factor = p_in->cgs_value / p_out->cgs_value;

        for (i=0; i < n; i++)
          x[i] *= factor;
     }
   else
     {
        factor = (PLANCK * CLIGHT / p_in->cgs_value) / p_out->cgs_value;

        for (i=0; i < n; i++)
          {
             if (x[i] <= 0.0)
               return -1;
             x[i] = factor / x[i];
          }
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ rebinning and interpolation */

/* Initial implementation of linear interpolation
 * taken from jdmath by John Davis <davis@space.mit.edu> */
double interpolate_d (double x, double *xp, double *yp, unsigned int n) /*{{{*/
{
   unsigned int n0, n1;
   double x0, x1;

   n0 = bsearch_d (x, xp, n);

   x0 = xp[n0];
   n1 = n0 + 1;

   if (x == x0)
     return yp[n0];
   if (n1 == n)
     {
        if (n0 == 0)
          return yp[n0];
        n1 = n0 - 1;
     }

   x1 = xp[n1];
   if (x1 == x0) return yp[n0];

   return yp[n0] + (yp[n1] - yp[n0]) / (x1 - x0) * (x - x0);
}

/*}}}*/

int interpolate_dvector (double *xp, double *yp, unsigned int n, /*{{{*/
                         double *oldxp, double *oldyp, unsigned int oldn)
{
   double y_0, y_1;
   double x_0, x_1, *xpmax, *oldxpmax;
   double x;

   if (oldn < 2)
     return -1;

   /* I think that this will win if log(oldn) < oldn/n */
   if (n * 10 < oldn)
     {
        unsigned int i;

        for (i = 0; i < n; i++)
          yp[i] = interpolate_d (xp[i], oldxp, oldyp, oldn);
        return 0;
     }

   xpmax = xp + n;
   oldxpmax = oldxp + (oldn - 1);               /* last value */

   oldxp++;
   oldyp++;

   /* Find out where to begin */
   x = *xp;
   while ((oldxp + 5 < oldxpmax)
          && (x > *(oldxp + 5)))
     {
        oldxp += 5;
        oldyp += 5;
     }

   while (xp < xpmax)
     {
        double *oldxp_save = oldxp;

        x = *xp++;

        /* Move along the old axis until x is between two values. */
        while ((oldxp < oldxpmax)
               && (x > *oldxp))
          {
             oldxp++;
          }

        oldyp += (oldxp - oldxp_save);

        y_0 = *(oldyp - 1);
        y_1 = *oldyp;
        x_0 = *(oldxp - 1);
        x_1 = *oldxp;

        /* linear interpolation -- more generally it may be better to do:
         * *yp++ = (*interp_fun) (x, x0, x1, y0, y1);
         */

        /* We have to form the test because the only thing that is assumed
         * is that the x values are ordered.  They may not be unique, */
        if (x_1 == x_0) *yp++ = y_0;
        else *yp++ = y_0 + (y_1 - y_0) * (x - x_0) / (x_1 - x_0);
     }
   return 0;
}

/*}}}*/

#if 0
static double overlap_frac (double flo, double fhi, double tlo, double thi) /*{{{*/
{
   double max_min, min_max;

   if (tlo > flo) max_min = tlo;
   else max_min = flo;

   if (thi < fhi) min_max = thi;
   else min_max = fhi;

   return (min_max - max_min) / (fhi - flo);
}

/*}}}*/
#endif

/* integrate from: (flo, fhi, fy)
               to: (tlo, thi) yielding => ty */
int rebin_histogram (double *fy, double *flo, double *fhi, int nf, /*{{{*/
                     double *ty, double *tlo, double *thi, int nt)
{
   int f, t;

   f = 0;
   for (t = 0; t < nt; t++)
     {
        double t0, t1, s;

        t0 = tlo[t];
        t1 = thi[t];

        /* Accumulate sum over 'from' bins
         * which overlap 'to' bin t
         */
        s = 0.0;
        for ( ;f < nf; f++)
          {
             double f0, f1;
             double min_max, max_min;

             f0 = flo[f];
             f1 = fhi[f];

             if (t0 > f1)
               continue;
             if (f0 > t1)
               break;

             if (t0 > f0) max_min = t0;
             else max_min = f0;

             if (t1 < f1) min_max = t1;
             else min_max = f1;

             /* prevent division by zero */
             if (f0 == f1)
               return -1;

             s += fy[f] * (min_max - max_min) / (f1 - f0);

             /* s += fy[f] * overlap_frac (f0, f1, t0, t1); */

             if (f1 > t1)
               break;
          }

        ty[t] = s;
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ grid validation */

static int examine_grid (double *x, unsigned int num, int *reversed) /*{{{*/
{
   unsigned int i = 0;
   double last_x;
   int order;

   if (num == 0)
     goto binning_error;

   last_x = x[0];
   order = 0;
   for (i = 0; i < num; i++)
     {
        double xi = x[i];

        if ((xi == DBL_MIN)
            || (xi < 0.0))
          goto binning_error;

        if (last_x < xi)
          {
             if (order < 0)
               goto binning_error;
             order++;
          }
        else if (last_x > xi)
          {
             if (order > 0)
               goto binning_error;
             order--;
          }

        last_x = xi;
     }

   *reversed = (order < 0);
   return 0;

   binning_error:
   isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "binning error in bin %u", i);
   return -1;
}

/*}}}*/

/* FIXME -- it would be nice to remove the duplicated functionality
 * between get_canonical_coordinates and _isis_fixup_lo_hi_grids
 */
int _isis_fixup_lo_hi_grids (double *lo, double *hi, unsigned int num, /*{{{*/
                             int *reversed, int *tweaked)
{
   int rev, tmp, num_tweaked;
   unsigned int i, imax;
   double last;

   num_tweaked = 0;

   if (-1 == examine_grid (lo, num, &rev))
     return -1;

   if (-1 == examine_grid (hi, num, &tmp))
     return -1;

   if (rev != tmp)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grid order mismatch");
        return -1;
     }

   *reversed = rev;

   if (rev)
     {
        reverse_d (lo, num);
        reverse_d (hi, num);
     }

#define BIN_OVERLAP_TOL 0.5e-3

   imax = num - 1;
   /* Now make sure lo[i] <= hi[i] <= lo[i+1] */
   for (i = 0; i < imax; i++)
     {
        if (lo[i] > hi[i])
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid grid, bin %d:  lo=%15.8e  hi=%15.8e",
                         i, lo[i], hi[i]);
             return -1;
          }

        if (hi[i] > lo[i+1])
          {
             double dxi, dxip, err;

             dxip = hi[i+1] - lo[i+1];
             dxi  = hi[i] - lo[i];
             err  = hi[i] - lo[i+1];

             if ((dxi <= 0.0)
                 || (dxip <= 0.0)
                 || (err > BIN_OVERLAP_TOL * (dxi + dxip)))
               {
                  isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid grid, overlap of bins %d and %d exceeds tolerance",
                              i, i+1);
                  return -1;
               }

             hi[i] = lo[i+1];
             num_tweaked++;
          }
     }

   if (lo[imax] > hi[imax])
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid grid, bin %d:  lo=%15.8e  hi=%15.8e",
                    imax, lo[imax], hi[imax]);
        return -1;
     }

   /* At this point we know that the grids are in ascending order.  Now
    * tweak them if necessary to get them in monotonically increasing order
    * with non-0 bin widths.
    */

   /* For simplicity, tweak the lo first then adjust the hi accordingly. */
   last = 0.0;
   i = 0;
   while (i < num)
     {
        unsigned int j;
        double next, delta;

        if (lo[i] > last)
          {
             last = lo[i];
             i++;
             continue;
          }

        next = last;
        j = i;
        while ((j < num) && (lo[j] == last))
          j++;                               /* loops at least once */

        if (j == num)
          {
             /* At end of array.  Damn. */
             next = last * 1.001;
             if (next == last)
               { /* punt */
                  isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid grid");
                  return -1;
               }
          }
        else next = lo[j];

        num_tweaked += (j-i);
        delta = ((next - last)/(j - i + 1)) * 0.001;
        while (i < j)
          {
             last = lo[i] = last + delta;
             i++;
          }
     }

   /* Now patch up the hi grid to sane values */
   for (i = 0; i < imax; i++)
     {
        if (hi[i] <= lo[i])
          {
             num_tweaked++;
             hi[i] = lo[i+1];
          }
     }
   if (hi[imax] <= lo[imax])
     {
        if (imax == 0)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "invalid grid");
             return -1;
          }

        hi[imax] = lo[imax] + (hi[imax-1] - lo[imax-1]);
        num_tweaked++;
     }

   if (num_tweaked)
     isis_vmesg (INFO, I_WARNING, __FILE__, __LINE__, "%u hi/lo grid values needed tweaking", num_tweaked);

   if (tweaked != NULL)
     *tweaked = (num_tweaked > 0);
   return 0;
}

/*}}}*/

/* FIXME -- it would be nice to remove the duplicated functionality
 * between get_canonical_coordinates and _isis_fixup_lo_hi_grids
 */
int get_canonical_coordinates (double *bin_lo, double *bin_hi, /*{{{*/
                               int *nbins, int input_units)
{
   int i, is_lambda, ret, closed_gap;

   is_lambda = unit_is_wavelength (input_units);

   /*  first, check for idiocy */

   if (bin_hi[0] == 0.0 && bin_lo[0] == 0.0)
     {
        *nbins = 0;
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grid has zeros");
        return BINS_INVALID;
     }
   if (bin_hi[*nbins-1] == 0.0 && bin_lo[*nbins-1] == 0.0)
     {
        *nbins -= 1;
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grid has zeros");
        if (!is_lambda)   /* is it the CXCDS? */
          return BINS_INVALID;
     }

   /* If we're handed an energy grid, one endpoint value might be zero
    * This is not physical, but people do it anyway.
    */

   if (!is_lambda)
     {
#define HACK_FRACTION 0.1

        if (bin_lo[0] == 0.0)
          bin_lo[0] = HACK_FRACTION * bin_hi[0];
        if (bin_hi[*nbins-1] == 0.0)
          bin_hi[*nbins-1] = bin_lo[*nbins-1] / HACK_FRACTION;
     }

   /*
    *   Grid coordinates must be positive and have non-zero width
    */

   for (i=0; i < *nbins; i++)
     {
        if (bin_lo[i] <= 0.0
            || bin_hi[i] <= 0.0
            || bin_lo[i] == bin_hi[i])
          {
             *nbins = i;
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "corrupted grid (bin %d)", i);
             if (*nbins == 0)
               return BINS_INVALID;
             break;
          }
     }

   /*
    *   Grid must be in Angstrom units
    */

   if (-1 == unit_convert_x (bin_lo, *nbins, U_ANGSTROM, input_units)
       ||-1 == unit_convert_x (bin_hi, *nbins, U_ANGSTROM, input_units))
     return BINS_INVALID;

   /*
    *     Must have bin_lo < bin_hi
    */

   for (i=0; i < *nbins; i++)
     {
        if (bin_lo[i] > bin_hi[i])
          {
             double tmp = bin_lo[i];
             bin_lo[i] = bin_hi[i];
             bin_hi[i] = tmp;
          }
     }

   /*
    *     Grid must be in monotonic increasing order
    */

   /* I'm testing bins in the middle of the array because
    * I've already seen at least one case where the *first* two
    * consecutive bins were screwed up.  Thinking that might
    * also happen with the last two bins, I opted for two in
    * the middle.  This isnt perfect, but there's a limit to
    * the sort of crap I'm going to sort through [JH].
    */

   i = (*nbins > 2) ? (*nbins/2) : 0;

   if ((*nbins == 1) || (bin_lo[i] < bin_lo[i+1]))
     {
        ret = BINS_OK;
     }
   else
     {
        if (-1 == reverse_d (bin_lo, *nbins)
            || -1 == reverse_d (bin_hi, *nbins))
          {
             return BINS_INVALID;
          }

        ret = BINS_REVERSED;
     }

   /*
    *   Gaps not allowed:  must have bin_hi[i-1] = bin_lo[i]
    */

   closed_gap = 0;

   for (i = 1; i < *nbins; i++)
     {
        /* Use cast to float since the bins may have originated from float
         * values.
         */
        if ((float) bin_hi[i-1] == (float) bin_lo[i])
          continue;

        if (Isis_Remove_Spectrum_Gaps)
          {
             closed_gap++;
             bin_hi[i-1] = bin_lo[i];
          }
        else
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                         "gap in grid at bin %d\n => set Remove_Spectrum_Gaps=1 to read it anyway",
                         i);
             return BINS_INVALID;
          }
     }

   if (closed_gap)
     isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "removed %d gaps from input grid", closed_gap);

   return ret;
}

/*}}}*/

int validate_wavelength_grid (int n, double *lo, double *hi) /*{{{*/
{
   int i;
   double last = 0.0;

   if (n <= 0 || lo == NULL || hi == NULL)
     return -1;

   for (i=0; i < n; i++)      /* validate grid */
     {
        if (lo[i] >= hi[i])
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "bin %d has lo >= hi", i);
             return -1;
          }
        if (lo[i] <= last)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "grid must be monotonic increasing");
             return -1;
          }
        last = lo[i];
     }

   return 0;
}

/*}}}*/

/*}}}*/

int transfer_notice (double *flo, double *fhi, int *f_notice_list, int f_num_notice, /*{{{*/
                     double *tlo, double *thi, int tn, int *t_notice)
{
   int i, k;

   /* I purposely do *not* zero out t_notice here.
    * This way, repeated calls can accumulate results
    * in a single array
    */

   i = 0;

   for (k = 0; k < f_num_notice; k++)
     {
        int j = f_notice_list[k];
        double flo_j, fhi_j;

        flo_j = flo[j];
        fhi_j = fhi[j];

        while ((i < tn) && (tlo[i] < fhi_j))
          {
             if (INTERVALS_OVERLAP(flo_j, fhi_j, tlo[i], thi[i]))
               {
                  t_notice[i] = 1;
               }
             i++;
          }
     }

   return 0;
}

/*}}}*/

int _update_notice_list (int *notice, int **notice_list, int *n_notice, int nbins) /*{{{*/
{
   int i, n;

   ISIS_FREE (*notice_list);
   *notice_list = NULL;

   n = 0;
   for (i=0; i < nbins; i++)
     {
        if (notice[i])
          n++;
     }
   *n_notice = n;

   if (n == 0)
     return 0;

   *notice_list = (int *) ISIS_MALLOC (n * sizeof(int));
   if (NULL == *notice_list)
     return -1;
   memset ((char *)*notice_list, 0, n * sizeof(int));

   n=0;
   for (i=0; i < nbins; i++)
     {
        if (notice[i])
          (*notice_list)[n++] = i;
     }

   return 0;
}

/*}}}*/

int unpack_noticed (double *packed, int *notice_list, int n_notice, /*{{{*/
                    int nbins, double *unpacked)
{
   int i, k;

   for (i = 0; i < n_notice; i++)
     {
        k = notice_list[i];
        if (k < nbins)
          unpacked[k] = packed[i];
        else return -1;
     }

   return 0;
}

/*}}}*/

int apply_rebin (double *x, unsigned int num_orig_data, /*{{{*/
                 int *rebin_flags, unsigned int num_rebin_data,
                 double *x_rebin)
{
   unsigned int k, n;

   k = 0;
   while ((k < num_orig_data) && (rebin_flags[k] == 0))
     k++;

   for (n = 0; n < num_rebin_data; n++)
     {
        x_rebin[n] = 0.0;

        if (k < num_orig_data)
          {
             int sign = rebin_flags[k];
             do
               {
                  if (rebin_flags[k] == 0)
                    continue;
                  if (rebin_flags[k] != sign)
                    break;

                  x_rebin[n] += x[k];

               } while (++k < num_orig_data);
          }
     }

   return 0;
}

/*}}}*/

int apply_rebin_weights (double *x, double *wt, unsigned int moment, /*{{{*/
                         unsigned int num_orig_data,
                         int *rebin_flags, unsigned int num_rebin_data,
                         double *x_rebin)
{
   unsigned int k, n, linear;
   double *xwt;

   linear = (moment == 2) ? 0 : 1;

   if (NULL == (xwt = (double *) ISIS_MALLOC (num_orig_data * sizeof(double))))
     return -1;

   if (linear)
     {
        for (k = 0; k < num_orig_data; k++)
          {
             xwt[k] = x[k] * wt[k];
          }
     }
   else
     {
        for (k = 0; k < num_orig_data; k++)
          {
             double xw = x[k] * wt[k];
             xwt[k] = xw * xw;
          }
     }

   k = 0;
   while ((k < num_orig_data) && (rebin_flags[k] == 0))
     k++;

   for (n = 0; n < num_rebin_data; n++)
     {
        double wt_sum = 0.0;
        x_rebin[n] = 0.0;

        if (k < num_orig_data)
          {
             int sign = rebin_flags[k];
             do
               {
                  if (rebin_flags[k] == 0)
                    continue;
                  if (rebin_flags[k] != sign)
                    break;

                  x_rebin[n] += xwt[k];
                  wt_sum += wt[k];

               } while (++k < num_orig_data);
          }

        if (wt_sum != 0.0)
          {
             x_rebin[n] /= linear ? wt_sum : (wt_sum * wt_sum);
          }
     }

   ISIS_FREE(xwt);

   return 0;
}

/*}}}*/

/* generic histogram type */

int Isis_Hist_has_grid (Isis_Hist_t *x) /*{{{*/
{
   return !((x == NULL)
            || ((x->bin_lo == NULL) || (x->bin_hi == NULL))
            || (x->notice_list == NULL));
}

/*}}}*/

void Isis_Hist_free (Isis_Hist_t *x) /*{{{*/
{
   if (x == NULL)
     return;

   ISIS_FREE (x->val);
   ISIS_FREE (x->val_err);
   ISIS_FREE (x->bin_lo);
   ISIS_FREE (x->bin_hi);
   ISIS_FREE (x->notice);
   ISIS_FREE (x->notice_list);
}

/*}}}*/

int Isis_Hist_allocate (unsigned int n, Isis_Hist_t *x) /*{{{*/
{
   unsigned int i;

   if (x == NULL)
     return -1;

   if (NULL == (x->val = (double *) ISIS_MALLOC (n * sizeof(double)))
       || NULL == (x->val_err = (double *) ISIS_MALLOC (n * sizeof(double)))
       || NULL == (x->bin_lo = (double *) ISIS_MALLOC (n * sizeof(double)))
       || NULL == (x->bin_hi = (double *) ISIS_MALLOC (n * sizeof(double)))
       || NULL == (x->notice = (int *) ISIS_MALLOC (n * sizeof(int)))
       || NULL == (x->notice_list = (int *) ISIS_MALLOC (n * sizeof(int))))
     {
        Isis_Hist_free (x);
        return -1;
     }

   x->nbins = n;
   x->n_notice = 0;

   for (i = 0; i < n; i++)
     {
        x->val[i] = 0.0;
        x->val_err[i] = 0.0;
        x->notice[i] = 0;
     }

   return 0;
}

/*}}}*/

int Isis_Hist_push_noticed_grid (Isis_Hist_t *g) /*{{{*/
{
   SLang_Array_Type *sl_binlo=NULL, *sl_binhi=NULL;
   double *flo, *fhi, *tlo, *thi;
   int i, num;
   int *list;

   if (g == NULL)
     return -1;

   if ((NULL == (sl_binlo = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &g->n_notice, 1)))
       ||(NULL == (sl_binhi = SLang_create_array (SLANG_DOUBLE_TYPE, 1, NULL, &g->n_notice, 1))))
     {
        SLang_free_array (sl_binlo);
        SLang_free_array (sl_binhi);
        return -1;
     }

   num = g->n_notice;
   list = g->notice_list;
   flo = g->bin_lo;   fhi = g->bin_hi;
   tlo = (double *) sl_binlo->data;
   thi = (double *) sl_binhi->data;

   for (i = 0; i < num; i++)
     {
        int n = list[i];
        tlo[i] = flo[n];
        thi[i] = fhi[n];
     }

   if ((-1 == SLang_push_array (sl_binlo, 1))
       || (-1 == SLang_push_array (sl_binhi, 1)))
     return -1;

   return 0;
}

/*}}}*/

int Isis_Hist_rebin_noticed (Isis_Hist_t *x, double *x_tmp, Isis_Hist_t *g, double *g_tmp) /*{{{*/
{
   int i;

   if (x == NULL || x_tmp == NULL || g == NULL || g_tmp == NULL)
     return -1;

   if (-1 == unpack_noticed (x->val, x->notice_list, x->n_notice,
                             x->nbins, x_tmp))
     return -1;

   if (-1 == rebin_histogram (x_tmp, x->bin_lo, x->bin_hi, x->nbins,
                              g_tmp, g->bin_lo, g->bin_hi, g->nbins))
     return -1;

   for (i = 0; i < g->n_notice; i++)
     {
        int k = g->notice_list[i];
        g->val[i] = g_tmp[k];
     }

   return 0;
}

/*}}}*/

int Isis_Hist_pop_valid_grid (Isis_Hist_t *g) /*{{{*/
{
   SLang_Array_Type *sl_lo=NULL, *sl_hi=NULL;
   unsigned int i, num;

   if (g == NULL)
     return -1;

   if (-1 == SLang_pop_array_of_type (&sl_hi, SLANG_DOUBLE_TYPE)
       || sl_hi == NULL
       || -1 == SLang_pop_array_of_type (&sl_lo, SLANG_DOUBLE_TYPE)
       || sl_lo == NULL)
     {
        SLang_free_array (sl_hi);
        SLang_free_array (sl_lo);
        return -1;
     }

   num = sl_lo->num_elements;

   if ((num != sl_hi->num_elements)
       || (0 != validate_wavelength_grid (num, (double *)sl_lo->data,
                                          (double *)sl_hi->data)))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "wavelength grid must be monotonic, increasing");
        SLang_free_array (sl_hi);
        SLang_free_array (sl_lo);
        return -1;
     }

   if (-1 == Isis_Hist_allocate (num, g))
     {
        SLang_free_array (sl_hi);
        SLang_free_array (sl_lo);
        return -1;
     }

   memcpy ((char *)g->bin_lo, (char *)sl_lo->data, num*sizeof(double));
   memcpy ((char *)g->bin_hi, (char *)sl_hi->data, num*sizeof(double));

   SLang_free_array (sl_hi);
   SLang_free_array (sl_lo);

   for (i = 0; i < num; i++)
     g->notice[i] = 1;

   if (-1 == _update_notice_list (g->notice, &g->notice_list, &g->n_notice, g->nbins))
     {
        Isis_Hist_free (g);
        return -1;
     }

   return 0;
}

/*}}}*/

int Isis_pop_double_array (double *x, SLindex_Type n) /*{{{*/
{
   SLang_Array_Type *s = NULL;

   if ((-1 == SLang_pop_array_of_type (&s, SLANG_DOUBLE_TYPE))
       || (s == NULL))
     {
        isis_throw_exception (Isis_Error);
        return -1;
     }

   if (s->num_elements != (unsigned int) n)
     {
        isis_vmesg (INTR, I_ERROR, __FILE__, __LINE__, "expected double array of size %d, got size %d\n",
                    n, s->num_elements);
        return -1;
     }

   memcpy ((char *)x, (char *)s->data, n * sizeof(double));
   SLang_free_array (s);

   return 0;
}

/*}}}*/

/* dynamic linking */

#define NO_DYNAMIC_LINKING \
   "Error:  ISIS was not compiled with support for dynamic linking\n"

static void handle_link_error (const char *path, const char *name) /*{{{*/
{
#ifndef HAVE_DLFCN_H
   (void) path;  (void) name;
#else
   const char * error = dlerror ();
   if (error == NULL)
     error = "(unknown)";

   fprintf (stderr, "Link error:  %s\n", error);

   if (name == NULL)
     name = "(null)";
   if (path == NULL)
     path  = "(null)";

   fprintf (stderr, "Link error:  failed loading %s from %s\n",
            name, path);
#endif
}

/*}}}*/

static int is_wrong_abi (void *handle, const char *path, const char *name) /*{{{*/
{
   unsigned int *pv;
   char *s;

#ifdef HAVE_DLFCN_H
   if (NULL == (s = isis_mkstrcat ("Isis_", name, "_api_version", NULL)))
     return -1;
   pv = (unsigned int *) dlsym (handle, s);
   ISIS_FREE(s);
#else
   fprintf (stderr, NO_DYNAMIC_LINKING);
   return -1;
#endif

   if ((pv == NULL) || (*pv != ISIS_API_VERSION))
     {
        fprintf (stderr, "*** Error: module is incompatible: %s\n", path);
        return 2;
     }

   return 0;
}

/*}}}*/

isis_fptr_type isis_load_function (const char *path, const char *name, const char *isis_type) /*{{{*/
{
#ifndef HAVE_DLFCN_H

   (void) path; (void) name; (void) isis_type;
   fprintf (stderr, NO_DYNAMIC_LINKING);
   return NULL;

#else

   isis_fptr_type f = NULL;
   void *handle = NULL;
   char *s = NULL;

   if (path == NULL || name == NULL)
     return NULL;

#ifndef RTLD_GLOBAL
# define RTLD_GLOBAL 0
#endif
#ifdef RTLD_NOW
# define DLOPEN_FLAG  (RTLD_NOW | RTLD_GLOBAL)
#else
# define DLOPEN_FLAG  (RTLD_LAZY | RTLD_GLOBAL)
#endif

   if (NULL == (handle = dlopen (path, DLOPEN_FLAG)))
     {
        handle_link_error (path, NULL);
        return NULL;
     }

   if (isis_type == NULL)
     s = isis_make_string (name);
   else
     {
        if (is_wrong_abi (handle, path, name))
          {
             dlclose (handle);
             return NULL;
          }
        s = isis_mkstrcat ("Isis_", name, "_", isis_type, NULL);
     }

   if (NULL == s)
     {
        dlclose(handle);
        return NULL;
     }

   if (NULL == (f = (isis_fptr_type) dlsym (handle, s)))
     {
        handle_link_error (path, s);
        ISIS_FREE(s);
        dlclose (handle);
        return NULL;
     }

   ISIS_FREE(s);

   return f;

#endif /* HAVE_DLFCN_H */
}
/*}}}*/

/* push_args/pop_args */

int pop_qualifiers_arg (SLang_Struct_Type **sp) /*{{{*/
{
   if (SLang_peek_at_stack() == SLANG_NULL_TYPE)
     {
        *sp = NULL;
        return SLang_pop_null ();
     }
   return SLang_pop_struct (sp);
}

/*}}}*/

void isis_free_args (Isis_Arg_Type *at) /*{{{*/
{
   while (at != NULL)
     {
        Isis_Arg_Type *n = at->next;
        SLang_free_anytype (at->arg);
        ISIS_FREE(at);
        at = n;
     }
}

/*}}}*/

static Isis_Arg_Type *new_arg (void) /*{{{*/
{
   Isis_Arg_Type *at;
   if (NULL == (at = (Isis_Arg_Type *) ISIS_MALLOC (sizeof *at)))
     return NULL;
   at->next = NULL;
   at->arg = NULL;
   at->type = 0;
   return at;
}

/*}}}*/

Isis_Arg_Type *isis_pop_args (int n) /*{{{*/
{
   Isis_Arg_Type *at = NULL;

   while (n-- > 0)
     {
        Isis_Arg_Type *t;
        if (NULL == (t = new_arg()))
          {
             isis_free_args (at);
             return NULL;
          }
        if (-1 == (t->type = SLang_peek_at_stack()))
          {
             isis_free_args (at);
             return NULL;
          }
        SLang_pop_anytype (&t->arg);
        t->next = at;
        at = t;
     }

   return at;
}

/*}}}*/

int isis_push_args (Isis_Arg_Type *at) /*{{{*/
{
   while (at != NULL)
     {
        Isis_Arg_Type *n = at->next;
        SLang_push_anytype (at->arg);
        at = n;
     }

   return 0;
}

/*}}}*/

int isis_coerce_array_to_type (SLang_Array_Type **at, int type) /*{{{*/
{
   SLang_Array_Type *tmp = NULL;
   if (at == NULL || *at == NULL)
     return -1;
   SLang_push_array (*at, 0);
   if (-1 == SLang_pop_array_of_type (&tmp, type))
     return -1;
   SLang_free_array(*at);
   *at = tmp;
   return 0;
}

/*}}}*/
