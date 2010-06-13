/*  Copyright (C) 1998-2004 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research.
    
    Author:  John E. Davis <davis@space.mit.edu>
             John C. Houck <houck@space.mit.edu>
    
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

/* The usual include files */
#include "config.h"
#include <stdio.h>

#include <signal.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <string.h>
/*
 * cpgplot.h contains prototypes with dummy arguments named y1, j1, etc..
 * This causes all sorts of warning about shadowing bessel functions with those
 * names.  Why does the C library use 'y1' instead of 'bessel_y1'??
 * 
 * just to avoid problems -- include math.h afterward -  sigh.
 */
#include <cpgplot.h>
#include <math.h>
#include <slang.h>

#ifdef __cplusplus
extern "C" 
{
#endif
SLANG_MODULE(pgplot);
#ifdef __cplusplus
}
#endif

#include "version.h"
 
/* Almost all of the functions implemented here have exactly the same number
 * of arguments as the pgplot library functions.  The main difference is that
 * the size of the array dimensions are not needed.  Another difference is
 * that temporary workspace arguments are handled transparantly.  For example,
 * the library function PGHI2D takes 12 arguments.  The corresponding slang
 * intrinsic _pghi2d takes only 9 arguments since the 2 of the 12 specify the
 * dimensions of a 2d array, and 1 is a temporary workspace argument.
 * 
 * Finally, functions that return values via the argument list are implemented
 * as slang functions that return multiple values.  An example of this is
 * _pgcurs which is prototyped in FORTRAN as:
 * 
 *      PGCURS (X, Y, CH);
 * 
 * but used in slang as:
 * 
 *      (x,y,ch) = _pgcurs ();
 * 
 * Also see _pgband for another example.
 */



#define DEFAULT_DEVICE	"/XWINDOW"

/* These are some utility routines that convert double arrays to float */
static int pop_float_vector (SLang_Array_Type **at)
{
   *at = NULL;

   if (-1 == SLclass_typecast (SLANG_FLOAT_TYPE, 1, 1))
     return -1;

   if (-1 == SLang_pop_array (at, 1))
     return -1;

   return 0;
}

static void free_arrays (SLang_Array_Type *a, SLang_Array_Type *b,
			 SLang_Array_Type *c, SLang_Array_Type *d)
{
   /* NULLs ok */
   SLang_free_array (a);
   SLang_free_array (b);
   SLang_free_array (c);
   SLang_free_array (d);
}

static int check_vectors (SLang_Array_Type *a, SLang_Array_Type *b)
{
   if (a->num_elements != b->num_elements)
     {
	SLang_verror (SL_TYPE_MISMATCH, "Arrays do not match in size");
	return -1;
     }
   return 0;
}

static int pop_two_float_vectors (SLang_Array_Type **a, SLang_Array_Type **b)
{
   *a = *b = NULL;
   
   if (-1 == pop_float_vector (b))
     return -1;
   if (-1 == pop_float_vector (a))
     {
	SLang_free_array (*b);
	*b = 0;
	return -1;
     }

   if (-1 == check_vectors (*a, *b))
     {
	free_arrays (*a, *b, NULL, NULL);
	*a = *b = NULL;
	return -1;
     }

   return 0;
}

static int pop_three_float_vectors (SLang_Array_Type **a, SLang_Array_Type **b, SLang_Array_Type **c)
{
   *a = *b = *c = NULL;

   if (-1 == pop_float_vector (c))
     return -1;

   if (-1 == pop_two_float_vectors (a, b))
     {
	SLang_free_array (*c);
	*c = NULL;
	return -1;
     }

   if (-1 == check_vectors (*a, *c))
     {
	free_arrays (*a, *b, *c, NULL);
	*a = *b = *c = NULL;
	return -1;
     }
   return 0;
}

static int pop_four_float_vectors (SLang_Array_Type **a, SLang_Array_Type **b,
				   SLang_Array_Type **c, SLang_Array_Type **d)
{
   *a = *b = *c = *d = NULL;

   if (-1 == pop_float_vector (d))
     return -1;

   if (-1 == pop_three_float_vectors (a, b, c))
     {
	SLang_free_array (*d);
	*d = NULL;
	return -1;
     }

   if (-1 == check_vectors (*a, *d))
     {
	free_arrays (*a, *b, *c, *d);
	*a = *b = *c = *d = NULL;
	return -1;
     }
   return 0;
}

/* Here nx corresponds to the fastest varying dimension and ny the slowest */
static SLang_Array_Type *pop_2d_float_array (float **data, unsigned int *ny, unsigned int *nx)
{
   SLang_Array_Type *at;

   *data = NULL;
   *nx = *ny = 0;

   if (-1 == SLclass_typecast (SLANG_FLOAT_TYPE, 1, 1))
     return NULL;

   if (-1 == SLang_pop_array (&at, 1))
     return NULL;

   if (at->num_dims > 2)
     {
	SLang_verror (SL_TYPE_MISMATCH, 
		      "A 2d numeric array is expected");
	SLang_free_array (at);
	return NULL;
     }
   
   *data = (float *)at->data;
   *ny = at->dims[0];
   if (at->num_dims == 1)
     *nx = 1;
   else
     *nx = at->dims[1];

   return at;
}


static int pop_5_floats (float *x1, float *x2, float *x3, float *x4, float *x5)
{
   if ((x5 != NULL)
       && (-1 == SLang_pop_float (x5)))
     return -1;
   if ((x4 != NULL)
       && (-1 == SLang_pop_float (x4)))
     return -1;
   if ((x3 != NULL)
       && (-1 == SLang_pop_float (x3)))
     return -1;
   if ((x2 != NULL)
       && (-1 == SLang_pop_float (x2)))
     return -1;
   if ((x1 != NULL)
       && (-1 == SLang_pop_float (x1)))
     return -1;
   
   return 0;
}

static int pop_4_ints (int *x1, int *x2, int *x3, int *x4)
{
   if ((x4 != NULL)
       && (-1 == SLang_pop_integer (x4)))
     return -1;
   if ((x3 != NULL)
       && (-1 == SLang_pop_integer (x3)))
     return -1;
   if ((x2 != NULL)
       && (-1 == SLang_pop_integer (x2)))
     return -1;
   if ((x1 != NULL)
       && (-1 == SLang_pop_integer (x1)))
     return -1;

   return 0;
}


static int push_2_doubles (double a, double b)
{
   if (-1 == SLang_push_double (a))
     return -1;

   return SLang_push_double (b);
}

static int push_3_doubles (double a, double b, double c)
{
   if (-1 == SLang_push_double (a))
     return -1;

   return push_2_doubles (b, c);
}

static int push_4_doubles (double a, double b, double c, double d)
{
   if (-1 == SLang_push_double (a))
     return -1;

   return push_3_doubles (b, c, d);
}

static int push_2_ints (int a, int b)
{
   if (-1 == SLang_push_integer (a))
     return -1;
   
   return SLang_push_integer (b);
}

static int push_int_2_doubles (int a, double b, double c)
{
   if (-1 == SLang_push_integer (a))
     return -1;
   
   return push_2_doubles (b, c);
}

   
/* The pgplot interface starts here */

static volatile sig_atomic_t Signal_In_Progress;
static void sig_segv (int signo)
{
   static char msg[] =
"\n**** PGPLOT is buggy:  Segmentation Fault occurred while in PGCLOSE\n";
   (void) signo;
   if (Signal_In_Progress)
     return;
   Signal_In_Progress = 1;
   write (STDERR_FILENO, msg, sizeof(msg));
   /* so more SEGVs won't interfere with exit() */
   SLsignal (SIGSEGV, SIG_DFL);
   SLang_exit_error ("Fatal error in PGPLOT--- exiting");
}

static void _pgclos (void)
{
   void (*sig_func)(int);
   
   sig_func = SLsignal (SIGSEGV, sig_segv);
   if (sig_func == SIG_ERR)
     fprintf (stderr, "warning:  failed initializing signal handler for SIGSEGV\n");
   
    cpgclos ();
   
   if (SLsignal (SIGSEGV, sig_func) == SIG_ERR)
     fprintf (stderr, "warning:  failed to re-set signal handler\n");
}


static int _pgopen (char *dev)
{
   int id;

   if (*dev == 0)
     dev = DEFAULT_DEVICE;
   
   id = cpgopen (dev);
   
   return id;			       /* <= 0 is error */
}

static void _pgask (int *i)
{
   cpgask (*i);
}


static void _pgbbuf (void)
{
   cpgbbuf ();
}

static void _pgebuf (void)
{
   cpgebuf ();
}

static void _pgend (void)
{
   cpgend ();
}

static void _pgldev (void)
{
    cpgldev();
}

static void _pgpage (void)
{
    cpgpage();
}

static void _pgupdt (void)
{
    cpgupdt();
}

static void _pgvstd (void)
{
    cpgvstd();
}


static void _pgsave (void)
{
    cpgsave();
}

static void _pgunsa (void)
{
    cpgunsa();
}


static void _pgiden (void)
{
    cpgiden();
}


static void _pgeras (void)
{
    cpgeras();
}

static void _pgetxt (void)
{
    cpgetxt();
}


static int _pgqitf (void)
{
   int i;

   cpgqitf(&i);
   return i;
}

static int _pgqls (void)
{
   int i;

   cpgqls(&i);
   return i;
}

static int _pgqlw (void)
{
   int i;

   cpgqlw(&i);
   return i;
}


static int _pgqcf (void)
{
   int i;

   cpgqcf(&i);
   return i;
}

static double _pgqch (void)
{
   float f;

   cpgqch (&f);
   return (double) f;
}


static int _pgqci (void)
{
   int i;

   cpgqci(&i);
   return i;
}

static int _pgqclp (void)
{
   int i;

   cpgqclp (&i);
   return i;
}

static void _pgqah (void)
{
   float angle, barb;
   int fs;
   
   cpgqah (&fs, &angle, &barb);
   push_int_2_doubles (fs, angle, barb);
}

static void _pgqcir (void)
{
   int lo, hi;
   cpgqcir (&lo, &hi);
   push_2_ints (lo, hi);
}

static void _pgqcol (void)
{
   int lo, hi;
   cpgqcol (&lo, &hi);
   push_2_ints (lo, hi);
}

static void _pgqcr (int *ci)
{
   float r, g, b;
   cpgqcr (*ci, &r, &g, &b);
   push_3_doubles (r, g, b);
}

#if 0
static void _pgqdt ()
{
}
#endif

static void _pgqhs (void)
{
   float ang, sep, phase;
   cpgqhs (&ang, &sep, &phase);
   push_3_doubles (ang, sep, phase);
}

static int _pgqndt (void)
{
   int n;
   cpgqndt (&n);
   return n;
}

static void _pgqpos (void)
{
   float x, y;
   cpgqpos (&x, &y);
   push_2_doubles (x, y);
}

static int _pgqtbg (void)
{
   int n;
   cpgqtbg (&n);
   return n;
}

#if 0
static void _pgqtxt ()
{
}
#endif

static void _pgqvsz (int *units)
{
   float a, b, c, d;
   
   cpgqvsz (*units, &a, &b, &c, &d);
   push_4_doubles (a, b, c, d);
}


static void _pgqcs (int *units)
{
   float xch, ych;
   
   cpgqcs (*units, &xch, &ych);
   
   (void) push_2_doubles (xch, ych);
}

static int _pgqfs (void)
{
   int i;

   cpgqfs(&i);
   return i;
}

static int _pgqid (void)
{
   int i;

   cpgqid(&i);
   return i;
}

static void _pgqinf (char *item)
{
   char buf[1024];
   int len;
   
   len = (int) sizeof(buf);

   cpgqinf (item, buf, &len);
   if (len < 0)
     len = 0;
   if (len >= (int) sizeof (buf))
     len = sizeof (buf) - 1;
   
   buf[len] = 0;
   (void) SLang_push_string (buf);
}

static void _pgqvp (int *units)
{
   float x1, x2, y_1, y2;
   
   cpgqvp (*units, &x1, &x2, &y_1, &y2);

   (void) push_4_doubles (x1, x2, y_1, y2);
}

static void _pgqwin (void)
{
   float x1, x2, y_1, y2;
   
   cpgqwin (&x1, &x2, &y_1, &y2);
   
   (void) push_4_doubles (x1, x2, y_1, y2);
}


/* select an open graphics device */
static void _pgslct (int *id)
{
   cpgslct (*id);
}

/* switch to a different panel on the view surface */
static void _pgpanl (int *ix, int *iy)
{
   cpgpanl (*ix, *iy);
}

/* subdivide view surface into panels */
static void _pgsubp (int *x, int *y)
{
   cpgsubp (*x, *y);
}

/* set window and viewport and draw labeled frame */
static void _pgenv (double *xmin, double *xmax, double *ymin, double *ymax,
		    int *just, int *axis)
{
   cpgenv ((float) *xmin, (float) *xmax, (float) *ymin, (float) *ymax,
	   *just, *axis);
}

/* set window */
static void _pgswin (double *xmin, double *xmax, double *ymin, double *ymax)
{
   cpgswin ((float) *xmin, (float) *xmax, (float) *ymin, (float) *ymax);
}

/*  set viewport (normalized device coordinates) */
static void _pgsvp (double *xmin, double *xmax, double *ymin, double *ymax)
{
   cpgsvp ((float) *xmin, (float) *xmax, (float) *ymin, (float) *ymax);
}

/* set viewport (inches) */
static void _pgvsiz (double *xmin, double *xmax, double *ymin, double *ymax)
{
   cpgvsiz ((float) *xmin, (float) *xmax, (float) *ymin, (float) *ymax);
}

/* set window and adjust viewport to same aspect ratio */
static void _pgwnad (double *xmin, double *xmax, double *ymin, double *ymax)
{
   cpgwnad ((float) *xmin, (float) *xmax, (float) *ymin, (float) *ymax);
}

/* draw labeled frame around viewport */
static void _pgbox (char *xopt, double *xtic, int *nx,
		    char *yopt, double *ytic, int *ny)
{
   cpgbox (xopt, (float) *xtic, *nx, yopt, (float) *ytic, *ny);
}

/* draw frame and write (DD) HH MM SS.S labelling */
static void _pgtbox (char *xopt, double *xtic, int *nx,
		     char *yopt, double *ytic, int *ny)
{
   cpgtbox (xopt, (float) *xtic, *nx, yopt, (float) *ytic, *ny);
}

/* change the size of the view surface */
static void _pgpap (double *w, double *a)
{
   cpgpap ((float) *w, (float) *a);
}

/* set color index */
static void _pgsci (int *i)
{
   cpgsci (*i);
}

/* set line clipping state */
static void _pgsclp (int *i)
{
   cpgsclp (*i);
}

/* set fill-area style */
static void _pgsfs (int *i)
{
   cpgsfs (*i);
}

/* set image transfer function */
static void _pgsitf (int *i)
{
   cpgsitf (*i);
}

/* set text background color index */
static void _pgstbg (int *i)
{
   cpgstbg (*i);
}

/* set line style */
static void _pgsls (int *i)
{
   cpgsls (*i);
}

/* set character font */
static void _pgscf (int *i)
{
   cpgscf (*i);
}

/* set line width */
static void _pgslw (int *i)
{
   cpgslw (*i);
}

/* set hatching style */
static void _pgshs (double *ang, double *s, double *p)
{
   cpgshs ((float) *ang, (float) *s, (float) *p);
}

/* set character height */
static void _pgsch (double *f)
{
   cpgsch((float) *f);
}

/* set arrow-head style */
static void _pgsah (int *i, double *x, double *y)
{
   cpgsah(*i, (float) *x, (float) *y);
}


/* write labels for x-axis, y-axis, and top of plot */
static void _pglab (char *a, char *b, char *c)
{
   cpglab (a, b, c);
}

/* write text at position relative to viewport */
static void _pgmtxt (char *s, double *f, double *g, double *h, char *t)
{
   cpgmtxt (s, (float)*f, (float)*g, (float)*h, t);
}

/* write text at arbitrary position and angle */
static void _pgptxt (double *x, double *y, double *a, double *j, char *s)
{
   cpgptxt ((float) *x, (float) *y, (float) *a, (float) *j, s);
}

/* draw an arrow */
static void _pgarro (double *x, double *y, double *a, double *b)
{
   cpgarro ((float) *x, (float) *y, (float) *a, (float) *b);
}

/* draw a circle, using fill-area attributes */
static void _pgcirc (double *x, double *y, double *a)
{
   cpgcirc ((float) *x, (float) *y, (float) *a);
}

/* move pen (change current pen position) */
static void _pgmove (double *x, double *y)
{
   cpgmove ((float) *x, (float) *y);
}

/* draw a line from the current pen position to a point */
static void _pgdraw (double *x, double *y)
{
   cpgdraw ((float) *x, (float) *y);
}

/* draw a rectangle, using fill-area attributes */
static void _pgrect (double *x, double *y, double *a, double *b)
{
   cpgrect ((float) *x, (float) *y, (float) *a, (float) *b);
}
 

/* draw several graph markers */
static void _pgpt (int *symbol)
{
   SLang_Array_Type *x, *y;
   
   if (-1 == pop_two_float_vectors (&x, &y))
     return;
   
   cpgpt ((int) x->num_elements, (float *)x->data, (float *)y->data, *symbol);
   
   free_arrays (x, y, NULL, NULL);
}

/* draw a polyline (curve defined by line-segments) */
static void _pgline (void)
{
   SLang_Array_Type *x, *y;

   if (-1 == pop_two_float_vectors (&x, &y))
     return;

   cpgline ((int) x->num_elements, (float *)x->data, (float *)y->data);

   free_arrays (x, y, NULL, NULL);
}

/* vertical error bar */
static void _pgerry (double *len)
{
   SLang_Array_Type *x, *y_1, *y_2;

   if (-1 == pop_three_float_vectors (&x, &y_1, &y_2))
     return;

   cpgerry ((int) x->num_elements, (float*)x->data, 
	    (float*)y_1->data, (float *)y_2->data, (float) *len);

   free_arrays (x, y_1, y_2, NULL);
}

/* horizontal error bar */
static void _pgerrx (double *len)
{
   SLang_Array_Type *x1, *x2, *y;

   if (-1 == pop_three_float_vectors (&x1, &x2, &y))
     return;

   cpgerrx ((int) x1->num_elements, 
	    (float*)x1->data, (float*)x2->data, (float*)y->data, (float) *len);

   free_arrays (x1, x2, y, NULL);
}

/* horizontal or vertical error bar */
static void _pgerrb (double *t)
{
   SLang_Array_Type *x, *y, *e;
   int dir;

   if (-1 == pop_three_float_vectors (&x, &y, &e))
     return;

   if (-1 != SLang_pop_integer (&dir))
     {
	cpgerrb (dir, (int)x->num_elements, 
		 (float*)x->data, (float*)y->data, (float*)e->data, 
		 (float) *t);
     }
   
   free_arrays (x, y, e, NULL);
}

/* horizontal or vertical error bar */
static void _pgerr1 (int *dir, double *x, double *y, double *e, double *t)
{
   cpgerr1 (*dir, (float) *x, (float) *y, (float) *e, (float) *t);
}

/* histogram of binned data */
static void _pgbin (int *center)
{
   SLang_Array_Type *x, *data;

   if (-1 == pop_two_float_vectors (&x, &data))
     return;
   
   cpgbin ((int) x->num_elements, (float*)x->data, (float*)data->data, *center);
   
   free_arrays (x, data, NULL, NULL);
}

/* histogram of unbinned data */
static void _pghist (double *dmin, double *dmax, int *nbins, int *flag)
{
   SLang_Array_Type *data;
   
   if (-1 == pop_float_vector (&data))
     return;
   
   cpghist ((int)data->num_elements, 
	    (float*)data->data, (float)*dmin, (float)*dmax, *nbins, *flag);

   SLang_free_array (data);
}

/* read cursor position */
static int _pgcurs (void)
{
   float x, y;
   char ch;
   
   x = y = 0;
   ch = 0;
   
   (void) cpgcurs (&x, &y, &ch);
   (void) push_2_doubles (x, y);
   
   return (int) (unsigned char) ch;
}


/* cross-sections through a 2D data array */
/* Prototype: _pghi2d (data[][], iy1, iy2, ix1, ix2, x[], ioff, bias, center) */
static void _pghi2d (int *ioff, double *bias, int *center)
{
   int ix1, ix2, iy1, iy2;
   SLang_Array_Type *x, *at;
   float *data;
   unsigned int nx, nxv, nyv;
   float *ylims;

   if (-1 == pop_float_vector (&x))
     return;
   
   nx = x->num_elements;

   if ((-1 == pop_4_ints (&iy1, &iy2, &ix1, &ix2))
       || (NULL == (at = pop_2d_float_array (&data, &nyv, &nxv))))
     {
	SLang_free_array (x);
	return;
     }

   if ((int) nx != (ix2 - ix1 + 1))
     {
	free_arrays (at, x, NULL, NULL);
	SLang_verror (SL_INVALID_PARM, 
		      "pghi2d: ix1 and ix2 are out of range");
	return;
     }
   
   if (NULL == (ylims = (float *) SLmalloc (sizeof (float) * (nx + 1))))
     {
	free_arrays (at, x, NULL, NULL);
	return;
     }

   /* Don't forget to add one to the array indices since FORTRAN arrays
    * are indexed from 1
    */
   if (nx && nxv && nyv)
     cpghi2d (data, nxv, nyv, ix1 + 1, ix2 + 1, iy1 + 1, iy2 + 1,
	      (float*)x->data, *ioff, (float) *bias, *center, ylims);
   
   SLfree ((char *) ylims);
   free_arrays (at, x, NULL, NULL);
}

static void _pgtext (double *x, double *y, char *s)
{
   cpgtext ((float) *x, (float) *y, s);
}

static void _pgaxis (void)
{
   char *opt;
   float x1, y_1, x2, y2, v1, v2, step;
   int nsub;
   float dmajl, dmajr, f_min, disp, orient;
   
   if (-1 == pop_5_floats (&dmajl, &dmajr, &f_min, &disp, &orient))
     return;
   if (-1 == SLang_pop_integer (&nsub))
     return;
   if (-1 == pop_5_floats (&x2, &y2, &v1, &v2, &step))
     return;
   if (-1 == pop_5_floats (&x1, &y_1, NULL, NULL, NULL))
     return;
   if (-1 == SLang_pop_slstring (&opt))
     return;
   cpgaxis (opt, x1, y_1, x2, y2, v1, v2, step, 
	    nsub, dmajl, dmajr, f_min, disp, orient);
   SLang_free_slstring (opt);
}

static int _pgband (int *mode, int *posn, double *xref, double *yref, 
		    SLang_Ref_Type *rx, SLang_Ref_Type *ry, SLang_Ref_Type *rc)
{
   float x, y;
   char c;
   int status;
   
   status = cpgband (*mode, *posn, *xref, *yref, &x, &y, &c);
   if (status == 1)
     {
	(void) SLang_assign_to_ref (rx, SLANG_FLOAT_TYPE, &x);
	(void) SLang_assign_to_ref (ry, SLANG_FLOAT_TYPE, &y);
	(void) SLang_assign_to_ref (rc, SLANG_CHAR_TYPE, &c);
     }
   return status;
}

static int pop_tr_vector (SLang_Array_Type **tr)
{
   if (-1 == pop_float_vector (tr))
     return -1;

   if ((*tr)->num_elements != 6)
     {
	SLang_verror (SL_INVALID_PARM, "_pgcon*: expecting a 6 element 'tr' vector");
	SLang_free_array (*tr);
	*tr = NULL;
	return -1;
     }
   return 0;
}

static void do_pgcon_bs (double *blank, int *nc_sign)
{
   unsigned int idim, jdim;
   float *a;
   int i_1, i_2, j_1, j_2;
   SLang_Array_Type *tr, *c, *at;

   at = NULL;

   if (-1 == pop_tr_vector (&tr))
     return;

   if (-1 == pop_float_vector (&c))
     goto return_error;

   if (-1 == pop_4_ints (&j_1, &j_2, &i_1, &i_2))
     goto return_error;

   if (NULL == (at = pop_2d_float_array (&a, &jdim, &idim)))
     goto return_error;

   /* Convert to FORTRAN indexing */
   i_1++; j_1++; i_2++; j_2++;

   if (blank != NULL)
     cpgconb (a, idim, jdim, i_1, i_2, j_1, j_2, 
	      (float *)c->data, c->num_elements,
	      (float *)tr->data, *blank);
   else if (nc_sign != NULL)
     cpgcont (a, idim, jdim, i_1, i_2, j_1, j_2, 
	      (float *)c->data, *nc_sign * (int)c->num_elements,
	      (float *)tr->data);
   else
     cpgcons (a, idim, jdim, i_1, i_2, j_1, j_2, 
	      (float *)c->data, c->num_elements,
	      (float *)tr->data);

   return_error:

   free_arrays (tr, c, at, NULL);
}

static void _pgconb (double *blank)
{
   do_pgcon_bs (blank, NULL);
}

static void _pgcons (void)
{
   do_pgcon_bs (NULL, NULL);
}

/* This differs from the API */
static void _pgcont (int *nc_sign)
{
   int s;
   
   if (*nc_sign >= 0)
     s = 1;
   else 
     s = -1;
   do_pgcon_bs (NULL, &s);
}

static void _pgconl (char *label, int *intval, int *minint)
{
   float *a;
   unsigned int idim, jdim;
   int i_1, i_2, j_1, j_2;
   float c;
   SLang_Array_Type *tr, *at;
   
   tr = at = NULL;

   if (-1 == pop_float_vector (&tr))
     return;

   if (-1 == SLang_pop_float (&c))
     goto return_error;

   if (-1 == pop_4_ints (&j_1, &j_2, &i_1, &i_2))
     goto return_error;

   if (NULL == (at = pop_2d_float_array (&a, &jdim, &idim)))
     goto return_error;
		
   /* Convert to FORTRAN indexing */
   i_1++; j_1++; i_2++; j_2++;

   cpgconl (a, idim, jdim, i_1, i_2, j_1, j_2, 
	    c, (float *)tr->data, 
	    label, *intval, *minint);
   
   return_error:
   
   free_arrays (at, tr, NULL, NULL);
}

static void do_pgconf_xxx (void (*f)(const float *, int, int, int, int, int, int,
				     float, float, const float *))
{
   float *a;
   unsigned int idim, jdim;
   int i_1, i_2, j_1, j_2;
   float c1, c2;
   SLang_Array_Type *at, *tr;
   
   at = tr = NULL;

   if (-1 == pop_float_vector (&tr))
     return;

   if (-1 == pop_5_floats (&c1, &c2, NULL, NULL, NULL))
     goto return_error;

   if (-1 == pop_4_ints (&j_1, &j_2, &i_1, &i_2))
     goto return_error;

   if (NULL == (at = pop_2d_float_array (&a, &jdim, &idim)))
     goto return_error;
   
   /* Convert to FORTRAN indexing */
   i_1++; j_1++; i_2++; j_2++;
   
   (*f) (a, idim, jdim, i_1, i_2, j_1, j_2, 
	 c1, c2, (float *)tr->data);
   
   return_error:
   free_arrays (at, tr, NULL, NULL);
}

static void _pgconf (void)
{
   do_pgconf_xxx (&cpgconf);
}

static void _pggray (void)
{
   do_pgconf_xxx (&cpggray);
}

static void _pgimag (void)
{
   do_pgconf_xxx (&cpgimag);
}

/* Warning: This routine differs from its pgplot counterpart. 
 * It does not allow the use of old arrays.  At most, 1024 points are allocated.
 */
static void _pglcur_pgncur_pgolin (SLang_Ref_Type *rx, SLang_Ref_Type *ry, 
				   int symbol, int what)
{
   SLang_Array_Type *a, *b;
   float x[1024];
   float y[1024];
   SLindex_Type n_it;
   int n;

   n = 0;

   switch (what)
     {
      case 1:
	cpglcur (1024, &n, x, y);
	break;
	
      case 2:
	cpgncur (1024, &n, x, y, symbol);
	break;
	
      case 3:
	cpgolin (1024, &n, x, y, symbol);
	break;
     }
   
   if (n < 0)
     n = 0;

   n_it = n;

   if (NULL == (a = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &n_it, 1)))
     return;
   if (NULL == (b = SLang_create_array (SLANG_FLOAT_TYPE, 0, NULL, &n_it, 1)))
     {
	SLang_free_array (a);
	return;
     }
   
   memcpy ((char *)a->data, (char *)x, n * sizeof (float));
   memcpy ((char *)b->data, (char *)y, n * sizeof (float));
   
   (void) SLang_assign_to_ref (rx, SLANG_ARRAY_TYPE, &a);
   (void) SLang_assign_to_ref (ry, SLANG_ARRAY_TYPE, &b);
   
   free_arrays (a, b, NULL, NULL);
}

static void _pglcur (SLang_Ref_Type *rx, SLang_Ref_Type *ry)
{
   _pglcur_pgncur_pgolin (rx, ry, -1, 1);
}

static void _pgncur (SLang_Ref_Type *rx, SLang_Ref_Type *ry, int *symbol)
{
   _pglcur_pgncur_pgolin (rx, ry, *symbol, 2);
}

static void _pgolin (SLang_Ref_Type *rx, SLang_Ref_Type *ry, int *symbol)
{
   _pglcur_pgncur_pgolin (rx, ry, *symbol, 3);
}

static void _pgpnts (void)
{
   SLang_Array_Type *x, *y;
   SLang_Array_Type *s;
   
   if (-1 == SLang_pop_array_of_type (&s, SLANG_INT_TYPE))
     return;
   if (0 == pop_two_float_vectors (&x, &y))
     {
	cpgpnts (x->num_elements, (float *) x->data, (float *) y->data, 
		 (int *) s->data, s->num_elements);
     }
   free_arrays (x, y, s, NULL);
}

static void _pgpoly (void)
{
   SLang_Array_Type *x, *y;

   if (-1 == pop_two_float_vectors (&x, &y))
     return;
   
   cpgpoly (x->num_elements, (float *) x->data, (float *) y->data);
   
   free_arrays (x, y, NULL, NULL);
}

static void _pgctab (double *contra, double *bright)
{
   SLang_Array_Type *l, *r, *g, *b;

   if (-1 == pop_four_float_vectors (&l, &r, &g, &b))
     return;
   
   cpgctab ((float *)l->data, (float *)r->data, (float *)g->data, (float *)b->data,
	    (int) l->num_elements, *contra, *bright);
   
   free_arrays (l, r, g, b);
}

/* assign RGB color to a color index */
static void _pgscr (int *ci, double *red, double *green, double *blue)
{
   cpgscr (*ci,  (float) *red, (float) *green, (float) *blue);
}

static void _pgscir (int *clo, int *chi)
{
   cpgscir (*clo, *chi);
}


static void _pgscrl (double *dx, double *dy)
{
   cpgscrl ((float) *dx, (float) *dy);
}

static int _pgscrn (int *ci, char *name)
{
   int er;
   cpgscrn (*ci, name, &er);
   return er;
}

/* assign HLS color to a color index (hue/lightness/saturation) */
static void _pgshls (int *ci, double *hue, double *lightness, double *saturation)
{
   cpgshls (*ci, (float) *hue, (float) *lightness,  (float) *saturation);
}

/*
 * Not implemented:
 *   Obsolete:
 *      pgbeg
 */
#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define D SLANG_DOUBLE_TYPE
#define R SLANG_REF_TYPE
static SLang_Intrin_Fun_Type Module_Intrinsics [] = 
{
   MAKE_INTRINSIC_4("_pgarro", _pgarro, V, D, D, D, D),
   MAKE_INTRINSIC_I("_pgask", _pgask, V),
   MAKE_INTRINSIC_0("_pgaxis", _pgaxis, V),
   MAKE_INTRINSIC_7("_pgband", _pgband, I, I, I, D, D, R, R, R),
   MAKE_INTRINSIC_0("_pgbbuf", _pgbbuf, V),
   MAKE_INTRINSIC_I("_pgbin", _pgbin, V),
   MAKE_INTRINSIC_6("_pgbox", _pgbox, V, S, D, I, S, D, I),
   MAKE_INTRINSIC_3("_pgcirc", _pgcirc, V, D, D, D),
   MAKE_INTRINSIC_0("_pgclos", _pgclos, V),
   MAKE_INTRINSIC_1("_pgconb", _pgconb, V, D),
   MAKE_INTRINSIC_3("_pgconl", _pgconl, V, S, I, I),
   MAKE_INTRINSIC_0("_pgconf", _pgconf, V),
   MAKE_INTRINSIC_0("_pgcons", _pgcons, V),
   MAKE_INTRINSIC_1("_pgcont", _pgcont, V, I),
   /* MAKE_INTRINSIC_0("_pgconx", _pgconx, V), */

   MAKE_INTRINSIC_2("_pgctab", _pgctab, V, D, D),
   MAKE_INTRINSIC_0("_pgcurs", _pgcurs, I),
   MAKE_INTRINSIC_2("_pgdraw", _pgdraw, V, D, D),
   MAKE_INTRINSIC_0("_pgebuf", _pgebuf, V),
   MAKE_INTRINSIC_0("_pgend", _pgend, V),
   MAKE_INTRINSIC_6("_pgenv", _pgenv, V, D,D,D,D,I,I),
   MAKE_INTRINSIC_0("_pgeras", _pgeras, V),
   MAKE_INTRINSIC_5("_pgerr1", _pgerr1, V, I, D, D, D, D),
   MAKE_INTRINSIC_1("_pgerrb", _pgerrb, V, D),
   MAKE_INTRINSIC_1("_pgerrx", _pgerrx, V, D),
   MAKE_INTRINSIC_1("_pgerry", _pgerry, V, D),
   MAKE_INTRINSIC_0("_pgetxt", _pgetxt, V),
   MAKE_INTRINSIC_0("_pggray", _pggray, V),

   MAKE_INTRINSIC_3("_pghi2d", _pghi2d, V, I, D, I),
   MAKE_INTRINSIC_4("_pghist", _pghist, V, D, D, I, I),
   MAKE_INTRINSIC_0("_pgiden", _pgiden, V),
   MAKE_INTRINSIC_0("_pgimag", _pgimag, V),
   MAKE_INTRINSIC_3("_pglab", _pglab, V, S, S, S),
   MAKE_INTRINSIC_2("_pglcur", _pglcur, V, R, R),

   MAKE_INTRINSIC_0("_pgldev", _pgldev, V),
   /* MAKE_INTRINSIC_0("_pglen", _pglen) */
   MAKE_INTRINSIC_0("_pgline", _pgline, V),
   MAKE_INTRINSIC_2("_pgmove", _pgmove, V, D, D),
   MAKE_INTRINSIC_5("_pgmtxt", _pgmtxt, V, S, D, D, D, S),
   MAKE_INTRINSIC_3("_pgncur", _pgncur, V, R, R, I),
   /* MAKE_INTRINSIC_0("_pgnumb", _pgnumb) */
   MAKE_INTRINSIC_3("_pgolin", _pgolin, V, R, R, I),
   MAKE_INTRINSIC_S("_pgopen", _pgopen, I),
   MAKE_INTRINSIC_0("_pgpage", _pgpage, V),
   MAKE_INTRINSIC_2("_pgpanl", _pgpanl, V, I, I),
   MAKE_INTRINSIC_2("_pgpap", _pgpap, V, D, D),
   /* MAKE_INTRINSIC_0("_pgpixl", _pgpixl), */
   MAKE_INTRINSIC_0("_pgpnts", _pgpnts, V),
   MAKE_INTRINSIC_0("_pgpoly", _pgpoly, V),

   MAKE_INTRINSIC_I("_pgpt", _pgpt, V),
   MAKE_INTRINSIC_5("_pgptxt", _pgptxt, V, D, D, D, D, S),
   MAKE_INTRINSIC_0("_pgqah", _pgqah, V),
   MAKE_INTRINSIC_0("_pgqcf", _pgqcf, I),
   MAKE_INTRINSIC_0("_pgqch", _pgqch, D),
   MAKE_INTRINSIC_0("_pgqci", _pgqci, I),
   MAKE_INTRINSIC_0("_pgqcir", _pgqcir, V),
   MAKE_INTRINSIC_0("_pgqclp", _pgqclp, I),
   MAKE_INTRINSIC_0("_pgqcol", _pgqcol, V),
   MAKE_INTRINSIC_1("_pgqcr", _pgqcr, V, I),
   MAKE_INTRINSIC_1("_pgqcs", _pgqcs, V, I),
   /* MAKE_INTRINSIC_0("_pgqdt", _pgqdt), */
   MAKE_INTRINSIC_0("_pgqfs", _pgqfs, I),
   MAKE_INTRINSIC_0("_pgqhs", _pgqhs, V),
   MAKE_INTRINSIC_0("_pgqid", _pgqid, I),
   MAKE_INTRINSIC_1("_pgqinf", _pgqinf, V, S),
   MAKE_INTRINSIC_0("_pgqitf", _pgqitf, I),
   MAKE_INTRINSIC_0("_pgqls", _pgqls, I),
   MAKE_INTRINSIC_0("_pgqlw", _pgqlw, I),
   MAKE_INTRINSIC_0("_pgqndt", _pgqndt, I),
   MAKE_INTRINSIC_0("_pgqpos", _pgqpos, V),
   MAKE_INTRINSIC_0("_pgqtbg", _pgqtbg, I),
   /* MAKE_INTRINSIC_0("_pgqtxt", _pgqtxt), */
   MAKE_INTRINSIC_1("_pgqvp", _pgqvp, V, I),
   MAKE_INTRINSIC_1("_pgqvsz", _pgqvsz, V, I),
   MAKE_INTRINSIC_0("_pgqwin", _pgqwin, V),
   MAKE_INTRINSIC_4("_pgrect", _pgrect, V, D, D, D, D),
   /* MAKE_INTRINSIC_0("_pgrnd", _pgrnd), */
   /* MAKE_INTRINSIC_0("_pgrnge", _pgrnge), */
   MAKE_INTRINSIC_3("_pgsah", _pgsah, V, I, D, D),
   MAKE_INTRINSIC_0("_pgsave", _pgsave, V),
   MAKE_INTRINSIC_I("_pgscf", _pgscf, V),
   MAKE_INTRINSIC_1("_pgsch", _pgsch, V, D),
   MAKE_INTRINSIC_I("_pgsci", _pgsci, V),
   MAKE_INTRINSIC_2("_pgscir", _pgscir, V, I, I),
   MAKE_INTRINSIC_I("_pgsclp", _pgsclp, V),
   MAKE_INTRINSIC_4("_pgscr", _pgscr, V, I, D, D, D),
   MAKE_INTRINSIC_2("_pgscrl", _pgscrl, V, D, D),
   MAKE_INTRINSIC_2("_pgscrn", _pgscrn, I, I, S),
   MAKE_INTRINSIC_I("_pgsfs", _pgsfs, V),
   MAKE_INTRINSIC_4("_pgshls", _pgshls, V, I, D, D, D),
   MAKE_INTRINSIC_3("_pgshs", _pgshs, V, D, D, D),
   MAKE_INTRINSIC_I("_pgsitf", _pgsitf, V),
   MAKE_INTRINSIC_I("_pgslct", _pgslct, V),
   MAKE_INTRINSIC_I("_pgsls", _pgsls, V),
   MAKE_INTRINSIC_I("_pgslw", _pgslw, V),
   MAKE_INTRINSIC_I("_pgstbg", _pgstbg, V),
   MAKE_INTRINSIC_2("_pgsubp", _pgsubp, V, I, I),
   MAKE_INTRINSIC_4("_pgsvp", _pgsvp, V, D, D, D, D),
   MAKE_INTRINSIC_4("_pgswin", _pgswin, V, D, D, D, D),
   MAKE_INTRINSIC_6("_pgtbox", _pgtbox, V, S, D, I, S, D, I),
   MAKE_INTRINSIC_3("_pgtext", _pgtext, V, D, D, S),
   /* MAKE_INTRINSIC_0("_pgtick", _pgtick), */
   MAKE_INTRINSIC_0("_pgunsa", _pgunsa, V),
   MAKE_INTRINSIC_0("_pgupdt", _pgupdt, V),
   MAKE_INTRINSIC_4("_pgvsiz", _pgvsiz, V, D, D, D, D),
   MAKE_INTRINSIC_0("_pgvstd", _pgvstd, V),
   /* MAKE_INTRINSIC_0("_pgwedg", _pgwedg), */
   MAKE_INTRINSIC_4("_pgwnad", _pgwnad, V, D, D, D, D),
   
   SLANG_END_INTRIN_FUN_TABLE
};

#undef V
#undef I
#undef S
#undef D

static SLang_Intrin_Var_Type Module_Variables [] =
{
   MAKE_VARIABLE("_pgplot_module_version_string", &Module_Version_String, SLANG_STRING_TYPE, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_IConstant_Type Module_IConstants [] =
{
   MAKE_ICONSTANT("_pgplot_module_version", MODULE_VERSION_NUMBER),
   SLANG_END_ICONST_TABLE
};

static int setenv_pgplot_dir (void)
{      
   static char *pgplot_dir;
    
   if (NULL != getenv ("PGPLOT_DIR"))
     return 0;
    
   pgplot_dir = "PGPLOT_DIR=" PGPLOTDIR;
   if (-1 == putenv (pgplot_dir))
     {
        fprintf (stderr, "Failed setting PGPLOT_DIR environment variable: %s\n", pgplot_dir);
        return -1;
     }
    
   return 0;
}

int init_pgplot_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns;
   
   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return -1;


   if ((-1 == SLns_add_intrin_fun_table (ns, Module_Intrinsics, "__PGPLOT__"))
       || (-1 == SLns_add_iconstant_table (ns, Module_IConstants, NULL))
       || (-1 == SLns_add_intrin_var_table (ns, Module_Variables, NULL)))
      return -1;
   
   return setenv_pgplot_dir ();
}

void deinit_pgplot_module (void)
{
}

#ifdef FC_DUMMY_MAIN
# ifdef __cplusplus
extern "C"
# endif
int FC_DUMMY_MAIN (void);
int FC_DUMMY_MAIN (void)
{
   SLang_vmessage ("*** pgplot-module: FC_DUMMY_MAIN called\n");
   return 1;
}
#endif

