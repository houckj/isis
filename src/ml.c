/* -*- mode: C; mode: fold -*- */

/*
    This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2019 Massachusetts Institute of Technology

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

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "isis.h"
#include "util.h"
#include "errors.h"

/* ML statistic */

static int ml_function (Isis_Fit_Statistic_Type *st, /*{{{*/
                          double *y, double *fx, double *w, unsigned int npts,
                          double *vec, double *stat)
{
   double sum, s, *val;
   unsigned int i;

   (void) w;
   (void) st;

   if (NULL == (val = (double *) ISIS_MALLOC (npts * sizeof(double))))
     return -1;

   for (i = 0; i < npts; i++)
     {
        double fxi = fx[i];
        double log_fxi = (fxi > 0) ? log(fxi) : (double) DBL_MIN_10_EXP;

        s = 2 * (fxi + lgamma (y[i] + 1) - y[i] * log_fxi);

        val[i] = s;
        vec[i] = s;
     }

   sum = isis_kahan_sum (val, npts);
   ISIS_FREE(val);

   if (0 == isfinite (sum))
     sum = DBL_MAX;

   *stat = sum;

   return 0;
}

/*}}}*/

static int ml_report (Isis_Fit_Statistic_Type *st, void *pfp, double stat, /*{{{*/
                        unsigned int npts, unsigned int nvpars)
{
   FILE *fp = (FILE *)pfp;
   (void) st;

   fprintf (fp, "           ML Statistic = %0.7g\n", stat);
   if (nvpars < npts)
     fprintf (fp, "   Reduced ML Statistic = %0.7g\n", stat/(npts - nvpars));

   return 0;
}

/*}}}*/

static void deallocate (Isis_Fit_Statistic_Type *s)
{
   ISIS_FREE (s->symbol);
   ISIS_FREE (s->option_string);
}

ISIS_FIT_STATISTIC_METHOD (ml)
{
   Isis_Fit_Statistic_Type *s;

   if (NULL == (s = (Isis_Fit_Statistic_Type *) ISIS_MALLOC (sizeof(Isis_Fit_Statistic_Type))))
     return NULL;
   memset ((char *)s, 0, sizeof (*s));

   s->compute_statistic = ml_function;
   s->report = ml_report;
   s->delta_is_chisqr_distributed = 1;
   s->uses_opt_data = 0;
   s->opt_data = NULL;
   s->deallocate = deallocate;
   if (NULL == (s->symbol = isis_make_string ("ML")))
     {
        deallocate (s);
        ISIS_FREE (s);
        return NULL;
     }

   s->option_string = isis_make_string ("ml");
   if (s->option_string == NULL)
     {
        deallocate (s);
        ISIS_FREE (s);
        return NULL;
     }

   return s;
}

