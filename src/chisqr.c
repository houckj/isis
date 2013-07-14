/* -*- mode: C; mode: fold -*- */

/*
    This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2012 Massachusetts Institute of Technology

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

/*$Id: chisqr.c,v 1.6 2004/02/09 11:14:17 houck Exp $*/

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#define ISIS_FIT_STATISTIC_PRIVATE_DATA \
   unsigned int sigma;      /* select chi-square variant */

#include "isis.h"
#include "util.h"

enum
{
     DATA_VARIANCE,
     GEHRELS_VARIANCE,
     MODEL_VARIANCE,
     LEAST_SQUARES
};

typedef struct Chisqr_Type
{
   char *name;
   unsigned int code;
}
Chisqr_Type;

static Chisqr_Type Chisqr_Types[] =
{
   {"data",    DATA_VARIANCE},
   {"gehrels", GEHRELS_VARIANCE},
   {"model",   MODEL_VARIANCE},
   {"lsq",     LEAST_SQUARES},
   {NULL,      DATA_VARIANCE}
};

/* Chi-square statistic */

static int chisqr_function (Isis_Fit_Statistic_Type *st,  /*{{{*/
                            double *y, double *fx, double *w,
                            unsigned int npts, double *vec, double *stat)
{
   double sum, dy, wt;
   unsigned int i, bad;

   bad = 0;
   sum = 0.0;

   *stat = -1.0;

   switch (st->sigma)
     {
      case DATA_VARIANCE:
        for (i = 0; i < npts; i++)
          {
             dy = y[i] - fx[i];
             /* sum += dy*dy * w[i]; */
             vec[i] = dy * sqrt(w[i]);
          }
        break;

      case GEHRELS_VARIANCE:
        for (i = 0; i < npts; i++)
          {
             dy = y[i] - fx[i];
             if (y[i] >= 0.0) dy /= 1.0 + sqrt(y[i] + 0.75);
             else bad++;
             /* sum += dy*dy; */
             vec[i] = dy;
          }
        if (bad && (st->message_string == NULL))
          {
             st->message_string =
               "*** Warning: chisqr used sigma=1.0 for some points with negative data values";
          }
        break;

      case MODEL_VARIANCE:
        for (i = 0; i < npts; i++)
          {
             dy = y[i] - fx[i];
             wt = 1.0;
             if (fx[i] != 0.0) wt = 1.0 / fabs(fx[i]);
             else bad++;
             /* sum += dy*dy * wt; */
             vec[i] = dy * sqrt(wt);
          }
        if (bad && (st->message_string == NULL))
          {
             st->message_string =
               "*** Warning: chisqr used sigma=1.0 for some points with model=0";
          }
        break;

      case LEAST_SQUARES:
        for (i = 0; i < npts; i++)
          {
             dy = y[i] - fx[i];
             /* sum += dy*dy; */
             vec[i] = dy;
          }
        break;

      default:
        fprintf (stderr, "*** chisqr_function:  internal error\n");
        return -1;
     }

   sum = isis_kahan_sum_squares (vec, npts);

   if (0 == isfinite(sum))
     sum = DBL_MAX;

   *stat = sum;

   return 0;
}

/*}}}*/

static int chisqr_report (Isis_Fit_Statistic_Type *st, void *pfp, double stat, /*{{{*/
                          unsigned int npts, unsigned int nvpars)
{
   FILE *fp = (FILE *)pfp;

   fprintf (fp, "           Chi-square = %0.7g", stat);
   if (st->sigma == DATA_VARIANCE)
     fputc ('\n', fp);
   else
     {
        Chisqr_Type *t = Chisqr_Types;
        char *sigma_name = NULL;
        while (t->name != NULL)
          {
             if (t->code == st->sigma)
               {
                  sigma_name = t->name;
                  break;
               }
             t++;
          }
        if (sigma_name != NULL)
          fprintf (fp, "    (sigma = %s)\n", sigma_name);
        else
          {
             fprintf (fp, "*** chisqr_report:  internal error\n");
             return -1;
          }
     }

   if (nvpars < npts)
     fprintf (fp, "   Reduced chi-square = %0.7g\n", stat/(npts - nvpars));

   return 0;
}

/*}}}*/

static int handle_sigma_option (char *subsystem, char *optname, char *value, void *clientdata) /*{{{*/
{
   Isis_Fit_Statistic_Type *s;
   Chisqr_Type *t = Chisqr_Types;
   s = (Isis_Fit_Statistic_Type *) clientdata;

   (void) subsystem;

   while (t->name != NULL)
     {
        if (0 == isis_strcasecmp (value, t->name))
          {
             s->sigma = t->code;
             return isis_update_option_string (&s->option_string, optname, value);
          }
        t++;
     }

   return -1;
}

/*}}}*/

static Isis_Option_Table_Type Option_Table [] =
{
     {"sigma", handle_sigma_option, ISIS_OPT_REQUIRES_VALUE, "data", "(data | gehrels | model | lsq)"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

static int set_options (Isis_Fit_Statistic_Type *s, Isis_Option_Type *opts)
{
   return isis_process_options (opts, Option_Table, (void *)s, 1);
}

static void deallocate (Isis_Fit_Statistic_Type *s)
{
   ISIS_FREE (s->symbol);
   ISIS_FREE (s->option_string);
}

ISIS_FIT_STATISTIC_METHOD (chisqr)
{
   Isis_Fit_Statistic_Type *s;

   if (NULL == (s = (Isis_Fit_Statistic_Type *) ISIS_MALLOC (sizeof(Isis_Fit_Statistic_Type))))
     return NULL;
   memset ((char *)s, 0, sizeof (*s));

   s->compute_statistic = chisqr_function;
   s->report = chisqr_report;
   s->delta_is_chisqr_distributed = 1;
   s->uses_opt_data = 0;
   s->opt_data = NULL;
   s->set_options = set_options;
   s->deallocate = deallocate;
   s->sigma = DATA_VARIANCE;
   if (NULL == (s->symbol = isis_make_string ("\\gD\\gx")))
     {
        deallocate (s);
        ISIS_FREE (s);
        return NULL;
     }

   s->option_string = isis_make_default_option_string ("chisqr", Option_Table);
   if (s->option_string == NULL)
     {
        deallocate (s);
        ISIS_FREE (s);
        return NULL;
     }

   return s;
}

