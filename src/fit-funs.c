/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2010 Massachusetts Institute of Technology

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

/* $Id: fit-funs.c,v 1.2 2004/02/09 11:14:20 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#include <slang.h>
#include "util.h"
#include "isis.h"
#include "fit.h"
#include "plot.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

/* integrated over bin width */
static int poly_b (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   int i;
   double xh, xl, dx, xm, xt;

   (void) npar;

   for (i=0; i < g->n_notice; i++)
     {
        int n = g->notice_list[i];
        xh = g->bin_hi[n];
        xl = g->bin_lo[n];
        xm = (xh + xl) / 2.0;
        xt = (xh*xh + xh*xl + xl*xl) / 3.0;
        dx =  xh - xl;
        val[i] = dx * (par[0] + xm * par[1] + xt * par[2]);
     }

   return 0;
}

/*}}}*/

static int poly_c (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   int i;
   (void) npar;

   for (i=0; i < g->npts; i++)
     {
        double x = g->x[i];
        val[i] = par[0] + x * (par[1] + x * par[2]);
     }

   return 0;
}

/*}}}*/

static int poly_p (unsigned int fun_id, double *par, double *par_min, double *par_max,  /*{{{*/
                   unsigned int npar)
{
   float c[3] = {0.0, 0.0, 0.0};
   float x[3] = {0.0, 0.0, 0.0};
   float y[3] = {0.0, 0.0, 0.0};
   float ca, cb, cc;
   int n = 0;
   int ret = 0;

   (void) npar;

   fprintf (stdout, "Initializing Polynomial[%u]:\n", fun_id);
   fprintf (stdout, "   Please click <= 3 points to define the polynomial\n");
   fprintf (stdout, "   Press 'x' to reposition the last point\n");
   fprintf (stdout, "   Press 'q' to quit\n");
   fflush (stdout);

   while (n < 3)
     {
        int reject;
        do
          {
             char ch = '\0';

             reject = 0;

             if (-1 == Plot_read_crosshair_cursor (&x[n], &y[n], &ch))
               {
                  isis_vmesg (WARN, I_FAILED, __FILE__, __LINE__, "reading cursor position");
                  ret = -1;
                  goto done;
               }

             if (-1 == change_xunits_to_angstrom (&x[n], &y[n]))
               {
                  isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "converting to wavelength units");
                  ret = -1;
                  goto done;
               }

             if (ch == 'q' || ch == 'Q')
               goto done;
             else if (ch == 'x' || ch == 'X')
               {
                  n = MAX(0, n-1);
                  reject = 1;
                  fprintf (stdout, "*** repositioning n=%d\n", n+1);
                  fflush (stdout);
               }
             else
               {                          /* reject duplicate x-values */
                  int j, k;
                  for (k=1; k < n+1; k++)
                    {
                       for (j=0; j < k; j++)
                         {
                            float dx = fabs(x[j]- x[k]);
                            float eps = 10*FLT_EPSILON * fabs(x[j] + x[k]);
                            if (dx > eps)
                              continue;
                            fprintf (stdout, "duplicate x value rejected ---  re-enter   (to quit, type 'q')\n");
                            fflush (stdout);
                            reject = 1;
                            goto dup;
                         }
                    }
               }

             dup:
             ;

          } while (reject);

        fprintf (stdout, "n=%d  x=%14.6e  y=%14.6e\n", n+1, x[n], y[n]);
        fflush (stdout);

        n++;
     }

   done:

   if (ret)
     return ret;

   switch (n)
     {
      case 1:                                        /* constant */
        c[0] = y[0];
        break;
      case 2:                                        /* linear */
        c[1] = (y[1] - y[0]) / (x[1] - x[0]);
        c[0] = y[0] - x[0] * c[1];
        break;
      case 3:                                        /* quadratic */
        ca = y[0] / (x[0] - x[1]) / (x[0] - x[2]);
        cb = y[1] / (x[1] - x[0]) / (x[1] - x[2]);
        cc = y[2] / (x[2] - x[0]) / (x[2] - x[1]);
        c[2] = ca + cb + cc;
        c[1] = - (ca * (x[1] + x[2]) +
                  cb * (x[0] + x[2]) +
                  cc * (x[0] + x[1]));
        c[0] = (ca * x[1] * x[2] +
                cb * x[0] * x[2] +
                cc * x[0] * x[1]);
        break;
      default:
        break;
     }

   par[0] = (double) c[0];
   par[1] = (double) c[1];
   par[2] = (double) c[2];

   par_min[0] = par_min[1] = par_min[2] = -DBL_MAX;
   par_max[0] = par_max[1] = par_max[2] = DBL_MAX;

   return ret;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(poly,p,options) /*{{{*/
{
   static char *names[] = {"a0", "a1", "a2", NULL};
   static double default_min []= {0.0, 0.0, 0.0};
   static double default_value []= {1.0, 1.0, 1.0};
   static double default_max []= {0.0, 0.0, 0.0};
   static unsigned int norm_indexes [] = {0,1,2};

   (void) options;

   if (p == NULL)
     return -1;

   p->binned = poly_b;
   p->function_exit = NULL;
   p->unbinned = poly_c;
   p->init_params_from_screen = poly_p;
   p->norm_indexes = norm_indexes;
   p->num_norms = 3;
   p->parameter_names = names;
   p->default_min = default_min;
   p->default_max = default_max;
   p->default_value = default_value;
   p->num_parameters = 3;

   return 0;
}

/*}}}*/

static void delta_b (double *val, Isis_Hist_t *g, double x0, double area) /*{{{*/
{
   int i;

   for (i = 0; i < g->n_notice; i++)
     {
        if (g->bin_lo[i] <= x0 && x0 < g->bin_hi[i])
          val[i] = area;
        else val[i] = 0.0;
     }
}

/*}}}*/

/* integrated over bin width */
static int lorentz_b (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double area = par[0];
   double x0   = par[1];
   double hfw  = par[2] / 2.0;           /* fwhm/2 */
   int i;

   (void) npar;

   if (hfw < 0.0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "lorentz: FWHM < 0");
        return -1;
     }

   if (hfw > 0.0)
     {
        for (i=0; i < g->n_notice; i++)
          {
             double dxh, dxl;
             int n = g->notice_list[i];

             dxh = (g->bin_hi[n] - x0) / hfw;
             dxl = (g->bin_lo[n] - x0) / hfw;

             val[i] = area * (atan (dxh) - atan (dxl)) / PI;
          }
     }
   else delta_b (val, g, x0, area);

   return 0;
}
/*}}}*/

static int lorentz_c (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   double area = par[0];
   double x0   = par[1];
   double hfw  = par[2] / 2.0;           /* fwhm/2 */
   double coef = area * hfw / PI;
   int i;

   (void) npar;

   if (hfw < 0.0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "lorentz: FWHM < 0");
        return -1;
     }

   for (i=0; i < g->npts; i++)
     {
        double x = g->x[i] - x0;
        val[i] = coef / ( x*x + hfw*hfw );
     }

   return 0;
}

/*}}}*/

static int lorentz_p (unsigned int fun_id, double *par, double *par_min, double *par_max, unsigned int npar) /*{{{*/
{
   float bx[2] = { 0.0, 0.0 };
   float by[2] = { 0.0, 0.0 };

   (void) npar;

   fprintf (stdout, "Initializing Lorentz[%u]:\n", fun_id);
   fprintf (stdout, "   Please draw a box to indicate the center, peak and fwhm\n");

   if (-1 == Plot_read_box_corners_from_cursor (bx,by))
     {
        isis_vmesg (WARN, I_FAILED, __FILE__, __LINE__, "reading cursor position");
        return -1;
     }

   if (-1 == change_xunits_to_angstrom (&bx[0], &by[0])
       || -1 == change_xunits_to_angstrom (&bx[1], &by[1]))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "converting to wavelength units");
        return -1;
     }
   else
     {
        float mx = MAX(bx[0], bx[1]);
        float mn = MIN(bx[0], bx[1]);
        bx[0] = mn;
        bx[1] = mx;
     }

   /*
    * Now estimate parameter values
    */

   par[0] = (double) (bx[1] - bx[0]) * 2.0 * by[1];      /* area */
   par[1] = (double) (bx[1] + bx[0]) / 2.0;              /* center */
   par[2] = (double) (bx[1] - bx[0]);                    /* fwhm */

   par_min[0] = par_min[1] = par_min[2] = -DBL_MAX;
   par_max[0] = par_max[1] = par_max[2] = DBL_MAX;

   return 0;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(Lorentz,p,options) /*{{{*/
{
   static char *names[] = {"area", "center", "fwhm", NULL};
   static char *units[] = {"photons/s/cm^2", "A", "A", NULL};
   static double default_max []= {0.0, 0.0, 1.0};
   static double default_value []= {1.0, 12.0, 0.025};
   static double default_min []= {0.0, 0.0, 1.e-6};
   static unsigned int norm_indexes[] = {0};

   (void) options;

   if (p == NULL)
     return -1;

   p->binned = lorentz_b;
   p->function_exit = NULL;
   p->unbinned = lorentz_c;
   p->init_params_from_screen = lorentz_p;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;
   p->parameter_names = names;
   p->parameter_units = units;
   p->default_min = default_min;
   p->default_max = default_max;
   p->default_value = default_value;
   p->num_parameters = 3;

   return 0;
}

/*}}}*/

/* integrated over bin width */
static int gauss_b (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double area = par[0];
   double x0   = par[1];
   double sigma= par[2];
   int i;

   (void) npar;

   if (sigma < 0.0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "gauss: sigma < 0");
        return -1;
     }

   if (sigma > 0.0)
     {
        for (i=0; i < g->n_notice; i++)
          {
             double dxh, dxl;
             int n = g->notice_list[i];

             dxh = (g->bin_hi[n] - x0) / sigma;
             dxl = (g->bin_lo[n] - x0) / sigma;

             val[i] = area * (isis_gpf (dxh) - isis_gpf (dxl));
          }
     }
   else delta_b (val, g, x0, area);

   return 0;
}

/*}}}*/

static int gauss_c (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   double area  = par[0];
   double x0    = par[1];
   double sigma = par[2];
   double coef = area / (sigma * sqrt(2*PI));
   int i;

   (void) npar;

   if (sigma <= 0.0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "gauss: sigma <= 0");
        return -1;
     }

   for (i=0; i < g->npts; i++)
     {
        double x = (g->x[i] - x0) / sigma;
        val[i] = coef * exp (-0.5 * x * x);
     }

   return 0;
}

/*}}}*/

static int gauss_p (unsigned int fun_id, double *par, double *par_min, double *par_max, unsigned int npar) /*{{{*/
{
   double factor = 2.354820044;     /* 2*sqrt(2*ln(2)) */
   float bx[2] = { 0.0, 0.0 };
   float by[2] = { 0.0, 0.0 };

   (void) npar;

   fprintf (stdout, "Initializing Gauss[%u]:\n", fun_id);
   fprintf (stdout, "   Please draw a box to indicate the center, peak and fwhm\n");

   if (-1 == Plot_read_box_corners_from_cursor (bx,by))
     {
        isis_vmesg (WARN, I_FAILED, __FILE__, __LINE__, "reading cursor position");
        return -1;
     }

   if (-1 == change_xunits_to_angstrom (&bx[0], &by[0])
       || -1 == change_xunits_to_angstrom (&bx[1], &by[1]))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "converting to wavelength units");
        return -1;
     }
   else
     {
        float mx = MAX(bx[0], bx[1]);
        float mn = MIN(bx[0], bx[1]);
        bx[0] = mn;
        bx[1] = mx;
     }

   /*
    * Now estimate parameter values
    */

   par[0] = (double) (bx[1] - bx[0]) * by[1];             /* area */
   par[1] = (double) (bx[1] + bx[0]) / 2.0;               /* center */
   par[2] = (double) (bx[1] - bx[0]) / factor;            /* sigma */

   par_min[0] = par_min[1] = par_min[2] = -DBL_MAX;
   par_max[0] = par_max[1] = par_max[2] = DBL_MAX;

   return 0;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(gauss,p,options) /*{{{*/
{
   static char *names[] = {"area", "center", "sigma", NULL};
   static char *units[] = {"photons/s/cm^2", "A", "A", NULL};
   static double default_max []= {0.0, 0.0, 1.0};
   static double default_value []= {1.0, 12.0, 0.025};
   static double default_min []= {0.0, 0.0, 1.e-6};
   static unsigned int norm_indexes[] = {0};

   (void) options;

   if (p == NULL)
     return -1;

   p->binned = gauss_b;
   p->function_exit = NULL;
   p->unbinned = gauss_c;
   p->init_params_from_screen = gauss_p;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;
   p->parameter_names = names;
   p->parameter_units = units;
   p->default_min = default_min;
   p->default_max = default_max;
   p->default_value = default_value;
   p->num_parameters = 3;

   return 0;
}

/*}}}*/

/* integrated over bin width */
static int egauss_b (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double area = par[0];
   double e0   = par[1];
   double sigma= par[2];
   int i;

   (void) npar;

   if (sigma < 0.0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "egauss: sigma < 0");
        return -1;
     }

   if (sigma > 0.0)
     {
        for (i=0; i < g->n_notice; i++)
          {
             double elo, ehi, dxh, dxl;
             int n = g->notice_list[i];

             elo = KEV_ANGSTROM / g->bin_hi[n];
             ehi = KEV_ANGSTROM / g->bin_lo[n];

             dxh = (ehi - e0) / sigma;
             dxl = (elo - e0) / sigma;
             val[i] = area * (isis_gpf (dxh) - isis_gpf (dxl));
          }
     }
   else delta_b (val, g, KEV_ANGSTROM/e0, area);

   return 0;
}

/*}}}*/

static int egauss_c (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   double area  = par[0];
   double e0    = par[1];
   double sigma = par[2];
   double coef = area / (sigma * sqrt(2*PI));
   int i;

   (void) npar;

   if (sigma <= 0.0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "egauss: sigma <= 0");
        return -1;
     }

   for (i=0; i < g->npts; i++)
     {
        double e, x;
        e = KEV_ANGSTROM / g->x[i];
        x = (e - e0) / sigma;
        val[i] = coef * exp (-0.5 * x * x);
     }

   return 0;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(egauss,p,options) /*{{{*/
{
   static char *names[] = {"area", "center", "sigma", NULL};
   static char *units[] = {"photons/s/cm^2", "keV", "keV", NULL};
   static double default_max[]   = {0.0, 0.0, 1.0};
   static double default_value[] = {1.0, 1.0, 0.002};
   static double default_min[]   = {0.0, 0.0, 1.e-6};
   static unsigned int norm_indexes[] = {0};

   (void) options;

   if (p == NULL)
     return -1;

   p->binned = egauss_b;
   p->function_exit = NULL;
   p->unbinned = egauss_c;
   p->init_params_from_screen = NULL;;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;
   p->parameter_names = names;
   p->parameter_units = units;
   p->default_min = default_min;
   p->default_max = default_max;
   p->default_value = default_value;
   p->num_parameters = 3;

   return 0;
}

/* integrated over bin width */
static int powr_b (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double norm = par[0];
   double alph = par[1];
   int i, n;

   (void) npar;

   if (fabs(alph+1.0) > 1.e4 * DBL_EPSILON)
     {
        double oma = 1.0 + alph;

        for (i=0; i < g->n_notice; i++)
          {
             n = g->notice_list[i];
             val[i] = (norm / oma) * (pow (KEV_ANGSTROM / g->bin_lo[n], oma)
                                      - pow (KEV_ANGSTROM / g->bin_hi[n], oma));
          }
     }
   else
     {
        for (i=0; i < g->n_notice; i++)
          {
             n = g->notice_list[i];
             val[i] = norm * log (g->bin_hi[n] / g->bin_lo[n]);
          }
     }

   return 0;
}

/*}}}*/

static int powr_c (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   double norm = par[0];
   double alph = par[1];
   int i;

   (void) npar;

   for (i=0; i < g->npts; i++)
     {
        val[i] = norm / pow (KEV_ANGSTROM / g->x[i], alph);
     }

   return 0;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(Powerlaw,p,options) /*{{{*/
{
   static char *names[] = {"norm", "alpha", NULL};
   static double default_max []= {1.e8,   0.0};
   static double default_value []= {1.0, -1.5};
   static double default_min []= {1.e-8, -10.0};
   static unsigned int norm_indexes[] = {0};

   (void) options;

   if (p == NULL)
     return -1;

   p->binned = powr_b;
   p->function_exit = NULL;
   p->unbinned = powr_c;
   p->init_params_from_screen = NULL;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;
   p->parameter_names = names;
   p->default_min = default_min;
   p->default_max = default_max;
   p->default_value = default_value;
   p->num_parameters = 2;

   return 0;
}

/*}}}*/

/* integrated over bin width */
static double planck (double e, double t) /*{{{*/
{
   double min_exp = 1.e-13;
   double max_exp = 500.0;
   double x = e / t;

   if (x < min_exp)
     return e * t;
   else if (x > max_exp)
     return 0.0;

   return e * e / (exp (x) - 1.0);
}

/*}}}*/

static int bbody_b (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double t, norm;
   int i;

   if (par[1] <= 0.0)
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "blackbody: temperature <= 0");
        return -1;
     }

   t    = par[1];
   norm = par[0] * 8.0525 / pow (t, 4.0);
   (void) npar;

   for (i=0; i < g->n_notice; i++)
     {
        double elo, ehi;
        int n = g->notice_list[i];
        elo = KEV_ANGSTROM / g->bin_hi[n];
        ehi = KEV_ANGSTROM / g->bin_lo[n];
        val[i] = (0.5 * (planck (ehi, t) + planck (elo, t))
                  * norm * (ehi - elo));
     }

   return 0;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(blackbody,p,options) /*{{{*/
{
   static char *names[] = {"norm", "kT", NULL};
   static char *units[] = {"", "keV", NULL};
   static double default_max[]   = {10.0, 10.0};
   static double default_value[] = {1.0,  1.0};
   static double default_min[]   = {0.0,  0.1};
   static unsigned int norm_indexes[] = {0};

   (void) options;

   if (p == NULL)
     return -1;

   p->binned = bbody_b;
   p->function_exit = NULL;
   p->unbinned = NULL;
   p->init_params_from_screen = NULL;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;
   p->parameter_names = names;
   p->parameter_units = units;
   p->default_min = default_min;
   p->default_max = default_max;
   p->default_value = default_value;
   p->num_parameters = 2;

   return 0;
}

/*}}}*/

/*}}}*/

#include "voigt.c"

typedef struct
{
   char *name;
   Isis_User_Source_Init_Fun_t *init;
}
Static_Function_Type;

#define STATIC_FUNCTION(name) {#name, Isis_##name##_function}

static Static_Function_Type Static_Function_Table [] =
{
  STATIC_FUNCTION(poly),
  STATIC_FUNCTION(Lorentz),
  STATIC_FUNCTION(gauss),
  STATIC_FUNCTION(egauss),
  STATIC_FUNCTION(Powerlaw),
  STATIC_FUNCTION(blackbody),
  STATIC_FUNCTION(voigt),
  {NULL,NULL}
};

int Fit_Append_Builtin_Functions (void)
{
   Static_Function_Type *t = Static_Function_Table;
   Static_Function_Type *u;

   for (u = t; u->name != NULL; u++)
     {
        if (-1 == Isis_Add_Static_Fun (u->init, u->name))
          return -1;
     }

   return 0;
}
