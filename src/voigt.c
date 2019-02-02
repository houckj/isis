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

#include <stdlib.h>
#include <math.h>
#include <slang.h>

#include "isis.h"

extern int Isis_Voigt_Is_Normalized;

#ifndef PI
#define      PI 3.14159265358979323846264338328
#endif

#ifndef SQRT_PI
#define SQRT_PI 1.77245385090551602729816748334
#endif

#define FOURPI  (4.0*PI)
#define C_KMS   299792.458

#ifndef KEV_ANGSTROM
#define KEV_ANGSTROM 12.3984185   /* hc */
#endif

/*{{{ w(z) = Faddeeva Function */

/*
 * Relation between the Voigt function and the Faddeeva function w(z):
 *
 *    see p. 276-279 Mihalas, Stellar Atmospheres
 *        p. 302, Abramowitz & Stegun, eq. 7.4.13
 *
 *  Natural resonance profile of emitted line:
 *  =========================================
 *
 *                                 (Gamma/(4*pi))
 *    Lorentz (nu) = (1/pi) -----------------------------
 *                          (nu-nu0)^2 + (Gamma/(4*pi))^2
 *
 *    where Gamma = resonance profile FWHM in frequency units (Hz)
 *
 *    In general, for a transition from some upper state U to
 *    some lower state L,
 *
 *    Gamma(U,L) = Sum[ i<U; A(U,i) ] + Sum[ i<L; A(L,i) ]
 *
 *           where A(k,i) = Einstein transition probability from
 *                          upper state k to lower state i (in sec^-1).
 *
 *  Maxwellian speed distribution along the line of sight through plasma:
 *  ====================================================================
 *
 *    Gauss (s) = (s0 * sqrt(pi))^(-1)  exp (-(s/s0)^2)
 *
 *                  where  s0 = sqrt (2*k*T/m)
 *                         m = ion mass
 *
 *  Convolution of Lorentzian and Maxwellian/Gaussian gives:
 *  =======================================================
 *
 *  Voigt function = H(x,y) = (1/Pi) Integral[ t, -Inf, Inf, f(t; x,y) ]
 *
 *                                       y * exp(-t^2)
 *                  where   f(t; x,y) = ---------------
 *                                       (x-t)^2 + y^2
 *
 *                          x = (nu - nu0) / Doppler_width
 *                          y = (Gamma/(4*pi)) / Doppler_width
 *
 *                          Doppler_width = s0 * nu0/c
 *
 *  A&S eq. 7.4.13 --->  H(x,y) = Re[w(z)];  w(z) = Faddeeva function
 *
 *                          z = x + i*y
 *                          i = sqrt(-1)
 *                          w (z) = exp (-z^2) * erfc (-i*z)
 */

/*
c      algorithm 680, collected algorithms from acm.
c      this work published in transactions on mathematical software,
c      vol. 16, no. 1, pp. 47.
c
c  given a complex number z = (xi,yi), this subroutine computes
c  the value of the Faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
c  where erfc is the complex complementary error-function and i
c  means sqrt(-1).
c  the accuracy of the algorithm for z in the 1st and 2nd quadrant
c  is 14 significant digits; in the 3rd and 4th it is 13 significant
c  digits outside a circular region with radius 0.126 around a zero
c  of the function.
c  all real variables in the program are double precision.
c
c
c  the code contains a few compiler-dependent parameters :
c     rmaxreal = the maximum value of rmaxreal equals the root of
c                rmax = the largest number which can still be
c                implemented on the computer in double precision
c                floating-point arithmetic
c     rmaxexp  = ln(rmax) - ln(2)
c     rmaxgoni = the largest possible argument of a double precision
c                goniometric function (dcos, dsin, ...)
c  the reason why these parameters are needed as they are defined will
c  be explained in the code by means of comments
c
c
c  parameter list
c     xi     = real      part of z
c     yi     = imaginary part of z
c     u      = real      part of w(z)
c     v      = imaginary part of w(z)
c     flag   = an error flag indicating whether overflow will
c              occur or not; type logical;
c              the values of this variable have the following
c              meaning :
c              flag=.false. : no error condition
c              flag=.true.  : overflow will occur, the routine
c                             becomes inactive
c  xi, yi      are the input-parameters
c  u, v, flag  are the output-parameters
c
c  furthermore the parameter factor equals 2/sqrt(pi)
c
c  the routine is not underflow-protected but any variable can be
c  put to 0 upon underflow;
c
c  reference - gpm poppe, cmj wijers; more efficient computation of
c  the complex error-function, acm trans. math. software.
c
*/

static int idnint(double x) /*{{{*/
{
   return (int)(x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x));
}

/*}}}*/

static double factor = 1.12837916709551257388;
static double rmaxreal = 0.5e+154;
static double rmaxexp  = 708.503061461606;
static double rmaxgoni = 3.53711887601422e+15;

static int wofz (double xi, double yi, double *u, double *v) /*{{{*/
{
   int a, b, n, i, j, kapn, nu, np1;
   double xabs, yabs, x, y, qrho, xabsq, xquad, yquad;
   double xsum, ysum, xaux, u1, v1, daux, u2, v2, h, h2;
   double qlambda, rx, ry, sx, sy, tx, ty, c, w1;

   u2 = v2 = h2 = qlambda = 0.0;

   xabs = fabs(xi);
   yabs = fabs(yi);
   x    = xabs/6.3;
   y    = yabs/4.4;

/*  the following if-statement protects
 *  qrho = (x**2 + y**2) against overflow
 */

   if ((xabs > rmaxreal) || (yabs > rmaxreal))
     return -1;

   qrho = x*x + y*y;

   xabsq = xabs*xabs;
   xquad = xabsq - yabs*yabs;
   yquad = 2*xabs*yabs;

/*
 * if (qrho.lt.0.085264d0) then the faddeeva-function is evaluated
 * using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
 * n is the minimum number of terms needed to obtain the required
 * accuracy
 */

   a = qrho < 0.085264;

   if (a)
     {
        qrho  = (1-0.85*y)*sqrt(qrho);
        n     = idnint(6 + 72*qrho);
        j     = 2*n+1;
        xsum  = 1.0/j;
        ysum  = 0.0;
        for (i=n; i >= 1; i--)
          {
             j    = j - 2;
             xaux = (xsum*xquad - ysum*yquad)/i;
             ysum = (xsum*yquad + ysum*xquad)/i;
             xsum = xaux + 1.0/j;
          }

        u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0;
        v1   =  factor*(xsum*xabs - ysum*yabs);
        daux =  exp(-xquad);
        u2   =  daux*cos(yquad);
        v2   = -daux*sin(yquad);

        *u    = u1*u2 - v1*v2;
        *v    = u1*v2 + v1*u2;
     }
   else
     {
        /* if (qrho.gt.1.o) then w(z) is evaluated using the laplace
         * continued fraction
         * nu is the minimum number of terms needed to obtain the required
         * accuracy
         *
         * if ((qrho.gt.0.085264d0).and.(qrho.lt.1.0)) then w(z) is evaluated
         * by a truncated taylor expansion, where the laplace continued fraction
         * is used to calculate the derivatives of w(z)
         * kapn is the minimum number of terms in the taylor expansion needed
         * to obtain the required accuracy
         * nu is the minimum number of terms of the continued fraction needed
         * to calculate the derivatives with the required accuracy
         */

        if (qrho > 1.0)
          {
             h    = 0.0;
             kapn = 0;
             qrho = sqrt(qrho);
             nu = (int) (1442 / (qrho * 26 + 77) + 3);
          }
        else
          {
             qrho = (1-y)*sqrt(1-qrho);
             h    = 1.88*qrho;
             h2   = 2*h;
             kapn = idnint(7  + 34*qrho);
             nu   = idnint(16 + 26*qrho);
          }

        b = h > 0.0;

        if (b)
          qlambda = pow(h2, (double) kapn);

        rx = 0.0;
        ry = 0.0;
        sx = 0.0;
        sy = 0.0;

        for (n=nu; n >= 0; n--)
          {
             np1 = n + 1;
             tx  = yabs + h + np1*rx;
             ty  = xabs - np1*ry;
             c   = 0.5/(tx*tx + ty*ty);
             rx  = c*tx;
             ry  = c*ty;
             if ( b && (n <= kapn))
               {
                  tx = qlambda + sx;
                  sx = rx*tx - ry*sy;
                  sy = ry*tx + rx*sy;
                  qlambda = qlambda/h2;
               }
          }

        if (h == 0.0)
          {
             *u = factor*rx;
             *v = factor*ry;
          }
        else
          {
             *u = factor*sx;
             *v = factor*sy;
          }

        if (yabs == 0.0)
          *u = exp(-xabs*xabs);

     }

/*
c  evaluation of w(z) in the other quadrants
c
*/
   if (yi < 0.0)
     {
        if (a)
          {
             u2    = 2*u2;
             v2    = 2*v2;
          }
        else
          {
             xquad =  -xquad;
/*
c         the following if-statement protects 2*exp(-z**2)
c         against overflow
*/
             if ((yquad > rmaxgoni) || (xquad > rmaxexp))
               return -1;

             w1 =  2*exp(xquad);
             u2  =  w1*cos(yquad);
             v2  = -w1*sin(yquad);
          }

        *u = u2 - *u;
        *v = v2 - *v;
        if (xi > 0.0)
          *v = -(*v);
     }
   else
     {
        if (xi < 0.0)
          *v = -(*v);
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*
 * f = Voigt function = H(x,y)   <dimensionless>
 *
 *     where  x = (nu - nu0) / Doppler_width
 *            y = (Gamma/(4*pi)) / Doppler_width
 *
 *            nu0 = line center frequency (Hz)
 *            Gamma = resonance profile FWHM (Hz)
 *            Doppler_width = s0 * nu0/c  (Hz)
 *            s0 = sqrt (2*k*T/m) = mean thermal velocity
 *            m = ion mass
 */

static int voigt (double x, double y, double *u) /*{{{*/
{
   double v;
   /* z = x + iy  =>   w(z) = u + iv */
   if (-1 == wofz (x, y, u, &v))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "evaluating voigt function for x=%g, y=%g",
                    x, y);
        return -1;
     }

   return 0;
}

/*}}}*/

static int gauss_quad4 (double lo, double hi, double y, double *q) /*{{{*/
{
   static double xi[] = {0.3399810436, 0.8611363116};
   static double wt[] = {0.6521451549, 0.3478548451};
   double a, b, v1, v2, v3, v4;

   a = 0.5 * (hi + lo);
   b = 0.5 * (hi - lo);

   if (   -1 == voigt (a - b*xi[1], y, &v1)
       || -1 == voigt (a - b*xi[0], y, &v2)
       || -1 == voigt (a + b*xi[0], y, &v3)
       || -1 == voigt (a + b*xi[1], y, &v4))
     return -1;

   *q = b * (wt[1]*v1 + wt[0]*v2 + wt[0]*v3 + wt[1]*v4);

   return 0;
}

/*}}}*/

static int binned_voigt (double *val, Isis_Hist_t *g, double *par, unsigned int npar) /*{{{*/
{
   double norm = par[0];
   double e0 = par[1];      /* center [keV] */
   double fwhm = par[2];    /* resonance FWHM [keV] */
   double vtherm = par[3];  /* thermal speed [km/s] */
   double *lo = g->bin_lo;
   double *hi = g->bin_hi;
   int *notice_list = g->notice_list;
   int num = g->n_notice;
   double y, width;
   int i, n;

   (void) npar;

   if ((e0 < 0) || (vtherm < 0) || (fwhm < 0))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "invalid voigt parameters: e0=%g, vtherm=%g fwhm=%g",
                    e0, vtherm, fwhm);
        return -1;
     }

   width = (e0 * vtherm / C_KMS);
   y = (fwhm / FOURPI) / width;

   if (Isis_Voigt_Is_Normalized)
     norm /= width * SQRT_PI;

   for (i=0; i < num; i++)
     {
        double xlo, xhi, v;

        n = notice_list[i];

        xlo = (KEV_ANGSTROM /hi[n] - e0) / width;
        xhi = (KEV_ANGSTROM /lo[n] - e0) / width;

        if (-1 == gauss_quad4 (xlo, xhi, y, &v))
          return -1;

        val[i] = norm * width * v;
     }

   return 0;
}

/*}}}*/

static int contin_voigt (double *val, Isis_User_Grid_t *g, double *par, unsigned int npar) /*{{{*/
{
   double norm = par[0];
   double e0 = par[1];      /* center [keV] */
   double fwhm = par[2];    /* resonance FWHM [keV] */
   double vtherm = par[3];  /* thermal speed [km/s] */
   double y, width;
   double *lam = g->x;
   int n = g->npts;
   int i;

   (void) npar;

   if ((e0 < 0) || (vtherm < 0) || (fwhm < 0))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "invalid voigt parameters: e0=%g, vtherm=%g fwhm=%g",
                    e0, vtherm, fwhm);
        return -1;
     }

   width = (e0 * vtherm / C_KMS);
   y = (fwhm / FOURPI) / width;

   if (Isis_Voigt_Is_Normalized)
     norm /= width * SQRT_PI;

   for (i=0; i < n; i++)
     {
        double x = (KEV_ANGSTROM / lam[i] - e0) / width;
        double v;
        if (-1 == voigt (x, y, &v))
          return -1;
        val[i] = norm * v;
     }

   return 0;
}

/*}}}*/

ISIS_USER_SOURCE_MODULE(voigt,p,options)
{
   static char *parameter_names[] = {"norm", "energy", "fwhm", "vtherm", NULL};
   static char *parameter_units[] = {"photons/s/cm^2", "keV", "keV", "km/s", NULL};
   static double default_min[] = { 0.0, 0.0, 1.e-3, 1.0};
   static double default_max[] = { 10.0, 10.0, 10.0, 1.e5};
   static double default_value[] = {1.0, 1.0, 0.02, 1.e3};
   static unsigned int default_freeze[] = {0, 0, 0, 0};
   static unsigned int norm_indexes[] = {1};

   (void) options;
   
   /* FIXME Very misleading parameter name!!
    * as it stands, fwhm=Gamma but the FWHM is Gamma/(2pi) */

   p->binned = binned_voigt;
   p->unbinned = contin_voigt;

   p->parameter_names = parameter_names;
   p->parameter_units = parameter_units;
   p->norm_indexes = norm_indexes;
   p->num_norms = 1;

   p->default_min = default_min;
   p->default_max = default_max;
   p->default_value = default_value;
   p->default_freeze = default_freeze;

   return 0;
}

