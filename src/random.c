/*
  Random Number Generators:

  History:
  June 22, 2001:  added to isismath by John Houck <houck@space.mit.edu>
                  minor edits to set up the interface.
  Nov  27, 2002:  included MT19937ar update with improved initialization
*/

/* $Id: random.c,v 1.3 2003/05/15 12:42:12 houck Exp $ */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include "isismath.h"

/* Mersenne Twister Random Number Generator.
 * According to the authors, the period of this generator
 * is 2^19937 - 1.
 */
#include "mt19937ar.c"

/* implementation copied from jdmath by John Davis <davis@space.mit.edu>*/
double grand (void)
{
   static int one_available = 0;
   double g1, g, s;
   static double g2;

   if (one_available)
     {
        one_available = 0;
        return g2;
     }

   one_available = 1;

   do
     {
        g1 = 2.0 * urand () - 1.0;
        g2 = 2.0 * urand () - 1.0;
        g = g1 * g1 + g2 * g2;
     }
   while ((g >= 1.0) || (g == 0.0));

   s = sqrt (-2.0 * log (g) / g);
   g2 = g2 * s;
   return g1 * s;
}

/* start _ptrs Poisson random generator, approximation for high rates
 * 
 * From W. H\"ormann "The Transformed Rejection Method
 * for Generating Poisson Random Variables"
 * http://statistik.wu-wien.ac.at/papers/92-04-13.wh.ps.gz
 */
#define LOG_SQRT_2PI  0.9189385332046727417803296  /* log(sqrt(2pi)) */

static double log_kfact (unsigned int k)
{
   static double fact[] = {1,1,    2,     6,     24,     120, 
                           720, 5040, 40320, 362880, 3628800};
   if (k <= 10)
     return log(fact[k]);
   
   return LOG_SQRT_2PI + (k + 0.5)*log(k) - k + (1.0/12 - 1.0/360/(k*k))/k;
}

static unsigned int _ptrs (double rate)
{
   double a, b, vr, u, v, us, lnmu, ra;
   int k;
   
   b = 0.931 + 2.53 * sqrt(rate);
   a = -0.059 + 0.02483 * b;
   vr = 0.9277 - 3.6224 / (b - 2);
   ra = 1.1239 + 1.1328 / (b - 3.4);
   lnmu = log(rate);
   
   for (;;)
     {
        double lhs, rhs;
        
        u = urand();
        v = urand();
        
        u -= 0.5;
        us = 0.5 - fabs(u);
        
        k = floor ((2 * a / us + b) * u + rate + 0.43);
        
        if (us >= 0.07 && v <= vr)
          return k;
        
        if (k < 0)
          continue;
        
        if (us < 0.013 && v > us)
          continue;
                
        lhs = log(v * ra / (a/(us*us) + b));
        rhs = -rate + k * lnmu - log_kfact ((unsigned int) k);

        if (lhs <= rhs)
          return k;
     }
}

/* end _ptrs */

unsigned int prand (double rate)
{
   double r, lgr, p, lg_nfact, cum;
   unsigned int n;

   if (rate <= 0.0)
     return 0;

   if (rate > 15.0)
     return _ptrs (rate);

   r = urand ();

   lgr = log(rate);
   p = 1.0;
   cum = 0.0;
   lg_nfact = 0.0;
   n = 0;

   while ((cum < r) && ((n <= rate) || (p > 0.0)))
     {
        double lgp = n * lgr - rate - lg_nfact;

        if (lgp > -30.0)
          p = exp (lgp);
        else
          p = 0.0;

        cum += p;
        n++;
        lg_nfact += log (n * 1.0);
     }

   return n - 1;
}

