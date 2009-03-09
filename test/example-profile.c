#include <stdlib.h>
#include <math.h>

#include "isis.h"

/* gcc -shared -fPIC -o example-profile.so example-profile.c */

ISIS_LINE_PROFILE_MODULE(square,g,flux,wl,atomic_weight,mid,
                         params,num_params,options)
{
   double width, lo_edge, hi_edge, frac, xl, xh;
   double *lo, *hi, *val;
   int i, n;

   (void) atomic_weight; (void) num_params; (void) options;

   /* params[0] = \Delta\lambda/\lambda */
   width = params[0] * wl;
   if (width <= 0.0)
     return 0;

   lo_edge = wl - 0.5*width;
   hi_edge = wl + 0.5*width;

   lo = g->bin_lo;
   hi = g->bin_hi;
   val = g->val;
   n = g->nbins;

#define INCREMENT_BIN(i) do { \
        if ((hi_edge < lo[i]) || (hi[i] < lo_edge)) \
          break; \
        xl = (lo_edge < lo[i]  ) ? lo[i] : lo_edge; \
        xh = (  hi[i] < hi_edge) ? hi[i] : hi_edge; \
        frac = (xh - xl) / width; \
        val[i] += flux * frac; } while (0)

   for (i = mid; i >= 0; i--)
     {
        INCREMENT_BIN(i);
     }

   for (i = mid+1; i < n; i++)
     {
        INCREMENT_BIN(i);
     }

   return 0;
}
