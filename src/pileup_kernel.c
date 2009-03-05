/*  
    Copyright (C) 2000 Massachusetts Institute of Technology 
 
    Author:  John E. Davis <davis@space.mit.edu>

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

/* $Id: pileup_kernel.c,v 1.3 2003/06/15 03:30:06 houck Exp $ */

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef ISIS_SRC         /* jch */
#include "isismath.h"
#include <limits.h>
#else
#include <jdmath.h>
#endif

#ifdef __linux__
# ifdef __GLIBC__
#  include <fpu_control.h>
# else
#  include <i386/fpu_control.h>
# endif
# ifndef _FPU_SETCW
#  define _FPU_SETCW __setfpucw
# endif
#endif

#define ISIS_KERNEL_PRIVATE_DATA \
   unsigned int num_terms; \
   unsigned int max_num_terms; \
   double num_frames; \
   int has_psf_fraction; \
   unsigned int num; \
   double *arf_time; 		       /* arf * frametime */ \
   double *energies; \
   double de; \
   double *arf_s;		       /* A(E)s(E) */ \
   double *arf_s1; \
   double *results; \
   double *arf_s_fft; \
   double *arf_s_tmp; \
   double *arf_s_fft_tmp; \
   double *pileup_fractions; \
   double integral_ae; \
   double arf_frac_exposure; \
   int verbose;


#include "isis.h"
#include "util.h"
#include "errors.h"

static unsigned int Num_Parameters;

#define PILEUP_VERSION "0.1.1"

/* Make this a power of 2 so that ffts are fast */
#define NUM_POINTS      (1024*4)
#define ANGSTROM_TO_KEV 12.3984185734
#define MAX_NUM_TERMS	30

#define XMALLOC(n,t) (t *)malloc((n)*sizeof(t))

static double Max_Probability_Cutoff = 0.99999;

#ifndef JDMATH_VERSION
/* THese functions borrowed from jdmath */
static unsigned int JDMbinary_search_d (double x, double *xp, unsigned int n)
{
   unsigned int n0, n1, n2;
   
   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
	n2 = (n0 + n1) / 2;
	if (xp[n2] >= x) 
	  {
	     if (xp[n2] == x) return n2;
	     n1 = n2;
	  }
	else n0 = n2;
     }
   return n0;
}

   
#if 0
static double JDMinterpolate_d (double x, double *xp, double *yp, unsigned int n)
{
   unsigned int n0, n1;
   double x0, x1;
   
   n0 = JDMbinary_search_d (x, xp, n);
   
   x0 = xp[n0];
   n1 = n0 + 1;
   
   if ((x == x0) || (n1 == n)) return yp[n0];
   
   x1 = xp[n1];
   if (x1 == x0) return yp[n0];
   
   return yp[n0] + (yp[n1] - yp[n0]) / (x1 - x0) * (x - x0);
}
#endif

#endif
#ifdef HAVE_DJBFFT
#include "djbfft.inc"
#else
/* Here fft_s is a complex array with num*2*2 complex elements */
static int setup_convolution_fft (double *s, unsigned int num, 
				  double *fft_s)
{
   double *right;
   unsigned int i;
   int num2;
   
   num2 = (int) num * 2;

   /* Set the left part to 0 --- this serves as a pad */
   memset ((char *) fft_s, 0, num2 * sizeof (double));
   right = fft_s + num2;

   for (i = 0; i < num; i++)
     {
	*right++ = s[i];	       /* real */
	*right++ = 0.0;	       /* imag */
     }

   return JDMfftn (1, &num2, fft_s, fft_s + 1, 2, -2);
}

static int do_convolution (double *fft_s, double *s, unsigned int num, 
			   double *fft_s_tmp)
{
   unsigned int num4;
   unsigned int i;
   int num2;

   (void) setup_convolution_fft (s, num, fft_s_tmp);

   /* Multiply transforms */
   num4 = 4 * num;
   i = 0;
   while (i < num4)
     {
	double re0, im0, re1, im1;
	unsigned int i1;
	
	i1 = i + 1;

	re0 = fft_s_tmp[i];
	re1 = fft_s[i];
	im0 = fft_s_tmp[i1];
	im1 = fft_s[i1];
	
	fft_s_tmp[i] = re0 * re1 - im0 * im1;
	fft_s_tmp[i1] = re0 * im1 + re1 * im0;

	i = i1 + 1;
     }

   num2 = (int) 2 * num;
   /* Perform inverse and return left half of real part */
   (void) JDMfftn (1, &num2, fft_s_tmp, fft_s_tmp + 1, -2, -2);
   for (i = 0; i < num; i++)
     s[i] = fft_s_tmp[2*i];
	
   return 0;
}
#endif
static int add_in_xspec (Isis_Kernel_t *k, double *arf_s, double psf_frac,
			 double *results)
{
   unsigned int i;
   unsigned int num = k->num;

   psf_frac *= k->num_frames;
   for (i = 0; i < num; i++)
     results[i] += arf_s[i] * psf_frac;
   
   return 0;
}

static int perform_pileup (Isis_Kernel_t *k,
			   double *arf_s,
			   double reg_size,
			   double g0,
			   double psf_frac,
			   double *gfactors,
			   double *pileup_dist,
			   double *results)
{
   unsigned int i;
   unsigned int num;
   unsigned int i_factorial;
   double fft_norm;
   double *fft_s;
   double *arf_s_fft_tmp;
   double *arf_s_tmp;
   double exp_factor;
   double integ_arf_s, integ_arf_s_n;
   double total_prob;

   num = k->num;
   fft_s = k->arf_s_fft;
   arf_s_tmp = k->arf_s_tmp;

   integ_arf_s = 0.0;
   psf_frac = psf_frac / reg_size;

   if (k->arf_frac_exposure > 0)
     psf_frac /= k->arf_frac_exposure;

   for (i = 0; i < num; i++)
     {
	arf_s_tmp[i] = results[i] = arf_s[i] * psf_frac;
	integ_arf_s += arf_s_tmp[i];
     }

   if (pileup_dist != NULL)
     {
	pileup_dist[0] = 0.0;
	pileup_dist[1] = integ_arf_s;
     }

   k->integral_ae = integ_arf_s/g0;
   k->num_terms = k->max_num_terms;

   if (integ_arf_s == 0.0)
     return 0;

   /* Normalize by integ_arf_s to avoid possible floating point overflow. 
    * This will be corrected below.
    */
   for (i = 0; i < num; i++)
     arf_s_tmp[i] /= integ_arf_s;

   exp_factor = exp (-k->integral_ae);

   (void) setup_convolution_fft (arf_s_tmp, num, fft_s);

   arf_s_fft_tmp = k->arf_s_fft_tmp;

   fft_norm = sqrt (2.0 * num);

   i_factorial = 1;
   integ_arf_s_n = integ_arf_s;
   total_prob = 1 + integ_arf_s;

   for (i = 2; i <= k->max_num_terms; i++)
     {
	unsigned int j;
	double norm_i;
	
	i_factorial *= i;
	integ_arf_s_n *= integ_arf_s;
	norm_i = integ_arf_s_n / i_factorial;
	total_prob += norm_i;

	(void) do_convolution (fft_s, arf_s_tmp, num, arf_s_fft_tmp);
	
	norm_i *= gfactors[i-2];

	for (j = 0; j < num; j++)
	  {
#ifndef HAVE_DJBFFT	     
	     arf_s_tmp[j] *= fft_norm;
#endif	     
	     results[j] += norm_i * arf_s_tmp[j];
	  }

	if (pileup_dist != NULL)
	  pileup_dist [i] = norm_i;

	if (total_prob * exp (-integ_arf_s) > Max_Probability_Cutoff)
	  {
	     k->num_terms = i;
	     break;
	  }
     }

   exp_factor *= k->num_frames * reg_size;

   /* Apply correction to account for the number of effective frames */
   exp_factor *= k->arf_frac_exposure;

   for (i = 0; i < num; i++)
     results [i] *= exp_factor;
   
   return 0;
}

static int evaluate_isis_model_on_grid (double *lamlo, double *lamhi,
					unsigned int num,
					double *val)
{
   Isis_Hist_t h;
   unsigned int i;
   int status;
	     
   memset ((char *)&h, 0, sizeof (h));
   
   h.bin_lo = lamlo;
   h.bin_hi = lamhi;
   h.nbins = num;
   h.n_notice = num;
   h.val = val;

   if (NULL == (h.notice_list = XMALLOC (num, int)))
     return -1;
   
   for (i = 0; i < num; i++)
     h.notice_list[i] = i;
   
   status = isis_eval_model_on_alt_grid (&h);
   
   free ((char *) h.notice_list);
   return status;
}


/* By definition, S(y)|dy| = s(E)|dE|, where E = hc/y.
 * So, s(E) = S(y) |dy/dE|
 *          = S(y) hc/E^2
 * We need the product, A(E)s(E)
 */
static int convert_spectrum (int (*fun)(Isis_Hist_t *), Isis_Kernel_t *k, Isis_Hist_t *g)
{
   double *s_dlam;
   double *arf_s, *energies, *arf;
   unsigned int i, num;
   double *lamlo, *lamhi;
   double de;

   /* FIXME???  Does isis depend up g->val having meaningful values? */
   (void) fun; (void) g;

   num = k->num;
   if (NULL == (lamlo = XMALLOC (3*num, double)))
     return -1;
   lamhi = lamlo + num;
   s_dlam = lamhi + num;

   de = k->de;
   energies = k->energies;

   if (de >= energies[0])
     lamhi[num-1] = ANGSTROM_TO_KEV/(energies[0]/2.0);
   else
     lamhi[num-1] = ANGSTROM_TO_KEV/(energies[0]);

   for (i = 1; i < num; i++)
     lamhi[num-1-i] = ANGSTROM_TO_KEV/(energies[i]);
   for (i = 1; i < num; i++)
     lamlo[i] = lamhi[i-1];
   lamlo[0] = ANGSTROM_TO_KEV/(energies[num-1] + de);

   if (-1 == evaluate_isis_model_on_grid (lamlo, lamhi, num, s_dlam))
     {
	free ((char *) lamlo);
	return -1;
     }

   arf_s = k->arf_s;
   arf = k->arf_time;

   for (i = 0; i < num; i++)
     {
	unsigned int j = num-1-i;
	arf_s[i] = arf[i] * s_dlam[j];

	if (arf_s[i] < 0)
	  arf_s[i] = 0;
     }
   free ((char *)lamlo);

   return 0;
}

static int convert_results (Isis_Kernel_t *k, Isis_Hist_t *g)
{
   unsigned int i;
   double *results;
   double *s_dlam;
   double *energies;
   unsigned int num;
   double *lamlo, *lamhi;
   double de;


   num = k->num;
   results = k->results;
   energies = k->energies;

   num = k->num;
   if (NULL == (lamlo = XMALLOC (3*num, double)))
     return -1;
   lamhi = lamlo + num;
   s_dlam = lamhi + num;

   de = k->de;
   energies = k->energies;

   if (de >= energies[0])
     lamhi[num-1] = ANGSTROM_TO_KEV/(energies[0]/2.0);
   else
     lamhi[num-1] = ANGSTROM_TO_KEV/(energies[0]);
   
   s_dlam[num-1] = results[0];
   for (i = 1; i < num; i++)
     {
	unsigned int j = num - 1 - i;
	lamhi[j] = ANGSTROM_TO_KEV/(energies[i]);
	s_dlam[j] = results[i];
     }

   for (i = 1; i < num; i++)
     lamlo[i] = lamhi[i-1];
   lamlo[0] = ANGSTROM_TO_KEV/(energies[num-1] + de);

   if (-1 == rebin_histogram (s_dlam, lamlo, lamhi, num,
			      g->val, g->bin_lo, g->bin_hi, g->nbins))
     {
	free ((char *)lamlo);
	return -1;
     }
   free ((char *)lamlo);

   return 0;
}

   
static void init_coeffs (double alpha,
			 double *coeffs, unsigned int npiles)
{
   unsigned int i;

   alpha = fabs(alpha);
   for (i = 2; i <= npiles; i++)
     {
	coeffs[i - 2] = pow (alpha, i-1);
     }
}

static int compute_kernel (Isis_Kernel_t *k, double *result, Isis_Hist_t *g, double *pars, unsigned int num_pars,
			   int (*fun)(Isis_Hist_t *))
{
   double coeffs[MAX_NUM_TERMS+1];
   double alpha, g0, num_regions;
   double psf_frac;

   if (k == NULL || g == NULL || fun == NULL)
     return -1;

   if (num_pars != Num_Parameters)
     {
	isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "pileup: Expecting %d pileup parameters", Num_Parameters);
	return -1;
     }

   if (g->n_notice != g->nbins)
     {
	isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "pileup: Pileup requires that all bins be noticed.");
	return -1;
     }
   
   if (-1 == convert_spectrum (fun, k, g))
     return -1;

   num_regions = pars[0];
   g0 = pars[1];
   alpha = pars[2];
   psf_frac = pars[3];

   init_coeffs (alpha, coeffs, k->max_num_terms);
   if (k->has_psf_fraction)
     psf_frac = pars[3];
   else
     psf_frac = 1.0;

   if (-1 == perform_pileup (k, k->arf_s, num_regions, g0, psf_frac,
			     coeffs, k->pileup_fractions,
			     k->results))
     return -1;
   
   psf_frac = 1.0 - psf_frac;
   if (psf_frac > 0.0)
     {
	if (-1 == add_in_xspec (k, k->arf_s, psf_frac, k->results))
	  return -1;
     }

   if (-1 == convert_results (k, g))
     return -1;

   return k->apply_rmf (k->rsp.rmf, result, k->num_orig_data, 
			g->val, g->notice_list, g->n_notice);
}

static int save_results (Isis_Kernel_t *k)
{
   unsigned int i;

   FILE *fp = fopen ("results.dat", "w");
   if (fp == NULL)
     return -1;

   for (i = 0; i < k->num; i++)
     {
	fprintf (fp, "%e\t%e\n", k->energies[i], k->results[i]);
     }
   
   fclose (fp);
   return 0;
}

static int print_kernel (Isis_Kernel_t *k)
{
   unsigned int i;
   double sum = 0.0;
   double pn, sum_piled;

   if (k == NULL)
     return -1;
   
   sum = sum_piled = 0.0;
   if (k->num_terms >= 1)
     sum = k->pileup_fractions[1];

   for (i = 2; i <= k->num_terms; i++)
     sum_piled += k->pileup_fractions[i];

   sum += sum_piled;

   pn = exp (-k->integral_ae);

   for (i = 1; i <= k->num_terms; i++)
     {
	pn *= k->integral_ae/i;
	fprintf (stdout, "%u: %g\t%g\n", i,
		 pn, k->pileup_fractions[i]/sum);
     }
   fprintf (stdout, "*** pileup fraction: %g\n", sum_piled/sum);
   fflush (stdout);

   /* FIXME? - requires mechanism for user-defined filename */
   if (0) save_results (k);

   return 0;
}

static int handle_nterms_option (char *subsystem, char *optname, char *value, void *clientdata)
{
   Isis_Kernel_t *k;
   k = (Isis_Kernel_t *) clientdata;
   if (1 != sscanf (value, "%u", &k->max_num_terms))
     {
	fprintf (stderr, "%s;%s option requires an unsigned int\n", subsystem, optname);
	return -1;
     }
   if (k->max_num_terms >= MAX_NUM_TERMS)
     k->max_num_terms = MAX_NUM_TERMS-1;

   return 0;
}

static int handle_fracexpo_option (char *subsystem, char *optname, char *value, void *clientdata)
{
   Isis_Kernel_t *k;
   k = (Isis_Kernel_t *) clientdata;
   if ((1 != sscanf (value, "%lf", &k->arf_frac_exposure))
       || (k->arf_frac_exposure <= 0.0))
     {
	fprintf (stderr, "%s;%s option requires a positive number\n", subsystem, optname);
	return -1;
     }
   return 0;
}

static int handle_verbose_option (char *subsystem, char *optname, char *value, void *clientdata)
{
   Isis_Kernel_t *k;
   k = (Isis_Kernel_t *) clientdata;
   if (1 != sscanf (value, "%d", &k->verbose))
     {
	fprintf (stderr, "%s;%s option requires an int\n", subsystem, optname);
	return -1;
     }
   return 0;
}


static Isis_Option_Table_Type Option_Table [] = 
{
     {"nterms", handle_nterms_option, ISIS_OPT_REQUIRES_VALUE, "30", "Max number of piled photons"},
     {"fracexpo", handle_fracexpo_option, ISIS_OPT_REQUIRES_VALUE, "1.0", "Fraction exposure for ARF"},
     {"verbose", handle_verbose_option, ISIS_OPT_REQUIRES_VALUE, "0", "verbose level"},
     ISIS_OPTION_TABLE_TYPE_NULL
};

static int process_options (Isis_Kernel_t *k, char *options)
{
   Isis_Option_Type *o;
   int status;

   k->max_num_terms = MAX_NUM_TERMS;
   k->arf_frac_exposure = 1.0;
   
   o = isis_parse_option_string (options);
   if (o == NULL)
     return -1;

   status = isis_process_options (o, Option_Table, (void *)k, 1);
   isis_free_options (o);
   return status;
}

	  
static void delete_kernel (Isis_Kernel_t *k)
{
   if (NULL == k)
     return;

   if (k->arf_time != NULL) free (k->arf_time);
   if (k->energies != NULL) free (k->energies);
   if (k->arf_s != NULL) free (k->arf_s);
   if (k->arf_s1 != NULL) free (k->arf_s1);
   if (k->results != NULL) free (k->results);
   if (k->arf_s_fft != NULL) free (k->arf_s_fft);
   if (k->arf_s_tmp != NULL) free (k->arf_s_tmp);
   if (k->arf_s_fft_tmp != NULL) free (k->arf_s_fft_tmp);
   if (k->pileup_fractions != NULL) free (k->pileup_fractions);

   free (k);
}

static Isis_Kernel_t *allocate_kernel (Isis_Obs_t *o, char *options)
{
   Isis_Kernel_t *k = NULL;
   double frame_time;
   unsigned int i;
   double *arf;
   double max_energy, delta_energy, min_energy;
   double *bin_lo, *bin_hi;
   double *eff_area;
   unsigned int num_bins;
   double *energies;
   unsigned int max_num_terms;
   double max_lambda, min_lambda;

   if (o->rsp.next != NULL)
     {
	isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "pileup:  Use with multiple responses is not supported");
	return NULL;
     }

   if (NULL == (eff_area = o->rsp.arf->arf))
     {
	isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "pileup:  An ARF is required.");
	return NULL;
     }

   if (o->rsp.arf->fracexpo_is_vector)
     {
	isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "pileup:  ARF inconsistent with pileup kernel (fracexpo should be a scalar)");
	return NULL;
     }

   k = isis_init_kernel (NULL, sizeof(*k), o);
   if (NULL == k) return NULL;

   /* options */
   if (-1 == process_options (k, options))
     {
	delete_kernel (k);
	return NULL;
     }
   max_num_terms = k->max_num_terms;
   
   /* public fields */
   k->delete_kernel = delete_kernel;
   k->compute_kernel = compute_kernel;
   k->print_kernel = print_kernel;
   
   k->max_num_terms = max_num_terms;
   
   k->arf_frac_exposure = o->rsp.arf->fracexpo.s;
   
   frame_time = o->frame_time;
   if (frame_time == 0.0)
     {
	isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "pileup: frame_time not set (use set_frame_time to set it)");
	frame_time = 3.2; /* 0.00285; */
	isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "pileup: assuming frame_time = %g sec", frame_time);
     }

   k->frame_time = frame_time;

   /* Private fields */
   k->num_frames = k->exposure_time / frame_time;
   k->num = NUM_POINTS;
   k->has_psf_fraction = 1;

   if (k->verbose)
     {
	isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "pileup: version %s\n", PILEUP_VERSION);
	isis_vmesg (WARN, I_INFO, __FILE__, __LINE__, "pileup: frame_time=%g, fracexpo=%g, num_frames=%g",
		    k->frame_time, k->arf_frac_exposure, k->num_frames);
     }

   if ((NULL == (k->arf_time = XMALLOC (NUM_POINTS, double)))
       || (NULL == (k->energies = XMALLOC (NUM_POINTS, double)))
       || (NULL == (k->arf_s = XMALLOC (NUM_POINTS, double)))
       || (NULL == (k->arf_s1 = XMALLOC (NUM_POINTS, double)))
       || (NULL == (k->results = XMALLOC (NUM_POINTS, double)))
       || (NULL == (k->arf_s_fft = XMALLOC (4*NUM_POINTS, double)))
       || (NULL == (k->arf_s_tmp = XMALLOC (NUM_POINTS, double)))
       || (NULL == (k->arf_s_fft_tmp = XMALLOC (4*NUM_POINTS, double)))
       || (NULL == (k->pileup_fractions = XMALLOC (max_num_terms + 1, double))))
     {
	delete_kernel (k);
	return NULL;
     }
   
   bin_lo = o->rsp.arf->bin_lo;
   bin_hi = o->rsp.arf->bin_hi;
   num_bins = o->rsp.arf->nbins;
   
   if (bin_lo[0] > 0.0)
     max_energy = ANGSTROM_TO_KEV / bin_lo[0];
   else
     max_energy = 12.0;

   max_energy = 15.0;
   min_energy = 0.0;
   /* max_energy = 1024 * 0.0146; */
   k->de = delta_energy = (max_energy - min_energy) / NUM_POINTS;
   
   arf = k->arf_time;
   energies = k->energies;

   arf[0] = 0.0;
   energies[0] = 0.0001;	       /* MUST be non-zero */
   
   min_lambda = bin_lo[0];
   max_lambda = bin_hi[num_bins-1];
   for (i = 1; i < NUM_POINTS; i++)
     {
	double lambda;
	
	energies[i] = min_energy + delta_energy * i;
	lambda = ANGSTROM_TO_KEV / energies[i];
	if ((lambda >= max_lambda) || (lambda < min_lambda))
	  {
	     arf[i] = 0;
	     continue;
	  }
	arf[i] = frame_time * eff_area[JDMbinary_search_d (lambda, bin_lo, num_bins)];
	if (arf[i] < 0)
	  arf[i] = 0;
     }
   
   return k;
}

static int parse_init_options (char *options)
{
   (void) options;
   
   return 0;
}

ISIS_USER_KERNEL_MODULE(pileup,def,options)
{
#define NUM_PILEUP_PARMS 4
   static char *kernel_parm_names[NUM_PILEUP_PARMS] = 
     {
	"nregions", "g0", "alpha", "psffrac"
     };
   static double default_min [] =
     {
	1, 0, 0, 0
     };
   static double default_max [] =
     {
	10, 1, 1, 1
     };
   static double default_value [] =
     {
	1, 1, 0.5, 0.95
     };
   static unsigned int default_freeze [] =
     {
	1, 1, 0, 1,	
     };
   
#if 0
#ifdef _FPU_SETCW
   unsigned int cw = 0x137f;
   cw = 0x1372;
   _FPU_SETCW (cw);
#endif
#endif
   
   if (-1 == parse_init_options (options))
     return -1;

   Num_Parameters = NUM_PILEUP_PARMS;

   def->kernel_name = "pileup";
   def->allocate_kernel = allocate_kernel;
   def->num_kernel_parms = Num_Parameters;
   def->kernel_parm_names = kernel_parm_names;
   
   def->default_min = default_min;
   def->default_max = default_max;   
   def->default_value = default_value;
   def->default_freeze = default_freeze;

   def->allows_ignoring_model_intervals = NULL;
   
   return 0;
}

#ifdef TEST_PILEUP
int main ()
{
   Isis_Hist_Obs_t o;
   double arf[10];
   double lo[10];
   double hi[10];
   double result[10];
   double val[10];
   Isis_Hist_t g;
   int i;
   double lam;
   Isis_Kernel_t *k;

   lam = .1;
   for (i = 0; i < 10; i++)
     {
	lo[i] = lam;
	lam += 0.01;
	hi[i] = lam;
	arf[i] = 1.0;
	val[i] = i;
     }
   o.frame_time = 0.0;
   o.exposure_time = 800.0;
   o.eff_area = arf;
   o.num = 10;
   o.arf_bin_lo = lo;
   o.arf_bin_hi = hi;
   
   g.n_notice = 10;
   g.nbins = 10;
   g.bin_lo = lo;
   g.bin_hi = lo;
   g.val = val;

   k = pileup_allocate_kernel (&o, 7);
   
   (void) k->compute_kernel (k, &g, result);
   (void) k->delete_kernel (k);
   
   return 0;
}

#endif
