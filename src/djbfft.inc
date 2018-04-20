/* -*- mode: C; mode: fold -*- */
/*
    This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2018 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

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

/*
 * These are alternate versions of the convolution routines
 * used by the pileup model.  They use the heavily optimized
 * FFT library 'djbfft' written by D. J. Bernstein, and
 * available from
 *
 *   http://cr.yp.to/djbfft.html
 *
 * Informal benchmarks show that the pileup model using djbfft
 * is about five times faster than with Singleton's mixed-radix
 * fftn.
 *
 * Excerpt from djbfft documentation:
 *
 * "fftr4_256 computes a 256-point real-to-complex discrete Fourier
 * transform. It evaluates the polynomial
 *
 *    a[0] + a[2] x + a[4] x^2 + ... + a[254] x^127
 *    + a[1] x^128 + a[3] x^129 + ... + a[255] x^255
 *
 * at all the 256th roots of 1 modulo conjugation, and puts the
 * values into a, overwriting the input. (Beware that the results
 * are stored in an unusual order; also note the input interleaving).
 * Each a[n] is a 4-byte real number.
 *
 * To compute the inverse transform, reconstructing a polynomial
 * from its values, call fftr4_scale256 and then fftr4_un256.
 * (fftr4_scale256 multiplies a[0] and a[1] by 1/256, and each
 * remaining a[n] by 1/128.)"
 *
 * The double-precision routines are analogous.
 */

/* $Id: djbfft.inc,v 1.6 2004/02/09 11:14:18 houck Exp $ */

#include <fftr8.h>

#define DJB_NUM_POINTS  (1024*4)

static void init_fft (double *a, unsigned int n, double *b) /*{{{*/
{
   unsigned int k;

   /* The following is equivalent to
    *  1) copying 'a' (of size n) into the upper half of
    *     a temporary array P (of size 2*n)
    *  2) padding the lower half of P with zeros
    *  3) interleaving the values of P into 'b'
    *     (interleaving is required for input to fftr8_8192).
    *
    * The interleaving permutation is:
    * for (k = 0; k < n; k++)
    * {
    *    b[2*k] = P[k]
    *    b[2*k+1] = P[n+k]
    * }
    * as described in the djbfft documentation excerpt above.
    * Because the lower half of P is all zeros, only the
    * upper half needs to be permuted.
    */

   memset ((char *)b, 0, 2*n*sizeof(double));

   for (k = 0; k < n; k++)
     {
        b[2*k+1] = a[k];
     }
}

/*}}}*/

/* Here fft_s is an array with num*2 elements */
static int setup_convolution_fft (double *s, unsigned int num, /*{{{*/
                                  double *fft_s)
{
   (void) s;

   if (num != DJB_NUM_POINTS)
     return -1;

   init_fft (s, num, fft_s);

   fftr8_8192 (fft_s);
   fftr8_scale8192 (fft_s);

   return 0;
}

/*}}}*/

static int do_convolution (double *fft_s, double *s, unsigned int num, /*{{{*/
                           double *fft_s_tmp)
{
   unsigned int i;

   if (num != DJB_NUM_POINTS)
     return -1;

   init_fft (s, num, fft_s_tmp);

   fftr8_8192 (fft_s_tmp);
   fftr8_mul8192 (fft_s_tmp, fft_s);
   fftr8_un8192 (fft_s_tmp);

   /* Packing Theorem:  The input was padded by a factor 2,
    * so the result is in the even numbered values
    */
   for (i = 0; i < num; i++)
     s[i] = fft_s_tmp[2*i];

   return 0;
}

/*}}}*/

