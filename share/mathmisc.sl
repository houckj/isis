% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2008 Massachusetts Institute of Technology
%
%    This software was developed by the MIT Center for Space Research under
%    contract SV1-61010 from the Smithsonian Institution.
%
%    Author:  John C. Houck  <houck@space.mit.edu>
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

% $Id: mathmisc.sl,v 1.14 2004/09/09 11:31:55 houck Exp $

#ifnexists cumsum
define cumsum () %{{{
{
   variable msg = "s[] = cumsum (x[] [, axis])";
   variable a, axis;

   axis = 0;

   if (_isis->get_varargs (&a, &axis, _NARGS, 1, msg))
     return;

   if (typeof(a) != Array_Type)
     a = [a];

   % pass a copy so the original is not overwritten
   return _isis->_cumsum (@a, axis);
}
%}}}
#endif

#ifnexists sum
define sum (a)
{
   variable tot = 0.0;

   foreach (a)
     {
	tot += ();
     }

   return tot;
}
#endif

% The slang reverse function doesn't support scalars so
% I'm replacing it.
%#ifnexists reverse
define reverse (a)
{
   variable i = length (a);
   if (i <= 1)
     return a;

   i--;
   __tmp(a)[[i:0:-1]];
}
%#endif

#ifnexists shift
define shift (x, n)
{
   variable len = length(x);
   variable i = [0:len-1];

   % allow n to be negative and large
   n = len + n mod len;
   return x[(i + n)mod len];
}
#endif

#ifnexists ones
define ones ()
{
   !if (_NARGS) return 1;
   variable a;

   a = __pop_args (_NARGS);
   return @Array_Type (Integer_Type, [__push_args (a)]) + 1;
}
#endif

#ifnexists howmany
define howmany (a)
{
   return length(where(a));
}
#endif

#ifnexists any
define any (a)
{
   return (0 != howmany(a));
}
#endif

define moment () %{{{
{
   variable msg = "stat_struct = moment (values)";
   variable values;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   values = ();

   return _isis->_moment (values);
}

%}}}

#ifnexists min
define min (a)
{
   a = _reshape(a,[length(a)]);

   return typecast (moment(a).min, _typeof(a));
}
#endif

#ifnexists max
define max (a)
{
   a = _reshape(a,[length(a)]);

   return typecast (moment(a).max, _typeof(a));
}
#endif

#ifnexists mean
define mean (a)
{
   a = _reshape(a,[length(a)]);

   return moment(a).ave;
}
#endif

define median () %{{{
{
   variable msg = "med = median (values[])";
   variable values;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   values = ();

   % operate on a copy to preserve the original
   return typecast (_isis->_median (@values), _typeof(values));
}

%}}}

#ifnexists hypot
define hypot () %{{{
{
   variable msg = "r[] = hypot (x[], y[])";
   variable x, y;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (x, y) = ();

   variable r = _isis->_hypot (x, y);

   if (r == NULL)
     return r;

   if (length(r) == 1)
     return r[0];

   return r;
}
#endif

%}}}

#ifnexists finite
define finite () %{{{
{
   variable msg = "k[] = finite (x[])";
   variable x;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   x = ();

   variable r = _isis->_finite(x);

   if (r == NULL)
     return r;

   if (length(r) == 1)
     return r[0];

   return r;
}
#endif

%}}}

define histogram () %{{{
{
   variable msg = "nx = histogram (x, lo [, hi [, &reverse]])";
   variable x, lo, hi, rev, flag;

   flag = 0;
   rev = NULL;
   hi = NULL;
   if (_isis->get_varargs (&x, &lo, &hi, &rev, _NARGS, 2, msg))
     return;

   if (rev != NULL)
     flag = 1;

   % generate the hi-bin edges with the last as an overflow bin.
   if (hi == NULL)
     {
	lo = typecast (lo, Double_Type);
	hi = [lo[[1:length(lo)-1]], _isis->DBL_MAX];
     }

   variable nx, r;
   (nx, r) = _isis->_make_1d_histogram (x, lo, hi, flag);

   if (r != NULL)
     @rev = r;

   return nx;
}

%}}}

define histogram2d () %{{{
{
   variable msg = "img[] = histogram2d (x, y, grid_x, grid_y [, &reverse])";
   variable x, y, grid_x, grid_y, rev, flag;

   flag = 0;
   rev = NULL;
   if (_isis->get_varargs (&x, &y, &grid_x, &grid_y, &rev, _NARGS, 4, msg))
     return;

   if (rev != NULL)
     flag = 1;

   variable num, r;
   (num, r) = _isis->_make_2d_histogram (x, y,
					 grid_x[array_sort(grid_x)],
					 grid_y[array_sort(grid_y)],
					 flag);
   if (rev != NULL)
     {
	variable dims;
	(dims,,) = array_info(num);
	reshape (r, dims);
	@rev = r;
     }

   return num;
}

%}}}

define fft1d () %{{{
{
   variable msg = "(re, im) = fft1d (real, imag, sign)";
   variable re, im, re_cpy, im_cpy, sgn;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (re, im, sgn) = ();

   if (typeof(re) != Array_Type
       or typeof(im) != Array_Type)
     {
	usage(msg);
	return;
     }
   
   variable num_re_dims, num_im_dims;
   (,num_re_dims,) = array_info (re);
   (,num_im_dims,) = array_info (im);
   if (num_re_dims != 1 or num_im_dims != 1)
     {
        message ("*** Warning:  fft supports 1-D arrays only");
        return;
     }   

   re_cpy = @re;
   im_cpy = @im;

   if (sgn > 0) sgn = 1;
   else sgn = -1;

   % -2.0 = normalization choice

   return _isis->_fft1d (re_cpy, im_cpy, sgn, -2.0);
}

%}}}

define fft () %{{{
{
   variable msg = "X[] = fft (x[], sign)";
   variable re, im, x, sgn;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;
   (x, sgn) = ();
   (re, im) = fft1d (Real(x), Imag(x), sgn);
   return re + im * 1i;
}

%}}}

define ks_diff (a,b) %{{{
{
   variable diff = _isis->_ks_difference (a,b);
   if (diff < 0.0)
     message ("ks_diff failed");
   return diff;
}

%}}}

define ks_prob (lam) %{{{
{
   return _isis->_ks_probability (lam);
}

%}}}

define seed_random () %{{{
{
   variable msg = "seed_random (int)";
   variable seed;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   seed = ();

   _isis->_seed_random (seed);
}

%}}}

private define _genrand (nargs, randfun)
{
   variable num, a = NULL;

   if (nargs == 0)
     num = 1;
   else if (nargs == 1)
     {
        a = ();
     }
   else
     {
        a = __pop_args (nargs);
        a = [__push_args(a)];
     }

   if (a != NULL)
     {
        num = 1;
        foreach (a)
          {
             variable x = ();
             num *= x;
          }
     }

   variable r = (@randfun)(num);

   if (a == NULL)
     return r;

   return _reshape(r,a);
}

define urand () %{{{
{
   return _genrand (_NARGS, &_isis->_urand_array);
}

%}}}

define grand () %{{{
{
   return _genrand (_NARGS, &_isis->_grand_array);
}

%}}}

define prand () %{{{
{
   variable msg = "x = prand (rate [,num])";
   variable r, n;

   n = 1;

   if (_isis->get_varargs (&r, &n, _NARGS, 1, msg))
     return;

   if (typeof(r) == Array_Type)
     return _isis->_prand_vec (@r);  % don't overwrite r

   return _isis->_prand_array (r, n);
}

%}}}

#ifnexists unique
define unique () %{{{
{
   variable i, j, len;
   variable a;

   if (_NARGS != 1)
     {
	_pop_n (_NARGS);
	usage ("i = unique (a); %% i = indices of unique elements of a");
     }

   a = ();
   len = length(a);
   if (len <= 1)
     return [0:len-1];

   a = _reshape (__tmp(a),[len]);

   i = array_sort(a);
   a = a[i];

   variable eps = 0.0;

   j = (abs(a[[1:]] - a[[0:len-2]]) > eps);
   % We need to get the last element
   j = [j, '\001'];
   %reshape (j, [len]);

   j = where(j);
   % Now, i contains the sorted indices, and j contains the indices into the
   % sorted array.  So, the unique elements are given by a[i][j] where a is
   % the original input array.  It seems amusing that the indices given by
   % [i][j] are also given by i[j].
   return i[j];
}

%}}}
#endif

% binary search a histogram grid
define bsearch_hist (x, lo, hi) %{{{
{
   variable n = length(lo);
   if (n == 0)
     return NULL;

   if (isnan(x))
     return NULL;

   if (x < lo[0])
     return -INT_MAX;
   else if (hi[n-1] < x)
     return INT_MAX;

   variable n0=0, n1=n, n2;

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

%}}}

% binary search for x in sorted array a,
% returning index i such that a[i] <= x < a[i+1]
define bsearch (x, a) %{{{
{
   variable n = length(a);
   if (n == 0)
     return NULL;

   if (isnan(x))
     return NULL;

   if (x < a[0])
     return -INT_MAX;
   else if (a[n-1] < x)
     return INT_MAX;

   variable n0=0, n1=n, n2;

   while (n1 > n0 + 1)
     {
        n2 = (n0 + n1) / 2;
        if (x < a[n2])
          {
             if (a[n2] == x)
               return n2;
             n1 = n2;
          }
        else n0 = n2;
     }

   return n0;
}

%}}}
