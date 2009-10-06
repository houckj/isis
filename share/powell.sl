% Powell's direction set method
% Copyright (C) 2009 John E. Davis
%
% This file is part of the S-Lang Optimization Package and may be
% distributed under the terms of the GNU General Public License.  See
% the file COPYING for more information.
%
% Powell's method is described in the paper:
%   "An efficient method for finding the minimum of a function of several
%    variables without calculating derivatives", The Computer Journal 1964
% Some of the implementation here makes use of:
%   W.I. Zangwill, "Minimizing a function without derivatives", Computer Journal 1967
%
require ("rand");

private variable _macheps = 1.0;
while (1+0.5*_macheps != 1.0) _macheps = 0.5*_macheps;

private variable Usage_Msg = "\
Usage:\n\
minF = obj.optimize(parms, min_parms, max_parms, func, func_data);\n\
Qualifiers:\n\
 maxnfe=500, reldiff=1e-4, absdiff=1e-9 verbose=0\n\
\n\
The value at the minimum is returned.  The parms array values are\n\
set to the parameters at the minimum.\n\
func should be defined as: define func(parms, func_data)\n\
";

private define powell ();
private variable Powell_Defaults = struct
{
   maxnfe = 500,
   nfe = 0,
   converged = 0,
   reldiff = 1e-4,
   absdiff = 1e-9,
   verbose = 0,
   new_min_hook = NULL,
   new_min_hook_args = {},
   optimize = &powell,
};

private define powell_usage ()
{
   usage (Usage_Msg);
}

private define override_defaults (sinfo, q)
{
   if (q == NULL)
     return;
   variable sinfo_fields = get_struct_field_names (sinfo);

   foreach (get_struct_field_names (q))
     {
	variable name = ();
	if (any (name == sinfo_fields))
	  set_struct_field (sinfo, name, get_struct_field (q, name));
     }
}

define powell_new ()
{
   variable s = @Powell_Defaults;
   override_defaults (s, __qualifiers());
   return s;
}

private define eval_func (p, func, func_data)
{
   variable y = (@func)(p, func_data);
   if (isnan(y))
     {
	() = fprintf (stderr, "NaN seen\n");
	y = _Inf;
     }
   return y;
}

private define call_func (lambda, cd)
{
   % func=cd[0], func_data=cd[1],
   variable p = cd.p + lambda * cd.xi;
   ifnot (all (cd.min_p <= p <= cd.max_p))
     return _Inf;
   return eval_func (p, cd.func, cd.func_data);
}

private define quadratic_min (a, b, c, fa, fb, fc, func, func_data)
{
   variable 
     apb=a+b, amb=a-b,
     bpc=b+c, bmc=b-c,
     cpa=c+a, cma=c-a;
   variable d=0.5*(bmc*bpc*fa + cma*cpa*fb + amb*apb*fc)
     /(bmc*fa + cma*fb + amb*fc);
  
   return d, (@func)(d,func_data);
}

% Here, p_norm = |p|, which is the norm of the current search point.
% xi_norm is |xi| and represents the value of the last search in the
% direction xi.
private define linemin (sinfo, f0, h, d2f, min_x, max_x, abserr, func, func_data, max_linemin_nfe)
{
   variable min_d2f = _macheps;
   
   variable x0 = 0.0, x1, x2, f1, f2;
   variable xmin = x0, fmin = f0;

   x1 = h;
   if (abs(x0-min_x) < abserr)
     {
	x1 = abs(h);
	if (x1 > max_x)
	  x1 = 0.5 * (min_x + max_x);
     }
   if (abs(x0-max_x) < abserr)
     {
	x1 = -abs(h);
	if (x1 < min_x)
	  x1 = 0.5 * (min_x + max_x);
     }

   if (x1 < min_x)
     x1 = 0.5 * (min_x + x0);
   if (x1 > max_x)
     x1 = 0.5 * (max_x + x0);

   f1 = (@func)(x1, func_data); sinfo.nfe++; max_linemin_nfe--;
   if (f1 <= fmin)
     (xmin, fmin) = (x1, f1);

   if (x1 < x0)
     (x0, f0, x1, f1) = (x1, f1, x0, f0);

   if (d2f < min_d2f)
     d2f = 0;

   while (d2f > min_d2f)
     {
	% f(x + dx) = f(x) + dx f'(x) + dx^2 f''(x)/2
	% Assume f'(x) = 0.  Then:
	%  f(x+dx) = f(x) + dx^2/2 f''(x)
	% Take: x+dx = x1 ==> dx = (x1-x)
	%    ==> f(x1) = f(x) + (x1-x)^2/2 f''
	%  And:  f(x0) = f(x) + (x0-x)^2/2 f''
	%    ==> 2*(f1-f0)/f'' = (x1^2 - 2x*x1 + x^2) - (x0^2 - 2x*x0 + x^2)
	%    ==>  2/f''*(f1-f0) = (x1^2-x0^2) -2x*(x1-x0)
	%                       = (x1-x0)*((x1+x0) - 2x)
	% ==> x = (x1+x0)/2 - (f1-f0)/(x1-x0)/f''
	x2 = 0.5*(x1+x0) - (f1-f0)/(x1-x0)/d2f;
	if (x2 < min_x) x2 = min_x;
	if (x2 > max_x) x2 = max_x;
	if ((abs (x2-x1) <= abserr) || (abs(x2 - x0) <= abserr))
	  {
	     if (xmin == 0.0) xmin = x2;
	     return xmin, fmin, d2f;
	  }

	f2 = (@func)(x2, func_data); sinfo.nfe++; max_linemin_nfe--;
	if (f2 <= fmin)
	  {
	     (xmin, fmin) = (x2, f2);
	     % Now use these to compute the new estimate of f''
	     % f(x) = a*x^2 + b*x + c;  f''(x) = 2*a = d2f
	     % f0 = a*x0^2 + b*x0 + c
	     % f1 = a*x1^2 + b*x1 + c
	     % f2 = a*x2^2 + b*x2 + c
	     % f1-f0 = a*(x1^2-x0^2) + b*(x1-x0)
	     %       = (x1-x0)*(a*(x1+x0) + b)
	     % f2-f1 = (x2-x1)*(a*(x2+x1) + b)
	     % ==>
	     %   b = (f2-f1)/(x2-x1) - a*(x2+x1) = (f1-f0)/(x1-x0) - a*(x1+x0)
	     % ==>
	     %      (f2-f1)/(x2-x1) - (f1-f0)/(x1-x0) = a*(x2-x0)
	     % ==>
	     %        a = ((f2-f1)/(x2-x1) - (f1-f0)/(x1-x0))/(x2-x0)
	     d2f = 2.0*((f2-f1)/(x2-x1) - (f1-f0)/(x1-x0))/(x2-x0);
	     if (isnan(d2f) || isinf(d2f) || (d2f < min_d2f))
	       d2f = 0;
	     return xmin, fmin, d2f;
	  }
	if (max_linemin_nfe == 0)
	  {
	     if (sinfo.verbose)
	       () = fprintf (stderr, "linemin: failed to converge\n");
	     return xmin, fmin,  d2f;
	  }

	% Failed to achieve a better min.  We have 3 possibilities:
	%  x0<=x1<=x2, x0<=x2<=x1, x2<=x0<=x1
	if (x2<x1)
	  {
	     if (x2<x0)
	       (x0,x1,x2,f0,f1,f2) = (x2,x0,x1,f2,f0,f1);
	     else
	       (x1,x2,f1,f2) = (x2,x1,f2,f1);
	  }
	if (f0 > f1 < f2)
	  {
	     d2f = 2.0*((f2-f1)/(x2-x1) - (f1-f0)/(x1-x0))/(x2-x0);
	     continue;
	  }
	% Discard the highest endpoint
	if (f2 < f0)
	  (x0, x1, f0, f1) = (x1, x2, f1, f2);
	break;
     }

   if (f1 < f0)
     {
	x2 = x1 + (x1-x0);
     }
   else
     {
	x2 = x0 - (x1-x0);
     }

   variable bisect_ok = 0;
   forever
     {
	if (x2 > max_x)
	  x2 = 0.5*(max_x + x1);
	if (x2 < min_x)
	  x2 = 0.5*(min_x + x0);

	if ((abs(x1-x2) < abserr) || (abs(x0-x2) < abserr))
	  {
	     if (xmin == 0.0) xmin = x2;
	     return xmin, fmin, 0.0;
	  }

	f2 = (@func)(x2, func_data); sinfo.nfe++; max_linemin_nfe--;
	if (f2 <= fmin)
	  {
	     (xmin, fmin) = (x2, f2);
	     if (d2f != 0) break;
	  }
	if (max_linemin_nfe == 0)
	  {
	     if (sinfo.verbose)
	       () = fprintf (stderr, "linemin: failed to converge\n");
	     break;
	  }

	% At this point, we will have one of the following orderings:
	%  x0<x1<x2,  x0<x2<x1,  x2<x0<x1

	if (x2 < x1)
	  {
	     if (x2 < x0)
	       {
		  % We get here if (x2<x0<x1). 
		  % If d2f!=0, then try to bisect if f2>f0<f1
		  (x0, x1, x2, f0, f1, f2) = (x2, x0, x1, f2, f0, f1);
		  if (bisect_ok && (d2f != 0) && (fmin == f1))
		    {
		       x2 = 0.382*(x0+x1); %  golden section
		       bisect_ok = 0;
		       continue;
		    }
	       }
	     else
	       (x1,x2,f1,f2) = (x2,x1,f2,f1);
	  }
	else if (bisect_ok && (d2f != 0) && (fmin == f1))
	  {
	     x2 = 0.382 * (x1 + x2);
	     bisect_ok = 0;
	     continue;
	  }
	   
	d2f = 2.0*((f2-f1)/(x2-x1) - (f1-f0)/(x1-x0))/(x2-x0);

	if (isnan(d2f) || (d2f < min_d2f) || isinf(d2f))
	  {
	     if (f2 < f0)
	       {
		  (x0, f0, x1, f1) = (x1, f1, x2, f2);
		  x2 = x1 + (x1-x0);
	       }
	     else
	       x2 = x0 - (x1-x0);
	     d2f = 0;
	     continue;
	  }
	
	% It looks like a parabola may work.  The min is at
	% x = -b/2a.  Let's find b:
	%   f0 = a*x0^2 + b*x0 + c
	%   f1 = a*x1^2 + b*x1 + c
	%   f2 = a*x2^2 + b*x2 + c
	% f1-f0 = a*(x1^2-x0^2) + b*(x1-x0) = (x1-x0)(a*(x1+x0)+b)
	% ==> b = (f1-f0)/(x1-x0) - a*(x1+x0)
	% ==> b = (f2-f1)/(x2-x1) - a*(x2+x1)
	% It is probably better to use a combination of both forms, but since
	% this is an approximate calculation, use just:
	%
	%   x = -b/2a = 0.5*(x0+x1) - (f1-f0)/(x1-x0)/2a)
	
	% We have: f0<f1<f2 or f0>f1>f2 or f0>f1<f2.  Throw away the largest
	% point.
	if (f0 > f2)
	  (x0,x1,f0,f1) = (x1,x2,f1,f2);

	x2 = 0.5*(x0+x1) - ((f1-f0)/(x1-x0))/d2f;
	bisect_ok = 1;
     }
   
   d2f = 2.0*((f2-f1)/(x2-x1) - (f1-f0)/(x1-x0))/(x2-x0);
   if (isnan(d2f) || isinf(d2f) || (d2f < min_d2f))
     d2f = 0;
   
   return xmin, fmin, d2f;
}

private define minimize_1d (sinfo, p, x, f, d2f, xi, cd, max_linemin_nfe)
{
   if (x == 0) x = 10*sqrt(_macheps);
   cd.xi = xi;
   cd.p = p;
   % min_p <= p+xi*lambda <= max_p
   % min_p - p <= xi*lambda <= max_p-p
   % Assume xi>0.  Then:
   %   (min_p-p)/xi <= lambda <= (max_p-p)/xi
   % If xi<0, then
   %   -xi*lambda <= (p-min_p), (p-max_p)<=(-xi)*lambda
   %   lambda <= (min_p-p)/xi, (max_p-p)/xi <= lambda
   % Or: (max_p-p)/xi <= lambda <= (min_p-p)/xi
   
   variable dp_min = cd.min_p - p;
   variable dp_max = cd.max_p - p;
   variable min_lambdas = 0.0*xi;
   variable max_lambdas = 0.0*xi;
   variable i, j;
   i = where (xi >= 0.0, &j);
   % Avoid problems with 1/(-0):
   xi[where(xi==0.0)] = 0.0;
   min_lambdas[i] = dp_min[i]/xi[i];
   max_lambdas[i] = dp_max[i]/xi[i];
   min_lambdas[j] = dp_max[j]/xi[j];
   max_lambdas[j] = dp_min[j]/xi[j];

   variable min_l = max(min_lambdas);
   variable max_l = min(max_lambdas);

   if (min_l >= max_l)
     {
	%() = fprintf (stderr, "min_lambda=%S>=max_lambda=%S in powell\n", min_l, max_l);
	return 0.0, f, 0.0;
     }

   % We will be evaluating f(p+lamda*xi), where xi is a unit vector.
   % Consider |p| + x.  Obviously, we want to compute the function
   % for values of x where the computed value |p|+x is different from
   % |p|.  We have: |p| (+) x = (|p|+x)(1+d), where |d|<=2eps.
   % Assuming x > 0, Consider  = (|p|(+)x) - |p|.  We would like this to
   % be as close to x as possible.  This means that p|d| << x.  Here,
   % I arbitrarily choose p|d|<x/2 ==> x > 2*p|d| = 4*p*eps
   variable norm_p = hypot (p);
   variable abserr = 4.0*(_macheps*norm_p);
   
   % Also, at the minimum x=x0, and 
   %   f(x) = f(x0) + (x-x0)^2/2 f''(x0) + ...
   % Assume f(x0) != 0.
   %        = f(x0) + Df
   %        = f(x0)(1 + Df/f0)
   % For f(x) to differ from f(x0), we must have |Df/f0|>eps.  
   % Or: (x-x0)^2/2 f'' > f0*eps ==> x-x0 > sqrt(eps)*sqrt(2*|f0|/f'').
   if (d2f > _macheps)
     abserr += sqrt(_macheps*2.0*f/d2f);

   if (10*abs(x) < abserr)
     x = 10*sign(x)*abserr;

   return linemin (sinfo, f, x, d2f, min_l, max_l, abserr, &call_func, cd, max_linemin_nfe);
}

private define powell ()
{
   if (_NARGS != 6)
     {
	_pop_n (_NARGS);
	powell_usage ();
     }

   variable sinfo, in_parms, min_parms, max_parms, func, func_data;
   (sinfo, in_parms, min_parms, max_parms, func, func_data) = ();

   min_parms = double(min_parms);
   max_parms = double(max_parms);

   ifnot (all(min_parms <= max_parms))
     throw InvalidParmError, "One or more minparms > maxparms";

   variable i = (min_parms <= in_parms <= max_parms);
   variable dsize = max_parms - min_parms;
   ifnot (all(i))
     {
	() = fprintf (stderr, "*** WARNING: min_powell called with initial parameters out of min/max range\n");
	() = fprintf (stderr, "*** WARNING: Randomly choosing a point within the allowed range\n");
	i = wherenot (i);
	in_parms[i] = min_parms[i] + rand_uniform(i) * dsize[i];
     }
   
   if (sinfo == NULL)
     sinfo = powell_new ();
   override_defaults (sinfo, __qualifiers());
   variable verbose = sinfo.verbose;
   
   variable cd = struct
     {
	p, xi, 
	func = func,
	func_data = func_data, 
	min_p = min_parms,
	max_p = max_parms,
     };

   variable num_parms = length(in_parms);
   variable p_0 = @in_parms;
   variable f_0 = (@func)(p_0, func_data); sinfo.nfe++;

   variable xi, lambda, xis = Array_Type[num_parms], 
     d2fs = Double_Type[num_parms], lambdas = Double_Type[num_parms];
   _for i (0, num_parms-1, 1)
     {
	xi = Double_Type[num_parms];
	xi[i] = 1.0;
	xis[i] = xi;
	lambda = 0.01*p_0[i];
	if (lambda == 0)
	  lambda = 0.01;
	lambdas[i] = lambda;
     }
   variable domin=0;
   variable delta = 1.0;
   variable max_linemin_nfe = 10;
   variable max_noprogress = 2;
   forever
     {
	variable f, d2f;
	variable p_n = @p_0;
	variable f_n = f_0;
	variable m = 0;

	if (sinfo.nfe >= sinfo.maxnfe)
	  {
	     sinfo.converged = 0;
	     break;
	  }

	if (verbose)
	  () = fprintf (stdout, "nfe=%d: f=%S @ [%S,..]\n",
			sinfo.nfe, f_0, p_0[0]);

	% Here I am using Zangwill's prescription.
	% Step (i): Calculate \lambda_r to minimize f(p_r) where
	%           p_r = p_{r-1} + \lambda_r\xi_r
	_for i (0, num_parms-1, 1)
	  {
	     xi = xis[i];

	     (lambda, f, d2f) = minimize_1d (sinfo, p_n, 0.5*lambdas[i], f_n, d2fs[i], xi, cd, max_linemin_nfe);
	     if (lambda != 0)
	       {
		  p_n += lambda * xi;
		  d2fs[i] = d2f;
	       }
	     lambdas[i] = lambda;
	     f_n = f;
	  }
	%max_linemin_nfe = 5;

	if ((f < f_0) && (NULL != sinfo.new_min_hook))
	  {
	     (@sinfo.new_min_hook)(p_n, f,
				   __push_list(sinfo.new_min_hook_args));
	  }

	if (feqs (f,f_0,sinfo.reldiff, sinfo.absdiff))
	  {
	     f_0 = f;
	     p_0 = p_n;
	     if (max_noprogress)
	       {
		  max_noprogress--;
		  continue;
	       }
	     sinfo.converged = 1;
	     break;
	  }
	max_noprogress = 2;

	  
	xi = p_n-p_0;
	variable alpha = hypot (xi);
	xi /= alpha;
	(lambda, f, d2f) = minimize_1d (sinfo, p_n, 0.5*alpha, f_n, 0, xi, cd, max_linemin_nfe);
	p_n = p_n + lambda * xi;

	variable tmp = abs(lambdas);
	i = wherefirst (tmp == max(tmp));
	variable lambda_max = tmp[i];
	tmp = delta * (lambda_max/alpha);
	if (tmp >= _macheps)
	  {
	     xis[i] = xi;
	     lambdas[i] = lambda;
	     d2fs[i] = d2f;
	     delta = tmp;
	  }

	p_0 = p_n;
	f_0 = f;
     }
   in_parms[*] = p_0;
   return f_0;
}

