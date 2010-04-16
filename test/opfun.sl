%-*- slang -*-
() = evalfile ("inc.sl");
msg ("testing operator functions.... ");

variable n = 100000;
variable x = grand(n);
x += abs(2*min(x));

variable lo, hi;
(lo,hi) = linear_grid (min(x), max(x), 100);
variable val = histogram (x, lo, hi);
() = define_counts (lo, hi, val, sqrt(val));
() = define_counts (lo, hi, val, sqrt(val));

% wlshift is an operator-function (somewhat like xspec's
% convolution models).  It wavelength-shifts the result of
% any model by an amount specific to each dataset.

define wlshift_fit (l,h,p,val)
{
   variable offset;

   switch (Isis_Active_Dataset)
     {
      % dataset 1 gets shifted by "offset1" = p[0]
      case 1:
	offset = p[0];
     }
     {
      % dataset 2 gets shifted by "offset2" = p[1]
      case 2:
	offset = p[1];
     }
     {
	% default
	error("*** wlshift: unsupported dataset index");
	return -1;
     }

   return rebin (l, h, l+offset, h+offset, val);
}
add_slang_function ("wlshift", ["offset1", "offset2"]);
set_function_category ("wlshift", ISIS_FUN_OPERATOR);

fit_fun ("gauss(1)");
set_par ("gauss(1).area", 1.e4);
set_par ("gauss(1).center", 8);
set_par ("gauss(1).sigma", 1);

() = renorm_counts();

variable info;
() = fit_counts (&info);
if (abs(1.0 - 207.941/info.statistic) > 1.e-3)
  failed ("opfun:  chisqr=%g  expected 207.941", info.statistic);

% apply the wavelength shift to a gaussian line

fit_fun ("wlshift (1,gauss(1))");
set_par ("wlshift(1).offset1", -0.1);
set_par ("wlshift(1).offset2", +0.1);

() = eval_counts (&info);
if (abs(1.0 - 2226.85/info.statistic) > 1.e-3)
  failed ("opfun:  chisqr=%g  expected 2226.85", info.statistic);

msg ("ok\n");
