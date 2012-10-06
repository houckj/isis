% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing rebin.... ");

%
%%  First, check the accuracy of the rebinning function
%

variable bin_lo, bin_hi, data, lo, dx;

data = [1:10]; data[*] = 1;
bin_lo = [0:9];
bin_hi = [1:10];

define _try_order (answer, to_lo, to_hi, from_lo, from_hi, from_y) %{{{
{
   variable to_y, s, err;

   to_y = rebin (to_lo, to_hi, from_lo, from_hi, from_y);
   s = sum(to_y);
   err = 1.0 - s / answer;
   if (abs(err) > 1.e-15)
     failed ("rebin");
}

%}}}

define rev() %{{{
{
   variable a, b, c;
   switch (_NARGS)
     {
      case 2:
	(a,b) = ();
	return (reverse(a), reverse(b));
     }
     {
      case 3:
	(a,b,c) = ();
	return (reverse(a), reverse(b), reverse(c));
     }
     {
	% default:
	error ("rev: unsupported number of args");
     }
}

%}}}

define _try_reversals (answer, lo, hi, bin_lo, bin_hi, data) %{{{
{
   _try_order (answer,     lo, hi,      bin_lo, bin_hi, data);
   _try_order (answer, rev(lo, hi),     bin_lo, bin_hi, data);
   _try_order (answer,     lo, hi,  rev(bin_lo, bin_hi, data));
   _try_order (answer, rev(lo, hi), rev(bin_lo, bin_hi, data));
}

%}}}

define try_rebin (lo, hi, answer) %{{{
{
   _try_reversals (answer, lo, hi, bin_lo, bin_hi, data);
   _try_reversals (answer, hi, lo, bin_lo, bin_hi, data);
   _try_reversals (answer, hi, lo, bin_hi, bin_lo, data);
}

%}}}

try_rebin (bin_lo, bin_hi, 10);
try_rebin ([2,3,4], [3,4,5], 3);
try_rebin ([2.5,3.5,4.5], [3.5,4.5,5.5], 3);

dx = 2;
lo = [0:10:dx];  try_rebin (lo, lo+dx, 10);
lo = [-2:12:dx];  try_rebin (lo, lo+dx, 10);
lo = [5:12:dx];  try_rebin (lo, lo+dx, 5);
lo = [-4:6:dx];  try_rebin (lo, lo+dx, 8);
lo = [5.3:12:dx];  try_rebin (lo, lo+dx, 4.7);
lo = [-4.3:7.7:dx];  try_rebin (lo, lo+dx, 7.7);

dx = 1;
lo = [0:10];  try_rebin (lo, lo+dx, 10);
lo = [-2:12];  try_rebin (lo, lo+dx, 10);
lo = [5:12];  try_rebin (lo, lo+dx, 5);
lo = [-4:6];  try_rebin (lo, lo+dx, 7);

try_rebin (0, 10, 10);
try_rebin (-4, 6, 6);
try_rebin (-2, 12, 10);
try_rebin (5, 12, 5);

%
%%  Now, check the rebin error-hook.
%

define init () %{{{
{
   fit_fun ("gauss(1)");
   set_par (1, 100);
   set_par (2, get_par(2), 0, 6, 24);
   set_par (3, 10*get_par(3), 0, 0.0025, 0.25);

   variable lo, hi, cts, i;
   (lo, hi) = linear_grid (10,14,256);
   cts = eval_fun (lo, hi);

   foreach ([0:length(cts)-1])
     {
	i = ();
	cts[i] = prand(cts[i]);
     }

   () = define_counts (lo, hi, cts, sqrt(cts));
   () = fit_counts();
}

%}}}

init ();

define hook (o_cts, o_err, grp) %{{{
{
   variable errs = sqrt(rebin_array (o_err^2, grp));
   errs[*] = 8.0;
   return errs;
}

%}}}

set_rebin_error_hook (1, "hook");
rebin_data (1,7);

variable i, d = get_data_counts (1);
if (length(d.err) != length(where (d.err == 8)))
  failed ("assigning user uncertainties");

set_rebin_error_hook (1, NULL);

define try_broken_hook () %{{{
{
   variable sv = Isis_Verbose;

   if (Isis_Verbose >= _isis->_INFO)
     message (" => Expect an error here:");
   ERROR_BLOCK
     {
	Isis_Verbose = sv;
	_clear_error();
	return;
     }

   Isis_Verbose = -5;
   rebin_data (1,2);
   Isis_Verbose = sv;
   failed ("catching intentional error");
}

%}}}

define returns_null (c,e,g) %{{{
{
   return NULL;
}

%}}}

define returns_zeros (c,e,g) %{{{
{
   variable errs = sqrt(rebin_array (e^2, g));
   errs[*] = 0.0;
   return errs;
}

%}}}

set_rebin_error_hook (1, "returns_null");
try_broken_hook();

% Here, isis will validate the uncertainties and reset
% them to min_stat_err
set_rebin_error_hook (1, "returns_zeros");
set_min_stat_err (1, 1.0);

rebin_data (1,0);
d = get_data_counts (1);
if (length(d.err) != length(where (d.err == 1)))
  failed ("applying min_stat_err");

define random_rebin ()
{
   variable mask = Integer_Type[256];
   variable u = urand(256);
   
   mask[where(u < 1.0/3)] = -1;
   mask[where(1.0/3 < u and u < 2.0/3)] =  0;
   mask[where(2.0/3 < u)] = 1;
   
   rebin_data (1,0);
   variable expect = sum (get_data_counts(1).value[where(mask != 0)]);
   
   rebin_data (1, mask);
   
   variable got = sum(get_data_counts(1).value);
   if (abs(1.0 - got/expect) > 1.e-15)
     failed ("rebin_data");   
}

loop (20)
{
   random_rebin();
}

msg ("ok\n");
