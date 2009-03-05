% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing par_fun.... ");

fit_fun ("gauss(1) + gauss(2) + gauss(3)");

set_par ("gauss(1).center", 12.0);
set_par ("gauss(2).center", 1.0);
set_par ("gauss(3).center", 12.0);

% Try simple operations
set_par_fun ("gauss(2).center", "0.5*(gauss(1).center + gauss(3).center)");
if (12.0 != get_par("gauss(2).center"))
  failed ("sum failed");

set_par_fun ("gauss(2).center", "gauss(1).center*gauss(3).center");
if (144.0 != get_par("gauss(2).center"))
  failed ("product failed");

% Try self-referential expression
set_par_fun ("gauss(2).center", "gauss(2).center*10");
set_par ("gauss(2).center", 12.0);
if (120.0 != get_par("gauss(2).center"))
  failed ("self-ref failed");

% Try kernel parameters
variable lo, hi, n;
(lo,hi) = linear_grid(1,20,200);
n = prand(20,200);
() = define_counts(lo,hi,n,sqrt(n));
set_kernel (1, "pileup");
variable _silly=1;
set_par_fun ("gauss(1).center", "pileup(_silly).g0 + 3.5");
if (4.5 != get_par("gauss(1).center"))
  failed ("kernel param failed");

% Try operator function parameters

define nullop_fit (lo,hi,par,fun)
{
   return fun;
}
add_slang_function ("nullop", "a");
set_function_category ("nullop", ISIS_FUN_OPERATOR);
fit_fun ("nullop(1,gauss(1)+gauss(2))");
set_par_fun ("gauss(1).center", NULL);
set_par ("gauss(1).center", 12.0);
set_par ("nullop(1).a", 24.1);
set_par_fun ("gauss(2).center", "gauss(1).center + nullop(1).a");
if (abs(1.0 - get_par("gauss(2).center")/36.1) > 1.e-6)
  failed ("operator param failed");

msg ("ok\n");
