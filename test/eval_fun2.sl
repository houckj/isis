% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing eval_fun2.... ");

variable lo, hi, y, num=100;
(lo, hi) = linear_grid (1,20,num);

% S-Lang user-defined function:
define foo_function (l,h,p)
{
   return p[0] * ones(length(l));
}
add_slang_function ("foo", &foo_function, "bar");
y = eval_fun2 ("foo", lo, hi, [20]);
if (abs(1.0 - sum(y)/(20*num)) > 1.e-4)
  failed ("S-Lang function failed");
y = eval_fun2 (fitfun_handle("foo"), lo, hi, [20]);
if (abs(1.0 - sum(y)/(20*num)) > 1.e-4)
  failed ("S-Lang function failed (handle)");

% Built-in C function:
variable x0=12.0, s=0.25;
variable gauss_pars = [1.0, x0, s];

(lo, hi) = linear_grid (x0-6*s, x0+6*s, 200);
y = eval_fun2 ("gauss", lo, hi, gauss_pars);
if (abs(1.0 - sum(y)) > 1.e-4)
  failed ("C function failed");
y = eval_fun2 (fitfun_handle ("gauss"), lo, hi, gauss_pars);
if (abs(1.0 - sum(y)) > 1.e-4)
  failed ("C function failed (handle)");

% S-Lang operator function:
define operator (l,h,p,f)
{
   return p[0] * mean(f) * ones(length(l)) + p[1];
}
add_slang_function ("op", &operator, ["fac", "offset"]);
set_function_category ("op", ISIS_FUN_OPERATOR);
y = eval_fun2 ("op", lo, hi, [12.3, -3], ones(length(lo)));
if (abs(1.0 - sum(y)/(9.3*length(lo))) > 1.e-4)
  failed ("operator function failed");
y = eval_fun2 (fitfun_handle("op"), lo, hi, [12.3, -3], ones(length(lo)));
if (abs(1.0 - sum(y)/(9.3*length(lo))) > 1.e-4)
  failed ("operator function failed (handle)");

msg ("ok\n");
