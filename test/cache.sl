% -*- mode: SLang; mode: fold =*=
() = evalfile ("inc.sl");
msg ("testing cacheing/aliasing.... ");

variable y0=1.0, y1=20.0;
variable lo, hi, nbins=2000;
(lo,hi) = linear_grid (y0, y1, nbins);

variable add = "cnst";
define cnst_fit(l,h,p){ return ones(length(l))*p[0]; }
add_slang_function (add, "a");
variable a_value = 1.0;
variable add_test_value = a_value * nbins;
variable
  mul_test_value = add_test_value,
  mul_cache_test_value = mul_test_value / 2.0;

variable add_cache = cache_fun (add, lo, hi; suffix="X");

fit_fun (add_cache);
set_par ("$add_cache(1).a"$, a_value);
variable e, excp_occurred = 0;
try (e)
{
   () = eval_fun (0.5, 2);
}
catch AnyError:
{
   excp_occurred = 1;
}
if (excp_occurred == 0)
   failed ("expected DomainError -- got but no exception was caught");
else if (e.error != DomainError)
   failed ("expected DomainError -- got %S", e.descr);

fit_fun (add);
set_par ("$add(1).a"$, a_value);
variable f_add = eval_fun (lo, hi);

fit_fun (add_cache);
set_par ("$add_cache(1).a"$, a_value);
variable f_add_cache = eval_fun (linear_grid(y0, y1, nbins/2));

if (sum(f_add) != add_test_value)
  failed ("add sum is wrong...");
if (sum(f_add_cache) != add_test_value)
  failed ("add_cache sum is wrong...");

variable mul = add;
variable mul_cache = cache_fun (mul, lo, hi; mult);

fit_fun (mul);
set_par ("$mul(1).a"$, a_value);
variable f_mul = eval_fun (lo,hi);

fit_fun (mul_cache);
set_par ("$mul_cache(1).a"$, a_value);
variable f_mul_cache = eval_fun (linear_grid(y0, y1, nbins/2));

if (sum(f_mul) != mul_test_value)
  failed ("mul sum is wrong...");
if (sum(f_mul_cache) != mul_cache_test_value)
  failed ("mul_cache sum is wrong...");

variable t = struct
{
   name = ["par1", "par2", "par3"],
   value = [100, 7.0, 0.02],
   min = [0, 1.0, 0.01],
   max = [1.e10, 10.0, 1.0],
   freeze = [0, 1, 0]
};

alias_fun ("gauss", "gauss_alias" ;
           names=t.name, values=t.value, min=t.min, max=t.max, freeze=t.freeze);
variable info = get_fun_info ("gauss_alias");
if (any(info.name != t.name)
    || any(info.value != t.value)
    || any(info.min != t.min)
    || any(info.max != t.max)
    || any(info.freeze != t.freeze))
  failed ("alias failed");

msg ("ok\n");
