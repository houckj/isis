% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing sys_err.... ");
%
private variable
  Sef_Test = 0.05,
  Value = 1.e5,
  Num = 100;

private define fake (n, value)
{
   variable lo, hi, c;
   (lo, hi) = linear_grid (1,20,n);
   c = prand (value * (1.0 + Sef_Test), n);
   return define_counts (lo, hi, c, sqrt(c));
}

private define check_sef (id)
{
   variable s = get_sys_err_frac (id);
   variable i = fneqs(s, Sef_Test, 1.e-4);
   if (any(i))
     {
        print(s[i]);
        throw IsisError, "sys_err_frac error, values should = $Sef_Test"$;
     }
}

variable id = fake (Num, Value);
variable s, sef_array = ones(Num)*Sef_Test;

% Exercise set/unset sys_err_frac:
s = get_sys_err_frac (id);
if (length(s) != 0)
  throw IsisError, "error retrieving sys_err_frac";

set_sys_err_frac (id, sef_array);
check_sef (id);

set_sys_err_frac (id, NULL);
s = get_sys_err_frac (id);
if (length(s) != 0)
  throw IsisError, "error retrieving sys_err_frac";

set_sys_err_frac (id, sef_array);

% Test grouping
group_data (1, 2);
check_sef (id);

group_data(1,1);
check_sef (id);

rebin_data (id, 100.0);
check_sef (id);

group_data (1,1);

% Test rebinning
variable new_lo, new_hi;
(new_lo, new_hi) = linear_grid (1,20, 2*Num);
rebin_dataset (1, new_lo, new_hi);
check_sef (id);

% Test fit-statistic
define const_fit(l,h,p)
{
   return p[0]*ones(length(l));
}
add_slang_function ("const", "value");
fit_fun ("const");
set_par ("const(1).value", Value);

define test_stat (id, sys_err)
{
   variable q;
   () = eval_counts (&q);

   variable
     d = get_data_counts(id),
     m = get_model_counts(id);

   variable sigma = sqrt(d.err^2 + (sys_err*d.value)^2);
   variable s = sum (sqr((d.value-m.value)/sigma));

   if (fneqs(q.statistic, s, 1.e-3))
     {
        variable msg = sprintf ("Error in statistic computation:  stat= %g, should be $s"$,
                                q.statistic);
        throw IsisError, msg;
     }
}

s = get_sys_err_frac (id);

set_sys_err_frac (id, NULL);
test_stat (id, 0);

set_sys_err_frac (id, s);
test_stat (id, s);

set_sys_err_frac (id, NULL);
test_stat (id, 0);

msg ("ok\n");
