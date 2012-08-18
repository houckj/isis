% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing assign_back.... ");

variable lo, hi, n=1024;
(lo,hi) = linear_grid (1,20,n);

variable t1=1.e5, a1=10.0;
variable t2=1.e6, a2=1.e3;
variable factor = 100.0;

define cnst_fit(l,h,p) { return p[0]*ones(length(l));}
add_slang_function ("cnst", ["factor"]);

fit_fun ("cnst");
set_par ("cnst(1).factor", factor);
variable f = eval_fun (lo,hi);

variable id1 = define_counts (lo,hi,f,sqrt(f));
set_data_exposure (id1, t1);
set_data_backscale (id1, a1);

variable id2 = define_counts (lo,hi,f,sqrt(f));
set_data_exposure (id2, t2);
set_data_backscale (id2, a2);

() = eval_counts;

assign_back (id2, id1);

variable b1 = get_back (id1);
variable bt1 = get_back_exposure (id1);
variable ba1 = get_back_backscale (id1);

if (bt1 != t2)
  throw ApplicationError, "*** assign_back: wrong background exposure time";

if (ba1 != a2)
   throw ApplicationError, "*** assign_back: wrong background area";

variable b2_model_value = (factor * t2);
variable expected_background_value = b2_model_value * (a1/a2) * (t1/t2);
if (any (b1 != expected_background_value))
   throw ApplicationError, "*** assign_back: wrong background value";

unassign_back (id1);
if (get_back (id1) != NULL)
  throw ApplicationError, "*** assign_back: unassign_back failed";

msg ("ok\n");
