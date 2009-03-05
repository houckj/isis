% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing fit constraints.... ");

define add_fake_dataset ()
{
   variable x = grand(100000);
   x += 2*abs(min(x));

   variable n, lo, hi;
   n = 1024;
   (lo, hi) = linear_grid (min(x), max(x), n);

   variable nx = histogram (x, lo, hi);
   return define_counts (_A(lo, hi), reverse(nx), sqrt(reverse(nx)));
}

define example_constraint (stat, pars)
{
   variable value = get_par ("constraint(1).a");
   return stat + value;
}

define no_params_constraint (stat, pars)
{
   return stat;
}

define  main ()
{
   variable info1, info2;
   () = add_fake_dataset();

   () = renorm_counts;
   fit_fun ("egauss(1)");
   set_par ("egauss(1).sigma", 1.0, 0, 0, 10);
   () = fit_counts (&info1);

   variable value = 1.e4;
   set_fit_constraint (&example_constraint, ["a"]);
   set_par ("constraint(1)", value);

   () = eval_counts (&info2);

   if (info2.statistic != info1.statistic + value)
     failed ("fit constraint value not applied");

   set_fit_constraint (NULL);
   () = eval_counts (&info2);

   if (info2.statistic != info1.statistic)
     failed ("fit constraint value wrongly applied");
   
   set_fit_constraint (&no_params_constraint);
   () = eval_counts (&info2);
}

main ();

msg ("ok\n");

