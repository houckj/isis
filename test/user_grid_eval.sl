%-*- slang -*-
() = evalfile ("inc.sl");
msg ("testing user-defined model evaluation grid.... ");

variable n = 100000;
variable x = grand(n);
x += abs(2*min(x));

variable lo, hi;
(lo,hi) = linear_grid (min(x), max(x), 100);
variable val = histogram (x, lo, hi);
() = define_counts (lo, hi, val, sqrt(val));

fit_fun ("gauss(1)");
set_par(1,1.e4);
set_par(2,8.2);
set_par(3,1);
() = renorm_counts;
() = fit_counts;

variable Num_Model_Eval_Bins;

define model_eval_grid_hook (id, s)
{
   (s.bin_lo, s.bin_hi) = linear_grid (0.5*lo[0], hi[-1]*2, Num_Model_Eval_Bins);
   return s;
}
set_eval_grid_method (USER_GRID, all_data(), &model_eval_grid_hook);

define do_test (num, chisqr)
{
   Num_Model_Eval_Bins = num;

   variable info;
   () = eval_counts (&info);
#iffalse
   oplot_model_flux;
   plot_pause;
#endif
   if (abs (1.0 - info.statistic/chisqr) > 0.001)
     failed ("user_grid_eval: num=%d  stat=%g  (should be %g)",
               num, info.statistic, chisqr);
}

variable sizes         = [12,      25,      50,   100,     200];
variable chisqr_values = [63155.7, 8189.29, 1522.61, 271.688, 103.378];

#iffalse
Num_Model_Eval_Bins = 200;
() = eval_counts;
plot_model_flux;
plot_pause;
#endif

array_map (Void_Type, &do_test, sizes, chisqr_values);

msg ("ok\n");
