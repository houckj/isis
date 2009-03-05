% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing stat.... ");

define foo_fit(l,h,p)
{
   return p[0]*ones(length(l));
}
add_slang_function ("foo", "a");
fit_fun ("foo");
set_par ("foo(1).a", 0);

variable Num_Bins = 2000;
variable Back = 2 * ones(Num_Bins);
variable Src = 8 * ones(Num_Bins);

define create_one_dataset (x, b)
{
   variable lo, hi, s;
   lo = [1:Num_Bins];
   hi = lo + 1;

   s = Src * (1.0 + 0.1*grand()) + prand(5, Num_Bins);

   variable id = define_counts (lo, hi, s, sqrt(s));
   if (b != NULL) () = _define_back (id, b);

   rebin_data (id, sum(s)/(Num_Bins/2));
}

define compute_sigma (n)
{
   variable sigma = Double_Type[length(n)];

   variable j = where (n < 1.0);
   if (length(j) > 0)
     sigma[j] = 1.0;

   variable i = where (n >= 1.0);
   if (length(i) > 0)
     sigma[i] = sqrt(n[i]);

   return sigma;
}

define compute_chi_square ()
{
   variable d, b, m, sigma, dchi2;

   variable list = all_data();
   variable chisqr = 0.0;

   foreach (list)
     {
        variable id = ();
        d = get_data_counts(id);
        b = get_back (id);
        m = get_model_counts(id);

        if (b == NULL)
          sigma = d.err;
        else
          sigma = hypot (d.err, compute_sigma(b));

        dchi2 = (d.value - m.value)^2 / sigma^2;
        chisqr += sum (dchi2);
     }

   return chisqr;
}

define compute_combined_chi_square ()
{
   variable d, m, e, chisqr;

   d = Double_Type[length(get_data_counts(1).value)];
   m = @d;
   e = @d;

   variable list = all_data();
   foreach (list)
     {
        variable id = ();
        variable _d, _m, _b, _e;
        _d = get_data_counts(id);
        _m = get_model_counts(id);
        _b = get_back (id);

        if (_b == NULL)
          _e = _d.err;
        else
          _e = hypot (_d.err, compute_sigma(_b));

        d += _d.value;
        m += _m.value;
        e += _e^2;             % add in quadrature
     }

   chisqr = sum ((d - m)^2 / e);

   return chisqr;
}

define compare_chisqr (fun, msg)
{
   variable info;
   Fit_Verbose=-1;
   () = eval_counts (&info);

   variable expected = (@fun)();

   if (abs(1.0 - info.statistic/expected) > 1.e-7)
     {
        failed ("%s\n expected stat=%15.7e, got stat=%15.7e",
                msg,
                expected,
                info.statistic);
     }
}

define main ()
{
   create_one_dataset (0, NULL);
   compare_chisqr(&compute_chi_square, "1 dataset, no background");
   delete_data(all_data);

   create_one_dataset (0, Back);
   compare_chisqr(&compute_chi_square, "1 dataset, with background");
   delete_data(all_data);

   create_one_dataset (0, Back);
   create_one_dataset (0, Back);
   create_one_dataset (0, Back);
   compare_chisqr(&compute_chi_square, "3 datasets, with background(s)   [not combined]");
   delete_data(all_data);

   create_one_dataset (0, Back);
   create_one_dataset (0, Back);
   create_one_dataset (0, Back);
   match_dataset_grids (all_data);
   () = combine_datasets (all_data);
   compare_chisqr(&compute_combined_chi_square, "3 datasets, with background(s)   [combined]");
   delete_data(all_data);
}

main();

msg ("ok\n");

