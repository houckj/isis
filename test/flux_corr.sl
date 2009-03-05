% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing flux_corr.... ");

fit_fun ("blackbody(1)");

assign_rsp (load_rmf ("data/acismeg1D1999-07-22rmfN0002.fits"),
            load_arf ("data/acisf01318_000N001MEG_-1_garf.fits"),
            load_data ("data/acisf01318N003_pha2.fits", 9));

variable xmin, xmax;
xmin = 1.0;
xmax = 25.0;

define doit (binfac, op)
{
   group_data (1, binfac);
   flux_corr(1);
   xnotice(1, xmin, xmax);
   () = (@op);

   variable df, cmf;
   df = get_data_flux(1);
   cmf = get_convolved_model_flux(1);
   
   variable i = where (xmin < df.bin_lo and df.bin_hi < xmax);

   variable tol = 1.e-3;
   variable df_sum = sum(df.value[i]);
   variable df_sum0 = 0.0527881;
   if (abs (df_sum/df_sum0 - 1.0) > tol)
       {
          failed ("[1]flux-corrected data sum = %g, should be %g",
                  df_sum, df_sum0);
       }

   variable cmf_sum = sum(cmf.value[i]);
   variable cmf_sum0 = 0.0200908;
   if (abs (cmf_sum/cmf_sum0 - 1.0) > tol)
     {
          failed ("[2]convolved model flux sum = %g, should be %g",
                  cmf_sum, cmf_sum0);
     }
   
   variable dc = get_data_counts(1);
   variable cerr, ferr;
   i = where (dc.value > 0);
   cerr = dc.err[i]/dc.value[i];
   ferr = df.err[i]/df.value[i];
   if (any (abs(ferr/cerr - 1) > 1.e-3))
     {
        failed ("[3] flux-corrected data uncertainties don't match counts data");
     }
}

doit(16, &fit_flux);
doit(4, &eval_flux);
doit(0, &eval_flux);

msg ("ok\n");
