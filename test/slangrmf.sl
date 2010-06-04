() = evalfile ("inc.sl");
msg ("testing slangrmf....");

private variable Tol = 1.e-6;

private define rmf_test (bin_lo, bin_hi, x, parms)
{
   % Depending on the 'grid' qualifier passed to load_slang_rmf
   % x, bin_lo, bin_hi may be either energy [keV] or wavelength [Angstrom].

   variable
     shape = parms[0],
     resolution = parms[1];

   variable rmf = Double_Type[length(bin_lo)];

   if (x < bin_lo[0] || bin_hi[-1] <= x)
     return rmf;

   variable i = where (x * (1.0 - resolution/2.0) <= bin_lo
                       and bin_hi < x * (1.0 + resolution/2.0));
   variable n = length(i);
   if (n == 0)
     {
        n = 1;
        i = bsearch_hist (x, bin_lo, bin_hi);
     }

   variable k;

   switch (shape)
     {
      case "square":
        rmf[i] = 1.0 / n;
     }
     {
      case "wedge":
        %     energy grid => high end of wedge at large energies
        % wavelength grid => high end of wedge at long wavelengths
        rmf[i[0]] = 1.0/n;
        foreach k (i[[1:]])
          {
             rmf[k] = rmf[k-1] + 1.0/n;
          }
     }
     {
        % default:
        throw ApplicationError, "Unsupported rmf shape = $shape"$;
     }

   return rmf/sum(rmf);
}

define test_norm ()
{
   delete_data (1);
   delete_rmf (1);

   () = load_slang_rmf (&rmf_test,
                        linear_grid (0.1, 10, 1024),
                        linear_grid (0.1, 10, 2048);; __qualifiers);

   assign_rmf (1, 1);
   () = eval_counts();

   variable m = _A(get_model_counts (1));

   if (any (m.value < 0))
     throw ApplicationError, "negative model values appeared ???";

   ifnot (feqs (sum(m.value), get_par ("delta(1).norm"), Tol))
     {
        throw ApplicationError, "norm error exceeds Tol=$Tol"$;
     }

   return m;
}

fit_fun ("delta(1)");
set_par (1, 10);
set_par (2, _A(6.4));

variable res = 0.2;

% Check profile norms

variable s_en = test_norm ( ; parms={"square", res});
variable s_wv = test_norm ( ; parms={"square", res}, grid="wv");

variable w_en = test_norm ( ; parms={"wedge", res});
variable w_wv = test_norm ( ; parms={"wedge", res}, grid="wv");

% Check profile shapes

variable high_en = w_en.value [ wherelast (w_en.value > 0) ];
variable   lo_en = w_en.value [ wherefirst (w_en.value > 0) ];
ifnot (feqs (high_en/lo_en, howmany(w_en.value > 0), Tol))
   throw ApplicationError, "profile error:  Tol=$Tol"$;

variable high_wv = w_wv.value [ wherefirst (w_wv.value > 0) ];
variable   lo_wv = w_wv.value [ wherelast (w_wv.value > 0) ];
ifnot (feqs (high_wv/lo_wv, howmany(w_wv.value > 0), Tol))
   throw ApplicationError, "profile error:  Tol=$Tol"$;

msg ("ok\n");
