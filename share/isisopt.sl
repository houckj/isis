% -*- mode: slang; mode: fold; -*-

require ("structfuns");

private define help_options (method, qualifiers)
{
   variable str = "Options for $method:\n"$;
   foreach (get_struct_field_names (qualifiers))
     {
	variable name = ();
	if (0 == strncmp (name, "help_", 5))
	  continue;

	variable val = get_struct_field (qualifiers, name);
	if (val == NULL)
	  str = sprintf ("%s %s\n", str, name);
	else
	  str = sprintf ("%s %s=%S\n", str, name, val);
	name = "help_" + name;
	if (struct_field_exists (qualifiers, name))
	  str = sprintf ("%s  [%s]\n", str, get_struct_field(qualifiers, name));
     }
   usage (str);
}

private define handle_options (options, method, qualifiers)
{
   if (struct_field_exists (options, "help"))
     help_options (method, qualifiers);

   foreach (get_struct_field_names (options))
     {
	variable name = ();
	if (struct_field_exists (qualifiers, name))
	  set_struct_field (qualifiers, name, get_struct_field (options, name));
     }
}

private define call_fun (p, cd)
{
   p = @p;
   variable pmin = cd[1], pmax = cd[2];
   variable i = where (p<pmin);
   p[i] = pmin[i];
   i = where (p>pmax);
   p[i] = pmax[i];
   __eval_stat (cd[0], p);
}

%{{{ Powell Interface
autoload ("powell_new", "powell");
private variable Powell_Qualifiers = struct
{
   verbose = 0,
   reldiff = 1e-6,
   absdiff = 1e-8,
   maxnfe = 500,
};

private define set_powell_options (options)
{
   handle_options (options, "powell", Powell_Qualifiers);
}
private define powellfit (mmt, p, pmin, pmax)
{
   variable cd = {mmt, pmin, pmax};
   variable obj = powell_new ();

   () = obj.optimize (p, pmin, pmax, &call_fun, cd ;; Powell_Qualifiers);
   return p;
}
register_slang_optimizer ("powell", &powellfit; set_options=&set_powell_options);

%}}}

%{{{ Simplex Interface
autoload ("simplex_new", "simplex");
private variable Simplex_Qualifiers = struct
{
   verbose = 0,
   reldiff = 1e-4,
   absdiff = 1e-8,
   maxnfe = 1000,
   reflect=1, contract=0.5, stretch=2.0, shrink=0.5, scale=0.1,
};

private define set_simplex_options (options)
{
   handle_options (options, "simplex", Simplex_Qualifiers);
}
private define simplexfit (mmt, p, pmin, pmax)
{
   variable cd = {mmt, pmin, pmax};
   variable obj = simplex_new ();

   () = obj.optimize (p, pmin, pmax, &call_fun, cd ;; Simplex_Qualifiers);
   return p;
}
register_slang_optimizer ("simplex", &simplexfit; set_options=&set_simplex_options);

%}}}

%{{{ Diffevol Interface
autoload ("diffevol_new", "diffevol");
private variable Diffevol_Qualifiers = struct
{
   verbose = 0,
   reldiff = 1e-4,
   absdiff = 1e-8,
   maxnfe = 1000,
   npop=10, help_npop="Population per parameter",
   diffscale=0.7, crossprob=0.5, max_gens=-1,
   strategy="best1bin",
   help_strategy="best1bin|best1exp|rand1exp|rand1bin|randtobest1exp|randtobest1bin",
   bc="rand", help_bc="Boundary Conditions: none|reject|reflect|truncate|rand",
};

private define set_diffevol_options (options)
{
   handle_options (options, "diffevol", Diffevol_Qualifiers);
}
private define diffevolfit (mmt, p, pmin, pmax)
{
   variable cd = {mmt, pmin, pmax};
   variable obj = diffevol_new ();
   variable q = @Diffevol_Qualifiers;
   q.npop *= length(p);
   () = obj.optimize (p, pmin, pmax, &call_fun, cd ;; q);
   return p;
}
register_slang_optimizer ("diffevol", &diffevolfit; set_options=&set_diffevol_options);

%}}}

%{{{ plm Interface
autoload ("plm_new", "plm");
private variable PLM_Qualifiers = struct
{
   tol = 1.e-4,
   lambda = 0.1,
   grow_factor = 10.0,
   shrink_factor = 0.1,
   max_loops = 100,
   delta = 1.e-4,
   verbose = 0,
   num_slaves = -1,
};

private define set_plm_options (options)
{
   handle_options (options, "plm", PLM_Qualifiers);
}
private define plmfit (mmt, p, pmin, pmax)
{
   variable obj = plm_new ();
   () = obj.optimize (p, pmin, pmax, mmt ;; PLM_Qualifiers);
   return p;
}
register_slang_optimizer ("plm", &plmfit; set_options=&set_plm_options);
%}}}

