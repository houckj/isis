% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing renorm.... ");

variable Norms = [1,4,5,6];
variable Num_Norms = length(Norms);
variable Test_Name;

variable Size = 256;
variable Noise = urand(Size);
variable Noise2 = urand(Size);
variable Noise_Scale = 0.1;

define init_data (def_ref, fit_ref) %{{{
{
   fit_fun ("gauss(1) + poly(1)");
   set_par (1, 100);
   set_par (2, get_par(2), 0, 6, 24);
   set_par (3, get_par(3), 0, 0.0025, 0.25);
   set_par (4, 10);
   set_par (5, 10);
   set_par (6, 10);

   variable lo, hi, x, i;
   (lo, hi) = linear_grid (10,14,Size);
   x = eval_fun (lo, hi);

   foreach ([0:length(x)-1])
     {
	i = ();
	x[i] += Noise_Scale * x[i] * Noise[i];
     }

   () = (@def_ref)(lo, hi, x, x * Noise_Scale * Noise2);
   () = (@fit_ref)();
}

%}}}

define init_counts () %{{{
{
   init_data (&define_counts, &fit_counts);
}

%}}}

define init_flux () %{{{
{
   init_data (&define_flux, &fit_flux);
}

%}}}

define apply_factor (f, pars) %{{{
{
   variable p;
   foreach (pars)
     {
	p = ();
	set_par (p, f*get_par(p));
     }
}

%}}}

define check_non_value (p, p0) %{{{
{
   if (p.freeze != p0.freeze
       or p.tie != p0.tie
       or p.min != p0.min
       or p.max != p0.max
       or p.fun != p0.fun)
     {
        print(p0);
        print(p);
	failed ("restoring initial parameter config: %s", Test_Name);
     }
}

%}}}

define check_value (p, p0) %{{{
{
   if ((p.is_a_norm == 0)
       or (p.is_a_norm == 1 and (p.freeze == 1 or p.tie != NULL)))
     {
	if (p.value != p0.value)
	  {
	     print(p);
	     print(p0);
	     failed ("restoring initial parameter values: %s", Test_Name);
	  }
	return;
     }

   % If we get to here, it must be a variable norm.

   if (abs(1.0 - p.value / p0.value) > 1.e-3)
     {
	print(p);
	print(p0);
	failed ("finding correct normalization: %s", Test_Name);
     }
}

%}}}

define apply_pars (fun, pars, pars0) %{{{
{
   variable i, ids = [1:get_num_pars()];
   variable p, p0;

   foreach (ids)
     {
	i = ();
	p = pars[i-1];
	p0 = pars0[i-1];

	(@fun) (p, p0);
     }
}

%}}}

define recover (pars0) %{{{
{
   variable pars = get_params();
   apply_pars (&check_non_value, pars, pars0);
   apply_pars (&check_value, pars, pars0);
}

%}}}

define check_renorm (init_ref, renorm_ref) %{{{
{
   (@init_ref) ();

   variable pars, pars0, initial_pars;
   initial_pars = get_params();

   Test_Name = "all norms vary";
   pars0 = get_params();
   apply_factor (0.5, Norms);
   () = (@renorm_ref)();
   recover (pars0);
   set_params (initial_pars);

   Test_Name = "some norms frozen";
   variable f, o, x;
   f = Norms[0];
   o = Norms[[1:Num_Norms-1]];
   freeze (f);
   pars0 = get_params();
   apply_factor (2.0, o);
   () = (@renorm_ref)();
   recover (pars0);
   set_params (initial_pars);

   Test_Name = "parameter function";
   variable fun = "_par(1) + 0.0*PI";
   set_par_fun (f, fun);
   () = eval_fun ([1,2], [2,3]);   % so param function is evaluated
   pars0 = get_params();
   apply_factor (0.75, o);
   () = (@renorm_ref)();
   recover (pars0);
   set_params (initial_pars);

   Test_Name = "everything frozen";
   freeze ([1:get_num_pars()]);
   pars0 = get_params();
   if (-1 != (@renorm_ref)())
     failed ("handling all params frozen");
   recover (pars0);
   set_params (initial_pars);
}

%}}}

check_renorm (&init_counts, &renorm_counts);
delete_data (all_data);
check_renorm (&init_flux, &renorm_flux);

msg ("ok\n");
