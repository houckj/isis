% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing conf_map.... ");

if (orelse
    {NULL == find_library_name ("cfitsio-module.so")}
    {NULL == find_script_name ("fits.sl")})
{
   msg ("skipped\n");
   exit(0);
}

variable Output_File = "calc_confmap.fits";

define init_demo () %{{{
{
   fit_fun ("gauss(1)");
   set_par (1, 100);
   set_par (2, get_par(2), 0, 6, 24);
   set_par (3, get_par(3), 0, 0.0025, 0.25);

   variable lo, hi, cts, i;
   (lo, hi) = linear_grid (10,14,256);
   cts = eval_fun (lo, hi);

   foreach ([0:length(cts)-1])
     {
	i = ();
	cts[i] = prand(cts[i]);
     }

   () = define_counts (lo, hi, cts, sqrt(cts));
   () = fit_counts();
}

%}}}

define map (scale, xo,yo) %{{{
{
   variable px, py;

   % X-param
   px = conf_grid (2, 11.96 + xo, 12.04 + xo, int(40 * scale));

   % Y-param
   py = conf_grid (3, 0.015 + yo, 0.035 + yo, int(30 * scale));

   return (px, py);
}

%}}}

init_demo ();

define mask_hook (p1, p2)
{
   return (p2 < 0.025);
}

variable info = struct
{
   fail, save, mask
};
info.mask = &mask_hook;

variable calc = conf_map_counts (map(2, 0.0, 0.0), info);
() = save_conf (calc, Output_File);

variable reload = load_conf (Output_File);
variable data_confmap = load_conf ("data/confmap.fits");

define diff_arrays (a,b)
{
   return sum(sqr(a-b));
}

variable epsilon = 1.e-8;
variable diff_a = diff_arrays (calc.chisqr, data_confmap.chisqr);
if (diff_a > epsilon)
{
   failed ("confidence map mismatch:  computed map appears incorrect; err = %S", diff_a);
}

variable diff_b = diff_arrays(reload.chisqr, data_confmap.chisqr);
if (diff_b > epsilon)
{
   failed ("confidence map mismatch:  reloaded map appears incorrect; err = %S", diff_b);
}

() = remove (Output_File);
msg ("ok\n");

