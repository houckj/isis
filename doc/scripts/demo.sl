
path = _isis_srcdir + "/../isis-examples/data";

load_data (path + "/acisf01318N003_pha2.fits.gz");

() = load_arf (path + "/acisf01318_000N001MEG_-1_garf.fits.gz");
assign_arf(1,9);

() = load_rmf (path + "/acismeg1D1999-07-22rmfN0002.fits.gz");
assign_rmf(1,9);

plot_data_counts (9);
alias ("plot_data_counts", "pdc");

xrange(12.0,13.0);
pdc (9);

ignore ([ [1:8], [10:12] ]);
xnotice (9, 12.05, 12.2);

fit_fun ("gauss(1)");
list_par;

set_par(1,100);
set_par(2,12.15);
set_par(3,0.03);

() = renorm_counts;
() = fit_counts;

alias ("oplot_model_counts", "opmc");
opmc (9);

save_par("gfit.p");

(amin, amax) = vconf(1);
a = get_par(1);
set_par(1, a, 0, 0.9*a, 1.1*a);
(amin, amax) = vconf(1);

fp = fopen ("gfit.txt", "w");
() = fprintf (fp, "Center  Sigma   Flux (F_lower, F_upper)\n");
() = fprintf (fp, "%10.7f %12.4e  %12.4e  (%9.4e, %9.4e)   \n",
	      get_par(2), get_par(3), get_par(1), amin, amax);
fclose(fp);




