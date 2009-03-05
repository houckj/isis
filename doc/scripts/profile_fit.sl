
variable use_lorentz = 0;

variable id = open_plot("profile_fit.ps/vcps");
resize(15);

plot_bin_density;

()=load_data(_isis_srcdir + "/doc/scripts/profile_fit.dat");

variable dx = 5.554e-3;              % bin width

if (use_lorentz)
{
   fit_fun("Lorentz(1) + Lorentz(2) + poly(1)");
   set_par(1, 300.0, 0, 50.0, 500.0);
   set_par(2, 16.0, 0, 15.95, 16.03);
   set_par(3, 0.02, 0, 0.001, 0.05);
   set_par(4, 120.0, 0, 50.0, 500.0);
   set_par(5, 16.03, 0, 15.95, 16.07);
   set_par(6, 0.01, 0, 0.001, 0.05);
   set_par(7, 10.0/dx);
   set_par(8, 0.0, 1, -1, 1);
   set_par(9, 0.0, 1, -1, 1);
}
else
{
   fit_fun("gauss(1) + gauss(2) + poly(1)");
   set_par(1, 2.7721065e+02,  0,  3.0000000e+01,    5.0000000e+02);
   set_par(2, 1.6026187e+01,  0,  1.5950000e+01,    1.6030000e+01);
   set_par(3, 1.0726321e-02,  0,  1.0000000e-03,    5.0000000e-02);
   set_par(4, 2.8011030e+02,  0,  5.0000000e+01,    5.0000000e+02);
   set_par(5, 1.6007506e+01,  0,  1.5950000e+01,    1.6070000e+01);
   set_par(6, 8.1140931e-03,  0,  1.0000000e-03,    5.0000000e-02);
   set_par(7, 10.0/dx);
   set_par(8, 0.0, 1, -1, 1);
   set_par(9, 0.0, 1, -1, 1);
}

list_par;
() = fit_counts;
list_par;

errorbars(1);
label("Wavelength (\\(2078))", "Counts/angstrom", "");
plot_data_counts(1);
oplot_model_counts(1);

variable x, y;
x = [15.95:16.1:7.5e-4];
y = get_cfun (x);
oplot(x,y);

close_plot(id);




