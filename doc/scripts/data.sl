
variable dir = _isis_srcdir + "/../isis-examples/data";

% load data file containing more than 9 spectra

() = load_data (dir + "/acisf01318N003_pha2.fits.gz");
() = load_arf (dir + "/acisf01318_000N001MEG_-1_garf.fits.gz");
assign_arf(1,9);
flux_corr(9,2);   % flux-correct down to 2-sigma significance level

variable id=open_plot ("data_cts.ps/vcps");
resize(15);
plot_data_counts (9);          % plot spectrum #9
close_plot (id);

id=open_plot ("data_flx.ps/vcps");
resize(15);
xrange(8,25);
plot_data_flux (9);          % plot spectrum #9
close_plot (id);



