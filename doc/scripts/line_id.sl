
variable dir = _isis_srcdir + "/../isis-examples/data";

% load data file containing more than 9 spectra
() = load_data (dir + "/acisf01318N003_pha2.fits.gz");

plasma(aped);
load_model(_isis_srcdir + "/doc/scripts/model.dat");

variable d, flux, g;

d = get_data_counts (9);
flux = model_spectrum (d.bin_lo, d.bin_hi);
g = brightest(10, where (wl(10,12)));    % get indices of bright lines

variable id=open_plot ("line_id.ps/vcps");
resize(15);

xrange(10,12);
yrange(0,50);
plot_data_counts (9);


% and modify the default line label style:
variable x = line_label_default_style();
x.angle = 45.0;    % angle the labels for easier reading
x.label_type = 1;  % use long labels

plot_group(g, 1, x);

close_plot(id);
