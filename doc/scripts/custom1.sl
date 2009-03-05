
variable pha2_file;
pha2_file = _isis_srcdir + "/../isis-examples/data/acisf01318N003_pha2.fits.gz";

variable id = 10;
() = load_data (pha2_file);
variable pid = open_plot ("custom1.ps/vcps");
resize(15);

xrange(10,20);
plot_data_counts (id);

% coordinates of region for inset plot:

variable xmin, xmax, ymin, ymax;
xmin = 12.35;
xmax = xmin + 0.2;
ymin = 0.0;
ymax = 30.0;

_pgsci (4);                       % Draw blue rectangle:
_pgsfs (2);
_pgrect (xmin, xmax, ymin, ymax);
_pgsfs (1);                       % restore default fill-area style

_pgsci (4);                       % Draw blue arrow:
_pgsch (0.7);
_pgsah (1, 45.0, 0.3);            % specify arrow-head shape
_pgarro (xmax, ymax, 15.9, 122);  % from (x1,y1) to (x2,y2)
_pgsci (1);  
                                  % New viewport in normalized
                                  % device coordinates:
_pgsvp (0.6, 0.8, 0.5, 0.8);      % xleft, xright, ybot, ytop
_pgswin (xmin, xmax, ymin, ymax);

                                  % Insert smaller plot
_pgsch (0.7);
_pgbox ("BCNST", 0.0, 1, "BCNST", 0.0, 1);
_pgsch (1.0);

variable d;
d = get_data(id);
_pgbin (d.bin_lo, d.value, 0);    % 0 means x value is bin left edge

plasma(aped);

charsize(0.7);
variable Fe = 26;
plot_group(where(wl(xmin, xmax) and el_ion(Fe,17)),4);
charsize(1.0);

close_plot(pid);
