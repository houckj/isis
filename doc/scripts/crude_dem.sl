
plasma(aped);

variable id=open_plot("crude_dem.ps/vcps");
resize(15);

variable k, t, em, i, flux_obs;

% Fe strong XVII line at 15.014 Angstrom
k = where(el_ion(Fe,17) and wl(15.0,15.02));   
flux_obs = 1.0;

t = 10.0^[5.8:7.5:0.01];             % define a Temperature grid
em = line_em (k,t);                  % get emissivity vs. temperature for
                                     % line k
i = where(em > 0.0);
limits;
variable y = flux_obs / em[i];
yrange (0.2*min(y), 10*min(y));
label ("Temperature (K)", "EM (arbitrary units)", "Fe XVII 15.014 Angstrom");
plot (t[i], y);      % plot the ratio
close_plot(id);
