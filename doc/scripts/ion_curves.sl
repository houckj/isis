plasma(aped);

variable Fe=26;
variable id=open_plot("ion_frac.ps/vcps");
resize(15);

variable t = 10.0^[6.0:7.2:0.05];          % define a log temperature grid (K)

label ("Temperature [K]", "Ion Fraction", "");
plot (t, ion_frac(Fe, 17, t));    % assuming Fe = 26 defined in ~/.isisrc
close_plot(id);

variable Ca = 20;
variable j_fe, f_fe;
variable j_ca, f_ca;
(j_fe, f_fe) = ion_bal (Fe, 1.e7);         % Fe at T = 1.e7 K
(j_ca, f_ca) = ion_bal (Ca, 1.e7);         % Ca at T = 1.e7 K

variable id=open_plot ("ion_bal.ps/vcps");
resize(15);
xrange(14,26);
label ("Ion Charge", "Ion Fraction", "");
hplot (j_ca-0.5, j_ca+0.5, f_ca);
ohplot (j_fe-0.5, j_fe+0.5, f_fe);
color(1);
close_plot(id);
