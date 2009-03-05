
variable id=open_plot ("elev.ps/vcps");
resize(15);

atoms(aped);
charsize(1.2);
label ("\\u(2S+1)\\dL", "Excitation Energy (eV)", "Ne IX");
plot_elev(Ne,9);
oplot_lines (Ne, 9, where(el_ion(Ne, 9) and wl(11,15)),4);

close_plot(id);

