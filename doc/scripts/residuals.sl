
public define plot_count_residuals (h, dev)
{
   variable id = open_plot (dev);
   resize(15);
   
   multiplot ([3,1]);

   errorbars (1);
   plot_data_counts (h);
   oplot_model_counts (h);
   
   variable d, m;
   d = get_data_counts (h);
   m = get_model_counts (h);

   errorbars (0);
   label ("Wavelength [Angstrom]", "Residuals", "");
   hplot (d.bin_lo, d.bin_hi, m.value - d.value);
   
   return id;
}

() = evalfile (_isis_srcdir + "/doc/scripts/profile_fit.sl");

message ("plotting residuals");
variable id = plot_count_residuals (1, "residuals.ps/vcps");

close_plot (id);

