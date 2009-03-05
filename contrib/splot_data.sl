%; Time-stamp: <2002-02-27 11:48:41 dph>
%; MIT Directory: ~dph/libisis
%; File: splot_data.sl
%; Author: D. Huenemoerder
%; Original version: 2001.10.26
%;====================================================================
% version: 0.1
%
% purpose: smooth histogram via convolution w/ gaussian before plotting.
%
%  History:
%  2002.05.29    jch   minimize duplicate code, support auto-color

require("gsmooth");
provide("splot_data_counts");
provide("osplot_data_counts");
provide("splot_data_flux");
provide("osplot_data_flux");

provide("splot_model_counts");
provide("osplot_model_counts");
provide("splot_model_flux");
provide("osplot_model_flux");

provide("shplot");
provide("oshplot");

% generic get-args and plot functions
private define get_splot_args (nargs, msg)
{
   variable id, sigma, color;

   switch (nargs)
     {
      case 2:
	(id, sigma) = ();
	color = NULL;
     }
     {
      case 3:
	(id, sigma, color) = ();
     }
     {
	% default:
	usage (msg);
     }

   return (id, sigma, color);
}

private define do_shplot (plot_ref, xlo, xhi, y, sigma, color)
{
   variable s = gsmooth (xlo, xhi, y, sigma);
   
   if (color != NULL)
     @plot_ref (xlo, xhi, s, color);
   else
     @plot_ref (xlo, xhi, s);
}

private define do_splot (id, sigma, color, get_ref, plot_ref)
{
   variable v = @get_ref (id);
   do_shplot (plot_ref, v.bin_lo, v.bin_hi, v.value, sigma, color);
}

% plot/oplot for data/model and counts/flux cases
define splot_data_counts ()
{
   variable id, sigma, color;
   variable msg = "splot_data_counts (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_data_counts, &hplot);
}

define osplot_data_counts ()
{
   variable id, sigma, color;
   variable msg = "osplot_data_counts (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_data_counts, &ohplot);
}

define splot_data_flux ()
{
   variable id, sigma, color;
   variable msg = "splot_data_flux (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_data_flux, &hplot);
}

define osplot_data_flux ()
{
   variable id, sigma, color;
   variable msg = "osplot_data_flux (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_data_flux, &ohplot);
}

define splot_model_counts ()
{
   variable id, sigma, color;
   variable msg = "splot_model_counts (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_model_counts, &hplot);
}

define osplot_model_counts ()
{
   variable id, sigma, color;
   variable msg = "osplot_model_counts (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_model_counts, &ohplot);
}

define splot_model_flux ()
{
   variable id, sigma, color;
   variable msg = "splot_model_flux (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_model_flux, &hplot);
}

define osplot_model_flux ()
{
   variable id, sigma, color;
   variable msg = "osplot_model_flux (id, sigma [, color])";

   (id, sigma, color) = get_splot_args (_NARGS, msg);
   do_splot (id, sigma, color, &get_model_flux, &ohplot);
}

% generic get-args
private define get_shplot_args (nargs, msg)
{
   variable xlo, xhi, y, sigma, color;

   switch (nargs)
     {
      case 4:
	(xlo, xhi, y, sigma) = ();
	color = NULL;
     }
     {
      case 5:
	(xlo, xhi, y, sigma, color) = ();
     }
     {
	% default:
	usage (msg);
     }

   return (xlo, xhi, y, sigma, color);
}

define shplot ()
{
   variable xlo, xhi, y, sigma, color;
   variable msg = "shplot( xlo, xhi, y, sigma[, color])";

   (xlo, xhi, y, sigma, color) = get_shplot_args (_NARGS, msg);
   do_shplot (&hplot, xlo, xhi, y, sigma, color);
}

define oshplot ()
{
   variable xlo, xhi, y, sigma, color;
   variable msg = "oshplot( xlo, xhi, y, sigma[, color])";

   (xlo, xhi, y, sigma, color) = get_shplot_args (_NARGS, msg);
   do_shplot (&ohplot, xlo, xhi, y, sigma, color);
}

