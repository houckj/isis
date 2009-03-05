% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing region_stats.... ");

public define contin_density (x, p) %{{{
{
   return p[0] + p[1] * (x - p[2]);
}

public define lev_fit(l,h,p)
{
   return (h-l)*contin_density (0.5 * (h + l), p);
}
add_slang_function ("lev", ["y0", "slope", "x0"]);

%}}}

public define square_fit(l,h,p) %{{{
{
   variable v = @l;
   v[*] = 0.0;
   variable mid = 0.5*(l+h);
   v[where (p[1] <= mid and mid <= p[2])] = p[0];
   return v * (h-l);
}
add_slang_function ("square", ["density", "lo_edge", "hi_edge"]);

%}}}

fit_fun ("lev(1) + square(1) + square(2)");

variable xa = 12.0 + [-0.1, 0.1];
variable xe = 14.0 + [-0.1, 0.1];
variable contin_params = [1000.0, -50, 12];

set_par ("lev(1)", contin_params);
set_par ("square(1)", [-500, xa[0], xa[1]]);
set_par ("square(2)", [ 500, xe[0], xe[1]]);

define check_values (name, got, expected) %{{{
{
   variable fail;

   % vmessage ("%s:  %g [%g]", name, got, expected);

   if (expected != 0)
     fail = (abs(got/expected-1.0) > 0.02);
   else
     fail = abs(got-expected) > 1.e-5;

   if (fail)
     failed ("got %0.4g, expected %0.4g", got, expected);
}

%}}}

define check_square_stats (which, s) %{{{
{
   variable name = sprintf ("square(%d)", which);

   variable lo, hi, dns;
   dns = get_par (name + ".density");
   lo = get_par (name + ".lo_edge");
   hi = get_par (name + ".hi_edge");

   variable sum, net, centroid, contin, slope, eqwidth;
   variable y1, y2;

   y1 = contin_density (lo, contin_params) + dns;
   y2 = contin_density (hi, contin_params) + dns;
   sum = 0.5 * (y1 + y2) * (hi - lo);
   net = dns * (hi - lo);
   centroid = 0.5 * (lo + hi);
   contin = contin_density (centroid, contin_params);
   slope = contin_params[1];
   eqwidth = abs(net) / contin;

   check_values ("net", s.net, net);
   check_values ("sum", s.sum, sum);
   check_values ("centroid", s.centroid, centroid);
   check_values ("contin", s.contin, contin);
   check_values ("slope", s.slope, slope);
   check_values ("eqwidth", s.eqwidth, eqwidth);
}

%}}}

variable lo, hi, f, binsize=0.001;
(lo, hi) = linear_grid (11, 15, (15-11.0)/binsize);

f = eval_fun (lo,hi);
() = define_flux (lo, hi, f, sqrt(f));

variable a, e, ya, ye;

ya = contin_density (xa, contin_params);
ye = contin_density (xe, contin_params);
a = region_flux (1, xa[0], xa[1], ya[0], ya[1]);
e = region_flux (1, xe[0], xe[1], ye[0], ye[1]);

check_square_stats (1, a);
check_square_stats (2, e);

#iffalse
define plot_line (x_end, b, s)
{
   variable x, y;
   x = [x_end[0]:x_end[1]:(x_end[1]-x_end[0])/100.0];
   y = b + s.slope * (x - x_end[0]);
   oplot (x, y);
}

xrange (11,15);
plot_unit ("a");
plot_data_flux;
plot_line (xa, ya[0], a);
plot_line (xe, ye[0], e);
#endif

msg ("ok\n");
