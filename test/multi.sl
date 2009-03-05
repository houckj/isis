% -*- mode: SLang; mode: fold -*- 
() = evalfile ("inc.sl");
msg ("testing multi-order.... ");

#ifndef __IMPORT__
  message (" *** isis not compiled with dynamic linking support ***");
  exit(0);
#endif

variable norm = 1.e3;
variable orders = [1:10];

% center wavelength
variable x0 = 2.0;       

variable num_orders = length(orders);

define init_data (orders)
{
   variable o, args;
   
   foreach (orders)
     {
	o = ();
	args = sprintf ("./rmf_user.so:gauss_rmf;order=%d;sigma=0.025", o);
	() = load_rmf (args);
     }

   fit_fun ("delta(1)");
   set_par (1, norm);
   set_par (2, x0);
   
   variable lo, hi, cts, id;
   (lo, hi) = linear_grid (1,24,2048);
   cts = @lo;
   cts[*] = 1.0;
   id = define_counts (lo, hi, cts, cts);
   
   assign_rsp ([1:num_orders]*0, [1:num_orders], id);
}

init_data (orders);
() = eval_counts();

variable m = get_model_counts (1);

% Check the overall total

variable err, tot;
tot = sum(m.value);
err = 1.0 - tot / (norm * num_orders);

if (abs(err) > 1.e-3)
  failed ("wrong total counts:  err = %e", err);

% Check the counts per dispersed order

define sum_roi (m, xmin, xmax)
{
   variable i = where (xmin <= m.bin_lo and m.bin_hi < xmax);
   
   return sum(m.value[i]);
}

variable o, tot_roi = Double_Type [num_orders];
foreach (orders)
{
   o = ();
   tot_roi[o-1] = sum_roi (m, x0*o - 1.0, x0*o + 1.0);
}

err = 1.0 - tot_roi/norm;
if (0 != length(where(abs(err) > 0.002)))
  failed ("wrong counts per order");

msg ("ok\n");


