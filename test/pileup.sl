%-*- slang -*-
() = evalfile ("inc.sl");
msg ("testing pileup.... ");

define cutoff_plaw_fit (_lo, _hi, par)
{
   variable lo = @_lo;
   variable hi = @_hi;

   variable s = Double_Type[length(lo)];

   (lo, hi) = (Const_hc/hi, Const_hc/lo);
   variable alpha = 1 - par[1];

   variable non_zero = (hi > 1.0);
   lo[where ((lo < 1.0) and non_zero)] = 1.0;
   variable i = where (non_zero);
   lo = lo[i];
   hi = hi[i];
   if (abs(alpha) < 0.00001)
     {
	variable log_b = log (hi);
	variable log_a = log (lo);
	s[i] = par[0] * (log_b - log_a)*(1.0 + 0.5*alpha*(log_b+log_a));
	return s;
     }
   s[i] = par[0]*(hi^alpha - lo^alpha)/alpha;
   return s;
}

add_slang_function ("cutoff_plaw", ["norm", "alpha"]);

define init_data ()
{
   % simulated data assumes delta-function RMF
   () = load_data ("data/pileup.dat");
   set_frame_time (1, 3.2);

   variable elo, ehi, arf;
   (elo, ehi, arf) = readcol ("data/arf.dat", [1,2,3]);
   elo[0] += 1.e-3*ehi[0];

   variable a = struct
     {
        bin_lo, bin_hi, value, err
     };
   (a.bin_lo, a.bin_hi) = _A(elo, ehi);
   a.value = @reverse(arf);
   a.err = Double_Type[length(arf)];

   () = define_arf (a);
   set_arf_exposure (1, 25600000);
   assign_arf (1, 1);
}

init_data();

fit_fun ("cutoff_plaw(1)");
set_par ("cutoff_plaw(1)", [0.005, 1],  [0,0], [0, 0], [0.01, 3]);
set_kernel (1, "pileup");

set_par ("pileup(1)", [1, 1, 0.5, 1], [1,1,0,1]);

static variable s;
() = eval_counts (&s);
if (s.statistic/s.num_bins > 1.0)
  failed ("pileup:  stat/num_bins = %g", s.statistic/s.num_bins);

msg ("ok\n");
