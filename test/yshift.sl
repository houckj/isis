% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing yshift.... ");

() = load_data ("data/acisf01318N003_pha2.fits", 9);
set_data_exposure(1,1);

variable D = get_data_counts(1);
define data_fit(l,h,p)
{
   return D.value;
}
add_slang_function ("data");
fit_fun ("data");
set_kernel (1, "yshift");

define do_shift (s, dy)
{
   return rebin (s.bin_lo-dy, s.bin_hi-dy, s.bin_lo, s.bin_hi, s.value);
}

define test_shift (dy)
{
   set_par ("yshift(1).offset", dy, 0, 0, 0);
   () = eval_counts;

   if (any (get_model_counts(1).value != do_shift (D, dy)))
     failed("shift=%g", dy);
}

test_shift (-0.2);
test_shift (0.2);
test_shift (1.0);
test_shift (10.0);
test_shift (20.0);

msg ("ok\n");
