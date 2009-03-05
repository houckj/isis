% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing rebin_dataset.... ");

fit_fun ("Powerlaw(1)");

assign_rsp (load_rmf ("data/acismeg1D1999-07-22rmfN0002.fits"),
            load_arf ("data/acisf01318_000N001MEG_-1_garf.fits"),
            load_data ("data/acisf01318N003_pha2.fits", 9));
variable d0 = get_data_counts(1);

() = eval_counts;
variable m0 = get_model_counts(1);

variable b = get_data_counts(1).value;
b[*] = 10.0;
() = _define_back (1, b);

() = eval_counts;
variable m1 = get_model_counts(1);

if (sum(m1.value) <= sum(m0.value))
{
   failed ("background wasn't added!");
}

variable lo, hi;
(lo,hi) = linear_grid(m0.bin_lo[0], m0.bin_hi[-1], 2000);
rebin_dataset (1, lo, hi);
variable d2 = get_data_counts(1);

variable d_ferr = 1 - sum(d2.value)/sum(d0.value);
if (abs(d_ferr) > 1.e-6)
{
   failed ("wrong count-total in rebinned dataset:  fractional error = %g",
           d_ferr);
}

() = eval_counts;
variable m2 = get_model_counts(1);

variable m_ferr = 1 - sum(m2.value)/sum(m1.value);
if (abs(m_ferr) > 1.e-6)
{
   failed ("wrong count-total in model for rebinned dataset:  fractional error = %g",
           m_ferr);
}

msg ("ok\n");
