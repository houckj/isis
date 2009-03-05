%-*- slang -*-
() = evalfile ("inc.sl");
msg ("testing dataset combination.... ");

variable elo, ehi, arf;
(elo, ehi, arf) = readcol ("arf.dat", [1,2,3]);
elo[0] += 1.e-3*ehi[0];

variable a = struct
{
   bin_lo, bin_hi, value, err
};
(a.bin_lo, a.bin_hi) = _A(elo, ehi);
a.value = @reverse(arf);
a.err = Double_Type[length(arf)];

() = define_arf (a);
set_arf_exposure (1, 1.46628223140590E+04);

fit_fun ("poly(1)");
set_par (1, 1, 0);
set_par (2, 1, 0);
set_par (3, 0, 1);

assign_arf (1, 1);
assign_arf (1, 2);

fakeit;

static variable c = get_data (1);
c.value *= 2.0; c.err = sqrt (c.value);
put_data (2, c);

set_fit_method ("subplex");

%message ("\n\n **** Fitting first data set ===>");
ignore (2); notice (1); ()=fit_counts();

%message ("\n\n **** Fitting second data set ===>");
ignore (1); notice (2); ()=fit_counts();

%message ("\n\n **** Performing a simultaneous fit ===>");
notice (1); notice (2); ()=fit_counts();

%message ("\n\n **** Fitting combined data sets ===>");
notice (1); notice (2);

()=combine_datasets ([1,2]);
variable s;
()=fit_counts(&s);

if (abs (1 - s.statistic/2691.301) > 1.e-3 or s.num_bins !=1024)
  failed ("ds_combine:  stat = %15.7e, num_bins = %d",
	  s.statistic, s.num_bins);

if (abs (1 - get_par(1)/1.5) > 2.e-3
    or abs (1 - get_par(2)/1.5) > 2.e-3)
{
   list_par;
   failed ("ds_combine:  ==> wrong parameter values");
}

% Now test rebin_combined:

delete_data (all_data);

variable s = struct
{
   bin_lo, bin_hi, value, err
};
variable n = 100;
variable value = 5.0;

(s.bin_lo, s.bin_hi) = linear_grid (1,20,n);
s.value = ones(n) * value;
s.err = sqrt(s.value);

() = array_map (Int_Type, &define_counts, [s,s,s,s]);

% rebin multiple combinations

variable gid1, gid2;
gid1 = combine_datasets (1,2);
gid2 = combine_datasets (3,4);
rebin_combined ([gid1, gid2], 20);
variable a1, a2;
a1 = get_combined (gid1, &get_data_counts);
a2 = get_combined (gid2, &get_data_counts);
if (any(a1.value != 20))
  failed ("ds_combine:  a1 failed");
if (any(a2.value != 20))
  failed ("ds_combine:  a2 failed");
rebin_data (all_data, 0);

% ok, we already know it works, but test anyway.

variable gid = combine_datasets (all_data);
variable c1, c2, c3;
c1 = get_combined (gid, &get_data_counts);
if (any(c1.value != 20))
  failed ("ds_combine:  c1 failed");

% test rebin_combined
rebin_combined (gid, 80);
c2 = get_combined (gid, &get_data_counts);
if (any(c2.value != 80))
  failed ("ds_combine:  c2 failed");

% test revert to original binning
rebin_combined (gid, 0);
c3 = get_combined (gid, &get_data_counts);
if (any(c3.value != 20))
  failed ("ds_combine:  c3 failed");

% test rebin_combined with rebin mask
s.value = urand(n)*25.0;
s.err = sqrt(s.value);
variable id = define_counts (s);
rebin_data (id, 20);
variable info = get_data_info (id);
delete_data (id);

rebin_combined (gid, info.rebin);
variable c4 = get_combined (gid, &get_data_counts);

if (sum(c1.value) != sum(c4.value))
  failed ("ds_combine:  c4 failed");

msg ("ok\n");
