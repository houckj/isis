% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing backscale.... ");

private variable Exposure = 1.e4;
private variable Test_Number = 0;

define fake_dataset (n)
{
   variable d = struct
     {
        bin_lo, bin_hi, value, err
     };

   (d.bin_lo, d.bin_hi) = linear_grid(1,20,n);
   d.value = prand (30, n);
   d.err = sqrt(d.value);

   variable i = define_counts (d);
   set_data_exposure (i, Exposure);

   return i;
}

define try_getset_data_backscale (i, a, expect_error)
{
   Test_Number++;
   ERROR_BLOCK
     {
        _clear_error();
        !if (expect_error)
          failed ("Test %d: unexpected error in set_data_backscale", Test_Number);
     }
   set_data_backscale (i, a);

   variable aa = get_data_backscale (i);

   variable err = 0;
   if (length(a) == length(aa))
     {
        if (any(a != aa) and expect_error == 0)
          err = 1;
     }
   else err = not expect_error;

   if (err)
     failed ("Test %d: get_data_backscale disagrees with set_data_backscale", Test_Number);
}

define try_define_back (i, area, expect_error)
{
   Test_Number++;
   ERROR_BLOCK
     {
        _clear_error();
        !if (expect_error)
          failed ("Test %d: unexpected error in _define_back", Test_Number);
     }

   variable n = length(get_data_counts(i).value);
   variable vec_back = 100*urand(n);

   variable status;
   status = _define_back (i, vec_back, area, Exposure);
   if (status != 0 and expect_error == 0)
     failed ("Test %d: _define_back failed", Test_Number);
}

variable n = 1024;
variable i = fake_dataset (n);

% turn off error messages.
Isis_Verbose = -2;

try_getset_data_backscale (i, 1.0, 0);
try_getset_data_backscale (i, 0.0, 0);
try_getset_data_backscale (i, -1.0, 1);
try_getset_data_backscale (i, urand(n), 0);
try_getset_data_backscale (i, urand(n)-0.5, 0);
try_getset_data_backscale (i, urand(n*2), 1);
try_getset_data_backscale (i, urand(n/2), 1);

try_define_back (i, 10.0, 0);
try_define_back (i, 0.0, 0);
try_define_back (i, -1.0, 1);
try_define_back (i, 10*urand(n), 0);
try_define_back (i, 10*urand(n) - 5, 1);
try_define_back (i, 10*urand(n*2), 1);
try_define_back (i, 10*urand(n/2), 1);

delete_data (all_data);

variable Expected_Value = 1000.0;

% Create a fake dataset.
variable n = 2000;
variable lo, hi, d;
(lo, hi) = linear_grid (1,20,n);
d = Double_Type[n];
d[*] = Expected_Value;
variable i = define_counts (lo, hi, d, sqrt(d));

% Assign a background
variable b = Double_Type[n];
b[*] = 200.0;
() = _define_back (i, b);

% Assign vector BACKSCAL values for data and background
variable da = Double_Type[n];
da[*] = 1.0;
set_data_backscale (i, da);

variable ba = Double_Type[n];
ba[*] = 10.0;
set_back_backscale (i, ba);

% Construct the model so that model+background
% fits the data exactly  (e.g. model = data - scaled_background)
variable model = d - b * da/ba;

define a_fit(l,h,p)
{
   return p[0] * ones(length(l));
}
add_slang_function ("a", "val");
fit_fun ("a");
set_par ("a(1).val", model[0]);

define tryit()
{
   Test_Number += 1;

   variable info;
   () = eval_counts (&info);

   % Compute chi-square.
   % By construction, chi-square should be zero.
   if (info.statistic != 0)
     {
        failed ("Test %d:  incorrect statistic", Test_Number);
     }

}

tryit();
   % By construction, rebinning should not affect
   % the fit of the model.

group_data(1, 5);
tryit();

msg ("ok\n");
