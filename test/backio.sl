% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing backio.... ");

fit_fun ("Powerlaw(1)");
set_par ("Powerlaw(1).norm", 1, 0, 0, 1.e5);

variable lo, hi, n = 1000;
(lo,hi) = linear_grid (1,20,n);
variable f = eval_fun(lo,hi);

% Create a dataset
() = define_counts (lo, hi, f, 0.1*f; min_stat_err = 1.e-10);
variable data_backscale = 2.0, data_exposure = 1.e4;
set_data_backscale (1, data_backscale);
set_data_exposure (1, data_exposure);
set_fake (1, 1);

% Define the background using a parametrized model
variable b_fun_name = "blackbody";
back_fun (1, b_fun_name);
set_par ("${b_fun_name}(1).norm"$, 1.e3, 0, 0, 0);
variable b_model = eval_fun2 (b_fun_name, lo, hi, get_par ("${b_fun_name}(1).*"$));

% Also define a background spectrum.  When back_backscale > 0,
% the uncertainties in this "data" background spectrum will be
% propagated to increase the uncertainty in each data bin.
variable b_data = eval_fun2 (b_fun_name, lo, hi, [1.e3, 2.0]);
variable back_backscale = PI, back_exposure = 1.e5;
() = _define_back (1, b_data, back_backscale, back_exposure);

% Compute the expected value of the background scale factor
variable b_scale = ((data_backscale * data_exposure)
                    / (back_backscale * back_exposure));

% generate fake data with Poisson errors and evaluate the model
fakeit();
() = eval_counts;

% Retrieve and check the various background components
variable b_back_fun_got = get_back_fun (1);
ifnot (0 == strcmp (b_back_fun_got, b_fun_name))
  failed ("get_back_fun failed");

variable b_data_got = get_back_data (1);
if (any(b_data_got != b_data))
  failed ("get_back_data failed");

variable b_model_got = get_back_model (1);
if (any(b_model_got != b_model))
  failed ("get_back_model failed");

variable b_scale_got =  get_back_data_scale_factor (1);
if (any(b_scale_got != b_scale))
  failed ("get_back_data_scale_factor failed");

variable b_back_got = get_back (1);
if (any (b_back_got != b_model))
  failed ("get_back failed");

variable b_backscale_got = get_back_backscale (1);
if (b_backscale_got != back_backscale)
  failed ("get_back_backscale failed");

variable b_exposure_got = get_back_exposure (1);
if (b_exposure_got != back_exposure)
  failed ("get_back_exposure failed");

msg ("ok\n");

