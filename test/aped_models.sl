% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");

variable atomdb = getenv ("ATOMDB");
if (NULL == atomdb
    || NULL == stat_file (atomdb))
{
   msg ("skipping aped models test -- no atomdb\n");
   exit (0);
}
msg ("testing aped models.... ");

% autoload aped.sl to define get_version_string function
() = aped();
variable version = get_version_string (atomdb);
variable supported_versions = ["1.3.1", "2.0.1"];
variable result_index = wherefirst (version == supported_versions);
if (result_index == NULL)
{
   msg ("skipping aped models test -- (test results for atomdb-$version are not available)"$);
   exit (0);
}

plasma(aped);
create_aped_fun ("xaped", default_plasma_state);

% load a line profile modifier function
variable name = "square";
variable lp = load_line_profile_function ("./example-profile.so", name);
create_aped_line_profile ("square_profile", lp, ["width"]);

% define a line emissivity modifier function
define line_emis_modifier (params, line_id, state, emis)
{
#iffalse
   print(params); print(line_id); print(state); print(emis);
   plot_pause;
#endif
   variable info = line_info(line_id);
   if (11.5 < info.lambda < 14) emis = 0.0;
   return emis;
}
create_aped_line_modifier ("modifier", &line_emis_modifier, ["a", "b", "c"]);

define check_value (test_name, value, expected_values)
{
   variable expected = expected_values[result_index];

   if ((typeof(value) == Int_Type && expected != value)
       || (abs (expected/value - 1.0) > 1.e-5))
     {
        failed ("$test_name:\n    got $value, expected $expected"$);
     }
}

variable lo, hi, f;
(lo, hi) = linear_grid (1,20,2000);

% test line profile modifier
fit_fun ("xaped (1, square_profile(1))");
set_par ("square_profile(1).width", 5.e-3);
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [1.3685056789748902,
              1.8739052675912387]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.011614329216854643,
              0.011900481798605351]);

% line profile modifier should now have no effect
fit_fun ("xaped (1)");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [1.3685056789748902,
              1.8739058698230815]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.07392995170953862,
              0.0713156235206039]);

% line profile modifier should work now
fit_fun ("xaped (1, square_profile(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [1.3685056789748902,
              1.8739052675912387]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.011614329216854643,
              0.011900481798605351]);

% test emissivity modifier
fit_fun ("xaped (1, modifier(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [0.8140366080700342,
              1.2492451849918764]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.07392995170953862,
              0.0713156235206039]);

% emissivity modifier should now be disabled
fit_fun ("xaped (1)");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [1.3685056789748902,
              1.8739058698230815]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.07392995170953862,
              0.0713156235206039]);

% emissivity modifier should work again
fit_fun ("xaped (1, modifier(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [0.8140366080700342,
              1.2492451849918764]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.07392995170953862,
              0.0713156235206039]);

% Test emissivity and line modifier together
fit_fun ("xaped (1, modifier(1), square_profile(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [0.8140353923673627,
              1.2492445827600334]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.0076028565585263075,
              0.011900481798605351]);

% Test emissivity and line modifier together, order should not matter
fit_fun ("xaped (1, square_profile(1), modifier(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [0.8140353923673627,
              1.2492445827600334]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.0076028565585263075,
              0.011900481798605351]);

% all modifiers should now be disabled
fit_fun ("xaped (1)");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [1.3685056789748902,
              1.8739058698230815]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.07392995170953862,
              0.0713156235206039]);

% Try introducing a contrib hook
variable Ne_Lines = where (el_ion(Ne,10));
variable Mg_Lines = where (el_ion(Mg,12));
variable Contrib_Flag = MODEL_LINES;

define xaped_hook (id)
{
   variable f = struct
     {
        contrib_flag, line_list
     };
   f.contrib_flag = Contrib_Flag;

   if (id == 1) f.line_list = Ne_Lines;
   else if (id == 2) f.line_list = Mg_Lines;
   else f = NULL;

   return f;
}

create_aped_fun ("xaped", default_plasma_state(), &xaped_hook);

% computed spectrum should now contain only the selected lines
fit_fun ("xaped (1)");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s with hook", get_fit_fun()),
             sum(f),
             [0.030758977978675755,
              0.0317477357507457]);
check_value (sprintf ("max fit_fun = %s with hook", get_fit_fun()),
             max(f),
             [0.016703498004343294,
              0.016849658261591658]);
check_value (sprintf ("number of non-zero bins with fit_fun = %s with hook", get_fit_fun()),
             length(where(f>0)),
             [29,
              34]);

% The hook and modifiers should work simultaneously
fit_fun ("xaped (1, modifier(1), square_profile(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s with hook", get_fit_fun()),
             sum(f),
             [0.005127893495442973,
              0.005851670512314665]);
check_value (sprintf ("max fit_fun = %s with hook", get_fit_fun()),
             max(f),
             [0.0006447144328433119,
              0.0006298649711275214]);
check_value (sprintf ("number of non-zero bins with fit_fun = %s with hook", get_fit_fun()),
             length(where(f>0)),
             [50,
              73]);

% All this machinery should work with the fit engine.
f = eval_fun (lo, hi);
() = define_counts (lo, hi, f, f);
set_fake(1, 1);
set_data_exposure(1, 1.e5);
fakeit;

set_par ("square_profile(1).width", 0.01);
thaw ("xaped(1).temperature");
thaw ("xaped(1).norm");
() = eval_counts;

variable m = get_model_counts(1);
check_value (sprintf ("sum fit_fun = %s with hook", get_fit_fun()),
             sum(m.value),
             [512.7893495443059,
              585.1670512314751]);
check_value (sprintf ("max fit_fun = %s with hook", get_fit_fun()),
             max(m.value),
             [32.23572164216559,
              31.493248556376074]);
check_value (sprintf ("number of non-zero bins with fit_fun = %s with hook", get_fit_fun()),
             length(where(m.value>0)),
             [80,
              108]);

delete_data (all_data);

% Test the eval_fun2 interface
create_aped_fun ("xaped", default_plasma_state);
variable p, mp, sp;
fit_fun ("xaped (1, square_profile(1), modifier(1))");
p = get_par(get_fun_params ("xaped(1)"));
mp = get_par(get_fun_params ("modifier(1)"));
sp = get_par(get_fun_params ("square_profile(1)"));
fit_fun ("bin_width");

% eval_fun2 with xaped alone
f = eval_fun2 ("xaped", lo, hi, p);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [1.3685068946775618,
              1.8739058698230815]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.07392995170953862,
              0.0713156235206039]);

% eval_fun2 with xaped + line emissivity modifier
f = eval_fun2 ("xaped", lo, hi, p,
               aped_line_modifier_args ("modifier", mp));
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [0.8140366080700342,
              1.2492451849918764]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.07392995170953862,
              0.0713156235206039]);

% eval_fun2 with xaped + line profile modifier
f = eval_fun2 ("xaped", lo, hi, p,
               aped_line_profile_args ("square_profile", sp));
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [1.368499481870153,
              1.8739000178550174]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.0069056109324198926,
              0.007765170850644935]);

% eval_fun2 with xaped + two modifiers
f = eval_fun2 ("xaped", lo, hi, p,
               aped_line_modifier_args ("modifier", mp),
               aped_line_profile_args ("square_profile", sp));
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             sum(f),
             [0.8140291952626265,
              1.2492393330238132]);
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             max(f),
             [0.0038878650885680738,
              0.006468087081832207]);

% Test aped_fun_details interface

% construct a model with two multi-temperature components
variable ps = default_plasma_state;
variable num_temps = 5;
variable temps = [1.0:10.0:#num_temps] * 1.e6;
ps.norm = ones(num_temps);
ps.temperature = temps;
create_aped_fun ("xaped", ps);
fit_fun ("xaped(1) + xaped(2)");
set_par ("xaped(2).temperature*", temps*2);

% evaluate the model
f = eval_fun (lo, hi);

% retrieve lists of line fluxes:
variable info_1, info_2;
info_1 = aped_fun_details ("xaped(1)");
info_2 = aped_fun_details ("xaped(2)");

% Find the contribution each temperature component
% makes to the brightest emission line, which just happens
% to be O VIII Ly alpha.

variable i, total = info_2[0].line_flux;
_for i (1, num_temps-1, 1)
{
   total += info_2[i].line_flux;
}
variable brightest_line_flux = max(total);
variable id = where (total == brightest_line_flux)[0];
variable o8lya = where (trans(O,8,4,1))[0];

check_value ("brightest line flux from aped_fun_details",
             brightest_line_flux,
             [0.2772986393051501,
              0.27665615]);
check_value ("brightest line flux from internal sum",
             line_info(id).flux,
             [0.2772986393051501,
              0.27665615025627893]);
check_value ("brightest line flux from O VIII Ly alpha lookup",
             line_info(o8lya).flux,
             [0.2772986393051501,
              0.27665615025627893]);

private define xaped_hook (id)
{
   return struct {contrib_flag = MODEL_CONTIN, line_list};
}
private define check_trace_element_contin_abundance_scaling (temp, dens, elements, abund)
{
   variable lo, hi, n = 8192;
   (lo,hi) = linear_grid (1,20,n);

   variable true_sum = Double_Type[n],
     pseudo_sum = Double_Type[n];

   variable i, ne = length(elements);
   _for i (0, ne-1, 1)
     {
        variable Z = elements[i];
        variable p = get_contin (lo, hi, temp, dens, Z);
        true_sum += p.true * abund[i];
        pseudo_sum += p.pseudo * abund[i];
     }

   variable c_sum = true_sum + pseudo_sum;

   variable s = default_plasma_state ();
   s.temperature = temp;
   s.density = dens;
   s.elem = elements,
   s.elem_abund = abund,
   create_aped_fun ("xaped", s, &xaped_hook);
   fit_fun ("xaped");
   variable c = 1.e-14 * eval_fun (lo, hi);

   if (any (c < 0 or c_sum < 0))
     throw ApplicationError, "negative continuum values!!??";

   i = where (c != 0);
   variable max_abs_rel_diff = max(abs(1.0 - c_sum[i] / c[i]));

#iffalse
   ylog;
   hplot (lo,hi,c);
   ohplot (lo, hi, c_sum);
   plot_pause;
#endif

   return max_abs_rel_diff;
}

private variable elements, abund, max_abs_rel_diff, max_Z,
  temp, dens = 1.0;

% low T test:
temp = 1.e7;

if (result_index == 1)
{
   max_Z = 28;
   elements = [1:max_Z];
   abund = Double_Type[max_Z];
   abund[ H-1] = 1.0;
   abund[He-1] = 1.0;
   abund[ C-1] = 10;
   abund[ O-1] = 10;
}
else if (result_index == 0)
{
   elements = [ H, He,  C,N, O,Ne,Mg,Al,Si,S,Ar,Ca,Fe,Ni];
   abund =    [1., 1., 10,0,10, 0, 0, 0, 0,0, 0, 0, 0, 0];
}

max_abs_rel_diff = check_trace_element_contin_abundance_scaling
                                     (temp, dens, elements, abund);
%vmessage ("max_abs_rel_diff = $max_abs_rel_diff"$);
if (max_abs_rel_diff > 1.e-6)
throw ApplicationError,
   "temp = $temp:  max_abs_rel_diff = $max_abs_rel_diff, expected <~5e-08"$;

% high T test:
temp = 7.e7;

if (result_index == 1)
{
   max_Z = 28;
   elements = [1:max_Z];
   abund = Double_Type[max_Z];
   abund[ H-1] = 1.0;
   abund[He-1] = 1.0;
   abund[ S-1] = 100;
   abund[Fe-1] = 100;
}
else if (result_index == 0)
{
   elements = [ H, He, C,N,O,Ne,Mg,Al,Si, S,Ar,Ca,Fe,Ni];
   abund =    [1., 1., 0,0,0, 0, 0, 0, 0,10, 0, 0,10, 0];
}

max_abs_rel_diff = check_trace_element_contin_abundance_scaling
                                     (temp, dens, elements, abund);
%vmessage ("max_abs_rel_diff = $max_abs_rel_diff"$);
if (max_abs_rel_diff > 1.e-6)
throw ApplicationError,
   "temp = $temp,  max_abs_rel_diff = $max_abs_rel_diff, expected <~5e-08"$;

msg ("ok\n");
