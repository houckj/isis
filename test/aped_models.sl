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

define check_value (test_name, expected_value, value)
{
   if ((typeof(value) == Int_Type && expected_value != value)
       || (abs (expected_value/value - 1.0) > 1.e-5))
     {
        failed ("$test_name:\n    got $value, expected $expected_value"$);
     }
}

variable lo, hi, f;
(lo, hi) = linear_grid (1,20,2000);

% test line profile modifier
fit_fun ("xaped (1, square_profile(1))");
set_par ("square_profile(1).width", 5.e-3);
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             1.3685056789748902, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.011614329216854643, max(f));

% line profile modifier should now have no effect
fit_fun ("xaped (1)");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             1.3685056789748902, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.07392995170953862, max(f));

% line profile modifier should work now
fit_fun ("xaped (1, square_profile(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             1.3685056789748902, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.011614329216854643, max(f));

% test emissivity modifier
fit_fun ("xaped (1, modifier(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             0.8140366080700342, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.07392995170953862, max(f));

% emissivity modifier should now be disabled
fit_fun ("xaped (1)");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             1.3685056789748902, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.07392995170953862, max(f));

% emissivity modifier should work again
fit_fun ("xaped (1, modifier(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             0.8140366080700342, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.07392995170953862, max(f));

% Test emissivity and line modifier together
fit_fun ("xaped (1, modifier(1), square_profile(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             0.8140353923673627, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.0076028565585263075, max(f));

% Test emissivity and line modifier together, order should not matter
fit_fun ("xaped (1, square_profile(1), modifier(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             0.8140353923673627, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.0076028565585263075, max(f));

% all modifiers should now be disabled
fit_fun ("xaped (1)");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             1.3685056789748902, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.07392995170953862, max(f));

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
             0.030758977978675755, sum(f));
check_value (sprintf ("max fit_fun = %s with hook", get_fit_fun()),
             0.016703498004343294, max(f));
check_value (sprintf ("number of non-zero bins with fit_fun = %s with hook", get_fit_fun()),
             29, length(where(f>0)));

% The hook and modifiers should work simultaneously
fit_fun ("xaped (1, modifier(1), square_profile(1))");
f = eval_fun (lo, hi);
check_value (sprintf ("sum fit_fun = %s with hook", get_fit_fun()),
             0.005127893495442973, sum(f));
check_value (sprintf ("max fit_fun = %s with hook", get_fit_fun()),
             0.0006447144328433119, max(f));
check_value (sprintf ("number of non-zero bins with fit_fun = %s with hook", get_fit_fun()),
             50,length(where(f>0)));

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
             512.7893495443059, sum(m.value));
check_value (sprintf ("max fit_fun = %s with hook", get_fit_fun()),
             32.23572164216559, max(m.value));
check_value (sprintf ("number of non-zero bins with fit_fun = %s with hook", get_fit_fun()),
             80,length(where(m.value>0)));

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
             1.3685068946775618, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.07392995170953862, max(f));

% eval_fun2 with xaped + line emissivity modifier
f = eval_fun2 ("xaped", lo, hi, p,
               aped_line_modifier_args ("modifier", mp));
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             0.8140366080700342, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.07392995170953862, max(f));

% eval_fun2 with xaped + line profile modifier
f = eval_fun2 ("xaped", lo, hi, p,
               aped_line_profile_args ("square_profile", sp));
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             1.368499481870153, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.0069056109324198926, max(f));

% eval_fun2 with xaped + two modifiers
f = eval_fun2 ("xaped", lo, hi, p,
               aped_line_modifier_args ("modifier", mp),
               aped_line_profile_args ("square_profile", sp));
check_value (sprintf ("sum fit_fun = %s", get_fit_fun()),
             0.8140291952626265, sum(f));
check_value (sprintf ("max fit_fun = %s", get_fit_fun()),
             0.0038878650885680738, max(f));

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
             0.2772986393051501, brightest_line_flux);
check_value ("brightest line flux from internal sum",
             0.2772986393051501, line_info(id).flux);
check_value ("brightest line flux from O VIII Ly alpha lookup",
             0.2772986393051501, line_info(o8lya).flux);

msg ("ok\n");
