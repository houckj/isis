% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing fit.... ");

%{{{ Test slang function definitions.

define testfun1_fit(l,h,p){return urand(length(l));}
add_slang_function ("testfun1", ["a", "b", "norm"]);
fit_fun ("testfun1(1)");
private variable _t = get_params();
!if (andelse
    {_t[0].is_a_norm == 0}
    {_t[1].is_a_norm == 0}
    {_t[2].is_a_norm == 1})
{
   failed ("[1]norm indexes don't match");
}
define testfun2_fit(l,h,p){return urand(length(l));}
add_slang_function ("testfun2", ["norm", "a", "b"], [0,2]);
fit_fun ("testfun2(1)");
private variable _t = get_params();
!if (andelse
    {_t[0].is_a_norm == 1}
    {_t[1].is_a_norm == 0}
    {_t[2].is_a_norm == 1})
{
   failed ("[2]norm indexes don't match");
}

%}}}

define check_num_pars (num_expected, id_string)
{
   if (num_expected != get_num_pars ())
     failed ("wrong number of params [%s]! => %d (expected %d)",
	     id_string,
	     get_num_pars(), num_expected);
}

define match_pars (pars) %{{{
{
   variable i, num_pars = get_num_pars ();
   for (i = 0; i < num_pars; i++)
     {
	if (get_par(i+1) != pars[i].value)
	  return -1;
     }

   return 0;
}

%}}}

% Test basic operations
define test_basic_ops (freeze_par, tie_par, test_name) %{{{
{
   variable old_pars, pars, num_pars;

   % copy/restore
   old_pars = get_params();
   if (old_pars == NULL) failed ("%s: %s", test_name, "get_params");
   randomize;
   pars = get_params();
   num_pars = get_num_pars();
   if (num_pars != length(pars)) failed ("%s: %s", test_name, "get_num_pars");
   if (0 != match_pars (pars)) failed ("%s: %s", test_name, "get_par");
   set_params (old_pars);
   if (0 != match_pars (old_pars)) failed ("%s: %s", test_name, "set_params");

   % freeze, thaw, tie, untie
   freeze(freeze_par);
   if (get_par_info(freeze_par).freeze != 1) failed ("%s: %s", test_name, "freeze");
   thaw(freeze_par);
   if (get_par_info(1).freeze != 0) failed ("%s: %s", test_name, "thaw");
   tie(tie_par,freeze_par);
   if (_get_index(get_par_info(freeze_par).tie) != _get_index(tie_par)) failed ("%s: %s", test_name, "tie");
   untie(freeze_par);
   if (_get_index(get_par_info(freeze_par).tie) != NULL) failed ("%s: %s", test_name, "untie");

   variable save_tie_par, fp, fun;

   % Define params as functions of other params
   %  tie_par <= f(freeze_par);
   save_tie_par = get_par (tie_par);
   fp = get_par (freeze_par);
   fun = sprintf ("_par(%d)*1.1", _get_index(freeze_par));
   set_par_fun (tie_par, fun);
   if (0 != strcmp (get_par_info (tie_par).fun, fun)
       or get_par(tie_par) != fp*1.1)
     failed ("%s: %s", test_name, "set_par_fun[1]");

   set_par_fun (tie_par, NULL);
   if (get_par_info (tie_par).fun != NULL)
     failed ("%s: %s", test_name, "set_par_fun[2]");
   set_par (tie_par, save_tie_par);

   % save and restore
   variable par_file = "/tmp/.isis." + string(getpid());
   () = remove (par_file);
   save_par (par_file);
   old_pars = get_params();
   fit_fun ("delta(1)");
   load_par (par_file);
   if (0 != match_pars (old_pars)) failed ("%s: %s", test_name, "save_par/load_par");
   () = remove (par_file);
}

%}}}

define add_fake_dataset () %{{{
{
   variable x = grand(100000);
   x += 2*abs(min(x));

   variable n, lo, hi;
   n = 1024;
   (lo, hi) = linear_grid (min(x), max(x), n);

   variable nx = histogram (x, lo, hi);
   nx += 10;
   variable id = define_counts (lo, hi, nx, sqrt(nx));
   if (id < 0) failed ("init_fake_data");

   set_frame_time (id, 3.2);
   variable info = get_data_info (id);
   info.order = 1;
   set_data_info (id, info);

   return id;
}

%}}}

define try_user_statistic () %{{{
{
   fit_fun ("Lorentz(1)");
   set_par (1,15000, 0, 1500, 1.e9);
   set_par (2,12);
   set_par (3,0.05);

   variable lo, hi, val;
   (lo,hi) = linear_grid (11.6,12.4,32);
   val = eval_fun (lo,hi);
   () = define_counts (lo,hi,val, prand(val));

   set_fit_statistic ("chisqr");
   variable chi_info;
   if (-1 == eval_counts(&chi_info))
     failed ("max_like [chisqr]");

   () = evalfile ("max_like.sl");
   variable max_info;
   if (-1 == eval_counts(&max_info))
     failed ("max_like [1]");

   if (max_info.statistic == chi_info.statistic)
     failed ("max_like [2]");
}

%}}}

% a simple fit function
fit_fun ("gauss(1)");
check_num_pars (3, "simple");
test_basic_ops (2, 3, "simple");
test_basic_ops ("gauss(1).center", "gauss(1).sigma", "simple");

% odd corner cases
fit_fun ("poly(1) + poly(1)");
check_num_pars (3, "corner");

% kernel parameters
variable id = add_fake_dataset ();
set_kernel (id, "pileup");
check_num_pars (7, "pileup");
test_basic_ops (1, 6, "pileup");

% A second dataset and kernel
id = add_fake_dataset ();
set_kernel (id, "pileup");
check_num_pars (11, "pileup");
test_basic_ops (2, 3, "two kernels");

% Add an instrumental background
back_fun (1, "Powerlaw(1)");
check_num_pars (13, "back[1]");
test_basic_ops (3, 9, "back");
back_fun (1, NULL);
check_num_pars (11, "back[2]");

% Now let each dataset get its own function....
define schizo () %{{{
{
   switch (Isis_Active_Dataset)
     {
      case 1:
	return poly(1) + blackbody(1);
     }
     {
      case 2:
	return Powerlaw(2);
     }
     {
	% default:
	message ("*** [schizo] too many datasets");
	exit(1);
     }
}

%}}}
fit_fun ("schizo");
back_fun (2, "poly(2)");
% poly(1) + blackbody(1) + powerlaw(2) => 7
%                      pileup + pileup => 5
%                              poly(2) => 3
%                             --------------
%                                        15 params
check_num_pars (18, "function per dataset");
test_basic_ops (14, 3, "function per dataset");
back_fun (2, NULL);
check_num_pars (15, "function per dataset");

% first dataset => standard kernel
% and lose pileup parameters
set_kernel (1, "std");
check_num_pars (11, "revert to std");
test_basic_ops (1, 4, "revert to std");

% Now delete a dataset
%  poly(1) + blackbody(1) => lose 5
delete_data (1);
check_num_pars (6, "delete dataset");
test_basic_ops (1, 2, "delete dataset");

% user-defined functions
() = evalfile("bpl.sl");
fit_fun ("bpl(1)");
bpl_set_default (0);
check_num_pars (8, "S-Lang fit-function");
test_basic_ops (3, 4, "S-Lang fit-function");

% user-defined statistic
delete_data (all_data);
check_num_pars (4, "UDF");
try_user_statistic ();

msg ("ok\n");
