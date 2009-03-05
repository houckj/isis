%-*- slang -*-
() = evalfile ("inc.sl");
msg ("testing post-model hook.... ");

Minimum_Stat_Err=1.e-10;

define create_data (n, include_background)
{
   variable lo, hi, val;
   (lo,hi) = linear_grid (1, 101, n);
   variable id = define_counts (lo, hi, hi-lo, ones(n));
   if (include_background)
     () = _define_back (id, hi-lo);
   return id;
}

define tweak_model (lo, hi, counts, bgd)
{
   return counts - bgd;
}

define do_dataset (id, answer)
{
   variable s = sum(get_model_counts(id).value);
   variable err = 0.001;
   variable msg = "post_model_hook";
   
   if (answer == 0.0)
     {
	if (abs(s) > err)
	  failed (msg);
     }
   else if (abs(1.0 - s/answer) > err)
     failed (msg);     
}

define do_check (answer)
{
   () = eval_counts;
   array_map (Void_Type, &do_dataset, all_data, answer);
}

variable Num_Bins = 100;

define try_hooks ()
{
   set_post_model_hook (all_data, NULL);
   do_check(2.0*Num_Bins);
   set_post_model_hook (all_data, &tweak_model);
   do_check(0.0);
   set_post_model_hook (all_data, NULL);
   do_check(2.0*Num_Bins);
}

%group_data(1,4);

fit_fun ("bin_width(1)");

% one dataset
() = create_data (Num_Bins, 1);
try_hooks();

% two datasets
() = create_data (Num_Bins, 1);
try_hooks();

% three datasets, 3rd with no background
() = create_data (Num_Bins, 0);
() = eval_counts;

msg("ok\n");
