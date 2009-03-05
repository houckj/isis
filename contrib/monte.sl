
%  Purpose:
%       Carry out a Monte-Carlo search of the fit-function 
%       parameter space to find a good set of starting parameters.
%       
%  Notes:
%       A Monte-Carlo search may be helpful when examining
%       a complex chi-square space such as that associated
%       with the pileup model.
%       
%       Each variable fit-parameter should be constrained
%       to lie within some reasonable min/max range.
%       If the allowed range is unconstrained, the Monte-Carlo
%       search will sample parameter values over the entire 
%       range of representable double-precision numbers.
%       
%  Usage:  
%          
%  e.g. to carry out 100 random trials with 'eval_counts':
%     monte_carlo (100, &eval_counts);
%

_debug_info = 1;
_traceback=1;

private define get_params ()
{
   variable i, n = get_num_pars ();   
   variable p = Struct_Type [n];
   
   for (i = 1; i <= n; i++)
     {
	p[i-1] = get_par_info (i);
     }
   
   return p;
}

private define set_params (p)
{
   variable i, n = get_num_pars ();
   
   for (i = 0; i < n; i++)
     {
	set_par (i+1, p[i].value, p[i].freeze, p[i].min, p[i].max);
     }
}

private define linear_random (mn, mx)
{
   return mn + (mx - mn) * urand();
}

private define select_at_random (mn, mx)
{
   return linear_random (mn, mx);
}

private define random_params (p)
{
   variable i, n = get_num_pars ();
   
   for (i = 0; i < n; i++)
     {
	if (p[i].freeze == 0 and p[i].tie == 0)
	  {
	     p[i].value = select_at_random (p[i].min, p[i].max);
	  }
     }
   
   return p;
}

private define test_params (p, evalfun)
{
   variable info;
   
   set_params (p);
   () = @evalfun (&info);
   
   return info.statistic;
}

public define monte_carlo (num, evalfun)
{
   variable bkp_verbose = Fit_Verbose;
   Fit_Verbose = -1;
   
   variable p, best_p, stat, best_stat;
   
   p = get_params ();
   best_stat = test_params (p, evalfun);
   best_p = get_params ();
         
   variable i = 0;
   
   loop (num)
     {
	if (bkp_verbose >= 0 and Isis_Batch_Mode == 0)
	  () = fprintf (stderr, "monte_carlo:  %3d/%3d\tbest=%7.3f\r", 
			i, num, best_stat);
	
	stat = test_params (random_params (p), evalfun);
	
	if (stat < best_stat)
	  {
	     best_stat = stat;
	     best_p = get_params ();
	  }
	
	i++;
     }
   
   set_params (best_p);
   
   Fit_Verbose = bkp_verbose;      
}
