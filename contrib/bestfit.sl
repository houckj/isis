% bestfit.sl
% Author: John E. Davis <davis@space.mit.edu

_debug_info = 1;

static define do_fit_counts ()
{
   variable s;
   if (-1 == fit_counts (&s))
     {
	set_fit_method ("marquardt");
	
	vmessage ("**** Subplex failed.  Trying marquardt\n");
	vmessage ("Here are the values so far:\n");
	list_par;
	() = fit_counts (&s);
	set_fit_method ("subplex;maxnfe=500");
     }

   if (s.statistic == 0.0)
     {
	vmessage ("Rejecting statisic==0\n");
	return 1e38, 1e38;
     }
   
   return (s.statistic, s.statistic/(s.num_bins - s.num_variable_params));
}


static define run_fits (dir, n)
{   
   variable best_chisqr, best_rchisqr;
   variable chisqr, rchisqr;
   variable i;
   variable j = 0;
   variable s;
   variable best_parms;

   if (dir != NULL) 
     () = mkdir (dir, 0777);

   randomize ();
   (best_chisqr, best_rchisqr) = do_fit_counts ();
   best_parms = get_params ();

   if (dir != NULL)
     {
	save_par (sprintf ("%s/best.%d", dir, j));  
	j++;
     }

   _for (1, n, 1)
     {
	i = ();
	randomize ();
	vmessage ("Working on number %d (best chisqr=%g [%g per dof])", i, best_chisqr, best_rchisqr);
	(chisqr, rchisqr) = do_fit_counts ();
	if (chisqr < best_chisqr)
	  {
	     if (dir != NULL)
	       {
		  save_par (sprintf ("%s/best.%d", dir, j));  j++;
	       }
	     
	     best_chisqr = chisqr;
	     best_rchisqr = rchisqr;
	     best_parms = get_params ();
	  }
     }
   vmessage ("best chisqr = %g", best_chisqr);
   set_params (best_parms);
   () = eval_counts;
   return best_chisqr;
}

public define find_best_fit ()
{
   if (_NARGS == 2)
     return run_fits ();
   else
     {
	variable n = ();
	return run_fits (NULL, n);
     }
}

