% Maximum-Likelihood fit-statistic

static variable Log_Factorial = Double_Type [100];
static define setup_factorial ()
{
   variable x = 1.0;
   Log_Factorial[0] = 0;
   _for (1, 99, 1)
     {
	variable i = ();
	x = x * i;
	Log_Factorial[i] = log(x);
     }
}

setup_factorial ();
define log_factorial (n)
{
   variable x = Double_Type [length(n)];
   variable i = where (n < 100);
   x[i] = Log_Factorial[n[i]];
   i = where (n >= 100);
   n = n[i];
   x[i] = n * log(n) - n;
   return x;
}


static define like_ml_function (y, fx, w)
{
   variable i = where (fx > 0);
   y = y[i];
   fx = fx[i];
   
   y = typecast (y+0.5, Int_Type);
   
   variable v = -(y * log(fx) - fx - log_factorial(y));
     
   return (v, sum(v));
}

static define like_ml_report (stat, npts, nvpars)
{
   variable s = sprintf ("   Likelihood = %0.4g\n", exp(-stat/npts));
   return s;
}

add_slang_statistic ("like", &like_ml_function, &like_ml_report);
set_fit_statistic ("like");

