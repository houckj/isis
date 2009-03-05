static variable MODULE_NAME = "template";
prepend_to_slang_load_path (".");
require (MODULE_NAME);

% Useful functions
static define print (x)
{
   x = string (x);
   () = fputs (x, stdout);
   () = fflush (stdout);
}

static define failed ()
{
   variable s = __pop_args (_NARGS);
   s = sprintf (__push_args(s));
   () = fprintf (stderr, "Failed: %s\n", s);
   exit (1);
}

static variable Random_Number = _time ();
static define urand_1 (x)
{
   Random_Number = typecast (Random_Number * 69069UL + 1013904243UL, UInt32_Type);
   return Random_Number/4294967296.0;
}
static define urand (n)
{
   if (n == 0)
     return Double_Type[0];

   return array_map (Double_Type, &urand_1, [1:n]);
}

vmessage ("testing %s routines\n", MODULE_NAME);
_debug_info = 1;

% Regression tests go here


% End of regression tests
message ("Ok\n");
exit (0);
