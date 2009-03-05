define msg (x)
{  
   x = string (x);
   () = fputs (x, stdout);
   () = fflush (stdout);
}

define failed ()
{
   variable s = __pop_args (_NARGS);
   s = sprintf (__push_args(s));
   () = fprintf (stderr, "Failed: %s\n", s);
   exit (1);
}

