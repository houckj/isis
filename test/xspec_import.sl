
private variable headas = getenv("HEADAS");
if (headas == NULL || NULL == stat_file (headas))
{
   () = fprintf (stderr, "skipping xspec module import\n");
   exit(0);
}

() = fprintf (stderr, "testing xspec module import.... ");
require ("xspec");
fit_fun ("mekal");
() = eval_fun(1,2);

variable c0 = [100.0, 0.1, 0.1], c1 = Double_Type[3];
xspec_set_cosmo (c0[0], c0[1], c0[2]);
(c1[0], c1[1], c1[2]) = xspec_get_cosmo();
c0 = typecast (c0, Float_Type);
c1 = typecast (c1, Float_Type);
if (any (c0 != c1))
   throw ApplicationError;

() = fprintf (stderr, "ok\n");
