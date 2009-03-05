
variable EPSILON = 1.e-12;

define plaw (alph, lo, hi)
{
   variable oma = 1.0 + alph;
   
   if (abs(alph + 1.0) > 1.e4 * EPSILON)
     return (1.0 / oma) * ((Const_hc/lo)^oma - (Const_hc/hi)^oma);
   else
     return log (lo / hi);
}

define bpl_fit (lo, hi, par)
{
   variable norm, ebrk, p1, p2, emid;
   
   norm = par[0];
   ebrk = par[1];
   p1   = par[2];
   p2   = par[3];

   emid = Const_hc / (0.5 * (lo + hi));
   
   variable i, j, val;
   i = where (emid <= ebrk);
   j = where (emid > ebrk);

   val = Float_Type [length(lo)];
   
   val[i] = norm * plaw (p1, lo[i], hi[i]);
   val[j] = norm * plaw (p2, lo[j], hi[j]) * (ebrk^(p2 - p1));
   
   return val;
}

add_slang_function ("bpl", ["norm", "ebrk", "a_lo", "a_hi"]);

define bpl_set_default (offset)
{
   set_par (offset+1, 0.01, 0, 1.e-6,  1.0);
   set_par (offset+2,  1.0, 0,   0.1, 10.0);
   set_par (offset+3, -0.5, 0,  -3.0,  3.0);
   set_par (offset+4, -1.7, 0,  -3.0,  3.0);
}
