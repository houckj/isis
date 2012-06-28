% -*- mode: SLang; mode: fold -*- 
() = evalfile ("inc.sl");
msg ("testing readcol.... ");

define write_file (file, nrows, ncols) %{{{
{
   variable r, c, x, fp;
   
   fp = fopen (file, "w");
   if (NULL == fp)
     failed ("opening file %s for writing", file);
   
   variable lastcol = Float_Type[nrows];
   
   _for (0, nrows-1, 1)
     {
	r = ();
	_for (0, ncols-1, 1)
	  {
	     c = ();
	     x = (c + 1)* 1.0;
	     () = fprintf (fp, "\t%f", x);
	     if (c == ncols-1)
	       lastcol[r] = x;
	  }
	() = fputs ("\n", fp);
     }
   
   () = fclose (fp);
}

%}}}

define mismatch (vec, val)
{
   return length (where (vec != val));
}

define try1_col (file, c)
{
   variable x = readcol (file, c);
   if (mismatch (x, c*1.0))
     failed ("read incorrect values");
}

define try2_cols (file, c1, c2)
{
   variable x, y;
   (x, y) = readcol (file, c1, c2);
   if (mismatch (x, c1*1.0) or mismatch (y, c2*1.0))
     failed ("read incorrect values");
}

variable file, nrows, ncols;
file = "isis_readcol.test";
nrows = 10;
ncols = 100;

write_file (file, nrows, ncols);
   
try1_col (file, 1);
try1_col (file, ncols);

try2_cols (file, ncols/2, ncols);
try2_cols (file, ncols, ncols/2);

variable x;

% this is a feature... 
x = readcol (file);
if (mismatch (x, 1.0))
  failed ("read incorrect values");

() = remove(file);

msg("ok\n");


