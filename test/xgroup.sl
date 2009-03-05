% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing xgroup.... ");

% First, test conversions between the isis
% and OGIP/xspec grouping conventions

define prt (g)
{
   variable ag = array_map (String_Type, &string, g);
   () = fprintf (stdout, "%s\n", strjoin(ag, " "));
}

define null_test (g)
{
   variable i2x = i2x_group(g);
   variable xig = x2i_group(i2x);

   if (all(xig == g))
     return;

   if (all(xig == -g))
     return;

   failed ("null_test");
   prt(g);
   prt(i2x);
   prt(xig);
}

null_test([ 1 ]);
null_test([ 1,  1]);
null_test([ 1, -1]);
null_test([-1, -1]);
null_test([ 1,  1,  1]);
null_test([-1,  1,  1]);
null_test([ 1, -1,  1]);
null_test([ 1,  1, -1]);
null_test([-1, -1,  1]);
null_test([-1,  1, -1]);
null_test([ 1, -1, -1]);
null_test([-1, -1, -1]);

null_test(ones(20));
null_test(-ones(20));

define random_grouping (len)
{
   variable r = urand(len);
   variable g = ones(len);
   g[where(r < 0.5)] = -1;
   return g;
}

loop (100)
{
   null_test (random_grouping (128));
}

% Now test FITS grouping I/O functions

#ifdef __CFITSIO__

define test_fitsio ()
{
   variable file = "data/pi.fits";
   variable file_cpy = "tmp_pi.fits";

   () = system (sprintf("/bin/cp %s %s; chmod u+w %s", file, file_cpy, file_cpy));

   try
     {
        () = load_data (file);

        group_data (1, 4);
        variable ig = get_data_info (1).rebin;

        regroup_file (1, file_cpy);
        variable xg = fits_read_col (file_cpy, "grouping");

        if (andelse
            {any (ig != x2i_group(xg))}
              {any (ig != -x2i_group(xg))})
          {
             failed ("MISMATCH (regroup_file)");
          }

        group_data (1,1);
        use_file_group (1, file_cpy);
        variable fg = get_data_info (1).rebin;

        if (andelse
            {any (ig != fg)}
              {any (ig != -fg)})
          {
             failed ("MISMATCH (use_file_group)");
          }
     }
   finally
     {
        () = remove (file_cpy);
     }
}

test_fitsio ();
#endif

msg("ok\n");

