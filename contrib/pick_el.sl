% Time-stamp: <2001-10-17 11:38:16 dph> 
% MIT Directory: ~dph/libisis
% File: pick_el.sl
% Author: D. Huenemoerder
% Original version: 2001.10.12
%
%%%%%%%%%%%%%%%%%%% pick ion from an index array %%%%%%%%%%%%%%%%%%%%%%
%+  
define pick_el()
  {
  % given array of line indices, a, return the subset of the array
  % which are lines of element, el, and ion, ion.
  % example:
  %     a = brightest(30, where(wl(6,7))); % brightest in a region
  %     a_Al = pick_el(a,Al);                % pick the iron lines in the list
  % If there aren't any, returns the array, [-1]

  % Note: this is different from:
  %      a = brightest(30, where( wl(6,7) and el_ion(Al) ) );
  %  which returns the 30 brightest set from only Al.    This would
  %  find any Al lines in the region, but they may not be among the
  %  brightest 30 features.

  variable a, ion, el;

  switch (_NARGS)

  {
  case 3:
  ion = ();
  el = ();
  a = ();
  }

  {
  case 2:
  el = ();
  a = ();
  ion = NULL;
  }

  {
  message("USAGE: g = pick_el(array, elem[, ion]);");
  }

  variable ii, s, n = length(a);
  variable result = -1;         % dummy value

  for (ii=0; ii<n; ii++)           % each line in array
    {
    s = line_info(a[ii]);

    switch (ion)

     {
     case NULL:
     if (s.Z == el)
       result = [result, a[ii]];
     }

     {
     % default case:
     if ( (s.Z == el) and (s.ion == ion) )
       result = [result, a[ii]];
     }

    } % end for each line.

  % were any lines found?
  if (length(result) > 1)
    {
    result = result[ [1:length(result)-1] ];  % truncate the -1
    }
  return result;
  }
%-  
