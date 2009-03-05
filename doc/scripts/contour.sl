%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PGPLOT has several contouring functions:
%
% PGCONB -- contour map of a 2D data array, with blanking
% PGCONF -- fill between two contours
% PGCONL -- label contour map of a 2D data array
% PGCONS -- contour map of a 2D data array (fast algorithm)
% PGCONT -- contour map of a 2D data array (contour-following)
% PGCONX -- contour map of a 2D data array (non rectangular)
%
%
% the corresponding ISIS wrappers are
% _pgconb, _pgconf, _pgconl, _pgcons, _pgcont
%
% At the moment, _pgconx isn't supported.

typedef struct
{
   min, max, num
} 
Param_Type;

public define plot_contours (a, px, py, levels)
{   
   % 1. define the plot coordinate system and draw the box
   % 2. plotting contours one at a time allows coloring
   %    them individually (see _pgsci()).
   
   _pgenv (px.min, px.max, py.min, py.max, 0, 0);
   _pgswin (1, px.num, 1, py.num);

   variable tr = [0, 1.0, 0, 0, 0, 1.0];   
   variable i, n = length(levels);

   _for (0, n-1, 1)
     {
        i = ();
	% negative number of contours means contour lines
	% are drawn with the current line attributes 
	% (color index, style, width)
	_pgcont (a, 1, 1, px.num, py.num, levels[i], tr, -1);
     }
}

public define fake_data (px, py)
{
   variable a = Float_Type [px.num, py.num];
   variable i, j, r, f;
   
   f = PI * py.max / px.max;   
   for (i = 0; i < px.num; i++)
     {
	for (j = 0; j < py.num; j++)
	  {
	     r = sqrt ( (i - px.num/2.0)^2 + (j - py.num/2.0)^2 );
	     a[i,j] = r - 0.1 * px.num * (sin(i*f/(j+1.0)))^2;
	  }
     }   
   return a;
}

public define do_demo ()
{
   variable px = @Param_Type;
   variable py = @Param_Type;

   % Define axis ranges and resolution
   px.min = 1.0;
   px.max = 175.0;
   px.num = 175;
   
   py.min = 1.0;
   py.max = 175.0;
   py.num = 175;
   
   % Define contour levels
   variable levels = [0.0:50.0:5.0];
   
   plot_contours (fake_data (px, py), px, py, levels);
   _pglab ("X", "Y", "Pretty Contours");
   
}

() = _pgopen ("contour.ps/vcps");
_pgpap (4.0, 1.0);      % 4 inches, aspect ratio = 1.
do_demo ();
_pgclos();

