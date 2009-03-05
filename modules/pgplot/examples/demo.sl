#!/usr/bin/env slsh
require ("pgplot");

static define demo1 ()
{
   variable xs = [0:60]*0.1;
   variable ys = xs^2;

   % Call cpgenv to specify the range of the axes and to draw a box, and
   % cpglab to label it. The x-axis runs from 0 to 10, and y from 0 to 20.
   _pgenv (0.0, 10.0, 0.0, 20.0, 0, 1);
   _pglab("(x)", "(y)", "PGPLOT Example 1: y = x\\u2\\d");
  
   % Mark integer-valued points using symbol number 9.
   variable i = [10:60:10];
   _pgpt(xs[i], ys[i], 9);
   
   % Draw smooth line through all points
   _pgline(xs, ys);
}

%--------------------------------------------------------------------
% Demonstration function for PGPLOT contouring routines.
%--------------------------------------------------------------------

static define demo2()
{
   variable nx = 40, ny = 40;  
   variable tr = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
  
   % Compute a suitable function.
   % First make a grid of points.  The maplib module has "meshgrid".  Here
   % for clarity, a loop will be used
   % Note: we want to create a 2d array in row major order, indexed [j,i]
   % so that i labels the x coord and j the y.
   % We want: I_ji = i for any j, and J_ji = j for any i
   % That way Img[J,I] has the proper indices
   variable I = Int_Type[ny, nx];
   variable J = Int_Type[ny, nx];
   variable i, j;
   i = [1:nx];
   _for (0, ny-1, 1)
     {
	j = ();
	I[j,*] = i;
     } 
   j = [1:ny];
   _for (0, nx-1, 1)
     {
	i = ();
	J[*,i] = j;
     }
   variable x = tr[0] + tr[1]*I + tr[2]*J;
   variable y = tr[3] + tr[4]*I + tr[5]*J;
   variable f = cos(0.3*sqrt(x*2)-0.13333*y)*cos(0.13333*x)+(x-y)/nx;
   variable fmin = min (f);
   variable fmax = max (f);

   % Clear the screen. Set up window and viewport.
  
   _pgpage();
   _pgsvp(0.05, 0.95, 0.05, 0.95);
   _pgswin(1.0, nx, 1.0, ny);
   _pgbox("bcts", 0.0, 0, "bcts", 0.0, 0);
   _pgmtxt("t", 1.0, 0.0, 0.0, "Contouring using cpgcont()");
  
   % Draw the map. _pgcont is called once for each contour, using
   % different line attributes to distinguish contour levels.
  
   variable lw, alev, ls, ci;
   _pgbbuf();
   _for (1, 20, 1)
     {
	i = ();
	alev = fmin + i*((fmax-fmin)/20.0);
	if (i mod 5) lw = 1; else lw = 3;
	if (i < 10) ci = 2; else ci = 3;
	if (i < 10) ls = 2; else ls = 2;
	_pgslw(lw);
	_pgsci(ci);
	_pgsls(ls);
	_pgcont(f, 0, ny-1, 0, nx-1, alev, tr, -1);
     }

   _pgslw(1);
   _pgsls(1);
   _pgsci(1);
   _pgebuf();
}

static define demo3 ()
{  
   variable n1 = [3, 4, 5, 5, 6, 8];
   variable n2 = [1, 1, 1, 2, 1, 3];

   variable lab = ["Fill style 1 (solid)",
		  "Fill style 2 (outline)",
		  "Fill style 3 (hatched)",
		  "Fill style 4 (cross-hatched)"];
  
   % Initialize the viewport and window.

   _pgbbuf();
   _pgsave();
   _pgpage();
   _pgsvp(0.0, 1.0, 0.0, 1.0);
   _pgwnad(0.0, 10.0, 0.0, 10.0);
  
   % Label the graph.

   _pgsci(1);
   _pgmtxt("T", -2.0, 0.5, 0.5, 
	   "PGPLOT fill area: routines cpgpoly(), cpgcirc(), cpgrect()");
  
   % Draw assorted polygons.
   variable i, j, k;
   
   for (k=1; k<5; k++) 
     {
	_pgsci(1);
	variable y0 = 10.0 -2.0*k;
	_pgtext(0.2, y0+0.6, lab[k-1]);
	_pgsfs(k);
	for (i=0; i<6; i++) 
	  {
	     _pgsci(i+1);
	     variable t = n2[i]*(2*PI*[0:n1[i]])/n1[i];
	     variable x = (i+1) + 0.5*cos(t);
	     variable y = y0 + 0.5*sin(t);
	     _pgpoly(n1[i], x, y);
	  }
	_pgsci(7);
	_pgshs(0.0, 1.0, 0.0);
	_pgcirc(7.0, y0, 0.5);
	_pgsci(8);
	_pgshs(-45.0, 1.0, 0.0);
	_pgrect(7.8, 9.5, y0-0.5, y0+0.5);
     }
  _pgunsa();
  _pgebuf();
}

static define main ()
{
   if (_pgopen ("?") <= 0)
     {
	() = fprintf (stderr, "Unable to open pgplot\n");
	exit (1);
     }
   _pgask(1);
   demo1();
   demo2();
   demo3();
   _pgend ();
   exit (0);
}

main ();
