static variable MODULE_NAME = "maplib";
prepend_to_slang_load_path (".");
set_import_module_path (".:" + get_import_module_path ());

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

% The idea behind these tests is that the specific cases will be compare 
% to a more general projection with parameters set to mimick to specific 
% case.

static define generate_lon_lats (n, lon0, lat0, dlon, dlat)
{
   variable lon, lat;
   
   lon = dlon * (2.0*urand (n)-1.0);
   lat = dlon * (2.0*urand (n)-1.0);
   
   variable m = maplib_new ("sphere");
   m.lon0 = 0;
   m.lat0 = 0;
   m.lon1 = lon0;
   m.lat1 = lat0;
   
   (lon, lat) = maplib_project (m, lat, lon);

   variable i, j;
   i = where (lon > 180.0);
   lon[i] -= 360.0;
   j = where (lon < -180);
   lon[j] += 360.0;
   if (length (i) or length (j))
     () = fprintf (stderr, "*** maplib_project sphere returned unnormalized lons\n");

   i = where (lat > 90.0);
   lat[i] = 180.0 - lat[i];
   j = where (lat < -90.0);
   lat[j] = -180.0 - lat[j];
   if (length (i) or length (j))
     () = fprintf (stderr, "*** maplib_project sphere returned unnormalized lats\n");

   return lon, lat;
}

static define check_a_diff (x, x1, tol, lon, lat, axis, name, fun)
{
   % FIXME!! I need to device a better check here.  The absolute difference
   % is not meaningful for large numbers.
   variable dx = abs(x-x1);
   variable i = where (isnan(x) != isnan(x1));
   if (length (i))
     {
	() = fprintf (stderr, "*** WARNING: NaN issue seen in by %s %s %d/%d times\n",
		      fun, name, length(i), length(x));
	i = i[0];
	() = fprintf (stderr, "             first at (%g,%g)\n", lon[i], lat[i]);
     }
   i = (where (dx > tol));

   if (0 == length(i))
     return;

   () = fprintf (stderr, "*** WARNING: %s tolerance for %s %s exceeded (%e>%g) %d/%d times\n",
		 axis, fun, name, max(dx), tol, length(i), length(x));
   i = wherefirst (max(dx) == dx);
   () = fprintf (stderr, "             %g vs %g at (%g,%g)\n", x[i], x1[i], lon[i], lat[i]);
}

static define check_diffs (x, y, x1, y1, tol_x, tol_y, lon, lat, name, fun)
{
   check_a_diff (x, x1, tol_x, lon, lat, "X", name, fun);
   check_a_diff (y, y1, tol_y, lon, lat, "Y", name, fun);
}

     
static define test_gnomic ()
{
   variable g = maplib_new ("gnomic");
   variable p = maplib_new ("plane");
   p.lon0 = g.lon0;
   p.lat0 = g.lat0;
   p.p_lon = 0.0;
   p.p_lat = 0.0;
   p.p_len = 0.0;
   
   variable lat, lon;
   (lon,lat) = generate_lon_lats (10000, p.lon0, p.lat0, 45, 45);
   variable x,y,x1,y1;
   
   variable tol_x = 1e-12;
   variable tol_y = 1e-12;

   loop (2)
     {
	(x,y) = maplib_project (p, lon, lat);
	(x1,y1) = maplib_project (g, lon, lat);
	check_diffs (x,y,x1,y1,tol_x,tol_y,lon,lat, "gnomic", "maplib_project");
	
	lon = typecast(lon,Float_Type);
	lat = typecast(lat,Float_Type);
     }
}

static define test_ortho ()
{
   variable g = maplib_new ("ortho");
   variable p = maplib_new ("plane");

   g.lon0 = 0;
   g.lat0 = 0;
   p.lon0 = g.lon0;
   p.lat0 = g.lat0;
   p.p_lat = p.lat0;
   p.p_lon = p.lon0;
#iffalse
   p.p_lon = 180 - p.lat0;
   if (p.p_lon > 180.0)
     p.p_lon -= 360.0;
#endif
   p.p_len = 1e30;	       %  Inf
   
   variable lat, lon;
   (lon,lat) = generate_lon_lats (10000, p.lon0, p.lat0, 89, 89);
   variable x,y,x1,y1;
   
   variable tol_x = 1e-12;
   variable tol_y = 1e-12;

   % Note: tests involving floats instead of doubles will give errors
   % because floats lack the necessary precision at larger angles.
   loop (2)
     {
	(x,y) = maplib_project (g, lon, lat);
	(x1,y1) = maplib_project (p, lon, lat);
	check_diffs (x,y,x1,y1,tol_x,tol_y,lon,lat, "ortho", "maplib_project");

	(x1, y1) = maplib_deproject (g, x, y);
	check_diffs (lon, lat, x1, y1, tol_x, tol_y, lon, lat, "ortho", "maplib_deproject");

%	(x1, y1) = maplib_project (g, typecast(lon,Float_Type), typecast(lat,Float_Type));
%	check_diffs (x, y, x1, y1, tol_x, tol_y, lon, lat, "ortho float vs double", "maplib_project");
	break;
	  
	lon = typecast(lon,Float_Type);
	lat = typecast(lat,Float_Type);
	tol_x = 1e-5;
	tol_y = 1e-5;
     }
}

static define test_stereo ()
{
   variable g = maplib_new ("stereo");
   variable p = maplib_new ("plane");
   p.lon0 = g.lon0;
   p.lat0 = g.lat0;
   p.p_lat = -p.lat0;
   p.p_lon = 180 - p.lon0;
   if (p.p_lon > 180.0)
     p.p_lon -= 360.0;
   p.p_len = 1.0;
   
   variable lat, lon;
   (lon,lat) = generate_lon_lats (10000, p.lon0, p.lat0, 45, 45);
   variable x,y,x1,y1;
   
   variable tol_x = 1e-12;
   variable tol_y = 1e-12;

   loop (2)
     {
	(x,y) = maplib_project (p, lon, lat);
	(x1,y1) = maplib_project (g, lon, lat);
	check_diffs (x,y,x1,y1,tol_x,tol_y,lon, lat, "stereo", "maplib_project");

	lon = typecast(lon,Float_Type);
	lat = typecast(lat,Float_Type);
     }
}

private define test_sphere ()
{
   variable m = maplib_new ("sphere");
   variable lon0, lat0, lon1, lat1, lon, lat;

   (lon,lat) = generate_lon_lats (10000, 0, 0, 180, 90);
   % Check unit transformation
   foreach lon0 ([0, 90, 180, -180, -90])
     {
	foreach lat0 ([90, 0, -90])
	  {
	     m.lon0 = lon0; m.lat0 = lat0;
	     m.beta = 31;
	     m.lon1 = 45; m.lat1 = 45;
	     loop (2)
	       {
		  variable tol = 1e-12;
		  (lon1, lat1) = maplib_reproject (m, m, lon, lat);
		  check_diffs (lon,lat,lon1,lat1,tol,tol,lon,lat, "sphere", "maplib_reproject");
		  lon = typecast (lon, Float_Type);
		  lat = typecast (lat, Float_Type);
		  tol = 1e-5;
	       }
	  }
     }

   m.lon0 = 0; m.lat0 = 90;
   m.lon1 = 90; m.lat1 = 90;
   m.beta = 45;
   (lon1, lat1) = maplib_project (m, 0, 90);
   check_diffs (90+m.beta, 90, lon1, lat1, 1e-12, 1e-12, 0, 90, "sphere", "maplib_project");
}

static define test_plane ()
{
   variable p = maplib_new ("plane");
   p.lon0 = 20;
   p.lat0 = 30;
   p.p_lon = 0;
   p.p_lat = 0;

   variable lat, lon;
   (lon,lat) = generate_lon_lats (10000, p.lon0, p.lat0, 19, 31);
   foreach ([0, 0.5, 1, 100])
     {
	p.p_len = ();
	variable x,y,lon1,lat1;
   
	variable tol_x = 1e-12;
	variable tol_y = 1e-12;

	(x,y) = maplib_project (p, lon, lat);
	(lon1, lat1) = maplib_deproject (p, x, y);
	check_diffs (lon,lat,lon1,lat1,tol_x,tol_y,lon,lat, "plane", 
		     sprintf ("maplib_project/deproject@plen=%g",p.p_len));
     }
}

define test_generic (name, lon0, lat0, dlon, dlat)
{
   variable m = maplib_new (name);
   m.lon0 = lon0;
   m.lat0 = lat0;

   variable lat, lon;
   (lon,lat) = generate_lon_lats (10000, lon0, lat0, 180, 90);

   variable tol_x = 1e-12;
   variable tol_y = 1e-12;

   variable lon1, lat1;
   (lon1, lat1) = maplib_deproject (m, maplib_project (m, lon, lat));

   check_diffs (lon,lat,lon1,lat1,tol_x,tol_y,lon,lat, m.name,
		"maplib_project/deproject");
}

test_generic ("sinusoidal", 20, 0, 180, 90);
test_generic ("bonne", 20, 45, 180, 90);
test_generic ("mercator", 20, 45, 180, 89);
   
test_plane ();
test_gnomic ();
test_ortho ();
test_stereo ();
test_sphere ();

% End of regression tests
message ("Ok\n");
exit (0);
