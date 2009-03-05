require ("xfig");
require ("maplib");
require ("structfuns");

_debug_info=1;
static variable Coastline_Data_Type = struct
{
   lon, lat, next
};

static define new_contour (last)
{
   variable next = @Coastline_Data_Type;
   if (last != NULL)
     last.next = next;
   return next;
}

define read_coastline_data (file)
{
   % The coastline data file has contours separated by
   % # -b
   variable list = new_contour (NULL);
   variable this = list;
   variable i = 0;
   variable n = 100;
   variable lon = Double_Type[n], lat = Double_Type[n];
   foreach (fopen (file, "r"))
     {
	variable line = ();
	variable x, y;

	if (2 != sscanf (line, "%g %g", &x, &y))
	  {
	     if (i > 0)
	       {
		  this.lat = lat[[0:i-1]];
		  this.lon = lon[[0:i-1]];
		  this = new_contour (this);
		  i = 0;
	       }
	     continue;
	  }
	
	if (i >= n)
	  {
	     lat = [lat, lat];
	     lon = [lon, lon];
	     n *= 2;
	  }
	
	lat[i] = y;
	lon[i] = x;
	
	i++;
     }
   
   if (i > 0)
     {
	this.lat = lat[[0:i-1]];
	this.lon = lon[[0:i-1]];
     }
   return list;
}

static define read_cities (file)
{
   variable lines = fgetslines (fopen (file,"r"));
   variable n = length (lines);
   variable city = String_Type[n];
   variable country = String_Type[n];
   variable lon = Float_Type[n];
   variable lat = Float_Type[n];
   
   variable i = 0;
   foreach (lines)
     {
	variable line = ();
	if (line[0] == '#')
	  continue;
	variable xcity, xcountry, xlon, xlat;
	if (4 != sscanf (line, "%[^,], %s\t%f\t%f", &xcity, &xcountry, &xlon, &xlat))
	  continue;
	
	city[i] = xcity;
	country[i] = xcountry;
	lon[i] = xlon;
	lat[i] = xlat;
	i++;
     }

   variable s = struct
     {
	city, country, lon, lat
     };
   i = [0:i-1];
   s.city = city[i];
   s.country = country[i];
   s.lon = lon[i];
   s.lat = lat[i];
   return s;
}

private define local_maplib_new (name)
{
   return maplib_new (name);

   if (name == "hammer")
     {
	variable m = struct 
	  {
	     name, lon0, lat0, sphere, lambert
	  };
	m.name = name;
	m.lon0 = 0.0;
	m.lat0 = 90.0;
	m.sphere = maplib_new ("sphere");
	m.lambert = maplib_new ("lambert");
	return m;
     }
   
   return maplib_new (name);
}

private define local_maplib_project (m, x, y)
{
   return maplib_project (m, x, y);

   if (m.name == "hammer")
     {
	variable sphere = m.sphere;
	sphere.lon0 = m.lon0;
	sphere.lat0 = m.lat0;
	(x,y) = maplib_project (sphere, x, y);
	(x,y) = maplib_project (m.lambert, 0.5*x, y);
	return 2.0*x, y;
     }
   return maplib_project (m, x, y);
}

define get_data_limits (c, m)
{
   variable x, y;
   variable min_x, max_x, min_y, max_y;

   min_x = 1e37;
   max_x = -1e37;
   min_y = 1e37;
   max_y = -1e37;

   variable cc;
   foreach (c) using ("next")
     {
	cc = ();
	(x, y) = local_maplib_project (m, cc.lon, cc.lat);
	variable i = where ((x==x) and (y==y));
	min_x = min([min_x, x[i]]);
	max_x = max([max_x, x[i]]);
	min_y = min([min_y, y[i]]);
	max_y = max([max_y, y[i]]);
     }
   return min_x, max_x, min_y, max_y;
}

define plot_coastline_data (c, fun, m, lon_min, lon_max, lat_min, lat_max)
{
   variable lon = c.lon, lat = c.lat;

   if (0 == length (where ((lon >= lon_min) and (lon < lon_max))))
     return;

   if (0 == length (where ((lat >= lat_min) and (lat < lat_max))))
     return;
   
   (@fun) (local_maplib_project (m, lon, lat));
}

#iffalse
static define add_latitude_longitude_tics (w, m, x0, x1, y0, y1)
{
   variable lat0, lon0, lat1, lon1, lat, lon, x, y;
   
   (lon0, lat0) = local_maplib_deproject (m, x0, y0);
   (lon1, lat1) = local_maplib_deproject (m, x0, y1);
   if (lon0 == NULL)
     {
	x = 0.5*(x0+x1);
	while ()
     }
}
#endif

static define make_ylabel (w, lon, x, y, x0, x1, y0, y1)
{
   variable i = where (x == x);
   x = x[i]; y = y[i];
   if (length(x) <= 2)
     return;
   variable sgn = (x >= x1);
   variable dsgn = shift(sgn,1)-sgn;
   dsgn[-1] = dsgn[-2];
   i = where (dsgn != 0);
   if (0 == length(i))
     return;

   i = i[0];
   y = 0.5*(y[i]+y[i+1]);
   if ((y < y0) or (y > y1))
     return;
   variable obj = xfig_new_text (sprintf ("%g", lon));
   xfig_plot_add_object (w, obj, x1, y, -0.6, 0);
}

private define plot_lines (w, x, y, x0, x1, y0, y1)
{
   variable dx = shift (x, 1) - x;
   dx[-1] = 0;
   variable dy = shift (x, 1) - x;
   dy[-1] = 0;
   variable dr = hypot (dx, dy);
   variable i = where (dr > 0.05 * (hypot (x1-x0, y1-y0))) + 1;
   x[i] = _NaN;
   y[i] = _NaN;
   xfig_plot_lines (w, x,y);
}

define plot_data (c, cities, m, dlon, dlat, x0, x1, y0, y1, theta, file)
{
   variable lon0 = m.lon0, lat0 = m.lat0;
   %variable lon0 = 0, lat0 = 0;
   theta *= PI/180.0;

   x0 = max([-350,x0]); x1 = min([350,x1]);
   y0 = max([-350,y0]); y1 = min([350,y1]);

   variable fp = xfig_create_file (file);
   variable w = xfig_plot_new (18, 18);
   xfig_plot_define_world (w, x0, x1, y0, y1);
   %xfig_plot_add_x_axis (w, 0, "X");
   %xfig_plot_add_y_axis (w, 0, "Y");
   %add_latitude_longitude_tics (w, m, x0, x1, y0, y1);
   xfig_plot_set_line_color (w, "blue4");
   xfig_plot_set_line_thickness (w, 1);
   variable lon, lat, lons, lats, x, y;
   variable i;
   
   lats = [-90:90:0.1];

   _for (-180, 180, 10)
     {
	lon = ();
	lons = lon + Double_Type[length(lats)];
	i = where ((abs(lons-lon0) <= dlon) and (abs(lats-lat0)<=dlat));
	(x,y) = local_maplib_project (m, lons[i], lats[i]);
	(x,y) = maplib_rotate2d (x, y, theta);
	%make_ylabel (w, lon, x, y, x0, x1, y0, y1);
	plot_lines (w, x, y, x0, x1, y0, y1);
     }
   
   lon = [-180:180.1:0.1];
   _for (-90, 90, 10)
     {
	lat = ();
	lat = lat + Double_Type[length(lon)];
	i = where ((abs(lon-lon0) <= dlon) and (abs(lat-lat0)<=dlat));
	(x,y) = local_maplib_project (m, lon[i], lat[i]);
	(x,y)=maplib_rotate2d (x, y, theta);
	plot_lines (w, x, y, x0, x1, y0, y1);
     }
   
   xfig_plot_set_line_color (w, "red");
   xfig_plot_set_line_thickness (w, 3);
   xfig_plot_inc_line_depth (w, -1);
   foreach (c) using ("next")
     {
	variable cc = ();
	variable cc_lon = cc.lon;
	variable cc_lat = cc.lat;

	i = where ((abs(cc_lon-lon0) <= dlon) and (abs(cc_lat-lat0)<=dlat));
	(x,y) = local_maplib_project (m, cc_lon[i], cc_lat[i]);
	(x,y)=maplib_rotate2d (x, y, theta);

	plot_lines (w, x, y, x0, x1, y0, y1);
	%xfig_plot_lines (w, x,y);
     }
   xfig_plot_set_line_color (w, "black");
   xfig_plot_set_point_size (w, 5);
   
   if (cities != NULL)
     {
	variable city_x, city_y;
	(city_x, city_y) = local_maplib_project (m, cities.lon, cities.lat);
	(city_x, city_y) = maplib_rotate2d (city_x, city_y, theta);
	i = where ((city_x > x0) and (city_x < x1) 
		   and (city_y > y0) and (city_y < y1));
	xfig_plot_points (w, city_x[i], city_y[i]);
	foreach (i)
	  {
	     i = ();
	     variable obj = xfig_new_text (sprintf ("\\Large\\bf %s", cities.city[i]));
	     xfig_plot_add_object (w, obj, city_x[i], city_y[i], -0.6, 0);
	  }
     }
   xfig_render_object (w, fp);
   xfig_close_file (fp);
}

define fixup_n (m, lon0, lat0)
{
   variable theta0 = lat0 * PI/180.0;
   variable phi0 = lon0 * PI/180.0;
   variable phi = m.p_lon * PI/180.0;
   variable theta = m.p_lat * PI/180.0;
   variable R = m.p_len;
   variable s0 = vector (cos(theta0)*cos(phi0), 
			 cos(theta0)*sin(phi0), 
			 sin(theta0));

   variable p = vector (R*cos(theta)*cos(phi),
			R*cos(theta)*sin(phi), R*sin(theta));
   variable n = p-s0;
   normalize_vector (n);
   theta = asin (n.z);
   phi = atan2 (n.y, n.x);
   m.n_lat = theta * 180.0/PI;
   m.n_lon = phi * 180.0/PI;
}

define cairo_map (c, cities, file)
{
   variable m = local_maplib_new ("plane");
   variable athens_lon = 23.44;
   variable athens_lat = 38.00;
   variable cairo_lon = 31.15;
   variable cairo_lat = 30.03;
   variable R = 1.35;

   m.p_lon = cairo_lon;
   m.p_lat = cairo_lat;
   m.p_len = R;
   m.lon0 = athens_lon;
   m.lat0 = athens_lat;
   
   variable i = where ((cities.city == "Cairo") 
		       or (cities.city == "Athens") or (cities.city == "Madrid"));

   struct_filter (cities, i);

   fixup_n (m, athens_lon, athens_lat);
   plot_data (c, cities, m, 180, 90, -12, 12, -10,14, -55, file);
   exit (0);
}

define make_plots ()
{
   variable c = read_coastline_data ("coastline.dat");
   variable cities = read_cities ("cities.dat");
   variable m;
   variable x0, x1, y0, y1;
#iffalse
   cairo_map (c, cities, "cairo.ps");

   m = local_maplib_new ("linear"); m.x0 = 0; m.y0 = 0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "linear.ps");
#endif


   cities = NULL;

   variable lon0 = -70, lat0 = 42;

   m = local_maplib_new ("azeqdist"); m.lon0 = lon0; m.lat0 = lat0; m.beta=0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   %c = NULL;
   %plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "lambert.ps");
   plot_data (c, cities, m, 3600, 1800, x0, x1, y0, y1, 0, "mercator.ps");
   return;

   m = local_maplib_new ("mercator"); m.lon0 = lon0; m.lat0 = lat0; m.beta=90;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   %c = NULL;
   %plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "lambert.ps");
   plot_data (c, cities, m, 3600, 1800, x0, x1, y0, y1, 0, "mercator.ps");
   return;

   m = local_maplib_new ("sinusoidal"); m.lon0 = lon0; m.lat0 = lat0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   %c = NULL;
   %plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "lambert.ps");
   plot_data (c, cities, m, 3600, 1800, x0, x1, y0, y1, 0, "bonne.ps");
   return;

   m = local_maplib_new ("bonne"); m.lon0 = lon0; m.lat0 = lat0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   %c = NULL;
   %plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "lambert.ps");
   plot_data (c, cities, m, 3600, 1800, x0, x1, y0, y1, 0, "bonne.ps");
   return;

   m = local_maplib_new ("hammer"); m.lon0 = lon0; m.lat0 = lat0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   %c = NULL;
   %plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "lambert.ps");
   plot_data (c, cities, m, 3600, 1800, x0, x1, y0, y1, 0, "hammer.ps");
   return;

   m = local_maplib_new ("lambert"); m.lon0 = lon0; m.lat0 = lat0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   %c = NULL;
   %plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "lambert.ps");
   plot_data (c, cities, m, 3600, 1800, x0, x1, y0, y1, 0, "lambert.ps");
   return;

   m = local_maplib_new ("stereo"); m.lon0 = 0; m.lat0 = 0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "stereo.ps");
   
   return;


   m = local_maplib_new ("hammer"); m.lon0 = 0; m.lat0 = 0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   plot_data (c, cities, m, 180, 90, x0, x1, y0, y1, 0, "hammer.ps");

   return;
   m = local_maplib_new ("gnomic"); m.lon0 = 0; m.lat0 = 0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   plot_data (c, cities, m, 80, 80, x0, x1, y0, y1, 0, "gnomic.ps");

   m = local_maplib_new ("ortho"); m.lon0 = 0; m.lat0 = 0;
   (x0, x1, y0, y1) = get_data_limits (c, m);
   x0 -= 10; x1 += 10; y0 -= 10; y1 += 10;
   plot_data (c, cities, m, 90, 90, x0, x1, y0, y1, 0, "ortho.ps");

}
make_plots ();
exit (0);
   
