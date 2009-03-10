private variable MODULE_NAME = "cfitsio";
prepend_to_slang_load_path (".");
set_import_module_path (".:" + get_import_module_path ());

require ("fitswcs");

private variable Failed = 0;
private define warn ()
{
   variable args = __pop_args (_NARGS);
   () = fprintf (stderr, "**** Warning: %s\n",
		 sprintf (__push_args (args)));
   Failed++;
}


private define test_wcs ()
{
   variable wcs = fitswcs_new (2);
   variable xmin = 100.0;
   variable ymin = 200.0;
   variable crpix = [ymin+10.0, xmin+20.0];

   wcs.cdelt = [1,1];
   wcs.crpix = crpix;
   wcs.crval = [7,8];
   wcs.pc=NULL;

   % Suppose that (x,y) coordinates are binned into an image of size
   % [M,M] using a grid [xmin, xmin+1, ..., xmin+(M-1)];
   variable M = 73;
   variable dx = 0.1, dy = 0.3;
   variable xgrid = xmin + dx*[0:M-1];
   variable ygrid = ymin + dy*[0:M-1];
   variable wcs_M = fitswcs_bin_wcs (wcs, ygrid, xgrid);
   
   % Now the reference pixel crpix will have the pixel coordinate i
   % given by
   %   crpix[1] = xmin + dx*(i-0.5)
   % or i = (crpix[1]-xmin)/dx + 0.5;
   
   variable ix = (crpix[1]-xmin)/dx + 0.5;
   if (fneqs (ix,wcs_M.crpix[1]))
     warn ("fitswcs_bin_wcs: CRPIX was improperly computed");
   if (fneqs (dx, wcs_M.cdelt[1]))
     warn ("fitswcs_bin_wcs: CDELT was improperly computed");

   % Now rebin the image to NxN
   % The image currently runs from xmin to M*dx.  We want it to go
   % from xmin to N*(M/N*dx)
   variable N = 93;
   dx = (M*dx)/N;
   dy = (M*dy)/N;
   xgrid = xmin + dx*[0:N-1];
   ygrid = ymin + dy*[0:N-1];
   
   variable wcs_N = fitswcs_rebin_wcs (wcs_M, [M,M], [N,N]);
   variable iy = (crpix[0]-ymin)/dy + 0.5;
   if (fneqs (iy,wcs_N.crpix[0]))
     warn ("fitswcs_rebin_wcs: CRPIX was improperly computed");
   if (fneqs (dy, wcs_N.cdelt[0]))
     warn ("fitswcs_rebin_wcs: CDELT was improperly computed");
}

test_wcs ();

if (Failed == 0)
  message ("Passed");
else
  message ("Failed");
