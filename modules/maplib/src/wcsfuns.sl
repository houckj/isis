% These routines perform FITS WCS transformations using maplib.
% Public routines:
%
%    (x',y') = wcsfuns_project (wcs, x, y);           % From WCS to image
%    (x',y') = wcsfuns_deproject (wcs, x, y);         % From image to WCS
%    (x',y') = wcsfuns_reproject (wcs', wcs, x, y);   % From image to image
%
% The wcs structure is assumed to be the same as given by the fitswcs.
%
% Author: John E. Davis <davis@space.mit.edu>
%---------------------------------------------------------------------------

$0 = 0;	 % Major
$1 = 2;  % Minor
$2 = 1;  % Patch-level

public variable _wcsfuns_version = $0*10000 + $1*100 + $2;
public variable _wcsfuns_version_string = sprintf ("%d.%d.%d", $0, $1, $2);

require ("maplib");
require ("fitswcs");
% Note: the fitswcs structure may not include the equinox and radsys
% fields.  For this reason, check the version number of the fitswcs
% file.

private define parse_ctype (ctype)
{
   ctype = strlow (ctype);

   variable sys, proj;
   if (strlen (ctype) < 6)
     return NULL, "linear";

   if (ctype[4] != '-')
     return NULL, "linear";

   sys = strtrim_end (substr (ctype, 1, 4), "-");
   proj = strtrim_end (substr (ctype, 6, -1), "- ");

   return sys, proj;
}

private define get_lon_lat_indices (sys_0, sys_1)
{
   variable lons = ["ra", "glon", "elon", "slon"];
   variable lats = ["dec", "glat", "elat", "slat"];
   % These correspond to
   %   RA/DEC:         Equatorial
   %   ELON/ELAT:      Ecliptic
   %   GLON/GLAT:      Galactic
   %   SLON/SLAT:      Supergalactic
   _for (0, length (lons)-1, 1)
     {
	variable i = ();
	variable lon_i = lons[i], lat_i = lats[i];
	if ((sys_0 == lon_i) and (sys_1 == lat_i))
	  return 0, 1;
	if ((sys_0 == lat_i) and (sys_1 == lon_i))
	  return 1, 0;
     }
   vmessage ("*** Warning: Assuming %s and %s specify lon and lat coordinates, resp",
	     strup(sys_0), strup(sys_1));
   return 0, 1;
}

% The FITSWCS-II paper defines the following conventions:
%
% The native system is the system where the coordinate mappings are
% defined.  For example, the FITS azimuthal projections use a native
% system in which the (native) pole is located at the center of the
% projection.  The native coordinates are denoted by (phi,theta) and
% the celestial ones by (alpha,delta).  The location of the celestial
% pole is give by (phi_p,theta_p) in native coordinates.  LATPOLEa and
% LONPOLEa keywords encode these values.  The location of the native
% pole is given by (alpha_p, delta_p) in celestial coordinates.
%
% A fiducial point is defined with coordinates (phi_0,theta_0) and
% (alpha_0,delta_0) in the native and celestial systems, resp.  The
% pair (alpha_0,delta_0) is given by CRVALia.  And (PVi_1a,PVi_2a)
% give the corresponding native coordinates.
%
% For azimuthal projections, (alpha_0, theta_0) correspond to maplib's
% lon0 and lat0 parameters.
%
% In FITS parlance, (lon0,lat0) specifies the celestial coordinate of the
% native pole, assumed to be located at the center of the map (given by crpix).
% The FITSWCS-II paper denotes the native coordinate by (phi,theta), and the
% celestial coordinate by (alpha,delta).
%
% The CRVAL keywords give the fiducial celestial coordinates
% (alpha0,delta0). The associated native coordinate is given by
% (phi0,theta0).  The azimuthal projections use (phi0,theta0)=(0,90).
% This means that for azimuthal projections, CRVAL provides the
% celestial coordinates of the native pole.
%
% The LONPOLEa keyword gives the native longitude of the celestial
% pole.  It corresponds to beta in the "sphere" projection.  Its
% default value is 0 for delta0>=theta0, and -180 otherwise.  Recall
% that we can always rotate about the native pole.  Specifying the
% longitude of the celestial pole fixes the rotation angle.  For the
% azimuthal projections, maplib automatically chooses the correct
% implicit value for beta by the constraint that (xhat,yhat) have the
% same directions as (lonhat,lathat) at the fiducial point.
%
% FIXME!!!
% Non-default values of LONPOLEa will require a rotation of the
% tangent plane.
%
%
%
private define init_simple_azimuthal_proj (proj, wcs)
{
   variable m = maplib_new (proj);
   m.lon0 = wcs.crval[0];
   m.lat0 = wcs.crval[1];
   m.x0 = wcs.crpix[0];
   m.y0 = wcs.crpix[1];
   m.xscale = wcs.cdelt[0];
   m.yscale = wcs.cdelt[1];

   return m;
}

% Fits-Name  Maplib-Name       Description        Keywords
%   AIT        hammer           Hammer-Aitoff
%   STG        stereo           Stereographic
%   TAN        gnomic           Gnomic
%   ZEA        lambert          Lambert EA
%   BON        bonne            Bonne EA             PVi_1a
%   SFL        sinusoidal       Sanson-Flamsteed
%   MER        mercator         Mercator
%   CAR        linear           Plate-Carree
%   AZP        plane*           Zenithal             PVi_1a, PVi_2a
%   SZP        plane*           Slant Zenithal       PVi_1a, PVi_2a, PVi_3a
%   SIN        plane*           Slant Orthographic   PVi_1a, PVi_2a
%   ARC        azeqdist         Zenith EqDist
%   ZPN                         Zenith Polynomial    PVi_0a,...,PVi_20a
%   AIR                         Airy                 PVi_1a
%   CYP                         Cylinder Perspective PVi_1a
%   CEA                         Cylinder EqArea      PVi_1a
%   PAR                         Parabolic
%   MOL                         Mollweide
%   COP                         Conic (Colles)
%   COE                         Conic EqArea
%   COD                         Conic EqDist
%   COO                         Conic (Lambert)
%   PCO                         Polyconic
%   TSC                         Tang Sphr Cube
%   CSC                         COBE Quad-Sphr Cube
%   QSC                         Quad-Sphr Cube
% *Note: These are special cases of the "plane" projection.
% Note: PVi_0a, PVi_1a, and PVi_2a specify (theta0,phi0)
%       PVi_3a = LONPOLEa, PVi_4a = LATPOLEa
%
private define init_tan_proj (proj, wcs)
{
   return init_simple_azimuthal_proj ("gnomic", wcs);
}
private define init_sin_proj (proj, wcs)
{
   return init_simple_azimuthal_proj ("ortho", wcs);
}

private define init_stg_proj (proj, wcs)
{
   return init_simple_azimuthal_proj ("stereo", wcs);
}

private define inverse_2x2 (a)
{
   variable a_00 = a[0,0], a_01 = a[0,1], a_10 = a[1,0], a_11 = a[1,1];
   variable det = (a_00*a_11 - a_01*a_10);
   if (det == 0.0)
     verror ("FITS wcs matrix has no inverse");
   a = Double_Type[2,2];
   a[0,0] = a_11;
   a[0,1] = -a_01;
   a[1,0] = -a_10;
   a[1,1] = a_00;
   return a/det;
}

private define init_linear_proj (proj, wcs)
{
   variable m = maplib_new ("linear");

   % In FITS parlance, projections take place from image planes to
   % the WCS system, which is the reverse of how maplib defines the
   % projections.  Hence, use the inverse matrix here, and reverse x0<->x1.
   m.x0 = wcs.crval[0];
   m.y0 = wcs.crval[1];
   m.x1 = wcs.crpix[0];
   m.y1 = wcs.crpix[1];

   if (wcs.pc != NULL)
     m.A = inverse_2x2 (wcs.pc);

   variable A = m.A;
   variable cdelt = wcs.cdelt;
   variable factor = 1.0/cdelt[0];
   A[0,0] *= factor;
   A[0,1] *= factor;
   factor = 1.0/cdelt[1];
   A[1,0] *= factor;
   A[1,1] *= factor;
   return m;
}

private define init_unsupported_proj (proj, wcs)
{
   verror ("FITS WCS map %s is not supported\n", proj);
}

static variable Fits_To_Maplib_Funs = Assoc_Type[Ref_Type, &init_unsupported_proj];
Fits_To_Maplib_Funs["tan"] = &init_tan_proj;
Fits_To_Maplib_Funs["sin"] = &init_sin_proj;
Fits_To_Maplib_Funs["stg"] = &init_stg_proj;
Fits_To_Maplib_Funs["linear"] = &init_linear_proj;
Fits_To_Maplib_Funs["car"] = &init_linear_proj;

% This routine converts a FITS WCS specification to a maplib specification.
% FITS WCS defines a transformation from "image" coordinates X_j to, e.g.,
% celestial coordinates, W_i via:
%
%     Q_i = \sum_j D_i PC_ij (X_j - X0_j)
%     W_i = F_i(Q_1, Q_2)
%
% where, generally speaking F_i are non-linear functions.  The precise meaning
% of the FITS crval keywords depend upon the function F_i.  Generally it appears
% that CRVAL_1 = F_1(0,0) and CRVAL_2 = F_2(0,0).  In other words, the CRVALs
% may be regarded as the WCS value of the reference pixel.  However, the FITS
% WCS documents do not force this interpretation.  For this reason, the last
% step will be written as
%
%     W_i = F_i(Q_1, Q_2; W0_1, W0_2)
%
% where the CRVALs have been represented by W0_i.
%
% The maplib projections from the tangent plane to the sphere are of the form
%
%    U_i = D_i (X_i - X0_i)
%    W_i = M_i(U_1, U_2, W0_1, W0_2)
%
% The maplib linear transformation is:
%
%    W_i = W0_i + \sum_j R_ij (X_j-X0_j)
%
% Hence, the FITS WCS transformation will have to be written as two maplib
% transformations: A linear transformation followed by a maplib projection.
% The linear transformation is required only when PC_ij is not diagonal.
%
% Case 1:
%
%   1.  F_i correspond to a linear transformation (ctypes are linear).
%
%      Want: W_i = W0_i + D_i \sum_j PC_ij (X_j - X0_j)
%      Use maplib linear:
%
%           W_i = W0_i + \sum_j R_ij (X_j - X0_j)
%           R_ij = D_i PC_ij
%
%   2. F_i correspond to a celestial mapping.
%   2a. PC_ij is present
%
%      Want:
%           Q_i = \sum_j D_i PC_ij (X_j - X0_j)
%           W_i = F_i(Q_1, Q_2; W0_1, W0_2)
%
%      Use: Maplib linear followed by Maplib projection
%
%       Linear:      U_i = 0 + \sum_j R_ij (X_j - X0_j) ; R_ij = D_i PC_ij
%       Projection:  W_i = M_i(U_1, U_2, W0_1, W0_2)
%
%   2b. PC_ij is not present
%
%       Want: Q_i = D_i (X_i - X0_i)
%             W_i = F_i(Q_1, Q_2, W0_1, W0_2)
%
%       Use single maplib projection.
%
% Finally note that the matrix for linear transformations should be specified
% as the inverse of PC_ij.  This is because FITS defines "projection" as being
% from the tangent plane to the WCS, whereas maplib treats it the other way.
define wcsfuns_fitswcs_to_maplib ()
{
   variable wcs, flipped_ref;

   if (_NARGS != 2)
     usage ("(m1, m2)=%s(wcs_struct, flipped_ref); %% if m2 is non-NULL apply it after m1 for maplib_project", _function_name());

   (wcs, flipped_ref) = ();
   @flipped_ref = 0;

   if (wcs.naxis != 2)
     verror ("%s: This interface currently supports only 2d maps", _function_name());

   variable m1, m2;

   % First use the FITS ctype to determine the projection.
   % The convention regarding the fits ctype is as follows:
   %
   %   If the ctypes are not in the so-called 4-3 format, then the projection
   %   is linear.  Otherwise, the name contains the coordinate system and
   %   projection
   variable sys_x, proj_x, sys_y, proj_y;
   variable ctype_x = wcs.ctype[0];
   variable ctype_y = wcs.ctype[1];

   (sys_x, proj_x) = parse_ctype (ctype_x);
   (sys_y, proj_y) = parse_ctype (ctype_y);

   if (proj_x != proj_y)
     verror ("%s: ctypes %s and %s are not compatible", _function_name(), ctype_x, ctype_y);

   if (proj_x != "linear")
     {
	variable lon_index = 0, lat_index = 1;
	(lon_index, lat_index) = get_lon_lat_indices (sys_x, sys_y);
	if (lon_index == 1)
	  {
	     % Flip the wcs
	     wcs = fitswcs_slice (wcs, [lon_index, lat_index]);
	     (sys_x, sys_y) = (sys_y, sys_x);
	     (proj_x, proj_y) = (proj_y, proj_x);
	     @flipped_ref = 1;
	  }
     }

   m2 = NULL;
   if ((proj_x != "linear")
       and (wcs.pc != NULL))
     {
	m2 = maplib_new ("linear");
	m2.x1 = wcs.crpix[0];
	m2.y1 = wcs.crpix[1];
	m2.x0 = 0.0;
	m2.y0 = 0.0;
	m2.A = inverse_2x2 (wcs.pc);
	%m2.A = @wcs.pc;

	% Adjust the WCS to conform to m2.
	wcs = @wcs;
	wcs.crpix = @wcs.crpix;
	wcs.crpix[0] = 0;
	wcs.crpix[1] = 0;
	wcs.pc = NULL;
     }

   m1 = (@Fits_To_Maplib_Funs[proj_x])(proj_x, wcs);

   return m1, m2;
}

private define maplib_xxxproject_helper (calling_func_name, nargs)
{
   variable wcs, x, y;

   if (nargs != 3)
     usage ("(x1,y1)=%s(wcs_struct,x,y)", calling_func_name);

   (wcs, x, y) = ();

   if (wcs.naxis != 2)
     verror ("%s: This interface currently supports only 2d maps", calling_func_name);

   variable m1, m2, flipped;
   (m1, m2) = wcsfuns_fitswcs_to_maplib (wcs, &flipped);
   return m1, m2, x, y, flipped;
}

define wcsfuns_project ()
{
   variable m1, m2, x, y, flipped;

   (m1, m2, x, y, flipped) = maplib_xxxproject_helper (_function_name(), _NARGS);

   (x,y) = maplib_project (m1, __tmp(x), __tmp(y));
   if (m2 != NULL)
     (x,y) = maplib_project (m2, __tmp(x), __tmp(y));

   if (flipped)
     return y, x;

   return x, y;
}

define wcsfuns_deproject ()
{
   variable m1, m2, x, y, flipped;
   (m1, m2, x, y, flipped) = maplib_xxxproject_helper (_function_name(), _NARGS);

   if (flipped)
     (x,y) = (y,x);

   if (m2 != NULL)
     (x,y) = maplib_deproject (m2, __tmp(x), __tmp(y));

   (x,y) = maplib_deproject (m1, __tmp(x), __tmp(y));
   return x, y;
}

define wcsfuns_reproject ()
{
   variable wcs_from, wcs_to, x, y;

   if (_NARGS != 4)
     usage ("(x1,y1)=%s(wcs_to, wcs_from, x, y)", _function_name());

   (wcs_to, wcs_from, x, y) = ();

   variable m1_from, m2_from, m1_to, m2_to, flipped_from, flipped_to;

   (m1_from, m2_from) = wcsfuns_fitswcs_to_maplib (wcs_from, &flipped_from);
   (m1_to, m2_to) = wcsfuns_fitswcs_to_maplib (wcs_to, &flipped_to);

   if (flipped_from) (x,y) = (y,x);

   if (m2_from != NULL)
     (x,y) = maplib_deproject (m2_from, __tmp(x), __tmp(y));

   (x, y) = maplib_reproject (m1_to, m1_from, __tmp(x), __tmp(y));

   if (m2_to != NULL)
     (x, y) = maplib_project (m2_to, __tmp(x), __tmp(y));

   if (flipped_to) return y,x;
   return x, y;
}

provide ("wcsfuns");
