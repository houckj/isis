require ("fits");

% FITS WCS may be attached to images or pixel-lists.  The images may be
% found in a standard image HDU, or as a cell of a binary table.
% Pixel lists are encoded as columns of a FITS binary table.  Hence, there
% are three forms of keywords that describe the wcs.

static variable CTYPE_INDX	= 0;
static variable CUNIT_INDX	= 1;
static variable CRVAL_INDX	= 2;
static variable CDELT_INDX	= 3;
static variable CRPIX_INDX	= 4;
static variable PC_INDX		= 5;
static variable CD_INDX		= 6;
static variable PV_INDX		= 7;
static variable PS_INDX		= 8;
static variable Image_Formats = 
  ["CTYPE%d%c", "CUNIT%d%c", "CRVAL%d%c", "CDELT%d%c", "CRPIX%d%c", 
   "PC%d_%d%c", "CD%d_%d%c", "PV%d_%d%c", "PS%d_%d%c"];
static variable Column_Formats =
  ["TCTYP%d%c", "TCUNI%d%c", "TCRVL%d%c", "TCDLT%d%c", "TCRPX%d%c", 
   "TP%d_%d%c", "TC%d_%d%c", "TV%d_%d%c", "TS%d_%d%c"];


% This structure defines a transformation of the following form:
% 
%    B_i = CDELT_i PC_ij (A_j - CRPIX_j)        (no sum on i)
%    {X} = PROJ_FUNC({CTYPE}, {PV}, {PS}, {B})
% 
% where {X} signifies the set of values.  Using this notation, the first 
% equation could have been written
% 
%    {B}_i = {CDELT}_i {PC}_ij ({A}_j - {CRPIX}_j)  (no sum on i)
%
static variable WCS_Type = struct
{
   naxis,			       %  number of axis to transform
   ctype,			       %  String_Type[naxis]
   cunit,			       %  String_Type[naxis]
   crval,			       %  Double_Type[naxis]
   crpix,			       %  Double_Type[naxis]
   cdelt,			       %  Double_Type[naxis]
   pc,				       %  Double_Type[naxis,naxis] or NULL
   pv,				       %  array of parameters 
   ps,				       %  array of string parms
   wcsname,			       %  String_Type[naxis]
};

static define read_simple_wcs_keywords ()
{
   variable indices = __pop_args (_NARGS-2);
   variable formats = ();
   variable fp = ();
   
   variable ctype, cunit, crval, crpix, cdelt;
   
   (ctype, cunit, crval, crpix, cdelt) = 
     fits_read_key (fp, sprintf (formats[CTYPE_INDX], __push_args(indices)),
		    sprintf (formats[CUNIT_INDX], __push_args(indices)),
		    sprintf (formats[CRVAL_INDX], __push_args(indices)),
		    sprintf (formats[CRPIX_INDX], __push_args(indices)),
		    sprintf (formats[CDELT_INDX], __push_args(indices)));

   if (ctype == NULL)
     ctype = "linear";
   if (crval == NULL)
     crval = 0.0;
   if (crpix == NULL)
     crpix = 0.0;
   if (cdelt == NULL)
     cdelt = 1.0;

   return (ctype, cunit, crval, crpix, cdelt);
}

static define make_diag_matrix (n, diag)
{
   variable d = Double_Type[n, n];
   d [[0:n*n-1:n]] = diag;
   return d;
}

static define read_matrix_wcs_keywords (fp, pc_fmt, is, js, a, diag)
{
   variable pc = NULL;
   variable n = length (is);
   foreach (is)
     {
	variable i = ();
	variable i1 = i-1;
	foreach (js)
	  {
	     variable j = ();
	     variable j1 = j-1;
	     variable pc_ij = sprintf (pc_fmt, i, j, a);
	     pc_ij = fits_read_key (fp, pc_ij);
	     if ((pc_ij == NULL) and (pc == NULL))
	       continue;
	     if (pc == NULL)
	       pc = make_diag_matrix (n, diag);
	     if (pc_ij != NULL)
	       pc [i1,j1] = pc_ij;
	  }
     }
   return pc;
}

static define get_wcs_naxis (fp)
{
   variable naxis = fits_read_key (fp, "WCSAXES");
   if (naxis == NULL)
     naxis = fits_read_key (fp, "NAXIS");

   return naxis;
}

static define open_interesting_hdu (file, type)
{
   variable fp = fits_open_file (file, "r");
   fits_move_to_interesting_hdu (fp, type);
   return fp;
}

static define allocate_wcs (naxis)
{
   variable wcs = @WCS_Type;
   wcs.naxis = naxis;
   wcs.ctype = String_Type[naxis];
   wcs.cunit = String_Type[naxis];
   wcs.crval = Double_Type[naxis];
   wcs.crpix = Double_Type[naxis];
   wcs.cdelt = Double_Type[naxis];
   wcs.pc = NULL;
   wcs.pv = NULL;
   wcs.ps = NULL;
   return wcs;
}

% This function will be used later when applying the wcs
static define simplify_wcs (wcs)
{
   variable pc = wcs.pc;
   variable n = wcs.naxis;

   if (pc != NULL)
     {
	% If pc is diagonal, then factor diagonal elements into the cdelts
	variable d = @pc;
	variable i = [0:n*n-1:n];      %  diagonal elements
	d[i] = 0.0;
	if (0 == length (where (d != pc)))
	  {
	     wcs.cdelts *= pc[i];
	     wcs.pc = NULL;
	  }
     }
   return wcs;
}

static define check_for_crota_hack (fp, wcs)
{
   if (wcs.naxis < 2)
     return NULL;

   variable rot = fits_read_key (fp, "CROTA2");
   if ((rot == NULL) or (rot == 0.0))
     {
	rot = fits_read_key (fp, "CROTA1");
	if ((rot == NULL) or (rot == 0.0))
	  return NULL;
     }
   rot *= PI/180.0;

   % Apparantly, this angle gets applied AFTER cdelts applied.  Hence, we have
   %
   %     B_i = R(rot)_ij CDELT_j (A_j - CRPIX_j)     (no sum on i)
   % 
   % instead of 
   %
   %    B_i = CDELT_i PC_ij (A_j - CRPIX_j)   (no sum on i)
   %
   % Hence, 
   %
   %   PC_ij = R(rot)_ij CDELT_j/CDELT_i  (no sum)
   variable pc = make_diag_matrix (wcs.naxis);
   variable cdelt0 = wcs.cdelts[0];
   variable cdelt1 = wcs.cdelts[0];
   variable c = cos(rot), s = sin(rot);
   pc[0,0] = c;               pc[0,1] = -s*cdelt1/cdelt0;
   pc[1,0] = s*cdelt0/cdelt1; pc[1,1] = c;

   return pc;
}

%!%+
%\function{fitswcs_get_img_wcs}
%\synopsis{Read a WCS for a FITS image}
%\usage{wcs = fitswcs_get_img_wcs (file [,alt_axis])}
%\description
%\example
%\notes
%\seealso{}
%!%-
public define fitswcs_get_img_wcs ()
{
   variable fp, a = 0;

   if (_NARGS == 1)
     fp = ();
   else if (_NARGS == 2)
     (fp, a) = ();
   else
     usage ("wcs = %s(file [,alt_axis_char])", _function_name);

   if (typeof (fp) == String_Type)
     fp = open_interesting_hdu (fp, _FITS_IMAGE_HDU);
   
   variable naxis = get_wcs_naxis (fp);
   variable wcs = allocate_wcs (naxis);
   variable ctype = wcs.ctype, cunit = wcs.cunit, crval = wcs.crval, 
     crpix = wcs.crpix, cdelt = wcs.cdelt;

   _for (0, naxis-1 , 1)
     {
	variable i = ();

	(ctype[i], cunit[i], crval[i], crpix[i], cdelt[i]) 
	  = read_simple_wcs_keywords (fp, Image_Formats, i+1, a);
     }

   variable pc;
   pc = read_matrix_wcs_keywords (fp, "PC%d_%d", [1:naxis], [1:naxis], a, 1.0);
   if (pc == NULL)
     pc = read_matrix_wcs_keywords (fp, "CD%d_%d", [1:naxis], [1:naxis], a, 0.0);

   if ((pc == NULL) and (a == 0))
     pc = check_for_crota_hack (fp, wcs);

   wcs.pc = pc;
   return wcs;
}

public define fitswcs_get_column_wcs ()
{
   variable fp, a = 0, column_names = NULL;

   if (_NARGS == 3)
     (fp, column_names, a) = ();
   else if (_NARGS == 2)
     (fp, column_names) = ();
   
   if ((column_names == NULL) or (typeof(a) == String_Type))
     usage ("wcs = %s(file, array_of_column_names [,alt_axis_char]);\n%s", 
	    _function_name,
	    sprintf ("Example: wcs = %s(\"evt1.fits\", [\"X\",\"Y\"],'A');\n",
		     _function_name));

   if (typeof (fp) == String_Type)
     fp = open_interesting_hdu (fp, _FITS_BINARY_TBL);
   
   variable naxis = length (column_names);

   variable wcs = allocate_wcs (naxis);
   variable ctype = wcs.ctype, cunit = wcs.cunit, crval = wcs.crval, 
     crpix = wcs.crpix, cdelt = wcs.cdelt;

   variable column_nums = Int_Type[naxis];
   _for (0, naxis-1 , 1)
     {
	variable i = ();
	variable col = fits_get_colnum (fp, column_names[i]);
	column_nums[i] = col;

	(ctype[i], cunit[i], crval[i], crpix[i], cdelt[i]) 
	  = read_simple_wcs_keywords (fp, Column_Formats, col, a);
     }

   variable pc;
   pc = read_matrix_wcs_keywords (fp, Column_Formats[PC_INDX], column_nums, column_nums, a, 1.0);
   if (pc == NULL)
     pc = read_matrix_wcs_keywords (fp, Column_Formats[CD_INDX], column_nums, column_nums, a, 0.0);

   wcs.pc = pc;
   return wcs;
}

provide("fitswcs");
