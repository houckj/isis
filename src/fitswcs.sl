%    Copyright (C) 2004-2005 Massachusetts Institute of Technology
%
%    Author:  John E. Davis <davis@space.mit.edu>
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

% TODO:
% 
%   * Add support for cross-references when dealing with vector columns.
% 
%   * Add support for writing binary table WCS info
%
% Notes:
% 
%   Some fits files do not use the WCSAXES keyword to specifiy the dimensionality
%   of the WCS.  Instead a common practice is to add so-called degenerate axes
%   where one increases the NAXIS keyword to the WCS dimensionality and adds
%   keywords NAXISj=1 for the new dimensions.  To deal with this hack, the 
%   user will be required to reshape the image to its proper dimensionality.
%   Then the projection routines defined here may be used on any subspace of 
%   the WCS by extracting that subspace via the fitswcs_slice routine.

require ("fits");
variable _fitswcs_version = 100;
variable _fitswcs_version_string = "0.1.1-0";

% FITS WCS may be attached to images or pixel-lists.  The images may be
% found in a standard image HDU, or as a cell of a binary table.
% Pixel lists are encoded as columns of a FITS binary table.  Hence, there
% are three forms of keywords that describe the wcs.

private variable CTYPE_INDX	= 0;
private variable CUNIT_INDX	= 1;
private variable CRVAL_INDX	= 2;
private variable CDELT_INDX	= 3;
private variable CRPIX_INDX	= 4;
private variable PC_INDX	= 5;
private variable CD_INDX	= 6;
private variable PV_INDX	= 7;
private variable PS_INDX	= 8;
private variable Image_Formats = 
  ["CTYPE%d", "CUNIT%d", "CRVAL%d", "CDELT%d", "CRPIX%d", 
   "PC%d_%d", "CD%d_%d", "PV%d_%d", "PS%d_%d"];
private variable Image_Formats_Alt = Image_Formats + "%c";

private variable Column_Formats =
  ["TCTYP%d", "TCUNI%d", "TCRVL%d", "TCDLT%d", "TCRPX%d", 
   "TP%d_%d", "TC%d_%d", "TV%d_%d", "TS%d_%d"];
private variable Column_Formats_Alt =
  ["TCTY%d%c", "TCUN%d%c", "TCRV%d%c", "TCDE%d%c", "TCRP%d%c", 
   "TP%d_%d%c", "TC%d_%d%c", "TV%d_%d%c", "TS%d_%d%c"];

private variable Vector_Formats = 
  ["%dCTYP%d", "%dCUNI%d", "%dCRVL%d", "%dCDLT%d", "%dCRPX%d", 
   "%dP%d_%d", "%dC%d_%d", "%dV%d_%d", "%dS%d_%d"];
private variable Vector_Formats_Alt =
  ["%dCTY%d%c", "%dCUN%d%c", "%dCRV%d%c", "%dCDE%d%c", "%dCRP%d%c", 
   "%dP%d_%d%c", "%dC%d_%d%c", "%dV%d_%d%c", "%dS%d_%d%c"];

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
private variable WCS_Type = struct
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
   wcsname,			       %  String_Type
   radsys,
   equinox
};

% Notes: For usage of the and default values for the radsys and
% equinox, see, e.g.,
% <http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?imcctran>

%!%+
%\function{fitswcs_new}
%\synopsis{Create a new-ndimensional linear WCS}
%\usage{wcs = fitswcs_new (Int_Type naxis)}
%\description
%  This function returns a new WCS structure of the specified dimensionality
%  that represents an identity (linear) transformation.
%\seealso{fitswcs_get_img_wcs, fitswcs_get_column_wcs, fitswcs_get_vector_wcs}
%!%-
define fitswcs_new (naxis)
{
   variable wcs = @WCS_Type;
   wcs.naxis = naxis;
   wcs.ctype = String_Type[naxis];  wcs.ctype[*] = "linear";
   wcs.cunit = String_Type[naxis];
   wcs.crval = Double_Type[naxis];  wcs.crval[*] = 0.0;
   wcs.crpix = Double_Type[naxis];  wcs.crpix[*] = 0.0;
   wcs.cdelt = Double_Type[naxis];  wcs.crpix[*] = 1.0;
   wcs.pc = NULL;
   wcs.pv = NULL;
   wcs.ps = NULL;
   return wcs;
}


%!%+
%\function{fitswcs_slice}
%\synopsis{Form a new wcs from one or more axes of another}
%\usage{new_wcs = fitswcs_slice (wcs, dims)}
%\description
%  This function may be used to construct a new wcs from another by rearranging 
%  its axes or by using a subset of them.  The \exmp{dims} argument specifies 
%  the dimensions to use.
%\example
%  Suppose that \exmp{wcs} represents a 4 dimensional WCS. Then
%#v+
%    wcs2 = fitswcs_slice (wcs, [0,1]);
%#v-
%  will result in a 2 dimensional WCS from the first 2 axis of the input WCS.
%  Similarly,
%#v+
%    wcs2 = fitswcs_slice (wcs, [1,0]);
%#v-
%  will produce a 2d WCS with the first two axes swapped.
%\seealso{fitswcs_get_img_wcs, fitswcs_get_column_wcs, fitswcs_get_vector_wcs}
%!%-
define fitswcs_slice ()
{
   if (_NARGS < 2)
     usage ("wcs1 = %s(wcs, dims-array)", _function_name ());
   
   variable wcs, dims;
   (wcs, dims) = ();
   variable new_wcs = @wcs;
   
   new_wcs.naxis = length (dims);
   new_wcs.ctype = wcs.ctype[dims];
   new_wcs.cunit = wcs.cunit[dims];
   new_wcs.crval = wcs.crval[dims];
   new_wcs.crpix = wcs.crpix[dims];
   new_wcs.cdelt = wcs.cdelt[dims];
   if (wcs.pc != NULL)
     new_wcs.pc = wcs.pc[dims, dims];
   return new_wcs;
}


%---------------------------------------------------------------------------
% Some simple utility functions
%---------------------------------------------------------------------------

private define make_diag_matrix (n, diag)
{
   variable d = Double_Type[n, n];
   d [[0:n*n-1:n+1]] = diag;
   return d;
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

private define dup_wcs (wcs)
{
   wcs = @wcs;
   
   foreach (get_struct_field_names (wcs))
     {
	variable field = ();
	variable value = get_struct_field (wcs, field);
	if (typeof (value) == Array_Type)
	  set_struct_field (wcs, field, @value);
     }
   return wcs;
}

% Convert to/from FORTRAN/C order
private define reverse_wcs (wcs)
{
   return fitswcs_slice (wcs, [wcs.naxis-1:0:-1]);
}

private define read_simple_wcs_keywords ()
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

private define read_matrix_wcs_keywords (fp, pc_fmt, is, js, a, diag)
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

private define get_wcs_naxis (fp)
{
   variable naxis = fits_read_key (fp, "WCSAXES");
   if (naxis == NULL)
     naxis = fits_read_key (fp, "NAXIS");

   return naxis;
}

private define open_interesting_hdu (file, type)
{
   variable fp = fits_open_file (file, "r");
   fits_move_to_interesting_hdu (fp, type);
   return fp;
}


private define check_for_crota_hack (fp, wcs)
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
%\usage{wcs = fitswcs_get_img_wcs (fp [,alt])}
%\description
%  The \sfun{fitswcs_get_img_wcs} returns a structure representing a WCS from
%  the specified file descriptor \exmp{fp} corresponding to an image HDU.
%  An optional parameter may be used to specified an alternate WCS.
%\example
%#v+
%   wcs = fitswcs_get_img_wcs ("img.fits[IMAGE]", 'P');
%#v-
%\seealso{fitswcs_put_img_wcs, fitswcs_get_column_wcs, fitswcs_get_vector_wcs}
%!%-
define fitswcs_get_img_wcs ()
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
   variable wcs = fitswcs_new (naxis);
   variable ctype = wcs.ctype, cunit = wcs.cunit, crval = wcs.crval, 
     crpix = wcs.crpix, cdelt = wcs.cdelt;

   wcs.wcsname = fits_read_key (fp, sprintf ("WCSNAME%c", a));
   variable formats = Image_Formats;
   if (a)
     formats = Image_Formats_Alt;
   _for (0, naxis-1 , 1)
     {
	variable i = ();

	(ctype[i], cunit[i], crval[i], crpix[i], cdelt[i]) 
	  = read_simple_wcs_keywords (fp, formats, i+1, a);
     }

   variable pc;
   pc = read_matrix_wcs_keywords (fp, formats[PC_INDX], [1:naxis], [1:naxis], a, 1.0);
   if (pc == NULL)
     pc = read_matrix_wcs_keywords (fp, formats[CD_INDX], [1:naxis], [1:naxis], a, 0.0);

   if ((pc == NULL) and (a == 0))
     pc = check_for_crota_hack (fp, wcs);

   wcs.pc = pc;
   
   return reverse_wcs (wcs);
}


%!%+
%\function{fitswcs_get_column_wcs}
%\synopsis{Get the WCS attached to one or more columns of a binary table}
%\usage{fitswcs_get_column_wcs (fp, columns-array [,alt]}
%\description
%  This function may be used to obtain the WCS associated with one or more
%  columns of a binary table.  The file descriptor \exmp{fp} must specify 
%  a binary table.  The \exmp{columns-array} parameter should be an array
%  of columns names.  The third parameter is optional and is used to specify
%  an alternate WCS.
%\example
%#v+
%   wcs = fitswcs_get_column_wcs ("evt1.fits[EVENTS]", ["X","Y"]);
%#v-
%\seealso{fitswcs_put_column_wcs, fitswcs_get_img_wcs, fitswcs_get_vector_wcs}
%!%-
define fitswcs_get_column_wcs ()
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

   variable wcs = fitswcs_new (naxis);
   variable ctype = wcs.ctype, cunit = wcs.cunit, crval = wcs.crval, 
     crpix = wcs.crpix, cdelt = wcs.cdelt;

   % The wcsname is not supported.  The specification of this is flawed for
   % "pixel-lists".  A wcsname applies to the WCS as a whole-- not just a 
   % column.
   wcs.wcsname = NULL;

   variable formats = Column_Formats;
   if (a)
     formats = Column_Formats_Alt;

   variable column_nums = Int_Type[naxis];
   _for (0, naxis-1 , 1)
     {
	variable i = ();
	variable col = fits_get_colnum (fp, column_names[i]);
	column_nums[i] = col;

	(ctype[i], cunit[i], crval[i], crpix[i], cdelt[i]) 
	  = read_simple_wcs_keywords (fp, formats, col, a);
     }

   variable pc;
   pc = read_matrix_wcs_keywords (fp, formats[PC_INDX], column_nums, column_nums, a, 1.0);
   if (pc == NULL)
     pc = read_matrix_wcs_keywords (fp, formats[CD_INDX], column_nums, column_nums, a, 0.0);

   wcs.pc = pc;
   return wcs;
}

private define read_key_or_col (fp, key, row)
{
   %vmessage ("Looking for %s in row %d", key, row);
   variable val = fits_read_key (fp, key);
   if (val != NULL)
     return val;
   
   % Perhaps it is the name of a column.
   if (fits_binary_table_column_exists (fp, key))
     val = fits_read_cell (fp, key, row);

   return val;
}

private define read_vector_matrix_wcs_keywords (fp, pc_fmt, is, js, col, row, a, diag)
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
	     variable pc_ij = sprintf (pc_fmt, i, j, col, a);
	     pc_ij = read_key_or_col (fp, pc_ij, row);
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


%!%+
%\function{fitswcs_get_vector_wcs}
%\synopsis{Get the WCS of an image in a specified table cell}
%\usage{wcs = fitswcs_get_vector_wcs (fp, column_name, row [,alt])}
%\description
%  This function reads the WCS of an image in a specified cell of a binary
%  table given by \exmp{fp} parameter.  The second and third parameters specify
%  the column name and row number of the cell.  An optional fourth parameter
%  may be used to obtain the corresponding alternate WCS.
%\example
%  This example reads the WCS associated with the image in the second row
%  of the QEU column of the binary table with HDUNAME equal to AXAF_QEU1
%  in the file "HRCQEU.fits":
%#v+
%    wcs = fitswcs_get_vector_wcs ("HRCQEU.fits[AXAF_QEU1], "QEU", 2);
%#v-
%\notes
%  The current implementation does not yet support references to the WCS
%  of other cells.
%\seealso{fitswcs_get_column_wcs, fitswcs_get_img_wcs}
%!%-
define fitswcs_get_vector_wcs ()
{
   variable fp, col, row, a = 0;

   switch (_NARGS)
     {
      case 3:
	(fp, col, row) = ();
     }
     {
      case 4:
	(fp, col, row, a) = ();
     }
     {
	usage ("wcs = %s(file, column, row [,alt])", _function_name);
     }
   
   if (typeof (fp) == String_Type)
     fp = open_interesting_hdu (fp, _FITS_BINARY_TBL);

   if (typeof (col) == String_Type)
     col = fits_get_colnum (fp, col);

   % Get the number of axes for the WCS.  This may be given explicitly by
   % the
   variable naxis = read_key_or_col (fp, sprintf ("WCAX%d%c", col, a), row);
   if (naxis == NULL)
     {
	% Read the image and get its dimensions
	variable img = fits_read_cell (fp, col, row);
	variable dims; (dims,,) = array_info (__tmp(img));
	naxis = length (dims);
     }
   variable wcs = fitswcs_new (naxis);
   variable ctype = wcs.ctype, cunit = wcs.cunit, crval = wcs.crval, 
     crpix = wcs.crpix, cdelt = wcs.cdelt;

   wcs.wcsname = read_key_or_col (fp, sprintf ("WCSN%d%c", col, a), row);

   variable formats = Vector_Formats;
   if (a)
     formats = Vector_Formats_Alt;

   _for (0, naxis-1 , 1)
     {
	variable i = ();
	variable i1 = i+1;
	variable val;
	
	val = read_key_or_col (fp, sprintf (formats[CTYPE_INDX], i1, col, a), row);
	if (val != NULL) ctype[i] = val;

	val = read_key_or_col (fp, sprintf (formats[CUNIT_INDX], i1, col, a), row);
	if (val != NULL) cunit[i] = val;

	val = read_key_or_col (fp, sprintf (formats[CRVAL_INDX], i1, col, a), row);
	if (val != NULL) crval[i] = val;

	val = read_key_or_col (fp, sprintf (formats[CRPIX_INDX], i1, col, a), row);
	if (val != NULL) crpix[i] = val;

	val = read_key_or_col (fp, sprintf (formats[CDELT_INDX], i1, col, a), row);
	if (val != NULL) cdelt[i] = val;
     }

   variable pc;
   dims = [1:naxis];
   pc = read_vector_matrix_wcs_keywords (fp, formats[PC_INDX], dims, dims, col, row, a, 1.0);
   if (pc == NULL)
     pc = read_vector_matrix_wcs_keywords (fp, formats[CD_INDX], dims, dims, col, row, a, 0.0);

   wcs.pc = pc;
   return reverse_wcs (wcs);
}


%!%+
%\function{fitswcs_new_img_wcs}
%\synopsis{Create a linear WCS for an image}
%\usage{wcs = fitswcs_new_img_wcs (grid0,grid1,...)}
%\description
%  This function may be used to construct a linear WCS for an image with the
%  specified grids.  The grids are assumed to be linear.
%\example
%  Use the histogram module's hist2d function to create an image from the X
%  and Y columns in a file, and the construct a corresponding WCS:
%#v+
%    (x,y) = fits_read_col ("table.fits", "X", "Y");
%    gridx = [min(x):max(x):0.5];
%    gridy = [min(y):max(y):0.5];
%    img = hist2d (y,x,gridy,gridx);
%    wcs = fitswcs_new_img_wcs (gridy, gridx);
%#v-
%\seealso{fitswcs_new, fitswcs_get_img_wcs}
%!%-
define fitswcs_new_img_wcs ()
{
   variable naxis = _NARGS;
   if (naxis == 0)
     usage ("wcs = %s(grid_dim0, grid_dim1, ..., grid_dimN)", _function_name());

   variable wcs = fitswcs_new (naxis);
   variable i = naxis;
   loop (naxis)
     {
	i--;
	variable grid = ();
	variable x0 = grid[0];
	wcs.crpix[i] = 0.5;
	wcs.crval[i] = x0;
	wcs.cdelt[i] = grid[1]-x0;
     }
   return wcs;
}

private define write_wcs_keyword (fp, format, index, axis, a, value, i)
{
   if (value == NULL)
     return;
   value = value[i];
   if (value == NULL)
     return;
   variable key = sprintf (format[index], axis, a);
   fits_update_key (fp, key, value, "");
}


%!%+
%\function{fitswcs_put_img_wcs}
%\synopsis{Write a WCS out to an image header}
%\usage{fitswcs_put_img_wcs (fp, wcs [,alt])}
%\description
%  The \sfun{fitswcs_put_img_wcs} may be used to write the specified wcs
%  out to the image HDU specified by the \exmp{fp} parameter.  An optional 
%  third parameter may be used to specify an alternate WCS.
%\example
%#v+
%    fp = fits_open_file ("img.fits", "w");
%      .
%      .
%      .
%    fits_put_img_wcs (fp, wcs, 'P');
%    fits_close_file (fp);
%#v-
%\seealso{fitswcs_put_column_wcs}
%!%-
define fitswcs_put_img_wcs ()
{
   variable fp, a = 0;
   variable wcs;

   if (_NARGS == 2)
     (fp, wcs) = ();
   else if (_NARGS == 3)
     (fp, wcs, a) = ();
   else
     usage ("%s(fp, wcs-struct [,alt_axis_char])", _function_name);

   if (typeof (fp) != Fits_File_Type)
     fp = fits_open_file (fp, "w");

   % Write the wcs out in FORTRAN order
   wcs = reverse_wcs (wcs);

   variable ctype = wcs.ctype, cunit = wcs.cunit, crval = wcs.crval, 
     crpix = wcs.crpix, cdelt = wcs.cdelt, naxis = wcs.naxis;

   fits_update_key (fp, sprintf ("WCSAXES%c", a), naxis);
   if (wcs.wcsname != NULL)
     fits_update_key (fp, sprintf ("WCSNAME%c", a), wcs.wcsname);

   variable i, j;
   variable formats = Image_Formats;
   if (a)
     formats = Image_Formats_Alt;

   variable axes = [1:naxis];
   _for (0, naxis-1 , 1)
     {
	i = ();
	j = axes[i];
	write_wcs_keyword (fp, formats, CTYPE_INDX, j, a, ctype, i);
	write_wcs_keyword (fp, formats, CUNIT_INDX, j, a, cunit, i);
	write_wcs_keyword (fp, formats, CRVAL_INDX, j, a, crval, i);
	write_wcs_keyword (fp, formats, CRPIX_INDX, j, a, crpix, i);
	write_wcs_keyword (fp, formats, CDELT_INDX, j, a, cdelt, i);
     }

   variable pc = wcs.pc;
   if (pc != NULL) _for (0, naxis-1, 1)
     {
	i = ();
	_for (0, naxis-1, 1)
	  {
	     j = ();
	     fits_update_key (fp, sprintf (formats[PC_INDX], axes[i], axes[j], a),
			      pc[i, j], NULL);
	  }
     }

   % FIXME: Write the rest of the wcs structure
}


%!%+
%\function{fitswcs_put_column_wcs}
%\synopsis{Write the WCS attached to one or more table columns}
%\usage{fitswcs_put_column_wcs (fp, wcs, columns-array [,alt])}
%\description
%  This function may be used to attach a WCS to one or more columns of a binary
%  table.  The dimensionality of the specified WCS must match the length of the
%  array specifying the column names.  The first parameter, \exmp{fp} must specify
%  a binary table extension.  The fourth parameter is optional and may be used
%  to specify an alternate WCS.
%\example
%#v+
%   fitswcs_put_column_wcs ("evt2.fits[EVENTS], wcs, ["X","Y"]);
%#v-
%\seealso{fitswcs_get_column_wcs, fitswcs_put_img_wcs, fitswcs_get_img_wcs}
%!%-
define fitswcs_put_column_wcs ()
{
   variable fp, wcs, columns, a = 0;

   switch (_NARGS)
     {
      case 3:
	(fp, wcs, columns) = ();
     }
     {
      case 4:
	(fp, wcs, columns, a) = ();
     }
     {
	usage ("%s(fp, wcs-struct, columns-array [,alt])", _function_name());
     }

   variable ctype = wcs.ctype, cunit = wcs.cunit, crval = wcs.crval, 
     crpix = wcs.crpix, cdelt = wcs.cdelt, naxis = wcs.naxis;

   if (length (columns) != naxis)
     verror ("The dimensionality of the specified WCS does not match the number of columns");
   
   if (typeof (fp) != Fits_File_Type)
     fp = fits_open_file (fp, "w");
   
   variable i, cols = Int_Type[naxis];
   _for (0, naxis-1, 1)
     {
	i = ();
	variable colname = columns[i];
	if (0 == fits_binary_table_column_exists (fp, colname))
	  verror ("Binary table column %s does not exist", colname);
	cols[i] = fits_get_colnum (fp, colname);
     }

   variable formats = Column_Formats;
   if (a)
     formats = Column_Formats_Alt;

   _for (0, naxis-1 , 1)
     {
	i = ();
	variable j = cols[i];
	write_wcs_keyword (fp, formats, CTYPE_INDX, j, a, ctype, i);
	write_wcs_keyword (fp, formats, CUNIT_INDX, j, a, cunit, i);
	write_wcs_keyword (fp, formats, CRVAL_INDX, j, a, crval, i);
	write_wcs_keyword (fp, formats, CRPIX_INDX, j, a, crpix, i);
	write_wcs_keyword (fp, formats, CDELT_INDX, j, a, cdelt, i);
     }

   variable pc = wcs.pc;
   if (pc != NULL) _for (0, naxis-1, 1)
     {
	i = ();
	_for (0, naxis-1, 1)
	  {
	     j = ();
	     fits_update_key (fp, sprintf (formats[PC_INDX], cols[i], cols[j], a),
			      pc[i, j], NULL);
	  }
     }

   % FIXME: Write the rest of the wcs structure
}



% This function will be used later when applying the wcs
private define simplify_wcs (wcs)
{
   variable pc = wcs.pc;
   variable n = wcs.naxis;

   if (pc != NULL)
     {
	% If pc is diagonal, then factor diagonal elements into the cdelts
	variable d = @pc;
	variable i = [0:n*n-1:n+1]*1;  %  *1 to force it from being a range
	d[i] = 0.0;
	if (0 == length (where (d != 0)))
	  {
	     wcs.cdelt *= pc[i];
	     wcs.pc = NULL;
	  }
#iffalse
	else
	  {
	     print (d);
	     print (pc);
	  }
#endif
     }
   return wcs;
}


%!%+
%\function{fitswcs_linear_transform_wcs}
%\synopsis{Apply a linear transformation to a WCS}
%\usage{wcs1 = fitswcs_linear_transform_wcs (wcs, X0, CD, I0)}
%#v+
%     wcs: The specified WCS to transform
%   X0,I0: 1-d arrays
%      CD: 2-d array
%#v-
%\description
%  This function may be used to create a new WCS by applying a linear 
%  transformation to an existing one.
%\notes
%  The dimensionality of the WCS is limited to 2 in the 
%  current implementation.
%\seealso{fitswcs_rebin_wcs}
%!%-

% The FITS linear transformation may be written:
% 
%     W_i = W0_i + \sum_j D_i PC_ij (X_j - X0_j)
%     
%  written as:  LX = W = W0 + CD#(X-X0)
%  
% Consider now another linear transformation:
% 
%    X = L'I = X0' + CD' # (I-I0)
%
% Then LX = LL'I
%
%      (LL')I = W0 + CD#(X0' + CD'#(I-I0) - X0)
%             = W0 + (CD#CD')#[I-I0 - invCD'#(X0-X0')]
%             = W0 + (CD#CD')#[I-I0']
%
% where I0' = I0 + invCD'#(X0-X0')
%
define fitswcs_linear_transform_wcs ()
{
   if (_NARGS != 4)
     usage ("new_wcs = %s(old_wcs, X0, CD, I0);\n%s\n%s",
	    _function_name,
	    "where X = X0 + CD#(I-I0) is the relationship between the new system I",
	    "and the old system X");
   
   variable wcs, x0, cd, i0;
   (wcs, x0, cd, i0) = ();

   variable wcs_naxis = wcs.naxis;
   if (wcs_naxis != 2)
     verror ("%s: Currently this routine supports only 2d wcs systems", _function_name);

   variable dims, num_dims;
   (dims, num_dims, ) = array_info (cd);
   variable naxis;
   
   if (num_dims == 1)
     {
	naxis = dims[0];
	variable tmp = Double_Type[naxis, naxis];
	tmp[[0:naxis*naxis-1:naxis+1]] = cd;
	cd = tmp;
     }
   else if (num_dims == 2)
     {
	naxis = dims[0];
	if (naxis != dims[1])
	  verror ("%s: CD matrix must be square", _function_name);
     }
   else verror ("%s: CD matrix must be square", _function_name);

   if (wcs_naxis != naxis)
     verror ("%s: wcs requires CD to be [%d,%d]", _function_name, wcs_naxis, wcs_naxis);
	
   wcs = dup_wcs (wcs);
   
   % Compute this: I0' = I0 + invCD'#(X0-X0')
   wcs.crpix = i0 + inverse_2x2(cd)#(wcs.crpix-x0);

   variable pc = wcs.pc;
   if (pc != NULL)
     pc = pc#cd;
   else
     pc = cd;
   
   wcs.pc = pc;
   return simplify_wcs (wcs);
}


%!%+
%\function{fitswcs_rebin_wcs}
%\synopsis{This function may be used to obtain the wcs for a rebinned image}
%\usage{wcs1 = fitswcs_rebin_wcs (wcs, grid0, grid1, ...)}
%\description
%  This function may be used to construct the WCS for a rebinned image from
%  the WCS of of the unbinned image.  The grid parameters specify the linear
%  grids the new image.
%\seealso{fitswcs_linear_transform_wcs, fitswcs_slice}
%!%-
define fitswcs_rebin_wcs ()
{
   if (_NARGS < 2)
     usage ("wcs = %s (wcs, grid_dim0, grid_dim1,...)", _function_name());

   variable naxis = _NARGS-1;
   variable X0 = Double_Type[naxis];
   variable dX = Double_Type[naxis];
   variable I0 = 0.5 + Double_Type[naxis];
   variable i = naxis;
   loop (naxis)
     {
	i--;
	variable grid = ();
	X0[i] = grid[0];
	dX[i] = grid[1]-grid[0];
     }
   variable wcs = ();
   return fitswcs_linear_transform_wcs (wcs, X0, dX, I0);
}

define fitswcs_translate_wcs ()
{
   if (_NARGS < 2)
     usage ("wcs = %s (wcs, dX_array)", _function_name());
   verror ("Not Implemented");
}

provide("fitswcs");
