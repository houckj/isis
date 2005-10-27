%    Copyright (C) 1998-2005 Massachusetts Institute of Technology
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

%if (current_namespace () != "")
%  import ("cfitsio", current_namespace ());
%else
import ("cfitsio");

variable _fits_sl_version = 402;
variable _fits_sl_version_string = "0.4.2-0";

% Forward declarations
private define do_fits_error ()
{
   variable status, file = "";
   
   if (_NARGS == 2)
     file = ();
   status = ();
   if (status == 0)
     return;
   
   if (strlen(file))
     file = strcat (file, ": ");
   verror ("%s%s", file, _fits_get_errstatus (status));
}

%!%+
%\function{fits_open_file}
%\synopsis{Open a fits file}
%\usage{Fits_File_Type fits_open_file (String_Type filename, String_Type mode)}
%\description
%  The \var{fits_open_file} function can be used to open and existing fits
%  file for reading or updating, or to create a new fits file, depending upon
%  the value of the \var{mode} parameter.  Specifically, if \var{mode} is 
%  \exmp{"r"}, the file will be opened for reading.  If \var{mode} is \exmp{"w"},
%  the file will be opened for updating (both reading and writing).  Otherwise, 
%  \var{mode} must be \var{"c"}, which indicates that a new file is to be created.
%  In the latter case, if a file already exists with the specified name, it will
%  get deleted and a new one created in its place.
% 
%  If the function fails, it will signal an error; otherwise an open file
%  pointer will be returned.
%\seealso{fits_close_file, fits_create_binary_table}
%!%-
define fits_open_file ()
{
   variable file, mode;

   if (_NARGS != 2)
     usage ("fp = fits_open_file (file, \"r|w|c\")");

   (file, mode) = ();
   variable fp;

   variable status = _fits_open_file (&fp, file, mode);
   if (status)
     do_fits_error (status, file);
   return fp;
}


%!%+
%\function{fits_close_file}
%\synopsis{Close a fits file}
%\usage{fits_close_file (Fits_File_Type f)}
%\description
%  The \var{fits_close_file} closes a previously opened fits file.  The function
%  will signal an error if the operation fails.
%\notes
%  This function could fail if it fails to write out any buffered data because
%  of filesystem errors (disk full, etc.).
%\seealso{fits_open_file}
%!%-
define fits_close_file (fp)
{
   do_fits_error (_fits_close_file (fp));
}

private define do_close_file (fp, needs_close)
{
   if (needs_close)
     fits_close_file (fp);
}

private define get_open_fp (fp, needs_close)
{
   @needs_close = 0;
   if (typeof (fp) != Fits_File_Type)
     {
	variable file = fp;
	do_fits_error (_fits_open_file (&fp, fp, "r"), file);
	@needs_close = 1;
     }
   return fp;
}


private define find_interesting_hdu (f, hdu_type, check_naxis)
{
   variable type;
   do
     {
	variable status;

	do_fits_error (_fits_get_hdu_type (f, &type));

	if ((hdu_type == NULL) or (type == hdu_type))
	  {
	     if (check_naxis == 0)
	       return 0;
	     
	     variable naxis;
	     do_fits_error (_fits_read_key (f, "NAXIS", &naxis, NULL));
	     if (naxis != 0)
	       return 0;
	  }

	status = _fits_movrel_hdu (f, 1);
     }
   while (not status);
   return -1;
}


private define get_open_hdu_of_type (f, hdu_type, needs_close, check_naxis)
{
   variable type;
   variable type_str;

   switch (hdu_type)
     {
      case _FITS_BINARY_TBL:
	type_str = "a binary table";
     }
     {
      case _FITS_IMAGE_HDU:
	type_str = "an image";
     }
     {
	% default
	type_str = "an interesting hdu";
     }
     
   @needs_close = 0;
   if (typeof (f) == Fits_File_Type)
     {
	do_fits_error (_fits_get_hdu_type (f, &type));
	if ((hdu_type != NULL) and (type != hdu_type))
	  verror ("Extension is not %s", type_str);
	return f;
     }

   f = get_open_fp (f, needs_close);
   
   if (0 == find_interesting_hdu (f, hdu_type, check_naxis))
     return f;

   verror ("Unable to locate %s", type_str);
}


%!%+
%\function{fits_move_to_interesting_hdu}
%\synopsis{Move to an extension that looks interesting}
%\usage{fits_move_to_interesting_hdu (fp [, hdu_type]}
%#v+
%  Fits_File_Type fp;
%  Int_Type hdu_type;
%#v-
%\description
%  The function move the fits file pointer \var{fp} forward to an HDU that looks 
%  interesting.  By definition, an interesting HDU is one in which NAXIS is 
%  non-zero.  The first parameter \var{fp} must be a pointer to an already open
%  fits file.  The second parameter, if present, may be used to specifiy the 
%  type of HDU, e.g., either an image (\exmp{hdu_type=_FITS_IMAGE_HDU}) or a 
%  binary table (\exmp{hdu_type=_FITS_BINARY_TBL}).
% 
%  If the function fails to find an interesting HDU of the appropriate type, 
%  an exception will be generated.
%\seealso{fits_open_file}
%!%-
define fits_move_to_interesting_hdu ()
{
   variable f, hdu_type = NULL;
   switch (_NARGS)
     {
      case 1:
	f = ();
     }
     {
      case 2:
	(f, hdu_type) = ();
     }
     {
	usage ("%s (f, hdu_type); % hdu_type = _FITS_IMAGE_HDU|_FITS_BINARY_TBL",
	       _function_name ());
     }

   if (-1 == find_interesting_hdu (f, hdu_type, 1))
     verror ("%s: Unable to find an interesting HDU", _function_name);
}

private define get_open_binary_table (f, needs_close)
{
   return get_open_hdu_of_type (f, _FITS_BINARY_TBL, needs_close, 1);
}

private define get_open_image_hdu (f, needs_close)
{
   return get_open_hdu_of_type (f, _FITS_IMAGE_HDU, needs_close, 1);
}

private define get_open_interesting_hdu (fp, needs_close)
{
   if (needs_close == NULL)
     {
	variable nc;
	needs_close = &nc;
     }

   @needs_close = 0;
   if (typeof (fp) == Fits_File_Type)
     return fp;

   return get_open_hdu_of_type (fp, NULL, needs_close, 1);
   %return get_open_binary_table (fp, needs_close);
}


%!%+
%\function{fits_key_exists}
%\synopsis{Check for the existence of a keyword}
%\usage{Int_Type fits_key_exists (fd, key)}
%#v+
%   Fits_File_Type or String_Type fd;
%   String_Type key;
%#v-
%\description
%  The \var{fits_key_exists} function checks for the existence of a specified 
%  keyword in the file specified by the descriptor \var{fd}, which must specify
%  the name of a file or an open file pointer.
% 
%  If the specified key exists, the function return \1, otherwise it returns \0.
%\seealso{fits_read_key, fits_read_header}
%!%-
define fits_key_exists ()
{
   if (_NARGS != 2)
     usage ("status = fits_key_exists (file, key)");

   variable fp, key;
   variable needs_close;

   (fp, key) = ();
   fp = get_open_interesting_hdu (fp, &needs_close);
   variable value;
   variable status = _fits_read_key (fp, key, &value, NULL);
   if (needs_close) fits_close_file (fp);
   if (status == 0)
     return 1;

   if (status == _FITS_KEY_NO_EXIST)
     return 0;

   do_fits_error (status);
}

private define get_fits_btable_info (fp)
{
   variable numrows, numcols, names, name, col;

   do_fits_error (_fits_get_num_rows (fp, &numrows));
   do_fits_error (_fits_get_num_cols (fp, &numcols));

   names = String_Type [numcols];
   _for (1, numcols, 1)
     {
	col = ();

	do_fits_error (_fits_read_key_string (fp, "TTYPE"+string(col), &name, NULL));
	names [col-1] = name;
     }

   return (numrows, names);
}


%!%+
%\function{fits_get_colnum}
%\synopsis{Get the column numbers of specified columns}
%\usage{column_num = fits_get_colnum (fd, column_name)}
%#v+
%   Fits_File_Type or String_Type fd;
%   String_Type column_name;
%#v-
%\description
%  This function returns the column number of the column with the specified name.
%  The file-descriptor \exmp{fd} must specify the name of a file, or an open
%  fits file pointer.
%\seealso{fits_binary_table_column_exists}
%!%-
define fits_get_colnum ()
{
   if (_NARGS != 2)
     usage ("colnum = %s (file, column_name)", _function_name ());

   variable f, column_names;  (f, column_names) = ();

   variable needs_close;
   f = get_open_binary_table (f, &needs_close);

   variable colnum;
   do_fits_error (_fits_get_colnum (f, column_names, &colnum));

   do_close_file (f, needs_close);
   return colnum;
}


%!%+
%\function{fits_binary_table_column_exists}
%\synopsis{Check for the existence of a binary table column}
%\usage{Int_Type fits_binary_table_column_exists (fd, col)}
%#v+
%   Fits_File_Type or String_Type fd;
%   String_Type col;
%#v-
%\description
%  This function may be used to determine whether or not a named column
%  exists in a binary table.  The table is specified via the \var{fd} 
%  parameter which must either be the name of a file containing the binary
%  table, or an file pointer.
% 
%  If the specified column exists, \1 will be returned; otherwise the function
%  will return \0.
%\seealso{fits_key_exists, fits_open_file}
%!%-
define fits_binary_table_column_exists ()
{
   if (_NARGS != 2)
     usage ("status = %s (file, column_name)", _function_name ());

   variable f, col;  (f, col) = ();

   variable needs_close;
   f = get_open_binary_table (f, &needs_close);

   variable names;
   (,names) = get_fits_btable_info (f);
   do_close_file (f, needs_close);

   col = strup (col);
   names = array_map (String_Type, &strup, names);
   return length (where (col == names));
}

private define get_tdim_string (fp, col)
{
   variable tdim = sprintf ("TDIM%d", col);
   !if (fits_key_exists (fp, tdim))
     return NULL;

   do_fits_error (_fits_read_key (fp, tdim, &tdim, NULL));
   return tdim;
}

private define make_tdim_string (dims)
{
   variable i;
   variable tdim = "(";

   dims = reverse (array_map (String_Type, &string, dims));
   return sprintf ("(%s)", strjoin (dims, ","));
}

private define convert_tdim_string (tdim, num_rows)
{
   tdim = strtrim (tdim, "()");
   tdim = reverse (strtok (tdim, " \t,"));

   tdim = array_map (Int_Type, &integer, tdim);
   if (num_rows == -1)
     return tdim;

   variable new_tdim = Int_Type[length(tdim)+1];
   new_tdim[0] = num_rows;
   new_tdim[[1:]] = tdim;
   return new_tdim;
}


private define check_vector_tdim (fp, first_row, col, data)
{
   variable tdim_col = sprintf ("TDIM%d", col);

   !if (fits_binary_table_column_exists (fp, tdim_col))
     return;
   
   do_fits_error (_fits_get_colnum (fp, tdim_col, &tdim_col));
   if (tdim_col == col)
     return;

   variable len = length (data);
   variable tdim;

   do_fits_error (_fits_read_col (fp, tdim_col, first_row, len, &tdim));
   if (_typeof (tdim) != String_Type)
     return;

   _for (0, len-1, 1)
     {
	variable i = ();
	reshape (data[i], convert_tdim_string (tdim[i], -1));
     }
}

% FITS column and keyword names can begin with a number or have dashes. 
% Bad Design.
private define normalize_names (names)
{
   names = @names;
   _for (0, length (names)-1, 1)
     {
	variable i = ();
	variable name = strlow (names[i]);
#iffalse
	(name,) = strreplace (name, "-", "_", strlen (name));
#else
	name = strtrans (name, "^a-z0-9", "_");
#endif
	names[i] = name;
	variable ch = name[0];
	if ((ch == '_')
	    or ((ch >= 'a') and (ch <= 'z')))
	  continue;
	names[i] = "_" + name;
     }
   return names;
}

private define get_column_number (fp, col)
{
   if (typeof (col) == String_Type)
     {
	do_fits_error (_fits_get_colnum (fp, col, &col));
	return col;
     }
   return int (col);
}

private define get_column_numbers (fp, args)
{
   variable column_nums = Int_Type[0];
   foreach (args)
     {
	variable arg = ();
	variable col = arg.value;
	if (typeof (col) == Array_Type)
	  col = array_map (Int_Type, &get_column_number, fp, col);
	else
	  col = get_column_number (fp, col);
	
	column_nums = [column_nums, col];
     }
   return column_nums;
}

% This function assumes that fp is an open pointer, and that columns is
% an array of column numbers.  The data are left on the stack.
private define read_cols (fp, columns, first_row, last_row)
{
   variable numrows;
   do_fits_error (_fits_get_num_rows (fp, &numrows));
   
   if (first_row < 0)
     first_row += (1+numrows);
   if (last_row < 0)
     last_row += (1+numrows);
   
   if ((first_row <= 0) or (last_row <= 0)
       or (first_row > numrows) or (last_row > numrows))
     verror ("Invalid first or last row parameters");

   variable numcols = length (columns);
   variable data_arrays;   
   do_fits_error (_fits_read_cols (fp, columns, first_row, last_row, &data_arrays));
   numrows = last_row - first_row + 1;
   _for (0, numcols-1, 1)
     {
	variable i = ();
	variable col = columns[i];
	variable data = data_arrays[i];
	variable tdim = get_tdim_string (fp, col);
	if (tdim != NULL)
	  {
	     tdim = convert_tdim_string (tdim, numrows);
	     reshape (data, tdim);
	  }
	else if (typeof (data) == Array_Type)
	  check_vector_tdim (fp, first_row, col, data);

	data;			       %  leave it on stack
     }
}

%!%+
%\function{fits_read_col}
%\synopsis{Read one or more columns from a FITS binary table}
%\usage{(x1, ...xN) = fits_read_col (file, c1, ... cN)}
%#v+
%   Fits_File_Type or String_Type file;
%   Int_Type or String_Type c1, ...cN;
%#v-
%\description
%  This function returns one or more vectors containing objects in the
%  specified columns of the binary table indicated by \var{file}.  If
%  \var{file} is a string, then the file will be opened via the virtual
%  file specification implied by \var{file}. Otherwise, \var{file}
%  should represent an already opened FITS file.  The column parameters
%  may either be strings denoting the column names, or integers
%  representing the column numbers.
%\seealso{fits_read_cell, fits_read_row, fits_read_table}
%!%-
define fits_read_col ()
{
   if (_NARGS < 2)
     usage ("(x1...xN) = fits_read_col (file, c1, ...cN)");

   variable fp, col;
   variable numrows;
   variable columns = __pop_args (_NARGS-1);
   fp = ();

   variable needs_close;
   fp = get_open_binary_table (fp, &needs_close);

   columns = get_column_numbers (fp, columns);

   read_cols (fp, columns, 1, -1);     %  data on stack

   do_close_file (fp, needs_close);
}

%!%+
%\function{fits_read_col_struct}
%\synopsis{Read one or more columns from a FITS binary table}
%\usage{struct = fits_read_col_struct (file, col1, ...)}
%#v+
%    Fits_File_Type or String_Type file;
%    String_Type col1, ...;
%#v-
%\description
%  This function works exactly like \var{fits_read_col} except it returns the
%  values in a structure.  See the documentation on that function for more
%  information.
% 
%\seealso{fits_read_col, fits_read_key_struct, fits_read_row, fits_read_header}
%!%-
define fits_read_col_struct ()
{
   !if (_NARGS)
     usage ("struct = fits_read_col_struct(file, COL1, ...)");
   
   variable cols = __pop_args (_NARGS - 1);
   variable file = ();
   variable fields = normalize_names ([__push_args(cols)]);
   variable s = @Struct_Type (fields);
   set_struct_fields (s, fits_read_col (file, __push_args (cols)));
   return s;
}

%!%+
%\function{fits_read_cell}
%\synopsis{Read a cell from a FITS binary table}
%\usage{X = fits_read_cell (file, c, r)}
%#v+
%   Fits_File_Type or String_Type file;
%   Int_Type r, c;
%#v-
%\description
%  This function returns the object in the column \var{c} and row
%  \var{r} of the binary table indicated by \var{file}.  If \var{file}
%  is a string, then the file will be opened via the virtual file
%  specification implied by \var{file}. Otherwise, \var{file} should
%  represent an already opened FITS file.
%\seealso{fits_read_col, fits_read_row}
%!%-
define fits_read_cell ()
{
   variable fp, r, c;
   variable needs_close;

   if (_NARGS != 3)
     usage ("x = fits_read_cell (file, c, r)");

   (fp, c, r) = ();
   fp = get_open_binary_table (fp, &needs_close);

   variable a = read_cols (fp, get_column_number (fp, c), r, r);
   variable dims, nd; (dims,nd,) = array_info (a);
   if (nd == 1)
     a = a[0];
   else
     reshape (a, dims[[1:]]);
   
   do_close_file (fp, needs_close);
   return a;
}

define fits_read_cells ()
{
   variable fp, r0, r1, columns;
   variable needs_close;

   if (_NARGS < 4)
     usage ("(x1,...xN) = %s (file, col1, ..., colN, r0, r1)", _function_name);

   (r0, r1) = ();
   columns = __pop_args (_NARGS-3);
   fp = ();
   fp = get_open_binary_table (fp, &needs_close);
   
   columns = get_column_numbers (fp, columns);
   read_cols (fp, columns, r0, r1);    %  on stack
   do_close_file (fp, needs_close);
}

%!%+
%\function{fits_read_row}
%\synopsis{Read a row from a FITS binary table}
%\usage{Struct_Type fits_read_cell (file, r)}
%#v+
%   Fits_File_Type or String_Type file;
%   Int_Type r;
%#v-
%\description
%  This function returns a structure containing the data in the columns
%  of the row \var{r} of the binary table indicated by \var{file}. If
%  \var{file} is a string, then the file will be opened via the virtual
%  file specification implied by \var{file}. Otherwise, \var{file}
%  should represent an already opened FITS file.
%\seealso{fits_read_col, fits_read_cell}
%!%-
define fits_read_row ()
{
   verror ("Not yet implemented");
}

			     
%!%+
%\function{fits_read_header}
%\synopsis{Read a FITS header}
%\usage{Struct_Type fits_read_header (file)}
%#v+
%    Fits_File_Type or String_Type file;
%#v-
%\description
%  This function reads the header of the fits file given by the
%  \var{file} parameter and returns it as a structure.  If \var{file} is
%  a string, then the file will be opened via the virtual file
%  specification implied by \var{file}. Otherwise, \var{file} should
%  represent an already opened FITS file.
%\seealso{fits_read_table}
%!%-
define fits_read_header ()
{
   !if (_NARGS)
     usage ("Struct_Type fits_read_header (file)");

   verror ("Not implemented");

   variable fp = ();
   variable needs_close;
   fp = get_open_fp (fp, &needs_close);
   do_close_file (fp, needs_close);
}


%!%+
%\function{fits_read_table}
%\synopsis{Read a FITS table}
%\usage{Struct_Type fits_read_table (file [,columns...])}
%#v+
%    Fits_File_Type or String_Type file;
%#v-
%\description
%  \var{fits_read_table} reads the data in a table of the FITS file
%  specified by \var{file} and returns it as a structure.  If the optional
%  column name parameters are specified, then only those columns will be read.
%  Otherwise, the entire table will be returned.
% 
%  If \var{file} is a string, then the file will be opened via the virtual file
%  specification implied by \var{file}. Otherwise, \var{file} should
%  represent an already opened FITS file.
%\seealso{fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
%!%-
define fits_read_table ()
{
   if (_NARGS == 0)
     usage ("S = fits_read_table (FILE [,columns,...])");

   variable f, names = NULL;
   if (_NARGS > 1)
     {
	names = __pop_args (_NARGS-1);
	names = [__push_args(names)];
     }
   f = ();
   variable needs_close;
   f = get_open_binary_table (f, &needs_close);

   if (names == NULL)
     (, names) = get_fits_btable_info (f);

   variable s = fits_read_col_struct (f, names);
   do_close_file (f, needs_close);
   return s;
}

define fits_info ()
{
   !if (_NARGS)
     %usage ("(numrows, numcols, colnames[]) = fits_info (file);");
     usage ("fits_info (file);");
   variable file = ();

   variable fp;
   variable numrows, numcols, names;
   variable needs_close;

   %do_fits_error (_fits_open_file (&fp, file, "r"));
   fp = get_open_interesting_hdu (file, &needs_close);

   (numrows, names) = get_fits_btable_info (fp);
   numcols = length (names);

   () = fprintf (stdout, "%s contains %d rows and %d columns:\n", file, numrows, numcols);
   _for (1, numcols, 1)
     {
	variable i = ();
	variable tform, name;

	name = names[i-1];
	do_fits_error (_fits_read_key_string (fp, "TFORM" + string(i), &tform, NULL));
	variable tdim = get_tdim_string (fp, i);
	if (tdim == NULL) tdim = "";
	else tdim = "TDIM=" + tdim;
	() = fprintf (stdout, "[%2d] %s %s %s\n", i, name, tform, tdim);
     }
   do_close_file (fp, needs_close);

   %return (numrows, numcols, names);
}


%!%+
%\function{fits_read_key}
%\synopsis{Read one or more keywords from a FITS file}
%\usage{(val1,...) = fits_read_key (file, key1, ...)}
%#v+
%    Fits_File_Type or String_Type file;
%    String_Type key1, ...;
%#v-
%\description
%  \var{fits_read_key} reads the values of one or more keywords in the fits
%  file specified by \var{file} and returns them.  If \var{file}
%  is a string, then the file will be opened via the virtual file
%  specification implied by \var{file}. Otherwise, \var{file} should
%  represent an already opened FITS file.  If any of the keywords do not exist,
%  a value of \NULL will be returned for the corresponding keyword.
%\seealso{fits_read_key_struct, fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
%!%-
define fits_read_key ()
{
   !if (_NARGS)
     usage ("(x,...) = fits_read_key (file, X_KEY, ...)");

   variable fp, keys;

   keys = __pop_args (_NARGS - 1);
   fp = ();
   
   variable needs_close;
   fp = get_open_interesting_hdu (fp, &needs_close);

   foreach (keys)
     {
	variable key = ().value;
	variable value, status;
	status = _fits_read_key (fp, key, &value, NULL);
	if (status == _FITS_KEY_NO_EXIST)
	  value = NULL;
	else if (status)
	  do_fits_error (status);

	value;
     }
   do_close_file (fp, needs_close);
}


%!%+
%\function{fits_read_key_struct}
%\synopsis{Read one or more keywords from a FITS file}
%\usage{struct = fits_read_key (file, key1, ...)}
%#v+
%    Fits_File_Type or String_Type file;
%    String_Type key1, ...;
%#v-
%\description
%  This function works exactly like \var{fits_read_key} excepts returns the
%  values in a structure.  See the documentation on that function for more
%  information.
%\seealso{fits_read_key, fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
%!%-
define fits_read_key_struct ()
{
   !if (_NARGS)
     usage ("struct = fits_read_key_struct(file, X_KEY, ...)");
   
   variable keys = __pop_args (_NARGS - 1);
   variable file = ();
   variable fields = normalize_names ([__push_args(keys)]);
   variable s = @Struct_Type (fields);
   set_struct_fields (s, fits_read_key (file, __push_args (keys)));
   return s;
}

		      
   
private define get_open_write_fp (fp, mode, needs_close)
{
   @needs_close = 0;
   if (typeof (fp) != Fits_File_Type)
     {
	@needs_close = 1;
	do_fits_error (_fits_open_file (&fp, fp, mode));
     }
   
   return fp;
}

%!%+
%\function{fits_create_binary_table}
%\synopsis{Prepare a binary table}
%\usage{fits_create_binary_table (file, extname, nrows, ttype, tform, tunit)}
%#v+
%    Fits_File_Type or String_Type file;
%    String_Type extname;
%    Int_Type nrows;
%    String_Type ttype[];
%    String_Type tform[];
%    String_Type tunit[];
%#v-
%\description
%  This creates a new binary table with the specified structure.  The parameters
%  \var{ttype}, \var{tform}, and \var{tunit} are string arrays that specify
%  the column names, column data type, and column units, respectively.
%  The binary table will be given the extension name \var{extname}.
%\seealso{fits_write_binary_table, fits_open_file}
%!%-
define fits_create_binary_table ()
{
   if (_NARGS != 6)
     usage ("fits_create_binary_table (file, extname, nrows, ttype[], tform[], tunit[])");

   variable fp, nrows, ttype, tform, tunit, extnam;

   (fp, extnam, nrows, ttype, tform, tunit) = ();

   variable needs_close;
   fp = get_open_write_fp (fp, "c", &needs_close);

   do_fits_error (_fits_create_binary_tbl (fp, nrows, ttype, tform, tunit, extnam));
   do_close_file (fp, needs_close);
}

%!%+
%\function{fits_write_binary_table}
%\synopsis{Write a binary table}
%\usage{fits_write_binary_table (file, extname, sdata, [skeys [,hist]])}
%#v+
%    Fits_File_Type or String_Type file;
%    String_Type extname;
%    Struct_Type sdata;
%    Struct_Type skeys;
%    Struct_Type hist;
%#v-
%\description
%  The \var{fits_write_binary_table} function creates a new binary table in
%  the specified file.  The parameter \var{file} specifies either a filename or
%  an open file pointer.  The \var{extname} parameter specifies the extension
%  name of the binary table.  The data written to table are specified in the 
%  \var{sdata} structure, where the name of the structure field specifies the 
%  column name.  If \var{skeys} is non-NULL, then it is a structure indicating
%  additional keywords to be written to the header of the binary table.  If the
%  optional parameter \var{hist} is present and non-NULL, then it is a structure
%  whose fields indicate either comment or history information to be written
%  to the header.
%\example
%  The following code
%#v+
%    variable data = struct { x, cosx, sinx };
%    data.x = [0:2*PI:0.01];
%    data.cosx = cos(data.x);
%    data.sinx = sin(data.x);
%
%    variable keys = struct { hduname, username};
%    keys.hduname = "COSXSINX";
%    keys.username = "John Doe";
%
%    variable hist = struct { history, comment};
%    hist.history = ["This is a history record", "This is another"];
%    hist.comment = ["This is a comment", "And this is another"];
%
%    fits_write_binary_table ("foo.fits", "COSXSINX", data, keys, hist);
%#v-
% produces a binary table with the header:
%#v+
%    XTENSION= 'BINTABLE' / binary table extension
%    BITPIX  =                   8 / 8-bit bytes
%    NAXIS   =                   2 / 2-dimensional binary table
%    NAXIS1  =                  24 / width of table in bytes
%    NAXIS2  =                 629 / number of rows in table
%    PCOUNT  =                   0 / size of special data area
%    GCOUNT  =                   1 / one data group (required keyword)
%    TFIELDS =                   3 / number of fields in each row
%    TTYPE1  = 'x       ' / label for field   1
%    TFORM1  = 'D       ' / data format of field: 8-byte DOUBLE
%    TTYPE2  = 'cosx    ' / label for field   2
%    TFORM2  = 'D       ' / data format of field: 8-byte DOUBLE
%    TTYPE3  = 'sinx    ' / label for field   3
%    TFORM3  = 'D       ' / data format of field: 8-byte DOUBLE
%    EXTNAME = 'COSXSINX' / name of this binary table extension
%    HDUNAME = 'COSXSINX'
%    USERNAME= 'John Doe'
%    HISTORY This is a history record
%    HISTORY This is another
%    COMMENT This is a comment
%    COMMENT And this is another
%#v-
%\notes
%  This function provides no mechanism to mix comments and keyword records.  As
%  the example shows, this function places the comment and history records at
%  the end of the table.
%\seealso{fits_create_binary_table, fits_open_file}
%!%-

private define add_keys_and_history_func (fp, keys, history)
{
   variable val;
   if (keys != NULL)
     {
	foreach (get_struct_field_names (keys))
	  {
	     variable keyword = ();
	     val = get_struct_field (keys, keyword);
	     do_fits_error (_fits_update_key (fp, keyword, val, NULL));
	  }
     }

   if (typeof (history) == String_Type)
     {
	history;
	history = struct {history}; 
	history.history = ();
     }

   if (history == NULL)
     return;

   foreach (get_struct_field_names (history))
     {
	keyword = ();
	val = get_struct_field (history, keyword);
	if (typeof (val) == String_Type)
	  val = [val];
	keyword = strlow (keyword);
	foreach (val)
	  {
	     val = ();
	     if (keyword == "history")
	       {
		  do_fits_error (_fits_write_history (fp, val));
		  continue;
	       }
	     if (keyword == "comment")
	       {
		  do_fits_error (_fits_write_comment (fp, val));
		  continue;
	       }
	     vmessage ("*** WARNING: history/comment record name '%s' is not supported",
		       history);
	  }
     }
}


define fits_write_binary_table ()
{
   variable fp, extname, s, keys, history;
   variable needs_close;
   variable keyfunc, keyfunc_args;

   variable usage_str = "\n"
     + "Form 1: fits_write_binary_table (file, extname, data_struct [,opt-keyword_struct [,opt-history]])\n"
     + "Form 2: fits_write_binary_table (file, extname, data_struct, &keyfunc [,opt-args...])";

   if (_NARGS < 3)
     usage (usage_str);

   keyfunc = NULL;
   if (_NARGS > 3)
     {
	_stk_reverse (_NARGS - 3);
	keyfunc = ();
	_stk_reverse (_NARGS - 4);
	
	if (typeof (keyfunc) != Ref_Type)
	  {
	     % keyfunc must be the keys struct
	     if (_NARGS == 4)
	       history = NULL;
	     else if (_NARGS == 5)
	       history = ();
	     else
	       {
		  _pop_n (_NARGS);
		  usage (usage_str);
	       }
	     (keyfunc, history);	       %  put back on stack
	     keyfunc_args = __pop_args (2);
	     keyfunc = &add_keys_and_history_func;
	  }
	else
	  keyfunc_args = __pop_args (_NARGS - 4);
     }

   (fp, extname, s) = ();

   fp = get_open_write_fp (fp, "c", &needs_close);

   variable ttype = get_struct_field_names (s);
   variable ncols = length (ttype);
   variable tform = String_Type [ncols];
   variable nrows = -1;
   variable tdim = String_Type[ncols];

   _for (0, ncols-1, 1)
     {
	variable i = ();
	variable colname = ttype[i];
	variable val = get_struct_field (s, colname);

	if (colname[0] == '_')	       %  unnormalize
	  colname = substr (colname, 2, -1);

	variable t = _typeof (val);
	variable ndims;
	
	switch (t)
	  {
	   case Int32_Type: t = "J";
	  }
	  {
	   case Float_Type: t = "E";
	  }
	  {
	   case Double_Type: t = "D";
	  }
	  {
	   case Int16_Type: t = "I";
	  }
	  {
	   case String_Type:
	     (,ndims,) = array_info (val);
	     if (ndims > 1)
	       verror ("This function does not support %d-d strings", ndims);

	     t = sprintf ("%dA", max (array_map (Int_Type, &strlen, val)));
	  }
	  {
	   case Char_Type or case UChar_Type:
	     t = "B";
	  }
	  {
	     verror ("%s: %s column: %S type not supported", _function_name, colname, t);
	  }

	variable nrows_i = length (val);
	if ((typeof (val) == Array_Type)
	    and nrows_i)
	  {
	     variable tdim_i;
	     (tdim_i,ndims,) = array_info (val);

	     if (ndims > 1)
	       {
		  t = string (nrows_i/tdim_i[0]) + t;
		  tdim[i] = make_tdim_string (tdim_i[[1:]]);
		  nrows_i = tdim_i[0];
	       }
	  }

	if (nrows != nrows_i)
	  {
	     if (nrows != -1)
	       verror ("Expecting field %s to have %d rows", ttype[i], nrows);
	     nrows = nrows_i;
	  }

	tform[i] = t;
	ttype[i] = colname;
     }

   fits_create_binary_table (fp, extname, nrows, ttype, tform, NULL);

   _for (0, ncols-1, 1)
     {
	i = ();
	if (NULL != tdim[i])
	  do_fits_error (_fits_update_key (fp, sprintf("TDIM%d", i+1), tdim[i], NULL));
     }

   if (keyfunc != NULL)
     (@keyfunc)(fp, __push_args(keyfunc_args));

   _for (0, ncols-1, 1)
     {
	i = ();
	val = get_struct_field (s, ttype[i]);
	do_fits_error (_fits_write_col (fp, i+1, 1, 1, val));
     }
   
   do_close_file (fp, needs_close);
}


private define do_write_xxx (func, nargs)
{
   variable args = __pop_args (nargs-1);
   variable fp = ();

   variable needs_close;
   fp = get_open_write_fp (fp, "w", &needs_close);

   if (nargs > 1)
     do_fits_error ((@func)(fp, __push_args(args)));
   else
     do_fits_error ((@func)(fp));

   do_close_file (fp, needs_close);
}

private define do_read_xxx (func, nargs)
{
   variable args = __pop_args (nargs-1);
   variable fp = ();

   variable needs_close;
   fp = get_open_fp (fp, &needs_close);

   if (nargs > 1)
     do_fits_error ((@func)(fp, __push_args(args)));
   else
     do_fits_error ((@func)(fp));

   do_close_file (fp, needs_close);
}


%!%+
%\function{fits_update_key}
%\synopsis{Update the value of a keyword}
%\usage{fits_update_key (fd, key, val [,comment])}
%#v+
%    String_Type or Fits_File_Type fd;
%    String_Type key;
%    Any type val;
%    String_Type comment;
%#v-
%\description
%  The \var{fits_update_key} function updates the value and comment fields
%  of an existing keyword with the specified name.  If the keyword does not 
%  exist, a new keyword will be appended to the end of the header.
%\seealso{fits_update_logical, fits_read_key}
%!%-
define fits_update_key ()
{
   variable nargs = _NARGS;
   if (nargs < 3)
     usage ("fits_update_key (fp, key, value, comment)");
   
   if (nargs == 3)
     {
	NULL;			       %  add comment
	nargs++;
     }

   do_write_xxx (&_fits_update_key, 4);
}

%!%+
%\function{fits_update_logical}
%\synopsis{Update the value of a logical (boolean) keyword}
%\usage{fits_update_logical (fd, key, val, comment)}
%#v+
%    String_Type or Fits_File_Type fd;
%    String_Type key;
%    Any type val;
%    String_Type comment;
%#v-
%\description
%  The \var{fits_update_logical} function updates the value and comment fields
%  of an existing keyword of the specified name with the specified boolean value.
%  If the keyword does not exist, a new keyword will be appended to the end of 
%  the header.
%\seealso{fits_update_key}
%!%-
define fits_update_logical ()
{
   if (_NARGS != 4)
     usage ("fits_update_key (fp, key, value, comment)");
   
   do_write_xxx (&_fits_update_logical, _NARGS);
}

%!%+
%\function{fits_write_comment}
%\synopsis{Write a comment to the header}
%\usage{fits_write_comment (fd, comment)}
%#v+
%  Fits_File_Type or String_Type fd;
%  String_Type comment;
%#v-
%\description
%  As the name indicates, this function writes a comment record to the specified
%  fits file.  The file-descriptor \exmp{fd} must either be the name of a fits
%  file or an open fits file pointer.
%\seealso{fits_update_key, fits_write_history}
%!%-
define fits_write_comment ()
{
   if (_NARGS != 2)
     usage ("fits_write_comment (fp, value)");

   do_write_xxx (&_fits_write_comment, _NARGS);
}

%!%+
%\function{fits_write_history}
%\synopsis{Write a history record to the header}
%\usage{fits_write_history (fd, history)}
%#v+
%  Fits_File_Type or String_Type fd;
%  String_Type history;
%#v-
%\description
%  As the name indicates, this function writes a history record to the specified
%  fits file.  The file-descriptor \exmp{fd} must either be the name of a fits
%  file or an open fits file pointer.
%\seealso{fits_update_key, fits_write_comment}
%!%-
define fits_write_history ()
{
   if (_NARGS != 2)
     usage ("fits_write_history (fp, value)");

   do_write_xxx (&_fits_write_history, _NARGS);
}

%!%+
%\function{fits_write_date}
%\synopsis{Write the DATE keyword to the current HDU}
%\usage{fits_write_date (fd)}
%#v+
%   Fits_File_Type or String_Type fd;
%#v-
%\description
%  The \sfun{fits_write_date} function calls \ifun{_fits_write_date} to write
%  the DATE to the header of the specified file descriptor, which  must either 
%  be the name of a fits file or an open fits file pointer.
%\seealso{fits_update_key}
%!%-
define fits_write_date ()
{
   if (_NARGS != 1)
     usage ("fits_write_date (fp)");
   do_write_xxx (&_fits_write_date, _NARGS);
}

%!%+
%\function{fits_write_chksum}
%\synopsis{Compute and write the DATASUM and CHECKSUM keywords}
%\usage{fits_write_chksum (fd)}
%#v+
%   Fits_File_Type or String_Type fd;
%#v-
%\description
%  The \sfun{fits_write_chksum} function calls \ifun{_fits_write_comment} to 
%  compute and write the DATASUM and CHECKSUM keywords to the 
%  header of the specified file descriptor, which  must either 
%  be the name of a fits file or an open fits file pointer.
%\seealso{fits_update_key, fits_verify_chksum}
%!%-
define fits_write_chksum ()
{
   if (_NARGS != 1)
     usage ("fits_write_chksum (fp)");
   do_write_xxx (&_fits_write_chksum, _NARGS);
}

%!%+
%\function{fits_verify_chksum}
%\synopsis{Verify the checksums for the current HDU}
%\usage{isok = fits_verify_chksum (fd [,dataok, hduok])}
%#v+
%   Fits_File_Type or String_Type fd;
%   Ref_Type dataok, hduok;
%#v-
%\description
%  The \sfun{fits_verify_chksum} function calls \ifun{_fits_verify_chksum} to 
%  verify the header and data checksums of the current HDU.  A non-zero return value
%  signifies that the checksums are ok, otherwise the function returns 0 to indicate
%  that the checksums are invalid.  The individual checksums of the HDU or data
%  can be checked through the use of the optional parameters.
%\seealso{fits_write_chksum}
%!%-
define fits_verify_chksum ()
{
   variable dataok_buf, hduok_buf;
   variable dataok = &dataok_buf, hduok = &dataok_buf;

   if (_NARGS == 3)
     (dataok, hduok) = ();
   else if (_NARGS != 1)
     usage ("ok = fits_verify_chksum (fp [,&dataok, &hduok])");

   if (dataok == NULL)
     dataok = &dataok_buf;
   if (hduok == NULL)
     hduok = &hduok_buf;

   &dataok, &hduok;		       %  push
   
   do_read_xxx (&_fits_verify_chksum, 3);
   
   return min([@dataok, @hduok]);
}

%!%+
%\function{fits_read_records}
%\synopsis{Read all the records in a fits header}
%\usage{String_Type[] fits_read_records (Fits_File_Type or String_Type fp)}
%\description
%  This function returns a list of all the header records associated with the
%  fits file descriptor as an array of strings.
%\seealso{fits_write_records, fits_read_key}
%!%-
define fits_read_records ()
{
   if (_NARGS != 1)
     usage ("String_Type[] fits_read_records (fp)");
   
   variable fp = ();
   fp = get_open_interesting_hdu (fp, NULL);

   variable nkeys;
   do_fits_error (_fits_get_num_keys (fp, &nkeys));

   variable recs = String_Type [nkeys];
   _for (0, nkeys-1, 1)
     {
	variable i = ();
	variable rec;

	do_fits_error (_fits_read_record (fp, i+1, &rec));
	recs[i] = rec;
     }
   return recs;
}

%!%+
%\function{fits_write_records}
%\synopsis{Write records to fits header}
%\usage{fits_write_records (fd, records)}
%#v+
%   Fits_File_Type or String_Type fd;
%   Array_Type records;
%#v-
%\description
%  This function uses the \ifun{_fits_write_record} function to write a series
%  of records to the current HDU.
%\seealso{fits_read_records}
%!%-
define fits_write_records ()
{
   if (_NARGS != 2)
     usage ("fits_read_records (fp, records[])");
   
   variable fp, records;
   (fp, records) = ();
   
   variable needs_close;
   fp = get_open_write_fp (fp, "w", &needs_close);

   if (String_Type == typeof (records))
     records = [records];

   foreach (records)
     {
	variable rec = ();
	do_fits_error (_fits_write_record (fp, rec));
     }
   do_close_file (fp, needs_close);
}


%!%+
%\function{fits_get_keyclass}
%\synopsis{Obtain the key classes for a set of cards}
%\usage{Int_Type[] = fits_get_keyclass (Array_Type cards)}
%\description
%  This function uses the \ifun{_fits_get_keyclass} function to obtain the 
%  key-classes associated with one or more cards.  The function returns an
%  integer-valued array of the same length as the \exmp{cards} array.
%\example
%  Obtain set of header cards to those that are not associated with the cards
%  describing the structure of the HDU:
%#v+
%    variable cards = fits_read_records ("evt2.fits[EVENTS]");
%    variable classes = fits_get_keyclass (cards);
%    cards = cards[where (classes != _FITS_TYP_STRUC_KEY)];
%#v-
%\seealso{fits_read_records, fits_read_key}
%!%-
define fits_get_keyclass ()
{
   if (_NARGS != 1)
     usage ("Int_Type[] = fits_get_keyclass (records)");
   
   variable records = ();
   if (String_Type == typeof (records))
     return _fits_get_keyclass (records);
   
   return array_map (Int_Type, &_fits_get_keyclass, records);
}

% Image routines

%!%+
%\function{fits_get_bitpix}
%\synopsis{Get the fits bitpix value for an array}
%\usage{Int_Type fits_get_bitpix (array)}
%\description
%  This function may be used to obtain the bitpix value for a specified image
%  array.  The array must be an integer or floating point type, otherwise
%  and error will be generated.  The bitpix value is returned.
%\seealso{fits_write_image_hdu, fits_read_img}
%!%-
define fits_get_bitpix (image)
{
   variable types = [Char_Type, UChar_Type, Int16_Type, UInt16_Type,
		     Int32_Type, UInt32_Type, Float32_Type, Float64_Type];
   variable bitpix = [8, 8, 16, 16, 32, 32, -32, -64];

   variable b;

   if (typeof (image) == DataType_Type)
     b = image;
   else
     b = _typeof (image);

   variable i = where (types == b);
   if (length (i) == 0)
     verror ("fits_get_bitpix: %S is not supported", b);
   
   return bitpix[i[0]];
}


%!%+
%\function{fits_read_img}
%\synopsis{Read image data from a fits file}
%\usage{Array_Type fits_read_img (fd)}
%#v+
%   Fits_File_Type or String_Type fd;
%#v-
%\description
%  This function reads an image from the specified file descriptor.  
%  The file descriptor must be either the name of an existing file, or an
%  open file pointer.  It returns the image upon sucess, or signals an error 
%  upon failure.
%\seealso{fits_read_table, fits_read_col, fits_open_file, fits_write_img}
%!%-
define fits_read_img ()
{
   !if (_NARGS)
     usage ("I=fits_read_img (file);");
   variable fp = ();
   
   variable needs_close;
   fp = get_open_image_hdu (fp, &needs_close);

   variable a;

   do_fits_error (_fits_read_img (fp, &a));
   do_close_file (fp, needs_close);

   return a;
}

%!%+
%\function{fits_create_image_hdu}
%\synopsis{Create a primary array or image extension}
%\usage{fits_create_image_hdu (fd, extname, type, dims)}
%#v+
%   Fits_File_Type or String_Type fd;
%   String_Type extname;
%   Array_Type dims;
%   DataType_Type type;
%#v-
%\description
%  This function make use of the \ifun{_fits_create_img} function to create an
%  image extension or primary array of the specified type and size.  If the
%  \exmp{extname} parameter is non-NULL, then an EXTNAME keyword will be 
%  written out with the value of the extname parameter.
%  The \exmp{dims} parameter must be a 1-d integer array that corresponds
%  to the dimensions of the array to be written.
%  
%  If \exmp{fd} is specified as a string, then a new file of that name will be 
%  created.  If a file by that name already exists, it will be deleted and
%  a new one created.  If this behavior is undesired, then explicitly open the
%  file and pass this routine the resulting file pointer.
%\seealso{fits_write_image_hdu}
%!%-
define fits_create_image_hdu ()
{
   if (_NARGS != 4)
     usage ("%s (file, extname, type, dims)", _function_name ());

   variable fp, extname, type, dims;

   (fp, extname, type, dims) = ();

   variable needs_close;
   fp = get_open_write_fp (fp, "c", &needs_close);

   do_fits_error (_fits_create_img (fp, fits_get_bitpix (type), dims));
   if (extname != NULL)
     do_fits_error (_fits_update_key (fp, "EXTNAME", extname, NULL));

   do_close_file (fp, needs_close);
}

%!%+
%\function{fits_write_image_hdu}
%\synopsis{Write an image extension}
%\usage{fits_write_image_hdu (file, extname, image [,skeys [,hist]])}
%#v+
%    Fits_File_Type or String_Type file;
%    String_Type extname;
%    Any_Type image
%    Struct_Type skeys;
%    Struct_Type hist;
%#v-
%\description
%  The \var{fits_write_image_hdu} function creates a new image extension in
%  the specified file.  The parameter \var{file} specifies either a filename or
%  an open file pointer.  The \var{extname} parameter specifies the extension
%  name of the image, or NULL for the primary image.  The image data written 
%  to the file are specified by the \var{image} parameter.
%  If the optional parameter \var{skeys} is non-NULL, then it is a 
%  structure indicating additional keywords to be written to the HDU. 
%  If the optional parameter \var{hist} is present and non-NULL, 
%  then it is a structure whose fields indicate either comment or history 
%  information to be written to the header.
%\example
%  The following code
%#v+
%     variable img = [1:128*128]; reshape (img, [128,128]);
%     variable keys = struct { hduname, username};
%     keys.hduname = "MY_IMAGE";
%     keys.username = "John Doe";
%     variable hist = struct { history, comment};
%     hist.history = ["This is a history record", "This is another"];
%     hist.comment = ["This is a comment", "And this is another"];
%     fits_write_image_hdu ("foo.fits", NULL, img, keys, hist);
%#v-
% produces an image HDU with the header:
%#v+
%     SIMPLE  =                   T / file does conform to FITS standard
%     BITPIX  =                  32 / number of bits per data pixel
%     NAXIS   =                   2 / number of data axes
%     NAXIS1  =                 128 / length of data axis 1
%     NAXIS2  =                 128 / length of data axis 2
%     EXTEND  =                   T / FITS dataset may contain extensions
%     COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
%     COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
%     HDUNAME = 'MY_IMAGE'
%     USERNAME= 'John Doe'
%     HISTORY This is a history record
%     HISTORY This is another
%     COMMENT This is a comment
%     COMMENT And this is another
%#v-
%\notes
%  This function provides no mechanism to mix comments and keyword records.  As
%  the example shows, this function places the comment and history records at
%  the end of the table.
%\seealso{fits_create_binary_table, fits_open_file}
%!%-
define fits_write_image_hdu ()
{
   variable fp, extname, image, keys = NULL, history = NULL;
   variable needs_close;

   switch (_NARGS)
     {
      case 4:
	keys = ();
     }
     {
      case 5:
	(keys, history) = ();
     }
     {
	if (_NARGS != 3)
	  {
	     _pop_n (_NARGS);
	     usage ("%s (file, extname, image [, keyword_struct [, history]]", _function_name ());
	  }
     }

   (fp, extname, image) = ();

   fp = get_open_write_fp (fp, "c", &needs_close);

   variable dims; (dims,,) = array_info (image);
   fits_create_image_hdu (fp, extname, _typeof (image), dims);
   
   if (keys != NULL)
     {
	foreach (get_struct_field_names (keys))
	  {
	     variable keyword = ();
	     variable val = get_struct_field (keys, keyword);
	     do_fits_error (_fits_update_key (fp, keyword, val, NULL));
	  }
     }

   if (typeof (history) == String_Type)
     {
	history;
	history = struct {history}; 
	history.history = ();
     }

   if (history != NULL)
     {
	foreach (get_struct_field_names (history))
	  {
	     keyword = ();
	     val = get_struct_field (history, keyword);
	     if (typeof (val) == String_Type)
	       val = [val];
	     keyword = strlow (keyword);
	     foreach (val)
	       {
		  val = ();
		  if (keyword == "history")
		    {
		       do_fits_error (_fits_write_history (fp, val));
		       continue;
		    }
		  if (keyword == "comment")
		    {
		       do_fits_error (_fits_write_comment (fp, val));
		       continue;
		    }
		  vmessage ("*** WARNING: history/comment record name '%s' is not supported",
			    history);
	       }
	  }
     }

   do_fits_error (_fits_write_img (fp, image));
   do_close_file (fp, needs_close);
}


%!%+
%\function{fits_write_img}
%\synopsis{Write the image data to an Image HDU}
%\usage{fits_write_img (Fits_File_Type fptr, Any_Type data)}
%\description
%  This function writes the image data out to current HDU, assumed to be 
%  an Image HDU.
%\seealso{fits_write_image_hdu, fits_create_image_hdu}
%!%-
% FIXME: Allow only a portion of the image to be written
define fits_write_img ()
{
   variable fp, data;

   switch (_NARGS)
     {
      case 2:
	(fp, data) = ();
     }
     {
	usage ("%s (fptr, img)", _function_name);
     }
   do_fits_error (_fits_write_img (fp, data));
}


#iffalse
define fits_iterate (fp, delta_rows, func, client_data, column_names)
{
   variable numrows;
   variable num_columns = length (column_names);
   variable col_nums = Int_Type[num_columns];
   variable data = Struct_Type[num_columns];
   variable i;
   
   for (i = 0; i < num_columns; i++)
     {
	variable col;
	do_fits_error (_fits_get_colnum (fp, column_names[i], &col));
	col_nums[i] = col;
	data[i] = struct { value };
     }

   do_fits_error (_fits_get_num_rows (fp, &numrows));

   variable current_row = 1;
   while (numrows)
     {
	variable value;

	if (numrows < delta_rows)
	  delta_rows = numrows;
	
	for (i = 0; i < num_columns; i++)
	  {
	     do_fits_error (_fits_read_col (fp, col_nums[i], current_row, 
					    delta_rows, &value));
	     data[i].value = value;
	  }
	@func (client_data, __push_args (data));
	
	current_row += delta_rows;
	numrows -= delta_rows;
     }
}

define test_func (info, x, y, ccdid, grade, status)
{
   variable i = where ((ccdid == info.ccdid)
		       and ((grade != 1) and (grade != 5) and (grade != 7))
		       and (status == 0));
   info.sum_x += sum (x[i]);
   info.sum_y += sum (y[i]);
   info.num += length (x[i]);
}


define test_fits_iterate ()
{
   variable file = "/tmp/test.fits[EVENTS]";
   
   variable info = struct 
     {
	sum_x, sum_y, ccdid, num
     };
   variable fp = fits_open_file (file, "r");
   
   variable delta = 1;
   while (delta < 100000000)
     {
	info.sum_x = 0;
	info.sum_y = 0;
	info.ccdid = 7;
	info.num = 0;

	tic ();
	fits_iterate (fp, 10000, &test_func, info, 
		      ["X", "Y", "CCD_ID", "GRADE", "STATUS"]);

	() = fprintf (stdout, "delta=%d, CPU=%g secs, mean([x,y]) is [%g,%g]\n",
		      delta, toc (), info.sum_x/info.num, info.sum_y/info.num);
	() = fflush (stdout);
	delta *= 10;
     }
   
}
#endif

% Obsolete functions

define fits_read_image ()
{
   () = fprintf (stderr, "*** Warning: fits_read_image is deprecated.\n");
   variable args = __pop_args (_NARGS);
   return fits_read_img (__push_args (args));
}

provide ("fits");

#ifexists add_doc_file
$1 = path_concat (path_concat (path_dirname (__FILE__), "help"),
		  "cfitsio.hlp");
if (NULL != stat_file ($1))
  add_doc_file ($1);
#endif

#iffalse
autoload ("fitswcs_get_img_wcs", "fitswcs.sl")
autoload ("fitswcs_get_column_wcs", "fitswcs.sl")
autoload ("fitswcs_write_img_wcs", "fitswcs.sl")
autoload ("fitswcs_slice", "fitswcs.sl")
#endif
