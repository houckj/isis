%    Copyright (C) 1998-2010 Massachusetts Institute of Technology
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

% Track the modules's version since it now follows the changes.txt file.
variable _fits_sl_version_string = _cfitsio_module_version_string;
variable _fits_sl_version = _cfitsio_module_version;

private variable Verbose = 1;
% Forward declarations

#ifexists new_exception
if (0 == is_defined ("FitsError"))
  new_exception ("FitsError", RunTimeError, "Fits Error");
#endif

private variable Last_Error_Messages = String_Type[0];

private define reverse (a)
{
#ifexists array_reverse
   a = @a;
   array_reverse (a);
   return a;
#else
   variable i = length (a);
   if (i <= 1)
     return a;

   i--;
   __tmp(a)[[i:0:-1]];
#endif
}

%!%+
%\function{fits_read_errmsgs}
%\synopsis{Retrieve all error messages from the CFITSIO error stack}
%\usage{String_Type[] fits_read_errmsgs ()}
%\description
% This function returns all the error messages from the CFITSIO error
% message stack as an array of strings.
%\notes
% Using this function will cause the error message stack to be cleared.
%\seealso{_fits_read_errmsg, fits_set_verbose_errors}
%!%-
define fits_read_errmsgs ()
{
   variable err;
   variable errlist = String_Type[0];
   while (err = _fits_read_errmsg (), err != NULL)
     {
	errlist = [errlist, err];
     }
   return errlist;
}

define fits_get_errmsgs ()
{
   return Last_Error_Messages;
}

%!%+
%\function{fits_set_verbose_errors}
%\synopsis{Set the verbosity level of the cfitsio error messages}
%\usage{fits_set_verbose_errors (Int_Type level)}
%\description
% When a call to a function in the high-level interface fails, a error
% message will get generated.  By default, all messages from the
% underlying cfitsio error stack are printed.  This behavior may be
% turned off by calling this function with \exmp{level} equal to 0.
%\seealso{fits_read_errmsgs}
%!%-
define fits_set_verbose_errors ()
{
   variable v = 1;
   if (_NARGS)
     v = ();

   Verbose = v;
}

define fits_check_error ()
{
   variable status, file = "";

   if (_NARGS == 2)
     file = ();
   status = ();
   if (status == 0)
     {
	_fits_clear_errmsg ();
	Last_Error_Messages = String_Type[0];
	return;
     }

   if (strlen(file))
     file = strcat (": ", file);
   variable errmsg = strcat (_fits_get_errstatus (status), file);
   Last_Error_Messages = fits_read_errmsgs ();
   if (Verbose)
     errmsg = strjoin ([Last_Error_Messages, errmsg], "\n");

#ifexists new_exception
   throw FitsError, errmsg, Last_Error_Messages;
#else
   error (errmsg);
#endif
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
     fits_check_error (status, file);
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
   fits_check_error (_fits_close_file (fp));
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
	fits_check_error (_fits_open_file (&fp, fp, "r"), file);
	@needs_close = 1;
     }
   return fp;
}

define fits_get_num_hdus ()
{
   if (_NARGS != 1)
     usage ("num = fits_get_num_hdus (file)");
   variable fp = ();
   variable needs_close, num;
   fp = get_open_fp (fp, &needs_close);
   variable status = _fits_get_num_hdus (fp, &num);
   do_close_file (fp, needs_close);
   if (status == 0)
     return num;

   fits_check_error (status);
}

define fits_movrel_hdu ()
{
   if (_NARGS != 2)
     usage ("num_actual = fits_movrel_hdu (fp, num_to_move)");
   variable fp, num_to_move, num_hdus, this_num;
   (fp, num_to_move) = ();
   this_num = _fits_get_hdu_num (fp);
   if (num_to_move > 0)
     {
	num_hdus = fits_get_num_hdus (fp);
	if (this_num + num_to_move > num_hdus)
	  num_to_move = num_hdus - this_num;
     }
   else if (num_to_move < 0)
     {
	if (this_num + num_to_move <= 0)
	  num_to_move = -this_num;
     }
   if (num_to_move == 0)
     return 0;

   fits_check_error (_fits_movrel_hdu (fp, num_to_move));
   return num_to_move;
}

define fits_movabs_hdu ()
{
   if (_NARGS != 2)
     usage ("fits_movabs_hdu (fp, n)");
   variable fp,n;
   (fp, n) = ();
   fits_check_error (_fits_movabs_hdu (fp, n));
}

private define find_interesting_hdu (f, hdu_type, check_naxis)
{
   variable type;
   do
     {
	variable status;

	fits_check_error (_fits_get_hdu_type (f, &type));
	if (type == _FITS_ASCII_TBL)
	  type = _FITS_BINARY_TBL;

	if ((hdu_type == NULL) or (type == hdu_type))
	  {
	     if (check_naxis == 0)
	       return 0;

	     variable naxis;
	     fits_check_error (_fits_read_key (f, "NAXIS", &naxis, NULL));
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
      case _FITS_ASCII_TBL:
	type_str = "an ascii table";
     }
     {
	% default
	type_str = "an interesting hdu";
     }

   @needs_close = 0;
   if (typeof (f) == Fits_File_Type)
     {
	fits_check_error (_fits_get_hdu_type (f, &type));
	if ((hdu_type == _FITS_BINARY_TBL)
	    and (type == _FITS_ASCII_TBL))
	  hdu_type = type;

	if ((hdu_type != NULL) and (type != hdu_type))
	  {
	     verror ("Extension is not %s", type_str);
	  }
	return f;
     }

   variable name_contains_extno
     = string_match (f, "\[.+\]$"R, 1) || string_match (f, "+\d+$"R, 1);

   f = get_open_fp (f, needs_close);

   if (name_contains_extno
       || (0 == find_interesting_hdu (f, hdu_type, check_naxis)))
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
   do_close_file (fp, needs_close);
   if (status == 0)
     return 1;

   if (status == _FITS_KEY_NO_EXIST)
     return 0;

   fits_check_error (status);
}

private define get_fits_btable_info (fp)
{
   variable numrows, numcols, names, name, col;

   fits_check_error (_fits_get_num_rows (fp, &numrows));
   fits_check_error (_fits_get_num_cols (fp, &numcols));

   names = String_Type [numcols];
   _for (1, numcols, 1)
     {
	col = ();

	fits_check_error (_fits_read_key_string (fp, "TTYPE"+string(col), &name, NULL));
	names [col-1] = name;
     }

   return (numrows, names);
}

private define get_casesens_qualifier ()
{
   % We allow the forms:
   %   func (...; casesen);
   %   func (...; casesen=0);
   %   func (...; casesen=1);
   ifnot (qualifier_exists ("casesen"))
     return 0;

   variable q = qualifier ("casesen");
   return (q == NULL) ? 1 : q;
}

private define _fits_get_colnum_maybe_casesen (ft, name, ref, casesen)
{
   if (casesen)
     return _fits_get_colnum_casesen (ft, name, ref);

   return _fits_get_colnum (ft, name, ref);
}

private define get_column_number (fp, col, casesen)
{
   if (typeof (col) == String_Type)
     {
	variable col_str = col;
	fits_check_error (_fits_get_colnum_maybe_casesen (fp, col_str, &col, casesen),
			  col_str);
	return col;
     }
   return int (col);
}

private define get_column_numbers (fp, col_list, casesen);
private define get_column_numbers (fp, col_list, casesen)
{
   variable column_nums = Int_Type[0];

   foreach (col_list)
     {
	variable col = ();
	if ((typeof (col) == Array_Type) || (typeof (col) == List_Type))
	  col = get_column_numbers (fp, col, casesen);
	else
	  col = get_column_number (fp, col, casesen);

	column_nums = [column_nums, col];
     }
   return column_nums;
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
%\qualifiers
%\qualifier{casesen}{use case-sensitive column names}
%\seealso{fits_binary_table_column_exists}
%!%-
define fits_get_colnum ()
{
   if (_NARGS != 2)
     usage ("colnum = %s (file, column_name)", _function_name ());

   variable f, column_names;  (f, column_names) = ();

   variable needs_close;
   f = get_open_binary_table (f, &needs_close);

   variable colnum, casesen = get_casesens_qualifier (;;__qualifiers);
   fits_check_error (_fits_get_colnum_maybe_casesen (f, column_names, &colnum, casesen));

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
%\qualifiers
%\qualifier{casesen}{use case-sensitive column names}
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

   ifnot (get_casesens_qualifier (;;__qualifiers))
     {
	col = strup (col);
	names = array_map (String_Type, &strup, names);
     }
   return length (where (col == names));
}

define fits_delete_col ()
{
   if (_NARGS != 2)
     usage ("fits_delete_col (file, col)");

   variable f, col;  (f, col) = ();
   variable needs_close;
   f = get_open_binary_table (f, &needs_close);

   col = get_column_number (f, col, get_casesens_qualifier(;;__qualifiers));
   fits_check_error (_fits_delete_col (f, col));
   do_close_file (f, needs_close);
}

private define get_tdim_string (fp, col)
{
   variable tdim = sprintf ("TDIM%d", col);
   !if (fits_key_exists (fp, tdim))
     return NULL;

   fits_check_error (_fits_read_key (fp, tdim, &tdim, NULL), tdim);
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

private define check_vector_tdim (fp, first_row, tdim_col, data)
{
   if (tdim_col == 0)
     return;

   variable len = length (data);
   variable tdim;

   fits_check_error (_fits_read_col (fp, tdim_col, first_row, len, &tdim));
   if (_typeof (tdim) != String_Type)
     return;

   _for (0, len-1, 1)
     {
	variable i = ();
	reshape (data[i], convert_tdim_string (tdim[i], -1));
     }
}

private define reshape_string_array (fp, col, data)
{
   variable tform, repeat, width;
   fits_check_error (_fits_read_key_string (fp, "TFORM" + string(col), &tform, NULL));
   % Look for rAw
   if (2 != sscanf (tform, "%dA%d", &repeat, &width))
     return data;

   variable num_substrs = repeat/width;
   if (num_substrs == 1)
     return data;

   variable len = length (data);
   variable new_data = String_Type[len, num_substrs];
   variable i, j;
   _for i (0, len-1, 1)
     {
	variable str = data[i];
	_for j (0, num_substrs-1, 1)
	  {
	     variable new_str = substrbytes (str, 1+j*width, width);
	     new_str = strtrim_end (new_str, " ");
	     if ((new_str == "") && width)
	       new_str = " ";
	     new_data[i,j] = new_str;
	  }
     }
   return new_data;
}

% FITS column and keyword names can begin with a number or have dashes.
% Bad Design.
private define normalize_names (names, casesen);
private define normalize_names (names, casesen)
{
   variable new_names = String_Type[0];
   _for (0, length (names)-1, 1)
     {
	variable i = ();
	variable name = names[i];
	variable t = typeof (name);
	if ((t == Array_Type) || (t == List_Type))
	  {
	     new_names = [new_names, normalize_names(name, casesen)];
	     continue;
	  }

	ifnot (casesen)
	  name = strlow (name);
	name = strtrans (name, "^a-zA-Z0-9", "_");
	if ('0' <= name[0] <= '9')
 	  name = "_" + name;

	new_names = [new_names, name];
     }
   return new_names;
}

private define open_read_cols (fp, columns)
{
   variable needs_close, numrows, numcols;
   fp = get_open_binary_table (fp, &needs_close);
   fits_check_error (_fits_get_num_rows (fp, &numrows));
   numcols = length(columns);

   variable s = struct
     {
	fp = fp,
	needs_close = needs_close,
	columns = get_column_numbers (fp, columns, get_casesens_qualifier(;;__qualifiers)),
	num_rows = numrows, num_cols = numcols,
	tdims = String_Type[numcols],
	tdim_cols = Int_Type[numcols],
     };

   _for (0, numcols-1, 1)
     {
	variable i = ();
	variable col = s.columns[i];
	variable tdim = get_tdim_string (fp, col);
	if (tdim != NULL)
	  s.tdims[i] = tdim;
	variable tdim_col = sprintf ("TDIM%d", col);
	if (fits_binary_table_column_exists (fp, tdim_col))
	  {
	     fits_check_error (_fits_get_colnum (fp, tdim_col, &tdim_col));
	     if (tdim_col != col)
	       s.tdim_cols[i] = tdim_col;
	  }
     }

   return s;
}

private define close_read_cols (s)
{
   do_close_file (s.fp, s.needs_close);
}

% This function assumes that fp is an open pointer, and that columns is
% an array of column numbers.  The data are left on the stack.
private define read_cols (fpinfo, first_row, last_row)
{
   variable
     fp = fpinfo.fp,
     numrows = fpinfo.num_rows,
     columns = fpinfo.columns,
     tdims = fpinfo.tdims,
     tdim_cols = fpinfo.tdim_cols;

   if (first_row < 0)
     first_row += (1+numrows);
   if (last_row < 0)
     last_row += (1+numrows);

   variable want_num_rows = last_row - first_row + 1;
   if ((first_row <= 0) or (last_row < 0)
       or (want_num_rows > numrows) or (want_num_rows < 0))
     throw FitsError, "Invalid first or last row parameters";

   variable data_arrays;
   fits_check_error (_fits_read_cols (fp, columns, first_row, want_num_rows, &data_arrays));
   _for (0, fpinfo.num_cols-1, 1)
     {
	variable i = ();
	variable col = columns[i];
	variable data = data_arrays[i];
	variable tdim = tdims[i];
	if (tdim != NULL)
	  {
	     tdim = convert_tdim_string (tdim, want_num_rows);
	     reshape (data, tdim);
	  }
	else if (typeof (data) == Array_Type)
	  {
	     if (_typeof (data) == String_Type)
	       data = reshape_string_array (fp, col, data);
	     if (tdim_cols[i]>0)
	       check_vector_tdim (fp, first_row, tdim_cols[i], data);
	  }
	data;			       %  leave it on stack
     }
}

private define pop_column_list (nargs)
{
   variable list = {};
   _stk_reverse (nargs);	       %  preserve the order
   loop (nargs)
     {
	variable arg = ();
	variable t = typeof(arg);
	if ((t == Array_Type) || (t == List_Type))
	  {
	     foreach arg (arg)
	       list_append (list, arg);
	     continue;
	  }
	list_append (list, arg);
     }
   return list;
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
%\qualifiers
%\qualifier{casesen}{use case-sensitive column names}
%\seealso{fits_read_cell, fits_read_row, fits_read_table}
%!%-
define fits_read_col ()
{
   if (_NARGS < 2)
     usage ("(x1...xN) = fits_read_col (file, c1, ...cN [;row=val, num=val])");

   variable cols = pop_column_list (_NARGS-1);
   variable fp = ();
   variable fpinfo = open_read_cols (fp, cols;; __qualifiers);

   variable first_row, last_row, num;
   first_row = qualifier ("row", 1);
   num = qualifier ("num");
   if (num == NULL)
     last_row = -1;
   else
     last_row = first_row + num - 1;

   read_cols (fpinfo, first_row, last_row);     %  data on stack

   close_read_cols (fpinfo);
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
%  Field names are converted to lowercase unless the \exmp{casesen} qualifier is set.
%\qualifiers
%\qualifier{casesen}{use case-sensitive column names}
%\seealso{fits_read_col, fits_read_key_struct, fits_read_row, fits_read_header}
%!%-
define fits_read_col_struct ()
{
   !if (_NARGS)
     usage ("struct = fits_read_col_struct(file, COL1, ...)");

   variable cols = pop_column_list (_NARGS - 1);
   variable file = ();
   variable fields = normalize_names (cols, get_casesens_qualifier (;;__qualifiers));
   variable s = @Struct_Type (fields);
   set_struct_fields (s, fits_read_col (file, cols;; __qualifiers));
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
   variable fpinfo = open_read_cols (fp, [c] ;; __qualifiers);

   variable a = read_cols (fpinfo, r, r);
   variable dims, nd; (dims,nd,) = array_info (a);
   if (nd == 1)
     a = a[0];
   else
     reshape (a, dims[[1:]]);

   close_read_cols (fpinfo);
   return a;
}

define fits_read_cells ()
{
   variable fp, r0, r1, columns;
   variable needs_close;

   if (_NARGS < 4)
     usage ("(x1,...xN) = %s (file, col1, ..., colN, r0, r1)", _function_name);

   (r0, r1) = ();
   columns = pop_column_list (_NARGS-3);
   fp = ();

   variable fpinfo = open_read_cols (fp, columns;; __qualifiers);
   read_cols (fpinfo, r0, r1);    %  on stack
   close_read_cols (fpinfo);
}

define fits_get_num_rows ()
{
   if (_NARGS != 1)
     usage ("nrows = fits_get_num_rows (file)");

   variable fp = ();
   variable needs_close, num_rows;
   fp = get_open_binary_table (fp, &needs_close);
   fits_check_error (_fits_get_num_rows (fp, &num_rows));
   do_close_file (fp, needs_close);
   return num_rows;
}

define fits_get_num_cols ()
{
   if (_NARGS != 1)
     usage ("ncols = fits_get_num_cols (file)");

   variable fp = ();
   variable needs_close, num_cols;
   fp = get_open_binary_table (fp, &needs_close);
   fits_check_error (_fits_get_num_cols (fp, &num_cols));
   do_close_file (fp, needs_close);
   return num_cols;
}

%!%+
%\function{fits_read_row}
%\synopsis{Read a row from a FITS binary table}
%\usage{Struct_Type fits_read_row (file, r)}
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
%  \sfun{fits_read_table} reads the data in a table of the FITS file
%  specified by \exmp{file} and returns it as a structure.  Field
%  names are converted to lowercase unless the \exmp{casesen} qualifier
%  is set.  If the optional column name parameters are specified,
%  then only those columns will be read.  Otherwise, the entire table
%  will be returned.
%
%  If \exmp{file} is a string, then the file will be opened via the virtual file
%  specification implied by \exmp{file}. Otherwise, \exmp{file} should
%  represent an already opened FITS file.
%\qualifiers
%\qualifier{casesen}{do not convert field names to lowercase}
%\seealso{fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
%!%-
define fits_read_table ()
{
   if (_NARGS == 0)
     usage ("S = fits_read_table (FILE [,columns,...] [;casesen])");

   variable f, names = NULL;
   if (_NARGS > 1)
     names = pop_column_list (_NARGS-1);

   f = ();
   variable needs_close;
   f = get_open_binary_table (f, &needs_close);

   if (names == NULL)
     (, names) = get_fits_btable_info (f);

   variable s = fits_read_col_struct (f, names;; __qualifiers);
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

   %fits_check_error (_fits_open_file (&fp, file, "r"));
   fp = get_open_interesting_hdu (file, &needs_close);

   (numrows, names) = get_fits_btable_info (fp);
   numcols = length (names);

   () = fprintf (stdout, "%S contains %d rows and %d columns:\n", file, numrows, numcols);
   _for (1, numcols, 1)
     {
	variable i = ();
	variable tform, name;

	name = names[i-1];
	fits_check_error (_fits_read_key_string (fp, "TFORM" + string(i), &tform, NULL));
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
	  {
	     value = NULL;
	     _fits_clear_errmsg ();
	  }
	else if (status)
	  fits_check_error (status, key);

	value;
     }
   do_close_file (fp, needs_close);
}

%!%+
%\function{fits_read_key_struct}
%\synopsis{Read one or more keywords from a FITS file}
%\usage{struct = fits_read_key_struct (file, key1, ...)}
%#v+
%    Fits_File_Type or String_Type file;
%    String_Type key1, ...;
%#v-
%\description
%  This function works exactly like \sfun{fits_read_key} except it returns the
%  values in a structure.  See the documentation on that function for more
%  information.
%
%  Field names are converted to lowercase unless the \exmp{casesen}
%  qualifier is set.
%\qualifiers
%\qualifier{casesen}{do not convert field names to lowercase}
%\seealso{fits_read_key, fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
%!%-
define fits_read_key_struct ()
{
   !if (_NARGS)
     usage ("struct = fits_read_key_struct(file, X_KEY, ...)");

   variable keys = __pop_args (_NARGS - 1);
   variable file = ();
   variable fields = normalize_names ([__push_args(keys)], get_casesens_qualifier(;;__qualifiers));
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
	fits_check_error (_fits_open_file (&fp, fp, mode));
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

   fits_check_error (_fits_create_binary_tbl (fp, nrows, ttype, tform, tunit, extnam));
   do_close_file (fp, needs_close);
}

%!%+
%\function{fits_write_binary_table}
%\synopsis{Write a binary table}
%\usage{fits_write_binary_table (file, extname, sdata, [skeys [,hist]])}
%#v+
%Fits_File_Type or String_Type file;
%String_Type extname;
%Struct_Type sdata;
%Struct_Type skeys;
%Struct_Type hist;
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
	     if (keyword[0] == '_') keyword = keyword[[1:]];   %  HACK!!! FIXME
	     fits_check_error (_fits_update_key (fp, keyword, val, NULL), keyword);
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
		  fits_check_error (_fits_write_history (fp, val));
		  continue;
	       }
	     if (keyword == "comment")
	       {
		  fits_check_error (_fits_write_comment (fp, val));
		  continue;
	       }
	     vmessage ("*** WARNING: history/comment record name '%s' is not supported",
		       history);
	  }
     }
}

private define shape_columns_before_write (s, ncols, ttype)
{
   variable reshapes_to = Array_Type[ncols];
   variable ndims, dims;

   _for (0, ncols-1, 1)
     {
	variable i = ();
	variable val = get_struct_field (s, ttype[i]);
	(dims,ndims,) = array_info (val);
	if (ndims > 1)
	  {
	     variable dim_0 = dims[0];
	     if (dim_0 != 0)
	       {
		  reshape (val, [dim_0, length(val)/dim_0]);
		  reshapes_to[i] = dims;
	       }
	  }
     }
   return reshapes_to;
}

private define unshape_columns_after_write (s, ncols, ttype, reshapes_to)
{
   _for (0, ncols-1, 1)
     {
	variable i = ();
	if (reshapes_to[i] == NULL)
	  continue;

	reshape (get_struct_field (s, ttype[i]), reshapes_to[i]);
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

   variable ttype;
   if (s == NULL)
     ttype = String_Type[0];
   else
     ttype = get_struct_field_names (s);
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
	   case UInt16_Type: t = "U";
	  }
	  {
	   case UInt32_Type: t = "V";
	  }
	  {
	   case String_Type:
	     (,ndims,) = array_info (val);
	     if (ndims > 1)
	       verror ("This function does not support %d-d strings", ndims);

	     t = sprintf ("%dA", max (array_map (Int_Type, &strlen, val)));
	     if (t == "0A") t = "1A";
	  }
	  {
	   case Char_Type or case UChar_Type:
	     t = "B";
	  }
	  {
	   case Int64_Type: t = "K";
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

   if (nrows == -1)		       %  ncols is 0
     nrows = 0;
   fits_create_binary_table (fp, extname, nrows, ttype, tform, NULL);

   _for (0, ncols-1, 1)
     {
	i = ();
	if (NULL != tdim[i])
	  fits_check_error (_fits_update_key (fp, sprintf("TDIM%d", i+1), tdim[i], NULL));
     }

   if (keyfunc != NULL)
     (@keyfunc)(fp, __push_args(keyfunc_args));

#iffalse
   _for (0, ncols-1, 1)
     {
	i = ();
	val = get_struct_field (s, ttype[i]);
	fits_check_error (_fits_write_col (fp, i+1, 1, 1, val));
     }
#else

   variable reshapes_to = shape_columns_before_write (s, ncols, ttype);

# ifeval (_slang_version >= 20000)
   try
     {
# else
	ERROR_BLOCK
	  {
	     unshape_columns_after_write (s, ncols, ttype, reshapes_to);
	  }
# endif
	variable r = 0;
	variable drows = 10;
	while (r < nrows)
	  {
	     variable r1 = r + nrows;
	     if (r1 > nrows)
	       r1 = nrows;

	     variable k = [r:r1-1];

	     _for (0, ncols-1, 1)
	       {
		  i = ();
		  val = get_struct_field (s, ttype[i]);
		  if (reshapes_to[i] == NULL)
		    fits_check_error (_fits_write_col (fp, i+1, r+1, 1, val[k]));
		  else
		    fits_check_error (_fits_write_col (fp, i+1, r+1, 1, val[k,*]));
	       }
	     r = r1;
	  }
#ifeval (_slang_version >= 20000)
     }
   finally
     {
	unshape_columns_after_write (s, ncols, ttype, reshapes_to);
     }
#endif

   do_close_file (fp, needs_close);
}

private define do_write_xxx (func, nargs)
{
   variable args = __pop_args (nargs-1);
   variable fp = ();

   variable needs_close;
   fp = get_open_write_fp (fp, "w", &needs_close);

   if (nargs > 1)
     fits_check_error ((@func)(fp, __push_args(args)));
   else
     fits_check_error ((@func)(fp));

   do_close_file (fp, needs_close);
}

private define do_read_xxx (func, nargs)
{
   variable args = __pop_args (nargs-1);
   variable fp = ();

   variable needs_close;
   fp = get_open_fp (fp, &needs_close);

   if (nargs > 1)
     fits_check_error ((@func)(fp, __push_args(args)));
   else
     fits_check_error ((@func)(fp));

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

   do_write_xxx (&_fits_update_key, nargs);
}

define fits_delete_key ()
{
   if (_NARGS != 2)
     usage ("fits_delete_key (fp, key)");

   do_write_xxx (&_fits_delete_key, _NARGS);
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
     usage ("fits_update_logical (fp, key, value, comment)");

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
   fits_check_error (_fits_get_num_keys (fp, &nkeys));

   variable recs = String_Type [nkeys];
   _for (0, nkeys-1, 1)
     {
	variable i = ();
	variable rec;

	fits_check_error (_fits_read_record (fp, i+1, &rec));
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
     usage ("fits_write_records (fp, records[])");

   variable fp, records;
   (fp, records) = ();

   variable needs_close;
   fp = get_open_write_fp (fp, "w", &needs_close);

   if (String_Type == typeof (records))
     records = [records];

   foreach (records)
     {
	variable rec = ();
	fits_check_error (_fits_write_record (fp, rec));
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
   variable bitpix = [10, 8, 16, 20, 32, 40, -32, -64];

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

   fits_check_error (_fits_read_img (fp, &a));
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

   fits_check_error (_fits_create_img (fp, fits_get_bitpix (type), dims));
   if (extname != NULL)
     fits_check_error (_fits_update_key (fp, "EXTNAME", extname, NULL));

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
	     fits_check_error (_fits_update_key (fp, keyword, val, NULL), keyword);
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
		       fits_check_error (_fits_write_history (fp, val));
		       continue;
		    }
		  if (keyword == "comment")
		    {
		       fits_check_error (_fits_write_comment (fp, val));
		       continue;
		    }
		  vmessage ("*** WARNING: history/comment record name '%s' is not supported",
			    history);
	       }
	  }
     }

   fits_check_error (_fits_write_img (fp, image));
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
   fits_check_error (_fits_write_img (fp, data));
}

define fits_iterate ()
{
   if (_NARGS != 4)
     {
	usage ("\
fits_iterate (fp, {col1, col2, ...}, &func, {arg1, arg2,...})\n\
  This function iterates over all rows in a binary table calling func,\n\
  which should be declared as:\n\
    define func(arg1, arg2, ..., data1, data2, ...)\n\
  where data1 is read form col1, etc.  The function must return 1 for\n\
  processing to continue.  Any other value will cause iteration to stop.\n\
\n\
  Qualifiers: drows=VAL\n\
    Use VAL rows for the number of rows to read at one time (default=4096)\n\
"
	      );
     }

   variable fp, col_list, func, func_list;
   (fp, col_list, func, func_list)=();

   variable delta_rows = qualifier ("drows", 4096);
   if (delta_rows <= 0)
     throw InvalidParmError, "drows must be >= 1";

   variable fpinfo = open_read_cols (fp, col_list;; __qualifiers);
   variable num_rows = fpinfo.num_rows;
   variable num_cols = fpinfo.num_cols;
   variable i;

   delta_rows--;
   variable r0 = 1;
   while (r0 <= num_rows)
     {
	variable r1 = r0 + delta_rows;
	if (r1 > num_rows)
	  r1 = num_rows;

	if (1 != (@func)(__push_list(func_list), read_cols (fpinfo, r0, r1)))
	  break;

	r0 = r1 + 1;
     }
}

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
