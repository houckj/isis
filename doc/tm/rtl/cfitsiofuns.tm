
\function{_fits_get_errstatus}
\synopsis{Retrieve a text string corresponding to an error code}
\usage{String_Type _fits_get_errstatus (Int_Type status)}
\description
  \xreferences{fits_get_errstatus}
\seealso{fits_read_errmsgs, _fits_read_errmsg}
\done

\function{_fits_read_errmsg}
\synopsis{Retrieve an error message from the cfitsio error message stack}
\usage{String_Type _fits_read_errmsg ()}
\description
  \xreferences{fits_read_errmsg}
\notes
  This function returns \NULL if there are no error messages available.
\seealso{fits_read_errmsgs, fits_set_verbose_errors}
\done

\function{_fits_open_file}
\synopsis{Open a fits file}
\usage{status = _fits_open_file (Ref_Type fptr, String_Type file, String_Type mode)}
\description
  \xreferences{fits_open_file}

  The main difference between this function and its cfitsio
  counterpart is that the \exmp{mode} argument is a string whose value
  must be one of the following:
#v+
    "r"     Open the file in read-only mode
    "w"     Open the file in read-write mode
    "c"     Create a new file.
#v-
  Note that if the mode argument is \exmp{"c"}, then the
  cfitsio \exmp{fits_create_file} function will be called.  An
  important difference between this intrinsic function and the
  \exmp{fits_create_file} library function is that if the file already
  exists, the library function will return an error, whereas
  \ifun{_fits_open_file} will simply delete the file before creating a
  new one.
\done

\function{_fits_delete_file}
\synopsis{Delete the file associated with a Fits_File_Type object}
\usage{status = _fits_delete_file (Fits_File_Type fptr)}
\description
  \xreferences{fits_delete_file}
\done

\function{_fits_close_file}
\synopsis{Close a fits file}
\usage{status = _fits_close_file (Fits_File_Type fptr)}
\description
  \xreferences{fits_close_file}
\done

\function{_fits_movabs_hdu}
\synopsis{Move to an absolute HDU number}
\usage{status = _fits_movabs_hdu (Fits_File_Type fptr, Int_Type hdunum)}
\description
  \xreferences{fits_movabs_hdu}
\notes
  The cfitsio counterpart also returns the HDU type of the resulting
  HDU.
\done

\function{_fits_movrel_hdu}
\synopsis{Move a relative number of HDUs}
\usage{status = _fits_movrel_hdu (Fits_File_Type fptr, Int_Type nmove)}
\description
  \xreferences{fits_movrel_hdu}
\notes
  The cfitsio counterpart also returns the HDU type of the resulting
  HDU.
\done

\function{_fits_movnam_hdu}
\synopsis{Move to a named HDU}
\usage{status = _fits_movnam_hdu (fptr, hdutype, extname, extvers)}
#v+
   Fits_File_Type fptr;
   Int_Type hdutype, extvers;
   String_Type extname;
#v-
\description
  \xreferences{fits_movnam_hdu}
\done

\function{_fits_get_num_hdus}
\synopsis{Return the number of HDUs in a FITS file}
\usage{status = _fits_get_num_hdus (Fits_File_Type fptr, Ref_Type hdunum)}
\description
  \xreferences{fits_get_num_hdus}
\done

\function{_fits_get_hdu_num}
\synopsis{Return the current HDU number}
\usage{hdunum = _fits_get_hdu_num (Fits_File_Type fptr)}
\description
  \xreferences{fits_get_hdu_num}
\done

\function{_fits_get_hdu_type}
\synopsis{Get the current HDU type}
\usage{status = _fits_get_hdu_type (Fits_File_Type fptr, Ref_Type hdutype)}
\description
  \xreferences{fits_get_hdu_type}

  Upon a sucessful return, the value of the variable associated with
  the \exmp{hdutype} reference will be set to one of the following
  constants:
#v+
   _FITS_IMAGE_HDU
   _FITS_ASCII_TBL
   _FITS_BINARY_TBL
#v-
\done

\function{_fits_copy_file}
\synopsis{Copy a fits file}
\usage{status = _fits_copy_file (infptr, outfptr, previous, current, following)}
#v+
   Fits_File_Type infptr, outfptr;
   Int_Type previous, current, following;
#v-
\description
  \xreferences{fits_copy_file}
\done

\function{_fits_copy_hdu}
\synopsis{Copy an HDU}
\usage{status = _fits_copy_hdu (infptr, outfptr, morekeys)}
#v+
   Fits_File_Type infptr, outfptr;
   Int_Type morekeys;
#v-
\description
  \xreferences{fits_copy_hdu}
\done

\function{_fits_copy_header}
\synopsis{Copy a fits header from one HDU to another}
\usage{status = _fits_copy_header (Fits_File_Type infptr, Fits_File_Type outfptr)}
\description
  \xreferences{fits_copy_header}
\done

\function{_fits_delete_hdu}
\synopsis{Delete the current HDU}
\usage{status = _fits_delete_hdu (Fits_File_Type fptr)}
\description
  \xreferences{fits_delete_hdu}
\notes
  The corresponding cfitsio function also returns the HDU type of the
  new HDU.  If that information is necessary, make a separate call to
  \ifun{_fits_get_hdu_type}.
\done

\function{_fits_create_img}
\synopsis{Create a new image extension}
\usage{status = _fits_create_img (fptr, bitpix, dims)}
#v+
   Fits_File_Type fptr;
   Int_Type bitpix;
   Array_Type dims;
#v-
\description
  \xreferences{fits_create_img}
\notes
  This function differs from the corresponding cfitsio function in
  that the \exmp{dims} array is assumed to be a 1-d integer array
  whose elements specify the number of axes and the size of each axis.
  In particular, the value of the NAXIS keyword will be given by
  \exmp{length(dims)}, the value of NAXIS1 will be given by
  \exmp{dims[0]}, and so forth.
\done

\function{_fits_write_img}
\synopsis{Write an image}
\usage{status = _fits_write_img (Fits_File_Type fptr, Array_Type img}
\description
  \xreferences{fits_write_img}
\notes
  This function differs from its cfitsio counterpart in that the whole
  image represented by the array \exmp{img} will be written out.
\done

\function{_fits_read_img}
\synopsis{Read an image}
\usage{status = _fits_read_img (Fits_File_Type fptr, Ref_Type img)}
\description
  \xreferences{fits_read_img}
\notes
  This function differs from the corresponding cfitsio routine in that
  the type of the image returned is automatically determined by the
  routine via a call to \exmp{fits_get_img_type}.  The dimensionality
  of the returned image is given by \exmp{[...,NAXIS2,NAXIS1]} such
  that \exmp{NAXIS1} corresponds to the fastest varying dimension,
  \exmp{NAXIS2} the next fastest varying, etc.
\done

\function{_fits_create_binary_tbl}
\synopsis{Create a binary table extension}
\usage{status = _fits_create_binary_tbl (fptr, naxis2, ttype, tform, tunit, extname)}
#v+
   Fits_File_Type fptr;
   Array_Type tunit, tform, ttype;
   Int_Type naxis2;
   String_Type extname;
#v-
\notes
  The \exmp{_fits_create_binary_tbl} function is a wrapper around the
  \cfitsioxref{fits_create_tbl} function explicitly creating a binary table.
  The input arrays \exmp{ttype}, \exmp{tform}, and \exmp{tunit} must
  be the same length, which determines the number of columns in the
  table.  The \exmp{tunit} and \exmp{extname} parameters may be set to
  \NULL.
\done

\function{_fits_update_key}
\synopsis{Update a keyword or append a new one}
\usage{status = _fits_update_key (fptr, keyname, value, comment)}
#v+
   Fits_File_Type fptr;
   String_Type keyname, comment;
   Any_Type value;
#v-
\description
  \xreferences{fits_update_key}
\notes
  The data-type for the \exmp{value} argument must be an approriate type for
  FITS keywords.  If \exmp{value} is a string, then the string will be
  written as a cfitsio long-string using the
  \cfitsioxref{fits_update_key_longstr} function.  If \exmp{value} is \NULL,
  then the \cfitsioxref{fits_update_key_null} function will be called to
  update the keyword.

  The \exmp{comment} parameter may be set to \NULL to if the comment
  field associated with the keyword is not to be modified.

  To update the value of a boolean (logical) keyword, use the
  \ifun{_fits_update_logical} function.
\done

\function{_fits_update_logical}
\synopsis{Update the value of a boolean keyword}
\usage{status = _fits_update_logical (fptr, keyname, value, comment)}
#v+
  Fits_File_Type fptr;
  String_Type keyname, comment;
  Int_Type value;
#v-
\description
  The \ifun{_fits_update_logical} function is a wrapper around the
  cfitsio \cfitsioxref{fits_update_key} function with TLOGICAL specified as
  the datatype argument.  If the \exmp{value} parameter is non-zero,
  then a value \exmp{T} (TRUE) will be given to the specified
  keyword.  Otherwise, the value of the keyword will be set to
  \exmp{F} (FALSE).
  
  If the \exmp{comment} parameter is NULL, then the keyword's comment
  field will not be modified.
\seealso{_fits_update_key}
\done

\function{_fits_write_comment}
\synopsis{Write a COMMENT header}
\usage{status = _fits_write_comment (Fits_File_Type fptr, String_Type comment)}
\description
  \xreferences{fits_write_comment}
\done

\function{_fits_write_history}
\synopsis{Write a HISTORY header}
\usage{status = _fits_write_history (Fits_File_Type fptr, String_Type history)}
\description
  \xreferences{fits_write_history}
\done

\function{_fits_write_date}
\synopsis{Write a DATE keyword}
\usage{status = _fits_write_date (Fits_File_Type fptr)}
\description
  \xreferences{fits_write_date}
\done

\function{_fits_write_record}
\synopsis{Write a keyword record}
\usage{status = _fits_write_record (Fits_File_Type fptr, String_Type card)}
\description
  \xreferences{fits_write_record}
\done

\function{_fits_modify_name}
\synopsis{Rename a keyword}
\usage{status = _fits_modify_name (fptr, oldname, newname)}
#v+
   Fits_File_Type fptr;
   String_Type oldname, newname;
#v-
\description
  \xreferences{fits_modify_name}
\done

\function{_fits_get_num_keys}
\synopsis{Get the number of keywords in the current HDU}
\usage{status _fits_get_num_keys (Fits_File_Type fptr, Ref_Type numkeys)}
\description
  This function is a wrapper around the cfitsio
  \cfitsioxref{fits_get_hdrspace} function.  It obtains the number of
  existing keywords in the current HDU (excluding the END keyword) and
  assigns that value to variable associated with the \exmp{numkeys}
  parameter.
\done

\function{_fits_read_key_integer}
\synopsis{Read the value of a keyword as an integer}
\usage{status = _fits_read_key_integer (fptr, keyname, value, comment)}
#v+
   Fits_File_Type fptr;
   String_Type keyname;
   Ref_Type value, comment;
#v-
\description
  This function uses the cfitsio \cfitsioxref{fits_read_key} function to read
  the value of the specifed keyword as an integer.  Its value is
  assigned to the variable referenced by the \exmp{value} parameter.
  If the comment parameter is non-NULL, then the value of the comment
  field will be assigned to it.
\seealso{_fits_read_key, _fits_read_key_string, _fits_read_key_double}
\done

\function{_fits_read_key_string}
\synopsis{Read the value of a keyword as a string}
\usage{status = _fits_read_key_string (fptr, keyname, value, comment)}
#v+
   Fits_File_Type fptr;
   String_Type keyname;
   Ref_Type value, comment;
#v-
\description
  This function uses the cfitsio \cfitsioxref{fits_read_key_longstr} function
  to read the value of the specifed keyword as a cfitsio long-string.
  The string is assigned to the variable referenced by the
  \exmp{value} parameter. If the comment parameter is non-NULL, then
  the value of the comment field will be assigned to it.
\seealso{_fits_read_key, _fits_read_key_integer, _fits_read_key_double}
\done

\function{_fits_read_key_double}
\synopsis{Read the value of a keyword as a double}
\usage{status = _fits_read_key_double (fptr, keyname, value, comment)}
#v+
   Fits_File_Type fptr;
   String_Type keyname;
   Ref_Type value, comment;
#v-
\description
  This function uses the cfitsio \cfitsioxref{fits_read_key} function
  to read the value of the specifed keyword as a double.
  The keyword's value is assigned to the variable referenced by the
  \exmp{value} parameter. If the comment parameter is non-NULL, then
  the value of the comment field will be assigned to it.
\seealso{_fits_read_key, _fits_read_key_integer, _fits_read_key_string}
\done

\function{_fits_read_key}
\synopsis{Read the value of a keyword}
\usage{status = _fits_read_key (fptr, keyname, value, comment)}
#v+
   Fits_File_Type fptr;
   String_Type keyname;
   Ref_Type value, comment;
#v-
\description
  This function uses the cfitsio \cfitsioxref{fits_read_key} function
  to read the value of the specifed keyword.  It first uses the cfitsio
  \cfitsioxref{fits_get_keytype} function to determine the data-type for the
  keyword and then calls \cfitsioxref{fits_read_key} using that data-type.
  The resulting value is assigned to the variable referenced by the
  \exmp{value} parameter. If the comment parameter is non-NULL, then
  the value of the comment field will be assigned to it.
\seealso{_fits_read_key_integer, _fits_read_key_string, _fits_read_key_double}
\done

\function{_fits_read_record}
\synopsis{Read a specified record from the current HDU}
\usage{status = _fits_read_record (fptr, keynum, card)}
#v+
   Fits_File_Type fptr;
   Int_Type keynum;
   Ref_Type card;
#v-
\description
  \xreferences{fits_read_record}
\done

\function{_fits_delete_key}
\synopsis{Delete a keyword from the header}
\usage{status = _fits_delete_key (Fits_File_Type fptr, String_Type keyname)}
\description
  \xreferences{fits_delete_key}
\done

\function{_fits_get_colnum}
\synopsis{Get the column number of a specfied table column}
\usage{status = _fits_get_colnum (fptr, colname, colnum)}
#v+
   Fits_File_Type fptr;
   String_Type colname;
   Ref_Type colnum;
#v-
\description
  \xreferences{fits_get_colnum}

\notes
  The corresponding cfitsio function permits a wildcard match to the
  \exmp{colname} parameter.  The current wrapping of this function
  does not support such matching.
  
  The \exmp{colname} parameter is treating in a case-insensitive manner.
\done

\function{_fits_insert_rows}
\synopsis{Insert rows into a table}
\usage{status = _fits_insert_rows (fptr, firstrow, nrows)}
#v+
   Fits_File_Type fptr;
   Int_Type firstrow, nrows;
#v-
\description
  \xreferences{fits_insert_rows}
\done

\function{_fits_delete_rows}
\synopsis{Delete rows from a table}
\usage{status = _fits_delete_rows (fptr, firstrow, nrows)}
#v+
   Fits_File_Type fptr;
   Int_Type firstrow, nrows;
#v-
\description
  \xreferences{fits_delete_rows}
\done

\function{_fits_insert_cols}
\synopsis{Insert columns into a table}
\usage{status = _fits_insert_cols (fptr, colnum, ttype, tform)}
#v+
   Fits_File_Type fptr;
   Int_Type colnum;
   Array_Type ttype, tform;
#v-
\description
  \xreferences{fits_insert_cols}
\notes
   The number of columns to be inserted is given by the length of
   the \exmp{ttype} and \exmp{tform} arrays, which must be of the same
   length.
\done

\function{_fits_delete_col}
\synopsis{Delete a column from a table}
\usage{status = _fits_delete_col (Fits_File_Type fptr, Int_Type colnum)}
\description
  \xreferences{fits_delete_col}
\done

\function{_fits_get_num_cols}
\synopsis{Get the number of table columns}
\usage{status = _fits_get_num_cols (Fits_File_Type fptr, Ref_Type ncols)}
\description
  \xreferences{fits_get_num_cols}
\done

\function{_fits_get_rowsize}
\synopsis{Get the number of rows to read or write for maximum efficiency}
\usage{status = _fits_get_rowsize (Fits_File_Type fptr, Ref_Type nrows)}
\description
  \xreferences{fits_get_rowsize}
\done

\function{_fits_get_num_rows}
\synopsis{Get the number of table rows}
\usage{status = _fits_get_num_cols (Fits_File_Type fptr, Ref_Type nrows)}
\description
  \xreferences{fits_get_num_rows}
\done

\function{_fits_write_col}
\synopsis{Write data to a table column}
\usage{status = _fits_write_col (fptr, colnum, firstrow, firstelem, array)}
#v+
   Fits_File_Type fptr;
   Int_Type colnum;
   Int_Type firstrow, firstelem;
   Array_Type array;
#v-
\description
  \xreferences{fits_write_col}
\notes
   The number of elements written out to the column by this function
   will be equal to the number of elements in the array.
\done

\function{_fits_read_col}
\synopsis{Read elements from a column}
\usage{status = _fits_read_col (fptr, colnum, firstrow, numrows, array}
#v+
   Fits_File_Type fptr;
   Int_Type colnum, firstrow, numrows;
   Ref_Type array;
#v-
\description
  This function is a complicated wrapper around a number of cfitsio
  functions to enable it to read any type of column, including vector
  and variable length columns.  The data read by the function is
  assigned as the appropriately typed array to the variable referenced
  by the \exmp{array} argument.

  For ordinary scalar columns, a 1-d array of size \exmp{numrows} will
  be produced.  For a vector column, a 2d array of size
  \exmp{[numrows,repeat]} will be generated.  Here \exmp{repeat} is
  given by the repeat value associated with the column.  For a
  variable length column, where data are stored in the heap of the
  HDU, the data will be read as a 1-d array of \exmp{numrows} arrays.
  
  If the column is a bit-valued column, then data will be returned as
  an array of integers of the appropriate size.  Currently only 8X,
  16X, and 32X bit columns are supported.
\seealso{_fits_read_cols, _fits_write_col}
\done

\function{_fits_get_keytype}
\synopsis{Get a keyword's data type}
\usage{status = _fits_get_keytype (fptr, keyname, type)}
#v+
   Fits_File_Type fptr;
   String_Type keyname;
   Ref_Type type;
#v-
\description
  \xreferences{fits_get_keytype}
\notes
  This function differs from its cfitsio counterpart in that instead
  of explicitly specifying the keyword's value string, this function
  uses the value of the specified keyword.  It also returns the type
  as a \slang \var{DataType_Type} object, e.g., \exmp{Int_Type},
  \exmp{Complex_Type}, etc.
\done

\function{_fits_get_keyclass}
\synopsis{Get the class of an input header record}
\usage{Int_Type _fits_get_keyclass (String_Type card)}
\description
  \xreferences{fits_get_keyclass}
\done

\function{_fits_read_cols}
\synopsis{Read one or more table columns}
\usage{status = _fits_read_cols (fptr, colnums, firstrow, nrows, arrays)}
#v+
   Fits_File_Type fptr;
   Array_Type colnums;
   Int_Type firstrow, numrows;
   Ref_Type arrays;
#v-
\description
  This function performs a similar task as the \exmp{_fits_read_col}.
  The main difference is that instead of reading a single column, it
  is capable of reading multiple columns specified by the \exmp{colnums}
  parameter.  It assigns the data as an array of arrays to the
  variable referenced by the \exmp{arrays} parameter.  See the
  documentation for the \ifun{_fits_read_col} function for more
  information.
\notes
  This function takes advantage of the cfitsio buffering mechanism to
  optimize the reads.
\done

\function{_fits_write_chksum}
\synopsis{Compute and write DATASUM and CHECKSUM keywords}
\usage{status = _fits_write_chksum (Fits_File_Type fptr)}
\description
  \xreferences{fits_write_chksum}
\done

\function{_fits_update_chksum}
\synopsis{Update the CHECKSUM keyword}
\usage{status = _fits_update_chksum (Fits_File_Type fptr)}
\description
  \xreferences{fits_update_chksum}
\done

\function{_fits_verify_chksum}
\synopsis{Verify the checksums for the current HDU}
\usage{status = _fits_verify_chksum (fptr, dataok, hduok)}
#v+
   Fits_File_Type fptr;
   Ref_Type dataok, hduok;
#v-
\description
  \xreferences{fits_verify_chksum}
\done

\function{_fits_get_chksum}
\synopsis{Get the checksums for the current HDU}
\usage{status = _fits_get_chksum (fptr, datasum, hdusum)}
#v+
   Fits_File_Type fptr;
   Ref_Type datasum, hdusum;
#v-
\description
  \xreferences{fits_get_chksum}
\done

\function{_fits_get_version}
\synopsis{Get the cfitsio library version number}
\usage{Float_Type _fits_get_version ()}
\description
  \xreferences{fits_get_version}
\done
