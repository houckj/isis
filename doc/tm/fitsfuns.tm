#c __FILE__: ../../src/fits.sl
#c __LINE__: 44
\function{fits_open_file}
\synopsis{Open a fits file}
\usage{Fits_File_Type fits_open_file (String_Type filename, String_Type mode)}
\description
  The \var{fits_open_file} function can be used to open and existing fits
  file for reading or updating, or to create a new fits file, depending upon
  the value of the \var{mode} parameter.  Specifically, if \var{mode} is 
  \exmp{"r"}, the file will be opened for reading.  If \var{mode} is \exmp{"w"},
  the file will be opened for updating (both reading and writing).  Otherwise, 
  \var{mode} must be \var{"c"}, which indicates that a new file is to be created.
  In the latter case, if a file already exists with the specified name, it will
  get deleted and a new one created in its place.
 
  If the function fails, it will signal an error; otherwise an open file
  pointer will be returned.
\seealso{fits_close_file, fits_create_binary_table}
\done
#c __LINE__: 79
\function{fits_close_file}
\synopsis{Close a fits file}
\usage{fits_close_file (Fits_File_Type f)}
\description
  The \var{fits_close_file} closes a previously opened fits file.  The function
  will signal an error if the operation fails.
\notes
  This function could fail if it fails to write out any buffered data because
  of filesystem errors (disk full, etc.).
\seealso{fits_open_file}
\done
#c __LINE__: 179
\function{fits_move_to_interesting_hdu}
\synopsis{Move to an extension that looks interesting}
\usage{fits_move_to_interesting_hdu (fp [, hdu_type]}
#v+
  Fits_File_Type fp;
  Int_Type hdu_type;
#v-
\description
  The function move the fits file pointer \var{fp} forward to an HDU that looks 
  interesting.  By definition, an interesting HDU is one in which NAXIS is 
  non-zero.  The first parameter \var{fp} must be a pointer to an already open
  fits file.  The second parameter, if present, may be used to specifiy the 
  type of HDU, e.g., either an image (\exmp{hdu_type=_FITS_IMAGE_HDU}) or a 
  binary table (\exmp{hdu_type=_FITS_BINARY_TBL}).
 
  If the function fails to find an interesting HDU of the appropriate type, 
  an exception will be generated.
\seealso{fits_open_file}
\done
#c __LINE__: 247
\function{fits_key_exists}
\synopsis{Check for the existence of a keyword}
\usage{Int_Type fits_key_exists (fd, key)}
#v+
   Fits_File_Type or String_Type fd;
   String_Type key;
#v-
\description
  The \var{fits_key_exists} function checks for the existence of a specified 
  keyword in the file specified by the descriptor \var{fd}, which must specify
  the name of a file or an open file pointer.
 
  If the specified key exists, the function return \1, otherwise it returns \0.
\seealso{fits_read_key, fits_read_header}
\done
#c __LINE__: 305
\function{fits_get_colnum}
\synopsis{Get the column numbers of specified columns}
\usage{column_num = fits_get_colnum (fd, column_name)}
#v+
   Fits_File_Type or String_Type fd;
   String_Type column_name;
#v-
\description
  This function returns the column number of the column with the specified name.
  The file-descriptor \exmp{fd} must specify the name of a file, or an open
  fits file pointer.
\seealso{fits_binary_table_column_exists}
\done
#c __LINE__: 337
\function{fits_binary_table_column_exists}
\synopsis{Check for the existence of a binary table column}
\usage{Int_Type fits_binary_table_column_exists (fd, col)}
#v+
   Fits_File_Type or String_Type fd;
   String_Type col;
#v-
\description
  This function may be used to determine whether or not a named column
  exists in a binary table.  The table is specified via the \var{fd} 
  parameter which must either be the name of a file containing the binary
  table, or an file pointer.
 
  If the specified column exists, \1 will be returned; otherwise the function
  will return \0.
\seealso{fits_key_exists, fits_open_file}
\done
#c __LINE__: 523
\function{fits_read_col}
\synopsis{Read one or more columns from a FITS binary table}
\usage{(x1, ...xN) = fits_read_col (file, c1, ... cN)}
v+
   Fits_File_Type or String_Type file;
   Int_Type or String_Type c1, ...cN;
v-
\description
  This function returns one or more vectors containing objects in the
  specified columns of the binary table indicated by \var{file}.  If
  \var{file} is a string, then the file will be opened via the virtual
  file specification implied by \var{file}. Otherwise, \var{file}
  should represent an already opened FITS file.  The column parameters
  may either be strings denoting the column names, or integers
  representing the column numbers.
\seealso{fits_read_cell, fits_read_row, fits_read_table}
\done
#c __LINE__: 561
\function{fits_read_col_struct}
\synopsis{Read one or more columns from a FITS binary table}
\usage{struct = fits_read_col_struct (file, col1, ...)}
#v+
    Fits_File_Type or String_Type file;
    String_Type col1, ...;
#v-
\description
  This function works exactly like \var{fits_read_col} except it returns the
  values in a structure.  See the documentation on that function for more
  information.
 
\seealso{fits_read_col, fits_read_key_struct, fits_read_row, fits_read_header}
\done
#c __LINE__: 589
\function{fits_read_cell}
\synopsis{Read a cell from a FITS binary table}
\usage{X = fits_read_cell (file, c, r)}
v+
   Fits_File_Type or String_Type file;
   Int_Type r, c;
v-
\description
  This function returns the object in the column \var{c} and row
  \var{r} of the binary table indicated by \var{file}.  If \var{file}
  is a string, then the file will be opened via the virtual file
  specification implied by \var{file}. Otherwise, \var{file} should
  represent an already opened FITS file.
\seealso{fits_read_col, fits_read_row}
\done
#c __LINE__: 645
\function{fits_read_row}
\synopsis{Read a row from a FITS binary table}
\usage{Struct_Type fits_read_cell (file, r)}
v+
   Fits_File_Type or String_Type file;
   Int_Type r;
v-
\description
  This function returns a structure containing the data in the columns
  of the row \var{r} of the binary table indicated by \var{file}. If
  \var{file} is a string, then the file will be opened via the virtual
  file specification implied by \var{file}. Otherwise, \var{file}
  should represent an already opened FITS file.
\seealso{fits_read_col, fits_read_cell}
\done
#c __LINE__: 667
\function{fits_read_header}
\synopsis{Read a FITS header}
\usage{Struct_Type fits_read_header (file)}
#v+
    Fits_File_Type or String_Type file;
#v-
\description
  This function reads the header of the fits file given by the
  \var{file} parameter and returns it as a structure.  If \var{file} is
  a string, then the file will be opened via the virtual file
  specification implied by \var{file}. Otherwise, \var{file} should
  represent an already opened FITS file.
\seealso{fits_read_table}
\done
#c __LINE__: 696
\function{fits_read_table}
\synopsis{Read a FITS table}
\usage{Struct_Type fits_read_table (file [,columns...])}
#v+
    Fits_File_Type or String_Type file;
#v-
\description
  \var{fits_read_table} reads the data in a table of the FITS file
  specified by \var{file} and returns it as a structure.  If the optional
  column name parameters are specified, then only those columns will be read.
  Otherwise, the entire table will be returned.
 
  If \var{file} is a string, then the file will be opened via the virtual file
  specification implied by \var{file}. Otherwise, \var{file} should
  represent an already opened FITS file.
\seealso{fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
\done
#c __LINE__: 773
\function{fits_read_key}
\synopsis{Read one or more keywords from a FITS file}
\usage{(val1,...) = fits_read_key (file, key1, ...)}
#v+
    Fits_File_Type or String_Type file;
    String_Type key1, ...;
#v-
\description
  \var{fits_read_key} reads the values of one or more keywords in the fits
  file specified by \var{file} and returns them.  If \var{file}
  is a string, then the file will be opened via the virtual file
  specification implied by \var{file}. Otherwise, \var{file} should
  represent an already opened FITS file.  If any of the keywords do not exist,
  a value of \NULL will be returned for the corresponding keyword.
\seealso{fits_read_key_struct, fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
\done
#c __LINE__: 819
\function{fits_read_key_struct}
\synopsis{Read one or more keywords from a FITS file}
\usage{struct = fits_read_key (file, key1, ...)}
#v+
    Fits_File_Type or String_Type file;
    String_Type key1, ...;
#v-
\description
  This function works exactly like \var{fits_read_key} excepts returns the
  values in a structure.  See the documentation on that function for more
  information.
\seealso{fits_read_key, fits_read_col, fits_read_cell, fits_read_row, fits_read_header}
\done
#c __LINE__: 860
\function{fits_create_binary_table}
\synopsis{Prepare a binary table}
\usage{fits_create_binary_table (file, extname, nrows, ttype, tform, tunit)}
#v+
    Fits_File_Type or String_Type file;
    String_Type extname;
    Int_Type nrows;
    String_Type ttype[];
    String_Type tform[];
    String_Type tunit[];
#v-
\description
  This creates a new binary table with the specified structure.  The parameters
  \var{ttype}, \var{tform}, and \var{tunit} are string arrays that specify
  the column names, column data type, and column units, respectively.
  The binary table will be given the extension name \var{extname}.
\seealso{fits_write_binary_table, fits_open_file}
\done
#c __LINE__: 895
\function{fits_write_binary_table}
\synopsis{Write a binary table}
\usage{fits_write_binary_table (file, extname, sdata, [skeys [,hist]])}
#v+
    Fits_File_Type or String_Type file;
    String_Type extname;
    Struct_Type sdata;
    Struct_Type skeys;
    Struct_Type hist;
#v-
\description
  The \var{fits_write_binary_table} function creates a new binary table in
  the specified file.  The parameter \var{file} specifies either a filename or
  an open file pointer.  The \var{extname} parameter specifies the extension
  name of the binary table.  The data written to table are specified in the 
  \var{sdata} structure, where the name of the structure field specifies the 
  column name.  If \var{skeys} is non-NULL, then it is a structure indicating
  additional keywords to be written to the header of the binary table.  If the
  optional parameter \var{hist} is present and non-NULL, then it is a structure
  whose fields indicate either comment or history information to be written
  to the header.
\example
  The following code
#v+
    variable data = struct { x, cosx, sinx };
    data.x = [0:2*PI:0.01];
    data.cosx = cos(data.x);
    data.sinx = sin(data.x);

    variable keys = struct { hduname, username};
    keys.hduname = "COSXSINX";
    keys.username = "John Doe";

    variable hist = struct { history, comment};
    hist.history = ["This is a history record", "This is another"];
    hist.comment = ["This is a comment", "And this is another"];

    fits_write_binary_table ("foo.fits", "COSXSINX", data, keys, hist);
#v-
 produces a binary table with the header:
#v+
    XTENSION= 'BINTABLE' / binary table extension
    BITPIX  =                   8 / 8-bit bytes
    NAXIS   =                   2 / 2-dimensional binary table
    NAXIS1  =                  24 / width of table in bytes
    NAXIS2  =                 629 / number of rows in table
    PCOUNT  =                   0 / size of special data area
    GCOUNT  =                   1 / one data group (required keyword)
    TFIELDS =                   3 / number of fields in each row
    TTYPE1  = 'x       ' / label for field   1
    TFORM1  = 'D       ' / data format of field: 8-byte DOUBLE
    TTYPE2  = 'cosx    ' / label for field   2
    TFORM2  = 'D       ' / data format of field: 8-byte DOUBLE
    TTYPE3  = 'sinx    ' / label for field   3
    TFORM3  = 'D       ' / data format of field: 8-byte DOUBLE
    EXTNAME = 'COSXSINX' / name of this binary table extension
    HDUNAME = 'COSXSINX'
    USERNAME= 'John Doe'
    HISTORY This is a history record
    HISTORY This is another
    COMMENT This is a comment
    COMMENT And this is another
#v-
\notes
  This function provides no mechanism to mix comments and keyword records.  As
  the example shows, this function places the comment and history records at
  the end of the table.
\seealso{fits_create_binary_table, fits_open_file}
\done
#c __LINE__: 1189
\function{fits_update_key}
\synopsis{Update the value of a keyword}
\usage{fits_update_key (fd, key, val [,comment])}
#v+
    String_Type or Fits_File_Type fd;
    String_Type key;
    Any type val;
    String_Type comment;
#v-
\description
  The \var{fits_update_key} function updates the value and comment fields
  of an existing keyword with the specified name.  If the keyword does not 
  exist, a new keyword will be appended to the end of the header.
\seealso{fits_update_logical, fits_read_key}
\done
#c __LINE__: 1220
\function{fits_update_logical}
\synopsis{Update the value of a logical (boolean) keyword}
\usage{fits_update_logical (fd, key, val, comment)}
#v+
    String_Type or Fits_File_Type fd;
    String_Type key;
    Any type val;
    String_Type comment;
#v-
\description
  The \var{fits_update_logical} function updates the value and comment fields
  of an existing keyword of the specified name with the specified boolean value.
  If the keyword does not exist, a new keyword will be appended to the end of 
  the header.
\seealso{fits_update_key}
\done
#c __LINE__: 1245
\function{fits_write_comment}
\synopsis{Write a comment to the header}
\usage{fits_write_comment (fd, comment)}
#v+
  Fits_File_Type or String_Type fd;
  String_Type comment;
#v-
\description
  As the name indicates, this function writes a comment record to the specified
  fits file.  The file-descriptor \exmp{fd} must either be the name of a fits
  file or an open fits file pointer.
\seealso{fits_update_key, fits_write_history}
\done
#c __LINE__: 1267
\function{fits_write_history}
\synopsis{Write a history record to the header}
\usage{fits_write_history (fd, history)}
#v+
  Fits_File_Type or String_Type fd;
  String_Type history;
#v-
\description
  As the name indicates, this function writes a history record to the specified
  fits file.  The file-descriptor \exmp{fd} must either be the name of a fits
  file or an open fits file pointer.
\seealso{fits_update_key, fits_write_comment}
\done
#c __LINE__: 1289
\function{fits_write_date}
\synopsis{Write the DATE keyword to the current HDU}
\usage{fits_write_date (fd)}
#v+
   Fits_File_Type or String_Type fd;
#v-
\description
  The \sfun{fits_write_date} function calls \ifun{_fits_write_date} to write
  the DATE to the header of the specified file descriptor, which  must either 
  be the name of a fits file or an open fits file pointer.
\seealso{fits_update_key}
\done
#c __LINE__: 1309
\function{fits_write_chksum}
\synopsis{Compute and write the DATASUM and CHECKSUM keywords}
\usage{fits_write_chksum (fd)}
#v+
   Fits_File_Type or String_Type fd;
#v-
\description
  The \sfun{fits_write_chksum} function calls \ifun{_fits_write_comment} to 
  compute and write the DATASUM and CHECKSUM keywords to the 
  header of the specified file descriptor, which  must either 
  be the name of a fits file or an open fits file pointer.
\seealso{fits_update_key, fits_verify_chksum}
\done
#c __LINE__: 1330
\function{fits_verify_chksum}
\synopsis{Verify the checksums for the current HDU}
\usage{isok = fits_verify_chksum (fd [,dataok, hduok])}
#v+
   Fits_File_Type or String_Type fd;
   Ref_Type dataok, hduok;
#v-
\description
  The \sfun{fits_verify_chksum} function calls \ifun{_fits_verify_chksum} to 
  verify the header and data checksums of the current HDU.  A non-zero return value
  signifies that the checksums are ok, otherwise the function returns 0 to indicate
  that the checksums are invalid.  The individual checksums of the HDU or data
  can be checked through the use of the optional parameters.
\seealso{fits_write_chksum}
\done
#c __LINE__: 1368
\function{fits_read_records}
\synopsis{Read all the records in a fits header}
\usage{String_Type[] fits_read_records (Fits_File_Type or String_Type fp)}
\description
  This function returns a list of all the header records associated with the
  fits file descriptor as an array of strings.
\seealso{fits_write_records, fits_read_key}
\done
#c __LINE__: 1400
\function{fits_write_records}
\synopsis{Write records to fits header}
\usage{fits_write_records (fd, records)}
#v+
   Fits_File_Type or String_Type fd;
   Array_Type records;
#v-
\description
  This function uses the \ifun{_fits_write_record} function to write a series
  of records to the current HDU.
\seealso{fits_read_records}
\done
#c __LINE__: 1436
\function{fits_get_keyclass}
\synopsis{Obtain the key classes for a set of cards}
\usage{Int_Type[] = fits_get_keyclass (Array_Type cards)}
\description
  This function uses the \ifun{_fits_get_keyclass} function to obtain the 
  key-classes associated with one or more cards.  The function returns an
  integer-valued array of the same length as the \exmp{cards} array.
\example
  Obtain set of header cards to those that are not associated with the cards
  describing the structure of the HDU:
#v+
    variable cards = fits_read_records ("evt2.fits[EVENTS]");
    variable classes = fits_get_keyclass (cards);
    cards = cards[where (classes != _FITS_TYP_STRUC_KEY)];
#v-
\seealso{fits_read_records, fits_read_key}
\done
#c __LINE__: 1468
\function{fits_get_bitpix}
\synopsis{Get the fits bitpix value for an array}
\usage{Int_Type fits_get_bitpix (array)}
\description
  This function may be used to obtain the bitpix value for a specified image
  array.  The array must be an integer or floating point type, otherwise
  and error will be generated.  The bitpix value is returned.
\seealso{fits_write_image_hdu, fits_read_img}
\done
#c __LINE__: 1499
\function{fits_read_img}
\synopsis{Read image data from a fits file}
\usage{Array_Type fits_read_img (fd)}
#v+
   Fits_File_Type or String_Type fd;
#v-
\description
  This function reads an image from the specified file descriptor.  
  The file descriptor must be either the name of an existing file, or an
  open file pointer.  It returns the image upon sucess, or signals an error 
  upon failure.
\seealso{fits_read_table, fits_read_col, fits_open_file, fits_write_img}
\done
#c __LINE__: 1530
\function{fits_create_image_hdu}
\synopsis{Create a primary array or image extension}
\usage{fits_create_image_hdu (fd, extname, type, dims)}
#v+
   Fits_File_Type or String_Type fd;
   String_Type extname;
   Array_Type dims;
   DataType_Type type;
#v-
\description
  This function make use of the \ifun{_fits_create_img} function to create an
  image extension or primary array of the specified type and size.  If the
  \exmp{extname} parameter is non-NULL, then an EXTNAME keyword will be 
  written out with the value of the extname parameter.
  The \exmp{dims} parameter must be a 1-d integer array that corresponds
  to the dimensions of the array to be written.
  
  If \exmp{fd} is specified as a string, then a new file of that name will be 
  created.  If a file by that name already exists, it will be deleted and
  a new one created.  If this behavior is undesired, then explicitly open the
  file and pass this routine the resulting file pointer.
\seealso{fits_write_image_hdu}
\done
#c __LINE__: 1573
\function{fits_write_image_hdu}
\synopsis{Write an image extension}
\usage{fits_write_image_hdu (file, extname, image [,skeys [,hist]])}
#v+
    Fits_File_Type or String_Type file;
    String_Type extname;
    Any_Type image
    Struct_Type skeys;
    Struct_Type hist;
#v-
\description
  The \var{fits_write_image_hdu} function creates a new image extension in
  the specified file.  The parameter \var{file} specifies either a filename or
  an open file pointer.  The \var{extname} parameter specifies the extension
  name of the image, or NULL for the primary image.  The image data written 
  to the file are specified by the \var{image} parameter.
  If the optional parameter \var{skeys} is non-NULL, then it is a 
  structure indicating additional keywords to be written to the header of the 
  binary table.  If the optional parameter \var{hist} is present and non-NULL, 
  then it is a structure whose fields indicate either comment or history 
  information to be written to the header.
\example
  The following code
#v+
    variable data = struct { x, cosx, sinx };
    data.x = [0:2*PI:0.01];
    data.cosx = cos(data.x);
    data.sinx = sin(data.x);

    variable keys = struct { hduname, username};
    keys.hduname = "COSXSINX";
    keys.username = "John Doe";

    variable hist = struct { history, comment};
    hist.history = ["This is a history record", "This is another"];
    hist.comment = ["This is a comment", "And this is another"];

    fits_write_image_hdu ("foo.fits", "COSXSINX", data, keys, hist);
#v-
 produces a binary table with the header:
#v+
    XTENSION= 'BINTABLE' / binary table extension
    BITPIX  =                   8 / 8-bit bytes
    NAXIS   =                   2 / 2-dimensional binary table
    NAXIS1  =                  24 / width of table in bytes
    NAXIS2  =                 629 / number of rows in table
    PCOUNT  =                   0 / size of special data area
    GCOUNT  =                   1 / one data group (required keyword)
    TFIELDS =                   3 / number of fields in each row
    TTYPE1  = 'x       ' / label for field   1
    TFORM1  = 'D       ' / data format of field: 8-byte DOUBLE
    TTYPE2  = 'cosx    ' / label for field   2
    TFORM2  = 'D       ' / data format of field: 8-byte DOUBLE
    TTYPE3  = 'sinx    ' / label for field   3
    TFORM3  = 'D       ' / data format of field: 8-byte DOUBLE
    EXTNAME = 'COSXSINX' / name of this binary table extension
    HDUNAME = 'COSXSINX'
    USERNAME= 'John Doe'
    HISTORY This is a history record
    HISTORY This is another
    COMMENT This is a comment
    COMMENT And this is another
#v-
\notes
  This function provides no mechanism to mix comments and keyword records.  As
  the example shows, this function places the comment and history records at
  the end of the table.
\seealso{fits_create_binary_table, fits_open_file}
\done
#c __LINE__: 1722
\function{fits_write_img}
\synopsis{Write the image data to an Image HDU}
\usage{fits_write_img (Fits_File_Type fptr, Any_Type data)}
\description
  This function writes the image data out to current HDU, assumed to be 
  an Image HDU.
\seealso{fits_write_image_hdu, fits_create_image_hdu}
\done
