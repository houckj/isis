#% -*- mode: tm; mode: fold -*-

#%{{{Macros

#i linuxdoc.tm
#d it#1 <it>$1</it>

#d slang \bf{S-lang}
#d exmp#1 \tt{$1}
#d var#1 \tt{$1}

#d ivar#1 \tt{$1}
#d ifun#1 \tt{$1}
#d cvar#1 \tt{$1}
#d cfun#1 \tt{$1}
#d svar#1 \tt{$1}
#d sfun#1 \tt{$1}
#d icon#1 \tt{$1}

#d chapter#1 <chapt>$1<p>
#d preface <preface>
#d tag#1 <tag>$1</tag>

#d function#1 \sect1{<bf>$1</bf>\label{$1}}<descrip>
#d variable#1 \sect1{<bf>$1</bf>\label{$1}}<descrip>
#d function_sect#1 \sect{$1}
#d begin_constant_sect#1 \sect{$1}<itemize>
#d constant#1 <item><tt>$1</tt>
#d end_constant_sect </itemize>

#d synopsis#1 <tag> Synopsis </tag> $1
#d keywords#1 <tag> Keywords </tag> $1
#d usage#1 <tag> Usage </tag> <tt>$1</tt>
#d description <tag> Description </tag>
#d qualifiers <tag> Qualifiers </tag>
#d qualifier#2:3 ; \tt{$1}: $2 \ifarg{$3}{(default: \tt{$3})}<newline>
#d example <tag> Example </tag>
#d notes <tag> Notes </tag>
#d seealso#1 <tag> See Also </tag> <tt>\linuxdoc_list_to_ref{$1}</tt>
#d done </descrip><p>
#d -1 <tt>-1</tt>
#d 0 <tt>0</tt>
#d 1 <tt>1</tt>
#d 2 <tt>2</tt>
#d 3 <tt>3</tt>
#d 4 <tt>4</tt>
#d 5 <tt>5</tt>
#d 6 <tt>6</tt>
#d 7 <tt>7</tt>
#d 8 <tt>8</tt>
#d 9 <tt>9</tt>
#d NULL <tt>NULL</tt>
#d file#1 <tt>$1</tt>

#d documentstyle book

#%}}}

#d module#1 \tt{$1}

\linuxdoc
\begin{\documentstyle}

\title S-Lang CFITSIO Module Reference
\author John E. Davis, \tt{davis@space.mit.edu}
\date \__today__

\toc

#d cfitsio_url http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
#d cfitsio_ref_url \
 http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/cfitsio.html
#d fits_url http://heasarc.gsfc.nasa.gov/docs/heasarc/fits.html
#d FITS FITS

#d CFITSIO CFITSIO

\chapter{Introduction} #%{{{

 \url{\fits_url}{FITS} (Flexible Image Transport System) is a data
 format that is in widespread use by the astronomical community.
 \url{\cfitsio_ref_url}{CFITSIO} is a popular C library that
 interfaces to such files and provides support for all features of the
 format.  Moreover \CFITSIO supports a number of unofficial or
 proposed \FITS conventions that are in widespread use.

 This package makes use of the \CFITSIO library allow one to
 manipulate \FITS files from the \slang interpreter.  The package
 consists of two interfaces: a high level interface and a low level
 one. The low level interface is implemented as a module and is more
 or less a straightforward wrapping of the functions of the \CFITSIO
 library.  Functions from this interface are prefixed with an
 underscore to indicate that they are part of the low-level interface.
 The high level interface is written in \slang and makes use of
 functions from the low level interface.  While there is some overlap
 with the semantics of the \CFITSIO library, the high level interface
 should be regarded as a separate interface to fits files.

 To illustrate the difference between the two interfaces, consider the
 low-level \ifun{_fits_read_col} function and its high level
 counterpart \sfun{fits_read_col}.  The low-level function reads a
 single column, specified via the column number, from a fits binary
 table and performs a minimal amount of error checking.  In contrast,
 \sfun{fits_read_col} is a high level function that reads one or more
 columns, specified either as column numbers or named columns, from a
 table and does so in a way that takes into account the way \CFITSIO
 performs buffering for maximum efficiency.  Moreover, the high level
 function checks for the presence of TDIM keywords or columns to give
 the arrays it returns the proper dimensionality.  If any errors
 occur, the function will throw an exception.

#%}}}

\chapter{The high-level interface}

\sect{Overview} #%{{{

 The high-level interface consists of a number of functions that are
 written in \slang and are designed to take some of the tedium out of
 performing standard operations on fits files.  To illustrate this
 point, consider the creation of a fits file with a binary table
 extension called called ``COSXSINX'':
#v+
    variable data = struct { x, cosx, sinx };
    data.x = [0:2*PI:0.01];
    data.cosx = cos(data.x);
    data.sinx = sin(data.x);
    fits_write_binary_table ("foo.fits", "COSXSINX", data);
#v-
 It can't get much easier than that!

#%}}}

\sect{Opening and Closing Files} #%{{{

 In general the high-level functions take an argument that represents
 the fits file to be manipulated.   It may be either an already open
 file pointer such as one returned by \sfun{fits_open_file}, or the
 name of a file to be opened.  In the documentation for the functions,
 this fact is indicated by
#v+
     Fits_File_Type or String_Type fd;
#v-
 showing that the file descriptor may be either an open file pointer
 or a string.

 If specified as a string, then a fits file of that name will be
 opened in a mode that is compatible with the operation being
 performed, with the current HDU (header-data unit) set to the first
 ``most interesting'' one.  Here, ``first most interesting''
 means the first HDU with a non-zero value for the NAXIS keyword.

 For example, sometimes one simply wants to read some keywords from a
 file.  In such a case it is not necessary to explicitly call
 \sfun{fits_open_file}.  Rather simply pass the name of the file to
 the appropriate function:
#v+
    (object, ra_targ, dec_targ)
       = fits_read_key ("foo.fits", "OBJECT", "RA_TARG", "DEC_TARG");
#v-
 It may be necessary to specify which HDU should be used if the
 ``first most interesting'' one is not the desired HDU.  The easiest way
 to do that is to specify the extension using \CFITSIO's virtual file
 syntax, e.g.,
#v+
    (object, ra_targ, dec_targ)
       = fits_read_key ("foo.fits+1", "OBJECT", "RA_TARG", "DEC_TARG");
    (object, ra_targ, dec_targ)
       = fits_read_key ("foo.fits[EVENTS]", "OBJECT", "RA_TARG", "DEC_TARG");
#v-

  If one is going to make a number of calls to functions in the
  high-level interface using the same file, then it is a good idea to
  explicitly open the file, e.g.,
#v+
    fptr = fits_open_file ("foo.fits[EVENTS]", "w");
#v-
  opens the file for both reading and writing and sets the current HDU
  to the ``EVENTS'' extension.

  The object returned by the \sfun{fits_open_file} function is an
  object of type \var{Fits_File_Type}.  It is automatically destroyed
  when it goes out of scope.  When this happens, the fits file
  attached to it will be silently closed.  Consider:
#v+
     define write_image_to_file (file, img)
     {
        variable fptr = fits_open_file ("new.fits", "c");
        fits_write_image (fptr, NULL, img, NULL, NULL);
        fits_write_date (fptr);
     }
#v-
  Here, the \exmp{write_image_to_file} function will create a new file
  called \file{new.fits} and write the specified image to the file.  In
  this example, the file pointer object was not explicitly closed.
  Since the \exmp{fptr} variable goes out of scope when the function
  returns, the cfitsio module will silently close the file.  While
  this is a convenient feature for many purposes, it is always better
  to explicitly close a file when it has been modified.  The reason
  for this is that \CFITSIO writes to internal buffers and then
  flushes those to the disk. Often some buffers will not get written
  to the disk until the file is closed.  If the disk is full, or
  something else goes wrong then the file will not be properly closed
  resulting in a incomplete or corrupt file. Hence it is strongly
  recommended that one explicitly close a file after writing to a
  file, i.e.,
#v+
     define write_image_to_file (file, img)
     {
        variable fptr = fits_open_file ("new.fits", "c");
        fits_write_image (fptr, NULL, img, NULL, NULL);
        fits_write_date (fptr);
        fits_close_file (fptr);
     }
#v-

#%}}}

\sect{Keywords}

  The high-level interface contains several functions for manipulating
  header keywords and records.
  
 \sect1{Reading Header Keywords} #%{{{

  The \sfun{fits_read_key} may be used to read the values of one or
  more keywords.  Suppose that the file \file{casA.fits} contains the
  following keywords in an extension called ``EVENTS'':
#v+
    TELESCOP= 'CHANDRA ' / Telescope
    INSTRUME= 'ACIS    ' / Instrument
    DETNAM  = 'ACIS-7  ' / Detector
    GRATING = 'NONE    ' / Grating
    OBJECT  = 'CAS A   ' / Source name
    RA_NOM  =     350.91781217089 / Nominal RA
    DEC_NOM =     58.792819089577 / Nominal Dec
    ROLL_NOM=     323.38710408328 / Nominal Roll
#v-    
  The \sfun{fits_read_key} function may be used to read, e.g, the
  OBJECT and RA_NOM keywords:
#v+
     (obj, ra) = fits_read_key ("casA.fits[EVENTS]", "OBJECT", "RA_NOM");
#v-
  After the function call, the variables \exmp{obj} and \exmp{ra} will
  have the data types \var{String_Type} and \var{Double_Type}, resp.
  If the requested keyword does not exist in the header, the function
  will return \NULL to signal its not existence:
#v+
     exptime = fits_read_key ("casA.fits[EVENTS]", "EXPTIME");
     if (exptime == NULL)
       {
          message ("*** Warning: EXPTIME does not exist.  Assuming 3.2);
          exptime = 3.2;
       }
#v-
  
  The \sfun{fits_read_key_struct} is an alternative to
  \sfun{fits_read_key} that returns a structure with field names that
  correspond to the keyword names.  In most cases, a field name will
  just be the lower case version of the keyword name.  However, if the
  keyword name does not start with an alphabetic character or contains
  a hyphen, then it will be normalized as follows:
\begin{enum}
   \item The keyword name will be lowercased.
   \item All non-alphanumeric characters will be changed to an underscore.
   \item If the first character of the resulting name is numeric, then
         the name will be prefixed with an underscore.
\end{enum}
  To illustrate the normalization process, consider:
#v+
     keys = fits_read_key_struct ("foo.fits", "OBJECT", "2CRVL3",
                                  "DATE-OBS");
#v-
  After the function call, \exmp{keys} will be a structure with the 3
  fields: \exmp{object}, \exmp{_2crvl3}, and \exmp{date_obs}.  If any
  of these keywords do not exist in the header, the value of the
  corresponding structure field will be NULL.

#%}}}

 \sect1{Writing  Header Keywords} #%{{{

  The \sfun{fits_update_key} function may be used to write or update
  the value of a keyword.   If a keyword of the specified name exists,
  then the value of the keyword will be updated to the new value.
  Otherwise a new keyword will be appended to the header.
  
  Other specialized keyword writing routines include
  \sfun{fits_write_date}, which write the current date in the required
  format, and \sfun{fits_write_chksum}, which computes and updates the
  checksum of the HDU.  Finally the \sfun{fits_write_comment} and
  \sfun{fits_write_history} functions may be used to write comments
  and history records to the header, respectively.  
  
#%}}}

\sect{Binary Tables} 

  \sect1{Reading Binary Tables} #%{{{

 There are a several functions for reading binary tables.  The
 simplest one, \sfun{fits_read_table} reads the entire binary table
 into a structure, whose fields correspond to the names of the columns
 in the table.  (If a column has a name that contains non-alphanumeric
 characters, or does not start with an alphabetic character, then the
 structure field name for the column will be undergo the normalization
 process described for keywords.) For example, consider a file called
 \exmp{foo.fits} with a binary table whose structure is defined by the
 FITS header:
#v+
   XTENSION= 'BINTABLE'
   BITPIX  =                   8
   NAXIS   =                   2
   NAXIS1  =                  34
   NAXIS2  =                 500
   PCOUNT  =                   0
   GCOUNT  =                   1
   TFIELDS =                   4
   TTYPE1  = 'TIME    '
   TFORM1  = 'D       '
   TTYPE2  = 'X       '
   TFORM2  = 'E       '
   TTYPE3  = 'Y       '
   TFORM3  = 'E       '
   TTYPE4  = 'PHAS    '
   TFORM4  = '9I      '
   EXTNAME = 'EXAMPLE '
   TDIM4   = '(3,3)   '
#v-
 This header shows that the binary table is 500 rows in length and
 contains 4 columns with names TIME, X, Y, and PHAS.  The table may be
 read via
#v+
   tbl = fits_read_table ("foo.fits[EXAMPLE]");
#v-
 assigning a structure to the \exmp{tbl} variable.  The structure has
 fields with names \exmp{time}, \exmp{x}, \exmp{y}, and \exmp{phas},
 which be displayed via
#v+
   vmessage ("tbl.time = %S", tbl.time);
   vmessage ("tbl.x = %S", tbl.x);
   vmessage ("tbl.y = %S", tbl.y);
   vmessage ("tbl.phas = %S", tbl.y);
#v-
 producing:
#v+
   tbl.time = Double_Type[500]
   tbl.x = Float_Type[500]
   tbl.y = Float_Type[500]
   tbl.z = Short_Type[500,3,3]
#v-
 Note that the \exmp{fits_read_table} function not only read the data
 in a way that preserved the data type, but it also correctly
 identified the \exmp{phas} column as one containing a 3x3 image in
 every row!

 Often one is interested in only a few columns of a table.  Instead of
 reading the entire table, which could use a lot of memory for a large
 table, the \exmp{fits_read_table} function may also be used to read
 just the specified columns, e.g.,
#v+
    tbl = fits_read_table ("foo.fits[EXAMPLE]", "x", "y");
#v-
 will read just the \exmp{X} and \exmp{Y} columns.

 An alternative interface with much the same functionality is provided
 by the \exmp{fits_read_col} function.  Instead of returning the data
 as a structure, it returns the data as multiple return values, e.g.,
#v+
   t = fits_read_col ("foo.fits[EXAMPLE]", "time");
   (x,y) = fits_read_col ("foo.fits[EXAMPLE]", "x", "y");
#v-

 The \sfun{fits_read_cell} function may be used to read a single cell
 in the table.  For example
#v+
   phas = fits_read_cell ("foo.fits[EXAMPLE]", "phas", 4);
#v-
 will return the 3x3 array of the ``phas'' column in the fourth row.
 Finally, the \exmp{fits_read_cells} function may be used to read a
 specified range of rows in the table.  For instance,
#v+
   (x,y) = fits_read_cells ("foo.fits[EXAMPLE]", "x", "y", 1, 1000);
#v-
 will read the first 1000 rows of the \exmp{X} and \exmp{Y} columns in
 the table.

#%}}}

  \sect1{Writing Binary Tables} #%{{{

 The high-level interface has several functions that are useful for
 the creation of a binary table.  Chief among them is the
 \sfun{fits_write_binary_table} function, which supports several
 methods of calling it.  The simplest use was illustrated earlier.
 Here more complicated uses will be considered.

 As a first step, suppose that the binary table is to contain data for
 the Lissajous curve constructed as follows:
#v+
    A_x = 10.0; omega_x = 3.0;  phi_x = 0.0;
    A_y = 20.0; omega_y = 7.0;  phi_y = 1.0;
    t = [0:100:0.01];
    x = A_x * cos (omega_x*t + phi_x);
    y = A_y * cos (omega_y*t + phi_y);
#v-
 The goal is to write out arrays \exmp{t}, \exmp{x}, and \exmp{y} to a
 binary table called LISSAJOUS, and with columns of the corresponding
 names.  The easiest way is to use:
#v+
    data = struct {t, x, y}; data.t = t; data.x = x; data.y = y;
    fits_write_binary_table ("foo.fits", "LISSAJOUS", data);
#v-
 Now suppose that it is desired to write the parameters defining the
 Lissajous pattern as keywords and to write a history record to the
 file.  One way to do this is via the \exmp{fits_update_key},
 \sfun{fits_write_history}, and \sfun{fits_write_comment} functions:
#v+
    fp = fits_open_file ("foo.fits[LISSAJOUS]");
    fits_write_comment (fp, "This table contains data for a Lissajous pattern");
    fits_update_key (fp, "A_X", A_x, "x(t) Amplitude");
    fits_update_key (fp, "A_Y", A_y, "y(t) Amplitude");
    fits_update_key (fp, "OMEGA_X", omega_x, "x(t) omega");
    fits_update_key (fp, "OMEGA_Y", omega_y, "y(t) omega");
    fits_update_key (fp, "PHI_X", phi_x, "x(t) phase");
    fits_update_key (fp, "PHI_Y", phi_y, "y(t) phase");
    fits_write_history (fp, "This was written as an example for the " +
                            "documentation of the slang cfitsio module");
    fits_close_file (fp);
#v-
 The advantage of using the \exmp{fits_update_key} is that it allows
 control over the comment associated with the keyword; however,
 repeated calls to \sfun{fits_update_key} can become tedious.

 A simpler mechanism to achieve this goal is to pass the keywords and
 history information to the \sfun{fits_write_binary_table} function as
 an optional arguments,
#v+
    keys = struct {A_x, omega_x, phi_x, A_y, omega_y, phi_y};
    set_struct_fields (keys, A_x, omega_x, phi_x, A_y, omega_y, phi_y);
    hist = struct {history, comment};
    hist.comment = "This table contains data for a Lissajous pattern";
    hist.history = "This was written as an example for the " +
                    "documentation of the slang cfitsio module";
    fits_write_binary_table ("foo.fits", "LISSAJOUS", data, keys, hist);
#v-
 to produce:
#v+
    XTENSION= 'BINTABLE' / binary table extension
    BITPIX  =                   8 / 8-bit bytes
    NAXIS   =                   2 / 2-dimensional binary table
    NAXIS1  =                  24 / width of table in bytes
    NAXIS2  =               10000 / number of rows in table
    PCOUNT  =                   0 / size of special data area
    GCOUNT  =                   1 / one data group (required keyword)
    TFIELDS =                   3 / number of fields in each row
    TTYPE1  = 't       ' / label for field   1
    TFORM1  = 'D       ' / data format of field: 8-byte DOUBLE
    TTYPE2  = 'x       ' / label for field   2
    TFORM2  = 'D       ' / data format of field: 8-byte DOUBLE
    TTYPE3  = 'y       ' / label for field   3
    TFORM3  = 'D       ' / data format of field: 8-byte DOUBLE
    EXTNAME = 'LISSAJOUS' / name of this binary table extension
    A_X     =                  10
    OMEGA_X =                   3
    PHI_X   =                   0
    A_Y     =                  20
    OMEGA_Y =                   7
    PHI_Y   =                   1
    COMMENT This table contains data for a Lissajous pattern
    HISTORY This was written as an example for the documentation of the slang cfitsi
    HISTORY o module
#v-

  It is important to note that in the the above examples, the name of
  the file to contain the binary table was explicitly passed to the
  \sfun{fits_write_binary_table} function.  This causes
  \sfun{fits_write_binary_table} to create a \em{new} file containing
  the binary table, and if a file of that name exists, \em{it will be
  deleted} before the new one is created. 
  
  To append a table to an existing file, first open it using the
  \sfun{fits_open_file} function, and then use the file pointer in
  place of the name:
#v+
     fp = fits_open_file ("foo.fits", "w");   % <<-- note the "w"
        .
        .
     fits_write_binary_table (fp, ....);
     fits_close_file (fp);
#v-
  This technique must also be used to create a file containing
  multiple binary tables:
#v+
     fp = fits_open_file ("foo.fits", "c");   % <<-- note the "c"
        .
        .
     fits_write_binary_table (fp, ....);      % first table
        .
        .
     fits_write_binary_table (fp, ....);      % second table
     fits_close_file (fp);
#v-

#%}}}

\sect{Images}

 \sect1{Preliminaries} #%{{{

 Dealing with FITS images in \slang is easy as long as one understands
 that FITS stores images in FORTRAN \em{column-major} order, whereas
 \slang utilizes a C \em{row-major} order.  That is, the first
 dimension of a FITS array varies fastest whereas it is the last
 dimension of a \slang array that varies fastest.  This difference is
 automatically accounted for by the underlying \module{cfitsio}
 module.  In other words, images may be used in \slang scope as
 ordinary \slang arrays where the first dimension varies slowest, and
 the \module{cfitsio} module will make the necessary translations when
 reading or writing an image from a file.  An easy way to remember the
 \slang or C ordering is that for a 2-d array, the first index is a
 row index and the second a column--- the same as matrices are indexed
 in linear algebra.  Do not fall into the trap of indexing a \slang array
 the same as you would of indexing a point in Cartesian space (x,y),
 instead think in terms of rows and columns.


#%}}}
   
 \sect1{Reading and Writing Images} #%{{{
 
 The \sfun{fits_read_img} may be used to read the image from the
 primary FITS HDU or a FITS image extension.  It is not designed to
 read images that are stored in binary tables--- there are other
 functions for that purpose.  The \sfun{fits_read_img} function
 simply returns the data as an array of the appropriate type and
 dimensions (in row-major order) with any scaling defined via the
 BZERO and BSCALE header keywords applied.

 Writing an image HDU is somewhat more involved than reading one because
 in addition to writing the image data, the header must first be set up to
 describe the image.  Fortunately, the high-level functions make this
 easy.
 
 Suppose that one has created an image array via, e.g., the
 \module{histogram} module's \ifun{hist2d} function from the X and Y
 columns of a binary table:
#v+
    (x,y) = fits_read_col ("evt2.fits[EVENTS]", "X","Y");
    xgrid = 3840.5 + [0:1023:2];
    ygrid = 3840.5 + [0:1023:2];
    img = hist2d (y, x, ygrid, xgrid);
#v-
 and that one wants to write this out to a fits file called
 \file{img.fits}.  The simplest way to do this is using
 \sfun{fits_write_image_hdu}:
#v+
    fits_write_image_hdu ("img.fits", NULL, img);
#v-
 
 Note that the data-type of the image array controls the type of image
 written to the fits file.  If the image array is an array of
 \exmp{Double_Type}s, then the image will be written with BITPIX set
 to -64.  To have such an array out as 16 bit integers, then the array
 must first be scaled to the range of a 16 bit integer and then
 typecast to \var{Int16_Type}:
#v+
    int16_image = typecast (img, Int16_Type);
#v-

 Additional keywords may be written to the image HDU using the
 \sfun{fits_update_key} function.  And like
 \sfun{fits_write_binary_table}, the \sfun{fits_write_image_hdu}
 function takes optional parameters that specify additional keywords,
 history, and comments to be written.  The reader is referred to the
 discussion of the \sfun{fits_write_binary_table} function for more
 information.
 
#%}}}

#% \sect{Copying Headers}

\sect{WCS Routines}

 \sect1{Introduction}

  The FITS package includes a set of routines for reading and
  writing \url{http://fits.gsfc.nasa.gov/fits_wcs.html}{WCS} keywords
  in the form proposed by
  \url{http://www.atnf.csiro.au/people/mcalabre/WCS/wcs.pdf}{Greisen
  and Calabretta}.
  Although they are part of the high-level interface, the routines are
  somewhat experimental and as such must be loaded separately via:
#v+
    require ("fitswcs");
#v-

  The routines in this interface deal with a structure that describes
  the WCS via the following fields:
\begin{descrip}
   \tag{naxis}
     The number of axes to transform (Int_Type)
   \tag{ctype}
     Specifies the WCS transformation (String_Type[naxis])
   \tag{cunit}
     Units (String_Type[naxis])
   \tag{crval}
     WCS values at the reference pixel (Double_Type[naxis])
   \tag{crpix}
     The coordinates of the reference pixel (Double_Type[naxis])
   \tag{cdelt}
     Species the gradient at the reference pixel (Double_Type[naxis])
   \tag{pc}
     An array used to linearly transform the WCS.
     (Double_Type[naxis,naxis] or NULL)
   \tag{pv} 
     An array of addition parameters used to specify the WCS (NULL in
     most cases)
   \tag{ps}
     An array of additional string parameters (NULL in most cases).
   \tag{wcsname}
     A name given to this coordinate system.
\end{descrip}

  While the user is encouraged to understand the FITS WCS conventions
  and the precise meanings of these fields, the \exmp{fitswcs} interface
  provides routines to make the use of this structure as transparent
  as possible for the most common uses.  A few examples will
  illustrate this.

 \sect1{Examples}

  Consider once again the example of creating a FITS image by binning
  two columns of a FITS binary table:
#v+
    (x,y) = fits_read_col ("evt2.fits[EVENTS]", "X","Y");
    xgrid = 3840.5 + [0:1023:2];
    ygrid = 3840.5 + [0:1023:2];
    img = hist2d (y, x, ygrid, xgrid);
    fits_write_image_hdu ("img.fits", NULL, img);
#v-
  Unfortunately the resulting file will contain none of the WCS
  information that was attached to the X and Y columns from which the
  image was constructed.  One might be tempted to simply copy that
  information to the output file with the aid of the \exmp{fitswcs}
  routines via
#v+
    wcs = fitswcs_get_column_wcs ("evt2.fits[EVENTS]", ["Y", "X"]);
    fitswcs_put_img_wcs ("img.fits", wcs);
#v-
  The problem with this approach is that the WCS read from the binary
  table does not describe the image created from it because it knows
  nothing about how the image was binned nor how the image pixel
  coordinates relate back to the X and Y columns.  That information is
  contained in the definition of the grids passed to the \ifun{hist2d}
  function:
#v+
    xgrid = 3840.5 + [0:1023:2];
    ygrid = 3840.5 + [0:1023:2];
#v-
  These grids describe a simple linear transformation from image pixel
  coordinates to the (X,Y) coordinates of the binary table.  Since the
  transformation is linear, the \exmp{fitswcs_bin_wcs} function may be 
  used to transform the WCS:
#v+
    wcs = fitswcs_bin_wcs (wcs, ygrid, xgrid);
#v-
  It is the transformed WCS that is to be written out:
#v+
    fitswcs_put_img_wcs ("img.fits", wcs);
#v-

  It is important to note the order in which the X and Y arguments
  were used.  Recall that FITS stores images in a FORTRAN column-major
  order whereas \slang uses a row-major order.  For this reason,
  ``row-like'' parameters come before ``column-like'' parameters in
  statements such as
#v+
    img = hist2d (y, x, ygrid, xgrid);
    wcs = fitswcs_get_column_wcs ("evt2.fits[EVENTS]", ["Y", "X"]);
    wcs = fitswcs_bin_wcs (wcs, ygrid, xgrid);
#v-

 \sect1{Alternate WCS}

  Sometimes it is useful to attach more than one coordinate system to
  an image.  For example, it is useful to have a coordinate system
  that maps the pixel coordinates back to (X,Y) coordinates from which
  they were derived.  The \sfun{fitswcs_new_img_wcs} may be used to
  construct a linear WCS corresponding to the linear coordinate grids
  of the image:
#v+
    wcsP = fitswcs_new_img_wcs (ygrid, xgrid);
    wcsP.wcsname = "PHYSICAL";
    fitswcs_put_img_wcs ("img.fits", wcsP, 'P');
#v-
  Note that the WCS was given the name ``PHYSICAL''.  While not
  required, this enables this alternate coordinate system to be
  displayed as the physical system by the DS9 image display program.

 \sect1{Degenerate Axes}

  Consider a FITS file \exmp{hydra.fits} containing an
  image HDU with the following FITS header: 
#v+
    SIMPLE  =                   T / file does conform to FITS standard
    BITPIX  =                 -32 / number of bits per data pixel
    NAXIS   =                   4 / number of data axes
    NAXIS1  =                 657 / length of data axis
    NAXIS2  =                 657 / length of data axis
    NAXIS3  =                   1 / length of data axis
    NAXIS4  =                   1 / length of data axis
    CTYPE1  = 'RA---SIN'
    CRVAL1  =         139.5235701
    CRPIX1  =                 330
    CDELT1  =    -0.0004166666768
    CTYPE2  = 'DEC--SIN'
    CRVAL2  =      -12.0955450949
    CRPIX2  =                 328
    CDELT2  =     0.0004166666768
    PC1_1   =                   1
    PC1_2   =-1.7318465835227e-09
    PC2_1   = 1.7318465835227e-09
    PC2_2   =                   1
    CTYPE3  = 'FREQ    '
    CRVAL3  =        332902343.75
    CRPIX3  =                   1
    CDELT3  =           2490234.5
    CTYPE4  = 'STOKES  '
    CRVAL4  =                   1
    CRPIX4  =                   1
    CDELT4  =                   1
#v-
 This particular image had so-called ``degenerate axes'' added, which
 had the effect of increasing its dimensionality from 2 to 4. As such,
 this image may be rejected by some image display programs that expect
 a 2-d image.  In fact,
#v+
    img = fits_read_img ("hydra.fits");
    wcs = fitswcs_get_img_wcs ("hydra.fits");
#v-
 will read the image as a \exmp{Float_Type[1,1,657,657]} object, and
 the WCS as a 4-d with wcs.ctype equal to
#v+
    ["STOKES", "FREQ", "DEC--SIN", "RA---SIN"]
#v-
 
 The degenerate dimensions may be removed from the image via
#v+
    img = img[0,0,*,*];
#v-
 producing a 2d image of type \exmp{Float_Type[657,657]}.  The
 corresponding wcs may be obtained using the \sfun{fitswcs_slice}
 function to extract the last two dimensions of the WCS:
#v+
   wcs = fitswcs_slice (wcs, [2,3]);
#v-

 Another use of the \sfun{fitswcs_slice} is to reorder the dimensions
 of the WCS.  For example, earlier it was pointed out that when
 constructing an image from columns in a table, that one read the WCS
 in a row-major order.  If the reverse order was used when obtaining
 the WCS from the columns of a binary table, e.g.,
#v+
    wcs = fitswcs_get_column_wcs ("evt2.fits[EVENTS]", ["X", "Y"]);
#v-  
 then it would have been necessary to reverse the order of the
 dimensions of the WCS structure.  The \sfun{fitswcs_slice} may be
 used to swap the dimensions of the WCS, e.g.,
#v+
    wcs = fitswcs_slice (wcs, [1,0]);
#v-

#%+
(Recall that the FITS convention is to assign the pixel coordinate
  (0.5,0.5) to the extreme corner of the first pixel, whereas many
  systems, including \slang, use (0,0) for that position and regard
  the center of the pixel as being at (0.5,0.5). 
#%-  

#% \sect{Examples}

\sect{High-level Function Reference}
#i fitsfuns.tm

\sect{WCS Function Reference}
#i fitswcsfuns.tm

\chapter{The low-level interface}

\sect{Overview} #%{{{

 Functions in the low-level module are usually needed when it is
 necessary to perform some task that is not readily achievable using
 the high-level interface.  This module may be loaded using
#v+
    require ("cfitsio");
#v-
 When mixing functions from both interfaces, it is not necessary to
 explicitly load the \module{cfitsio} module in this manner since it
 is loaded automatically by the high-level interface.

 For the most part, for those functions that have been wrapped, the
 \module{cfitsio} module represents a 1-1 mapping between the
 functions of the cfitsio library and those of the module.  For this
 reason a detailed description of the functions in the
 \module{cfitsio} module will not be given here; the reader is
 referred to the documentation for the \CFITSIO library itself for the
 details.  Here only the semantic differences between the functions in
 the module and those of the library, and how the functions are
 documented.

 Most \CFITSIO functions adhere to a so-called ``inherited status''
 convention via a ``status'' argument as the last parameter.  In
 addition functions also return the error status as a return value.
 For simplicity the wrapping by the module does not respect this
 convention.  That is, none of the module's functions take a status
 argument.  For example, the \CFITSIO documentation for the
 \exmp{fits_get_num_hdus} specifies that it is to be called from C
 via:
#v+
    status = fits_get_num_hdus (fptr, &hdunum, &status);
#v-
 This function has been wrapped such that it is to be called from
 \slang via
#v+
   status = _fits_get_num_hdus (fptr, &hdunum);
#v-

#%}}}

\sect{Low-level Function Reference}


#s+
#i mkindex.sl
#s-
#d iflatex#2 <#if output=latex2e>$1</#if><#unless output=latex2e>$2</#unless>
#d ifhtml#2 <#if output=html>$1</#if><#unless output=html>$2</#unless>

#d xreferences#1 This function is a wrapper around the cfitsio library \
 function \exmp{$1}. \
 See \ifhtml{\url{\cfitsio_fun_url{$1}}{its documentation}}{its documentation}\
 for additional information.

#d cfitsioxref#1 \ifhtml{\url{\cfitsio_fun_url{$1}}{$1}}{\exmp{$1}}

#i rtl/cfitsiofuns.tm

\end{\documentstyle}
