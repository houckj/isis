#c __FILE__: ../../src/fitswcs.sl
#c __LINE__: 96
\function{fitswcs_new}
\synopsis{Create a new-ndimensional linear WCS}
\usage{wcs = fitswcs_new (Int_Type naxis)}
\description
  This function returns a new WCS structure of the specified dimensionality
  that represents an identity (linear) transformation.
\seealso{fitswcs_get_img_wcs, fitswcs_get_column_wcs, fitswcs_get_vector_wcs}
\done
#c __LINE__: 121
\function{fitswcs_slice}
\synopsis{Form a new wcs from one or more axes of another}
\usage{new_wcs = fitswcs_slice (wcs, dims)}
\description
  This function may be used to construct a new wcs from another by rearranging 
  its axes or by using a subset of them.  The \exmp{dims} argument specifies 
  the dimensions to use.
\example
  Suppose that \exmp{wcs} represents a 4 dimensional WCS. Then
#v+
    wcs2 = fitswcs_slice (wcs, [0,1]);
#v-
  will result in a 2 dimensional WCS from the first 2 axis of the input WCS.
  Similarly,
#v+
    wcs2 = fitswcs_slice (wcs, [1,0]);
#v-
  will produce a 2d WCS with the first two axes swapped.
\seealso{fitswcs_get_img_wcs, fitswcs_get_column_wcs, fitswcs_get_vector_wcs}
\done
#c __LINE__: 312
\function{fitswcs_get_img_wcs}
\synopsis{Read a WCS for a FITS image}
\usage{wcs = fitswcs_get_img_wcs (fp [,alt])}
\description
  The \sfun{fitswcs_get_img_wcs} returns a structure representing a WCS from
  the specified file descriptor \exmp{fp} corresponding to an image HDU.
  An optional parameter may be used to specified an alternate WCS.
\example
#v+
   wcs = fitswcs_get_img_wcs ("img.fits[IMAGE]", 'P');
#v-
\seealso{fitswcs_put_img_wcs, fitswcs_get_column_wcs, fitswcs_get_vector_wcs}
\done
#c __LINE__: 371
\function{fitswcs_get_column_wcs}
\synopsis{Get the WCS attached to one or more columns of a binary table}
\usage{fitswcs_get_column_wcs (fp, columns-array [,alt]}
\description
  This function may be used to obtain the WCS associated with one or more
  columns of a binary table.  The file descriptor \exmp{fp} must specify 
  a binary table.  The \exmp{columns-array} parameter should be an array
  of columns names.  The third parameter is optional and is used to specify
  an alternate WCS.
\example
#v+
   wcs = fitswcs_get_column_wcs ("evt1.fits[EVENTS]", ["X","Y"]);
#v-
\seealso{fitswcs_put_column_wcs, fitswcs_get_img_wcs, fitswcs_get_vector_wcs}
\done
#c __LINE__: 480
\function{fitswcs_get_vector_wcs}
\synopsis{Get the WCS of an image in a specified table cell}
\usage{wcs = fitswcs_get_vector_wcs (fp, column_name, row [,alt])}
\description
  This function reads the WCS of an image in a specified cell of a binary
  table given by \exmp{fp} parameter.  The second and third parameters specify
  the column name and row number of the cell.  An optional fourth parameter
  may be used to obtain the corresponding alternate WCS.
\example
  This example reads the WCS associated with the image in the second row
  of the QEU column of the binary table with HDUNAME equal to AXAF_QEU1
  in the file "HRCQEU.fits":
#v+
    wcs = fitswcs_get_vector_wcs ("HRCQEU.fits[AXAF_QEU1], "QEU", 2);
#v-
\notes
  The current implementation does not yet support references to the WCS
  of other cells.
\seealso{fitswcs_get_column_wcs, fitswcs_get_img_wcs}
\done
#c __LINE__: 577
\function{fitswcs_new_img_wcs}
\synopsis{Create a linear WCS for an image}
\usage{wcs = fitswcs_new_img_wcs (grid0,grid1,...)}
\description
  This function may be used to construct a linear WCS for an image with the
  specified grids.  The grids are assumed to be linear.
\example
  Use the histogram module's hist2d function to create an image from the X
  and Y columns in a file, and the construct a corresponding WCS:
#v+
    (x,y) = fits_read_col ("table.fits", "X", "Y");
    gridx = [min(x):max(x):0.5];
    gridy = [min(y):max(y):0.5];
    img = hist2d (y,x,gridy,gridx);
    wcs = fitswcs_new_img_wcs (gridy, gridx);
#v-
\seealso{fitswcs_new, fitswcs_get_img_wcs}
\done
#c __LINE__: 628
\function{fitswcs_put_img_wcs}
\synopsis{Write a WCS out to an image header}
\usage{fitswcs_put_img_wcs (fp, wcs [,alt])}
\description
  The \sfun{fitswcs_put_img_wcs} may be used to write the specified wcs
  out to the image HDU specified by the \exmp{fp} parameter.  An optional 
  third parameter may be used to specify an alternate WCS.
\example
#v+
    fp = fits_open_file ("img.fits", "w");
      .
      .
      .
    fits_put_img_wcs (fp, wcs, 'P');
    fits_close_file (fp);
#v-
\seealso{fitswcs_put_column_wcs}
\done
#c __LINE__: 705
\function{fitswcs_put_column_wcs}
\synopsis{Write the WCS attached to one or more table columns}
\usage{fitswcs_put_column_wcs (fp, wcs, columns-array [,alt])}
\description
  This function may be used to attach a WCS to one or more columns of a binary
  table.  The dimensionality of the specified WCS must match the length of the
  array specifying the column names.  The first parameter, \exmp{fp} must specify
  a binary table extension.  The fourth parameter is optional and may be used
  to specify an alternate WCS.
\example
#v+
   fitswcs_put_column_wcs ("evt2.fits[EVENTS], wcs, ["X","Y"]);
#v-
\seealso{fitswcs_get_column_wcs, fitswcs_put_img_wcs, fitswcs_get_img_wcs}
\done
#c __LINE__: 818
\function{fitswcs_linear_transform_wcs}
\synopsis{Apply a linear transformation to a WCS}
\usage{wcs1 = fitswcs_linear_transform_wcs (wcs, X0, CD, I0)}
#v+
     wcs: The specified WCS to transform
   X0,I0: 1-d arrays
      CD: 2-d array
#v-
\description
  This function may be used to create a new WCS by applying a linear 
  transformation to an existing one.
\notes
  The dimensionality of the WCS is limited to 2 in the 
  current implementation.
\seealso{fitswcs_rebin_wcs}
\done
#c __LINE__: 907
\function{fitswcs_rebin_wcs}
\synopsis{This function may be used to obtain the wcs for a rebinned image}
\usage{wcs1 = fitswcs_rebin_wcs (wcs, grid0, grid1, ...)}
\description
  This function may be used to construct the WCS for a rebinned image from
  the WCS of of the unbinned image.  The grid parameters specify the linear
  grids the new image.
\seealso{fitswcs_linear_transform_wcs, fitswcs_slice}
\done
