
\function{maplib_new}
\synopsis{Create a new mapping instance}
\usage{m = maplib_new (String_Type mapname)}
\description
 This function creates a new instance of a mapping of the specified
 name.  If the mapping exists, the function returns a structure whose
 fields depend upon the specified mapping.  Otherwise the function
 will generate an exception.
 
 The following mappings are supported:
#v+
    mapname          description
      gnomic           Gnomic
      ortho            Orthographic
      stereo           Stereographic
      lambert          Lambert Azimuthal Equal-Area
      hammer           Hammer-Aitoff
      sphere           Rotation of a sphere
      linear           Linear Transformation
#v-
 The fields of the structure returned by the \ifun{maplib_new}
 function vary with the name of the mapping.  Unless otherwise
 specified, most of the maps contain the following fields:
#v+
   field-name   default-value             description
     name         <varies>         The name of the transformation
     lon0            0             The longitude of the map center
     lat0            90            The latitude of the map center
     x0              0             (x0,y0) are the coordinates of the
     y0              0               plane at the map center
     xscale          1             Scale factor in the x direction
     yscale          1             Scale factor in the y direction
#v-
 For tangent plane projections, (lon0,lat0) specify the position of
 the tangent plane on the sphere.
 
 The following transformations return a structure whose fields differ
 from the above:
\begin{descrip}
  \tag{sphere}
#v+
   field-name   default-value             description
     name         sphere         Rotation of the unit sphere
     lon0,lat0    (0,90)         The rotation takes the point (lon0,lat0)
     lon1,lat1    (0,90)           into the point (lon1,lat1)
     beta           0            Additional rotation about (lon1,lat1)
#v-
  \tag{linear}
#v+
   field-name   default-value             description
     name         linear         A linear transformation
     x0,y0        (0,0)            (x',y') = (x1,y1) + A#(x-x0,y-y0)
     x1,y1        (0,0)
     A           [1,0,0,1]
#v-
\end{descrip}
\example
#v+
    m = maplib_new ("gnomic");
    m.lon0 = 20;
    m.lat0 = 35;
    (x,y) = maplib_project (m, 45, 22);
#v-
\seealso{maplib_project, maplib_deproject, maplib_reproject}
\done

\function{maplib_project}
\synopsis{Project}
\usage{(x,y) = maplib_project (Struct_Type m, Double_Type lon, Double_Type lat)}
\description
 This function applies the mapping specified by the first parameter
 \exmp{m} to the longitude and latitude pairs \exmp{(lon,lat)}.  The
 coordinates of the projected points will be returned.
\notes
  For some transformations, the parameters \exmp{lon} and \exmp{lat}
  need not represent longitude and latitude coordinates.  The
  \exmp{"linear"} transformtion is one example.

 This function is the inverse of \ifun{maplib_deproject}.
\seealso{maplib_new, maplib_deproject, maplib_reproject}
\done

\function{maplib_deproject}
\synopsis{Project}
\usage{(lon,lat) = maplib_project (Struct_Type m, Double_Type x, Double_Type y)}
\description
 This function applies the inverse of the mapping specified by the
 first parameter \exmp{m} to the projected coordinates \exmp{(x,y)}
 and returns the result.
\notes
 This function is the inverse of \ifun{maplib_project}.
\seealso{maplib_new, maplib_project, maplib_reproject}
\done


\function{maplib_reproject}
\synopsis{Reproject the coordinates to another projection}
\usage{(x1,y1)=maplib_reproject (m_to, m_from, x, y)}
#v+
   Struct_Type m_to, m_from;
   Double_Type x, Double_Type y;
#v-
\description
 The \ifun{maplib_reproject} function reprojects the coordinates
 \exmp{(x,y)} associated with the mapping \exmp{m_from} to coordinates
 described by the mapping \exmp{m_to}.  This function is equivalent to
#v+
     (lon, lat) = maplib_deproject (m_from, x, y);
     (x1, y1) = maplib_reproject (m_to, lon, lat);
#v-
\notes
 It is strongly recommended that this function be used for
 reprojections rather than the two step process outlined above.  The
 main reason is that if \exmp{(x,y)} are expressed in single
 precision, then the intermediate coordinates \exmp{(lon, lat)} will
 be returned as single precision values.  Hence, there could be some
 loss in precision when using the two step process.  Using
 \ifun{maplib_reproject} does not have this problem because double
 precision is used throughout.
\seealso{maplib_new, maplib_project, maplib_deproject}
\done
