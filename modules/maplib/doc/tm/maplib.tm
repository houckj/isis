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

#d function#1 \sect{<bf>$1</bf>\label{$1}}<descrip>
#d variable#1 \sect{<bf>$1</bf>\label{$1}}<descrip>
#d function_sect#1 \sect{$1}
#d begin_constant_sect#1 \sect{$1}<itemize>
#d constant#1 <item><tt>$1</tt>
#d end_constant_sect </itemize>

#d synopsis#1 <tag> Synopsis </tag> $1
#d keywords#1 <tag> Keywords </tag> $1
#d usage#1 <tag> Usage </tag> <tt>$1</tt>
#d description <tag> Description </tag>
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
#d documentstyle book

#%}}}

#d module#1 \tt{$1}
#d file#1 \tt{$1}
#d slang-documentation \
 \url{http://www.jedsoft.org/slang/doc/html/slang.html}{S-Lang documentation}

\linuxdoc
\begin{\documentstyle}

\title S-Lang Maplib Module Reference
\author John E. Davis, \tt{<davis@space.mit.edu>}
\date \__today__

\toc

\chapter{Introduction to the Maplib Module}
The routines in this module implement various mapping projections from
a unit sphere to a surface such as a plane or cylinder.

In all the projections, points on the sphere are specified using
latitude and longitude coordinates expressed in degrees. Values of
latitude run from -90 degrees at the south pole to +90 degrees at the
north pole.  Longitude values run from -180 to 180 in an "easterly"
direction.

At the moment, the following map mappings are supported by the module:
\begin{descrip}
   \tag{Gnomic}
     A perspective projection from the center of a sphere to a plane
     tangent to some point on the surface of the sphere.  Great
     circles on the sphere are mapped to straight lines on the plane.
   \tag{Orthographic}
     A perspective projection of the unit sphere to a tangent plane
     on the unit sphere from an infinite distance.
   \tag{Stereographic}
     A perspective projection to a tangent plane from the point on the
     unit sphere that is antipodal to the point of tangency.
   \tag{Lambert}
     An azimuthal equal area projection.
   \tag{Hammer-Aitoff}
     An equal area projection of the sphere onto an ellipse.
\end{descrip}
 In addition, the module supports the additional transformations:
\begin{descrip}
   \tag{Linear}
     A simple 1-1 mapping from a point on one plane to a point on
     another.  The mapping may involve translation, scaling, or
     rotation.
   \tag{Sphere}
     A rotation of a sphere onto itself.
\end{descrip}

Addition map projections will be added in future versions of the module.

\chapter{Usage}

 In order to a particular map projection, it must first be created via
 the \ifun{maplib_new} function.  It takes a single string parameter
 that specifies the particular mapping and returns a structure that
 defines the map.  For example, a Gnomic transformation may be
 allocated using
#v+
     p = maplib_new ("gnomic");
#v-
  The structure returned by the \ifun{maplib_new} function contains
  fields that are specified to the mapping.  In the case of the Gnomic
  transformation, the structure has the fields:
#v+
   field-name   default-value             description
     name          gnomic          The name of the transformation
     lon0            0             The longitude of the tangent point
     lat0            90            The latitude of the tangent point
     x0              0             (x0,y0) are the coordinates of the
     y0              0               tangent point in the tangent plane
     xscale          1             Scale factor in the x direction
     yscale          1             Scale factor in the y direction
#v-
  As the table shows, the default location of the tangent plane is at
  the north pole of the unit sphere.  The location of the tangent
  plane may be changed to an arbitrary point on the sphere by setting
  the \exmp{lon0} and \exmp{lat0} parameters to the desired values.

  The \ifun{maplib_project} function is used to carry out the mapping
  of a point on the unit sphere.  It takes three arguments: the
  structure that defines the transformation, and the longitude and
  latitude of the point to be mapped.  For example, 
#v+
     (x,y) = maplib_project (p, 10, 70);
#v-
  computes the coordinates of point with a longitude of 10 degrees and
  a latitude of 70 degrees under the mapping specified by \exmp{p}.

  Similarly, the reverse transformation back to the sphere may be
  carried out using the \ifun{maplib_deproject} function:
#v+
    (lon,lat) = maplib_deproject (p, x, y);
#v-
  
  For some purposes one is interested in transformation from one image
  surface to another.  For example, consider two tangent planes: one at
  the north pole and one at a longitude of 20 and a latitude of 80
  degrees.  To transform points from the latter plane to the one at the
  north pole under the gnomic mapping, it is necessary to first
  allocate the appropriate transformations:
#v+
    A = maplib_new ("gnomic"); A.lon0 = 0; A.lat0 = 90;
    B = maplib_new ("gnomic"); B.lon0 = 20; B.lat0 = 80;
#v-
  One way to transform a point \exmp{(x_B,y_B)} on the plane defined
  by the transformation \exmp{B} to the plane at the north pole
  defined by \exmp{A} is via:
#v+
    (lon, lat) = maplib_deproject (B, x_B, y_B);
    (x_A, y_A) = maplib_project (A, lon, lat);
#v-
  The problem with this approach is that there may be some loss of
  accuracy by computing the intermediate coordinates.  This is because
  if the coordinates given to either \exmp{maplib_project} or
  \exmp{maplib_deproject} are single precision values, the functions
  will return the results as single precision values.  For this
  reason, the \module{maplib} module includes a function called
  \ifun{maplib_reproject} that avoids this problem:
#v+
    (x_A,y_A) = maplib_reproject (A, B, x_A, x_B);
#v-

\chapter{Maplib Module Function Reference}
#i maplibfuns.tm

\end{\documentstyle}
