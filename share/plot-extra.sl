% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2012 Massachusetts Institute of Technology
%
%    This software was developed by the MIT Center for Space Research under
%    contract SV1-61010 from the Smithsonian Institution.
%
%    Author:  John C. Houck  <houck@space.mit.edu>
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
%

% FIXME - this is a hack.
private variable Saved_Transform;

define isis_plot_library_interface ();

define cursor() %{{{
{
   variable x_ref, y_ref, ch_ref;
   variable x,y,ch;

   variable pli = isis_plot_library_interface ();

   variable ci = (@pli.get_color)();
   () = (@pli.set_color)(1);
   EXIT_BLOCK
     {
	() = (@pli.set_color)(ci);
     }

   switch (_NARGS)
     {
      case 3:
	message("Click for cursor coordinates\n");
	(x, y, ch) = _isis->_cursor(0);
	(x_ref, y_ref, ch_ref) = ();
	@x_ref = x;
	@y_ref = y;
	@ch_ref = char(ch);
     }
     {
      case 2:
	message("Click for cursor coordinates\n");
	(x, y, ch) = _isis->_cursor(0);
	(x_ref, y_ref) = ();
	@x_ref = x;
	@y_ref = y;
     }
     {
      case 0:
	message("Click for cursor coordinates, type 'q' in plot window to quit\n");
	(x, y, ch) = _isis->_cursor(1);
     }
     {
	% default:
	_pop_n (_NARGS);
	message("cursor( [ &x, &y [, &ch]] );");
     }
}

%}}}

define plot_pause () %{{{
{
   variable p, t;

   t = -1;
   p = NULL;

   switch (_NARGS)
     {
      case 2:
        (t,p) = ();
     }
     {
      case 1:
        t = ();
        if (typeof (t) == String_Type)
          {
             p = t;
             t = -1;
          }
     }

   variable pli = isis_plot_library_interface ();

   if ((@pli.is_hardcopy)())
     return;

   if (p != NULL)
     {
        () = fputs (p, stdout);
        () = fflush (stdout);
     }

   if (t >= 0)
     {
        sleep (t);
        return;
     }

   () = fputs ("Press RETURN to continue\n", stdout);
   () = fflush (stdout);
   () = fgets (&p, stdin);
}

%}}}

% simple image plotting

private define img_transform (gx, gy) %{{{
{
   variable tr = Float_Type[6];
   variable dx, dy;

   variable nx, ny;
   nx = length(gx);
   ny = length(gy);

   dx = (gx[nx-1] - gx[0]) / (nx - 1.0);
   dy = (gy[ny-1] - gy[0]) / (ny - 1.0);

   tr[0] = gx[0] - dx;
   tr[1] = dx;
   tr[2] = 0.0;
   tr[3] = gy[0] - dy;
   tr[4] = 0.0;
   tr[5] = dy;

   return tr;
}
%}}}

private define restore_axes (r) %{{{
{
   xrange(r.xmin,r.xmax);
   yrange(r.ymin,r.ymax);
}

%}}}

private define set_axes (a) %{{{
{
   variable dims;
   (dims,,) = array_info (a);

   variable gx, gy;

   variable xr = struct {min, max};
   (xr.min, xr.max) = _isis->_get_xrange();
   if (0 <= xr.min and xr.max < dims[1]-1)
     gx = [int(xr.min) : int(xr.max)];
   else gx = [0:dims[1]-1];

   variable yr = struct {min, max};
   (yr.min, yr.max) = _isis->_get_yrange();
   if (0 <= yr.min and yr.max < dims[0]-1)
     gy = [int(yr.min) : int(yr.max)];
   else gy = [0:dims[0]-1];

   return (gx, gy);
}

%}}}

% FIXME -- this is now a no-op.  if nothing broke, get rid of it!
private define apply_axis_ranges (a, gx, gy) %{{{
{
#iffalse
   variable xr = struct {min, max};
   (xr.min, xr.max) = _isis->_get_xrange();
   variable ix = where (xr.min <= gx and gx <= xr.max);

   variable yr = struct {min, max};
   (yr.min, yr.max) = _isis->_get_yrange();
   variable iy = where (yr.min <= gy and gy <= yr.max);
#endif

   return (a, gx, gy); %(a[iy,ix], gx[ix], gy[iy]);
}

%}}}

private define validate_array_sizes (a, gx, gy) %{{{
{
   variable
     a_dims = array_shape (a),
     gx_dims = array_shape (gx),
     gy_dims = array_shape (gy);

   % a[row,col]  dim[0] = row = "y-axis"
   %             dim[1] = col = "x-axis"

   variable ok = ((length(a_dims) == 2 && length(gx_dims) == 1 && length(gy_dims) == 1)
                 && (gy_dims[0] == a_dims[0])
                 && (gx_dims[0] == a_dims[1]));

   ifnot (ok)
     {
        vmessage ("*** Error:  Inconsistent array shapes:");
        vmessage ("   2D array shape is [nrow, ncol] = %S", a);
        vmessage ("    axis array shapes are X[ncol] = %S", gx);
        vmessage ("                          Y[nrow] = %S", gy);
        throw UsageError;
     }
}

%}}}

private define plot_image_args (nargs, msg) %{{{
{
   variable a, aspect, gx, gy, amin, amax;

   switch (nargs)
     {
      case 1:
	a = ();
	(gx, gy) = set_axes (a);
	a = a[[gy[0]:gy[-1]], [gx[0]:gx[-1]]];
	amin = min(a);
	amax = max(a);
	aspect = 0;
     }
     {
      case 2:
	(a, aspect) = ();
	(gx, gy) = set_axes (a);
	a = a[[gy[0]:gy[-1]],[gx[0]:gx[-1]]];
	amin = min(a);
	amax = max(a);
     }
     {
      case 4:
	(a, aspect, gx, gy) = ();
	(a, gx, gy) = apply_axis_ranges (a, gx, gy);
	amin = min(a);
	amax = max(a);
     }
     {
      case 6:
	(a, aspect, gx, gy, amin, amax) = ();
	(a, gx, gy) = apply_axis_ranges (a, gx, gy);
     }
     {
	% default:
	usage (msg);
	return;
     }

   validate_array_sizes (a, gx, gy);

   return (a, aspect, gx, gy, amin, amax);
}

%}}}

private define plot_contour_args (nargs, msg) %{{{
{
   variable a, aspect, gx, gy, c, amin, amax;

   c = [1:10]/10.0;

   switch (nargs)
     {
      case 1:
	a = ();
	(gx, gy) = set_axes (a);
	a = a[[gy[0]:gy[-1]],[gx[0]:gx[-1]]];
	amin = min(a);
	amax = max(a);
	aspect = 0;

	c = amin + (amax-amin) * c;
     }
     {
      case 2:
	(a, aspect) = ();
	(gx, gy) = set_axes (a);
	a = a[[gy[0]:gy[-1]],[gx[0]:gx[-1]]];
	amin = min(a);
	amax = max(a);

	c = amin + (amax-amin) * c;
     }
     {
      case 4:
	(a, aspect, gx, gy) = ();
	if (gx == NULL or gy == NULL) (gx, gy) = set_axes (a);
	(a, gx, gy) = apply_axis_ranges (a, gx, gy);
	amin = min(a);
	amax = max(a);
	c = amin + (amax-amin) * c;
     }
     {
      case 5:
	(a, aspect, gx, gy, c) = ();
	if (gx == NULL or gy == NULL) (gx, gy) = set_axes (a);
	(a, gx, gy) = apply_axis_ranges (a, gx, gy);
     }
     {
	% default:
	usage (msg);
	return;
     }

   validate_array_sizes (a, gx, gy);

   return (a, aspect, gx, gy, c);
}

%}}}

% plot confidence contours

% The world coordinates of the array point A(I,J)
% are given by:
%     X = TR(0) + TR(1)*I + TR(2)*J
%     Y = TR(3) + TR(4)*I + TR(5)*J

private define get_transform (s) %{{{
{
   variable tr = Float_Type[6];
   variable dx, dy;

   dx = (s.px.max - s.px.min) / (s.px.num - 1.0);
   dy = (s.py.max - s.py.min) / (s.py.num - 1.0);

   tr[0] = s.px.min - dx;
   tr[1] = dx;
   tr[2] = 0.0;
   tr[3] = s.py.min - dy;
   tr[4] = 0.0;
   tr[5] = dy;

   return tr;
}

%}}}

private define make_axes (gx, gy, aspect) %{{{
{
   variable pli = isis_plot_library_interface ();

   variable r = struct
     {
	xmin,xmax,ymin,ymax
     };
   (r.xmin,r.xmax) = _isis->_get_xrange();
   (r.ymin,r.ymax) = _isis->_get_yrange();

   variable
     _xmin = (r.xmin == -DOUBLE_MAX) ? gx[ 0] : r.xmin,
     _xmax = (r.xmax == +DOUBLE_MAX) ? gx[-1] : r.xmax,
     _ymin = (r.ymin == -DOUBLE_MAX) ? gy[ 0] : r.ymin,
     _ymax = (r.ymax == +DOUBLE_MAX) ? gy[-1] : r.ymax;

   xrange(_xmin, _xmax);
   yrange(_ymin, _ymax);

   if (aspect != 0)
     () = (@pli.match_viewport_to_window)(gx[0], gx[-1], gy[0], gy[-1]);
   else if ((@pli.device_id)() != 0)
     () = (@pli.set_standard_viewport)();

   plot_box();

   return r;
}

%}}}

define plot_image () %{{{
{
   variable msg = "plot_image (img [, aspect [, x, y [, img_min, img_max]]])";
   variable a, aspect, gx, gy, amin, amax;
   variable tr;

   if (_NARGS == 0)
     usage (msg);

   (a, aspect, gx, gy, amin, amax) = plot_image_args (_NARGS, msg);

   variable r = make_axes (gx, gy, aspect);
   tr = img_transform (gx, gy);

   Saved_Transform = @tr;

   variable pli = isis_plot_library_interface ();
   () = (@pli.plot_image)(a, amin, amax, tr);

   restore_axes (r);
}

%}}}

define plot_contour () %{{{
{
   variable msg = "plot_contour (img[, aspect [, x, y [, c]]])";
   variable a, aspect, gx, gy, c;
   variable tr;

   if (_NARGS == 0)
     usage (msg);

   (a, aspect, gx, gy, c) = plot_contour_args (_NARGS, msg);

   variable r = make_axes (gx, gy, aspect);
   tr = img_transform (gx, gy);
#iffalse
   % 10/2013:  Apparently this bug got fixed in the
   %           HEASARC pgplot, so this work-around is
   %           no longer needed.
   % This seems to be a PGPLOT bug:
   tr[0] -= 0.5*tr[1];
   tr[3] -= 0.5*tr[5];
#endif
   Saved_Transform = @tr;

   % If attribute arg is
   % negative => use current line attributes
   % positive => solid lines for (+) contours, dashed for (-)
   variable pli = isis_plot_library_interface ();
   () = (@pli.plot_contour)(a, tr, c, 1);

   restore_axes (r);
}

%}}}

define oplot_contour () %{{{
{
   variable msg = "oplot_contour (img [, c])";
   variable a, c, tr;

   c = NULL;

   if (_isis->get_varargs (&a, &c, _NARGS, 1, msg))
     return;

   if (c == NULL)
     {
	variable amin, amax;
	amin = min(a);
	amax = max(a);
	c = amin + (amax - amin) * [1:10]/10.0;
     }

   if (Saved_Transform != NULL)
     tr = @Saved_Transform;
   else
     {
	tr = [0.0, 1.0, 0.0,
	      0.0, 0.0, 1.0];
     }

   % If number of contours is
   % negative => use current line attributes
   % positive => solid lines for (+) contours, dashed for (-)
   variable pli = isis_plot_library_interface ();
   () = (@pli.plot_contour)(a, tr, c, 1);
}

%}}}

define set_palette () %{{{
{
   variable msg = "set_palette (id [, contrast [,brightness]])";
   variable id, contrast, brightness;
   brightness = 0.5;
   contrast = 1.0;

   if (_isis->get_varargs (&id, &contrast, &brightness, _NARGS, 1, msg))
     return;

   variable pli = isis_plot_library_interface ();
   () = (@pli.select_color_table)(id, contrast, brightness);
}

%}}}

define plot_image_ctrl () %{{{
{
   variable pli = isis_plot_library_interface ();

   variable msg = "Use cursor to adjust color table:\n"
     + " Keys 1,2,3,4,5 select different palettes\n"
     + " Key P cycles through available palettes\n"
     + " Key F adjusts contrast and brightness, with\n"
     + "  cursor X position setting brightness [0.0 - 1.0]\n"
     + "   and Y position setting contrast [0.0 - 10.0]\n"
     + "  (Hold down F key while moving cursor to change\n"
     + "  contrast and brightness continuously)\n"
     + " Key C resets contrast=1.0, brightness=0.5\n"
     + " Key - reverses color palette\n"
     + " Key X or right mouse button exits this mode";

   message (msg);

   variable p, contra, bright;
   variable x1, x2, y1, y2;
   variable x, y, ch, sgn;

   p = 2;
   contra = 1.0;
   bright = 0.5;
   x = 0.5;
   y = 1.0;
   sgn = 1.0;

   (x1, x2, y1, y2) = (@pli.query_plot_limits)();

   variable b1, b2, c1, c2;
   b1 = 0.0;   b2 = 1.0;
   c1 = 0.0;   c2 = 10.0;

   () = (@pli.set_plot_limits)((b1, b2, c1, c2));

   forever
     {
        variable cc = struct {xanchor, yanchor, position, mode};
        cc.xanchor = 0.0;
        cc.yanchor = 0.0;
        cc.mode = 0;
        cc.position = 0;
        (x, y, ch) = (@pli.read_cursor)(cc);

	switch (ch)
	  {
	   case 'x' or case 'X':
	     %() = (@pli.set_plot_limits)(x1, x2, y1, y2);
	     xrange(x1,x2);
	     yrange(y1,y2);
	     return;
	  }
	  {
	   case 'f' or case 'F':
	     bright = max([b1, min([b2,x])]);
	     contra = max([c1, min([c2,y])]);
	  }
	  {
	   case 'c' or case 'C':
	     contra = 1.0;
	     y = 1.0;
	     bright = 0.5;
	     x = 0.5;
	  }
	  {
	   case '-':
	     sgn = -sgn;
	  }
	  {
	     case '1':
	     p = 1;
	  }
	  {
	   case '2':
	     p = 2;
	  }
	  {
	   case '3':
	     p = 3;
	  }
	  {
	   case '4':
	     p = 4;
	  }
	  {
	   case '5':
	     p = 5;
	  }
	  {
	   case 'p' or case 'P':
	     p = 1 + p mod 5;
	  }

	set_palette (p, sgn*contra, bright);
     }
}

%}}}

private define plot_best_point (s) %{{{
{
   variable pli = isis_plot_library_interface ();

   % put a crosshair at the best-fit point
   () = (@pli.set_char_size)(2);
   () = (@pli.set_color)(1);
   () = (@pli.plot_points)(s.px_best, s.py_best, 2);
   () = (@pli.set_char_size)(1);
}

%}}}

private define plot_conf_contours (s, line, c) %{{{
{
   variable pli = isis_plot_library_interface ();

   variable i, nc, tr;
   tr = get_transform (s);
   nc = length(c);
#iffalse
   % 10/2013:  Apparently this bug got fixed in the
   %           HEASARC pgplot, so this work-around is
   %           no longer needed.
   % This seems to be a PGPLOT bug:
   tr[0] -= 0.5*tr[1];
   tr[3] -= 0.5*tr[5];
#endif
   % default color sequence:
   % Red=2, Green=3, Blue=4 (PGPLOT specific!)
   variable color = [2, 3, 4];
   if (line != NULL && struct_field_exists (line, "color"))
     color = line.color;

   _for i (0, nc-1, 1)
     {
	if (line != NULL)
	  {
             if (typeof(line.width) == Array_Type)
               () = (@pli.set_line_width)(line.width[i]);
             else
               () = (@pli.set_line_width)(line.width);

             if (typeof(line.type) == Array_Type)
               () = (@pli.set_line_style)(line.type[i]);
             else
               () = (@pli.set_line_style)(line.type);
	  }

        if (typeof(color) == Array_Type)
          () = (@pli.set_color)(color[i]);
        else
          () = (@pli.set_color)(color);

        variable cpy = @s.chisqr;
        variable nan_indices = where(isnan(cpy));
        if (length(nan_indices) > 0)
          {
             cpy[nan_indices] = FLOAT_MAX;
          }

        () = (@pli.plot_contour)(cpy, tr, c[i], -1);
     }

   plot_best_point (s);
}

%}}}

private define set_axis_ranges (s) %{{{
{
   variable r = struct
     {
	xmin, xmax, ymin, ymax
     };
   (r.xmin, r.xmax) = _isis->_get_xrange();
   (r.ymin, r.ymax) = _isis->_get_yrange();

   variable xmin, xmax, ymin, ymax;
   xmin = r.xmin;
   xmax = r.xmax;
   ymin = r.ymin;
   ymax = r.ymax;

   if (xmin == -_isis->DBL_MAX) xmin = s.px.min;
   if (xmax ==  _isis->DBL_MAX) xmax = s.px.max;
   if (ymin == -_isis->DBL_MAX) ymin = s.py.min;
   if (ymax ==  _isis->DBL_MAX) ymax = s.py.max;

   xrange (xmin, xmax);
   yrange (ymin, ymax);

   return r;
}

%}}}

% Confidence limits for 2 degrees of freedom:
% 1-sigma (68.3%), 2-sigma (90%) and 3-sigma (99%)
private variable std_conf_limits = [2.30, 4.61, 9.21];

define plot_conf () %{{{
{
   variable msg = "plot_conf (Struct_Type[, line [, dchisqr_array]])";
   variable s, c, line;

   line = NULL;
   c = @std_conf_limits;

   if (_isis->get_varargs (&s, &line, &c, _NARGS, 1, msg))
     return;

   variable r = set_axis_ranges (s);
   plot_box ();
   plot_conf_contours (s, line, c);
   restore_axes (r);
}

%}}}

define oplot_conf () %{{{
{
   variable msg = "oplot_conf (Struct_Type [, line [, dchisqr_array]])";
   variable s, c, line;

   line = NULL;

   c = @std_conf_limits;

   if (_isis->get_varargs (&s, &line, &c, _NARGS, 1, msg))
     return;

   plot_conf_contours (s, line, c);
}

%}}}

% back-compatibility
alias ("plot_conf", "conf_plot");
alias ("oplot_conf", "conf_oplot");

%}}}
