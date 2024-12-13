% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2025 Massachusetts Institute of Technology
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

% $Id: plot-cmds.sl,v 1.12 2004/09/09 17:50:51 houck Exp $

define plot_unit () %{{{
{
   variable msg = "plot_unit (x_unit_name)";
   variable name;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   name = ();
   _isis->_set_plot_units (name);
}

%}}}

define plot_bin_integral ()
{
   _isis->_plot_bin_integral ();
}

define plot_bin_density ()
{
   _isis->_plot_bin_density ();
}

%{{{ open/ resize

define plot_device ()
{
   variable msg = "plot_device (device_name)";
   variable device_string;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   device_string = ();
   _isis->_set_plot_device (device_string);
}

define plot_auto_color ()
{
   variable msg = "plot_auto_color (on=1/off=0)";
   variable flag;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   flag = ();
   _isis->_set_auto_increment (flag);
}

define label()
{
   variable msg = "label (\"x-label\", \"y-label\", \"title\")";
   variable xlabel, ylabel, title;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (xlabel, ylabel, title) = ();
   _isis->_set_labels (xlabel, ylabel, title);
}

define title()
{
   variable msg = "title (\"the title\")";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   _isis->_set_title ();
}

define xlabel()
{
   variable msg = "xlabel (\"the label\")";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   _isis->_set_xlabel ();
}

define ylabel()
{
   variable msg = "ylabel (\"the label\")";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   _isis->_set_ylabel ();
}

define xrange ()
{
   variable args = __pop_args (_NARGS);
   _isis->_set_xrange(__push_args(args));
}

define yrange ()
{
   variable args = __pop_args (_NARGS);
   _isis->_set_yrange(__push_args(args));
}

define limits ()
{
   variable args = __pop_args (_NARGS);
   _isis->_reset_plot_ranges(__push_args(args));
}

define xlog ()
{
   variable args = __pop_args (_NARGS);
   _isis->_logx_on(__push_args(args));
}

define xlin ()
{
   variable args = __pop_args (_NARGS);
   _isis->_logx_off(__push_args(args));
}

define ylog ()
{
   variable args = __pop_args (_NARGS);
   _isis->_logy_on(__push_args(args));
}

define ylin ()
{
   variable args = __pop_args (_NARGS);
  _isis->_logy_off(__push_args(args));
}

define close_plot ()
{
   variable args = __pop_args (_NARGS);
   _isis->_close_plot(__push_args(args));
}
alias ("close_plot", "plot_close");

define plot_quit ()
{
   variable args = __pop_args (_NARGS);
   _isis->_plot_quit(__push_args(args));
}

define erase ()
{
   variable args = __pop_args (_NARGS);
   _isis->_erase_plot(__push_args(args));
}

define clear ()
{
   variable args = __pop_args (_NARGS);
   _isis->_erase_plot(__push_args(args));
}

define cursor_box ()
{
   variable args = __pop_args (_NARGS);
   _isis->_get_cursor_box(__push_args(args));
}

define get_color ()
{
   variable args = __pop_args (_NARGS);
   _isis->_get_color(__push_args(args));
}

define window()
{
   variable msg = "window (pid)";
   variable pid;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   pid = ();
   _isis->_set_window (pid);
}

define mpane()
{
   variable msg = "mpane (pane)";
   variable pid;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   pid = ();
   pid -= 1;
   if (pid < 0) pid = 0;

   _isis->_set_pane (pid);
}

define line_or_color()
{
   variable msg = "line_or_color (1=colors/0=line-styles)";
   variable style;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   style = ();
   _isis->_set_style_meaning (style);
}

define color()
{
   variable msg = "color (color_index)";
   variable idx;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   idx = ();
   _isis->_set_color (idx);
}

define multiplot ()
{
   variable msg = "multiplot (relative_size_list)";
   variable ysizes = _isis->pop_list (_NARGS, msg);
   if (ysizes == NULL)
     return;

   _isis->_multiplot (ysizes);
}

define get_outer_viewport ()
{
   _isis->get_outer_viewport_intrin();
}

define set_outer_viewport ()
{
   variable msg = "set_outer_viewport (Struct_Type);";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   _isis->set_outer_viewport_intrin ();
}

define open_plot ()
{
   variable msg = "open_plot (device [, nxpanes, nypanes])";
   variable dev, nx, ny;

   dev = "";     % empty string implies use current ISIS default plot device
   nx = 1;
   ny = 1;

   if (_isis->get_varargs (&dev, &nx, &ny, _NARGS, 0, msg))
     return -1;

   if (String_Type != typeof(dev))
     {
	usage(msg);
	return -1;
     }

   return _isis->_open_plot (dev, nx, ny);
}
alias ("open_plot", "plot_open");

define dup_plot ()
{
   variable msg = "dup_plot ([device], [plot_id]);";
   variable dev, id;

   dev = "";
   id = -1;    % interpreted as current plot window

   if (_isis->get_varargs (&dev, &id, _NARGS, 0, msg))
     return -1;

   return _isis->_copy_plot (dev, id);
}
alias ("dup_plot", "plot_dup");

define resize ()
{
   variable width, aspect;
   variable msg = "resize ([width-cm], [aspect])";

   width = 0.0;
   aspect = 0.618;

   if (_isis->get_varargs (&width, &aspect, _NARGS, 0, msg))
     return;

   if (NULL == width)
     width = 0.0;

   if (NULL == aspect)
     aspect = 0.618;

   _isis->_resize_window (width, aspect);
}

define get_plot_info () %{{{
{
   variable p = struct
     {
	xmin, xmax,
	ymin, ymax,
        xlog, ylog,
	line_width,
	line_color
     };

   (p.xmin, p.xmax) = _isis->_get_xrange ();
   (p.ymin, p.ymax) = _isis->_get_yrange ();
   p.line_width = _isis->_get_line_width ();
   p.line_color = _isis->_get_color ();
   p.xlog = _isis->_xaxis_is_log ();
   p.ylog = _isis->_yaxis_is_log ();

   return p;
}

%}}}

%}}}

% plot diffraction orders
define lambda_mth_order()%{{{
{
   variable msg = "lambda_mth_order(m, lambda_n, [n])";
   variable m, n, lambda_n;

   n = 1.0;

   if (_isis->get_varargs (&m, &lambda_n, &n, _NARGS, 2, msg))
     return;

   _isis->_lambda_mth_order (m, lambda_n, n);
}

%}}}

%{{{ xylabel/ cursor

define xylabel ()
{
   variable msg = "xylabel (x, y, text [, angle [, justify]])";
   variable x,y,text, angle, justify;

   angle = 0.0;
   justify = 0.0;

   if (_isis->get_varargs (&x, &y, &text, &angle, &justify, _NARGS, 3, msg))
     return;

   _isis->_put_text( x, y, text, angle, justify );
}

define xinterval ()
{
   return _isis->_xinterval ();      % returns (xmin, xmax)
}

define yinterval ()
{
   return _isis->_yinterval ();      % returns (ymin, ymax)
}

define charsize ()
{
   variable msg = "charsize (height)";
   variable height;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   height = ();
   _isis->_set_charsize (height);
}

%}}}

%{{{ pointstyle/ linestyle/ points

define pointstyle () %{{{
{
   variable idx;

   if (_NARGS != 1)
     {
        message("pointstyle (idx);\n");
        message("Where idx is one of the following:");
        message("   -1, -2    a single dot");
        message("  -3..-31    a regular polygon with abs(idx) edges");
        message("    0..31    standard marker symbols");
        message("  32..127    ASCII characters (in current font)");
        message("    > 127    a Hershey symbol number\n");
        _pop_n(_NARGS);
        return;
     }

     idx = ();

     _isis->_set_pointstyle (idx);
}
alias ("pointstyle", "point_style");
%}}}

define point_size ()
{
   variable msg = "point_size (size)";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;
   
   _isis->_set_point_size ();
}

private define match_linestyle (s) %{{{
{
   switch (s)
     {
      case "-" or   % dashed
      case "--":
	return 2;
     }
     {
      case ".-" or  % dot-dash
      case "-.":
	return 3;
     }
     {
      case "." or   % dotted
      case ".." or
      case "...":
	return 4;
     }
     {
      case "..-" or  % dash-dot-dot-dot
      case "-.." or
      case "...-" or
      case "-...":
	return 5;
     }
     {
	% default  (full line)
	return 1;
     }
}

%}}}

define linestyle () %{{{
{
   variable idx;

   if (_NARGS != 1)
     {
        message ("linestyle (idx);\n");
        message ("Where idx is one of:");
        message (" 1 = full line");
	message (" 2 = dashed");
	message (" 3 = dot-dash");
	message (" 4 = dotted");
	message (" 5 = dash-dot-dot-dot\n");
        _pop_n(_NARGS);
        return;
     }

   idx = ();

   if (typeof(idx) == String_Type)
     idx = match_linestyle (idx);

   _isis->_set_linestyle (idx);
}
alias ("linestyle", "line_style");
%}}}

define set_frame_line_width () %{{{
{
   variable msg = "set_frame_line_width (width)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   _isis->_set_frame_line_width ();
}

%}}}

define set_line_width () %{{{
{
   variable msg = "set_line_width (width)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   _isis->_set_line_width ();
}

%}}}

private variable X_Unit = Assoc_Type[Int_Type];
X_Unit["eV"] = U_EV;
X_Unit["keV"] = U_KEV;
X_Unit["MeV"] = U_MEV;
X_Unit["GeV"] = U_GEV;
X_Unit["TeV"] = U_TEV;
X_Unit["Hz"] = U_HZ;
X_Unit["kHz"] = U_KHZ;
X_Unit["MHz"] = U_MHZ;
X_Unit["GHz"] = U_GHZ;
X_Unit["A"] = U_ANGSTROM;
X_Unit["nm"] = U_NANOMETER;
X_Unit["mm"] = U_MILLIMETER;
X_Unit["cm"] = U_CENTIMETER;
X_Unit["m"] = U_METER;

define get_plot_options ()
{
   variable o = _isis->_get_plot_options_struct();
   % x_unit always returned as an integer.
   variable i = where (assoc_get_values (X_Unit) == o.x_unit)[0];
   o.x_unit = assoc_get_keys (X_Unit)[i];
   return o;
}

define set_plot_options () %{{{
{
   variable msg = "set_plot_options (Struct_Type)";
   variable s;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   s = ();

   if (typeof(s.x_unit) == String_Type)
     {
        !if (assoc_key_exists (X_Unit, s.x_unit))
          {
             verror ("unsupported unit %s", s.x_unit);
             return;
          }
        s.x_unit = X_Unit[s.x_unit];
     }   
   
   _isis->_set_plot_options_struct (s);
   
}

%}}}

define connect_points () %{{{
{
   variable msg = "connect_points (on=1/off=0)";
   variable flag;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   flag = ();

   _isis->_set_connect_points (flag);
}

%}}}

%}}}

%{{{ [o]plot (x,y)

private define get_plot_args (num)
{
   variable x, y, style, symbol;

   x = NULL;
   y = NULL;
   style = _isis->_INT_MAX;      % use default style
   symbol = [0];

   switch (num)
     {
      case 3:
	(x, y, style) = ();
     }
     {
      case 2:
	(x, y) = ();
     }
     {
	% default:
	_pop_n(num);
     }

   if (Array_Type == typeof(style))
     {
	symbol = style;
	style = _isis->_INT_MAX;
     }

   return (x, y, style, symbol);
}

private define __plot_xy (x, y, style, symbol, overlay, msg)
{
   if (NULL == x or NULL == y)
     {
	usage(msg);
	return;
     }

   _isis->_plot_xy (x,y, style, overlay, symbol);
}

define plot_box ()
{
   _isis->_plot_box();
}

define plot()
{
   variable msg = "plot (x[], y[] [, linestyle or symbol[] ])";
   variable x, y, style, symbol;

   (x, y, style, symbol) = get_plot_args(_NARGS);

   __plot_xy (x, y, style, symbol, 0, msg);  % no overlay
}

define oplot()
{
   variable msg = "oplot (x[], y[] [, linestyle or symbol[] ])";
   variable x, y, style, symbol;

   (x, y, style, symbol) = get_plot_args(_NARGS);

   __plot_xy (x, y, style, symbol, 1, msg);  % overlay
}

%}}}

%{{{ [o]hplot (lo,hi,y)

private define get_hplot_args (num)
{
   variable lo, hi, y, style, s;

   s = NULL;
   lo = NULL;
   hi = NULL;
   y  = NULL;
   style = _isis->_INT_MAX;      % use default style

   switch (num)
     {
      case 4:
	(lo, hi, y, style) = ();
     }
     {
      case 3:
	(lo, hi, y) = ();
     }
     {
      case 2:
	(s, style) = ();
     }
     {
      case 1:
	s = ();
     }
     {
	% default:
	_pop_n(num);
        return (lo, hi, y, style);
     }

   if (typeof(s) == Struct_Type)
     {
	lo = s.bin_lo;
	hi = s.bin_hi;
	y = s.value;
     }

   % Sort into ascending order works around a bug in pgbin.
   % When the X coordinates are in descending order, pgbin
   % gets slightly confused.
   variable i = array_sort (lo);   
   return (lo[i], hi[i], y[i], style);
}

private define __hplot (lo, hi, y, style, overlay, msg)
{
   if (NULL == lo or NULL == hi or NULL == y)
     {
	usage(msg);
	return;
     }

   _isis->_hplot (lo, hi, y, style, overlay);
}

define hplot()
{
   variable msg = "hplot (lo, hi, y [, linestyle])";
   variable lo, hi, y, style;

   (lo, hi, y, style) = get_hplot_args(_NARGS);

   __hplot (lo, hi, y, style, 0, msg);  % no overlay
}

define ohplot()
{
   variable msg = "ohplot (lo, hi, y [, linestyle])";
   variable lo, hi, y, style;

   (lo, hi, y, style) = get_hplot_args(_NARGS);

   __hplot (lo, hi, y, style, 1, msg);  % overlay
}

%}}}

% transform latex string format to PGPLOT format

%{{{ support for PGPLOT Greek characters

private variable pg_greek = Assoc_Type [];
pg_greek ["alpha"]   = "a";
pg_greek ["beta"]    = "b";
pg_greek ["gamma"]   = "g";
pg_greek ["delta"]   = "d";
pg_greek ["epsilon"] = "e";
pg_greek ["zeta"]    = "z";
pg_greek ["eta"]     = "y";
pg_greek ["theta"]   = "h";
pg_greek ["iota"]    = "i";
pg_greek ["kappa"]   = "k";
pg_greek ["lambda"]  = "l";
pg_greek ["mu"]      = "m";
pg_greek ["nu"]      = "n";
pg_greek ["xi"]      = "c";
pg_greek ["omicron"] = "o";
pg_greek ["pi"]      = "p";
pg_greek ["rho"]     = "r";
pg_greek ["sigma"]   = "s";
pg_greek ["tau"]     = "t";
pg_greek ["upsilon"] = "u";
pg_greek ["phi"]     = "f";
pg_greek ["chi"]     = "x";
pg_greek ["psi"]     = "q";
pg_greek ["omega"]   = "w";

private define strup1 (s) %{{{
{
   return strsub (s, 1, toupper(s[0]));
}

%}}}

private define greek_letters () %{{{
{
   variable k, keys = assoc_get_keys (pg_greek);
   variable g = Assoc_Type[];

   foreach (keys)
     {
	k = ();
	g["\\" +        k ] = "\\g" +       pg_greek[k];
	g["\\" + strup1(k)] = "\\g" + strup(pg_greek[k]);
     }

   return g;
}

%}}}

private variable greek = greek_letters ();

private define greek_substring (s) %{{{
{
   variable matched;
   matched = string_match(s, "\\\\[a-zA-Z]+", 1);
   if (matched == 0)
     return NULL;

   variable len, g;
   (, len) = string_match_nth (0);
   g = s[[0:len-1]];

   if (0 == assoc_key_exists (greek, g))
     return NULL;

   return g;
}

%}}}

private define replace_greek_letter (s, g) %{{{
{
   variable repl = greek [g];
   (s.string, ) = strreplace (s.string, g, repl, 1);
   s.imax = strlen (s.string);
   s.i += strlen(repl) - 1;
   return repl;
}

%}}}

private define next_char (s);
private define parse_greek (s) %{{{
{
   variable g = greek_substring (s.string[[s.i-1:]]);
   if (g == NULL)
     return char(s.string[s.i-1]);
   return replace_greek_letter (s, g);
}

%}}}

%}}}

%{{{ latex string parser

%  This is a recursive decent parser which makes the
%  following string transformations:
%
%     '~'  --> ' '
%     ^X   --> \\uX\\d
%     _{Y} --> \\dY\\u
%     ^{Y} --> \\uY\\d
% \\alpha  --> \\ga
%      .. etc. for other Greek characters
%
%  It is intended to transform pseudo-LaTeX-formatted
%  label strings from the APED database into plottable
%  strings understood by PGPLOT.
%
%  The original version written by John Davis did
%  not contain support for Greek characters.
%  JH added that.
%
private define next_char (s) %{{{
{
   if (s.i == s.imax)
     return 0;
   variable ch = s.string[s.i];

   s.i++;
   return ch;
}

%}}}

private define parse_block ();
private define parse_substring ();
private define parse_one_char (ch, s) %{{{
{
   switch (ch)
     {
	case '^':
	  return "\\u" + parse_block (s) + "\\d";
     }
     {
	case '_':
	  return "\\d" + parse_block (s) + "\\u";
     }
     {
        case '\\':
   	  return parse_greek (s);
     }
     {
	  return char (ch);
     }
}

%}}}

private define parse_block (s) %{{{
{
   variable ch = next_char (s);
   if (ch == '{')
     return parse_substring (s, '}');

   return parse_one_char (ch, s);
}

%}}}

private define parse_substring (s, end_char) %{{{
{
   variable t = "";

   forever
     {
	  variable ch = next_char (s);
	  if ((ch == 0) or (ch == end_char))
	    return t;

	  t += parse_one_char (ch, s);
     }
}

%}}}

define latex2pg (str) %{{{
{
   (str, ) = strreplace (str, "~", " ", strlen(str));

   variable s = struct
     {
	  string, imax, i
     };
   s.string = str;
   s.imax = strlen (str);
   s.i = 0;

   return parse_substring (s, 0);
}

%}}}

%}}}

