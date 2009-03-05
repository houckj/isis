% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2008 Massachusetts Institute of Technology
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

% $Id: atom-cmds.sl,v 1.4 2004/02/09 11:14:14 houck Exp $

define atoms () %{{{
{
   variable msg = "atoms ([db_struct])";
   variable dir = NULL;

   if (_isis->get_varargs (&dir, _NARGS, 0, msg))
     return;
   
   if (dir == NULL)
     {
	usage (msg);
	return;
     }
   
   _isis->Dbase = dir;
   
   variable filemap = "";
   
   if (_isis->Dbase.atomic_data_filemap != NULL)
     filemap = path_concat (_isis->Dbase.dir, _isis->Dbase.atomic_data_filemap);
   
   _isis->_db_start (filemap);
}

%}}}

%{{{ energy levels

define list_elev ()
{
   variable msg = "list_elev (Z, ion)";
   variable Z, ion;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (Z, ion) = ();
   _isis->_list_elev (Z,ion-1);
}

define plot_elev ()
{
   variable msg ="plot_elev (Z, ion)";
   variable Z, ion;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (Z, ion) = ();

   _isis->_plot_elev (Z, ion-1,
	             -1,       % no line list
   	              0,       % plot all levels
	              0);      % not an overlay
}

define plot_elev_subset ()
{
   variable msg ="plot_elev_subset (Z, ion, line-list)";
   variable Z, ion, lines;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (Z, ion, lines) = ();

   _isis->_plot_elev (Z, ion-1, lines,
	              1,              % plot subset of levels
	              0);             % not an overlay
}

define oplot_elev_subset ()
{
   variable msg ="oplot_elev_subset (Z, ion, line-list)";
   variable Z, ion, lines;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (Z, ion, lines) = ();

   _isis->_plot_elev (Z, ion-1, lines,
	              1,               % plot subset of levels
	              1);              % overlay
}

define oplot_lines ()
{
   variable Z, ion, lines, style;
   variable msg ="oplot_lines (Z, ion, line-list, [style])";

   style = 0;

   if (_isis->get_varargs (&Z, &ion, &lines, &style, _NARGS, 3, msg))
     return;

   _isis->_oplot_transitions (Z, ion-1, lines, style);
}

%}}}

%{{{ line groups

define list_group ()
{
   _isis->_list_group ();
}

define name_group ()
{
   variable msg = "name_group (group_index, \"name\")";
   variable g, name;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (g, name) = ();
   _isis->_name_group (g,name);
}

define delete_group ()
{
   variable msg = "delete_group (group_index)";
   variable g;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   g = ();

   _isis->_delete_group (g);
}

define page_group ()
{
   variable msg = "page_group (group_index or line_list)";
   variable g;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   g = ();

   if (Array_Type == typeof (g))
     {
	_isis->_page_line_list (g);
     }
   else
     {
	_isis->_page_group (g);
     }
}

define brightest()
{
   variable k, list;
   variable msg = "brightest(k, line-list)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (k, list) = ();

   _isis->get_k_brightest_lines (k, list);
}

define unblended ()
{
   variable frac, wl_sep, type, list;
   variable msg =
     "unblended (frac, wl_sep, type, line-list)";

   if (_isis->chk_num_args (_NARGS, 4, msg))
     return;

   (frac, wl_sep, type, list) = ();

   _isis->get_unblended (frac, wl_sep, type, list);
}

define el_ion ()
{
   variable msg = "el_ion ([Z-list], [ion-list])";
   variable Z, ion;

   ion = -1;

   if (_isis->get_varargs (&Z, &ion, _NARGS, 1, msg))
     return;

   if (NULL == Z)
     Z = -1;
   if (NULL == ion)
     ion = -1;
   else
     ion = ion - 1;

   _isis->filter_line_by_el_ion (Z, ion);
}

define trans ()
{
   variable msg = "trans ([Z] [, ion[, upper[, lower]]])";
   variable Z, ion, upper, lower;

   Z = NULL;
   ion = NULL;
   upper = NULL;
   lower = NULL;

   if (_isis->get_varargs (&Z, &ion, &upper, &lower, _NARGS, 0, msg))
     return;
   
   if (Z == NULL)
     Z = -1;
   if (ion == NULL)
     ion = -1;
   else ion = ion -1;

   if (NULL == upper) upper = -1;
   else upper = upper[where(upper > 0)];
   
   if (NULL == lower) lower = -1;
   else lower = lower[where (lower > 0)];
   
   _isis->filter_line_by_trans (upper, lower, Z, ion);
}

private define filter_by_float_range (num, fcn, msg)
{
   variable fmin, fmax;
   variable HUGE = -1.0;   % interpreted as either FLT_MAX or DBL_MAX
   variable TINY = 0.0;

   fmin = TINY;
   fmax = HUGE;

   if (_isis->get_varargs (&fmin, &fmax, num, 1, msg))
     return;

   if (NULL == fmin)
     fmin = TINY;

   if (NULL == fmax)
     fmax = HUGE;

   @fcn (fmin, fmax);
}

define flx ()
{
   variable msg = "flx ([flux_min], [flux_max])";
   filter_by_float_range (_NARGS, &_isis->filter_line_by_flux, msg);
}

define wl ()
{
   variable msg = "wl ([lambda_min], [lambda_max])";
   filter_by_float_range (_NARGS, &_isis->filter_line_by_wavelength, msg);
}

define define_group ()
{
   variable msg = "define_group (group, index-list)";
   variable group, list;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (group, list) = ();
   _isis->make_group_from_list (group, list);
}

define save_group ()
{
   variable filename, group;
   variable msg = "save_group (group or line_list, \"filename\")";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (group, filename) = ();

   if (String_Type != typeof(filename))
     {
	usage(msg);
	return;
     }

   if (Array_Type == typeof (group))
     {
	_isis->_save_line_list (group, filename);
     }
   else
     {
	_isis->_save_group (group, filename);
     }
}

define line_label_default_style ()
{
   variable style = struct
     {
	angle, justify, top_frac, bottom_frac, offset, char_height, label_type
     };

   style.label_type = 0;     % short label
   style.justify = 0.0;      % left justified
   style.angle = 90.0;       % ccw degrees
   style.top_frac = 0.7;       % as fraction of plot Y limits
   style.bottom_frac = 0.6;
   style.offset = 0.25;        % text offset from indicator line
                               % as fraction of indicator line length
   style.char_height = 1.0;

   return style;
}

private define do_plot_line_group (group, redshift, s, gr_color, msg)
{
   if (Array_Type == typeof (group))
     {
	_isis->_plot_line_list (group, redshift, gr_color, 
			 s.angle, s.justify, s.top_frac, s.bottom_frac,
			 s.offset, s.char_height, s.label_type);
     }
   else 
     {
	_isis->plot_line_group (group, redshift, gr_color, 
			 s.angle, s.justify, s.top_frac, s.bottom_frac,
			 s.offset, s.char_height, s.label_type);
     }
}

define plot_group ()
{
   variable msg = "plot_group (group_id | line_list [, color, [, label_style [, redshift]]])";
   variable group, s, style, redshift, gr_color;

   redshift = 0.0;
   gr_color = 0;
   s = 0;

   if (_isis->get_varargs (&group, &gr_color, &s, &redshift, _NARGS, 1, msg))
     return;

   if (is_struct_type(s))
     {
	style = s;
     }
   else
     {
	style = line_label_default_style ();
	style.label_type = s;                    % back-compatibility
     }

   do_plot_line_group (group, redshift, style, gr_color, msg);
}

% default linelabel hook:
define isis_linelabel_hook (s)
{
   return latex2pg (s);
}

define plot_linelist ()
{
   variable msg = "plot_linelist (lambdas, labels [, color, [, label_style [, redshift]]])";
   variable lambdas, labels, label_style, s, redshift, gr_color;

   redshift = 0.0;
   gr_color = 0;
   label_style = 0;

   if (_isis->get_varargs (&lambdas, &labels, &gr_color, 
			   &label_style, &redshift, _NARGS, 2, msg))
     {
	return;
     }

   if (is_struct_type(label_style))
     {
	s = label_style;
     }
   else
     {
	s = line_label_default_style ();
	s.label_type = label_style;                    % back-compatibility
     }
   
   variable m_labels = array_map (String_Type, &isis_linelabel_hook, labels);

   _isis->_plot_line_list2 (lambdas, m_labels, redshift, gr_color, 
			    s.angle, s.justify, s.top_frac, s.bottom_frac,
			    s.offset, s.char_height, s.label_type);
}

define lines_in_group ()
{
   variable msg = "index[] = lines_in_group (group)";
   variable g;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   g = ();

   return _isis->_group_members (g);
}

%}}}

define change_wl () %{{{
{
   variable msg = "change_wl (\"filename\")";
   variable file;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   file = ();
   _isis->_change_wl (file);
}

%}}}

define line_info () %{{{
{
   variable msg = "Struct_Type s = line_info (line_index)";
   variable s, id;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;
   
   id = ();

   s = _isis->_get_line_info (id);
   if (s.Z == 0)
     return NULL;
   s.ion = s.ion + 1;
   
   return s;
}

%}}}

define list_branch () %{{{
{
   variable msg = "list_branch (proton_number, ion)";
   variable Z, ion;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (Z, ion) = ();
   _isis->_list_branching (Z, ion-1);
}

%}}}

