% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2017 Massachusetts Institute of Technology
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

% $Id: model-cmds.sl,v 1.4 2004/02/09 11:14:16 houck Exp $

define edit_model () %{{{
{
   variable msg = "edit_model ([\"file\"])";
   variable file="";

   if (_isis->get_varargs (&file, _NARGS, 0, msg))
     return;

   _isis->_edit_model (file);
}

%}}}

define list_model () %{{{
{
   _isis->_list_model ();
}

%}}}

private define do_model_file (nargs, fun) %{{{
{
   variable filename;

   if (nargs != 1)
     {
	_pop_n (nargs);
	return -1;
     }

   filename = ();

   if (String_Type != typeof(filename))
     return -1;

   return @fun(filename);
}

%}}}

define load_model () %{{{
{
   variable msg = "load_model (\"filename\")";
   variable ret;

   ret = do_model_file (_NARGS, &_isis->_load_ascii_model);

   if (ret)
     usage (msg);
}

%}}}

define save_model () %{{{
{
   variable msg = "save_model (\"filename\")";
   variable ret;

   ret = do_model_file (_NARGS, &_isis->_save_ascii_model);

   if (ret)
     usage (msg);
}

%}}}

define use_thermal_profile () %{{{
{
   _isis->_set_model_profile (1);
}

%}}}

define use_delta_profile () %{{{
{
   _isis->_set_model_profile (0);
}

%}}}

define model_spectrum () %{{{
{
   variable msg = "emis = model_spectrum (binlo, binhi [,contrib_flag [, line_list]])";
   variable lo, hi, list, flag;

   list = NULL;
   flag = MODEL_LINES_AND_CONTINUUM;

   if (_isis->get_varargs (&lo, &hi, &flag, &list, _NARGS, 2, msg))
     return;

   if (length(lo) != length(hi))
     {
	message ("Invalid grid:  grid arrays should be the same size");
	return Null_Type;
     }

   if (list == NULL)
     list = [-1];

   if (length(where((lo-hi) >= 0.0)) > 0)
     {
	message ("Invalid grid:  grid should be in ascending order, with lo[i] < hi[i]");
	return Null_Type;
     }

   variable s = struct
     {
        contrib_flag, line_list
     };
   s.line_list = list;
   s.contrib_flag = flag;

   return _isis->_calc_model_list (s, "aped_hook", lo, hi);
}

%}}}

define default_plasma_state () %{{{
{
   variable state = struct
     {
	norm, temperature, density, metal_abund,
	  elem_abund, elem,
	  vturb, redshift
     };

   state.norm = 1.0;           % (Emission measure) / (4*PI*D^2)
   state.temperature = 1.0e7;  % Kelvin
   state.density = 1.0;        % cm^-3
   state.metal_abund = 1.0;    % relative to solar
   state.elem = NULL;          % element Z
   state.elem_abund = NULL;    % element abundance relative to solar
   state.vturb = 0.0;          % cm/s
   state.redshift = 0.0;

   return state;
}

%}}}

define define_model () %{{{
{
   variable msg = "define_model (state[])";
   variable s, state;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   state = ();

   variable init = 1;
   foreach (state)
     {
	s = ();
	_isis->_add_component (s, s.elem_abund, s.elem, init);
	init = 0;
     }
}

%}}}

define append_model () %{{{
{
   variable msg = "append_model (state[])";
   variable s, state;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   state = ();

   foreach (state)
     {
	s = ();
	_isis->_add_component (s, s.elem_abund, s.elem, 0);
     }
}

%}}}

define make_hi_grid () %{{{
{
   variable msg = "hi = make_hi_grid (lo)";
   variable x, n, x_last;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   x = ();
   n = length(x);
   if (n < 2)
     return NULL;

   x_last = 2*x[n-1] - x[n-2];

   return [x[[1:n-1]], x_last];
}

%}}}

define linear_grid () %{{{
{
   variable msg = "(lo, hi) = linear_grid (min, max, nbins)";
   variable n, min, max, lo, hi, dx;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (min, max, n) = ();

   dx = (max - min) / double(n);
   n = int(n);

   if (n <= 0)
     {
	usage (msg);
	return;
     }

   hi = Double_Type [n];
   lo = Double_Type [n];

   lo = min + dx * [0:n-1];
   hi[[0:n-2]] = lo[[1:n-1]];
   hi[n-1]   = max;

   return (lo, hi);
}

%}}}

define mt_def_model () %{{{
{
   variable msg = "Model_Type = mt_def_model (state[])";
   variable m, s, state;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   state = ();

   m = NULL;
   foreach (state)
     {
	s = ();
	m = _isis->_cl_add_component (s, s.elem_abund, s.elem, m);
     }

   return m;
}

%}}}

private define make_model_details_list () %{{{
{
   variable 
     args = __pop_args (_NARGS),
     num_args = length(args);

   variable l = list_new ();

   variable i = 0;
   while (i < num_args)
     {
        variable s = args[i].value;
        s = struct_combine (s, "line_flux");
        s.line_flux = args[i+1].value;
        list_append (l, s);
        i += 2;
     }

   return l;
}

%}}}

define mt_get_model_details (m) %{{{
{
   return make_model_details_list (_isis->_cl_get_model_details (m));
}

%}}}

define mt_calc_model () %{{{
{
   variable msg =
      "\nemis = mt_calc_model (Model_Type, lo, hi[, args])\n"
    + "   OR\n"
    + "emis = mt_calc_model (Model_Type, lo, hi[, contrib_flag [, line_list]])";
   variable model, lo, hi;

   variable args = __pop_args (_NARGS);

   if (_NARGS < 3)
     {
        usage(msg);
        return;
     }

   model = args[0].value;
   lo = args[1].value;
   hi = args[2].value;

   if (typeof (model) != Model_Type)
     {
	usage (msg);
	return Null_Type;
     }

   if (length (lo) != length (hi))
     {
	message ("Invalid grid:  grid arrays should be the same size");
	return Null_Type;
     }

   if (length(where((lo-hi) >= 0.0)) > 0)
     {
	message ("Invalid grid:  grid should be in ascending order, with lo[i] < hi[i]");
	return Null_Type;
     }

   % Support old interface for back-compatibility:
   if (_NARGS == 4 || _NARGS == 5)
     {
        if (typeof(args[3].value) == Int_Type)
          {
             variable s = struct {contrib_flag, line_list};
             s.contrib_flag = args[3].value;
             if (_NARGS == 5) s.line_list = args[4].value;
             s, "aped_hook";
          }
        else __push_args (args[[3:]]);
     }

   return _isis->_cl_calc_model (lo, hi, model);
}

%}}}

define mt_load_model () %{{{
{
   variable msg = "m = mt_load_model (file)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable file = ();
   _isis->_cl_load_ascii_model (file);
}

%}}}

define mt_save_model () %{{{
{
   variable msg = "mt_save_model (Model_Type, file)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable mt, file;
   (mt, file) = ();

   _isis->_cl_save_ascii_model (mt, file);
}

%}}}

define mt_list_model () %{{{
{
   variable msg = "mt_list_model (Model_Type)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable mt = ();

   _isis->_cl_list_model (mt);
}

%}}}

define load_line_profile_function () %{{{
{
   variable msg = "lp = load_line_profile_function (file, name)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   _isis->load_line_profile_function_intrin ();
}

%}}}
