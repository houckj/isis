% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2010 Massachusetts Institute of Technology
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
% $Id: fit-cmds.sl,v 1.35 2004/09/09 17:50:51 houck Exp $

%{{{ verbose hook

define chisqr_report_hook ()
{
}

define open_fit_verbose_hook ()
{
%   message ("called open_fit_verbose_hook\n");
}

define close_fit_verbose_hook ()
{
%   message ("called close_fit_verbose_hook\n");
}

define fit_verbose_warn_hook ()
{
   variable s = ();
   () = fputs (s, stderr);
}

define fit_verbose_info_hook () %{{{
{
   variable chisqr, params, param_names;
   (chisqr, params, param_names) = ();

   () = fprintf (stdout, "%s=%11.4e:\n", Fit_Statistic, chisqr);

   variable n = length(params);
   variable i;

   _for (0, n-1, 1)
     {
	i = ();
	() = fprintf (stdout, " %s= %15.8e", param_names[i], params[i]);
	if (0 == ((i+1) mod 2))
	  () = fputs ("\n", stdout);
     }

   () = fputs ("\n", stdout);
}

%}}}

%}}}

%{{{ user-defined fit functions

define list_functions ()
{
   _isis->list_functions_hook (_isis->_function_list ());
}

define del_function () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "del_function (\"fun_name\")";
   variable f_name;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   f_name = ();
   _isis->_del_function (f_name);
}

%}}}

private define parse_parameter_units (names)
{
   variable i, n, t, s;

   if (typeof(names) == String_Type)
     {
	names = [names];
     }

   n = length(names);
   if (n == 0)
     return names;

   variable units = String_Type[n];
   foreach ([0:n-1])
     {
	i = ();
	s = names[i];

	!if (is_substr (s, "["))
	  {
	     units[i] = "";
	  }
	else
	  {
	     t = strtok (s, "[]");
	     if (length(t) != 2)
	       return NULL;

	     names[i] = strtrim(t[0]);
	     units[i] = strtrim(t[1]);
	  }
     }

   return units;
}

define add_compiled_function () %{{{
{
   variable msg = "add_compiled_function (lib_name, function_name [, option_string])";
   variable lib, f_name, option_string = "", dummy_units = "";

   if (_isis->get_varargs (&lib, &f_name, &option_string, _NARGS, 2, msg))
     return;

   if (String_Type != typeof(lib)
       or String_Type != typeof(f_name))
     {
	message ("Invalid arguments");
	return;
     }

   variable path = find_library_name (lib);
   if (path == NULL)
     {
	vmessage ("File not found:  %s", lib);
	return;
     }

   _isis->_add_cfun_intrin (path, f_name, option_string, dummy_units, 0);
}

%}}}

define add_slang_function () %{{{
{
   variable msg =
     "add_slang_function (name [, parameter_names [, norm_indices]])\n"
    +"       add_slang_function (name, ref, [, parameter_names [, norm_indices]])";
   variable f_name, f = NULL;
   variable norm_indices = NULL;
   variable parameter_names = "";
   variable parameter_units = "";
   variable a1, a2;

   switch (_NARGS)
     {
      case 1:
        f_name = ();
     }
     {
      case 2:
        (f_name, a1) = ();
        if (_typeof(a1) == Ref_Type)
          {
             f = a1;
          }
        else
          {
             parameter_names = a1;
          }
     }
     {
      case 3:
        (f_name, a1, a2) = ();
        if (_typeof (a1) == Ref_Type)
          {
             f = a1;
             parameter_names = a2;
          }
        else
          {
             parameter_names = a1;
             norm_indices = a2;
          }
     }
     {
      case 4:
        (f_name, f, parameter_names, norm_indices) = ();
     }
     {
        % default:
        _pop_n (_NARGS);
        usage (msg);
        return;
     }

   if (String_Type != _typeof(parameter_names)
       or String_Type != typeof(f_name))
     {
	message ("Invalid arguments");
	return;
     }

   if (f == NULL)
     {
        variable fit_func = f_name + "_fit";
        f = __get_reference (fit_func);
        if (f == NULL)
          {
             vmessage ("Function '%s' not found --> Required by new fit-function '%s'",
                       fit_func, f_name);
             return;
          }
     }

   parameter_units = parse_parameter_units (parameter_names);
   if (parameter_units == NULL)
     {
	vmessage ("Error parsing parameter units:  %s", f_name);
	return;
     }

   if (norm_indices == NULL)
     {
        % out-of-range norm index indicates no norms were provided
        norm_indices = [length(parameter_names)];
     }
   else if (any(norm_indices < 0 or length(parameter_names) <= norm_indices))
     {
	vmessage ("Norm indices must be numbered from 0 <= n < N-1");
	return;
     }

   if (typeof(f) == Array_Type)
     {
        _isis->_add_slangfun_intrin (f[0], f[1],
                                     f_name, parameter_names, parameter_units,
                                     norm_indices);
     }
   else
     {
        _isis->_add_slangfun_intrin (f, NULL,
                                     f_name, parameter_names, parameter_units,
                                     norm_indices);
     }
}

%}}}

define set_param_default_hook () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_param_default_hook (slangfun_name, hook [, args])";
   variable name, hook, args=NULL;

   if (_NARGS < 2)
     {
        usage(msg);
        return;
     }
   else if (_NARGS > 2)
     {
        args = __pop_args (_NARGS-2);
     }

   (name, hook) = ();

   if (typeof(hook) == String_Type)
     hook = __get_reference (hook);

   if (args == NULL)
     _isis->_set_slangfun_param_default_hook (hook, name);
   else
     _isis->_set_slangfun_param_default_hook (__push_args(args), hook, name);
}

%}}}

define set_function_category () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_function_category (name, category);\n   e.g. category = ISIS_FUN_ADDMUL | ISIS_FUN_OPERATOR";
   variable name, category;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (name, category) = ();
   _isis->_set_function_category (name, category);
}

%}}}

%}}}

% freeze / thaw / tie / untie fit params

define get_params();
define set_params();

%{{{ expand parameter name expressions

% Parameter names can be specified by
% - integer indices
% - strings
% - slang regular expressions
% - globbing expressions

private define par_match (pat) %{{{
{
   variable parms = get_params ();
   variable i, n = length (parms);
   variable ok = Int_Type[n]-1;
   _for i (0, n-1, 1)
     {
        variable p = parms[i];
        if (string_match (p.name, pat, 1))
          ok[i] = p.index;
     }
   return ok[where(ok != -1)];
}

%}}}

define get_num_pars();
define get_par_info();

define _get_index (pat) %{{{
{
   if (typeof (pat) != String_Type)
     {
        % integer?
        return pat;
     }

   pat = strtrans (pat, " \t\n", "");

   variable n;
   if (pat[0] == '^' or pat[-1] == '$')
     {
        % regular expression?
        n = par_match (pat);
        if (length(n) != 0)
          return n;
     }
   else if (string_match (pat, "[[^*?]", 1))
     {
        % globbing expression?
#ifexists glob_to_regexp
        n = par_match (glob_to_regexp (strtrim(pat, "$^")));
        if (length(n) != 0)
          return n;
#endif
     }

   % string?
   return _isis->_get_index_for_param_name (pat);
}

%}}}

private define unary_param_op (fun, num, msg) %{{{
{
   variable a = _isis->pop_list (num, msg);
   if (a == NULL)
     return;

   variable a1, i = Integer_Type[0];
   foreach a1 ([a])
     {
        i = [i, _get_index(a1)];
     }

   foreach (i)
     {
	(@fun) ();
     }
}

%}}}

define freeze()
{
   _isis->error_if_fit_in_progress (_function_name);
   unary_param_op (&_isis->_freeze, _NARGS, "freeze (par-list)");
}

define thaw()
{
   _isis->error_if_fit_in_progress (_function_name);
   unary_param_op (&_isis->_thaw, _NARGS, "thaw (par-list)");
}

define untie()
{
   _isis->error_if_fit_in_progress (_function_name);
   unary_param_op (&_isis->_untie, _NARGS, "untie (par-list)");
}

define tie() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "tie (par, par-list)";
   variable to, args;

   if (_NARGS > 2)
     {
	args = __pop_args (_NARGS-1);
	args = [__push_args (args)];
	args = array_map (Integer_Type, &_get_index, args);
     }
   else
     {
	args = ();
	args = _get_index (args);
     }

   to = ();

   if (typeof(to) == Array_Type)
     {
	usage (msg);
	return;
     }

   to = _get_index(to);

   variable i;
   foreach (args)
     {
	i = ();
	_isis->_tie (i, to);
     }
}

%}}}

%}}}

%{{{ set / save / load fit params

define get_num_pars ()
{
   return _isis->_get_num_params ();
}

private define do_get_param_info (idx)
{
   variable s = _isis->_get_param_info (idx);
   if (s.fun == "")
     s.fun = NULL;
   return s;
}

define get_par_info () %{{{
{
   variable msg = "struct = get_par_info (idx)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable idx = ();

   idx = _get_index (idx);

   if (typeof (idx) == Array_Type)
     {
        return array_map (Struct_Type, &do_get_param_info, idx);
     }
   else if (idx > 0)
     {
        return do_get_param_info (idx);
     }

   return NULL;
}

%}}}

define get_par () %{{{
{
   variable msg = "value = get_par (idx)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable idx = ();

   idx = _get_index (idx);

   if (typeof (idx) == Array_Type)
     return array_map (Double_Type, &_isis->_get_par, idx);
   else
     return _isis->_get_par (idx);
}

%}}}

define get_fit_fun () %{{{
{
   return _isis->_get_fit_fun();
}

%}}}

define __set_hard_limits ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg =
"__set_hard_limits (fun_name, par_name, hard_min, hard_max);\n\
  or\n\
__set_hard_limits (index, hard_min, hard_max);\n\
          Qualifier:  parm=[min, max, value]";

   variable fun_name = "", par_name = "", idx = -1,
     hard_min, hard_max;

   switch (_NARGS)
     {
      case 4:
        (fun_name, par_name, hard_min, hard_max) = ();
     }
     {
      case 3:
        (idx, hard_min, hard_max) = ();
        idx = _get_index (idx);
     }
     {
        %default:
        _pop_n(_NARGS);
        message(msg);
        return;
     }

   variable p = struct
     {
        hard_min = hard_min,
        hard_max = hard_max,
        min = _NaN,
        max = _NaN,
        value = _NaN,
        step = _NaN, % unused
        freeze = 0   % unused
     };

   variable
     hard_limits_only,
     parm = qualifier ("parm", NULL);

   if (parm == NULL)
     hard_limits_only = 1;
   else if (typeof(parm) == Array_Type && length(parm) == 3)
     {
        hard_limits_only = 0;
        variable i = array_sort (parm);
        p.min = parm[i[0]];
        p.value = parm[i[1]];
        p.max = parm[i[2]];
     }
   else
     {
        message(msg);
        return;
     }

   if (_isis->set_hard_limits (p, fun_name, par_name, idx, hard_limits_only) != 0)
     throw IsisError;
}

define get_fun_info () %{{{
{
   variable msg = "Struct_Type = get_fun_info (\"fun_name\")";
   variable fun;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   fun = ();

   variable info = _isis->get_fitfun_info (fun);
   if (info.name == NULL)
     return NULL;

   return info;
}

%}}}

private define uniq_strings (a) %{{{
{
   variable last, n, num, uniq;

   a = a[ array_sort (a, &strcmp) ];

   num = length(a);
   uniq = Int_Type[num];

   last = 0;
   uniq[0] = 0;
   n = 1;

   foreach ([1:num-1])
     {
	variable k = ();
	if (a[k] != a[last])
	  {
	     uniq[n] = k;
	     n++;
	     last = k;
	  }
     }

   a[uniq[[0:n-1]]];
}

%}}}

define get_fun_components () %{{{
{
   variable dummy_result, num, x;
   variable msg = "names[] = get_fun_components ()";

   if (NULL == get_fit_fun())
     return NULL;

   _isis->Model_Component_Names = NULL;
   _isis->_do_name_eval_mode ();

   num = 0;
   foreach (_isis->Model_Component_Names) using ("next")
     {
	x = ();
	x.name;
	num++;
     }

   if (num <= 1)
     return;

   uniq_strings ( _isis->pop_list (num, msg));
}

%}}}

define get_fun_params () %{{{
{
   variable msg = "id = get_fun_params (\"fun_name\")";
   variable fun, n, idx, len;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   fun = ();

   fun = strtrans (fun, " \t\n", "");

   len = strlen(fun);
   n = get_num_pars();
   idx = Integer_Type[n];

   variable i, numf = 0;

   _for (1, n, 1)
     {
	i = ();
	if (0 == strncmp (get_par_info (i).name, fun, len))
	  {
	     idx[numf] = i;
	     numf++;
	  }
     }

   if (numf > 0)
     return idx[[0:numf-1]];
   else
     return NULL;
}

%}}}

define _par (id) %{{{
{
   id = _get_index (id);
   return get_par (id);
}

%}}}

define set_par_fun () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_par_fun (id, function_string)";
   variable id, s;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (id, s) = ();

   if (s == NULL)
     s = "";

   id = _get_index (id);
   _isis->_set_par_fun (id, s);
}

%}}}

define set_par() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_par (idx, value [,freeze, [ min, max [, step]]])";
   variable idx, freeze, value, update_minmax;
   variable par_min, par_max, preserve_tie, step;

   freeze = -1;        % default means dont change current setting
   preserve_tie = -1;  % default means dont change current setting
   update_minmax = 0;
   par_min = 0;
   par_max = 0;
   step = -1.0;       % default means don't change current setting

   switch (_NARGS)
     {
      case 2:
	(idx, value) = ();
     }
     {
      case 3:
	(idx, value, freeze) = ();
     }
     {
      case 5:
	(idx, value, freeze, par_min, par_max) = ();
	update_minmax = 1;
     }
     {
      case 6:
	(idx, value, freeze, par_min, par_max, step) = ();
	update_minmax = 1;
     }
     {
	% default:
	_pop_n(_NARGS);
	usage(msg);
	return;
     }

   variable id = _get_index (idx);

   if (andelse
       {typeof(idx) == String_Type}
       {typeof (value) == Array_Type}
       {is_substr(idx, ".")}
       {length(id) != length(value)}
      )
     {
        vmessage ("*** parameter name/value mismatch");
        return;
     }

   if (andelse
       {typeof(id) == Integer_Type}
       {id < 0})
     {
	id = get_fun_params(idx);
	if (id == NULL)
	  {
	     verror ("*** invalid parameter list %S", idx);
	     return;
	  }
     }

   if (orelse
       {typeof (value) == String_Type}
       {value == NULL})
     {
	set_par_fun (id, value);
	_isis->_set_par (id, preserve_tie, freeze, get_par(id), par_min, par_max, update_minmax, step);
     }
   else
     {
	array_map (Void_Type, &_isis->_set_par, id, preserve_tie,
		   freeze, value, par_min, par_max,
		   update_minmax, step);
     }
}

%}}}

define load_par() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "load_par (\"filename\")";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;
   variable fname = ();
   () = _isis->_load_par (fname);
}

%}}}

define edit_par () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "edit_par ([\"file\"])";
   variable file="";

   if (_isis->get_varargs (&file, _NARGS, 0, msg))
     return;

   _isis->_edit_par (file);
}

%}}}

private variable Fmt_Colheads, Fmt_Params;
Fmt_Colheads = " idx%-*stie-to  freeze         value         min         max";
Fmt_Params   = "%3u  %-*s  %2d    %2d   %14.7g  %10.7g  %10.7g  %s%s";

private define param_string (p, len) %{{{
{
   return sprintf (Fmt_Params, p.index, len, p.name,
                   p.tie, p.freeze, p.value, p.min, p.max, p.units, p.fun);
}

%}}}

private define massage_param_structs (pars) %{{{
{
   variable len = 0;

   foreach (pars)
     {
        variable p = ();

        len = max([len, strlen(p.name)]);

        if (p.tie != NULL)
          {
             p.tie = _get_index(p.tie);
          }
        else p.tie = 0;

        if (p.fun != NULL)
          {
             p.fun = sprintf ("\n#=>  %s", strtrim(p.fun));
          }
        else p.fun = "";

        if (andelse
            {p.min == -_isis->DBL_MAX}
            {p.max ==  _isis->DBL_MAX})
          {
             p.min = 0;
             p.max = 0;
          }
     }

   return len;
}

%}}}

private define is_free (p) %{{{
{
   return (p.tie == NULL and p.freeze == 0 and p.fun == NULL);
}

%}}}

private define back_fun_string (i) %{{{
{
   variable s = _isis->_get_instrumental_background_hook_name(i);
   if (s == NULL) return "";
   return sprintf ("# [bgd %d] = %s", i, strtrim(s));
}

%}}}

private define list_par_string (free_only) %{{{
{
   variable pars = get_params();
   if (pars == NULL)
     return NULL;

   if (free_only)
     {
        variable free_pars = array_map (Int_Type, &is_free, pars);
        pars = pars[where(free_pars)];
     }
   variable len;
   len = massage_param_structs (pars);

   variable s = get_fit_fun();
   s = (s == NULL) ? String_Type[0] : [strtrim(s)];

   variable bs, datasets = all_data();
   if (datasets != NULL)
     {
        bs = array_map (String_Type, &back_fun_string, datasets);
        bs = bs[where(bs != "")];
        s = [s, bs];
     }

   variable colheads;
   colheads = sprintf (Fmt_Colheads, len, "  param");
   variable ps;
   ps = array_map (String_Type, &param_string, pars, len);

   return strjoin ([s, colheads, ps], "\n");
}

%}}}

define list_par () %{{{
{
   variable dest = NULL;
   if (_NARGS) dest = ();
   _isis->put_string (dest, _function_name, list_par_string (0));
}

%}}}

define list_free () %{{{
{
   variable dest = NULL;
   if (_NARGS) dest = ();
   _isis->put_string (dest, _function_name, list_par_string (1));
}

%}}}

private variable __isis_save_par_hook = NULL;
define save_par () %{{{
{
   if (_NARGS == 0)
     {
        usage ("%s (filename)", _function_name);
     }
   variable s, file = ();
   variable fp = fopen (file, "w");
   if (fp == NULL)
     throw IsisError, "Failed opening $file for writing"$;
   _isis->put_string (fp, _function_name, list_par_string (0));
   if (2 == is_defined ("isis_save_par_hook"))
     {
        eval ("isis_save_par_hook (\"$file\")"$);
        s = ();
        () = fputs (s, fp);
     }
   else if (__isis_save_par_hook != NULL)
     {
        s = (@__isis_save_par_hook)();
        () = fputs (s, fp);
     }
   if (-1 == fclose (fp))
     {
        variable err = "Failed closing $file"$;
        if (errno) err += sprintf (" (%s)", errno_string(errno));
        throw IsisError, err;
     }
}

%}}}

% randomize fit variables within specified min/max ranges

define randomize () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable list, val;

   if (_NARGS == 0)
     list = [1:get_num_pars ()];
   else
     list = _isis->pop_list (_NARGS, "randomize ([params])");

   foreach (list)
     {
        variable idx = ();
        variable s = get_par_info(idx);

        if (s.freeze) continue;

	if ((s.min == -_isis->DBL_MAX) or (s.max == _isis->DBL_MAX))
	  val = 2*(urand () - 0.5) * _isis->DBL_MAX;
	else
	  val = s.min + urand () * (s.max - s.min);

	set_par (idx, val);
     }
}

%}}}

%}}}

%{{{ set fit function

private define do_fit_fun (num, interactive, msg)
{
   variable def;

   if (_isis->chk_num_args (num, 1, msg))
     return;

   def = ();
   if (def == NULL) def = "";

   _isis->_set_fit_fun (interactive, def);
}

define ifit_fun()
{
   _isis->error_if_fit_in_progress (_function_name);
   do_fit_fun(_NARGS, 1, "ifit_fun (\"function\")");
}

define fit_fun()
{
   _isis->error_if_fit_in_progress (_function_name);
   do_fit_fun(_NARGS, 0, "fit_fun (\"function\")");
}

%}}}

%{{{ fit model to data

private define _do_eval_fit (_do_eval, data_type, msg, nargs) %{{{
{
   variable tmp, ref = &tmp,
     response_type = qualifier ("response", Assigned_ARFRMF);

   switch (nargs)
     {
      case 0:
        % do nothing
     }
     {
      case 1:
        variable x = ();
	if (typeof (x) == Ref_Type)
	  {
             ref = x;
	  }
	else
	  {
	     response_type = x;
	  }
     }
     {
      case 2:
	(response_type, ref) = ();
     }
     {
	% default:
	_pop_n (nargs);
	usage (msg);
     }

   _isis->_set_fit_type (response_type, data_type);

   variable saved_fit_verbose = Fit_Verbose,
     fit_verbose = qualifier ("fit_verbose", NULL);

   variable status;
   try
     {
        if (fit_verbose != NULL)
          {
             Fit_Verbose = fit_verbose;
          }
        status = (@_do_eval)(ref);
     }
   finally
     {
        Fit_Verbose = saved_fit_verbose;
     }
   return status;
}

%}}}

define fit_counts ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = fit_counts ([response_type][, &info_struct])";
   return _do_eval_fit (&_isis->_fit, 0, msg, _NARGS ;; __qualifiers);
}

define fit_flux ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = fit_flux ([response_type][, &info_struct])";
   return _do_eval_fit (&_isis->_fit, 1, msg, _NARGS ;; __qualifiers);
}

define eval_counts ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = eval_counts ([response_type][, &info_struct])";
   return _do_eval_fit (&_isis->_eval_model, 0, msg, _NARGS ;; __qualifiers);
}

define eval_flux ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = eval_flux ([response_type][, &info_struct])";
   return _do_eval_fit (&_isis->_eval_model, 1, msg, _NARGS ;; __qualifiers);
}

private define eval_statistic (data_type)
{
   _isis->_set_fit_type (Assigned_ARFRMF, data_type);
   variable s, status;
   status = _isis->_eval_statistic_only(&s);
   if (status == 0)
     return s;
   else return NULL;
}

define eval_stat_counts ()
{
   _isis->error_if_fit_in_progress (_function_name);
   return eval_statistic (0);
}

define eval_stat_flux ()
{
   _isis->error_if_fit_in_progress (_function_name);
   return eval_statistic (1);
}

%}}}

% differential (unbinned) function

define fitfun_handle () %{{{
{
   variable msg = "handle = fitfun_handle (\"name\")";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;
   return _isis->get_fitfun_handle_intrin ();
}

%}}}

private define resolve_fitfun_handle (handle) %{{{
{
   switch (typeof(handle))
     {
      case Fitfun_Type:
        return handle;
     }
     {
      case String_Type:
        return _isis->get_fitfun_handle_intrin (handle);
     }
     {
      case Ref_Type:
        foreach (_isis->_function_list())
          {
             variable name = ();
             if (handle == __get_reference(name))
               {
                  return  _isis->get_fitfun_handle_intrin (name);
               }
          }
     }

   return NULL;
}

%}}}

define get_cfun2 () %{{{
{
   variable msg = "y = get_fun2 (fitfun, x, [, pars])";
   variable handle, x;
   variable pars=NULL, args=NULL;

   variable nargs = _NARGS;

   if (nargs < 2)
     {
        usage(msg);
        return;
     }

   if (nargs >= 3)
     {
        args = __pop_args (nargs-3);
        (handle, x, pars) = ();
     }
   else
     {
        (handle, x) = ();
     }

   handle = resolve_fitfun_handle (handle);

   if (pars == NULL)
     pars = Double_Type[0];

   if (args == NULL)
     {
        return _isis->eval_diff_fitfun_using_handle_intrin (pars, x, handle);
     }
   else
     {
        return _isis->eval_diff_fitfun_using_handle_intrin (__push_args(args), pars, x, handle);
     }
}

%}}}

define get_cfun () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "y = get_cfun (x)";

   if (_NARGS == 0)
     {
        usage(msg);
        return;
     }
   else if (_NARGS > 1)
     {
        _stk_roll(-_NARGS);
     }

   return _isis->_get_differential_model ();
}

%}}}

define eval_fun2 () %{{{
{
   variable msg = "y = eval_fun2 (fitfun, binlo, binhi [, pars [,args..]])";
   variable handle, lo, hi;
   variable pars=NULL, args=NULL;

   variable nargs = _NARGS;

   if (nargs < 3)
     {
        usage(msg);
        return;
     }

   if (nargs >= 4)
     {
        args = __pop_args (nargs-4);
        (handle, lo, hi, pars) = ();
     }
   else
     {
        (handle, lo, hi) = ();
     }

   handle = resolve_fitfun_handle (handle);

   if (pars == NULL)
     pars = Double_Type[0];

   if (args == NULL)
     {
        return _isis->eval_fitfun_using_handle_intrin (pars, lo, hi, handle);
     }
   else
     {
        return _isis->eval_fitfun_using_handle_intrin (__push_args(args), pars, lo, hi, handle);
     }
}

%}}}

define eval_fun () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "y = eval_fun (lo, hi)";
   variable lo, hi;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (lo, hi) = ();

   variable f = _isis->_get_model_on_user_grid (lo, hi);

   if (length(f) == 1)
     return f[0];

   return f;
}

%}}}

define get_back () %{{{
{
   variable msg = "bgd = get_back (hist_index)";
   variable id, bgd, save_active_dataset;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   id = ();

   save_active_dataset = Isis_Active_Dataset;
   try
     {
        Isis_Active_Dataset = id;
        bgd = _isis->_get_instrumental_background (id);
     }
   finally
     {
        Isis_Active_Dataset = save_active_dataset;
     }

   return bgd;
}
define get_iback ()
{
   vmessage ("*** %s will soon be renamed to get_back", _function_name);
   variable args = __pop_args (_NARGS);
   get_back (__push_args(args));
}
%}}}

%{{{ function caching and cloning

% Model cacheing, and model renaming contributed by
% Mike Noble, Mike Nowak, Manfred Hanke

%!%+
%\function{cache_fun}
%\usage{caching_name = cache_fun (name, lo, hi)}
%\qualifiers{
%\qualifier{suffix}{: identifier string for the caching model (to label different instances)}
%\qualifier{mult}{: the model is not additive, but multiplicative}
%}
%!%-
define cache_fun ()%{{{
{
   variable msg = "caching_name = cache_fun (name, lo, hi);\n\t qualifiers = suffix,mult";

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   variable name, lo, hi;
   (name, lo, hi) = ();

   variable info = get_fun_info (name);
   if (info == NULL)
     throw ApplicationError, "*** Error:  nonexistent fit function '$name'"$;

   if (length(info.name) == 0)
     throw ApplicationError,  "*** Error: no support for caching functions without parameters";

   variable caching_name
     = sprintf ("%s%s_cache", name, qualifier("suffix", ""));

   variable s = struct
     {
        handle = fitfun_handle (name),
        bin_lo = lo,
        bin_hi = hi,
        pars = NULL,
        value = NULL
     };

   _isis->Model_Cache[caching_name] = s;

   variable m = "define ${caching_name}_fit(lo,hi,pars){"$
     + "variable s = _isis->Model_Cache[\"${caching_name}\"];"$
     + "if (lo[0] < s.bin_lo[0] || hi[-1] > s.bin_hi[-1])"
     +  sprintf ("   throw DomainError, \"*** Error: %s: no support for extrapolating beyond cached grid [%g, %g]\";",
                 caching_name, s.bin_lo[0], s.bin_hi[-1])
     + "if(s.pars==NULL || any(pars!=s.pars))"
     + "{"
     +   "s.value = eval_fun2 (s.handle, s.bin_lo, s.bin_hi, pars);"
     +   "s.pars = pars;"
     + "}";

   % multiplicative models should be bin-averaged,
   % additive models are bin-integrated.
   if(qualifier_exists("mult"))
     {
        m+="return rebin(lo,hi, s.bin_lo,s.bin_hi, s.value*(s.bin_hi-s.bin_lo))/(hi-lo);}";
     }
   else
     {
        m+="return rebin(lo,hi, s.bin_lo,s.bin_hi, s.value);}";
     }

   variable i, md = "define ${caching_name}_defaults(i){switch(i)"$;
   _for i (0, length(info.name)-1, 1)
     {
        md += sprintf("{case %d: return(%S, %S, %S, %S);}",
                      i, info.value[i], info.freeze[i], info.min[i], info.max[i]);
     }
   md+="}";

   variable asf = "add_slang_function(\"${caching_name}\", ["$
     + strjoin("\"" + info.name + "\"", ", ")
       + "]);";

   variable spdh = "set_param_default_hook(\"${caching_name}\", &${caching_name}_defaults);"$;

   array_map (Void_Type, &eval, [m, md, asf, spdh]);

   return caching_name;
}
%}}}

%!%+
%\function{alias_fun}
%\synopsis{create a copy of a fit-function with a new name}
%\usage{alias_fun (String_Type name, String_Type alias_name);}
%\qualifiers{
%\qualifier{names}{: array of parameter names}
%\qualifier{values}{: array of default values}
%\qualifier{freeze}{: array of flags whether the parameter is frozen by default}
%\qualifier{min}{: array of default minimal parameter values}
%\qualifier{max}{: array of default maximal parameter values}
%\qualifier{params}{: array of { name, value, freeze, min, max } lists}
%}
%\examples
%     alias_fun("gaussian", "FeKa";
%                  names= ["area [ph/s/cm^2]", "E [keV]", "sigma [keV]"],
%                  values=[ 0.01,               6.4,       0.5         ],
%                  freeze=[ 0,                  1,         0           ],
%                  min=   [ 0,                  5.8,       1e-6        ],
%                  max=   [ 1,                  7,         2           ]);
%
%     % If all properties are specified, the same can be written as follows:
%     alias_fun("gaussian", "FeKa";
%                  params=[ { "area [ph/s/cm^2]", 0.01, 0, 0,    1 },
%                           {          "E [keV]", 6.4,  1, 5.8,  7 },
%                           { "     sigma [keV]", 0.5,  0, 1e-6, 2 } ]);
%!%-
define alias_fun ()%{{{
{
   variable msg = "alias_fun (name, alias_name);\n\t qualifiers = names,values,freeze,min,max,params";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable name, alias_name;
   (name, alias_name) = ();

   variable info = get_fun_info (name);
   if (info == NULL)
     throw ApplicationError, "*** Error:  fit function '$name' does not exist"$;

   variable i, num_pars = length(info.name);

   variable params = qualifier ("params", NULL);
   if (params != NULL)
     {
        _for i (0, num_pars-1, 1)
          {
             (info.name[i], info.value[i], info.freeze[i], info.min[i], info.max[i]) =
               (params[i][0], params[i][1], params[i][2], params[i][3], params[i][4]);
          }
     }
   else
     {
        info.name = qualifier ("names", info.name);
        info.value = qualifier ("values", info.value);
        info.freeze = qualifier ("freeze", info.freeze);
        info.min = qualifier ("min", info.min);
        info.max = qualifier ("max", info.max);
     }

   variable m = "define ${alias_name}_fit(lo,hi,par){return eval_fun2(\"$name\",lo,hi,par);}"$;

   variable md = "define ${alias_name}_defaults(i){"$;
   if (num_pars > 0)
     {
        md += "switch(i)";
        _for i (0, num_pars-1, 1)
          {
             md += sprintf("{case %d: return(%S, %S, %S, %S);}",
                      i, info.value[i], info.freeze[i], info.min[i], info.max[i]);
          }
     }
   md += "}";

   variable asf = "add_slang_function(\"$alias_name\""$;
   if (num_pars > 0)
     {
        asf += ",["$
          + strjoin("\"" + info.name + "\"", ",")
            + "]";
     }
   asf += ")";

   variable spdh = "set_param_default_hook(\"$alias_name\",&${alias_name}_defaults);"$;

   array_map (Void_Type, &eval, [m, md, asf, spdh]);
}
%}}}

%}}}

%{{{ confidence limits

private variable Default_Chisqr_Tolerance = 1.e-3;
private variable Delta_Chisqr =
[
   1.00, % 1-sigma = 68% confidence
   2.71, % 2-sigma = 90% confidence
   6.63  % 3-sigma = 99% confidence
];

private define _conf (num, msg, verbose) %{{{
{
   variable idx, lev, tolerance;

   lev = 1;
   tolerance = Default_Chisqr_Tolerance;

   if (_isis->get_varargs (&idx, &lev, &tolerance, num, 1, msg))
     return;

   _isis->_set_fit_type (qualifier ("response", Assigned_ARFRMF),
                         qualifier_exists ("flux"));

   idx = _get_index (idx);
   return _isis->_confidlev (idx, Delta_Chisqr[lev], verbose, tolerance);
}

%}}}

private define _fconf (num, msg, verbose) %{{{
{
   variable idx, flev, tolerance;

   flev = Delta_Chisqr[1];  % 90% confidence
   tolerance = Default_Chisqr_Tolerance;

   if (_isis->get_varargs (&idx, &flev, &tolerance, num, 1, msg))
     return;

   _isis->_set_fit_type (qualifier ("response", Assigned_ARFRMF),
                         qualifier_exists ("flux"));

   idx = _get_index (idx);
   return _isis->_confidlev (idx, flev, verbose, tolerance);
}

%}}}

define conf ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "(low, high) = conf (param_index [, level [, tolerance]])";
   return _conf (_NARGS, msg, 0 ;; __qualifiers);         % quiet
}

define vconf ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "(low, high) = vconf (param_index [, level [, tolerance]])";
   return _conf (_NARGS, msg, 1 ;; __qualifiers);         % verbose
}

define fconf ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "(low, high) = fconf (param_index [, dchisqr [, tolerance]])";
   return _fconf (_NARGS, msg, 0 ;; __qualifiers);         % quiet
}

define vfconf ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "(low, high) = vfconf (param_index [, dchisqr [, tolerance]])";
   return _fconf (_NARGS, msg, 1 ;; __qualifiers);         % verbose
}

%}}}

%{{{ renorm

private define _renorm (fit_ref, response_type) %{{{
{
   variable i, p, pars, ids, num_pars;
   variable vnorms, num_vnorms, i_scale, norm_scale;

   % Save initial parameter configuration
   pars = get_params();
   if (pars == NULL)
     return -1;

   ERROR_BLOCK
     {
	set_params (pars);
     }

   num_pars = get_num_pars();
   ids = [1:num_pars];

   % Make a list of variable normalization parameters
   % and pick a non-zero one to use to scale the rest.

   vnorms = Int_Type[num_pars];
   num_vnorms = 0;
   i_scale = 0;

   foreach (ids)
     {
	i = ();
	p = pars[i-1];
	if ((p.is_a_norm == 1) and (p.freeze == 0) and (p.tie == NULL))
	  {
	     vnorms[num_vnorms] = i;
	     num_vnorms++;
	     if (p.value != 0.0 and i_scale == 0)
	       i_scale = i;
	  }
     }

   if (num_vnorms == 0)
     {
	if (Isis_Verbose >= _isis->_WARN)
	  message ("*** renorm failed:  no variable norm parameters");
	return -1;
     }

   vnorms = vnorms[[0:num_vnorms-1]];

   if (i_scale == 0)
     {
	% if *all* the variable norms were zero,
	% we don't have anything to go on.....
        throw IsisError, "*** renorm failed:  all variable norms are zero";
     }

   % Guaranteed to be non-zero
   norm_scale = get_par (i_scale);

   if (num_vnorms > 1)
     {
	foreach (vnorms[where(vnorms != i_scale)])
	  {
	     i = ();
	     set_par_fun (i,
			  sprintf ("%e * (_par(%d)/%e)",
				   get_par(i), i_scale, norm_scale));
	  }
     }

   freeze (ids);
   thaw (i_scale);

   % At this point, there should be exactly one free parameter

   variable status, new_pars;
   status = (@fit_ref) (response_type ;; __qualifiers);
   if (-1 == status)
     {
	if (Isis_Verbose >= _isis->_WARN)
	  message ("*** renorm failed");
	set_params (pars);
	return -1;
     }

   new_pars = get_params ();
   set_params (pars);

   foreach (vnorms)
     {
	i = ();
	p = new_pars[i-1];
	set_par (i, p.value);
     }

   return status;
}
%}}}

define renorm_counts () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = renorm_counts ([response_type])";
   variable response_type = Assigned_ARFRMF;

   if (_isis->get_varargs (&response_type, _NARGS, 0, msg))
     return;

   _renorm (&fit_counts, response_type ;; __qualifiers);
}

%}}}

define renorm_flux () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = renorm_flux ([response_type])";
   variable response_type = Assigned_ARFRMF;   % won't apply ARF twice

   if (_isis->get_varargs (&response_type, _NARGS, 0, msg))
     return;

   _renorm (&fit_flux, response_type ;; __qualifiers);
}

%}}}

%}}}

%{{{ aux fit kernels

define load_kernel () %{{{
{
   variable msg = "ret = load_kernel (file, name[, args])";
   variable file, name, args = NULL;

   if (_isis->get_varargs (&file, &name, &args, _NARGS, 2, msg))
     return;

   variable path = find_library_name (file);
   if (path == NULL)
     {
	vmessage ("File not found:  %s", file);
	return;
     }

   return _isis->_load_kernel (args, path, name);
}

%}}}

define set_kernel () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_kernel (data_list, kernel_name)";
   variable i, hist_index, kernel_name;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (hist_index, kernel_name) = ();

   if (kernel_name == NULL)
     kernel_name = "std";

   foreach (hist_index)
     {
	i = ();
	_isis->_set_kernel (kernel_name, i);
     }
}

%}}}

define print_kernel () %{{{
{
   variable msg = "print_kernel (hist_index[])";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable i, list = ();

   foreach (list)
     {
	i = ();
	_isis->_print_kernel (i);
     }
}

%}}}

define list_kernels () %{{{
{
   variable msg = "list_kernels ()";

   if (_isis->chk_num_args (_NARGS, 0, msg))
     return;

   variable list = [_isis->_list_kernels ()];

   variable i = array_sort (list, &strcmp);
   writecol (stdout, list[i]);
}

%}}}

%}}}

%{{{ alternate minimizer and fit-statistic

define load_fit_method ()
{
   variable msg = "status = load_fit_method (file, name)";
   variable file, name;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return -1;

   (file, name) = ();

   variable path = find_library_name (file);
   if (path == NULL)
     {
	vmessage ("File not found:  %s", file);
	return;
     }

   return _isis->_load_fit_method (path, name);
}

define load_fit_statistic ()
{
   variable msg = "load_fit_statistic (file, name)";
   variable file, name;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (file, name) = ();

   variable path = find_library_name (file);
   if (path == NULL)
     {
	vmessage ("File not found:  %s", file);
	return;
     }

   return _isis->_load_fit_statistic (path, name);
}

define add_slang_statistic ()
{
   variable msg = "add_slang_statistic (name, &fun, &report)";
   variable name, fun, report;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (name, fun, report) = ();

   _isis->_add_slang_statistic (fun, report, name);
}

define list_fit_methods ()
{
   _isis->list_statistics_and_engines();
}

define get_fit_method ()
{
   _isis->_get_fit_method_name ();
}

define get_fit_statistic ()
{
   _isis->_get_fit_statistic_name ();
}

define set_fit_method ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_fit_method (name)";
   variable name;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   name = ();

   _isis->_set_fit_method_name (name);
}

define set_fit_statistic ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_fit_statistic (name)";
   variable name;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   name = ();

   _isis->_set_fit_statistic_name (name);
}

define set_fit_constraint ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "Usage:  set_fit_constraint (&function, [, param_names[] ])";
   variable ref, names = NULL;

   if (_NARGS == 0)
     {
        vmessage (msg);
        return;
     }

   if (_NARGS > 1)
     {
        names = __pop_args (_NARGS-1);
        names = [__push_args(names)];
     }

   ref = ();

   if (typeof(ref) == String_Type)
     ref = __get_reference(ref);

   if (ref != NULL)
     {
        variable constraint_name = "constraint";
        eval(sprintf("public define %s_fit(l,h,p){return 0*l;}",
                     constraint_name));
        if (andelse
            {_NARGS > 1}
            {names[0] != NULL})
          add_slang_function(constraint_name, names);
        else
          add_slang_function(constraint_name);
     }

   _isis->_set_fit_constraint_fun (ref);
}

define set_fit_range_hook ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_fit_range_hook (&func)";
   variable ref;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   ref = ();

   _isis->_set_fit_range_hook (ref);
}

define set_post_model_hook ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_post_model_hook (id[], &func)";
   variable ids, hook;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (ids, hook) = ();

   if (typeof(hook) != Array_Type)
     hook = [hook];

   array_map (Void_Type, &_isis->_set_post_model_hook, ids, hook);
}

define set_pre_combine_hook () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_pre_combine_hook (hist_index[], &hook)";
   variable id, hook;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (id, hook) = ();

   if (typeof(hook) == String_Type)
     hook = __get_reference (hook);

   array_map (Void_Type, &_isis->_set_pre_combine_hook, id, hook);
}

%}}}

define set_define_model_hook ()
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_define_model_hook (&func)";
   variable ref;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   ref = ();

   _isis->set_define_model_hook_intrin (ref);
}

%}}}

%{{{ Some utility fit-functions

define delta_fit (lo, hi, par)
{
   variable result = Double_Type [length(lo)];
   result[ where(lo <= par[1] and par[1] < hi) ] = par[0];
   return result;
}

define null_fit (lo, hi, par)
{
   return Double_Type[length(lo)];
}

define bin_width_fit (lo, hi, par)
{
   return (hi - lo);
}

define bin_width_en_fit (lo, hi, par)
{
   variable hc = _A(1);
   return hc/lo - hc/hi;
}

define bin_center_fit (lo, hi, par)
{
   return 0.5 * (hi + lo);
}

define bin_center_en_fit (lo, hi, par)
{
   variable hc = _A(1);
   return 0.5 * (hc/hi + hc/lo);
}

define init_static_slang_functions () %{{{
{
   add_slang_function ("delta", ["norm", "lambda"]);
   add_slang_function ("null", [""]);
   add_slang_function ("bin_width", [""]);
   add_slang_function ("bin_width_en");
   add_slang_function ("bin_center", [""]);
   add_slang_function ("bin_center_en", [""]);
}

%}}}

init_static_slang_functions();

%}}}

%{{{ Generate and plot 2D confidence contours

private variable Contour_Type = struct {index, min, max, num};
private variable Param_Type = struct {value, conf_min, conf_max};

define conf_grid () %{{{
{
   variable msg = "Struct_Type = conf_grid (index, min, max, num)";
   variable p = @Contour_Type;

   if (_isis->chk_num_args (_NARGS, 4, msg))
     return;

   (p.index, p.min, p.max, p.num) = ();

   if (_get_index (p.index) < 0)
     return NULL;

   return p;
}

%}}}

private define generate_grid (pmin, pmax, num) %{{{
{
   return pmin + ((pmax - pmin)/double(num)) * [0.0:num:1.0];
}

%}}}

define get_params () %{{{
{
   variable msg = "p = get_params ([list])";
   variable n = get_num_pars ();
   if (n == 0)
     {
        _pop_n (_NARGS);
        return NULL;
     }

   variable k, pars;

   if (_NARGS == 0)
     pars = [1:n];
   else
     {
	pars = _isis->pop_list (_NARGS, msg);
	if (pars == NULL)
	  return;

        % User might mistake this for get_fun_params().
        % Maybe both functions should use a single interface?
        if (_NARGS == 1 and typeof(pars) == String_Type)
          usage(msg);

	n = length(pars);
     }

   variable i, par_info = Struct_Type [n];

   _for (0, n-1, 1)
     {
	i = ();
	k = pars[i];
	par_info[i] = get_par_info (k);
     }

   return par_info;
}

%}}}

define set_params () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_params (pars)";
   variable par_info, p;

   if (_isis->chk_num_args(_NARGS, 1, msg))
     return;

   par_info = ();

   foreach (par_info)
     {
	p = ();
	set_par_fun (p.index, p.fun);
	set_par (p.index, p.value, p.freeze, p.min, p.max, p.step);
     }
}

%}}}

define num_free_params () %{{{
{
   variable p = get_params();
   if (p == NULL)
     return -1;

   variable x, n = 0;
   foreach (p)
     {
	x = ();
	if (x.freeze == 0)
	  n++;
     }

   return n;
}

%}}}

private define map_chisqr (ip1, p1, ip2, p2, info) %{{{
{
   variable save_fit_verbose = Fit_Verbose;
   ERROR_BLOCK
     {
	Fit_Verbose = save_fit_verbose;
     }
   Fit_Verbose = -1;

   variable fit_ref, eval_ref, fail_ref, save_ref, mask_ref;

   fit_ref = info.fit_ref;
   eval_ref = info.eval_ref;
   fail_ref = info.fail_ref;
   save_ref = info.save_ref;
   mask_ref = info.mask_ref;

   variable k, num_free, pars;
   variable p, best_pars, try_pars;

   pars = [1:get_num_pars ()];
   best_pars = get_params ();
   num_free = num_free_params();

   variable i1, num_p1 = length(p1);
   variable i2, num_p2 = length(p2);

   variable chisqr = Float_Type [num_p2, num_p1];
   chisqr[*,*] = -1.0;
   variable fit_info;

   variable print_status;
   print_status = (save_fit_verbose >= 0) or (Isis_Verbose >= _isis->_WARN);

   if (print_status)
     () = fprintf (stderr, "Mapping %s:\n", Fit_Statistic);

   _for (0, num_p1 - 1, 1)
     {
	i1 = ();

	_for (0, num_p2 - 1, 1)
	  {
	     i2 = ();

	     foreach (pars)
	       {
		  k = ();
		  p = best_pars[k-1];
		  set_par (k, p.value, p.freeze, p.min, p.max);
	       }

	     variable pa = p1[i1]*[1.0, 1.0];
	     variable pb = p2[i2]*[1.0, 1.0];

	     set_par (ip1, p1[i1], 1, min(pa), max(pa));
	     set_par (ip2, p2[i2], 1, min(pb), max(pb));

             variable masked_out =
               ((mask_ref != NULL)
                && (0 == (@mask_ref) (p1[i1], p2[i2] ;; __qualifiers)));

	     ifnot (masked_out)
	       {
		  if (num_free <= 2)
		    () = (@eval_ref) (&fit_info ;; __qualifiers);
		  else
		    {
		       try_pars = get_params();
                       if (-1 == (@fit_ref) (&fit_info ;; __qualifiers))
                         {
                            if (fail_ref != NULL)
                              (@fail_ref) (ip1, ip2, @best_pars, try_pars, fit_info ;; __qualifiers);
                         }
		    }
                  chisqr[i2, i1] = fit_info.statistic;
	       }
             else fit_info = NULL;

	     if (save_ref != NULL)
	       (@save_ref) (fit_info ;; __qualifiers);

             if (print_status)
               {
                  () = fprintf (stderr, "Completed %d/%d\r",
                                1 + i2 + i1 * num_p2,
                                num_p1*num_p2);
               }
	  }
     }

   if (print_status)
     () = fputs ("\n", stderr);

   Fit_Verbose = save_fit_verbose;

   set_params (best_pars);
   () = (@eval_ref) ( ;; __qualifiers);

   return chisqr;
}

%}}}

require ("fork_socket");

#ifexists fork_slave

private variable Parallel_Map_Info;

private define send_next_task (s) %{{{
{
   if (write_array (s.fp, Parallel_Map_Info.next_task))
     throw IOError, "*** master:  write failed";
   Parallel_Map_Info.next_task++;
}

%}}}

private define subarray_indices (task) %{{{
{
   variable q = Parallel_Map_Info;

   variable
     ix = task / q.num_ysub,
     iy = task mod q.num_ysub;

   return q.xsub[ix], q.ysub[iy];
}

%}}}

private define slave_has_result (s) %{{{
{
   variable task = read_array (s.fp, 1, Int_Type)[0];

   variable xsub, ysub;
   (xsub, ysub) = subarray_indices (task);

   variable
     nx = length(xsub),
     ny = length(ysub);

   variable subarray = read_array (s.fp, nx * ny, Float_Type);
   Parallel_Map_Info.map[ysub, xsub] = _reshape (subarray, [ny, nx]);
}

%}}}

private define message_handler (s, msg) %{{{
{
   switch (msg.type)
     {
      case SLAVE_RESULT:
        slave_has_result (s);
        send_next_task (s);
     }
     {
        print(msg);
        throw ApplicationError, "Unsupported message type";
     }
}

%}}}

private define slave_process (s, task, ip1, pxs, ip2, pys, info) %{{{
{
   while (task < Parallel_Map_Info.num_tasks)
     {
        variable xsub, ysub;
        (xsub, ysub) = subarray_indices (task);

        variable e;
        try (e)
          {
             variable map = map_chisqr (ip1, pxs[xsub], ip2, pys[ysub], info);
          }
        catch AnyError:
          {
             pid_vmessage ("Caught exception calling map_chisqr!!!");
             print(e);
             _exit(1);
          }

        send_msg (s, SLAVE_RESULT);
        if (0 != write_array (s.fp, task)
            or 0 != write_array (s.fp, map))
          throw IOError, "*** slave: write failed";

        task = read_array (s.fp, 1, Int_Type)[0];
     }

   return 0;
}

%}}}

private define partition_indices (num, num_sub) %{{{
{
   variable
     sub = Array_Type[num_sub],
     block = num / num_sub,
     num_left = num - block * num_sub;

   variable i, mx, mn = 0;

   _for i (0, num_sub-1, 1)
     {
        mx = mn + block;
        if (num_left > 0)
          {
             mx += 1;
             num_left--;
          }

        sub[i] = [mn:mx-1];
        mn = mx;
     }

   return sub;
}

%}}}

private define parallel_map_chisqr (num_slaves, px, py, %{{{
                                    ix, pxs, iy, pys, info)
{
   variable num_xsub, num_ysub;
   num_xsub = qualifier ("num_xsub", 1);
   num_ysub = qualifier ("num_ysub", num_slaves);

   Parallel_Map_Info = struct
     {
        num_ysub = num_ysub,
        num_xsub = num_xsub,
        num_tasks = num_xsub * num_ysub,
        next_task = 0,
        xsub = partition_indices (px.num, num_xsub),
        ysub = partition_indices (py.num, num_ysub),
        map = Float_Type[py.num, px.num],
     };

   variable slaves = new_slave_list (;; __qualifiers);

   variable i, s;
   _for i (0, num_slaves-1, 1)
     {
        variable task = Parallel_Map_Info.next_task;
        s = fork_slave (&slave_process, task, ix, pxs, iy, pys, info ;; __qualifiers);
        append_slave (slaves, s);
        Parallel_Map_Info.next_task++;
     }

   manage_slaves (slaves, &message_handler ;; __qualifiers);

   return Parallel_Map_Info.map;
}

%}}}

#endif

private define generate_contours (px, py, info) %{{{
{
   variable
     pxs = generate_grid (px.min, px.max, px.num),
     pys = generate_grid (py.min, py.max, py.num);

   variable
     ix = _get_index(px.index),
     iy = _get_index(py.index);

#ifexists fork_slave
   variable serial = qualifier_exists ("serial");
   variable num_slaves = qualifier ("num_slaves", _num_cpus());
   if ((serial == 0) && (num_slaves > 1))
     {
        return parallel_map_chisqr (num_slaves, px, py,
                                    ix, pxs, iy, pys, info ;; __qualifiers);
     }
#endif

   return map_chisqr (ix, pxs, iy, pys, info);

}

%}}}

private define _conf_map (px, py, info) %{{{
{
   variable save_fit_verbose = Fit_Verbose;
   variable save_isis_verbose = Isis_Verbose;
   _isis->_set_conf_limit_search (1);

   ERROR_BLOCK
     {
	_isis->_set_conf_limit_search (0);
	Fit_Verbose = save_fit_verbose;
        Isis_Verbose = save_isis_verbose;
     }
   if (Fit_Verbose == 0)
     Fit_Verbose = -1;
   if (Isis_Verbose == 0)
     Isis_Verbose = -1;

   variable s = struct {chisqr, px, py, best, px_best, py_best};
   variable fit_info, best_pars, eval_ref;

   eval_ref = info.eval_ref;

   () = (@eval_ref) (&fit_info;; __qualifiers);
   best_pars = get_params ();

   s.px_best = get_par (px.index);
   s.py_best = get_par (py.index);
   s.best = fit_info.statistic;
   s.px = @px;
   s.py = @py;

   s.chisqr = generate_contours (px, py, info ;; __qualifiers);
   s.chisqr -= s.best;

   set_params (best_pars);
   () = (@eval_ref) ( ;; __qualifiers);

   Isis_Verbose = save_isis_verbose;
   Fit_Verbose = save_fit_verbose;
   _isis->_set_conf_limit_search (0);

   return s;
}

%}}}

private define merge_user_info_struct (info, user_info) %{{{
{
   if (user_info == NULL)
     {
        info.mask_ref = NULL;
        info.fail_ref = NULL;
        info.save_ref = NULL;
        return;
     }

   if (struct_field_exists (user_info, "mask"))
     info.mask_ref = user_info.mask;
   else
     info.mask_ref = NULL;

   if (struct_field_exists (user_info, "fail"))
     info.fail_ref = user_info.fail;
   else
     info.fail_ref = NULL;

   if (struct_field_exists (user_info, "save"))
     info.save_ref = user_info.save;
   else
     info.save_ref = NULL;
}

%}}}

define conf_map_counts () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "Struct_Type = conf_map_counts (Struct_Type, Struct_Type [, info])";
   variable px, py, user_info;

   user_info = NULL;

   if (_isis->get_varargs (&px, &py, &user_info, _NARGS, 2, msg))
     return;

   variable info = struct
     {
	fit_ref, eval_ref, fail_ref, save_ref, mask_ref
     };

   info.fit_ref = &fit_counts;
   info.eval_ref = &eval_counts;
   merge_user_info_struct (info, user_info);

   return _conf_map (px, py, info ;; __qualifiers);
}

%}}}

define conf_map_flux () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "Struct_Type = conf_map_flux (Struct_Type, Struct_Type [, info])";
   variable px, py, user_info;

   user_info = NULL;

   if (_isis->get_varargs (&px, &py, &user_info, _NARGS, 2, msg))
     return;

   variable info = struct
     {
	fit_ref, eval_ref, fail_ref, save_ref, mask_ref
     };

   info.fit_ref = &fit_flux;
   info.eval_ref = &eval_flux;
   merge_user_info_struct (info, user_info);

   return _conf_map (px, py, info ;; __qualifiers);
}

%}}}

private define x_interpolate_joint (chisqr, xs, c, p) %{{{
{
   variable clo, chi, xlo, xhi;;

   xlo = xs[p[0] - 1];
   xhi = xs[p[0]    ];

   clo = chisqr[p[1], p[0] - 1];
   chi = chisqr[p[1], p[0]    ];

   return xlo + (xhi - xlo) * (c - clo) / (chi - clo);
}

%}}}

private define y_interpolate_joint (chisqr, ys, c, p) %{{{
{
   variable clo, chi, ylo, yhi;;

   ylo = ys[p[1] - 1];
   yhi = ys[p[1]    ];

   clo = chisqr[p[1] - 1, p[0]];
   chi = chisqr[p[1]    , p[0]];

   return ylo + (yhi - ylo) * (c - clo) / (chi - clo);
}

%}}}

define conf_joint () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "Struct_Type = conf_joint (Struct_Type, [, delta_chisqr])";
   variable s, px, py;

   % default to 1-sigma (68.3% conf.) limit for 2 d.o.f
   variable c = 2.30;

   if (_isis->get_varargs (&s, &c, _NARGS, 1, msg))
     return;

   px = s.px;
   py = s.py;

   variable dims;
   (dims,,) = array_info (s.chisqr);

   variable xs, ys;
   xs = generate_grid (px.min, px.max, px.num);
   ys = generate_grid (py.min, py.max, py.num);

   variable pxmin, pxmax, pymin, pymax;
   pxmin = [dims[1]-1, 0];
   pxmax = [0, 0];
   pymin = [0, dims[0]-1];
   pymax = [0, 0];

   variable i, j;

   _for (1, dims[0]-1, 1)
     {
        j = ();
        _for (1, dims[1]-1, 1)
          {
             i = ();
	     if (i <= pxmin[0])
	       {
		  if (s.chisqr[j, i-1] > c and c >= s.chisqr[j,i])
		    pxmin = [i,j];
	       }
	     else if (i >= pxmax[0])
	       {
		  if (s.chisqr[j, i-1] < c and c <= s.chisqr[j,i])
		    pxmax = [i,j];
	       }

	     if (j <= pymin[1])
	       {
		  if (s.chisqr[j-1, i] > c and c >= s.chisqr[j,i])
		    pymin = [i,j];
	       }
	     else if (j >= pymax[1])
	       {
		  if (s.chisqr[j-1, i] < c and c <= s.chisqr[j,i])
		    pymax = [i,j];
	       }
	  }
     }

   variable jxy = struct
     {
	xmin, xmax, ymin, ymax
     };

   jxy.xmin = x_interpolate_joint (s.chisqr, xs, c, pxmin);
   jxy.xmax = x_interpolate_joint (s.chisqr, xs, c, pxmax);
   jxy.ymin = y_interpolate_joint (s.chisqr, ys, c, pymin);
   jxy.ymax = y_interpolate_joint (s.chisqr, ys, c, pymax);

   return jxy;
}

%}}}

define array_fit () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "(pars, stat) = array_fit (x, y, wt, pars, par_min, par_max, fun_ref)";
   variable x, y, wt, pars, par_min, par_max, fun;

   !if (_NARGS == 7)
     usage (msg);

   (x, y, wt, pars, par_min, par_max, fun) = ();

   if (x == NULL or y == NULL or pars == NULL or fun == NULL)
     usage (msg);

   if (typeof(pars) != Array_Type)
     pars = [pars];

   if (length(x) != length(y)
       or (wt != NULL and length(x) != length(wt)))
     usage (msg);

   if (wt == NULL)
     {
	wt = @x;
	wt[*] = 1.0;
     }
   else if (length(where(wt <= 0.0)) > 0)
     {
	verror ("*** array_fit:  weights must be positive");
     }

   if (par_min == NULL)
     {
	par_min = @pars;
	par_min[*] = -_isis->DBL_MAX;
     }
   else if (length(par_min) != length(pars))
     {
	verror ("*** array_fit:  length(par_min) != length(pars)");
     }

   if (par_max == NULL)
     {
	par_max = @pars;
	par_max[*] = _isis->DBL_MAX;
     }
   else if (length(par_max) != length(pars))
     {
	verror ("*** array_fit:  length(par_max) != length(pars)");
     }

   _isis->_array_fit (x, y, wt, @pars, par_min, par_max, fun);
}

%}}}

%}}}

%{{{ Monte-Carlo search
%
%{{{  Purpose:
%       Carry out a Monte-Carlo search of the fit-function
%       parameter space to find a good set of starting parameters.
%
%  Notes:
%       A Monte-Carlo search may be helpful when examining
%       a complex chi-square space such as that associated
%       with the pileup model.
%
%       Each variable fit-parameter should be constrained
%       to lie within some reasonable min/max range.
%       If the allowed range is unconstrained, the Monte-Carlo
%       search will sample parameter values over the entire
%       range of representable double-precision numbers.
%
%  Usage:  best = fit_search (num, &ref [,"dir" [, save_flag]);
%
%  e.g. to carry out 100 random trials with 'eval_counts':
%     best = fit_search (100, &eval_counts);
%}}}

private define linear_random (mn, mx) %{{{
{
   return mn + (mx - mn) * urand();
}

%}}}

private define random_params (p) %{{{
{
   variable i, n = get_num_pars ();

   for (i = 0; i < n; i++)
     {
	if (p[i].freeze == 0 and p[i].tie == NULL)
	  {
	     p[i].value = linear_random (p[i].min, p[i].max);
	  }
     }

   return p;
}

%}}}

private variable Fit_Search_Ctrl = NULL;
private variable _Fit_Search_Info = NULL;
define fit_search_info ()
{
   return _Fit_Search_Info;
}

private define fit_search_save_par_hook () %{{{
{
   variable info = _Fit_Search_Info;
   if (info == NULL)
     return "";

   return sprintf ("#  statistic = %20.15e  num_vary = %d  num_bins = %d",
                   info.statistic, info.num_variable_params, info.num_bins);
}

%}}}

define save_search_pars (dir, prefix, suffix) %{{{
{
   variable file = sprintf ("%s.%d.%S", prefix, getpid(), suffix);
   save_par (path_concat (dir, file));
}

%}}}

private define test_params (p, trial) %{{{
{
   variable
     func_ref = Fit_Search_Ctrl.method,
     dir = Fit_Search_Ctrl.dir,
     save_all = Fit_Search_Ctrl.save_all;

   variable info;

   set_params (p);
   () = (@func_ref)(&info ;; __qualifiers);

   _Fit_Search_Info = @info;

   if (save_all)
     {
        save_search_pars (dir, "all/fit", trial);
     }

   return info.statistic;
}

%}}}

private define search_loop (trials, seed, best) %{{{
{
   variable
     forked = (getpid() != Fit_Search_Ctrl.master_pid),
     verbose = qualifier_exists ("verbose") && not forked;

   variable
     dir = qualifier ("dir", NULL),
     save_all = qualifier_exists ("save_all");

   seed_random (seed);

   variable k, i = 0,
     p = best.p,
     num = length(trials);

   foreach k (trials)
     {
	if (verbose)
	  () = fprintf (stderr, "fit_search:  %3d/%3d\tbest=%7.3f\r",
			i, num, best.stat);

	variable stat = test_params (random_params (p), i ;; __qualifiers);

	if (stat < best.stat)
	  {
	     best.stat = stat;
	     best.p = get_params ();
	     if (dir != NULL)
	       {
                  k++;
                  save_search_pars (dir, "best", k);
	       }
	  }

	i++;
     }

   return best;
}

%}}}

private define fit_search_slave (s, trials, seed, best) %{{{
{
   best = search_loop (trials, seed, best ;; __qualifiers);

   send_msg (s, SLAVE_RESULT);
   send_objs (s, best);

   return 0;
}

%}}}

private define fit_search_handler (s, msg) %{{{
{
   switch (msg.type)
     {
      case SLAVE_RESULT:
        s.data = recv_objs (s);
     }
}

%}}}

private define iterate_func_ref (num, func_ref) %{{{
{
   variable
     dir = qualifier ("dir", sprintf ("fit_search_%d", getpid())),
     save_all = qualifier_exists ("save_all"),
     seed = qualifier ("seed", _time);

   Fit_Search_Ctrl = struct
     {
        method = func_ref,
        dir = dir,
        save_all = save_all,
        master_pid = getpid()
     };

   if (dir != NULL)
     {
	() = mkdir (dir, 0777);
	if (save_all != 0)
	  () = mkdir (dir + "/all", 0777);
     }

   variable i = 0, best = struct {p, stat};

   best.p = get_params ();
   best.stat = test_params (best.p, i ;; __qualifiers);
   best.p = get_params ();

   variable trials = [0:num-1];
   if (dir != NULL)
     {
        save_search_pars (dir, "best", trials[0]);
     }

   variable serial = qualifier_exists ("serial");
   variable num_slaves = qualifier ("num_slaves", _num_cpus());
   variable parallel = (serial == 0) && (num_slaves > 1);

   ifnot (parallel)
     {
        best = search_loop (trials, seed, best ;; __qualifiers);
     }
   else
     {
        variable slave_trials = partition_indices (num, num_slaves);
        variable ii = 0, s, slaves = new_slave_list ();
        loop (num_slaves)
          {
             variable slave_seed = int (urand() * INT_MAX);
             s = fork_slave (&fit_search_slave, slave_trials[ii], slave_seed, best ;; __qualifiers);
             append_slave (slaves, s);
             ii++;
          }

        manage_slaves (slaves, &fit_search_handler);
        foreach s (slaves)
          {
             variable b = s.data;
             if (b == NULL)
               continue;
             if (b.stat < best.stat)
               {
                  best.stat = b.stat;
                  best.p = @b.p;
               }
          }
     }

   if (qualifier_exists ("verbose"))
     vmessage ("\nbest statistic = %g", best.stat);

   set_params (best.p);

   if (func_ref == __get_reference ("fit_counts"))
     () = eval_counts ( ;; __qualifiers);
   else if (func_ref == __get_reference ("fit_flux"))
     () = eval_flux ( ;; __qualifiers);

   return best;
}

%}}}

define fit_search () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable num, func_ref;
   variable msg = "best = fit_search (num, &func [; qualifiers])";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (num, func_ref) = ();

   variable stat, bkp_fit_verbose;
   try
     {
        bkp_fit_verbose = Fit_Verbose;
        if (Fit_Verbose == 0) Fit_Verbose = -1;
        _Fit_Search_Info = NULL;
        __isis_save_par_hook = &fit_search_save_par_hook;
        stat = iterate_func_ref (num, func_ref ;; __qualifiers);
     }
   finally
     {
        _Fit_Search_Info = NULL;
        __isis_save_par_hook = NULL;
        Fit_Verbose = bkp_fit_verbose;
     }

   return stat;
}

%}}}

%}}}

% dataset combinations

define combine_datasets () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "gid = combine_datasets (indices [, weights[] ])";
   variable n, args, indices, weights = NULL;

   if (_NARGS == 0)
     {
        usage(msg);
        return;
     }

   ERROR_BLOCK
     {
	_clear_error();
	return -1;
     }

   n = _NARGS;
   args = __pop_args(n);

   % combine_datasets (1, 2, [0.5, 0.5])
   % combine_datasets ([1,2], [0.5,0.5])
   % combine_datasets (1, 2)
   % combine_datasets ([1,2])

   if (andelse
       {n > 2}
       {length(args[n-1].value) == n-1})
     {
	weights = args[n-1].value;
	indices = [__push_args(args[[0:n-2]])];
     }
   else if (andelse
	    {n == 2}
	    {typeof(args[0].value) == Array_Type}
	    {length(args[1].value) == length(args[0].value)})
     {
	__push_args (args);
	(indices, weights) = ();
     }
   else
     {
	__push_args (args);
	indices = _isis->pop_list (n, msg);
	if (indices == NULL)
	  return -1;
	weights = ones(length(indices));
     }

   return _isis->_define_hist_combination (indices, weights);
}

%}}}

define uncombine_datasets () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "uncombine_datasets (gid[])";
   variable x, gid = _isis->pop_list (_NARGS, msg);
   if (gid == NULL)
     return;

   foreach (gid)
     {
	x = ();
	_isis->_undefine_hist_combination (x);
     }
}

%}}}

define combination_members () %{{{
{
   variable msg= "combination_members (gid)";
   variable gid;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   gid = ();

   variable num_members = 0;

   variable all, n, d;
   all = all_data();
   n = length(all);

   for (d = 0; d < n; d++)
     {
	if (0 == _isis->_hist_is_combined (all[d], gid))
	  continue;
	all[d];
	num_members++;
     }

   if (num_members == 0)
     return NULL;
   else if (num_members > 1)
     {
	variable list = __pop_args (num_members);
	return [__push_args (list)];
     }
}

%}}}

define get_combined (gid, get_func) %{{{
{
   variable msg = "Struct_Type = get_combined (gid, &get_function)";
   variable members, s, h, m, weight, use_weight;

   weight = 1.0;
   use_weight = 0;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   members = combination_members (gid);
   if (members == NULL)
     return NULL;

   m = members[0];

   if (orelse
       {get_func == &get_data_counts}
       {get_func == &get_data_flux})
     {
	use_weight = 1;
	weight = get_data_info (m).combo_weight;
     }

   s = @get_func(m);
   s.value *= weight;

   if (s.err != NULL)
     s.err = (weight * s.err)^2;

   foreach (members[[1:]])
     {
	m = ();
	h = @get_func (m);
	if (use_weight)
	  weight = get_data_info (m).combo_weight;

	s.value += weight * h.value;

	if (s.err != NULL)
	  s.err += (weight * h.err)^2;
     }

   if (s.err != NULL)
     s.err = sqrt(s.err);

   return s;
}

%}}}

% support storing most recent combined data/model
private variable Combined_Datasets;

private define combination_struct (combo_id, indices, data, weights)
{
   variable s = struct
     {
        combo_id = combo_id,
        indices = indices,
        data = data,
        err = sqrt(1./weights),
        model
     };
   return s;
}

private define find_combo (gid) %{{{
{
   ifnot (__is_initialized (&Combined_Datasets))
     return NULL;

   variable x;
   foreach x (Combined_Datasets)
     {
        if (x.combo_id == gid)
          return x;
     }

   return NULL;
}

%}}}

private define store_combined_data (combo_ids, offsets, nbins, data, weights) %{{{
{
   variable uniq_indices = unique(combo_ids),
     n = length(uniq_indices);

   variable j, lst = {};
   _for j (0, n-1, 1)
     {
        variable i = uniq_indices[j];
        variable k = offsets[i] + [0:nbins[i]-1];
        variable s = combination_struct (combo_ids[i], k, data[k], weights[k]);
        list_append (lst, s);
     }

   Combined_Datasets = lst;
}

%}}}

private define store_combined_models (combo_ids, models) %{{{
{
   variable uniq_indices = unique(combo_ids),
     n = length(uniq_indices);

   variable j;
   _for j (0, n-1, 1)
     {
        variable i = uniq_indices[j];
        variable s = find_combo (combo_ids[i]);
        s.model = models[s.indices];
     }
}

%}}}

_isis->set_combined_store_data_ref (&store_combined_data);
_isis->set_combined_store_models_ref (&store_combined_models);

private define get_combo_result (gid) %{{{
{
   ifnot (__is_initialized (&Combined_Datasets))
     return NULL;

   if (gid == NULL)
     {
        return Combined_Datasets;
     }

   variable g, a = Struct_Type[0];
   foreach g (gid)
     {
        variable x = find_combo(g);
        if (x == NULL)
          continue;
        a = [a, x];
     }

   return length(a) ? a : NULL;
}

%}}}

define get_combined2 () %{{{
{
   variable msg = "get_combined2 (gid_list)";
   variable gid;

   if (_NARGS == 0)
     gid = NULL;
   else if (_NARGS == 1)
     gid = ();
   else
     {
        gid = __pop_args (_NARGS);
        gid = [__push_args (gid)];
     }

   get_combo_result (gid);
}
%}}}

define rebin_combined () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "rebin_combined (gid[], min_counts | rebin_mask);";
   variable gids, rebin_info;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (gids, rebin_info) = ();

   foreach (gids)
     {
        variable gid = ();

        variable m = combination_members (gid);

        if (typeof(rebin_info) != Array_Type)
          {
             rebin_data (m, 0);
             match_dataset_grids (m);
             variable sid = define_counts (get_combined (gid, &get_data_counts));
             rebin_data (sid, rebin_info);
             rebin_info = get_data_info (sid).rebin;
             delete_data (sid);
          }

        foreach (m)
          {
             variable id = ();
             rebin_data (id, rebin_info);
          }
     }
}

%}}}

define set_eval_grid_method () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_eval_grid_method (method, datasets[][, &hook[, cache_model_values]]);\n     method = SEPARATE_GRID | MERGED_GRID | USER_GRID";
   variable method, ids, ref = NULL, cache_model_values = 0;

   if (_isis->get_varargs (&method, &ids, &ref, &cache_model_values, _NARGS, 2, msg))
     return;

   _isis->_set_eval_grid_method (ref, ids, method, cache_model_values);
}

%}}}

% Support using optimization methods implemented in S-Lang

define __eval_residuals () %{{{
{
   variable msg = "res[] = __eval_residuals (obj, params[]);";
   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable enable_copying = 1;
   if (qualifier_exists ("nocopy"))
     enable_copying = 0;

   variable obj, params;
   (obj, params) = ();

   return _isis->fobj_eval_residuals (params, obj, enable_copying);
}

%}}}

define __eval_stat () %{{{
{
   variable msg = "s = __eval_stat (obj, params[]);";
   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable enable_copying = 1;
   if (qualifier_exists ("nocopy"))
     enable_copying = 0;

   variable obj, params;
   (obj, params) = ();

   variable err, stat;
   (err, stat,,) = _isis->fobj_eval_statistic (params, obj, enable_copying);

   return stat;
}

%}}}

define __data_weights () %{{{
{
   variable msg = "w[] = __data_weights (obj);";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable obj = ();
   return _isis->fobj_get_data_weights (obj);
}

%}}}

private define fit_object_eval_statistic (s, pars) %{{{
{
   variable enable_copying = 1;
   if (qualifier_exists ("nocopy"))
     enable_copying = 0;

   (s.status, s.statistic, s.num_vary, s.num_points) =
     _isis->fobj_eval_statistic (pars, s.object, enable_copying);
   return s.statistic;
}

%}}}

private define fit_object_close (s) %{{{
{
   s.object = 0;
}

%}}}

define open_fit () %{{{
{
   variable s = struct
     {
        object, close, eval_statistic,
        status, statistic, num_vary, num_points,
        response_type, data_type
     };

   s.data_type = qualifier_exists ("flux");
   % (data_type == 0) ? counts : flux

   s.response_type = qualifier ("response", Assigned_ARFRMF);
   % response = Ideal_ARF
   % response = Ideal_RMF
   % response = Ideal_ARF | Ideal_RMF
   % response = Ideal_ARFRMF

   s.close = &fit_object_close;
   s.eval_statistic = &fit_object_eval_statistic;

   _isis->_set_fit_type (s.response_type, s.data_type);
   s.object = _isis->open_fit_object_mmt_intrin ();

   return s;
}

%}}}

define register_slang_optimizer () %{{{
{
   variable msg = "register_slang_optimizer (name, &method [;qualifiers])\n"
                + "   qual:  set_options=&ref, stat_name=name";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable name, method, set_options, stat_name;
   (name, method) = ();
   set_options = qualifier ("set_options", NULL);
   stat_name = qualifier ("stat_name", "chisqr");

   _isis->_add_slang_optimizer (set_options, method, name, stat_name);
}

%}}}

