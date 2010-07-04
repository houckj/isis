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

% $Id: aped_fun.sl,v 1.8 2004/06/14 14:38:31 houck Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:  Automatically generate a S-Lang fit-function
%           to compute a multi-component APED model with a
%           parameter list defined by a single Struct_Type.
%
% Public interface:
%    create_aped_fun (name, Struct_Type [, hook_ref])
%    name = function name
%
% The Struct_Type has the form:
%
% variable x = struct
%              {
%                norm, temperature, density,
%                elem, elem_abund
%                metal_abund,
%                vturb, redshift
%              }
%
% where    norm[] = float array of size Num_Components
%   temperature[] = float array of size 1 or Num_Components
%       density[] = float array of size 1 or Num_Components
%          elem[] = NULL or int array of size Num_Elements
%    elem_abund[] = NULL or float array of size Num_Elements
%   metal_abund = float
%         vturb = float
%      redshift = float
%
% All components in the resulting model have the same
% values of metal_abund, vturb, redshift, elem, elem_abund.
%
% Because the element abundances are referenced as separate
% scalar parameters in the fitting engine but as two vectors
% in the model evaluator, some considerable juggling of
% the abundances is required.
%
% The low-level problem solved is in mapping between
% these two forms:
%
% (I) Struct_Type {
%     norm = Double_Type[3]
%     temperature = Double_Type[3]
%     density = 1
%     metal_abund = 0.5
%     elem_abund = Double_Type[3] = [1.0, 3.1, 0.0];
%     elem = Integer_Type[3] = [Ne, Mg, Fe]
%     vturb = 100
%     redshift = 0
% }
%
% (II) Parameter List:
%
% xaped(1)
%  idx  param              tie-to  freeze  value    min   max
%   1  xaped(1).norm1          0     1          1      0     0
%   2  xaped(1).norm2          0     1          1      0     0
%   3  xaped(1).norm3          0     1          1      0     0
%   4  xaped(1).temperature1   0     1    4333333      0     0
%   5  xaped(1).temperature2   0     1    7666667      0     0
%   6  xaped(1).temperature3   0     1    1.1e+07      0     0
%   7  xaped(1).density        0     1          1      0     0
%   8  xaped(1).vturb          0     1        100      0     0
%   9  xaped(1).redshift       0     1          0      0     0
%  10  xaped(1).metal_abund    0     1        0.5      0     0
%  11  xaped(1).Ne             0     1          1      0     0
%  12  xaped(1).Mg             0     1        3.1      0     0
%  13  xaped(1).Fe             0     1          0      0     0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private define array_to_struct (fields) %{{{
{
   eval (sprintf ("define __atos__(){return struct {%s};}",
                  strjoin (fields, ",")));
   return eval ("__atos__()");
}

%}}}

private define concat_structs () %{{{
{
   variable a, x, n = _NARGS;

   if (n == 0)
     return NULL;

   a = __pop_args (n);

   if (n == 1)
     return a;

   a = [__push_args (a)];

   variable names, s;

   names = get_struct_field_names (a[0]);
   foreach (a[[1:n-1]])
     {
	x = ();
	names = [names, get_struct_field_names(x)];
     }

   s = array_to_struct (names);
   foreach (a)
     {
	x = ();
	foreach (get_struct_field_names (x))
	  {
	     n = ();
	     set_struct_field (s, n, get_struct_field (x, n));
	  }
     }

   return s;
}

%}}}

% Internal global data structures

% must match default_plasma_state() struct, sans elem/elem_abund.
private variable Struct_Field_Names =
 ["norm", "temperature", "density", "vturb", "redshift", "metal_abund"];

variable _Isis_Aped_Model_Table = Assoc_Type[Struct_Type];
variable _Isis_Aped_Profile_Table = Assoc_Type[Line_Profile_Type];
variable _Isis_Aped_Line_Modifier_Table = Assoc_Type[Ref_Type];

private variable Element_Names = %{{{
[
   "",
   "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
   "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca",
   "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni"
];

%}}}

private variable Max_Element_Name_Index = length(Element_Names) - 1;
private variable Elements = Assoc_Type[Int_Type];

private define init_element_names () %{{{
{
   foreach ([0:length(Element_Names)-1])
     {
	variable i = ();
	Elements [Element_Names[i]] = i;
     }
}

%}}}
init_element_names();

private variable Abund_Prefix;

private define abund_name (Z) %{{{
{
   if (Max_Element_Name_Index < Z)
     verror ("Unknown element, Z = %d", Z);

   variable s = Element_Names[Z];
   return (Abund_Prefix == NULL) ? s : Abund_Prefix + s;
}

%}}}

private define parse_abund_name (s) %{{{
{
   return (Abund_Prefix == NULL) ? s : strreplace (s, Abund_Prefix, "");
}

%}}}

private define get_abund_names (x) %{{{
{
   if (x.elem == NULL)
     return NULL;

   return array_map (String_Type, &abund_name, x.elem);
}

%}}}

private define abund_arrays_to_fields (x) %{{{
{
   variable i, n, names, s;

   names = get_abund_names (x);
   if (names == NULL)
     return NULL;

   s = array_to_struct (names);
   n = length(names);

   _for (0, n-1, 1)
     {
	i = ();
	set_struct_field (s, names[i], x.elem_abund[i]);
     }

   return s;
}

%}}}

private define abund_fields_to_arrays (y, names) %{{{
{
   if (names == NULL)
     return NULL;

   variable n = length(names);
   variable x = struct
     {
	elem, elem_abund
     };

   x.elem = Int_Type [n];
   x.elem_abund = Double_Type[n];

   variable i, Z;

   _for (0, n-1, 1)
     {
	i = ();
	Z = Elements[parse_abund_name(names[i])];

	x.elem[i] = Z;
	x.elem_abund[i] = get_struct_field (y, names[i])[0];
     }

   return x;
}

%}}}

private define sf_length (x, f) %{{{
{
   return length (get_struct_field (x, f));
}

%}}}

private define check_model_size (x) %{{{
{
   variable n, len;

   n = length(x.norm);
   if (n == 0) return -1;

   len = array_map (Int_Type, &sf_length, x, Struct_Field_Names);

   variable num_var, num_vec, num_scalar;

   num_var = length(len);
   num_scalar = length(where(len == 1));

   if (n == 1) num_vec = 0;
   else num_vec = length(where(len == n));

   if (num_var != num_vec + num_scalar)
     {
	error ("Invalid model specification");
	return -1;
     }

   variable n_elem, n_abund;

   if (x.elem == NULL) n_elem = 0;
   else n_elem = length(x.elem);

   if (x.elem_abund == NULL) n_abund = 0;
   else n_abund = length(x.elem_abund);

   if (n_elem != n_abund)
     {
	error ("elemental abundance arrays must have equal length");
	return -1;
     }

   return n;
}

%}}}

private define append_num (s, n) %{{{
{
   if (n <= 1)
     return s;

   return s + array_map (String_Type, &string, [1:n]);
}

%}}}

private define expanded_struct_field_name_list (x) %{{{
{
   variable abund_names = get_abund_names (x);

   if (abund_names == NULL)
     return Struct_Field_Names;

   return [Struct_Field_Names, abund_names];
}

%}}}

private define make_param_names (x) %{{{
{
   variable f, n, names, fields, abund_names;

   names = NULL;

   fields = expanded_struct_field_name_list (x);

   foreach (fields)
     {
	f = ();
	n = append_num (f, length(get_struct_field (x, f)));
	if (names == NULL)
	  names = n;
	else names = [names, n];
     }

   return names;
}

%}}}

private define array_assign_struct_field (structs, field, a) %{{{
{
   variable i, n = length(structs);
   variable val;

   if (length(a) == 1)
     {
	val = Double_Type[n];
	val[*] = a[0];
     }
   else val = a;

   _for (0, n-1, 1)
     {
	i = ();
	set_struct_field (structs[i], field, val[i]);
     }
}

%}}}

private define set_defaults () %{{{
{
   variable k = ();
   return default_plasma_state();
}

%}}}

define mt_create_from_struct () %{{{
{
   variable msg = "Model_Type = mt_create_from_struct (Struct_Type)";
   variable n, x;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   x = ();

   n = check_model_size (x);
   if (n < 1)
     return;

   variable structs = Struct_Type[n];
   structs = array_map (Struct_Type, &set_defaults, [1:n]);

   variable f;

   foreach (Struct_Field_Names)
     {
	f = ();
	array_assign_struct_field (structs, f, get_struct_field (x, f));
     }

   if (x.elem != NULL)
     {
	foreach (structs)
	  {
	     variable s = ();
	     s.elem = [x.elem];
	     s.elem_abund = [x.elem_abund];
	  }
     }

   return mt_def_model (structs);
}

%}}}

private define append_abund_fields (x) %{{{
{
   variable abund_struct = abund_arrays_to_fields (x);
   if (abund_struct == NULL)
     return x;

   return concat_structs (x, abund_struct);
}

%}}}

private define param_array_to_struct (p, template) %{{{
{
   variable i, n;
   n = length(p);
   i = 0;

   % if necessary, add a struct field for each element
   variable x, fields;
   x = append_abund_fields (template);
   fields = expanded_struct_field_name_list (template);

   foreach (fields)
     {
	variable f = ();

	variable old, len;
	old = get_struct_field (x, f);
	len = length(old);

	set_struct_field(x, f, p[ i+[0:len-1] ]);
	i += len;
     }

   % if necessary, update the elemental abundance value array
   variable t = abund_fields_to_arrays (x, get_abund_names (template));
   if (t != NULL)
     x.elem_abund = [t.elem_abund];

   return x;
}

%}}}

private define struct_to_param_array (x) %{{{
{
   variable fields = expanded_struct_field_name_list (x);
   variable p = NULL;

   foreach (fields)
     {
	variable f = ();

	variable v, len;
	v = get_struct_field (x, f);

	if (v == NULL)
	  v = 0.0;

	if (p != NULL)
	  {
	     p = [p, v];
	  }
	else
	  {
	     if (typeof(v) == Array_Type)
	       p = (@v) * 1.0;
	     else
	       p = v * 1.0;
	  }
     }

   return p;
}

%}}}

% Generate line modifier args
define aped_line_modifier_args () %{{{
{
   variable msg = "args = aped_line_modifier_args (name, p [,num_extra_args])";
   variable name, p, num_extra_args = 0;

   switch (_NARGS)
     {case 2:  (name, p) = ();}
     {case 3:  (name, p, num_extra_args) = ();}
     { usage(msg); }

   variable ref = _Isis_Aped_Line_Modifier_Table[name];

   if (num_extra_args == 0)
     return 0, ref, p, "aped_line_modifier";

   variable extra_args = __pop_args (num_extra_args);
   return __push_args (extra_args), num_extra_args, ref, p, "aped_line_modifier";
}

%}}}

% Generate line profile args
define aped_line_profile_args (name, p) %{{{
{
   return _Isis_Aped_Profile_Table[name], p, "aped_profile";
}

%}}}

% Generate hook args
define aped_hook_args (hook_ref, which_aped_fun_instance) %{{{
{
   return @hook_ref (which_aped_fun_instance), "aped_hook";
}

%}}}

define _isis_calc_model_using_template (l,h,p,info) %{{{
{
   variable
     template = info.param_struct,
     hook_ref = info.hook_ref;

   variable m, s;
   s = param_array_to_struct (p, template);
   m = mt_create_from_struct (s);

   if (hook_ref != NULL)
     {
        % leave hook args on the stack.
        aped_hook_args (hook_ref, Isis_Active_Function_Id);
     }

   variable v = mt_calc_model (m, l, h);

   variable this_model =
     sprintf ("%s(%d)", info.model_name, Isis_Active_Function_Id);
   info.model_details[this_model] = mt_get_model_details (m);

   return v;
}
%}}}

define norm_indices (names) %{{{
{
   variable indices = [0:length(names)-1];
   foreach (indices)
     {
	variable i = ();
	if (0 == is_substr (names[i], "norm"))
	  indices[i] = -1;
     }
   return indices[where(indices >= 0)];
}

%}}}

define create_aped_fun () %{{{
{
   variable msg =
`create_aped_fun (name, Struct_Type [,hook_ref])
      qualifiers:    no_abund_prefix`;
   variable fun_name, param_struct, hook_ref = NULL;

   Abund_Prefix = qualifier_exists ("no_abund_prefix") ? NULL : "abund_";

   if (_isis->get_varargs (&fun_name, &param_struct, &hook_ref, _NARGS, 2, msg))
     return;

   if (check_model_size (param_struct) < 0)
     return;

   variable fun_fmt =
   "define %s_fit (l,h,p)"
  +"{"
  +    "variable s = _Isis_Aped_Model_Table[\"%s\"];"
  +    "return _isis_calc_model_using_template (l,h,p,s);"
  +"}";

   variable defaults_fmt =
   "define %s_defaults (i)"
  +"{"
  +   "variable v = _Isis_Aped_Model_Table[\"%s\"].param_values;"
  +   "variable d = struct {value, freeze, min, max};"
  +   "d.value = v[i];"
  +   "d.freeze = 1;"
  +   "d.min = 0.0;"
  +   "d.max = 0.0;"
  +   "return (d.value, d.freeze, d.min, d.max);"
  +"}";

   eval (sprintf (fun_fmt, fun_name, fun_name));
   eval (sprintf (defaults_fmt, fun_name, fun_name));

   variable xs, names;

   xs = append_abund_fields (param_struct);
   names = make_param_names (xs);
   add_slang_function (fun_name, names, norm_indices (names));

   variable set_hook_fmt;
   set_hook_fmt = "set_param_default_hook(\"%s\", \"%s_defaults\");";
   eval (sprintf(set_hook_fmt, fun_name, fun_name));

   % Save model specs in a global associative array,
   % indexed by the function name.
   variable s = struct
     {
	param_values, param_struct, hook_ref,
        model_name, model_details
     };

   s.param_struct = @param_struct;
   s.param_values = struct_to_param_array (xs);
   s.hook_ref = hook_ref;
   s.model_name = strtrim (fun_name);
   s.model_details = Assoc_Type[List_Type];

   _Isis_Aped_Model_Table[s.model_name] = s;
}

%}}}

define create_aped_line_profile () %{{{
{
   variable msg = "create_aped_line_profile (name, &Line_Profile_Type [, param_names[]])";
   variable profile_name, mmt, param_names = NULL;

   if (_isis->get_varargs (&profile_name, &mmt, &param_names, _NARGS, 2, msg))
     return;

   variable s;
   % returns id_string, param_array, arg, ...
   s = "define ${profile_name}_fit(l,h,p){"$
     + "  return aped_line_profile_args (\"${profile_name}\", p);"$
     + "}";
   eval(s);
   if (param_names != NULL)
     add_slang_function (profile_name, param_names);
   else
     add_slang_function (profile_name);
   _Isis_Aped_Profile_Table[profile_name] = mmt;
}

%}}}

define create_aped_line_modifier () %{{{
{
   variable msg = "create_aped_line_modifier (name, &ref, param_names[] [,num_extra_args])";
   variable modifier_name, modifier_ref, param_names = NULL;
   variable num_extra_args = 0;

   if (_isis->get_varargs (&modifier_name, &modifier_ref, &param_names,
                           &num_extra_args, _NARGS, 2, msg))
     return;

   variable s;
   % returns id_string, param_array, arg, ...
   s = "define ${modifier_name}_fit(l,h,p){"$
     + "     return aped_line_modifier_args (\"${modifier_name}\", p, ${num_extra_args});"$
     + "}";
   eval(s);

   if (param_names != NULL)
     add_slang_function (modifier_name, param_names);
   else
     add_slang_function (modifier_name);

   _Isis_Aped_Line_Modifier_Table[modifier_name] = modifier_ref;
}

%}}}

% Provide access to line emissivities from each (n,T) component
define aped_fun_details () %{{{
{
   variable msg = "List_Type = aped_fun_details (\"model(k)\")";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable name = ();

   name = strcompress(name, " ");
   variable model = strtok (name, "(")[0];
   variable s = _Isis_Aped_Model_Table[model];

   if (assoc_key_exists (s.model_details, name))
     return s.model_details[name];

   return NULL;
}

%}}}
