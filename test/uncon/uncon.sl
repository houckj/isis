% -*- mode: SLang; mode: fold -*-

%   This file is part of ISIS, the Interactive Spectral Interpretation System
%   Copyright (C) 1998-2010 Massachusetts Institute of Technology
%
%   This software was developed by the MIT Center for Space Research under
%   contract SV1-61010 from the Smithsonian Institution.
%
%   Author:  John C. Houck  <houck@space.mit.edu>
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
%--------------------------------------------------------------------

% Purpose:  Provide a S-Lang interface to the 'uncon'
%           suite of optimizer tests from netlib.org
%  Author:  John C. Houck <houck@space.mit.edu>
%    Date:  12/2004

import ("uncon");

public variable Uncon_Functions;
public variable _Uncon_Symbol = Assoc_Type[];
public variable _Uncon_Struct_Def = Assoc_Type[];

static define load_names (file) %{{{
{
   variable fp = fopen (file, "r");
   if (fp == NULL)
     {
        vmessage ("*** Error: uncon: failed opening file %s", file);
        return NULL;
     }
   variable names = fgetslines (fp);
   names = array_map (String_Type, &strtrim, names);
   () = fclose (fp);

   return names;
}

%}}}

static define load_symbol (name, file) %{{{
{
   variable symbol;
   % FIXME: the underscore is system-dependent
   if (-1 == _load_uncon_fun (&symbol, file, name + "_"))
     return NULL;
   _Uncon_Symbol[name] = symbol;
}

%}}}

define uncon_struct () %{{{
{
   variable u = struct
     {
        x, n, f, m, ftf, mode
     };

   variable big_enough = Double_Type[100];

   u.x = @big_enough;
   u.f = @big_enough;
   u.ftf = 0.0;
   u.n = 0;
   u.m = 0;
   u.mode = 0;

   return u;
}

%}}}

static define trim_uncon_struct (_u) %{{{
{
   variable u = @_u;
   u.x = u.x[[0:u.n-1]];
   u.f = u.f[[0:u.m-1]];
   return u;
}

%}}}

define _uncon_massage_hook (u, ref) %{{{
{
   return trim_uncon_struct (_uncon_hook (u, ref));
}

%}}}

define _uncon_init (name) %{{{
{
   variable ref = __get_reference (sprintf ("_uncon_%s", name));
   variable u = uncon_struct();
   variable modes = [UNCON_SET_SIZE, UNCON_SET_INITPT];
   foreach (modes)
     {
        u.mode = ();
        u = @ref(u);
     }
   return u;
}

%}}}

static variable Num_Function_Evaluations;

define uncon_setup_fit (name) %{{{
{
   variable u = _uncon_init (name);
   variable d = struct
     {
        bin_lo, bin_hi, value, err
     };

   delete_data (all_data);
   d.bin_lo = [1:u.m];
   d.bin_hi = make_hi_grid (d.bin_lo);
   d.value = Double_Type[u.m];
   d.err = ones(u.m);

   () = define_counts (d);

   fit_fun (name + "(1)");
   set_par ([1:u.n], u.x);

   Num_Function_Evaluations = 0;
}

%}}}

define _uncon_eval (p, name) %{{{
{
   variable ref = __get_reference (sprintf ("_uncon_%s", name));
   variable u = _Uncon_Struct_Def[name];
   u.mode = UNCON_COMPUTE;
   u.x = @p;
   u = @ref(u);
   Num_Function_Evaluations++;
   return u.f;
}

%}}}

define uncon_best (name) %{{{
{
   variable u = _uncon_init (name);
   variable ref = __get_reference (sprintf ("_uncon_%s", name));
   u.mode = UNCON_SET_BESTF;
   % When available, the initpt values for u.x will be
   % overwritten by by the best-fit values.
   u = @ref (u);
   return u;
}

%}}}

define uncon_num_function_evaluations () %{{{
{
   return Num_Function_Evaluations;
}

%}}}

static define init_uncon_fun (name) %{{{
{
   eval(sprintf ("public define _uncon_%s(u){return _uncon_massage_hook(u,_Uncon_Symbol[\"%s\"]);}",
                  name, name));
   eval(sprintf("_Uncon_Struct_Def[\"%s\"] = _uncon_init (\"%s\")",
                name, name));
   eval(sprintf ("public define %s_defaults(i){variable u=_Uncon_Struct_Def[\"%s\"]; return (u.x[i],0,0,0);}",
                 name, name));
   eval(sprintf ("public define %s_fit(l,h,p){return _uncon_eval(p,\"%s\");}",
                 name, name));
   variable u = _Uncon_Struct_Def[name];
   variable param_names = "x" + array_map (String_Type, &string, [0:u.n-1]);
   param_names = array_map (String_Type, &make_printable_string, param_names);
   eval(sprintf ("add_slang_function (\"%s\", [%s])",
                 name, strjoin (param_names, ",")));
}

%}}}

static define init_uncon (file) %{{{
{
   variable names;

   names = load_names (file);
   if (names == NULL)
     return;

   array_map (Void_Type, &load_symbol, names, "src/libuncon.so");
   array_map (Void_Type, &init_uncon_fun, names);

   return names;
}

%}}}

Uncon_Functions = init_uncon ("src/modules.lis");

provide ("uncon");
