% -*- mode: SLang; mode: fold -*-

%   This file is part of ISIS, the Interactive Spectral Interpretation System
%   Copyright (C) 1998-2018 Massachusetts Institute of Technology
%
%   Author:  John C. Houck <houck@space.mit.edu>
%
%   This software was developed by the MIT Center for Space Research under
%   contract SV1-61010 from the Smithsonian Institution.
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

#ifndef __XSPEC_STATIC_LINKED__
#ifdef __XSPEC__
#stop
#endif
import ("xspec");
#endif

private variable Model_File;    % path to .so file

private variable shared_lib_ext = "so";
#ifdef __APPLE__
  shared_lib_ext = "dylib";
#endif

require ("readascii");
require ("parse_model_dat");

private variable _Temp_Xspec_Link_Errors, _Xspec_Link_Errors;
private define register_link_error (assoc, errmsg, name)
{
   if (assoc_key_exists (assoc, errmsg))
     list_append (assoc[errmsg], name);
   else assoc[errmsg] = {name};
}
private define xspec_keep_link_errors ()
{
   variable errmsg, namelist;
   foreach errmsg, namelist (_Temp_Xspec_Link_Errors) using ("keys", "values")
     {
        variable n;
        foreach n (namelist)
          {
             register_link_error (_Xspec_Link_Errors, errmsg, n);
          }
     }
}
private define xspec_clear_link_errors ()
{
   _Temp_Xspec_Link_Errors = Assoc_Type[];
}
private define xspec_track_link_errors ()
{
   _Xspec_Link_Errors = Assoc_Type[];
}
private define xspec_print_link_errors ()
{
   variable errmsg, namelist;
   foreach errmsg, namelist (_Xspec_Link_Errors) using ("keys", "values")
     {
        vmessage ("Link error: $errmsg"$);
     }
}

%=======================================================================
%
% Put some important symbols in the _isis namespace so
% we always know where to find them:
%
$1 = current_namespace();
if ($1 == "") $1= "Global";

use_namespace("_isis");

static variable __Xspec_Symbol = Assoc_Type[];

static define _lmodel_set_default (value, freeze_val, hard_min, min, max, hard_max, step) %{{{
{
   variable t = struct {value, freeze, hard_min, hard_max, min, max, step, relstep};
   t.value = value;
   t.freeze = freeze_val;
   t.hard_min = hard_min;
   t.hard_max = hard_max;
   t.min = min;
   t.max = max;
   t.step = step;
   t.relstep = Isis_Default_Relstep;
   return t;
}

%}}}

static define xspec_register_link_error (name, path, errmsg)
{
   register_link_error (_Temp_Xspec_Link_Errors, errmsg, name);
}

use_namespace($1);

static define load_xspec_symbol (name) %{{{
{
   return load_xspec_fun (Model_File, name);
}

%}}}

static define parse_model_table ();
private define load_shared_model_table (buf)
{
   variable names;
   names = parse_model_table (buf, 1, &load_xspec_symbol, 0);
   return names;
}

private define create_param_defaults_array (m) %{{{
{
   variable name = sprintf ("__%s_Defaults", m.model_name);

   eval(sprintf ("variable %s = Struct_Type[%d];",
		name, m.num_pars));

   variable s, freeze_state, step;
   variable pi = m.par_info;
   variable indices = [0:m.num_pars-1];

   foreach (indices)
     {
	variable k = ();
        step = typecast (eval(pi[7,k]), Double_Type);
        if (step <= 0)
          {
             freeze_state = 1;
             step = abs(step);
          }
        else freeze_state = 0;

	s = sprintf ("%s[%d] = _isis->_lmodel_set_default(%s, %d, %s, %s, %s, %s, %e);",
		     name, k,
		     pi[2,k], freeze_state, pi[3,k], pi[4,k], pi[5,k], pi[6,k], step);
	eval(s);
     }

   return name;
}

%}}}

private define attach_unit_string (name, units) %{{{
{
   if (strlen(strtrim(units)) > 0)
     return sprintf ("%s [%s]", name, units);
   else
     return name;
}

%}}}

private variable __model_renamings = Assoc_Type[String_Type];

private define create_wrapper (m) %{{{
{
   if (m.model_type != "add"
       and m.model_type != "mul"
       and m.model_type != "con")
     return;

   variable arg;
   if (m.model_type == "con")
     arg = ",arg";
   else
     arg = "";

   variable set_init_string = "";
   if (m.init_string != NULL)
     {
        set_init_string =
          sprintf ("_xspec_model_init_string (_isis->__Xspec_Symbol[\"%s\"], \"%s\");",
                   m.model_name, m.init_string);
     }

   variable isis_model_name = assoc_key_exists (__model_renamings, m.model_name) ? __model_renamings[m.model_name] : m.model_name;
   variable def;
   def = sprintf ("define %s_fit(l,h,p%s){%s return %s(l,h,p%s, _isis->__Xspec_Symbol[\"%s\"]);}",
		  isis_model_name, arg, set_init_string,
		  m.exec_symbol_hook, arg,
		  m.model_name);

   variable param_names, units, s;

   param_names = m.par_info[0,*];
   units = m.par_info[1,*];

   param_names = array_map (String_Type, &attach_unit_string, param_names, units);
   s = array_map (String_Type, &make_printable_string, param_names);
   param_names = strjoin (s, ", ");

   variable decl;
   decl = sprintf ("add_slang_function (\"%s\", [%s]);",
		   isis_model_name, param_names);

   variable pdef_array_name = create_param_defaults_array (m);
   variable pdef;
   pdef = sprintf ("define %s_default(i){return %s[i];}",
                   isis_model_name, pdef_array_name);

   variable psetdef;
   psetdef = sprintf ("set_param_default_hook (\"%s\", &%s_default)",
                      isis_model_name, isis_model_name);
   eval(def);
   eval(decl);
   eval(pdef);
   eval(psetdef);

   if (m.model_type == "con")
     set_function_category (isis_model_name, ISIS_FUN_OPERATOR);

   return isis_model_name;
}

%}}}

static define choose_exec_symbol_hook (m) %{{{
{
   if (m.exec_symbol_hook != NULL)
     return;
   variable fmt = "_xspec_%s_f_hook";
   m.has_fortran_linkage = 1;

   if (strlen(m.routine_name) > 3)
     {
        variable prefix = m.routine_name[[0:1]];
        switch (prefix)
          {
           case "C_" or case "c_":
             m.routine_name = m.routine_name[[2:]];
             fmt = "_xspec_%s_C_hook";
             m.has_fortran_linkage = 0;
          }
          {
           case "F_":
             m.routine_name = m.routine_name[[2:]];
             fmt = "_xspec_%s_F_hook";
             m.has_fortran_linkage = 1;
          }
     }

   m.exec_symbol_hook = sprintf (fmt, m.model_type);
}

%}}}

static define mangle_name (m, name) %{{{
{
   variable s = name;

   if (m.has_fortran_linkage)
     {
        s += __FC_MANGLE_SUFFIX;
        s = __FC_MANGLE_UPCASE ? strup(s) : strlow(s);
     }

   return s;
}

%}}}

static define parse_model_table (t, mangle, load_symbol_ref, lookup_by_model_name) %{{{
{
   variable buf = t.buf;

   if (buf == NULL)
     return NULL;

   if (buf[length(buf)-1] != "")
     buf = [buf, ""];

   variable blanks = where (buf == "");
   if (length(blanks) == 0) blanks = [0];

   variable e, m, s;
   variable names = {};

   xspec_track_link_errors ();

   s = 0;

   foreach (blanks)
     {
	e = ();

	m = parse_function (buf[[s:e-1]]);
	s = e+1;

	if (m == NULL)
	  continue;

        !if (assoc_key_exists (Supported_Model_Types, m.model_type))
          continue;

	if (m.num_pars == 0)
	  continue;

        choose_exec_symbol_hook (m);
        xspec_clear_link_errors ();

        variable lookup_name = (lookup_by_model_name ?
                                m.model_name : m.routine_name);
        if (mangle)
          {
             lookup_name = mangle_name (m, lookup_name);
          }

        m.exec_symbol = (@load_symbol_ref) (lookup_name);
	if (m.exec_symbol == NULL)
          {
             % save recent link errors only if load ultimately failed
             xspec_keep_link_errors();
             vmessage ("Error loading symbol %s for model %s",
                       lookup_name, m.model_name);
             continue;
          }

	_isis->__Xspec_Symbol[m.model_name] = m.exec_symbol;
	variable isis_model_name = create_wrapper (m);
	list_append (names, isis_model_name);
     }

   if (Isis_Verbose > 0)
     {
	variable list = length (names) ? strjoin (list_to_array (names), ", ") : "(null)";
	vmessage ("loaded:  %S", list);
     }

   xspec_print_link_errors ();

   return names;
}

%}}}

#ifdef __HAVE_XSPEC_12__
private variable Lmoddir;

private define set_local_model_dir (dir)
{
   Lmoddir = dir;
}

private define get_local_model_dir ()
{
   variable dir;
   if (__is_initialized (&Lmoddir))
     dir = Lmoddir;
   else
     {
        dir = getenv ("LMODDIR");
     }
   return dir;
}

define build_xspec_local_models () %{{{
{
   variable msg = "build_xspec_local_models ([dir [, pkg_name]][;lmodel=filename]);";
   variable lmoddir = get_local_model_dir ();
   variable pkg_name = "xspeclocal";

   switch (_NARGS)
     {
      case 0:
        % do nothing
     }
     {
      case 1:
        lmoddir = ();
     }
     {
      case 2:
        (lmoddir, pkg_name) = ();
     }
     {
        % default:
        _pop_n (_NARGS);
        usage(msg);
        return;
     }

   if (-1 == access (lmoddir, F_OK | R_OK | W_OK))
     {
        throw ApplicationError, "*** Error: No read/write access to $lmoddir"$;
     }

   variable lmodeldat = qualifier ("lmodel", "lmodel.dat");

   variable cmd;

   cmd = ". $HEADAS/headas-init.sh ; initpackage $pkg_name $lmoddir/$lmodeldat $lmoddir -udmget"$;
   if (0 != system (cmd))
     {
        vmessage ("*** Error while running 'initpackage'");
        vmessage ("    The failed command line was:");
        vmessage (cmd);
        throw ApplicationError;
     }

   cmd = ". $HEADAS/headas-init.sh ; cd $lmoddir ; hmake"$;
   if (0 != system (cmd))
     {
        throw ApplicationError,
          "*** Error while running 'hmake' in $lmoddir"$;
     }
}

%}}}

private define make_lib_path (dir, name) %{{{
  {
   variable path = path_concat (dir, "lib${name}.${shared_lib_ext}"$);

   if (NULL == stat_file (path))
     {
        vmessage ("*** Error:  cannot stat library $path"$);
        return NULL;
     }

   return path;
}

%}}}

% Try to find the shared library in each of the following cases:
%  - given {dir,name}, look for $dir/lib$name.$ext
%  - given dir, look for $dir/lib$name.$ext
%  - given name, look for $LMODDIR/lib$name.$ext
%  - given <nothing>, look for $LMODDIR/lib$name.$ext
%
private define guess_lib_path (nargs) %{{{
{
   variable s, dir = get_local_model_dir();

   if (nargs == 2)
    {
        (dir, s) = ();

	return make_lib_path (dir, s);
     }
   else if (nargs == 1)
     {
        s = ();
        variable st = stat_file (s);
        if (st != NULL)
          {
             if (stat_is ("reg", st.st_mode))
               return s;
             else if (stat_is ("dir", st.st_mode))
               dir = s;
          }
        else if (dir != NULL)
          {
             return make_lib_path (dir, s);
          }
     }

   if (dir == NULL)
     return NULL;

   variable names = glob ("${dir}/*.${shared_lib_ext}"$);
   if (length(names) == 1)
     return names[0];

   return NULL;
}

%}}}

define load_xspec_local_models () %{{{
{
   if (_NARGS > 2)
     {
        _pop_n (_NARGS);
        usage ("load_xspec_local_models ([dir [, pkg_name]][;lmodel=filename]);");
        return;
     }

   variable lmodeldat = qualifier ("lmodel", "lmodel.dat");

   variable lib_path = guess_lib_path (_NARGS);
   if (lib_path == NULL)
     return;

   Model_File = lib_path;
   vmessage ("Loading %s with parameter file %s", Model_File, lmodeldat);

   variable dir = path_dirname (lib_path);
   variable buf = load_buf (path_concat (dir, lmodeldat));
   () = load_shared_model_table (buf);
}

%}}}

private define try_loading_local_models ()
{
   load_xspec_local_models ();
}

#else

private define get_model_lib_path (dir, ext) %{{{
{
   if (NULL == stat_file (dir))
     {
	vmessage ("*** xspec:  Directory %S not found", dir);
	return NULL;
     }

   % Since the release version number is embedded in the name
   % of the shared library, we can't search for a fixed library
   % name.  Instead, search for the shared library using a
   % regular expression for the name.  Thank you HEASARC.

   variable names = listdir (dir);
   if (names == NULL)
     {
	vmessage ("*** xspec:  Directory %s is empty", dir);
	return NULL;
     }

   variable root_name = "libxspec_lfn";
   variable pat = sprintf ("%s.*%s", root_name, ext);
   variable i = array_map (Integer_Type, &string_match, names, pat, 1);

   variable len = length(i);
   if (len == 0)
     return NULL;

   names = names[where(i)];
   if (length(names) == 0)
     return NULL;

   if (len > 1)
     {
	variable lib = root_name + ext;
	if (any(names == lib))
	  return path_concat (dir, lib);
     }

   return path_concat (dir, names[0]);
}

%}}}

private define load_lmodels_from_dir (dir, file) %{{{
{
   Model_Dat = path_concat (dir, file);

   Model_File = get_model_lib_path (dir, ".${shared_lib_ext}"$);
   if (Model_File == NULL)
     {
        vmessage ("*** xspec: Did not find a shared library in directory %s", dir);
        return NULL;
     }
   vmessage ("Loading %s", Model_File);

   return load_shared_model_table (load_buf(Model_Dat));
}

%}}}

define load_xspec_local_models () %{{{
{
   variable dirs, env;

   switch (_NARGS)
     {
      case 0:
	env = getenv ("LMODDIR");
	if (env == NULL)
	  return;
	dirs = strtok(env, ":");
     }
     {
	% default:
	dirs = __pop_args (_NARGS);
	dirs = [__push_args (dirs)];
     }

   foreach (dirs)
     {
	variable d = ();
	() = load_lmodels_from_dir (d, "lmodel.dat");
     }
}

%}}}

define _load_xspec_local_models () %{{{
{
   variable msg = "_load_xspec_local_models (libso_path, lmodeldat_path)";
   variable libso_path, lmodeldat_path;

   if (_NARGS != 2)
     {
	message (msg);
	_pop_n (_NARGS);
	return;
     }

   (libso_path, lmodeldat_path) = ();

   Model_Dat = lmodeldat_path;
   Model_File = libso_path;

   () = load_shared_model_table (load_buf(Model_Dat));
}

%}}}

private define try_loading_local_models () %{{{
{
   % check for local models
   variable env = getenv ("LMODDIR");
   if (env == NULL)
     return;
   variable dirs = strtok (env, ":");

   if (length(dirs) == 0)
     return;

   if (length(dirs) > 1)
     {
        load_xspec_local_models ();
        return;
     }

   % A single local model directory might contain either a
   % shared or a static library.

   if (NULL != get_model_lib_path (dirs[0], ".${shared_lib_ext}"$))
     {
        load_xspec_local_models ();
        return;
     }

   if (NULL == get_model_lib_path (dirs[0], ".a"))
     return;

   % If it's a static library, we'll assume it has already
   % been linked.
   Model_Dat = path_concat (dirs[0], "lmodel.dat");
   if (NULL != stat_file (Model_Dat))
     {
        () = parse_model_table (load_buf (Model_Dat), 0, &find_xspec_fun, 0);
     }
}

%}}}
#endif

private define load_xspec_models () %{{{
{
   % load static models first
   Model_Dat = find_model_dat_file ("$HEADAS"$, Xspec_Version);
   () = parse_model_table (load_buf (Model_Dat), 0, &find_xspec_fun, 1);

   try_loading_local_models ();
}

%}}}

% In order to rename an xspec model, define the following function: %{{{
%
% define xspec_rename_model_hook (isis_model_name)
% {
%    isis_model_name[ "old-xspec-name" ] = "new-isis-name";
% }
if (2 == is_defined ("xspec_rename_model_hook"))
{
   __model_renamings;  % left on stack
   eval ("xspec_rename_model_hook");
}
%}}}

load_xspec_models ();

% define lower-case aliases for back-compatibility
private define make_lcase_aliases () %{{{
{
   variable names = ["compLS", "compPS", "compST", "compTT", "nthComp",
                     "SSS_ice", "TBabs", "TBgrain", "TBvarabs", "zTBabs"];
   foreach (names)
     {
        variable n = ();
        if (2 == is_defined (n)) alias_fun (n, strlow(n));
     }
}
make_lcase_aliases ();
%}}}

private define table_param_default (file, name, is_norm) %{{{
{
   variable fcn_fmt =
`   define ${name}_default_hook(i,is_norm,pval,pmin,pmax)
    {
       variable t = struct
          {
             value, min, max,
             freeze = 0,
             hard_min = -_Inf, hard_max = _Inf,
             step = 0, relstep = Isis_Default_Relstep
          };
       if (is_norm) i--;
       if (i < 0)
         {
            % the norm parameter, eh?
            t.value = 1.0;  t.min = 0.0;   t.max = 1.e10;
         }
       else if (i >= length(pval))
         {
            % the redshift, I presume.
            t.value = 0.0;  t.min = -10.0; t.max = 10.0;
         }
       else
         {
            t.value = pval[i]; t.min = pmin[i]; t.max = pmax[i];
         }
       return t;
     }`$;
   variable reg_fmt =
`    variable t = fits_read_table ("${file}");
     set_param_default_hook ("${name}", &${name}_default_hook, ${is_norm},
                             t.initial, t.minimum, t.maximum);`$;
   eval(fcn_fmt);
   eval(reg_fmt);
}

%}}}

private define _add_table_model (nargs, msg, fmt, args_hook, is_norm) %{{{
{
   variable file, name, args;

   if (nargs != 2)
     {
	_pop_n (nargs);
	usage (msg);
	return;
     }

   (file, name) = ();

   args = fits_read_col (file, "name");
   if (args == NULL)
     {
	verror ("failed reading parameter names from %S", file);
	return;
     }

   if (args_hook != NULL)
     {
	if (-1 == @args_hook (file, &args))
	  {
	     _pop_n (nargs);
	     usage (msg);
	     return;
	  }
	if (typeof (args) == String_Type) args = [args];
	args = array_map (String_Type, &str_delete_chars, args, " ");
     }

   eval (sprintf (fmt, name, file));
   add_slang_function (name, args);
   table_param_default (file, name, is_norm);
}

%}}}

private define redshift_hook (file, args_ref) %{{{
{
   variable z, args;

   variable fp = fits_open_file (file, "r");
   if (fp == NULL)
     return -1;

   % read logical keyword value
   if (0 == _fits_read_key_integer (fp, "REDSHIFT", &z, NULL))
     {
	fp = NULL;   % close the file

	args = @args_ref;
        if (z == 1) args = [args, "redshift"];

	variable dims, num_dims;
	(dims, num_dims, ) = array_info(args);
	if (num_dims == 2)
	  {
	     reshape (args, dims[0]*dims[1]);
	     @args_ref = args;
	  }
	else if (num_dims > 2)
	  {
	     return -1;
	  }
     }

   fp = NULL;  % close the file

   @args_ref = args;

   return 0;
}

%}}}

private define atable_hook (file, args_ref) %{{{
{
   variable args;

   if (-1 == redshift_hook (file, args_ref))
     return -1;

   args = ["norm", @args_ref];
   @args_ref = args;

   return 0;
}

%}}}

define add_atable_model () %{{{
{
   variable msg = "add_atable_model (file, name);";
   variable fmt = "define %s_fit(l,h,p){_set_table_model_filename(\"%s\"); variable n = length(p); if (n > 1) return p[0]*_atbl(l,h,p[[1:n-1]]); else return p[0] * _atbl(l,h,NULL);}";

   _add_table_model (_NARGS, msg, fmt, &atable_hook, 1);
}

%}}}

define add_mtable_model () %{{{
{
   variable msg = "add_mtable_model (file, name);";
   variable fmt = "define %s_fit(l,h,p){_set_table_model_filename(\"%s\");return _mtbl(l,h,p);}";

   _add_table_model (_NARGS, msg, fmt, &redshift_hook, 0);
}

%}}}

define add_etable_model () %{{{
{
   variable msg = "add_etable_model (file, name);";
   variable fmt = "define %s_fit(l,h,p){_set_table_model_filename(\"%s\");return exp(-_mtbl(l,h,p));}";

   _add_table_model (_NARGS, msg, fmt, &redshift_hook, 0);
}

%}}}

define xspec_abund () %{{{
{
   variable msg =
`Usage:
  table_name = xspec_abund();
  status = xspec_abund (name)
            where name is an xspec table name or
            where name has the form "file ascii_filename"
  status = xspec_abund (Float_Type[30])
`;

   if (_NARGS == 0)
     return _xs_get_abundances ();

   variable arg = ();

   if (typeof(arg) == Array_Type)
     return _xs_set_abundances (typecast (arg, Float_Type), "file");

   if (typeof(arg) != String_Type)
     usage (msg);

   if (0 != strncmp (arg, "file", 4))
     return _xs_set_abundances (NULL, arg);

   variable abund, n, file = arg[[5:]];
   if (readascii (file, &abund; format="%f") < 0)
     throw IOError, "*** error reading $file"$;

   return _xs_set_abundances (abund, "file");
}

%}}}

define xspec_xsect () %{{{
{
   if (_NARGS == 1)
     _xs_set_xsections ();
   else
     _xs_get_xsections ();
}

%}}}

define xspec_photo () %{{{
{
   variable msg = "xsec[] = xspec_photo (kev1[], kev2[], Z [,version])";
   variable kev1, kev2, z, version=3;

   switch (_NARGS)
     {case 4: (kev1,kev2,z,version) = ();}
     {case 3: (kev1,kev2,z) = ();}
     {%default
        usage(msg);
     }
   if (typeof(kev1) == Array_Type)
     {
        array_map (Double_Type, &_xs_photo,
                   typecast(kev1,Float_Type),
                   typecast(kev2,Float_Type), z, version);
     }
   else
     {
        _xs_photo (typecast(kev1,Float_Type),
                   typecast(kev2,Float_Type), z, version);
     }
}

%}}}

define xspec_gphoto () %{{{
{
   variable msg = "xsec[] = xspec_gphoto (kev1[], kev2[], Z)";
   variable kev1, kev2, z;

   if (_NARGS != 3)
     usage(msg);

   (kev1,kev2,z) = ();

   if (typeof(kev1) == Array_Type)
     {
        array_map (Double_Type, &_xs_gphoto,
                   typecast(kev1,Float_Type), typecast(kev2,Float_Type), z);
     }
   else
     {
        _xs_gphoto (typecast(kev1,Float_Type),typecast(kev2,Float_Type), z);
     }
}

%}}}

define xspec_phfit2 () %{{{
{
   variable msg = "xsec_Mb = xspec_phfit2 (Z, Ne, shell, E_eV[])";
   variable z, ne, shell, e;

   if (_NARGS != 4)
     usage(msg);

   (z, ne, shell, e) = ();

   if (typeof(e) == Array_Type)
     {
        array_map (Double_Type, &_xs_phfit2, z, ne, shell,
                   typecast(e,Float_Type));
     }
   else
     {
        _xs_phfit2 (z, ne, shell, typecast(e,Float_Type));
     }
}

%}}}

#ifexists _xs_ionsneqr
define xspec_ionsneqr () %{{{
{
   variable msg = "Struct_Type = xspec_ionsneqr(T[], tau[]);  %% T [K], tau [s/cm^3]";
   variable temp, tau;

   if (_NARGS != 2)
     usage(msg);

   (temp, tau) = ();

   variable s = struct
     {
        fout, ionel, ionstage
     };

   (s.fout, s.ionel, s.ionstage) = _xs_ionsneqr (temp, tau);
   return s;
}

%}}}
#endif

#ifexists _xs_noneq
define xspec_noneq () %{{{
{
   variable msg =
`Struct_Type = xspec_noneq(tempr, tau, weight);  %% T [K], tau [s/cm^3]
               tempr = Float_Type[ntmp]
         tau, weight = Float_Type[ntau]
 Note: the 'weight' array will be normalized so that sum(weight)=1
`;
   variable temp, tau, weight;

   if (_NARGS != 3)
     usage(msg);

   (temp, tau, weight) = ();

   % normalize weights.
   weight /= sum(weight);

   % Stupid routine *requires* two temperatures
   % or it doesn't work at all.
   if (length (temp) == 1)
     {
        temp = [temp, temp];
     }

   variable s = struct
     {
        fout, ionel, ionstage
     };

   (s.fout, s.ionel, s.ionstage) = _xs_noneq (temp, tau, weight);
   return s;
}

%}}}
#endif

define xspec_elabund () %{{{
{
   if (_NARGS == 0)
     usage ("X = xspec_elabund (\"element_name\")");
   variable elname = ();
   _xs_get_element_solar_abundance (elname);
}

%}}}

define xspec_xset () %{{{
{
   variable msg = "xspec_xset (\"symbol\", \"value\")";

   if (_NARGS != 2)
     usage (msg);

   _xs_fpmstr ();
}

%}}}

define xspec_set_cosmo () %{{{
{
   variable msg = "xspec_set_cosmo (H0, q0, lambda0)";
   variable h0, q0, l0;

   if (_NARGS != 3)
     usage (msg);

   (h0, q0, l0) = ();

   _xs_set_cosmo_hubble (h0);
   _xs_set_cosmo_decel (q0);
   _xs_set_cosmo_lambda (l0);
}

%}}}

define xspec_get_cosmo () %{{{
{
   variable msg = "(h0, q0, lambda0) = xspec_get_cosmo ()";

   return (_xs_get_cosmo_hubble (),
	   _xs_get_cosmo_decel (),
	   _xs_get_cosmo_lambda ());
}

%}}}

xspec_set_cosmo (70.0, 0.0, 0.73);  % H0, q0, lambda0

define xspec_get_chatter () %{{{
{
   _xs_gchat ();
}

%}}}

define xspec_set_chatter () %{{{
{
   variable msg = "xspec_set_chatter (level);";

   if (_NARGS != 1)
     usage(msg);

   _xs_pchat ();
}

%}}}

% xspec 12's default is too chatty
xspec_set_chatter(5);

%
% -------- Parse XSPEC help file
%

#ifdef __HAVE_XSPEC_12__
variable Xspec_Module_Help;

if (NULL != getenv ("DISPLAY"))
  Xspec_Module_Help = "pdf ; xpdf %s";
else
  Xspec_Module_Help = "html ; lynx %s";

define xspec_help () %{{{
{
   if (_NARGS != 1)
     {
	_pop_n (_NARGS);
	usage ("xspec_help (\"model\")");
        return NULL;
     }

   variable model = ();

   variable headas = getenv ("HEADAS");
   if (headas == NULL)
     return NULL;

   variable help_env = getenv ("XSPEC_MODULE_HELP");
   if (help_env != NULL)
     Xspec_Module_Help = help_env;

   variable s = strtok (qualifier ("method", Xspec_Module_Help), ";");
   variable format = strtrim(s[0], " "), cmd = s[1];

   variable file =
     sprintf ("XSmodel%c%s.%s",
              toupper(model[0]), model[[1:]], format);

   variable dir;
   dir = path_dirname (strtrim_end (headas, "/"));
   dir = path_concat (dir, "spectral/help");

   switch (format)
     { case "pdf":  }
     { case "html": dir += "/html"; }
     {
        % default:
        throw ApplicationError, "Unrecognized file format:  '$format'"$;
     }

   variable path = path_concat (dir, file);
   if (NULL == stat_file (path))
     return NULL;

   if (0 != system (sprintf (cmd, path)))
     {
        vmessage ("*** Unable to display help file '$file'"$);
        vmessage ("    See the documentation for 'xspec_help' for details");
        return NULL;
     }

   return 0;
}

%}}}

#else
private variable Xspec_Names = Assoc_Type [];

private define build_lookup_table () %{{{
{
   variable file = path_concat (path_dirname(__FILE__), Xspec_Model_Names_File);
   variable fp = fopen (file, "r");
   variable names = fgetslines (fp);
   () = fclose (fp);

   names = array_map (String_Type, &strtrim, names);

   foreach (names)
     {
        variable n = ();
        Xspec_Names[n] = n;
     }
}

%}}}

build_lookup_table ();

private define read_help_file (path) %{{{
{
   variable fp, buf;

   fp = fopen (path, "r");
   if (fp == NULL)
     return NULL;
   buf = fgetslines (fp);
   () = fclose (fp);

   return array_map (String_Type, &strtrim, buf);
}

%}}}

private define find_entry_end (buf, i) %{{{
{
   variable pat = "^3 [a-zA-Z0-9]+$";
   variable k, m, n = length (buf);

   _for (i+1, n-1, 1)
     {
	k = ();
	m = string_match (buf[k], pat, 1);
	if (m != 0 or k == n-1)
	  return buf[[i:k-1]];
     }

   return NULL;
}

%}}}

private define find_entry (buf, model_name) %{{{
{
   variable pat = sprintf ("^3 %s", model_name);
   variable i, m, n = length (buf);

   _for (0, n-1, 1)
     {
	i = ();
	m = string_match (buf[i], pat, 1);
	if (m != 0)
	  return find_entry_end (buf, i);
     }

   return NULL;
}

%}}}

private define default_help_path () %{{{
{
   variable xanadu = getenv ("XANADU");
   if (xanadu == NULL)
     {
	vmessage ("*** XANADU environment variable not set");
	return NULL;
     }

   variable help = "/spectral/xspec/help/xspec.hlp";
   variable files = xanadu + [help,
			      "/src"+help];

   foreach (files)
     {
	variable f = ();
	if (NULL != stat_file (f))
	  return f;
     }

   vmessage ("*** File not found:  xspec.hlp");

   return NULL;
}

%}}}

private variable Xspec_Help_Path = default_help_path ();

define set_xspec_help_path (path) %{{{
{
   Xspec_Help_Path = path;
}

%}}}

define xspec_help () %{{{
{
   if (_NARGS != 1)
     {
	_pop_n (_NARGS);
	usage ("xspec_help (\"model\")");
	return NULL;
     }

   variable model_name = ();

   if (0 == assoc_key_exists (Xspec_Names, model_name))
     {
	vmessage ("'%s' is not among the available XSPEC models", model_name);
	return NULL;
     }

   variable path = Xspec_Help_Path;

   if (path == NULL)
     {
	vmessage ("*** Couldn't find help file 'xspec.hlp'");
	vmessage ("*** If you know where it is:  set_xspec_help_path (\"path\")");
        vmessage ("*** Otherwise, see http://heasarc.gsfc.nasa.gov/lheasoft/xanadu/xspec/manual/");
	return NULL;
     }

   variable buf, s;

   buf = read_help_file (path);
   if (buf == NULL)
     {
	vmessage ("Failed reading %s", path);
	return NULL;
     }

   s = find_entry (buf, model_name);
   if (s == NULL)
     {
	vmessage ("No help for %s", model_name);
	return NULL;
     }

   % blank out the '3' in column 1 on the first line
   s[0] = strsub (s[0], 1, ' ');
   return strjoin (s, "\n");
}

%}}}

#endif

add_help_hook ("xspec_help");

public define __unsupported (name,nargs)
{
   _pop_n(nargs);
   vmessage ("*** Error:  %s is obsolete -- please use eval_fun2", name);
}

private define alias_unsupported (name)
{
   eval(sprintf ("define %s (){ return __unsupported(_function_name, _NARGS);}",
                 name));
}

array_map (Void_Type, &alias_unsupported,
           ["_wabs", "_mekal", "_raym", "_brems", "_vaped",
            "_shock", "_nei", "_sedov", "_bmc"]);

% Support local customization.  For example, one may
% wish to modify the default hard limits of some xspec functions
% by calling __set_hard_limits -- this must be done after the
% functions are defined but before they are used.
if (2 == is_defined ("xspec_config_hook"))
{
   eval ("xspec_config_hook");
}
