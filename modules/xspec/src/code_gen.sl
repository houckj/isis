#! /usr/bin/env slsh
% -*- mode: SLang; mode: fold -*-

%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2017 Massachusetts Institute of Technology
%
%     This software was developed by the MIT Center for Space Research under
%     contract SV1-61010 from the Smithsonian Institution.
%
%          Author:  John Houck <houck@space.mit.edu>
%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%

$1 = path_concat (path_dirname (__FILE__), "parse_model_dat.sl");
() = evalfile ($1);

% These file names are specified
% incorrectly in model.dat:
private variable Filename_Alias = Assoc_Type[];
Filename_Alias["xneqs"] = "xsneqs";
Filename_Alias["xshock"] = "xsshock";
Filename_Alias["xssedov"] = "xsedov";
Filename_Alias["C_zpowerLaw"] = "powerLaw";

typedef struct
{
     model_name, routine_name, model_type,
     interface, num_pars, par_info
}
Model_Fun_Type;

typedef struct
{
   name, col, type, print
}
Field_Type;

#ifnexists any
define any (a)
{
   return (0 != length(where(a)));
}
#endif

private define set_field (f, col, name, type, print) %{{{
{
   f.col = col;
   f.name = name;
   f.type = type;
   f.print = print;
}

%}}}

private define name_split (name) %{{{
{
   variable true_name = name;
   variable ext = ".f";

   if (name[1] == '_')
     {
        switch (name[0])
          {
           case 'C':
             ext = ".cxx";
             true_name = name[[2:]];
          }
          {
           case 'c':
             ext = ".c";
             true_name = name[[2:]];
          }
          {
           case 'f':
             ext = ".f";
             true_name = name[[2:]];
          }
          {
           case 'F':
             ext = ".f";
             true_name = name[[2:]];
          }
     }

   return (name[0], true_name, ext);
}

%}}}

private define emit_name_string (fp, m) %{{{
{
   % Newer functions have the 'photer' argument
   variable photer;
   if (m.interface != NULL)
     photer = (0 != is_substr (strlow(m.interface), "photer"));
   else
     photer = 1;

   variable interface_type = "f";
   if (m.routine_name[1] == '_')
     interface_type = m.routine_name[[0]];

   !if (photer)
     interface_type += "n";

   variable fields;
   fields = [Supported_Model_Types[m.model_type],
             interface_type];

   variable interface_name = strjoin(fields, "_");

   variable symbol_name;
   (,symbol_name,) = name_split (m.routine_name);

   if (m.routine_name[[0:1]] == "C_")
     symbol_name = m.routine_name;
   else if (m.routine_name[[0:1]] == "c_")
     symbol_name = m.routine_name[[2:]];
   else
     {
        variable s = symbol_name;
        symbol_name = sprintf ("FC_FUNC(%s,%s)", s, strup(s));
     }

   () = fprintf (fp.names, "%s\n", m.model_name);

   () = fprintf (fp.struct_fields, "{\"%s\", (fptr_type *) %s, \"_xspec_%s_hook\", NULL, 0},\n",
                 m.model_name, symbol_name, interface_name);

   () = fprintf (fp.externs, "extern Fcn_%s_Type %s;\n",
                 interface_type, symbol_name);
}

%}}}

private define scan_buf (b, fp) %{{{
{
   variable buf = b.buf;

   variable n_add, n_mul, n_con;
   variable e, m, s = 0;

   variable blanks = where (buf == "");

   n_add = 0;
   n_mul = 0;
   n_con = 0;

   foreach (blanks)
     {
	e = ();
	m = parse_function (buf[[s:e-1]]);
	s = e+1;

	if (m == NULL)
	  continue;

	if (m.num_pars == 0)
	  continue;

	switch (m.model_type)
	  {
	   case "add":
	     n_add++;
	  }
	  {
	   case "mul":
	     n_mul++;
	  }
	  {
	   case "con":
	     n_con++;
	  }
	  {
	     % default
	     continue;
	  }

	emit_name_string (fp, m);
     }

   vmessage ("Got %d models (%d additive, %d multiplicative, %d convolution)",
	     n_add + n_mul + n_con,
	     n_add, n_mul, n_con);
}

%}}}

private define do_scan (model_dat, fp) %{{{
{
   variable buf = load_buf (model_dat);
   scan_buf (buf, fp);
}

%}}}

private define open_output_files (xspec_version) %{{{
{
   variable fp = struct
     {
        names, struct_fields, externs
     };

   fp.names = fopen ("_names_${xspec_version}.dat"$, "w");

   fp.struct_fields = fopen ("_model_table_${xspec_version}.inc"$, "w");
   () = fprintf (fp.struct_fields, "/* -- machine generated:  do not edit -- */\n");

   fp.externs = fopen ("_model_externs_${xspec_version}.inc"$, "w");
   () = fprintf (fp.externs, "/* -- machine generated:  do not edit -- */\n");

   return fp;
}

%}}}

private define close_output_files (fp) %{{{
{
   () = fclose(fp.struct_fields);
   () = fclose(fp.externs);
}

%}}}

private define make_interface (model_dat_files, xspec_version) %{{{
{
   variable fp = open_output_files (xspec_version);
   array_map (Void_Type, &do_scan, model_dat_files, fp);
   close_output_files (fp);
}

%}}}

private define xspec_model_dat (dir, xspec_version) %{{{
{
   variable model_dat = find_model_dat_file (dir, xspec_version);
   if (model_dat == NULL)
     {
        vmessage ("*** Could not find model.dat file for %s", dir);
        exit(1);
     }

   variable local_dir = getenv ("LMODDIR");
   if (local_dir != NULL)
     {
	variable lmodel_dat = path_concat (local_dir, "lmodel.dat");

	model_dat = [model_dat, lmodel_dat];
     }

   return model_dat;
}

%}}}

define slsh_main ()
{
   variable msg = "Usage:  code_gen.sl <xspec-version> [<headas-source-dir>]";
   variable root_dir, xspec_version;

   switch (__argc)
     {
      case 3:
        xspec_version = eval(__argv[1]);
        root_dir = __argv[2];
     }
     {
      case 2:
        xspec_version = eval(__argv[1]);
        root_dir = getenv ("HEADAS");
     }
     {
        % default
        vmessage (msg);
        return;
     }

   ifnot (xspec_version == 11 or xspec_version == 12)
     {
        verror (" xspec version must be either 11 or 12");
        return;
     }

   if (root_dir == NULL)
     {
        vmessage (msg);
        return;
     }
   else if (NULL == stat_file (root_dir))
     {
        vmessage ("Cannot access %s", root_dir);
        return;
     }

   vmessage ("Searching for XSPEC model.dat in: %s", root_dir);
   make_interface (xspec_model_dat (root_dir, xspec_version),
                   sprintf ("xspec%d", xspec_version));
}

