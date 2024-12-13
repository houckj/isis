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

require ("arrayfuns");

define __eval();  % prototype
% eval_result performs the same function as 'eval' except that it
% explicitly returns a result.  Use eval_result instead of eval
% to avoid warnings from the stack checker.
define eval_result (s)
{
   eval(sprintf ("define __eval (){return %s;}", s));
   return __eval();
}

$1 =
[
   "physconst"
   ,"_isis"
   ,"mathmisc"
   ,"cmds"
   ,"plot-cmds"
   ,"plot-extra"
   ,"pgplot_interface"
   ,"hist-cmds"
   ,"atom-cmds"
   ,"emis-cmds"
   ,"fit-cmds"
   ,"isisopt"
   ,"group"
   ,"conf_loop"
   ,"parallel_map"
   ,"model-cmds"
   ,"aped_fun"
];

private define isis_byte_compile_file (path, f) %{{{
{
   variable file = path_concat (path, f);

   if (NULL != stat_file (file))
     {
	%vmessage ("Processing %s", file);
	byte_compile_file (file, 0);
     }
   else vmessage ("%s not found", f);
}

%}}}

private define _isis_preparse () %{{{
{
   variable path = path_concat (_isis_srcdir, "share");

   foreach ($1)
     {
	variable f = ();
	isis_byte_compile_file (path, f + ".sl");
     }

   _isis->_quit (0);
}

%}}}

if (_isis->Isis_Preparse_Only)
  _isis_preparse();

use_namespace (_isis->Isis_Public_Namespace_Name);

private define load_isis_intrinsics() %{{{
{
   ERROR_BLOCK
     {
	vmessage ("*** isis initialization failed; using ISIS_SRCDIR = %s\n", _isis_srcdir);
     }

   variable ns_pub = _isis->Isis_Public_Namespace_Name;

   foreach ($1)
     {
	variable m = ();
	() = evalfile (m, ns_pub);
     }
}

%}}}
load_isis_intrinsics();
private define load_isis_intrinsics();

% 'make check' must work before 'make install'
private variable Installed =
  (NULL == stat_file (path_concat (_isis_srcdir, "configure")));

private define load_modules () %{{{
{
   ifnot (Installed)
     {
        variable cfitsio_dir = path_concat (_isis_srcdir, "modules/cfitsio/src");
        prepend_to_isis_load_path (cfitsio_dir);
        prepend_to_isis_module_path (cfitsio_dir);
     }

   require ("fits");
}
%}}}
load_modules();
private define load_modules();

private define load_xspec_model_name_file (file) %{{{
{
   variable fp = fopen (file, "r");
   if (fp == NULL)
     return NULL;

   variable names = fgetslines (fp);
   () = fclose (fp);
   if (length(names) == 0)
     return NULL;

   return array_map (String_Type, &strtrim, names);
}

%}}}

private define init_xspec_autoload () %{{{
{
   ifnot (Installed)
     {
        variable p = path_concat (_isis_srcdir, "modules/xspec/src");
        prepend_to_isis_load_path (p);
        variable d = "elfobjs";
        prepend_to_isis_module_path (path_concat (p, d));
     }

#ifdef __XSPEC_IS_STATIC__
   () = evalfile ("xspec", _isis->Isis_Public_Namespace_Name);
   return;
#endif

#ifdef __IMPORT__
   variable dir;
   if (Installed)
     dir = path_concat (_isis_srcdir, "share");
   else
     dir = path_concat (_isis_srcdir, "modules/xspec/src");

   variable version_file = path_concat (dir, "config-xspec-version");
   if (NULL == stat_file (version_file))
     return;
   variable fp = fopen (version_file, "r");
   if (fp == NULL)
     {
        throw ApplicationError, "Failed opening $version_file"$;
     }
   variable xspec_version = fgetslines (fp, 1);
   xspec_version = strtrim (xspec_version[0]);
   () = fclose (fp);

   variable name_file;
   if (xspec_version == "HAVE_XSPEC_11")
     name_file = path_concat (dir, "_names_xspec11.dat");
   else if (xspec_version == "HAVE_XSPEC_12")
     name_file = path_concat (dir, "_names_xspec12.dat");
   else
     {
        throw ApplicationError, "$version_file specifies an unsupported version of xspec"$;
     }
   variable names = load_xspec_model_name_file (name_file);

   variable other_names =
     [
      "add_atable_model", "add_etable_model", "add_mtable_model",
      "xspec_abund", "xspec_xsect", "xspec_elabund", 
      "xspec_photo", "xspec_gphoto", "xspec_phfit2",
      "xspec_ionsneqr", "xspec_xset",
      "xspec_set_cosmo", "xspec_get_cosmo"
      ];

   foreach ([names, other_names])
     {
	variable n = ();
	autoload (n, "xspec");
     }
#endif
}
%}}}

private define init_autoloads () %{{{
{
   % cfitsio module dependence
   array_map (Void_Type, &autoload,
              ["save_conf", "load_conf", "use_file_group", "regroup_file",
               "aped_bib", "aped_bib_query_string"
              ],
              "fits_module_dep");

   autoload ("sldb", "isisdb");

   init_xspec_autoload ();
}

%}}}

init_autoloads();
private define init_autoloads();

private define init_search_paths () %{{{
{
#ifdef __ISIS_MODULE__
   prepend_to_isis_load_path (getenv ("ISIS_LOAD_PATH"));
   prepend_to_isis_module_path (getenv ("ISIS_MODULE_PATH"));
   return;
#else
   
   % we might be running 'make check' before anything has been installed:
   variable cfitsio_module_build_path, cfitsio_module_install_path;
   cfitsio_module_build_path = path_concat (_isis_srcdir, "modules/cfitsio/src/cfitsio-module.so");
   cfitsio_module_install_path = path_concat (_isis_srcdir, "lib/modules/cfitsio-module.so");
   if (NULL == stat_file (cfitsio_module_install_path))
     {
        prepend_to_isis_module_path (_isis_srcdir + "/modules/cfitsio/src");
     }

   variable opt = path_concat (_isis_srcdir, "opt");
   if (NULL != stat_file (opt))
     {
        prepend_to_isis_load_path (path_concat (opt, "share/slsh"));
        prepend_to_isis_load_path (path_concat (opt, "share/slsh/local-packages"));
        prepend_to_isis_module_path (path_concat (opt, "lib/slang/v2/modules"));
     }

   % isis-specific modules are here:
   if (NULL != stat_file (_isis_install_prefix_sans_subdir))
     {
        prepend_to_isis_load_path (path_concat (_isis_install_prefix_sans_subdir, "share/isis"));
        prepend_to_isis_module_path (path_concat (_isis_install_prefix_sans_subdir, "lib/isis/modules"));
     }

   variable ascds_path_file = path_concat (_isis_srcdir, "ascds_install_path");
   if (NULL != stat_file (ascds_path_file))
     {
        variable fp = fopen (ascds_path_file, "r");
        if (NULL == fp)
          {
             verror ("Failed opening %s for reading", ascds_path_file);
             return;
          }
        variable ascds_path = strtrim (fgetslines (fp)[0], "\n\s");
        () = fclose (fp);
        prepend_to_isis_load_path (path_concat (ascds_path, "ots/share/slsh"));
        prepend_to_isis_module_path (path_concat (ascds_path, "lib/slang/v2/modules"));        
     }   

   prepend_to_isis_load_path (getenv ("SLANG_LOAD_PATH"));
   prepend_to_isis_load_path (getenv ("ISIS_LOAD_PATH"));
   prepend_to_isis_load_path (".");

   prepend_to_isis_module_path (getenv ("SLANG_MODULE_PATH"));
   prepend_to_isis_module_path (getenv ("ISIS_MODULE_PATH"));
#endif
}

%}}}
init_search_paths ();
private define init_search_paths ();

require ("glob");
  
private define load_site_config() %{{{
{
   if (NULL != getenv ("ISIS_SKIP_SITE_CONFIG"))
     return;

   variable path = path_concat (_isis_srcdir, "etc/local.sl");
   if (NULL == stat_file(path))
     return;

   ERROR_BLOCK
     {
        vmessage ("*** failed loading %s", path);
     }

   () = evalfile (path, current_namespace());
}

%}}}
load_site_config();
private define load_site_config ();

