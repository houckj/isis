% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2018 Massachusetts Institute of Technology
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

% $Id: cmds.sl,v 1.28 2004/02/09 11:14:14 houck Exp $

%_debug_info = 1;           % compile line number information into functions

define isis_get_pager ();

#ifnexists source
define source ()
{
   variable append_semicolon = Isis_Append_Semicolon;
   try
     {
        Isis_Append_Semicolon = 0;
        _isis->_source ();
     }
   finally
     {
        Isis_Append_Semicolon = append_semicolon;
     }
}
#endif

define stop_log () %{{{
{
   return _isis->_stop_log ();
}

%}}}

define start_log () %{{{
{
   variable args = __pop_args (_NARGS);
   return _isis->_start_log (__push_args (args));
}

%}}}

private variable Isis_Namespace = _isis->Isis_Public_Namespace_Name;

private define make_apropos_list (list, ns) %{{{
{
   variable n = length(list);
   variable s = {};

   list = list[array_sort(list, &strcmp)];

   variable i = 0;
   loop (n / 3)
     {
        list_append (s, sprintf ("%-26s %-26s %s\n", list[i], list[i+1], list[i+2]));
	i += 3;
     }
   n = n mod 3;
   loop (n)
     {
        list_append (s, sprintf ("%-26s ", list[i]));
	i++;
     }

   if (i mod 3 != 0)
     list_append (s, "\n");

   return s;
}

%}}}

define apropos () %{{{
{
   !if (_NARGS)
     usage ("apropos (\"keyword\")");

   variable what = ();
   what = str_delete_chars (what, ";\"");

   % default to case sensitive
   if (0 == is_substr (what, "\\"))
     what = "\\C" + what;

   variable namespaces = [Isis_Namespace, current_namespace()];

   variable masks = Assoc_Type[];
   masks["function"] = 1|2;
   masks["variable"] = 4|8;

   variable ns, type, list, len = 0,
     s = {};

   foreach type (["function", "variable"])
     {
        foreach ns (namespaces)
          {
             list = _apropos (ns, what, masks[type]);
             variable n = length(list);
             if (n == 0)
               continue;
             len += n;
             list_append (s, sprintf ("\nFound %d %s match%s in namespace %S:\n",
                                      n, type,
                                      n>1 ? "es" : "",
                                      ns));
             list_append (s, make_apropos_list (list, ns));
          }
     }

   if (len == 0)
     {
        vmessage ("No matches");
        return;
     }

   variable fp = stdout;
   if (sum(len) > 24*3)  % 3 per line, 24 lines
     {
        variable p = isis_get_pager ();
        if (p != NULL)
          fp = popen (p, "w");
     }

   variable ss;
   foreach ss (s)
     {
        if (typeof(ss) != String_Type)
          ss = strjoin (list_to_array (ss), "");
        () = fprintf (fp, "%s", ss);
     }

   if (fp != stdout)
     () = pclose (fp);
}

%}}}

private define send_to_pager (doc) %{{{
{
   variable p, fp = stdout;

   p = isis_get_pager ();
   if (p != NULL)
     fp = popen (p, "w");

   () = fputs (doc, fp);

   if (fp != stdout)
     () = pclose (fp);
}

%}}}

private define show_doc (help_file, topic) %{{{
{
   variable doc;

   if (help_file != NULL)
     doc = get_doc_string_from_file (help_file, topic);
   else
     {
        doc = get_doc_string_from_file (topic);
     }

   if (doc != NULL)
     {
	send_to_pager (doc);
	return 0;
     }

   if ((help_file != NULL)
       and (Isis_Verbose >= _isis->_WARN))
     {
        () = printf ("Searched: %s\n", help_file);
     }

   return -1;
}

%}}}

private variable Help_Hooks = Assoc_Type[];

private define try_help_hooks (topic) %{{{
{
   foreach (Help_Hooks) using ("values")
     {
	variable ref = ();
        variable s = (@ref) (topic);
        if (typeof(s) == String_Type)
          {
             send_to_pager (s);
             return 0;
          }
        else if (s == 0)
          return 0;
     }

   return -1;
}

%}}}

define add_help_hook () %{{{
{
   variable msg = "add_help_hook (name)";
   variable name, ref;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   name = ();

   ref = __get_reference (name);
   if (ref == NULL)
     {
	verror ("*** function '%s' not found", name);
	return;
     }

   Help_Hooks[name] = ref;
}

%}}}

require ("slshhelp");
add_help_hook ("slsh_get_doc_string");

define add_help_file () %{{{
{
   variable msg = "add_help_file (file)";
   variable f;

   if (_isis->get_varargs (&f, _NARGS, 1, msg))
     return;

   !if (path_is_absolute (f))
     f = path_concat (getcwd(), f);

   add_doc_file(f);
}

%}}}

private define init_help_file_list () %{{{
{
   variable files = ["maplib.hlp", "cfitsio.hlp", "help.txt", "local_help.txt"];
   variable doc_dir = path_concat (_isis_srcdir, "doc");

   if (__is_initialized (&_slang_doc_dir))
     add_help_file (path_concat (_slang_doc_dir, "slangfun.txt"));

   foreach (files)
     {
	variable f = ();
	add_help_file (path_concat (doc_dir, f));
     }
}

%}}}

init_help_file_list ();

define help () %{{{
{
   variable help_file, topic;
   variable msg =
`Several sources of help are available:
1. To get help on a particular topic, try one of the following:
        help topic
       .help topic
        help ("topic");
2. To find symbols containing a particular substring, try one of:
        apropos substring
       .apropos substring
        apropos ("substring");
3. Many functions provide a brief usage message when called with
   no parameters.
4. For help on isis command-line options, try:
        > isis --help
5. ISIS documentation is available at:
         http://space.mit.edu/cxc/isis/docs.html
   S-Lang documentation is available at:
         http://www.jedsoft.org/slang/docs.html
6. The ISIS FAQ is available at:
         http://space.mit.edu/cxc/isis/faq.html
7. If none of these sources provides an answer to your question,
   feel free to ask isis-users@space.mit.edu`;

   if (_NARGS != 1)
     {
        message (msg);
        return;
     }

   topic = ();
   topic = str_delete_chars (topic, ";");

   if (topic == "isis")
     {
        message (msg);
        return;
     }
   else if (topic == "changes")
     {
        variable change_log = path_concat (_isis_srcdir, "Changes.txt");
        variable fp = fopen (change_log, "r");
        if (fp == NULL)
          throw ApplicationError, "*** Error:  cannot open $change_log for reading"$;
        variable s = fgetslines (fp);
        () = fclose (fp);
        send_to_pager (strjoin (s, ""));
        return;
     }

   if (0 == show_doc(NULL, topic))
     return;

   if (0 == try_help_hooks (topic))
     return;

   ()=printf ("No help on '%s' is available\n", topic);

   apropos (topic);
}

%}}}

define isis_set_pager () %{{{
{
   variable msg = "isis_set_pager (\"pager\")";
   variable p;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   p = ();

   if (p == NULL) p = "";
   _isis->_isis_set_pager (p);
}

%}}}

define isis_get_pager () %{{{
{
   return _isis->_isis_get_pager();
}

%}}}

define alias () %{{{
{
   variable msg = "alias (from-cmd, to-cmd)";
   variable from, to;
   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;
   (from, to) = ();

   variable forbidden_chars = "[)(\]\[]+"R;
   if (string_matches (from, forbidden_chars) != NULL)
     {
        throw UsageError, "*** unsupported syntax:  '$from'"$;
     }
   if (string_matches (to, forbidden_chars) != NULL)
     {
        throw UsageError, "*** unsupported syntax:  '$to'"$;
     }

   eval(sprintf("define %s() { variable args = __pop_args (_NARGS); %s(__push_args(args) ;; __qualifiers);}", to, from),
	Isis_Namespace);
}

%}}}

define reset () %{{{
{
#iftrue
   _pop_n (_NARGS);
   vmessage ("***WARNING: 'reset' has been disabled and will soon be removed");
   return;
#else
   variable force = 0;

   switch (_NARGS)
     {
      case 0:
	variable buf;
	() = fprintf (stderr, "Confirm reset (y/n):  ");
	() = fgets (&buf, stdin);
	if ((buf[0] == 'y') or (buf[0] == 'Y'))
	  force = 1;
	else message ("ok, no reset");
     }
     {
      case 1:
	force = ();
     }
     {
	% default:
	_pop_n (_NARGS);
	message ("reset ([force])");
	return;
     }

   if (force)
     _isis->_reset;
#endif
}

%}}}

define array_struct_field () %{{{
{
   variable msg = "a[] = array_struct_field (Struct_Type[], \"field_name\")";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable s, f;
   (s, f) = ();

   if (typeof (s) == Struct_Type)
     return get_struct_field (s, f);

   variable type = typeof (get_struct_field (s[0], f));
   return array_map (type, &get_struct_field, s, f);
}

%}}}

private define print_struct () %{{{
{
   variable msg = "print_struct (s)";
   variable name, value, s;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   s = ();

   foreach (get_struct_field_names (s))
     {
	name = ();
	value = get_struct_field (s, name);
	vmessage ("    %s = %s", name, string(value));
     }
}

%}}}

require ("print");
#ifnexists print
% borrowed from John Davis's jdl code
define print () %{{{
{
   variable msg = "print (struct or array)";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable i;
   variable args = __pop_args (_NARGS);

   for (i = 0; i < _NARGS; i++)
     {
	variable arg = args[i].value;

	if (is_struct_type (arg))
	  arg = typecast (arg, Struct_Type);

	switch (typeof (arg))
	  {
	   case Struct_Type:
	     print_struct (arg);
	  }
	  {
	   case Array_Type:
	     _isis->print_array (arg);
	  }
	  {
	     ()=printf ("%S\n", args[i].value);
	  }
     }
   () = fflush (stdout);
}

%}}}
#endif

define who () %{{{
{
   variable ns = current_namespace();
   variable a_var, a_fun, name, obj;
   variable pat, str;

   if (_NARGS == 0)
     "";
   pat = ();

   a_var = _apropos (ns, pat, 8);
   a_fun = _apropos (ns, pat, 2);

   variable fp = stdout;

   variable num = length(a_var) + length(a_fun);
   if (num > 24)
     {
        variable p = isis_get_pager ();
        if (p != NULL)
          fp = popen (p, "w");
     }

   foreach (a_var[array_sort(a_var)])
     {
	name = ();
	obj = __get_reference (name);
	if (__is_initialized (obj))
	  str = string (@obj);
	else
	  str = "*** Not Initialized ***";

	() = fprintf (fp, "%s: %s\n", name, str);
     }
   foreach (a_fun[array_sort(a_fun)])
     {
	name = ();
	() = fprintf (fp, "%s()\n", name);
     }

   if (fp != stdout)
     () = pclose (fp);
}

%}}}

define delete () %{{{
{
   variable pat;

   if (_NARGS == 0)
     "";
   pat = ();

   foreach (_apropos (current_namespace(), pat, 8))
     __uninitialize (__get_reference ());
}

%}}}

define save_input () %{{{
{
   variable file = NULL;

   if (_NARGS)
     file = ();

   _isis->do_save_input (file);
}

%}}}

define find_library_name () %{{{
{
   variable msg = "full_path = find_library_name (\"libname\")";

   if (_isis->chk_num_args(_NARGS, 1, msg))
     return;

   variable libname = ();
   return _isis->_find_file_in_path (get_isis_module_path(), libname);
}

%}}}

define find_script_name () %{{{
{
   variable msg = "full_path = find_script_name (\"name\")";

   if (_isis->chk_num_args(_NARGS, 1, msg))
     return;

   variable name = ();
   return _isis->_find_file_in_path (get_isis_load_path(), name);
}

%}}}

#ifnexists struct_field_exists
define struct_field_exists (s, name) %{{{
{
   variable fields = get_struct_field_names (s);
   return length(where(fields == name));
}

%}}}
#endif

define _A ();  % declare for recursive call
define _A () %{{{
{
   variable n = _NARGS;
   loop (n)
     {
	variable e = ();
	if (typeof (e) == Struct_Type)
	  {
	     variable cpy = @e;
	     (cpy.bin_lo, cpy.bin_hi) = _A(e.bin_lo, e.bin_hi);

	     if (struct_field_exists (e, "value"))
	       cpy.value = reverse (e.value);

	     if (struct_field_exists (e, "err"))
	       cpy.err = reverse (e.err);

	     cpy;
	  }
	else
	  reverse (Const_keV_A / e);
        _stk_roll (n);
     }
   _stk_reverse (n);
}

%}}}

private variable Units = Assoc_Type[];

define unit_add ()
{
   variable msg = "unit_add (String_Type name, Int_Type is_energy, scale)";
   variable name, is_energy, scale;

   if (_NARGS != 3)
     usage (msg);
   (name, is_energy, scale) = ();

   ifnot (typeof(name) == String_Type
          && typeof (is_energy) == Integer_Type
          && __is_numeric (scale))
     {
        usage (msg);
     }

   variable u = struct
     {
        name = name,
        is_energy = is_energy,
        scale = scale
        % Multiplying by 'scale' converts energy units to keV
        % and wavelength units to Angstrom.
     };

   Units[strlow(name)] = u;
}

unit_add ("Angstrom",  0, 1.0);
unit_add ("A",  0, 1.0);
unit_add ("nm", 0, 1.e1);
unit_add ("um", 0, 1.e4);
unit_add ("mm", 0, 1.e7);
unit_add ("cm", 0, 1.e8);
unit_add ("m",  0, 1.e10);

unit_add ("eV",  1, 1.e-3);
unit_add ("keV", 1, 1.0);
unit_add ("MeV", 1, 1.e3);
unit_add ("GeV", 1, 1.e6);
unit_add ("TeV", 1, 1.e9);

private variable hz_to_kev = Const_h / (Const_eV * 1.e3);
unit_add ("Hz",  1, hz_to_kev);
unit_add ("kHz", 1, hz_to_kev * 1.e3);
unit_add ("MHz", 1, hz_to_kev * 1.e6);
unit_add ("GHz", 1, hz_to_kev * 1.e9);

define unit_exists ()
{
   if (_NARGS != 1) usage ("unit_exists (unit)");
   variable unit = ();
   return assoc_key_exists (Units, strlow(unit));
}

private variable Default_Unit;

private define unit_default_init (u)
{
   ifnot (unit_exists (u))
     {
        vmessage ("Unsupported unit: %s", u);
        return;
     }
   Default_Unit = list_new();
   list_insert (Default_Unit, u);
}
unit_default_init ("Angstrom");

define unit_push ()
{
   if (_NARGS != 1)
     usage ("unit_push (unit)");

   variable unit = ();
   if (unit_exists (unit))
     list_insert (Default_Unit, unit);
   else vmessage ("Unsupported unit: %s", unit);
}

define unit_pop ()
{
   if (length(Default_Unit) == 1)
     return Default_Unit[0];
   return list_pop (Default_Unit);
}

define unit_default ()
{
   switch (_NARGS)
     {case 0: return Default_Unit[0];}
     {case 1: unit_default_init ();}
     {
        usage ("unit = unit_default()  OR  unit_default (unit)");
     }
}

define unit_info ()
{
   variable s = unit_default ();
   if (_NARGS == 1)
     s = ();
   else if (_NARGS > 1)
     usage ("Struct_Type = unit_info ([unit])");

   if (unit_exists (s))
     return Units[strlow(s)];
   else vmessage ("Unsupported unit: %s", s);
}

define unit_is_energy ()
{
   variable s = unit_default ();
   if (_NARGS == 1)
     s = ();
   else if (_NARGS > 1)
     usage ("Int_Type = unit_is_energy ([unit])");

   return unit_info(s).is_energy;
}

define to_angstrom ()
{
   variable msg = "x_angstrom = to_angstrom (x)";
   if (_NARGS != 1) usage(msg);

   variable x = ();

   variable units, u, xa;
   units = qualifier ("unit", unit_default());

   u = unit_info (units);
   xa = x * u.scale;
   if (u.is_energy) xa = _A(xa);

   return xa;
}

define readcol () %{{{
{
   variable msg = "(a, b, ..) = readcol (file, [c1, c2, ..]);";
   variable file, cols;

   if (_NARGS == 0)
     {
	_pop_n (_NARGS);
	usage (msg);
	return;
     }
   else if (_NARGS == 1)
     {
	file = ();
	cols = 1;
     }
   else
     {
	cols = _isis->pop_list (_NARGS-1, msg);
	file = ();
     }

   % make sure returned columns are returned in the
   % correct order

   variable i, n;
   i = array_sort (cols);
   n = _stkdepth();

   _isis->_readcol (cols[i], file);

   variable args = __pop_args (_stkdepth() - n);
   __push_args (args[array_sort(i)]);
}

%}}}

define writecol () %{{{
{
   if (_NARGS < 2)
     {
        _pop_n (_NARGS);
        vmessage ("%s (fp, a, b, ....);", _function_name);
        return;
     }
   variable nargs = _NARGS-1;
   variable args = __pop_args (nargs);
   variable fp = ();

   if (String_Type == typeof (fp))
     {
        variable _fp = fopen (fp, "w");
        if (_fp == NULL)
          verror ("writecol:  unable to open %s", fp);
        fp = _fp;
     }

   variable fmt = "";
   loop (nargs)
     fmt = strcat ("  %12S", fmt);
   fmt = strcat (strtrim (fmt), "\n");

   variable len = length(args[0].value);
   foreach (args)
     {
	variable x = ();
	if (len != length(x.value))
	  verror ("writecol:  arrays have unequal length");
     }

   variable status = array_map (Int_Type, &fprintf,
				fp, fmt, __push_args (args));
   if (length (where (status == -1)))
     verror ("writecol:  error writing to file");

   () = fflush (fp);
   %if (fp != stdout and fp != stderr)
   %  _isis->close_file (fp);
}

%}}}

%  Useful mnemonic variables:

variable  % Elements (Z):
  H = 1,
  He = 2,
  C  = 6,
  N  = 7,
  O  = 8,
  Ne = 10,
  Na = 11,
  Mg = 12,
  Al = 13,
  Si = 14,
  S  = 16,
  Ar = 18,
  Ca = 20,
  Fe = 26,
  Ni = 28;

%                                            Plot Colors:
variable white  = 1;   % or black...
variable red    = 2;
variable green  = 3;
variable blue   = 4;
variable ltblue = 5;
variable purple = 6;
variable yellow = 7;
variable orange = 8;
variable grey   = 15;

private define check_file (f) %{{{
{
   if (NULL != stat_file (f))
     return;
   verror ("*** Error: file not found: %s", f);
}

%}}}

define backtrace_filter () %{{{
{
   variable msg = "backtrace_filter (\"trace_file.txt\")";
   variable trace_file;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   trace_file = ();

   variable isis_binary = path_concat (_isis_srcdir, "bin/isis");
   variable btscript = path_concat (_isis_srcdir, "share/backtrace.sh");

   array_map (Void_Type, &check_file, [trace_file, isis_binary, btscript]);

   variable cmd = sprintf ("%s %s %s", btscript, trace_file, isis_binary);
   () = system (cmd);
}

%}}}

% to turn this off, just re-define motd_print() to do nothing.
define motd_print ()
{
   variable fp = fopen (path_concat (_isis_srcdir, "etc/motd.txt"), "r");
   if (fp == NULL)
     return;
   variable motd = fgetslines (fp);
   () = fclose (fp);
   motd = array_map (String_Type, &strtrim_end, motd, "\n");
   array_map (Void_Type, &message, motd);
}
