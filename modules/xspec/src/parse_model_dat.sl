% -*- mode: SLang; mode: fold -*-
%
%   This file is part of ISIS, the Interactive Spectral Interpretation System
%   Copyright (C) 1998-2009 Massachusetts Institute of Technology
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

%  These routines are used to parse the XSPEC model.dat file.
%  They are shared by code_gen.sl and xspec.sl.
%

variable Model_Dat;     % path to model.dat file
variable Model_Srcdir;  % path to source code directory

variable Supported_Model_Types = Assoc_Type[];
Supported_Model_Types["add"] = "add";
Supported_Model_Types["mul"] = "mul";
Supported_Model_Types["con"] = "con";

define set_model_srcdir (s) %{{{
{
   Model_Srcdir = s;
}

%}}}

define get_model_srcdir () %{{{
{
   return Model_Srcdir;
}

%}}}

private define filter (s) %{{{
{
   % The model.dat file contains some quoted strings
   % that complicate parsing.
   % Since I don't care about those, I just replace
   % quoted text with a series of X's

   variable m, pos, len;
   %m = string_match (s, "\"[ \ta-zA-Z0-9+/-*\.]*\"", 1);
   m = string_match (s, "\"[^\"]*\"", 1);
   if (m == 0)
     return s;
   (pos, len) = string_match_nth (0);

   variable i;
   for (i = 1; i <= len; i++)
     {
	s = strsub (s, pos + i, 'X');
     }

   return s;
}

%}}}

define load_buf (model_dat, model_srcdir) %{{{
{
   set_model_srcdir (model_srcdir);

   variable fp = fopen (model_dat, "r");
   if (fp == NULL)
     {
	vmessage ("Failed to open %s for reading", model_dat);
	return NULL;
     }

   variable buf = fgetslines (fp);
   () = fclose(fp);

   buf = array_map (String_Type, &filter, buf);

   variable i, lines = [0:length(buf)-1];

   foreach (lines)
     {
	i = ();
	buf[i] = strtrim_end (buf[i]);
     }

   % append a blank line to simplify later parsing
   buf = [buf, ""];

   variable s = struct
     {
	buf, model_dat, model_srcdir
     };

   s.buf = buf;
   s.model_dat = model_dat;
   s.model_srcdir = model_srcdir;

   return s;
}

%}}}

private define delete_substrings (a, b) %{{{
{
   (a, ) = strreplace (a, b, "", strlen (a));
   return a;
}

%}}}

private define parse_parameter_line (line) %{{{
{
   variable t = strtok (line);

   switch (length(t))
     {
      case 8:
        return t;
     }
     {
      case 9:
        % sometimes there's a space in the parameter name.. (grrr)
        return [sprintf ("%s_%s", t[0],t[1]),
                t[[2:]]];
     }
     {
        % default:
        vmessage ("*** Failed parsing line: '%s'", line);
        exit(1);
     }
}

%}}}

private define xspec12_filter (buf) %{{{
{
   variable i,
     n = length(buf),
     drop = Int_Type[n];

   _for i (0, n-1, 1)
     {
        variable ts = strtrim_beg(buf[i]);

        if (length(ts) == 0)
          {
             drop[i] = 1;
             continue;
          }

        variable tok, num_tok;

        variable flag_char = NULL;
        if ((ts[0] == '$') || (ts[0] == '*'))
          flag_char = ts[[0]];

        if (flag_char != NULL)
          {
             tok = strtok (ts[[1:]]);
             num_tok = length(tok);

             if (num_tok >= 8)
               {
                  (buf[i], ) = strreplace (buf[i], flag_char, "", 1);
               }
             else
               {
                  variable name, value;
                  name = tok[0];

                  if (num_tok > 1 && (0 == is_substr(tok[1], "X")))
                    value = tok[1];
                  else if (num_tok > 2 && (0 == is_substr(tok[2], "X")))
                    value = tok[2];
                  else
                    value = "0";

                  buf[i] = "$name  XXX  $value  0   0  1    10   -1"$;
               }
          }
     }

   return buf[where (drop == 0)];
}

%}}}

define parse_function (buf) %{{{
{
   if (0 == length(buf))
     return NULL;

   buf = xspec12_filter (buf);

   variable hdr = strtok (buf[0]);

   % header line contains:
   % model name, num_pars, [], [], subroutine name, model_type, [], [init_string]
   variable m = struct
     {
          model_name, routine_name, model_type,
	  num_pars, par_info, init_string,
	  exec_symbol, exec_symbol_hook,
          interface
     };

   m.model_name = hdr[0];
   m.routine_name = hdr[4];
   m.model_type = strlow(hdr[5]);

   if (m.model_type == "acn")
     return NULL;

   % this should follow the test for "acn" models.
   if (length(hdr) > 7)
     {
        m.init_string = strjoin (hdr[[7:]], " ");
     }
   else m.init_string = NULL;

   variable num_pars;
   if (1 != sscanf (hdr[1], "%d", &num_pars))
     return NULL;

   if (length(buf) != num_pars+1)
     return NULL;

   m.num_pars = num_pars;

   if (num_pars == 0)
     {
	m.par_info = NULL;
	return m;
     }

   variable line, pars, i = 0;

   if (0 != strcmp (m.model_type, "add"))
	pars = String_Type[8, num_pars];
   else
     {
	pars = String_Type[8, num_pars+1];
	% add the norm using this line entry
	line = "norm  XXX  1. 0.0 0.0  1.0e10  1.e38  0.01";
	pars[*,i] = parse_parameter_line (line);
	pars[1,i] = delete_substrings (pars[1,i], "X");
	i++;
	m.num_pars++;
     }

   foreach (buf[[1:length(buf)-1]])
     {
	line = ();
	pars[*, i] = parse_parameter_line (line);
	pars[1,i] = delete_substrings (pars[1,i], "X");
	i++;
     }

   m.par_info = pars;

   return m;
}

%}}}

private define handle_xspec_object_not_found (object, xspec_version) %{{{
{
   variable s = "*** Error:  xspec $object not found."$;

#ifexists Xspec_Compiled_Headas_Path
   variable headas = getenv ("HEADAS");
   if (Xspec_Compiled_Headas_Path != headas)
     {
        s += "  This module was compiled for xspec-${xspec_version}\n"$
          +  "    with HEADAS=$Xspec_Compiled_Headas_Path\n"$
          +  "     Now HEADAS=$HEADAS"$;
     }
   else
     {
        s += "\n  Using HEADAS=$headas"$;
     }
#endif

   throw ApplicationError, s;
}

%}}}

define find_model_srcdir (dir, xspec_version) %{{{
{
   variable model_srcdirs;
   switch (xspec_version)
     {
      case 11:
        model_srcdirs = ["src/spectral/xspec/src/functions/",
                         "ftools/spectral/xspec/src/functions/"];
     }
     {
      case 12:
        model_srcdirs = ["Xspec/src/XSFunctions/"];
     }
     {
        % default:
        throw ApplicationError, "xspec version ${xspec_version} is unsupported"$;
     }

   variable i, n = length(model_srcdirs);

   for (i = 0 ; i < n; i++)
     {
        variable path = path_concat (dir, model_srcdirs[i]);
        if (NULL != stat_file (path))
          return path;
     }

   handle_xspec_object_not_found ("source code", xspec_version);

   return NULL;
}

%}}}

define find_model_dat_file (dir, xspec_version) %{{{
{
   variable model_dat_paths;
   switch (xspec_version)
     {
      case 11:
        model_dat_paths = ["src/spectral/xspec/src/MODEL_COMP/model.dat"
                           , "spectral/xspec/src/local_mod/model.dat"];
     }
     {
      case 12:
        model_dat_paths = ["Xspec/src/manager/model.dat"
                           , "spectral/manager/model.dat"];
     }
     {
        % default:
        throw ApplicationError, "xspec version ${xspec_version} is unsupported"$;
     }

   foreach (model_dat_paths)
     {
        variable x = ();
        variable path = path_concat (dir, x);
        if (NULL != stat_file (path))
          return path;
        path = path_concat (dir + "/..", x);
        if (NULL != stat_file (path))
          return path;
     }

   handle_xspec_object_not_found ("model.dat", xspec_version);

   return NULL;
}

%}}}

provide ("parse_model_dat");
