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

% $Id: _isis.sl,v 1.10 2004/02/09 11:14:14 houck Exp $

use_namespace ("_isis");

define pop_list (num, msg) %{{{
{
   if (num == 0)
     {
	usage (msg);
	return NULL;
     }

   variable a;

   if (num > 1)
     {
	a = __pop_args (num);
	a = [__push_args (a)];
     }
   else a = ();

   return a;
}

%}}}

define get_varargs (nargs, min_num_args, msg) %{{{
{
   variable max_num_args;
   variable arg_ptr;
   variable addrs, i;

   max_num_args = _NARGS - 3;

   if ((nargs < min_num_args) or (nargs > max_num_args))
     {
	_pop_n (max_num_args + nargs);
	usage (msg);
	return 1;
     }

   addrs = __pop_args (max_num_args);
   _stk_reverse (nargs);

   for (i = 0; i < nargs; i++)
     {
     % Use temp variable since @addrs[i].value = (); is not yet supported
	variable addr = addrs[i].value;
	@addr = ();
     }

   return 0;
}

%}}}

define chk_num_args (nargs, num, msg) %{{{
{
   if (nargs == num)
     return 0;

   _pop_n (nargs);
   usage (msg);
   return 1;
}

%}}}

define _sys_shell_cmd (cmd) %{{{
{
   variable status;

   status = system (cmd);
   !if (status)
     return;

   () = fprintf (stdout, "shell command returned %d\n", status);
}

%}}}

define _sys_chdir_cmd (dir) %{{{
{
   if (-1 == chdir (dir))
     verror ("chdir(%s) failed: %s", dir, errno_string (errno));
}

%}}}

define my_str_replace (a, b, c) %{{{
{
   () = str_replace (a, b, c);
   !if (_stkdepth) return a;
}

%}}}

% Support interactive mode which doesn't require the ';' EOL mark.
% Note that this no-semicolon-required mode causes surprising
% behavior in the '.source' intrinsic.  I'm taking the view
% that that's the price users must pay for the convenience of
% not typing the semicolons.
public variable Isis_Append_Semicolon = 0;
define maybe_append_semicolon (input) %{{{
{
   if (Isis_Append_Semicolon != 0 and input[-1] != ';')
     input += ";";
   return input;
}

%}}}

define isis_massage_input_hook (input) %{{{
{
   variable s, n;
   variable s0;
   variable ch;

   input = strtrim (input);
   !if (strlen (input))
     return input;

   % If cut-and-paste entered a copy of the prompt, ignore it.
   variable prompt = get_prompt();
   variable prompt_len = strlen(prompt);

   if (0 == strncmp (prompt, input, prompt_len))
     input = input[[prompt_len:]];

   ch = input[0];

   if (ch == '!')
     {
	% shell command, unlesss !if
	if (string_match (input, "^!if[ \t]*(", 1))
	  return input;

	input = strtrim (input[[1:]]);
        if (strlen(input) == 0)
          return input;

        input = str_quote_string (input, "\"", '\\');

	% cd is special.  On unix systems we do not want to execute it in
	% a shell.
	if (strtok (input)[0] == "cd")
	  {
	     return sprintf ("_isis->_sys_chdir_cmd(\"%s\");", strtrim (input[[2:]]));
	  }

	return sprintf ("_isis->_sys_shell_cmd(\"%s\");", input);
     }

   if (ch != '.')
     {
	s = strtok (input);
	s0 = s[0];
	if ((s0 != "help") and (s0 != "apropos") and (s0 != "quit"))
	  return maybe_append_semicolon (input);

	if (length(s) > 1)
	  {
	     if (s[1][0] == '(')
	       return input;
	  }
     }
   else
     {
	if (String_Type != _slang_guess_type (strtok (input, "-+*/<>&|; \t")[0]))
	  return sprintf ("eval(\" %s\");", input);

	input = input[[1:]];

	s = strtok (input);
	s0 = s[0];
     }

   if (length (s) < 2)
     {
	if (s0 == "help")
	  return "help();";

	if (s0 == "quit")
	  return "quit();";

	return input;
     }

   if (length (s) == 2)
     {
	variable s1 = s[1];

	if (s1[0] != '(')
	  {
	     if ((s0 == "apropos")
		 or (s0 == "help"))
	       {
		  if (0 == is_substr (s1, "\""))
		    return sprintf ("%s(\"%s\");", s0, s1);
		  else
		    return sprintf ("%s(%s);", s0, s1);
	       }

	     if (s0 == "cd")
	       return sprintf ("_isis->_sys_chdir_cmd(\"%s\");", s1);

	     %  Assume the next arg is a slang script.
	     if ((strlen(path_extname (s1)) == 0)
		 and (NULL != stat_file (s1 + ".sl")))
	       s1 = s1 + ".sl";

	     if (s0 == "load")
	       return sprintf ("()=evalfile(\"%s\");", s1);

	     if (s0 == "source")
	       return sprintf ("%s(\"%s\");", s0, s1);
	  }
     }

   % Anything else that does not look like a slang command will be converted
   % to a function call
   input = sprintf ("%s(%s);", s0, input[[strlen (s0):]]);

   return input;
}

%}}}

%{{{ GNU readline symbol completion hooks:
private variable gnu_rl_completions;
private variable list_index;
define gnu_rl_completion_hook (text, state) %{{{
{
   if (state == 0)
     {
        variable names = _apropos ("Global", "^" + text, 0xF);
        variable ns = current_namespace ();
        if ("Global" != ns)
          {
             names = [names, _apropos (ns, "^" + text, 0xF)];
          }
        gnu_rl_completions = names[array_sort(names)];
        list_index = 0;
     }

   if (gnu_rl_completions == NULL)
     return NULL;

   if (list_index >= length(gnu_rl_completions))
     {
        gnu_rl_completions = NULL;
        return NULL;
     }

   variable len = strlen(text);

   foreach (gnu_rl_completions[[list_index:]])
     {
        variable s = ();
        list_index += 1;
        if (0 == strncmp (s, text, len))
          return s;
     }

   return NULL;
}

%}}}

%}}}

private variable Isis_Log_File = "isis.log";
private variable Isis_Log_File_Fp = NULL;
private variable Input_Line_List = NULL;
private variable Last_Input_Line = NULL;

define warn_fprintf () %{{{
{
   variable status, fp, fmt, args;
   variable msg = "status = warn_fprintf (fp, fmt [,args]);";

   _isis->_isis_set_errno (0);

   switch (_NARGS)
     {
      case 0 or case 1:
	usage (msg);
	return 0;
     }
     {
      case 2:
	(fp, fmt) = ();
	() = fprintf (fp, fmt);
     }
     {
	% default:
	args = __pop_args (_NARGS-2);
	(fp, fmt) = ();
	() = fprintf (fp, fmt, __push_args (args));
     }

   () = fflush (fp);

   if (errno != 0)
     {
	verror ("*** Error:  write failed (%s)", errno_string(errno));
	return -1;
     }

   return 0;
}

%}}}

define close_file (fp) %{{{
{
   _isis->_isis_set_errno (0);

   ()=fflush(fp);

   if (0 != fclose (fp)
	or errno != 0)
     {
	verror ("*** Error:  closing file (%s)", errno_string(errno));
     }
}

%}}}

define _stop_log () %{{{
{
   ERROR_BLOCK
     {
	Isis_Log_File_Fp = NULL;
	_isis->_isis_set_errno (0);
	_clear_error;
	vmessage ("Stopped logging input to %S", Isis_Log_File);
	return;
     }

   if (Isis_Log_File_Fp == NULL)
     return;

   () = warn_fprintf (Isis_Log_File_Fp, "%% Log stopped %s\n", time ());
   close_file (Isis_Log_File_Fp);

   EXECUTE_ERROR_BLOCK;
}

%}}}

define open_log_file (file) %{{{
{
   variable fp = fopen (file, "a");
   if (fp == NULL)
     vmessage ("***Warning: Unable to log to %s\n", file);
   return fp;
}

%}}}

define _start_log () %{{{
{
   if (Isis_Log_File_Fp != NULL)
     _stop_log();

   if (_NARGS > 0)
     Isis_Log_File = ();

   Isis_Log_File_Fp = open_log_file (Isis_Log_File);
   if (Isis_Log_File_Fp == NULL)
     return;

   () = warn_fprintf (Isis_Log_File_Fp, "%% Log started %s\n", time ());
   vmessage ("Logging input to %s", Isis_Log_File);
}

%}}}

define do_save_input (file) %{{{
{
   if (file == NULL)
     file = Isis_Log_File;

   if (Input_Line_List == NULL)
     return;

   variable fp = open_log_file (file);
   if (fp == NULL)
     return;

   variable l = Input_Line_List;
   while (l != NULL)
     {
	() = warn_fprintf (fp, "%s\n", l.line);
	l = l.next;
     }

   close_file (fp);
}

%}}}

define log_input (buf) %{{{
{
   buf = strtrim (buf);
   !if (strlen (buf))
     return;

   variable l = struct
     {
	line, next
     };
   l.line = buf;
   if (Input_Line_List == NULL)
     {
	Input_Line_List = l;
     }
   else Last_Input_Line.next = l;
   Last_Input_Line = l;

   if (Isis_Log_File_Fp != NULL)
     {
	() = fputs (buf, Isis_Log_File_Fp);
	() = fputs ("\n", Isis_Log_File_Fp);
	() = fflush (Isis_Log_File_Fp);
     }
}

%}}}

define isis_readline_hook (buf) %{{{
{
   log_input (buf);
}

%}}}

define take_input_hook () %{{{
{
   while (_stkdepth ())
     message(string());
}

%}}}

variable Dbase = struct
{
   dir, atomic_data_filemap,
     abundance, ion_balance,
     line_emissivity, continuum_emissivity
};

define atomdb_start_hook (filemap_filename) %{{{
{
   variable fp = fopen (filemap_filename, "r");
   if (NULL == fp)
     {
        vmessage ("Failed opening file %s for reading", filemap_filename);
        return;
     }

   variable s = fgetslines (fp);
   () = fclose(fp);

   if (NULL == s)
     {
	vmessage ("Failed reading file %s", filemap_filename);
	return;
     }

   variable ss = array_map (String_Type, &strcompress, s, " \t\n");
   variable type = array_map (String_Type, &extract_element, ss, 0, ' ');

   variable e = where (type == "2");
   variable w = where (type == "3");

   variable elev_files = array_map (String_Type, &extract_element, ss[e], 3, ' ');
   variable wavelen_files = array_map (String_Type, &extract_element, ss[w], 3, ' ');

   variable dir = _isis->Dbase.dir;

   variable v, env_variables = ["$ATOMDB", "$APEC_DIR"];
   variable env = NULL;
   foreach (env_variables)
     {
	v = ();
	if (is_substr (elev_files[0], v))
	  {
	     env = v;
	     break;
	  }
     }

   % FIXME - do something better here?
   if (env == NULL)
     env = "XXX";

   variable elev_paths, wavelen_paths;
   elev_paths = array_map (String_Type, &my_str_replace, elev_files, env, dir);
   wavelen_paths = array_map (String_Type, &my_str_replace, wavelen_files, env, dir);

   return (wavelen_paths, elev_paths);
}

%}}}

define emisdb_start_hook () %{{{
{
   variable s = struct
     {
	line_emis, contin_emis, ionization, abundance, filemap
     };

   s.line_emis = ""; s.contin_emis = ""; s.ionization = "";
   s.abundance = ""; s.filemap = "";

   if (_isis->Dbase.atomic_data_filemap != NULL)
     s.filemap = path_concat (_isis->Dbase.dir, _isis->Dbase.atomic_data_filemap);

   if (_isis->Dbase.abundance != NULL)
     s.abundance = path_concat (_isis->Dbase.dir, _isis->Dbase.abundance);

   if (_isis->Dbase.ion_balance != NULL)
     s.ionization = path_concat (_isis->Dbase.dir, _isis->Dbase.ion_balance);

   if (_isis->Dbase.line_emissivity != NULL)
     s.line_emis = path_concat (_isis->Dbase.dir, _isis->Dbase.line_emissivity);

   if (_isis->Dbase.continuum_emissivity != NULL)
     s.contin_emis = path_concat (_isis->Dbase.dir, _isis->Dbase.continuum_emissivity);

   return s;
}

%}}}

define kernel_init_hook (name, arg_string) %{{{
{
   message (arg_string);
   variable pnames = strchop (arg_string, ',', 0);
   _isis->_add_slangfun_intrin (__get_reference (name + "_fit"), pnames, 0);
}

%}}}

define list_functions_hook () %{{{
{
   variable num_cols = 4;
   variable names = ();

   if (names == NULL)
     return;

   variable n = length(names);
   variable rows = n / num_cols;
   if ((n mod num_cols) != 0) rows += 1;

   variable indx = array_sort (names, &strcmp);

   variable p, i, r;

   for (r = 0; r < rows; r++)
     {
	p = r;
	for (i = 0; i < num_cols and p < n; i++)
	  {
	     () = fprintf (stdout, "%-15s ", names[indx[p]]);
	     p += rows;
	  }
	() = fputs ("\n", stdout);
     }

   %() = fputs ("\n", stdout);
}

%}}}

variable Model_Component_Names = NULL;

define name_mode_eval_hook (name) %{{{
{
   variable s = struct {next, name};

   % When starting a new list, Model_Component_Names=NULL
   s.name = name;
   s.next = Model_Component_Names;
   Model_Component_Names = @s;
}

%}}}

define error_if_fit_in_progress (s) %{{{
{
   if (Isis_Fit_In_Progress == 0)
     return;
   verror ("*** Warning: %s may not be used while a fit is in progress.", s);
}

%}}}

define put_string (dest, fname, s) %{{{
{
   if (s == NULL)
     return;

   if (dest == NULL)
     dest = stdout;

   switch (typeof(dest))
     {
      case Ref_Type:
        @dest = s;
     }
     {
      case File_Type:
        if (-1 == fprintf (dest, "%s\n", s))
          throw IsisError, "$fname: Failed writing list"$;
     }
     {
      case String_Type:
        if (is_substr (dest, "\n"))
          throw IsisError, "$fname: Error output file name contains a newline character"$;
        variable fp = fopen (dest, "w");
        if (fp == NULL)
          throw IsisError, "$fname: Failed opening file $dest"$;
        if (-1 == fprintf (fp, "%s\n", s))
          throw IsisError, "$fname: Error writing to file $dest"$;
        if (fclose (fp))
          throw IsisError, sprintf("$fname: Error closing file $dest (%s) "$, errno_string(errno));
     }
     {
        % default:
        throw IsisError, "$fname: unsupported destination type"$;
     }
}

%}}}

% Supports spectral model caching (cache_fun)
variable Model_Cache = Assoc_Type[];

private define array_to_struct (fields) %{{{
{
   eval (sprintf ("define __atos__(){return struct {%s};}",
                  strjoin (fields, ",")));
   return eval ("__atos__()");
}
%}}}

% Supports passing options to slang optimizers
define options_to_struct (options) %{{{
{
   % options is an array of strings of the form 'name=value'
   variable n = length(options);

   if (n == 0)
     return 0;

   variable i, names = String_Type[n];
   _for i (0, n-1, 1)
     {
        names[i] = strtok(options[i], "=")[0];
     }
   variable s = array_to_struct (names);

   _for i (0, n-1, 1)
     {
        variable t = strtok (options[i], "=");
        if (length(t) < 2)
          continue;
        variable v = t[1];
        switch (__is_datatype_numeric(_slang_guess_type (v)))
          { case 0:             }
          { case 1: v = atoi(v);}
          { case 2: v = atof(v);}
          {
             % default:
             throw ApplicationError, "Unsupported numeric type in option";
          }
        set_struct_field (s, names[i], v);
     }

   return s;
}

%}}}

% Partial support for changing hard limits on parameters of
% fit-functions implemented in slang.  This works only for
% functions that keep their defaults in a global structure
% named __${fname}_Defaults (and the only functions known
% to do that are those defined by the xspec module...)
%
define set_slangfun_param_hard_limits (fun_name, par_name, index, hard_min, hard_max) %{{{
{
   variable msg;

   variable default_struct_array_name = "__${fun_name}_Defaults"$;
   variable rs = __get_reference (default_struct_array_name);

   if (rs == NULL || 0 == _is_struct_type (@rs))
     {
        msg = "*** Sorry, I can't change the hard limits for function '$fun_name' (struct $default_struct_array_name is undefined)"$;
        vmessage(msg);
        return -1;
     }

   variable s = (@rs)[index];

   if (0 == isnan(hard_min))
     {
        if (s.min < hard_min)
          {
             msg = sprintf ("parameter %s.%s: hard minimum not changed.  (hard_min=%g is inconsistent with default minimum=%g)",
                            fun_name, par_name,
                            hard_min,
                            s.min);
             vmessage(msg);
             return -1;
          }
        s.hard_min = hard_min;
     }

   if (0 == isnan(hard_max))
     {
        if (s.max > hard_max)
          {
             msg = sprintf ("parameter %s.%s: hard maximum not changed.  (hard_max=%g is inconsistent with default maximum=%g)",
                            fun_name, par_name,
                            hard_max,
                            s.max);
             vmessage(msg);
             return -1;
          }
        s.hard_max = hard_max;
     }

   return 0;
}

%}}}
