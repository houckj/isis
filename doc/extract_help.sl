% -*- mode: SLang; mode: fold -*-

% Usage:  jed -batch -l extract_help.sl

CASE_SEARCH=1;
_debug_info = 1;
_traceback = 1;

typedef struct
{
   beg, end
}
Region_Type;

typedef struct
{
   name, purpose, usage, seealso, description
}
Entry_Fields;

static define find_region (entry) %{{{
{
   !if (bol_fsearch (entry.beg))
     return 0;

   push_mark ();

   !if (bol_fsearch (entry.end))
     return 0;

   go_down_1();

   return 1;
}

%}}}

static define strip_bol_comments (sym) %{{{
{
   bob();
   while (not(eobp()))
     {
	if (0 == bol_fsearch (sym))
	  break;
	del_eol();
     }
   bob();
}

%}}}

static define regex_replace (old, new) %{{{
{
   bob ();

   forever
     {
	!if (re_fsearch (old))
	  break;
	!if (replace_match (new, 1))
	  break;
     }

   bob();
}

%}}}

static define indent_buffer () %{{{
{
   eob ();
   variable num = what_line();
   bob ();

   loop (num)
     {
	insert_spaces (4);
	go_down_1();
     }
}

%}}}

static define reflow_paragraphs () %{{{
{
   bob();
   WRAP=65;
   while (not(eobp()))
     {
	skip_chars ("\t \n");
	call ("format_paragraph");
	forward_paragraph ();
     }
   bob();
}

%}}}

static define filter_text () %{{{
{
   check_region (0);
   narrow ();
   bob();

   replace ("\\isisx", "ISIS");
   replace ("\\slang", "S-Lang");
   replace ("\\cfitsio", "CFITSIO");
   replace ("\\pgplot", "PGPLOT");
   replace ("\\xspec", "XSPEC");
   replace ("\\HRMA", "HRMA");
   replace ("\\verb", "");
   replace ("\\left", "");
   replace ("\\right", "");
   replace ("\\leq", "<=");
   replace ("\\times", " x ");
   replace ("\\ldots", "...");
   replace ("\\tt", "");
   replace ("\\sc", "");
   replace ("\\bf", "");
   replace ("\\it", "");
   replace ("\\rm", "");
   replace ("|", "");
      replace ("'$'", "__DOLLAR__");
      replace ("$", "");
      replace ("__DOLLAR__", "'$'");
   replace ("~", "");
   replace ("{", "");
   replace ("}", "");
   replace ("``", "\"");
   replace ("''", "\"");
   replace ("\\AA", "A");
   replace ("\\S\\", "*section*");
   replace ("\\S ", "*section* ");
   replace ("\\", "");

   strip_bol_comments ("%");

   reflow_paragraphs ();

   widen ();
}

%}}}

static define massage_value () %{{{
{
   variable cbuf = whatbuf ();
   check_region (0);
   narrow ();

   setbuf ("tmp");
   erase_buffer ();
   insbuf (cbuf);

   bob();

   replace ("\\vspace*{\\baselineskip}", "");
   regex_replace ("\\\\index{[-0-9a-zA-Z _?\\\\]*[!@|]?[-0-9a-zA-Z (){},_?\\\\]*}", "");
   regex_replace ("\\\\ref{[-0-9a-zA-Z_:]*}", "");

  while (not(eobp()))
     {
	push_mark ();
	if (0 == bol_fsearch ("\\begin"))
	  {
	     eob();
	     filter_text();
	     break;
	  }
	del_eol();
	insert_char ('\n');

	filter_text();

	() = bol_fsearch ("\\end");
	del_eol();
     }

   indent_buffer();
   trim_buffer();

   bob ();
   push_mark ();
   eob ();
   variable s = bufsubstr ();

   setbuf (cbuf);
   widen ();

   return s;
}

%}}}

static define get_massaged_text (entry) %{{{
{
   go_down_1();
   if (looking_at (entry.end))
     return "";

   push_mark ();
   () = bol_fsearch (entry.end);
   go_up_1();

   return massage_value ();
}

%}}}

static define get_delim_value () %{{{
{
   bol ();
   go_down_1 ();

   !if (fsearch_char ('{'))
     return "";

   push_mark ();

   if (1 != find_matching_delimiter ('{'))
     {
	pop_mark (1);
	return "";
     }

   variable value = bufsubstr();
   pop_mark (0);

   if (value == NULL)
     return "";

   value = str_replace_all (value, "\\sc ", "");
   value = str_replace_all (value, "\\tt ", "");
   value = str_replace_all (value, "\\isisx", "ISIS");
   value = str_replace_all (value, "{", "");
   value = str_replace_all (value, "}", "");
   value = str_replace_all (value, "$", "");
   value = str_replace_all (value, "~", "");
   value = str_replace_all (value, "\\", "");

   return strtrim (value);
}

%}}}

static define indent_wrap (s, pad, line_width) %{{{
{
   line_width -= strlen(pad);

   variable toks = strtok (s, ",");

   variable _t,
     len = 0,
     new_s = "",
     delim = "";

   foreach _t (toks)
     {
        variable t = strtrim(_t);
        variable nt = strlen(t) + 2;

        if (len + nt > line_width)
          {
             new_s += ",\n$pad"$;
             delim = "";
             len = 0;
          }

        len += nt;
        new_s += "${delim}${t}"$;
        delim = ", ";
     }

   return new_s;
}

%}}}

static define parse_region (entry) %{{{
{
   bob();
   !if (looking_at (entry.beg))
     return;

   variable f = @Entry_Fields;

   f.name = get_delim_value ();
   f.purpose = get_delim_value ();
   f.usage = get_delim_value ();
   f.seealso = indent_wrap (get_delim_value (), "    ", 72);
   f.description = get_massaged_text (entry);

   return f;
}

%}}}

static define insert_fields (f, buf) %{{{
{
   variable cbuf = whatbuf();
   setbuf (buf);
   insert (sprintf ("%s\n\n SYNOPSIS\n    %s\n\n USAGE\n    %s\n\n",
		    f.name,
		    f.purpose,
		    f.usage));
   insert (sprintf (" DESCRIPTION\n%s\n\n", f.description));
   insert (sprintf (" SEE ALSO\n    %s\n\n", f.seealso));
   insert ("------------------------------------------------------------------------\n");
   setbuf (cbuf);
}

%}}}

static define parse_file (entry, buf) %{{{
{
   variable f;

   forever
     {
	!if (find_region (entry))
	  return;

	check_region (0);
	narrow ();
	f = parse_region (entry);
	insert_fields (f, buf);
	widen ();
     }
}

%}}}

variable Syn = "whatever";
create_syntax_table (Syn);
define_syntax ("{", "}", '(', Syn);

static define extract_help (entry, infile, outfile) %{{{
{
   () = read_file (infile);
   use_syntax_table (Syn);

   parse_file (entry, outfile);
   setbuf (outfile);
   write_buffer (outfile);
}

%}}}

variable entry = @Region_Type;
entry.beg = "\\begin{isisfunction}";
entry.end = "\\end{isisfunction}";

extract_help (entry, "manual.tex", "help.txt");

exit_jed (0);
