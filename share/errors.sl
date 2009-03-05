#! /usr/bin/env jed-script

if (__argc != 3)
{
   vmessage ("Usage:  %s <infile> <outfile>", __FILE__);
   exit_jed (0);
}

private variable Input = __argv[1];
private variable Output = __argv[2];

TAB = 0;
_traceback = 1;
_debug_info = 1;

private define narrow_to_enums ()
{
   !if (bol_fsearch ("/*!+!"))
     return -1;
   
   () = down (1);
   push_mark ();

   !if (bol_fsearch ("/*!-!"))
     {
	pop_mark ();
	return -1;
     }
   
   () = up (1);
   narrow ();
   bob();
   
   return 0;
}

private define parse_symbol ()
{
   !if (re_fsearch ("[ \\t]*\\([a-zA-Z_]+\\)[,]?"))
     return NULL;

   return regexp_nth_match (1);
}

private define parse_comment ()
{
   variable value;
   
   !if (fsearch ("/*"))
     return NULL;
   !if (fsearch_char ('"'))
     return NULL;
   
   push_mark();
   
   !if (fsearch ("*/"))
     {
	pop_mark (1);
        return NULL;
     }

   value = bufsubstr();
   pop_mark (0);
   
   return value;
}

private define split ()
{
   variable s = struct 
     {
	symbol, comment
     };
   
   s.symbol = parse_symbol ();
   s.comment = parse_comment ();
   
   return s;
}

private variable Separator = ' ';

private define emit (s, buf)
{
   variable cbuf, c, sy, n, width = 80;
   
   if (s.comment == NULL or s.symbol == NULL)
     return 0;
   
   cbuf = whatbuf ();
   setbuf (buf);
   
   c = sprintf ("%c\"%s\"",
		   Separator, 
		   strtrim(strtrans (s.comment, "\"", "")));
   Separator = ',';
   
   sy = sprintf ("/* %s */\n",s.symbol);
   n = width - (strlen(c) + strlen(sy));
   
   insert (c);
   if (n > 0) whitespace (n);
   insert(sy);
   
   setbuf (cbuf);
   
   return 1;
}

private define parse (buf)
{
   variable s, cbuf;
   
   cbuf = whatbuf();
   setbuf (buf);
   insert ("/* automatically generated - do not edit */\n");
   setbuf (cbuf);
   
   narrow_to_enums ();
   
   forever
     {
	s = split ();
	!if (emit (s, buf))
	  break;
     }

   widen ();
   
   return 0;
}

private define extract_msg_strings (infile, outfile) 
{
   () = read_file (infile);
   
   if (parse (outfile))
     error ("*** parse failed");
   
   setbuf (outfile);
   write_buffer (outfile);
}

extract_msg_strings (Input, Output);

exit_jed (0);

