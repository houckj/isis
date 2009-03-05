%quit_jed;
!if (is_defined ("__argv"))
{
   message ("You need a newer version of jed to run this script");
   quit_jed ();
}

if (__argc != 4)
{
   message ("Usage: jed -script fixtex.sl <filename>");
   quit_jed ();
}

variable file = __argv[3];
() = read_file (file);

% Patch up the >,< signs
bob ();
replace ("$<$", "<");
replace ("$>$", ">");

% It appears that sgml2tex screws up _for in section titles, producing \_{for}.
replace ("ion\\_{", "ion{\\_");

% Make the first chapter a preface
bob ();
if (bol_fsearch ("\\chapter{Preface}"))
{
   push_spot ();
   push_mark ();
   go_right (8); insert ("*");	       %  \chapter{ --> \chapter*{
   () = bol_fsearch ("\\chapter{");
   push_spot ();

   insert("\\tableofcontents\n");
   eol ();
   insert ("\n\\pagenumbering{arabic}");

   pop_spot ();
   narrow ();
   bob ();
   replace ("\\section{", "\\section*{");
   widen ();

   if (bol_bsearch ("\\tableofcontents"))
     delete_line ();
   
   pop_spot ();
   if (bol_bsearch ("\\maketitle"))
     insert ("\\pagenumbering{roman}\n");

}

static define fixup_urldefs ()
{
   % pdflatex cannot grok urldef
   bob ();
   while (bol_fsearch("\\urldef{") and ffind ("\\url{"))
     {
	variable line = line_as_string ();
	bol ();
	insert ("\\ifpdf\n");
	
	deln (7); insert ("\\newcommand");
	push_mark ();
	()=ffind ("}");
	variable macro = bufsubstr ();
	() = ffind ("\\url");
	go_left (1);
	trim ();
	insert("{");eol(); insert("}");
	insert ("\n\\else\n");
	insert (line); newline ();
	%insert ("\\newcommand"); insert(macro); insert("}{}\n");
	insert ("\\fi\n");
     }
}

fixup_urldefs ();
save_buffer ();
quit_jed ();


