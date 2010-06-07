% The code here constructs an index for the low-level routines.  It uses
% the index page at listed off
%  http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/cfitsio.html
% Each function is indexed using:
%  <TR><TD ALIGN="LEFT">fits_clear_errmark</TD>
%  <TD ALIGN="RIGHT"><A HREF="node33.html#ffpmrk"><img src="crossref.png" ALIGN="BOTTOM" BORDER="1" ALT="[*]"></A></TD>
%  </TR>

static variable Index_Page = "node122.html";
static variable CFitsio_Root_URL 
  = "http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user";

static variable Wget_Cmd 
  = "wget --quiet";

static define retrieve_url (url, outfile)
{
   vmessage ("Obtaining the cfitsio index from heasarc...");
   () = system (sprintf ("%s -O '%s' '%s'", Wget_Cmd, outfile, url));
}

static define mkindex ()
{
   variable cbuf = whatbuf ();
   variable index_url = strcat (CFitsio_Root_URL, "/", Index_Page);
   variable tmpfile = sprintf ("/tmp/cfitsio_index_%d_%d", getpid(),_time());
   retrieve_url (index_url, tmpfile);
   () = read_file (tmpfile);
   () = remove (tmpfile);

   variable index = Assoc_Type[String_Type];

   bob ();
   while (bol_fsearch ("<TR><TD ALIGN=\"LEFT\">fits_"))
     {
	() = ffind ("fits_");
	push_mark ();
	() = ffind ("</TD>");
	variable fun = bufsubstr ();
	go_down(1);
	if (0 == ffind ("<A HREF=\""))
	  continue;
	go_right (9);
	push_mark ();
	() = ffind_char ('"');
	variable node = bufsubstr ();
	index[fun] = node;
     }
   delbuf (whatbuf());
   setbuf (cbuf);
   return index;
}

static variable CFitsio_Function_Index = NULL;
static define cfitsio_fun_url (fun)
{
   if (CFitsio_Function_Index == NULL)
     CFitsio_Function_Index = mkindex ();
   variable url;
   if (assoc_key_exists (CFitsio_Function_Index, fun))
     url = CFitsio_Function_Index[fun];
   else
     {
	() = fprintf (stderr, "Warning: Unable to find a URL for %s\n", fun);
	url = "";
     }
   vinsert ("%s/%s", CFitsio_Root_URL, url);
}
tm_add_macro ("cfitsio_fun_url", &cfitsio_fun_url, 1, 1);
