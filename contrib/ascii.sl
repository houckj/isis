
%
%  (Note that ISIS already has this functionality built-in)
%

_debug_info = 1;
	
define writecol ()
{
   if (_NARGS < 2)
     {
	_pop_n (_NARGS);
	vmessage ("Usage:  %s (fp, a, b, ....);", _function_name);
	return;
     }
   variable args = __pop_args (_NARGS - 1);
   variable fp = ();
   
   if (String_Type == typeof (fp))
     {
	variable _fp = fopen (fp, "w");
	if (_fp == NULL)
	  verror ("Unable to open %s", fp);
	fp = _fp;
     }

   variable fmt = "";
   loop (_NARGS-1)
     fmt = strcat ("%S\t", fmt);
   fmt = strcat (strtrim (fmt), "\n");

   variable status = array_map (Int_Type, &fprintf, fp, fmt, __push_args (args));
   if (length (where (status == -1)))
     verror ("Error writing to file");
}

private define chop (str)
{
   strchop (strcompress(str, "\t "), '\t', 0);
}

define readcol ()
{
   variable fp, cols =NULL;
   
   switch (_NARGS)
     {
      case 2:
	(fp, cols) = ();
     }
     {
      case 1:
	fp = ();
     }
     {
	% default:
	
	message ("% Usage:  (a, b, c, ...) = readcol (fp, [cols] )");
	return;
     }
   
   if (String_Type == typeof (fp))
     {
	variable _fp = fopen (fp, "r");
	if (_fp == NULL)
	  verror ("Unable to open %s", fp);
	fp = _fp;
     }

   variable lines = fgetslines (fp);
   
   variable not_a_comment = array_map (Integer_Type, &strncmp, lines, "#", 1);
   
   lines = lines[where(not_a_comment)];
   variable nlines = length(lines);

   variable ncols;
   
   if (NULL != cols)
     ncols = length(cols);
   else
     {
	ncols = length(chop(lines[0]));
	cols = [0:ncols-1];
     }
   
   variable values = Double_Type [nlines, ncols];

   variable col, fields, nfields, k, line;
   
   line = 0;
   _for (0, nlines-1, 1)
     {
	k = ();
	
	fields = chop (lines[k]);
	nfields = length(fields);
	
	_for (0, ncols-1, 1)
	  {
	     col = ();
	     if (cols[col] >= nfields)
		  break;
	     values[line, col] = atof(fields[cols[col]]);
	  }
	
	line++;
     }
   
   _for (0, ncols-1, 1)            % return as 1-D arrays
     {
	col = ();
	values [[0:line-1], col];
     }
}

