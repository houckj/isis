/* This file borrowed from the JED distribution and slightly modified
 * by John C. Houck 1999/09/30
 * 
 * $Id: chkslang.c,v 1.3 2003/04/18 21:29:46 houck Exp $
 */

/* Copyright (c) 1992, 1998 John E. Davis
 * This file is part of JED editor library source.
 *
 * You may distribute this file under the terms the GNU General Public
 * License.  See the file COPYING for more information.
 */
/* It is too bad that this cannot be done at the preprocessor level.
 * Unfortunately, C is not completely portable yet.  Basically the #error
 * directive is the problem.
 */

#include "config.h"

#include <stdio.h>
#ifdef VMS
#include <ssdef.h>
#endif
#include <slang.h>

#ifndef ISIS_SRC                         
#include "jdmacros.h"
#endif

static char *make_version (unsigned int v)
{
   static char v_string[16];
   unsigned int a, b, c;
   
   a = v/10000;
   b = (v - a * 10000) / 100;
   c = v - (a * 10000) - (b * 100);
   sprintf (v_string, "%u.%u.%u", a, b, c);
   return v_string;
}

int main (int argc, char **argv)
{
   unsigned int min_version, sl_version;
   unsigned int sug_version;
   int ret;
   
   if ((argc < 3) || (argc > 4))
     {
	fprintf (stderr, "Usage: %s <PGM> <SLANG-VERSION> <SUGG VERSION>\n", argv[0]);
	return 1;
     }
#ifndef SLANG_VERSION
   sl_version = 0;
#else
   sl_version = SLANG_VERSION;
/* ISIS_SRC   
#ifdef REAL_UNIX_SYSTEM
*/
   if (SLang_Version != SLANG_VERSION)
     {
	fprintf (stderr, "\n\n******\n\
slang.h does not match slang library version.  Did you install slang as\n\
as a shared library?  Did you run ldconfig?  You have an installation problem\n\
and you will need to check the SLANG include and library paths and properly\n\
set them.  Also try: make clean; make\n******\n\n");
	fprintf (stderr, "\t*****       slang.h is version = %d *****\n", SLANG_VERSION);
	fprintf (stderr, "\t***** slang library is version = %d *****\n\n", SLang_Version);
	return 1;
     }
/* ISIS_SRC   
#endif
*/
#endif
   
   sscanf (argv[2], "%u", &min_version);
   if (argc == 4) sscanf (argv[3], "%u", &sug_version);
   else sug_version = sl_version;

   ret = 0;
   if (sl_version < min_version)
     {
	fprintf (stderr, "This version of %s requires slang version %s.\n",
		 argv[1], make_version(min_version));
#ifdef VMS
	ret = SS$_ABORT;
#else
	ret = 1;
#endif
     }
   
   if (sl_version < sug_version)
     {
	fprintf (stderr, "Your slang version is %s.\n", make_version(sl_version));
	fprintf (stderr, "To fully utilize this program, you should upgrade the slang library to\n");
	fprintf (stderr, "  version %s\n", make_version(sug_version));
	fprintf (stderr, "This library is available via anonymous ftp from\n\
space.mit.edu in pub/davis/slang.\n");
     }
   
   return ret;
}




