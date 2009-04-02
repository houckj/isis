/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Author:  John C. Houck  <houck@space.mit.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/* $Id: isis.c,v 1.128 2004/09/10 11:13:58 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#include <slang.h>

#include "isis.h"
#include "util.h"
#include "_isis.h"
#include "errors.h"

/*}}}*/

#define SYSTEM_ISISRC "etc/isis.rc"
#define USER_ISISRC ".isisrc"

extern int Isis_Preparse_Only;
extern int Isis_Debug_Mode;

unsigned int Isis_Version = ISIS_VERSION;

static char *Isis_Private_Namespace = "_isis";
static char Isis_Copyright [] =
    "Copyright (C) 1998-2008 Massachusetts Institute of Technology";

static int load_path_file (char *path, char *file) /*{{{*/
{
#ifdef HAVE_STAT
   struct stat st;
#endif
   int status = -1;

   if (file == NULL)
     return -1;

   if (path == NULL)
     path = SLpath_find_file_in_path (".", file);
   else
     path = SLpath_dircat (path, file);

   if (path == NULL)
     return -1;

#ifdef HAVE_STAT
   if ((-1 == stat (path, &st))
       || (0 == (st.st_mode & S_IRUSR)))
     {
        SLfree(path);
        return -1;
     }
#endif

   status = SLang_load_file (path);

   if ((status != 0) || (0 != SLang_get_error()))
     status = -1;

   SLfree (path);
   return status;
}

/*}}}*/

static int load_home_file (char *file) /*{{{*/
{
   char *home;

   if (NULL == (home = getenv ("ISIS_HOME")))
     {
        if (NULL == (home = getenv ("HOME")))
          return 0;
     }

   return load_path_file (home, file);
}

/*}}}*/

static int file_is_script (char *file) /*{{{*/
{
   FILE *fp;
   char buf[3];
   int is;

   if (file == NULL)
     return 0;

   if (NULL == (fp = fopen (file, "r")))
     return 0;

   is = ((NULL != fgets (buf, sizeof(buf), fp))
         && (buf[0] == '#') && (buf[1] == '!'));

   fclose (fp);
   return is;
}

/*}}}*/

/* environment setup */

static unsigned int version_from_version_file (char *srcdir) /*{{{*/
{
   static char *version_file = "isis-version";
   unsigned int version;
   FILE *fp;
   char *s;

   if ((srcdir == NULL)
       || (NULL == (s = SLpath_dircat (srcdir, version_file))))
     {
        return 0;
     }

   if (NULL == (fp = fopen (s, "r")))
     {
        SLfree (s);
        return 0;
     }

   if (1 != fscanf (fp, "%u", &version))
     {
        SLfree (s);
        fclose (fp);
        return 0;
     }

   fclose (fp);
   SLfree (s);

   return version;
}

/*}}}*/

static int setenv_isis_srcdir (void) /*{{{*/
{
   char *isis_srcdir;
   unsigned int version;

   isis_srcdir = getenv ("ISIS_SRCDIR");

   version = version_from_version_file (isis_srcdir);
   if (version == Isis_Version)
     {
        Isis_Srcdir = isis_srcdir;
        return 0;
     }

   if (-1 == putenv ("ISIS_SRCDIR=" SRCDIR))
     {
        fprintf (stderr, "*** Failed setting environment variable ISIS_SRCDIR= " SRCDIR "\n");
        return -1;
     }

   Isis_Srcdir = SRCDIR;

   version = version_from_version_file (Isis_Srcdir);
   if (version == Isis_Version)
     return 0;

   fprintf (stderr, "*** isis-%s initialization problem\n", ISIS_VERSION_STRING);
   fprintf (stderr, "Can't find the install directory matching this version of isis.\n");
   fprintf (stderr, "Either unset ISIS_SRCDIR, or set it to the correct path.\n");

   return -1;
}

/*}}}*/

static int setenv_atomdb_dir (void) /*{{{*/
{
   static char *atomdb_dir, *p;

   if (NULL != getenv ("ATOMDB"))
     return 0;

   p = ATOMDBDIR;
   if (p == NULL || *p == 0)
     return 0;

   atomdb_dir = "ATOMDB=" ATOMDBDIR;
   if (-1 == putenv (atomdb_dir))
     {
        fprintf (stderr, "Failed setting ATOMDB environment variable: %s\n", atomdb_dir);
        return -1;
     }

   return 0;
}

/*}}}*/

static int set_load_paths (char *prefix) /*{{{*/
{
   char *s;

   if (prefix == NULL)
     return -1;

   if (NULL == (s = SLpath_dircat (prefix, "share")))
     return -1;

   if (-1 == SLpath_set_load_path (s))
     {
        SLfree(s);
        return -1;
     }
   SLfree(s);

   if (NULL == (s = SLpath_dircat (prefix, "lib/modules")))
     return -1;

   if (-1 == SLang_set_module_load_path (s))
     {
        SLfree(s);
        return -1;
     }
   SLfree(s);

   return 0;
}

/*}}}*/

static int env_setup (void) /*{{{*/
{
   if (-1 == setenv_isis_srcdir ()
       || -1 == setenv_atomdb_dir ())
     return -1;

   return 0;
}

/*}}}*/

/* version, copyright, usage */

static void output_version (void) /*{{{*/
{
   fprintf (stdout, "ISIS version %s; S-Lang version %s\r\n",
            ISIS_VERSION_STRING, SLang_Version_String);
   fflush (stdout);
}

/*}}}*/

static int output_copyright (void) /*{{{*/
{
   output_version ();
   fprintf (stdout, "%s\r\n", Isis_Copyright);
   fprintf (stdout, "This is free software with ABSOLUTELY NO WARRANTY.\r\n");
   fprintf (stdout, "\n");
   fflush (stdout);

   return 0;
}

/*}}}*/

static void usage (char *pgm) /*{{{*/
{
   fprintf (stdout, "\n");
   output_copyright ();
   fprintf (stdout, "ISIS usage forms:\n\n");
   fprintf (stdout, "1.  %s [options] [FILE [args....]]\n\n", pgm);
   fprintf (stdout, "    Options:\n");
   fprintf (stdout, "     --help         Display this message\n");
   fprintf (stdout, "     --sldb [FILE]  Invoke S-Lang debugger\n");
   fprintf (stdout, "     --sldb-isis    Invoke S-Lang debugger on isis internals\n");
   fprintf (stdout, "     --prof [FILE]  Invoke S-Lang profiler on FILE\n");
   fprintf (stdout, "     --batch        Exit after loading FILE\n");
   fprintf (stdout, "     --script       Script mode\n");
   fprintf (stdout, "     --force-interactive   Prompt after loading FILE\n");
   fprintf (stdout, "     --secure       Secure mode\n");
   fprintf (stdout, "     --stdin        Read commands from stdin\n");
   fprintf (stdout, "     -              Read commands from stdin\n");
   fprintf (stdout, "     --version      Show version\n");
   fprintf (stdout, "     -Dsymbol       Define 'symbol' for #ifdef\n");
   fprintf (stdout, "     -a             Prompt has anonymous namespace\n");
   fprintf (stdout, "     -g             Compile with debugging info\n");
   fprintf (stdout, "     -i FILE        Use alternative init FILE (.isisrc is default)\n");
   fprintf (stdout, "     -n             Do not load .isisrc file\n");
   fprintf (stdout, "     -t             If isis_main or slsh_main exists, do not execute it\n");
   fprintf (stdout, "     -q n           Set verbosity (n < 0 means quieter)\n");
   fprintf (stdout, "     -v             Show verbose loading messages\n");
   fprintf (stdout, "\n");
   fprintf (stdout, "2.  %s-script [options] [FILE [arg, ....]]\n\n", pgm);
   fprintf (stdout, "    This form sets __argv[0] = FILE, __argv[1] = arg, ..\n\n");
   fprintf (stdout, "For more information see\n");
   fprintf (stdout, "       http://space.mit.edu/cxc/isis/\n");
   fprintf (stdout, "\n");
   fflush (stdout);
   exit_isis (0);
}

/*}}}*/

/* errno/user break */

static void isis_errno_hook (int err) /*{{{*/
{
   if ((err == (int) I_MALLOC_FAILED)
       || (Isis_Errno == (int) I_MALLOC_FAILED))
     isis_throw_exception (Isis_Error);
}

/*}}}*/

static int isis_user_break_hook (void) /*{{{*/
{
   return (SLang_get_error() == SL_USER_BREAK);
}

/*}}}*/

/* initialization */

static int check_for_alternate_name (char *pgm_basename) /*{{{*/
{
   char *buf;

   if (0 == strncmp (pgm_basename, "isis", 4))
     return 0;

   /* Look for companion startup script with .sl extension:
    * failure to find is not fatal, but any failure during load is */

   if (NULL == (buf = isis_mkstrcat (pgm_basename, ".sl", NULL)))
     return -1;

   if (-1 == SLang_load_file (buf))
     {
        if (SLang_get_error() != SL_IO_OPEN_ERROR)
          {
             ISIS_FREE(buf);
             return -1;
          }
        SLang_set_error (0);
     }

   sprintf (buf, "%s> ", pgm_basename);
   (void) SLang_run_hooks("set_prompt", 1, buf);

   ISIS_FREE(buf);

   return 0;
}

/*}}}*/

static int run_default_main (int testing) /*{{{*/
{
   char *mains[] = {"isis_main", "slsh_main", NULL};
   char **m;

   if (testing)
     return 0;

   for (m = mains; *m != NULL; m++)
     {
        if (2 == SLang_is_defined (*m))
          break;
     }

   if (*m == NULL)
     return 0;

   if (0 == SLang_load_string (*m))
     return 1;

   return 0;
}

/*}}}*/

static int initialize (int argc, char **argv) /*{{{*/
{
   static char *init_file = NULL;
   char *readline_namespace = "isis";
   char *script_name = NULL;
   char *pgm, *pgm_basename, *p_hyphen;
   int no_init_file = 0;
   int is_script = 0;
   int set_verbose = 0;
   int test_mode = 0;
   int force_interactive = 0;
   int ran_main = 0;
   int use_profiler = 0;

   pgm = argv[0];
   pgm_basename = SLpath_basename (pgm);
   argv++;
   argc--;

   p_hyphen = strrchr (pgm_basename, '-');
   is_script = ((p_hyphen != NULL) && (0 == strcmp (p_hyphen, "-script")));
   if (is_script)
     {
        struct stat st;
        if ((argc != 0)
            && (0 == stat (argv[0], &st)))
          {
             script_name = argv[0];
             argv++;
             argc--;
          }
     }

   Isis_Public_Namespace_Name = "Global";
   Isis_Errno_Hook = isis_errno_hook;
   Isis_User_Break_Hook = isis_user_break_hook;
   SLang_Error_Hook = error_hook;

   isis_set_pager (NULL);

   /* some args must be processed before initializing the intrinsics */
   if (is_script == 0 && argc == 1)
     {
        if (0 == strcmp (argv[0], "-V"))
          {
             fprintf (stdout, "%d\n", ISIS_VERSION);
             exit (0);
          }
        if (0 == strcmp (argv[0], "--version"))
          {
             output_version ();
             exit (0);
          }
        if (0 == strcmp (argv[0], "--preparse"))
          {
             Isis_Preparse_Only = 1;
             argc--;
             argv++;
          }
        else if (0 == strcmp (argv[0], "--secure"))
          {
             Isis_Secure_Mode=1;
             argc--;
             argv++;
          }
     }

   if ((-1 == SLang_init_all ())
       || (-1 == SLang_init_import ())
       || (-1 == SLang_init_array_extra ()))
     {
        fprintf (stderr, "*** S-Lang initialization failed\n");
        exit (EXIT_FAILURE);
     }

   SLang_Traceback = 4;
   SLang_generate_debug_info (1);
   (void) SLdefine_for_ifdef ("__ISIS__");

#ifdef __APPLE__
   (void) SLdefine_for_ifdef ("__APPLE__");
#endif

#ifdef HAVE_SLXPA_MODULE
   (void) SLdefine_for_ifdef ("__HAVE_SLXPA_MODULE__");
#endif

#ifdef HAVE_ISIS_EXTRAS
   (void) SLdefine_for_ifdef ("__HAVE_ISIS_EXTRAS__");
#endif

   Isis_Error = SLerr_new_exception (SL_Usage_Error, "IsisError", "Isis Error");

   if ((argc > 0) && (0 == strcmp (argv[0], "-T")))
     {
        Isis_Trace = 1;
        argc--;
        argv++;
     }

   /* check for verbose loading */
   if (((is_script == 0) || (script_name == NULL))
       && (argc > 0)
       && (0 == strcmp (argv[0], "-v")))
     {
        (void) SLang_load_file_verbose (1);
        argc--;
        argv++;
     }

   if ((argc > 0) && (0 == strcmp (argv[0], "--sldb")))
     {
        Isis_Debug_Mode = 1;
        force_interactive = 1;
        argc--;
        argv++;
     }

   if ((argc > 0) && (0 == strcmp (argv[0], "--sldb-isis")))
     {
        Isis_Debug_Mode = 2;
        force_interactive = 1;
        argc--;
        argv++;
     }

   if ((-1 == env_setup ())
       || (-1 == set_load_paths (Isis_Srcdir))
       || (-1 == init_isis_intrinsics (Isis_Private_Namespace)))
     {
        fprintf (stderr, "*** isis initialization failed\n");
        exit (EXIT_FAILURE);
     }

   if (is_script)
     {
        Isis_Verbose = -1;
        Isis_Batch_Mode = 1;
        /* remove '-script' from pgm_basename */
        *p_hyphen = 0;
     }

   if (is_script && script_name)
     goto skip_option_processing;

   while (argc > 0)
     {
        char *arg = argv[0];

        if (0 == strcmp (arg, "--script"))
          {
             is_script = 1;
             script_name=NULL;
             Isis_Verbose = -1;
             Isis_Batch_Mode = 1;
             argc--;
             argv++;
             continue;
          }

        if (0 == strcmp (arg, "--batch"))
          {
             Isis_Batch_Mode = 1;
             argc--;
             argv++;
             continue;
          }

        if (0 == strcmp (arg, "--force-interactive"))
          {
             force_interactive = 1;
             argc--;
             argv++;
             continue;
          }

        if (0 == strcmp (arg, "--prof"))
          {
             Isis_Verbose = -1;
             Isis_Batch_Mode = 1;
             use_profiler = 1;
             argc--;
             argv++;
             break;
          }

        if (0 == strcmp (arg, "--sldb"))
          {
             Isis_Debug_Mode = 1;
             force_interactive = 1;
             argc--;
             argv++;
             continue;
          }

        if (0 == strcmp (arg, "--help"))
          usage (pgm);

        if (0 == strncmp (arg, "-D", 2))
          {
             char *prep = arg + 2;
             if (*prep != 0)
               (void) SLdefine_for_ifdef (prep);
             argc--;
             argv++;
             continue;
          }

         if (0 == strcmp (arg, "-g"))
          {
             argc--;
             argv++;
             SLang_Traceback = 1;
             continue;
          }

        if (0 == strcmp (arg, "-a"))
          {
             readline_namespace = NULL;
             argc -= 1;
             argv += 1;
             continue;
          }

        if (0 == strcmp (arg, "-n"))
          {
             argc--;
             argv++;
             no_init_file = 1;
             continue;
          }

        if (0 == strcmp (arg, "-t"))
          {
             argc--;
             argv++;
             test_mode = 1;
             continue;
          }

        if (0 == strcmp (arg, "-v"))
          {
             argc--;
             argv++;
             (void) SLang_load_file_verbose (1);
             continue;
          }

        if ((0 == strcmp (arg, "--stdin"))
            || (0 == strcmp (arg, "-")))
          {
             argc--;
             argv++;
             input_from_stdin();
             continue;
          }

        if ((0 == strcmp (arg, "-i"))
            && (argc > 1))
          {
             init_file = argv[1];
             argc -= 2;
             argv += 2;
             continue;
          }

        if ((0 == strcmp (arg, "-q"))
            && (argc > 1))
          {
             Isis_Verbose = strtol (argv[1], NULL, 0);
             argc -= 2;
             argv += 2;
             set_verbose = 1;
             continue;
          }

        if (is_script
            && (script_name == NULL))
          {
             script_name = arg;
             argv++;
             argc--;
             break;
          }

        break;
     }

   skip_option_processing:

   /* At this point, argc is either 0, or argv points to the first
    * valid argument for the program.
    */

   if (is_script)
     {
        if (script_name == NULL)
          usage (pgm);
        argv--;
        argc++;
     }

   if ((argc > 0)
       && file_is_script (script_name))
     {
        is_script = 1;
        Isis_Batch_Mode = 1;
     }

   if (0 == isatty (fileno(stdout)))
     {
        input_from_stdin();
     }
   else if (0 == isatty (fileno(stdin)))
     {
        input_from_stdin();
     }

   if (Isis_Batch_Mode)
     {
        /* unless CL args say otherwise, run a bit quieter in batch mode */
        if ((Isis_Verbose == 0) && (set_verbose == 0))
          Isis_Verbose = FAIL;
     }
   else (void) SLdefine_for_ifdef ("__INTERACTIVE__");

   if (-1 == SLang_set_argc_argv (argc, argv))
     {
        fprintf (stderr, "*** Failed handling command-line arguments\n");
        exit_isis (1);
     }

   init_signals ();

   if (no_init_file == 0)
     {
        load_path_file (Isis_Srcdir, SYSTEM_ISISRC);
        if (init_file)
          load_path_file (NULL, init_file);
        else
          load_home_file (USER_ISISRC);
     }

   if (-1 == check_for_alternate_name (pgm_basename))
     exit_isis (1);

   if (use_profiler)
     {
        if (-1 == SLang_load_file ("isisprof"))
          exit_isis (1);
     }
   else if (Isis_Debug_Mode == 0)
     {
        if (is_script)
          {
             (void) SLang_load_file (argv[0]);
          }
        else if (argc > 0)
          {
             if (-1 == load_path_file (NULL, argv[0]))
               fprintf (stderr, "Failed loading file %s\n", argv[0]);
          }
        ran_main = run_default_main (test_mode);
     }

   if (SLang_get_error())
     exit_isis (1);

   if (force_interactive == 0)
     {
        if (Isis_Batch_Mode || ran_main)
          exit_isis (0);
     }

   if (isatty(fileno(stdin)) && (is_script == 0))
     {
        fprintf (stdout, "\nWelcome to ISIS Version %s\n%s\n\n",
                 ISIS_VERSION_STRING, Isis_Copyright);
        fflush(stdout);
     }

   if (-1 == open_readline (readline_namespace))
     exit_isis (1);

   if (Isis_Debug_Mode && (argc > 0))
     {
        (void) SLang_load_file ("isisdb");
        if (is_script)
          exit_isis (0);
     }

   return 0;
}

/*}}}*/

#ifdef MODULE

SLANG_MODULE(isis);
int init_isis_module_ns (char *ns_name) /*{{{*/
{
   SLang_NameSpace_Type *ns;

   /* silence compiler warning about functions defined but not used */
   if (0) {(void) initialize;}

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return -1;

   Isis_Public_Namespace_Name = SLang_create_slstring (ns_name);

   (void) SLdefine_for_ifdef ("__ISIS_MODULE__");

   Isis_Error = SLerr_new_exception (SL_Usage_Error, "IsisError", "Isis Error");

   if ((-1 == env_setup ())
       || (-1 == init_isis_intrinsics (Isis_Private_Namespace)))
     {
        fprintf (stderr,  "*** Failed loading isis module\n");
        SLang_free_slstring (Isis_Public_Namespace_Name);
        Isis_Public_Namespace_Name = NULL;
        isis_throw_exception (Isis_Error);
        return -1;
     }

   return 0;
}

/*}}}*/

void deinit_isis_module (void) /*{{{*/
{
   SLang_free_slstring (Isis_Public_Namespace_Name);
   quit_isis (0);
}

/*}}}*/

#else   /* MODULE */

static void check_slang_version (void) /*{{{*/
{
   if (SLang_Version >= SLANG_VERSION)
     return;

   fprintf (stderr, "***Warning: Executable compiled against S-Lang %d but linked to %d\n",
            SLang_Version, SLANG_VERSION);
   fflush (stderr);
   sleep (2);
}

/*}}}*/

int main (int argc, char **argv) /*{{{*/
{
   check_slang_version ();

   if (-1 == initialize (argc, argv))
     exit (EXIT_FAILURE);

   if (-1 == SLang_set_argc_argv (argc, argv))
     {
        fprintf (stderr, "*** Failed handling command-line arguments\n");
        exit_isis (1);
     }

   interpreter_loop ();

   return 0;
}

/*}}}*/

#endif   /* MODULE */

/*{{{ fortran entry */

/* This is to allow the use of Fortran compilers that insist
 * on having MAIN__ (or something like that) as the internal
 * entry point for the main program.
 */

#ifdef FC_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
int FC_DUMMY_MAIN (void);
int FC_DUMMY_MAIN (void)
{
   static char msg[] = "*** isis: FC_DUMMY_MAIN: this should never happen\n";
   write (STDERR_FILENO, msg, sizeof(msg));
   return 1;
}
#endif

/*}}}*/
