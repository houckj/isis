/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2008  Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Authors:  John C. Houck  <houck@space.mit.edu>
              John E. Davis  <davis@space.mit.edu>
              Michael S. Noble <mnoble@space.mit.edu>

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

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <signal.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef __linux__
# ifdef __GLIBC__
#  include <fpu_control.h>
# else
#  include <i386/fpu_control.h>
# endif
# ifndef _FPU_SETCW
#  define _FPU_SETCW __setfpucw
# endif
#endif

#include <slang.h>
#include "isis.h"
#include "_isis.h"
#include "util.h"
#include "errors.h"

#ifdef HAVE_GNU_READLINE
static char *Gnu_Rl_History_File;
# include <readline/readline.h>
# include <readline/history.h>
#endif

#ifdef HAVE_GNU_READLINE
# define USE_SLANG_READLINE 0
#else
# define USE_SLANG_READLINE 1
#endif

static int Use_SLang_Readline = USE_SLANG_READLINE;
static char *isis_prompt_default = "isis> ";
static char *isis_prompt;

#ifndef SLANGLIBDIR
# define SLANGLIBDIR NULL
#endif
static char *slang_libdir = SLANGLIBDIR;

/*}}}*/

int Isis_Silent_Errors;
int Isis_Preparse_Only;
int Isis_Debug_Mode;
char *Isis_Public_Namespace_Name;

static pid_t Isis_Initial_Pid;

static int Input_From_Stdin;
static int TTY_Inited;
static char *Isis_Readline_Namespace_Name;

static SLang_Load_Type *Readline_Load_Object = NULL;

static void close_readline (void);
static void sig_sigtstp (int sig);
static void (*last_sig_sigtstp) (int);

static int enable_fpe (int);
static void check_fpe (void);

/*{{{ slang2 compatibility */

static SLang_RLine_Info_Type *RLI;

/*}}}*/

/*{{{ tty configuration */

#ifdef HAVE_GNU_READLINE
static void (*last_sig_sigint) (int);
static void gnu_rl_sigint_handler (int sig) /*{{{*/
{
   (void) sig;
   rl_delete_text (0, rl_end);
   rl_point = rl_end = 0;
   fprintf (stdout, "\n");
   rl_forced_update_display ();
}

/*}}}*/
#endif

/* This hook if a signal occurs while waiting for input. */
static int getkey_intr_hook (void)
{
   return SLang_handle_interrupt ();
}

static void init_tty (void) /*{{{*/
{
#ifdef HAVE_GNU_READLINE
   if (Use_SLang_Readline == 0)
     {
        SLsig_block_signals ();
        last_sig_sigint = SLsignal (SIGINT, gnu_rl_sigint_handler);
        last_sig_sigtstp = SLsignal (SIGTSTP, sig_sigtstp);
        SLsig_unblock_signals ();
        return;
     }
#endif

   if (TTY_Inited)
     return;
   TTY_Inited++;

   SLsig_block_signals ();
   SLang_TT_Read_FD = fileno (stdin);
   last_sig_sigtstp = SLsignal (SIGTSTP, sig_sigtstp);

   if (-1 == SLang_init_tty (-1, 1, 0))
     {
        SLsignal (SIGTSTP, last_sig_sigtstp);
        SLsig_unblock_signals ();
        fprintf (stderr, "Error initializing terminal.");
        exit (EXIT_FAILURE);
     }

   SLang_getkey_intr_hook = getkey_intr_hook;

   SLtty_set_suspend_state (1);
   SLsig_unblock_signals ();
}

/*}}}*/

static void reset_tty (void) /*{{{*/
{
#ifdef HAVE_GNU_READLINE
   if (Use_SLang_Readline == 0)
     {
        SLsig_block_signals ();
        SLsignal (SIGINT, last_sig_sigint);
        SLsignal (SIGTSTP, last_sig_sigtstp);
        SLsig_unblock_signals ();
        return;
     }
#endif

   if (TTY_Inited == 0)
     return;
   TTY_Inited = 0;

   SLsig_block_signals ();
   SLsignal (SIGTSTP, last_sig_sigtstp);
   SLang_reset_tty ();
   SLsig_unblock_signals ();
   fputs ("\r\n", stdout);   /* need this with buggy Linux rlogin? */
   fflush (stdout);
}

/*}}}*/

static void sig_sigtstp (int sig) /*{{{*/
{
   (void) sig;
   SLsig_block_signals ();
   reset_tty ();
   kill(getpid(),SIGSTOP);
   init_tty ();

   if (Use_SLang_Readline == 0)
     {
#ifdef HAVE_GNU_READLINE
        rl_refresh_line (0,0);
#endif
     }
   else SLrline_redraw (RLI);

   SLsig_unblock_signals ();
}

/*}}}*/

/*}}}*/

/*{{{ init/exit functions */

#ifdef WITH_XSPEC_STATIC_LINKED
extern int init_xspec_module_ns (char *);
extern void deinit_xspec_module (void);
#endif

static int init_optional_modules (char *ns_name) /*{{{*/
{
   (void) ns_name;

#ifdef WITH_XSPEC_STATIC_LINKED
   if (-1 == init_xspec_module_ns (Isis_Public_Namespace_Name))
     return isis_trace_return(-1);
   (void) SLdefine_for_ifdef ("__XSPEC_IS_STATIC__");
#endif

   return 0;
}

/*}}}*/

static void deinit_optional_modules (void) /*{{{*/
{
#ifdef WITH_XSPEC_STATIC_LINKED
   deinit_xspec_module ();
#endif
}

/*}}}*/

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
}
#endif
extern int init_readline_module_ns (char *);
extern int init_miscio_module_ns (char *);
extern int init_math_module_ns (char *);
extern int init_plot_module_ns (char *);
extern int init_hist_module_ns (char *);
extern int init_atom_module_ns (char *);
extern int init_emis_module_ns (char *);
extern int init_fit_module_ns (char *);
extern int init_model_module_ns (char *);
#if 0
{
#endif
#ifdef __cplusplus
}
#endif

static int sldb_load_hook (char *file) /*{{{*/
{
   if (file == NULL)
     return -1;

   if (1 != SLang_run_hooks ("sldb", 1, file))
     return -1;

   return 0;
}

/*}}}*/

int init_isis_intrinsics (char *ns_name) /*{{{*/
{
   enable_fpe (1);
   Isis_Initial_Pid = getpid();

   /* init_readline_ must come first
    * fits module must come before init.sl
    */
   if (-1 == init_readline_module_ns (ns_name)
       || -1 == init_miscio_module_ns (ns_name)
       || -1 == init_math_module_ns (ns_name)
       || -1 == init_plot_module_ns (ns_name)
       || -1 == init_hist_module_ns (ns_name)
       || -1 == init_atom_module_ns (ns_name)
       || -1 == init_emis_module_ns (ns_name)
       || -1 == init_fit_module_ns (ns_name)
       || -1 == init_model_module_ns (ns_name)
       || -1 == init_optional_modules (ns_name)
       || -1 == SLang_load_file ("setup"))
     {
        return isis_trace_return(-1);
     }

   if (Isis_Debug_Mode == 2)
     {
        if (-1 == SLang_load_file ("isisdb"))
          return isis_trace_return(-1);
        SLang_Load_File_Hook = &sldb_load_hook;
     }

   if (-1 == SLang_load_file ("init"))
     return isis_trace_return(-1);

   if (Isis_Debug_Mode == 2)
     {
        SLang_Load_File_Hook = NULL;
        Isis_Debug_Mode = 0;
     }

   return 0;
}

/*}}}*/

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
}
#endif
extern void deinit_readline_module (void);
extern void deinit_miscio_module (void);
extern void deinit_math_module (void);
extern void deinit_plot_module (void);
extern void deinit_hist_module (void);
extern void deinit_atom_module (void);
extern void deinit_emis_module (void);
extern void deinit_fit_module (void);
extern void deinit_model_module (void);
#if 0
{
#endif
#ifdef __cplusplus
}
#endif

static int Called_Quit_Isis = 0;

void quit_isis (int reset) /*{{{*/
{
   /* If a signal arrives while we're exiting, it's possible
    * that this routine could be called again, before it has
    * finished.  Set a global variable to prevent that.
    */
   if (reset == 0)
     {
        if (Called_Quit_Isis)
          return;
        Called_Quit_Isis = 1;
     }

   deinit_miscio_module ();
   deinit_math_module ();
   deinit_plot_module ();
   deinit_hist_module ();
   deinit_atom_module ();
   deinit_emis_module ();
   deinit_model_module ();
   deinit_fit_module ();
   deinit_optional_modules ();

   if (reset)
     init_fit_module_internals ();
}

/*}}}*/

typedef struct _AtExit_Type
{
   SLang_Name_Type *nt;
   struct _AtExit_Type *next;
}
AtExit_Type;

static AtExit_Type *AtExit_Hooks;

/* at_exit() shamelessly stolen from slsh */
static void at_exit (SLang_Ref_Type *ref) /*{{{*/
{
   SLang_Name_Type *nt;
   AtExit_Type *a;

   if (NULL == (nt = SLang_get_fun_from_ref (ref)))
     return;

   a = (AtExit_Type *) SLmalloc (sizeof (AtExit_Type));
   if (a == NULL)
     return;

   a->nt = nt;
   a->next = AtExit_Hooks;
   AtExit_Hooks = a;
}

/*}}}*/

static void call_at_exit_hooks (void) /*{{{*/
{
   while (AtExit_Hooks != NULL)
     {
	AtExit_Type *next = AtExit_Hooks->next;
        if (SLang_get_error() == 0)
          (void) SLexecute_function (AtExit_Hooks->nt);
	SLfree ((char *) AtExit_Hooks);
	AtExit_Hooks = next;
     }
}

/*}}}*/

void exit_isis (int err) /*{{{*/
{
   int status = 0;

   /* Use _exit() if this process was created by fork() */
   if (Isis_Initial_Pid != getpid())
     {
        if (err || SLang_get_error())
          status = EXIT_FAILURE;
        _exit (status);
     }

   SLang_run_hooks ("stop_log", 1, NULL);
   call_at_exit_hooks ();
   quit_isis (0);
   deinit_readline_module ();

   if (err || SLang_get_error())
     status = EXIT_FAILURE;

   exit (status);
}

/*}}}*/

static void _quit (int *err) /*{{{*/
{
   int s = 0;
   if (err) s = *err;
   exit_isis (s);
}

/*}}}*/

static void reset_isis (void) /*{{{*/
{
   quit_isis (1);
}

/*}}}*/

static volatile sig_atomic_t Exit_Signal_In_Progress = 0;
static void sig_exit_isis (int sig) /*{{{*/
{
   if (Isis_Initial_Pid != getpid ())
     _exit (EXIT_FAILURE);

   if (Exit_Signal_In_Progress)
     return;

   Exit_Signal_In_Progress = 1;
   SLsig_block_signals ();

   fprintf (stderr, "isis[%d]:  Killed by signal %d.\n", getpid (), sig);

   exit_isis (sig);
}

/*}}}*/

#ifdef __GLIBC__
# include <execinfo.h>
#endif
static void print_backtrace (void)
{
#ifdef __GLIBC__
   static char msg1[] = "\n**** isis:  Backtrace start\n\n";
   static char msg2[] = "\n**** isis:  Backtrace end\n****        Please save this backtrace for debugging.\n";
   void *bt[64];
   int size = backtrace (bt, sizeof(bt) / sizeof(*bt));
   write (STDERR_FILENO, msg1, sizeof(msg1)-1);
   backtrace_symbols_fd (bt, size, STDERR_FILENO);
   write (STDERR_FILENO, msg2, sizeof(msg2)-1);
#endif
}

static volatile sig_atomic_t Abort_Signal_In_Progress = 0;
static void sig_abort_isis (int sig) /*{{{*/
{
   static char msg[] =
"\n**** isis:  Fatal Error:  Please report this bug to houck@space.mit.edu.\n\n";

   (void) sig;

   if (Abort_Signal_In_Progress)
     abort();

   Abort_Signal_In_Progress = 1;

   SLsignal (SIGSEGV, SIG_DFL);
   SLsignal (SIGBUS, SIG_DFL);
   SLsignal (SIGILL, SIG_DFL);

   SLsig_block_signals ();

   print_backtrace ();
   write (STDERR_FILENO, msg, sizeof(msg)-1);

   abort ();
}

/*}}}*/

static void enable_keyboard_interrupt (void) /*{{{*/
{
   static int is_enabled = 0;

   if (is_enabled == 0)
     {
	(void) SLang_set_abort_signal (NULL);
	is_enabled = 1;
     }
}

/*}}}*/

void init_signals (void) /*{{{*/
{
   enable_keyboard_interrupt ();
   SLsignal (SIGPIPE, SIG_IGN);

   SLsignal (SIGSEGV, sig_abort_isis);
   SLsignal (SIGBUS, sig_abort_isis);
   SLsignal (SIGILL, sig_abort_isis);

   SLsignal (SIGHUP, sig_exit_isis);
   SLsignal (SIGTERM, sig_exit_isis);
}

/*}}}*/

/*}}}*/

/*{{{ 'source' input */

/* Note: In the following section the word 'source' is used in the csh
 * sense.  Sourcing a file means to evaluate it as if the user typed it in.
 */
typedef struct
{
   char *file;
   FILE *fp;
   char buf[1024];
   char *b;
}
Sourced_File_Type;

static void delete_sourced_file_type (Sourced_File_Type *sf) /*{{{*/
{
   if (sf == NULL)
     return;

   SLang_free_slstring (sf->file);     /* NULL ok */
   if (sf->fp != NULL)
     (void) fclose (sf->fp);
   SLang_free_slstring (sf->b);               /* NULL ok */
   SLfree ((char *) sf);
}

/*}}}*/

static char *read_source_line (SLang_Load_Type *lt) /*{{{*/
{
   Sourced_File_Type *sf;

   if (lt == NULL)
     return NULL;

   sf = (Sourced_File_Type *) lt->client_data;
   if (NULL == fgets (sf->buf, sizeof (sf->buf), sf->fp))
     return NULL;

   if (sf->b != NULL)
     {
        SLang_free_slstring (sf->b);
        sf->b = NULL;
     }

   if (1 == SLang_run_hooks ("_isis->isis_massage_input_hook", 1, sf->buf))
     {
        if (-1 == SLang_pop_slstring (&sf->b))
          return NULL;

        return sf->b;
     }

   return sf->buf;
}

/*}}}*/

static int source_file (char *file) /*{{{*/
{
   FILE *fp;
   Sourced_File_Type *sf;
   SLang_Load_Type *lt;
   int status;

   if (NULL == (fp = fopen (file, "r")))
     {
        SLang_verror (SL_OBJ_NOPEN, "source %s: unable to open file", file);
        return -1;
     }

   sf = (Sourced_File_Type *) SLmalloc (sizeof(Sourced_File_Type));
   if (sf == NULL)
     {
        (void) fclose (fp);
        return -1;
     }
   memset ((char *) sf, 0, sizeof (Sourced_File_Type));
   if (NULL == (sf->file = SLang_create_slstring (file)))
     {
        (void) fclose (fp);
        delete_sourced_file_type (sf);
        return -1;
     }
   sf->fp = fp;

   if (NULL == (lt = SLns_allocate_load_type ("<stdin>", Isis_Readline_Namespace_Name)))
     {
        delete_sourced_file_type (sf);
        return -1;
     }
   lt->read = read_source_line;
   lt->auto_declare_globals = 1;
   lt->client_data = sf;

   status = SLang_load_object (lt);
   if (SLang_get_error())
     {
        fprintf (stderr, "Error occured on line %d of %s\n",
                 lt->line_num, file);
     }
   SLdeallocate_load_type (lt);

   delete_sourced_file_type (sf);
   return status;
}

/*}}}*/

static void source_file_cmd (char *file) /*{{{*/
{
   (void) source_file (file);
}

/*}}}*/

/*}}}*/

/*{{{ readline input */

void input_from_stdin (void) /*{{{*/
{
   Input_From_Stdin = 1;
   isis_set_pager ("cat");
}

/*}}}*/

static void use_slang_readline (void) /*{{{*/
{
   Use_SLang_Readline = 1;
}

/*}}}*/

/* adapted from jdl_getkey() by John Davis */
unsigned int isis_getkey (void) /*{{{*/
{
   int reset;
   unsigned int ch;

   reset = (TTY_Inited == 0);

   if (SLANG_GETKEY_ERROR == (ch = SLang_getkey ()))
     {
        if (reset)
          reset_tty ();
        (void) fprintf (stderr, "isis_getkey: read failed\n");
        exit_isis (1);
     }

   if (reset)
     reset_tty ();
   return ch;
}

/*}}}*/

static int buf_has_text (char *buf) /*{{{*/
{
   if (buf == NULL)
     return 0;

   while (*buf == ' ' || *buf == '\t' || *buf == '\n')
     buf++;

   return (*buf != 0 && *buf != ';');
}

/*}}}*/

static int handle_fgets_error (char *line) /*{{{*/
{
#ifdef EINTR
   if (errno == EINTR)
     {
        (void) SLang_handle_interrupt ();
        line[0] = 0;
        return 0;
     }
#endif

   fflush (stdout);
   SLfree (line);

   return -1;
}

/*}}}*/

static char *read_line_stdin (SLang_RLine_Info_Type *rl, char *prompt, unsigned int *n) /*{{{*/
{
   char *line;
   unsigned int bufsize;

   (void) rl;
   *n = 0;

   if (isatty (fileno(stdin)))
     {
        fputs (prompt, stdout);
        fflush (stdout);
     }

   if (NULL == (line = SLmalloc (BUFSIZE*sizeof(char))))
     return NULL;
   bufsize = BUFSIZE*sizeof(char);

   if (NULL == fgets (line, bufsize, stdin))
     {
        if (-1 == handle_fgets_error (line))
          return NULL;
     }

   *n = strlen (line);
   return line;
}

/*}}}*/

static char *read_input_line (SLang_RLine_Info_Type *rline, char *prompt) /*{{{*/
{
   unsigned int n;
   char *line = NULL;

   if (Input_From_Stdin)
     {
        line = read_line_stdin (rline, prompt, &n);
     }
   else
     {
        SLtt_get_screen_size ();
        SLrline_set_display_width (rline, SLtt_Screen_Cols);

        init_tty ();

        if (Use_SLang_Readline == 0)
          {
#ifdef HAVE_GNU_READLINE
             line = readline (prompt);
#endif
          }
        else
          {
             line = SLrline_read_line (rline, prompt, &n);
          }

        reset_tty ();
     }

   return line;
}

/*}}}*/

static int save_input_line (SLang_RLine_Info_Type *rline, char *line) /*{{{*/
{
   if ((Input_From_Stdin == 0) && (0 != buf_has_text (line)))
     {
        if (Use_SLang_Readline == 0)
          {
#ifdef HAVE_GNU_READLINE
             (void) rline;
             add_history (line);
#endif
          }
        else SLrline_save_line (rline);
     }

   return 0;
}

/*}}}*/

/*{{{ readline support for debugger (borrowed from slsh) */

static SLang_RLine_Info_Type *Intrinsic_Rline_Info;
#if USE_SLANG_READLINE
static void close_intrinsic_readline (void) /*{{{*/
{
   if (Intrinsic_Rline_Info != NULL)
     {
	SLrline_close (Intrinsic_Rline_Info);
	Intrinsic_Rline_Info = NULL;
     }
}
/*}}}*/
#endif

static int readline_intrinsic_internal (SLang_RLine_Info_Type *rli, char *prompt) /*{{{*/
{
   char *line;

   if (rli == NULL)
     rli = Intrinsic_Rline_Info;

#if USE_SLANG_READLINE
   if ((rli == NULL) && Input_From_Stdin)
     {
	Intrinsic_Rline_Info = SLrline_open (SLtt_Screen_Cols, SL_RLINE_BLINK_MATCH);
	if (Intrinsic_Rline_Info == NULL)
	  return -1;
	(void) SLang_add_cleanup_function (close_intrinsic_readline);
	rli = Intrinsic_Rline_Info;
     }
#endif
   enable_keyboard_interrupt ();
   line = read_input_line (rli, prompt);
   (void) save_input_line (rli, line);
   (void) SLang_push_malloced_string (line);
   return 0;
}

/*}}}*/

static int Rline_Type_Id = 0;

static SLang_MMT_Type *pop_rli_type (SLang_RLine_Info_Type **rlip) /*{{{*/
{
   SLang_MMT_Type *mmt;

   if (NULL == (mmt = SLang_pop_mmt (Rline_Type_Id)))
     return NULL;
   if (NULL == (*rlip = (SLang_RLine_Info_Type *)SLang_object_from_mmt (mmt)))
     {
	SLang_free_mmt (mmt);
	return NULL;
     }
   return mmt;
}

/*}}}*/

static void readline_intrinsic (char *prompt) /*{{{*/
{
   SLang_RLine_Info_Type *rli = NULL;
   SLang_MMT_Type *mmt = NULL;

   if (SLang_Num_Function_Args == 2)
     {
	if (NULL == (mmt = pop_rli_type (&rli)))
	  return;
     }
   (void) readline_intrinsic_internal (rli, prompt);

   if (mmt != NULL)
     SLang_free_mmt (mmt);
}

/*}}}*/

static void new_slrline_intrinsic (char *name) /*{{{*/
{
   SLang_RLine_Info_Type *rli;
   SLang_MMT_Type *mmt;

   if (NULL == (rli = SLrline_open2 (name, SLtt_Screen_Cols, SL_RLINE_BLINK_MATCH)))
     return;

   if (NULL == (mmt = SLang_create_mmt (Rline_Type_Id, (VOID_STAR) rli)))
     {
	SLrline_close (rli);
	return;
     }

   if (-1 == SLang_push_mmt (mmt))
     SLang_free_mmt (mmt);
}

/*}}}*/

static int init_readline (char *appname) /*{{{*/
{
   static int inited = 0;

   (void) appname;

   if (inited)
     return 0;

#if USE_SLANG_READLINE
   if (Input_From_Stdin == 0)
     {
	inited = 1;
	return 0;
     }
   if (-1 == SLrline_init (appname, NULL, NULL))
     return -1;
#endif

   inited = 1;
   return 0;
}

/*}}}*/

static void init_readline_intrinsic (char *appname) /*{{{*/
{
   (void) init_readline (appname);
}

/*}}}*/

/*}}}*/

/*{{{ Prompt manipulation */
static int set_prompt (char *prompt) /*{{{*/
{
   char *s;

   if (prompt == NULL || *prompt == 0)
     return 0;

   if (NULL == (s = isis_make_string (prompt)))
     return -1;

   ISIS_FREE(isis_prompt);
   isis_prompt = s;

   return 0;
}

/*}}}*/

static void set_prompt_intrin (void) /*{{{*/
{
   char *prompt;

   if (SLang_Num_Function_Args == 0)
     {
        (void) set_prompt (isis_prompt_default);
     }
   else if (SLANG_NULL_TYPE == SLang_peek_at_stack())
     {
        (void) SLdo_pop();
        (void) set_prompt (isis_prompt_default);
     }
   else if (0 == SLang_pop_slstring (&prompt))
     {
        (void) set_prompt (prompt);
        SLang_free_slstring (prompt);
     }
}

/*}}}*/

/*}}}*/

static void get_prompt_intrin (void) /*{{{*/
{
   SLang_push_string (isis_prompt);
}

/*}}}*/

/* slang2:  returns a malloced value */
static char *get_input_line (SLang_Load_Type *x) /*{{{*/
{
   char *prompt;
   char *line;

   if (x->parse_level == 0)
     {
        prompt = isis_prompt;
        if (-1 == SLang_run_hooks ("_isis->take_input_hook", 0))
          return NULL;
     }
   else prompt = "      ";

   line = read_input_line (RLI, prompt);

   if (line == NULL)
     {
        int err = SLang_get_error ();
        if (err == 0)
          exit_isis (0);   /* EOF */
        if (err == SL_UserBreak_Error)
          return NULL;     /* Ctrl-c */
        (void) fprintf (stderr, "Error Occurred: %s\n", SLerr_strerror (err));
        return NULL;
     }

   (void) SLang_run_hooks ("_isis->isis_readline_hook", 1, line);

   (void) save_input_line (RLI, line);

   return line;
}

/*}}}*/

static char *read_using_readline (SLang_Load_Type *x) /*{{{*/
{
   char *s;
   static char *last_s;

   if (last_s != NULL)
     {
        SLfree (last_s);
        last_s = NULL;
     }

   check_fpe ();

   if (SLang_get_error())
     return NULL;

   SLKeyBoard_Quit = 0;

   s = get_input_line (x);

   if (s == NULL)
     return NULL;

   if ((x->parse_level == 0)
       && (1 == SLang_run_hooks ("_isis->isis_massage_input_hook", 1, s)))
     {
        SLfree (s);
        if (-1 == SLpop_string (&s))
          return NULL;
     }

   if (SLang_get_error())
     {
        SLfree (s);
        return NULL;
     }

   last_s = s;

   return s;
}

/*}}}*/

#ifdef HAVE_GNU_READLINE
static void gnu_rl_close_history (void) /*{{{*/
{
   if (Gnu_Rl_History_File == NULL)
     return;

   write_history (Gnu_Rl_History_File);
   SLfree (Gnu_Rl_History_File);
}

/*}}}*/
#endif

#ifdef HAVE_GNU_READLINE
static int gnu_rl_open_history (void) /*{{{*/
{
   char *history_file_name = getenv ("ISIS_HISTORY_FILE");

   if ((Gnu_Rl_History_File != NULL)
       || (history_file_name == NULL))
     return 0;

   Gnu_Rl_History_File = SLpath_dircat (NULL, history_file_name);
   if ((Gnu_Rl_History_File == NULL)
       || ((read_history (Gnu_Rl_History_File) != 0) && (errno != ENOENT)))
     {
        fprintf (stderr, "Couldn't open/read readline history file: %s\n",
                 Gnu_Rl_History_File ? Gnu_Rl_History_File : history_file_name);
        fprintf (stderr, "Is the ISIS_HISTORY_FILE environment variable set correctly?\n");
        return -1;
     }

   return 0;
}

/*}}}*/
#endif

#ifdef HAVE_GNU_READLINE
/* Return the next name which partially matches from the symbol list. */
static char *gnu_rl_symbol_generator (const char *text, int state) /*{{{*/
{
   char *name = NULL;

   SLang_start_arg_list ();
   SLang_push_string ((char *)text);
   SLang_push_integer (state);
   SLang_end_arg_list ();

   SLang_execute_function ("_isis->gnu_rl_completion_hook");
   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
        SLdo_pop();
        return NULL;
     }

   if (-1 == SLpop_string (&name))
     return NULL;

   return name;
}

/*}}}*/
#endif

#ifdef HAVE_GNU_READLINE
/* Attempt to complete on the contents of TEXT.  START and END bound the
   region of rl_line_buffer that contains the word to complete.  TEXT is
   the word to complete.  We can use the entire contents of rl_line_buffer
   in case we want to do some simple parsing.  Return the array of matches,
   or NULL if there aren't any. */
static char **gnu_rl_completion_function (const char *text, int start, int end) /*{{{*/
{
   char **matches = NULL;
   int num_quotes, do_symbol_completion;
   char *p;

   (void) start; (void) end;

   /* Do symbol completion when we've seen an even number of
    * double-quotes, otherwise do filename completion */
   num_quotes = 0;
   for (p = rl_line_buffer; *p != 0; p++)
     {
        if (*p == '"')
          num_quotes++;
     }
   do_symbol_completion = ((num_quotes/2)*2 == num_quotes);

   if (do_symbol_completion)
     matches = rl_completion_matches (text, gnu_rl_symbol_generator);

   return matches;
}

/*}}}*/
#endif

#ifdef HAVE_GNU_READLINE
static int gnu_rl_open (void) /*{{{*/
{
   /* Allow conditional parsing of the ~/.inputrc file. */
   rl_readline_name = "isis";
   rl_attempted_completion_function = gnu_rl_completion_function;

   if (-1 == gnu_rl_open_history ())
     return -1;

   return 0;
}

/*}}}*/
#endif

static int set_readline_namespace_name (char *namespace_name) /*{{{*/
{
   char *s;
   if (namespace_name == NULL)
     namespace_name = "";
   if (NULL == (s = isis_make_string (namespace_name)))
     return -1;
   ISIS_FREE (Isis_Readline_Namespace_Name);
   Isis_Readline_Namespace_Name = s;
   return 0;
}

/*}}}*/

int open_readline (char *namespace_name) /*{{{*/
{
   if (RLI != NULL)
     SLrline_close (RLI);

   if (-1 == set_readline_namespace_name (namespace_name))
     return -1;

   if (Use_SLang_Readline == 0)
     {
#ifdef HAVE_GNU_READLINE
        if (-1 == gnu_rl_open ())
          return -1;
#endif
     }
   else
     {
        if (-1 == SLrline_init ("ISIS", NULL, NULL))
          return -1;
        if (NULL == (RLI = SLrline_open2 ("isis", SLtt_Screen_Cols, SL_RLINE_BLINK_MATCH)))
          return -1;
     }

   if (*Isis_Readline_Namespace_Name == 0)
     Readline_Load_Object = SLallocate_load_type ("<stdin>");
   else
     Readline_Load_Object = SLns_allocate_load_type ("<stdin>", Isis_Readline_Namespace_Name);

   if (NULL == Readline_Load_Object)
     {
        close_readline ();
        fputs ("Couldn't open readline\n", stderr);
        return -1;
     }

   Readline_Load_Object->read = read_using_readline;
   Readline_Load_Object->auto_declare_globals = 1;

   return 0;
}

/*}}}*/

static void close_readline (void) /*{{{*/
{
   if (Use_SLang_Readline == 0)
     {
#ifdef HAVE_GNU_READLINE
        gnu_rl_close_history ();
#endif
     }
   else if (RLI != NULL)
     {
        SLrline_close (RLI);
     }

   if (Readline_Load_Object == NULL)
     return;

   SLdeallocate_load_type (Readline_Load_Object);
   Readline_Load_Object = NULL;
}

/*}}}*/

static void set_readline_method_intrinsic (char *method) /*{{{*/
{
   if (method == NULL)
     return;

   if (Readline_Load_Object)
     {
        fprintf (stderr, "*** no support for interactive changes -- edit .isisrc instead\n");
        return;
     }

   if (0 == strcmp (method, "stdin"))
     {
        input_from_stdin();
     }
   else if (0 == strcmp (method, "slang"))
     {
        use_slang_readline();
     }
#ifdef HAVE_GNU_READLINE
   else if (0 == strcmp (method, "gnu"))
     {
        /* for completeness */
     }
#endif
   else
     {
        fprintf (stderr, "readline method '%s' is not supported\n", method);
     }
}

/*}}}*/

/*}}}*/

/*{{{ initialize, clear_errors */

static void print_error_message (char *err) /*{{{*/
{
   fputs (err, stderr);
   fputs ("\r\n", stderr);
   fflush (stderr);
}

/*}}}*/

static int Err_Message_Count = 0;

void error_hook (char *err) /*{{{*/
{
   if (NULL == err)
     return;

   if (SL_USAGE_ERROR == SLang_get_error ())
     {
        if (Isis_Batch_Mode == 0)
          {
             if (0 == strncmp (err, "<stdin>:", 8))
               return;
          }
        print_error_message (err);
        return;
     }

   if ((SLang_Traceback == 0) && (Err_Message_Count > 0))
     return;

   Err_Message_Count++;

   if (Isis_Silent_Errors == 0)
     print_error_message (err);
}

/*}}}*/

#ifdef SIGFPE
static int Float_Exception;
static void floating_point_exception (int sig) /*{{{*/
{
   (void) sig;
   enable_fpe (0);
   SLsignal (SIGFPE, SIG_IGN);
   Float_Exception = 1;
}

/*}}}*/
#endif

static int enable_fpe (int on) /*{{{*/
{
#ifdef __linux__
   unsigned int cw = 0x137f;
   /* It appears that this is flawed.  I will wait for glibc's fenv
    * and gcc support for add-fwait.  Sigh.
    */
   on = 0;
   if (on)
     {
	cw = 0x1372;
	_FPU_SETCW (cw);
	(void) SLsignal (SIGFPE, floating_point_exception);
     }
   else _FPU_SETCW (cw);
#else
   (void) on;
#endif

   /* Force a floating point calculation to ensure that the control
    * word gets loaded now.  I discovered this by looking at the linux
    * kernel source code.  Another point for open software.
    */
   (void) SLmath_hypot (3.0, 4.0);
   return 0;
}

/*}}}*/

static void check_fpe (void) /*{{{*/
{
#ifdef SIGFPE
   if (Float_Exception)
     {
	fprintf (stderr, "*** Floating Point Exception Occurred ***\r\n");
	Float_Exception = 0;
	enable_fpe (1);
     }
#endif
}

/*}}}*/

static void clear_errors (void) /*{{{*/
{
#ifndef __linux__
   if (Float_Exception)
     {
        fprintf (stderr, "*** Floating point exception\n");
        Float_Exception = 0;
     }
   (void) SLsignal (SIGFPE, floating_point_exception);
#endif

   check_fpe ();

   /* Workaround temporary S-Lang problem:
    * Make sure stack underflow error messages appear
    * when _traceback=0
    */
   if ((SLang_Version < 20103)
       && SLang_get_error () && (SLang_Traceback == 0))
     {
        SLang_verror (0, "%s", SLerr_strerror (0));
     }

   if (SLang_get_error())
     SLang_restart (1);

   isis_set_errno (0);
   SLang_set_error (0);
   SLKeyBoard_Quit = 0;
   Err_Message_Count = 0;
}

/*}}}*/

/*}}}*/

void interpreter_loop (void) /*{{{*/
{
   SLang_Traceback = 0;
   (void) SLdefine_for_ifdef ("__ISIS_PROMPT__");
   (void) SLang_run_hooks ("isis_interactive_hook", 0, NULL);

   for (;;)
     {
        clear_errors ();
        (void) SLang_load_object (Readline_Load_Object);
     }
}

/*}}}*/

static int safe_system (char *cmd) /*{{{*/
{
   (void) cmd;
   (void) isis_secure_mode ();
   return -1;
}

/*}}}*/

static int safe_popen (char *cmd, char *mode) /*{{{*/
{
   (void) cmd;  (void) mode;
   (void) isis_secure_mode ();
   return -1;
}

/*}}}*/

/*{{{ SLang intrinsics */

static SLang_Intrin_Var_Type Misc_Intrin_Vars [] =
{
   MAKE_VARIABLE("Isis_Public_Namespace_Name", &Isis_Public_Namespace_Name, SLANG_STRING_TYPE, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_Intrin_Var_Type Private_Intrin_Vars [] =
{
   MAKE_VARIABLE("Isis_Preparse_Only", &Isis_Preparse_Only, SLANG_INT_TYPE, 1),
   MAKE_VARIABLE("Isis_Debug_Mode", &Isis_Debug_Mode, SLANG_INT_TYPE, 1),
   MAKE_VARIABLE("_isis_config_slang_libdir", &slang_libdir, SLANG_STRING_TYPE, 0),
   SLANG_END_INTRIN_VAR_TABLE
};

#define V SLANG_VOID_TYPE
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE

static SLang_Intrin_Fun_Type Private_Intrinsics[] =
{
   MAKE_INTRINSIC_I("_quit", _quit, V),
   MAKE_INTRINSIC("_reset", reset_isis, V, 0),
   MAKE_INTRINSIC_S("_source", source_file_cmd, V),
   SLANG_END_INTRIN_FUN_TABLE
};

static SLang_Intrin_Fun_Type Readline_Intrinsics[] =
{
   MAKE_INTRINSIC_S("set_readline_method", set_readline_method_intrinsic, V),
   MAKE_INTRINSIC_0("set_prompt", set_prompt_intrin, V),
   MAKE_INTRINSIC_0("get_prompt", get_prompt_intrin, V),
   MAKE_INTRINSIC_1("atexit", at_exit, VOID_TYPE, SLANG_REF_TYPE),
   /* readline intrinsic is used mainly by the slang debugger */
   MAKE_INTRINSIC_S("slsh_readline_init", init_readline_intrinsic, VOID_TYPE),
   MAKE_INTRINSIC_S("slsh_readline_new", new_slrline_intrinsic, VOID_TYPE),
   MAKE_INTRINSIC_S("slsh_readline", readline_intrinsic, VOID_TYPE),
   SLANG_END_INTRIN_FUN_TABLE
};

/* Use these definitions to disable
 *   o command-line shell escapes
 *   o S-Lang system(), popen()
 *
 * Perhaps more functionality should be disabled?
 */

static SLang_Intrin_Fun_Type Secure_Intrinsics [] =
{
   MAKE_INTRINSIC_1("system", safe_system, I, S),
   MAKE_INTRINSIC_2("popen", safe_popen, I, S, S),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef V
#undef I
#undef S

/*}}}*/

void deinit_readline_module (void) /*{{{*/
{
   close_readline ();
   fflush (stdout);
   ISIS_FREE (Isis_Readline_Namespace_Name);
   ISIS_FREE (isis_prompt);
}

/*}}}*/

static void destroy_rline (SLtype type, VOID_STAR f) /*{{{*/
{
   SLang_RLine_Info_Type *rli;
   (void) type;

   rli = (SLang_RLine_Info_Type *) f;
   if (rli != NULL)
     SLrline_close (rli);
}

/*}}}*/

static int register_rline_type (void) /*{{{*/
{
   SLang_Class_Type *cl;

   if (Rline_Type_Id != 0)
     return 0;

   if (NULL == (cl = SLclass_allocate_class ("Isis_RLine_Type")))
     return -1;

   if (-1 == SLclass_set_destroy_function (cl, destroy_rline))
     return -1;

   /* By registering as SLANG_VOID_TYPE, slang will dynamically allocate a
    * type.
    */
   if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE, sizeof (SLang_RLine_Info_Type*),
                                     SLANG_CLASS_TYPE_MMT))
     return -1;

   Rline_Type_Id = SLclass_get_class_id (cl);

   return 0;
}

/*}}}*/

SLANG_MODULE(readline);
int init_readline_module_ns (char *ns_name) /*{{{*/
{
   SLang_NameSpace_Type *ns;
   SLang_NameSpace_Type *pub_ns;

   if (NULL == (ns = SLns_create_namespace (ns_name)))
     return isis_trace_return(-1);

   if (NULL == (pub_ns = SLns_create_namespace (Isis_Public_Namespace_Name)))
     return isis_trace_return(-1);

   if (Isis_Secure_Mode)
     {
        /* redefine insecure functions before any slang code is interpreted */
        if (-1 == SLns_add_intrin_fun_table (NULL, Secure_Intrinsics, NULL))
          return isis_trace_return(-1);
     }

   if (-1 == register_rline_type ())
     return -1;

   if ((-1 == SLns_add_intrin_var_table (ns, Misc_Intrin_Vars, NULL))
       || (-1 == SLns_add_intrin_var_table (ns, Private_Intrin_Vars, NULL))
       || (-1 == SLns_add_intrin_fun_table (ns, Private_Intrinsics, NULL))
       || (-1 == SLns_add_intrin_fun_table (pub_ns, Readline_Intrinsics, NULL)))
     return isis_trace_return(-1);

   if (-1 == set_prompt (isis_prompt_default))
     return isis_trace_return(-1);

   return 0;
}

/*}}}*/

