% -*- mode: SLang; mode: fold -*-

#ifndef __HAVE_SLXPA_MODULE__
# ifndef __HAVE_ISIS_EXTRAS__
# stop
# endif
$1 = path_dirname (__FILE__);
if (NULL != stat_file (path_concat ($1, "extras.sl")))
  require ("extras");
else
{
   foreach (["../extras/slxpa/src", "../extras"])
     {
        $2 = ();
        prepend_to_isis_load_path (path_concat ($1, $2));
        prepend_to_isis_module_path (path_concat ($1, $2));
     }
}
#endif   % __HAVE_SLXPA_MODULE__

require ("xpa");
require ("fork");

private variable Server_Info,
  Server_Name = "isis_xpa_server";

private variable Master_Mask;
define init_slaves () %{{{
{
   % Avoid zombies:  to make sure we wait for all
   % child processes, block SIGINT while forking.
   sigprocmask (SIG_BLOCK, SIGINT, &Master_Mask);

   return list_new();
}

%}}}

define server_name () %{{{
{
   return Server_Name;
}

%}}}

define pid_vmessage () %{{{
{
   variable args = __pop_args (_NARGS);
   () = fprintf (stderr, "pid %d: ", getpid());
   vmessage (__push_args(args));
}

%}}}

define slave_lookup () %{{{
{
   variable msg = "Struct_Type = slave_lookup (pid);";
   variable pid;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   pid = ();

   variable s;
   foreach s (Server_Info)
     {
        if (s.pid == pid)
          return s;
     }

   return NULL;
}

%}}}

private define wait_for_server (n) %{{{
{
   variable name = server_name();

   % Note that this fails if called from the server process...
   loop (n)
     {
        if (xpaaccess (name))
          return 0;
        sleep (1);
     }

   return -1;
}

%}}}

private variable Slave_Counter = 0;
private define process_struct ()
{
   variable s = struct {pid, status, done, seq_num, slave_info};
   s.seq_num = Slave_Counter;   Slave_Counter++;
   s.pid = -1;
   s.status = -1;
   s.done = 0;
   s.slave_info = NULL;
   return s;
}

define fork_slave () %{{{
{
   variable msg = "Struct_Type = fork_slave (&task [, args])";

   if (_NARGS == 0)
     usage(msg);

   variable ref, args = NULL;
   if (_NARGS > 1)
     args = __pop_args (_NARGS - 1);

   ref = ();

   variable verbose = qualifier_exists ("verbose");
   variable s = process_struct ();

   s.pid = fork();

   if (s.pid < 0)
     exit(1);
   else if (s.pid > 0) % parent
     return s;

   % child (s.pid == 0)
   if (verbose)
     pid_vmessage ("child started");

   % slave should exit when master gets keyboard interrupt.
   sigprocmask (SIG_BLOCK, SIGINT);
   signal (SIGINT, SIG_DFL);
   sigprocmask (SIG_UNBLOCK, SIGINT);

   if (wait_for_server (10))
     {
        pid_vmessage ("Exiting -- failed to contact server %s"$,
                      server_name());
        s.status = -1;
        return;
     }

   (@ref)(__push_args(args));
   exit(s.status);
}

%}}}

private define perform_wait (slaves) %{{{
{
   variable verbose = qualifier_exists ("verbose");

   variable s, num_slaves = 0;
   foreach s (slaves)
     {
        if (s.done == 0)
          num_slaves++;
     }

   while (num_slaves > 0)
     {
        variable msec = 1;
        () = _XPAPoll (msec, 0);

        foreach s (slaves)
          {
             if (s.done)
               continue;

             variable w = waitpid (s.pid, WNOHANG | WUNTRACED);

             if (w == NULL)
               {
                  if (errno == ECHILD)
                    {
                       % The process nolonger exits-- something else collected it!
                       num_slaves--;
                       s.done = 1;
                       continue;
                    }

                  continue;
               }

             if (w.pid == 0)
               continue;

             if (w.exited && verbose)
               vmessage ("child pid=%d exited, status = %d", w.pid, w.exit_status);
             else if (w.signal)
               vmessage ("child pid=%d killed by signal %d", w.pid, w.signal);

             if (w.exited || w.signal)
               {
                  num_slaves--;
                  s.done = 1;
               }
          }
     }
}

%}}}

private define null_send (paramlist){return pack ("i",getpid());}
private define null_recv (paramlist,buf){return 0;}
private variable _XPAServer;

define start_server () %{{{
{
   variable msg = "status = start_server (&send, &recv [; qualifiers])";
   variable send, recv;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (send, recv) = ();

   if (send == NULL) send = &null_send;
   if (recv == NULL) recv = &null_recv;

   variable name = server_name();

   if (xpaaccess (name))
     return 0;

   variable s, dir = path_dirname (__FILE__);
   foreach s (["../bin", "../extras/xpa"])
     {
        variable xpans_dir = path_concat (dir, s);
        if (NULL != stat_file (path_concat (xpans_dir, "xpans")))
          {
             putenv("PATH=${xpans_dir}:${PATH}"$);
          }
     }

   % XPANew's return value is useless, but to hold the server
   % open, the value must be kept someplace for as long as the
   % server is needed.
   _XPAServer = XPANew (name, send, recv ;; __qualifiers);

   if (_XPAServer == NULL)
     return -1;

   return 0;
}

%}}}

private define master_sigint_handler (sig) %{{{
{
   throw UserBreakError;
}

%}}}

define wait_for_slaves () %{{{
{
   variable msg = "wait_for_slaves (List_Type slaves)";
   variable slaves;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   slaves = ();

   if (typeof(slaves) != List_Type)
     usage(msg);

   if (length(slaves) == 0)
     return;

   Server_Info = slaves;

   try
     {
        % ignore SIGINT while we wait for the slaves to exit
        sigprocmask (SIG_BLOCK, SIGINT);
        signal (SIGINT, SIG_IGN);
        sigprocmask (SIG_UNBLOCK, SIGINT);

        perform_wait (slaves ;; __qualifiers);
     }
   finally
     {
        % restore default signal handling
        sigprocmask (SIG_BLOCK, SIGINT);
        signal (SIGINT, &master_sigint_handler);
        sigprocmask (SIG_SETMASK, Master_Mask);
     }
}

%}}}

public variable Isis_Num_Slaves;
#iftrue
require ("glob");
define guess_num_slaves () %{{{
{
   if (__is_initialized (&Isis_Num_Slaves) && Isis_Num_Slaves >= 0)
     return Isis_Num_Slaves;

   % linux-specific guess...
   variable num_slaves = 0;
   variable dir = "/sys/devices/system/cpu";
   if (NULL != stat_file (dir))
     {
        num_slaves = length(glob("$dir/cpu?"$));
     }

   return num_slaves;
}

%}}}
#endif
