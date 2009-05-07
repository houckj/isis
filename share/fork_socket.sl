% -*- mode: SLang; mode: fold -*-
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2008 Massachusetts Institute of Technology
%
%    Author:  John C. Houck <houck@space.mit.edu>
%
%    This software was developed by the MIT Center for Space Research under
%    contract SV1-61010 from the Smithsonian Institution.
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
#ifndef __HAVE_ISIS_EXTRAS__
#stop
#endif

$1 = path_dirname (__FILE__);
if (NULL != stat_file (path_concat ($1, "extras.sl")))
  require ("extras");
else
{
   $2 = "../extras";
   prepend_to_isis_load_path (path_concat ($1, $2));
   prepend_to_isis_module_path (path_concat ($1, $2));
}

require ("fork");
require ("socket");
require ("select");

#ifnexists socketpair
#stop
#endif

%  Public functions:
%     slaves = new_slave_list ( [; qualifiers]);
%     s = fork_slave (&task [,args [; qualifiers]])
%     append_slave (slaves, s);
%     manage_slaves (slaves, &message_handler [; qualifiers]);
%     send_msg (fp, type)
%     msg = recv_msg (fp)
%     array = read_n_array_vals (fp, n, type);
%     status = write_array (fp, array);
%     pid_vmessage (...)

private variable User_Message_Handler;
private variable Slaves_Running;
private variable Sigchld_Received;
private variable Verbose = 0;

variable
   SLAVE_EXITED  = -101,
   SLAVE_STARTED =  100,
   SLAVE_RUNNING =  101,
   SLAVE_READY   =  102,
   SLAVE_RESULT  =  103,
   SLAVE_EXITING =  104;

private define slave_is_active (s)
{
   return s.status > 0;
}

private define process_struct (pid, sock)
{
   variable s = struct {pid, sock, fp, status, data};
   s.pid = pid;
   s.sock = sock;
   s.status = SLAVE_STARTED;
   s.data = NULL;
   s.fp = fdopen (s.sock, "w+");
   if (s.fp == NULL)
     {
        vmessage ("%s: fdopen failed (%s)", _function_name, errno_string());
        return NULL;
     }
   return s;
}

private define sigchld_handler (sig);
private define sigchld_handler (sig)
{
   Sigchld_Received++;
   signal (SIGCHLD, &sigchld_handler);
}

define do_close (fd)
{
   variable status;
   do
     {
        status = close (fd);
     }
   while (status < 0 && errno == EINTR);

   return status;
}

define do_fflush (fp)
{
   variable status;
   do
     {
        status = fflush (fp);
     }
   while (status < 0 && errno == EINTR);

   return status;
}

private define cleanup_slave (s)
{
   Slaves_Running--;
   s.status = SLAVE_EXITED;

   if ((s.sock != NULL)
       && (do_close (s.sock) != 0))
     {
        variable msg = sprintf ("%s:  close failed (%s)", _function_name, errno_string());
        throw IOError, msg;
     }

   return 0;
}

private define find_slave (slaves, pid)
{
   variable s;
   foreach s (slaves)
     {
        if (s.pid == pid)
          return s;
     }

   throw ApplicationError, "slave pid=$pid not found!"$;
   %return NULL;
}

private define call_waitpid_for_slave (s)
{
   variable w;
   do
     {
        w = waitpid (s.pid, WNOHANG|WUNTRACED);
        if (w == NULL)
          {
             if (errno == ECHILD)
               {
                  () = cleanup_slave (s);
                  return;
               }
             continue;
          }
     }
   while (w.pid == 0);

   if (w.exited && Verbose)
     vmessage ("child pid=%d exited, status = %d", w.pid, w.exit_status);
   if (w.signal)
     vmessage ("child pid=%d killed by signal %d", w.pid, w.signal);

   if (w.exited || w.signal)
     {
        () = cleanup_slave (s);
     }
}

define read_n_array_vals (fp, n, type)
{
   variable array = type[0];
   while (n > 0)
     {
        variable darray;
        variable num_read = fread (&darray, type, n, fp);
        if (num_read <= 0)
          {
             if (errno == EINTR)
               continue;
             throw IOError, sprintf ("-%d- %s:  %s", getpid(), _function_name, errno_string());
          }
        array = [array, darray];
        n -= num_read;
     }
   if (n > 0)
     throw IOError, sprintf ("-%d- %s:  %s", getpid(), _function_name, errno_string());
   return array;
}

define write_array (fp, array)
{
   variable i = 0, n = length(array);
   while (n > 0)
     {
        variable num_written = fwrite (array[[i:]], fp);
        if (num_written <= 0)
          {
             if (errno == EINTR)
               continue;
             throw IOError, sprintf ("-%d- %s:  %s", getpid(), _function_name, errno_string());
          }
        n -= num_written;
        i += num_written;
     }
   if (n > 0)
     throw IOError, sprintf ("-%d- %s:  %s", getpid(), _function_name, errno_string());

   return do_fflush (fp);
}

define recv_msg (fp)
{
   if (feof(fp))
     return NULL;

   variable msg = read_n_array_vals (fp, 2, Int_Type);

   variable s = struct {from_pid, type};
   s.from_pid = msg[0];
   s.type = msg[1];

   return s;
}

define send_msg (fp, type)
{
   if (write_array (fp, [getpid(), type]))
     throw IOError;
}

private define handle_message (s, msg)
{
   switch (msg.type)
     {
      case SLAVE_EXITING:
        % handshake to avoid race condition -- we don't want the
        % slave to exit before we've finished reading its results
        send_msg (s.fp, SLAVE_EXITING);
        call_waitpid_for_slave (s);
     }
     {
        % default:
        if (User_Message_Handler == NULL)
          {
             print(msg);
             throw ApplicationError, "Unsupported message type";
          }
        (@User_Message_Handler)(s, msg);
     }

   % return non-zero if the number of slaves changed,
   % otherwise, return zero.

   return (msg.type == SLAVE_EXITING);
}

private define slave_sockets (slaves)
{
   variable s, fds = FD_Type[0], fps = File_Type[0];

   foreach s (slaves)
     {
        if (slave_is_active(s))
          {
             fds = [fds, s.sock];
             fps = [fps, s.fp];
          }
     }

   if (length(fds) == 0)
     return NULL;

   variable ss = struct {fd, fp};
   ss.fd = fds;
   ss.fp = fps;

   return ss;
}

private define pending_messages (slaves)
{
   variable socks = slave_sockets (slaves);
   if (socks == NULL)
     return NULL;

   variable ss = select (socks.fd, NULL, NULL, -1);
   if (ss == NULL)
     throw IOError, "pending_messages:  error on socket";

   if (ss.nready == 0)
     return NULL;

   return socks.fp[ss.iread];
}

private variable Do_Sigtest = 0;
private variable Num_Sigusr1 = 0;
private define sigusr1_handler (sig)
{
   Num_Sigusr1++;
}
private define catch_sigusr1()
{
   sigprocmask (SIG_BLOCK, SIGUSR1);
   signal (SIGUSR1, &sigusr1_handler);
   sigprocmask (SIG_UNBLOCK, SIGUSR1);
}

define fork_slave ()
{
   variable args = NULL;
   if (_NARGS > 1)
     args = __pop_args (_NARGS-1);

   variable func_ref = ();

   variable s1, s2, e;
   try (e)
     {
        (s1, s2) = socketpair (PF_UNIX, SOCK_STREAM, 0);
     }
   catch AnyError:
     {
        vmessage ("%s: socketpair failed (%s)", _function_name, errno_string());
        return NULL;
     }

   signal (SIGCHLD, &sigchld_handler);

   variable pid = fork ();
   if (pid < 0)
     {
        vmessage ("%s: fork failed (%s)", _function_name, errno_string());
        return NULL;
     }

   if (pid == 0)
     {
        % child
        signal (SIGINT, SIG_DFL);
        signal (SIGCHLD, SIG_DFL);
        if (Do_Sigtest) catch_sigusr1();
        variable status = 1;
        try (e)
          {
             if (do_close (s1) != 0)
               throw IOError, errno_string();
             variable p = process_struct (0, s2);
             if (args != NULL)
               status = (@func_ref)(p, __push_args(args));
             else
               status = (@func_ref)(p);
          }
        catch AnyError:
          {
             vmessage ("*** caught exception from slave process pid=%d", getpid());
             print(e);
          }
        finally
          {
             % handshake to avoid race condition -- we don't want the
             % slave to exit before we've finished reading its results.
             send_msg (p.fp, SLAVE_EXITING);
             () = recv_msg (p.fp);
             _exit (status);
          }
     }

   % parent
   Slaves_Running++;
   if (do_close (s2) != 0)
     throw IOError, errno_string();
   return process_struct (pid, s1);
}

private define shuffle (array)
{
   % Fisher-Yates-Durstenfeld in-place shuffle
   variable n = length(array);
   while (n > 1)
     {
        variable k = int(n * urand());
        n--;
        variable tmp = array[n];
        array[n] = array[k];
        array[k] = tmp;
     }
}

private define handle_pending_messages (slaves)
{
   variable fp_set = pending_messages (slaves);
   % The fp_set must be updated whenever the number of slaves
   % changes (e.g. when one or more slaves exit).

   variable fp, num_slaves_changed = 0;

   do
     {
        foreach fp (fp_set)
          {
             variable msg = recv_msg (fp);
             if (msg == NULL)
               continue;

             variable s = find_slave (slaves, msg.from_pid);

             if (handle_message (s, msg))
               num_slaves_changed = 1;
          }

        shuffle (fp_set);
     }
   while (num_slaves_changed == 0);
}

private define master_sigint_handler (sig)
{
   throw UserBreakError;
}

private define kill_slaves (slaves)
{
   signal (SIGCHLD, SIG_DFL);

   variable mask;
   try
     {
        sigprocmask (SIG_BLOCK, SIGINT, &mask);
        signal (SIGINT, SIG_IGN);
        sigprocmask (SIG_UNBLOCK, SIGINT);

        variable s;
        foreach s (slaves)
          {
             if (kill (s.pid, 0) == 0)
               {
                  if (kill (s.pid, SIGTERM) == 0)
                    call_waitpid_for_slave (s);
               }
          }
     }
   finally
     {
        sigprocmask (SIG_BLOCK, SIGINT);
        signal (SIGINT, &master_sigint_handler);
        sigprocmask (SIG_SETMASK, mask);
     }
}

define sigtest_slave (s, pids)
{
   variable n = length(pids), slaves_are_running;
   do
     {
        variable i;
        slaves_are_running = 0;
        _for i (0, n-1, 1)
          {
             variable pid = pids[i];
             if (pid < 0)
               continue;
             if (kill (pid, 0) != 0)
               {
                  pids[i] = -1;
                  continue;
               }
             slaves_are_running++;
             () = kill (pid, SIGUSR1);
          }
     }
   while (slaves_are_running);

   return 0;
}

define append_slave (slaves, s);

private define start_sigtest_slave (slaves)
{
   vmessage ("###### RUNNING SIGTEST ######");
   variable s, n=0, pids = Integer_Type[length(slaves)];
   foreach s (slaves)
     {
        pids[n] = s.pid;
        n++;
     }
   append_slave (slaves, fork_slave (&sigtest_slave, pids));
}

define manage_slaves (slaves, mesg_handler)
{
   Verbose = qualifier_exists ("verbose");
   User_Message_Handler = mesg_handler;

   if (length(slaves) != Slaves_Running)
     {
        throw ApplicationError,
          "manage_slaves:  slave list doesn't match number of running slaves";
     }

   if (Sigchld_Received)
     {
        kill_slaves (slaves);
        throw ApplicationError,
          "*** At least $Sigchld_Received slaves exited before manager could start!"$;
     }

   if (Do_Sigtest) start_sigtest_slave (slaves);

   variable e;
   try (e)
     {
        while (Slaves_Running)
          {
             handle_pending_messages (slaves);
          }
     }
   catch AnyError:
     {
        vmessage ("*** manage_slaves: caught exception!");
        print(e);
     }
   finally
     {
        kill_slaves (slaves);
        if (Do_Sigtest) () = list_pop (slaves, -1);
     }
}

define new_slave_list ()
{
   Sigchld_Received = 0;
   Slaves_Running = 0;
   Do_Sigtest = qualifier_exists ("sigtest");
   return list_new();
}

define append_slave (slaves, s)
{
   list_append (slaves, s);
}

define pid_vmessage ()
{
   variable args = __pop_args (_NARGS);
   () = fprintf (stderr, "%d: ", getpid());
   vmessage (__push_args(args));
}

public variable Isis_Num_Slaves;

require ("glob");
define guess_num_slaves ()
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

#ifdef FORK_SOCKET_TEST %{{{

private variable M = 2;

define task (s, which, num_loops)
{
   %pid_vmessage ("started");
   seed_random (_time - getpid());

   variable i;
   _for i (0, num_loops-1, 1)
     {
        %pid_vmessage ("ready");
        send_msg (s.fp, SLAVE_READY);
        variable x = read_n_array_vals (s.fp, 2, Double_Type);

        send_msg (s.fp, SLAVE_RESULT);

        if (write_array (s.fp, urand(M,M)))
          throw IOError, "*** slave: write failed";
     }

   %pid_vmessage ("finished");

   return 0;
}

define slave_is_ready (s)
{
   %vmessage ("slave %d is ready", s.pid);
   if (write_array (s.fp, urand(2)))
     throw IOError;
}

define slave_has_result (s)
{
   %vmessage ("slave %d has result", s.pid);
   s.data = read_n_array_vals (s.fp, M*M, Double_Type);
}

define message_handler (s, msg)
{
   switch (msg.type)
     {
      case SLAVE_READY:
        slave_is_ready (s);
     }
     {
      case SLAVE_RESULT:
        slave_has_result (s);
     }
     {
        % default:
        print(msg);
        throw ApplicationError, "Unsupported message type";
     }
}

define isis_main()
{
   seed_random (_time);

   variable i, s, slaves,
     num_slaves = 4,
     args = Int_Type[num_slaves];

   args[*] = 100;

   loop (5)
     {
        slaves = new_slave_list();
        _for i (0, num_slaves-1, 1)
          {
             s = fork_slave (&task, i, args[i]);
             append_slave (slaves, s);
          }

        manage_slaves (slaves, &message_handler; verbose);

        foreach s (slaves)
          {
             vmessage ("pid=%d sum_result=%g",
                       s.pid, sum(s.data));
          }
     }
}
%}}}
#endif
