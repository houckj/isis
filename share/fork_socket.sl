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
%     slaves = new_slave_list ();
%     s = fork_slave (&task [,args])
%     append_slave (slaves, s);
%     manage_slaves (slaves, &message_handler);
%     send_msg (fp, type)
%     msg = recv_msg (fp)
%     pid_vmessage (...)

private variable User_Message_Handler;
private variable Slaves_Running;
private variable Sigchld_Received = 0;
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

private define cleanup_slave (s)
{
   Sigchld_Received--;
   Slaves_Running--;

   s.status = SLAVE_EXITED;

   if ((s.sock != NULL)
       && (close (s.sock) != 0))
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

   return NULL;
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

define recv_msg (fp)
{
   if (feof(fp))
     return NULL;

   variable n, msg, mask;
   sigprocmask (SIG_BLOCK, SIGCHLD, &mask);
   n = fread (&msg, Int_Type, 2, fp);
   sigprocmask (SIG_SETMASK, mask);

   variable err_msg;
   if (n < 0)
     {
        err_msg = sprintf ("%s: read failed: (%s)", _function_name, errno_string());
        throw ApplicationError, err_msg;
     }
   else if (n != 2)
     {
        err_msg = sprintf ("%s: invalid message buffer (fread returned %d) %s",
                           _function_name, n, errno_string());
        throw ApplicationError, err_msg;
     }

   variable s = struct {from_pid, type};
   s.from_pid = msg[0];
   s.type = msg[1];

   return s;
}

define send_msg (fp, type)
{
   variable n, mask;
   sigprocmask (SIG_BLOCK, SIGCHLD, &mask);
   n = fwrite ([getpid(), type], fp);
   sigprocmask (SIG_SETMASK, mask);

   variable err_msg;
   if (n != 2)
     {
        err_msg = sprintf ("%s: fwrite failed %s", _function_name, errno_string());
        throw ApplicationError, err_msg;
     }
   if (0 != fflush (fp))
     {
        err_msg = sprintf ("%s: fflush failed %s", _function_name, errno_string());
        throw ApplicationError, err_msg;
     }
}

private define handle_message (s, msg)
{
   switch (msg.type)
     {
      case SLAVE_EXITING:
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

private define messages_pending (slaves)
{
   variable socks = slave_sockets (slaves);
   if (socks == NULL)
     return NULL;

   variable ss = select (socks.fd, NULL, NULL, -1);
   if (ss == NULL)
     throw ApplicationError, "messages_pending:  error on socket";

   if (ss.nready == 0)
     return NULL;

   return socks.fp[ss.iread];
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
        if (close (s1) != 0)
          throw IOError, errno_string();
        variable p = process_struct (0, s2);
        variable status;
        if (args != NULL)
          status = (@func_ref)(p, __push_args(args));
        else
          status = (@func_ref)(p);
        send_msg (p.fp, SLAVE_EXITING);
        exit (status);
     }

   % parent
   if (close (s2) != 0)
     throw IOError, errno_string();
   return process_struct (pid, s1);
}

private define random_indices (n)
{
   variable i, o = [0:n-1];

   _for i (0, n-1, 1)
     {
        if (urand() < 0.5)
          {
             variable tmp = o[n-i-1];
             o[n-i-1] = o[i];
             o[i] = tmp;
          }
     }

   return o;
}

private define round_robin (slaves, fp_set)
{
   % The fp_set must be updated whenever the number of slaves
   % changes (e.g. when one or more slaves exit).

   % visit the slaves in random order.
   variable
     o = random_indices (length(fp_set)),
     num_slaves_changed = 0;

   variable fp, s, msg;

   foreach fp (fp_set[o])
     {
        msg = recv_msg (fp);
        if (msg == NULL)
          continue;

        s = find_slave (slaves, msg.from_pid);

        if (handle_message (s, msg))
          num_slaves_changed = 1;
     }

   return num_slaves_changed;
}

private define sigchld_handler (sig);
private define sigchld_handler (sig)
{
   Sigchld_Received++;
   signal (SIGCHLD, &sigchld_handler);
}

define manage_slaves (slaves, mesg_handler)
{
   Verbose = qualifier_exists ("verbose");

   User_Message_Handler = mesg_handler;

   try
     {
        Sigchld_Received = 0;
        signal (SIGCHLD, &sigchld_handler);

        while (Slaves_Running)
          {
             variable fp_set = messages_pending (slaves);
             while (round_robin (slaves, fp_set) == 0)
               {}
          }
     }
   finally
     {
        signal (SIGCHLD, SIG_DFL);
     }
}

define new_slave_list ()
{
   Slaves_Running = 0;
   return list_new();
}

define append_slave (slaves, s)
{
   Slaves_Running++;
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

define task (s, num_loops)
{
   pid_vmessage ("started");
   seed_random (_time - getpid());

   variable i;
   _for i (0, num_loops-1, 1)
     {
        send_msg (s.fp, SLAVE_READY);
        variable x, n;
        n = fread (&x, Double_Type, 2, s.fp);
        if (n != 2)
          {
             vmessage ("task:  read failed");
             throw ApplicationError;
          }
        send_msg (s.fp, SLAVE_RESULT);

        n = fwrite (urand(100,100), s.fp);
        if (n != 10000)
          throw ApplicationError, "*** slave: fwrite failed";
     }

   return 0;
}

define slave_is_ready (s)
{
   %vmessage ("slave %d is ready", s.pid);
   variable n, mask;
   sigprocmask (SIG_BLOCK, SIGCHLD, &mask);
   n = fwrite (urand(2), s.fp);
   sigprocmask (SIG_SETMASK, mask);

   variable err_msg;
   if (n != 2)
     {
        err_msg = sprintf ("%s: fwrite failed %s",
                           _function_name, errno_string());
        throw ApplicationError, err_msg;
     }
   if (0 != fflush (s.fp))
     {
        err_msg = sprintf ("%s: fflush failed %s",
                           _function_name, errno_string());
        throw ApplicationError, err_msg;
     }
}

define slave_has_result (s)
{
   %vmessage ("slave %d has result", s.pid);
   variable n, x, mask;
   sigprocmask (SIG_BLOCK, SIGCHLD, &mask);
   n = fread (&x, Double_Type, 10000, s.fp);
   sigprocmask (SIG_SETMASK, mask);

   if (n != 10000)
     throw ApplicationError, "*** master: fread failed";

   s.data = x;
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

   slaves = new_slave_list();
   _for i (0, num_slaves-1, 1)
     {
        s = fork_slave (&task, args[i]);
        append_slave (slaves, s);
     }

   manage_slaves (slaves, &message_handler);

   foreach s (slaves)
     {
        vmessage ("pid=%d sum_result=%g",
                  s.pid, sum(s.data));
     }
}
%}}}
#endif
