% -*- mode: SLang; mode: fold -*-
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2011 Massachusetts Institute of Technology
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

require ("fork");
require ("socket");
require ("select");

%  Public functions:
%     slaves = new_slave_list ( [; qualifiers]);
%     s = fork_slave (&task [,args [; qualifiers]]);
%     append_slave (slaves, s);
%     manage_slaves (slaves, &message_handler [; qualifiers]);
%     s = replace_slave (&s, &task [, args [; qualifiers]]);
%     send_msg (s, type)
%     msg = recv_msg (s)
%     array = read_array (fp, n, type);
%     status = write_array (fp, array);
%     pid_vmessage (...)

variable Isis_Slaves = struct
{
   num_slaves,            % number of slaves I can create
   num_slave_slaves = 0,  % number of slaves each slave can create
   nice = 0               % slaves run at this nice level
};

variable
   SLAVE_EXITED  = -101,
   SLAVE_EXITING = -100,
   SLAVE_STARTED =  100,
   SLAVE_RUNNING =  101,
   SLAVE_READY   =  102,
   SLAVE_RESULT  =  103;

private variable User_Message_Handler;
private variable Slaves_Running = 0;
private variable Slaves_Were_Replaced;
private variable Sigchld_Received = 0;
private variable Verbose = 0;
private variable Ctrl_C_Kills_Slaves;

private define slave_is_active (s)
{
   return s.status > 0;
}

private variable Check_For_Unread_Data = 0;
% If the master checks for unread data in the message handler
% loop, it can end up blocking in a read if replace_slave is
% being used.  If replace_slave is *not* being used, then
% it doesn't seem to matter whether Check_For_Unread_Data is zero
% or not.

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
   % With the default stdio buffering, deadlocks sometimes occur
   % when signals interrupt fflush (the sigtest deadlock disappears
   % when delivery of SIGUSR1 is blocked during the fflush call).
   % Using unbuffered streams seems to solve this problem.
   if (setvbuf (s.fp, _IONBF, 0) != 0)
     {
        vmessage ("%s: unable to turn off stream buffering! (%s)",
                  _function_name, errno_string());
     }
   return s;
}

private define sigchld_handler (sig);
private define sigchld_handler (sig)
{
   Sigchld_Received++;
   signal (SIGCHLD, &sigchld_handler);
}

private define cleanup_slave (s)
{
   Slaves_Running--;
   s.status = SLAVE_EXITED;

   if (s.fp != NULL)
     {
        if (fclose (s.fp) != 0)
          {
             throw IOError, sprintf ("%s:  fclose failed (%s)",
                                     _function_name, errno_string());
          }
     }
   else if (s.sock != NULL)
     {
        if (close (s.sock) != 0)
          {
             throw IOError, sprintf ("%s:  close failed (%s)",
                                     _function_name, errno_string());
          }
     }

   % Be sure slang knows these descriptors are no longer valid.
   s.sock = NULL;
   s.fp = NULL;

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

define read_array (fp, n, type)
{
   variable array;
   if (fread (&array, type, n, fp) < 0)
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
             if (errno == EINTR || errno == 0)
               continue;
             throw IOError, sprintf ("-%d- %s:  errno=%d (%s)",
                                     getpid(), _function_name,
                                     errno, errno_string());
          }
        n -= num_written;
        i += num_written;
     }
   if (n > 0)
     throw IOError, sprintf ("-%d- %s:  %s",
                             getpid(), _function_name, errno_string());

   if (-1 == fflush(fp))
     {
        throw IOError, sprintf ("-%d- %s: fflush failed: errno=%d (%s)",
                                getpid(), _function_name,
                                errno, errno_string());
     }

   return 0;
}

define _recv_msg (fp)
{
   if (feof(fp))
     return NULL;

   variable msg = read_array (fp, 2, Int_Type);

   variable s = struct {from_pid, type};
   s.from_pid = msg[0];
   s.type = msg[1];

   return s;
}

define recv_msg ()
{
   if (_NARGS != 1)
     usage ("Struct_Type = recv_msg (Struct_Type slave)");

   variable s = ();
   return _recv_msg (s.fp);
}

define _send_msg (fp, type)
{
   if (write_array (fp, [getpid(), type]))
     throw IOError;
}

define send_msg ()
{
   if (_NARGS != 2)
     usage ("send_msg (Struct_Type slave, Integer_Type msg_type)");

   variable s, type;
   (s, type) = ();
   _send_msg (s.fp, type);
}

private define handle_message (s, msg)
{
   switch (msg.type)
     {
      case SLAVE_EXITING:
        % handshake to avoid race condition -- we don't want the
        % slave to exit before we've finished reading its results
        _send_msg (s.fp, SLAVE_EXITING);
        call_waitpid_for_slave (s);
     }
     {
      case SLAVE_RUNNING:
        %vmessage ("%d is running", s.pid);
     }
     {
        % default:
        if (User_Message_Handler != NULL)
          (@User_Message_Handler)(s, msg ;; __qualifiers);
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

   if (Check_For_Unread_Data)
     {
        % handle any descriptors with unread data,
        variable i,
          status = array_map (Int_Type, &feof, socks.fp);
        i = where (status == 0);
        if (length(i) > 0)
          return socks.fp[i];
     }

   variable timeout = qualifier ("timeout", -1);

   variable ss = select (socks.fd, NULL, NULL, timeout);
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

private define renice (dp)
{
   if (dp == 0)
     return 0;
   variable p = getpriority (PRIO_PROCESS, 0);
   if (p == NULL)
     return -1;
   variable s = setpriority (PRIO_PROCESS, 0, p + dp);
   if (s == -1)
     return -1;
   return getpriority (PRIO_PROCESS, 0);
}

private variable Initial_Pid = getpid();
private define check_fork_permission (pid)
{
   % is this a child process?
   if (pid == Initial_Pid)
     return;

   % we're child process:  can we fork?
   if (Slaves_Running < Isis_Slaves.num_slave_slaves)
     return;

   vmessage ("\n*** Error: slave %d cannot fork (Isis_Slaves.num_slave_slaves = %d)\n",
             pid, Isis_Slaves.num_slave_slaves);
   throw ApplicationError;
}

define fork_slave ()
{
   variable args = NULL;
   if (_NARGS > 1)
     args = __pop_args (_NARGS-1);

   variable func_ref = ();

   Ctrl_C_Kills_Slaves = qualifier_exists ("ctrl_c_kills_slaves");

   variable parent_pid = getpid();
   check_fork_permission (parent_pid);

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
        ifnot (Ctrl_C_Kills_Slaves)
          {
             if (parent_pid == Initial_Pid)
               () = setpgid (0, 0);
             else () = setpgid (0, parent_pid);
          }
        signal (SIGINT, SIG_DFL);
        signal (SIGCHLD, SIG_DFL);
        if (Do_Sigtest) catch_sigusr1();
        variable status = 1;
        try (e)
          {
             if (close (s1) != 0)
               throw IOError, errno_string();
             variable p = process_struct (0, s2);
             () = renice (typecast (qualifier ("nice", Isis_Slaves.nice), Integer_Type));
             if (args != NULL)
               status = (@func_ref)(p, __push_args(args) ;; __qualifiers);
             else
               status = (@func_ref)(p ;; __qualifiers);
          }
        catch AnyError:
          {
             vmessage ("*** caught exception from slave process pid=%d", getpid());
             variable traceback = e.traceback;
             e.traceback = NULL;
             print(e);
             message (traceback);
          }
        finally
          {
             % handshake to avoid race condition -- we don't want the
             % slave to exit before we've finished reading its results.
             if (Do_Sigtest) vmessage ("%d: Num_Sigusr1 = %d", getpid(), Num_Sigusr1);
             _send_msg (p.fp, SLAVE_EXITING);
             () = _recv_msg (p.fp);
             () = close (s2);
             _exit (status);
          }
     }

   % parent
   Slaves_Running++;
   if (close (s2) != 0)
     throw IOError, errno_string();
   return process_struct (pid, s1);
}

private define handle_pending_messages (slaves)
{
   variable fp_set = pending_messages (slaves ;; __qualifiers);
   if (fp_set == NULL)
     return;

   variable i, num_fp = length(fp_set);

   Slaves_Were_Replaced = 0;

   _for i (0, num_fp-1, 1)
     {
        variable msg = _recv_msg (fp_set[i]);
        if (msg == NULL)
          continue;

        variable s = find_slave (slaves, msg.from_pid);
        () = handle_message (s, msg ;; __qualifiers);

        if (Slaves_Were_Replaced)
          break;
     }
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
        if (Check_For_Unread_Data)
          {
             % If the master checks for unread data before calling
             % select, it will end up reading from this slave because
             % feof will return 0 from this slave's descriptor even
             % if the slave has written nothing to it.
             % For this reason, this slave must send something once
             % in a while to keep that read from blocking.
             _send_msg (s.fp, SLAVE_RUNNING);
          }

        variable i;
        slaves_are_running = 0;
        loop (10000){
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

private variable Current_Slave_List;

private define handler_condition ()
{
   variable condition = qualifier ("while_", NULL);

   if (condition == NULL)
     return Slaves_Running;

   return Slaves_Running && (@condition)( ;; __qualifiers);
}

define manage_slaves (slaves, mesg_handler)
{
   Current_Slave_List = slaves;
   User_Message_Handler = mesg_handler;

   Verbose = qualifier_exists ("verbose");

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
        while (handler_condition( ;; __qualifiers))
          {
             handle_pending_messages (slaves ;; __qualifiers);
          }
     }
   catch AnyError:
     {
        vmessage ("*** manage_slaves: caught exception!");
        variable traceback = e.traceback;
        e.traceback = NULL;
        print(e);
        message (traceback);
     }
   finally
     {
        ifnot (qualifier_exists ("while_"))
          kill_slaves (slaves);
        if (Do_Sigtest) () = list_pop (slaves, -1);
     }
}

define new_slave_list ()
{
   Sigchld_Received = 0;
   Slaves_Running = 0;
   Slaves_Were_Replaced = 0;
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

define replace_slave ()
{
   variable s, task_ref, args = NULL;
   variable msg = "Struct_Type = replace_slave (s, &task [, args] [; qualifiers])";

   if (_NARGS == 0)
     usage (msg);
   else if (_NARGS > 2)
     args = __pop_list (_NARGS-2);

   (s, task_ref) = ();

   variable i,
     num_slaves = length(Current_Slave_List);

   _for i (0, num_slaves-1, 1)
     {
        if (Current_Slave_List[i].pid == s.pid)
          break;
     }
   if (i == num_slaves)
     {
        vmessage ("*** %s:  unregistered slave:", _function_name);
        print(s);
        return NULL;
     }

   if (kill (s.pid, 0) != 0)
     {
        vmessage ("*** %s: slave pid=%d does not exist", _function_name, s.pid);
        return NULL;
     }

   % check for a process group in case the slave forked.
   if (kill (-s.pid, 0) == 0)
     {
        if (killpg (s.pid, SIGTERM) != 0)
          {
             vmessage ("*** %s: error sending TERM to pid=%d (%s)",
                       _function_name, s.pid, errno_string());
             return NULL;
          }
     }
   else if (kill (s.pid, SIGTERM) != 0)
     {
        vmessage ("*** %s: error sending TERM to pid=%d (%s)",
                  _function_name, s.pid, errno_string());
        return NULL;
     }

   % Call waitpid here while we know who to wait for.
   call_waitpid_for_slave (s);
   Current_Slave_List[i] = NULL;

   variable new_s;
   if (args != NULL)
     new_s = fork_slave (task_ref, __push_list (args) ;; __qualifiers);
   else
     new_s = fork_slave (task_ref ;; __qualifiers);

   Current_Slave_List[i] = new_s;
   Slaves_Were_Replaced++;

   return new_s;
}

require ("sysconf");
require ("glob");
define _num_cpus ()
{
   if (Isis_Slaves.num_slaves != NULL && Isis_Slaves.num_slaves >= 0)
     return Isis_Slaves.num_slaves;

   variable num_slaves = sysconf ("_SC_NPROCESSORS_ONLN");
   if (num_slaves != NULL)
     return num_slaves;

   % Sigh. Try a linux-specific guess.
   variable dir = "/sys/devices/system/cpu";
   if (NULL != stat_file (dir))
     return length(glob("$dir/cpu?"$));

   % We have at least 1 CPU
   return 1;
}

% Communications interface:
%
%   send_objs (s, obj [, obj, ...]);
%   List_Type = recv_objs (s [, num]);
%   buf = __fs_initsend();
%   __fs_pack (buf, item [, item, ..]);
%   x = __fs_unpack (buf [, num]);
%   status = __fs_send_buffer (fp, buf);
%   buf = __fs_recv_buffer (fp);
%
%  TODO:
%    - add support for:
%        array of List_Type,
%        array of Assoc_Type,
%        multi-dimensional arrays of "non-basic" types
%        Ref_Types (for linked lists, etc),

define __fs_initsend ()
{
   return list_new();
}

private define do_pack (buf, item)
{
   if (typeof(item) == Array_Type)
     {
        variable t = _typeof(item);
        if (t == List_Type || t == Assoc_Type)
          {
             foreach (item)
               {
                  variable x = ();
                  list_append (buf, x);
               }
             return;
          }
        % fall-through
     }
   list_append (buf, item);
}

define __fs_pack ()
{
   variable msg = "__fs_pack (buf, item [, item, ...])";

   if (_NARGS < 2) usage(msg);

   variable item,
     args = __pop_args (_NARGS-1),
     buf = ();

   foreach item (args)
     {
        do_pack (buf, item.value);
     }
}

define __fs_unpack ()
{
   variable msg = "x = __fs_unpack (buf [, num])";
   variable buf, num = 1;

   switch (_NARGS)
     {case 1: buf = ();}
     {case 2: (buf, num) = ();}
     {
        usage (msg);
     }

   ifnot (__is_initialized (&buf))
     return NULL;

   if (num < 0)
     {
        num = length(buf);
     }

   if (num == 1)
     return list_pop (buf);

   variable x = list_new();

   while (num > 0)
     {
        list_append (x, list_pop (buf));
        num--;
     }

   return x;
}

private variable Types =
{
   Char_Type, UChar_Type,
   Integer_Type, UInteger_Type, Long_Type, ULong_Type,
   Float_Type, Double_Type,
   String_Type, BString_Type,
   Null_Type, Struct_Type, Assoc_Type, List_Type
};

private define datatype (i)
{
   return Types[i];
}

private define datatype_index (object)
{
   variable i, n = length(Types);

   _for i (0, n-1, 1)
     {
        if (_typeof (object) == Types[i])
          return i;
     }

   return NULL;
}

private define recv_basic (fp, type)
{
   variable num_dims, dims;
   num_dims = read_array (fp, 1, Integer_Type)[0];
   dims = read_array (fp, num_dims, Integer_Type);

   variable len, array;
   len = int(prod(dims));
   array = read_array (fp, len, type);
   reshape (array, dims);

   return array;
}

private define send_basic (fp, x)
{
   variable dims = array_shape(x);
   if (write_array (fp, length(dims)) < 0)
     return -1;
   if (write_array (fp, dims) < 0)
     return -1;

   reshape (x, length(x));
   variable status = write_array (fp, x);
   reshape (x, dims);

   return status;
}

private define recv_string (fp)
{
   variable num = read_array (fp, 1, Integer_Type)[0];

   variable s = String_Type[num],
     len = read_array (fp, num, Integer_Type);

   variable i;

   _for i (0, num-1, 1)
     {
        variable u = read_array (fp, len[i], UChar_Type);
        s[i] = typecast (array_to_bstring(u) , String_Type);
     }

   if (num == 1)
     return s[0];

   return s;
}

private define write_string (fp, s)
{
   variable n = strbytelen(s) + 1;
   variable b = pack ("s$n"$, s);
   return write_array (fp, bstring_to_array (b));
}

private define send_string (fp, s)
{
   variable num = length(s);

   if (write_array (fp, num) < 0)
     return -1;

   variable i, len = Integer_Type[num];

   if (typeof(s) == String_Type)
     {
        len[0] = strbytelen(s) + 1;
        if (write_array (fp, len) < 0)
          return -1;
        return write_string (fp, s);
     }

   _for i (0, num-1, 1)
     {
        len[i] = strbytelen(s[i]) + 1;
     }

   if (write_array (fp, len) < 0)
     return -1;

   _for i (0, num-1, 1)
     {
        if (write_string (fp, s[i])<0)
          return -1;
     }

   return 0;
}

private define recv_null (fp)
{
   variable num = read_array (fp, 1, Integer_Type)[0];
   if (num == 1)
     return NULL;
   return Null_Type[num];
}

private define send_null (fp, s)
{
   variable num = length(s);
   if (write_array (fp, num) < 0)
     return -1;
   return 0;
}

private define recv_struct();
private define recv_assoc();
private define recv_list();

private define recv_item (fp)
{
   variable type_index, type, item;

   type_index = read_array (fp, 1, Integer_Type)[0];
   type = datatype(type_index);

   switch (type)
     {case Struct_Type:  item = recv_struct (fp);}
     {case Assoc_Type:   item = recv_assoc (fp);}
     {case List_Type:    item = recv_list (fp);}
     {case String_Type:  item = recv_string (fp);}
     {case Null_Type:    item = recv_null (fp);}
     {
      case Array_Type:
        switch (_typeof(item))
          {case Struct_Type:  item = recv_struct (fp);}
          {case String_Type:  item = recv_string (fp);}
          {case Null_Type:    item = recv_null (fp);}
          {
             item = recv_basic (fp);
          }
     }
     {
        % default
        item = recv_basic (fp, type);
     }

   return item;
}

private define send_struct();
private define send_assoc();
private define send_list();

private define send_item (fp, item)
{
   variable type_index = datatype_index (item);

   if (write_array (fp, type_index) < 0)
     return -1;

   variable status;

   switch (typeof(item))
     {case Struct_Type:  status = send_struct (fp, item);}
     {case Assoc_Type:   status = send_assoc (fp, item);}
     {case List_Type:    status = send_list (fp, item);}
     {case String_Type:  status = send_string (fp, item);}
     {case Null_Type:    status = send_null (fp, item);}
     {
      case Array_Type:
        switch (_typeof(item))
          {case Struct_Type:   status = send_struct (fp, item);}
          {case String_Type:   status = send_string (fp, item);}
          {case Null_Type:     status = send_null (fp, item);}
          {
             status = send_basic (fp, item);
          }
     }
     {
        % default
        status = send_basic (fp, item);
     }

   return status;
}

% private define array_to_struct (fields)
% {
%    eval (sprintf ("define __atos__(){return struct {%s};}",
%                   strjoin (fields, ",")));
%    return eval ("__atos__()");
% }

private define recv_one_struct (fp)
{
   variable names = recv_string (fp);
   variable s = @Struct_Type(names);
   %variable s = array_to_struct (names);

   variable i, n = length(names);
   _for i (0, n-1, 1)
     {
        variable value = recv_item (fp);
        set_struct_field (s, names[i], value);
     }

   return s;
}

private define send_one_struct (fp, s)
{
   variable names = get_struct_field_names (s);

   if (-1 == send_string (fp, names))
     return -1;

   variable i, n = length(names);
   _for i (0, n-1, 1)
     {
        variable v = get_struct_field (s, names[i]);
        if (-1 == send_item (fp, v))
          return -1;
     }

   return 0;
}

private define recv_struct (fp)
{
   variable num = read_array (fp, 1, Integer_Type)[0];
   variable s = Struct_Type[num];

   if (num == 0)
     return s;

   variable i;
   _for i (0, num-1, 1)
     {
        s[i] = recv_one_struct (fp);
     }

   if (length(s) == 1)
     s = s[0];

   return s;
}

private define send_struct (fp, array)
{
   variable num = length (array);

   if (write_array (fp, num) < 0)
     return -1;

   variable i;
   _for i (0, num-1, 1)
     {
        if (send_one_struct (fp, array[i]) < 0)
          return -1;
     }

   return 0;
}

private define recv_assoc (fp)
{
   variable keys = recv_string (fp);
   variable a = Assoc_Type[];

   variable i, n = length(keys);
   _for i (0, n-1, 1)
     {
        variable value = recv_item (fp);
        a[keys[i]] = value;
     }

   return a;
}

private define send_assoc (fp, a)
{
   variable keys = assoc_get_keys (a);

   if (-1 == send_string (fp, keys))
     return -1;

   variable k;

   foreach k (keys)
     {
        if (-1 == send_item (fp, a[k]))
          return -1;
     }

   return 0;
}

private define recv_list (fp)
{
   variable num = read_array (fp, 1, Integer_Type)[0];

   variable x = {};

   loop (num)
     list_append (x, recv_item (fp));

   return x;
}

private define send_list (fp, list)
{
   if (write_array (fp, length(list)) < 0)
     return -1;

   variable x;
   foreach x (list)
     {
        if (-1 == send_item (fp, x))
          return -1;
     }

   return 0;
}

private variable
  BUFFER_START = 201,
  BUFFER_END   = 202;

define __fs_recv_buffer ()
{
   variable umsg = "List_Type = __fs_recv_buffer (File_Type)";

   if (_NARGS != 1)
     usage(umsg);

   variable fp = ();

   variable msg = _recv_msg (fp);
   if (msg.type != BUFFER_START)
     throw ApplicationError, "expected BUFFER_START message";

   variable num_objects = read_array (fp, 1, Integer_Type)[0];
   variable buf = list_new();

   loop (num_objects)
     {
        variable item = recv_item (fp);
        list_append (buf, item);
     }

   msg = _recv_msg (fp);
   if (msg.type != BUFFER_END)
     throw ApplicationError, "expected BUFFER_END message";

   return buf;
}

define __fs_send_buffer ()
{
   variable umsg = "status = __fs_send_buffer (File_Type, buf)";

   if (_NARGS != 2)
     usage(umsg);

   variable fp, buf;
   (fp, buf) = ();

   if (buf == NULL)
     throw ApplicationError, "invalid buffer";

   _send_msg (fp, BUFFER_START);
   if (write_array (fp, length(buf)) < 0)
     {
        vmessage ("failed writing buffer length");
        return -1;
     }

   variable item;
   foreach item (buf)
     {
        if (send_item (fp, item) < 0)
          {
             vmessage ("failed writing item of type %S", typeof(item));
             return -1;
          }
     }

   _send_msg (fp, BUFFER_END);

   return 0;
}

define send_objs ()
{
   variable msg = "send_objs (Struct_Type slv, obj [, obj, ...])";
   variable s, args = NULL;

   if (_NARGS < 2)
     usage (msg);

   args = __pop_list (_NARGS-1);
   s = ();

   variable x, buf = __fs_initsend ();
   foreach x (args)
     {
        __fs_pack (buf, x);
     }

   if (__fs_send_buffer (s.fp, buf))
     throw IOError;
}

define recv_objs ()
{
   variable msg = "List_Type = recv_objs (Struct_Type slv [, num])";
   variable s, num = -1;

   if (_NARGS > 1) num = ();
   s = ();

   variable buf = __fs_recv_buffer (s.fp);
   return __fs_unpack (buf, num);
}

#ifdef FORK_SOCKET_TEST %{{{

#ifndef urand
require ("rand");
private define urand ()
{
   if (_NARGS == 0)
     return rand_uniform ();
   variable a;
   if (_NARGS == 1)
     {
	a = ();
	return rand_uniform (a);
     }
   a = __pop_args (_NARGS);
   return rand_uniform (Char_Type[__push_args(a)]);
}
private define seed_random ()
{
   srand();
}
#endif

private variable M = 2;

define task (s, which, num_loops)
{
   %pid_vmessage ("started");
   seed_random (_time - getpid());

   variable i;
   _for i (0, num_loops-1, 1)
     {
        %pid_vmessage ("ready");
        _send_msg (s.fp, SLAVE_READY);
        variable x = read_array (s.fp, 2, Double_Type);

        _send_msg (s.fp, SLAVE_RESULT);

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
   s.data = read_array (s.fp, M*M, Double_Type);
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

define slsh_main()
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
