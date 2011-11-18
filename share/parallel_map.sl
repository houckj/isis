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

private variable Slave_Task;
private variable Parallel_Map_Info;

private define parallel_map_slave (s)
{
   send_msg (s, SLAVE_READY);

   forever
     {
        variable objs = recv_objs (s);

        variable index = objs[0];
        if (index < 0)
          break;

        send_msg (s, SLAVE_RESULT);
        send_objs (s, index, (@Slave_Task)(__push_list (objs[1]);;
                                           __qualifiers));
     }

   return 0;
}

private define maybe_finished (slaves)
{
   variable s;
   foreach s (slaves)
     {
        if (s.status != SLAVE_READY)
          return;
     }

   foreach s (slaves)
     {
        send_objs (s, -1);
     }
}

private define make_arglist (index)
{
   variable x = Parallel_Map_Info;
   variable i, args = {},
     num_args = length(x.arg_lists);

   _for i (0, num_args-1, 1)
     {
        variable lst = x.arg_lists[i];
        list_append (args,
                     (length(lst) == 1) ? lst[0] : lst[index]);
     }
   return args;
}

private define send_next_task (slv)
{
   variable x = Parallel_Map_Info;

   if (x.next_task < length(x.arg_indices))
     {
        variable index = x.arg_indices[x.next_task];
        send_objs (slv, index, make_arglist (index));
        x.next_task++;
        slv.status = SLAVE_RUNNING;
     }
   else maybe_finished (x.slaves);
}

private define recv_slave_result (slv)
{
   % objs = {index [, result1, result2, ...]}
   variable objs = recv_objs (slv);

   slv.status = SLAVE_READY;

   if (length (objs) > 1)
     {
        variable x = Parallel_Map_Info;
        variable index = objs[0];
        variable names = get_struct_field_names (x.results);
        variable i, n = length(names);
        _for i (0, n-1, 1)
          {
             variable v = get_struct_field (x.results, names[i]);
             variable o = objs[i+1];
             if (x.return_types[i] == Array_Type
                 && typeof(o) != Array_Type)
               {
                  o = [o];
               }
             v[index] = o;
          }
     }
}

private define parallel_map_handler (s, msg)
{
   switch (msg.type)
     {
      case SLAVE_READY:
        send_next_task (s);
     }
     {
      case SLAVE_RESULT:
        recv_slave_result (s ;; __qualifiers);
        send_next_task (s);
     }
}

private define prepare_args (arg_lists)
{
   variable num_args = length(arg_lists);
   variable i, max_arg_len = 0,
     arg_len = Int_Type[num_args];

   _for i (0, num_args-1, 1)
     {
        if ((typeof(arg_lists[i]) != Array_Type)
            && (typeof(arg_lists[i]) != List_Type))
          {
             variable a = list_pop (arg_lists, i);
             list_insert (arg_lists, [a], i);
          }
        variable len = length(arg_lists[i]);
        arg_len[i] = len;
        if (len > max_arg_len)
          max_arg_len = len;
     }

   variable vec, scl = where (arg_len != max_arg_len, &vec);
   if (((length(vec) > 0) && (any(arg_len[vec] != max_arg_len)))
       || ((length(scl) > 0) && (any(arg_len[scl] != 1))))
     {
        throw UsageError, "parallel_map:  vector args must have equal length";
     }

   return max_arg_len;
}

private define result_struct (return_types, max_arg_len)
{
   variable num_return_values = length(return_types);
   if (num_return_values == 0)
     return NULL;
   variable field_names = "result"
     + array_map (String_Type, &string, [0:num_return_values-1]);

   variable s = struct_combine (field_names);

   variable i;
   _for i (0, num_return_values-1, 1)
     {
        set_struct_field (s, field_names[i],
                          @Array_Type(return_types[i], max_arg_len));
     }

   return s;
}

define parallel_map ()
{
   variable msg = "parallel_map ([Return_Types...,] &func, args, ... [; qualifiers]);";

   if (_NARGS < 2)
     {
        _pop_n (_NARGS);
        usage(msg);
     }

   variable return_types, arg_lists;

   variable in_args = __pop_list (_NARGS);
   variable i, num_in_args = length(in_args);
   _for i (0, num_in_args-1, 1)
     {
        if (typeof(in_args[i]) == Ref_Type)
          {
             return_types = (i > 0) ? in_args[[0:i-1]] : [Void_Type];
             Slave_Task = in_args[i];
             arg_lists = in_args[[i+1:]];
             break;
          }
     }

   variable max_arg_len = prepare_args (arg_lists);
   variable slaves = new_slave_list ( ;;__qualifiers);

   Parallel_Map_Info = struct
     {
        next_task = 0,
        arg_lists = arg_lists,
        arg_indices = [0:max_arg_len-1],
        return_types = return_types,
        results = result_struct (return_types, max_arg_len),
        slaves = slaves,
        qualifiers = __qualifiers
     };

   variable num_slaves =
     qualifier ("num_slaves", min ([_num_cpus(), max_arg_len]));

   loop (num_slaves)
     {
        variable s = fork_slave (&parallel_map_slave;; __qualifiers);
        s.status = SLAVE_READY;
        append_slave (slaves, s);
     }

   manage_slaves (slaves, &parallel_map_handler ;; __qualifiers);

   variable num_return_values = length(return_types);
   if ((num_return_values > 1)
       || (num_return_values == 1) && (return_types[0] != Void_Type))
     {
        () = _push_struct_field_values (Parallel_Map_Info.results);
        _stk_reverse (num_return_values);
     }
}
