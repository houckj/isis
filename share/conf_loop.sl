% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2008 Massachusetts Institute of Technology
%
%    This software was developed by the MIT Center for Space Research under
%    contract SV1-61010 from the Smithsonian Institution.
%
%    Author:  John C. Houck  <houck@space.mit.edu>
%             Mike Nowak <mnowak@space.mit.edu>
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

public define save_par_limits()
{
   variable msg = "save_par_limits (indices, pmin, pmax, file)";
   variable indices, pmin, pmax, file;

   if (_NARGS != 4) usage (msg);
   (indices, pmin, pmax, file) = ();

   if (any (length (indices) != [length(pmin), length(pmax)]))
     throw UsageError, "*** Inconsistent array lengths";

   if (typeof(file) != String_Type)
     throw UsageError, "*** Must input a String_Type file name";

   variable current_params = get_params (indices);
   array_map (Void_Type, &set_par,
              indices, get_par(indices), 0, pmin, pmax);
   save_par (file);
   set_params (current_params);

   return;
}

private define valid_param_indices (index_array)
{
   variable s, pi = get_par_info (index_array),
     indices = list_new();

   foreach s (pi)
     {
        if (s.freeze != 0 || s.tie != NULL)
          throw UsageError, sprintf ("*** %s is not a free parameter", s.name);
        list_append (indices, s.index);
     }

   return indices;
}

private define try_conf (index, ctrl)
{
   variable save = qualifier_exists ("save");

   variable prefix = qualifier ("prefix", "conf_loop"),
     file_conf = prefix + ".conf",
     file_parm = prefix + ".par";

   if (save && not ctrl.serial)
     {
        variable pid = sprintf (".%d", getpid());
        file_conf += pid;
        file_parm += pid;
     }

   variable fp;
   if (save)
     {
        save_par (file_parm);
        if (NULL == stat_file (file_conf))
          {
             fp = fopen (file_conf, "w");
             () = fprintf (fp, "\n%s\n\n", get_fit_fun());
             () = fprintf (fp, "TRIES   PAR#      P_VALUE        P_MIN        P_MAX  TIME \n");
          }
        else
          {
             fp = fopen (file_conf, "a");
          }
     }

   variable s = struct {pmin, pmax, num_retries = 0};

   do
     {
        (s.pmin, s.pmax) = conf (index, ctrl.level, ctrl.tol);
        if (s.pmin == s.pmax)
          {
             % conf failed -- refit
             variable info;
             () = fit_counts (&info);
             s.num_retries++;

             if (save)
               {
                  save_par (file_parm);
                  () = fprintf (fp, "\n chi^2 =%11.4e,  DoF =%6i, time = %s \n \n",
                                info.statistic,
                                (info.num_bins-info.num_variable_params),
                                time());
               }
          }

     } while (s.pmin == s.pmax and s.num_retries <= ctrl.max_param_retries);

   if (save != 0 && s.pmin != s.pmax)
     {
        variable pi = get_par_info (index);
        () = fprintf (fp, "%5i   %4i  %11.4e  %11.4e  %11.4e  %s\n",
                      1 + s.num_retries, pi.index, pi.value,
                      s.pmin, s.pmax, time());
        if (fclose (fp))
          throw IOError, "failed closing " + file_conf;
     }

   return s;
}

private define conf_slave (s, ctrl)
{
   send_msg (s, SLAVE_READY);

   forever
     {
        variable objs = recv_objs (s);

        variable
          index = objs[0],
          param_version = objs[1],
          params = objs[2];

        if (index < 0)
          break;

        set_params (params);
        variable result = try_conf (index, ctrl ;; __qualifiers);

        variable info;
        () = eval_counts (&info);

        send_msg (s, SLAVE_RESULT);
        send_objs (s, index, param_version,
                   result, info, get_params());
     }

   return 0;
}

private variable Parallel_Conf_Loop_Info;

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
        send_objs (s, -1, NULL, NULL);
     }
}

private define send_next_task (slv)
{
   variable x = Parallel_Conf_Loop_Info;

   if (x.next_task < length(x.indices))
     {
        variable index = x.indices[x.next_task];
        send_objs (slv, index, x.param_version, get_params());
        x.next_task++;
        slv.status = SLAVE_RUNNING;
     }
   else maybe_finished (x.slaves);
}

private define recv_slave_result (slv)
{
   variable x = Parallel_Conf_Loop_Info;

   variable objs = recv_objs (slv);

   variable
     index = objs[0],
     param_version = objs[1],
     result = objs[2],
     info = objs[3],
     slave_params = objs[4];

   x.num_retries += result.num_retries;
   slv.status = SLAVE_READY;

   variable keep, restart;

   restart = (info.statistic < x.best_stat
              || (result.pmin == result.pmax
                  && result.num_retries < x.max_param_retries));

   keep = (restart ||
           (param_version == x.param_version
            && result.pmin != result.pmax));

   ifnot (keep || restart)
     return 0;

   if (keep)
     {
        % indices[] may have been shuffled
        variable k = where (x.ordered_indices == index)[0];
        x.pmin_final[k] = result.pmin;
        x.pmax_final[k] = result.pmax;
     }

   % If an index failed to converge, it's likely to
   % happen again -- move it to the front of the list
   % so that failures happen as early as possible.
   if (result.num_retries != 0 || result.pmin == result.pmax)
     {
        variable i = where (list_to_array (x.indices) == index)[0];
        list_insert (x.indices, list_pop (x.indices, i), 0);
        x.next_task = 0;
     }

   if (restart)
     {
        x.best_stat = info.statistic;
        set_params (slave_params);
        x.param_version++;
        x.next_task = 0;
     }

   return restart;
}

private define restart_slaves (slaves)
{
   variable s;
   foreach s (slaves)
     {
        if (s.status == SLAVE_READY)
          send_next_task (s);
     }
}

private define conf_handler (s, msg)
{
   variable x = Parallel_Conf_Loop_Info;

   switch (msg.type)
     {
      case SLAVE_READY:
        send_next_task (s);
     }
     {
      case SLAVE_RESULT:
        if (recv_slave_result (s))
          restart_slaves (x.slaves);
        else send_next_task (s);
     }
}

public define conf_loop()
{
   variable ctrl = struct
     {
        level = 1,
        tol = 1.e-3,
        max_param_retries = qualifier ("max_param_retries", 10),
        serial
     };
   variable indices;

   switch(_NARGS)
     {case 1: indices = ();}
     {case 2: (indices, ctrl.level) = ();}
     {case 3: (indices, ctrl.level, ctrl.tol) = ();}
     {
        usage (" (pmin,pmax) = conf_loop(indices[], [,level [,tolerance]] ;qualifiers);");
     }

   indices = valid_param_indices (indices);
   variable ordered_indices = list_to_array (indices);

   variable num_indices = length(indices),
     pmin_final = Double_Type[num_indices],
     pmax_final = Double_Type[num_indices];

   variable i, s, num_retries = 0,
     max_num_retries = qualifier ("max_num_retries", 100);

   variable save = qualifier_exists ("save"),
     prefix = qualifier ("prefix", "conf_loop");

   if (save)
     {
        variable dir = path_dirname (prefix);
        if (dir != ".")
          {
             if (NULL == stat_file (dir))
               {
                  if (mkdir (dir, 0777) != 0)
                    throw IOError, "failed opening directory $dir"$;
               }
          }
     }

   variable slaves,
     num_slaves = qualifier ("num_slaves", _num_cpus());
   ctrl.serial = qualifier_exists ("serial") || num_slaves < 2;

   ifnot (ctrl.serial)
     {
        slaves = new_slave_list ( ;;__qualifiers);

        Parallel_Conf_Loop_Info = struct
          {
             next_task = 0,
             param_version = 0,
             indices = indices,
             ordered_indices = ordered_indices,
             num_retries = num_retries,
             max_num_retries = max_num_retries,
             pmin_final = pmin_final,
             pmax_final = pmax_final,
             best_stat = _Inf,
             slaves = slaves
          };

        loop (num_slaves)
          {
             s = fork_slave (&conf_slave, ctrl ;; __qualifiers);
             s.status = SLAVE_READY;
             append_slave (slaves, s);
          }

        manage_slaves (slaves, &conf_handler);
     }

   while (ctrl.serial && num_retries <= max_num_retries)
     {
        variable starting_over = 0;

        _for i (0, num_indices-1, 1)
          {
             s = try_conf (indices[i], ctrl ;; __qualifiers);
             num_retries += s.num_retries;

             % indices[] may have been shuffled
             variable k = where (ordered_indices == indices[i])[0];
             pmin_final[k] = s.pmin;
             pmax_final[k] = s.pmax;

             % If an index failed to converge, it's likely to
             % happen again -- move it to the front of the list
             % so that failures happen as early as possible.
             if (s.num_retries != 0 || s.pmin == s.pmax)
               {
                  list_insert (indices, list_pop (indices, i), 0);
                  starting_over = 1;
                  break;
               }
          }

        ifnot (starting_over)
          break;
     }

   if (save != 0  && num_retries <= max_num_retries)
     {
        save_par_limits (ordered_indices,
                         pmin_final, pmax_final, prefix + ".save");
     }

   return pmin_final, pmax_final;
}
