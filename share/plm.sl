% This file is part of ISIS, the Interactive Spectral Interpretation System
% Copyright (C) 1998-2020 Massachusetts Institute of Technology
%
% This software was developed by the MIT Center for Space Research under
% contract SV1-61010 from the Smithsonian Institution.
%
%    Author:  John C. Houck  <houck@space.mit.edu>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
%%%%%%%%%%%%%%%
%
% This is a parallel implementation of the Levenberg-Marquardt
% optimization algorithm, written June 26-27 2009.
%
% The general approach used to parallelize the LM algorithm
% was taken from Cao et al (2009) ``A Parallel Levenberg-Marquardt
% Algorithm'', in the proceedings of the 23rd International Conference
% on Supercomputing, Yorktown Heights, NY, pp 450-459.

private define plmfit();
private variable PLM_Defaults = struct
{
   tol = 1.e-4,
   lambda = 0.1,
   grow_factor = 10.0,
   shrink_factor = 0.1,
   max_loops = 100,
   delta = 1.e-6,
   verbose = 0,
   num_slaves = -1,
   optimize = &plmfit
};

private define override_defaults (sinfo, q)
{
   if (q == NULL)
     return;
   variable sinfo_fields = get_struct_field_names (sinfo);

   foreach (get_struct_field_names (q))
     {
	variable name = ();
	if (any (name == sinfo_fields))
	  set_struct_field (sinfo, name, get_struct_field (q, name));
     }
}

define plm_new ()
{
   variable s = @PLM_Defaults;
   override_defaults (s, __qualifiers ());
   return s;
}

private define svd_solve (A, b)
{
   return _isis->svd_solve_intrin (A, b);
}

private define lu_solve (A, b)
{
   return _isis->lu_solve_intrin (A, b);
}

private define enforce_param_limits (p, pmin, pmax)
{
   variable i = where (p < pmin or p > pmax);
   variable n = length(i);
   if (n == 0)
     return;

   p[i] = pmin[i] + (pmax[i] - pmin[i]) * urand(n);
}

private variable
   MATRIX_TASK = 1,
   SEARCH_TASK = 2;

private variable State = struct
{
   num_pars, params, pmin, pmax, weight,
   r0, chisqr, matrix, search
};

private define perform_matrix_task (s, mmt)
{
   variable objs = recv_objs (s);

   variable id = objs[0],
     p = objs[1];

   variable r = __eval_residuals (mmt, p);

   send_msg (s, SLAVE_RESULT);
   send_objs (s, id, r);
}

private define perform_search_task (s, mmt)
{
   variable objs = recv_objs (s);
   variable
     p = objs[0],
     lambda = objs[1],
     alpha = objs[2],
     beta = objs[3];

   if (typeof(alpha) != Array_Type)
     alpha = alpha[[0],[0]];

   variable i, n = length(p);
   _for i (0, n-1, 1)
     {
        alpha[i,i] *= 1.0 + lambda;
     }

   variable dp = lu_solve (alpha, beta);
   if (dp == NULL)
     {
        dp = svd_solve (alpha, beta);
     }

   variable p1 = p + dp;

   enforce_param_limits (p1, State.pmin, State.pmax);

   variable r = __eval_residuals (mmt, p1);

   variable chisqr1 = _Inf;
   if (r != NULL)
     {
        chisqr1 = sum (sqr(__tmp(r)) * State.weight);
     }

   send_msg (s, SLAVE_RESULT);
   send_objs (s, chisqr1, lambda, p1);
}

private define lm_slave (s, mmt)
{
   forever
     {
        variable msg = recv_msg (s);

        switch (msg.type)
          {case MATRIX_TASK: perform_matrix_task (s, mmt);}
          {case SEARCH_TASK: perform_search_task (s, mmt);}
          {
             % default
             return 0;
          }

        send_msg (s, SLAVE_READY);
     }

   return 0;
}

private define matrix_queue_setup (mmt, qual)
{
   variable sqrt_delta = sqrt(qual.delta);

   variable i,
     num_pars = State.num_pars,
     weight = State.weight,
     queue = Struct_Type[num_pars],
     p = @State.params;

   State.r0 = __eval_residuals (mmt, p);
   State.chisqr = sum (sqr(State.r0) * weight);

   _for i (0, num_pars-1, 1)
     {
        variable
          pi = p[i],
          dp = (abs(pi) + sqrt_delta) * sqrt_delta;

        variable p_test = pi + dp;
        if (p_test < State.pmin[i])
          p_test = State.pmin[i];
        else if (State.pmax[i] < p_test)
          p_test = State.pmax[i];

        p[i] = p_test;
        variable q = struct {pars = @p, residuals = NULL};
        p[i] = pi;

        queue[i] = q;
     }

   State.matrix = struct
     {
        queue = queue,
        next = 0,
        num_finished = 0,
        done = 0,
        alpha, beta
     };
}

private define matrix_send_next (s)
{
   variable m = State.matrix;
   if (m.next >= State.num_pars)
     return;

   send_msg (s, MATRIX_TASK);
   variable q = m.queue[m.next];
   send_objs (s, m.next, q.pars);
   m.next++;
}

private define matrix_recv_result (s)
{
   variable m = State.matrix;
   variable objs = recv_objs (s);
   variable i = objs[0];
   m.queue[i].residuals = objs[1];
   m.num_finished++;
}

private define compute_alpha_beta ()
{
   variable i, j,
     n = State.num_pars,
     dyda = Array_Type[n],
     beta = Double_Type[n],
     alpha = Double_Type[n,n];

   _for i (0, n-1, 1)
     {
        variable q = State.matrix.queue[i];
        dyda[i] = (q.residuals - State.r0) / (q.pars[i] - State.params[i]);
        beta[i] = -sum (State.weight * State.r0 * dyda[i]);
     }

   _for i (0, n-1, 1)
     {
        variable wt_dyda = State.weight * dyda[i];
        _for j (i, n-1, 1)
          {
             alpha[i,j] = sum(wt_dyda * dyda[j]);
             alpha[j,i] = alpha[i,j];
          }
     }

   State.matrix.alpha = alpha;
   State.matrix.beta = beta;
}

private define matrix_handler (s, msg)
{
   switch (msg.type)
     {
      case SLAVE_READY:
        matrix_send_next (s);
     }
     {
      case SLAVE_RESULT:
        variable m = State.matrix;
        matrix_recv_result (s);
        if (m.num_finished == State.num_pars)
          {
             compute_alpha_beta ();
             State.matrix.done = 1;
          }
        else matrix_send_next (s);
     }
}

private define search_queue_setup (lambdas, chisqr, p, alpha, beta, sinfo)
{
   State.search = struct
     {
        params = p, alpha = alpha, beta = beta,
        lambdas = lambdas,
        num_lambdas = length(lambdas),
        next = 0,
        best = struct {chisqr=chisqr, params=p, lambda=NULL},
        num_finished = 0,
        done = 0,
        verbose = sinfo.verbose
     };
}

private define search_send_next (s)
{
   variable x = State.search;
   if (x.next == x.num_lambdas)
     return;

   variable m = State.matrix;
   send_msg (s, SEARCH_TASK);
   send_objs (s, x.params, x.lambdas[x.next], m.alpha, m.beta);
   x.next++;
}

private define search_recv_result (s)
{
   variable x = State.search;

   variable objs = recv_objs (s);
   variable
     chisqr1 = objs[0],
     lambda = objs[1],
     p1 = objs[2];

   x.num_finished++;

   variable best_label = " ";

   if (chisqr1 < x.best.chisqr)
     {
        best_label = "B";
        x.best.chisqr = chisqr1;
        x.best.lambda = lambda;
        % force typeof(params)=Array_Type
        x.best.params = [p1];
     }

   if (x.verbose > 0)
     {
        variable ps;
        print(p1, &ps);
        ps = strtrans (ps, "\n", " ");
        vmessage ("%s (%0.6g) %s", best_label, chisqr1, ps);
     }
}

private define search_handler (s, msg)
{
   switch (msg.type)
     {
      case SLAVE_READY:
        search_send_next (s);
     }
     {
      case SLAVE_RESULT:
        variable x = State.search;
        search_recv_result (s);
        if (x.num_finished == x.num_lambdas)
          {
             x.done = 1;
          }
        else search_send_next (s);
     }
}

private define computing_matrix ()
{
   return not State.matrix.done;
}

private define performing_search ()
{
   return not State.search.done;
}

private define halt_slaves (slaves)
{
   variable s;
   foreach s (slaves)
     {
        send_msg (s, SLAVE_EXITING);
     }
   manage_slaves (slaves, NULL);
}

private define plmfit (obj, in_parms, pmin, pmax, mmt)
{
   variable p = @in_parms;

   override_defaults (obj, __qualifiers());

   variable
     tol = obj.tol,
     lambda = obj.lambda,
     grow_factor = obj.grow_factor,
     shrink_factor = obj.shrink_factor,
     max_loops = obj.max_loops,
     verbose = obj.verbose;

   variable new_pars = @p;
   variable weight = __data_weights (mmt);

   State.num_pars = length(p);
   State.pmin = @pmin;
   State.pmax = @pmax;
   State.weight = @weight;

   variable s, slaves = new_slave_list();

   variable num_slaves = min([State.num_pars, _num_cpus()]);
   if (obj.num_slaves > 0)
     num_slaves = obj.num_slaves;

   loop (num_slaves)
     {
        s = fork_slave (&lm_slave, mmt);
        append_slave (slaves, s);
     }

   loop (max_loops)
     {
        State.params = @p;

        matrix_queue_setup (mmt, obj);
        foreach s (slaves) {matrix_send_next (s);}
        manage_slaves (slaves, &matrix_handler; while_=&computing_matrix);

        variable
          chisqr = State.chisqr,
          alpha = State.matrix.alpha,
          beta = State.matrix.beta;

        variable powers = [0:num_slaves-1];
        forever
          {
             variable p1, chisqr1,
               lambdas = lambda * grow_factor^powers;

             search_queue_setup (lambdas, chisqr, p, alpha, beta, obj);
             foreach s (slaves) {search_send_next (s);}
             manage_slaves (slaves, &search_handler; while_=&performing_search);

             variable b = State.search.best;
             p1 = b.params;
             chisqr1 = b.chisqr;

             if (b.lambda != NULL)
               lambda = b.lambda;
             else lambda = lambdas[-1];

             if (chisqr1 < chisqr)
               break;

             if (lambda > 1.e10)
               {
                  if (verbose > 0) vmessage ("too many iterations");
                  halt_slaves (slaves);
                  in_parms[*] = p1;
                  return -1;
               }

             powers += powers[-1] + 1;
          }

        lambda *= shrink_factor;
        if (lambda < DOUBLE_EPSILON)
          lambda = DOUBLE_EPSILON;

        p = @p1;

        if (chisqr - chisqr1 <= tol * chisqr)
          {
             halt_slaves (slaves);
             in_parms[*] = p;
             return 0;
          }

        chisqr = chisqr1;
     }

   halt_slaves (slaves);
   in_parms[*] = p;
   return 0;
}

