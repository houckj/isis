% Nelder-Mead Simplex
% Copyright (C) 2009 John E. Davis
%
% This file is part of the S-Lang Optimization Package and may be
% distributed under the terms of the GNU General Public License.  See
% the file COPYING for more information.

require ("rand");

private define simplex_usage ()
{
   variable str = "Usage:\n";
   str += "minF = obj.optimize(parms, min_parms, max_parms, func, func_data);\n";
   str += "Qualifiers:\n";
   str += " reflect=1.0, contact=0.5, stretch=2.0, shrink=0.5,\n";
   str += " maxnfe=1000, reldiff=1e-4, absdiff=1e-9, scale=0.5\n";
   str += " verbose=0, maxsize=0\n",
   str += "\n";
   str += "The value at the minimum is returned.  The parms array values are\n";
   str += "set to the parameters at the minimum.\n";
   str += "func should be defined as: define func(parms, func_data)\n";
   
   usage (str);
}

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

private define simplex ();
define simplex_new ()
{
   variable s = struct
     {
	% Required of all Optimization Class objects
	maxnfe = 1000,
	nfe = 0,
	converged = 0,
	reldiff = 1e-4,
	absdiff = 1e-9,
	verbose = 0,
	new_min_hook = NULL,
	new_min_hook_args = {},
	optimize = &simplex,
	
	% Information specific to this class
	reflect = 1.0,
	contract = 0.5,
	stretch = 2.0,
	shrink = 0.5,
	max_size = 0,
	scale = 0.5,
	% if non-zero projected points are not placed at the boundary
	avoid_bc = 0,

     };
   
   override_defaults (s, __qualifiers());
   return s;
}

private define print_stats (parm_sets, com)
{
   variable n = length(parm_sets)-1;
   variable ssum = 0.0, smax = 0.0;
   _for (0, n-1, 1)
     {
	variable i = ();
	variable dp = (parm_sets[i] - com);
	variable s = sum (dp*dp);
	if (s > smax)
	  smax = s;
	ssum += s;
     }
   ssum = sqrt (ssum/n);
   vmessage ("Simplex RMS size: %g, (max dist=%g)", ssum, sqrt(smax));
}

private define sort_simplex (parm_sets, func_vals)
{
   variable i = array_sort (func_vals);
   variable parm_sets_orig = @parm_sets;
   variable func_vals_orig = @func_vals;
   _for (0, length(parm_sets)-1, 1)
     {
	variable j = ();
	variable i_j = i[j];
	parm_sets[j] = parm_sets_orig[i_j];
	func_vals[j] = func_vals_orig[i_j];
     }
}

% Compute COM from first n sets
private define compute_com (parm_sets, func_vals, n, sizep)
{
   variable i, j;
   variable p0 = parm_sets[0];
   variable com = p0;
   variable num = 1;

   _for i (1, n-1, 1)
     {
	variable p = parm_sets[i];
	com += p;
	num++;
     }
   com/=num;

   variable size = 0;
   _for i (0, n-1, 1)
     {
	size = _max (size, maxabs (parm_sets[i]-p0));
     }
   @sizep = size;
   return com;
}

% contract all points toward the nth one by a factor of lambda.
private define contract (parm_sets, nth, lambda)
{
   variable i;
   variable parms_nth = parm_sets[nth];

   _for i (0, length(parm_sets)-1, 1)
     {
	if (i == nth)
	  continue;
	
	parm_sets[i] = parms_nth + lambda * (parm_sets[i] - parms_nth);
     }
}

private define project_parms (info, p0, p1, t)
{
   variable ileft, iright;
   variable p2 = p0 + t*(p1-p0);   % p0 + (p0-p1)*(-t)
   ifnot (any((p2 < 0.0) or (p2 > 1.0)))
     return p2;

   ileft = where (p2 < 0.0);
   iright = where (p2 > 1.0);
   variable dp = p1-p0;
   variable tleft = (-p0/dp)[ileft];
   variable tright = ((1.0-p0)/dp)[iright];

   t = sign(t) * minabs([tleft,tright]);
   if (info.avoid_bc)
     t *= 0.25 + 0.5*rand_uniform ();
   p2 = p0 + t*dp;
   % Avoid roundoff error
   p2[where(p2<0.0)] = 0.0;
   p2[where(p2>1.0)] = 1.0;
   return p2;
}

private define call_func (func, parms, func_data)
{
   variable val = (@func)(parms, func_data);
   if (isnan (val))
     {
	() = fprintf (stderr, "NaN value seen-- using Inf\n");
	val = _Inf;
     }
   return val;
}

private define call_normalized_func (func, all_p, i, p, func_data, p0, dp)
{
   all_p[i] = p0+p*dp;
   return call_func (func, all_p, func_data);
}

private define simplex ()
{
   if (_NARGS == 0)
     simplex_usage ();

   variable sinfo, in_parms, min_parms, max_parms, func, func_data;
   
   (sinfo, in_parms, min_parms, max_parms, func, func_data) = ();

   min_parms = double(min_parms);
   max_parms = double(max_parms);

   ifnot (all(min_parms <= max_parms))
     throw InvalidParmError, "One or more minparms > maxparms";

   variable i = (min_parms <= in_parms <= max_parms);
   ifnot (all(i))
     {
	() = fprintf (stderr, "*** WARNING: simplex called with initial parameters out of min/max range\n");
	() = fprintf (stderr, "*** WARNING: Randomly choosing a point within the allowed range\n");
	i = wherenot (i);
	in_parms[i] = min_parms[i] + rand_uniform(i) * (max_parms[i]-min_parms[i]);
     }

   if (sinfo == NULL)
     sinfo = simplex_new ();
   else
     sinfo = @sinfo;

   override_defaults (sinfo, __qualifiers());

   variable vindices = where (min_parms < max_parms);
   variable vmin_parms = min_parms[vindices];
   variable vmax_parms = max_parms[vindices];
   variable nparms = length (vindices);
   variable nsets = nparms + 1;

   variable parm_sets=Array_Type[nsets];
   variable scale = sinfo.scale;
   variable verbose = sinfo.verbose;
   variable max_size = sinfo.max_size;

   variable dvparms = vmax_parms - vmin_parms;

   variable func_vals = Double_Type[nsets];
   variable nfe = sinfo.nfe;
   
   if (verbose)
     vmessage ("************* Computing simplex point %d/%d ****", 0, nsets);

   variable val;
   loop (50)
     {
	val = (@func)(in_parms, func_data);
	ifnot (isnan (val))
	  break;

	() = fprintf (stderr, "Warning: NaN found, choosing a random parameter set\n");
	in_parms[vindices] = vmin_parms + rand_uniform(dvparms)*dvparms;
     }
   then 
     {
	() = fprintf (stderr, "Unable to generate non-NaN values\n");
	return _NaN;
     }
   func_vals[0] = val;

   variable reldiff = sinfo.reldiff;
   variable absdiff = sinfo.absdiff;

   % From here on, normalized parameter sets will be used such that:
   %   min_p <= p <= max_p
   %     0  <= p' <= 1
   %   p' = (p-min_p)/(max_p-min_p)
   %   p = min_p + (max_p-min_p)*p'
   variable vparms = in_parms[vindices];
   parm_sets[0] = (vparms-vmin_parms)/dvparms;

   i = 1;
   variable j = 0;
   variable parm_sets_0 = parm_sets[0];
   while (i < nsets)
     {
	variable e_j;
#iffalse
	e_j = Double_Type[nparms]; e_j[j] = 1.0;
	if (parm_sets_0[j] > 0.5)
	  e_j[j] = -1.0;
#else
	e_j = (0.5-rand_uniform (nparms));
	e_j /= hypot (e_j);
#endif
	parm_sets[i] = project_parms (sinfo, parm_sets_0, parm_sets_0+e_j, scale);

	if (verbose)
	  vmessage ("************* Computing simplex point %d/%d ****", i, nsets);

	func_vals [i] 
	  = call_normalized_func (func, in_parms, vindices, parm_sets[i], func_data,
				  vmin_parms, dvparms);
	nfe++;
	i++;
	j++;
     }

   variable maxnfe = sinfo.maxnfe;
   variable stretch_factor = sinfo.stretch;
   variable contract_factor = sinfo.contract;
   variable reflect_factor = sinfo.reflect;
   variable shrink_factor = sinfo.shrink;

   %stretch_factor = 1.0 + 3.0/nparms;

   variable indices = [0:nsets-1];
   variable last_lo = NULL;

   while (nfe < maxnfe)
     {
	%sort_simplex (parm_sets, func_vals, indices);
	i = array_sort (func_vals);
	indices = indices[i];
	func_vals = func_vals[i];
	parm_sets = parm_sets[i];

	variable jth = nsets-1;
	%indices = sort_simplex1 (func_vals);
	variable hi = func_vals[-1];
	variable lo = func_vals[0];
	variable hi2 = func_vals[-2];

	variable this_vertex = indices[jth];
	variable try_parms, try_y, save_try_y, save_try_parms;

	if ((last_lo == NULL) || (lo < last_lo))
	  {
	     if (sinfo.new_min_hook != NULL)
	       {
		  in_parms[vindices] = vmin_parms + dvparms*parm_sets[0];
		  (@sinfo.new_min_hook)(in_parms, lo,
					__push_list(sinfo.new_min_hook_args));
	       }
	     last_lo = lo;
	  }


	variable hi_parms = parm_sets[jth];
	variable size;
	variable com = compute_com (parm_sets, func_vals, nsets-1, &size);

	if (hi == lo)
	  {
	     if (verbose)
	       vmessage ("%s", "********* hi == lo ************, checking com");
	     try_parms = com;
	     try_y = call_normalized_func (func, in_parms, vindices, try_parms, func_data, vmin_parms, dvparms);
	     nfe++;
	     if (try_y == lo)
	       break;
	     if (try_y < lo)
	       {
		  func_vals[jth] = try_y;
		  parm_sets[jth] = try_parms;
		  continue;
	       }
	     break;
	  }
	  
	if (all (feqs (parm_sets[0], hi_parms, reldiff, absdiff)))
	  {
	     %vmessage ("hi=%g, lo=%g", hi, lo);
	     %print (func_vals);
	     %print_stats (parm_sets, com);
	     break;
	  }

	if (verbose)
	  vmessage ("nfe=%d:  hi=%g, hi2=%g, lo=%g, size=%g", nfe, hi, hi2, lo, size);

	if (verbose)
	  print_stats (parm_sets, com);
	
	% Try to reflect highest about COM
	% Tweak the COM to avoid a degenerate simplex
	%com += 0.001*(hi_parms-com)*(rand_uniform(com)-0.5);

	if (verbose)
	  vmessage ("********* reflecting %d ************", this_vertex);
	try_parms = project_parms (sinfo, com, hi_parms, -reflect_factor);
	try_y = call_normalized_func (func, in_parms, vindices, try_parms, func_data, vmin_parms, dvparms);
	nfe++;

	if (try_y <= lo)
	  {
	     if (nfe >= maxnfe)
	       {
		  func_vals[jth] = try_y;
		  parm_sets[jth] = try_parms;
		  break;
	       }

	     save_try_y = try_y;
	     save_try_parms = try_parms;
	     % try extra stretch
	     if (verbose)
	       vmessage ("********** stretching %d **********", this_vertex);
	     try_parms = project_parms (sinfo, com, try_parms, stretch_factor);
	     try_y = call_normalized_func (func, in_parms, vindices, try_parms, func_data, vmin_parms, dvparms);
	     nfe++;

	     if (try_y > save_try_y)
	       {
		  try_y = save_try_y;
		  try_parms = save_try_parms;
	       }
	     func_vals[jth] = try_y;
	     parm_sets[jth] = try_parms;
	     continue;
	  }

	if (try_y < hi2)
	  {
	     func_vals[jth] = try_y;
	     parm_sets[jth] = try_parms;
	     continue;
	  }

	% we know it is >= hi2, so try to contract it.  But do it in a smart
	% way.  Consider:
	%
	%                          P_R>=hi
	%    P_hi                  P_R>=hi2
	%          P_hi2
	%                  COM
	%
	% If P_R>=hi, then contract to point betweeen P_hi and COM.  Otherwise
	% contract to point between COM and P_R.
	
	variable ref_parms = try_parms;
	variable ref_y = try_y;

	if (try_y > hi)
	  {
	     try_parms = hi_parms;
	     try_y = hi;
	     % This is not in the simplex algorithm.
	     %com = compute_com(parm_sets, nsets/2);
	     if (verbose)
	       vmessage ("********** left contracting %d **********", this_vertex);
	  }
	else if (verbose)
	  vmessage ("********** right contracting %d **********", this_vertex);


	save_try_parms = try_parms;
	save_try_y = try_y;

	variable contract_parms, contract_y;
	contract_parms = project_parms (sinfo, com, try_parms, contract_factor);
	contract_y = call_normalized_func (func, in_parms, vindices, contract_parms, func_data, vmin_parms, dvparms);
	nfe++;
	
	% The test below used to use `<='.  But I have seen cases where this
	% led to a repeated reflect-contract cycle that went no where.
	% So it is better to simply shrink.
	if (contract_y < save_try_y)
	  {
	     func_vals[jth] = contract_y;
	     parm_sets[jth] = contract_parms;
	     continue;
	  }
	if (nfe >= maxnfe)
	  break;

#ifntrue
	if (verbose)
	  vmessage ("********** shrinking %d **********", jth, nsets);
	parm_sets[jth] = parm_sets[0] + shrink_factor * (parm_sets[jth] - parm_sets[0]);
	func_vals[jth] = call_normalized_func (func, in_parms, vindices, parm_sets[jth], func_data, vmin_parms, dvparms);
	nfe++;
#else	
	% contraction failed.  So shrink.
	contract (parm_sets, 0, shrink_factor);
	_for i (1, nsets-1, 1)
	  {
	     if (verbose)
	       vmessage ("********** shrinking %d/%d**********", i, nsets);
	     func_vals[i] = call_normalized_func (func, in_parms, vindices, parm_sets[i], func_data, vmin_parms, dvparms);
	     nfe++;
	     if (nfe >= maxnfe)
	       break;
	  }
#endif
     }
   
   sinfo.nfe = nfe;

   sort_simplex (parm_sets, func_vals);
   in_parms[vindices] = vmin_parms + dvparms*parm_sets[0];

   % The hook may not have been called if a lower value was hit during
   % the shrink step, but the number of function evaluations was
   % exhausted.
   if ((sinfo.new_min_hook != NULL)
       && ((last_lo == NULL) || (func_vals[0] < last_lo)))
     (@sinfo.new_min_hook)(in_parms, func_vals[0],
			   __push_list(sinfo.new_min_hook_args));
   return func_vals[0];
}
