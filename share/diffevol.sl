% Differential Evolution
% Copyright (C) 2009 John E. Davis
%
% This file is part of the S-Lang Optimization Package and may be
% distributed under the terms of the GNU General Public License.  See
% the file COPYING for more information.
%
require ("rand");

private variable Sort_Samples=0;
private define sample_population (s, current_index, n)
{
   variable npop = s.npop;
   variable samples = Int_Type[n+1] - 1;
   samples[n] = current_index;

   _for (0, n-1, 1)
     {
	variable i = ();
	do
	  {
	     variable sample = rand_int (0, npop-1);
	  }
	while (any(sample == samples));
	samples[i] = sample;
     }

   samples = samples[[0:n-1]];
   if (Sort_Samples)
     {
	variable vals = s.population_energies[samples];
	samples = samples[array_sort (vals)];
     }
   return samples;
}

private define binary_crossover (parms, mutant_parms, prob)
{
   variable i = where(rand_uniform(parms)<prob);
   ifnot (length (i))
     i = rand_int (0, length(parms)-1);
   parms[i] = mutant_parms[i];
   return parms;
}

private define exp_crossover (parms, mutant_parms, prob)
{
   variable n = length(parms);
   variable i = rand_int (0, n-1);
   loop (n)
     {
	parms[i] = mutant_parms[i];
	if (rand_uniform ()>prob)
	  break;
	i = (i + 1) mod n;
     }
   return parms;
}

private variable Strategies = Assoc_Type[Ref_Type];

private define trig_mutant (s, p1, p2, p3, w1, w2, w3)
{
   variable w0 = s.best_energy;
   w1 -= w0;
   w2 -= w0;
   w3 -= w0;

   variable norm = (w1 + w2 + w3)/3.0;
   if (isinf (norm))
     return NULL;
   
   variable mutant = (p1+p2+p3)/3.0;
   if (norm != 0.0)
     {
	w1 = w1/norm;
	w2 = w2/norm;
	w3 = w3/norm;
	mutant -= w1*(p1-mutant) + w2*(p2-mutant) + w3*(p3-mutant);
     }
   return mutant;
}

private define rand1_mutant (s, current_index)
{
   variable samples = sample_population (s, current_index, 3);
   variable parms0 = s.population_parms[samples[0]];
   variable parms1 = s.population_parms[samples[1]];
   variable parms2 = s.population_parms[samples[2]];

   if (rand_uniform () < s.trigprob)
     {
	variable mutant;
	mutant = trig_mutant (s, parms0, parms1, parms2, 
			      s.population_energies[samples[0]],
			      s.population_energies[samples[1]],
			      s.population_energies[samples[2]]);
	if (mutant != NULL)
	  return mutant;
     }

   return parms0 + s.diffscale * (parms1 - parms2);
}
   
private define rand1exp (s, current_index)
{
   variable cparms = @s.population_parms[current_index];
   variable mutant = rand1_mutant (s, current_index);
   return exp_crossover (cparms, mutant, s.crossprob);
}
Strategies["rand1exp"] = &rand1exp;

private define rand1bin (s, current_index)
{
   variable cparms = @s.population_parms[current_index];
   variable mutant = rand1_mutant (s, current_index);
   return binary_crossover (cparms, mutant, s.crossprob);
}
Strategies["rand1bin"] = &rand1bin;

private define best1_mutant (s, current_index)
{
   variable samples = sample_population (s, current_index, 2);
   variable parms0 = s.best_parms;
   variable parms1 = s.population_parms[samples[0]];
   variable parms2 = s.population_parms[samples[1]];

   if (rand_uniform () < s.trigprob)
     {
	variable mutant;
	mutant = trig_mutant (s, parms0, parms1, parms2, 
			      s.best_energy,
			      s.population_energies[samples[0]],
			      s.population_energies[samples[1]]);
	if (mutant != NULL)
	  return mutant;
     }

   return parms0 + s.diffscale * (parms1 - parms2);
}

private define rand_to_best_mutant (s, current_index)
{
   variable samples = sample_population (s, current_index, 2);
   variable cparms = s.population_parms[current_index];
   variable parms0 = s.best_parms;
   variable parms1 = s.population_parms[samples[0]];
   variable parms2 = s.population_parms[samples[1]];

   if (rand_uniform () < s.trigprob)
     {
	variable mutant;
	mutant = trig_mutant (s, parms0, parms1, parms2, 
			      s.best_energy,
			      s.population_energies[samples[0]],
			      s.population_energies[samples[1]]);
	if (mutant != NULL)
	  return mutant;
     }

   return cparms + s.diffscale * ((parms0-cparms) + (parms1-parms2));
}

private define best1exp (s, current_index)
{
   variable cparms = @s.population_parms[current_index];
   variable mutant = best1_mutant (s, current_index);
   return exp_crossover (cparms, mutant, s.crossprob);
}
Strategies["best1exp"] = &best1exp;

private define randtobest1exp (s, current_index)
{
   variable cparms = @s.population_parms[current_index];
   variable mutant = rand_to_best_mutant (s, current_index);
   return exp_crossover (cparms, mutant, s.crossprob);
}
Strategies["randtobest1exp"] = &randtobest1exp;

private define best1bin (s, current_index)
{
   variable cparms = @s.population_parms[current_index];
   variable mutant = best1_mutant (s, current_index);
   return binary_crossover (cparms, mutant, s.crossprob);
}
Strategies["best1bin"] = &best1bin;

private define randtobest1bin (s, current_index)
{
   variable cparms = @s.population_parms[current_index];
   variable mutant = rand_to_best_mutant (s, current_index);
   return binary_crossover (cparms, mutant, s.crossprob);
}
Strategies["randtobest1bin"] = &randtobest1bin;

private define mixedbin (s, current_index)
{
   variable cutoff_en = s.best_energy + 0.3*(s.worst_energy - s.best_energy);
   if (s.population_energies[current_index] > cutoff_en)
     return best1bin (s, current_index);
   
   return rand1bin (s, current_index);
}
Strategies["mixedbin"] = &mixedbin;

private define mixedexp (s, current_index)
{
   variable cutoff_en = s.best_energy + 0.3*(s.worst_energy - s.best_energy);
   if (s.population_energies[current_index] > cutoff_en)
     return randtobest1exp (s, current_index);
   
   return rand1exp (s, current_index);
}
Strategies["mixedexp"] = &mixedexp;

private define lookup_strategy (strategy)
{
   if (assoc_key_exists (Strategies, strategy))
     return Strategies[strategy];
   
   throw InvalidParmError, "Unknown strategy: $strategy"$;
}

private define de_usage ()
{
   variable str = "Usage:\n";
   str += "minF = diffevol(sinfo, parms, min_parms, max_parms, func, func_data);\n";
   str += "Qualifiers:\n";
   str += " npop=10*length(parms), maxnfe=50000, reldiff=1e-4, absdiff=1e-9,\n";
   str += " diffscale=0.7, crossprob=0.5 max_gens=-1\n";
   str += " strategy=\"best1bin|best1exp|rand1exp|rand1bin\"\n";
   str += " bc=\"none|reject|reflect|truncate|rand\"\n",
   str += " verbose=0\n",
   str += "\n";
   str += "sinfo may be NULL or it may be created via a call to diffevol_new.\n";
   str += "It contains fields that correspond to the qualifiers.\n";
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

private define diffevol ();
define diffevol_new ()
{
   variable s = struct
     {
	npop,
	maxnfe = 50000,
	max_gens = -1,
	diffscale = 0.7,
	crossprob = 0.5,
	strategy = "best1bin",
	trigprob = 0.0,
	ode_rate = 0.0,
	verbose = 0,
	reldiff = 1e-4,
	absdiff = 1e-9,
	bc = "truncate",
	population_parms,
	population_energies,
	best_energy,
	best_parms,
	worst_energy,
	nfe,
	
	optimize = &diffevol,
	converged = 0,
     };
   override_defaults (s, __qualifiers());
   return s;
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

private define bc_truncate (parms, min_parms, max_parms)
{
   variable i = where (parms < min_parms);
   parms[i] = min_parms[i];
   i = where (parms > max_parms);
   parms[i] = max_parms[i];
   return 0;
}

private define bc_rand (parms, min_parms, max_parms)
{
   variable i = wherenot (min_parms <= parms <= max_parms);
   if (length (i))
     parms[i] = min_parms[i] + (max_parms[i]-min_parms[i])*rand_uniform(i);

   return 0;
}

private define bc_reflect (parms, min_parms, max_parms)
{
   variable i;
   variable again = 1;

   % This is flawed when parms<<min_parms or parms>>max_parms.  For
   % such a case, "mod" has to be used.
   variable nloops = 10;
   while (again && nloops)
     {
	again = 0;
	i = where (parms < min_parms);
	if (length (i))
	  {
	     parms[i] = (2*min_parms - parms)[i];
	     again++;
	  }
	i = where (parms > max_parms);
	if (length (i))
	  {
	     parms[i] = (2*max_parms - parms)[i];
	     again++;
	  }
	nloops--;
     }
   if (nloops == 0)
     {
	return bc_rand (parms, min_parms, max_parms);
     }
   return 0;
}

private define bc_reject (parms, min_parms, max_parms)
{
   if (all (min_parms <= parms <= max_parms))
     return 0;

   return -1;
}

private define bc_none (parms, min_parms, max_parms)
{
   return 0;
}

private define lookup_bcfunc (bc)
{
   if (bc == NULL) bc = "none";
   variable bcs = ["none", "reject", "reflect", "truncate", "rand"];
   variable funs = [&bc_none, &bc_reject, &bc_reflect, &bc_truncate, &bc_rand];
   variable i = wherefirst (bc == bcs);
   if (i == NULL)
     throw InvalidParmError, "Invalid bc qualifier: $bc"$;
   return funs[i];
}

private define initialize_population (npop, in_parms, min_parms, max_parms,
				      func, funcdata,
				      verbose, nfep, maxnfe)
{
   variable i = 0;
   variable pop_parms = Array_Type[npop];
   variable pop_energies = Double_Type[npop];
   variable nfe = 0;
   variable dparms = max_parms - min_parms;
   variable parms = @in_parms;
   
   pop_energies[*] = _Inf;
   while ((i < npop) && (nfe < maxnfe))
     {
	variable en = call_func (func, parms, funcdata); nfe++;
	ifnot (isinf (en))
	  {
	     if (verbose)
	       vmessage ("Found initial parameters for individual %d: en=%S, nfe=%d",
			 i, en, nfe);
	     pop_parms[i] = parms;
	     pop_energies[i] = en;
	     i++;
	  }
	parms = min_parms + rand_uniform(dparms)*dparms;
     }
   @nfep = nfe;
   return pop_parms, pop_energies;
}

private define compute_population_size (parm_sets)
{
   variable nsets = length (parm_sets);
   variable min_parms = @parm_sets[0];
   variable max_parms = @min_parms;
   variable nparms = length (min_parms);
   variable i, j;
   _for i (1, nsets-1, 1)
     {
	variable parms = parm_sets[i];
	_for j (0, nparms-1, 1)
	  {
	     if (parms[j] < min_parms[j]) min_parms[j] = parms[j];
	     if (parms[j] > max_parms[j]) max_parms[j] = parms[j];
	  }
     }
   return min_parms, max_parms;
}

private define merge_opposite_population (next_population_energies, next_population_parms, op_pop_energies, op_pop_parms)
{
   variable i = wherenot (isnan (op_pop_energies));
   variable pop_parms = [next_population_parms, op_pop_parms[i]];
   variable pop_energies = [next_population_energies, op_pop_energies[i]];
   i = array_sort (pop_energies);
   _for (0, length(next_population_energies)-1, 1)
     {
	variable k = ();
	variable i_k = i[k];
	next_population_energies[k] = pop_energies[i_k];
	next_population_parms[k] = pop_parms[i_k];
     }
}

private define diffevol ()
{
   if (_NARGS != 6)
     {
	_pop_n (_NARGS);
	de_usage ();
     }

   variable sinfo, in_parms, min_parms, max_parms, func, func_data;
   (sinfo, in_parms, min_parms, max_parms, func, func_data) = ();

   min_parms = double(min_parms);
   max_parms = double(max_parms);

   ifnot (all(min_parms <= max_parms))
     throw InvalidParmError, "One or more minparms > maxparms";

   variable i = (min_parms <= in_parms <= max_parms);
   variable dsize = max_parms - min_parms;
   ifnot (all(i))
     {
	() = fprintf (stderr, "*** WARNING: diffevol called with initial parameters out of min/max range\n");
	() = fprintf (stderr, "*** WARNING: Randomly choosing a point within the allowed range\n");
	i = wherenot (i);
	in_parms[i] = min_parms[i] + rand_uniform(i) * dsize[i];
     }
   
   if (sinfo == NULL)
     sinfo = diffevol_new ();
   override_defaults (sinfo, __qualifiers());

   variable bcfun = lookup_bcfunc (sinfo.bc);
   variable strategy = lookup_strategy (sinfo.strategy);

   variable population_energies, population_parms, npop;

   variable nfe = sinfo.nfe;
   variable maxnfe = sinfo.maxnfe;
   variable verbose = sinfo.verbose;
   variable energy;

   if (sinfo.population_parms == NULL)
     {
	npop = sinfo.npop;
	if (npop == NULL)
	  npop = 10*length(in_parms);

	(population_parms, population_energies)
	    = initialize_population (npop, in_parms, min_parms, max_parms,
				     func, func_data,
				     verbose, &nfe, maxnfe);

	sinfo.population_parms = population_parms;	
	sinfo.population_energies = population_energies;
     }
   else
     {
	population_parms = sinfo.population_parms;
	population_energies = sinfo.population_energies;
	npop = length (population_parms);

	if (population_energies == NULL)
	  {
	     population_energies = Double_Type[npop];
	     _for i (0, npop-1, 1)
	       {
		  energy = call_func (func, population_parms[i], func_data);
		  population_energies[i] = energy;
	       }
	     sinfo.population_energies = population_energies;
	     nfe = npop;
	  }
     }
   sinfo.npop = npop;

   i = wherefirst (population_energies == min(population_energies));
   sinfo.best_parms = population_parms[i];
   sinfo.best_energy = population_energies[i];
   sinfo.converged = 0;

   variable best_energy = sinfo.best_energy;
   variable reldiff = sinfo.reldiff;
   variable absdiff = sinfo.absdiff;
   variable max_gens = sinfo.max_gens;

   variable next_population_parms = @population_parms;
   variable next_population_energies = @population_energies;
   while (nfe < maxnfe)
     {
	if (max_gens == 0)
	  break;
	if (max_gens > 0)
	  max_gens--;

	variable min_energy, max_energy;

	max_energy = max(population_energies);
	min_energy = min(population_energies);
	if (feqs (min_energy, max_energy, reldiff, absdiff))
	  {
	     sinfo.converged = 1;
	     break;
	  }
	if (verbose)
	  vmessage ("nfe=%d: min_energy=%g, max_energy=%g", nfe, min_energy, max_energy);
	
	sinfo.worst_energy = max_energy;

	variable changed = 0;
	_for i (0, npop-1, 1)
	  {
	     variable try_parms;

	     do
	       {
		  do
		    {
		       try_parms = (@strategy)(sinfo, i);
		    }
		  while (-1 == (@bcfun)(try_parms, min_parms, max_parms));
		  energy = call_func (func, try_parms, func_data);
		  nfe++;
	       }
	     while (isinf(energy));

	     if (verbose>2)
	       {
		  vmessage ("nfe=%d: changed: %d, candidate %d/%d(%g) : %g",
			    nfe, changed, i, npop, sinfo.population_energies[i], energy);
	       }
	     % Accept a value only if it lowers the energy.  A tie
	     % could mean that the parameter's value has no effect.
	     if (energy < sinfo.population_energies[i])
	       {
		  changed++;

		  if (verbose>1)
		    vmessage ("Found better candidate %d: %g ->%g", i, sinfo.population_energies[i], energy);
		  next_population_energies[i] = energy;
		  next_population_parms[i] = try_parms;

		  if (energy < best_energy)
		    {
		       if (verbose)
			 () = fprintf (stdout, "*** Lower energy found: %g @ nfe=%d\n",
				       energy, nfe);

		       best_energy = energy;
		       sinfo.best_energy = best_energy;
		       sinfo.best_parms = try_parms;
		    }
	       }

	     if (nfe >= maxnfe)
	       break;

	  }

	if (nfe >= maxnfe)
	  break;

	if ((sinfo.ode_rate > 0) && (rand_uniform () < sinfo.ode_rate))
	  {
	     if (verbose)
	       () = fprintf (stdout, "*** Computing opposite population\n");
	     variable op_min_parms, op_max_parms;
	     
	     (op_min_parms, op_max_parms) = compute_population_size (next_population_parms);
	     variable sum_op_min_max_parms = op_min_parms + op_max_parms;
	     variable op_pop_energies = Double_Type[npop];
	     variable op_pop_parms = Array_Type[npop];
	     _for i (0, npop-1, 1)
	       {
		  try_parms = sum_op_min_max_parms - next_population_parms[i];
		  energy = call_func (func, try_parms, func_data);
		  nfe++;
		  if (verbose>2)
		    {
		       vmessage ("nfe=%d: opposite candidate %d/%d(%g) : %g",
				 nfe, i, npop, next_population_energies[i], energy);
		    }
		  op_pop_parms[i] = try_parms;
		  op_pop_energies[i] = energy;
	       }
	     merge_opposite_population (next_population_energies, next_population_parms, op_pop_energies, op_pop_parms);
	  }
	
	population_energies[*] = next_population_energies;
	_for i (0, npop-1, 1)
	  population_parms[i] = next_population_parms[i];
     }
   
   sinfo.nfe = nfe;
   in_parms[*] = sinfo.best_parms;
   return best_energy;
}
