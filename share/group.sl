% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2010 Massachusetts Institute of Technology
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
%

private define mask_intervals ()
{
   variable nargs = _NARGS;
   _stk_roll(-nargs);
   variable s = ();
   nargs--;

   if (nargs mod 2 != 0)
     {
        _pop_n(nargs);
        throw ApplicationError;
     }

   variable len = length(s.bin_lo),
     mask = Int_Type[len];

   variable num_pairs = nargs / 2;
   while (num_pairs > 0)
     {
        variable lo, hi;
        (lo, hi) = ();

        variable xa = to_angstrom ([lo,hi];; __qualifiers);
        (lo, hi) = (xa[0], xa[1]);

        variable i = where (lo <= s.bin_lo and s.bin_hi <= hi);
        mask[i] = 1;

        num_pairs--;
     }

   return mask;
}

private define mask_values (datasets)
{
   variable d = get_data_counts(datasets[0]),
     len = length(d.value),
     veto = Int_Type[len];

   variable
     min_sum = qualifier ("min_sum", NULL),
     min_val = qualifier ("min_val", NULL);

   if (min_sum == NULL && min_val == NULL)
     return not veto;

   variable i, s = Double_Type[len];

   foreach i (datasets)
     {
        variable vi = get_data_counts(i).value;
        s += vi;
        if (min_val != NULL)
          {
             veto[where (vi < min_val)] = 1;
          }
     }

   if (min_sum != NULL)
     {
        veto[where (s < min_sum)] = 1;
     }

   return not veto;
}

private define mismatched_grids (datasets)
{
   if (length(datasets) == 1)
     return 0;

   variable i, s = get_data_counts (datasets[0]);
   foreach i (datasets[[1:]])
     {
        variable ss = get_data_counts (i);
        if (any (ss.bin_lo != s.bin_lo or ss.bin_hi != s.bin_hi))
          return 1;
     }
   return 0;
}

private define perform_select ()
{
   variable nargs = _NARGS;

   variable select_method = ();
   nargs--;

   variable exclusive = ();
   nargs--;

   _stk_roll (-nargs);
   variable datasets = ();
   nargs--;

   variable args = __pop_args (nargs);

   if (mismatched_grids (datasets))
     throw ApplicationError, "*** Error:  mismatched dataset grids";

   if (exclusive) ignore (datasets);

   variable value_mask, interval_mask;
   value_mask = mask_values (datasets ;; __qualifiers);

   variable grid = get_data_counts (datasets[0]);
   interval_mask = mask_intervals (grid, __push_args(args) ;; __qualifiers);

   (@select_method)(datasets, where (value_mask and interval_mask));
}

define ignore_values ()
{
   variable msg = "ignore_values (datasets, lo1, hi1 [, lo2, hi2...] ; qualifiers)\n"
                + "   qualifiers:  min_sum, min_val, unit";
   if (_NARGS < 3) usage (msg);
   variable args = __pop_args (_NARGS);
   perform_select (__push_args (args), 0, &ignore_list ;; __qualifiers);
}

define notice_values ()
{
   variable msg = "notice_values (datasets, lo1, hi1 [, lo2, hi2...] ; qualifiers)\n"
                + "   qualifiers:  min_sum, min_val, unit";
   if (_NARGS < 3) usage (msg);
   variable args = __pop_args (_NARGS);
   perform_select (__push_args (args), 1, &notice_list ;; __qualifiers);
}

private define init_data_struct (k)
{
   rebin_data (k, 0);

   variable s = struct {value, back, scaled_back, bin_lo, bin_hi};

   variable d = get_data_counts(k);
   s.bin_lo = d.bin_lo;
   s.bin_hi = d.bin_hi;
   s.value = d.value;
   s.back = get_back (k);

   variable scale = 1.0;
   if (s.back != NULL)
     {
        scale = [((get_data_exposure (k) * get_data_backscale(k))
                  /(get_back_exposure (k) * get_back_backscale(k)))];
        scale[where(isnan(scale))] = 1.0;
        if (length(scale) == 1) scale = scale[0];
     }
   else s.back = Double_Type[length(s.value)];

   s.scaled_back = scale * s.back;

   return s;
}

private define force_array_length (type, n, q)
{
   variable value = type[n];

   if (q != NULL)
     {
        variable nq = length(q);

        switch (nq)
          {case 1:  value[*] = q;}
          {(nq < n):
             value[[0:nq-1]] = q;
             value[[nq:n-1]] = q[-1];
          }
          {(nq >= n):
            value = q[[0:n-1]];
          }
     }

   return value;

}

public define group ()
{
   variable msg = "group (datasets [; qualifiers])\n"
     + "   qualifiers:  min_sn, min_chan, bounds, unit, sn_data";

   if (_NARGS == 0)
     usage(msg);

   variable datasets = ();

   variable
     data = array_map (Struct_Type, &init_data_struct,
                       qualifier ("sn_data", datasets));

   variable
     num_bins = length (data[0].value),
     rebin_flags = Integer_Type[num_bins];

   variable
     bounds = qualifier ("bounds", [0]),
     units = qualifier ("unit", unit_default()),
     energy_ordered_bounds = unit_is_energy (units),
     length_bounds = length (bounds);

   bounds = to_angstrom (bounds ;; __qualifiers);

   variable min_sn, min_chan;
   min_sn = force_array_length (Double_Type, length_bounds,
                                qualifier("min_sn", 0.0));

   min_chan = force_array_length (Integer_Type, length_bounds,
                                  qualifier("min_chan",1));

   if (energy_ordered_bounds)
     {
        min_sn = reverse (min_sn);
        min_chan = reverse (min_chan);
     }

   % IMPORTANT -- the output grouping will be different
   % depending on whether an energy unit or a wavelength
   % unit has been provided!
   %
   % Wavelength units --> start grouping at short-wavelength end
   % Energy units --> start grouping at low-energy end

   variable i_inc,
     i_start = 0,
     i_end = num_bins-1,
     bin_edges = data[0].bin_lo;

   if (energy_ordered_bounds)
     {
        bin_edges = data[0].bin_hi;
        (i_start, i_end) = (i_end, i_start);
     }

   i_inc = sign (i_end - i_start);

   if (i_inc == 0)
     {
        vmessage ("*** group:  null grouping specified (no bins are affected)");
        return;
     }

   variable i, j, d, sn, _sign = 1, bin = NULL;

   _for i (i_start, i_end, i_inc)
     {
        rebin_flags[i] = _sign;

        if (bin == NULL)
          bin = struct {ssum=0, bsum=0, bsum_scaled=0, num_chan=0};

        foreach d (data)
          {
             bin.ssum += d.value[i];
             bin.bsum += d.back[i];
             bin.bsum_scaled += d.scaled_back[i];
          }

        bin.num_chan++;

        sn = ((bin.ssum - bin.bsum)
              / sqrt(bin.ssum + bin.bsum_scaled));

        % What min_sn, min_chan applies to this bin?
        variable this_edge = bin_edges[i];
        if (energy_ordered_bounds)
          j = min([where(bounds >= this_edge),length_bounds-1]);
        else
          j = max([where(bounds <= this_edge),0]);

        variable force_new_bin =
          ((energy_ordered_bounds==0 && this_edge < bounds[0])
           || (energy_ordered_bounds==1 && this_edge > bounds[-1]));

        variable new_bin_meets_criteria =
          ((bin.num_chan >= min_chan[j])
           && (sn >= min_sn[j] || min_sn[j] <= 0));

        if (force_new_bin || new_bin_meets_criteria)
          {
             _sign *= -1;
             bin = NULL;
          }
     }

   foreach i (datasets)
     {
        rebin_data(i, rebin_flags);
     }
}

private define matching_grids (datasets)
{
   variable i, n = length(datasets);

   variable k,
     d = Struct_Type[n],
     rebin_flags = Array_Type[n];

   _for i (0, n-1, 1)
     {
        k = datasets[i];
        rebin_flags[i] = get_data_info(k).rebin;
        rebin_data (k, 0);
        d[i] = get_data_counts(k);
     }

   variable mismatch = 0;

   _for i (1, n-1, 1)
     {
        if (any (d[i].bin_lo != d[0].bin_lo or d[i].bin_hi != d[0].bin_hi))
          {
             mismatch = 1;
             break;
          }
     }

   if (mismatch)
     {
        _for i (0, n-1, 1)
          {
             k = datasets[i];
             rebin_data (k, rebin_flags[i]);
          }
        throw ApplicationError, "mismatched dataset grids";
     }

   return d, rebin_flags;
}

%  x and grid are input floating point arrays sorted in ascending order.
%  Their values presumably overlap, but this isn't necessary.
%  They need not be the same size.
%
%  Return x_status = Integer_Type[n], n = length(x) such that:
%
%    x_status[i] = j  to indicate that grid[j] is closest to x[i].
%                  If x[i] falls exactly at the midpoint
%                  between grid[k] and grid[k+1], then j=k.
%    x_status[i] = -1  otherwise.
%
private define closest_alignments (x, grid)
{
   variable n = length(x), num_grid = length(grid),
     grid_index = Integer_Type[n],
     min_offset = Double_Type[n];

   min_offset[*] = _Inf;
   grid_index[*] = -1;

   variable i, j;

   _for i (0, n-1, 1)
     {
        if (x[i] < grid[0])
          {
             grid_index[i] = 0;
             continue;
          }
        else if (x[i] >= grid[-1])
          {
             grid_index[i] = num_grid-1;
             continue;
          }

        j = bsearch (x[i], grid);

        if ((j+1 < num_grid) && (grid[j+1]-x[i] < x[i] - grid[j]))
          j++;

        variable offset = abs(x[i] - grid[j]);

        if (offset < min_offset[i])
          {
             grid_index[i] = j;
             min_offset[i] = offset;
          }
     }

   return grid_index;
}

define group_bin ()
{
   variable msg = "group_bin (datasets, lo, hi [; qualifiers])\n"
     + "   qualifiers:  unit";

   if (_NARGS != 3)
     usage(msg);

   variable datasets, _lo, _hi;
   (datasets, _lo, _hi) = ();

   variable
     lo = to_angstrom (_lo ;; __qualifiers),
     hi = to_angstrom (_hi ;; __qualifiers);

   variable n = length(lo);
   if (length(hi) != length(lo))
     throw ApplicationError, "bin edge arrays have unequal length";

   variable d, old_flags;
   (d, old_flags)= matching_grids (datasets);

   variable alo = closest_alignments (lo, d[0].bin_lo);
   variable ahi = closest_alignments (hi, d[0].bin_hi);
   variable ia = where (alo >= 0 and ahi >= 0);

   alo = alo[ia];
   ahi = ahi[ia];

   variable num_a = length(ia),
     num_bins = length(d[0].bin_lo),
     flags = Integer_Type[num_bins];

   variable k, _sign=1;

   _for k (0, num_a-1, 1)
     {
        if (alo[k] <= ahi[k])
          {
             flags[[ alo[k] : ahi[k] ]] = _sign;
             _sign *= -1;
          }
     }

   variable num_datasets = length(datasets),
     indices = [alo[0]:ahi[-1]];

   _for k (0, num_datasets-1, 1)
     {
        variable new_flags = old_flags[k];
        new_flags[indices] = flags[indices];
        rebin_data (datasets[k], new_flags);
     }
}
