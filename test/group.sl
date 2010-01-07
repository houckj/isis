define generate_fake_data (num_datasets, num_bins, cvalue, bvalue)
{
   variable lo, hi, n=num_bins,
     emin = 1.0, emax=10.0,
     c = Double_Type[n],
     b = Double_Type[n];

   (lo, hi) = linear_grid (emin, emax, n);
   (lo, hi) = _A(lo, hi);

   c[*] = cvalue;
   b[*] = bvalue;

   loop (num_datasets)
     {
        variable id = define_counts (lo, hi, c, sqrt(c));
        () = _define_back (id, b);
     }
}

define compute_sn (c, b, scale)
{
   % b already contains one factor of 'scale'
   return (c-b) / sqrt (c + b * scale);
}

define compute_dataset_sn (i)
{
   variable
     d = get_data_counts(i),
     b = get_back (i),
     scale = ((get_data_exposure(i) * get_data_backscale(i))
              /(get_back_exposure(i) * get_back_backscale(i)));

   return compute_sn (d.value, b, scale);
}

define check_sn_bounds (k, bounds, is_egy, sn_values)
{
   variable d = get_data_counts (k);
   if (is_egy) d = _A(d);

   variable i, num_values = length(sn_values),
     expected_sn = Double_Type[length(d.value)];

   expected_sn[where (d.bin_lo <= bounds[0])] = sn_values[0];
   _for i (1, num_values-2, 1)
     {
        variable r = where (bounds[i-1] <= d.bin_hi and d.bin_lo <= bounds[i]);
        expected_sn[r] = sn_values[i];
     }
   expected_sn[where (bounds[-1] <= d.bin_lo)] = sn_values[-1];
   variable sn = compute_dataset_sn (1);
   if (is_egy) sn = reverse(sn);

   i = where (sn != expected_sn);

   % the expected values may not exactly match at sub-region
   % boundaries because the expected value calculation doesn't
   % include the associated "edge" effects.
   if (length(i) > length(bounds))
     throw ApplicationError;
}

define isis_main()
{
   () = fputs ("testing group... ", stderr);
   variable d, expected_sn;

   %
   % Construct fake data sets that yield simple S/N tests
   %

   variable num_bins = 32;
   generate_fake_data (1, num_bins, 9, 0);
   % test min_sn alone, no background
   group (1; min_sn=3);
   if (any (compute_dataset_sn (1) != 3))
     throw ApplicationError;
   % test min_chan alone, no background
   group (1; min_chan=4);
   if (length(get_data_counts(1).value) != num_bins/4)
     throw ApplicationError;

   delete_data (all_data);
   generate_fake_data (1, num_bins, 10, 1);
   % test min_sn alone, with background
   group (1; min_sn=5);
   expected_sn = compute_sn (4*10, 4*1, 1);
   if (any(expected_sn != compute_dataset_sn (1)))
     throw ApplicationError;

   delete_data (all_data);
   generate_fake_data (4, num_bins, 2, 0);
   % test min_sn alone, no background, multiple datasets
   group (all_data; min_sn=3);
   if (length(get_data_counts(1).value) != num_bins/2)
     throw ApplicationError;

   delete_data (all_data);
   generate_fake_data (1, 1024, 1.5^2, 0);
   % test min_sn with energy-bounded subregions, no background
   variable bounds, sn_values, sn, i;
   bounds = [4, 6, 8];
   sn_values = [1.5, 1.5*sqrt(2), 3];
   group (1; bounds=bounds, unit="keV", min_sn=sn_values);
   check_sn_bounds (1, bounds, 1, [1.5, sn_values]);

   % test min_sn with wavelength-bounded subregions, no background
   bounds = _A(bounds);
   sn_values = reverse(sn_values);
   group (1; bounds=bounds, unit="a", min_sn=sn_values);
   check_sn_bounds (1, bounds, 0, [1.5, sn_values]);

   vmessage ("ok");
   () = fputs ("testing group_bin... ", stderr);

   delete_data (all_data);
   variable lo, hi, n=1000, c = Double_Type[n];
   (lo, hi) = linear_grid (1, 10, n);
   c[*] = prand(100, n) + 25 * sin (lo);
   () = define_counts (lo, hi, c, sqrt(c));

   variable k, num= 2 + int(log(n)/log(2));
   _for k (0, num-1, 1)
     {
        variable new_lo, new_hi, new_num_bins = 1 shl k;
        (new_lo, new_hi) = linear_grid (1, 10, new_num_bins);
        group_bin (1, new_lo, new_hi ; unit="a");

        variable new = get_data_counts (1);

        % if the new grid is finer than the old one, there
        % may not be any grouping at all.
        if (new_num_bins >= n && length(new.bin_lo) == n)
          break;

        % verify that the grouped bin edges really are
        % the ones closest to the requested bin edges.
        variable ll;
        foreach ll (new_lo)
          {
             variable correct_min_diff = min(abs(lo - ll));
             variable min_diff = min(abs(new.bin_lo - ll));
             if (min_diff != correct_min_diff)
               throw ApplicationError;
          }
     }

   vmessage ("ok");
}
