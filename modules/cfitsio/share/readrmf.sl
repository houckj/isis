require ("histogram");
require ("fits");
static define convert_to_array_type (a)
{
   if (_typeof (a) == Array_Type)
     return;
   
   variable dims, ndims, type;
   (dims, ndims, type) = array_info (a);
   variable nrows = dims[0];
   variable aa = Array_Type[nrows];
   _for (0, nrows-1, 1)
     {
	variable i = ();
	aa[i] = a[i, *];
     }
   return aa;
}


public define fits_read_rmf (file)
{
   variable rmf = struct
     {
        energ_lo, energ_hi, n_grp, f_chan, n_chan, matrix,
	e_min, e_max, channel
     };
   
   variable t = fits_read_table (file + "[MATRIX]");
   rmf.energ_lo = t.energ_lo;
   rmf.energ_hi = t.energ_hi;
   rmf.n_grp = t.n_grp;
   
   variable f_chan = t.f_chan;
   variable n_chan = t.n_chan;
   variable matrix = t.matrix;
   variable nrows = length (rmf.energ_lo);
   
   f_chan = convert_to_array_type (f_chan);
   n_chan = convert_to_array_type (n_chan);

   _for (0, nrows-1, 1)
     {
	variable i = ();
	variable a;
	a = f_chan[i];
	reshape (a, [length (a)]);
	a = n_chan[i];
	reshape (a, [length (a)]);
	a = matrix[i];
	reshape (a, [length (a)]);
     }
   rmf.f_chan = f_chan;
   rmf.n_chan = n_chan;
   rmf.matrix = matrix;

   (rmf.e_min, rmf.e_max, rmf.channel) = fits_read_col (file+"[EBOUNDS]", "e_min","e_max","channel");
   rmf.channel = int (rmf.channel + 0.5);
   return rmf;
}

public define eval_rmf_E(rmf, e)
{
   variable channels = rmf.channel;
   variable num_channels = max (channels);
   variable num_energies = length (e);
   variable i_list = hist_bsearch (e, rmf.energ_lo);

   variable matrix = rmf.matrix;
   variable r = @Array_Type(_typeof(matrix[0]), [num_energies, num_channels]);

   variable min_chan = min (channels);
   _for (0, num_energies-1, 1)
     {
	variable ii = ();
	variable i = i_list[ii];
	variable n_chan = rmf.n_chan[i];
	variable f_chan = rmf.f_chan[i];
	variable matrix_i = matrix[i];

	variable matrix_start = 0;
	_for (0, rmf.n_grp[i]-1, 1)
	  {
	     variable j = ();
	     variable n = n_chan[j];
	     variable start_chan = f_chan[j];
	     variable stop_chan = start_chan + (n-1);
	     variable chans = [start_chan:stop_chan];
	     variable matrix_stop = matrix_start + (n-1);
	     r[ii,chans] = matrix_i[[matrix_start:matrix_stop]];
	     matrix_start = matrix_stop+1;
	  }
     }
   
   if (typeof(e) != Array_Type)
     reshape (r, [num_channels]);

   return r;
}
