% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2010 Massachusetts Institute of Technology
%
%    This software was developed by the MIT Center for Space Research under
%    contract SV1-61010 from the Smithsonian Institution.
%
%    Author:  John C. Houck  <houck@space.mit.edu>
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

%  $Id: emis-cmds.sl,v 1.6 2004/06/06 04:25:23 houck Exp $

define plasma () %{{{
{
   variable msg = "plasma ([db_struct [, density_range[2] [, temp_range[2] ]]] )";
   variable trange = Float_Type[2];
   variable drange = Float_Type[2];
   variable s = NULL;

   trange = [0.0, 1.e9];        % Kelvin
   drange = [0.0, 1.e30];       % cm^-3

   if (_isis->get_varargs (&s, &drange, &trange, _NARGS, 0, msg))
     return;

   if (s == NULL)
     {
	usage(msg);
	return;
     }

   _isis->Dbase = s;

   if (NULL == _isis->get_atomic_db_pointer())
     {
        atoms (s);
        if (NULL == _isis->get_atomic_db_pointer())
          throw ApplicationError, "Failed loading atomic database";
     }

   variable em = _isis->_em_start (min(trange), max(trange), min(drange), max(drange));
   _isis->set_emis_db_pointer (em);
}
%}}}

define db_push () %{{{
{
   variable msg = "db_push (Struct_Type)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable s = ();

   if (_isis->db_push(s))
     {
        atoms(s);
        plasma(s);
     }
}

%}}}

define db_pop () %{{{
{
   variable msg = "Struct_Type = db_pop ([k])";
   variable k = NULL;
   if (_NARGS == 1)
     k = ();
   else if (_NARGS > 1)
     usage(msg);
   _isis->db_pop (k);
}

%}}}

define db_indices () %{{{
{
   _isis->db_indices();
}

%}}}

define db_list () %{{{
{
   variable args = __pop_list(_NARGS);
   _isis->db_list (__push_list(args));
}

%}}}

define db_select () %{{{
{
   variable msg = "db_select (k)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable k = ();
   _isis->db_select (k);
}

%}}}

private define path_to_emis_file (x) %{{{
{
   variable file;
   if (x == NULL and _isis->Dbase.dir != NULL)
     {
	file = "";
     }
   else if (typeof(x) == Struct_Type)
     {
	file = path_concat (x.dir, x.line_emissivity);
     }
   else file = x;

   return file;
}

%}}}

define db_grid () %{{{
{
   variable msg = "Struct_Type = db_grid ([file | Struct_Type])";
   variable x = NULL;

   if (_isis->get_varargs (&x, _NARGS, 0, msg))
     return;

   variable file = path_to_emis_file (x);

   variable s = struct
     {
	temp, dens
     };

   (s.temp, s.dens) = _isis->_get_filemap (file);

   return s;
}

%}}}

define list_db () %{{{
{
   variable msg = "list_db ([file | Struct_Type])";
   variable x = NULL;

   if (_isis->get_varargs (&x, _NARGS, 0, msg))
     return;

   variable file = path_to_emis_file (x);

   if (file == NULL)
     {
	usage (msg);
	return;
     }

   variable s = struct
     {
	temp, dens
     };

   (s.temp, s.dens) = _isis->_get_filemap (file);

   variable tm, dm;
   tm = moment(s.temp);
   dm = moment(s.dens);

   vmessage ("Temperature/Density Coverage: ");

   if (tm.max > tm.min)
     vmessage ("  kT:  %0.4g -> %0.4g K", tm.min, tm.max);
   else
     vmessage ("  kT:  %0.4g K", tm.max);

   if (dm.max > dm.min)
     vmessage ("   n:  %0.4g -> %0.4g cm^-3", dm.min, dm.max);
   else
     vmessage ("   n:  %0.4g cm^-3", dm.max);
}

%}}}

define vlist_db () %{{{
{
   variable msg = "vlist_db ([file | Struct_Type])";
   variable x = NULL;

   if (_isis->get_varargs (&x, _NARGS, 0, msg))
     return;

   variable file = path_to_emis_file (x);

   if (file == NULL)
     {
	usage (msg);
	return;
     }

   variable s = struct
     {
	temp, dens
     };

   (s.temp, s.dens) = _isis->_get_filemap (file);

   variable p, fp = stdout;
   p = isis_get_pager ();
   if (p != NULL)
     fp = popen (p, "w");

   () = fprintf (fp, "  hdu      kT (K)    n (cm^-3)\n");
   foreach ([0:length(s.temp)-1])
     {
	variable k = ();
	() = fprintf (fp, " %4d  %10.4g   %10.4g\n", k, s.temp[k], s.dens[k]);
     }

   if (fp != stdout)
     () = pclose (fp);
}

%}}}

%{{{ ionization

define ion_bal()
{
   variable msg = "(ion, frac) = ion_bal (proton_number, temp, [ioniz_table_id])";
   variable proton_number, temp, ioniz_table_id;
   variable q, frac;

   ioniz_table_id = 0;

   if (_isis->get_varargs (&proton_number, &temp, &ioniz_table_id, _NARGS, 2, msg))
     return;

   (q, frac) = _isis->_ioniz_bal (proton_number, temp, 0.0, ioniz_table_id);  % FIXME density not used

   return (q + 1, frac);
}

define ion_frac()
{
   variable msg = "frac[] = ion_frac (proton_number, ion, temp[], [ioniz_table_id])";
   variable proton_number, ion, temp, ioniz_table_id;

   ioniz_table_id = 0;

   if (_isis->get_varargs (&proton_number, &ion, &temp, &ioniz_table_id,
		    _NARGS, 3, msg))
     return;

   return _isis->_ioniz_fcn (proton_number, ion-1, temp, ioniz_table_id);
}

%}}}

%{{{ line emissivities and ratios

define line_em1 ()
{
   variable msg = "(emis[], temp[], dens[]) = line_em1 (line_index)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable line_index = ();

   return _isis->_line_emissivities_on_datafile_grid (line_index);
}

define line_em()
{
   variable msg = "em[] = line_em (line[], temp[] [, dens[] ])";
   variable line_list, temp, dens;

   dens = [0.0];

   if (_isis->get_varargs (&line_list, &temp, &dens, _NARGS, 2, msg))
     return;

   if (Array_Type == typeof (line_list))
     {
	variable n = length (line_list);
	variable nd = length (dens);
	variable nt = length (temp);
	variable em;
	variable i;

	if (nd == 1)
	  em = Float_Type [ nt ];
	else if (nt == 1)
	  em = Float_Type [ nd ];
	else if (nd > 1 and nt > 1)
	  em = Float_Type [ nd, nt ];
	else
	  {
	     message ("invalid array dimension");
	     return NULL;
	  }

	em = _isis->_line_em (line_list, temp, dens);

	return em;
     }
   else
     {
	return _isis->_line_em (line_list, temp, dens);
     }
}

define ratio_em ()
{
   variable msg = "ratio_array = ratio_em (l1_array, l2_array, temp_array, [dens_array])";
   variable top, bot, ratio, i, nd, nt;
   variable t, b, temp, dens;

   dens=[0.0];

   if (_isis->get_varargs (&t, &b, &temp, &dens, _NARGS, 3, msg))
     return;

   top = line_em (t, temp, dens);
   bot = line_em (b, temp, dens);

   if (NULL == top or NULL == bot)
     {
	message("Failed to compute emissivity ratio");
	return NULL;
     }

   nd = length (dens);
   nt = length (temp);

   if (nt == 1)
     ratio = Float_Type [ nd ];
   else if (nd == 1)
     ratio = Float_Type [ nt ];
   else
     ratio = Float_Type [ nd, nt ];

   i = where ( bot > 0.0 );

   ratio[i] = top[i] / bot[i];

   if (nt == 1)
     reshape (ratio, [ nd ]);
   else if (nd == 1)
     reshape (ratio, [ nt ]);
   else
     reshape (ratio, [ nd, nt ]);

   return ratio;
}

%}}}

define list_abund () %{{{
{
   variable msg = "list_abund ([verbose])";
   variable verbose = 0;

   if (_isis->get_varargs (&verbose, _NARGS, 0, msg))
     return;

   _isis->_list_abund_tables (verbose);
}

%}}}

private define _abund_table (id) %{{{
{
   if (typeof(id) == String_Type)
     id = _isis->_index_for_abundance_table (id);

   return id;
}

%}}}

define set_abund () %{{{
{
   variable msg = "set_abund (index | name)";
   variable id;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   id = ();

   _isis->_choose_abs_abundance_table (_abund_table(id));
}

%}}}

define get_abundances () %{{{
{
   variable msg = "Struct_Type = get_abundances ([index|name])";
   variable id = -1;

   if (_isis->get_varargs (&id, _NARGS, 0, msg))
     return;

   variable t = _isis->_get_abundance_table (_abund_table(id));

   if (t.abun == NULL)
     t = NULL;
   return t;
}

%}}}

define add_abundances () %{{{
{
   variable msg = "id = add_abundances (Struct_Type | name, abun[], z[])";
   variable name, abun, z, t;

   if (_isis->get_varargs (&name, &abun, &z, _NARGS, 1, msg))
     return -1;

   if (Struct_Type == typeof(name))
     t = name;
   else
     {
	t = struct {abun, z, name};
	t.name = name;
	t.abun = abun;
	t.z = z;
     }

   return _isis->_add_abund_table (t);
}

%}}}

define free_alt_ioniz () %{{{
{
   _isis->_free_alt_ionization_table ();
}

%}}}

define load_alt_ioniz () %{{{
{
   variable msg = "load_alt_ioniz (file)";
   variable file;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   file = ();

   _isis->_load_alt_ionization_table (file);
}

%}}}

%{{{ continuum emissivity

define get_contin ()
{
   variable msg = "Struct_Type p = get_contin (lo, hi, temp, [dens], [Z], [ion])";
   variable lo, hi, temp, dens, Z, ion;

   dens = 1.0;     % low density by default.
   Z = 0;          % default is summed over elements and ions.
   ion = -1;

   if (_isis->get_varargs (&lo, &hi, &temp, &dens, &Z, &ion, _NARGS, 3, msg))
     return;

   variable p = struct { true, pseudo };

   (p.true, p.pseudo)
     = _isis->_get_continuum (lo, hi, temp, dens, Z, ion-1);

   return p;
}

%}}}
