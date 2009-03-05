% -*- mode: SLang; mode: fold -*-
%
%    This file is part of ISIS, the Interactive Spectral Interpretation System
%    Copyright (C) 1998-2008 Massachusetts Institute of Technology
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

% $Id: fits_module_dep.sl,v 1.1 2004/09/09 09:54:21 houck Exp $

require ("fits");

private define do_fits_error (status) %{{{
{
   if (status)
     error (_fits_get_errstatus (status));
}

%}}}

private define add_fits_header_blurb (fp) %{{{
{
   variable blurb, user, host;
   user = getenv ("USER");
   if (user == NULL) user = "unknown";
   host = getenv ("HOST");
   if (host == NULL) host = "unknown";

   blurb = sprintf ("Image created with ISIS version %s", _isis_version_string);
   do_fits_error (_fits_write_history (fp, blurb));
   blurb = sprintf ("by %s@%s, %s", user, host, time);
   do_fits_error (_fits_write_history (fp, blurb));
   do_fits_error (_fits_write_history (fp, "See http://space.mit.edu/cxc/isis for more on ISIS"));
}

%}}}

define save_conf () %{{{
{
   variable msg = "status = save_conf (map, \"file.fits\")";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable file, m;
   (m, file) = ();

   variable img = m.chisqr;

   variable dims, num_dims, data_type;
   (dims, num_dims, data_type) = array_info (img);

   variable fp;
   do_fits_error (_fits_open_file (&fp, file, "c"));
   do_fits_error (_fits_create_img (fp, -32, dims));
   do_fits_error (_fits_write_img (fp, img));

   variable xname, yname, s;

   % this hack allows saving a struct when the data is available
   % but the function isnt defined.
   s = get_par_info (m.px.index);
   if (s != NULL) xname = s.name;
   else xname = "";

   s = get_par_info (m.py.index);
   if (s != NULL) yname = s.name;
   else yname = "";

   variable best_label = sprintf ("Best-fit %s value", Fit_Statistic);

   do_fits_error (_fits_update_key (fp, "BESTSTAT",  m.best, best_label));
   do_fits_error (_fits_update_key (fp, "BEST_X",  m.px_best, "X-coordinate, " + best_label));
   do_fits_error (_fits_update_key (fp, "BEST_Y",  m.py_best, "Y-coordinate, " + best_label));
   do_fits_error (_fits_update_key (fp, "PXNAME", xname, ""));
   do_fits_error (_fits_update_key (fp, "PXMIN", m.px.min, ""));
   do_fits_error (_fits_update_key (fp, "PXMAX", m.px.max, ""));
   do_fits_error (_fits_update_key (fp, "PXNUM", m.px.num, ""));
   do_fits_error (_fits_update_key (fp, "PYNAME", yname, ""));
   do_fits_error (_fits_update_key (fp, "PYMIN", m.py.min, ""));
   do_fits_error (_fits_update_key (fp, "PYMAX", m.py.max, ""));
   do_fits_error (_fits_update_key (fp, "PYNUM", m.py.num, ""));

   variable fun = get_fit_fun();
   if (fun == NULL) fun = "NONE";
   variable fun_comment = sprintf ("fit_fun (\"%s\")", fun);
   % I suspect its safer to write function as a comment
   % because it might sometimes be too long for a keyword.
   do_fits_error (_fits_write_comment (fp, fun_comment));

   variable dx, dy;
   dx = (m.px.max-m.px.min)/double(m.px.num);
   dy = (m.py.max-m.py.min)/double(m.py.num);

   do_fits_error (_fits_update_key (fp, "CTYPE1P", xname, ""));
   do_fits_error (_fits_update_key (fp, "CRVAL1P", m.px.min, ""));
   do_fits_error (_fits_update_key (fp, "CRPIX1P", 0.5, ""));
   do_fits_error (_fits_update_key (fp, "CDELT1P", dx, ""));
   do_fits_error (_fits_update_key (fp, "WCSTY1P", "PHYSICAL", ""));
   do_fits_error (_fits_update_key (fp, "CUNIT1P", "", ""));
   if (dx != 0.0 and dy != 0.0)
     {
        do_fits_error (_fits_update_key (fp, "LTV1",   1.0 - m.px.min/dx, ""));
        do_fits_error (_fits_update_key (fp, "LTM1_1", 1.0/dx, ""));
     }

   do_fits_error (_fits_update_key (fp, "CTYPE2P", yname, ""));
   do_fits_error (_fits_update_key (fp, "CRVAL2P", m.py.min, ""));
   do_fits_error (_fits_update_key (fp, "CRPIX2P", 0.5, ""));
   do_fits_error (_fits_update_key (fp, "CDELT2P", dy, ""));
   do_fits_error (_fits_update_key (fp, "WCSTY2P", "PHYSICAL", ""));
   do_fits_error (_fits_update_key (fp, "CUNIT2P", "", ""));

   if (dx != 0.0 and dy != 0.0)
     {
        do_fits_error (_fits_update_key (fp, "LTV2",   1.0 -m.py.min/dy, ""));
        do_fits_error (_fits_update_key (fp, "LTM2_2", 1.0/dy, ""));
     }

   add_fits_header_blurb (fp);

   _fits_close_file (fp);
}

%}}}

define load_conf () %{{{
{
   variable msg = "m = load_conf (\"file.fits\")";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable file = ();

   variable fp;
   do_fits_error (_fits_open_file (&fp, file, "r"));

   variable m = struct
     {
	chisqr, px, py, best, px_best, py_best
     };

   m.px = conf_grid (0,0,0,0);
   m.py = conf_grid (0,0,0,0);

   variable tmp, status;

   do_fits_error (_fits_read_key (fp, "BESTSTAT", &tmp, NULL));
   m.best = tmp;
   do_fits_error (_fits_read_key (fp, "BEST_X", &tmp, NULL));
   m.px_best = tmp;
   do_fits_error (_fits_read_key (fp, "BEST_Y", &tmp, NULL));
   m.py_best = tmp;

   status = _fits_read_key (fp, "PXNAME", &tmp, NULL);
   if (status == _FITS_KEY_NO_EXIST)
     tmp = "";
   m.px.index = tmp;

   do_fits_error (_fits_read_key (fp, "PXMIN", &tmp, NULL));
   m.px.min = tmp;
   do_fits_error (_fits_read_key (fp, "PXMAX", &tmp, NULL));
   m.px.max = tmp;
   do_fits_error (_fits_read_key (fp, "PXNUM", &tmp, NULL));
   m.px.num = tmp;

   status = _fits_read_key (fp, "PYNAME", &tmp, NULL);
   if (status == _FITS_KEY_NO_EXIST)
     tmp = "";
   m.py.index = tmp;

   do_fits_error (_fits_read_key (fp, "PYMIN", &tmp, NULL));
   m.py.min = tmp;
   do_fits_error (_fits_read_key (fp, "PYMAX", &tmp, NULL));
   m.py.max = tmp;
   do_fits_error (_fits_read_key (fp, "PYNUM", &tmp, NULL));
   m.py.num = tmp;

   do_fits_error (_fits_read_img (fp, &tmp));
   m.chisqr = tmp;

   do_fits_error(_fits_close_file (fp));

   return m;
}

%}}}

% Contributed by Mike Nowak <mnowak@space.mit.edu>
define regroup_file() %{{{
{
   variable msg = "regroup_file (id [, file])  OR  regroup_file (grp[], file)";
   variable id, file = NULL;

   if (_isis->get_varargs (&id, &file, _NARGS, 1, msg))
     return;

   variable grp;
   if (typeof(id) == Int_Type)
     {
        variable info = get_data_info (id);
        grp = i2x_group (info.rebin);
        if (file == NULL)
          {
             file = info.file;
          }
     }
   else
     {
        grp = id;
        if (file == NULL)
          {
             vmessage ("Usage:  %s", msg);
             throw RunTimeError;
             return;
          }
     }

   variable hdu = "[SPECTRUM]";
   if (andelse {0 == is_substr (file, hdu)}
         {0 == is_substr (file, strlow(hdu))})
     {
        file += hdu;
     }

   variable fp = fits_open_file (file, "w");
   if (fp == NULL)
     {
        vmessage ("failed opening %s", file);
        throw OpenError;
     }

   variable colnum, colname = "grouping";
   if (fits_binary_table_column_exists (file, colname))
     {
        do_fits_error(_fits_get_colnum(fp, colname, &colnum));
     }
   else
     {
        variable num_cols;
        do_fits_error(_fits_get_num_cols (fp, &num_cols));
        colnum = num_cols+1;
        do_fits_error(_fits_insert_cols (fp, colnum, colname, "1I"));
     }

   do_fits_error(_fits_write_col(fp, colnum, 1, 1, grp));
   fits_close_file(fp);
}

%}}}

% Contributed by Mike Nowak <mnowak@space.mit.edu>
define use_file_group() %{{{
{
   variable msg = "use_file_group (id [, file]);";
   variable id, file=NULL;

   if (_isis->get_varargs (&id, &file, _NARGS, 1, msg))
     return;

   if (file == NULL)
     {
        file = get_data_info(id).file;
     }

   variable colname = "grouping";
   if (fits_binary_table_column_exists (file, colname))
     {
        variable grp = fits_read_col (file, colname);
        rebin_data (id, x2i_group(grp));
     }
}

%}}}

% APED database bibliography retrieval

require ("readascii");

private define expand_string (s) %{{{
{
   return eval(sprintf ("\"%s\"$", s));
}

%}}}

private define load_filemap (s) %{{{
{
   variable f = struct
     {
        type, Z, ion, file
     };

   variable filemap_file = path_concat (s.dir, s.atomic_data_filemap);
   variable nlines = readascii (filemap_file,
                                &f.type, &f.Z, &f.ion, &f.file;
                                format="%d %d %d %s");
   f.file = array_map (String_Type, &expand_string, f.file);

   return f;
}

%}}}

private variable Warn_Duplicates = 0;
private define warn_dup (who, num, up, lo, file) %{{{
{
   if (Warn_Duplicates)
     {
        vmessage ("$who:  $num entries for $up -> $lo  $file"$);
     }
}

%}}}

private define load_dr_refs1 (upper, lower, dr_file) %{{{
{
   variable r = fits_read_col_struct (dr_file,
                                      ["upper_lev", "lower_lev",
                                       "drrate_ref", "wave_ref", "wv_obs_ref"]);

   variable refs = struct {dr_drrate_ref, dr_wave_ref, dr_wv_obs_ref};
   variable num_tr = length(upper);
   refs.dr_drrate_ref = String_Type[num_tr];
   refs.dr_wave_ref = String_Type[num_tr];
   refs.dr_wv_obs_ref = String_Type[num_tr];

   refs.dr_drrate_ref[*] = "";
   refs.dr_wave_ref[*] = "";
   refs.dr_wv_obs_ref[*] = "";

   _for (0, num_tr-1, 1)
     {
        variable j = ();
        variable k = where (r.upper_lev == upper[j]
                            and r.lower_lev == lower[j]);
        if (length(k) == 0)
          continue;
        if (length(k) != 1)
          {
             warn_dup ("DR1", length(k), upper[j], lower[j], dr_file);
          }

        variable ik = k[0];
        refs.dr_drrate_ref[j] = r.drrate_ref[ik];
        refs.dr_wave_ref[j] = r.wave_ref[ik];
        refs.dr_wv_obs_ref[j] = r.wv_obs_ref[ik];
     }

   return refs;
}

%}}}

private define load_dr_refs2 (upper, lower, dr_file) %{{{
{
   variable r = fits_read_col_struct (dr_file, ["reference"]);

   variable levs = [upper, lower];
   levs = levs[unique(levs)];

   variable refs = struct {dr_reference};
   variable num_levs = length(levs);
   refs.dr_reference = String_Type[num_levs];
   num_levs.dr_reference[*] = "";

   variable j;
   _for j (0, num_levs-1, 0)
     {
        refs.dr_reference[j] = r.reference[levs[j]-1];
     }

   return refs;
}

%}}}

private define load_dr_refs (upper, lower, dr_file) %{{{
{
   variable num_tr = length(upper);

   if (fits_binary_table_column_exists (dr_file, "wv_obs_ref"))
     {
        return load_dr_refs1 (upper, lower, dr_file);
     }
   else
     {
        return load_dr_refs2 (upper, lower, dr_file);
     }
}

%}}}

private define load_ec_refs (upper, lower, ec_file) %{{{
{
   variable r = fits_read_col_struct (ec_file,
                                      ["upper_lev", "lower_lev", "reference"]);

   variable num_tr = length(upper);
   variable refs = struct {ec_reference};
   refs.ec_reference = String_Type[num_tr];
   refs.ec_reference[*] = "";

   _for (0, num_tr-1, 1)
     {
        variable j = ();
        variable k = where (r.upper_lev == upper[j]
                            and r.lower_lev == lower[j]);
        if (length(k) == 0)
          continue;
        if (length(k) > 1)
          {
             warn_dup ("EC", length(k), upper[j], lower[j], ec_file);
          }

        variable s = "";
        if (length(k) == 1) s = r.reference[k[0]];
        refs.ec_reference[j] = s;
     }

   return refs;
}

%}}}

private define load_lv_refs1 (levs, elev_file) %{{{
{
   variable r = fits_read_col_struct (elev_file,
                                      ["energy_ref", "phot_ref"]);

   variable num_table_levs = length(r.energy_ref);
   variable i = where (levs < num_table_levs);
   levs = levs[i];

   variable num_levs = length(levs);

   variable refs = struct {lv_energy_ref, lv_phot_ref};
   refs.lv_energy_ref = String_Type[num_levs];
   refs.lv_phot_ref = String_Type[num_levs];

   variable j;
   _for j (0, num_levs-1, 1)
     {
        variable k = levs[j]-1;
        refs.lv_energy_ref[j] = r.energy_ref[k];
        refs.lv_phot_ref[j] = r.phot_ref[k];
     }

   return refs;
}

%}}}

private define load_lv_refs2 (levs, elev_file) %{{{
{
   variable r = fits_read_col_struct (elev_file, ["reference"]);

   variable num_table_levs = length(r.reference);
   variable i = where (levs < num_table_levs);
   levs = levs[i];

   variable num_levs = length(levs);

   variable refs = struct {lv_reference};
   refs.lv_reference = String_Type[num_levs];

   variable j;
   _for j (0, num_levs-1, 1)
     {
        variable k = levs[j]-1;
        refs.lv_reference[j] = r.reference[k];
     }

   return refs;
}

%}}}

private define load_lv_refs (levs, elev_file) %{{{
{
   if (fits_binary_table_column_exists (elev_file, "energy_ref"))
     {
        return load_lv_refs1 (levs, elev_file);
     }
   else
     {
        return load_lv_refs2 (levs, elev_file);
     }
}

%}}}

private define load_la_refs1 (upper, lower, wave_file) %{{{
{
   variable x = fits_read_col_struct (wave_file,
                                      ["upper_lev", "lower_lev",
                                       "wave_ref", "wv_obs_ref", "ein_a_ref"]);

   variable refs = struct
     {
        la_wave_ref, la_wv_obs_ref, la_ein_a_ref
     };
   variable num_tr = length(upper);
   refs.la_wave_ref = String_Type[num_tr];
   refs.la_wv_obs_ref = String_Type[num_tr];
   refs.la_ein_a_ref = String_Type[num_tr];

   refs.la_wave_ref[*] = "";
   refs.la_wv_obs_ref[*] = "";
   refs.la_ein_a_ref[*] = "";

   variable j;
   _for j (0, num_tr-1, 1)
     {
        variable i_ref = where (x.upper_lev == upper[j]
                                and x.lower_lev == lower[j]);
        if (length(i_ref) != 1)
          {
             warn_dup("LA1", length(i_ref), upper[j], lower[j], wave_file);
          }
        i_ref = i_ref[0];

        refs.la_wave_ref[j] = x.wave_ref[i_ref];
        refs.la_wv_obs_ref[j] = x.wv_obs_ref[i_ref];
        refs.la_ein_a_ref[j] = x.ein_a_ref[i_ref];
     }

   return refs;
}

%}}}

private define load_la_refs2 (upper, lower, wave_file) %{{{
{
   variable x = fits_read_col_struct (wave_file,
                                      ["upper_lev", "lower_lev", "reference"]);

   variable refs = struct {la_reference};
   variable num_tr = length(upper);
   refs.la_reference = String_Type[num_tr];
   refs.la_reference[*] = "";

   variable j;
   _for j (0, num_tr-1, 1)
     {
        variable i_ref = where (x.upper_lev == upper[j]
                                and x.lower_lev == lower[j]);
        if (length(i_ref) != 1)
          {
             warn_dup ("LA2", length(i_ref), upper[j], lower[j], wave_file);
          }
        i_ref = i_ref[0];

        refs.la_reference[j] = x.reference[i_ref];
     }

   return refs;
}

%}}}

private define load_la_refs (upper, lower, wave_file) %{{{
{
   if (fits_binary_table_column_exists (wave_file, "wave_ref"))
     {
        return load_la_refs1 (upper, lower, wave_file);
     }
   else
     {
        return load_la_refs2 (upper, lower, wave_file);
     }
}

%}}}

private define load_pc_refs (upper, lower, pc_file) %{{{
{
   variable r = fits_read_col_struct (pc_file,
                                      ["upper_lev", "lower_lev", "reference"]);

   variable num_tr = length(upper);
   variable refs = struct {pc_reference};
   refs.pc_reference = String_Type[num_tr];
   refs.pc_reference[*] = "";

   variable j;
   _for j (0, num_tr-1, 1)
     {
        variable k = where (r.upper_lev == upper[j]
                            and r.lower_lev == lower[j]);
        if (length(k) == 0)
          continue;
        if (length(k) != 1)
          {
             warn_dup ("PC", length(k), upper[j], lower[j], pc_file);
          }

        refs.pc_reference[j] = r.reference[k[0]];
     }

   return refs;
}

%}}}

private define struct_append_strings (a, b) %{{{
{
   variable na, nb;

   foreach nb (get_struct_field_names (b))
     {
        variable fb = get_struct_field (b, nb);
        if (a == NULL || 0 == struct_field_exists (a, nb))
          {
             a = struct_combine (a, nb);
             set_struct_field (a, nb, fb);
          }
        else
          {
             variable fa = get_struct_field (a, nb);
             set_struct_field (a, nb, [fa, fb]);
          }
     }

   return a;
}

%}}}

private define unique_strings (ss) %{{{
{
   if (typeof(ss) == String_Type)
     return ss;

   variable i = array_sort (ss);
   variable s = ss[i];
   s = array_map (String_Type, &strtrim, s);

   variable num_strings = length(s);

   variable k = 0;
   while (k < num_strings)
     {
        if (s[k] != "")
          break;
        k += 1;
     }

   if (k == num_strings)
     return NULL;

   variable n=0, u=[ s[k] ];

   _for i (k, num_strings-1, 1)
     {
        variable t = s[i];
        if (t != "" && t != u[n])
          {
             u = [u, t];
             n += 1;
          }
     }

   return u;
}

%}}}

private define filter_unique (refs) %{{{
{
   foreach (get_struct_field_names(refs))
     {
        variable n = ();
        variable s = unique_strings (get_struct_field (refs, n));
        set_struct_field (refs, n, s);
     }
}

%}}}

define aped_bib () %{{{
{
   variable msg =
     "Struct_Type = aped_bib (Struct_Type aped, Integer_Type[] lines)\n"
     + " qualifiers:  verbose | unique | warn_dup";
   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable db_struct, lines;
   (db_struct, lines) = ();

   if (length(lines) == 0)
     return;

   Warn_Duplicates = qualifier_exists ("warn_dup");

   if (qualifier_exists ("verbose"))
     () = fprintf (stderr, "looking up references...");

   variable fm = load_filemap (db_struct);
   variable info = array_map (Struct_Type, &line_info, lines);

   variable n = length(lines);
   variable refs = NULL;

   variable
     z_list = array_struct_field (info, "Z"),
     ion_list = array_struct_field (info, "ion"),
     upper = array_struct_field (info, "upper"),
     lower = array_struct_field (info, "lower");

   variable species = z_list * 100 + ion_list;
   variable unique_species = species[unique(species)];

   variable us;
   foreach us (unique_species)
     {
        variable
          z = us / 100,
          i = us - 100 * z;

        variable izi = where (fm.Z == z and fm.ion == i);

        if (length(izi) == 0)
          continue;

        variable
          us_lines = where (z_list == z and ion_list == i),
          up = upper[us_lines],
          lo = lower[us_lines];

        variable ifile, num_files = length(izi);

        _for ifile (0, num_files-1, 1)
          {
             variable
               k = izi[ifile],
               type = fm.type[k],
               file = fm.file[k];

             variable r;
             switch (type)
               {
                case 2:
                  variable levs = [up, lo];
                  levs = levs[unique(levs)];
                  r = load_lv_refs (levs, file);
               }
               {case 3: r = load_la_refs (up, lo, file);}
               {case 4: r = load_ec_refs (up, lo, file);}
               {case 5: r = load_pc_refs (up, lo, file);}
               {case 6: r = load_dr_refs (up, lo, file);}
               {
                  % default:
                  vmessage ("unimplemented type: $type file: $file"$);
               }

             refs = struct_append_strings (refs, r);
          }
     }

   if (qualifier_exists ("verbose"))
     () = fprintf (stderr, "done\n");

   if (qualifier_exists ("unique"))
     filter_unique (refs);

   return refs;
}

%}}}

define aped_bib_query_string () %{{{
{
   variable msg =
     "String_Type = aped_bib_query_string (Struct_Type refs)\n"
     + "  qualifiers:   host=, type=BIBTEX|ABSTRACT..";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable refs = ();

   variable r = @refs;
   filter_unique (r);

   variable
     url = qualifier ("host", "adsabs.harvard.edu"),
     type = qualifier ("type", "BIBTEX");

   variable v_array = String_Type[0];

   foreach (get_struct_field_names (r))
     {
        variable n = ();
        variable v = get_struct_field (r, n);
        if (v != NULL)
          v_array = [v_array, v];
     }

   if (length(v_array) == 0)
     return NULL;

   variable bibcodes = strjoin ("bibcode=" + v_array, "&");

   return "http://${url}/cgi-bin/nph-bib_query/?${bibcodes}&data_type=${type}&db_key=AST&nocookieset=1"$;
}

%}}}

% - end - APED database bibliography retrieval code
