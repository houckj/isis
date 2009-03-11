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

% $Id: hist-cmds.sl,v 1.34 2004/06/06 22:23:16 houck Exp $
require ("structfuns");

private variable _V = struct 
{
   data_counts, 
   data_flux, 
   model_counts,
   model_flux,
   convolved_model_flux
};

% Setting various bits.
_V.model_counts         = 0;
_V.data_counts          = _isis->_Data;
_V.data_flux            = _isis->_Data | _isis->_Flux;
_V.model_flux           = _isis->_Flux;
_V.convolved_model_flux = _isis->_Flux | _isis->_Convolved;

define get_data_counts ();
define get_data_exposure ();
define get_data_info ();
define all_data ();
private define massage_data_info_struct (s) %{{{
{
   if (s.grating != "")
     {
        s.instrument = sprintf ("%s-%s", s.grating, s.instrument);
     }

   return strlen (s.instrument);
}

%}}}

private define dataset_string (s, i, width) %{{{
{
   variable fmt =
     "%3u %c %*s %3d %2d  %2d   %5d/%5d  %s  %s  %10.4e   %7.3f  %s";

   variable d = get_data_counts(i);

   variable arf_id;
   if (length(s.arfs) > 1)
     arf_id = " L";
   else if (s.arfs[0] == 0)
     arf_id = " -";
   else arf_id = sprintf ("%2d", s.arfs[0]);

   variable rmf_id;
   if (length(s.rmfs) > 1)
     rmf_id = " L";
   else if (s.rmfs[0] == 0)
     rmf_id = " -";
   else rmf_id = sprintf ("%2d", s.rmfs[0]);

   variable exclude_char;
   if (s.exclude == 0)
     exclude_char = ' ';
   else exclude_char = 'x';

   variable str;
   str = sprintf (fmt, i, exclude_char, width, s.instrument,
                  s.order, s.part, s.srcid,
                  length(s.notice_list), length(s.notice),
                  arf_id, rmf_id,
                  sum(d.value), get_data_exposure(i)/1.e3,
                  s.target);

   if (Isis_List_Filenames)
     {
        str += sprintf ("\nfile:  %s", s.file);

        variable back = _isis->_get_instrumental_background_hook_name(i);
        if (back == NULL) back = s.bgd_file;
        if (back != "")
          {
             str += sprintf ("\nback:  %s", back);
          }
     }

   if (2 == is_defined ("list_data_hook"))
     {
        Isis_Active_Dataset = i;
        eval ("list_data_hook");
        variable sh = ();
        str += sprintf ("\n%s", sh);
     }

   return str;
}

%}}}

private define list_data_string () %{{{
{
   variable ids = all_data();
   if (ids == NULL)
     return NULL;

   variable info;
   info = array_map (Struct_Type, &get_data_info, ids);

   variable width, min_width = 11;
   width = array_map (Int_Type, &massage_data_info_struct, info);
   width = max([width, min_width]);

   variable title = "Current Spectrum List:";
   variable blanks = "";
   if (width > min_width)
     {
        blanks = String_Type[width-min_width];
        blanks[*] = ' ';
     }
   variable colheads = " id " + blanks +
     "   instrument   m prt src    use/nbins   A   R     totcts   exp(ksec)  target";

   variable s;
   s = array_map (String_Type, &dataset_string, info, ids, width);

   return strjoin ([title, colheads, s], "\n");
}

%}}}

define list_data () %{{{
{
   variable dest = NULL;
   if (_NARGS) dest = ();
   _isis->put_string (dest, _function_name, list_data_string());
}

%}}}

define get_arf ();
define get_arf_info ();
define all_arfs ();
private define arf_string (i) %{{{
{
   variable fmt = "%3u %7s %8s %3d %2d  %2d    %5d  %7.2f  %s";

   variable arf = get_arf(i);
   variable a = get_arf_info(i);

   variable str;
   str = sprintf (fmt, i, a.grating, a.instrument,
                  a.order, a.part, a.srcid,
                  length(arf.value), a.exposure/1.e3, a.object);

   if (Isis_List_Filenames)
     {
        str += sprintf ("\nfile:  %s", a.file);
     }

   return str;
}

%}}}

private define list_arf_string () %{{{
{
   variable ids = all_arfs();
   if (ids == NULL)
     return NULL;
   variable title = "Current ARF List:";
   variable colheads = 
     " id grating detector   m prt src   nbins  exp(ksec)  target";
   variable s = array_map (String_Type, &arf_string, ids);
   return strjoin ([title, colheads, s], "\n");
}

%}}}

define list_arf () %{{{
{
   variable dest = NULL;
   if (_NARGS) dest = ();
   _isis->put_string (dest, _function_name, list_arf_string());
}

%}}}

define get_rmf_info ();
define all_rmfs ();
private define rmf_string (i) %{{{
{
   variable fmt = "%3u %7s %8s  %3d %s  %s";
   variable r = get_rmf_info(i);

   variable type;
   switch (r.method)
     { case _isis->RMF_FILE: type = "file:"; }
     { case _isis->RMF_USER: type = "user:"; }
     { case _isis->RMF_DELTA: type = "delta:"; }
     {
        % default
        type = "???:";
     }

   variable str;
   str = sprintf (fmt, i, r.grating, r.instrument, r.order,
                  type, r.arg_string);

   return str;
}

%}}}

private define list_rmf_string () %{{{
{
   variable ids = all_rmfs();
   if (ids == NULL)
     return NULL;
   variable title = "Current RMF List:";
   variable colheads =
     " id grating detector    m type   file";
   variable s = array_map (String_Type, &rmf_string, ids);
   return strjoin ([title, colheads, s], "\n");
}

%}}}

define list_rmf () %{{{
{
   variable dest = NULL;
   if (_NARGS) dest = ();
   _isis->put_string (dest, _function_name, list_rmf_string());
}

%}}}

% data/arf/rmf input

private define load_data_file (load, num, msg) %{{{
{
   variable file;

   if (_isis->chk_num_args (num, 1, msg))
     return;

   file = ();
   return (@load) (file);
}

%}}}

% For maximum consistency with previous behavior, 
% turn this off by default:
variable Isis_Use_PHA_Grouping = 0;

define load_data () %{{{
{
   variable msg = "id = load_data (\"pha_filename\", [row_num])";
   variable file, row = 0;

   if (_isis->get_varargs(&file, &row, _NARGS, 1, msg))
     return;

   variable id;

   if (typeof(row) == Integer_Type)
     {
	id = _isis->_load_data (file, row);
        if (Isis_Use_PHA_Grouping)
          {
             variable fun = "use_file_group";
             variable ref = __get_reference (fun);
             if (NULL != ref)(@ref)(id, file);
          }
     }
   else
     {
        variable n = length(row);
        id = Integer_Type[n];
        variable k;

        _for (0, n-1, 1)
          {
             k = ();
             id[k] = _isis->_load_data (file, row[k]);
          }
     }
   
   variable hook = "load_data_hook";
   if (2 == is_defined(hook))
     {
        % is there a cleaner way to do this?
        file, id;  
        eval (hook);
     }   

   return id;
}

%}}}

define load_arf () %{{{
{
   variable msg = "id = load_arf (\"filename\")";
   return load_data_file (&_isis->_load_arf, _NARGS, msg);
}

%}}}

define get_data_exposure () %{{{
{
   variable msg = "time = get_data_exposure (hist_index)";
   variable index;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   index = ();
   return _isis->_get_histogram_exposure_time (index);
}

%}}}

define get_back_exposure () %{{{
{
   variable msg = "time = get_back_exposure (hist_index)";
   variable index;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   index = ();
   return _isis->_get_back_exposure (index);
}

%}}}

define get_arf_exposure () %{{{
{
   variable msg = "time = get_arf_exposure (arf_index)";
   variable index;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   index = ();
   return _isis->_get_arf_exposure_time (index);
}

%}}}

define get_data_backscale () %{{{
{
   variable msg = "area = get_data_backscale (hist_index)";
   variable index;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   index = ();
   
   variable b = _isis->_get_data_region_area (index);
   if (length(b) == 1)
     return b[0];

   return b;
}

%}}}

define set_data_backscale () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_data_backscale (hist_index, area)";
   variable index, area;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (index, area) = ();
   _isis->_set_data_region_area (area, index);
}

%}}}

define get_back_backscale () %{{{
{
   variable msg = "area = get_back_backscale (hist_index)";
   variable index;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   index = ();
   variable b = _isis->_get_back_region_area (index);
   if (length(b) == 1)
     return b[0];

   return b;
}

%}}}

define set_back_backscale () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_back_backscale (hist_index, area)";
   variable index, area;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (index, area) = ();
   _isis->_set_back_region_area (area, index);
}

%}}}

define have_data() %{{{
{
   variable msg = "flag = have_data(index);";
   _isis->have_data();
}

%}}}

define all_data () %{{{
{
   variable msg = "list = all_data ([noticed])";
   variable noticed = 0;

   if (_isis->get_varargs (&noticed, _NARGS, 0, msg))
     return;

   return _isis->_all_data(noticed);
}

%}}}

define all_arfs () %{{{
{
   variable msg = "list = all_arfs ()";

   if (_isis->chk_num_args (_NARGS, 0, msg))
     return;

   return _isis->_all_arfs ();
}

%}}}

define all_rmfs () %{{{
{
   variable msg = "list = all_rmfs ()";

   if (_isis->chk_num_args (_NARGS, 0, msg))
     return;

   return _isis->_all_rmfs ();
}

%}}}

typedef struct
{
   spec_num, order, part, srcid, exclude, combo_id, combo_weight, target,
     instrument, grating, 
     tstart, frame_time, 
     arfs, rmfs,
     notice, notice_list, rebin,
     file, bgd_file
}
Isis_Data_Info_Type;

private define _copy_struct_fields (to, from) %{{{
{
   variable name, value;

   foreach (get_struct_field_names (from))
     {
	name = ();
	value = get_struct_field (from, name);
	set_struct_field (to, name, value);
     }
}

%}}}

private define _get_data_info_struct (id, data_type) %{{{
{
   if (have_data(id) == 0)
     return NULL;

   variable s = @Isis_Data_Info_Type;
   variable info, status;

   (info, s.target, s.instrument, s.grating, s.arfs, s.rmfs, status) 
     = _isis->_get_hist_info (id);
   if (status < 0)
     return NULL;
   
   s.frame_time = _isis->_get_hist_frame_time (id);

   _copy_struct_fields (s, info);
   (s.notice, s.notice_list, s.rebin) = _isis->_get_hist_notice_info (id, data_type);

   return s;
}

%}}}

define get_data_info () %{{{
{
   variable msg = "Struct_Type[] = get_data_info (index[])";
   variable list, data_type, n;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   list = ();

   n = length(list);
   if (n == 0)
     return NULL;

   data_type = _V.data_counts;

   if (n == 1)
     return _get_data_info_struct (list, data_type);

   variable k, s = Isis_Data_Info_Type[n];
   _for (0, n-1, 1)
     {
	k = ();
	s[k] = _get_data_info_struct (list[k], data_type);
     }

   return s;
}

%}}}

define get_kernel () %{{{
{
   variable msg = "get_kernel (index)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   return _isis->_get_kernel ();
}

%}}}

define set_data_info () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_data_info (index, Struct_Type)";
   variable index, s, info;

   if (_isis->get_varargs (&index, &s, _NARGS, 2, msg))
     return;

   info = get_data_info (index);
   if (info == NULL)
     return;

   _copy_struct_fields (info, s);
   _isis->_set_hist_info (info, index);
   
   _isis->_set_hist_frame_time (index, s.frame_time);

   if (NULL != get_struct_field (s, "target"))
     _isis->_set_hist_object_name (index, s.target);
}

%}}}

% data/arf/rmf list maintenance

define assign_arf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "assign_arf (arf_index[], hist_index[])";
   variable arf_index, hist_index, i;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (arf_index, hist_index) = ();

   if (length(arf_index) == 1)
     arf_index = ones(length(arf_index))*arf_index;
   
   array_map (Void_Type, &_isis->_assign_arf_to_hist, arf_index, hist_index);
}

%}}}

define unassign_arf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "unassign_arf (hist_index[])";
   variable idx;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   idx = ();

   assign_arf (0, idx);
}

%}}}

define load_rmf () %{{{
{
   variable msg = "status = load_rmf (\"filename[;args]\")";
   variable arg, method;

   if (_isis->get_varargs (&arg, _NARGS, 1, msg))
     return;

   arg = strtrans(arg, " ", "");
   
   if (is_substr (arg, ".so:"))
     method = _isis->RMF_USER;
   else
     method = _isis->RMF_FILE;

   if (method == _isis->RMF_USER)
     {
	variable n, libname, libpath;
	
	n = is_substr (arg, ":");
	libname = arg[[0:n-2]];
	
	libpath = find_library_name (libname);
	if (libpath == NULL)
	  {
	     vmessage ("File not found:  %s", libname);
	     return -1;
	  }
	
	(arg, ) = strreplace (arg, libname, libpath, 1);
     }
   
   return _isis->_load_rmf (arg, method);
}

%}}}

define rebin_data();
private define do_rmf_assign (rmf, data) %{{{
{
   variable is_grouped = _isis->_is_grouped_data (data);
   variable di = NULL;

   ERROR_BLOCK
     {
	if (is_grouped != 0)
	  {
	     rebin_data (data, di.rebin);
	     _isis->_set_notice_using_mask (di.notice, data);
	  }
     }

   if (is_grouped != 0)
     {
	di = get_data_info (data);
	rebin_data (data, 0);
     }

   _isis->_assign_rmf_to_hist (rmf, data);

   EXECUTE_ERROR_BLOCK;
}

%}}}

define assign_rmf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "assign_rmf (rmf_index[], hist_index[])";
   variable rmf_index, hist_index, i;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (rmf_index, hist_index) = ();
   
   if (length(rmf_index) == 1)
     rmf_index = ones(length(rmf_index))*rmf_index;
   
   array_map (Void_Type, &do_rmf_assign, rmf_index, hist_index);
}

%}}}

define unassign_rmf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "unassign_rmf (hist_index[])";
   variable hist_index, i;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   hist_index = ();

   foreach (hist_index)
     {
	i = ();
	do_rmf_assign (0, i);
     }
}

%}}}

define assign_rsp () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "assign_rsp (arf[], rmf[], hist_index[])";
   variable arfs, rmfs, hist_index, i;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (arfs, rmfs, hist_index) = ();

   if ((length(arfs) == 1) and (length(rmfs) == 1))
     {
	if (have_data(hist_index) != 0)
	  {
	     unassign_rmf(hist_index);
	     unassign_arf(hist_index);
	  }

	assign_rmf (rmfs[0], hist_index);
	assign_arf (arfs[0], hist_index);

	return;
     }
   
   foreach (hist_index)
     {
	i = ();
        
        variable is_grouped = _isis->_is_grouped_data (i);
        variable di = NULL;

        try
          {
             if (is_grouped)
               {
                  di = get_data_info (i);
                  rebin_data (i, 0);
               }
             _isis->_assign_rsp_list (arfs, rmfs, i);
          }
        finally
          {
             if (is_grouped)
               {
                  rebin_data (i, di.rebin);
                  _isis->_set_notice_using_mask (di.notice, i);
               }
          }
     }
}

%}}}

private define _delete_list_members (delete_fun, num, msg) %{{{
{
   variable a = _isis->pop_list (num, msg);
   if (a == NULL)
     return;

   foreach (a) 
     {
	(@delete_fun) ();
     }   
}

%}}}

define delete_data () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "delete_data (hist_index[])";
   return _delete_list_members (&_isis->_delete_hist, _NARGS, msg);
}

%}}}

define delete_arf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "delete_arf (index[])";
   return _delete_list_members (&_isis->_delete_arf, _NARGS, msg);
}

%}}}

define delete_rmf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "delete_rmf (index[])";
   return _delete_list_members (&_isis->_delete_rmf, _NARGS, msg);
}

%}}}

define load_dataset () %{{{
{
   variable msg = "load_dataset (data-file [, rmf-file [, arf-file]])";
   variable p, a, r;

   a = NULL;
   r = NULL;
   
   if (_isis->get_varargs (&p, &r, &a, _NARGS, 1, msg))
     return;

   p = load_data (p);
   if (p == -1)
     verror ("load_dataset: failed to load data file");

   if (r != NULL)
     {
        r = load_rmf (r);
        if (r == -1)
          verror ("load_dataset: failed to load rmf file");
        assign_rmf (r, p);
     }

   if (a != NULL)
     {
        a = load_arf (a);
        if (a == -1)
          verror ("load_dataset: failed to load arf file");
        assign_arf (a, p);
     }
}

%}}}

% rebin/rebin

define set_rebin_error_hook () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_rebin_error_hook (hist_index, &hook)";
   variable id, hook;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (id, hook) = ();
   
   if (typeof(hook) == String_Type)
     hook = __get_reference (hook);
   
   _isis->_set_stat_error_hook (hook, id);
}

%}}}

define rebin_data () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "rebin_data (hist_index[], min_counts | rebin_index[])";
   variable hist_index, arg2, i, fun;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (hist_index, arg2) = ();

   if (typeof(arg2) == Array_Type and length(arg2) > 1)
     fun = &_isis->_rebin_index;
   else
     fun = &_isis->_rebin_min_counts;

   foreach (hist_index)
     {
	i = ();
	(@fun) (i, arg2);
     }
}

%}}}

% This function rebins a data set by an integer factor
% e.g.  group_data (3, 4);
% means to group data set 3 by a factor of 4.

define get_data_counts();
define group_data () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "group_data (hist_index[], factor)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable ds, num;
   (ds, num) = ();

   rebin_data (ds, 0);     % unbin data
   if (num == 0 or num == 1)
     return;

   variable d, len, idx, i;

   foreach (ds)
     {
	i = ();
	d = get_data_counts (i);
	len = length (d.bin_lo);
	idx = ([0:len-1:1]/num mod 2) * 2 - 1;
	rebin_data (i, idx);
     }
}

%}}}

define rebin_array () %{{{
{
   variable msg = "result[] = rebin_array (array[], rebin_flags[])";
   variable a, rebin;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (a, rebin) = ();

   return _isis->_rebin_array (a, rebin);
}

%}}}

private define is_decreasing (x) %{{{
{
   if (length(x) < 2)
     return 0;

   % checking all bins would be slower.
   return (x[0] > x[1]);
}

%}}}

private define order_grid (l, h) %{{{
{
   variable o = is_decreasing(l);

   if (o != is_decreasing(h))
     return (NULL, NULL, -1);

   % checking all bins would be slower.
   if (h[0] < l[0])
     (l, h) = (h, l);

   if (o != 0)
     return (reverse(l), reverse(h), 1);

   return (l, h, 0);
}

%}}}

private define _rebin (l, h, ol, oh, oy) %{{{
{
   variable new_is_reversed, old_is_reversed;

   (l, h, new_is_reversed) = order_grid (l, h);
   if (new_is_reversed == -1)
     return NULL;

   (ol, oh, old_is_reversed) = order_grid (ol, oh);
   if (old_is_reversed == 1)
     oy = reverse(oy);
   else if (old_is_reversed == -1)
     return NULL;

   variable y = _isis->_rebin_histogram (l, h, ol, oh, oy);

   if (new_is_reversed == 1 and y != NULL)
     return reverse(y);

   return y;
}

%}}}

define rebin () %{{{
{
   variable msg = "new_value[] = rebin (new_lo[], new_hi[], lo[], hi[], value[]);";
   variable inlo, inhi, inval, outlo, outhi;

   if (_isis->chk_num_args (_NARGS, 5, msg))
     return;

   (outlo, outhi, inlo, inhi, inval) = ();

   return _rebin (outlo, outhi, inlo, inhi, inval);
}

%}}}

define interpol () %{{{
{
   variable msg = "new_y[] = interpol (new_x[], old_x[], old_y[]);";
   variable newx, oldx, oldy;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (newx, oldx, oldy) = ();
   
   variable newy = _isis->_array_interp (newx, oldx, oldy);
   
   if (length(newy) == 1)
     newy = newy[0];
   
   return newy;   
}

%}}}

define interpol_points () %{{{
{
   variable msg = "new_y[] = interpol_points (new_x[], old_x[], old_y[]);";
   variable newx, oldx, oldy;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (newx, oldx, oldy) = ();

   variable newy = _isis->_array_interp_points (newx, oldx, oldy);
   
   if (length(newy) == 1)
     newy = newy[0];
   
   return newy;
}

%}}}

define rebin_rmf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "rebin_rmf (index, new_lo, new_hi)";
   variable i, lo, hi;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (i, lo, hi) = ();

   _isis->_rebin_rmf (i, lo, hi);
}

%}}}

define rebin_dataset () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "rebin_dataset (hist_index, lo, hi)";
   variable id, lo, hi;

   if (_isis->chk_num_args (_NARGS, 3, msg))
     return;

   (id, lo, hi) = ();

   _isis->_rebin_dataset (lo, hi, id);
}

%}}}

define factor_rsp () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "arf_id = factor_rsp (rmf_id[])";
   variable arf_id, rmf_id;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   rmf_id = ();
   
   arf_id = array_map (Int_Type, &_isis->_factor_rsp, rmf_id);
   
   if (length(arf_id) == 1)
     return arf_id[0];
   else
     return arf_id;
}

%}}}

% get/ put

private define _region_sum (num, msg, version) %{{{
{
   variable hist_index, xmin, xmax, ymin, ymax;

   switch (num)
     {
      case 3:
	(hist_index, xmin, xmax) = ();
	ymin = 0.0;
	ymax = 0.0;
     }
     {
      case 5:
	(hist_index, xmin, xmax, ymin, ymax) = ();
     }
     {
	%default:
	_pop_n (num);
	usage (msg);
	return;
     }

   return _isis->_get_region_sum (hist_index, version, xmin, xmax, ymin, ymax);
}

%}}}

define region_sum () %{{{
{
   variable msg =
     " Struct_Type = region_sum (hist_index, rmin, rmax [, ymin, ymax])";

   return _region_sum (_NARGS, msg, _V.data_counts);
}

%}}}

define region_counts () %{{{
{
   variable msg =
     " Struct_Type = region_counts (hist_index, rmin, rmax [, ymin, ymax])";

   return _region_sum (_NARGS, msg, _V.data_counts);
}

%}}}

define region_flux () %{{{
{
   variable msg =
     " Struct_Type = region_flux (hist_index, rmin, rmax [, ymin, ymax])";

   return _region_sum (_NARGS, msg, _V.data_flux);
}

%}}}

private define _cursor_stats (num, msg, version) %{{{
{
   variable hist_index, file, subtract;

   file = NULL;
   subtract = 1;     % continuum subtraction on by default

   if (_isis->get_varargs (&hist_index, &file, &subtract, num, 1, msg))
     return;

   if (NULL == file)
     file = "region_stats_" + string(hist_index);

   _isis->_cursor_region_stats (hist_index, version, subtract, file);
}

%}}}

define cursor_sum () %{{{
{
   variable msg = "cursor_sum (hist_index [, out-file [, subtract_flag]])";
   _cursor_stats (_NARGS, msg, _V.data_counts);
}

%}}}

define cursor_counts () %{{{
{
   variable msg = "cursor_counts (hist_index [, out-file[, subtract_flag]])";
   _cursor_stats (_NARGS, msg, _V.data_counts);
}

%}}}

define cursor_flux () %{{{
{
   variable msg = "cursor_flux (hist_index [, out-file[, subtract_flag]])";
   _cursor_stats (_NARGS, msg, _V.data_flux);
}

%}}}

private define __get_hist (num, version, msg) %{{{
{
   variable s = struct
     {
	bin_lo, bin_hi, value, err
     };

   if (_isis->chk_num_args (num, 1, msg))
     return;

   variable hist_index = ();

   variable depth = _stkdepth ();
   _isis->_get_hist (hist_index, version);
   depth = _stkdepth () - depth;

   switch (depth)
     {
      case 3:
	(s.bin_lo, s.bin_hi, s.value) = ();
	if (s.bin_lo == NULL) return NULL;
	return s;
     }
     {
      case 4:
	(s.bin_lo, s.bin_hi, s.value, s.err) = ();
	if (s.bin_lo == NULL) return NULL;
	return s;
     }
     {
	% default
	_pop_n (depth);
	vmessage ("Couldn't get data set %d\n", hist_index);
	return NULL;
     }
}

%}}}

define get_data () %{{{
{
   variable msg = "Struct_Type = get_data (hist_index)";

   __get_hist (_NARGS, _V.data_counts, msg);
}

%}}}

define get_data_counts () %{{{
{
   variable msg = "Struct_Type = get_data_counts (hist_index)";

   __get_hist (_NARGS, _V.data_counts, msg);
}

%}}}

define get_data_flux () %{{{
{
   variable msg = "Struct_Type = get_data_flux (hist_index)";

   __get_hist (_NARGS, _V.data_flux, msg);
}

%}}}

private define get_args_for_put_hist (nargs, msg) %{{{
{
   variable hist_index;

   switch (nargs)
     {
      case 2:
	variable s;
	(hist_index, s) = ();
	return (hist_index, s.bin_lo, s.bin_hi, s.value, s.err);
     }
     {
      case 5:
	variable bin_lo, bin_hi, value, err;
	(hist_index, bin_lo, bin_hi, value, err) = ();
	return (hist_index, bin_lo, bin_hi, value, err);
     }
     {
	_pop_n (nargs);
	usage(msg);
     }
}

%}}}

define put_data () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "put_data (hist_index, struct | binlo, binhi, value, uncert)";
   variable hist_index, binlo, binhi, value, uncert;

   (hist_index, binlo, binhi, value, uncert) = get_args_for_put_hist (_NARGS, msg);
   _isis->_put_hist (hist_index, _V.data_counts, binlo, binhi, value, uncert);
}

%}}}

define put_data_counts () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "put_data_counts (hist_index, struct | binlo, binhi, counts, uncert)";
   variable hist_index, binlo, binhi, counts, uncert;

   (hist_index, binlo, binhi, counts, uncert) = get_args_for_put_hist (_NARGS, msg);
   _isis->_put_hist (hist_index, _V.data_counts, binlo, binhi, counts, uncert);
}

%}}}

define put_data_flux () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "put_data_flux (hist_index, struct | binlo, binhi, flux, uncert)";
   variable hist_index, binlo, binhi, flux, uncert;

   (hist_index, binlo, binhi, flux, uncert) = get_args_for_put_hist (_NARGS, msg);
   _isis->_put_hist (hist_index, _V.data_flux, binlo, binhi, flux, uncert);
}

%}}}

private define get_args_for_define_hist (num, msg) %{{{
{
   variable s = struct
     {
	bin_lo, bin_hi, value, err
     };

   variable v;

   switch (num)
     {
      case 1:
	v = ();
	if (typeof(v) == Struct_Type)
	  s = v;
	else
	  {
	     s.value = v;
	     s.bin_lo = 0;
	     s.bin_hi = 0;
	     s.err = 0;
	  }
     }
     {
      case 4:
	(s.bin_lo, s.bin_hi, s.value, s.err) = ();
     }
     {
	% default:
	_pop_n (num);
	usage (msg);
	return NULL;
     }

   return s;
}

%}}}

define define_counts () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = define_counts (struct | bins | lo, hi, n, err);";
   variable s;

   s = get_args_for_define_hist (_NARGS, msg);
   if (s == NULL)
     return;

   return _isis->_define_hist (_V.data_counts, s.bin_lo, s.bin_hi, s.value, s.err);
}

%}}}

define define_flux () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "s = define_flux (struct | lo, hi, f, err);";
   variable s;

   s = get_args_for_define_hist (_NARGS, msg);
   if (s == NULL) return;

   return _isis->_define_hist (_V.data_flux, s.bin_lo, s.bin_hi, s.value, s.err);
}

%}}}

define back_fun () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "back_fun (hist_index, fun_string)";
   variable index, s;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (index, s) = ();
   if (s == NULL) s = "";
   _isis->_set_instrument_background (index, s);
}
define iback_fun ()
{
   vmessage ("*** %s will soon be renamed to back_fun", _function_name);
   variable args = __pop_args (_NARGS);
   back_fun (__push_args(args));
}
%}}}

define _define_back () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "status = _define_back (index, bgd [, area [, exposure]]])";
   variable i, bgd, area, t;

   area = -1.0;
   t = -2.0;

   if (_isis->get_varargs (&i, &bgd, &area, &t, _NARGS, 2, msg))
     return;

   return _isis->_define_bgd (area, bgd, i, t);
}
define _define_bgd ()
{
   vmessage ("*** %s will soon be renamed to _define_back", _function_name);
   variable args = __pop_args (_NARGS);
   _define_back (__push_args(args));
}
%}}}

private define do_define_bgd (i, file) %{{{
{
   if (file != NULL)
     return _isis->_define_bgd_file (i, file);
   else
     return _isis->_define_bgd (NULL, i, -1.0);
}

%}}}

define define_back () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "status = define_back (index, file)";
   variable n, i, indexes, files;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (indexes, files) = ();

   n = length (indexes);

   switch (length(files))
     {
      case 1:
	foreach (indexes)
	  {
	     i = ();
	     if (-1 == do_define_bgd (i, files))
	       return -1;
	  }
     }
     {
      case n:
	_for (0, n-1, 1)
	  {
	     i = ();
	     if (-1 == do_define_bgd (indexes[i], files[i]))
	       return -1;
	  }
     }
     {
	% default
	usage(msg);
     }

   return 0;
}
define define_bgd ()
{
   vmessage ("*** %s will soon be renamed to define_back", _function_name);
   variable args = __pop_args (_NARGS);
   define_back (__push_args(args));
}
%}}}

define get_model () %{{{
{
   variable msg = "Struct_Type = get_model (hist_index)";

   __get_hist (_NARGS, _V.model_counts, msg);
}

%}}}

define get_model_counts () %{{{
{
   variable msg = "Struct_Type = get_model_counts (hist_index)";

   __get_hist (_NARGS, 0, msg);
}
%}}}

define get_model_flux () %{{{
{
   variable msg = "Struct_Type = get_model_flux (hist_index)";

   __get_hist (_NARGS, _V.model_flux, msg);
}

%}}}

define get_convolved_model_flux () %{{{
{
   variable msg = "Struct_Type = get_convolved_model_flux (hist_index)";

   __get_hist (_NARGS, _V.convolved_model_flux, msg);
}

%}}}

define put_model_counts () %{{{
{
   variable msg = "put_model_counts (hist_index, counts[])";
   variable index, value;
   if (_isis->chk_num_args(_NARGS, 2, msg))
     return;
   (index, value) = ();
   _isis->put_model_intrin (index, _V.model_counts, value);
}
%}}}

define put_model_flux () %{{{
{
   variable msg = "put_model_flux (hist_index, flux[])";
   variable index, value;
   if (_isis->chk_num_args(_NARGS, 2, msg))
     return;
   (index, value) = ();
   _isis->put_model_intrin (index, _V.model_flux, value);
}
%}}}

define put_convolved_model_flux () %{{{
{
   variable msg = "put_convolved_model_flux (hist_index, value[])";
   variable index, value;
   if (_isis->chk_num_args(_NARGS, 2, msg))
     return;
   (index, value) = ();
   _isis->put_model_intrin (index, _V.convolved_model_flux, value);
}
%}}}

define get_rmf_data_grid () %{{{
{
   variable msg = "Struct_Type = get_rmf_data_grid (index)";   
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;
   variable g = _isis->_get_rmf_data_grid ();
   if (g.bin_lo == NULL)
     g = NULL;
   return g;
}

%}}}

define get_rmf_arf_grid () %{{{
{
   variable msg = "Struct_Type = get_rmf_arf_grid (index)";   
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;
   variable g = _isis->_get_rmf_arf_grid ();
   if (g.bin_lo == NULL)
     g = NULL;
   return g;
}

%}}}

define find_rmf_peaks() %{{{
{
   variable msg = "h_P = find_rmf_peaks (index)";
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;
   return _isis->_find_rmf_peaks ();
}

%}}}

define get_arf () %{{{
{
   variable msg = "Struct_Type = get_arf (index)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable index = ();
   variable a = struct
     {
	bin_lo, bin_hi, value, err
     };

   (a.bin_lo, a.bin_hi, a.value, a.err) = _isis->_get_arf (index);

   if (a.bin_lo == NULL)
     return NULL;

   return a;
}

%}}}

define put_arf () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "put_arf (arf_index, struct | bin_lo, bin_hi, value, err)";
   variable arf_index, bin_lo, bin_hi, value, err;

   switch (_NARGS)
     {
      case 2:
	variable s;
	(arf_index, s) = ();
	bin_lo = s.bin_lo;
	bin_hi = s.bin_hi;
	value = s.value;
	err = s.err;
     }
     {
      case 5:
	(arf_index, bin_lo, bin_hi, value, err) = ();
     }
     {
	_pop_n (_NARGS);
	usage(msg);
	return;
     }

   _isis->_put_arf (arf_index, bin_lo, bin_hi, value, err);
}

%}}}

define define_arf () %{{{
{
   variable msg = "id = define_arf (Struct_Type | bin_lo, bin_hi, value, err)";
   variable s = struct {bin_lo, bin_hi, value, err};
   variable v;

   switch (_NARGS)
     {
      case 1:
	v = ();
	if (typeof(v) == Struct_Type)
	  s = v;
	else
	  {
	     s.value = v;
	     s.bin_lo = 0;
	     s.bin_hi = 0;
	     s.err = 0;
	  }
     }
     {
      case 4:
	(s.bin_lo, s.bin_hi, s.value, s.err) = ();
     }
     {
	% default:
	usage (msg);
	_pop_n (_NARGS);
	return -1;
     }

   return _isis->_define_arf (s.bin_lo, s.bin_hi, s.value, s.err);
}

%}}}

define get_arf_info () %{{{
{
   variable msg = "Struct_Type = get_arf_info (arf_index)";
   variable id;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   id = ();

   variable t, fracexpo, status;
   (t, fracexpo, status) = _isis->_get_arf_info (id);

   if (length(fracexpo) == 1)
     fracexpo = fracexpo[0];

   if (status < 0)
     return NULL;

   variable s = struct_combine (t, "fracexpo");
   s.fracexpo = fracexpo;

   return s;
}

%}}}

define set_arf_info () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_arf_info (index, Struct_Type)";
   variable index, s;

   if (_isis->get_varargs (&index, &s, _NARGS, 2, msg))
     return;

   if (typeof (s) != Struct_Type)
     {
	usage(msg);
	return;
     }

   _isis->_set_arf_info (s.fracexpo, s, index);
}

%}}}

define get_rmf_info () %{{{
{
   variable msg = "Struct_Type = get_rmf_info (rmf_index)";
   variable id;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   id = ();

   variable t, status;
   (t, status) = _isis->_get_rmf_info (id);

   if (status < 0)
     return NULL;

   return t;
}

%}}}

define set_rmf_info () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_rmf_info (index, Struct_Type)";
   variable index, s;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (index, s) = ();

   if (typeof (s) != Struct_Type)
     {
	usage(msg);
	return;
     }

   _isis->_set_rmf_info (s, index);
}

%}}}

define get_flux_corr_weights () %{{{
{
   variable msg = "wt[] = get_flux_corr_weights (hist_index)";   
   variable id;
   
   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   id = ();
   return _isis->get_flux_corr_weights_intrin (id);
}

%}}}

private define do_flux_corr (nargs, msg, flux_corr_ref) %{{{
{
   variable hist_index = 1;
   variable threshold = 0.0;
   variable options = NULL;

   if (_isis->get_varargs (&hist_index, &threshold, &options, nargs, 1, msg))
     return;

   % last arg is 0 because we want the ARF & RMF that
   % would be used for fitting counts data (ie for forward-folding)
   _isis->_set_fit_type (Assigned_ARFRMF, 0);

   % Instantiate the kernel (if necessary)
   % and sync model config with current dataset config.
   % For efficiency, do this outside the loop over datasets
   if (-1 == _isis->_sync_model_with_data ())
     return;

   variable k;
   foreach (hist_index)
     {
	k = ();
	(@flux_corr_ref) (options, k, threshold);
     }
}

%}}}

define flux_corr () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "flux_corr (hist[] [, threshold[, options]])";
   do_flux_corr (_NARGS, msg, &_isis->_flux_correct);
}

%}}}

define flux_corr_model_counts () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "flux_corr_model_counts (hist[] [, threshold[, options]])";
   do_flux_corr (_NARGS, msg, &_isis->_flux_correct_model_counts);
}

%}}}

define set_data_exposure () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_data_exposure (hist_index, exposure_sec)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable hist_index, exposure_sec;
   (hist_index, exposure_sec) = ();
   _isis->_set_histogram_exposure_time (hist_index, exposure_sec);
}

%}}}

define set_arf_exposure () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_arf_exposure (arf_index, exposure_sec)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable arf_index, exposure_sec;
   (arf_index, exposure_sec) = ();

   variable k;
   foreach (arf_index)
     {
	k = ();
	_isis->_set_arf_exposure_time (k, exposure_sec);
     }
}

%}}}

define set_frame_time () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_frame_time (hist_index, frame_time_sec)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable hist_index, frame_time_sec;
   (hist_index, frame_time_sec) = ();
   _isis->_set_hist_frame_time (hist_index, frame_time_sec);
}

%}}}

define get_frame_time () %{{{
{
   variable msg = "frame_time_sec = get_frame_time (hist_index)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable hist_index = ();
   return _isis->_get_hist_frame_time (hist_index);
}

define copy_data_keywords () %{{{
{
   variable msg = "copy_data_keywords (dest_index, src_index)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable d, s;
   (d, s) = ();
   _isis->_copy_hist_keywords (d, s);
}

%}}}

%}}}

define set_dataset_metadata () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "set_dataset_metadata (hist_index, Any_Type)";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable hist_index, any;
   (hist_index, any) = ();
   _isis->_set_dataset_metadata (any, hist_index);
}

%}}}

define get_dataset_metadata () %{{{
{
   variable msg = "Any_Type = get_dataset_metadata (hist_index)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable hist_index = ();
   _isis->_get_dataset_metadata (hist_index);
}

%}}}

define set_sys_err_frac () %{{{
{
   variable msg = "set_sys_err_frac (hist_index, sys_err_frac[])";

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   variable hist_index, sef;
   (hist_index, sef) = ();

   _isis->_set_sys_err_frac_intrin (sef, hist_index);
}

%}}}

define get_sys_err_frac () %{{{
{
   variable msg = "sys_err_frac = set_sys_err_frac (hist_index)";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable id = ();
   _isis->_get_sys_err_frac_intrin (id);
}

%}}}

% plot style

define errorbars () %{{{
{
   variable msg = "errorbars (on_off [, terminal_length])";
   variable flag, length;

   length = -1;  % means keep the internal default

   if (_isis->get_varargs (&flag, &length, _NARGS, 1, msg))
     return;

   _isis->_set_errorbar (flag, length);
}

%}}}

define set_data_color () %{{{
{
   variable msg = "set_data_color (hist-index-list, color)";
   variable index, color, i;

   if (_isis->chk_num_args (_NARGS, 2, msg))
     return;

   (index, color) = ();

   foreach (index)
     {
	i = ();
	_isis->_set_hist_color (i, color);
     }
}

%}}}

define unset_data_color () %{{{
{
   variable msg = "unset_data_color (hist-index)";
   variable index, i, n;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   index = ();

   foreach (index)
     {
	i = ();
	_isis->_unset_hist_color (i);
     }
}

%}}}

% plotting

private define get_args_hist_style (num) %{{{
{
   variable hist_index = NULL;
   variable style = NULL;

   switch (num)
     {
      case 2:
	(hist_index, style) = ();
     }
     {
      case 1:
	hist_index = ();
	style = -1;
     }
     {
      case 0:
	if (all_data != NULL)
	  {
	     if (length(all_data) == 1)
	       hist_index = all_data[0];
	     style = -1;
	  }
     }
     {
	% default:
	_pop_n (num);
     }

   return (hist_index, style);
}

%}}}

private define __plot_hist(nargs, msg, version, overlay, residuals) %{{{
{
   variable hist_index, style;

   (hist_index, style) = get_args_hist_style (nargs);
   if (hist_index == NULL)
     {
	usage (msg);
	return;
     }

   _isis->_plot_hist (hist_index, version, style, overlay, residuals);
}

%}}}

define _rplot_counts () %{{{
{
   variable msg = "_rplot_counts (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.model_counts, 0, 1);    % no overlay
}

%}}}

define _orplot_counts () %{{{
{
   variable msg = "_orplot_counts (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.model_counts, 1, 1);    % overlay
}

%}}}

define _rplot_flux () %{{{
{
   variable msg = "_rplot_flux (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.convolved_model_flux, 0, 1);    % no overlay
}

%}}}

define _orplot_flux () %{{{
{
   variable msg = "_orplot_flux (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.convolved_model_flux, 1, 1);    % overlay
}

%}}}

define plot_data () %{{{
{
   variable msg = "plot_data (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.data_counts, 0, 0);    % no overlay
}

%}}}

define oplot_data () %{{{
{
   variable msg = "oplot_data (hist-index [, style])";
   __plot_hist(_NARGS, msg, _V.data_counts, 1, 0);     % overlay
}

%}}}

define plot_data_counts () %{{{
{
   variable msg = "plot_data_counts (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.data_counts, 0, 0);    % no overlay
}

%}}}

define oplot_data_counts () %{{{
{
   variable msg = "oplot_data_counts (hist-index [, style])";
   __plot_hist(_NARGS, msg, _V.data_counts, 1, 0);     % overlay
}

%}}}

define plot_data_flux () %{{{
{
   variable msg = "plot_data_flux (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.data_flux, 0, 0);    % no overlay
}

%}}}

define oplot_data_flux () %{{{
{
   variable msg = "oplot_data_flux (hist-index [, style])";
   __plot_hist(_NARGS, msg, _V.data_flux, 1, 0);     % overlay
}

%}}}

define plot_model() %{{{
{
   variable msg = "plot_model (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.model_counts, 0, 0);  % no overlay
}

%}}}

define oplot_model() %{{{
{
   variable msg = "oplot_model(hist-index [, style])";
   __plot_hist(_NARGS, msg, _V.model_counts, 1, 0);  % overlay
}

%}}}

define plot_model_counts () %{{{
{
   variable msg = "plot_model_counts (hist-index [, style])";
   __plot_hist (_NARGS, msg, 0, 0, 0);    % no overlay
}

%}}}

define oplot_model_counts () %{{{
{
   variable msg = "oplot_model_counts (hist-index [, style])";
   __plot_hist (_NARGS, msg, 0, 1, 0);    % overlay
}

%}}}

define plot_model_flux () %{{{
{
   variable msg = "plot_model_flux (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.model_flux, 0, 0);    % no overlay
}
%}}}

define oplot_model_flux () %{{{
{
   variable msg = "oplot_model_flux (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.model_flux, 1, 0);    % overlay
}

%}}}

define plot_convolved_model_flux () %{{{
{
   variable msg = "plot_convolved_model_flux (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.convolved_model_flux, 0, 0);    % no overlay
}
%}}}

define oplot_convolved_model_flux () %{{{
{
   variable msg = "oplot_convolved_model_flux (hist-index [, style])";
   __plot_hist (_NARGS, msg, _V.convolved_model_flux, 1, 0);    % overlay
}

%}}}

private define __get_index_to_plot (num, msg) %{{{
{
   variable id = -1;

   switch (num)
     {
      case 1:
	id = ();
     }
     {
      case 0:
	if (length(all_data) == 1)
	  id = all_data[0];
	else
	  usage (msg);
     }
     {
	%default:
	usage(msg);
     }

   return id;
}

%}}}

define rplot_counts () %{{{
{
   variable msg = "rplot_counts (hist_index)";

   if (_NARGS == 0 and (all_data == NULL or length(all_data) > 1))
     {
	usage (msg);
	return;
     }

   variable id = __get_index_to_plot (_NARGS, msg);
   if (id == -1)
     return;
   
   ERROR_BLOCK
     {
	multiplot (1);
     }

   multiplot ([3,1]);

   variable ebar_state = _isis->_errorbar_state();
   variable new_ebar_state;
   if (ebar_state > 0)
     new_ebar_state = ebar_state;
   else
     new_ebar_state = 1;
   errorbars (new_ebar_state);

   variable info = get_plot_info();
   plot_data_counts(id);
   set_line_width (2*info.line_width);
   oplot_model_counts(id);
   set_line_width (info.line_width);

%   variable ymin, ymax;
%   (ymin, ymax) = _isis->_get_yrange();
   yrange;
   errorbars(0);
   _rplot_counts (id);
   if (ebar_state >= 0) errorbars (ebar_state);

   if (info.ylog == 0)
     yrange (info.ymin, info.ymax);
   else
     {
        ylin;
        yrange(10.0^info.ymin, 10.0^info.ymax);
        ylog;
     }

   multiplot (1);
}

%}}}

define rplot_flux () %{{{
{
   variable msg = "rplot_flux (hist_index)";

   if (_NARGS == 0 and (all_data == NULL or length(all_data) > 1))
     {
	usage (msg);
	return;
     }

   variable id = __get_index_to_plot (_NARGS, msg);
   if (id == -1)
     return;

   ERROR_BLOCK
     {
	multiplot (1);
     }
   
   multiplot ([3,1]);

   variable ebar_state = _isis->_errorbar_state();
   variable new_ebar_state;
   if (ebar_state > 0)
     new_ebar_state = ebar_state;
   else
     new_ebar_state = 1;
   errorbars (new_ebar_state);

   variable info = get_plot_info();
   plot_data_flux(id);
   set_line_width (2*info.line_width);
   oplot_convolved_model_flux(id);
   set_line_width (info.line_width);

%   variable ymin, ymax;
%   (ymin, ymax) = _isis->_get_yrange();
   yrange;
   errorbars(0);
   _rplot_flux (id);
   if (ebar_state >= 0) errorbars (ebar_state);

   if (info.ylog == 0)
     yrange (info.ymin, info.ymax);
   else
     {
        ylin;
        yrange(10.0^info.ymin, 10.0^info.ymax);
        ylog;
     }

   multiplot (1);
}

%}}}

define eval_counts ();
define fakeit () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "fakeit ([&noise_fun])";
   variable noticed = 1;

   % default to Poisson noise
   variable noise_fun = &prand;

   if (_isis->get_varargs (&noise_fun, _NARGS, 0, msg))
     return;

   variable datasets = all_data (noticed);
   if (datasets == NULL)
     return;

   % Don't overwrite real data !
   variable is_fake = array_map (Integer_Type, &_isis->_is_fake_data, datasets);
   datasets = datasets[where (is_fake)];
   if (length(datasets) == 0)
     return;

   % save the current rebin state
   variable di = array_map (Struct_Type, &get_data_info, datasets);
   array_map (Void_Type, &rebin_data, datasets, 0);

   try
     {
        variable verbose = Fit_Verbose;
        Fit_Verbose = -1;
        variable ret = eval_counts ();
        Fit_Verbose = verbose;
        if (ret < 0)
          error ("*** Failed evaluating model function");

        variable id, m, stat_err;

        foreach id (datasets)
          {
             m = get_model_counts (id);
             if (noise_fun != NULL)
               m.value = array_map (Integer_Type, noise_fun, m.value);
             stat_err = sqrt(abs(m.value));
             put_data_counts (id, m.bin_lo, m.bin_hi, m.value, stat_err);
          }
     }
   finally
     {
        % restore the grouping state
        variable i;
        _for i (0, length(datasets)-1, 1)
          {
             rebin_data (datasets[i], di[i].rebin);
          }
     }
}

%}}}

define set_fake () %{{{
{
   variable msg = "set_fake (id, 0|1)";
   variable id, value;
   
   if (_isis->get_varargs (&id, &value, _NARGS, 2, msg))
     return;
   
   () = _isis->_set_fake (id, value);
}

%}}}

% notice/ ignore

private define get_notice_args (num) %{{{
{
   variable hist_index, lo, hi;
   lo = _isis->DBL_MIN;
   hi = _isis->DBL_MAX;
   hist_index = NULL;

   switch (num)
     {
      case 3:
	(hist_index, lo, hi) = ();
     }
     {
      case 2:
	(lo, hi) = ();
	if (length(all_data) == 1)
	  hist_index = all_data[0];
     }
     {
      case 1:
	hist_index = ();
     }
     {
	% default:
	_pop_n(num);
     }

   if (hist_index != NULL
       and _typeof(hist_index) != Integer_Type
       and _typeof(hist_index) != UInteger_Type)
     {
	error ("*** dataset index should be an integer");
     }

   if (NULL == lo)
     lo = _isis->DBL_MIN;
   if (NULL == hi)
     hi = _isis->DBL_MAX;

   return (lo, hi, hist_index);
}

%}}}

private define apply_notice (val, lo, hi, datasets) %{{{
{
   variable id;
   foreach (datasets)
     {
	id = ();
	_isis->_set_notice (val, lo, hi, id);
     }
}

%}}}

private define set_notice_list (nargs, value, msg) %{{{
{
   variable id, list;

   if (_isis->chk_num_args (nargs, 2, msg))
     return;

   (id, list) = ();

   _isis->_set_notice_using_list (list, id, value);
}

%}}}

define notice_list () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   set_notice_list (_NARGS, 1, "notice_list (id, where)");
}

%}}}

define ignore_list () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   set_notice_list (_NARGS, 0, "ignore_list (id, where)");
}

%}}}

define notice() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "notice (hist_index[] [, lo, hi])";
   variable lo, hi, hist_index;

   (lo, hi, hist_index) = get_notice_args(_NARGS);
   if (NULL == hist_index)
     {
	usage(msg);
	return;
     }

   apply_notice (1, lo, hi, hist_index);
}

%}}}

define notice_en () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "notice_en (hist_index[] [, lo, hi])";
   variable lo, hi, hist_index;

   (lo, hi, hist_index) = get_notice_args(_NARGS);
   if (NULL == hist_index)
     {
	usage(msg);
	return;
     }

   apply_notice (1, _A(lo, hi), hist_index);
}

%}}}

define ignore_bad() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg = "ignore_bad (hist_index[])";
   variable index, i;

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   index = ();

   foreach (index)
     {
	i = ();
	_isis->_ignore_bad (i);
     }
}

%}}}

define ignore() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "ignore (hist_index[] [, lo, hi])";
   variable lo, hi, hist_index;

   (lo, hi, hist_index) = get_notice_args(_NARGS);
   if (NULL == hist_index)
     {
	usage(msg);
	return;
     }

   apply_notice (0, lo, hi, hist_index);
}

%}}}

define ignore_en() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "ignore_en (hist_index[] [, lo, hi])";
   variable lo, hi, hist_index;

   (lo, hi, hist_index) = get_notice_args(_NARGS);
   if (NULL == hist_index)
     {
	usage(msg);
	return;
     }

   apply_notice (0, _A(lo, hi), hist_index);
}

%}}}

define xnotice() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "xnotice (hist_index[] [, lo, hi])";
   variable lo, hi, hist_index;

   (lo, hi, hist_index) = get_notice_args (_NARGS);
   if (NULL == hist_index)
     {
	usage (msg);
	return;
     }

   apply_notice (0, -_isis->DBL_MAX, _isis->DBL_MAX, hist_index);
   apply_notice (1, lo, hi, hist_index);
}

%}}}

define xnotice_en() %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "xnotice_en (hist_index[] [, lo, hi])";
   variable lo, hi, hist_index;

   (lo, hi, hist_index) = get_notice_args (_NARGS);
   if (NULL == hist_index)
     {
	usage (msg);
	return;
     }

   apply_notice (0, -_isis->DBL_MAX, _isis->DBL_MAX, hist_index);
   apply_notice (1, _A(lo, hi), hist_index);
}

%}}}

define exclude () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "exclude (data[])";
   variable x, list = _isis->pop_list (_NARGS, msg);
   if (list == NULL)
     return;

   foreach (list)
     {
	x = ();
	_isis->_set_exclude_flag (x, 1);
     }
}

%}}}

define include () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "include (data[])";
   variable x, list = _isis->pop_list (_NARGS, msg);
   if (list == NULL)
     return;

   foreach (list)
     {
	x = ();
	_isis->_set_exclude_flag (x, 0);
     }
}

%}}}

% dataset combinations

define match_dataset_grids () %{{{
{
   _isis->error_if_fit_in_progress (_function_name);
   variable msg= "match_dataset_grids (indices[])";
   variable list = _isis->pop_list (_NARGS, msg);
   if (list == NULL)
     return;

   variable n = length(list);

   if (n < 2)
     return;

   % Undo any grouping so rebin_dataset wont complain.
   % This also notices all bins.
   rebin_data (list, 0);

   % Match all the datasets to the grid of the first.
   variable d = get_data_counts (list[0]);

   foreach (list[[1:n-1]])
     {
	variable k = ();
	rebin_dataset (k, d.bin_lo, d.bin_hi);
     }
}

%}}}

% Contributed by Mike Nowak <mnowak@space.mit.edu>
define i2x_group () %{{{
{
   variable msg = "isis_group[] = i2x_group (xspec_group[])";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable gdi = ();

   variable i, lb, ub, lgdi, iw, liw, grp;
   gdi = [reverse(gdi),5];
   lgdi = length(gdi)-1;
   grp = Integer_Type[lgdi];

   iw = [where( shift(gdi,-1) != gdi )];
   liw = length(iw);
   iw[liw-1] = lgdi;

   grp[iw[[0:liw-2]]] = 1;

   foreach i ([0:liw-2])
     {
        lb = iw[i];
        ub = iw[i+1]-1;
        if(ub>lb)
          {
             grp[[lb+1:ub]] = -1;
          }
     }

   return grp;
}

%}}}

% Contributed by Mike Nowak <mnowak@space.mit.edu>
define x2i_group () %{{{
{
   variable msg = "xspec_group[] = x2i_group (isis_group[])";

   if (_isis->chk_num_args (_NARGS, 1, msg))
     return;

   variable grp = ();

   variable lgrp, iw, liw;
   lgrp = length(grp);
   iw = where(grp == 1);
   liw = length(iw);
   iw = [iw,lgrp];

   variable i, sn = 1;
   foreach i ([0:liw-1])
     {
        grp[[iw[i]:max([iw[i],iw[i+1]-1])]] = sn;
        sn = -1*sn;
     }

   return reverse(grp);
}

%}}}

