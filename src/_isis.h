#ifndef ISIS_ISIS_H
#define ISIS_ISIS_H

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2011 Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.
    
    Author:  John C. Houck  <houck@space.mit.edu>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/* $Id: _isis.h,v 1.10 2004/09/10 03:10:39 houck Exp $ */

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

extern int init_isis_intrinsics (char *ns_name);
extern void init_signals (void);
extern void quit_isis (int reset);
extern void exit_isis (int err);
extern void error_hook (char *err);

extern void input_from_stdin (void);

extern void interpreter_loop (void);
extern int open_readline (char *);
extern unsigned int isis_getkey (void);

extern void isis_set_pager (char *pager);

#ifdef ISIS_DBATOMIC_H
extern DB_t *ptr_to_atomic_db (void);
extern DB_t *_ptr_to_atomic_db (void);  /* silent */
#endif

#ifdef ISIS_DBEM_H
extern EM_t *ptr_to_emissivity_db (void);
extern EM_t *_ptr_to_emissivity_db (void);  /* silent */
#endif

#ifdef ISIS_PLOT_H
extern Plot_t *current_format (void);
extern void reset_current_pane (Plot_t *fmt);
#endif

extern int current_plot_shows_bin_density (void);
extern int force_open_plot_device (void);
extern int change_xunits_to_angstrom (float *x, float *y);

extern int sync_model_with_data (void);

extern int update_user_model (void);

#ifdef ISIS_HISTOGRAM_H
extern Hist_t *find_hist (int hist_index);
extern int get_kernel_params_for_hist (Hist_t *, double **, unsigned int *);
extern int pop_instrumental_background (SLang_Array_Type **bgd, Hist_t *h);
#endif

extern void set_hook_from_stack (SLang_Name_Type **hook);

extern void add_slangfun_intrin (void);
extern void add_cfun_intrin (void);
extern void del_function (char *fun_name);
extern void set_function_category (char *fun_name, unsigned int *category);
extern void function_list (void);
extern void _add_slang_statistic (char *);
extern void list_statistics_and_engines (void);

extern int init_fit_functions (void);
extern void deinit_fit_functions (void);
extern int init_fit_module_internals (void);

extern Isis_Fit_Statistic_Type *current_fit_statistic (void);
extern Isis_Fit_Engine_Type *current_fit_method (void);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
