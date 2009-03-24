/* -*- mode: C; mode: fold -*- */

/*{{{ Includes  */

#include "config.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdarg.h>
#include <string.h>

#include <setjmp.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include <slang.h>

/*}}}*/

#define ISIS_FIT_ENGINE_PRIVATE_DATA \
   SLang_Name_Type *sl_optimize; \
   SLang_Name_Type *sl_set_options;

#include "isis.h"
#include "util.h"
#include "fit.h"
#include "errors.h"

static SLang_MMT_Type *Current_Fit_Object_MMT;

static void slfe_deallocate (Isis_Fit_Engine_Type *e) /*{{{*/
{
   ISIS_FREE (e->engine_name);
   ISIS_FREE (e->default_statistic_name);
   ISIS_FREE (e->option_string);
   SLang_free_function (e->sl_optimize);
   SLang_free_function (e->sl_set_options);
}

/*}}}*/

static int slfe_optimize (Isis_Fit_Type *ift, void *clientdata, /*{{{*/
                          double *x, double *y, double *weights, unsigned int npts,
                          double *pars, unsigned int npars)
{
   Isis_Fit_Engine_Type *e;
   SLang_Array_Type *sl_pars=NULL, *sl_pars_min=NULL, *sl_pars_max=NULL;
   SLang_Array_Type *sl_new_pars=NULL;
   int n, status = -1;

   (void) clientdata; (void) x; (void) y; (void) weights; (void) npts;

   if ((ift == NULL) || (pars == NULL) || (npars <= 0)
       || (Current_Fit_Object_MMT == NULL))
     return -1;

   e = ift->engine;

   n = (int) npars;
   sl_pars = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1);
   sl_pars_min = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1);
   sl_pars_max = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1);

   if ((NULL == sl_pars) || (NULL == sl_pars_min) || (NULL == sl_pars_max))
     return -1;

   memcpy ((char *)sl_pars->data, (char *)pars, npars * sizeof(double));
   memcpy ((char *)sl_pars_min->data, (char *)e->par_min, npars * sizeof(double));
   memcpy ((char *)sl_pars_max->data, (char *)e->par_max, npars * sizeof(double));

   /* Increment the reference count to prevent a segv
    * if the user deletes their copy.  Is there a better
    * way to handle this? */
   SLang_inc_mmt (Current_Fit_Object_MMT);

   SLang_start_arg_list ();
   if ((-1 == SLang_push_mmt (Current_Fit_Object_MMT))
       || (-1 == SLang_push_array (sl_pars, 1))
       || (-1 == SLang_push_array (sl_pars_min, 1))
       || (-1 == SLang_push_array (sl_pars_max, 1)))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "calling user-defined optimization method '%s'",
                    e->engine_name);
        goto return_error;
     }
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (e->sl_optimize))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "executing optimization method '%s'",
                    e->engine_name);
        goto return_error;
     }

   if (-1 == SLang_pop_array_of_type (&sl_new_pars, SLANG_DOUBLE_TYPE))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "returning results from optimization method '%s'",
                    e->engine_name);
        goto return_error;
     }

   if ((sl_new_pars == NULL) || (sl_new_pars->num_elements != npars))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "corrupted parameter array returned from optimization method '%s'",
                    e->engine_name);
        goto return_error;
     }

   memcpy ((char *)pars, (char *)sl_new_pars->data, npars * sizeof(double));

   status = 0;
return_error:
   SLang_free_array (sl_new_pars);

   if (SLang_get_error())
     {
        isis_throw_exception (SLang_get_error());
        return -1;
     }

   return status;
}

/*}}}*/

static int slfe_set_options (Isis_Fit_Engine_Type *e, Isis_Option_Type *opts) /*{{{*/
{
   SLang_Array_Type *sl_opts;
   int i, n, status;

   if (opts == NULL)
     return -1;

   n = opts->num_options;
   if (n == 0)
     return 0;

   if (NULL == (sl_opts = SLang_create_array (SLANG_STRING_TYPE, 1, NULL, &n, 1)))
     return -1;

   for (i = 0; i < n; i++)
     {
        int have_value = (opts->option_values[i] != 0);
        char *s;

        if (have_value)
          {
             s = isis_mkstrcat (opts->option_names[i], "=",
                                opts->option_values[i], NULL);
          }
        else s = opts->option_names[i];

        if ((s == NULL)
            || (-1 == SLang_set_array_element (sl_opts, &i, &s)))
          {
             SLang_free_array (sl_opts);
             if (have_value) ISIS_FREE(s);
          }

        if (have_value) ISIS_FREE(s);
     }

   SLang_start_arg_list();
   (void) SLang_push_array (sl_opts, 1);
   SLang_end_arg_list();

   /* converts options array to a struct */
   SLang_execute_function ("_isis->options_to_struct");

   /* this function then pops the struct off the stack */
   if (-1 == SLexecute_function (e->sl_set_options))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "setting options for fit method '%s'",
                    e->engine_name);
        return -1;
     }

   if (-1 == SLang_pop_integer (&status))
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__,
                    "status return failed after setting options for fit method '%s'",
                    e->engine_name);
        return -1;
     }

   return status;
}

/*}}}*/

static Isis_Fit_Engine_Type *add_slang_fit_engine (char *eng_name, char *stat_name) /*{{{*/
{
   Isis_Fit_Engine_Type *e;

   if (NULL == (e = ISIS_MALLOC (sizeof(Isis_Fit_Engine_Type))))
     return NULL;
   memset ((char *)e, 0, sizeof (*e));

   if ((NULL == (e->engine_name = isis_make_string (eng_name)))
       || (NULL == (e->default_statistic_name = isis_make_string (stat_name))))
     {
        slfe_deallocate (e);
        ISIS_FREE (e);
        return NULL;
     }

   e->method = &slfe_optimize;
   e->deallocate = &slfe_deallocate;
   e->set_options = &slfe_set_options;
   e->set_range_hook = NULL;
   e->range_hook = NULL;
   e->verbose_hook = NULL;
   e->warn_hook = NULL;

   if (NULL == (e->sl_optimize = SLang_pop_function ()))
     {
        slfe_deallocate (e);
        return NULL;
     }

   if (SLANG_NULL_TYPE == SLang_peek_at_stack())
     SLdo_pop();
   else if (NULL == (e->sl_set_options = SLang_pop_function ()))
     {
        slfe_deallocate (e);
        return NULL;
     }

   if (NULL == (e->option_string = isis_make_string (eng_name)))
     {
        slfe_deallocate (e);
        return NULL;
     }

   return e;
}

/*}}}*/

void add_slang_fit_engine_intrin (char *eng_name, char *stat_name) /*{{{*/
{
   (void) isis_fit_add_engine (eng_name, stat_name, add_slang_fit_engine);
}

/*}}}*/

/* FIXME!!! This is really ugly.... */
void set_slopt_fit_object (SLang_MMT_Type *mmt)
{
   Current_Fit_Object_MMT = mmt;
}

/* FIXME!!! This is really ugly.... */
int is_slang_optimizer (Isis_Fit_Type *ft) /*{{{*/
{
   return (ft->engine->method == &slfe_optimize);
}

/*}}}*/
