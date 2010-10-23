/* -*- mode: C; mode: fold -*- */

/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2010 Massachusetts Institute of Technology

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

/* $Id: model.c,v 1.10 2004/02/09 11:14:23 houck Exp $ */

/*{{{ includes */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>

#include <slang.h>

#include "isis.h"
#include "db-atomic.h"
#include "db-em.h"
#include "db-cie.h"
#include "util.h"
#include "model.h"
#include "errors.h"

/*}}}*/

/*{{{ data types */

static Isis_Line_Profile_Type *Model_Profile = NULL;
static Isis_Line_Profile_Type map_thermal_profile;

/*}}}*/

/*{{{ new/free start/end */

static void free_model_node (Model_t *m) /*{{{*/
{
   if (m != NULL)
     {
        SLang_free_array (m->line_flux);
        m->line_flux = NULL;
        SLang_free_array (m->last_ionpop);
        m->last_ionpop = NULL;
     }
   ISIS_FREE (m);
}

/*}}}*/

static Model_t *new_model_node (void) /*{{{*/
{
   Model_t *m;

   if (NULL == (m = (Model_t *) ISIS_MALLOC (sizeof(Model_t))))
     return NULL;

   memset ((char *) m, 0, sizeof(*m));

   m->id = 1;
   m->next = NULL;
   m->metal_abund = 1.0;
   m->redshift = 0.0;
   m->line_flux = NULL;
   m->last_ionpop = NULL;

   memset ((char *) m->rel_abund, 0, (ISIS_MAX_PROTON_NUMBER+1)*sizeof(float));

   return m;
}

/*}}}*/

static void free_model_list (Model_t *head) /*{{{*/
{
   if (head == NULL)
     return;

   while (head != NULL)
     {
        Model_t *h = head->next;
        free_model_node (head);
        head = h;
     }
}

/*}}}*/

int Model_start (void) /*{{{*/
{
   return 0;
}

/*}}}*/

void Model_end (Model_t *m) /*{{{*/
{
   free_model_list (m);
}

/*}}}*/

static int append_node (Model_t *m, Model_t *h) /*{{{*/
{
   Model_t *p;

   if (h == NULL || m == NULL)
     return -1;

   for (p = h; p->next != NULL; p = p->next)
     {
     }

   p->next = m;
   m->next = NULL;
   m->id = p->id + 1;

   return m->id;
}

/*}}}*/

/*}}}*/

/*{{{ variable abundances */

static int set_rel_abund (Model_t *m, float metal_abund, float *rel_abund,  /*{{{*/
                          unsigned int *element, unsigned int num_rel_abund)
{
   unsigned int i, Z;

   if (NULL == m)
     return -1;

   m->rel_abund[1] = 1.0;   /* H is cosmic */
   m->rel_abund[2] = 1.0;   /* He is cosmic */

   for (Z = 3; Z <= ISIS_MAX_PROTON_NUMBER; Z++)
     m->rel_abund[Z] = metal_abund;

   if ((rel_abund == NULL) || (element == NULL))
     return 0;

   for (i = 0; i < num_rel_abund; i++)
     {
        Z = element[i];
        if (Z <= ISIS_MAX_PROTON_NUMBER)
          m->rel_abund[Z] = rel_abund[i];
        else isis_vmesg (FAIL, I_INVALID, __FILE__, __LINE__, "Z = %d", Z);
     }

   return 0;
}

/*}}}*/

Model_t *Model_add_component (Model_t *head, Model_t *x,  /*{{{*/
                              unsigned int *elem, float *elem_abund, unsigned int num_elems)
{
   Model_t *m = new_model_node ();

   if ((m == NULL) || (x == NULL))
     return NULL;

   m->temperature = x->temperature;
   m->density = x->density;
   m->vturb = x->vturb * 1.e5;    /* km/s -> cm/s */
   m->norm = x->norm * 1.e14;
   m->redshift = x->redshift;
   m->metal_abund = x->metal_abund;

   if (-1 == set_rel_abund (m, x->metal_abund, elem_abund, elem, num_elems))
     {
        free_model_node (m);
        return NULL;
     }

   if (head == NULL)
     head = m;
   else if (-1 == append_node (m, head))
     {
        free_model_node (m);
        return NULL;
     }

   return head;
}

/*}}}*/

static int scan_for_specific_abundances (FILE *fp, char *buf, int size, Model_t *m) /*{{{*/
{
   if (m == NULL || size < BUFSIZE)
     return -1;

   /* Scan lines looking for keyword-value pairs of the form
    *         cc = float
    *     OR   c = float
    *
    * Skip blank lines and lines with '#' in column 1.
    *
    * Return 0 if a line's first non-blank char is a digit.
    *
    * Return 1 if corrupt input on a line (telling caller to
    *          read another line).
    *
    * Return -1 on input error
    *
    */

   while (!feof(fp))
     {
        char *p;

        if (NULL == fgets(buf, size, fp))
          {
             if (feof (fp)) return 0;
             isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "abundances");
             return -1;
          }

        if (buf[0] == COMMENT_CHAR || buf[0] == '\n')
          continue;

        p = buf;

        while (*p != '\0' && (*p == ' ' || *p == '\t' || *p == '\n'))
             p++;

        /* skip blank lines */
        if (*p == '\0')
          continue;

           if (isdigit((int) *p))
          return 0;

        /* skip lines with garbage characters */
        if (!isalpha((int) *p))
          continue;

        do
          {
             char el[3];
             float abund;
             int off, Z, bogus;

             el[0] = *p++;
             el[1] = isalpha((int) *p) ? *p++ : '\0';
             el[2] = '\0';

             bogus = _DB_get_element_Z (&Z, el);
             if (bogus)
               isis_vmesg (WARN, I_WARNING, __FILE__, __LINE__, "ignoring unrecognized element:  `%s'", el);

             while (*p != '\0' &&
                    (*p == ' ' || *p == '\t' || *p == '\n' || *p == '='))
               p++;

             if (*p == '\0')
               break;

             if (!isdigit((int) *p))
                  return 1;
             else if (1 != sscanf (p, "%f%n", &abund, &off))
               {
                  isis_vmesg (FAIL, I_READ_FAILED, __FILE__, __LINE__, "%s abundance", el);
                  return 1;
               }

             if (!bogus)
               m->rel_abund[Z] = abund;

             p += off;

             while (*p != '\0' && (*p == ' ' || *p == '\t' || *p == '\n'))
               p++;

          } while (*p != '\0');
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ load/print */

Model_t *Model_load_ascii_file (char *filename) /*{{{*/
{
   FILE *fp;
   Model_t *head, *m;
   int need_new_buf = 1;
   int failed = -1;

   if ((0 == is_regular_file (filename))
       || (NULL == (fp = fopen (filename, "r"))))
     {
        isis_vmesg (FAIL, I_READ_OPEN_FAILED, __FILE__, __LINE__, "%s", filename);
        return NULL;
     }

   head = NULL;
   while (!feof(fp))
     {
        float temp, dens, metal_abund, redshift, vturb;
        double norm;
        char buf[BUFSIZE];

        if (need_new_buf
            && NULL == fgets(buf, BUFSIZE, fp))
          goto fail;

        if (buf[0] == COMMENT_CHAR)
          continue;

        if (6 != sscanf (buf, "%*d  %e %e %e %le %e %e",
                         &temp, &dens, &metal_abund, &norm, &vturb, &redshift))
          continue;

        m = new_model_node ();
        if (NULL == m)
          goto fail;

        m->temperature = temp;
        m->density = dens;
        m->vturb = vturb * 1.e5;            /* km/s -> cm/s */
        m->norm = norm * 1.e14;
        m->redshift = redshift;
        m->metal_abund = metal_abund;

        if (-1 == set_rel_abund (m, metal_abund, NULL, NULL, 0))
          goto fail;

        need_new_buf = scan_for_specific_abundances (fp, buf, (int) sizeof(buf), m);
        if (-1 == need_new_buf)
          goto fail;

        if (head == NULL)
          head = m;
        else
          {
             if (append_node (m, head) < 0)
               goto fail;
          }
     }

   failed = 0;

   fail:

   if (failed)
     {
        isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "loading model from %s", filename);
        Model_end (head);
        head = NULL;
     }

   fclose (fp);
   return head;
}

/*}}}*/

static int print_model (FILE *fp, Model_t *m) /*{{{*/
{
   enum {N_PER_LINE = 5};
   int Z, k;

   fprintf (fp, "%4d  %11.4e  %11.4e  %11.4e  %11.4e %11.4e %11.4e\n",
            m->id,
            m->temperature,
            m->density,
            m->metal_abund,
            m->norm * 1.e-14 ,
            m->vturb * 1.e-5,           /* cm/s -> km/s */
            m->redshift);

   k = 0;
   for (Z = 3; Z <= ISIS_MAX_PROTON_NUMBER; Z++)
     {
        char el[3];
        if (fabs(m->rel_abund[Z] - m->metal_abund) < 100 * FLT_EPSILON * m->metal_abund)
          continue;
        if (-1 == _DB_get_element_name (el, Z))
          return -1;
        k++;
        fprintf (fp, "  %2s=%6.4f %c",
                 el, m->rel_abund[Z], (k % N_PER_LINE) ? ' ' : '\n');
     }

   if (k % N_PER_LINE != 0)
     fputc ('\n', fp);

   fflush(fp);

   return 0;
}

/*}}}*/

int Model_print_model (FILE *fp, Model_t *h) /*{{{*/
{
   Model_t *m;

   if (NULL == fp)
     return -1;

   fprintf (fp, "# id    Temp         Dens         Abund        Norm        Vturb       redshift\n");
   fprintf (fp, "#       (K)          (cm^-3)                               (km/s)              \n");

   if (NULL == h)
     return 0;

   for (m = h; m != NULL; m = m->next)
     {
        if (-1 == print_model (fp, m))
          return -1;
     }

   return 0;
}

/*}}}*/

/*}}}*/

/*{{{ spectrum calculation */

void Model_set_profile_function (int idx) /*{{{*/
{
   const char *msg;
   switch (idx)
     {
      case 0:
        Model_Profile = NULL;
        msg = "delta-function";
        break;
      default:
        Model_Profile = map_thermal_profile;
        msg = "thermal";
        break;
     }
   isis_vmesg (INFO, I_INFO, __FILE__, __LINE__, "Use %s line profiles", msg);
}

/*}}}*/

static int profile_thermal (double *f, double wllo, double wlhi, /*{{{*/
                            double wl0, double atwt, double *params, int num_params,
                            void *cl)
{
   double temperature = params[0];
   double vturb = params[1];
   double sigma;
   double dxl, dxh;

   (void) num_params; (void) cl;

   /* use wl0 as the line center rather than the value listed
    * in the DB_line_t struct -- might be redshifted.
    */

   sigma = (wl0 / CLIGHT) * sqrt (BOLTZ * temperature / atwt / AMU
                                  + 0.5 * vturb * vturb);

   dxh = (wlhi - wl0) / sigma;
   dxl = (wllo - wl0) / sigma;

   *f = isis_gpf (dxh) - isis_gpf (dxl);

   return 0;
}

/*}}}*/

#define FAINT_TOL 1.e-4
/* Line wing faintness limit.  Larger values should make the code run
 * faster but may underestimate the total contribution of faint
 * line wings.   Smaller will take more time, but should be more accurate.
 */

static int map_thermal_profile (Isis_Hist_t *g, double flux, double wl, double atomic_weight, int mid, /*{{{*/
                                double *profile_params, int num_profile_params, void *profile_options)
{
   double f, av, de, faint;
   double *wllo, *wlhi, *val;
   int i, nbins;

   if (g == NULL)
     return -1;

   wllo = g->bin_lo;
   wlhi = g->bin_hi;
   val = g->val;
   nbins = g->nbins;

   for (i = mid; i >= 0; i--)
     {
        if (-1 == profile_thermal (&f, wllo[i], wlhi[i], wl, atomic_weight,
                                   profile_params, num_profile_params,
                                   profile_options))
          return -1;
        de = f * flux;
        av = fabs (val[i]);
        faint = (av > 0.0) && (fabs(de) < FAINT_TOL * av);
        val[i] += de;
        if (faint)
          break;
     }

   for (i = mid+1; i < nbins; i++)
     {
        if (-1 == profile_thermal (&f, wllo[i], wlhi[i], wl, atomic_weight,
                                   profile_params, num_profile_params,
                                   profile_options))
          return -1;
        de = f * flux;
        av = fabs (val[i]);
        faint = (av > 0.0) && (fabs(de) < FAINT_TOL * av);
        val[i] += de;
        if (faint)
          break;
     }

   return 0;
}

/*}}}*/

typedef struct
{
   double temperature;
   double ndensity;
}
Plasma_State_Type;

static SLang_CStruct_Field_Type Plasma_State_Layout [] =
{
   MAKE_CSTRUCT_FIELD (Plasma_State_Type, temperature, "temperature", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD (Plasma_State_Type, ndensity, "ndensity", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int apply_line_modifier (Model_t *m, Model_Info_Type *info, DB_line_t *line, double *emis) /*{{{*/
{
   Plasma_State_Type lm_state;

   lm_state.temperature = m->temperature;
   lm_state.ndensity = m->density;

   SLang_start_arg_list ();
   SLang_push_array (info->line_emis_modifier_params, 0);
   SLang_push_int (line->indx);
   if (-1 == SLang_push_cstruct ((VOID_STAR)&lm_state, Plasma_State_Layout))
     return -1;
   SLang_push_double (*emis);
   if (info->line_emis_modifier_args != NULL)
     isis_push_args (info->line_emis_modifier_args);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (info->line_emis_modifier))
     return -1;

   if (-1 == SLang_pop_double (emis))
     return -1;

   return 0;
}

/*}}}*/

static int call_ionpop_modifier (Model_t *m, Model_Info_Type *info, float *ionpop_new) /*{{{*/
{
   Plasma_State_Type s;
   SLang_Array_Type *sl_ionpop = NULL;
   int Z, q, status = -1;
   int n = ISIS_MAX_PROTON_NUMBER;

   s.temperature = m->temperature;
   s.ndensity = m->density;

   /* Float_Type[n,n] = ionpop_modifier (params, state, last_ionpop, [,args]) */
   SLang_start_arg_list ();
   SLang_push_array (info->ionpop_params, 0);
   if (-1 == SLang_push_cstruct ((VOID_STAR)&s, Plasma_State_Layout))
     {
        SLang_end_arg_list ();
        return -1;
     }
   SLang_push_array (m->last_ionpop, 0);
   if (info->ionpop_args != NULL)
     isis_push_args (info->ionpop_args);
   SLang_end_arg_list ();

   if (-1 == SLexecute_function (info->ionpop_modifier))
     return -1;

   if (-1 == SLang_pop_array_of_type (&sl_ionpop, SLANG_FLOAT_TYPE))
     return -1;

   if ((sl_ionpop == NULL)
       || (sl_ionpop->num_dims != 2)
       || ((sl_ionpop->dims[0] != n+1) || (sl_ionpop->dims[1] != n+1)))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__,
                    "ionpop_modifier: invalid return value, expecting: Float_Type[n,n] with n=%d", n+1);
        goto return_status;
     }

   for (Z = 1; Z <= n; Z++)
     {
        for (q = 0; q <= Z; q++)
          {
             int i[2];
             i[0] = Z;  i[1] = q;
             if (-1 == SLang_get_array_element (sl_ionpop, i, &ionpop_new[Z*(n+1)+q]))
               goto return_status;
          }
     }

   SLang_free_array (m->last_ionpop);
   m->last_ionpop = sl_ionpop;

   status = 0;
return_status:
   if (status != 0)
     {
       SLang_free_array (sl_ionpop);
     }

   return status;
}

/*}}}*/

static int add_spread_lines (double *val, double *wllo, double *wlhi, int nbins, /*{{{*/
                             EM_line_emis_t *t, Model_t *m, Model_Info_Type *info)
{
   Isis_Line_Profile_Type *map_profile;
   double thermal_profile_params[2];
   double *profile_params = NULL;
   int num_profile_params = 0;
   void *profile_options = NULL;
   int k, nlines;
   float atwt;
   Isis_Hist_t g;

   if (NULL == wllo || NULL == wlhi || NULL == val || NULL == t
       || nbins < 1)
     return -1;

   memset ((char *)&g, 0, sizeof g);
   g.bin_lo = wllo;
   g.bin_hi = wlhi;
   g.val = val;
   g.nbins = nbins;

   /* model-specific line profile over-rides the "global" setting */
   if (info->profile)
     {
        SLang_Array_Type *pa = info->profile_params;
        if (pa != NULL)
          {
             profile_params = (double *) pa->data;
             num_profile_params = pa->num_elements;
          }
        map_profile = info->profile;
        profile_options = info->profile_options;
     }
   else
     {
        map_profile = Model_Profile;
        if (map_profile != NULL)
          {
             thermal_profile_params[0] = (double) m->temperature;
             thermal_profile_params[1] = (double) m->vturb;
             profile_params = thermal_profile_params;
             num_profile_params = 2;
          }
        else
          {
             profile_params = NULL;
             num_profile_params = 0;
          }
        profile_options = NULL;
     }

   nlines = EM_get_nlines (t);
   if (-1 == nlines)
     return -1;

   for (k = 0; k < nlines; k++)
     {
        DB_line_t *line;
        double flux, emis;
        float emis_f, wl, flux_f;
        int mid, Z, q;

        if (-1 == _EM_get_line_emis_wl (&line, &emis_f, &wl, k, t))
          return -1;
        emis = (double) emis_f;

        if (wl < wllo[0] || wlhi[nbins-1] < wl)
          continue;

        if (info->line_emis_modifier != NULL)
          {
             if (-1 == apply_line_modifier (m, info, line, &emis))
               return -1;
          }

        if (emis <= 0.0)
          continue;

        mid = find_bin ((double) wl, wllo, wlhi, nbins);
        if (mid < 0)
          continue;

        flux = m->norm * emis;

        if (-1 == DB_get_line_ion (&Z, &q, line))
          return -1;

        flux *= m->rel_abund[Z];

        /* Side-effect: increment line fluxes stored in wavelength tables: */
        line->flux += flux;

        /* each model component remembers its contribution to the line flux */
        flux_f = (float) flux;
        SLang_set_array_element (m->line_flux, &line->indx, &flux_f);

        /* profile = NULL imples a delta function. */
        if (NULL == map_profile)
          {
             val[mid] += flux;
             continue;
          }
        else
          {
             if (-1 == DB_get_atomic_weight_amu (&atwt, line))
               return -1;

             if (-1 == (*map_profile)(&g, flux, (double) wl, (double) atwt, mid,
                                      profile_params, num_profile_params, profile_options))
               return -1;
          }
     }

   return 0;
}

/*}}}*/

static float lorentz_gamma (float z) /*{{{*/
{
   float x = z + 1.0;
   float beta = (x*x - 1.0) / (x*x + 1.0);
   return (1.0 / sqrt (1.0 - beta*beta));
}

/*}}}*/

static int shift_grid_to_emitter_frame (EM_cont_type_t *cont, double *wllo, double *wlhi, float redshift) /*{{{*/
{
   float fac;
   int i;

   if ((cont == NULL) || (redshift == -1.0))
     return -1;

   /* shift from lab to (moving) emitter frame */
   fac = 1.0 / (1.0 + redshift);

   for (i = 0; i < cont->nbins; i++)
     {
        cont->wllo[i] = fac * wllo[i];
        cont->wlhi[i] = fac * wlhi[i];
     }

   return 0;
}

/*}}}*/

int Model_spectrum (Model_t *h, Model_Info_Type *info, /*{{{*/
                    double *wllo, double *wlhi, int nbins, double *val)
{
   Model_t *m;
   EM_cont_type_t *cont = NULL;
   double *tmp_val = NULL;
   float *ionpop_new = NULL;
   char *flag = NULL;
   int i, cont_nbins, include_lines, include_contin, ret = -1;
   SLindex_Type db_nlines;

   if (NULL == h || NULL == info)
     return -1;

   if (NULL == val || NULL == wllo || NULL == wlhi || nbins < 1)
     return -1;

   if (NULL == info->em || NULL == info->db)
     return -1;

   switch (info->contrib_flag)
     {
      case MODEL_LINES_AND_CONTINUUM:
        include_lines = 1;
        include_contin = 1;
        break;

      case MODEL_LINES:
        include_lines = 1;
        include_contin = 0;
        break;

      case MODEL_CONTIN:
      case MODEL_CONTIN_PSEUDO:
      case MODEL_CONTIN_TRUE:
        include_lines = 0;
        include_contin = 1;
        break;

      default:
        isis_vmesg (WARN, I_INVALID, __FILE__, __LINE__, "contrib_flag=%d; using default",
                    info->contrib_flag);
        info->contrib_flag = MODEL_LINES_AND_CONTINUUM;
        include_lines = 1;
        include_contin = 1;
        break;
     }

   memset ((char *)val, 0, nbins * sizeof(double));

   /* add_spread_lines has side effect of incrementing these
    * fluxes -- so zero them out first.
    */

   if (-1 == DB_zero_line_flux (info->db))
     return -1;

   if (-1 == (db_nlines = DB_get_nlines (info->db)))
     return -1;

   if ((info->line_list != NULL) && (info->line_list->data_type != SLANG_NULL_TYPE))
     {
        int *line_list = (int *) info->line_list->data;
        int num_lines = info->line_list->num_elements;
        if (num_lines > 0)
          {
             flag = DB_flag_array_from_list (line_list, num_lines, info->db);
             if (NULL == flag)
               return -1;
          }
     }

   if (info->ionpop_modifier != NULL)
     {
        int n = ISIS_MAX_PROTON_NUMBER+1;
        if (NULL == (ionpop_new = ISIS_MALLOC (n*n*sizeof(float))))
          goto finish;
        memset ((char *)ionpop_new, 0, n*n*sizeof(float));
     }

   if ((NULL == (cont = EM_new_continuum (nbins)))
       || (NULL == (tmp_val = (double *) ISIS_MALLOC (nbins * sizeof(double)))))
     goto finish;

   cont_nbins = cont->nbins;

   /* 1) Shift the input observer frame grid into the source frame
    * 2) Compute the emissivity in each rest-frame bin
    * 3) Shift the emissivity in each bin back to the observer frame,
    *    (just the time-dilation factor, since I'm assuming that
    *     the model values are photons/sec/whatever
    *                      NOT ergs/sec/whatever!).
    */

   for (m = h; m != NULL; m = m->next)
     {
        EM_line_emis_t *emis_list = NULL;
        double *val_ptr;
        float par[2];
        par[0] = m->temperature;
        par[1] = m->density;

        if (isis_user_break())
          {
             ret = 0;
             goto finish;
          }

        if (m->line_flux == NULL)
          {
             if (NULL == (m->line_flux = SLang_create_array (SLANG_FLOAT_TYPE, 1, NULL, &db_nlines, 1)))
               goto finish;
          }
        memset ((char *)m->line_flux->data, 0, db_nlines * sizeof(float));

        if (m->norm == 0)
          continue;

        if (m->redshift != 0.0)
          {
             memset ((char *)tmp_val, 0, nbins * sizeof(double));
             val_ptr = tmp_val;
          }
        else val_ptr = val;

        if (-1 == shift_grid_to_emitter_frame (cont, wllo, wlhi, m->redshift))
          goto finish;

        if (info->ionpop_modifier != NULL)
          {
             if (-1 == call_ionpop_modifier (m, info, ionpop_new))
               goto finish;
          }

        if (include_lines)
          {
             emis_list = EM_get_line_spectrum (flag, par, ionpop_new, info->em);
             if (NULL == emis_list)
               goto finish;

             if (-1 == add_spread_lines (val_ptr, cont->wllo, cont->wlhi, cont->nbins,
                                         emis_list, m, info))
               goto finish;

             EM_free_line_emis_list (emis_list);
             emis_list = NULL;
          }

        if (include_contin)
          {
             EM_cont_select_t s;
             double m_norm;
             double *c_val, *p_val;

             s.Z = 0; s.q = -1;  s.rel_abun = m->rel_abund;
             if (-1 == EM_get_continuum (cont, par, ionpop_new, &s, info->em))
               goto finish;

             m_norm = m->norm;

             switch (info->contrib_flag)
               {
                case MODEL_CONTIN_TRUE:
                  c_val = cont->true_contin;
                  for (i = 0; i < cont_nbins; i++)
                    val_ptr[i] += m_norm * c_val[i];
                  break;

                case MODEL_CONTIN_PSEUDO:
                  c_val = cont->pseudo;
                  for (i = 0; i < cont_nbins; i++)
                    val_ptr[i] += m_norm * c_val[i];
                  break;

                default:
                  c_val = cont->true_contin;
                  p_val = cont->pseudo;
                  for (i = 0; i < cont_nbins; i++)
                    val_ptr[i] += m_norm * (c_val[i] + p_val[i]);
                  break;
               }
          }

        if (m->redshift != 0.0)
          {
             float gm = lorentz_gamma (m->redshift);
             for (i = 0; i < cont_nbins; i++)
               val[i] += tmp_val[i] / gm;
          }
     }

   ret = 0;

   finish:

   if (ret)
     isis_vmesg (FAIL, I_FAILED, __FILE__, __LINE__, "computing model spectrum");

   ISIS_FREE (ionpop_new);
   ISIS_FREE (flag);
   ISIS_FREE (tmp_val);
   EM_free_continuum (cont);

   return ret;
}

/*}}}*/

/*}}}*/
