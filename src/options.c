/*  This file is part of ISIS, the Interactive Spectral Interpretation System
    Copyright (C) 1998-2012  Massachusetts Institute of Technology

    This software was developed by the MIT Center for Space Research under
    contract SV1-61010 from the Smithsonian Institution.

    Author:  John E. Davis  <davis@space.mit.edu>
             John C. Houck  <houck@space.mit.edu>

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

/* $Id: options.c,v 1.5 2004/02/09 11:14:23 houck Exp $ */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "isis.h"
#include "util.h"
#include "errors.h"

#define FIELD_SEP_CHAR        ';'

static char *isis_skip_whitespace (char *s)
{
   if (s == NULL)
     return s;

   while (isspace ((unsigned char)*s))
     s++;

   return s;
}

static char *skip_and_truncate_field (char *str, char delim)
{
   while (1)
     {
        char ch = *str;

        if (ch == 0)
          return str;

        if (ch == delim)
          {
             *str++ = 0;
             return str;
          }

        if (ch == '\\')
          {
             str++;
             if (*str == 0)
               {
                  *(str - 1) = 0;      /* kill trailing \ */
                  return str;
               }
          }
        str++;
     }
}

static void unescape_and_trim_string (char *s)
{
   char *p, *p1;
   char ch;

   if (s == NULL)
     return;

   /* First, trim off trailing whitespace.  Be careful to take \\ into account
    */
   p = s;
   p1 = NULL;

   while ((ch = *p) != 0)
     {
        p++;
        if (ch == '\\')
          {
             p1 = p;
             if (*p == 0)
               break;
             p++;
          }
     }
   while (ch != 0)
     {
     }

   /* Now p is at the end of the string, and p1 marks the last escaped
    * character (or NULL).  So, trim back to it.
    */
   while (p > s)
     {
        p--;
        if (0 == isspace ((unsigned char) *p))
          break;
        if (p == p1)
          break;
        *p = 0;
     }

   /* Now unescape */
   p = s;
   do
     {
        ch = *s++;
        if (ch == '\\')
          ch = *s++;

        *p++ = ch;
     }
   while (ch != 0);
}

void isis_free_options (Isis_Option_Type *t)
{
   if (t == NULL)
     return;

   /* NULLs ok */
   ISIS_FREE (t->subsystem);
   ISIS_FREE (t->option_names);
   ISIS_FREE (t->option_values);
   ISIS_FREE (t);
}

Isis_Option_Type *isis_parse_option_string (char *str)
{
   unsigned int num;
   Isis_Option_Type *t;
   char *estr;

   if (NULL == (t = (Isis_Option_Type *) ISIS_MALLOC (sizeof(Isis_Option_Type))))
     return NULL;

   memset ((char *) t, 0, sizeof (Isis_Option_Type));

   if (str == NULL)
     str = "";

   str = isis_skip_whitespace (str);
   if (NULL == (str = isis_make_string (str)))
     {
        ISIS_FREE (t);
        return NULL;
     }

   t->subsystem = str;
   /* Now skip to the first field */
   str = skip_and_truncate_field (str, FIELD_SEP_CHAR);
   unescape_and_trim_string (t->subsystem);

   /* Count the number of commas as an estimate on the number of fields */
   /* Allow for a final NULL. */
   num = 2;
   estr = str;
   while (NULL != (estr = strchr (estr, FIELD_SEP_CHAR)))
     {
        num++;
        estr++;
     }

   if (NULL == (t->option_names = (char **) ISIS_MALLOC (num * sizeof(char *))))
     {
        isis_free_options (t);
        return NULL;
     }

   if (NULL == (t->option_values = (char **) ISIS_MALLOC (num * sizeof(char *))))
     {
        isis_free_options (t);
        return NULL;
     }

   /* Ok, now get the options */
   num = 0;
   while (1)
     {
        str = isis_skip_whitespace (str);
        if (*str == 0)
          break;

        t->option_names[num] = str;
        estr = skip_and_truncate_field (str, FIELD_SEP_CHAR);
        /* Now look for = form */
        str = skip_and_truncate_field (str, '=');
        str = isis_skip_whitespace (str);
        if (*str == 0)
          t->option_values[num] = NULL;
        else
          t->option_values[num] = str;

        unescape_and_trim_string (t->option_names[num]);
        unescape_and_trim_string (t->option_values[num]);

        str = estr;
        num++;
     }

   t->option_values [num] = NULL;
   t->option_names [num] = NULL;
   t->num_options = num;
   return t;
}

static Isis_Option_Table_Type *find_option (char *name, Isis_Option_Table_Type *t)
{
   if (t == NULL)
     return NULL;

   while (t->optname != NULL)
     {
        if (0 == isis_strcasecmp (t->optname, name))
          return t;

        t++;
     }

   return NULL;
}

static int check_for_help (Isis_Option_Type *opt, Isis_Option_Table_Type *t)
{
   unsigned int n, i;
   char **names;

   names = opt->option_names;
   n = opt->num_options;

   for (i = 0; i < n; i++)
     {
        if (0 == isis_strcasecmp (names[i], "HELP"))
          {
             fprintf (stdout, "Specify \"%s;default\" to restore defaults.\n", opt->subsystem);
             fprintf (stdout, "Valid options for subsystem \"%s\" include:\n", opt->subsystem);
             while (t->optname != NULL)
               {
                  char *v;

                  if (t->value_flags & ISIS_OPT_REQUIRES_VALUE)
                    v = "=value";
                  else if (t->value_flags & ISIS_OPT_NO_VALUE)
                    v = "";
                  else
                    v = "[=value]";

                  fprintf (stdout, "%8s%-10s %s\n", t->optname, v,
                           (t->help_string ? t->help_string : ""));
                  if (t->default_value_string)
                    {
                       fprintf (stdout, "         default:  %-s='%-s'\n",
                                t->optname, t->default_value_string);
                    }

                  t++;
               }

             return 1;
          }
     }
   return 0;
}

static int set_defaults (char *subsystem, Isis_Option_Table_Type *t, void *client_data)
{
   while (t->optname != NULL)
     {
        char *value = t->default_value_string;
        char *name = t->optname;

        if (-1 == t->fun (subsystem, name, value, client_data))
          {
             static char *fmt = "processing \"%s\" option \"%s=%s\".";
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, fmt, subsystem, name,
                         (value ? value : "(null)"));
             return -1;
          }

        t++;
     }

   return 0;
}

static int check_for_keyword (Isis_Option_Type *opt, Isis_Option_Table_Type *t,
                              void *client_data,
                              char *keyword,
                              int (*action)(char *, Isis_Option_Table_Type *, void *))
{
   unsigned int n, i;
   char **names;

   (void) t;

   names = opt->option_names;
   n = opt->num_options;

   for (i = 0; i < n; i++)
     {
        if (0 == isis_strcasecmp (names[i], keyword))
          {
             names[i][0] = 0;
             return (*action)(opt->subsystem, t, client_data);
          }
     }

   return 0;
}

int isis_process_options (Isis_Option_Type *opt,
                          Isis_Option_Table_Type *table,
                          void *client_data,
                          int err_on_unsupported)
{
   char **names;
   char **values;
   char *subsystem;
   unsigned int i, n;

   if (opt == NULL)
     return 0;

   subsystem = opt->subsystem;
   names = opt->option_names;
   values = opt->option_values;
   n = opt->num_options;

   if ((table == NULL)
       && (opt->num_options != 0))
     {
        isis_vmesg (WARN, I_ERROR, __FILE__, __LINE__, "Subsystem \"%s\" does not permit options.", subsystem);
        if (err_on_unsupported)
          return -1;
        return 0;
     }

   if (check_for_help (opt, table))
     return -1;

   if (check_for_keyword (opt, table, client_data, "DEFAULT", &set_defaults))
     return 0;

   for (i = 0; i < n; i++)
     {
        Isis_Option_Table_Type *t;
        char *name;
        char *value;

        name = names[i];
        if (*name == 0)
          continue;

        if (NULL == (t = find_option (name, table)))
          {
             isis_vmesg (WARN, I_ERROR, __FILE__, __LINE__, "%s option `%s' is unknown.", subsystem, name);
             if (err_on_unsupported)
               return -1;
             continue;
          }

        if (t->fun == NULL)
          {
             isis_vmesg (WARN, I_NOT_IMPLEMENTED, __FILE__, __LINE__, "%s option `%s'.", subsystem, name);
             if (err_on_unsupported)
               return -1;
             continue;
          }

        value = values[i];

        if (value == NULL)
          {
             if (t->value_flags & ISIS_OPT_REQUIRES_VALUE)
               {
                  isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "\"%s\" option \"%s\" requires a value.",
                              subsystem, name);
                  return -1;
               }
          }
        else if (t->value_flags & ISIS_OPT_NO_VALUE)
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "\"%s\" option \"%s\" does not permit a value.",
                         subsystem, name);
             return -1;
          }

        if (-1 == t->fun (subsystem, name, value, client_data))
          {
             static char *fmt = "processing \"%s\" option \"%s=%s\".";
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, fmt, subsystem, name,
                         (value ? value : "(null)"));
             return -1;
          }
#if 0
        if (value != NULL)
          _ardlib_add_to_history ("%s used option %s=%s", subsystem, name, value);
        else
          _ardlib_add_to_history ("%s used option %s", subsystem, name);
#endif
     }

   return 0;
}

int isis_is_option_present (char *options, char *option, char **value)
{
   unsigned int i, num;
   Isis_Option_Type *opt;
   int status;
   char **option_names;

   if (options == NULL)
     return 0;

   if (option == NULL)
     return 0;

   opt = isis_parse_option_string (options);
   if (opt == NULL)
     return -1;

   option_names = opt->option_names;
   num = opt->num_options;

   status = 0;

   for (i = 0; i < num; i++)
     {
        if (0 != isis_strcasecmp (option, option_names[i]))
          continue;

        status = 1;

        if (value != NULL)
          {
             char *val = opt->option_values[i];

             if (val != NULL)
               {
                  if (NULL == (val = isis_make_string (val)))
                    status = -1;
               }
             *value = val;
          }
        break;
     }

   isis_free_options (opt);
   return status;
}

char *isis_add_option (char *subsystem, char *option)
{
   unsigned int len;
   char *s;

   if (option == NULL)
     option = "";
   if (subsystem == NULL)
     subsystem = "";

   len = strlen (subsystem) + strlen (option) + 2;
   if (NULL == (s = (char *) ISIS_MALLOC (len * sizeof(char))))
     return s;

   sprintf (s, "%s;%s", subsystem, option);
   return s;
}

static int get_opt_value (Isis_Option_Type *opt, char *option, char **vp)
{
   unsigned int i, num;
   char **option_names;

   if ((opt == NULL) || (option == NULL))
     return 0;

   option_names = opt->option_names;
   num = opt->num_options;

   for (i = 0; i < num; i++)
     {
        if (0 != isis_strcasecmp (option, option_names[i]))
          continue;

        if (NULL == (*vp = opt->option_values[i]))
          {
             isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "%s option %s requires a value",
                         opt->subsystem, option);
             return -1;
          }

        return 1;
     }
   return 0;
}

static int get_scalar_option (Isis_Option_Type *opt, char *option,
                       void *v, char *fmt, char *type)
{
   int status;
   char *s;

   if (1 != (status = get_opt_value (opt, option, &s)))
     return status;

   if (1 != sscanf (s, fmt, v))
     {
        isis_vmesg (FAIL, I_ERROR, __FILE__, __LINE__, "Unable to parse %s;%s option as %s.",
                    opt->subsystem, option, type);
        return -1;
     }

   return 1;
}

int isis_get_option_d (Isis_Option_Type *opt, char *option, double *d)
{
   return get_scalar_option (opt, option, (void *) d, "%lf", "a double");
}

int isis_get_option_i (Isis_Option_Type *opt, char *option, int *d)
{
   return get_scalar_option (opt, option, (void *) d, "%d", "an int");
}

int isis_get_option_u (Isis_Option_Type *opt, char *option, unsigned int *d)
{
   return get_scalar_option (opt, option, (void *) d, "%d", "an unsigned int");
}

char *isis_make_default_option_string (const char *subsystem, Isis_Option_Table_Type *table)
{
   Isis_Option_Table_Type *t;
   char *s, *ps;
   int len, o;

   /* start with subsystem name */
   len = strlen(subsystem);

   for (t = table; t->optname != NULL; t++)
     {
        /* semicolon + option name */
        len += strlen(t->optname) + 1;
        if (t->default_value_string != NULL)
          {
             /* equals-sign plus option value */
             len += strlen(t->default_value_string) + 1;
          }
     }

   /* trailing null */
   len += 1;

   if (NULL == (s = (char *) ISIS_MALLOC(len*sizeof *s)))
     return NULL;

   o = sprintf (s, "%s", subsystem);
   ps = s + o;

   for (t = table; t->optname != NULL; t++)
     {
        char *name = t->optname;
        char *value = t->default_value_string;

        if (value == NULL)
          o = sprintf (ps, ";%s", name);
        else
          o = sprintf (ps, ";%s=%s", name, value);

        ps += o;
     }

   return s;
}

int isis_update_option_string (char **optstring, char *optname, char *optvalue)
{
   char *start, *end, *newstring, *s;
   int len1, len2, len3, len4;

   if (optstring == NULL || optname == NULL)
     return -1;

   /* find option */
   if (NULL == (start = strstr (*optstring, optname)))
     return -1;
   /* NULL means this option appears last */
   end = strchr (start, ';');

   len1 = start - *optstring;
   len2 = strlen(optname) + 1;
   len3 = optvalue ? strlen(optvalue) + 1 : 0;
   len4 = end ? strlen (end) + 1: 0;

   if (NULL == (newstring = (char *) ISIS_MALLOC (len1 + len2 + len3 + len4 + 1)))
     return -1;

   s = isis_strcpy (newstring, *optstring, len1);
   s = isis_strcpy (s, optname, len2);
   if (len3 > 0)
     {
        sprintf (s, "=%s", optvalue);
        s += len3;
     }
   if (len4 > 0)
     {
        s = isis_strcpy (s, end, len4);
     }

   ISIS_FREE (*optstring);
   *optstring = newstring;

   return 0;
}
