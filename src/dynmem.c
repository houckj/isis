
/* This file is a subset of routines from the JDMath library by
 * John E. Davis, modified by John C. Houck (1/28/99)
 * to work apart from JDMath
 *
 * Copyright (c) 1994, 1998 John E. Davis 
 *
 * You may distribute this file under the terms the GNU General Public
 * License.  See the file COPYING for more information.
 */

/* $Id: dynmem.c,v 1.2 2003/02/10 00:46:33 houck Exp $ */

#include "config.h"
#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "isis.h"
#include "isismath.h"

double *JDMdouble_vector (unsigned int n)
{
   return (double *) ISIS_MALLOC (n * sizeof(double));
}

void JDMfree_double_matrix (double **matrix, unsigned int n)
{
   unsigned int i;
   if (matrix == NULL) return;

   for (i = 0; i < n; i++)
     {
        ISIS_FREE (matrix[i]);
     }
   ISIS_FREE (matrix);
}

double **JDMdouble_matrix (unsigned int n, unsigned int m)
{
   double **matrix;
   unsigned int i;

   if (NULL == (matrix = (double **) ISIS_MALLOC (n * sizeof (double *))))
     return NULL;

   /* initialize everything to NULL */
   for (i = 0; i < n; i++) matrix[i] = NULL;

   for (i = 0; i < n; i++)
     {
        double *ptr;
        ptr = JDMdouble_vector (m);
        if (ptr == NULL)
          {
             JDMfree_double_matrix (matrix, n);
             return NULL;
          }
        matrix[i] = ptr;
     }

   return matrix;
}
