/* cfitsio interface */
/*
  Copyright (c) 2003-2009 Massachusetts Institute of Technology

  This software was developed by the MIT Center for Space Research
  under contract SV1-61010 from the Smithsonian Institution.

  Permission to use, copy, modify, distribute, and sell this software
  and its documentation for any purpose is hereby granted without fee,
  provided that the above copyright notice appear in all copies and
  that both that copyright notice and this permission notice appear in
  the supporting documentation, and that the name of the Massachusetts
  Institute of Technology not be used in advertising or publicity
  pertaining to distribution of the software without specific, written
  prior permission.  The Massachusetts Institute of Technology makes
  no representations about the suitability of this software for any
  purpose.  It is provided "as is" without express or implied warranty.

  THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS.  IN NO EVENT SHALL THE MASSACHUSETTS
  INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY SPECIAL, INDIRECT OR
  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
  OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
  NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
  WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <slang.h>

#include <errno.h>

#include "cfitsio.h"

#ifdef __cplusplus
extern "C"
{
#endif
SLANG_MODULE(cfitsio);
#ifdef __cplusplus
}
#endif

#include "version.h"

/* This is a hack that works for all 32 and 64 bit systems that I know */
/* The CFITSIO_INT*_TYPE objects must refer to the corresponding C type
 * and have the approriate length.
 */
#if SIZEOF_SHORT == 2
# define SLANG_UINT16_TYPE SLANG_USHORT_TYPE
# define SLANG_INT16_TYPE SLANG_SHORT_TYPE
# define CFITSIO_INT16_TYPE	TSHORT
# define CFITSIO_UINT16_TYPE	TUSHORT
typedef short int16_type;
#endif

#if SIZEOF_INT == 4
# define SLANG_INT32_TYPE	SLANG_INT_TYPE
# define SLANG_UINT32_TYPE	SLANG_UINT_TYPE
# define CFITSIO_INT32_TYPE	TINT
# define CFITSIO_UINT32_TYPE	TUINT
typedef int int32_type;
#endif

#ifdef TLONGLONG
# if (SIZEOF_LONG == 8)
#  define SLANG_INT64_TYPE	SLANG_LONG_TYPE
#  define CFITSIO_INT64_TYPE	TLONG
# else
#  define SLANG_INT64_TYPE	SLANG_LLONG_TYPE
#  define CFITSIO_INT64_TYPE	TLONGLONG
# endif
#endif

typedef struct
{
   fitsfile *fptr;
}
FitsFile_Type;

static SLtype Fits_Type_Id = 0;

/* This routine is used for binary tables --- not keywords.  For a binary table,
 * TLONG always specifies a 32 bit integer, but for a keyword is simply means
 * a long integer.
 */
static int map_fitsio_type_to_slang (int *typep, long *repeat, SLtype *stype)
{
   int type = *typep;
   int sgn = (type < 0) ? -1 : 1;

   /* Variable length objects have negative type values */
   if (sgn == -1)
     type = -type;

   switch (type)
     {
      case TSHORT:		       /* cfitsio 16 bit type */
	*stype = SLANG_INT16_TYPE;
	*typep = sgn * CFITSIO_INT16_TYPE;
	break;

      case TUSHORT:		       /* cfitsio 16 bit type */
	*stype = SLANG_UINT16_TYPE;
	*typep = sgn * CFITSIO_UINT16_TYPE;
	break;

      case TINT:
	*stype = SLANG_INT_TYPE;
	break;

#ifdef TUINT
      case TUINT:
	*stype = SLANG_UINT_TYPE;
	break;
#endif

      case TLONG:		       /* cfitsio 32 bit value */
	*stype = SLANG_INT32_TYPE;
	*typep = sgn * CFITSIO_INT32_TYPE;
	break;

      case TULONG:
	*stype = SLANG_UINT32_TYPE;
	*typep = sgn * CFITSIO_UINT32_TYPE;
	break;

#ifdef TLONGLONG
      case TLONGLONG:
	*stype = SLANG_INT64_TYPE;
	break;
#endif
      case TDOUBLE:
	*stype = SLANG_DOUBLE_TYPE;
	break;

      case TFLOAT:
	*stype = SLANG_FLOAT_TYPE;
	break;

      case TLOGICAL:
      case TBYTE:
	*stype = SLANG_UCHAR_TYPE;
	break;
      case TBIT:
	/* Make sure these all map to SIGNED types -- not unsigned.  This
	 * way they will be written out as SIGNED types and avoid the
	 * bit corruption that takes place when cfitsio adds, e.g., 0x7FFFFFFF
	 * to push the value into the unsigned range.
	 */
	if (*repeat <= 8)
	  {
	     *repeat = 1;
	     *stype = SLANG_CHAR_TYPE;
	     break;
	  }
	if (*repeat <= 16)
	  {
	     *repeat = 1;
	     *stype = SLANG_INT16_TYPE;
	     break;
	  }
	if (*repeat <= 32)
	  {
	     *repeat = 1;
	     *stype = SLANG_INT32_TYPE;
	     break;
	  }
	SLang_verror (SL_NOT_IMPLEMENTED, "bit type %ldX is not supported", *repeat);
	return -1;

      case TSTRING:
	*stype = SLANG_STRING_TYPE;
	break;

      default:
	SLang_verror (SL_NOT_IMPLEMENTED, "Fits column type %d is not supported",
		      type);
	return -1;
     }

   return 0;
}

static int open_file (SLang_Ref_Type *ref, char *filename, char *mode)
{
   fitsfile *fptr;
   int status;
   FitsFile_Type *ft;
   SLang_MMT_Type *mmt;

   if (-1 == SLang_assign_to_ref (ref, SLANG_NULL_TYPE, NULL))
     return -1;

   status = 0;
   fptr = NULL;
   switch (*mode)
     {
      case 'r':
	(void) fits_open_file (&fptr, filename, READONLY, &status);
	break;

      case 'w':
	(void) fits_open_file (&fptr, filename, READWRITE, &status);
	break;

      case 'c':
	if ((-1 == remove (filename))
	    && (errno != ENOENT))
	  {
	     SLang_verror (SL_OBJ_NOPEN, "Unable to create a new version of %s--- check permissions", filename);
	     return -1;
	  }
	(void) fits_create_file (&fptr, filename, &status);
	break;

      default:
	SLang_verror (SL_INVALID_PARM, "fits_open_file: iomode \"%s\" is invalid", mode);
	return -1;
     }

   if (status)
     return status;

   if (fptr == NULL)
     return -1;

   ft = (FitsFile_Type *) SLmalloc (sizeof (FitsFile_Type));
   if (ft == NULL)
     {
	fits_close_file (fptr, &status);
	return -1;
     }
   memset ((char *) ft, 0, sizeof (FitsFile_Type));

   ft->fptr = fptr;

   if (NULL == (mmt = SLang_create_mmt (Fits_Type_Id, (VOID_STAR) ft)))
     {
	fits_close_file (fptr, &status);
	SLfree ((char *) fptr);
	return -1;
     }

   if (-1 == SLang_assign_to_ref (ref, Fits_Type_Id, &mmt))
     {
	SLang_free_mmt (mmt);	       /* This will close the file */
	return -1;
     }

   return status;
}

static int delete_file (FitsFile_Type *ft)
{
   int status = 0;

   if (ft->fptr != NULL)
     fits_delete_file (ft->fptr, &status);
   ft->fptr = NULL;
   return status;
}

static int close_file (FitsFile_Type *ft)
{
   int status = 0;

   status = 0;
   if (ft->fptr != NULL)
     {
	(void) fits_close_file (ft->fptr, &status);
	ft->fptr = NULL;
     }
   return status;
}

static int movnam_hdu (FitsFile_Type *ft, int *hdutype, char *extname, int *extvers)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_movnam_hdu (ft->fptr, *hdutype, extname, *extvers, &status);
}

static int movabs_hdu (FitsFile_Type *ft, int *n)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_movabs_hdu (ft->fptr, *n, NULL, &status);
}

static int movrel_hdu (FitsFile_Type *ft, int *n)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_movrel_hdu (ft->fptr, *n, NULL, &status);
}

static int get_num_hdus (FitsFile_Type *ft, SLang_Ref_Type *ref)
{
   int status = 0;
   int num;

   if (ft->fptr == NULL)
     return -1;

   if (0 == fits_get_num_hdus (ft->fptr, &num, &status))
     {
	if (-1 == SLang_assign_to_ref (ref, SLANG_INT_TYPE, &num))
	  return -1;
     }

   return status;
}

static int get_hdu_num (FitsFile_Type *ft)
{
   int num;
   if (ft->fptr == NULL)
     return -1;
   return fits_get_hdu_num (ft->fptr, &num);
}

static int get_hdu_type (FitsFile_Type *ft, SLang_Ref_Type *ref)
{
   int hdutype;
   int status = 0;

   if (ft->fptr == NULL)
     return -1;
   if (0 == fits_get_hdu_type (ft->fptr, &hdutype, &status))
     {
	if (-1 == SLang_assign_to_ref (ref, SLANG_INT_TYPE, &hdutype))
	  return -1;
     }
   return status;
}

static int copy_file (FitsFile_Type *ft, FitsFile_Type *gt,
		      int *prev, int *cur, int *next)
{
   int status = 0;

   if ((ft->fptr == NULL) || (gt->fptr == NULL))
     return -1;
#ifndef fits_copy_file
   (void) status; (void) prev; (void) cur; (void) next;
   SLang_verror (SL_NOT_IMPLEMENTED, "Not supported by this version of cfitsio");
   return -1;
#else
   return fits_copy_file (ft->fptr, gt->fptr, *prev, *cur, *next, &status);
#endif
}

static int copy_hdu (FitsFile_Type *ft, FitsFile_Type *gt, int *morekeys)
{
   int status = 0;

   if ((ft->fptr == NULL) || (gt->fptr == NULL))
     return -1;

   return fits_copy_hdu (ft->fptr, gt->fptr, *morekeys, &status);
}

static int copy_header (FitsFile_Type *ft, FitsFile_Type *gt)
{
   int status = 0;

   if ((ft->fptr == NULL) || (gt->fptr == NULL))
     return -1;

   return fits_copy_header (ft->fptr, gt->fptr, &status);
}

static int delete_hdu (FitsFile_Type *ft)
{
   int status = 0;

   if (ft->fptr == NULL)
     return -1;

   return fits_delete_hdu (ft->fptr, NULL, &status);
}

static int pop_string_or_null (char **s)
{
   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
	*s = NULL;
	return SLang_pop_null ();
     }

   return SLang_pop_slstring (s);
}

static int pop_array_or_null (SLang_Array_Type **a)
{
   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
	*a = NULL;
	return SLang_pop_null ();
     }
   return SLang_pop_array (a, 1);
}

static FitsFile_Type *pop_fits_type (SLang_MMT_Type **mmt)
{
   FitsFile_Type *ft;

   if (NULL == (*mmt = SLang_pop_mmt (Fits_Type_Id)))
     return NULL;

   if (NULL == (ft = (FitsFile_Type *) SLang_object_from_mmt (*mmt)))
     {
	SLang_free_mmt (*mmt);
	*mmt = NULL;
     }
   return ft;
}

static int create_img (FitsFile_Type *ft, int *bitpix,
		       SLang_Array_Type *at_naxes)
{
   long *axes;
   unsigned int i, imax;
   int status = 0;

   if (ft->fptr == NULL)
     return -1;

   if (at_naxes->data_type != SLANG_INT_TYPE)
     {
	SLang_verror (SL_TYPE_MISMATCH,
		      "fits_create_img: naxis must be an integer array");
	return -1;
     }

   imax = at_naxes->num_elements;
   axes = (long *) SLmalloc ((imax+1) * sizeof (long));
   if (axes == NULL)
     return -1;

   /* Transpose to FORTRAN order */
   for (i = 0; i < imax; i++)
     axes[i] = ((int *) at_naxes->data)[imax-(i+1)];

   (void) fits_create_img (ft->fptr, *bitpix, imax, axes, &status);
   SLfree ((char *) axes);
   return status;
}

static int write_img (FitsFile_Type *ft, SLang_Array_Type *at)
{
   int type;
   int status = 0;

   if (ft->fptr == NULL)
     return -1;

   switch (at->data_type)
     {
      case SLANG_STRING_TYPE:
	type = TSTRING;
	break;

      case SLANG_DOUBLE_TYPE:
	type = TDOUBLE;
	break;

      case SLANG_FLOAT_TYPE:
	type = TFLOAT;
	break;

      case SLANG_INT16_TYPE:
	type = CFITSIO_INT16_TYPE;
	break;

      case SLANG_UINT16_TYPE:
	type = CFITSIO_UINT16_TYPE;
	break;

      case SLANG_INT32_TYPE:
	type = CFITSIO_INT32_TYPE;
	break;

      case SLANG_UINT32_TYPE:
	type = CFITSIO_UINT32_TYPE;
	break;

      case SLANG_CHAR_TYPE:
      case SLANG_UCHAR_TYPE:
	type = TBYTE;
	break;

#ifdef SLANG_INT64_TYPE
      case SLANG_INT64_TYPE:
	type = CFITSIO_INT64_TYPE;
	break;
#endif

      default:
	SLang_verror (SL_NOT_IMPLEMENTED,
		      "fits_write_img: %s not supported",
		      SLclass_get_datatype_name (at->data_type));
	return -1;
     }

   return fits_write_img (ft->fptr, type, 1, at->num_elements,
			  at->data, &status);
}

static int read_img (FitsFile_Type *ft, SLang_Ref_Type *ref)
{
   int status = 0;
   int anynul = 0;
   int type, stype;
   int num_dims, i;
   long ldims[SLARRAY_MAX_DIMS];
   int dims[SLARRAY_MAX_DIMS];
   SLang_Array_Type *at;

   if (ft->fptr == NULL)
     return -1;

#ifdef fits_get_img_equivtype
   status = fits_get_img_equivtype (ft->fptr, &type, &status);
#else
   status = fits_get_img_type (ft->fptr, &type, &status);
#endif
   if (status)
     return status;

   switch (type)
     {
      case BYTE_IMG:
	stype = SLANG_UCHAR_TYPE;
	type = TBYTE;
	break;

      case SHORT_IMG:		       /* 16 bit image */
	stype = SLANG_INT16_TYPE;
	type = CFITSIO_INT16_TYPE;
	break;

      case USHORT_IMG:		       /* 16 bit image */
	stype = SLANG_UINT16_TYPE;
	type = CFITSIO_UINT16_TYPE;
	break;

      case LONG_IMG:		       /* 32 bit image */
	stype = SLANG_INT32_TYPE;
	type = CFITSIO_INT32_TYPE;
	break;

      case ULONG_IMG:		       /* 32 bit image */
	stype = SLANG_UINT32_TYPE;
	type = CFITSIO_UINT32_TYPE;
	break;

#ifdef TLONGLONG
      case LONGLONG_IMG:	       /* 64 bit image */
	stype = SLANG_INT64_TYPE;
	type = CFITSIO_INT64_TYPE;
	break;
#endif

      case DOUBLE_IMG:
	stype = SLANG_DOUBLE_TYPE;
	type = TDOUBLE;
	break;

      case FLOAT_IMG:
      default:
	stype = SLANG_FLOAT_TYPE;
	type = TFLOAT;
	break;
     }

   if (fits_get_img_dim (ft->fptr, &num_dims, &status))
     return status;

   if ((num_dims > SLARRAY_MAX_DIMS) || (num_dims < 0))
     {
	SLang_verror (SL_NOT_IMPLEMENTED, "Image dimensionality is not supported");
	return -1;
     }

   if (fits_get_img_size (ft->fptr, num_dims, ldims, &status))
     return status;

#if 0
   for (i = 0; i < num_dims; i++) dims[i] = (int) ldims[i];
#else
   for (i = 0; i < num_dims; i++) dims[num_dims-1-i] = (int) ldims[i];
#endif

   if (NULL == (at = SLang_create_array (stype, 0, NULL, dims, num_dims)))
     return -1;

   status = fits_read_img (ft->fptr, type, 1, at->num_elements, NULL,
			   at->data, &anynul, &status);

   if (status)
     {
	SLang_free_array (at);
	return status;
     }

   if (-1 == SLang_assign_to_ref (ref, SLANG_ARRAY_TYPE, (VOID_STAR)&at))
     status = -1;

   SLang_free_array (at);
   return status;
}

static int create_binary_tbl (void)
{
   SLang_MMT_Type *mmt;
   FitsFile_Type *ft;
   SLang_Array_Type *at_ttype, *at_tform, *at_tunit;
   char *extname;
   int tfields, nrows;
   int status;

   status = -1;
   at_ttype = at_tform = at_tunit = NULL;
   mmt = NULL;
   ft = NULL;

   if (-1 == pop_string_or_null (&extname))
     return -1;

   if (-1 == pop_array_or_null (&at_tunit))
     goto free_and_return;

   if (-1 == SLang_pop_array (&at_tform, 1))
     goto free_and_return;

   if (-1 == SLang_pop_array (&at_ttype, 1))
     goto free_and_return;

   if (-1 == SLang_pop_integer (&nrows))
     goto free_and_return;

   if (NULL == (ft = pop_fits_type (&mmt)))
     goto free_and_return;

   if (ft->fptr == NULL)
     goto free_and_return;

   tfields = (int) at_ttype->num_elements;

   if (at_ttype->data_type != SLANG_STRING_TYPE)
     {
	SLang_verror (SL_TYPE_MISMATCH,
		      "fits_create_binary_tbl: ttype must be String_Type[%d]",
		      tfields);
	goto free_and_return;
     }

   if ((tfields != (int) at_tform->num_elements)
       || (at_tform->data_type != SLANG_STRING_TYPE))
     {
	SLang_verror (SL_TYPE_MISMATCH,
		      "fits_create_binary_tbl: tform must be String_Type[%d]",
		      tfields);
	goto free_and_return;
     }

   if ((at_tunit != NULL)
       && ((tfields != (int) at_tunit->num_elements)
	   || (at_tunit->data_type != SLANG_STRING_TYPE)))
     {
	SLang_verror (SL_TYPE_MISMATCH,
		      "fits_create_binary_tbl: tunit must be String_Type[%d]",
		      tfields);
	goto free_and_return;
     }

   status = 0;
   fits_create_tbl (ft->fptr, BINARY_TBL, nrows, tfields,
		    (char **) at_ttype->data,
		    (char **) at_tform->data,
		    ((at_tunit == NULL) ? NULL : (char **) at_tunit->data),
		    extname, &status);

   /* drop */

   free_and_return:
   SLang_free_array (at_ttype);
   SLang_free_array (at_tform);
   SLang_free_array (at_tunit);
   SLang_free_mmt (mmt);
   SLang_free_slstring (extname);

   return status;
}

static int update_key (void)
{
   SLang_MMT_Type *mmt;
   FitsFile_Type *ft;
   char *comment;
   int i;
   unsigned int ui;
   double d;
   long l;
   unsigned long ul;
   char *s;
   char *key;
   int type;
   VOID_STAR v;
   int status;

   if (-1 == pop_string_or_null (&comment))
     return -1;

   key = s = NULL;
   mmt = NULL;
   status = -1;

   type = SLang_peek_at_stack ();
   switch (type)
     {
      case SLANG_STRING_TYPE:
	type = TSTRING;
	if (-1 == SLang_pop_slstring (&s))
	  goto free_and_return;
	v = (VOID_STAR) s;
	break;

      case SLANG_CHAR_TYPE:
      case SLANG_SHORT_TYPE:
      case SLANG_INT_TYPE:
	type = TINT;
	if (-1 == SLang_pop_integer (&i))
	  goto free_and_return;
	v = (VOID_STAR) &i;
	break;

      case SLANG_UCHAR_TYPE:
      case SLANG_USHORT_TYPE:
      case SLANG_UINT_TYPE:
	type = TUINT;
	if (-1 == SLang_pop_uint (&ui))
	  goto free_and_return;
	v = (VOID_STAR) &ui;
	break;

      case SLANG_LONG_TYPE:
	type = TLONG;
	if (-1 == SLang_pop_long (&l))
	  goto free_and_return;
	v = (VOID_STAR) &l;
	break;

      case SLANG_ULONG_TYPE:
	type = TULONG;
	if (-1 == SLang_pop_long (&l))
	  goto free_and_return;
	v = (VOID_STAR) &ul;
	break;

      case SLANG_NULL_TYPE:
	if (-1 == SLang_pop_null ())
	  goto free_and_return;
	v = NULL;
	break;

      case -1:			       /* stack underflow */
	goto free_and_return;

      case SLANG_DOUBLE_TYPE:
      default:
	type = TDOUBLE;
#if SLANG_VERSION < 20000
	if (-1 == SLang_pop_double (&d, NULL, NULL))
	  goto free_and_return;
#else
	if (-1 == SLang_pop_double (&d))
	  goto free_and_return;
#endif
	v = (VOID_STAR) &d;
	break;
     }

   if (-1 == SLang_pop_slstring (&key))
     goto free_and_return;

   if (NULL == (ft = pop_fits_type (&mmt)))
     goto free_and_return;

   if (ft->fptr == NULL)
     goto free_and_return;

   status = 0;
   if (v != NULL)
     {
	if (type == TSTRING)
	  fits_update_key_longstr (ft->fptr, key, (char *)v, comment, &status);
	else
	  fits_update_key (ft->fptr, type, key, v, comment, &status);
     }
   else
     fits_update_key_null (ft->fptr, key, comment, &status);

   free_and_return:

   SLang_free_mmt (mmt);
   SLang_free_slstring (key);
   SLang_free_slstring (comment);
   SLang_free_slstring (s);

   return status;
}

static int update_logical (void)
{
   SLang_MMT_Type *mmt;
   FitsFile_Type *ft;
   char *comment;
   int i;
   char *key;
   int status;

   if (-1 == pop_string_or_null (&comment))
     return -1;

   key = NULL;
   mmt = NULL;
   status = -1;

   if ((0 == SLang_pop_integer (&i))
       && (0 == SLang_pop_slstring (&key))
       && (NULL != (ft = pop_fits_type (&mmt)))
       && (ft->fptr != NULL))
     {
	status = 0;
	fits_update_key (ft->fptr, TLOGICAL, key,
			 (VOID_STAR) &i, comment, &status);
     }

   SLang_free_mmt (mmt);
   SLang_free_slstring (key);
   SLang_free_slstring (comment);

   return status;
}

static int write_comment (FitsFile_Type *ft, char *comment)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_write_comment (ft->fptr, comment, &status);
}

static int write_history (FitsFile_Type *ft, char *comment)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_write_history (ft->fptr, comment, &status);
}

static int write_date (FitsFile_Type *ft)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_write_date (ft->fptr, &status);
}

static int write_record (FitsFile_Type *ft, char *card)
{
   int status = 0;

   if (ft->fptr == NULL)
     return -1;

   /* How robust is fits_write_record to cards that are not 80 characters long? */
   return fits_write_record (ft->fptr, card, &status);
}

static int insert_record (FitsFile_Type *ft, int *keynum, char *card)
{
   int status = 0;

   if (ft->fptr == NULL)
     return -1;

   return fits_insert_record (ft->fptr, *keynum, card, &status);
}

static int modify_name (FitsFile_Type *ft, char *oldname, char *newname)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_modify_name (ft->fptr, oldname, newname, &status);
}

static int do_get_keytype (fitsfile *f, char *name, int *stype)
{
   int status = 0;
   char type;
   int s;
   char card [FLEN_CARD + 1];
   char value [FLEN_CARD + 1];

   if (f == NULL)
     return -1;

   if (0 != fits_read_card (f, name, card, &status))
     return status;

   if (0 != fits_parse_value (card, value, NULL, &status))
     return status;

   if (0 != fits_get_keytype (value, &type, &status))
     return status;

   switch (type)
     {
      case 'C':
	s = SLANG_STRING_TYPE;
	break;

      case 'L':
	s = SLANG_INT_TYPE;
	break;

      case 'F':
	s = SLANG_DOUBLE_TYPE;
	break;

      case 'X':
	s = SLANG_COMPLEX_TYPE;
	break;

      case 'I':
      default:
	s = SLANG_INT_TYPE;
	break;
     }
   *stype = s;

   return 0;
}

static int read_key (int type)
{
   SLang_MMT_Type *mmt;
   FitsFile_Type *ft;
   int status;
   char *name;
   char comment_buf [FLEN_COMMENT];
   SLang_Ref_Type *comment_ref, *v_ref;
   int ival;
   long lval;
   double dval;
   char *sval = NULL;
   int ftype;
   VOID_STAR v;

   v_ref = comment_ref = NULL;
   name = NULL;
   mmt = NULL;
   status = -1;

   if (SLANG_NULL_TYPE == SLang_peek_at_stack ())
     {
	if (-1 == SLang_pop_null ())
	  return -1;
     }
   else if (-1 == SLang_pop_ref (&comment_ref))
     return -1;

   if (-1 == SLang_pop_ref (&v_ref))
     goto free_and_return;

   if (-1 == SLang_pop_slstring (&name))
     goto free_and_return;

   if (NULL == (ft = pop_fits_type (&mmt)))
     goto free_and_return;

   if (ft->fptr == NULL)
     goto free_and_return;

   if (type == SLANG_VOID_TYPE)
     {
	if (0 != (status = do_get_keytype (ft->fptr, name, &type)))
	  goto free_and_return;
	/* if (type == SLANG_INT_TYPE) type = SLANG_LONG_TYPE; */

	status = -1;
     }

   switch (type)
     {
      case SLANG_INT_TYPE:
	v = (VOID_STAR) &ival;
	ftype = TINT;
	ival = 0;
	break;

      case SLANG_LONG_TYPE:
	v = (VOID_STAR) &lval;
	ftype = TLONG;
	ival = 0;
	break;

      case SLANG_DOUBLE_TYPE:
	ftype = TDOUBLE;
	v = (VOID_STAR) &dval;
	dval = 0.0;
	break;

      case SLANG_STRING_TYPE:
	ftype = TSTRING;
	v = (VOID_STAR) &sval;
	break;

      default:
	SLang_verror (SL_INVALID_PARM,
		      "fits_read_key: type %s not supported",
		      SLclass_get_datatype_name (type));
	goto free_and_return;
     }

   status = 0;
   if (ftype == TSTRING)
     fits_read_key_longstr (ft->fptr, name, &sval, comment_buf, &status);
   else
     fits_read_key (ft->fptr, ftype, name, v, comment_buf, &status);

   if (status == 0)
     {
	if (comment_ref != NULL)
	  {
	     char *cptr = comment_buf;
	     if (-1 == SLang_assign_to_ref (comment_ref,
					    SLANG_STRING_TYPE, (VOID_STAR) &cptr))
	       {
		  status = -1;
		  goto free_and_return;
	       }
	  }
	if (-1 == SLang_assign_to_ref (v_ref, type, v))
	  {
	     status = -1;
	     goto free_and_return;
	  }
     }

   free_and_return:
   SLfree (sval);
   SLang_free_ref (comment_ref);
   SLang_free_ref (v_ref);
   SLang_free_slstring (name);
   SLang_free_mmt (mmt);
   return status;
}

static int read_key_integer (void)
{
   return read_key (SLANG_INT_TYPE);
}

static int read_key_double (void)
{
   return read_key (SLANG_DOUBLE_TYPE);
}

static int read_key_string (void)
{
   return read_key (SLANG_STRING_TYPE);
}

static int read_generic_key (void)
{
   return read_key (SLANG_VOID_TYPE);
}

static int read_record (FitsFile_Type *ft, int *keynum, SLang_Ref_Type *ref)
{
   int status = 0;
   char card[FLEN_CARD+1];

   if (ft->fptr == NULL)
     return -1;

   if (0 == fits_read_record (ft->fptr, *keynum, card, &status))
     {
	char *c = card;
	if (-1 == SLang_assign_to_ref (ref, SLANG_STRING_TYPE, &c))
	  return -1;
     }

   return status;
}

static int delete_key (FitsFile_Type *ft, char *key)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_delete_key (ft->fptr, key, &status);
}

static int get_colnum_internal (FitsFile_Type *ft, char *name, SLang_Ref_Type *ref, int casesen)
{
   int status = 0;
   int col;

   if (ft->fptr == NULL)
     return -1;
   /* FIXME: fits_get_colnum may be used to get columns matching a pattern */
   col = 1;
   fits_get_colnum (ft->fptr, casesen, name, &col, &status);

   if (-1 == SLang_assign_to_ref (ref, SLANG_INT_TYPE, (VOID_STAR) &col))
     status = -1;

   return status;
}

static int get_colnum (FitsFile_Type *ft, char *name, SLang_Ref_Type *ref)
{
   return get_colnum_internal (ft, name, ref, CASEINSEN);
}

static int get_colnum_casesen (FitsFile_Type *ft, char *name, SLang_Ref_Type *ref)
{
   return get_colnum_internal (ft, name, ref, CASESEN);
}


static int insert_rows (FitsFile_Type *ft, int *first, int *num)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   if ((*first < 0) || (*num < 0))
     {
	SLang_verror (SL_INVALID_PARM, "fits_insert_rows: first and num must be non-negative");
	return -1;
     }

   return fits_insert_rows (ft->fptr, *first, *num, &status);
}

static int delete_rows (FitsFile_Type *ft, int *first, int *num)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   if ((*first <= 0) || (*num < 0))
     {
	SLang_verror (SL_INVALID_PARM, "fits_delete_rows: first and num must be positive");
	return -1;
     }

   return fits_delete_rows (ft->fptr, *first, *num, &status);
}

static int insert_cols (FitsFile_Type *ft, int *colnum,
			SLang_Array_Type *at_ttype,
			SLang_Array_Type *at_tform)
{
   int ncols;
   char **ttype, **tform;
   int i;
   int status = 0;

   if (ft->fptr == NULL)
     return -1;
   ncols = at_ttype->num_elements;
   if ((ncols < 0) || (ncols != (int) at_tform->num_elements)
       || (at_ttype->data_type != SLANG_STRING_TYPE)
       || (at_tform->data_type != SLANG_STRING_TYPE))
     {
	SLang_verror (SL_INVALID_PARM,
		      "fits_insert_cols: ttype and tform must be string arrays of same size");
	return -1;
     }

   if (*colnum <= 0)
     {
	SLang_verror (SL_INVALID_PARM, "fits_insert_cols: colnum must be positive");
	return -1;
     }

   tform = (char **)at_tform->data;
   ttype = (char **)at_ttype->data;

   for (i = 0; i < ncols; i++)
     {
	if ((tform[i] == NULL) || (ttype[i] == NULL))
	  {
	     SLang_verror (SL_INVALID_PARM,
			   "fits_insert_cols: ttype and tform elements muts be non NULL");
	     return -1;
	  }
     }
   return fits_insert_cols (ft->fptr, *colnum, ncols, ttype, tform, &status);
}

static int delete_col (FitsFile_Type *ft, int *col)
{
   int status = 0;
   if (ft->fptr == NULL)
     return -1;
   return fits_delete_col (ft->fptr, *col, &status);
}

static int get_num_rows (FitsFile_Type *ft, SLang_Ref_Type *ref)
{
   long nrows;
   int status = 0;

   if (ft->fptr == NULL)
     return -1;
   if (0 == fits_get_num_rows (ft->fptr, &nrows, &status))
     {
	int inrows = (int) nrows;
	if (-1 == SLang_assign_to_ref (ref, SLANG_INT_TYPE, (VOID_STAR) &inrows))
	  return -1;
     }
   return status;
}

static int get_rowsize (FitsFile_Type *ft, SLang_Ref_Type *ref)
{
   long nrows;
   int status = 0;

   if (ft->fptr == NULL)
     return -1;
   if (0 == fits_get_rowsize (ft->fptr, &nrows, &status))
     {
	int inrows = (int) nrows;
	if (-1 == SLang_assign_to_ref (ref, SLANG_INT_TYPE, (VOID_STAR) &inrows))
	  return -1;
     }
   return status;
}

static int get_num_cols (FitsFile_Type *ft, SLang_Ref_Type *ref)
{
   int ncols;
   int status = 0;

   if (ft->fptr == NULL)
     return -1;
   if (0 == fits_get_num_cols (ft->fptr, &ncols, &status))
     {
	if (-1 == SLang_assign_to_ref (ref, SLANG_INT_TYPE, (VOID_STAR) &ncols))
	  return -1;
     }
   return status;
}

static void byte_swap32 (unsigned char *ss, unsigned int n)
{
   unsigned char *p, *pmax, ch;

   p = (unsigned char *) ss;
   pmax = p + 4 * n;
   while (p < pmax)
     {
	ch = *p;
	*p = *(p + 3);
	*(p + 3) = ch;

	ch = *(p + 1);
	*(p + 1) = *(p + 2);
	*(p + 2) = ch;
	p += 4;
     }
}

static void byte_swap16 (unsigned char *p, unsigned int nread)
{
   unsigned char *pmax, ch;

   pmax = p + 2 * nread;
   while (p < pmax)
     {
	ch = *p;
	*p = *(p + 1);
	*(p + 1) = ch;
	p += 2;
     }
}

/* MAJOR HACK!!!! */
static int hack_write_bit_col (fitsfile *f, unsigned int col,
			       unsigned int row, unsigned int firstelem,
			       unsigned int sizeof_type, unsigned int num_elements,
			       unsigned char *bytes)
{
   static int status = 0;
   tcolumn *colptr;
   long trepeat;
   int tcode;

   if ((f == NULL) || (f->Fptr == NULL)
       || (NULL == (colptr  = (f->Fptr)->tableptr)))
     return WRITE_ERROR;

   colptr += (col - 1);     /* offset to correct column structure */
   trepeat = colptr->trepeat;
   tcode = colptr->tdatatype;

   colptr->tdatatype = TBYTE;
   colptr->trepeat = sizeof_type;

   (void) fits_write_col (f, TBYTE, col, row, firstelem,
			  num_elements*sizeof_type, bytes, &status);

   colptr->tdatatype = tcode;
   colptr->trepeat = trepeat;

   return status;
}

static int write_tbit_col (fitsfile *f, unsigned int col, unsigned int row,
			   unsigned int firstelem, unsigned int repeat,
			   unsigned int width, SLang_Array_Type *at)
{
   int status = 0;
   unsigned int num_elements;
   unsigned int sizeof_type;
   unsigned char *data, *buf;
   unsigned short s;
   void (*bs) (unsigned char *, unsigned int);

   (void) width;
   num_elements = at->num_elements;
   sizeof_type = at->sizeof_type;
   data = (unsigned char *) at->data;

   if (8 * sizeof_type != repeat)
     {
	SLang_verror (SL_NOT_IMPLEMENTED, "Writing a %dX bit column requires the appropriately sized integer",
		      repeat);
	return -1;
     }

   s = 0x1234;
   if ((*(unsigned char *) &s == 0x12)
       || (sizeof_type == 1))
     {
	return hack_write_bit_col (f, col, row, firstelem,
				   sizeof_type, num_elements, data);
     }

   /* Sigh.  Need to byteswap */
   switch (sizeof_type)
     {
      case 2:
	bs = byte_swap16;
	break;

      case 4:
	bs = byte_swap32;
	break;

      default:
	SLang_verror (SL_NOT_IMPLEMENTED, "writing to a %dX column is not supported", repeat);
	return -1;
     }

   buf = (unsigned char *) SLmalloc (num_elements * sizeof_type);
   if (buf == NULL)
     return -1;

   memcpy (buf, data, num_elements * sizeof_type);
   (*bs) (buf, num_elements);

   status = hack_write_bit_col (f, col, row, firstelem,
				sizeof_type, num_elements, buf);

   SLfree ((char *) buf);
   return status;
}

#ifdef fits_get_eqcoltype
# define GET_COL_TYPE fits_get_eqcoltype
#else
# define GET_COL_TYPE my_fits_get_coltype
static int my_fits_get_coltype (fitsfile *fptr, int col, int *type,
				long *repeat, long *width, int *statusp)
{
   int status = *statusp;
   char tscaln[32];
   char tzeron[32];
   double tscal, tzero;
   double min_val, max_val;

   if (ft->fptr == NULL)
     return -1;

   if (0 != fits_get_coltype (fptr, col, type, repeat, width, &status))
     {
	*statusp = status;
	return status;
     }

   sprintf (tscaln, "TSCAL%d", col);
   sprintf (tzeron, "TZERO%d", col);

   if (0 != fits_read_key (fptr, TDOUBLE, tscaln, &tscal, NULL, &status))
     {
	fits_clear_errmsg ();
	tscal = 1.0;
     }

   if (0 != fits_read_key (fptr, TDOUBLE, tzeron, &tzero, NULL, &status))
     {
	fits_clear_errmsg ();
	tzero = 0;
     }

   switch (*type)
     {
      case TSHORT:
	min_val = -32768.0;
	max_val = 32767.0;
	break;

      case TLONG:
	min_val = -2147483648.0;
	max_val = 2147483647.0;
	break;

      default:
	return 0;
     }

   min_val = tzero + tscal * min_val;
   max_val = tzero + tscal * max_val;

   if (min_val > max_val)
     {
	double tmp = max_val;
	max_val = min_val;
	min_val = tmp;
     }

   if ((min_val >= -32768.0) && (max_val <= 32767.0))
     {
	*type = TSHORT;
	return 0;
     }

   if ((min_val >= 0) && (max_val <= 65535.0))
     {
	*type = TUSHORT;
	return 0;
     }

   if ((min_val >= -2147483648.0) && (max_val <= 2147483647.0))
     {
	*type = TLONG;
	return 0;
     }

   if ((min_val >= 0.0) && (max_val <= 4294967295.0))
     {
	*type = TULONG;
	return 0;
     }

   *type = TDOUBLE;
   return 0;
}
#endif

static int write_col (FitsFile_Type *ft, int *colnum,
		      int *firstrow, int *firstelem, SLang_Array_Type *at)
{
   int type;
   int status = 0;
   int col;
   long repeat;
   long width;

   if (ft->fptr == NULL)
     return -1;

   col = *colnum;

   if (0 != GET_COL_TYPE (ft->fptr, col, &type, &repeat, &width, &status))
     return status;

   if (type == TBIT)
     return write_tbit_col (ft->fptr, col, *firstrow, *firstelem,
			    repeat, width, at);

   switch (at->data_type)
     {
      case SLANG_INT16_TYPE:
	type = CFITSIO_INT16_TYPE;
	break;
      case SLANG_UINT16_TYPE:
	type = CFITSIO_UINT16_TYPE;
	break;
      case SLANG_INT32_TYPE:
	type = CFITSIO_INT32_TYPE;
	break;
      case SLANG_UINT32_TYPE:
	type = CFITSIO_UINT32_TYPE;
	break;
#ifdef CFITSIO_INT64_TYPE
      case SLANG_INT64_TYPE:
	type = CFITSIO_INT64_TYPE;
	break;
#endif

      case SLANG_DOUBLE_TYPE:
	type = TDOUBLE;
	break;

      case SLANG_FLOAT_TYPE:
	type = TFLOAT;
	break;

      case SLANG_STRING_TYPE:
	type = TSTRING;
	break;

      case SLANG_CHAR_TYPE:
      case SLANG_UCHAR_TYPE:
	type = TBYTE;
	break;

      default:
	SLang_verror (SL_NOT_IMPLEMENTED,
		      "fits_write_col: %s not suppported",
		      SLclass_get_datatype_name (at->data_type));
	return -1;
     }

   return fits_write_col (ft->fptr, type, *colnum, *firstrow, *firstelem,
			  at->num_elements, at->data, &status);
}

static int read_string_cell (fitsfile *f, unsigned int row, unsigned int col,
			     unsigned int len, unsigned int num_substrs, char **sp)
{
   char *s, *ss;
   char *sls;
   int status = 0;
   int anynul;
   unsigned int i;
   unsigned int width;

   *sp = NULL;

   if (f == NULL)
     return -1;

   if (NULL == (s = SLmalloc (len + 1)))
     return -1;
   memset (s, ' ', len);

   width = len/num_substrs;
   ss = s;
   for (i = 0; i < num_substrs; i++)
     {
	/* Note that the number of elements passed to this must be 1 */
	if (0 != fits_read_col (f, TSTRING, col, row, i+1, 1, NULL,
				&ss, &anynul, &status))
	  {
	     SLfree (s);
	     return status;
	  }
	/* If there is more than one substring, append them together.  Only the
	 * last substring will have trailing whitespace removed.  This is
	 * probably ok because fits does not like trailing whitespace.
	 */
	if (i + 1 < num_substrs)
	  {
	     unsigned int len1 = strlen (ss);
	     ss[len1] = ' ';
	  }
	ss += width;
     }
   sls = SLang_create_slstring (s);
   SLfree (s);
   if (sls == NULL)
     return -1;

   *sp = sls;
   return 0;
}

static int read_string_column (fitsfile *f, int is_var, long repeat, unsigned int num_substrs,
			       int col, unsigned int firstrow, unsigned int numrows,
			       SLang_Array_Type **atp)
{
   int num_elements;
   char **ats;
   unsigned int i;
   int status = 0;
   SLang_Array_Type *at;

   *atp = NULL;

   if (f == NULL)
     return -1;

   num_elements = (int) numrows;
   at = SLang_create_array (SLANG_STRING_TYPE, 0, NULL, &num_elements, 1);
   if (at == NULL)
     return -1;

   ats = (char **) at->data;

   for (i = 0; i < numrows; i++)
     {
	long offset;
	long row;

	row = firstrow + i;
	if (is_var)
	  {
	     if (0 != fits_read_descript (f, col, row, &repeat, &offset, &status))
	       {
		  SLang_free_array (at);
		  return status;
	       }
	  }

	status = read_string_cell (f, row, col, repeat, num_substrs, ats+i);
	if (status != 0)
	  {
	     SLang_free_array (at);
	     return status;
	  }
     }

   *atp = at;
   return 0;
}

static int read_bit_column (fitsfile *f, unsigned int col, unsigned int row,
			    unsigned int firstelem, unsigned int num_elements,
			    unsigned char *data, unsigned int bytes_per_elem,
			    unsigned int bits_per_elem)
{
   int status, anynul;
   unsigned short s;

   if (f == NULL)
     return -1;

   status = 0;
   if (0 != fits_read_col (f, TBYTE, col, row,
			   firstelem, num_elements*bytes_per_elem,
			   NULL, data, &anynul, &status))
     return status;

   s = 0x1234;
   if (*(unsigned char *) &s == 0x12)
     return status;

   /* Otherwise, byteswap */
   switch (bytes_per_elem)
     {
	SLuindex_Type i;
	int shift;
	int16_type *data16;
	int32_type *data32;

      case 1:
	shift = 8*bytes_per_elem - bits_per_elem;
	if (shift) for (i = 0; i < num_elements; i++)
	  {
	     data[i] = (data[i] >> shift);
	  }
	break;

      case 2:
	byte_swap16 (data, num_elements);
	shift = 8*bytes_per_elem - bits_per_elem;
	data16 = (int16_type *)data;
	if (shift) for (i = 0; i < num_elements; i++)
	  {
	     data16[i] = (data16[i] >> shift);
	  }
	break;

      case 4:
	byte_swap32 (data, num_elements);
	shift = 8*bytes_per_elem - bits_per_elem;
	data32 = (int32_type *)data;
	if (shift) for (i = 0; i < num_elements; i++)
	  {
	     data32[i] = (data32[i] >> shift);
	  }
	break;

      default:
	SLang_verror (SL_NOT_IMPLEMENTED, "%u byte integers are unsupported",
		      bytes_per_elem);
	return -1;
     }

   return 0;
}

static int read_column_values (fitsfile *f, int type, SLtype datatype,
			       unsigned int row, unsigned int col, unsigned int num_rows,
			       int repeat, int repeat_orig, SLang_Array_Type **atp)
{
   int num_elements;
   int status = 0;
   int anynul;
   SLang_Array_Type *at;
   int dims[2];
   int num_dims;

   *atp = NULL;

   if (f == NULL)
     return -1;

   num_elements = num_rows * repeat;
   if (num_rows <= 1)
     {
	dims[0] = num_elements;
	num_dims = 1;
     }
   else				       /* was repeat>1 */
     {
	dims[0] = num_rows;
	dims[1] = repeat;
	num_dims = 2;
     }

   at = SLang_create_array (datatype, 0, NULL, dims, num_dims);
   if (at == NULL)
     return -1;

   if (num_elements)
     {
	if (type == TBIT)
	  status = read_bit_column (f, col, row, 1, num_elements, (unsigned char *)at->data, at->sizeof_type, repeat_orig);
	else
	  (void) fits_read_col (f, type, col, row, 1, num_elements, NULL,
				at->data, &anynul, &status);
     }

   if (status)
     {
	SLang_free_array (at);
	return status;
     }

   *atp = at;
   return 0;
}

static int read_var_column (fitsfile *f, int ftype, SLtype datatype,
			    int col, unsigned int firstrow, unsigned int num_rows,
			    SLang_Array_Type **atp)
{
   int num_elements;
   unsigned int i;
   SLang_Array_Type **ati;
   SLang_Array_Type *at;

   *atp = NULL;
   if (f == NULL)
     return -1;

   num_elements = (int) num_rows;

   at = SLang_create_array (SLANG_ARRAY_TYPE, 0, NULL, &num_elements, 1);
   if (at == NULL)
     return -1;

   ati = (SLang_Array_Type **) at->data;
   for (i = 0; i < num_rows; i++)
     {
	long offset;
	long repeat;
	unsigned int row;
	int status = 0;

	row = firstrow + i;
	if (0 != fits_read_descript (f, col, row, &repeat, &offset, &status))
	  {
	     SLang_free_array (at);
	     return status;
	  }

	status = read_column_values (f, ftype, datatype, row, col, 1, repeat, repeat, ati+i);
	if (status)
	  {
	     SLang_free_array (at);
	     return status;
	  }
     }

   *atp = at;
   return 0;
}

static int read_col (FitsFile_Type *ft, int *colnum, int *firstrowp,
		     int *num_rowsp, SLang_Ref_Type *ref)
{
   SLang_Array_Type *at;
   int type;
   long num_rows;
   long width;
   int status;
   SLtype datatype;
   int num_columns;
   int firstrow;
   long repeat, save_repeat;
   int col;

   if (ft->fptr == NULL)
     return -1;

   status = 0;
   if (0 != fits_get_num_cols (ft->fptr, &num_columns, &status))
     return status;

   if (0 != fits_get_num_rows (ft->fptr, &num_rows, &status))
     return status;

   if (*num_rowsp <= 0)
     {
	SLang_verror (SL_INVALID_PARM, "Number of rows must positive");
	return -1;
     }

   col = *colnum;

   if ((col <= 0) || (col > num_columns))
     {
	SLang_verror (SL_INVALID_PARM, "Column number out of range");
	return -1;
     }
   firstrow = *firstrowp;
   if ((firstrow <= 0) || (firstrow > num_rows))
     {
	SLang_verror (SL_INVALID_PARM, "Row number out of range");
	return -1;
     }

   if (firstrow + *num_rowsp > num_rows + 1)
     num_rows = num_rows - (firstrow - 1);
   else
     num_rows = *num_rowsp;

   if (0 != GET_COL_TYPE (ft->fptr, col, &type, &repeat, &width, &status))
     return status;

   save_repeat = repeat;
   if (-1 == map_fitsio_type_to_slang (&type, &repeat, &datatype))
     return -1;

   if (datatype == SLANG_STRING_TYPE)
     {
	unsigned int num_substrs;
	/* This assumes an ASCII_TBL, which will always have a
	 * repeat of 1, and the number of bytes is given by the
	 * width field.  In contrast, a BINARY_TBL will have
	 * repeat = number of bytes, and width represents the
	 * number of bytes in a substring.
	 */
	if ((repeat == 1) && (width != 1))
	  {
	     num_substrs = 1;
	     repeat = width;
	  }
	else
	  {
	     if (width > 0)
	       num_substrs = repeat/width;
	     else
	       num_substrs = 0;
	  }
	status = read_string_column (ft->fptr, (type < 0), repeat, num_substrs, col, firstrow, num_rows, &at);
     }
   else if (type < 0)
     status = read_var_column (ft->fptr, -type, datatype, col, firstrow, num_rows, &at);
   else
     status = read_column_values (ft->fptr, type, datatype, firstrow, col, num_rows, repeat, save_repeat, &at);

   if (status)
     return status;

   if (-1 == SLang_assign_to_ref (ref, SLANG_ARRAY_TYPE, (VOID_STAR)&at))
     status = -1;

   SLang_free_array (at);
   return status;
}

typedef struct
{
   int type;
   long repeat, width;
   long repeat_orig;		       /* used for tbit columns */
   SLtype datatype;
   unsigned int data_offset;
}
Column_Info_Type;

static int read_var_column_data (fitsfile *f, int ftype, SLtype datatype,
				 int col, unsigned int firstrow, unsigned int num_rows,
				 SLang_Array_Type **at_data)
{
   unsigned int i;

   for (i = 0; i < num_rows; i++)
     {
	long offset;
	long repeat;
	unsigned int row;
	int status = 0;

	row = firstrow + i;
	if (0 != fits_read_descript (f, col, row, &repeat, &offset, &status))
	  return status;

	status = read_column_values (f, ftype, datatype, row, col, 1, repeat, repeat, at_data+i);
	if (status)
	  return status;
     }
   return 0;
}

static int read_string_column_data (fitsfile *f, int is_var, long repeat, unsigned int num_substrs, int col,
				    long firstrow, unsigned int num_rows,
				    char **strs)
{
   unsigned int i;
   int status = 0;

   for (i = 0; i < num_rows; i++)
     {
	long offset;
	long row;

	row = firstrow + i;
	if (is_var)
	  {
	     if (0 != fits_read_descript (f, col, row, &repeat, &offset, &status))
	       return status;
	  }

	status = read_string_cell (f, row, col, repeat, num_substrs, strs + i);
	if (status != 0)
	  return status;
     }
   return 0;
}

/* Usage: read_cols (ft, [columns...], firstrow, nrows, &ref) */
/* TODO: Add support for the following calling convention:
 *    read_cols (ft, [columns...], [rows], &ref)
 * Here rows is an integer-valued array that specifies what rows to be
 * read.
 */
static int read_cols (void)
{
   SLang_MMT_Type *mmt;
   FitsFile_Type *ft;
   fitsfile *f;
   int status;
   int num_columns_in_table;
   long num_rows_in_table, delta_rows;
   int num_rows;
   int firstrow;
   int *cols;
   int num_cols;
   int i;
   SLang_Ref_Type *ref;
   Column_Info_Type *ci = NULL;
   SLang_Array_Type *data_arrays_at = NULL;
   SLang_Array_Type **data_arrays = NULL;
   SLang_Array_Type *columns_at = NULL;

   if (-1 == SLang_pop_ref (&ref))
     return -1;
   if ((-1 == SLang_pop_integer (&num_rows))
       || (-1 == SLang_pop_integer (&firstrow))
       || (-1 == SLang_pop_array (&columns_at, 1)))
     {
	SLang_free_ref (ref);
	return -1;
     }
   if (NULL == (ft = pop_fits_type (&mmt)))
     {
	SLang_free_array (columns_at);
	SLang_free_ref (ref);
	return -1;
     }

   status = -1;
   f = ft->fptr;
   if (f == NULL)
     goto free_and_return_status;

   status = 0;
   if ((0 != fits_get_num_cols (f, &num_columns_in_table, &status))
       || (0 != fits_get_num_rows (f, &num_rows_in_table, &status)))
     goto free_and_return_status;

   if (num_rows < 0)
     {
	SLang_verror (SL_INVALID_PARM, "Number of rows must be non-negative");
	status = -1;
	goto free_and_return_status;
     }

   if ((firstrow <= 0)
       || ((firstrow > num_rows_in_table) && (num_rows > 0)))
     {
	SLang_verror (SL_INVALID_PARM, "Row number out of range");
	return -1;
     }

   if (firstrow + num_rows > num_rows_in_table + 1)
     num_rows = num_rows_in_table - (firstrow - 1);

   cols = (int *)columns_at->data;
   num_cols = columns_at->num_elements;

   if (NULL == (ci = (Column_Info_Type *) SLmalloc (num_cols*sizeof (Column_Info_Type))))
     {
	status = -1;
	goto free_and_return_status;
     }

   data_arrays_at = SLang_create_array (SLANG_ARRAY_TYPE, 0, NULL, &num_cols, 1);
   if (data_arrays_at == NULL)
     {
	status = -1;
	goto free_and_return_status;
     }
   data_arrays = (SLang_Array_Type **)data_arrays_at->data;

   for (i = 0; i < num_cols; i++)
     {
	SLang_Array_Type *at;
	SLtype datatype;
	long repeat;
	int type;
	int col;

	col = cols[i];
	if ((col <= 0) || (col > num_columns_in_table))
	  {
	     SLang_verror (SL_INVALID_PARM, "Column number out of range");
	     status = -1;
	     goto free_and_return_status;
	  }

	if (0 != GET_COL_TYPE (f, col, &type, &repeat, &ci[i].width, &status))
	  goto free_and_return_status;

	ci[i].repeat_orig = repeat;
	if (-1 == map_fitsio_type_to_slang (&type, &repeat, &datatype))
	  {
	     status = -1;
	     goto free_and_return_status;
	  }
	ci[i].repeat = repeat;
	ci[i].type = type;
	ci[i].datatype = datatype;
	ci[i].data_offset = 0;

	if (datatype == SLANG_STRING_TYPE)
	  {
	     at = SLang_create_array (SLANG_STRING_TYPE, 0, NULL, &num_rows, 1);
	  }
	else if (type < 0)	       /* variable length */
	  {
	     if (-type == TBIT)
	       {
		  SLang_verror (SL_NOT_IMPLEMENTED, "Read bit-data from the heap is not supported.  Please report this problem");
		  status = -1;
		  goto free_and_return_status;
	       }
	     at = SLang_create_array (SLANG_ARRAY_TYPE, 0, NULL, &num_rows, 1);
	  }
	else
	  {
	     int dims[2];
	     int num_dims = 1;
	     dims[0] = num_rows;
	     if (repeat > 1)
	       {
		  dims[1] = repeat;
		  num_dims++;
	       }
	     at = SLang_create_array (datatype, 0, NULL, dims, num_dims);
	  }

	if (at == NULL)
	  {
	     status = -1;
	     goto free_and_return_status;
	  }
	data_arrays[i] = at;
     }

   if (fits_get_rowsize (f, &delta_rows, &status))
     goto free_and_return_status;

   if (delta_rows < 1)
     delta_rows = 1;

   while (num_rows)
     {
	if (num_rows < delta_rows)
	  delta_rows = num_rows;

	for (i = 0; i < num_cols; i++)
	  {
	     SLtype datatype = ci[i].datatype;
	     int type = ci[i].type;
	     long repeat = ci[i].repeat;
	     int col = cols[i];
	     SLang_Array_Type *at = data_arrays[i];
	     unsigned int data_offset = ci[i].data_offset;

	     if (datatype == SLANG_STRING_TYPE)
	       {
		  unsigned int num_substrs;
		  /* This assumes an ASCII_TBL, which will always have a
		   * repeat of 1, and the number of bytes is given by the
		   * width field.  In contrast, a BINARY_TBL will have
		   * repeat = number of bytes, and width represents the
		   * number of bytes in a substring.
		   */
		  if ((repeat == 1) && (ci[i].width != 1))
		    {
		       repeat = ci[i].width;
		       num_substrs = 1;
		    }
		  else
		    {
		       if (ci[i].width > 0)
			 num_substrs = repeat / ci[i].width;
		       else
			 num_substrs = 0;
		    }

		  status = read_string_column_data (f, (type < 0), repeat, num_substrs, col, firstrow, delta_rows,
						    (char **)at->data + data_offset);
		  data_offset += delta_rows;
	       }
	     else if (type < 0)
	       {
		  status = read_var_column_data (f, -type, datatype, col, firstrow, delta_rows,
						 (SLang_Array_Type **)at->data + data_offset);
		  data_offset += delta_rows;
	       }
	     else
	       {
		  unsigned int num_elements = repeat * delta_rows;
		  unsigned char *data = (unsigned char *)at->data + data_offset;
		  /* int anynul; */

		  if (type == TBIT)
		    status = read_bit_column (f, col, firstrow, 1, num_elements, data, at->sizeof_type, ci[i].repeat_orig);
		  else
		    (void) fits_read_col (f, type, col, firstrow, 1, num_elements, NULL,
					  data, NULL, &status);

		  data_offset += num_elements * at->sizeof_type;
	       }
	     ci[i].data_offset = data_offset;

	     if (status)
	       goto free_and_return_status;
	  }
	firstrow += delta_rows;
	num_rows -= delta_rows;
     }

   if (status)
     return status;

   if (-1 == SLang_assign_to_ref (ref, SLANG_ARRAY_TYPE, (VOID_STAR)&data_arrays_at))
     status = -1;

   /* drop */

   free_and_return_status:
   SLfree ((char *)ci);
   SLang_free_mmt (mmt);
   SLang_free_array (columns_at);
   SLang_free_ref (ref);
   SLang_free_array (data_arrays_at);

   return status;
}

static void clear_errmsg (void)
{
   fits_clear_errmsg ();
}

static void get_errstatus (int *status)
{
   char errbuf [FLEN_ERRMSG];

   *errbuf = 0;
   fits_get_errstatus (*status, errbuf);
   (void) SLang_push_string (errbuf);
}

static void read_errmsg (void)
{
   char errbuf [FLEN_ERRMSG];

   if (0 == fits_read_errmsg (errbuf))
     (void) SLang_push_null ();
   else
     (void) SLang_push_string (errbuf);
}

static int get_num_keys (FitsFile_Type *f, SLang_Ref_Type *ref)
{
   int status = 0;
   int nkeys;

   if (f->fptr == NULL)
     return -1;
   if (0 == fits_get_hdrspace (f->fptr, &nkeys, NULL, &status))
     return SLang_assign_to_ref (ref, SLANG_INT_TYPE, (VOID_STAR) &nkeys);

   return status;
}

static int get_keytype (FitsFile_Type *f, char *name, SLang_Ref_Type *v)
{
   int status = 0;
   int type;

   if (f->fptr == NULL)
     return -1;
   if (0 == (status = do_get_keytype (f->fptr, name, &type)))
     return SLang_assign_to_ref (v, SLANG_DATATYPE_TYPE, (VOID_STAR) &type);

   return status;
}

#if 0
static int _read_key_n (FitsFile_Type *f, SLang_Ref_Type *v, SLang_Ref_Type *c)
{
   int status = 0;

   if (fits_get_keytype (card, &type, &status))
     return status;
   switch (type)
     {
      case 'C':			       /* string */
	break;

      case 'L':			       /* logical */
	break;

      case 'I':			       /* integer */
	break;

      case 'F':			       /* floating point */
	break;

      case 'X':			       /* complex */
	break;
     }
}
#endif

static void get_version (void)
{
   float v;
   (void) fits_get_version (&v);
   (void) SLang_push_float (v);
}

static int get_keyclass (char *card)
{
   return fits_get_keyclass (card);
}

static int do_fits_fun_f(int (*fun)(fitsfile *, int *), FitsFile_Type *f)
{
   int status = 0;

   if (f->fptr == NULL)
     return -1;

   (void) (*fun) (f->fptr, &status);
   return status;
}

static int write_chksum (FitsFile_Type *f)
{
   return do_fits_fun_f (fits_write_chksum, f);
}
static int update_chksum (FitsFile_Type *f)
{
   return do_fits_fun_f (fits_update_chksum, f);
}

static int verify_chksum (FitsFile_Type *f, SLang_Ref_Type *dataok, SLang_Ref_Type *hduok)
{
   int status = 0;
   int dok=0, hok=0;

   if (f->fptr == NULL)
     return -1;

   if (0 == fits_verify_chksum (f->fptr, &dok, &hok, &status))
     {
	if ((-1 == SLang_assign_to_ref (dataok, SLANG_INT_TYPE, (VOID_STAR)&dok))
	    || (-1 == SLang_assign_to_ref (hduok, SLANG_INT_TYPE, (VOID_STAR)&hok)))
	  status = -1;
     }
   return status;
}

static int get_chksum (FitsFile_Type *f, SLang_Ref_Type *datasum, SLang_Ref_Type *hdusum)
{
   int status = 0;
   unsigned long dsum, hsum;

   if (f->fptr == NULL)
     return -1;

   if (0 == fits_get_chksum (f->fptr, &dsum, &hsum, &status))
     {
	if ((-1 == SLang_assign_to_ref (datasum, SLANG_ULONG_TYPE, (VOID_STAR)&dsum))
	    || (-1 == SLang_assign_to_ref (hdusum, SLANG_ULONG_TYPE, (VOID_STAR)&hsum)))
	  status = -1;
     }
   return status;
}

static int set_bscale (FitsFile_Type *f, double *scale, double *zero)
{
   int status = 0;

   if (f->fptr == NULL)
     return -1;

   return fits_set_bscale (f->fptr, *scale, *zero, &status);
}

static int set_tscale (FitsFile_Type *f, int *colp, double *scale, double *zero)
{
   int status = 0;

   if (f->fptr == NULL)
     return -1;

   return fits_set_tscale (f->fptr, *colp, *scale, *zero, &status);
}

/* DUMMY_FITS_FILE_TYPE is a temporary hack that will be modified to the true
 * id once the interpreter provides it when the class is registered.  See below
 * for details.  The reason for this is simple: for a module, the type-id
 * must be assigned dynamically.
 */
#define DUMMY_FITS_FILE_TYPE	255
#define I SLANG_INT_TYPE
#define S SLANG_STRING_TYPE
#define F DUMMY_FITS_FILE_TYPE
#define R SLANG_REF_TYPE
#define A SLANG_ARRAY_TYPE
#define T SLANG_DATATYPE_TYPE
#define D SLANG_DOUBLE_TYPE

static SLang_Intrin_Fun_Type Fits_Intrinsics [] =
{
   MAKE_INTRINSIC_0("_fits_clear_errmsg", clear_errmsg, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("_fits_read_errmsg", read_errmsg, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_I("_fits_get_errstatus", get_errstatus, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_3("_fits_open_file", open_file, I, R, S, S),
   MAKE_INTRINSIC_1("_fits_delete_file", delete_file, I, F),
   MAKE_INTRINSIC_1("_fits_close_file", close_file, SLANG_INT_TYPE, F),

   /* HDU Access Routines */
   MAKE_INTRINSIC_2("_fits_movabs_hdu", movabs_hdu, I, F, I),
   MAKE_INTRINSIC_2("_fits_movrel_hdu", movrel_hdu, I, F, I),
   MAKE_INTRINSIC_4("_fits_movnam_hdu", movnam_hdu, I, F, I, S, I),
   MAKE_INTRINSIC_2("_fits_get_num_hdus", get_num_hdus, I, F, R),
   MAKE_INTRINSIC_1("_fits_get_hdu_num", get_hdu_num, I, F),
   MAKE_INTRINSIC_2("_fits_get_hdu_type", get_hdu_type, I, F, R),
   MAKE_INTRINSIC_5("_fits_copy_file", copy_file, I, F, F, I, I, I),
   MAKE_INTRINSIC_3("_fits_copy_hdu", copy_hdu, I, F, F, I),
   MAKE_INTRINSIC_2("_fits_copy_header", copy_header, I, F, F),
   MAKE_INTRINSIC_1("_fits_delete_hdu", delete_hdu, I, F),

   MAKE_INTRINSIC_3("_fits_create_img", create_img, I, F, I, A),
   MAKE_INTRINSIC_2("_fits_write_img", write_img, I, F, A),
   MAKE_INTRINSIC_2("_fits_read_img", read_img, I, F, R),

   /* Keword Writing Routines */
   MAKE_INTRINSIC_0("_fits_create_binary_tbl", create_binary_tbl, I),
   MAKE_INTRINSIC_0("_fits_update_key", update_key, I),
   MAKE_INTRINSIC_0("_fits_update_logical", update_logical, I),
   MAKE_INTRINSIC_2("_fits_write_comment", write_comment, I, F, S),
   MAKE_INTRINSIC_2("_fits_write_history", write_history, I, F, S),
   MAKE_INTRINSIC_1("_fits_write_date", write_date, I, F),
   MAKE_INTRINSIC_2("_fits_write_record", &write_record, I, F, S),
   MAKE_INTRINSIC_3("_fits_insert_record", &insert_record, I, F, I, S),

   MAKE_INTRINSIC_3("_fits_modify_name", modify_name, I, F, S, S),
   MAKE_INTRINSIC_2("_fits_get_num_keys", get_num_keys, I, F, R),
   MAKE_INTRINSIC_0("_fits_read_key_integer", read_key_integer, I),
   MAKE_INTRINSIC_0("_fits_read_key_string", read_key_string, I),
   MAKE_INTRINSIC_0("_fits_read_key_double", read_key_double, I),
   MAKE_INTRINSIC_0("_fits_read_key", read_generic_key, I),
   MAKE_INTRINSIC_3("_fits_read_record", read_record, I, F, I, R),

   MAKE_INTRINSIC_2("_fits_delete_key", delete_key, I, F, S),
   MAKE_INTRINSIC_3("_fits_get_colnum", get_colnum, I, F, S, R),
   MAKE_INTRINSIC_3("_fits_get_colnum_casesen", get_colnum_casesen, I, F, S, R),

   MAKE_INTRINSIC_3("_fits_insert_rows", insert_rows, I, F, I, I),
   MAKE_INTRINSIC_3("_fits_delete_rows", delete_rows, I, F, I, I),

   MAKE_INTRINSIC_4("_fits_insert_cols", insert_cols, I, F, I, A, A),
   MAKE_INTRINSIC_2("_fits_delete_col", delete_col, I, F, I),

   MAKE_INTRINSIC_2("_fits_get_num_cols", get_num_cols, I, F, R),
   MAKE_INTRINSIC_2("_fits_get_rowsize", get_rowsize, I, F, R),
   MAKE_INTRINSIC_2("_fits_get_num_rows", get_num_rows, I, F, R),
   MAKE_INTRINSIC_5("_fits_write_col", write_col, I, F, I, I, I, A),
   MAKE_INTRINSIC_5("_fits_read_col", read_col, I, F, I, I, I, R),
   MAKE_INTRINSIC_3("_fits_get_keytype", get_keytype, I, F, S, R),
   MAKE_INTRINSIC_1("_fits_get_keyclass", get_keyclass, I, S),

   MAKE_INTRINSIC_0("_fits_read_cols", read_cols, I),
   MAKE_INTRINSIC_3("_fits_set_bscale", set_bscale, I, F, D, D),
   MAKE_INTRINSIC_4("_fits_set_tscale", set_tscale, I, F, I, D, D),

   /* checksum routines */
   MAKE_INTRINSIC_1("_fits_write_chksum", write_chksum, I, F),
   MAKE_INTRINSIC_1("_fits_update_chksum", update_chksum, I, F),
   MAKE_INTRINSIC_3("_fits_verify_chksum", verify_chksum, I, F, R, R),
   MAKE_INTRINSIC_3("_fits_get_chksum", get_chksum, I, F, R, R),

   MAKE_INTRINSIC_0("_fits_get_version", get_version, SLANG_VOID_TYPE),
   SLANG_END_INTRIN_FUN_TABLE
};

static SLang_IConstant_Type IConst_Table [] =
{
   MAKE_ICONSTANT("_FITS_BINARY_TBL", BINARY_TBL),
   MAKE_ICONSTANT("_FITS_ASCII_TBL", ASCII_TBL),
   MAKE_ICONSTANT("_FITS_IMAGE_HDU", IMAGE_HDU),

   MAKE_ICONSTANT("_FITS_SAME_FILE",	SAME_FILE),
   MAKE_ICONSTANT("_FITS_TOO_MANY_FILES",	TOO_MANY_FILES),
   MAKE_ICONSTANT("_FITS_FILE_NOT_OPENED",	FILE_NOT_OPENED),
   MAKE_ICONSTANT("_FITS_FILE_NOT_CREATED",	FILE_NOT_CREATED),
   MAKE_ICONSTANT("_FITS_WRITE_ERROR",	WRITE_ERROR),
   MAKE_ICONSTANT("_FITS_END_OF_FILE",	END_OF_FILE),
   MAKE_ICONSTANT("_FITS_READ_ERROR",	READ_ERROR),
   MAKE_ICONSTANT("_FITS_FILE_NOT_CLOSED",	FILE_NOT_CLOSED),
   MAKE_ICONSTANT("_FITS_ARRAY_TOO_BIG",	ARRAY_TOO_BIG),
   MAKE_ICONSTANT("_FITS_READONLY_FILE",	READONLY_FILE),
   MAKE_ICONSTANT("_FITS_MEMORY_ALLOCATION",	MEMORY_ALLOCATION),
   MAKE_ICONSTANT("_FITS_BAD_FILEPTR",	BAD_FILEPTR),
   MAKE_ICONSTANT("_FITS_NULL_INPUT_PTR",	NULL_INPUT_PTR),
   MAKE_ICONSTANT("_FITS_SEEK_ERROR",	SEEK_ERROR),
   MAKE_ICONSTANT("_FITS_BAD_URL_PREFIX",	BAD_URL_PREFIX),
   MAKE_ICONSTANT("_FITS_TOO_MANY_DRIVERS",	TOO_MANY_DRIVERS),
   MAKE_ICONSTANT("_FITS_DRIVER_INIT_FAILED",	DRIVER_INIT_FAILED),
   MAKE_ICONSTANT("_FITS_NO_MATCHING_DRIVER",	NO_MATCHING_DRIVER),
   MAKE_ICONSTANT("_FITS_URL_PARSE_ERROR",	URL_PARSE_ERROR),
   MAKE_ICONSTANT("_FITS_HEADER_NOT_EMPTY",	HEADER_NOT_EMPTY),
   MAKE_ICONSTANT("_FITS_KEY_NO_EXIST",	KEY_NO_EXIST),
   MAKE_ICONSTANT("_FITS_KEY_OUT_BOUNDS",	KEY_OUT_BOUNDS),
   MAKE_ICONSTANT("_FITS_VALUE_UNDEFINED",	VALUE_UNDEFINED),
   MAKE_ICONSTANT("_FITS_NO_QUOTE",	NO_QUOTE),
   MAKE_ICONSTANT("_FITS_BAD_KEYCHAR",	BAD_KEYCHAR),
   MAKE_ICONSTANT("_FITS_BAD_ORDER",	BAD_ORDER),
   MAKE_ICONSTANT("_FITS_NOT_POS_INT",	NOT_POS_INT),
   MAKE_ICONSTANT("_FITS_NO_END",	NO_END),
   MAKE_ICONSTANT("_FITS_BAD_BITPIX",	BAD_BITPIX),
   MAKE_ICONSTANT("_FITS_BAD_NAXIS",	BAD_NAXIS),
   MAKE_ICONSTANT("_FITS_BAD_NAXES",	BAD_NAXES),
   MAKE_ICONSTANT("_FITS_BAD_PCOUNT",	BAD_PCOUNT),
   MAKE_ICONSTANT("_FITS_BAD_GCOUNT",	BAD_GCOUNT),
   MAKE_ICONSTANT("_FITS_BAD_TFIELDS",	BAD_TFIELDS),
   MAKE_ICONSTANT("_FITS_NEG_WIDTH",	NEG_WIDTH),
   MAKE_ICONSTANT("_FITS_NEG_ROWS",	NEG_ROWS),
   MAKE_ICONSTANT("_FITS_COL_NOT_FOUND",	COL_NOT_FOUND),
   MAKE_ICONSTANT("_FITS_BAD_SIMPLE",	BAD_SIMPLE),
   MAKE_ICONSTANT("_FITS_NO_SIMPLE",	NO_SIMPLE),
   MAKE_ICONSTANT("_FITS_NO_BITPIX",	NO_BITPIX),
   MAKE_ICONSTANT("_FITS_NO_NAXIS",	NO_NAXIS),
   MAKE_ICONSTANT("_FITS_NO_NAXES",	NO_NAXES),
   MAKE_ICONSTANT("_FITS_NO_XTENSION",	NO_XTENSION),
   MAKE_ICONSTANT("_FITS_NOT_ATABLE",	NOT_ATABLE),
   MAKE_ICONSTANT("_FITS_NOT_BTABLE",	NOT_BTABLE),
   MAKE_ICONSTANT("_FITS_NO_PCOUNT",	NO_PCOUNT),
   MAKE_ICONSTANT("_FITS_NO_GCOUNT",	NO_GCOUNT),
   MAKE_ICONSTANT("_FITS_NO_TFIELDS",	NO_TFIELDS),
   MAKE_ICONSTANT("_FITS_NO_TBCOL",	NO_TBCOL),
   MAKE_ICONSTANT("_FITS_NO_TFORM",	NO_TFORM),
   MAKE_ICONSTANT("_FITS_NOT_IMAGE",	NOT_IMAGE),
   MAKE_ICONSTANT("_FITS_BAD_TBCOL",	BAD_TBCOL),
   MAKE_ICONSTANT("_FITS_NOT_TABLE",	NOT_TABLE),
   MAKE_ICONSTANT("_FITS_COL_TOO_WIDE",	COL_TOO_WIDE),
   MAKE_ICONSTANT("_FITS_COL_NOT_UNIQUE",	COL_NOT_UNIQUE),
   MAKE_ICONSTANT("_FITS_BAD_ROW_WIDTH",	BAD_ROW_WIDTH),
   MAKE_ICONSTANT("_FITS_UNKNOWN_EXT",	UNKNOWN_EXT),
   MAKE_ICONSTANT("_FITS_UNKNOWN_REC",	UNKNOWN_REC),
   MAKE_ICONSTANT("_FITS_END_JUNK",	END_JUNK),
   MAKE_ICONSTANT("_FITS_BAD_HEADER_FILL",	BAD_HEADER_FILL),
   MAKE_ICONSTANT("_FITS_BAD_DATA_FILL",	BAD_DATA_FILL),
   MAKE_ICONSTANT("_FITS_BAD_TFORM",	BAD_TFORM),
   MAKE_ICONSTANT("_FITS_BAD_TFORM_DTYPE",	BAD_TFORM_DTYPE),
   MAKE_ICONSTANT("_FITS_BAD_TDIM",	BAD_TDIM),
   MAKE_ICONSTANT("_FITS_BAD_HDU_NUM",	BAD_HDU_NUM),
   MAKE_ICONSTANT("_FITS_BAD_COL_NUM",	BAD_COL_NUM),
   MAKE_ICONSTANT("_FITS_NEG_FILE_POS",	NEG_FILE_POS),
   MAKE_ICONSTANT("_FITS_NEG_BYTES",	NEG_BYTES),
   MAKE_ICONSTANT("_FITS_BAD_ROW_NUM",	BAD_ROW_NUM),
   MAKE_ICONSTANT("_FITS_BAD_ELEM_NUM",	BAD_ELEM_NUM),
   MAKE_ICONSTANT("_FITS_NOT_ASCII_COL",	NOT_ASCII_COL),
   MAKE_ICONSTANT("_FITS_NOT_LOGICAL_COL",	NOT_LOGICAL_COL),
   MAKE_ICONSTANT("_FITS_BAD_ATABLE_FORMAT",	BAD_ATABLE_FORMAT),
   MAKE_ICONSTANT("_FITS_BAD_BTABLE_FORMAT",	BAD_BTABLE_FORMAT),
   MAKE_ICONSTANT("_FITS_NO_NULL",	NO_NULL),
   MAKE_ICONSTANT("_FITS_NOT_VARI_LEN",	NOT_VARI_LEN),
   MAKE_ICONSTANT("_FITS_BAD_DIMEN",	BAD_DIMEN),
   MAKE_ICONSTANT("_FITS_BAD_PIX_NUM",	BAD_PIX_NUM),
   MAKE_ICONSTANT("_FITS_ZERO_SCALE",	ZERO_SCALE),
   MAKE_ICONSTANT("_FITS_NEG_AXIS",	NEG_AXIS),

   /* get_keyclass return value */
   MAKE_ICONSTANT("_FITS_TYP_STRUC_KEY",TYP_STRUC_KEY),
   MAKE_ICONSTANT("_FITS_TYP_CMPRS_KEY",TYP_CMPRS_KEY),
   MAKE_ICONSTANT("_FITS_TYP_SCAL_KEY",	TYP_SCAL_KEY),
   MAKE_ICONSTANT("_FITS_TYP_NULL_KEY",	TYP_NULL_KEY),
   MAKE_ICONSTANT("_FITS_TYP_DIM_KEY",	TYP_DIM_KEY),
   MAKE_ICONSTANT("_FITS_TYP_RANG_KEY",	TYP_RANG_KEY),
   MAKE_ICONSTANT("_FITS_TYP_UNIT_KEY",	TYP_UNIT_KEY),
   MAKE_ICONSTANT("_FITS_TYP_DISP_KEY",	TYP_DISP_KEY),
   MAKE_ICONSTANT("_FITS_TYP_HDUID_KEY",TYP_HDUID_KEY),
   MAKE_ICONSTANT("_FITS_TYP_CKSUM_KEY",TYP_CKSUM_KEY),
   MAKE_ICONSTANT("_FITS_TYP_WCS_KEY",	TYP_WCS_KEY),
   MAKE_ICONSTANT("_FITS_TYP_REFSYS_KEY",TYP_REFSYS_KEY),
   MAKE_ICONSTANT("_FITS_TYP_COMM_KEY",	TYP_COMM_KEY),
   MAKE_ICONSTANT("_FITS_TYP_CONT_KEY",	TYP_CONT_KEY),
   MAKE_ICONSTANT("_FITS_TYP_USER_KEY",	TYP_USER_KEY),

   MAKE_ICONSTANT("_cfitsio_module_version", MODULE_VERSION_NUMBER),

   SLANG_END_ICONST_TABLE
};

static char *Module_Version_String = MODULE_VERSION_STRING;
static SLang_Intrin_Var_Type Intrin_Vars[] =
{
   MAKE_VARIABLE("_cfitsio_module_version_string", &Module_Version_String, SLANG_STRING_TYPE, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

static void patchup_intrinsic_table (void)
{
   SLang_Intrin_Fun_Type *f;

   f = Fits_Intrinsics;
   while (f->name != NULL)
     {
	unsigned int i, nargs;
	SLtype *args;

	nargs = f->num_args;
	args = f->arg_types;
	for (i = 0; i < nargs; i++)
	  {
	     if (args[i] == DUMMY_FITS_FILE_TYPE)
	       args[i] = Fits_Type_Id;
	  }

	/* For completeness */
	if (f->return_type == DUMMY_FITS_FILE_TYPE)
	  f->return_type = Fits_Type_Id;

	f++;
     }
}

static void free_fits_file_type (SLtype type, VOID_STAR f)
{
   FitsFile_Type *ft;
   int status = 0;

   (void) type;

   ft = (FitsFile_Type *) f;
   if (ft->fptr != NULL)
     fits_close_file (ft->fptr, &status);

   SLfree ((char *) ft);
}

static int check_version (void)
{
   float compiled_version = 0;
   float linked_version = 0;
   float tol = 0.0001;

#ifdef CFITSIO_VERSION
   compiled_version = CFITSIO_VERSION;
   (void) fits_get_version (&linked_version);
#endif

   if (fabs (linked_version - compiled_version) <= tol)
     return 0;

   fprintf (stderr, "\n\
***WARNING: The version of CFITSIO that this module is linked against (%g)\n\
   is not the same as the version it was compiled against (%g).\n\
   As the CFITSIO developers make no guarantees of binary compatibility,\n\
   you may experience problems with this module.  You are stongly urged to\n\
   recompile the module.\n\n\
", linked_version, compiled_version);

   return -1;
}

int init_cfitsio_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns;

   ns = SLns_create_namespace (ns_name);
   if (ns == NULL)
     return -1;

   if (Fits_Type_Id == 0)
     {
	SLang_Class_Type *cl;

	(void) check_version ();

	cl = SLclass_allocate_class ("Fits_File_Type");
	if (cl == NULL) return -1;
	(void) SLclass_set_destroy_function (cl, free_fits_file_type);

	/* By registering as SLANG_VOID_TYPE, slang will dynamically allocate a
	 * type.
	 */
	if (-1 == SLclass_register_class (cl, SLANG_VOID_TYPE,
					  sizeof (FitsFile_Type),
					  SLANG_CLASS_TYPE_MMT))
	  return -1;

	Fits_Type_Id = SLclass_get_class_id (cl);
	patchup_intrinsic_table ();
     }

   if (-1 == SLns_add_intrin_fun_table (ns, Fits_Intrinsics, "__CFITSIO__"))
     return -1;

   if (-1 == SLns_add_iconstant_table (ns, IConst_Table, NULL))
     return -1;

   if (-1 == SLns_add_intrin_var_table (ns, Intrin_Vars, NULL))
     return -1;

   return 0;
}

