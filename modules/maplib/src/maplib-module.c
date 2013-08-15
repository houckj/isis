/* -*- mode: C; mode: fold; -*- */
/*
  Copyright (c) 2004-2008 Massachusetts Institute of Technology

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

/* Author: John E. Davis (davis@space.mit.edu) */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <slang.h>
#if SLANG_VERSION < 20000
# define POP_DOUBLE(x,y,z) SLang_pop_double(x,y,z)
#else
# define POP_DOUBLE(x,y,z) SLang_pop_double(x)
#endif

#ifdef __cplusplus
extern "C" 
{
#endif
SLANG_MODULE(maplib);
#ifdef __cplusplus
}
#endif

#include "version.h"

#include <math.h>
#ifndef PI
# define PI 3.14159265358979323846264338327950288
#endif

#ifdef DEGREES
# undef DEGREES
#endif
#define DEGREES(x) ((180.0/PI)*(x))
#ifdef RADIANS
# undef RADIANS
#endif
#define RADIANS(x) ((PI/180.0)*(x))

/* In this file, theta and phi specify standard spherical coordinates.
 * whose values are given in radians.  In contrast, lat and lon represent
 * latitude and longitude coordinates specified in degrees, unless otherwise
 * specified.
 * They are related as follows:
 * 
 *    theta = PI/180*(90-lat);  theta_hat = -lat_hat
 *    phi = PI/180*lon;         phi_hat = lon_hat
 */
#define THETA_TO_LAT(theta)	(90.0 - DEGREES(theta))
#define LAT_TO_THETA(lat)	RADIANS(90.0 - (lat))
#define PHI_TO_LON(phi)		DEGREES(phi)
#define LON_TO_PHI(lon)		RADIANS(lon)

#define VECTOR_DIFF(c,a,b) \
   c[0]=a[0]-b[0]; \
   c[1]=a[1]-b[1]; \
   c[2]=a[2]-b[2]
#define VECTOR_SUM(c,a,b) \
   c[0]=a[0]+b[0]; \
   c[1]=a[1]+b[1]; \
   c[2]=a[2]+b[2]
#define VECTOR_COMBINE_SUM(c,x,s,a,t,b) \
   c[0]=x[0]+(s)*a[0]+(t)*b[0]; \
   c[1]=x[1]+(s)*a[1]+(t)*b[1]; \
   c[2]=x[2]+(s)*a[2]+(t)*b[2];
#define VECTOR_A_PLUS_BT(c,a,b,t) \
   c[0]=a[0]+(t)*b[0]; \
   c[1]=a[1]+(t)*b[1]; \
   c[2]=a[2]+(t)*b[2]
#define VECTOR_MUL(b,a,t) \
   b[0] = (t)*a[0]; \
   b[1] = (t)*a[1]; \
   b[2] = (t)*a[2];
#define VECTOR_ASSIGN(b,a) \
   b[0] = a[0]; \
   b[1] = a[1]; \
   b[2] = a[2]

#define DOTPROD(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

static double Maplib_NaN;

/*{{{ Spherical Coordinate Routines */

static double hypot3d (double x, double y, double z)
{
   double tmp;

   x = fabs(x);
   y = fabs(y);
   z = fabs(z);

   /* Want x <= y <= z */
   if (y < x)
     {
	tmp = x;
	x = y;
	y = tmp;
     }
   
   if (z < y)
     {
	tmp = z;
	if (z < x)
	  {
	     z = x;
	     x = tmp;
	  }
	else
	  {
	     z = y;
	     y = tmp;
	  }
     }
   
   if (z == 0.0)
     return 0.0;
   
   x /= z;
   y /= z;
   return z * sqrt (1.0 + x*x + y*y);
}


static void vector_to_spherical (double x, double y, double z, 
				 double *rp, double *thetap, double *phip)
{
   double r = hypot (x, y);

   if (r == 0.0)
     {
	if (z < 0)
	  {
	     *thetap = PI;
	     *rp = -z;
	  }
	else
	  {
	     *thetap = 0.0;
	     *rp = z;
	  }
	*phip = 0.0;
	return;
     }

   *thetap = atan2 (r, z);
   *phip = atan2 (y, x);
   
   *rp = hypot (r, z);
}

static void spherical_to_vector (double r, double theta, double phi, double v[3])
{
   double ct = r*cos(theta);
   double st = r*sin(theta);
   double cp = cos(phi);
   double sp = sin(phi);
   
   v[0] = st*cp;
   v[1] = st*sp;
   v[2] = ct;
}

static void get_spherical_basis (double theta, double phi, double *rhat, double *theta_hat, double *phi_hat)
{
   double ct = cos(theta);
   double st = sin(theta);
   double cp = cos(phi);
   double sp = sin(phi);

   if (rhat != NULL)
     {
	rhat[0] = st*cp; rhat[1] = st*sp; rhat[2] = ct;
     }
   if (theta_hat != NULL)
     {
	theta_hat[0] = ct*cp; theta_hat[1] = ct*sp; theta_hat[2] = -st;
     }
   if (phi_hat != NULL)
     {
	phi_hat[0] = -sp; phi_hat[1] = cp; phi_hat[2] = 0.0;
     }
}

/* Returns -r <= x <= r */
static double normalize_to_range (double x, double r)
{
   double r2;

   if ((x >= -r) && (x <= r))
     return x;
   
   r2 = 2.0 * r;
   x = fmod (x, r2);
   if (x < -r)
     x += r2;
   else if (x > r)
     x -= r2;
   return x;
}

static void normalize_lon_lat_degrees (double *lonp, double *latp)
{
   double lon, lat;
   
   lon = normalize_to_range (*lonp, 180.0);
   lat = normalize_to_range (*latp, 180.0);

   if (lat > 90.0)
     {
	lat = 180.0 - lat;
	lon = normalize_to_range (lon-180.0, 180.0);
     }
   else if (lat < -90)
     {
	lat = 180.0 + lat;
	lon = normalize_to_range (lon-180.0, 180.0);
     }
   
   *lonp = lon;
   *latp = lat;
}

static int normalize_vector (double v[3])
{
   double len = hypot3d (v[0],v[1],v[2]);
   if (len == 0.0)
     return -1;
   len = 1.0/len;
   VECTOR_MUL(v,v,len);
   return 0;
}

/* Given two sets of othonormal systems, (x1,x2,x3) and (y1,y2,y3), find the
 * rotation angle and direction that rotates (x1,x2,x3) into (y1,y2,y3).
 */
static int find_rotation (double x_1[3], double x_2[3], double x_3[3],
			  double y_1[3], double y_2[3], double y_3[3],
			  double axis[3], double *cos_thetap, double *sin_thetap)
{
   double cos_theta, sin_theta;
   double omega[3];

   cos_theta = 0.5 * (DOTPROD(y_1,x_1) + DOTPROD(y_2,x_2) + DOTPROD(y_3,x_3) - 1.0);
   if (cos_theta > 1.0) cos_theta = 1.0; 
   else if (cos_theta < -1.0) cos_theta = -1.0;

   sin_theta = sqrt (1.0 - cos_theta*cos_theta);
   if (sin_theta == 0.0)
     {
	if (cos_theta > 0.0)
	  {
	     /* Identity transform */
	     omega[0] = omega[1] = omega[2] = 0.0;
	  }
	else
	  {
	     omega[0] = sqrt (0.5*(1.0 + DOTPROD(y_1,x_1)));
	     omega[1] = sqrt (0.5*(1.0 + DOTPROD(y_2,x_2)));
	     omega[2] = sqrt (0.5*(1.0 + DOTPROD(y_3,x_3)));
	     /* Now determine the signs. */
	     if ((omega[0] >= omega[1]) && (omega[0] >= omega[2]))
	       {
		  if (DOTPROD(y_1,x_2) < 0) omega[1] = -omega[1];
		  if (DOTPROD(y_1,x_3) < 0) omega[2] = -omega[2];
	       }
	     else if ((omega[1] >= omega[0]) && (omega[1] >= omega[2]))
	       {
		  if (DOTPROD(y_2,x_1) < 0) omega[0] = -omega[0];
		  if (DOTPROD(y_2,x_3) < 0) omega[2] = -omega[2];
	       }
	     else
	       {
		  if (DOTPROD(y_3,x_1) < 0) omega[0] = -omega[0];
		  if (DOTPROD(y_3,x_2) < 0) omega[1] = -omega[1];
	       }
	     (void) normalize_vector (omega);
	  }
     }
   else
     {
	double factor = 0.5/sin_theta;
	omega[0] = factor * (DOTPROD(y_2,x_3) - DOTPROD(y_3,x_2));
	omega[1] = factor * (DOTPROD(y_3,x_1) - DOTPROD(y_1,x_3));
	omega[2] = factor * (DOTPROD(y_1,x_2) - DOTPROD(y_2,x_1));
	(void) normalize_vector (omega);
     }
   axis[0] = omega[0]*x_1[0] + omega[1]*x_2[0] + omega[2]*x_3[0];
   axis[1] = omega[0]*x_1[1] + omega[1]*x_2[1] + omega[2]*x_3[1];
   axis[2] = omega[0]*x_1[2] + omega[1]*x_2[2] + omega[2]*x_3[2];
   
   *cos_thetap = cos_theta;
   *sin_thetap = sin_theta;
   return 0;
}

static void create_lonlat_basis (double lon, double lat,
				 double n_hat[3], double lon_hat[3], double lat_hat[3])
{
   double cos_lat = cos (lat);
   double sin_lat = sin (lat);
   double cos_lon = cos (lon);
   double sin_lon = sin (lon);
   
   n_hat[0] = cos_lat*cos_lon; n_hat[1] = cos_lat*sin_lon; n_hat[2] = sin_lat;
   lon_hat[0] = -sin_lon; lon_hat[1] = cos_lon; lon_hat[2] = 0.0;
   lat_hat[0] = -sin_lat*cos_lon; lat_hat[1] = -sin_lat*sin_lon; lat_hat[2] = cos_lat;
}

/* Rotate a vector v about unit vector n by theta */
static void rotate_vector (double v[3], double n[3], double cos_theta, double sin_theta)
{
   double v_dot_n;
   double v0, v1, v2, n0, n1, n2;
   
   v0 = v[0]; v1 = v[1]; v2 = v[2];
   n0 = n[0]; n1 = n[1]; n2 = n[2];

   v_dot_n = v0*n0 + v1*n1 + v2*n2;
   v_dot_n *= (1.0 - cos_theta);
   v[0] = v0*cos_theta + v_dot_n*n0 + sin_theta*(n1*v2 - n2*v1);
   v[1] = v1*cos_theta + v_dot_n*n1 + sin_theta*(n2*v0 - n0*v2);
   v[2] = v2*cos_theta + v_dot_n*n2 + sin_theta*(n0*v1 - n1*v0);
}
   
typedef struct 
{
   double theta;
   double cos_theta, sin_theta;
   double omega[3];
}
Sphere_Rotate_Type;
/* Initialize a transform that maps (lon0,lat0) --> (lon1,lat1) and applies
 * an additional rotation about the final point by beta.
 */
static int init_sphere_xform (Sphere_Rotate_Type *s, 
			      double lon0, double lat0, 
			      double lon1, double lat1, 
			      double beta)
{
   double lon0_hat[3], lat0_hat[3], n0_hat[3];
   double lon1_hat[3], lat1_hat[3], n1_hat[3];

   create_lonlat_basis (lon0, lat0, n0_hat, lon0_hat, lat0_hat);
   create_lonlat_basis (lon1, lat1, n1_hat, lon1_hat, lat1_hat);
   
   if (beta != 0.0)
     {
	/* Rotate lon and lat vectors about the normal by beta */
	double cos_beta = cos (beta);
	double sin_beta = sin (beta);
	rotate_vector (lon1_hat, n1_hat, cos_beta, sin_beta);
	rotate_vector (lat1_hat, n1_hat, cos_beta, sin_beta);
     }
   find_rotation (n0_hat, lon0_hat, lat0_hat, n1_hat, lon1_hat, lat1_hat, 
		  s->omega, &s->cos_theta, &s->sin_theta);
   s->theta = atan2 (s->sin_theta, s->cos_theta);

   return 0;
}

static void rotate_sphere (Sphere_Rotate_Type *s, int dir, 
			   double lon, double lat, 
			   double *lonp, double *latp)
{
   double p[3], pz;
   double cos_lat = cos (lat);

   p[0] = cos_lat * cos(lon);
   p[1] = cos_lat * sin(lon);
   p[2] = sin(lat);
   
   /* FIXME: : If rotating about the pole, and if lon,lat corresponds to the
    * pole, then the returned lon should be lon+rotation angle
    */
   if (((lat == PI/2.0) || (lat == -PI/2.0)) 
       && ((s->omega[2] == 1.0) || (s->omega[2] == -1.0)))
     {
	*lonp = normalize_to_range (lon + s->omega[2] * s->theta, PI);
	*latp = lat;
	return;
     }
   
   rotate_vector (p, s->omega, s->cos_theta, dir*s->sin_theta);
   lon = atan2 (p[1], p[0]);
   pz = p[2];
   if (pz > 1.0) pz = 1.0; else if (pz < -1.0) pz = -1.0;
   lat = asin (pz);
   
   *lonp = lon;
   *latp = lat;
}

/*}}}*/

/*{{{ Generic Tangent Plane projection routines */

/* These routines assume that the plane is tangent to the unit sphere at some point
 * s0.  The normal to the plane is also s0.
 * 
 * The location of the observer from where rays eminate is given by p.
 * 
 * The vector x represents a point on the plane.
 * 
 * The basis vectors at the point s0 are aligned with the local 
 * lon,lat vectors.
 */

/* If solution_type is 1, then we are seeking the primary soln. */
static int project_to_plane (double s0[3], double p, double phat[3], 
			     double s[3], 
			     /* int solution_type, */ 
			     double x[3])
{
   double t, t1;
   double phat_dot_s0, s_dot_s0;
   double a, b;
   double pt1[3];
   int solution_type = 1;

   if (p > 1.0)
     {				       /* outside sphere */
	if (DOTPROD(phat,s) < 1.0/p)
	  {
	     if (solution_type == 1)
	       {
		  x[0] = x[1] = x[2] = Maplib_NaN;
		  return -1;
	       }
	  }
	else if (solution_type != 1)
	  {
	     x[0] = x[1] = x[2] = Maplib_NaN;
	     return -1;
	  }
     }

   /* See pg 9 of maplib notebook */
   phat_dot_s0 = DOTPROD(phat,s0);
   s_dot_s0 = DOTPROD(s,s0);
   if (p >= 1.0)
     {
	a = phat_dot_s0 - 1.0/p;
	b = phat_dot_s0 - s_dot_s0/p;
	t = a/b;
	t1 = (1.0 - s_dot_s0)/b;       /* this is p(1-t) */
     }
   else
     {
	a = 1.0 - p*phat_dot_s0;
	b = s_dot_s0 - p*phat_dot_s0;
	t=a/b;
	t1 = p*(1-t);
     }
   VECTOR_MUL(pt1, phat, t1);
   VECTOR_A_PLUS_BT(x,pt1,s,t);
   return 0;
}

static int project_to_plane_coords (double s[3], double s0[3], double p, double phat[3],
				    double txhat[3], double tyhat[3],
				    double *tx, double *ty)
{
   int ret;
   double x[3];

   ret = project_to_plane (s0, p, phat, s, x);
   /* We have:
    *   x - s0 = tx*tx_hat + ty*ty_hat
    */
   
   VECTOR_DIFF(x,x,s0);
   *tx = DOTPROD(x,txhat);
   *ty = DOTPROD(x,tyhat);
   return ret;
}

				    

static int project_from_plane (double s0[3], double p, double phat[3], double x[3], double s[3])
{
   double t, t1;
   double a, b, c;
   double phat_dot_x, p2, x2, p_dot_x;
   double alpha, beta;
   int solution_type = 1;

   phat_dot_x = DOTPROD(phat,x);
   p_dot_x = p*phat_dot_x;
   x2 = DOTPROD(x,x);
   p2 = p*p;
   if (p2 > 1.0)
     {
	int sgn = 1;

	/* Pages 18,19 of maplib notebook */
	a = 1.0 - 1.0/p2;
	b = (phat_dot_x/p - 1.0)/a;
	c = (1 + x2/p2 - 2.0*phat_dot_x/p)/a;
	alpha = -b;
	beta = sqrt (b*b-c);
	if (alpha >= 0)
	  {
	     if (solution_type != 1)
	       sgn = -1;
	  }
	else if (solution_type == 1)
	  sgn = -1;
	t = alpha + sgn * beta;

	/* For large p, 1-t goes like 1/p. */
	alpha = (-1.0/p + phat_dot_x)/a;
	beta = sqrt (-x2*a + 1.0 + phat_dot_x*(phat_dot_x-2.0/p));
	beta = beta/a;
	t1 = alpha - sgn*beta;	       /* (1-t)*p */
     }
   else 
     {
	/* Pg 17 of maplib notebook */
	int sgn = 1;
	b = 2.0 * (p_dot_x - p2);
	c = (x2 + p2 - 2.0*p_dot_x);
	if (p2 == 1.0)
	  t = -c/b;
	else
	  {
	     a = p2 - 1.0;
	     b = b/(2.0*a);
	     c = c / a;
	     alpha = -b;
	     beta = sqrt (b*b - c);
	     if (DOTPROD(s0,x) >= p*DOTPROD(s0,phat))
	       {
		  if (solution_type != 1)
		    sgn = -1;
	       }
	     else if (solution_type == 1)
	       sgn = -1;
	     t = alpha + sgn * beta;
	  }
	t1 = (1.0 - t)*p;
     }
   t1 = -t1;
   VECTOR_A_PLUS_BT(s,x,phat,t1);
   t = 1.0/t;
   VECTOR_MUL(s,s,t);
   return 0;
}

static int project_from_plane_coords (double tx, double ty,
				      double s0[3], double p, double phat[3], 
				      double txhat[3], double tyhat[3],
				      double s[3])
{
   int ret;
   double x[3];
   
   VECTOR_COMBINE_SUM(x,s0,tx,txhat,ty,tyhat);
   ret = project_from_plane (s0, p, phat, x, s);

   return ret;
}

static int compute_plane_vectors (double theta_0, double phi_0,   /* s0 */
				  double xhat[3], double yhat[3], double zhat[3])
{
   double theta0_hat[3], phi0_hat[3];

   get_spherical_basis (theta_0, phi_0, zhat, theta0_hat, phi0_hat);

   VECTOR_ASSIGN (xhat, phi0_hat);
   VECTOR_ASSIGN (yhat, theta0_hat);
   VECTOR_MUL(yhat,yhat,-1.0);	       /* lat_hat = -theta_hat */
   
   return 0;
}


/*}}}*/

/*{{{ Tangent Plane Projection Routines */

/* This routine returns a vector from the center of the sphere to the point
 * (tx,ty) on the tangent plane.  This vector does not depend upon the
 * specific projection.
 */
static void tangent_plane_to_vector (double ex[3], double ey[3], double ez[3],
				     double tx, double ty, 
				     double v[3])
{
   VECTOR_COMBINE_SUM(v,ez,tx,ex,ty,ey);
}

/* In this routine, v is a vector from the center of the sphere to the 
 * tangent plane at v0.  The tangent plane coordinates are returned.
 * 
 * v*t = v0 + tx*ex+ty*ey
 * t = 1.0/(v.v0)
 * 
 */
static int vector_to_tangent_plane (double v[3],
				    double ex[3], double ey[3], double ez[3],
				    double *tx, double *ty)
{
   double v0, v1, v2;
   double t = DOTPROD(v,ez);
   int ret;
   
   /* If t == 0, then we will get infinity.  Let's hope IEEE Inf is available */
   if (t == 0.0)
     ret = -1;
   else
     ret = 0;

   t = 1.0/t;

   v0 = v[0]*t;
   v1 = v[1]*t;
   v2 = v[2]*t;

   *tx = v0*ex[0] + v1*ex[1] + v2*ex[2];
   *ty = v0*ey[0] + v1*ey[1] + v2*ey[2];

   return ret;
}


/* In this routine, p is a unit vector on the sphere to be projected from 
 * infinity to the tangent plane located at ez.
 * The tangent plane coordinates are returned.
 * 
 * v = p + (1-p.ez)ez;
 */
static int ortho_project_vector (double p[3],
				 double ex[3], double ey[3], double ez[3],
				 double *tx, double *ty)
{
   double v[3];
   double c;
   c = (1.0 - DOTPROD(p,ez));
   VECTOR_A_PLUS_BT(v,p,ez,c);
   
   return vector_to_tangent_plane (v, ex, ey, ez, tx, ty);
}

/* This routine is the inverse of ortho_project_vector */
static int ortho_deproject_vector (double ex[3], double ey[3], double ez[3],
				   double tx, double ty, 
				   double p[3])
{
   double x[3];
   double c;
   int ret = 0;

   tangent_plane_to_vector (ex, ey, ez, tx, ty, x);
   c = 2.0 - DOTPROD(x,x);
   if (c < 0.0)
     ret = -1;
   
   c = sqrt (c) - 1.0;
   VECTOR_A_PLUS_BT(p, x, ez, c);
   return ret;
}


/*}}}*/

/*{{{ push/pop functions */

static int pop_2_floating_arrays_1 (SLang_Array_Type **atp, SLang_Array_Type **btp,
				    int *is_scalarp)
{
   int btype, atype, ctype;
   SLang_Array_Type *at, *bt;
   int is_scalar;

   if (-1 == (btype = SLang_peek_at_stack1 ()))
     return -1;
   
   is_scalar = (SLANG_ARRAY_TYPE != SLang_peek_at_stack ());

   if (-1 == SLroll_stack (2))
     return -1;
   
   atype = SLang_peek_at_stack1 ();
   is_scalar = is_scalar && (SLANG_ARRAY_TYPE != SLang_peek_at_stack ());

   if ((atype == SLANG_FLOAT_TYPE) && (btype == SLANG_FLOAT_TYPE))
     ctype = SLANG_FLOAT_TYPE;
   else 
     ctype = SLANG_DOUBLE_TYPE;

   if (-1 == SLang_pop_array_of_type (&at, ctype))
     return -1;
   if (-1 == SLang_pop_array_of_type (&bt, ctype))
     {
	SLang_free_array (at);
	return -1;
     }
   *atp = at;
   *btp = bt;
   *is_scalarp = is_scalar;
   return 0;
}

static int pop_2_floating_arrays (SLang_Array_Type **atp, SLang_Array_Type **btp, 
				  int *is_scalarp)
{
   SLang_Array_Type *at, *bt;
   int is_scalar;

   if (-1 == pop_2_floating_arrays_1 (&at, &bt, &is_scalar))
     return -1;

   if (at->num_elements != bt->num_elements)
     {
	SLang_verror (SL_TYPE_MISMATCH, "Expecting two arrays of the same size");
	SLang_free_array (at);
	SLang_free_array (bt);
	return -1;
     }
   
   *atp = at;
   *btp = bt;
   *is_scalarp = is_scalar;
   return 0;
}

static SLang_Array_Type *create_or_reuse_array (SLang_Array_Type *at)
{
   if ((at->num_refs == 1) && (0 == (at->flags & SLARR_DATA_VALUE_IS_READ_ONLY)))
     {
	at->num_refs++;
	return at;
     }

   return SLang_create_array (at->data_type, 0, NULL, at->dims, at->num_dims);
}

static int create_or_reuse_2_arrays (SLang_Array_Type *at0, SLang_Array_Type *at1,
				     SLang_Array_Type **bt0p, SLang_Array_Type **bt1p)
{
   SLang_Array_Type *bt0, *bt1;

   if (NULL == (bt0 = create_or_reuse_array (at0)))
     return -1;

   if (NULL == (bt1 = create_or_reuse_array (at1)))
     {
	SLang_free_array (bt0);
	return -1;
     }

   *bt0p = bt0;
   *bt1p = bt1;
   return 0;
}

static int pop_reusable_arrays (SLang_Array_Type **xp, SLang_Array_Type **yp,
				SLang_Array_Type **x1p, SLang_Array_Type **y1p,
				int *is_scalarp)
{
   SLang_Array_Type *x, *y;

   if (-1 == pop_2_floating_arrays (&x, &y, is_scalarp))
     return -1;

   if (-1 == create_or_reuse_2_arrays (x, y, x1p, y1p))
     {
	SLang_free_array (y);
	SLang_free_array (x);
	return -1;
     }
   
   *xp = x;
   *yp = y;
   return 0;
}


   
static int push_array_maybe_scalar (SLang_Array_Type *a, int is_scalar)
{
   if (is_scalar)
     return SLang_push_value (a->data_type, a->data);
   
   return SLang_push_array (a, 0);
}

static void get_tp_basis (double theta, double phi, 
			  double tx_hat[3], double ty_hat[3], double tz_hat[3])
{
   get_spherical_basis (theta, phi, tz_hat, ty_hat, tx_hat);
   ty_hat[0] = -ty_hat[0];
   ty_hat[1] = -ty_hat[1];
   ty_hat[2] = -ty_hat[2];
}

static int pop_mxm_matrix (SLang_Array_Type **ap, unsigned int m)
{
   SLang_Array_Type *a;

   if (-1 == SLang_pop_array_of_type (&a, SLANG_DOUBLE_TYPE))
     return -1;
   
   if ((a->num_dims != 2)
       || (a->dims[0] != (int)m)
       || (a->dims[1] != (int)m))
     {
	SLang_verror (SL_TYPE_MISMATCH, "Expecting a %u x %u array", m, m);
	SLang_free_array (a);
	return -1;
     }
   *ap = a;
   return 0;
}

/*}}}*/

/* Projection Methods */

#define REQUIRED_FIELDS \
   char *name; \
   double lon0, lat0		       /* map center */

/*{{{ Linear Projection */

/* The linear projection transforms (x,y) via
 *   (x',y') = (x1,y1) + A#(x-x0,y-y0)
 */
typedef struct
{
   REQUIRED_FIELDS;
   double x_0, y_0;
   double x_1, y_1;
   SLang_Array_Type *at;
   double a[4];
   double a_inv[4];
}
Linear_Projection_Type;

static SLang_CStruct_Field_Type Linear_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Linear_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Linear_Projection_Type, x_0, "x0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Linear_Projection_Type, y_0, "y0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Linear_Projection_Type, x_1, "x1", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Linear_Projection_Type, y_1, "y1", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Linear_Projection_Type, at, "A", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int linear_new (char *name)
{
   Linear_Projection_Type g;
   SLang_Array_Type *at;
   int dims[2];
   int ret;

   dims[0] = 2;
   dims[1] = 2;

   if (NULL == (at = SLang_create_array (SLANG_DOUBLE_TYPE, 0, NULL, dims, 2)))
     return -1;
   ((double *)at->data)[0] = 1.0;
   ((double *)at->data)[1] = 0.0;
   ((double *)at->data)[2] = 0.0;
   ((double *)at->data)[3] = 1.0;
   g.at = at;

   g.name = name;
   g.x_0 = g.y_0 = 0.0;
   g.x_1 = g.y_1 = 0.0;
   
   ret = SLang_push_cstruct (&g, Linear_Projection_Struct);
   SLang_free_array (at);
   return ret;
}

static void linear_free (void *g)
{
   if (g == NULL)
     return;
   
   SLang_free_cstruct (g, Linear_Projection_Struct);
   SLfree ((char *)g);
}

static int make_2x2_inverse (double a[4], double ainv[4])
{
   double det = a[0]*a[3] - a[1]*a[2];

   if (det == 0.0)
     {
	SLang_verror (SL_INVALID_PARM, "Transformation matrix must have an inverse");
	return -1;
     }
   ainv[0] = a[3]/det;
   ainv[1] = -a[1]/det;
   ainv[2] = -a[2]/det;
   ainv[3] = a[0]/det;
   
   return 0;
}

static void *linear_pop (char *name)
{
   Linear_Projection_Type *g;
   SLang_Array_Type *at;
   double *a;
   
   (void) name;

   if (NULL == (g = (Linear_Projection_Type *)SLmalloc (sizeof (Linear_Projection_Type))))
     return NULL;
   
   if (-1 == SLang_pop_cstruct (g, Linear_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }

   at = g->at;
   a = g->a;
   if (at == NULL)
     {
	a[0] = 1.0;
	a[1] = 0.0;
	a[2] = 0.0;
	a[3] = 1.0;
     }
   else 
     {
	if (-1 == SLang_push_array (at, 0))
	  {
	     linear_free (g);
	     return NULL;
	  }

	if (-1 == pop_mxm_matrix (&at, 2))
	  {
	     linear_free (g);
	     return NULL;
	  }

	a[0] = ((double *)at->data)[0];
	a[1] = ((double *)at->data)[1];
	a[2] = ((double *)at->data)[2];
	a[3] = ((double *)at->data)[3];
	SLang_free_array (at);
     }
   
   if (-1 == make_2x2_inverse (a, g->a_inv))
     {
	linear_free (g);
	return NULL;
     }

   return g;
}

static int linear_project_f_internal (double *a,
				      double x_0, double y_0,
				      double x_1, double y_1,
				      float *x, float *y,
				      float *newx, float *newy,
				      unsigned int num)
{
   double a_00 = a[0], a_01 = a[1];
   double a_10 = a[2], a_11 = a[3];
   unsigned int i;

   for (i = 0; i < num; i++)
     {
	double xx = x[i] - x_0;
	double yy = y[i] - y_0;
	newx[i] = x_1 + (a_00*xx + a_01*yy);
	newy[i] = y_1 + (a_10*xx + a_11*yy);
     }
   return 0;
}

static int linear_project_f (void *vg,
			     float *x, float *y,
			     float *newx, float *newy,
			     unsigned int num)
{
   Linear_Projection_Type *g = (Linear_Projection_Type *)vg;
   return linear_project_f_internal (g->a, g->x_0, g->y_0, g->x_1, g->y_1, x, y, newx, newy, num);
}

static int linear_deproject_f (void *vg,
			       float *x, float *y,
			       float *newx, float *newy,
			       unsigned int num)
{
   Linear_Projection_Type *g = (Linear_Projection_Type *)vg;
   return linear_project_f_internal (g->a_inv, g->x_1, g->y_1, g->x_0, g->y_0, x, y, newx, newy, num);
}

static int linear_project_d_internal (double *a, 
				      double x_0, double y_0,
				      double x_1, double y_1,
				      double *x, double *y,
				      double *newx, double *newy,
				      unsigned int num)
{
   double a_00 = a[0], a_01 = a[1];
   double a_10 = a[2], a_11 = a[3];
   unsigned int i;

   for (i = 0; i < num; i++)
     {
	double xx = x[i] - x_0;
	double yy = y[i] - y_0;
	newx[i] = x_1 + (a_00*xx + a_01*yy);
	newy[i] = y_1 + (a_10*xx + a_11*yy);
     }
   return 0;
}

static int linear_project_d (void *vg,
			     double *x, double *y,
			     double *newx, double *newy,
			     unsigned int num)
{
   Linear_Projection_Type *g = (Linear_Projection_Type *)vg;
   return linear_project_d_internal (g->a, g->x_0, g->y_0, g->x_1, g->y_1, x, y, newx, newy, num);
}

static int linear_deproject_d (void *vg,
			       double *x, double *y,
			       double *newx, double *newy,
			       unsigned int num)
{
   Linear_Projection_Type *g = (Linear_Projection_Type *)vg;
   return linear_project_d_internal (g->a_inv, g->x_1, g->y_1, g->x_0, g->y_0, x, y, newx, newy, num);
}


	
/*}}}*/

/*{{{ Sphere Projection (sphere to sphere) */

/* This mapping maps a sphere onto itself: lon,lat --> lon,lat.  It takes the
 * position of the pole, and a final rotation angle (beta) about that.
 */
typedef struct
{
   REQUIRED_FIELDS;
   double lon1, lat1;
   double beta;
   Sphere_Rotate_Type sphere_xform;
}
Sphere_Projection_Type;

static SLang_CStruct_Field_Type Sphere_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Sphere_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Sphere_Projection_Type, lon0, "lon0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Sphere_Projection_Type, lat0, "lat0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Sphere_Projection_Type, lon1, "lon1", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Sphere_Projection_Type, lat1, "lat1", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Sphere_Projection_Type, beta, "beta", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int sphere_new (char *name)
{
   Sphere_Projection_Type g;

   g.name = name;
   g.lon0 = 0.0;
   g.lat0 = 90.0;
   g.lon1 = 0.0;
   g.lat1 = 90.0;
   g.beta = 0.0;
   
   return SLang_push_cstruct (&g, Sphere_Projection_Struct);
}

static void sphere_free (void *g)
{
   if (g == NULL)
     return;
   
   SLang_free_cstruct (g, Sphere_Projection_Struct);
   SLfree ((char *)g);
}

static void *sphere_pop (char *name)
{
   Sphere_Projection_Type *g;
   
   (void) name;

   if (NULL == (g = (Sphere_Projection_Type *)SLmalloc (sizeof (Sphere_Projection_Type))))
     return NULL;

   if (-1 == SLang_pop_cstruct (g, Sphere_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }
   if (-1 == init_sphere_xform (&g->sphere_xform, RADIANS(g->lon0), RADIANS(g->lat0),
				RADIANS(g->lon1), RADIANS(g->lat1),
				RADIANS(g->beta)))
     {
	SLfree ((char *)g);
	return NULL;
     }
   return g;
}

static int sphere_project_f (void *vg,
			     float *x, float *y,
			     float *newx, float *newy,
			     unsigned int num)
{
   Sphere_Projection_Type *g = (Sphere_Projection_Type *)vg;
   Sphere_Rotate_Type *s;
   unsigned int i;
   
   s = &g->sphere_xform;
   for (i = 0; i < num; i++)
     {
	double lon1, lat1;
	rotate_sphere (s, 1, RADIANS(x[i]), RADIANS(y[i]), &lon1, &lat1);
	newx[i] = (float)DEGREES(lon1);
	newy[i] = (float)DEGREES(lat1);
     }
   return 0;
}

static int sphere_deproject_f (void *vg,
			       float *x, float *y,
			       float *newx, float *newy,
			       unsigned int num)
{
   Sphere_Projection_Type *g = (Sphere_Projection_Type *)vg;
   Sphere_Rotate_Type *s;
   unsigned int i;
   
   s = &g->sphere_xform;
   for (i = 0; i < num; i++)
     {
	double lon1, lat1;
	rotate_sphere (s, -1, RADIANS(x[i]), RADIANS(y[i]), &lon1, &lat1);
	newx[i] = (float)DEGREES(lon1);
	newy[i] = (float)DEGREES(lat1);
     }
   return 0;
}

static int sphere_project_d (void *vg,
			     double *x, double *y,
			     double *newx, double *newy,
			     unsigned int num)
{
   Sphere_Projection_Type *g = (Sphere_Projection_Type *)vg;
   Sphere_Rotate_Type *s;
   unsigned int i;
   
   s = &g->sphere_xform;
   for (i = 0; i < num; i++)
     {
	double lon1, lat1;
	rotate_sphere (s, 1, RADIANS(x[i]), RADIANS(y[i]), &lon1, &lat1);
	newx[i] = DEGREES(lon1);
	newy[i] = DEGREES(lat1);
     }
   return 0;
}

static int sphere_deproject_d (void *vg,
			       double *x, double *y,
			       double *newx, double *newy,
			       unsigned int num)
{
   Sphere_Projection_Type *g = (Sphere_Projection_Type *)vg;
   Sphere_Rotate_Type *s;
   unsigned int i;
   
   s = &g->sphere_xform;
   for (i = 0; i < num; i++)
     {
	double lon1, lat1;
	rotate_sphere (s, -1, RADIANS(x[i]), RADIANS(y[i]), &lon1, &lat1);
	newx[i] = DEGREES(lon1);
	newy[i] = DEGREES(lat1);
     }
   return 0;
}
	
/*}}}*/

/*{{{ Gnomic Projection */

typedef struct
{
   REQUIRED_FIELDS;
   double tx0, ty0, txscale, tyscale;
   double tx_hat[3], ty_hat[3], tz_hat[3];
}
Gnomic_Projection_Type;

static SLang_CStruct_Field_Type Gnomic_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Gnomic_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Gnomic_Projection_Type, lon0, "lon0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Gnomic_Projection_Type, lat0, "lat0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Gnomic_Projection_Type, tx0, "x0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Gnomic_Projection_Type, ty0, "y0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Gnomic_Projection_Type, txscale, "xscale", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Gnomic_Projection_Type, tyscale, "yscale", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int gnomic_new (char *name)
{
   Gnomic_Projection_Type g;
   
   g.name = name;
   g.lon0 = 0.0;
   g.lat0 = 90.0;
   g.tx0 = g.ty0 = 0;
   g.txscale = 1.0;
   g.tyscale = 1.0;
   
   return SLang_push_cstruct (&g, Gnomic_Projection_Struct);
}

static void *gnomic_pop (char *name)
{
   Gnomic_Projection_Type *g;
   double phi0, theta0;

   if (NULL == (g = (Gnomic_Projection_Type *)SLmalloc (sizeof (Gnomic_Projection_Type))))
     return NULL;
   
   if (-1 == SLang_pop_cstruct (g, Gnomic_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }
   
   if ((0.0 == g->txscale) || (0.0 == g->tyscale))
     {
	SLang_verror (SL_INVALID_PARM, "%s txscale and tyscale must be non-zero", name);
	SLang_free_cstruct (g, Gnomic_Projection_Struct);
	return NULL;
     }

   g->txscale = RADIANS(g->txscale);
   g->tyscale = RADIANS(g->tyscale);
   phi0 = LON_TO_PHI(g->lon0);
   theta0 = LAT_TO_THETA(g->lat0);
   get_tp_basis (theta0, phi0, g->tx_hat, g->ty_hat, g->tz_hat);
   return g;
}

static void gnomic_free (void *g)
{
   if (g == NULL)
     return;
   
   SLang_free_cstruct (g, Gnomic_Projection_Struct);
   SLfree ((char *)g);
}

static int gnomic_project_f (void *vg,
			     float *lon, float *lat,
			     float *tx, float *ty,
			     unsigned int num)
{
   Gnomic_Projection_Type *g = (Gnomic_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   double inv_txscale, inv_tyscale;
   unsigned int i;

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double theta, phi;

	theta = LAT_TO_THETA(lat[i]);
	phi = LON_TO_PHI(lon[i]);
	spherical_to_vector (1.0, theta, phi, v);
	(void) vector_to_tangent_plane (v, tx_hat, ty_hat, tz_hat, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}

static int gnomic_project_d (void *vg,
			     double *lon, double *lat,
			     double *tx, double *ty,
			     unsigned int num)
{
   Gnomic_Projection_Type *g = (Gnomic_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   double inv_txscale, inv_tyscale;
   unsigned int i;

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double theta, phi;

	theta = LAT_TO_THETA(lat[i]);
	phi = LON_TO_PHI(lon[i]);
	spherical_to_vector (1.0, theta, phi, v);
	(void) vector_to_tangent_plane (v, tx_hat, ty_hat, tz_hat, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}

static int gnomic_deproject_f (void *vg,
			       float *tx, float *ty,
			       float *lon, float *lat,
			       unsigned int num)
{
   Gnomic_Projection_Type *g = (Gnomic_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   unsigned int i;
   
   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double r, theta, phi;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	tangent_plane_to_vector (tx_hat, ty_hat, tz_hat, tx_i, ty_i, v);
	vector_to_spherical (v[0], v[1], v[2], &r, &theta, &phi);
	lon[i] = PHI_TO_LON(phi);
	lat[i] = THETA_TO_LAT(theta);
     }
   return 0;
}

static int gnomic_deproject_d (void *vg,
			       double *tx, double *ty,
			       double *lon, double *lat,
			       unsigned int num)
{
   Gnomic_Projection_Type *g = (Gnomic_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   unsigned int i;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double r, theta, phi;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	tangent_plane_to_vector (tx_hat, ty_hat, tz_hat, tx_i, ty_i, v);
	vector_to_spherical (v[0], v[1], v[2], &r, &theta, &phi);
	lon[i] = PHI_TO_LON(phi);
	lat[i] = THETA_TO_LAT(theta);
     }
   return 0;
}

static int gnomic_reproject_f (void *vga, float *txa, float *tya, 
			       void *vgb, float *txb, float *tyb,
			       unsigned int num)
{
   Gnomic_Projection_Type *ga = (Gnomic_Projection_Type *)vga;
   Gnomic_Projection_Type *gb = (Gnomic_Projection_Type *)vgb;
   double txa0 = ga->tx0, tya0 = ga->ty0;
   double txb0 = gb->tx0, tyb0 = gb->ty0;
   double *txa_hat = ga->tx_hat, *tya_hat = ga->ty_hat, *tza_hat = ga->tz_hat;
   double *txb_hat = gb->tx_hat, *tyb_hat = gb->ty_hat, *tzb_hat = gb->tz_hat;
   double txascale, tyascale, inv_txbscale, inv_tybscale;
   unsigned int i;

   txascale = ga->txscale;
   tyascale = ga->tyscale;
   inv_txbscale = 1.0/gb->txscale;
   inv_tybscale = 1.0/gb->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];

	tx_i = (txa[i] - txa0) * txascale;
	ty_i = (tya[i] - tya0) * tyascale;

	tangent_plane_to_vector (txa_hat, tya_hat, tza_hat, tx_i, ty_i, v);
	(void) vector_to_tangent_plane (v, txb_hat, tyb_hat, tzb_hat, &tx_i, &ty_i);
	txb[i] = txb0 + tx_i*inv_txbscale;
	tyb[i] = tyb0 + ty_i*inv_tybscale;
     }
   return 0;
}

static int gnomic_reproject_d (void *vga, double *txa, double *tya, 
			       void *vgb, double *txb, double *tyb,
			       unsigned int num)
{
   Gnomic_Projection_Type *ga = (Gnomic_Projection_Type *)vga;
   Gnomic_Projection_Type *gb = (Gnomic_Projection_Type *)vgb;
   double txa0 = ga->tx0, tya0 = ga->ty0;
   double txb0 = gb->tx0, tyb0 = gb->ty0;
   double *txa_hat = ga->tx_hat, *tya_hat = ga->ty_hat, *tza_hat = ga->tz_hat;
   double *txb_hat = gb->tx_hat, *tyb_hat = gb->ty_hat, *tzb_hat = gb->tz_hat;
   double txascale, tyascale, inv_txbscale, inv_tybscale;
   unsigned int i;

   txascale = ga->txscale;
   tyascale = ga->tyscale;
   inv_txbscale = 1.0/gb->txscale;
   inv_tybscale = 1.0/gb->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];

	tx_i = (txa[i] - txa0) * txascale;
	ty_i = (tya[i] - tya0) * tyascale;

	tangent_plane_to_vector (txa_hat, tya_hat, tza_hat, tx_i, ty_i, v);
	(void) vector_to_tangent_plane (v, txb_hat, tyb_hat, tzb_hat, &tx_i, &ty_i);
	txb[i] = txb0 + tx_i*inv_txbscale;
	tyb[i] = tyb0 + ty_i*inv_tybscale;
     }
   return 0;
}

	
/*}}}*/

/*{{{ Ortho Projection */

typedef struct
{
   REQUIRED_FIELDS;
   double tx0, ty0, txscale, tyscale;
   double tx_hat[3], ty_hat[3], tz_hat[3];
}
Ortho_Projection_Type;

static SLang_CStruct_Field_Type Ortho_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Ortho_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Ortho_Projection_Type, lon0, "lon0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Ortho_Projection_Type, lat0, "lat0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Ortho_Projection_Type, tx0, "x0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Ortho_Projection_Type, ty0, "y0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Ortho_Projection_Type, txscale, "xscale", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Ortho_Projection_Type, tyscale, "yscale", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int ortho_new (char *name)
{
   Ortho_Projection_Type g;
   
   g.name = name;
   g.lon0 = 0.0;
   g.lat0 = 90.0;
   g.tx0 = g.ty0 = 0;
   g.txscale = 1.0;
   g.tyscale = 1.0;
   
   return SLang_push_cstruct (&g, Ortho_Projection_Struct);
}

static void *ortho_pop (char *name)
{
   Ortho_Projection_Type *g;
   double phi0, theta0;

   if (NULL == (g = (Ortho_Projection_Type *)SLmalloc (sizeof (Ortho_Projection_Type))))
     return NULL;
   
   if (-1 == SLang_pop_cstruct (g, Ortho_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }
   
   if ((0.0 == g->txscale) || (0.0 == g->tyscale))
     {
	SLang_verror (SL_INVALID_PARM, "%s txscale and tyscale must be non-zero", name);
	SLang_free_cstruct (g, Ortho_Projection_Struct);
	return NULL;
     }

   g->txscale = RADIANS(g->txscale);
   g->tyscale = RADIANS(g->tyscale);
   phi0 = LON_TO_PHI(g->lon0);
   theta0 = LAT_TO_THETA(g->lat0);
   get_tp_basis (theta0, phi0, g->tx_hat, g->ty_hat, g->tz_hat);
   return g;
}

static void ortho_free (void *g)
{
   if (g == NULL)
     return;
   
   SLang_free_cstruct (g, Ortho_Projection_Struct);
   SLfree ((char *)g);
}

static int ortho_project_f (void *vg,
			    float *lon, float *lat,
			    float *tx, float *ty,
			    unsigned int num)
{
   Ortho_Projection_Type *g = (Ortho_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   double inv_txscale, inv_tyscale;
   unsigned int i;

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double theta, phi;

	theta = LAT_TO_THETA(lat[i]);
	phi = LON_TO_PHI(lon[i]);
	spherical_to_vector (1.0, theta, phi, v);
	(void) ortho_project_vector (v, tx_hat, ty_hat, tz_hat, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}

static int ortho_project_d (void *vg,
			     double *lon, double *lat,
			     double *tx, double *ty,
			     unsigned int num)
{
   Ortho_Projection_Type *g = (Ortho_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   double inv_txscale, inv_tyscale;
   unsigned int i;

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double theta, phi;

	theta = LAT_TO_THETA(lat[i]);
	phi = LON_TO_PHI(lon[i]);
	spherical_to_vector (1.0, theta, phi, v);
	(void) ortho_project_vector (v, tx_hat, ty_hat, tz_hat, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}

static int ortho_deproject_f (void *vg,
			       float *tx, float *ty,
			       float *lon, float *lat,
			       unsigned int num)
{
   Ortho_Projection_Type *g = (Ortho_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   unsigned int i;
   
   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double r, theta, phi;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) ortho_deproject_vector (tx_hat, ty_hat, tz_hat, tx_i, ty_i, v);
	vector_to_spherical (v[0], v[1], v[2], &r, &theta, &phi);
	lon[i] = PHI_TO_LON(phi);
	lat[i] = THETA_TO_LAT(theta);
     }
   return 0;
}

static int ortho_deproject_d (void *vg,
			       double *tx, double *ty,
			       double *lon, double *lat,
			       unsigned int num)
{
   Ortho_Projection_Type *g = (Ortho_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat, *tz_hat = g->tz_hat;
   unsigned int i;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double v[3];
	double r, theta, phi;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) ortho_deproject_vector (tx_hat, ty_hat, tz_hat, tx_i, ty_i, v);
	vector_to_spherical (v[0], v[1], v[2], &r, &theta, &phi);
	lon[i] = PHI_TO_LON(phi);
	lat[i] = THETA_TO_LAT(theta);
     }
   return 0;
}
	
/*}}}*/

/*{{{ Generic Plane Projection */

typedef struct
{
   REQUIRED_FIELDS;
   double tx0, ty0, txscale, tyscale;
   double p_lon, p_lat, p_len;
   double tx_hat[3], ty_hat[3];
   double s0[3];		       /* pos of tangent plane on sphere -- also its normal */
   double phat[3];		       /* unit vector to observer position */
}
Plane_Projection_Type;

static SLang_CStruct_Field_Type Plane_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, lon0, "lon0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, lat0, "lat0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, tx0, "x0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, ty0, "y0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, txscale, "xscale", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, tyscale, "yscale", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, p_lon, "p_lon", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, p_lat, "p_lat", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Plane_Projection_Type, p_len, "p_len", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int plane_new (char *name)
{
   Plane_Projection_Type g;
   
   g.name = name;
   g.lon0 = 0.0;
   g.lat0 = 90.0;
   g.tx0 = g.ty0 = 0;
   g.txscale = 1.0;
   g.tyscale = 1.0;
   
   g.p_lon = 0.0;
   g.p_lat = -90.0;
   g.p_len = 0.0;		       /* gnomic */

   return SLang_push_cstruct (&g, Plane_Projection_Struct);
}

static void *plane_pop (char *name)
{
   Plane_Projection_Type *g;

   if (NULL == (g = (Plane_Projection_Type *)SLmalloc (sizeof (Plane_Projection_Type))))
     return NULL;
   
   if (-1 == SLang_pop_cstruct (g, Plane_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }
   
   if ((0.0 == g->txscale) || (0.0 == g->tyscale))
     {
	SLang_verror (SL_INVALID_PARM, "%s txscale and tyscale must be non-zero", name);
	SLang_free_cstruct (g, Plane_Projection_Struct);
	return NULL;
     }

   g->txscale = RADIANS(g->txscale);
   g->tyscale = RADIANS(g->tyscale);

   g->p_len = fabs (g->p_len);
   (void) compute_plane_vectors (LAT_TO_THETA(g->lat0), LON_TO_PHI(g->lon0),
				 g->tx_hat, g->ty_hat, g->s0);
   spherical_to_vector (1.0, LAT_TO_THETA(g->p_lat), LON_TO_PHI(g->p_lon), g->phat);
   return g;
}

static void plane_free (void *g)
{
   if (g == NULL)
     return;
   
   SLang_free_cstruct (g, Plane_Projection_Struct);
   SLfree ((char *)g);
}

static int plane_project_f (void *vg,
			    float *lon, float *lat,
			    float *tx, float *ty,
			    unsigned int num)
{
   Plane_Projection_Type *g = (Plane_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat;
   double *phat = g->phat, *s0 = g->s0;
   double plen = g->p_len;
   double inv_txscale, inv_tyscale;
   unsigned int i;
   
   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double s[3];
	double theta, phi;

	theta = LAT_TO_THETA(lat[i]);
	phi = LON_TO_PHI(lon[i]);
	spherical_to_vector (1.0, theta, phi, s);
	(void) project_to_plane_coords (s, s0, plen, phat, tx_hat, ty_hat, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}

static int plane_project_d (void *vg,
			    double *lon, double *lat,
			    double *tx, double *ty,
			    unsigned int num)
{
   Plane_Projection_Type *g = (Plane_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat;
   double *phat = g->phat, *s0 = g->s0;
   double plen = g->p_len;
   double inv_txscale, inv_tyscale;
   unsigned int i;

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double s[3];
	double theta, phi;

	theta = LAT_TO_THETA(lat[i]);
	phi = LON_TO_PHI(lon[i]);
	spherical_to_vector (1.0, theta, phi, s);
	(void) project_to_plane_coords (s, s0, plen, phat, tx_hat, ty_hat, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}


static int plane_deproject_f (void *vg,
			      float *tx, float *ty,
			      float *lon, float *lat,
			      unsigned int num)
{
   Plane_Projection_Type *g = (Plane_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat;
   double *phat = g->phat, *s0 = g->s0;
   double plen = g->p_len;
   unsigned int i;
   
   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double s[3];
	double r, theta, phi;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) project_from_plane_coords (tx_i, ty_i, s0, plen, phat, tx_hat, ty_hat, s);
	vector_to_spherical (s[0], s[1], s[2], &r, &theta, &phi);
	lon[i] = PHI_TO_LON(phi);
	lat[i] = THETA_TO_LAT(theta);
     }
   return 0;
}

static int plane_deproject_d (void *vg,
			      double *tx, double *ty,
			      double *lon, double *lat,
			      unsigned int num)
{
   Plane_Projection_Type *g = (Plane_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double *tx_hat = g->tx_hat, *ty_hat = g->ty_hat;
   double *phat = g->phat, *s0 = g->s0;
   double plen = g->p_len;
   unsigned int i;
   
   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double s[3];
	double r, theta, phi;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) project_from_plane_coords (tx_i, ty_i, s0, plen, phat, tx_hat, ty_hat, s);

	vector_to_spherical (s[0], s[1], s[2], &r, &theta, &phi);
	lon[i] = PHI_TO_LON(phi);
	lat[i] = THETA_TO_LAT(theta);
     }
   return 0;
}
	
/*}}}*/

/*{{{ Stereographic */

typedef Plane_Projection_Type Stereo_Projection_Type;

static SLang_CStruct_Field_Type Stereo_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Stereo_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Stereo_Projection_Type, lon0, "lon0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Stereo_Projection_Type, lat0, "lat0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Stereo_Projection_Type, tx0, "x0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Stereo_Projection_Type, ty0, "y0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Stereo_Projection_Type, txscale, "xscale", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Stereo_Projection_Type, tyscale, "yscale", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int stereo_new (char *name)
{
   Stereo_Projection_Type g;
   
   g.name = name;
   g.lon0 = 0.0;
   g.lat0 = 90.0;
   g.tx0 = g.ty0 = 0;
   g.txscale = 1.0;
   g.tyscale = 1.0;
   
   return SLang_push_cstruct (&g, Stereo_Projection_Struct);
}

static void *stereo_pop (char *name)
{
   Stereo_Projection_Type *g;

   if (NULL == (g = (Stereo_Projection_Type *)SLmalloc (sizeof (Stereo_Projection_Type))))
     return NULL;
   
   if (-1 == SLang_pop_cstruct (g, Stereo_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }
   
   if ((0.0 == g->txscale) || (0.0 == g->tyscale))
     {
	SLang_verror (SL_INVALID_PARM, "%s txscale and tyscale must be non-zero", name);
	SLang_free_cstruct (g, Stereo_Projection_Struct);
	return NULL;
     }

   g->txscale = RADIANS(g->txscale);
   g->tyscale = RADIANS(g->tyscale);

   normalize_lon_lat_degrees (&g->lon0, &g->lat0);

   g->p_lon = g->lon0;
   g->p_lat = -g->lat0;
   g->p_len = 1.0;

   (void) compute_plane_vectors (LAT_TO_THETA(g->lat0), LON_TO_PHI(g->lon0),
				 g->tx_hat, g->ty_hat, g->s0);
   spherical_to_vector (1.0, LAT_TO_THETA(g->p_lat), LON_TO_PHI(g->p_lon), g->phat);
   return g;
}

static void stereo_free (void *g)
{
   if (g == NULL)
     return;

   SLang_free_cstruct (g, Stereo_Projection_Struct);
   SLfree ((char *)g);
}
/*}}}*/

/*{{{ Lambert */

typedef struct
{
   REQUIRED_FIELDS;
   double tx0, ty0;
   double txscale, tyscale;
   double cos_lat0, sin_lat0;
}
Lambert_Projection_Type;

static SLang_CStruct_Field_Type Lambert_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Lambert_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Lambert_Projection_Type, lon0, "lon0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Lambert_Projection_Type, lat0, "lat0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Lambert_Projection_Type, tx0, "x0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Lambert_Projection_Type, ty0, "y0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Lambert_Projection_Type, txscale, "xscale", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Lambert_Projection_Type, tyscale, "yscale", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int lambert_new (char *name)
{
   Lambert_Projection_Type g;
   
   g.name = name;
   g.lon0 = 0.0;
   g.lat0 = 90.0;
   g.tx0 = g.ty0 = 0;
   g.txscale = 1.0;
   g.tyscale = 1.0;
   
   return SLang_push_cstruct (&g, Lambert_Projection_Struct);
}

static void *lambert_pop (char *name)
{
   Lambert_Projection_Type *g;
   double lat0;

   if (NULL == (g = (Lambert_Projection_Type *)SLmalloc (sizeof (Lambert_Projection_Type))))
     return NULL;
   
   if (-1 == SLang_pop_cstruct (g, Lambert_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }
   
   if ((0.0 == g->txscale) || (0.0 == g->tyscale))
     {
	SLang_verror (SL_INVALID_PARM, "%s txscale and tyscale must be non-zero", name);
	SLang_free_cstruct (g, Lambert_Projection_Struct);
	return NULL;
     }

   g->txscale = RADIANS(g->txscale);
   g->tyscale = RADIANS(g->tyscale);
   normalize_lon_lat_degrees (&g->lon0, &g->lat0);
   lat0 = RADIANS(g->lat0);
   g->cos_lat0 = cos (lat0);
   g->sin_lat0 = sin (lat0);
   return g;
}

static void lambert_free (void *g)
{
   if (g == NULL)
     return;

   SLang_free_cstruct (g, Lambert_Projection_Struct);
   SLfree ((char *)g);
}

static int lambert_lonlat_to_xy (double lon, double lat, double lon0,
				 double cos_lat0, double sin_lat0,
				 double *xp, double *yp)
{
   double r, x, y;
   double cos_lat = cos (lat);
   double sin_lat = sin (lat);
   double cos_lon, sin_lon;

   lon = fmod (lon, 2*PI);	       /* lon could be huge */
   lon -= lon0;

   cos_lon = cos (lon);
   sin_lon = sin (lon);

   r = 1.0 + sin_lat0*sin_lat + cos_lat0*cos_lat*cos_lon;
   if (r <= 1e-9)
     {
#if 0
	x = 2.0*sin_lon;
	y = -2.0*cos_lon;
	if (sin_lat0 < 0)
	  y = -y;
	*xp = x;
	*yp = y;
	return 0;
#endif
     }
   r = sqrt (2.0/r);
   
   x = r * (cos_lat*sin_lon);
   y = r * (cos_lat0 * sin_lat - sin_lat0*cos_lat*cos_lon);
   
   *xp = x;
   *yp = y;
   return 0;
}

static int lambert_xy_to_lonlat (double x, double y, 
				 double lon0,
				 double cos_lat0, double sin_lat0,
				 double *lonp, double *latp)
{
   double r, y_r, x_r;
   double theta, c, s;
   double arg, lon, lat;

   r = hypot (x, y);
   if (r > 2.0)
     {
	*latp = *lonp = Maplib_NaN;
	return -1;
     }
     
   theta = 2.0 * asin(0.5*r);
   c = cos (theta);
   s = sin (theta);

   if (r == 0.0)
     y_r = x_r = 1.0;
   else 
     {
	y_r = y/r;
	x_r = x/r;
     }

   arg = c * sin_lat0 + y_r * s * cos_lat0;
   if (fabs(arg) > 1.0)
     {
	*lonp = *latp = Maplib_NaN;
	return -1;
     }
   lat = asin (arg);
   lon = lon0 + atan2 (s*x_r, c*cos_lat0 - y_r*s*sin_lat0);
   
   lon = fmod (lon, 2*PI);
   if (lon < -PI)
     lon += 2*PI;

   *lonp = lon;
   *latp = lat;
   return 0;
}

static int lambert_project_d (void *vg,
			      double *lon, double *lat,
			      double *tx, double *ty,
			      unsigned int num)
{
   Lambert_Projection_Type *g = (Lambert_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double inv_txscale, inv_tyscale;
   double cos_lat0, sin_lat0, lon0;
   unsigned int i;
   
   lon0 = RADIANS(g->lon0);
   cos_lat0 = g->cos_lat0;
   sin_lat0 = g->sin_lat0;

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	lon_i = RADIANS(lon[i]);
	lat_i = RADIANS(lat[i]);
	(void) lambert_lonlat_to_xy (lon_i, lat_i, lon0, cos_lat0, sin_lat0, 
				     &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}


static int lambert_deproject_d (void *vg,
				double *tx, double *ty,
				double *lon, double *lat,
				unsigned int num)
{
   Lambert_Projection_Type *g = (Lambert_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double cos_lat0, sin_lat0, lon0;
   unsigned int i;
   
   lon0 = RADIANS(g->lon0);
   cos_lat0 = g->cos_lat0;
   sin_lat0 = g->sin_lat0;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) lambert_xy_to_lonlat (tx_i, ty_i, lon0, cos_lat0, sin_lat0, 
				      &lon_i, &lat_i);
	lon[i] = DEGREES(lon_i);
	lat[i] = DEGREES(lat_i);
     }
   return 0;
}

static int lambert_project_f (void *vg,
			      float *lon, float *lat,
			      float *tx, float *ty,
			      unsigned int num)
{
   Lambert_Projection_Type *g = (Lambert_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double inv_txscale, inv_tyscale;
   double cos_lat0, sin_lat0, lon0;
   unsigned int i;
   
   lon0 = RADIANS(g->lon0);
   cos_lat0 = g->cos_lat0;
   sin_lat0 = g->sin_lat0;

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	lon_i = RADIANS(lon[i]);
	lat_i = RADIANS(lat[i]);
	(void) lambert_lonlat_to_xy (lon_i, lat_i, lon0, cos_lat0, sin_lat0, 
				     &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}


static int lambert_deproject_f (void *vg,
				float *tx, float *ty,
				float *lon, float *lat,
				unsigned int num)
{
   Lambert_Projection_Type *g = (Lambert_Projection_Type *)vg;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   double cos_lat0, sin_lat0, lon0;
   unsigned int i;
   
   lon0 = RADIANS(g->lon0);
   cos_lat0 = g->cos_lat0;
   sin_lat0 = g->sin_lat0;

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) lambert_xy_to_lonlat (tx_i, ty_i, lon0, cos_lat0, sin_lat0, 
				      &lon_i, &lat_i);
	lon[i] = DEGREES(lon_i);
	lat[i] = DEGREES(lat_i);
     }
   return 0;
}


/*}}}*/

/*{{{ Generic */

typedef struct
{
   REQUIRED_FIELDS;  /* name, lon0, lat0 (degrees) */
   double beta;
   double tx0, ty0;
   double txscale, tyscale;
   double ref_lon, ref_lat;	       /* sphere will be rotated to this position */
   int apply_sphere_rotation;
   Sphere_Rotate_Type sphere_xform;
   void *client_data;
   struct Generic_Info_Type *info;
}
Generic_Projection_Type;

typedef struct Generic_Info_Type
{
   char *name;
   int (*init)(Generic_Projection_Type *);
   void (*deinit)(void *);
   int (*project)(void *, double, double, double *, double *);
   int (*deproject)(void *, double, double, double *, double *);
   double lon0, lat0;		       /* degrees */
}
Generic_Info_Type;

/*{{{ Hammer */

static int hammer_project (void *cd, double lon, double lat, double *xp, double *yp)
{
   double x;
   int ret = lambert_lonlat_to_xy (0.5*lon, lat, 0.0, 1.0, 0.0, &x, yp);
   (void) cd;
   *xp = 2.0 * x;
   return ret;
}

static int hammer_deproject (void *cd, double x, double y, double *lonp, double *latp)
{
   double lon;
   int ret = lambert_xy_to_lonlat (0.5*x, y, 0.0, 1.0, 0.0, &lon, latp);
   (void) cd;
   *lonp = 2.0 * lon;
   return ret;
}

/*}}}*/

/*{{{ Mercator */

static int mercator_project (void *cd, double lon, double lat, double *xp, double *yp)
{
   (void) cd;
   *xp = lon;
   *yp = log (tan (PI/4 + 0.5*lat));
   return 0;
}

static int mercator_deproject (void *cd, double x, double y, double *lonp, double *latp)
{
   (void) cd;
   *lonp = x;
   *latp = 2.0*(atan(exp(y)) - (PI/4));
   return 0;
}

/*}}}*/

/*{{{ Azimuthal Equidist (azeqdist) */

static int azeqdist_project (void *cd, double lon, double lat, double *xp, double *yp)
{
   (void) cd;
   lat = PI/2 - lat;
   *xp = lat * sin(lon);
   *yp = -lat * cos(lon);
   return 0;
}

static int azeqdist_deproject (void *cd, double x, double y, double *lonp, double *latp)
{
   (void) cd;
   
   *lonp = atan2(x,-y);
   *latp = PI/2 - hypot (x,y);
   return 0;
}

/*}}}*/

/*{{{ Bonne */

/* The equations here are taken from Snyder 1987 */
typedef struct
{
   double lon0, lat0;		       /* radians */
   double cot_lat0;
   double cot_lat0_plus_lat0;
}
Bonne_Type;

static void bonne_deinit (void *p)
{
   if (p != NULL)
     SLfree ((char *)p);
}

static int bonne_init (Generic_Projection_Type *g)
{
   Bonne_Type *b;

   if (NULL == (b = (Bonne_Type *)SLmalloc (sizeof (Bonne_Type))))
     return -1;
   g->client_data = b;
   /* Transformation is oblique */
   g->ref_lon = g->lon0;
   g->ref_lat = g->lat0;
   
   b->lon0 = RADIANS(g->lon0);
   b->lat0 = RADIANS(g->lat0);
   b->cot_lat0 = 1.0/tan(b->lat0);
   b->cot_lat0_plus_lat0 = b->cot_lat0 + b->lat0;
   return 0;
}

static int bonne_project (void *cd, double lon, double lat, double *xp, double *yp)
{
   double rho, e, dlon;
   Bonne_Type *b = (Bonne_Type *)cd;
   
   if ((lat == b->lat0)
       && (lat == PI/2))
     {
	*xp = *yp = 0.0;
	return 0;
     }

   rho = (b->cot_lat0_plus_lat0 - lat);
   dlon = lon - b->lon0;
   if ((dlon < -PI) || (dlon > PI))
     dlon = normalize_to_range (dlon, PI);
   
   e = dlon*cos(lat)/rho;
   *xp = rho*sin(e);
   *yp = b->cot_lat0 - rho*cos(e);
   return 0;
}

static int bonne_deproject (void *cd, double x, double y, double *lonp, double *latp)
{
   double rho, dy;
   Bonne_Type *b = (Bonne_Type *)cd;
   double lon, lat;

   dy = b->cot_lat0 - y;
   rho = hypot (x, dy);
   if (b->cot_lat0 < 0)
     rho = -rho;

   lat = b->cot_lat0_plus_lat0 - rho;
   if (b->lat0 < 0)
     {
	x = -x;
	dy = -dy;
     }
   lon = b->lon0 + rho * atan2 (x, dy)/cos(lat);
   if ((lon < -PI) || (lon > PI))
     lon = normalize_to_range (lon, PI);
   
   *lonp = lon;
   *latp = lat;
   return 0;
}

/*}}}*/

/*{{{ Sinusoidal */

/* The equations here are taken from Snyder 1987 */

typedef struct
{
   double lon0;		       /* radians */
}
Sinusoidal_Type;

static void sinusoidal_deinit (void *p)
{
   if (p != NULL)
     SLfree ((char *)p);
}

static int sinusoidal_init (Generic_Projection_Type *g)
{
   Sinusoidal_Type *b;

   if (NULL == (b = (Sinusoidal_Type *)SLmalloc (sizeof (Sinusoidal_Type))))
     return -1;
   g->client_data = b;

   /* Transformation is oblique */
   g->ref_lon = g->lon0;
   g->ref_lat = g->lat0;
   
   b->lon0 = RADIANS(g->lon0);
   return 0;
}

static int sinusoidal_project (void *cd, double lon, double lat, double *xp, double *yp)
{
   Sinusoidal_Type *b = (Sinusoidal_Type *)cd;

   lon = lon - b->lon0;
   if ((lon < -PI) || (lon > PI))
     lon = normalize_to_range (lon, PI);

   *xp = lon * cos (lat);
   *yp = lat;
   return 0;
}

static int sinusoidal_deproject (void *cd, double x, double y, double *lonp, double *latp)
{
   Sinusoidal_Type *b = (Sinusoidal_Type *)cd;
   double lon;

   *latp = y;
   lon = b->lon0 + x/cos(y);

   if ((lon < -PI) || (lon > PI))
     lon = normalize_to_range (lon, PI);
   
   *lonp = lon;
   return 0;
}

/*}}}*/

static Generic_Info_Type Generic_Info_Table[] = 
{
     {"hammer", NULL, NULL, hammer_project, hammer_deproject, 0.0, 0.0},
     {"mercator", NULL, NULL, mercator_project, mercator_deproject, 0.0, 0.0},
     {"azeqdist", NULL, NULL, azeqdist_project, azeqdist_deproject, 0.0, 90.0},
     {"bonne", bonne_init, bonne_deinit, bonne_project, bonne_deproject, 0.0, 45.0},
     {"sinusoidal", sinusoidal_init, sinusoidal_deinit, sinusoidal_project, sinusoidal_deproject, 0.0, 0.0},
     {NULL, NULL, NULL, NULL, NULL, 0.0, 0.0}
};

static Generic_Info_Type *lookup_generic_xform (char *name)
{
   Generic_Info_Type *info;
   
   info = Generic_Info_Table;
   while (info->name != NULL)
     {
	if (0 == strcmp (name, info->name))
	  return info;
	info++;
     }
   SLang_verror (SL_INTERNAL_ERROR, "Generic projection %s not found", name);
   return NULL;
}

static SLang_CStruct_Field_Type Generic_Projection_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, lon0, "lon0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, lat0, "lat0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, beta, "beta", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, tx0, "x0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, ty0, "y0", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, txscale, "xscale", SLANG_DOUBLE_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Generic_Projection_Type, tyscale, "yscale", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int generic_new (char *name)
{
   Generic_Projection_Type g;
   Generic_Info_Type *info;

   if (NULL == (info = lookup_generic_xform (name)))
     return -1;

   g.name = name;
   g.lon0 = info->lon0;
   g.lat0 = info->lat0;
   g.beta = 0.0;
   g.tx0 = g.ty0 = 0;
   g.txscale = 1.0;
   g.tyscale = 1.0;

   return SLang_push_cstruct (&g, Generic_Projection_Struct);
}

static void generic_free (void *vg)
{
   Generic_Projection_Type *g = (Generic_Projection_Type *)vg;
   if (g == NULL)
     return;
   
   if ((g->info != NULL) && (g->info->deinit != NULL))
     (*g->info->deinit)(g->client_data);

   SLang_free_cstruct (g, Generic_Projection_Struct);
   SLfree ((char *)g);
}

static void *generic_pop (char *name)
{
   Generic_Projection_Type *g;
   Generic_Info_Type *info;

   if (NULL == (info = lookup_generic_xform (name)))
     return NULL;

   if (NULL == (g = (Generic_Projection_Type *)SLmalloc (sizeof (Generic_Projection_Type))))
     return NULL;
   
   if (-1 == SLang_pop_cstruct (g, Generic_Projection_Struct))
     {
	SLfree ((char *)g);
	return NULL;
     }

   if ((0.0 == g->txscale) || (0.0 == g->tyscale))
     {
	SLang_verror (SL_INVALID_PARM, "%s txscale and tyscale must be non-zero", name);
	SLfree ((char *)g);
	return NULL;
     }

   g->txscale = RADIANS(g->txscale);
   g->tyscale = RADIANS(g->tyscale);
   normalize_lon_lat_degrees (&g->lon0, &g->lat0);
   g->info = info;
   g->ref_lon = info->lon0;
   g->ref_lat = info->lat0;
   g->client_data = NULL;

   if ((info->init != NULL) && (-1 == (*info->init)(g)))
     {
	SLfree ((char *)g);
	return NULL;
     }
   /* From here on, we can call generic_free to free g */

   g->apply_sphere_rotation = 0;
   if ((g->lon0 != g->ref_lon) || (g->lat0 != g->ref_lat) || (g->beta != 0.0))
     {
	if (-1 == init_sphere_xform (&g->sphere_xform, RADIANS(g->lon0), RADIANS(g->lat0), RADIANS(g->ref_lon), RADIANS(g->ref_lat), RADIANS(g->beta)))
	  {
	     generic_free (g);
	     return NULL;
	  }
	g->apply_sphere_rotation = 1;
     }
   return g;
}

static int generic_project_d (void *vg,
			      double *lon, double *lat,
			      double *tx, double *ty,
			      unsigned int num)
{
   Generic_Projection_Type *g = (Generic_Projection_Type *)vg;
   int (*project)(void *, double, double, double *, double *);
   void *client_data;
   double tx0 = g->tx0, ty0 = g->ty0;
   double inv_txscale, inv_tyscale;
   unsigned int i;

   client_data = g->client_data;
   project = g->info->project;
   if (project == NULL)
     {
	SLang_verror (SL_NOT_IMPLEMENTED, "Projection %s does not support projecting", g->name);
	return -1;
     }
   
   if (g->apply_sphere_rotation)
     {
	Sphere_Rotate_Type *s = &g->sphere_xform;
	for (i = 0; i < num; i++)
	  {
	     double tx1, ty1;
	     rotate_sphere (s, 1, RADIANS(lon[i]), RADIANS(lat[i]), &tx1, &ty1);
	     tx[i] = DEGREES(tx1);
	     ty[i] = DEGREES(ty1);
	  }
	lon = tx;
	lat = ty;
     }

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;
   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	lon_i = RADIANS(lon[i]);
	lat_i = RADIANS(lat[i]);
	(void) (*project)(client_data, lon_i, lat_i, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}


static int generic_deproject_d (void *vg,
				double *tx, double *ty,
				double *lon, double *lat,
				unsigned int num)
{
   Generic_Projection_Type *g = (Generic_Projection_Type *)vg;
   int (*deproject)(void *, double, double, double *, double *);
   void *client_data;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   unsigned int i;

   client_data = g->client_data;
   deproject = g->info->deproject;
   if (deproject == NULL)
     {
	SLang_verror (SL_NOT_IMPLEMENTED, "Projection %s does not support deprojecting", g->name);
	return -1;
     }

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) (*deproject)(client_data, tx_i, ty_i, &lon_i, &lat_i);
	lon[i] = DEGREES(lon_i);
	lat[i] = DEGREES(lat_i);
     }

   if (g->apply_sphere_rotation)
     {
	Sphere_Rotate_Type *s = &g->sphere_xform;
	for (i = 0; i < num; i++)
	  {
	     double lon1, lat1;
	     rotate_sphere (s, -1, RADIANS(lon[i]), RADIANS(lat[i]), &lon1, &lat1);
	     lon[i] = DEGREES(lon1);
	     lat[i] = DEGREES(lat1);
	  }
     }
   return 0;
}

static int generic_project_f (void *vg,
			      float *lon, float *lat,
			      float *tx, float *ty,
			      unsigned int num)
{
   Generic_Projection_Type *g = (Generic_Projection_Type *)vg;
   int (*project)(void *, double, double, double *, double *);
   void *client_data;
   double tx0 = g->tx0, ty0 = g->ty0;
   double inv_txscale, inv_tyscale;
   unsigned int i;

   client_data = g->client_data;
   project = g->info->project;
   if (project == NULL)
     {
	SLang_verror (SL_NOT_IMPLEMENTED, "Projection %s does not support projecting", g->name);
	return -1;
     }
   
   if (g->apply_sphere_rotation)
     {
	Sphere_Rotate_Type *s = &g->sphere_xform;
	for (i = 0; i < num; i++)
	  {
	     double tx1, ty1;
	     rotate_sphere (s, 1, RADIANS(lon[i]), RADIANS(lat[i]), &tx1, &ty1);
	     tx[i] = DEGREES(tx1);
	     ty[i] = DEGREES(ty1);
	  }
	lon = tx;
	lat = ty;
     }

   inv_txscale = 1.0/g->txscale;
   inv_tyscale = 1.0/g->tyscale;
   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	lon_i = RADIANS(lon[i]);
	lat_i = RADIANS(lat[i]);
	(void) (*project)(client_data, lon_i, lat_i, &tx_i, &ty_i);
	tx[i] = tx0 + tx_i*inv_txscale;
	ty[i] = ty0 + ty_i*inv_tyscale;
     }
   return 0;
}


static int generic_deproject_f (void *vg,
				float *tx, float *ty,
				float *lon, float *lat,
				unsigned int num)
{
   Generic_Projection_Type *g = (Generic_Projection_Type *)vg;
   int (*deproject)(void *, double, double, double *, double *);
   void *client_data;
   double tx0 = g->tx0, ty0 = g->ty0;
   double txscale = g->txscale, tyscale = g->tyscale;
   unsigned int i;

   client_data = g->client_data;
   deproject = g->info->deproject;
   if (deproject == NULL)
     {
	SLang_verror (SL_NOT_IMPLEMENTED, "Projection %s does not support deprojecting", g->name);
	return -1;
     }

   for (i = 0; i < num; i++)
     {
	double tx_i, ty_i;
	double lon_i, lat_i;

	tx_i = (tx[i] - tx0) * txscale;
	ty_i = (ty[i] - ty0) * tyscale;
	(void) (*deproject)(client_data, tx_i, ty_i, &lon_i, &lat_i);
	lon[i] = DEGREES(lon_i);
	lat[i] = DEGREES(lat_i);
     }

   if (g->apply_sphere_rotation)
     {
	Sphere_Rotate_Type *s = &g->sphere_xform;
	for (i = 0; i < num; i++)
	  {
	     double lon1, lat1;
	     rotate_sphere (s, -1, RADIANS(lon[i]), RADIANS(lat[i]), &lon1, &lat1);
	     lon[i] = DEGREES(lon1);
	     lat[i] = DEGREES(lat1);
	  }
     }
   return 0;
}

/*}}}*/

/* Intrinsic functions */

typedef struct Projection_Info_Type
{
   char *name;
   int (*new_projection) (char *);
   void *(*pop_projection) (char *);
   void (*free_projection) (void *);
   int (*project_fun_d) (void *, double *, double *, double *, double *, unsigned int);
   int (*deproject_fun_d) (void *, double *, double *, double *, double *, unsigned int);
   int (*reproject_fun_d) (void *, double *, double *, void *, double *, double *, unsigned int);
   int (*project_fun_f) (void *, float *, float *, float *, float *, unsigned int);
   int (*deproject_fun_f) (void *, float *, float *, float *, float *, unsigned int);
   int (*reproject_fun_f) (void *, float *, float *, void *, float *, float *, unsigned int);
}
Projection_Info_Type;

static Projection_Info_Type Projection_List[] = 
{
     {"linear",
	  linear_new, linear_pop, linear_free,
	  linear_project_d, linear_deproject_d, NULL,
	  linear_project_f, linear_deproject_f, NULL,
     },
     {"sphere",
	  sphere_new, sphere_pop, sphere_free,
	  sphere_project_d, sphere_deproject_d, NULL,
	  sphere_project_f, sphere_deproject_f, NULL,
     },
     {"gnomic",
	  gnomic_new, gnomic_pop, gnomic_free,
	  gnomic_project_d, gnomic_deproject_d, gnomic_reproject_d,
	  gnomic_project_f, gnomic_deproject_f, gnomic_reproject_f
     },
     {"plane",
	  plane_new, plane_pop, plane_free,
	  plane_project_d, plane_deproject_d, NULL,
	  plane_project_f, plane_deproject_f, NULL,
     },
     {"ortho",
	  ortho_new, ortho_pop, ortho_free,
	  ortho_project_d, ortho_deproject_d, NULL,
	  ortho_project_f, ortho_deproject_f, NULL,
     },
     {"stereo",
	  stereo_new, stereo_pop, stereo_free,
	  plane_project_d, plane_deproject_d, NULL,
	  plane_project_f, plane_deproject_f, NULL,
     },
     {"lambert",
	  lambert_new, lambert_pop, lambert_free,
	  lambert_project_d, lambert_deproject_d, NULL,
	  lambert_project_f, lambert_deproject_f, NULL,
     },
     {"hammer",
	  generic_new, generic_pop, generic_free,
	  generic_project_d, generic_deproject_d, NULL,
	  generic_project_f, generic_deproject_f, NULL,
     },
     {"mercator",
	  generic_new, generic_pop, generic_free,
	  generic_project_d, generic_deproject_d, NULL,
	  generic_project_f, generic_deproject_f, NULL,
     },
     {"azeqdist",
	  generic_new, generic_pop, generic_free,
	  generic_project_d, generic_deproject_d, NULL,
	  generic_project_f, generic_deproject_f, NULL,
     },
     {"bonne",
	  generic_new, generic_pop, generic_free,
	  generic_project_d, generic_deproject_d, NULL,
	  generic_project_f, generic_deproject_f, NULL,
     },
     {"sinusoidal",
	  generic_new, generic_pop, generic_free,
	  generic_project_d, generic_deproject_d, NULL,
	  generic_project_f, generic_deproject_f, NULL,
     },
     {
	NULL,
	  NULL, NULL, NULL,
 	  NULL, NULL, NULL, 
	  NULL, NULL, NULL
     }
};

static void get_projections (void)
{
   SLang_Array_Type *at;
   Projection_Info_Type *p;
   int i, num;
   char **names;
   
   num = 0;
   p = Projection_List;
   while (p->name != NULL)
     {
	num++;
	p++;
     }
   
   if (NULL == (at = SLang_create_array (SLANG_STRING_TYPE, 1, NULL, &num, 1)))
     return;

   names = (char **) at->data;
   p = Projection_List;
   for (i = 0; i < num; i++)
     {
	if (NULL == (names[i] = SLang_create_slstring (p->name)))
	  {
	     SLang_free_array (at);
	     return;
	  }
	p++;
     }
   (void) SLang_push_array (at, 1);
}


typedef struct
{
   Projection_Info_Type *pinfo;
   void *projection_data;
}
Projection_Type;

static Projection_Info_Type *lookup_projection (char *name)
{
   Projection_Info_Type *p;

   p = Projection_List;
   while (p->name != NULL)
     {
	if ((p->name[0] == name[0])
	    && (0 == strcmp (p->name, name)))
	  return p;
	p++;
     }
   SLang_verror (SL_INVALID_PARM, "Projection %s is unknown or unsupported", name);
   return NULL;
}

static void project_new (char *name)
{
   Projection_Info_Type *pinfo;

   if (NULL == (pinfo = lookup_projection (name)))
     return;
   
   (void) (*pinfo->new_projection)(name);
}

typedef struct
{
   char *name;
}
Any_Projection_Type;

static SLang_CStruct_Field_Type Any_Projection_Struct[] = 
{
   MAKE_CSTRUCT_FIELD(Any_Projection_Type, name, "name", SLANG_STRING_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static int pop_projection (Projection_Info_Type **pinfop,
			   void **projection_datap)
{
   Any_Projection_Type p;
   Projection_Info_Type *pinfo;
   int ret;

   if (-1 == SLdup_n (1))
     return -1;
   
   if (-1 == SLang_pop_cstruct (&p, Any_Projection_Struct))
     return -1;
   
   if ((NULL == (pinfo = lookup_projection (p.name)))
       || (NULL == (*projection_datap = (pinfo->pop_projection (p.name)))))
     ret = -1;
   else
     ret = 0;

   SLang_free_cstruct (&p, Any_Projection_Struct);
   *pinfop = pinfo;
   return ret;
}


static int do_reprojection (Projection_Info_Type *pinfo_to, void *projection_data_to,
			    Projection_Info_Type *pinfo_from, void *projection_data_from,
			    SLang_Array_Type *at_xfrom, SLang_Array_Type *at_yfrom,
			    SLang_Array_Type *at_xto, SLang_Array_Type *at_yto)
{
   unsigned int num = at_xfrom->num_elements;
   unsigned int i;
#define MAX_REPROJECT 1024

   if (at_xfrom->data_type == SLANG_FLOAT_TYPE)
     {
	float *xa = (float *)at_xfrom->data;
	float *ya = (float *)at_yfrom->data;
	float *xb = (float *)at_xto->data;
	float *yb = (float *)at_yto->data;
	double tmpx[MAX_REPROJECT], tmpy[MAX_REPROJECT];
#if 1
	if ((pinfo_to == pinfo_from)
	    && (pinfo_to->reproject_fun_f != NULL))
	  return (*pinfo_to->reproject_fun_f)(projection_data_from, xa, ya,
					      projection_data_to, xb, yb, num);
#endif
	i = 0;
	while (i < num)
	  {
	     unsigned int j;
	     unsigned int dn = num - i;
	     if (dn > MAX_REPROJECT)
	       dn = MAX_REPROJECT;
	     
	     for (j = 0; j < dn; j++)
	       {
		  tmpx[j] = xa[j];
		  tmpy[j] = ya[j];
	       }
	     if (-1 == (*pinfo_from->deproject_fun_d)(projection_data_from,
						      tmpx, tmpy, tmpx, tmpy, 
						      dn))
	       return -1;

	     if (-1 == (*pinfo_to->project_fun_d)(projection_data_to,
						  tmpx, tmpy, tmpx, tmpy, 
						  dn))
	       return -1;

	     for (j = 0; j < dn; j++)
	       {
		  xb[j] = tmpx[j];
		  yb[j] = tmpy[j];
	       }
	     xa += dn;
	     xb += dn;
	     ya += dn;
	     yb += dn;
	     i += dn;
	  }
	return 0;
     }
   else
     {
	
	double *xa = (double *)at_xfrom->data;
	double *ya = (double *)at_yfrom->data;
	double *xb = (double *)at_xto->data;
	double *yb = (double *)at_yto->data;

	if ((pinfo_to == pinfo_from)
	    && (pinfo_to->reproject_fun_d != NULL))
	  return (*pinfo_to->reproject_fun_d)(projection_data_from, xa, ya,
					      projection_data_to, xb, yb, num);

	if (-1 == (*pinfo_from->deproject_fun_d)(projection_data_from, xa, ya,
						 xb, yb, num))
	  return -1;
	if (-1 == (*pinfo_to->project_fun_d)(projection_data_to, xb, yb,
					       xb, yb, num))
	  return -1;
	
	return 0;
     }
}

static void project_intrin_internal (int direction)
{
   SLang_Array_Type *at_xfrom, *at_yfrom;
   SLang_Array_Type *at_xto, *at_yto;
   int is_scalar;
   Projection_Info_Type *pinfo;
   void *projection_data;
   int (*fun_d) (void *, double *, double *, double *, double *, unsigned int);
   int (*fun_f) (void *, float *, float *, float *, float *, unsigned int);

   if (-1 == pop_reusable_arrays (&at_xfrom, &at_yfrom, &at_xto, &at_yto, &is_scalar))
     return;

   if (-1 == pop_projection (&pinfo, &projection_data))
     {
	SLang_free_array (at_yto);
	SLang_free_array (at_xto);
	SLang_free_array (at_yfrom);
	SLang_free_array (at_xfrom);
	return;
     }

   if (direction == 0)
     {
	void *projection_data_to;
	Projection_Info_Type *pinfo_to;
	
	if (-1 == pop_projection (&pinfo_to, &projection_data_to))
	  goto free_return;

	if (-1 == do_reprojection (pinfo_to, projection_data_to,
				   pinfo, projection_data,
				   at_xfrom, at_yfrom, at_xto, at_yto))
	  {
	     (*pinfo_to->free_projection)(projection_data_to);
	     goto free_return;
	  }
	(*pinfo_to->free_projection)(projection_data_to);
     }
   else
     {
	if (direction == 1)
	  {
	     fun_d = pinfo->project_fun_d;
	     fun_f = pinfo->project_fun_f;
	  }
	else
	  {
	     fun_d = pinfo->deproject_fun_d;
	     fun_f = pinfo->deproject_fun_f;
	  }

	if (at_xfrom->data_type == SLANG_FLOAT_TYPE)
	  {
	     if (-1 == (*fun_f)(projection_data, 
				(float*)at_xfrom->data, (float*)at_yfrom->data, 
				(float*)at_xto->data, (float*)at_yto->data, 
				at_xfrom->num_elements))
	       goto free_return;
	  }
	else
	  {
	     if (-1 == (*fun_d)(projection_data,
				(double*)at_xfrom->data, (double*)at_yfrom->data, 
				(double*)at_xto->data, (double*)at_yto->data, 
				at_xfrom->num_elements))
	       goto free_return;
	  }
     }


   (void) push_array_maybe_scalar (at_xto, is_scalar);
   (void) push_array_maybe_scalar (at_yto, is_scalar);

   free_return:
   
   (*pinfo->free_projection)(projection_data);

   SLang_free_array (at_yto);
   SLang_free_array (at_xto);
   SLang_free_array (at_yfrom);
   SLang_free_array (at_xfrom);
}

static void project_intrin (void)
{
   if (SLang_Num_Function_Args != 3)
     {
	SLang_verror (SL_USAGE_ERROR, "Usage: (x,y)=maplib_project(map_info,lon,lat)");
	return;
     }
   project_intrin_internal (1);
}

static void deproject_intrin (void)
{
   if (SLang_Num_Function_Args != 3)
     {
	SLang_verror (SL_USAGE_ERROR, "Usage: (x,y)=maplib_deproject(map_info,lon,lat)");
	return;
     }
   project_intrin_internal (-1);
}

static void reproject_intrin (void)
{
   if (SLang_Num_Function_Args != 4)
     {
	SLang_verror (SL_USAGE_ERROR, "Usage: (x1,y1)=maplib_reproject(map_info1,map_info0,x0,y0)");
	return;
     }
   project_intrin_internal (0);
}

static void rotate2d_d (double *x, double *y, unsigned int n,
			double theta, double x_0, double y_0, 
			double *x_1, double *y_1)
{
   unsigned int i;
   double c = cos (theta);
   double s = sin (theta);
   
   for (i = 0; i < n; i++)
     {
	double x_i = x[i]-x_0;
	double y_i = y[i]-y_0;
	
	x_1[i] = x_0 + c*x_i - s*y_i;
	y_1[i] = y_0 + s*x_i + c*y_i;
     }
}

static void rotate2d_f (float *x, float *y, unsigned int n,
			double theta, double x_0, double y_0, 
			float *x_1, float *y_1)
{
   unsigned int i;
   double c = cos (theta);
   double s = sin (theta);
   
   for (i = 0; i < n; i++)
     {
	double x_i = x[i]-x_0;
	double y_i = y[i]-y_0;
	
	x_1[i] = x_0 + c*x_i - s*y_i;
	y_1[i] = y_0 + s*x_i + c*y_i;
     }
}

static void rotate2d_intrin (void)
{
   double theta;
   SLang_Array_Type *x, *y, *x_1, *y_1;
   int is_scalar;
   double x_0 = 0.0, y_0 = 0.0;

   if (SLang_Num_Function_Args == 5)
     {
	if (-1 == POP_DOUBLE (&y_0,NULL,NULL))
	  return;
	if (-1 == POP_DOUBLE (&x_0,NULL,NULL))
	  return;
     }
   else if (SLang_Num_Function_Args != 3)
     {
	SLang_verror (SL_USAGE_ERROR, "Usage: (x1,y1)=maplib_rotate2d(x,y,theta_radians[,x0,y0])");
	return;
     }
   if (-1 == POP_DOUBLE (&theta, NULL, NULL))
     return;

   if (-1 == pop_reusable_arrays (&x, &y, &x_1, &y_1, &is_scalar))
     return;
   
   if (x->data_type == SLANG_FLOAT_TYPE)
     rotate2d_f ((float *)x->data, (float *)y->data, x->num_elements,
		 theta, x_0, y_0, (float *)x_1->data, (float *)y_1->data);
   else
     rotate2d_d ((double *)x->data, (double *)y->data, x->num_elements,
		 theta, x_0, y_0, (double *)x_1->data, (double *)y_1->data);
   
   (void) push_array_maybe_scalar (x_1, is_scalar);
   (void) push_array_maybe_scalar (y_1, is_scalar);

   SLang_free_array (y_1);
   SLang_free_array (x_1);
   SLang_free_array (y);
   SLang_free_array (x);
}


static void meshgrid_intrin (void)
{
   SLang_Array_Type *at_x, *at_y, *at_xx, *at_yy;
   int data_type;
   int is_scalar;
   int dims[2];
   unsigned int nx, ny;

   if (SLang_Num_Function_Args != 2)
     {
	SLang_verror (SL_USAGE_ERROR, "(X,Y)=maplib_meshgrid (x,y)");
	return;
     }
   
   if (-1 == pop_2_floating_arrays_1 (&at_x, &at_y, &is_scalar))
     return;
   
   data_type = at_x->data_type;
   nx = (unsigned int) at_x->num_elements;
   ny = (unsigned int) at_y->num_elements;
   dims[0] = (int) nx;
   dims[1] = (int) ny;

   if (NULL == (at_xx = SLang_create_array (data_type, 0, NULL, dims, 2)))
     {
	SLang_free_array (at_y);
	SLang_free_array (at_x);
	return;
     }
   
   if (NULL == (at_yy = SLang_create_array (data_type, 0, NULL, dims, 2)))
     {
	SLang_free_array (at_xx);
	SLang_free_array (at_y);
	SLang_free_array (at_x);
	return;
     }
   
   if (data_type == SLANG_FLOAT_TYPE)
     {
	float *x = (float *)at_x->data;
	float *y = (float *)at_y->data;
	float *xx = (float *)at_xx->data;
	float *yy = (float *)at_yy->data;
	unsigned int i, j;
	
	/* xx_ij = x_i; yy_ij = y_j */
	for (i = 0; i < nx; i++)
	  {
	     float x_i = x[i];
	     for (j = 0; j < ny; j++)
	       {
		  xx[j] = x_i;
		  yy[j] = y[j];
	       }
	     xx += ny;
	     yy += ny;
	  }
     }
   else
     {
	double *x = (double *)at_x->data;
	double *y = (double *)at_y->data;
	double *xx = (double *)at_xx->data;
	double *yy = (double *)at_yy->data;
	unsigned int i, j;
	
	/* xx_ij = x_i; yy_ij = y_j */
	for (i = 0; i < nx; i++)
	  {
	     double x_i = x[i];
	     for (j = 0; j < ny; j++)
	       {
		  xx[j] = x_i;
		  yy[j] = y[j];
	       }
	     xx += ny;
	     yy += ny;
	  }
     }
   
   SLang_free_array (at_y);
   SLang_free_array (at_x);
   (void) SLang_push_array (at_xx, 1);
   (void) SLang_push_array (at_yy, 1);
}

static SLang_Intrin_Fun_Type Module_Intrinsics [] =
{
   MAKE_INTRINSIC_0 ("maplib_get_projections", get_projections, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0 ("maplib_project", project_intrin, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0 ("maplib_deproject", deproject_intrin, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0 ("maplib_reproject", reproject_intrin, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_1 ("maplib_new", project_new, SLANG_VOID_TYPE, SLANG_STRING_TYPE),
   MAKE_INTRINSIC_0 ("maplib_meshgrid", meshgrid_intrin, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0 ("maplib_rotate2d", rotate2d_intrin, SLANG_VOID_TYPE),

   SLANG_END_INTRIN_FUN_TABLE
};

static SLang_Intrin_Var_Type Module_Variables [] =
{
   MAKE_VARIABLE("_maplib_module_version_string", &Module_Version_String, SLANG_STRING_TYPE, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_IConstant_Type Module_IConstants [] =
{
   MAKE_ICONSTANT("_maplib_module_version", MODULE_VERSION_NUMBER),
   SLANG_END_ICONST_TABLE
};


int init_maplib_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns = SLns_create_namespace (ns_name);
   if (ns == NULL)
     return -1;

   if (
       (-1 == SLns_add_intrin_var_table (ns, Module_Variables, NULL))
       || (-1 == SLns_add_intrin_fun_table (ns, Module_Intrinsics, NULL))
       || (-1 == SLns_add_iconstant_table (ns, Module_IConstants, NULL)))
     return -1;

   /* Hack */
   Maplib_NaN = sqrt (-1.0);

   return 0;
}

/* This function is optional */
void deinit_maplib_module (void)
{
}
