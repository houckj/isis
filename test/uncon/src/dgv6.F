*************************************************************************
*                 dgv6 
* dennis, gay, and vu, tech. rep. 83-16, rice university
*************************************************************************

      subroutine getfun( x, n, f, m, ftf, fj, lfj, g, mode)

      implicit double precision (a-h,o-z)

      integer            n, m, lfj, mode

      double precision   x(n), f(m), ftf, fj(lfj,n), g(n)

      integer            nprob, nprobs, nstart, nstrts
      common /PROBLM/    nprob, nprobs, nstart, nstrts

      integer            nout
      common /IOUNIT/    nout

      logical            lf, lj

      double precision   summx,summy,suma,sumb,sumc,sumd,sume,sumf
      common /PARAM1/    summx,summy,suma,sumb,sumc,sumd,sume,sumf
      save   /PARAM1/    

      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0, six = 6.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .lt. -1)  goto    30

      na = mode / 1000
      nt = mode - na*1000
      nb = nt / 100
      nh = nt - nb*100
      nc = nh / 10
      nd = nh - nc*10

      lf = (na .ne. 0) .or. (nb .ne. 0) .or. (nd .ne. 0)
      lj = (nc .ne. 0) .or. (nd .ne. 0)

      a = x(1)
      b = summx - a
      c = x(2)
      d = summy - c
      t = x(3)
      u = x(4)
      v = x(5)
      w = x(6)

      tv    = t*v
      tt    = t*t
      vv    = v*v
      tsvs  = tt - vv
      ts3vs = tt - three*vv
      vs3ts = vv - three*tt
      uw    = u*w
      uu    = u*u
      ww    = w*w
      usws  = uu - ww
      us3ws = uu - three*ww
      ws3us = ww - three*uu

      if (lf .and. lj)  goto 300
      if (lf)           goto 100
      if (lj)           goto 200

*-----------------------------------------------------------------------

   10 continue

      n      = 6
      m      = n

      nprobs = 5
      nstrts = 1

      if (nout .gt. 0)  write( nout, 9999) nprob, n, m

      goto ( 11, 12, 13, 14, 15), nprob

   11 continue

      summx =  0.485d0
      summy = -0.0019d0
      suma  = -0.0581d0
      sumb  =  0.015d0
      sumc  =  0.105d0
      sumd  =  0.0406d0
      sume  =  0.167d0
      sumf  = -0.399d0

      return

   12 continue

      summx = - 0.69d0
      summy = - 0.044d0
      suma  = - 1.57d0
      sumb  = - 1.31d0
      sumc  = - 2.65d0
      sumd  =   2.0d0
      sume  = -12.6d0
      sumf  =   9.48d0

      return

   13 continue

      summx = - 0.816d0
      summy = - 0.017d0
      suma  = - 1.826d0
      sumb  = - 0.754d0
      sumc  = - 4.839d0
      sumd  = - 3.259d0
      sume  = -14.023d0
      sumf  =  15.467d0

      return

   14 continue

      summx = - 0.809d0
      summy = - 0.021d0
      suma  = - 2.04d0
      sumb  = - 0.614d0
      sumc  = - 6.903d0
      sumd  = - 2.934d0
      sume  = -26.328d0
      sumf  =  18.639d0

      return

   15 continue

      summx = - 0.807d0
      summy = - 0.021d0
      suma  = - 2.379d0
      sumb  = - 0.364d0
      sumc  = -10.541d0
      sumd  = - 1.961d0
      sume  = -51.551d0
      sumf  =  21.053d0

      return

*-----------------------------------------------------------------------

   20 continue

      goto ( 21, 22, 23, 24, 25), nprob

   21 continue

      x(1)  =  0.299d0
c     x(2)  =  0.186d0
      x(2)  = -0.0273d0
c     x(4)  =  0.0254d0
      x(3)  = -0.474d0
      x(4)  =  0.474d0
      x(5)  = -0.0892d0
      x(6)  =  0.0892d0

      return

   22 continue

      x(1)  = -0.3d0
c     x(2)  = -0.39d0
      x(2)  =  0.3d0
c     x(4)  = -0.344d0
      x(3)  = -1.2d0
      x(4)  =  2.69d0
      x(5)  =  1.59d0
      x(6)  = -1.5d0

      return

   23 continue

      x(1)  = -0.041d0
c     x(2)  = -0.775d0
      x(2)  =  0.03d0
c     x(4)  = -0.047d0
      x(3)  = -2.565d0
      x(4)  =  2.565d0
      x(5)  = -0.754d0
      x(6)  =  0.754d0

      return

   24 continue

      x(1)  = -0.056d0
c     x(2)  = -0.753d0
      x(2)  =  0.026d0
c     x(4)  = -0.047d0
      x(3)  = -2.991d0
      x(4)  =  2.991d0
      x(5)  = -0.568d0
      x(6)  =  0.568d0

      return

   25 continue

      x(1)  = -0.074d0
c     x(2)  = -0.733d0
      x(2)  =  0.013d0
c     x(4)  = -0.034d0
      x(3)  = -3.632d0
      x(4)  =  3.632d0
      x(5)  = -0.289d0
      x(6)  =  0.289d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

c     f(1) = a + b - summx
c     f(2) = c + d - summy
      f(1) = t*a + u*b - v*c - w*d - suma
      f(2) = v*a + w*b + t*c + u*d - sumb
      f(3) = a*tsvs - two*c*t*v + b*usws - two*d*u*w - sumc
      f(4) = c*tsvs + two*a*t*v + d*usws + two*b*u*w - sumd
      f(5) = a*t*ts3vs + c*v*vs3ts + b*u*us3ws + d*w*ws3us - sume
      f(6) = c*t*ts3vs - a*v*vs3ts + d*u*us3ws - b*w*ws3us - sumf

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 210  continue

      fj( 1, 1) =  t - u
      fj( 1, 2) = -v + w
      fj( 1, 3) =  a
      fj( 1, 4) =  b
      fj( 1, 5) = -c
      fj( 1, 6) = -d

      fj( 2, 1) =  v - w
      fj( 2, 2) =  t - u
      fj( 2, 3) =  c
      fj( 2, 4) =  d
      fj( 2, 5) =  a
      fj( 2, 6) =  b

      fj( 3, 1) =  tsvs - usws
      fj( 3, 2) = -two*(tv - uw)
      fj( 3, 3) =  two*(a*t - c*v)
      fj( 3, 4) =  two*(b*u - d*w)
      fj( 3, 5) = -two*(a*v + c*t)
      fj( 3, 6) = -two*(b*w + d*u)

      fj( 4, 1) =  two*(tv - uw)
      fj( 4, 2) =  tsvs - usws
      fj( 4, 3) =  two*(c*t + a*v)
      fj( 4, 4) =  two*(d*u + b*w)
      fj( 4, 5) =  two*(a*t - c*v)
      fj( 4, 6) =  two*(b*u - d*w)

      fj( 5, 1) =  t*ts3vs - u*us3ws
      fj( 5, 2) =  v*vs3ts - w*ws3us
      fj( 5, 3) =  three*(a*tsvs - two*c*tv)
      fj( 5, 4) =  three*(b*usws - two*d*uw)
      fj( 5, 5) = -three*(c*tsvs + two*a*tv)
      fj( 5, 6) = -three*(d*usws + two*b*uw)

      fj( 6, 1) = -v*vs3ts + w*ws3us
      fj( 6, 2) =  t*ts3vs - u*us3ws
      fj( 6, 3) =  three*(c*tsvs + two*a*tv)
      fj( 6, 4) =  three*(d*usws + two*b*uw)
      fj( 6, 5) =  three*(a*tsvs - two*c*tv)
      fj( 6, 6) =  three*(b*usws - two*d*uw)

      return

 300  continue

c     f(1) = a + b - summx
c     f(2) = c + d - summy
      f(1) = t*a + u*b - v*c - w*d - suma
      f(2) = v*a + w*b + t*c + u*d - sumb
      f(3) = a*tsvs - two*c*t*v + b*usws - two*d*u*w - sumc
      f(4) = c*tsvs + two*a*t*v + d*usws + two*b*u*w - sumd
      f(5) = a*t*ts3vs + c*v*vs3ts + b*u*us3ws + d*w*ws3us - sume
      f(6) = c*t*ts3vs - a*v*vs3ts + d*u*us3ws - b*w*ws3us - sumf

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      do 310 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 310  continue

      fj( 1, 1) =  t - u
      fj( 1, 2) = -v + w
      fj( 1, 3) =  a
      fj( 1, 4) =  b
      fj( 1, 5) = -c
      fj( 1, 6) = -d

      fj( 2, 1) =  v - w
      fj( 2, 2) =  t - u
      fj( 2, 3) =  c
      fj( 2, 4) =  d
      fj( 2, 5) =  a
      fj( 2, 6) =  b

      fj( 3, 1) =  tsvs - usws
      fj( 3, 2) = -two*(tv - uw)
      fj( 3, 3) =  two*(a*t - c*v)
      fj( 3, 4) =  two*(b*u - d*w)
      fj( 3, 5) = -two*(a*v + c*t)
      fj( 3, 6) = -two*(b*w + d*u)

      fj( 4, 1) =  two*(tv - uw)
      fj( 4, 2) =  tsvs - usws
      fj( 4, 3) =  two*(c*t + a*v)
      fj( 4, 4) =  two*(d*u + b*w)
      fj( 4, 5) =  two*(a*t - c*v)
      fj( 4, 6) =  two*(b*u - d*w)

      fj( 5, 1) =  t*ts3vs - u*us3ws
      fj( 5, 2) =  v*vs3ts - w*ws3us
      fj( 5, 3) =  three*(a*tsvs - two*c*tv)
      fj( 5, 4) =  three*(b*usws - two*d*uw)
      fj( 5, 5) = -three*(c*tsvs + two*a*tv)
      fj( 5, 6) = -three*(d*usws + two*b*uw)

      fj( 6, 1) = -v*vs3ts + w*ws3us
      fj( 6, 2) =  t*ts3vs - u*us3ws
      fj( 6, 3) =  three*(c*tsvs + two*a*tv)
      fj( 6, 4) =  three*(d*usws + two*b*uw)
      fj( 6, 5) =  three*(a*tsvs - two*c*tv)
      fj( 6, 6) =  three*(b*usws - two*d*uw)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' dennis, gay, and vu function                  data set #', i3//,
     *'        number of variables =', i4,'  (variable)'/,
     *'        number of functions =', i4,'  (  >= n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   summx,summy,suma,sumb,sumc,sumd,sume,sumf
      common /PARAM1/    summx,summy,suma,sumb,sumc,sumd,sume,sumf
      save   /PARAM1/   

      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0, six = 6.d0)

*=======================================================================

      a = x(1)
      b = summx - a
      c = x(2)
      d = summy - c
      t = x(3)
      u = x(4)
      v = x(5)
      w = x(6)

      tv    = t*v
      tt    = t*t
      vv    = v*v
      tsvs  = tt - vv
      ts3vs = tt - three*vv
      vs3ts = vv - three*tt
      uw    = u*w
      uu    = u*u
      ww    = w*w
      usws  = uu - ww
      us3ws = uu - three*ww
      ws3us = ww - three*uu

      do 10 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
   10 continue

      linear = .false.
      goto ( 100, 200, 300, 400, 500, 600), k

 100  continue

      hess( 1, 3) =  one
      hess( 3, 1) =  one
      hess( 1, 4) = -one
      hess( 4, 1) = -one
      hess( 2, 5) = -one
      hess( 5, 2) = -one
      hess( 2, 6) =  one
      hess( 6, 2) =  one

      return

 200  continue

      hess( 1, 5) =  one
      hess( 5, 1) =  one
      hess( 1, 6) = -one
      hess( 6, 1) = -one
      hess( 2, 3) =  one
      hess( 3, 2) =  one
      hess( 2, 4) = -one
      hess( 4, 2) = -one

      return

 300  continue

      hess( 1, 3) =  two*t
      hess( 3, 1) =  hess( 1, 3)
      hess( 1, 4) = -two*u
      hess( 4, 1) =  hess( 1, 4)
      hess( 1, 5) = -two*v
      hess( 5, 1) =  hess( 1, 5)
      hess( 1, 6) =  two*w
      hess( 6, 1) =  hess( 1, 6)
      hess( 2, 3) = -two*v
      hess( 3, 2) =  hess( 2, 3)
      hess( 2, 4) =  two*w
      hess( 4, 2) =  hess( 2, 4)
      hess( 2, 5) = -two*t
      hess( 5, 2) =  hess( 2, 5)
      hess( 2, 6) =  two*u
      hess( 6, 2) =  hess( 2, 6)
      hess( 3, 3) =  two*a
      hess( 3, 5) = -two*c
      hess( 5, 3) =  hess( 3, 5)
      hess( 4, 4) =  two*b
      hess( 4, 6) = -two*d
      hess( 6, 4) =  hess( 4, 6)
      hess( 5, 5) = -two*a
      hess( 6, 6) = -two*b

      return

 400  continue

      hess( 1, 3) =  two*v
      hess( 3, 1) =  hess( 1, 3)
      hess( 1, 4) = -two*w
      hess( 4, 1) =  hess( 1, 4)
      hess( 1, 5) =  two*t
      hess( 5, 1) =  hess( 1, 5)
      hess( 1, 6) = -two*u
      hess( 6, 1) =  hess( 1, 6)
      hess( 2, 3) =  two*t
      hess( 3, 2) =  hess( 2, 3)
      hess( 2, 4) = -two*u
      hess( 4, 2) =  hess( 2, 4)
      hess( 2, 5) = -two*v
      hess( 5, 2) =  hess( 2, 5)
      hess( 2, 6) =  two*w
      hess( 6, 2) =  hess( 2, 6)
      hess( 3, 3) =  two*c
      hess( 3, 5) =  two*a
      hess( 5, 3) =  hess( 3, 5)
      hess( 4, 4) =  two*d
      hess( 4, 6) =  two*b
      hess( 6, 4) =  hess( 4, 6)
      hess( 5, 5) = -two*c
      hess( 6, 6) = -two*d

      return

 500  continue

      hess( 1, 3) =  three*tsvs
      hess( 3, 1) =  hess( 1, 3)
      hess( 1, 4) = -three*usws
      hess( 4, 1) =  hess( 1, 4)
      hess( 1, 5) = -six*tv
      hess( 5, 1) =  hess( 1, 5)
      hess( 1, 6) =  six*uw
      hess( 6, 1) =  hess( 1, 6)
      hess( 2, 3) = -six*tv
      hess( 3, 2) =  hess( 2, 3)
      hess( 2, 4) =  six*uw
      hess( 4, 2) =  hess( 2, 4)
      hess( 2, 5) = -three*tsvs
      hess( 5, 2) =  hess( 2, 5)
      hess( 2, 6) =  three*usws
      hess( 6, 2) =  hess( 2, 6)
      hess( 3, 3) =  six*(a*t - c*v)
      hess( 3, 5) = -six*(a*v + c*t)
      hess( 5, 3) =  hess( 3, 5)
      hess( 4, 4) =  six*(b*u - d*w)
      hess( 4, 6) = -six*(b*w + d*u)
      hess( 6, 4) =  hess( 4, 6)
      hess( 5, 5) =  six*(c*v - a*t)
      hess( 6, 6) =  six*(d*w - b*u)

      return

 600  continue

      hess( 1, 3) =  six*tv
      hess( 3, 1) =  hess( 1, 3)
      hess( 1, 4) = -six*uw
      hess( 4, 1) =  hess( 1, 4)
      hess( 1, 5) =  three*tsvs
      hess( 5, 1) =  hess( 1, 5)
      hess( 1, 6) = -three*usws
      hess( 6, 1) =  hess( 1, 6)
      hess( 2, 3) =  three*tsvs
      hess( 3, 2) =  hess( 2, 3)
      hess( 2, 4) = -three*usws
      hess( 4, 2) =  hess( 2, 4)
      hess( 2, 5) = -six*tv
      hess( 5, 2) =  hess( 2, 5)
      hess( 2, 6) =  six*uw
      hess( 6, 2) =  hess( 2, 6)
      hess( 3, 3) =  six*(c*t + a*v)
      hess( 3, 5) = -six*(c*v - a*t)
      hess( 5, 3) =  hess( 3, 5)
      hess( 4, 4) =  six*(d*u + b*w)
      hess( 4, 6) = -six*(d*w - b*u)
      hess( 6, 4) =  hess( 4, 6)
      hess( 5, 5) = -six*(a*v + c*t)
      hess( 6, 6) = -six*(b*w + d*u)

      return

      end

