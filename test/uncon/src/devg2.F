***************************************************************************
*              devilliers and glasser function 2
*   siam j. n. a., vol. 18, no. 6 (december 1981) 1139-1154
***************************************************************************

      subroutine getfun( x, n, f, m, ftf, fj, lfj, g, mode)

      implicit double precision (a-h,o-z)

      integer            n, m, lfj, mode

      double precision   x(n), f(m), ftf, fj(lfj,n), g(n)

      integer            nprob, nprobs, nstart, nstrts
      common /PROBLM/    nprob, nprobs, nstart, nstrts

      integer            nout
      common /IOUNIT/    nout

      logical            lf, lj

      integer            na, nb, nc, nd, nt, nh

      integer            i, j

      double precision   x1, x2, x3, x4, x5
      double precision   ac, as, ch, ci, e5, si, ti, t2, yi
      double precision   d2, d3, d4, d5

      double precision   ddot

      intrinsic          dble, exp, cos, sin, cosh, tanh

      double precision   one, two, point1
      parameter         (one = 1.d0, two = 2.d0, point1 = .1d0)

      double precision   e0
      parameter         (e0 = 1.66030d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .lt. -1)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)
      x5 = x(5)

      if (x2 .lt. 0.d0) then
c        write( 6, *) '*** getfun : function not defined : x2 < 0'
c        stop
         mode=999
         return        
      end if

      na = mode / 1000
      nt = mode - na*1000
      nb = nt / 100
      nh = nt - nb*100
      nc = nh / 10
      nd = nh - nc*10

      lf = (na .ne. 0) .or. (nb .ne. 0) .or. (nd .ne. 0)
      lj = (nc .ne. 0) .or. (nd .ne. 0)

      if (lf .and. lj)  goto 300
      if (lf)           goto 100
      if (lj)           goto 200

*-----------------------------------------------------------------------

  10  continue

      nprobs =  1
      nstrts =  5

      n      =  5
      m      = 16

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      goto ( 21, 22, 23, 24, 25), nstart

   21 continue

      x(1) = 45.0d0
      x(2) =  2.0d0
      x(3) =  2.5d0
      x(4) =  1.5d0
      x(5) =  0.9d0

      return

   22 continue

      x(1) = 42.0d0
      x(2) =  0.8d0
      x(3) =  1.4d0
      x(4) =  1.8d0
      x(5) =  1.0d0

      return

   23 continue

      x(1) = 45.0d0
      x(2) =  2.0d0
      x(3) =  2.1d0
      x(4) =  2.0d0
      x(5) =  0.9d0

      return

   24 continue

      x(1) = 45.0d0
      x(2) =  2.5d0
      x(3) =  1.7d0
      x(4) =  1.0d0
      x(5) =  1.0d0

      return

   25 continue

      x(1) = 35.0d0
      x(2) =  2.5d0
      x(3) =  1.7d0
      x(4) =  1.0d0
      x(5) =  1.0d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf  = 0.d0

      x(1) =  53.81d0
      x(2) =   1.27d0
      x(3) =   3.012d0
      x(4) =   2.13d0
      x(5) =   0.507d0

      return

*-----------------------------------------------------------------------

 100  continue

c     e0 = exp(0.507d0)
      e5 = exp(x5)
      do 110 i = 1, m
        ti   = point1*dble(i-1)
        yi   = 53.81d0*(1.27d0**ti)*tanh(3.012d0*ti+sin(2.13d0*ti))
     *         *cos(ti*e0)
        f(i) = x1*(x2**ti)*tanh(x3*ti+sin(x4*ti))*cos(ti*e5) - yi
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      e5 = exp(x5)
      do 210 i = 1, m
        ti = point1*dble(i-1)
        t2 = x2**ti
        as = x4*ti
        si = sin(as)
        ah = x3*ti + si
        hi = tanh(ah)
        ch = one / cosh(ah)
        ac = ti*e5
        ci = cos(ac)
        d2   =  ti * (x2**(ti - one))
        d3   =  ch * ch * ti
        d4   =  d3 * cos(as)
        d5   = -ac * sin(ac)
        fj( i, 1) =  t2 * hi * ci
        fj( i, 2) =  x1 * d2 * hi * ci
        fj( i, 3) =  x1 * t2 * d3 * ci
        fj( i, 4) =  x1 * t2 * d4 * ci
        fj( i, 5) =  x1 * t2 * hi * d5
 210  continue

      return

 300  continue

c     e0 = exp(0.507d0)
      e5 = exp(x5)
      do 310 i = 1, m
        ti = point1*dble(i-1)
        t2 = x2**ti
        as = x4*ti
        si = sin(as)
        ah = x3*ti + si
        hi = tanh(ah)
        ch = one / cosh(ah)
        ac = ti*e5
        ci = cos(ac)
        yi   = 53.81d0*(1.27d0**ti)*tanh(3.012d0*ti+sin(2.13d0*ti))
     *        *cos(ti*e0)
        f(i) = x1 * t2 * hi * ci - yi
        d2   =  ti * (x2**(ti - one))
        d3   =  ch * ch * ti
        d4   =  d3 * cos(as)
        d5   = -ac * sin(ac)
        fj( i, 1) =  t2 * hi * ci
        fj( i, 2) =  x1 * d2 * hi * ci
        fj( i, 3) =  x1 * t2 * d3 * ci
        fj( i, 4) =  x1 * t2 * d4 * ci
        fj( i, 5) =  x1 * t2 * hi * d5
 310  continue
c
      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)
c
      if (nd .eq. 0)  return
c
      do 320 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 320  continue
c
      return
c
9999  format(/'1',70('=')//,
     *' devilliers and glasser function 2       starting value no. 1'//,
     *'        number of variables =', i4/,
     *'        number of functions =', i4//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear
     
      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   x1, x2, x3, x4
      double precision   bk, dd
      double precision   a4, a5, c4, c5, d2, d3, d4, d5, s4, e5, t2
      double precision   a34, s34, t34, d22, d33, d34, d44, d55

      intrinsic          dble, exp, cos, sin, cosh, tanh

      double precision   zero, one, two, point1
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = .1d0)

*=======================================================================

      linear = .false.

      x1  = x(1)
      x2  = x(2)
      x3  = x(3)
      x4  = x(4)

      e5  = exp(x(5))
      bk  = point1*dble(k-1)
      t2  = x2**bk
      d2  = bk * (x2**(bk - one))
      d22 = bk * (bk - one) * (x2 ** (bk - two))
      a4  = bk * x4
      s4  = sin(a4)
      c4  = cos(a4)
      a34 = x3*bk + s4
      t34 = tanh(a34)
      s34 = one / (cosh(a34)**2)
      d3  =  s34 * bk
      d4  =  d3 * c4
      dd  = two * s34 * s34 * t34
      d33 =  bk * bk * dd
      d34 =  bk * c4 * dd
      d44 =  c4 * d34 - d3 * bk * s4
      a5  = bk * e5
      c5  = cos(a5)
      t5  = c5
      s5  = sin(a5)
      d5  = -a5 * s5
      d55 = -a5 * (a5 * c5 + s5)

      hess(1,1) =  zero
      hess(1,2) =  d2 * t34 * c5
      hess(1,3) =  t2 * d3 * c5
      hess(1,4) =  t2 * d4 * c5
      hess(1,5) =  t2 * t34 * d5
      hess(2,1) =  hess(1,2)
      hess(2,2) =  x1 * d22 * t34 * c5
      hess(2,3) =  x1 * d2 * d3 * c5
      hess(2,4) =  x1 * d2 * d4 * c5
      hess(2,5) =  x1 * d2 * t34 * d5
      hess(3,1) =  hess(1,3)
      hess(3,2) =  hess(2,3)
      hess(3,3) =  x1 * t2 * d33 * c5
      hess(3,4) =  x1 * t2 * d34 * c5
      hess(3,5) =  x1 * t2 * d3 * d5
      hess(4,1) =  hess(1,4)
      hess(4,2) =  hess(2,4)
      hess(4,3) =  hess(3,4)
      hess(4,4) =  x1 * t2 * d44 * c5
      hess(4,5) =  x1 * t2 * d4 * d5
      hess(5,1) =  hess(1,5)
      hess(5,2) =  hess(2,5)
      hess(5,3) =  hess(3,5)
      hess(5,4) =  hess(4,5)
      hess(5,5) =  x1 * t2 * t34 * d55

      return
      end
