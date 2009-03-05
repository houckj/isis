*******************************************************************************
*              devilliers and glasser function 1
*   siam j. n. a., vol. 18, no. 6 (december 1981) 1139-1154
*******************************************************************************

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

      double precision   x1, x2, x3, x4
      double precision   ai, ci, si, ti, t2, yi

      double precision   ddot

      intrinsic          dble, cos, sin

      double precision   one, point1
      parameter         (one = 1.d0, point1 = .1d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .lt. -1)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)

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
      nstrts =  4

      n      =  4
      m      = 24

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      goto ( 21, 22, 23, 24), nstart

   21 continue

      x(1) =  1.d0
      x(2) =  8.d0
      x(3) =  4.d0
      x(4) =  4.412d0

   22 continue

      x(1) =  1.d0
      x(2) =  8.d0
      x(3) =  8.d0
      x(4) =  1.d0

      return

   23 continue

      x(1) =  1.d0
      x(2) =  8.d0
      x(3) =  1.d0
      x(4) =  4.412d0

      return

   24 continue

      x(1) =  1.d0
      x(2) =  8.d0
      x(3) =  4.d0
      x(4) =  1.d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf  = 0.d0

      x(1) =  60.137d0
      x(2) =   1.371d0
      x(3) =   3.112d0
      x(4) =   1.761d0

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ti   = point1*dble(i-1)
        yi   = 60.137d0 * (1.371d0 ** ti) * sin( 3.112d0*ti + 1.761d0)
        f(i) = x1 * (x2 ** ti) * sin( x3*ti + x4 ) - yi
 110  continue

      if (nb .gt. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        ti = point1*dble(i-1)
        t2 = x2**ti
        ai = x3*ti + x4
        si = sin(ai)
        ci = cos(ai)
        fj( i, 1) = t2 * si
        fj( i, 2) = x1 * (ti*(x2**(ti - one))) * si
        fj( i, 3) = x1 * t2 * ti * ci
        fj( i, 4) = x1 * t2 * ci
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ti = point1*dble(i-1)
        t2 = x2**ti
        ai = x3*ti + x4
        si = sin(ai)
        ci = cos(ai)
        yi   = 60.137d0 * (1.371d0 ** ti) * sin( 3.112d0*ti + 1.761d0)
        f(i) = x1 * t2 * si - yi
        fj( i, 1) = t2 * si
        fj( i, 2) = x1 * (ti*(x2**(ti - one))) * si
        fj( i, 3) = x1 * t2 * ti * ci
        fj( i, 4) = x1 * t2 * ci
 310  continue

      if (nb .gt. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 320 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 320  continue

      return

9999  format(/'1',70('=')//,
     *' devilliers and glasser function 1       starting value no. 1'//,
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
      double precision   ak, d2, ck, sk, tk, t2

      intrinsic          dble, cos, sin

      double precision   zero, one, two, point1
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = .1d0)

*=======================================================================

      linear = .false.

      x1  = x(1)
      x2  = x(2)
      x3  = x(3)
      x4  = x(4)

      tk  = point1*dble(k-1)
      ak  = x3*tk + x4
      sk =  sin(ak)
      ck =  cos(ak)
      t2 = x2**tk
      d2 = tk * (x2 ** (tk - one))

      hess(1,1) =  zero
      hess(1,2) =  d2 * sk
      hess(1,3) =  t2 * tk * ck
      hess(1,4) =  t2 * ck
      hess(2,1) =  hess(1,2)
      hess(2,2) =  x1 * tk * (tk - one) * (x2 ** (tk - two)) * sk
      hess(2,3) =  x1 * d2 * tk * ck
      hess(2,4) =  x1 * d2 * ck
      hess(3,1) =  hess(1,3)
      hess(3,2) =  hess(2,3)
      hess(3,3) = -x1 * t2 * tk * tk * sk
      hess(3,4) = -x1 * t2 * tk * sk
      hess(4,1) =  hess(1,4)
      hess(4,2) =  hess(2,4)
      hess(4,3) =  hess(3,4)
      hess(4,4) = -x1 * t2 * sk

      return
      end
