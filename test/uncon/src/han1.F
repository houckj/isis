******************************************************************************
*              hanson 1
* from salane, siam j. sci. stat. comp. 1987
******************************************************************************

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

      intrinsic          cos, sin

      double precision   a, b
      common /PARAM1/    a(16), b(16)
      save   /PARAM1/ 

      double precision   zero, one, two, point1  
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = .1d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)

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

   10 continue

      nprobs = 1
      nstrts = 1

      n      = 2
      m      = 16

      a( 1)  =  0.d0
      a( 2)  =  1.d0
      a( 3)  =  2.d0
      a( 4)  =  3.d0
      a( 5)  =  4.d0
      a( 6)  =  5.d0
      a( 7)  =  6.d0
      a( 8)  =  8.d0
      a( 9)  = 10.d0
      a(10)  = 12.d0
      a(11)  = 15.d0
      a(12)  = 20.d0
      a(13)  = 25.d0
      a(14)  = 30.d0
      a(15)  = 40.d0
      a(16)  = 50.d0

      b( 1)  =  0.d0
      b( 2)  =  1.d0
      b( 3)  =  2.d0
      b( 4)  =  3.d0
      b( 5)  =  5.d0
      b( 6)  =  6.d0
      b( 7)  =  8.d0
      b( 8)  = 11.d0
      b( 9)  = 13.d0
      b(10)  = 12.d0
      b(11)  =  9.d0
      b(12)  =  6.d0
      b(13)  =  3.d0
      b(14)  =  2.d0
      b(15)  =  1.5d0
      b(16)  =  1.d0

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  10.d0
      x(2) =  0.15d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 2.08700d2

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        f(i) = x1 * sin( x2*a(i) ) - b(i)
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        fj( i, 1) = sin( x2*a(i) )
        fj( i, 2) = x1 * a(i) * cos( x2*a(i) )
 210  continue

      return

 300  continue

      do 310 i = 1, m
        f(i) = x1 * sin( x2*a(i) ) - b(i)
        fj( i, 1) = sin( x2*a(i) )
        fj( i, 2) = x1 * a(i) * cos( x2*a(i) )
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 320 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 320  continue

      return

9999  format(/'1',70('=')//,
     *' hansons function 1 (salane) '//,
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

      double precision   x1, x2, ak

      intrinsic          cos, sin

      double precision   a, b
      common /PARAM1/    a(16), b(16)
      save   /PARAM1/

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      linear = .false.

      x1  = x(1)
      x2  = x(2)

      ak  = a(k)

      hess(1,1) =  zero
      hess(1,2) =  ak * cos( x2*ak )
      hess(2,1) =  hess(1,2)
      hess(2,2) = -x1 * ak * ak * sin( x2*ak )

      return
      end

