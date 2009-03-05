***************************************************************************
*               broyden tridiagonal function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
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

      integer            n1, i, ip1, im1, j

      double precision   xi, xip1, xim1

      double precision   ddot

      double precision   zero, one, two, three, four
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0, four = 4.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

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

      nprobs = 1
      nstrts = 1

      n      = 10
      m      = n

      n1  = n + 1

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      call dcopy( n, (-one), 0, x, 1)

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      i  =   1
      xi = x(1)
      do 110 ip1 = 2, n1
        if (i .lt. n)  xip1 = x(ip1)
        f(i) = (three - two*xi)*xi + one
        if (i .gt. 1)  f(i) = f(i) - xim1
        if (i .lt. n)  f(i) = f(i) - two*xip1
        im1  = i
        i    = ip1
        xim1 = xi
        xi   = xip1
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 210  continue

      i = 1
      do 220 ip1 = 2, n1
        xi      = x(i)
        fj(i,i) = three - four*xi
        if (i .gt. 1)  fj(i,im1) = -one
        if (i .lt. n)  fj(i,ip1) = -two
        im1 = i
        i   = ip1
 220  continue

      return

 300  continue

      do 310 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 310  continue

      i  =   1
      xi = x(1)
      do 320 ip1 = 2, n1
        if (i .lt. n)  xip1 = x(ip1)
        f(i)    = (three - two*xi)*xi + one
        if (i .gt. 1)  f(i) = f(i) - xim1
        if (i .lt. n)  f(i) = f(i) - two*xip1
        fj(i,i) = three - four*xi
        if (i .gt. 1)  fj(i,im1) = -one
        if (i .lt. n)  fj(i,ip1) = -two
        im1  = i
        i    = ip1
        xim1 = xi
        xi   = xip1
 320  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 330 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 330  continue

      return

9999  format(/'1',70('=')//,
     *' broyden tridiagonal function (more et al.)'//,
     *'        number of variables =', i4,'  (variable)'/,
     *'        number of functions =', i4,'  (   = n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      double precision   zero, four
      parameter         (zero = 0.d0)
      parameter         (four = 4.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      nonzro(k) = 1

      dfj(k,k)  = -four

      return
      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            j

      double precision   zero, four
      parameter         (zero = 0.d0)
      parameter         (four = 4.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      hess(k,k) = -four

      return
      end
