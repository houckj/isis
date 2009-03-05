******************************************************************************
*               penalty function i
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
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

      double precision   a, xi

      double precision   ddot

      intrinsic          sqrt, dble

      double precision   a2
      common /PARAM1/    a2
      save   /PARAM1/

      double precision   zero, one, two, qtr
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (qtr = .25d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      na = mode / 1000
      nh = mode - na*1000
      nb = nh / 100
      nt = nh - nb*100
      nc = nt / 10
      nd = nt - nc*10

      lf = (na .ne. 0) .or. (nb .ne. 0) .or. (nd .ne. 0)
      lj = (nc .ne. 0) .or. (nd .ne. 0)

      if (lf .and. lj)  goto 300
      if (lf)           goto 100
      if (lj)           goto 200

*-----------------------------------------------------------------------

   10 continue

      nprobs = 1
      nstrts = 1

      n      = 4
      m      = n + 1

      a      = 1.d-5

      a2     = sqrt(a)

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      do 21 j = 1, n
        x(j) = dble(j)
   21 continue

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 2.24997d-5

      return

*-----------------------------------------------------------------------

 100  continue

      sum = - qtr
      do 110 i = 1, n
        xi   = x(i)
        sum  = sum + xi*xi
        f(i) = a2*(xi - one)
 110  continue
      f(n+1) = sum

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        call dcopy( n, zero, 0, fj( 1, j), 1)
        fj(j,j) = a2
        fj(m,j) = two*x(j)
 210  continue

      return

 300  continue

      sum = - qtr
      do 310 j = 1, n
        xj   = x(j)
        sum  = sum + xj*xj
        f(j) = a2*(xj - one)
        call dcopy( n, zero, 0, fj( 1, j), 1)
        fj(j,j) = a2
        fj(m,j) = two*xj
 310  continue
      f(n+1) = sum

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' penalty function i (more et al.) '//,
     *'        number of variables =', i4, '  (variable)'/,
     *'        number of functions =', i4, '  (  = n+1 )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( n, zero, 0, dfj( 1, j), 1)
  100 continue

      nonzro(k) = 1
      dfj(n+1,k) = two

      return
      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, lhess, hess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .true.
      if (k .ne. n+1)  return
      linear = .false.

      do 200 i = 1, n
        hess(i,i) = two
 200  continue

      return
      end
