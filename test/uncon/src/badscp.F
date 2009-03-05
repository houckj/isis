**************************************************************************
*               powell badly scaled function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
**************************************************************************

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

      double precision   e1, e2, x1, x2

      double precision   ddot

      intrinsic          exp

      double precision   zero, two, ten4th
      parameter         (zero = 0.d0, two = 2.d0, ten4th = 1.d4)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      e1 = exp(-x1)
      e2 = exp(-x2)

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
      m      = 2

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  0.d0
      x(2) =  1.d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  1.098d-5
      x(2) =  9.106d0

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      f(1) =  ten4th*x1*x2 - one
      f(2) =  e1 + e2 - 1.0001d0

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      fj( 1, 1) = ten4th*x2
      fj( 2, 1) = -e1
      fj( 1, 2) = ten4th*x1
      fj( 2, 2) = -e2

      return

 300  continue

      f(1) =  ten4th*x1*x2 - one
      f(2) =  e1 + e2 - 1.0001d0

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      fj( 1, 1) = ten4th*x2
      fj( 2, 1) = -e1
      fj( 1, 2) = ten4th*x1
      fj( 2, 2) = -e2

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' powell badly scaled function (more et al.)'//,
     *'        number of variables =', i4, '  (2)'/,
     *'        number of functions =', i4, '  (2)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      intrinsic          exp

      double precision   zero, ten4th
      parameter         (zero = 0.d0, ten4th = 1.d4)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 1
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      goto ( 210, 220), k

  210 continue

      dfj( 1, 2) = ten4th
      dfj( 2, 1) = exp(-x(1))

      return

  220 continue

      dfj( 1, 1) = ten4th
      dfj( 2, 2) = exp(-x(2))

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

      double precision   zero, ten4th
      parameter         (zero = 0.d0, ten4th = 1.d4)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      goto ( 210, 220), k

  210 continue

      hess( 1, 2) = ten4th
      hess( 2, 1) = ten4th

      return

  220 continue

      hess( 1, 1) = exp(-x(1))
      hess( 2, 2) = exp(-x(2))

      return

      end
