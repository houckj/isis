*****************************************************************************
*               brown badly scaled function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
*****************************************************************************

      subroutine getfun( x, n, f, m, ftf, fj, lfj, g, mode)

      implicit double precision (a-h,o-z)

      integer            n, m, lfj, mode

      double precision   x(n), f(m), ftf, fj(lfj,n), g(m)

      integer            nprob, nprobs, nstart, nstrts
      common /PROBLM/    nprob, nprobs, nstart, nstrts

      integer            nout
      common /IOUNIT/    nout

      logical            lf, lj

      integer            na, nb, nc, nd, nt, nh

      integer            j

      double precision   x1, x2

      double precision   ddot

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

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
      m      = 3

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  1.d0
      x(2) =  1.d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  1.d6
      x(2) =  2.d-6

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      f(1) =  x1 - 1.d6
      f(2) =  x2 - 2.d-6
      f(3) =  x1*x2 - two

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      fj( 1, 1) = one
      fj( 1, 2) = zero
      fj( 2, 1) = one
      fj( 2, 2) = one
      fj( 3, 1) = x2
      fj( 3, 2) = x1

      return

 300  continue

      f(1) =  x1 - 1.d6
      f(2) =  x2 - 2.d-6
      f(3) =  x1*x2 - two

      fj( 1, 1) = one
      fj( 1, 2) = zero
      fj( 2, 1) = zero
      fj( 2, 2) = one
      fj( 3, 1) = x2
      fj( 3, 2) = x1

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 310 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 310  continue

      return

9999  format(/'1',70('=')//,
     *' brown badly scaled function (more et al.)'//,
     *'        number of variables =', i4, '  (2)'/,
     *'        number of functions =', i4, '  (3)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

*=======================================================================

      goto ( 210, 220), k

 210  continue

      nonzro(1) = 0
      nonzro(2) = 1

      dfj( 1, 1) = zero
      dfj( 1, 2) = zero
      dfj( 2, 1) = zero
      dfj( 2, 2) = zero
      dfj( 3, 1) = zero
      dfj( 3, 2) = one

      return

 220  continue

      nonzro(1) = 1
      nonzro(2) = 0

      dfj( 1, 1) = zero
      dfj( 1, 2) = zero
      dfj( 2, 1) = zero
      dfj( 2, 2) = zero
      dfj( 3, 1) = one
      dfj( 3, 2) = zero

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

*=======================================================================

      goto ( 210, 220, 230), k

 210  continue

      linear = .true.

      hess( 1, 1) = zero
      hess( 1, 2) = zero
      hess( 2, 1) = zero
      hess( 2, 2) = zero

      return

 220  continue

      linear = .true.

      hess( 1, 1) = zero
      hess( 1, 2) = zero
      hess( 2, 1) = zero
      hess( 2, 2) = zero

      return

 230  continue

      linear = .false.

      hess( 1, 1) = zero
      hess( 1, 2) = one
      hess( 2, 1) = one
      hess( 2, 2) = zero

      return

      end
