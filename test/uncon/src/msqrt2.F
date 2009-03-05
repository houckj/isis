*******************************************************************************
*          2 by 2  matrix square root problem (hammarling)
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

      double precision   x1, x2, x3, x4

      double precision   ddot

      double precision   a
      common /PARAM1/    a(2,2)
      save   /PARAM1/

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .lt. -1)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)

      na = mode / 1000
      nh = mode - na*1000
      nb = nh / 100
      nt = nh - nb*100
      nc = nt / 10
      nd = nt - nc*10

      lf = (na .ne. 0) .or. (nb .ne. 0) .or. (nd .ne. 0)
      lj = (nc .ne. 0) .or. (nd .ne. 0)

      if (lf .and. lj)  goto 100
      if (lf)           goto 100
      if (lj)           goto 200

*-----------------------------------------------------------------------

  10  continue

      numbr  = 1
      nprobs = 1
      nstrts = 1

      na    = 2
      n     = 4
      m     = 4

      if (nout .gt. 0)  write( nout, 9999)  numbr, na, n

      a(1,1) = 1.0d-4
      a(1,2) = one
      a(2,1) = zero
      a(2,2) = 1.0d-4

      return

*-----------------------------------------------------------------------

  20  continue

      x(1) = 1.0d+0
      x(2) = 0.0d+0
      x(3) = 0.0d+0
      x(4) = 1.0d+0

      return

*-----------------------------------------------------------------------

  30  continue

      ftf = zero

      x(1) = 1.0d-2
      x(2) = 5.0d+1
      x(3) = 0.0d+0
      x(4) = 1.0d-2

      return

*-----------------------------------------------------------------------

 100  continue

      f(1) = (x1*x1 + x2*x3) - a(1,1)
      f(2) = (x1*x2 + x2*x4) - a(1,2)
      f(3) = (x3*x1 + x4*x3) - a(2,1)
      f(4) = (x3*x2 + x4*x4) - a(2,2)

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (.not. lj)   return

 200  continue

      fj( 1, 1) = two*x1
      fj( 1, 2) = x3
      fj( 1, 3) = x2
      fj( 1, 4) = zero
      fj( 2, 1) = x2
      fj( 2, 2) = x1 + x4
      fj( 2, 3) = zero
      fj( 2, 4) = x2
      fj( 3, 1) = x3
      fj( 3, 2) = zero
      fj( 3, 3) = x1 + x4
      fj( 3, 4) = x3
      fj( 4, 1) = zero
      fj( 4, 2) = x3
      fj( 4, 3) = x2
      fj( 4, 4) = two*x4

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' matrix square root problem (hammarling) ', i4//,
     *'                  rows in matrix = ', i4/,
     *'  # variables = # matrix entries = # constraints =', i4//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      do 10 j = 1, n
        nonzro(j) = 1
        call dcopy( m, zero, 0, dfj( 1, j), 1)
   10 continue

      goto ( 100, 200, 300, 400), k

 100  continue

      nonzro(4) = 0

      dfj( 1, 1) = two
      dfj( 2, 2) = one
      dfj( 3, 3) = one

      return

 200  continue

      nonzro(2) = 0

      dfj( 1, 3) = one
      dfj( 2, 1) = one
      dfj( 2, 4) = one
      dfj( 4, 3) = one

      return

 300  continue

      nonzro(3) = 0

      dfj( 1, 2) = one
      dfj( 3, 1) = one
      dfj( 3, 4) = one
      dfj( 4, 2) = one

      return

 400  continue

      nonzro(1) = 0

      dfj( 2, 2) = one
      dfj( 3, 3) = one
      dfj( 4, 4) = two

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            j

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      do  10 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
   10 continue

      linear = .false.

      goto ( 100, 200, 300, 400), k

 100  continue

      hess(1,1) = two
      hess(2,3) = one
      hess(3,2) = one

      return

 200  continue

      hess(1,2) = one
      hess(2,1) = one
      hess(2,4) = one
      hess(4,2) = one

      return

 300  continue

      hess(1,3) = one
      hess(3,1) = one
      hess(3,4) = one
      hess(4,3) = one

      return

 400  continue

      hess(2,3) = one
      hess(3,2) = one
      hess(4,4) = two

      return

      end
