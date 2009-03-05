*******************************************************************************
*          3 by 3 matrix square root problem (hammarling)
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

      double precision   x1, x2, x3, x4, x5, x6, x7, x8, x9

      double precision   ddot

      double precision   a
      common /PARAM1/    a(3,3)
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
      x5 = x(5)
      x6 = x(6)
      x7 = x(7)
      x8 = x(8)
      x9 = x(9)

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

   10 continue

      numbr  = 2
      nprobs = 3
      nstrts = 1

      na     = 3
      n      = 9
      m      = 9

      if (nout .gt. 0)  write( nout, 9999)  numbr, na, n

      goto ( 11, 12, 13), nprob

  11  continue

      a(1,1) = 1.0d-4
      a(1,2) = one
      a(1,3) = zero
      a(2,1) = zero
      a(2,2) = 1.0d-4
      a(2,3) = zero
      a(3,1) = zero
      a(3,2) = zero
      a(3,3) = 1.0d-4

      return

  12  continue

      a(1,1) = zero
      a(1,2) = one
      a(1,3) = zero
      a(2,1) = zero
      a(2,2) = zero
      a(2,3) = zero
      a(3,1) = zero
      a(3,2) = zero
      a(3,3) = zero

      return

  13  continue

      a(1,1) = one
      a(1,2) = one
      a(1,3) = one
      a(2,1) = zero
      a(2,2) = zero
      a(2,3) = zero
      a(3,1) = zero
      a(3,2) = zero
      a(3,3) = zero

      return

*-----------------------------------------------------------------------

  20  continue

      goto ( 21, 22, 23), nprob

  21  continue

      x(1) =  1.0d+0
      x(2) =  0.0d+0
      x(3) =  0.0d+0
      x(4) =  0.0d+0
      x(5) =  1.0d+0
      x(6) =  0.0d+0
      x(7) =  0.0d+0
      x(8) =  0.0d+0
      x(9) =  1.0d+0

      return

  22  continue

      x(1) = 1.0d+0
      x(2) = 0.0d+0
      x(3) = 0.0d+0
      x(4) = 0.0d+0
      x(5) = 1.0d+0
      x(6) = 0.0d+0
      x(7) = 0.0d+0
      x(8) = 0.0d+0
      x(9) = 1.0d+0

      return

  23  continue

      x(1) = 1.0d+0
      x(2) = 0.0d+0
      x(3) = 0.0d+0
      x(4) = 0.0d+0
      x(5) = 1.0d+0
      x(6) = 0.0d+0
      x(7) = 0.0d+0
      x(8) = 0.0d+0
      x(9) = 1.0d+0

      return

*-----------------------------------------------------------------------

  30  continue

      ftf = zero

      goto ( 31, 32, 33), nprob

  31  continue

      x(1) =  0.01d+0
      x(2) = 50.0d+0
      x(3) =  0.0d+0
      x(4) =  0.0d+0
      x(5) =  0.01d+0
      x(6) =  0.0d+0
      x(7) =  0.0d+0
      x(8) =  0.0d+0
      x(9) =  0.01d+0

      return

  32  continue

      x(1) = 0.0d+0
      x(2) = 0.0d+0
      x(3) = 1.0d+0
      x(4) = 0.0d+0
      x(5) = 0.0d+0
      x(6) = 0.0d+0
      x(7) = 0.0d+0
      x(8) = 1.0d+0
      x(9) = 0.0d+0

      return

  33  continue

      x(1) =  1.d+0
      x(2) =  1.d+0
      x(3) =  1.d+0
      x(4) =  0.d+0
      x(5) =  0.d+0
      x(6) =  0.d+0
      x(7) =  0.d+0
      x(8) =  0.d+0
      x(9) =  0.d+0

      return

*-----------------------------------------------------------------------

 100  continue

      f(1) = (x1*x1 + x2*x4 + x3*x7) - a(1,1)
      f(2) = (x1*x2 + x2*x5 + x3*x8) - a(1,2)
      f(3) = (x1*x3 + x2*x6 + x3*x9) - a(1,3)
      f(4) = (x4*x1 + x5*x4 + x6*x7) - a(2,1)
      f(5) = (x4*x2 + x5*x5 + x6*x8) - a(2,2)
      f(6) = (x4*x3 + x5*x6 + x6*x9) - a(2,3)
      f(7) = (x7*x1 + x8*x4 + x9*x7) - a(3,1)
      f(8) = (x7*x2 + x8*x5 + x9*x8) - a(3,2)
      f(9) = (x7*x3 + x8*x6 + x9*x9) - a(3,3)

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (.not. lj)  return

 200  continue

      fj( 1, 1) = two * x1
      fj( 1, 2) = x4
      fj( 1, 3) = x7
      fj( 1, 4) = x2
      fj( 1, 5) = zero
      fj( 1, 6) = zero
      fj( 1, 7) = x3
      fj( 1, 8) = zero
      fj( 1, 9) = zero
      fj( 2, 1) = x2
      fj( 2, 2) = x1 + x5
      fj( 2, 3) = x8
      fj( 2, 4) = zero
      fj( 2, 5) = x2
      fj( 2, 6) = zero
      fj( 2, 7) = zero
      fj( 2, 8) = x3
      fj( 2, 9) = zero
      fj( 3, 1) = x3
      fj( 3, 2) = x6
      fj( 3, 3) = x1 + x9
      fj( 3, 4) = zero
      fj( 3, 5) = zero
      fj( 3, 6) = x2
      fj( 3, 7) = zero
      fj( 3, 8) = zero
      fj( 3, 9) = x3
      fj( 4, 1) = x4
      fj( 4, 2) = zero
      fj( 4, 3) = zero
      fj( 4, 4) = x1 + x5
      fj( 4, 5) = x4
      fj( 4, 6) = x7
      fj( 4, 7) = x6
      fj( 4, 8) = zero
      fj( 4, 9) = zero
      fj( 5, 1) = zero
      fj( 5, 2) = x4
      fj( 5, 3) = zero
      fj( 5, 4) = x2
      fj( 5, 5) = two*x5
      fj( 5, 6) = x8
      fj( 5, 7) = zero
      fj( 5, 8) = x6
      fj( 5, 9) = zero
      fj( 6, 1) = zero
      fj( 6, 2) = zero
      fj( 6, 3) = x4
      fj( 6, 4) = x3
      fj( 6, 5) = x6
      fj( 6, 6) = x5 + x9
      fj( 6, 7) = zero
      fj( 6, 8) = zero
      fj( 6, 9) = x6
      fj( 7, 1) = x7
      fj( 7, 2) = zero
      fj( 7, 3) = zero
      fj( 7, 4) = x8
      fj( 7, 5) = zero
      fj( 7, 6) = zero
      fj( 7, 7) = x1 + x9
      fj( 7, 8) = x4
      fj( 7, 9) = x7
      fj( 8, 1) = zero
      fj( 8, 2) = x7
      fj( 8, 3) = zero
      fj( 8, 4) = zero
      fj( 8, 5) = x8
      fj( 8, 6) = zero
      fj( 8, 7) = x2
      fj( 8, 8) = x5 + x9
      fj( 8, 9) = x8
      fj( 9, 1) = zero
      fj( 9, 2) = zero
      fj( 9, 3) = x7
      fj( 9, 4) = zero
      fj( 9, 5) = zero
      fj( 9, 6) = x8
      fj( 9, 7) = x3
      fj( 9, 8) = x6
      fj( 9, 9) = two*x9

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

      goto ( 100, 200, 300, 400, 500, 600, 700, 800, 900), k

 100  continue

      nonzro(5) = 0
      nonzro(6) = 0
      nonzro(8) = 0
      nonzro(9) = 0

      dfj( 1, 1) = two
      dfj( 2, 2) = one
      dfj( 3, 3) = one
      dfj( 4, 4) = one
      dfj( 7, 7) = one

      return

 200  continue

      nonzro(2) = 0
      nonzro(3) = 0
      nonzro(8) = 0
      nonzro(9) = 0

      dfj( 1, 4) = one
      dfj( 2, 1) = one
      dfj( 2, 5) = one
      dfj( 3, 6) = one
      dfj( 5, 4) = one
      dfj( 8, 7) = one

      return

 300  continue

      nonzro(2) = 0
      nonzro(3) = 0
      nonzro(5) = 0
      nonzro(6) = 0

      dfj( 1, 7) = one
      dfj( 2, 8) = one
      dfj( 3, 1) = one
      dfj( 3, 9) = one
      dfj( 6, 4) = one
      dfj( 9, 7) = one

      return

 400  continue

      nonzro(4) = 0
      nonzro(6) = 0
      nonzro(7) = 0
      nonzro(9) = 0

      dfj( 1, 2) = one
      dfj( 4, 1) = one
      dfj( 4, 5) = one
      dfj( 5, 2) = one
      dfj( 6, 3) = one
      dfj( 7, 8) = one

      return

 500  continue

      nonzro(1) = 0
      nonzro(3) = 0
      nonzro(7) = 0
      nonzro(9) = 0

      dfj( 2, 2) = one
      dfj( 4, 4) = one
      dfj( 5, 5) = two
      dfj( 6, 6) = one
      dfj( 8, 8) = one

      return

 600  continue

      nonzro(1) = 0
      nonzro(3) = 0
      nonzro(4) = 0
      nonzro(6) = 0

      dfj( 3, 2) = one
      dfj( 4, 7) = one
      dfj( 5, 8) = one
      dfj( 6, 5) = one
      dfj( 6, 9) = one
      dfj( 9, 8) = one

      return

 700  continue

      nonzro(4) = 0
      nonzro(5) = 0
      nonzro(7) = 0
      nonzro(8) = 0

      dfj( 1, 3) = one
      dfj( 4, 6) = one
      dfj( 7, 1) = one
      dfj( 7, 9) = one
      dfj( 8, 2) = one
      dfj( 9, 3) = one

      return

 800  continue

      nonzro(1) = 0
      nonzro(2) = 0
      nonzro(7) = 0
      nonzro(8) = 0

      dfj( 2, 3) = one
      dfj( 5, 6) = one
      dfj( 7, 4) = one
      dfj( 8, 5) = one
      dfj( 8, 9) = one
      dfj( 9, 6) = one

      return

 900  continue

      nonzro(1) = 0
      nonzro(2) = 0
      nonzro(4) = 0
      nonzro(5) = 0

      dfj( 3, 3) = one
      dfj( 6, 6) = one
      dfj( 7, 7) = one
      dfj( 8, 8) = one
      dfj( 9, 9) = two

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

      do 10 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
   10 continue

      linear = .false.

      goto ( 100, 200, 300, 400, 500, 600, 700, 800, 900), k

 100  continue

      hess(1,1) = two
      hess(2,4) = one
      hess(4,2) = one
      hess(3,7) = one
      hess(7,3) = one

      return

 200  continue

      hess(1,2) = one
      hess(2,1) = one
      hess(2,5) = one
      hess(5,2) = one
      hess(3,8) = one
      hess(8,3) = one

      return

 300  continue

      hess(1,3) = one
      hess(3,1) = one
      hess(2,6) = one
      hess(6,2) = one
      hess(3,9) = one
      hess(9,3) = one

      return

 400  continue

      hess(1,4) = one
      hess(4,1) = one
      hess(4,5) = one
      hess(5,4) = one
      hess(6,7) = one
      hess(7,6) = one

      return

 500  continue

      hess(2,4) = one
      hess(4,2) = one
      hess(5,5) = two
      hess(6,8) = one
      hess(8,6) = one

      return

 600  continue

      hess(3,4) = one
      hess(4,3) = one
      hess(5,6) = one
      hess(6,5) = one
      hess(6,9) = one
      hess(9,6) = one

      return

 700  continue

      hess(1,7) = one
      hess(7,1) = one
      hess(4,8) = one
      hess(8,4) = one
      hess(7,9) = one
      hess(9,7) = one

      return

 800  continue

      hess(2,7) = one
      hess(7,2) = one
      hess(5,8) = one
      hess(8,5) = one
      hess(8,9) = one
      hess(9,8) = one

      return

 900  continue

      hess(3,7) = one
      hess(7,3) = one
      hess(6,8) = one
      hess(8,6) = one
      hess(9,9) = two

      return

      end


