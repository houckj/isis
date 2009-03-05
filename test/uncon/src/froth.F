****************************************************************************
*               freudenstein and roth function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
****************************************************************************

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

      double precision   x1, x2

      double precision   ddot

      double precision   zero, one, two, three, five, ten
      double precision   thirtn, fourtn, twnine
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0, five = 5.d0, ten = 10.d0)
      parameter         (thirtn = 13.d0, fourtn = 14.d0, twnine = 29.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)

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

      n      = 2
      m      = 2

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  0.5d0
      x(2) = -2.0d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  5.d0
      x(2) =  4.d0

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      f(1) = -thirtn + x1 + ((five - x2)*x2 -    two)*x2
      f(2) = -twnine + x1 + (( x2 + one)*x2 - fourtn)*x2

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      fj( 1, 1) = one
      fj( 2, 1) = one
      fj( 1, 2) = (ten - three*x2)*x2 - two
      fj( 2, 2) = (three*x2 + two)*x2 - fourtn

      return

 300  continue

      f(1) = -thirtn + x1 + ((five - x2)*x2 -    two)*x2
      f(2) = -twnine + x1 + (( x2 + one)*x2 - fourtn)*x2

      fj( 1, 1) = one
      fj( 2, 1) = one
      fj( 1, 2) = (ten - three*x2)*x2 - two
      fj( 2, 2) = (three*x2 + two)*x2 - fourtn

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' freudenstein and roth function (more et al.)'//,
     *'        number of variables =', i4, '  (2)'/,
     *'        number of functions =', i4, '  (2)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      double precision   s2

      double precision   zero, two, six, ten
      parameter         (zero = 0.d0, two = 2.d0)
      parameter         (six = 6.d0, ten = 1.d1)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      if (k .eq. 1)  return

      s2 = six*x(2)

      nonzro(2) = 1

      dfj( 1, 2) = ten - s2
      dfj( 2, 2) = s2 + two

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, lhess, hess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            j

      double precision   zero, two, six, ten
      parameter         (zero = 0.d0, two = 2.d0)
      parameter         (six = 6.d0, ten = 10.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      x2 = x(2)
      goto ( 210, 220), k

  210 continue

      hess( 2, 2) = ten - six*x2

      return

  220 continue

      hess( 2, 2) = six*x2 + two

      return

      end


