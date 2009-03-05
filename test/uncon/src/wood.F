***************************************************************************
*               wood function
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

      integer            j 

      double precision   x1, x2, x3, x4

      double precision   ddot

      double precision   zero, one, two, ten, twenty
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (ten = 10.d0, twenty = 20.d0)

      double precision   rt10, r10inv, rt90
      common /PARAM1/    rt10, r10inv, rt90
      save   /PARAM1/

c     parameter         (rt10   = 3.1622776601683795d0)
c     parameter         (r10inv = .31622776601683794d0)
c     parameter         (rt90   = 9.4868329805051381d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)

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

      n      = 4
      m      = 6

      rt10    = sqrt(10.d0)
      r10inv  = 1.d0 / rt10
      rt90    = sqrt(90.d0)

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) = -3.d0
      x(2) = -1.d0
      x(3) = -3.d0
      x(4) = -1.d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  1.d0
      x(2) =  1.d0
      x(3) =  1.d0
      x(4) =  1.d0

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      f(1) = ten*(x2 - x1*x1)
      f(2) = one - x1
      f(3) = rt90*(x4 - x3*x3)
      f(4) = one - x3
      f(5) = rt10*(x2 + x4 - two)
      f(6) = r10inv*(x2 - x4)

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      fj( 1, 1) = -twenty*x1
      fj( 1, 2) =  ten
      fj( 1, 3) =  zero
      fj( 1, 4) =  zero
      fj( 2, 1) = -one
      fj( 2, 2) =  zero
      fj( 2, 3) =  zero
      fj( 2, 4) =  zero
      fj( 3, 1) =  zero
      fj( 3, 2) =  zero
      fj( 3, 3) = -two*rt90*x3
      fj( 3, 4) =  rt90
      fj( 4, 1) =  zero
      fj( 4, 2) =  zero
      fj( 4, 3) = -one
      fj( 4, 4) =  zero
      fj( 5, 1) =  zero
      fj( 5, 2) =  rt10
      fj( 5, 3) =  zero
      fj( 5, 4) =  rt10
      fj( 6, 1) =  zero
      fj( 6, 2) =  r10inv
      fj( 6, 3) =  zero
      fj( 6, 4) = -r10inv

      return

 300  continue

      f(1) = ten*(x2 - x1*x1)
      f(2) = one - x1
      f(3) = rt90*(x4 - x3*x3)
      f(4) = one - x3
      f(5) = rt10*(x2 + x4 - two)
      f(6) = r10inv*(x2 - x4)

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      fj( 1, 1) = -twenty*x1
      fj( 1, 2) =  ten
      fj( 1, 3) =  zero
      fj( 1, 4) =  zero
      fj( 2, 1) = -one
      fj( 2, 2) =  zero
      fj( 2, 3) =  zero
      fj( 2, 4) =  zero
      fj( 3, 1) =  zero
      fj( 3, 2) =  zero
      fj( 3, 3) = -two*rt90*x3
      fj( 3, 4) =  rt90
      fj( 4, 1) =  zero
      fj( 4, 2) =  zero
      fj( 4, 3) = -one
      fj( 4, 4) =  zero
      fj( 5, 1) =  zero
      fj( 5, 2) =  rt10
      fj( 5, 3) =  zero
      fj( 5, 4) =  rt10
      fj( 6, 1) =  zero
      fj( 6, 2) =  r10inv
      fj( 6, 3) =  zero
      fj( 6, 4) = -r10inv

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' wood function (more et al.)'//,
     *'        number of variables =', i4,'  (4)'/,
     *'        number of functions =', i4,'  (6)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      double precision   rt10, r10inv, rt90
      common /PARAM1/    rt10, r10inv, rt90
      save   /PARAM1/

c     double precision   zero, rt90
c     parameter         (zero = 0.d0, rt90 = 9.4868329805051d0)
 
      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      if (k .ne. 0)  goto 210
      if (k .eq. 3)  goto 230

      return

  210 continue

      nonzro(1) = 1
      dfj( 1, 1) = -20.d0

      return

  230 continue

      nonzro(3) = 1
      dfj( 3, 3) = 2.d0*rt90

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

      double precision   rt10, r10inv, rt90
      common /PARAM1/    rt10, r10inv, rt90
      save   /PARAM1/

c     double precision   zero, rt90
c     parameter         (zero = 0.d0, rt90 = 9.4868329805051d0)
 
      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue
c
      goto ( 210, 220, 230, 220, 220, 220), k

  210 continue

      linear = .false.
      hess( 1, 1) = -20.d0

      return

  220 continue

      linear = .true.

      return

  230 continue

      linear = .false.
      hess( 3, 3) = -2.d0*rt90

      return

      end
