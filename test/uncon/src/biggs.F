*****************************************************************************
*               biggs exp6 function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
*****************************************************************************

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

      double precision   x1, x2, x3, x4, x5, x6
      double precision   e1, e2, e5, ti, ym, yp

      double precision   ddot

      intrinsic          dble, exp

      double precision   zero, one, point1, three, four, five, ten
      parameter         (zero = 0.d0, one = 1.d0, point1 = .1d0)
      parameter         (three = 3.d0, four = 4.d0, five = 5.d0)
      parameter         (ten = 10.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)
      x5 = x(5)
      x6 = x(6)

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

      n      =  6
      m      = 13

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) = 1.d0
      x(2) = 2.d0
      x(3) = 1.d0
      x(4) = 1.d0
      x(5) = 1.d0
      x(6) = 1.d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  1.d0
      x(2) =  10.d0
      x(3) =  1.d0
      x(4) =  5.d0
      x(5) =  4.d0
      x(6) =  3.d0

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ti = point1*dble(i)
        yp = exp(-ti) + three*exp(-four*ti)
        ym = five*exp(-ten*ti)
        e1 = exp(-ti*x1)
        e2 = exp(-ti*x2)
        e5 = exp(-ti*x5)
        f(i) = (x3*e1 + x6*e5  + ym) - (x4*e2 + yp)
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        ti = point1*dble(i)
        e1 = exp(-ti*x1)
        e2 = exp(-ti*x2)
        e5 = exp(-ti*x5)
        fj( i, 1) = -ti*x3*e1
        fj( i, 2) =  ti*x4*e2
        fj( i, 3) =  e1
        fj( i, 4) = -e2
        fj( i, 5) = -ti*x6*e5
        fj( i, 6) =  e5
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ti = point1*dble(i)
        yp = exp(-ti) + three*exp(-four*ti)
        ym = five*exp(-ten*ti)
        e1 = exp(-ti*x1)
        e2 = exp(-ti*x2)
        e5 = exp(-ti*x5)
        f(i) = (x3*e1 + x6*e5 + ym) - (x4*e2 + yp)
        fj( i, 1) = -ti*x3*e1
        fj( i, 2) =  ti*x4*e2
        fj( i, 3) =  e1
        fj( i, 4) = -e2
        fj( i, 5) = -ti*x6*e5
        fj( i, 6) =  e5
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 320 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 320  continue

      return

9999  format(/'1',70('=')//,
     *' biggs exp6 function (more et al.) '//,
     *'        number of variables =', i4, '  (    6   )'/,
     *'        number of functions =', i4, '  (  >= 6  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, j

      double precision   x1, x2, x3, x4, x5, x6
      double precision   e1, e2, e5, ti, t1, t2, t5

      intrinsic          dble, exp

      double precision   zero, point1
      parameter         (zero = 0.d0, point1 = .1d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      goto ( 210, 220, 230, 240, 250, 260 ), k

 230  continue

      x1 = x(1)

      nonzro(1) = 1

      do 235 i = 1, m
        ti = point1*dble(i)
        dfj( i, 1) = -ti*exp(-ti*x1)
 235  continue

      return

 240  continue

      x2 = x(2)

      nonzro(2) = 1

      do 245 i = 1, m
        ti = point1*dble(i)
        dfj( i, 2) =  ti*exp(-ti*x2)
 245  continue

      return

 260  continue

      x5 = x(5)

      nonzro(5) = 1

      do 265 i = 1, m
        ti = point1*dble(i)
        dfj( i, 5) = -ti*exp(-ti*x5)
 265  continue

      return

 210  continue

      x1 = x(1)
      x3 = x(3)

      nonzro(1) = 1
      nonzro(3) = 1

      do 215 i = 1, m
        ti = point1*dble(i)
        e1 = exp(-ti*x1)
        t1 = ti*e1
        dfj( i, 1) =  ti*t1*x3
        dfj( i, 3) = -t1
 215  continue

      return

 220  continue

      x2 = x(2)
      x4 = x(4)

      nonzro(2) = 1
      nonzro(4) = 1

      do 225 i = 1, m
        ti = point1*dble(i)
        e2 = exp(-ti*x2)
        t2 = ti*e2
        dfj( i, 2) = -ti*t2*x4
        dfj( i, 4) =  t2
 225  continue

      return

 250  continue

      x5 = x(5)
      x6 = x(6)

      nonzro(5) = 1
      nonzro(6) = 1

      do 255 i = 1, m
        ti = point1*dble(i)
        e5 = exp(-ti*x5)
        t5 = ti*e5
        dfj( i, 5) =  ti*t5*x6
        dfj( i, 6) = -t5
 255  continue

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

      double precision   e1, e2, e5, tk, t1, t2, t5
       
      intrinsic          dble, exp

      double precision   zero, point1
      parameter         (zero = 0.d0, point1 = .1d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      tk = point1*dble(k)
      e1 = exp(-tk*x(1))
      e2 = exp(-tk*x(2))
      e5 = exp(-tk*x(5))
      t1 = tk * e1
      t2 = tk * e2
      t5 = tk * e5
      hess(1,3) = -t1
      hess(3,1) = hess(3,1)
      hess(2,4) =  t2
      hess(4,2) = hess(2,4)
      hess(5,6) = -t5
      hess(6,5) = hess(5,6)
      hess(1,1) =  tk * x(3) * t1
      hess(2,2) = -tk * x(4) * t2
      hess(5,5) =  tk * x(6) * t5

      return
      end
