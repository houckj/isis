*************************************************************************
*               gaussian function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
*************************************************************************

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

      double precision   x1, x2, x3
      double precision   ei, si, ti, tt, t3

      double precision   ddot

      intrinsic          dble, exp

      double precision   y
      common /PARAM1/    y(15)
      save   /PARAM1/

      double precision   two, eight
      parameter         (two = 2.d0, eight = 8.d0)

*========================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

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

      n      =  3
      m      = 15

      y( 1) = 0.0009d0
      y( 2) = 0.0044d0
      y( 3) = 0.0175d0
      y( 4) = 0.0540d0
      y( 5) = 0.1295d0
      y( 6) = 0.2420d0
      y( 7) = 0.3521d0
      y( 8) = 0.3989d0
      y( 9) = 0.3521d0
      y(10) = 0.2420d0
      y(11) = 0.1295d0
      y(12) = 0.0540d0
      y(13) = 0.0175d0
      y(14) = 0.0044d0
      y(15) = 0.0009d0

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  0.4d0
      x(2) =  1.d0
      x(3) =  0.d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 1.12793d-8

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ti   = (eight - dble(i)) / two
        f(i) = x1*exp(-(x2*(ti-x3)**2)/two) - y(i)
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        ti = (eight - dble(i)) / two
        t3 = ti - x3
        tt = (t3*t3)/two
        ei = exp(-x2*tt)
        si = x1*ei
        fj(i,1) =  ei
        fj(i,2) = -tt*si
        fj(i,3) =  x2*t3*si
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ti = (eight - dble(i)) / two
        t3 = ti - x3
        tt = (t3*t3)/two
        ei = exp(-x2*tt)
        si = x1*ei
        f(i) = x1*exp((-x2*(ti-x3)**2)/two) - y(i)
        fj(i,1) =  ei
        fj(i,2) = -tt*si
        fj(i,3) =  x2*t3*si
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' gaussian function (more et al.)'//,
     *'        number of variables =', i4, '  ( 3)'/,
     *'        number of functions =', i4, '  (15)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i

      double precision   x1, x2, x3
      double precision   ei, si, ti, tt, t2, t3

      intrinsic          dble, exp

      double precision   zero, one, two, eight
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (eight = 8.d0)

*=======================================================================

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      goto ( 100, 200, 300), k

 100  continue

      nonzro(1) = 0
      nonzro(2) = 1
      nonzro(3) = 1

      do 110 i = 1, m
        ti = (eight - dble(i)) / two
        t3 = ti - x3
        tt = (t3*t3)/two
        ei = exp(-x2*tt)
        dfj(i,1) =  zero
        dfj(i,2) = -tt*ei
        dfj(i,3) =  x2*t3*ei
 110  continue

      return

 200  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 1

      do 210 i = 1, m
        ti = (eight - dble(i)) / two
        t3 =  ti - x3
        tt =  (t3*t3)/two
        t2 = -x2*tt
        ei =  exp(t2)
        si =  x1*ei
        dfj(i,1) = -tt*ei
        dfj(i,2) =  tt*tt*si
        dfj(i,3) =  t3*(one + t2)*si
 210  continue

      return

 300  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 1

      do 310 i = 1, m
        ti = (eight - dble(i)) / two
        t3 =  ti - x3
        tt =  (t3*t3)/two
        t2 = -x2*tt
        ei =  exp(t2)
        si =  x1*ei
        dfj(i,1) =  x2*t3*ei
        dfj(i,2) =  t3*(one + t2)*si
        dfj(i,3) = -x2*(two*t2 + one)*si
 310  continue

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear
 
      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   x1, x2, x3
      double precision   ek, s1, s2, tk, tt, t2, t3

      intrinsic          dble, exp

      double precision   zero, one, two, eight
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (eight = 8.d0)

*=======================================================================

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      linear = .false.

      tk = (eight - dble(k))/two
      t3 = tk - x3
      tt = (t3*t3)/two
      t2 = -x2*tt
      ek = exp(t2)
      s1 = x1*ek
      s2 = x2*ek

      hess( 1, 1) =  zero
      hess( 1, 2) = -tt*ek
      hess( 2, 1) =  hess( 1, 2)
      hess( 2, 2) =  (tt*tt)*s1
      hess( 1, 3) =  t3*s2
      hess( 3, 1) =  hess( 1, 3)
      hess( 2, 3) =  t3*(one + t2)*s1
      hess( 3, 2) =  hess( 2, 3)
      hess( 3, 3) = -x1*(two*t2 + one)*s2

      return
      end

