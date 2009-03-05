****************************************************************************
*               osborne 1 function
*  more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
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

      integer            i, im1, j

      double precision   x1, x2, x3, x4, x5
      double precision   e4, e5, t2, t3, ti

      double precision   ddot

      intrinsic          dble, exp

      double precision   y
      common /PARAM1/    y(33)
      save   /PARAM1/

      double precision   zero, one, ten
      parameter         (zero = 0.d0, one = 1.d0, ten = 10.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)
      x5 = x(5)

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

      n      =  5
      m      = 33

      y( 1) = 0.844d0
      y( 2) = 0.908d0
      y( 3) = 0.932d0
      y( 4) = 0.936d0
      y( 5) = 0.925d0
      y( 6) = 0.908d0
      y( 7) = 0.881d0
      y( 8) = 0.850d0
      y( 9) = 0.818d0
      y(10) = 0.784d0
      y(11) = 0.751d0
      y(12) = 0.718d0
      y(13) = 0.685d0
      y(14) = 0.658d0
      y(15) = 0.628d0
      y(16) = 0.603d0
      y(17) = 0.580d0
      y(18) = 0.558d0
      y(19) = 0.538d0
      y(20) = 0.522d0
      y(21) = 0.506d0
      y(22) = 0.490d0
      y(23) = 0.478d0
      y(24) = 0.467d0
      y(25) = 0.457d0
      y(26) = 0.448d0
      y(27) = 0.438d0
      y(28) = 0.431d0
      y(29) = 0.424d0
      y(30) = 0.420d0
      y(31) = 0.414d0
      y(32) = 0.411d0
      y(33) = 0.406d0

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  0.5d0
      x(2) =  1.5d0
      x(3) = -1.0d0
      x(4) =  0.01d0
      x(5) =  0.02d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  0.3754d0
      x(2) =  1.9358d0
      x(3) = -1.4647d0
      x(4) =  0.01287d0
      x(5) =  0.02212d0

      ftf = 5.46489d-5

      return

*-----------------------------------------------------------------------

 100  continue

      im1 = 0
      do 110 i = 1, m
        ti   =  dble(im1)*ten
        e4   =  exp(-ti*x4)
        e5   =  exp(-ti*x5)
        f(i) = (x1 + x2*e4 + x3*e5) - y(i)
        im1  =  i
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      im1 = 0
      do 210 i = 1, m
        ti = dble(im1)*ten
        e4 = exp(-ti*x4)
        e5 = exp(-ti*x5)
        fj( i, 1) = one
        fj( i, 2) =  e4
        fj( i, 3) =  e5
        fj( i, 4) = -ti*x2*e4
        fj( i, 5) = -ti*x3*e5
        im1 = i
 210  continue

      return

 300  continue

      im1 = 0
      do 310 i = 1, m
        ti = dble(im1)*ten
        e4 = exp(-ti*x4)
        e5 = exp(-ti*x5)
        t2 = x2*e4
        t3 = x3*e5
        f(i) = (x1 + t2 + t3) - y(i)
        fj( i, 1) = one
        fj( i, 2) =  e4
        fj( i, 3) =  e5
        fj( i, 4) = -ti*t2
        fj( i, 5) = -ti*t3
        im1 = i
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' osborne 1 function (more et al.) - exponential fitting'//,
     *'        number of variables =', i4, '  ( 5)'/,
     *'        number of functions =', i4, '  (33)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, im1, j

      double precision   x2, x3, x4, x5, ti

      intrinsic          dble, exp

      double precision   zero, ten
      parameter         (zero = 0.d0, ten = 10.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      goto ( 210, 220, 230, 240, 250 ), k

 210  continue

      return

 220  continue

      x4 = x(4)

      nonzro(4) = 1

      im1 = 0
      do 225 i = 1, m
        ti       = dble(im1)*ten
        dfj(i,4) = -ti*exp(-ti*x4)
        im1      = i
 225  continue

      return

 230  continue

      x5 = x(5)

      nonzro(5) = 1

      im1  = 0
      do 235 i = 1, m
        ti       = dble(im1)*ten
        dfj(i,5) = -ti*exp(-ti*x5)
        im1      = i
 235  continue

      return

 240  continue

      x2 = x(2)
      x4 = x(4)

      nonzro(2) = 1
      nonzro(4) = 1

      im1 = 0
      do 245 i = 1, m
        ti       = dble(im1)*ten
        t4       = ti*exp(-ti*x4)
        dfj(i,2) = -t4
        dfj(i,4) =  ti*x2*t4
        im1      = i
 245  continue

      return

 250  continue

      x3 = x(3)
      x5 = x(5)

      nonzro(3) = 1
      nonzro(5) = 1

      im1 = 0
      do 255 i = 1, m
        ti       = dble(im1)*ten
        t5       = ti*exp(-ti*x5)
        dfj(i,3) = -t5
        dfj(i,5) =  ti*x3*t5
        im1      = i
 255  continue

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

      double precision   tk, t4, t5

      intrinsic          dble, exp

      double precision   zero, ten
      parameter         (zero = 0.d0, ten = 10.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      tk = ten*dble(k-1)
      t4 = tk*exp(-tk*x(4))
      t5 = tk*exp(-tk*x(5))

      hess(2,4) = -t4
      hess(4,2) = hess(2,4)
      hess(3,5) = -t5
      hess(5,3) = hess(3,5)
      hess(4,4) = tk*x(2)*t4
      hess(5,5) = tk*x(3)*t5

      return
      end

