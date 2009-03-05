****************************************************************************
*               osborne 2 function
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
    
      integer            i, im1, j

      double precision   x01, x02, x03, x04, x05, x06
      double precision   x07, x08, x09, x10, x11

      double precision   ti, t09, t10, t11, s09, s10, s11
      double precision   e1, e2, e3, e4, r2, r3, r4

      double precision   ddot

      intrinsic          dble, exp

      double precision   y
      common /PARAM1/    y(65)
      save   /PARAM1/    

      double precision   zero, one, two, point1
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = 0.1d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x01 = x( 1)
      x02 = x( 2)
      x03 = x( 3)
      x04 = x( 4)
      x05 = x( 5)
      x06 = x( 6)
      x07 = x( 7)
      x08 = x( 8)
      x09 = x( 9)
      x10 = x(10)
      x11 = x(11)

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

      n      = 11
      m      = 65

      y( 1) = 1.366d0
      y( 2) = 1.191d0
      y( 3) = 1.112d0
      y( 4) = 1.013d0
      y( 5) = 0.991d0
      y( 6) = 0.885d0
      y( 7) = 0.831d0
      y( 8) = 0.847d0
      y( 9) = 0.786d0
      y(10) = 0.725d0
      y(11) = 0.746d0
      y(12) = 0.679d0
      y(13) = 0.608d0
      y(14) = 0.655d0
      y(15) = 0.615d0
      y(16) = 0.606d0
      y(17) = 0.602d0
      y(18) = 0.626d0
      y(19) = 0.651d0
      y(20) = 0.724d0
      y(21) = 0.649d0
      y(22) = 0.649d0
      y(23) = 0.694d0
      y(24) = 0.644d0
      y(25) = 0.624d0
      y(26) = 0.661d0
      y(27) = 0.612d0
      y(28) = 0.558d0
      y(29) = 0.533d0
      y(30) = 0.495d0
      y(31) = 0.500d0
      y(32) = 0.423d0
      y(33) = 0.395d0
      y(34) = 0.375d0
      y(35) = 0.372d0
      y(36) = 0.391d0
      y(37) = 0.396d0
      y(38) = 0.405d0
      y(39) = 0.428d0
      y(40) = 0.429d0
      y(41) = 0.523d0
      y(42) = 0.562d0
      y(43) = 0.607d0
      y(44) = 0.653d0
      y(45) = 0.672d0
      y(46) = 0.708d0
      y(47) = 0.633d0
      y(48) = 0.668d0
      y(49) = 0.645d0
      y(50) = 0.632d0
      y(51) = 0.591d0
      y(52) = 0.559d0
      y(53) = 0.597d0
      y(54) = 0.625d0
      y(55) = 0.739d0
      y(56) = 0.710d0
      y(57) = 0.729d0
      y(58) = 0.720d0
      y(59) = 0.636d0
      y(60) = 0.581d0
      y(61) = 0.428d0
      y(62) = 0.292d0
      y(63) = 0.162d0
      y(64) = 0.098d0
      y(65) = 0.054d0

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x( 1) =  1.3d0
      x( 2) =  0.65d0
      x( 3) =  0.65d0
      x( 4) =  0.7d0
      x( 5) =  0.6d0
      x( 6) =  3.d0
      x( 7) =  5.d0
      x( 8) =  7.d0
      x( 9) =  2.d0
      x(10) =  4.5d0
      x(11) =  5.5d0

      return

*-----------------------------------------------------------------------

   30 continue

      x( 1) =  1.3100d0
      x( 2) =  0.4315d0
      x( 3) =  0.6336d0
      x( 4) =  0.5993d0
      x( 5) =  0.7539d0
      x( 6) =  0.9056d0
      x( 7) =  1.3651d0
      x( 8) =  4.8248d0
      x( 9) =  2.3988d0
      x(10) =  4.5689d0
      x(11) =  5.6754d0

      ftf = 4.01683d-2

      return

*-----------------------------------------------------------------------

 100  continue

      im1 = 0
      do 110 i = 1, m
        ti  = dble(im1)*point1
        t09 = ti - x09
        t10 = ti - x10
        t11 = ti - x11
        s09 = t09 * t09
        s10 = t10 * t10
        s11 = t11 * t11
        e1  = exp(-ti*x05)
        e2  = exp(-s09*x06)
        e3  = exp(-s10*x07)
        e4  = exp(-s11*x08)
        f(i) = (x01*e1 + x02*e2 + x03*e3 + x04*e4) - y(i)
        im1 = i
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      im1 = 0
      do 210 i = 1, m
        ti  = dble(im1)*point1
        t09 = ti - x09
        t10 = ti - x10
        t11 = ti - x11
        s09 = t09 * t09
        s10 = t10 * t10
        s11 = t11 * t11
        e1  = exp(-ti*x05)
        e2  = exp(-s09*x06)
        e3  = exp(-s10*x07)
        e4  = exp(-s11*x08)
        r2  = x02*e2
        r3  = x03*e3
        r4  = x04*e4
        fj( i, 1) = e1
        fj( i, 2) = e2
        fj( i, 3) = e3
        fj( i, 4) = e4
        fj( i, 5) = -ti*x01*e1
        fj( i, 6) = -s09*r2
        fj( i, 7) = -s10*r3
        fj( i, 8) = -s11*r4
        fj( i, 9) =  two*t09*x06*r2
        fj( i,10) =  two*t10*x07*r3
        fj( i,11) =  two*t11*x08*r4
        im1 = i
 210  continue

      return

 300  continue

      im1 = 0
      do 310 i = 1, m
        ti  = dble(im1)*point1
        t09 = ti - x09
        t10 = ti - x10
        t11 = ti - x11
        s09 = t09 * t09
        s10 = t10 * t10
        s11 = t11 * t11
        e1  = exp(-ti*x05)
        e2  = exp(-s09*x06)
        e3  = exp(-s10*x07)
        e4  = exp(-s11*x08)
        r1  = x01*e1
        r2  = x02*e2
        r3  = x03*e3
        r4  = x04*e4
        f(i) = (r1 + r2 + r3 + r4) - y(i)
        fj( i, 1) = e1
        fj( i, 2) = e2
        fj( i, 3) = e3
        fj( i, 4) = e4
        fj( i, 5) = -ti*r1
        fj( i, 6) = -s09*r2
        fj( i, 7) = -s10*r3
        fj( i, 8) = -s11*r4
        fj( i, 9) =  two*t09*x06*r2
        fj( i,10) =  two*t10*x07*r3
        fj( i,11) =  two*t11*x08*r4
        im1 = i
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' osborne 2 function (more et al.) : gaussian + exponential'//,
     *'        number of variables =', i4, '  (11)'/,
     *'        number of functions =', i4, '  (65)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, im1, j

      double precision   x01, x02, x03, x04, x05, x06
      double precision   x07, x08, x09, x10, x11

      double precision   e1, e2, e3, e4, q, r, s09, s10, s11
      double precision   ti, t2, t3, t4, t09, t10, t11

      intrinsic          dble, exp

      double precision   zero, one, two, point1
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = 0.1d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      goto(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100,2110), k

2010  continue

      x05 = x( 5)

      nonzro( 5) = 1

      im1 = 0
      do 2015 i = 1, m
        ti       =  dble(im1)*point1
        e1       =  exp(-ti*x05)
        dfj(i,5) = -ti*e1
        im1      =  i
2015  continue

      return

2020  continue

      x06 = x( 6)
      x09 = x( 9)

      nonzro( 6) = 1
      nonzro( 9) = 1

      im1 = 0
      do 2025 i = 1, m
        ti       =  dble(im1)*point1
        t09      =  ti - x09
        s09      =  t09 * t09
        e2       =  exp(-s09*x06)
        dfj(i,6) = -s09 * e2
        dfj(i,9) =  two * t09 * x06 * e2
        im1      =  i
2025  continue

      return

2030  continue

      x07 = x( 7)
      x10 = x(10)

      nonzro( 7) = 1
      nonzro(10) = 1

      im1 = 0
      do 2035 i = 1, m
        ti        =  dble(im1)*point1
        t10       =  ti - x10
        s10       =  t10 * t10
        e3        =  exp(-s10*x07)
        dfj(i, 7) = -s10 * e3
        dfj(i,10) =  two * t10 * x07 * e3
        im1       =  i
2035  continue

      return

2040  continue

      x08 = x( 8)
      x11 = x(11)

      nonzro( 8) = 1
      nonzro(11) = 1

      im1 = 0
      do 2045 i = 1, m
        ti        =  dble(im1)*point1
        t11       =  ti - x11
        s11       =  t11 * t11
        e4        =  exp(-s11*x08)
        dfj(i, 8) = -s11 * e4
        dfj(i,11) =  two * t11 * x08 * e4
        im1 = i
2045  continue

      return

2050  continue

      x01 = x( 1)
      x05 = x( 5)

      nonzro( 1) = 1
      nonzro( 5) = 1

      im1 = 0
      do 2055 i = 1, m
        ti       =  dble(im1)*point1
        e1       =  exp(-ti*x05)
        dfj(i,1) = -ti * e1
        dfj(i,5) = (ti * ti) * x01 * e1
        im1      =  i
2055  continue

      return

2060  continue

      x02 = x( 2)
      x06 = x( 6)
      x09 = x( 9)

      nonzro( 2) = 1
      nonzro( 6) = 1
      nonzro( 9) = 1

      im1 = 0
      do 2065 i = 1, m
        ti       =  dble(im1)*point1
        t09      =  ti - x09
        s09      =  t09 * t09
        r        =  s09*x06
        e2       =  exp(-r)
        q        =  x02*e2
        dfj(i,2) = -s09 * e2
        dfj(i,6) =  s09 * s09 * q
        dfj(i,9) =  two * t09 * q * (one - r)
        im1      =  i
2065  continue

      return

2070  continue

      x03 = x( 3)
      x07 = x( 7)
      x10 = x(10)

      nonzro( 3) = 1
      nonzro( 7) = 1
      nonzro(10) = 1

      im1 = 0
      do 2075 i = 1, m
        ti        = dble(im1)*point1
        t10       = ti - x10
        s10       = t10 * t10
        r         = s10*x07
        e3        = exp(-r)
        q         = x03*e3
        dfj(i, 3) = -s10 * e3
        dfj(i, 7) =  s10 * s10 * q
        dfj(i,10) =  two * t10 * q * (one - r)
        im1       = i
2075  continue

      return

2080  continue

      x04 = x( 4)
      x08 = x( 8)
      x11 = x(11)

      nonzro( 4) = 1
      nonzro( 8) = 1
      nonzro(11) = 1

      im1 = 0
      do 2085 i = 1, m
        ti        =  dble(im1)*point1
        t11       =  ti - x11
        s11       =  t11 * t11
        r         =  s11*x08
        e4        =  exp(-r)
        q         =  x04*e4
        dfj(i, 4) = -s11 * e4
        dfj(i, 8) =  s11 * s11 * q
        dfj(i,11) =  two * t11 * q * (one - r)
        im1       =  i
2085  continue

      return

2090  continue

      x02 = x( 2)
      x06 = x( 6)
      x09 = x( 9)

      nonzro( 2) = 1
      nonzro( 6) = 1
      nonzro( 9) = 1

      im1 = 0
      do 2095 i = 1, m
        ti       = dble(im1)*point1
        t09      = ti - x09
        r        = x06 * (t09 * t09) 
        t2       = two * exp(-r)
        q        = x02 * t2
        dfj(i,2) = t09 * x06 * t2
        dfj(i,6) = t09 * q * (one - r)
        dfj(i,9) = x06 * q * (two*r - one)
        im1      = i
2095  continue

      return

2100  continue

      x03 = x( 3)
      x07 = x( 7)
      x10 = x(10)

      nonzro( 3) = 1
      nonzro( 7) = 1
      nonzro(10) = 1

      im1 = 0
      do 2105 i = 1, m
        ti        = dble(im1)*point1
        t10       = ti - x10
        r         = x07 * (t10 * t10)
        t3        = two * exp(-r)
        q         = x03 * t3
        dfj(i, 3) = t10 * x07 * t3
        dfj(i, 7) = t10 * q * (one - r)
        dfj(i,10) = x07 * q * (two*r - one)
        im1 = i
2105  continue

      return

2110  continue

      x04 = x( 4)
      x08 = x( 8)
      x11 = x(11)

      nonzro( 4) = 1
      nonzro( 8) = 1
      nonzro(11) = 1

      im1 = 0
      do 2115 i = 1, m
        ti        = dble(im1)*point1
        t11       = ti - x11
        r         = x08 * (t11 * t11)
        t4        = two  * exp(-r)
        q         = x04 * t4
        dfj(i,4)  = t11 * x08 * t4
        dfj(i,8)  = t11 * q * (one - r)
        dfj(i,11) = x08 * q * (two*r - one)
        im1       = i
2115  continue

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

      double precision   d2, d3, d4, e1, e2, e3, e4
      double precision   q2, q3, q4, r2, r3, r4, s09, s10, s11
      double precision   tk, t09, t10, t11, v, u2, u3, u4

      intrinsic          dble, exp

      double precision   zero, one, two, point1
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = 0.1d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      x01 = x( 1)
      x02 = x( 2)
      x03 = x( 3)
      x04 = x( 4)
      x05 = x( 5)
      x06 = x( 6)
      x07 = x( 7)
      x08 = x( 8)
      x09 = x( 9)
      x10 = x(10)
      x11 = x(11)

      tk  = point1*dble(k-1)
      t09 = tk - x09
      t10 = tk - x10
      t11 = tk - x11
      s09 = t09 * t09
      s10 = t10 * t10
      s11 = t11 * t11
      r2  = s09 * x06
      r3  = s10 * x07
      r4  = s11 * x08
      e1  = exp(-tk*x05)
      e2  = exp(-r2)
      e3  = exp(-r3)
      e4  = exp(-r4)
      v   = tk * e1
      d2  = two * e2
      d3  = two * e3
      d4  = two * e4
      q2  = x02 * d2
      q3  = x03 * d3
      q4  = x04 * d4
      u2  = s09 * e2
      u3  = s10 * e3
      u4  = s11 * e4

      hess( 1, 5) = - v
      hess( 5, 1) =  hess( 1, 5)

      hess( 2, 6) =  - u2
      hess( 6, 2) =  hess( 2, 6)
      hess( 2, 9) =  t09 * x06 * d2
      hess( 9, 2) =  hess( 2, 9)

      hess( 3, 7) =  - u3
      hess( 7, 3) =  hess( 3, 7)
      hess( 3,10) =  t10 * x07 * d3
      hess(10, 3) =  hess( 3,10)

      hess( 4, 8) =  - u4
      hess( 8, 4) =  hess( 4, 8)
      hess( 4,11) =  t11 * x08 * d4
      hess(11, 4) =  hess( 4,11)

      hess( 5, 5) = tk * x01 * v

      hess( 6, 6) =  x02 * s09 * u2
      hess( 6, 9) =  t09 * q2 * (one - r2)
      hess( 9, 6) =  hess( 6, 9)

      hess( 7, 7) =  x03 * s10 * u3
      hess( 7,10) =  t10 * q3 * (one - r3)
      hess(10, 7) =  hess(10, 7)

      hess( 8, 8) =  x04 * s11 * u4
      hess( 8,11) =  t11 * q4 * (one - r4)
      hess(11, 8) =  hess( 8,11)

      hess( 9, 9) =  x06 * q2 * (two*r2 - one)

      hess(10,10) =  x07 * q3 * (two*r3 - one)

      hess(11,11) =  x08 * q4 * (two*r4 - one)

      return
      end
