*************************************************************************
*               extended powell singular function
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

      double precision   ddot

      intrinsic          sqrt

      double precision   rtfive, rtten, rtten2
      common /PARAM1/    rtfive, rtten, rtten2
      save   /PARAM1/

      double precision   zero, one, two, three, four, five, ten
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0, four = 4.d0, five = 5.d0)
      parameter         (ten = 10.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

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

  10  continue

      nprobs = 1
      nstrts = 1

      nqrtr  = 3
      n      = 4*nqrtr
      m      = n

      rtten  = dsqrt( ten)
      rtten2 = two * rtten
      rtfive = dsqrt(five)

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      do 21 i = 1, nqrtr
        i4 = 4*i
        x(i4-3) =  three
        x(i4-2) = -one
        x(i4-1) =  zero
        x(i4  ) =  one
   21 continue

      return

*-----------------------------------------------------------------------

   30 continue

      do 31 i = 1, n
        x(i) = zero
   31 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, nqrtr
        i4   = 4*i
        i4m1 = i4 - 1
        i4m2 = i4 - 2
        i4m3 = i4 - 3
        xi4   = x(i4)
        xi4m1 = x(i4m1)
        xi4m2 = x(i4m2)
        xi4m3 = x(i4m3)
        t1 = xi4m2 - two*xi4m1
        t0 = xi4m3 - xi4
        f(i4m3) = xi4m3 + ten*xi4m2
        f(i4m2) = rtfive * (xi4m1 - xi4)
        f(i4m1) = t1 * t1
        f(i4  ) = rtten * t0 * t0
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 210  continue

      do 220 i = 1, nqrtr
        i4   = 4*i
        i4m1 = i4 - 1
        i4m2 = i4 - 2
        i4m3 = i4 - 3
        xi4   = x(i4)
        xi4m1 = x(i4m1)
        xi4m2 = x(i4m2)
        xi4m3 = x(i4m3)
        s1 = two*(xi4m2 - two*xi4m1)
        s0 = rtten2*(xi4m3 - xi4)
        fj(i4m3,i4m3) =  one
        fj(i4m3,i4m2) =  ten
        fj(i4m2,i4m1) =  rtfive
        fj(i4m2,  i4) = -rtfive
        fj(i4m1,i4m2) =  s1
        fj(i4m1,i4m1) = -two*s1
        fj(  i4,i4m3) =  s0
        fj(  i4,  i4) = -s0
 220  continue

      return

 300  continue

      do 310 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 310  continue

      do 320 i = 1, nqrtr
        i4   = 4*i
        i4m1 = i4 - 1
        i4m2 = i4 - 2
        i4m3 = i4 - 3
        xi4   = x(i4)
        xi4m1 = x(i4m1)
        xi4m2 = x(i4m2)
        xi4m3 = x(i4m3)
        t1 = xi4m2 - two*xi4m1
        t0 = xi4m3 - xi4
        s0 = rtten2 * t0
        f(i4m3) = xi4m3 + ten*xi4m2
        f(i4m2) = rtfive * (xi4m1 - xi4)
        f(i4m1) = t1 * t1
        f(  i4) = rtten * t0 * t0
        fj(i4m3,i4m3) =  one
        fj(i4m3,i4m2) =  ten
        fj(i4m2,i4m1) =  rtfive
        fj(i4m2,  i4) = -rtfive
        fj(i4m1,i4m2) =  two * t1
        fj(i4m1,i4m1) = -four * t1
        fj(  i4,i4m3) =  s0
        fj(  i4,  i4) = -s0
 320  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' extended powell singular function (more et al.)'//,
     *'        number of variables =', i4,'  (multiple of 4)'/,
     *'        number of functions =', i4,'  (     =  n    )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      double precision   rtfive, rtten, rtten2
      common /PARAM1/    rtfive, rtten, rtten2
      save   /PARAM1/

      double precision   zero, two, four, eight
      parameter         (zero = 0.d0, two = 2.d0)
      parameter         (four = 4.d0, eight = 8.d0)

*=======================================================================

      do 50 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
   50 continue

      goto ( 100, 200, 300) mod(k,4)

      if (mod(k,4) .eq. 1)  goto 100
      if (mod(k,4) .eq. 2)  goto 200
      if (mod(k,4) .eq. 3)  goto 300

      nonzro(  k) = 1
      nonzro(k-3) = 1

      dfj(k,k-3) = -rtten2
      dfj(k,  k) =  rtten2

      return

 100  continue

      nonzro(  k) = 1
      nonzro(k+3) = 1

      dfj(k+3,  k) =  rtten2
      dfj(k+3,k+3) = -rtten2

      return

 200  continue

      nonzro(  k) = 1
      nonzro(k+1) = 1

      dfj(k+1,  k) =  two
      dfj(k+1,k+1) = -four

      return

 300  continue

      nonzro(  k) = 1
      nonzro(k-1) = 1

      dfj(k,k-1) = -four
      dfj(k,  k) =  eight

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

      double precision   rtfive, rtten, rtten2
      common /PARAM1/    rtfive, rtten, rtten2
      save   /PARAM1/

      double precision   zero, two, four, eight
      parameter         (zero = 0.d0, two = 2.d0)
      parameter         (four = 4.d0, eight = 8.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      mk4 = mod(k,4)

      if (mk4 .eq. 1 .or. mk4 .eq. 2)  then
        linear = .true.
        return
      end if

      linear = .false.

      if (mk4 .eq. 0)  goto 400

      hess(k-1,k-1) =  two
      hess(  k,k-1) = -four
      hess(k-1,  k) = -four
      hess(  k,  k) =  eight

      return

 400  continue

      hess(k-3,k-3) =  rtten2
      hess(  k,k-3) = -rtten2
      hess(k-3,  k) = -rtten2
      hess(  k,  k) =  rtten2

      return
      end
