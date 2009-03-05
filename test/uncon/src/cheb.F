*************************************************************************
*               chebyquad function
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

      integer            m2
      double precision   coeff
      common /PARAM1/    coeff, m2
      save   /PARAM1/

      logical            lf, lj

      integer            na, nb, nc, nd, nt, nh

      integer            i, i2, j

      double precision   dn1
      double precision   s, s0, s1, t, t0, t1, tj, u

      double precision   ddot

      intrinsic          dble

      double precision   zero, one, two, four
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (four = 4.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

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

      n      = 10
      m      = n

      m2     = m/2

      coeff  = one / dble(n)

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

  20  continue

      dn1 = one / dble(n+1)
      do 21 i = 1, n
        x(i) = dble(i) * dn1
  21  continue

      return

*-----------------------------------------------------------------------

  30  continue

      ftf = 4.77271d-3

      return

*-----------------------------------------------------------------------

 100  continue

      call dcopy( m, zero, 0, f, 1)
      do 110 i = 1, m2
        i2 = 2*i
        f(i2) = one / dble((i2*i2) - 1)
 110  continue

      do 125 j = 1, n
        t0 = one
        t1 = two * x(j) - one
        tj = t1
        f(1) = f(1) + coeff*t1
        do 120 i = 2, m
          t = two*tj*t1 - t0
          f(i) = f(i) + coeff*t
          t0 = t1
          t1 = t
 120    continue
 125  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 225 j = 1, n
        t0 = one
        s0 = zero
        s1 = two
        t1 = two * x(j) - one
        tj = t1
        fj( 1, j) = coeff*s1
        do 220 i = 2, m
          u  = two * tj
          t  = u*t1 - t0
          s  = four*t1 + u*s1 - s0
          fj( i, j) = coeff*s
          t0 = t1
          t1 = t
          s0 = s1
          s1 = s
 220    continue
 225  continue

      return

 300  continue

      call dcopy( m, zero, 0, f, 1)

      do 310 i = 1, m2
        i2 = 2*i
        f(i2) = one / dble((i2*i2) - 1)
 310  continue

      do 325 j = 1, n
        t0 = one
        s0 = zero
        s1 = two
        t1 = two * x(j) - one
        tj = t1
        f(1) = f(1) + coeff*t1
        fj( 1, j) = coeff*s1
        do 320 i = 2, m
          u  = two*tj
          t  = u*t1 - t0
          s  = four*t1 + u*s1 - s0
          f(i) = f(i) + coeff*t
          fj( i, j) = coeff*s
          t0 = t1
          t1 = t
          s0 = s1
          s1 = s
 320    continue
 325  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' chebyquad function (more et al.)'//,
     *'        number of variables =', i4,'  (variable)'/,
     *'        number of functions =', i4,'  (  >= n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, j

      double precision   r, r0, r1, s, s0, s1, t, t0, t1, tk, u

      integer            m2
      double precision   coeff
      common /PARAM1/    coeff, m2
      save   /PARAM1/

      double precision   zero, one, two, four, eight
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (four = 4.d0, eight = 8.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      nonzro(k) = 1

      t0 = one
      s0 = zero
      s1 = two
      r0 = zero
      r1 = zero
      t1 = two * x(k) - one
      tk = t1
      dfj( 1, k) = zero
      do 200 i = 2, m
        u = two*tk
        t = u*t1 - t0
        s = four*t1 + u*s1 - s0
        r = eight*s1 + u*r1 - r0
        dfj( i, k) = coeff*r
        t0 = t1
        t1 = t
        s0 = s1
        s1 = s
        r0 = r1
        r1 = r
 200  continue

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

      double precision   r, r0, r1, s, s0, s1, t, t0, t1, tj, u

      integer            m2
      double precision   coeff
      common /PARAM1/    coeff, m2
      save   /PARAM1/

      double precision   zero, one, two, four, eight
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (four = 4.d0, eight = 8.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      if (k .eq. 1)  linear = .true.
      if (k .eq. 1)  return
      linear = .false.

      do 220 j = 1, n
        t0 = one
        s0 = zero
        s1 = two
        r0 = zero
        r1 = zero
        t1 = two * x(j) - one
        tj = t1
        do 210 i = 2, k
          u = two*tj
          t = u*t1 - t0
          s = four*t1 + u*s1 - s0
          r = eight*s1 + u*r1 - r0
          t0 = t1
          t1 = t
          s0 = s1
          s1 = s
          r0 = r1
          r1 = r
  210   continue
        hess( j, j) = coeff*r
 220  continue

      return
      end

