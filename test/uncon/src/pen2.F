*****************************************************************************
*               penalty function ii
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

      integer            j, jm1

      double precision   a, vj, x1, xj
      double precision   sj, sjm1, tj, tjm1, uj, ujm1

      double precision   ddot

      intrinsic          sqrt, dble, exp

      double precision   ap, aq, ar, d0, scale
      common /PARAM1/    ap, aq, ar, d0, scale
      save   /PARAM1/

      double precision   zero, one, two, point2
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point2 = .2d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1  = x(1)

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

      n      = 4
      m      = 2*n

      scale = 0.1d0
      a     = 1.d-5

      ar = sqrt(a)
      ap = ar*scale
      aq = ap*scale
      d0 = ar*exp(-scale)

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      do 21 j = 1, n
        x(j) = 0.5d0
   21 continue

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 9.37626d-6

      return

*-----------------------------------------------------------------------

 100  continue

      sum  = dble(n)*(x1*x1)
      f(1) = x1 - point2
      sj   = ar*exp(scale)
      tj   = ar*exp(x1*scale)
      jm1  = 1
      do 110 j = 2, n
        sjm1     = sj
        sj       = ar*exp(dble(j)*scale)
        xj       = x(j)
        tjm1     = tj
        tj       = ar*exp(xj*scale)
        f(j)     = (tj + tjm1) - (sj + sjm1)
        f(n+jm1) = tj - d0
        sum      = sum + dble(n-jm1)*(xj*xj)
        jm1      = j
 110  continue
      f(m) = sum - one

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 210  continue

      fj(1,1) = one
      fj(m,1) = two*dble(n)*x1
      uj      = ap*exp(x1*scale)
      jm1     = 1
      do 220 j = 2, n
        xj   = x(j)
        ujm1 = uj
        uj   = ap*exp(xj*scale)
        fj(j    , jm1) = ujm1
        fj(j    , j  ) = uj
        fj(n+jm1, j  ) = uj
        fj(m    , j  ) = two*dble(n-jm1)*xj
        jm1 = j
 220  continue

      return

 300  continue

      do 310 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 310  continue

      sum     = dble(n)*x1*x1
      f(1)    = x1 - point2
      fj(1,1) = one
      fj(m,1) = two*dble(n)*x1
      sj  = ar*exp(scale)
      tj  = ar*exp(x1*scale)
      uj  = scale*tj
      jm1 = 1
      do 320 j = 2, n
        sjm1 = sj
        sj   = ar*exp(dble(j)*scale)
        xj   = x(j)
        tjm1 = tj
        ujm1 = uj
        tj   = ar*exp(xj*scale)
        uj   = scale*tj
        f(j) = (tj + tjm1) - (sj + sjm1)
        f(n+jm1) = tj - d0
        vj   = dble(n-jm1)*xj
        sum  = sum + vj*xj
        fj(    j, jm1) = ujm1
        fj(    j,   j) = uj
        fj(n+jm1,   j) = uj
        fj(    m,   j) = two*vj
        jm1 = j
 320  continue
      f(m) = sum - one

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' penalty function ii (more et al.) '//,
     *'        number of variables =', i4, '  (variable)'/,
     *'        number of functions =', i4, '  (  = 2n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      intrinsic          dble, exp

      double precision   ap, aq, ar, d0, scale
      common /PARAM1/    ap, aq, ar, d0, scale
      save   /PARAM1/

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 1
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      t = aq*exp(x(k)*scale)

      if (k+1 .ge. 2 .and. k+1 .le. n)  dfj(k+1  ,k) = t
      if (k   .ge. 2 .and. k   .le. n)  dfj(k    ,k) = t
      if (k   .ge. 2 .and. k   .le. n)  dfj(k+n-1,k) = t

      dfj(m,k) = two*dble(n-k+1)

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, lhess, hess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      intrinsic          dble, exp

      double precision   ap, aq, ar, d0, scale
      common /PARAM1/    ap, aq, ar, d0, scale
      save   /PARAM1/

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      if (k .eq. 1)  linear = .true.
      if (k .eq. 1)  return
      linear = .false.

      if ( k .ge.   2  .and. k .le.  n)
     *  hess(k   ,k    ) = aq*exp(x(k)*scale)
      if ( k .ge.   2  .and. k .le.  n)
     *  hess(k-1 ,k-1  ) = aq*exp(x(k-1)*scale)
      if ( k .gt. n .and. k .lt. 2*n)
     * hess(k-n+1,k-n+1) = aq*exp(x(k-n+1)*scale)

      if (k .ne. 2*n)  return

      do 200 i = 1, n
        hess(i,i) = two*dble(n-i+1)
 200  continue

      return
      end
