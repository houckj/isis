**************************************************************************
*               trigonometric function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
**************************************************************************

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

      double precision   dn, sum, cxi, sxj, xi, xj

      double precision   ddot

      intrinsic          dble, cos, sin

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

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

      dn  = dble(n)

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      call dcopy( n, (one/dn), 0, x, 1)

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      sum = zero
      do 110 i = 1, n
        xi   = x(i)
        cxi  = cos(xi)
        sum  = sum + cxi
        f(i) = dn + dble(i)*(one - cxi) - sin(xi)
 110  continue

      do 120 i = 1, n
        f(i) = f(i) - sum
 120  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 220 j = 1, n
        xj  = x(j)
        sxj = sin(xj)
        do 210 i = 1, n
          fj( i, j) = sxj
 210  continue
      fj(j, j) =  dble(j+1)*sxj - cos(xj)
 220  continue

      return

 300  continue

      sum = zero
      do 310 i = 1, n
        xi   = x(i)
        cxi  = cos(xi)
        sum  = sum + cxi
        f(i) = dn + dble(i)*(one - cxi) - sin(xi)
 310  continue

      do 320 i = 1, n
        f(i) = f(i) - sum
 320  continue

      do 340 j = 1, n
        xj  = x(j)
        sxj = sin(xj)
        do 330 i = 1, n
          fj( i, j) = sxj
 330  continue
      fj( j, j) = dble(j+1)*sxj - cos(xj)
 340  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' trigonometric function (more et al.) '//,
     *'        number of variables =', i4, '  (variable)'/,
     *'        number of functions =', i4, '  (   = n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, j
    
      double precision   cxk, xk

      intrinsic          dble, cos, sin

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      nonzro(k) = 1

      xk  = x(k)
      cxk = cos(xk)
      do 200 i = 1, n
        dfj( i, k) = cxk
  200 continue
      dfj(k,k) = dble(k+1)*cxk + sin(xk)

      return
      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            i, j

      double precision   cxi, xi

      intrinsic          dble, cos, sin

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      do 220 i = 1, n
        xi  = x(i)
        cxi = cos(xi)
        if (i .ne. k) hess( i, i) = cxi
        if (i .eq. k) hess( i, i) = dble(i+1)*cxi + sin(xi)
 220  continue

      return
      end
