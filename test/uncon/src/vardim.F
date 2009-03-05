****************************************************************************
*               variably dimensioned function
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

      integer            i, j

      double precision   dn, dj, sum, xim1, xjm1

      double precision   ddot

      intrinsic          dble

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

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

   10 continue

      nprobs = 1
      nstrts = 1

      n      = 10
      m      = n + 2

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      dn = dble(n)
      do 21 j = 1, n
        x(j) = one - (dble(j)/dn)
   21 continue

      return

*-----------------------------------------------------------------------

   30 continue

      do 31 i = 1, n
        x(i) = one
   31 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      sum = zero
      do 110 i = 1, n
        xim1 = x(i) - one
        sum  = sum + dble(i)*xim1
        f(i) = xim1
 110  continue
      f(n+1) = sum
      f(n+2) = sum*sum

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      sum = zero
      do 210 j = 1, n
        sum = sum + dble(j)*(x(j) - one)
        call dcopy( m, zero, 0, fj( 1, j), 1)
 210  continue

      do 220 j = 1, n
        dj        = dble(j)
        fj(j  ,j) = one
        fj(n+1,j) = dj
        fj(n+2,j) = two*sum*dj
 220  continue

      return

 300  continue

      sum = zero
      do 310 j = 1, n
        xjm1 = x(j) - one
        sum  = sum + dble(j)*xjm1
        call dcopy( m, zero, 0, fj( 1, j), 1)
        f(j) = xjm1
 310  continue
      f(n+1) = sum
      f(n+2) = sum*sum

      do 320 j = 1, n
        dj        = dble(j)
        fj(j  ,j) = one
        fj(n+1,j) = dj
        fj(n+2,j) = two*sum*dj
 320  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' variably dimensioned function (more et al.) '//,
     *'        number of variables =', i4, '  (variable)'/,
     *'        number of functions =', i4, '  (  = n+2 )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            n1, j

      double precision   tk

      intrinsic          dble

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      n1 = n+1
      do 100 j = 1, n
        nonzro(j) = 1
        call dcopy( n1, zero, 0, dfj(1, j), 1)
  100 continue

      tk = two*dble(k)
      do 200 j = 1, n
        dfj(n+2,j) = tk*dble(j)
 200  continue

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            i, im1, j

      double precision   di, ti

      intrinsic          dble

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .true.
      if (k .ne. n+2)  return
      linear = .false.

      hess(1,1) = two
      do 200 i = 2, n
        im1 = i - 1
        di  = dble(i)
        ti  = two*dfi
        hess(i,i) = ti*di
        do 200 j = 1, im1
          hess(i,j) = ti*dble(j)
          hess(j,i) = hess(i,j)
 200  continue

      return
      end

