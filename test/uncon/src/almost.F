***************************************************************************
*               brown almost linear function
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

      integer            i, j

      double precision   xi, sum, prd, t

      double precision   ddot

      intrinsic          dble

      integer            nvm1
      common /PARAM1/    nvm1
      save   /PARAM1/

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

  10  continue

      nprobs = 1
      nstrts = 1

      n      = 10
      nvm1   = n - 1
      m      = n

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      call dcopy( n, (.5d0), 0, x, 1)

      return

*-----------------------------------------------------------------------

   30 continue

      call dcopy( n, one, 0, x, 1)

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      sum = zero
      prd =  one
      do 110 i = 1, n
        xi  = x(i)
        sum = sum + xi
        prd = prd * xi
        f(i) = xi - dble(n+1)
 110  continue

      do 120 i = 1, nvm1
        f(i) = f(i) + sum
 120  continue
      f(n) = prd - one

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 220 j = 1, n
        do 210 i = 1, nvm1
          fj( i, j) = one
 210  continue
      fj(    j, j) = two
 220  continue

      do 240 j = 1, n
        t = one
        do 230 i = 1, n
          if (i .ne. j)  t = t * x(i)
 230  continue
      fj( n, j) = t
 240  continue

      return

 300  continue

      sum = zero
      prd =  one
      do 310 i = 1, n
        xi  = x(i)
        sum = sum + xi
        prd = prd * xi
        f(i) = xi - dble(n+1)
 310  continue

      do 320 i = 1, nvm1
        f(i) = f(i) + sum
 320  continue
      f(n) = prd - one

      do 340 j = 1, n
        do 330 i = 1, nvm1
          fj( i, j) = one
 330  continue
      fj(    j, j) = two
 340  continue

      do 360 j = 1, n
        t = one
        do 350 i = 1, n
          if (i .ne. j)  t = t * x(i)
 350  continue
      fj( n, j) = t
 360  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 370 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 370  continue

      return

9999  format(/'1',70('=')//,
     *' brown almost-linear function (more et al.) '//,
     *'        number of variables =', i4, '  (variable)'/,
     *'        number of functions =', i4, '  (  >= n  )'//,
     *        ' ', 70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

c  derivative of the jacobian with respect to the kth variable

      implicit double precision (a-h,o-z)

      integer            k, n, m, ldfj, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, j

      double precision   t

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 1
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      do 220 j = 1, n
        t = one
        do 210 i = 1, n
          if (i .ne. j .and. i .ne. k)  t = t * x(i)
 210    continue
        dfj( n, j) = t
 220  continue

      nonzro(k) = 0
      dfj( n, k) = zero

      return
      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear
        
      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            i, im1, j

      double precision   t

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      if (k .ne. n)  linear = .true.
      if (k .ne. n)  return
      linear = .false.

      do 220 i = 2, n
        im1 = i - 1
        do 220 j = 1, im1
          t = one
          do 210 m = 1, n
            if (m .ne. i .and. m .ne. j)  t = t * x(m)
 210      continue
        hess(i,j) = t
        hess(j,i) = t
 220  continue

      return
      end

