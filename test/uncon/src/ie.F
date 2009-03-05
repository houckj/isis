*****************************************************************************
*               discrete integral function
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

      double precision   h2, bj, cj, oj, sj, ti, tj, xj 

      double precision   ddot

      intrinsic          dble

      double precision   h
      common /PARAM1/    h
      save   /PARAM1/

      double precision   zero, half, one, two, three
      parameter         (zero = 0.d0, half = .5d0, one = 1.d0)
      parameter         (two = 2.d0, three = 3.d0)

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

   10 continue

      nprobs = 1
      nstrts = 1

      n      = 10
      m      = n

      h      = one / dble(n+1)
      h2     = half*h

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      do 21 j = 1, n
        tj   = dble(j)*h
        x(j) = tj*(tj-one)
  21  continue

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      call dcopy( n, x, 1, f, 1)
      do 110 j = 1, n
        tj = dble(j)*h
        oj = one - tj
        cj = h2*(x(j) + tj + one)**3
        do 110 i = 1, n
          ti = dble(i)*h
          oi = one - ti
          if (j .le. i)  f(i) = f(i) + oi*tj*cj
          if (j .gt. i)  f(i) = f(i) + ti*oj*cj
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        tj = dble(j)*h
        oj = one - tj
        sj = three*h2*(x(j) + tj + one)**2
        fj(j,j) =  one + oj*tj*sj
        do 210 i = 1, n
          ti = dble(i)*h
          oi = one - ti
          if (j .lt. i)  fj(i,j) =  oi*tj*sj
          if (j .gt. i)  fj(i,j) =  ti*oj*sj
 210  continue

      return

 300  continue

      call dcopy( n, x, 1, f, 1)
      do 310 j = 1, n
        xj = x(j)
        tj = dble(j)*h
        oj = one - tj
        bj = (xj + tj + one)
        sj = h2*bj*bj
        cj = sj*bj
        sj = sj*three
        fj(j,j) = one + oj*tj*sj
        do 310 i = 1, n
          ti = dble(i)*h
          oi = one - ti
          if (j .le. i)  f(i) = f(i) + oi*tj*cj
          if (j .gt. i)  f(i) = f(i) + ti*oj*cj
          if (j .lt. i)  fj(i,j) = oi*tj*sj
          if (j .gt. i)  fj(i,j) = ti*oj*sj
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' discrete integral function (more et al.)'//,
     *'        number of variables =', i4,'  (variable)'/,
     *'        number of functions =', i4,'  (   = n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, j

      double precision   ti, tk, xk

      intrinsic          dble

      double precision   h
      common /PARAM1/    h
      save   /PARAM1/

      double precision   zero, one, three
      parameter         (zero = 0.d0, one = 1.d0)
      parameter         (three = 3.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
 100  continue

      nonzro(k) = 1

      tk = dble(k)*h
      xk = x(k)
      do 200 i = 1, n
        ti = dble(i)*h
        if (k .le. i)  dfj(i,k) =
     *    three*h*(one-ti)*tk*(xk+tk+one)**2
        if (k .gt. i)  dfj(i,k) =
     *    three*h*ti*(one-tk)*(xk+tk+one)**2
  200 continue

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

      double precision   tk, tj

      intrinsic          dble

      double precision   h
      common /PARAM1/    h
      save   /PARAM1/

      double precision   zero, one, three
      parameter         (zero = 0.d0, one = 1.d0)
      parameter         (three = 3.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy(  n, zero, 0, hess( 1, j), 1)
 100  continue

      linear = .false.

      tk = dble(k)*h
      do 200 j = 1, n
        tj = dble(j)*h
        if (j .le. k)  hess(j,j) = three*h*(one-tk)*tj*(x(j)+tj+one)**2
        if (j .gt. k)  hess(j,j) = three*h*tk*(one-tj)*(x(j)+tj+one)**2
 200  continue

      return
      end
