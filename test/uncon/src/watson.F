****************************************************************************** 
*               watson function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
******************************************************************************

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

      integer            i, j, jm1, n1, nj, nj1

      double precision   r1, r2, s1, s2, t2, ti, t2i, tijm2, x1, xnj

      double precision   ddot

      intrinsic          dble

      double precision   zero, one, two, oneneg, point1, twnty9
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (oneneg = -1.d0, point1 = .1d0, twnty9 = 29.d0)

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
      nstart = 1

      n      = 20
      m      = 31

      n1   = n  - 1
      n2   = n1 - 1

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

  20  continue

      call dcopy( n, zero, 0, x, 1)

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 2.48631d-20

      return

*-----------------------------------------------------------------------

 100  continue

      x1 = x(1)
      r2 = x(n)
      r1 = n1*r2
      do 120 i = 1, 29
        ti  = dble(i) / twnty9
        s1  = r1
        s2  = r2
        nj  = n1
        nj1 = n2
        do 110 j = 2, n1
          xnj = x(nj)
          s1  = dble(nj1)*xnj + ti*s1
          s2  = xnj + ti*s2
          nj  = nj1
          nj1 = nj1 - 1
 110    continue
        s2   = x1 + ti*s2
        f(i) = s1 - s2*s2 - one
 120  continue

      f(30) = x1
      f(31) = x(2) - x1*x1 - one

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      r2 = x(n)
      do 220 i = 1, 29
        ti = dble(i) / twnty9
        s2 = r2
        nj = n1
        do 210 j = 1, n1
          s2 = x(nj) + ti*s2
          nj = nj - 1
 210    continue
        t2      =  two*s2
        fj(i,1) = -t2
        t2i     =  t2*ti
        fj(i,2) =  one - t2i
        jm1     =  2
        tijm2   =  ti
        do 220 j = 3, n
          fj(i,j) = (dble(jm1) - t2i) * tijm2
          jm1     = j
          tijm2   = tijm2 * ti
 220  continue

      fj(30,1) = one
      fj(31,1) = -two*x(1)
      fj(30,2) = zero
      fj(31,2) = one

      do 230 j = 3, n
        fj(30,j) = zero
        fj(31,j) = zero
 230  continue

      return

 300  continue

      x1 = x(1)
      r2 = x(n)
      r1 = n1*r2
      do 320 i = 1, 29
        ti = dble(i) / twnty9
        s1  = r1
        s2  = r2
        nj  = n1
        nj1 = n2
        do 310 j = 2, n1
          xnj = x(nj)
          s1  = dble(nj1)*xnj + ti*s1
          s2  = xnj + ti*s2
          nj  = nj1
          nj1 = nj1 - 1
 310    continue
        s2      =  x1 + ti*s2
        f(i)    =  s1 - s2*s2 - one
        t2      =  two*s2
        fj(i,1) = -t2
        t2i     =  t2*ti
        fj(i,2) =  one - t2i
        jm1     = 2
        tijm2   = ti
        do 320 j = 3, n
          fj(i,j) = (dble(jm1) - t2i) * tijm2
          jm1     = j
          tijm2   = tijm2 * ti
 320  continue

      f(30) = x1
      f(31) = x(2) - x1*x1 - one

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      fj(30,1) =  one
      fj(31,1) = -two*x1
      fj(30,2) =  zero
      fj(31,2) =  one

      do 330 j = 3, n
        fj(30,j) = zero
        fj(31,j) = zero
 330  continue

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' watson function (more et al.)'//,
     *'        number of variables =', i4, ' ( 2 <= n <= 31 )'/,
     *'        number of functions =', i4, ' (      31      )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, j

      double precision   ti

      double precision   zero, two, twnty9
      parameter         (zero = 0.d0, two = 2.d0, twnty9 = 29.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 1
        dfj( 30, j) = zero
        dfj( 31, j) = zero
  100 continue

      if (k .eq. 1)  dfj(31,1) = -two

      do 210 i = 1, 29
        ti = dble(i) / twnty9
        do 210 j = 1, n
          dfj(i,j) = -two*ti**(k+j-2)
 210  continue

      return
      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, lhess, hess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            i, j

      intrinsic          dble

      double precision   tk

      double precision   zero, two, twnty9
      parameter         (zero = 0.d0, two = 2.d0, twnty9 = 29.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      if (k .eq. 30)  linear = .true.
      if (k .eq. 30)  return

      linear = .false.

      if (k .eq. 31)  hess( 1, 1) = -two
      if (k .eq. 31)  return

      tk = dble(k)/twnty9
      do 210 i = 1, n
        do 210 j = 1, n
          hess( i, j) = -two*tk**(i+j-2)
 210  continue

      return
      end
