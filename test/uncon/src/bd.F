******************************************************************************
*               brown and dennis function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
******************************************************************************

      subroutine getfun( x, n, f, m, ftf, fj, lfj, g, mode)

      implicit double precision (a-h,o-z)

      integer            n, m, lfj, mode

      double precision   x(n), f(m), ftf, fj(lfj,n), g(n)

      common /PROBLM/    nprob, nprobs, nstart, nstrts

      common /IOUNIT/    nout

      logical            lf, lj

      integer            na, nb, nc, nd, nt, nh

      integer            i, j

      double precision   x1, x2, x3, x4
      double precision   f1, f3, ei, ci, si, ti
      
      double precision   ddot

      intrinsic          dble, exp, cos, sin

      double precision   zero, one, two, point2
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point2 = .2d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)

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

      n      =  4
      m      = 20

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) = 25.d0
      x(2) =  5.d0
      x(3) = -5.d0
      x(4) = -1.d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 8.58222d+4

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ti   = dble(i)*point2
        ei   = exp(ti)
        si   = sin(ti)
        ci   = cos(ti)
        f(i) = (x1 + ti*x2 - ei)**2 + (x3 + x4*si - ci)**2
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        ti = dble(i)*point2
        ei = exp(ti)
        si = sin(ti)
        ci = cos(ti)
        f1 = two*(x1 + ti*x2 - ei)
        f3 = two*(x3 + x4*si - ci)
        fj( i, 1) = f1
        fj( i, 2) = f1 * ti
        fj( i, 3) = f3
        fj( i, 4) = f3 * si
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ti = dble(i)*point2
        ei = exp(ti)
        si = sin(ti)
        ci = cos(ti)
        f1 = two*(x1 + ti*x2 - ei)
        f3 = two*(x3 + x4*si - ci)
        f(i) = (x1 + ti*x2 - ei)**2 + (x3 + x4*si - ci)**2
        fj( i, 1) = f1
        fj( i, 2) = f1 * ti
        fj( i, 3) = f3
        fj( i, 4) = f3 * si
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 320 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 320  continue

      return

9999  format(/'1',70('=')//,
     *' brown and dennis function (more et al.) '//,
     *'        number of variables =', i4, '  (    4   )'/,
     *'        number of functions =', i4, '  (  >= 4  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)
     
      integer            i

      double precision   ti, t2, si, s2

      intrinsic          dble, sin

      double precision   zero, two, point2
      parameter         (zero = 0.d0, two = 2.d0, point2 = .2d0)

*=======================================================================

      goto ( 210, 220, 230, 240), k

 210  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 0
      nonzro(4) = 0

      do 215 i = 1, m
        ti = point2*dble(i)
        dfj( i, 1) = two
        dfj( i, 2) = two*ti
        dfj( i, 3) = zero
        dfj( i, 4) = zero
 215  continue

      return

 220  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 0
      nonzro(4) = 0

      do 225 i = 1, m
        ti = point2*dble(i)
        t2 = two*ti
        dfj( i, 1) = t2
        dfj( i, 2) = t2*ti
        dfj( i, 3) = zero
        dfj( i, 4) = zero
 225  continue

      return

 230  continue

      nonzro(1) = 0
      nonzro(2) = 0
      nonzro(3) = 1
      nonzro(4) = 1

      do 235 i = 1, m
        ti = point2*dble(i)
        si = sin(ti)
        dfj( i, 1) = zero
        dfj( i, 2) = zero
        dfj( i, 3) = two
        dfj( i, 4) = two*si
 235  continue

      return

 240  continue

      nonzro(1) = 0
      nonzro(2) = 0
      nonzro(3) = 1
      nonzro(4) = 1

      do 245 i = 1, m
        ti = point2*dble(i)
        si = sin(ti)
        s2 = two*si
        dfj( i, 1) = zero
        dfj( i, 2) = zero
        dfj( i, 3) = s2
        dfj( i, 4) = s2*si
 245  continue

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   tk, t2, sk, s2

      intrinsic          dble, sin

      double precision   zero, two, point2
      parameter         (zero = 0.d0, two = 2.d0, point2 = .2d0)

*=======================================================================

      linear = .false.

      tk = point2*dble(k)
      sk = sin(tk)
      t2 = two*tk
      s2 = two*sk

      hess(1,1) = two
      hess(1,2) = t2
      hess(1,3) = zero
      hess(1,4) = zero
      hess(2,1) = hess(1,2)
      hess(2,2) = t2*tk
      hess(2,3) = zero
      hess(2,4) = zero
      hess(3,1) = zero
      hess(3,2) = zero
      hess(3,3) = two
      hess(3,4) = s2
      hess(4,1) = zero
      hess(4,2) = zero
      hess(4,3) = hess(3,4)
      hess(4,4) = s2*sk

      return
      end

