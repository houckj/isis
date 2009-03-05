************************************************************************
*               beale function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
************************************************************************

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

      double precision   x1, x2
      double precision   s, t

      double precision   ddot

      double precision   y
      common /PARAM1/    y(3)
      save   /PARAM1/

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .lt. -1)  goto    30

      x1 = x(1)
      x2 = x(2)

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

      n      = 2
      m      = 3

      y(1)   = 1.5d0
      y(2)   = 2.25d0
      y(3)   = 2.625d0

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  1.d0
      x(2) =  1.d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  3.d0
      x(2) =   .5d0

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        f(i) = y(i) - x1*(one - x2**i)
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      im1 = 0
      do 210 i = 1, m
        t = x2**im1
        fj(i,1) = t*x2 - one
        fj(i,2) = x1*dble(i)*t
        im1 = i
 210  continue

      return

 300  continue

      im1 = 0
      do 310 i = 1, m
        t = x2**im1
        s = t*x2
        f(i) = y(i) - x1*(one - s)
        fj(i,1) = s - one
        fj(i,2) = x1*dble(i)*t
        im1 = i
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 320 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 320  continue

      return

9999  format(/'1',70('=')//,
     *' beale function (more et al.)'//,
     *'        number of variables =', i4, '  (2)'/,
     *'        number of functions =', i4, '  (3)'//,
     *        ' ',70('=')/ )
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      double precision   zero, one, two, three
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0)

*======================================================================

      goto ( 210, 220), k

 210  continue

      nonzro(1) = 0
      nonzro(2) = 1

      x2 = x(2)
      dfj( 1, 1) = zero
      dfj( 2, 1) = zero
      dfj( 3, 1) = zero
      dfj( 1, 2) = one
      dfj( 2, 2) = two*x2
      dfj( 3, 2) = three*x2*x2

      return

 220  continue

      nonzro(1) = 1
      nonzro(2) = 1

      x2 = x(2)
      t1 = two*x(1)
      t2 = three*x2
      dfj( 1, 1) = one
      dfj( 2, 1) = two*x2
      dfj( 3, 1) = t2*x2
      dfj( 1, 2) = zero
      dfj( 2, 2) = t1
      dfj( 3, 2) = t1*t2

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear
     
      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   t2, x2

      double precision   zero, one, two, three
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0)

*======================================================================

      linear = .false.

      goto ( 210, 220, 230), k

 210  continue

      hess( 1, 1) = zero
      hess( 1, 2) = one
      hess( 2, 1) = hess( 1, 2)
      hess( 2, 2) = zero

      return

 220  continue

      hess( 1, 1) = zero
      hess( 1, 2) = two*x(2)
      hess( 2, 1) = hess( 1, 2)
      hess( 2, 2) = two*x(1)

      return

 230  continue

      x2 = x(2)
      t2 = three*x2
      hess( 1, 1) = zero
      hess( 1, 2) = t2*x2
      hess( 2, 1) = hess( 1, 2)
      hess( 2, 2) = two*x(1)*t2

      return

      end
