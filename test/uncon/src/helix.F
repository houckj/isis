***************************************************************************
*               helical valley function
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

      integer            j

      double precision   quo, theta

      intrinsic          atan, atan2, sqrt

      double precision   piinv
      common /PARAM1/    piinv
      save   /PARAM1/    

      double precision   zero, half, one, ten, fifty
      parameter         (zero = 0.d0, half = .5d0, one = 1.d0)
      parameter         (ten = 10.d0, fifty = 50.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

c      write( nu1, 2) x1
c   2  format(/' x1 = ', 1pe14.5)
      if (x1 .eq. zero)  then
c         stop
          mode=999
          return        
      endif

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

      n      = 3
      m      = 3

c     piinv  = one / acos(-one)
      piinv  = .25d0 / atan(one)

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) = -one
      x(2) =  zero
      x(3) =  zero

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  one
      x(2) =  zero
      x(3) =  zero

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      quo   = x2 / x1
c     theta = half*piinv*atan2(x2/x1)
      theta = half*piinv*atan(quo)
      if (x1 .lt. zero)  theta = theta + half

      f(1) = ten*(x3 -ten*theta)
      f(2) = ten*(sqrt(x1*x1 + x2*x2) - one)
      f(3) = x3

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      q = one / (x1*x1 + x2*x2)
      r = sqrt(q)

      fpq = fifty*piinv*q
      tr  = ten*r

      fj(1,1) =  fpq*x2
      fj(1,2) = -fpq*x1
      fj(1,3) =  ten
      fj(2,1) =  tr*x1
      fj(2,2) =  tr*x2
      fj(2,3) =  zero
      fj(3,1) =  zero
      fj(3,2) =  zero
      fj(3,3) =  one

      return

 300  continue

      quo   = x2 / x1
c     theta = half*piinv*atan2(x2,x1)
      theta = half*piinv*atan(quo)
      if (x1 .lt. zero)  theta = theta + half

      t = x1*x1 + x2*x2
      s = dsqrt(t)
      q = one / t
      r = dsqrt(q)

      f(1) = ten*(x3 -ten*theta)
      f(2) = ten*(s - one)
      f(3) = x3

      fpq = fifty*piinv*q
      tr  = ten*r

      fj(1,1) =  fpq*x2
      fj(1,2) = -fpq*x1
      fj(1,3) =  ten
      fj(2,1) =  tr*x1
      fj(2,2) =  tr*x2
      fj(2,3) =  zero
      fj(3,1) =  zero
      fj(3,2) =  zero
      fj(3,3) =  one

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' helical valley function (more et al.)'//,
     *'        number of variables =', i4, '  (3)'/,
     *'        number of functions =', i4, '  (3)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      double precision   piinv
      common /PARAM1/    piinv
      save   /PARAM1/    

      double precision   zero, one, two, ten, fifty
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (ten = 10.d0, fifty = 50.d0)

*=======================================================================

      nonzro(3) = 0

      dfj( 1, 3) =  zero
      dfj( 2, 3) =  zero
      dfj( 3, 1) =  zero
      dfj( 3, 2) =  zero
      dfj( 3, 3) =  zero

      if (k .eq. 3)  goto 230

      x1 = x(1)
      x2 = x(2)

      x1sq = x1*x1
      x2sq = x2*x2
      x12  = x1*x2

      nonzro(1) = 1
      nonzro(2) = 1

      q =  one / (x1sq + x2sq)
      d = q*q

      fpq = fifty*piinv*d
      td  = ten*dsqrt(q)*q

      goto ( 210, 220), k

 210  continue

      dfj( 1, 1) = -two*fpq*x12
      dfj( 1, 2) =  fpq*(x1sq - x2sq)
      dfj( 2, 1) =  td*x2sq
      dfj( 2, 2) = -td*x12

      return

 220  continue

      dfj( 1, 1) =  fpq*(x1sq - x2sq)
      dfj( 1, 2) =  two*fpq*x12
      dfj( 2, 1) = -td*x12
      dfj( 2, 2) =  td*x1sq

      return

 230  continue

      nonzro(1) = 0
      nonzro(2) = 0

      dfj( 1, 1) =  zero
      dfj( 1, 2) =  zero
      dfj( 2, 1) =  zero
      dfj( 2, 2) =  zero

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   piinv
      common /PARAM1/    piinv
      save   /PARAM1/    

      double precision   zero, one, two, ten, fifty
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (ten = 10.d0, fifty = 50.d0)

*=======================================================================

      goto ( 210, 220, 230), k

 210  continue

      linear = .false.

      x1 = x(1)
      x2 = x(2)

      x1sq = x1*x1
      x2sq = x2*x2

      q  = one / (x1sq + x2sq)
      d  = q*q

      fpq = fifty*piinv*d
      tmp = two*fpq*x1*x2

      hess( 1, 1) = -tmp
      hess( 1, 2) =  fpq*(x1sq - x2sq)
      hess( 2, 1) =  hess(1,2)
      hess( 1, 3) =  zero
      hess( 3, 1) =  zero
      hess( 2, 2) =  tmp
      hess( 2, 3) =  zero
      hess( 3, 2) =  zero
      hess( 3, 3) =  zero

      return

 220  continue

      linear = .false.

      x1 = x(1)
      x2 = x(2)

      x1sq = x1*x1
      x2sq = x2*x2

      q  = one / (x1sq + x2sq)
      d  = q*q

      td = ten*dsqrt(q)*q

      hess( 1, 1) =  td * x2sq
      hess( 1, 2) = -td * x1 * x2
      hess( 2, 1) =  hess(1,2)
      hess( 1, 3) =  zero
      hess( 3, 1) =  zero
      hess( 2, 2) =  td * x1sq
      hess( 2, 3) =  zero
      hess( 3, 2) =  zero
      hess( 3, 3) =  zero

      return

 230  continue

      linear = .true.

      hess( 1, 1) =  zero
      hess( 1, 2) =  zero
      hess( 2, 1) =  zero
      hess( 1, 3) =  zero
      hess( 3, 1) =  zero
      hess( 2, 2) =  zero
      hess( 2, 3) =  zero
      hess( 3, 2) =  zero
      hess( 3, 3) =  zero

      return

      end
