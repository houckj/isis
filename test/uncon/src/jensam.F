*************************************************************************
*               jennrich and sampson function
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

      logical            lf, lj

      integer            na, nb, nc, nd, nt, nh

      integer            i, j

      double precision   di, e1, e2, x1, x2
    
      double precision   ddot

      intrinsic          dble, exp

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

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

      n      =  2
      m      = 10

      if (nu1 .gt. 0)  write( nu1, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  0.3d0
      x(2) =  0.4d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  0.2578d0
      x(2) =  0.2578d0

      ftf = 1.24362d+2

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        di   = dble(i)
        f(i) = two + two*di - (exp(di*x1) + exp(di*x2))
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        di = dble(i)
        fj(i,1) = -di*(exp(di*x1))
        fj(i,2) = -di*(exp(di*x2))
 210  continue

      return

 300  continue

      do 310 i = 1, m
        di = dble(i)
        e1 = exp(di*x1)
        e2 = exp(di*x2)
        f(i) = two + two*di - (e1 + e2)
        fj(i,1) = -di*e1
        fj(i,2) = -di*e2
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' jennrich and sampson function (more et al.)'//,
     *'        number of variables =', i4, '  (    2   )'/,
     *'        number of functions =', i4, '  (  >= 2  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      double precision   di 

      intrinsic          dble, exp     

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      j = 3 - k
      nonzro(j) = 0
      nonzro(k) = 1

      do 200 i = 1, m
        di = dble(i)
        dfj( i, j) =  zero
        dfj( i, k) = -di*di*exp(di*x(k))
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

      double precision   t, dk

      intrinsic          dble, exp

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      linear = .false.

      dk =  dble(k)
      t  = -dk*dk
      hess( 1, 1) = t*exp(dk*x(1))
      hess( 1, 2) = zero
      hess( 2, 1) = zero
      hess( 2, 2) = t*exp(dk*x(2))

      return

      end
