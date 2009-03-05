******************************************************************************
*               linear function - rank one
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

      integer            i, j
    
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
      m      = 20

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

  20  continue

      call dcopy( n, one, 0, x, 1)

      return

*-----------------------------------------------------------------------

  30  continue

c     ftf = ( dble(m*(m-1))/dble(2*(2*m+1)) )
      ftf = 4.63414d0

      return

*-----------------------------------------------------------------------

 100  continue

      sum = zero
      do 110 j = 1, n
        sum = sum + x(j)*dble(j)
 110  continue

      do 120 i = 1, m
        f(i) = dble(i)*sum - one
 120  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        di = dble(i)
        do 210 j = 1, n
          fj(i,j) = di*dble(j)
 210  continue

      return

 300  continue

      sum = zero
      do 310 j = 1, n
        sum = sum + x(j)*dble(j)
 310  continue

      do 320 i = 1, m
        di = dble(i)
        f(i) = di*sum - one
        do 320 j = 1, n
          fj(i,j) = di*dble(j)
 320  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' linear function - rank one (more et al.)'//,
     *'        number of variables =', i4,'  (variable)'/,
     *'        number of functions =', i4,'  (  >= n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      return
      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            j

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .true.

      return
      end
