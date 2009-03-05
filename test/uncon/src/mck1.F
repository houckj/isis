**************************************************************************
*           mckeown problem 1     ( math. prog. 9 (1975) 57-68 )
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

      integer            i, j

      double precision   rhomin, rhomax
      
      double precision   ddot

      double precision   hmu, a0, g0, bx
      common /PARAM1/    hmu, a0(3), g0(3,2), bx(3)
      save   /PARAM1/

      double precision   rmu, b0, d0
      common /PARAM2/    rmu, b0(2,2), d0(3)
      save   /PARAM2/

      double precision   rlamin, rlamax
      parameter         (rlamin = 0.606151d0, rlamax = 7.18396d0)

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

      if (lf .and. lj)  goto 100
      if (lf)           goto 100
      if (lj)           goto 200

*-----------------------------------------------------------------------

   10 continue

      numbr  = 1
      nprobs = 1
      nstrts = 1

      n      = 2
      m      = 3

      rmu    = 1.d0
      hmu    = .5d0 * rmu

      rhomax = rmu / (rlamin * rlamin)
      rhomin = rmu / (rlamax * rlamax)

      a0(2) = -0.244378d0
      a0(3) =  0.325895d0

      d0(1) =  2.5074d0
      d0(2) = -1.36401d0
      d0(3) =  1.02282d0

      g0(1,1) = -0.564255d0
      g0(1,2) =  0.392417d0
      g0(2,1) = -0.404979d0
      g0(2,2) =  0.927589d0
      g0(3,1) = -0.0735084d0
      g0(3,2) =  0.535493d0

      b0(1,1) =  5.66598d0
      b0(1,2) =  2.77141d0
      b0(2,1) =  2.77141d0
      b0(2,2) =  2.12413d0

      if (nout .gt. 0)  write( nout, 9999)  rmu,rhomin,rhomax,n,m

      return

*-----------------------------------------------------------------------

  20  continue

      x(1) = 0.1d0
      x(2) = 0.1d0

      return

*-----------------------------------------------------------------------

  30  continue

      ftf = 0.183601d0

      return

*-----------------------------------------------------------------------

 100  continue

      call dcopy( m, a0, 1, f, 1)

      do 110 i = 1, n
        call daxpy( m, x(i), g0( 1, i), 1, f, 1)
        bx(i) = ddot( n, b0( 1, i), 1, x, 1)
 110  continue

      call daxpy( m, (hmu*ddot( n, x, 1, bx, 1)), d0, 1, f, 1)

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (.not. lj)   return

 200  continue

      do 210 i = 1, m
        call dcopy( n, g0( i, 1), m, fj( i, 1), lfj)
        call daxpy( n, (rmu*d0(i)), bx, 1, fj( i, 1), lfj)
 210  continue

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *        '  mckeown problem 1 :   rmu = ', 1pe9.0/,
     *        '   rhomin =', 1pe11.2, '   rhomax =', 1pe11.2//,
     *        '  number of variables =', i4/,
     *        '  number of functions =', i4//
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            j

      double precision   rmu, b0, d0
      common /PARAM2/    rmu, b0(2,2), d0(3)
      save   /PARAM2/

      double precision   zero
      parameter         (zero = 0.d0)

*=======================================================================

      linear = .false.

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
        call daxpy( n, (rmu*d0(k)), b0( 1, j), 1, hess( 1, j), 1)
  100 continue

      return
      end
