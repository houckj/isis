****************************************************************************
*           mckeown problem 3     ( math. prog. 9 (1975) 57-68 )
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

      double precision   rhomin, rhomax

      double precision   ddot

      double precision   hmu, a0, g0, bx
      common /PARAM1/    hmu, a0(10), g0(10,5), bx(5)
      save   /PARAM1/

      double precision   rmu, b0, d0
      common /PARAM2/    rmu, b0(5,5), d0(10)
      save   /PARAM2/

      double precision   rlamin, rlamax
      parameter         (rlamin = 0.126d0, rlamax = 1.190d0)

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

  10  continue

      numbr  = 1
      nprobs = 1
      nstrts = 1

      n     =  5
      m     = 10

      rmu    = 1.d0
      hmu    = .5d0 * rmu

      rhomax = rmu / (rlamin * rlamin)
      rhomin = rmu / (rlamax * rlamax)

      a0( 1) =  0.0426149d0
      a0( 2) =  0.0352053d0
      a0( 3) =  0.0878058d0
      a0( 4) =  0.0330812d0
      a0( 5) =  0.0580924d0
      a0( 6) =  0.649704d0
      a0( 7) =  0.344144d0
      a0( 8) = -0.627443d0
      a0( 9) =  0.001828d0
      a0(10) = -0.224783d0

      d0( 1) =  2.34659d0
      d0( 2) =  2.84048d0
      d0( 3) =  1.13888d0
      d0( 4) =  3.02286d0
      d0( 5) =  1.72139d0
      d0( 6) =  0.153917d0
      d0( 7) =  0.290577d0
      d0( 8) = -0.159378d0
      d0( 9) = 54.6910d0
      d0(10) = -0.444873d0

      g0( 1,1) = -0.564255d0
      g0( 1,2) =  0.392417d0
      g0( 1,3) = -0.404979d0
      g0( 1,4) =  0.927589d0
      g0( 1,5) = -0.0735083d0
      g0( 2,1) =  0.535493d0
      g0( 2,2) =  0.658799d0
      g0( 2,3) = -0.636666d0
      g0( 2,4) = -0.681091d0
      g0( 2,5) = -0.869487d0
      g0( 3,1) =  0.586387d0
      g0( 3,2) =  0.289826d0
      g0( 3,3) =  0.854402d0
      g0( 3,4) =  0.789312d0
      g0( 3,5) =  0.949721d0
      g0( 4,1) =  0.608734d0
      g0( 4,2) =  0.984915d0
      g0( 4,3) =  0.375699d0
      g0( 4,4) =  0.239547d0
      g0( 4,5) =  0.463136d0
      g0( 5,1) =  0.774227d0
      g0( 5,2) =  0.325421d0
      g0( 5,3) = -0.151719d0
      g0( 5,4) =  0.448051d0
      g0( 5,5) =  0.149926d0
      g0( 6,1) = -0.435033d0
      g0( 6,2) = -0.688583d0
      g0( 6,3) =  0.222278d0
      g0( 6,4) = -0.524653d0
      g0( 6,5) =  0.413248d0
      g0( 7,1) =  0.759468d0
      g0( 7,2) = -0.627795d0
      g0( 7,3) =  0.0403142d0
      g0( 7,4) =  0.724666d0
      g0( 7,5) = -0.0182537d0
      g0( 8,1) = -0.152448d0
      g0( 8,2) = -0.546437d0
      g0( 8,3) =  0.484134d0
      g0( 8,4) =  0.353951d0
      g0( 8,5) =  0.887866d0
      g0( 9,1) = -0.821772d0
      g0( 9,2) = -0.53412d0
      g0( 9,3) = -0.798498d0
      g0( 9,4) = -0.658572d0
      g0( 9,5) =  0.662362d0
      g0(10,1) =  0.819831d0
      g0(10,2) = -0.910632d0
      g0(10,3) = -0.480344d0
      g0(10,4) = -0.871758d0
      g0(10,5) = -0.978666d0

      b0(1,1) =  0.354033d0
      b0(1,2) = -0.0230349d0
      b0(1,3) = -0.211938d0
      b0(1,4) = -0.0554288d0
      b0(1,5) =  0.220429d0
      b0(2,1) = -0.0230349d0
      b0(2,2) =  0.29135d0
      b0(2,3) = -0.00180333d0
      b0(2,4) = -0.111141d0
      b0(2,5) =  0.0485461d0
      b0(3,1) = -0.211938d0
      b0(3,2) = -0.00180333d0
      b0(3,3) =  0.815808d0
      b0(3,4) = -0.133538d0
      b0(3,5) = -0.38067d0
      b0(4,1) = -0.0554288d0
      b0(4,2) = -0.111141d0
      b0(4,3) = -0.133538d0
      b0(4,4) =  0.389198d0
      b0(4,5) = -0.131586d0
      b0(5,1) =  0.220429d0
      b0(5,2) =  0.0485461d0
      b0(5,3) = -0.38067d0
      b0(5,4) = -0.131586d0
      b0(5,5) =  0.534706d0

      if (nout .gt. 0)  write( nout, 9999)  rmu,rhomin,rhomax,n,m

      return

*-----------------------------------------------------------------------

  20  continue

      x(1) = 0.1d0
      x(2) = 0.1d0
      x(3) = 0.1d0
      x(4) = 0.1d0
      x(5) = 0.1d0

      return

*-----------------------------------------------------------------------

  30  continue

      ftf = 1.0000d0

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
        call dcopy( n, g0( i, 1), m, fj( i, 1), m)
        call daxpy( n, (rmu*d0(i)), bx, 1, fj( i, 1), m)
 210  continue

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *        '  mckeown problem 3 :   rmu = ', 1pe9.0/,
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
      common /PARAM2/    rmu, b0(5,5), d0(10)
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
