*************************************************************************
*               rosenbrock function
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

      double precision   ddot

      double precision   coeff
      common /PARAM1/    coeff
      save   /PARAM1/

      double precision   zero, one, two, ten, hundrd
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (ten = 10.d0, hundrd = 100.d0)

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
      nstrts = 1

      neven  = 5
      n      = 2*neven
      m      = n

      coeff = ten

      if (nout .gt. 0)  write( nout, 9999)  coeff, n, m

      return

*-----------------------------------------------------------------------

   20 continue

      do 21 i = 1, neven
        i2      =   2*i
        x(i2-1) = - 1.2d0
        x(i2)   =   one
   21 continue

      return

*-----------------------------------------------------------------------

   30 continue

      do 31 i = 1, n
        x(i) = one
   31 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, neven
        i2    = 2*i
        i2m1  = i2 - 1
        xi2   = x(i2)
        xi2m1 = x(i2m1)
        f(i2m1) = coeff * (xi2 - xi2m1*xi2m1)
        f(  i2) = one - xi2m1
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 210  continue

      do 220 i = 1, neven
        i2    = 2*i
        i2m1  = i2 - 1
        xi2   = x(i2)
        xi2m1 = x(i2m1)
        fi2m1 = coeff * (xi2 - xi2m1*xi2m1)
        fi2   = one - xi2m1
        fj(i2m1,i2m1) = -coeff * two * xi2m1
        fj(i2m1,  i2) =  coeff
        fj(  i2,i2m1) = -one
 220  continue

      return

 300  continue

      do 310 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 310  continue

      do 320 i = 1, neven
        i2    = 2*i
        i2m1  = i2 - 1
        xi2   = x(i2)
        xi2m1 = x(i2m1)
        f(i2m1) = coeff * (xi2 - xi2m1*xi2m1)
        f(  i2) = one - xi2m1
        fj(i2m1,i2m1) = -coeff * two * xi2m1
        fj(i2m1,  i2) =  coeff
        fj(  i2,i2m1) = -one
 320  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *'  extended rosenbrock function (more et al.);   coefficient = ',
     *  2pe9.1//,
     *'        number of variables =', i4,'  (multiple of 2)'/,
     *'        number of functions =', i4,'  (     =  n    )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            j

      common /PARAM1/    coeff
      save   /PARAM1/

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      if (mod(k,2) .eq. 0)  return

      nonzro(k) = 1

      dfj(k,k) = - coeff * two

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

      double precision   coeff
      common /PARAM1/    coeff
      save   /PARAM1/

      double precision   zero, two
      parameter         (zero = 0.d0, two = 2.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      if (mod(k,2) .eq. 0) then
        linear = .true.
        return
      end if

      hess(k,k) = - coeff * two

      return
      end

