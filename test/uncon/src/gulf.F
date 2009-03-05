*************************************************************************
*               gulf research and development function
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

      double precision   ai, av, bi, ci, d1, d2, d3, ei, ti, yi
      double precision   x1, x2, x3, x1inv

      intrinsic          abs, dble, exp, log

      double precision   ddot

      double precision   two3rd
      common /PARAM1/    two3rd
      save   /PARAM1/    

      double precision   zero, one, point1, twnty5, fifty
      parameter         (zero = 0.d0, one = 1.d0)
      parameter         (point1 = .01d0, twnty5 = 25.d0, fifty = 50.d0)

*=======================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      if (x1 .eq. zero) then
c        write( 6, 1)
c 1      format(/' +++ singularity in gulf function evaluation'/)
c        stop
        mode=999
        return
      endif
      x1inv = one / x1

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

      n     =  3
      m     = 10

      two3rd = 2.d0 / 3.d0

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =  5.d0
      x(2) =  2.5d0
      x(3) =  0.15d0

      return

*-----------------------------------------------------------------------

   30 continue

      x(1) =  50.d0
      x(2) =  25.d0
      x(3) =  1.5d0

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ti       = dble(i)*point1
        yi       = twnty5 + (-fifty*log(ti)) ** two3rd
        ai       = yi - x2
        f(i)     = exp(-((abs(ai)**x3)/x1)) - ti
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return
c
 200  continue
c
      do 210 i = 1, m
        ti = dble(i)*point1
        yi = twnty5 + (-fifty*log(ti)) ** two3rd
        ai = yi - x2
        av = abs(ai)
        bi = av ** x3
        ci = bi*x1inv
        ei = exp(-ci)
        d1 =   ci * x1inv
        d2 =   x3 * x1inv * av**(x3 - one)
        d3 = - log(av) * ci
        fj(i,1) = d1 * ei
        if (ai .ge. zero)  fj(i,2) =  d2 * ei
        if (ai .lt. zero)  fj(i,2) = -d2 * ei
        fj(i,3) = d3 * ei
 210  continue
      return

 300  continue

      do 310 i = 1, m
        ti = dble(i)*point1
        yi = twnty5 + (-fifty*log(ti)) ** two3rd
        ai = yi - x2
        av = abs(ai)
        bi = av**x3
        ci = bi*x1inv
        ei = exp(-ci)
        d1 =   ci * x1inv
        d2 =   x3 * x1inv * av**(x3 - one)
        d3 = - log(av) * ci
        f(i) = ei - ti
        fj(i,1) = d1 * ei
        if (ai .ge. zero)  fj(i,2) =  d2 * ei
        if (ai .lt. zero)  fj(i,2) = -d2 * ei
        fj(i,3) = d3 * ei
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue
      return

9999  format(/'1',70('=')//,
     *' gulf research and development function (more et al.)'//,
     *'        number of variables =', i4, '  (      3       )'/,
     *'        number of functions =', i4, '  ( >= 3, <= 100 )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      intrinsic          abs, dble, exp, log

      double precision   two3rd
      common /PARAM1/    two3rd
      save   /PARAM1/    

      double precision   zero, one, point1, twnty5, fifty
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = .01d0, twnty5 = 25.d0, fifty = 50.d0)

*=======================================================================

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      x1inv = one/x1

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 1

      goto ( 100, 200, 300), k

 100  continue

      do 110 i = 1, m
        ti = dble(i)*point1
        yi = twnty5 + (-fifty*log(ti)) ** two3rd
        ai = yi - x2
        av = abs(ai)
        b1 = av**(x3-one)
        bi = av**x3
        gi = log(av)
        ci = bi * x1inv
        ei = exp(-ci)
        d1 =   ci * x1inv
        d2 =   x3 * x1inv * b1
        d3 = - log(av) * ci
        d11 = - two * ci * x1inv * x1inv
        d12 = - x3 * b1 * x1inv * x1inv
        d13 =   gi * ci * x1inv
        dfj(i,1) = (d11 + d1*d1) * ei
        if (ai .ge. zero)  dfj(i,2) =  (d12 + d2*d2) * ei
        if (ai .lt. zero)  dfj(i,2) = -(d12 + d1*d2) * ei
        dfj(i,3) = (d13 + d1*d3) * ei
 110  continue

      return

 200  continue

      do 210 i = 1, m
        ti = dble(i)*point1
        yi = twnty5 + (-fifty*log(ti)) ** two3rd
        ai = yi - x2
        av = abs(ai)
        b2 = av**(x3-two)
        b1 = av**(x3-one)
        bi = av**x3
        gi = log(av)
        ci = bi * x1inv
        ei = exp(-ci)
        d1 =   ci * x1inv
        d2 =   x3 * x1inv * b1
        d3 = - log(av) * ci
        d21 = - x3 * b1 * x1inv * x1inv
        d22 = - x3 * (x3 - one) * b2 * x1inv
        d23 =   b1 * x1inv * (one + x3 * gi)
        if (ai .ge. zero)  dfj(i,1) =  (d21 + d2*d1) * ei
        if (ai .lt. zero)  dfj(i,1) = -(d21 + d2*d1) * ei
        dfj(i,2) =  (d22 + d2*d2) * ei
        if (ai .ge. zero)  dfj(i,3) =  (d23 + d2*d3) * ei
        if (ai .lt. zero)  dfj(i,3) = -(d23 + d2*d3) * ei
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ti = dble(i)*point1
        yi = twnty5 + (-fifty*log(ti)) ** two3rd
        ai = yi - x2
        av = abs(ai)
        b1 = av**(x3-one)
        bi = av**x3
        gi = log(av)
        ci = bi * x1inv
        ei = exp(-ci)
        d1 =   ci * x1inv
        d2 =   x3 * x1inv * b1
        d3 = - log(av) * ci
        d31 =   gi * ci * x1inv
        d32 =   b1 * x1inv * (one + x3 * gi)
        d33 = - gi * gi * ci
        dfj(i,1) =  (d31 + d3*d1) * ei
        if (ai .ge. zero)  dfj(i,2) =  (d32 + d3*d2) * ei
        if (ai .lt. zero)  dfj(i,2) = -(d32 + d3*d2) * ei
        dfj(i,3) =  (d33 + d3*d3) * ei
 310  continue

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      intrinsic          abs, dble, exp, log

      double precision   two3rd
      common /PARAM1/    two3rd
      save   /PARAM1/    

      double precision   zero, one, point1, twnty5, fifty
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (point1 = .01d0, twnty5 = 25.d0, fifty = 50.d0)

*=======================================================================

      linear = .false.

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      tk = dble(k)*point1
      yk = twnty5 + (-fifty*log(tk)) ** two3rd
      ak = yk - x2
      av = abs(ak)
      b2 = av**(x3 - two)
      b1 = av**(x3 - one)
      bk = av**x3
      if (x1 .eq. zero)  then
        write( 6, 1)
   1    format( ' second derivative does not exist')
        stop
      endif
      x1inv = one / x1
      ck = bk * x1inv
      ek = exp(-ck)
      gk = log(av)

      d1 =   ck * x1inv
      d2 =   x3 * x1inv * b1
      d3 = - gk * ck

      d11 = - two * ck * x1inv * x1inv
      d12 = - x3 * b1 * x1inv * x1inv
      d13 =   gk * ck * x1inv
      d22 = - x3 * (x3 - one) * b2 * x1inv
      d23 =   b1 * x1inv * (one + x3 * gk)
      d33 = - gk * gk * ck
      if (ak .ge. zero)  goto 100
        d2  = -d2
        d12 = -d12
        d23 = -d23

 100  continue

      hess( 1, 1) = (d11 + d1*d1) * ek
      hess( 1, 2) = (d12 + d1*d2) * ek
      hess( 1, 3) = (d13 + d1*d3) * ek
      hess( 2, 2) = (d22 + d2*d2) * ek
      hess( 2, 3) = (d23 + d2*d3) * ek
      hess( 3, 3) = (d33 + d3*d3) * ek
      hess( 2, 1) = hess( 1, 2)
      hess( 3, 1) = hess( 1, 3)
      hess( 3, 2) = hess( 2, 3)

      return
      end
