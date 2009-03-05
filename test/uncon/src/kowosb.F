*************************************************************************
*               kowalik and osborne function
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

      double precision   x1, x2, x3, x4

      double precision   ddot

      double precision   u, y
      common /PARAM1/    u(11), y(11)
      save   /PARAM1/

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

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

  10  continue

      nprobs = 1
      nstrts = 1

      n      =  4
      m      = 11

      y( 1) = 0.1957d0
      y( 2) = 0.1947d0
      y( 3) = 0.1735d0
      y( 4) = 0.1600d0
      y( 5) = 0.0844d0
      y( 6) = 0.0627d0
      y( 7) = 0.0456d0
      y( 8) = 0.0342d0
      y( 9) = 0.0323d0
      y(10) = 0.0235d0
      y(11) = 0.0246d0

      u( 1) = 4.0000d0
      u( 2) = 2.0000d0
      u( 3) = 1.0000d0
      u( 4) = 0.5000d0
      u( 5) = 0.2500d0
      u( 6) = 0.1670d0
      u( 7) = 0.1250d0
      u( 8) = 0.1000d0
      u( 9) = 0.0833d0
      u(10) = 0.0714d0
      u(11) = 0.0625d0

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) = 0.25d0
      x(2) = 0.39d0
      x(3) = 0.415d0
      x(4) = 0.39d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 3.07505d-4

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ui = u(i)
        t2 = ui*(ui + x2)
        di = ui*(ui + x3) + x4
        f(i) = (x1*t2)/di - y(i)
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue
      do 210 i = 1, m
        ui = u(i)
        t2 = ui*(ui + x2)
        di = ui*(ui + x3) + x4
        qi = one / di
        r  =  t2*qi
        s  =  x1*qi
        v  = -r*s
        fj(i,1) = r
        fj(i,2) = s*ui
        fj(i,3) = v*ui
        fj(i,4) = v
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ui = u(i)
        t2 = ui*(ui + x2)
        di = ui*(ui + x3) + x4
        qi = one / di
        r  =  t2*qi
        s  =  x1*qi
        v  = -r*s
        f(i) = x1*r - y(i)
        fj(i,1) = r
        fj(i,2) = s*ui
        fj(i,3) = v*ui
        fj(i,4) = v
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' kowalik and osborne function (more et al.)'//,
     *'        number of variables =', i4,'  ( 4)'/,
     *'        number of functions =', i4,'  (11)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i

      double precision   x1, x2, x3, x4

      double precision   u, y
      common /PARAM1/    u(11), y(11)
      save   /PARAM1/

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)

      goto ( 100, 200, 300, 400), k

 100  continue

      nonzro(1) = 0
      nonzro(2) = 1
      nonzro(3) = 1
      nonzro(4) = 1

      do 110 i = 1, m
        ui = u(i)
        t2 = ui*(ui + x2)
        di = ui*(ui + x3) + x4
        qi = one / di
        s  = -t2*qi*qi
        dfj(i,1) = zero
        dfj(i,2) = ui*qi
        dfj(i,3) = s*ui
        dfj(i,4) = s
 110  continue

      return

 200  continue

      nonzro(1) = 1
      nonzro(2) = 0
      nonzro(3) = 1
      nonzro(4) = 1

      do 210 i = 1, m
        ui = u(i)
        t2 = ui*(ui + x2)
        di = ui*(ui + x3) + x4
        qi = one / di
        r  =  ui*qi
        s  = -x1*r*qi
        dfj(i,1) = r
        dfj(i,2) = zero
        dfj(i,3) = s*ui
        dfj(i,4) = s
 210  continue

      return

 300  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 1
      nonzro(4) = 1

      do 310 i = 1, m
        ui = u(i)
        t2 = ui*(ui + x2)
        di = ui*(ui + x3) + x4
        qi = one / di
        q  = qi*qi*ui
        r  = x1*ui
        s  = t2*q
        w  = x1*qi*s*two
        dfj(i,1) = -s
        dfj(i,2) = -q*r
        dfj(i,3) =  w*ui
        dfj(i,4) =  w
 310  continue

      return

 400  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 1
      nonzro(4) = 1

      do 410 i = 1, m
        ui = u(i)
        t2 = ui*(ui + x2)
        di = ui*(ui + x3) + x4
        qi = one / di
        q  = t2*qi*qi
        r  = x1*qi
        s  = ui*qi
        v  = q*r*two
        dfj(i,1) = -q
        dfj(i,2) = -r*s
        dfj(i,3) =  v*ui
        dfj(i,4) =  v
 410  continue

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij ( k, x, n, hess, lhess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   u, y
      common /PARAM1/    u(11), y(11)
      save   /PARAM1/

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      linear = .false.

      uk = u(k)
      t2 = uk*(uk + x(2))
      dk = uk*(uk + x(3)) + x(4)
      qk = one / dk
      q  = uk*qk
      r  = t2*qk
      s  = x(1)*qk
      v  = -r*qk
      w  = -q*s
      z  = two*r*s*qk
      t  = z*uk

      hess( 1, 1) = zero
      hess( 1, 2) = q
      hess( 2, 1) = hess( 1, 2)
      hess( 1, 3) = v*uk
      hess( 3, 1) = hess( 1, 3)
      hess( 1, 4) = v
      hess( 4, 1) = hess( 1, 4)
      hess( 2, 2) = zero
      hess( 2, 3) = w*uk
      hess( 3, 2) = hess( 2, 3)
      hess( 2, 4) = w
      hess( 4, 2) = hess( 2, 4)
      hess( 3, 3) = t*uk
      hess( 3, 4) = t
      hess( 4, 3) = hess( 3, 4)
      hess( 4, 4) = z

      return
      end
