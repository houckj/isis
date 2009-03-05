*************************************************************************
*               meyer function
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

      double precision   x1, x2, x3
      double precision   di, ei, qi, si, ti

      double precision   ddot

      intrinsic          dble, exp

      double precision   y
      common /PARAM1/    y(16)
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

   10 continue

      nprobs = 1
      nstrts = 1

      n     =  3
      m     = 16

      y( 1) = 34780.d0
      y( 2) = 28610.d0
      y( 3) = 23650.d0
      y( 4) = 19630.d0
      y( 5) = 16370.d0
      y( 6) = 13720.d0
      y( 7) = 11540.d0
      y( 8) =  9744.d0
      y( 9) =  8261.d0
      y(10) =  7030.d0
      y(11) =  6005.d0
      y(12) =  5147.d0
      y(13) =  4427.d0
      y(14) =  3820.d0
      y(15) =  3307.d0
      y(16) =  2872.d0

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) =    0.02d0
      x(2) = 4000.d0
      x(3) =  250.d0

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 8.79458d+1

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ti = dble(45+5*i)
        di = ti + x3
        ei = exp(x2/di)
        f(i) = (x1 * ei) - y(i)
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)
      return

 200  continue

      do 210 i = 1, m
        ti = dble(45+5*i)
        di = ti + x3
        qi = one / di
        ei = exp(x2*qi)
        si = x1*qi*ei
        fj(i,1) =  ei
        fj(i,2) =  si
        fj(i,3) = -x2*qi*si
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ti = dble(45+5*i)
        di = ti + x3
        qi = one / di
        ei = exp(x2*qi)
        si = x1*qi*ei
        f(i) = (x1*ei) - y(i)
        fj(i,1) =  ei
        fj(i,2) =  si
        fj(i,3) = -x2*qi*si
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return
9999  format(/'1',70('=')//,
     *' meyer function (more et al.)'//,
     *'        number of variables =', i4, '  ( 3)'/,
     *'        number of functions =', i4, '  (16)'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)
  
      integer            i
    
      double precision   x1, x2, x3
      double precision   di, ei, ri, si, ti

      intrinsic          dble, exp

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      goto ( 100, 200, 300), k

 100  continue

      nonzro(1) = 0
      nonzro(2) = 1
      nonzro(3) = 1

      do 110 i = 1, m
        ti = dble(45+5*i)
        di = ti + x3
        qi = one / di
        ei = exp(x2*qi)
        si = ei*qi
        dfj(i,1) = zero
        dfj(i,2) = si
        dfj(i,3) = -x2*qi*si
 110  continue

      return

 200  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 1

      do 210 i = 1, m
        ti = dble(45+5*i)
        di = ti + x3
        qi = one / di
        ei = exp(x2*qi)
        si = ei*qi
        ri = x1*qi*si
        dfj(i,1) = si
        dfj(i,2) = ri
        dfj(i,3) = -ri*(one+x2*qi)
 210  continue

      return

 300  continue

      nonzro(1) = 1
      nonzro(2) = 1
      nonzro(3) = 1

      do 310 i = 1, m
        ti = dble(45+5*i)
        di = ti + x3
        qi = one / di
        ei =  exp(x2*qi)
        si = ei*(qi*qi)
        ri = x2*si
        dfj(i,1) = -ri
        dfj(i,2) = -x1*si*(one+x2*qi)
        dfj(i,3) =  x1*ri*qi*(two+x2*qi)
 310  continue

      return

      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, lhess, hess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      double precision   x1, x2, x3
      double precision   dk, ek, q2, qk, r1, r2, sk, tk

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      linear = .false.

      tk = dble(45+5*k)
      dk = tk + x3
      qk = one / dk
      q2 = qk * qk
      ek = exp(x2*qk)
      sk = ek*q2
      r1 = x1*sk
      r2 = x2*sk

      hess( 1, 1) = zero
      hess( 1, 2) = ek*qk
      hess( 2, 1) = hess( 1, 2)
      hess( 2, 2) =  r1
      hess( 1, 3) = -r2
      hess( 3, 1) = hess( 1, 3)
      hess( 2, 3) = -r1*(one+x2*qk)
      hess( 3, 2) = hess( 2, 3)
      hess( 3, 3) = x1*rk*qk*(two+x2*qk)

      return
      end

