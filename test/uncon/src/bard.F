*************************************************************************
*                 bard function
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
      double precision   di, qi, ui, vi, wi
     
      double precision   ddot
   
      intrinsic          dble, min

      double precision   y
      common /PARAM1/    y(15)
      save   /PARAM1/

      double precision   zero, one
      parameter         (zero = 0.d0, one = 1.d0)

*=============================================================================

      if (mode .eq.  0)  goto    20
      if (mode .eq. -1)  goto    10
      if (mode .eq. -2)  goto    30

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

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

      nprobs =  1
      nstrts =  1

      n      =  3
      m      = 15

      y( 1) = 0.14d0
      y( 2) = 0.18d0
      y( 3) = 0.22d0
      y( 4) = 0.25d0
      y( 5) = 0.29d0
      y( 6) = 0.32d0
      y( 7) = 0.35d0
      y( 8) = 0.39d0
      y( 9) = 0.37d0
      y(10) = 0.58d0
      y(11) = 0.73d0
      y(12) = 0.96d0
      y(13) = 1.34d0
      y(14) = 2.10d0
      y(15) = 4.39d0

      if (nout .gt. 0)  write( nout, 9999)  n, m

      return

*-----------------------------------------------------------------------

   20 continue

      x(1) = one
      x(2) = one
      x(3) = one

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = 8.21487d-3

      return

*-----------------------------------------------------------------------

 100  continue

      do 110 i = 1, m
        ui = dble(i)
        vi = dble(16-i)
        wi = min( ui, vi)
        di = vi*x2 + wi*x3
        if (di .eq. zero) then
c          write(6,*) '+++ getfun : attempt to divide by zero'
c          stop
          mode=999
          return
        end if
        f(i) = x1 + (ui/di) - y(i)
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 i = 1, m
        ui = dble(i)
        vi = dble(16-i)
        wi = min( ui, vi )
        di = vi*x2 + wi*x3
        if (di .eq. zero) then
c          write(6,*) '+++ getfun : attempt to divide by zero'
c          stop
          mode=999
          return          
        end if
        qi = -ui / (di*di)
        fj(i,1) = one
        fj(i,2) = qi * vi
        fj(i,3) = qi * wi
 210  continue

      return

 300  continue

      do 310 i = 1, m
        ui = dble(i)
        vi = dble(16-i)
        wi = min( ui, vi )
        di = vi*x2 + wi*x3
        if (di .eq. zero) then
c          write(6,*) '+++ getfun : attempt to divide by zero'
c          stop
          mode=999
          return          
        end if
        qi = -ui / (di*di)
        f(i) = x1 + (ui/di) - y(i)
        fj(i,1) = one
        fj(i,2) = qi * vi
        fj(i,3) = qi * wi
 310  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 320 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 320  continue

      return

9999  format(/'1',70('=')//,
     *' bard function (more et al.)'//,
     *'        number of variables =', i4, '  ( 3) '/,
     *'        number of functions =', i4, '  (15) '//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk ( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)

      integer            i, j

      double precision   x1, x2, x3
      double precision   di, qi, ti, ui, vi, wi, tv

      intrinsic          dble, min

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      goto ( 100, 200, 300), k

 100  continue

      do 110 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
 110  continue

      return

 200  continue

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      nonzro(1) = 0
      nonzro(2) = 1
      nonzro(3) = 1

      do 210 i = 1, m
        ui = dble(i)
        vi = dble(16-i)
        wi = min( ui, vi)
        di = vi*x2 + wi*x3
        if (di .eq. zero) then
c          write(6,*) '+++ getfun : attempt to divide by zero'
c          stop
          mode=999
          return          
        end if
        qi = one / (di*di)
        ti = two * di * qi * qi * ui
        tv = ti * vi
        dfj(i,1) = zero
        dfj(i,2) = tv * vi
        dfj(i,3) = tv * wi
 210  continue

      return

 300  continue

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      nonzro(1) = 0
      nonzro(2) = 1
      nonzro(3) = 1

      do 310 i = 1, m
        ui = dble(i)
        vi = dble(16-i)
        wi = min( ui, vi)
        di = vi*x2 + wi*x3
        if (di .eq. zero) then
c          write(6,*) '+++ getfun : attempt to divide by zero'
c          stop
          mode=999
          return          
        end if
        qi = one / (di*di)
        ti = two * di * qi * qi * ui
        tw = ti * wi
        dfj(i,1) = zero
        dfj(i,2) = tw * vi
        dfj(i,3) = tw * wi
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

      double precision   x1, x2, x3
      double precision   dk, qk, tk, uk, vk, wk, tw

      intrinsic          dble, min

      double precision   zero, one, two
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

*=======================================================================

      x1 = x(1)
      x2 = x(2)
      x3 = x(3)

      linear = .false.

      uk = dble(k)
      vk = dble(16-k)
      wk = min( uk, vk)
      dk = vk*x2 + wk*x3
      qk = one / (dk*dk*dk)
      tk = two * qk * uk
      tw = tk * wk

      hess( 1, 1) = zero
      hess( 1, 2) = zero
      hess( 2, 1) = zero
      hess( 2, 2) = tk * vk * vk
      hess( 1, 3) = zero
      hess( 3, 1) = zero
      hess( 2, 3) = tw * vk
      hess( 3, 2) = hess( 2, 3)
      hess( 3, 3) = tw * wk

      return
      end
