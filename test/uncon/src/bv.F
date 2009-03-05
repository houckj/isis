*****************************************************************************
*               discrete boundary value function
* more, garbow, and hillstrom, acm toms vol. 7 no. 1 (march 1981) 17-41
*****************************************************************************

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

      integer            i, ip1, im1, j

      double precision   ri, si, ti, tj, t, xi, xip1, xim1 

      double precision   ddot

      intrinsic          dble

      double precision   h
      common /PARAM1/    h
      save   /PARAM1/

      double precision   zero, one, two, three
      parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter         (three = 3.d0)

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

      n      = 10
      m      = n

      n1 = n + 1

      h = one / dble(n1)

      if (nout .gt. 0)  write( nout, 9999) n, m

      return

*-----------------------------------------------------------------------

   20 continue

      do 21 j = 1, n
        tj   = dble(j)*h
        x(j) = tj*(tj-one)
  21  continue

      return

*-----------------------------------------------------------------------

   30 continue

      ftf = zero

      return

*-----------------------------------------------------------------------

 100  continue

      i  =   1
      xi = x(1)
      do 110 ip1 = 2, n1
        if (i .ne. n)  xip1 = x(ip1)
        t = two*xi + (h*h*(xi + dble(i)*h + one)**3)/two
        if (i .ne. 1)  t = t - xim1
        if (i .ne. n)  t = t - xip1
        f(i) = t
        im1  = i
        i    = ip1
        xim1 = xi
        xi   = xip1
 110  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      return

 200  continue

      do 210 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 210  continue

      i  =   1
      do 220 ip1 = 2, n1
        fj(i,i) = two + (three*h*h*(x(i) + dble(i)*h + one)**2)/two
        if (i .ne. 1)  fj(i,im1) = -one
        if (i .ne. n)  fj(i,ip1) = -one
        im1 = i
        i   = ip1
 220  continue

      return

 300  continue

      do 310 j = 1, n
        call dcopy( m, zero, 0, fj( 1, j), 1)
 310  continue

      i  =   1
      xi = x(1)
      do 320 ip1 = 2, n1
        if (i .ne. n)  xip1 = x(ip1)
        ti = dble(i)*h
        ri = (xi + ti + one)
        si = h*h*ri*ri/two
        t  = two*xi + si*ri
        if (i .ne.    1)  t = t - xim1
        if (i .ne. n)  t = t - xip1
        f(i)    = t
        fj(i,i) = two + three*si
        if (i .ne. 1)  fj(i,im1) = -one
        if (i .ne. n)  fj(i,ip1) = -one
        im1  = i
        i    = ip1
        xim1 = xi
        xi   =  xip1
 320  continue

      if (nb .ne. 0)  ftf = ddot( m, f, 1, f, 1)

      if (nd .eq. 0)  return

      do 400 j = 1, n
        g(j) = ddot( m, fj( 1, j), 1, f, 1)
 400  continue

      return

9999  format(/'1',70('=')//,
     *' discrete boundary value function (more et al.)'//,
     *'        number of variables =', i4,'  (variable)'/,
     *'        number of functions =', i4,'  (   = n  )'//,
     *        ' ',70('=')/)
      end

************************************************************************
************************************************************************

      subroutine dfjdxk( k, x, n, dfj, ldfj, m, nonzro)

      implicit double precision  (a-h,o-z)

      integer            k, n, ldfj, m, nonzro(n)

      double precision   x(n), dfj(ldfj,n)
      
      integer            j

      intrinsic          dble

      double precision   h
      common /PARAM1/    h
      save   /PARAM1/

      double precision   zero, one, three
      parameter         (zero = 0.d0, one = 1.d0)
      parameter         (three = 3.d0)

*=======================================================================

      do 100 j = 1, n
        nonzro(j) = 0
        call dcopy( m, zero, 0, dfj( 1, j), 1)
  100 continue

      nonzro(k) = 1

      dfj(k,k) = three*h*h*(x(k) + dble(k)*h + one)

      return
      end

************************************************************************
************************************************************************

      subroutine dfkdij( k, x, n, lhess, hess, linear)

      implicit double precision (a-h,o-z)

      logical            linear

      integer            k, n, lhess

      double precision   x(n), hess(lhess,n)

      integer            j

      intrinsic          dble

      double precision   h
      common /PARAM1/    h
      save   /PARAM1/

      double precision   zero, one, three
      parameter         (zero = 0.d0, one = 1.d0)
      parameter         (three = 3.d0)

*=======================================================================

      do 100 j = 1, n
        call dcopy( n, zero, 0, hess( 1, j), 1)
  100 continue

      linear = .false.

      hess(k,k) = three*h*h*(x(k) + dble(k)*h + one)

      return
      end
