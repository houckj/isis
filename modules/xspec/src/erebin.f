      SUBROUTINE erebin(nin, inar, nout, istart, istop, fstart, fstop,
     &                  outar)
c		rashafer 1 april 1987
c	subroutine to perform a rebinning of an array, using simple
c	block interpolation.
c	note that many of the routines must be initialized by INIBIN
c
      INTEGER*4 nin
c                                I: The dimension of the input array
      REAL*4 inar(nin)
c                                I: input array
      INTEGER*4 nout
c                                I: The dimension of the output array
      INTEGER*4 istart(nout)
c                                I: The first bin
      INTEGER*4 istop(nout)
c                                I: The last bin
      REAL*4 fstart(nout)
c                                I: The fraction of the first bin
      REAL*4 fstop(nout)
c                                I: The fraction of the last bin
      REAL*4 outar(nout)
c                                R: The rebinned values.
c
      INTEGER*4 i, j
c
      DO i = 1, nout
         IF (istart(i).GT.0) THEN
            outar(i) = fstart(i)*inar(istart(i))
            IF (istop(i).GT.istart(i)) THEN
               DO j = istart(i) + 1, istop(i) - 1
                  outar(i) = outar(i) + inar(j)
               ENDDO
               outar(i) = outar(i) + fstop(i)*inar(istop(i))
            ENDIF
         ELSE
            outar(i) = 0.
         ENDIF
      ENDDO
      RETURN
      END
