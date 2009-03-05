
      SUBROUTINE xsmtbl(ear, ne, param, filenm, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(*), photar(ne), photer(ne)
      CHARACTER*(*) filenm

c	xsmtbl	kaa		29 Mar 1990
c		XSPEC model subroutine
c		Multiplicative table model

c  10/2/94  kaa    modified for FITS format table models
C   9/9/96  kaa    added dynamic memory
c  10/18/97 han    added model error 

c Interpolate on the table. Set the value for any energies outside
c those tabulated to 1.0.

      LOGICAL qmerr

      CALL mtbint(ear, ne, param, filenm, ifl, photar, photer, 1.0,
     &            qmerr)

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE mtbint(ear, ne, param, filenm, ifl, photar, photer, 
     &                  defval, moderr)

      INTEGER ne, ifl
      REAL ear(0:ne), param(*), photar(ne), photer(ne)
      REAL defval
      LOGICAL moderr
      CHARACTER*(*) filenm


c Local variables

      INCLUDE 'xspec.inc'

      REAL zshift
      INTEGER ierr, npars, nbins, nesave, ivin

      LOGICAL qfirst, qnew, qrshft

      CHARACTER contxt*255, wrtstr*255, sfilenm*255

      INTEGER lenact
      EXTERNAL lenact

C pointers for dynamic arrays

      INTEGER istart, iend, ifstart, ifend
      INTEGER ibin, iefil

      SAVE istart, iend, ifstart, ifend, nesave, qfirst
      SAVE ibin, iefil, sfilenm

      DATA qfirst /.TRUE./

      ierr = 0

C If necessary (re)allocate memory for rebinning work arrays

      IF ( qfirst .OR. (ne .NE. nesave) ) THEN
         CALL udmget(ne, 4, istart, ierr)
         contxt = 'Failed to get memory for xsmtbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         CALL udmget(ne, 4, iend, ierr)
         contxt = 'Failed to get memory for xsmtbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         CALL udmget(ne, 6, ifstart, ierr)
         contxt = 'Failed to get memory for xsmtbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         CALL udmget(ne, 6, ifend, ierr)
         contxt = 'Failed to get memory for xsmtbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         nesave = ne
      ENDIF

C  check to see whether we are still dealing with the same table model

      qnew = .FALSE.
      IF ( qfirst ) qnew = .TRUE.
      IF ( filenm(:lenact(filenm)) .NE. sfilenm(:lenact(sfilenm)) ) 
     &      qnew = .TRUE.

c  interpolate on parameters to get model spectrum - the pointers ibin
c  and iefil are set within tblint.

      CALL tblint(1, param, filenm, qnew, npars, nbins, qrshft, ibin,
     &            ivin, iefil, moderr, ierr)

      IF ( ierr .NE. 0 ) RETURN

c  now rebin onto required energies

      zshift = 0.
      IF ( qrshft ) zshift = param(npars+1)

      CALL inibin(nbins, MEMR(iefil), ne, ear(0), MEMI(istart), 
     &            MEMI(iend), MEMR(ifstart), MEMR(ifend), zshift)

c  interpolate (NB - do not use erebin since interpolating not rebinning)

      CALL eintrp(MEMR(ibin), ne, MEMI(istart), MEMI(iend), 
     &            MEMR(ifstart), MEMR(ifend), photar, defval)

c if model has associated error, simply interpolate the variance,
c dont think a precise error propagation is worth the mess here 
c if there is no error then return -1 as the first element of the
c photer array.

      IF(moderr)THEN
         CALL eintrp(MEMR(ivin), ne, MEMI(istart), MEMI(iend), 
     &            MEMR(ifstart), MEMR(ifend), photer, 0.0)
      ELSE
         photer(1)  = -1.
      ENDIF


      qfirst = .FALSE.
      sfilenm = filenm(:min(len(filenm),len(sfilenm)))

 100  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         WRITE(wrtstr,'(''Error '', i3, '' in XSMTBL'')') ierr 
         CALL xwrite(wrtstr, 10)
         CALL xwrite(contxt, 10)
      ENDIF

      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE eintrp(bin, ne, start, end, fstart, fend, photar, 
     &                  defval)

      INTEGER ne, start(*), end(*)
      REAL bin(*), photar(*), fstart(*), fend(*)
      REAL defval


      INTEGER i, j
      REAL total


      DO i = 1, ne
         total = 0.
         IF (start(i).GT.0) THEN
            photar(i) = fstart(i)*bin(start(i))
            total = total + fstart(i)
            IF (end(i).GT.start(i)) THEN
               DO j = start(i) + 1, end(i) - 1
                  photar(i) = photar(i) + bin(j)
                  total = total + 1.
               ENDDO
               photar(i) = photar(i) + fend(i)*bin(end(i))
               total = total + fend(i)
            ENDIF
         ELSE
            photar(i) = defval
         ENDIF
         IF (total .GT. 0) THEN
            photar(i) = photar(i)/total
         ELSE
            photar(i) = defval
         ENDIF

      ENDDO

      RETURN
      END



