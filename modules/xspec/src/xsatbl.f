
      SUBROUTINE xsatbl(ear, ne, param, filenm, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(*), photar(ne), photer(*)
      CHARACTER*(*) filenm

c	xsatbl	kaa		29 Mar 1990
c		XSPEC model subroutine
c		Additive table model

c  10/2/94  kaa    modified for FITS format table models
C   9/9/96  kaa    added dynamic memory

c Local variables

      INCLUDE 'xspec.inc'

      INTEGER i, ierr, npars, nbins, nesave
      REAL zfac, zshift

      LOGICAL qfirst, qnew, qrshft, moderr

      CHARACTER contxt*255, wrtstr*255, sfilenm*255

      INTEGER lenact
      EXTERNAL lenact

C pointers for dynamic arrays

      INTEGER istart, iend, ifstart, ifend
      INTEGER ibin, ivin, iefil

      SAVE istart, iend, ifstart, ifend, nesave, qfirst
      SAVE ibin, ivin, iefil, sfilenm

      DATA qfirst /.TRUE./

      ierr = 0

C If necessary (re)allocate memory for rebinning work arrays

      IF ( qfirst .OR. (ne .NE. nesave) ) THEN
         CALL udmget(ne, 4, istart, ierr)
         contxt = 'Failed to get memory for xsatbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         CALL udmget(ne, 4, iend, ierr)
         contxt = 'Failed to get memory for xsatbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         CALL udmget(ne, 6, ifstart, ierr)
         contxt = 'Failed to get memory for xsatbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         CALL udmget(ne, 6, ifend, ierr)
         contxt = 'Failed to get memory for xsatbl work array'
         IF ( ierr .NE. 0 ) GOTO 100
         nesave = ne
      ENDIF

C  check to see whether we are still dealing with the same table model

      qnew = .FALSE.
      IF ( qfirst ) qnew = .TRUE.
      IF ( filenm(:lenact(filenm)) .NE. sfilenm(:lenact(sfilenm)) ) 
     &      qnew = .TRUE.

C  if the first time through then make sure that the iefil and ibin pointers
C  are set to something that will not be interpreted as a real address

      IF ( qfirst ) THEN
         iefil = -1
         ibin  = -1
         ivin  = -1
      ENDIF

c  interpolate on parameters to get model spectrum - the pointers ibin
c  and iefil are set within tblint.


      CALL tblint(2, param, filenm, qnew, npars, nbins, qrshft, ibin,
     &            ivin, iefil, moderr, ierr)
      IF ( ierr .NE. 0 ) RETURN

c  now rebin onto the correct energies

      zshift = 0.
      IF ( qrshft ) zshift = param(npars+1)

      CALL inibin(nbins, MEMR(iefil), ne, ear(0), MEMI(istart), 
     &            MEMI(iend), MEMR(ifstart), MEMR(ifend), zshift)

      CALL erebin(nbins, MEMR(ibin), ne, MEMI(istart), MEMI(iend), 
     &            MEMR(ifstart), MEMR(ifend), photar)

c if model has errors associated, then bin the variance vector

      IF ( moderr ) THEN
         DO i = 0, ne-1
            MEMR(ifstart+i) = MEMR(ifstart+i)**2 
            MEMR(ifend+i)   = MEMR(ifend+i)**2 
         ENDDO
         CALL erebin(nbins, MEMR(ivin), ne, MEMI(istart), MEMI(iend), 
     &               MEMR(ifstart), MEMR(ifend), photer)
      ENDIF

c  If redshift is non-zero then put in time dilation

      IF ( qrshft ) THEN
         zfac = zshift + 1.0
         IF (zfac .GT. 1.) THEN
            DO i = 1, ne
               photar(i) = photar(i)/zfac
            ENDDO
            IF ( moderr ) THEN
               zfac=zfac**2
               DO i = 1, ne
                  photer(i) = photer(i)/zfac
               ENDDO
            ENDIF
         ENDIF
      ENDIF

      qfirst = .FALSE.
      sfilenm = filenm(:min(len(filenm),len(sfilenm)))

 100  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         WRITE(wrtstr,'(''Error '', i3, '' in XSATBL'')') ierr 
         CALL xwrite(wrtstr, 10)
         CALL xwrite(contxt, 10)
      ENDIF


      RETURN
      END


