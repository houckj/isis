
      SUBROUTINE tblint(itype, param, filenm, qhead, nparam, nbins, 
     &                  qshift, ibino, ivino, iefilo, qmerr, ierr)

      INCLUDE 'xspec.inc'

      REAL param(*)
      CHARACTER*(*) filenm
      INTEGER itype, ierr, nparam, nbins, ibino, ivino, iefilo
      LOGICAL qhead, qshift, qmerr

c  tblint		kaa  4/3/90  - algorithm by Mike Corcoran
c				     - generalized and improved by kaa
c                       kaa 6/30/94  - trapped out case of only one tabulated
c                                      value for a parameter, mainly of use
c                                      when the user wants to input a single
c                                      spectrum
c                       kaa 10/2/94  - FITS format version
C                       kaa  9/9/96  - added dynamic memory
C                       kaa 4/25/97  - fixed so works correctly when both
C                                      mtable and atable are in use
c
c  routine to generate interpolated model spectrum for a table model.

c       i       itype           i:1=mtable, 2=atable
c	r	param(*)	i:model parameter values
c	c*	filenm		i:table model filename
c       l       qhead           i:if true then load the header information
c       i       nparam          r:the total number of parameters (int+add)
c       i       nbins           r:the number of spectral bins
c       l       qshift          r:whether a redshift parameter should be included
c       i       ibino           r:pointer for the interpolated model spectrum
c       i       ivino           r:pointer for the interpolated model variance
c                                 negative, if no error
C       i       iefilo          r:pointer for the tabulated energies
c	i	ierr		r:error flag:
c				 0 - successful
c				-1 - a parameter lies below tabulated values
c				+1 - a parameter lies above tabulated values
c                               other - a FITSIO error

C  pointers for local dynamic arrays
C     imethod       interpolation method for interpolation parameter
C     inump         number of values for interpolation parameter
C     ivalues       tabulated values for interpolation parameters
C     inrec         number of records in each parameter block
C     inbmin        work array
C     iintpar       parameters that need to be interpolated
C     iexact        interpolation parameters that are exact
C     irec          array of record numbers to be used in interpolation
C     iworkphot     work array storing interpolation steps
C     ibin          array of interpolated spectrum
C     iefil         array of tabulated energies
C     imerr         model error avaiable

      INTEGER imethod(2), inump(2), ivalues(2), inrec(2), inbmin(2)
      INTEGER iintpar(2), iexact(2), irec(2), iworkphot
      INTEGER ibin(2), ivin(2), iefil(2), imerr(2)

C Saved local parameters
C     npar          number of tabulated interpolation parameters
C     cntpar        number of parameters used in interpolation
C     naddpr        number of additional parameters
C     nval          number of tabulated energy bins
C     qrshft        whether redshift is to be included as a parameter
C     ilun          i/o channel used to read the table file 

      INTEGER npar(2), cntpar(2), naddpr(2), nval(2), ilun(2)
      LOGICAL qrshft(2)

      INTEGER i, j, index, icol, ntot, ioff, nsize
      
      LOGICAL qanyf, terr(2)

      CHARACTER*255 wrtstr
      CHARACTER*72 contxt, comment

      INTEGER lenact
      EXTERNAL lenact

      SAVE imethod, inump, ivalues, inrec, inbmin, iintpar, iexact
      SAVE irec, ibin, ivin, iefil, imerr
      SAVE npar, cntpar, naddpr, nval, ilun, qrshft, terr

      DATA ilun /2*0/
      DATA imethod, inump, ivalues, inrec, inbmin, iintpar /12*-1/
      DATA iexact, irec, ibin, iefil, imerr /10*-1/

      ierr = 0

c If required then read the header information and open the file

      IF ( qhead ) THEN

c If a file is already open then close it

         IF (ilun(itype) .NE. 0) THEN

            CALL ftclos(ilun(itype), ierr)
            contxt = 'TBLINT: Failed to close FITS file'
            IF ( ierr .NE. 0 ) GOTO 999
            CALL frelun(ilun(itype))

         ENDIF

c Open the FITS file

         CALL getlun(ilun(itype))
         CALL xsftop(ilun(itype), filenm, 0, index, ierr)
         contxt = 'TBLINT: Failed to open FITS file '//
     &               filenm(:MIN(lenact(filenm),47))
         IF ( ierr .NE. 0 ) GOTO 999

         wrtstr = 'Opening '//filenm(:MIN(lenact(filenm),247))
         CALL xwrite(wrtstr, 25)

c  read the redshift flag

         CALL ftgkyl(ilun(itype), 'REDSHIFT', qrshft(itype), comment, 
     &               ierr)
         contxt = 'TBLINT: Failed to read REDSHIFT keyword'
         IF ( ierr .NE. 0 ) GOTO 999

c Go to the parameter extension

         CALL ftmahd(ilun(itype), 2, index, ierr)
         contxt = 'TBLINT: Failed to go to second extension'
         IF ( ierr .NE. 0 ) GOTO 999

c Get the number of interpolated and additive parameters

         CALL ftgkyj(ilun(itype), 'NINTPARM', npar(itype), comment, 
     &               ierr)
         contxt = 'TBLINT: Failed to read NINTPARM keyword'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL ftgkyj(ilun(itype), 'NADDPARM', naddpr(itype), comment, 
     &               ierr)
         contxt = 'TBLINT: Failed to read NADDPARM keyword'
         IF ( ierr .NE. 0 ) GOTO 999

C get the memory for the interpolation methods

         CALL udmget(npar(itype), 4, imethod(itype), ierr)
         contxt = 
     &      'TBLINT: Failed to get memory for interpolation methods'
         IF ( ierr .NE. 0 ) GOTO 999

c Get the interpolation methods

         CALL ftgcno(ilun(itype), .true., 'METHOD', icol, ierr)
         contxt = 'Failed to find METHOD column'
         IF ( ierr .NE. 0 ) GOTO 999
         CALL ftgcvj(ilun(itype), icol, 1, 1, npar(itype), 0, 
     &               MEMI(imethod(itype)), qanyf, ierr)
         contxt = 
     &    'TBLINT: Failed to read the METHOD column in extension 2'
         IF ( ierr .NE. 0 ) GOTO 999

C get the memory for the number of parameter values

         CALL udmget(npar(itype), 4, inump(itype), ierr)
         contxt = 
     &    'TBLINT: Failed to get memory for number of parameter values'
         IF ( ierr .NE. 0 ) GOTO 999

c Get the number of parameter values

         CALL ftgcno(ilun(itype), .true., 'NUMBVALS', icol, ierr)
         contxt = 'Failed to find NUMBVALS column'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL ftgcvj(ilun(itype), icol, 1, 1, npar(itype), 0, 
     &               MEMI(inump(itype)), qanyf, ierr)
         contxt = 
     &    'TBLINT: Failed to read the NUMBVALS column in extension 2'
         IF ( ierr .NE. 0 ) GOTO 999

C get the memory for the parameter values

         ntot = 0
         DO i = 0, npar(itype)-1
            ntot = ntot + MEMI(inump(itype)+i)
         ENDDO

         CALL udmget(ntot, 6, ivalues(itype), ierr)
         contxt = 
     &    'TBLINT: Failed to get memory for parameter values'
         IF ( ierr .NE. 0 ) GOTO 999

c and the values themselves

         CALL ftgcno(ilun(itype), .true., 'VALUE', icol, ierr)
         contxt = 'Failed to find VALUE column'
         IF ( ierr .NE. 0 ) GOTO 999

         ioff = 0
         DO i = 1, npar(itype)

            CALL ftgcve(ilun(itype), icol, i, 1, 
     &                  MEMI(inump(itype)+i-1), 1.e-31, 
     &                  MEMR(ivalues(itype)+ioff), qanyf, ierr)
            contxt = 
     &       'TBLINT: Failed to read the VALUE column in extension 2'
            IF ( ierr .NE. 0 ) GOTO 999

            ioff = ioff + MEMI(inump(itype)+i-1)

         ENDDO

C Get the memory for the number of records in each parameter block

         CALL udmget(npar(itype), 4, inrec(itype), ierr)
         contxt = 
     &    'TBLINT: Failed to get memory for nrec array'
         IF ( ierr .NE. 0 ) GOTO 999

c Calculate the number of records in each parameter block

         DO j = 0, npar(itype)-1
            MEMI(inrec(itype)+j) = 1
         ENDDO
         DO j = 1, npar(itype)
            DO i = j + 1, npar(itype)
               MEMI(inrec(itype)+j-1) = MEMI(inump(itype)+i-1)
     &                                 *MEMI(inrec(itype)+j-1)
            ENDDO
         ENDDO

C Get the memory for other work arrays for interpolated parameters

         CALL udmget(npar(itype), 4, inbmin(itype), ierr)
         contxt = 
     &    'TBLINT: Failed to get memory for nbmin array'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL udmget(npar(itype), 4, iintpar(itype), ierr)
         contxt = 
     &    'TBLINT: Failed to get memory for intpar array'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL udmget(npar(itype), 1, iexact(itype), ierr)
         contxt = 
     &    'TBLINT: Failed to get memory for exact array'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL udmget(2**npar(itype), 4, irec(itype), ierr)
         contxt = 
     &      'TBLINT: Failed to get memory for record number array'
         IF ( ierr .NE. 0 ) GOTO 999

c Get the number of bins in the spectra

         CALL ftmahd(ilun(itype), 3, index, ierr)
         contxt = 'TBLINT: Failed to go to third extension'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL ftgkyj(ilun(itype), 'NAXIS2', nval(itype), comment, ierr)
         contxt = 'TBLINT: Failed to read NAXIS2 keyword in extension 3'
         IF ( ierr .NE. 0 ) GOTO 999

C Get the memory for the energies

         CALL udmget((nval(itype)+1), 6, iefil(itype), ierr)
         contxt = 'TBLINT: Failed to get memory for tabulated energies'
         IF ( ierr .NE. 0 ) GOTO 999

c and load the energies into the efil array

         CALL ftgcve(ilun(itype), 1, 1, 1, nval(itype), 0., 
     &               MEMR(iefil(itype)), qanyf, ierr) 
         contxt = 'TBLINT: Failed to read ENERG_LO column'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL ftgcve(ilun(itype), 2, nval(itype), 1, 1, 0., 
     &               MEMR(iefil(itype)+nval(itype)), qanyf, ierr) 
         contxt = 'TBLINT: Failed to read ENERG_HI column'
         IF ( ierr .NE. 0 ) GOTO 999

C Get the memory for the total interpolated spectrum

         CALL udmget(nval(itype), 6, ibin(itype), ierr)
         contxt = 
     &      'TBLINT: Failed to get memory for interpolated spectrum'
         IF ( ierr .NE. 0 ) GOTO 999

c Move to the fourth extension to be ready to read the spectra

         CALL ftmahd(ilun(itype), 4, index, ierr)
         contxt = 'TBLINT: Failed to go to fourth extension'
         IF ( ierr .NE. 0 ) GOTO 999

C Get the memory for the error indicator 

         CALL udmget(naddpr(itype)+1, 1, imerr(itype), ierr)
         contxt = 'TBLINT: Failed to get memory for model error'
         IF ( ierr .NE. 0 ) GOTO 999

         CALL merr_info(ilun(itype), nval(itype), naddpr(itype), 
     &                  MEMB(imerr(itype)), terr(itype), ierr) 
         contxt = 'TBLINT: Failed to get error information'
         IF ( ierr .NE. 0 ) GOTO 999

         IF ( terr(itype) ) THEN
           CALL udmget(nval(itype), 6, ivin(itype), ierr)
           contxt = 
     &      'TBLINT: Failed to get memory for interpolated spectrum'
           IF ( ierr .NE. 0 ) GOTO 999
         ELSE
           ivin(itype)=ibin(itype)
         ENDIF

      ENDIF

C Set up the arrays used to determine which records are used for the
C interpolation. Also checks the input parameters are within the bounds
C of the tabulated values.

      CALL intstp(npar(itype), MEMI(inump(itype)), param, 
     &            MEMR(ivalues(itype)), MEMI(inrec(itype)), 
     &            MEMB(iexact(itype)), MEMI(inbmin(itype)), 
     &            MEMI(iintpar(itype)), cntpar(itype), 
     &            MEMI(irec(itype)), ierr)
      contxt = 'TBLINT: Error detected in INTSTP'
      IF ( ierr .NE. 0 ) GOTO 999

C Get the memory for the workphot array - really a 2-D array with
C first axis being parameter values followed by spectrum and variance 
C if any.

      nsize = npar(itype)+nval(itype)
      if ( terr(itype) ) nsize = nsize + nval(itype) 
      iworkphot = -1
      CALL udmget(nsize*(cntpar(itype)+1), 6, 
     &            iworkphot, ierr)
      contxt = 
     &   'TBLINT: Failed to get memory for workphot array'
      IF ( ierr .NE. 0 ) GOTO 999

C Now do the actual calculation of the interpolated spectrum

      CALL grdint(ilun(itype), cntpar(itype), npar(itype), nval(itype),
     &            naddpr(itype), param, MEMI(iintpar(itype)), 
     &            MEMI(imethod(itype)), MEMI(irec(itype)), 
     &            nsize, MEMB(imerr(itype)), terr(itype), 
     &            MEMR(iworkphot), MEMR(ibin(itype)), 
     &            MEMR(ivin(itype)), contxt, ierr)

      CALL udmfre(iworkphot, 6, ierr)
      contxt = 
     &   'TBLINT: Failed to free memory for workphot array'
      IF ( ierr .NE. 0 ) GOTO 999
 
C Set the output variables from those saved internally

      nbins = nval(itype)
      nparam = npar(itype)+naddpr(itype)
      qshift = qrshft(itype)
      ibino = ibin(itype)
      ivino = ivin(itype)
      iefilo = iefil(itype)
      qmerr = terr(itype)

 999  CONTINUE

      IF ( ierr .GT. 1 ) THEN
         WRITE(wrtstr, '(a,i4)') 'TBLINT: FITSIO error ', ierr
         CALL xwrite(wrtstr, 10)
         CALL xwrite(contxt, 10)
      ENDIF

      RETURN
      END

C *******************************************************************

      SUBROUTINE intstp(npar, nump, param, values, nrec, exact, nbmin,
     &                  intpar, cntpar, rec, ierr)

      REAL values(*), param(*)
      INTEGER nump(*), nrec(*), nbmin(*), intpar(*), rec(*)
      INTEGER npar, cntpar, ierr
      LOGICAL exact(*)

C Subroutine to set up for the interpolation grid. Works out which parameters
C need to be interpolated and stores the pointers to records in the table
C file
C Parameters :
C     npar     i        i: number of interpolation parameters
C     nump     i        i: number of tabulated values for each parameter
C     param    r        i: current parameter values
C     values   r        i: tabulated parameter values
C     nrec     i        i: number of records in parameter block
C     exact    l        r: parameters that exactly match a tabulated value
C     nbmin    i        r: first entry required for each parameter
C     intpar   i        r: parameters that require interpolation
C     cntpar   i        r: number of parameters to be interpolated
C     rec      i        r: record numbers of interpolation

      INTEGER i, j, k, ioff

      CHARACTER wrtstr*255

c    TEST THAT ALL INPUT PARAMETERS ARE WITHIN BOUNDS OF TABULATED VALUES

      ioff = 0
      DO i = 1, npar
         IF ( nump(i) .GT. 1 ) THEN
            IF ( param(i) .LT. values(ioff+1) ) THEN
               WRITE (wrtstr, '(a,i2,a)') ' Parameter #', i,
     &                   ' is less than the Tabulated Value'
               CALL xwrite(wrtstr, 10)
               WRITE (wrtstr, *) i, 1, param(i), values(ioff+1)
               CALL xwrite(wrtstr, 10)
               ierr = -1
               GOTO 999
            ENDIF
            IF ( param(i) .GT. values(ioff+nump(i)) ) THEN
               WRITE (wrtstr, '(a,i2,a)') ' Parameter #', i,
     &                ' is greater than the Tabulated Value'
               CALL xwrite(wrtstr, 10)
               WRITE (wrtstr, *) i, nump(i), param(i), 
     &                           values(ioff+nump(i))
               CALL xwrite(wrtstr, 10)
               ierr = 1
               GOTO 999
            ENDIF
         ENDIF
         ioff = ioff + nump(i)
      ENDDO

c     DETERMINE THE # OF PARAMETER BLOCKS TO SKIP, AND TEST FOR EQUALITY
c         OF INPUT AND TABULATED PARAMETER VALUES

      DO i = 1, npar
         exact(i) = .FALSE.
      ENDDO
      ioff = 0
      DO i = 1, npar
         j = 2
         DO WHILE ( (param(i) .GE. values(ioff+j) )
     &             .AND. (j .LE. nump(i)) )
            j = j + 1
         ENDDO
         nbmin(i) = j - 2
         IF ( ABS(param(i)-values(ioff+j-1))
     &         .LT. 0.0001*ABS(param(i)) ) THEN
            exact(i) = .TRUE.
         ENDIF
         ioff = ioff + nump(i)
      ENDDO
      DO i = 1, npar
         IF ( nump(i) .LE. 1 ) exact(i) = .TRUE.
      ENDDO

c     DETERMINE WHICH PARAMETERS NEED TO BE INTERPOLATED

      cntpar = 0
      DO i = npar, 1, -1
         IF (.NOT. exact(i)) THEN
            cntpar = cntpar + 1
            intpar(cntpar) = i
         ENDIF
      ENDDO

c     CALCULATE THE NECESSARY RECORD NUMBERS AND STORE IN REC ARRAY

      rec(1) = 1
      DO i = 1, npar
         rec(1) = rec(1) + nbmin(i)*nrec(i)
      ENDDO

      DO i = 2, 2**cntpar
         j = 0
         k = i - 1
         DO WHILE (mod(k,2).EQ.0)
            k = k/2
            j = j + 1
         ENDDO
         rec(i) = rec(i-2**j) + nrec(intpar(1+j))
      ENDDO

 999  CONTINUE

      RETURN
      END

C *******************************************************************

      SUBROUTINE grdint(ilun, cntpar, npar, nval, naddpr, param,
     &                  intpar, method, rec, n1, qmerr, terr,
     &                  workphot, bin, vin, contxt, ierr)

      INTEGER n1
      REAL param(*), bin(*), vin(*), workphot(n1,*)
      INTEGER intpar(*), method(*), rec(*)
      INTEGER ilun, cntpar, npar, nval, naddpr, ierr
      LOGICAL qmerr(*), terr  
      CHARACTER contxt*(*)

C Do the actual interpolation on the grid
C Parameters :
C       ilun     i        i: I/O unit for table file
C       cntpar   i        i: number of parameters to be interpolated
C       npar     i        i: total number of interpolation parameters
C       nval     i        i: number of energy ranges
C       naddpr   i        i: number of addition parameters
C       param    r        i: current parameter values
C       intpar   i        i: array of parameters to be interpolated
C       method   i        i: interpolation methods
C       rec      i        i: records for parameter interpolation
C       n1       i        i: size of first axis of workphot (npar+nval)
c       qmerr    l        i: array: true if model error 
c       terr     l        i: true if one of qmerr element true 
C       workphot r        w: work array to accumulate model spectra
C       bin      r        r: output array of model spectrum
C       contxt   c        r: error description
C       ierr     i        r: error flag
C
C Note: if model has associated errors (terr=true)
C       the variance of each spectrum is appended to the end of first dim   
c       ie,  variance of i spectrum -> workphot(npar+nval:npar+2*nval, i)
c       the variance is interpolated in the same way as spectrum.
c       When it comes to interpolation, this might be a better approach
c       to handle the error than the usual quadratic thing. 
c 

      INTEGER i, j, ipoint, jpoint, ipar, jcount

c Handle case where the input parameter values exactly match those
c stored in a record in the file

      IF (cntpar .EQ. 0) THEN
         CALL rdmdsp(ilun, rec(1), npar, nval, naddpr, param(npar+1), 
     &               workphot(1,1), qmerr, bin, contxt, ierr)
         IF ( ierr .NE. 0 ) GOTO 999
         GOTO 10
      ENDIF

c     READ IN AND INTERPOLATE FIRST TWO RECORDS - A SPECIAL CASE

      DO i = 1, 2
         CALL rdmdsp(ilun, rec(i), npar, nval, naddpr, param(npar+1), 
     &              workphot(1,i), qmerr, bin, contxt, ierr)
         IF ( ierr .NE. 0 ) GOTO 999
      ENDDO

      CALL interp(n1-npar, npar, intpar(1), param(intpar(1)),
     &            workphot(1,1), workphot(1,2), method(intpar(1)))

c     NOW READ IN AND INTERPOLATE REST OF GRID
c     pointers track current position in workspace array and list
c     of records to be read. ipar is the current parameter. This
c     somewhat complex algorithm goes as cntpar+1 in work space
c     required as opposed to 2**cntpar for the obvious binary
c     algorithm.

      ipoint = 1
      jpoint = 2

      DO i = 2, cntpar
         DO j = 1, 2**(i-2)
            CALL rdmdsp(ilun, rec(jpoint+1), npar, nval, naddpr, 
     &                  param(npar+1), workphot(1,ipoint+1), 
     &                  qmerr,  bin, contxt, ierr)
            IF ( ierr .NE. 0 ) GOTO 999
            CALL rdmdsp(ilun, rec(jpoint+2), npar, nval, naddpr, 
     &                  param(npar+1), workphot(1,ipoint+2), 
     &                  qmerr, bin, contxt, ierr)
            IF ( ierr .NE. 0 ) GOTO 999
            ipar = 1
            ipoint = ipoint + 2
            jpoint = jpoint + 2
            jcount = jpoint
            DO WHILE(mod(jcount,2).EQ.0)
               CALL interp(n1-npar, npar, intpar(ipar),
     &                     param(intpar(ipar)), workphot(1,ipoint-1),
     &                     workphot(1,ipoint), method(intpar(ipar)))
               ipoint = ipoint - 1
               ipar = ipar + 1
               jcount = jcount / 2
            ENDDO
         ENDDO
      ENDDO

 10   CONTINUE
      DO j = 1, nval
         bin(j) = workphot(j+npar, 1)
      ENDDO

      IF(terr)THEN
        DO j = 1, nval
           vin(j) = workphot(j+npar+nval, 1)
        ENDDO
      ENDIF

 999  CONTINUE

      RETURN
      END


c *******************************************************************

      SUBROUTINE interp(nval, npar, ipar, parint, val1, val2, method)

c        PERFORMS  1-D INTERPOLATION
c        writes results into val1 array
c


      INTEGER j, nval, method, npar, ipar
      REAL val1(*), val2(*), par1, par2, parint
      REAL fact

c  set the parameter values for the two arrays

      par1 = val1(ipar)
      par2 = val2(ipar)


      IF ( par1 .NE. par2 ) THEN

c  linear interpolation

         IF (method.EQ.0) THEN
            fact = (parint-par1)/(par2-par1)

c  logarithmic interpolation

         ELSEIF (method.EQ.1) THEN
            fact = (log10(parint)-log10(par1))
     &             /(log10(par2)-log10(par1))
         ENDIF

         DO j = 1+npar, nval+npar
            val1(j) = val1(j) + (val2(j)-val1(j))*fact
         ENDDO

      ENDIF

      val1(ipar) = parint

      RETURN
      END

c *******************************************************************

      SUBROUTINE rdmdsp(ilun, rec, npar, nval, naddpr, addpar, phot, 
     &                  qmerr, work, contxt, ierr)

      INTEGER ilun, rec, npar, nval, naddpr, ierr
      REAL addpar(*), phot(*), work(*)
      LOGICAL qmerr(0:*)
      CHARACTER*(*) contxt

c        Reads the model spectrum for a given row and adds in the
c        additional parameter contributions, including variance if any

      INTEGER i, k
      LOGICAL qanyf
      CHARACTER*255 wrtstr

c Read the parameter values

      CALL ftgcve(ilun, 1, rec, 1, npar, 1.e-31, phot, qanyf, ierr)
      WRITE(contxt,'(a,i6,a)') 'RDMDSP: Failed to read record ', rec, 
     &                         ' column 1 in extension 4'
      IF ( ierr .NE. 0 ) GOTO 999

      WRITE(wrtstr,'(a,i5,20(1x,1pe9.2))') 'Interp ', rec,
     &                               (phot(i),i=1,MIN(20,npar-naddpr))
      CALL xwrite(wrtstr, 25)

c Read the interpolation parameters spectrum

      CALL ftgcve(ilun, 2, rec, 1, nval, 1.e-31, phot(1+npar), qanyf, 
     &            ierr)
      WRITE(contxt,'(a,i6,a)') 'RDMDSP: Failed to read record ', rec, 
     &                         ' column 2 in extension 4'
      IF ( ierr .NE. 0 ) GOTO 999
c
c read model error if any
c
      IF(qmerr(0))THEN 
          CALL ftgcve(ilun, 2, rec, nval+1, nval, 1.e-31, 
     &                phot(1+nval+npar), qanyf, ierr)
          WRITE(contxt,'(a,i6,a)') 
     &     'RDMDSP: Failed to read error in record ', rec, 
     &                         ' column 2 in extension 4'
          IF ( ierr .NE. 0 ) GOTO 999
c
c    convert to variance
c
         DO k = nval+npar+1, 2*nval+npar
            phot(k) = phot(k)**2 
         ENDDO
      END IF

c Read the additional parameters spectra and sum them in

      DO i = 1, naddpr

         CALL ftgcve(ilun, i+2, rec, 1, nval, 1.e-31, work, qanyf, 
     &               ierr)
         WRITE(contxt,'(a,i6,a,i6,a)') 'RDMDSP: Failed to read record ',
     &           rec, ' column ', i+2, ' in extension 4'
         IF ( ierr .NE. 0 ) GOTO 999

         DO k = 1, nval
            phot(k+npar) = phot(k+npar) + addpar(i)*work(k)
         ENDDO

         IF(qmerr(i))THEN
             CALL ftgcve(ilun, i+2, rec, nval+1, nval, 1.e-31, 
     &                   work, qanyf, ierr)
             WRITE(contxt,'(a,i6,a,i6,a)') 
     &           'RDMDSP: Failed to read error in record ',
     &           rec, ' column ', i+2, ' in extension 4'
             IF ( ierr .NE. 0 ) GOTO 999
c
c accumulate variance
c 
             DO k = 1, nval 
               phot(npar+nval+k) = 
     &           phot(npar+nval+k)+(addpar(i)*work(k))**2
             ENDDO
         ENDIF

      ENDDO 

 999  CONTINUE
      RETURN
      END

c *******************************************************************
c find out which columns have associated errors in SPECTRA extension

      SUBROUTINE  merr_info (unit, ne, nadd, qmerr, terr, ierr)
      INTEGER unit, ne, nadd, ierr
      LOGICAL qmerr(*), terr

c Arguments :
c    unit       i        i: i/o channel
c    ne         i        i: number of energy bins
c    nadd       i        i: number of additional parameters
c    qmerr      l        r: set to true if model errors for this additional
c                           parameter
c    terr       l        r: set to true if any model errors

      INTEGER i, j, ndim
      CHARACTER  key*8, value*20, comment*80, contxt*72

      ierr = 0
      terr = .false.
      DO i = 1, nadd+1
         qmerr(i) = .false.
      ENDDO

      DO i = 2, nadd+2 
         WRITE (key, '(a5,i1)') 'TFORM', i
         IF ( i .GT. 9  ) WRITE (key, '(a5,i2)') 'TFORM', i
         IF ( i .GT. 99 ) WRITE (key, '(a5,i3)') 'TFORM', i

         CALL FTGKEY (unit, key, value, comment, ierr)
         contxt = 'Failed to read '//key
         IF ( ierr .NE. 0 ) GOTO 999

         j = index(value(2:),'E')
         IF ( j .EQ. 2 ) READ (value(2:j), '(i1)', iostat=ierr) ndim
         IF ( j .EQ. 3 ) READ (value(2:j), '(i2)', iostat=ierr) ndim
         IF ( j .EQ. 4 ) READ (value(2:j), '(i3)', iostat=ierr) ndim
         IF ( j .EQ. 5 ) READ (value(2:j), '(i4)', iostat=ierr) ndim
         IF ( j .EQ. 6 ) READ (value(2:j), '(i5)', iostat=ierr) ndim
         IF ( j .EQ. 7 ) READ (value(2:j), '(i6)', iostat=ierr) ndim
         IF ( j .EQ. 8 ) READ (value(2:j), '(i7)', iostat=ierr) ndim
         contxt = 'Failed to extract size from '//value
         IF ( ierr .NE. 0 ) GOTO 999

         IF ( ndim .EQ. 2*ne ) THEN
            qmerr(i-1) = .true.
            terr = .true.
         ENDIF
      ENDDO 

 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 5)
         WRITE(contxt, '(a, i3)') 'Error in MERR_INFO = ', ierr
         CALL xwrite(contxt, 5)
      ENDIF      

      RETURN
      END 
