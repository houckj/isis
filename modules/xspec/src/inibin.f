
      SUBROUTINE inibin(nin, ein, nout, eout, start, end, fstart, fend,
     &                  z)

      INTEGER nin, nout
      INTEGER start(nout), end(nout)
      REAL    ein(0:nin), eout(0:nout), fstart(nout), fend(nout), z

c		rashafer 1 april 1987
c	subroutine to initialize the variable used by EREBIN to rebin
c	arrays, assuming that the samples are flat topped histograms
c	from the lower limit to the  upper limit of the bin (e.g. photon
c	arrays)
c 		version 2 6 April 1987 --- include redshift parameter
c
c Arguments :
c     nin     i                 I: number of unbinned energy ranges
c     ein     r                 I: The unbinned boundaries
c     nout    i                 I: The number of binned ranges
c     eout    r                 I: The binned boundaries
c     start   i                 R: The first unbinned point
c     end     i                 R: The last unbinned point (if <=start(i)
c                                    then there is only one point used)
c     fstart  r                 R: The weighting for the first(or only) point
c     fend    r                 R: The weighting for the last point.  If
c                                  end(i)>start(i)+1, then all points between
c                                  have a weighting of 1.
c     z       r                 I: Redshift:  The input model is redshifted
c                                  by the factor Z, i.e. all input energies
c                                  are treated as coming from ranges
c                                  divided by (1+z)

      INTEGER in, out, i
      REAL    zf
      CHARACTER wrtstr*255

c Initialize the output arrays to give no contributions

      DO i = 1, nout
         start(i) = 1
         end(i) = 0
         fstart(i) = 0.
         fend(i) = 0.
      ENDDO

c check whether the redshift is reasonable - if not write out a warning
c and return

      IF ( z .NE. -1. ) THEN
         zf = 1./(1.+z)
      ELSE
         CALL xwrite(
     &     ' Redshift is -1 which will lead to a divide by zero', 2)
         RETURN
      ENDIF

c find the first unbinned boundary below the first binned boundary

      out = 0
      DO WHILE ((out.LT.nout) .AND. (eout(out)*1.000001.LT.zf*ein(0)))
         out = out + 1
      ENDDO
      out = out + 1

      IF ( out .GT. 1 ) THEN

         CALL xwrite (
     & ' No model information is available for observed frame energies',
     &                15)
         WRITE (wrtstr, '(a,1pg14.4,a)')
     & ' below ', eout(out-1), ' so the model is assumed to be zero'
         CALL xwrite(wrtstr, 15)

c  No contribution


c  N.B.  In the current implementation, all the output bins
c        that have any part below the first input bin are set to
c        zero flux.  However, if a bin is partially above the
c        last input bin, then it is allowed a partial accumulation
c        (It just was easier to do that one at the time).

      ENDIF

      in = 1
      DO WHILE ((in.LE.nin) .AND. (zf*ein(in).LT.eout(out-1)))
         in = in + 1
      ENDDO

      IF ( in .GT. nin ) THEN
         CALL xwrite(
     & ' There is no model information available in the requested ',
     &               5)
         WRITE(wrtstr, '(a,1pg14.4,a,1pg14.4)')
     & ' observed frame energy range of ', eout(0), ' to ', eout(nout)
         CALL xwrite(wrtstr, 5)
         WRITE(wrtstr, '(a,1pg14.4,a,1pg14.4)')
     & ' The model rest frame energy range is ', ein(0), ' to ', 
     &                                           ein(nin)
         CALL xwrite(wrtstr, 5)
         WRITE(wrtstr, '(a,1pg14.4)') ' and the redshift is ', z
         CALL xwrite(wrtstr, 5)
         RETURN
      ENDIF

      DO WHILE ( out .LE. nout )

c the output bin is completely within a single input range

         IF ( eout(out) .LE. zf*ein(in) ) THEN

            end(out) = 0
            start(out) = in
            fstart(out) = (eout(out)-eout(out-1))
     &                    /(zf*(ein(in)-ein(in-1)))

c Exactly finished the bin

            IF ( eout(out) .EQ. zf*ein(in) ) in = in + 1

            out = out + 1

c Now do the case where the output bin overlaps input ranges

         ELSE

            start(out) = in
            fstart(out) = (zf*ein(in)-eout(out-1))
     &                    /(zf*(ein(in)-ein(in-1)))

            DO WHILE ((in.LE.nin) .AND. (eout(out).GT.zf*ein(in)))
               in = in + 1
            ENDDO

            IF ( in .GT. nin ) THEN

               end(out) = nin
               IF (end(out).NE.start(out)) THEN
                  fend(out) = 1.
               ELSE
                  end(out) = 0
               ENDIF
               out = out + 1

            ELSE

               end(out) = in
               fend(out) = (eout(out)-zf*ein(in-1))
     &                     /(zf*(ein(in)-ein(in-1)))
               IF (zf*ein(in).EQ.eout(out)) in = in + 1
               out = out + 1

            ENDIF

         ENDIF

c Now check for case that we have run out of input bins

         IF ( in .GT. nin ) THEN

            CALL xwrite (
     & ' No model information is available for observed frame energies',
     &                   15)
            WRITE (wrtstr, '(a,1pg14.4,a)')
     & ' above ', eout(out-2), ' so the model is assumed to be zero'
            CALL xwrite(wrtstr, 15)

c Signal the dowhile loop that we are done

            out = nout + 1

         ENDIF

      ENDDO

      RETURN
      END


