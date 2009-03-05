
      SUBROUTINE xsetbl(ear, ne, param, filenm, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(*), photar(ne), photer(*)
      CHARACTER*(*) filenm

c  XSPEC model subroutine. Multiplicative table model which takes 
c  EXP(-(table input)). Mainly useful for complicated variable
c  abundance absorption where the abundances can then be treated
c  as additive parameters

c  6/21/95  kaa

      INTEGER ie
      LOGICAL moderr 

c Get the photar as if this were a standard multiplicative model but with
c the value for energies outside those tabulated set to 0.

      CALL mtbint(ear, ne, param, filenm, ifl, photar, photer, 0.0, 
     &            moderr)

c Now do the exponentiation

      DO ie = 1, ne
         IF ( photar(ie) .GT. 100. ) THEN
            photar(ie) = 0.
         ELSE
            photar(ie) = EXP(-photar(ie))
         ENDIF
      ENDDO
c
c variance
c
      IF(moderr)THEN
         DO ie = 1, ne
             photer(ie) = photar(ie)**2 * photer(ie)
         ENDDO
      ENDIF 

      RETURN
      END


