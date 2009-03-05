
c  Include file for XSPEC internal parameters


c      MAXMOD          Number of model components
c      MAXMPR          Number of model PARAMETERs
c      MAXFPR          Number of fit PARAMETERs
c      MAXVPR          Number of variable PARAMETERs
c      MAXFIL          Number of PHA files
c      MAXFLT          Maximum number of filters from XFLTxxxx keywords


      INTEGER MAXMOD, MAXMPR, MAXFPR, MAXVPR, MAXFIL, MAXFLT

      PARAMETER (MAXMOD=250,MAXMPR=1000,MAXFPR=1000,MAXVPR=100)
      PARAMETER (MAXFIL=500,MAXFLT=10)

      CHARACTER*12 VERSION
      PARAMETER (VERSION='11.3.2d')


c Common block for dynamic memory using udmget

      LOGICAL          MEMB(1)
      INTEGER*2        MEMS(1)
      INTEGER*4        MEMI(1)
      INTEGER*4        MEML(1)
      REAL             MEMR(1)
      DOUBLE PRECISION MEMD(1)
      COMPLEX          MEMX(1)
      CHARACTER*1      MEMC(1)
      EQUIVALENCE (MEMB, MEMS, MEMI, MEML, MEMR, MEMD, MEMX, MEMC)
      COMMON /MEM/ MEMD

