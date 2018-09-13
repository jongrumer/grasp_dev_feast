************************************************************************
*                                                                      *
      SUBROUTINE HFS (host)
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the  hyperfine structure parameters.               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, GRACAH1, ISPAR, ITJPO, *
*                        ONEPARTICLEJJ.                                *
*               [HFS92]: MATELT, RINT, RINTHF.                         *
*                                                                      *
*   Written by Per Jonsson and Farid A. Parpia                         *
*                                                                      *
*                                         Last revision: 10 Nov 1995   *
*                                                                      *
*   Modified by G. Gaigalas by includeing the new spin-angular         *
*   libraries.                                      Vilnius May 2012   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
c
cbieron include 'parameters.def'
c
      include 'parameters.def'
c
c      PARAMETER (NNNW = 120)
c
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      INTEGER FFMIN,FFMAX,FF
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION TSHELL(NNNW),RINTME(2,NNNW,NNNW),AMELT(2,NNNW,NNNW)
*
      POINTER (PNTHFC,HFC(5,*))
*P    POINTER (PNTHFCtmp,HFCtmp(5,*))
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF3/EMPAM,RBCM
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /DEF11/FMTOAU,B1
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /NSMDAT/HFSI,HFSD,HFSQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /OPT6/NTC(10)
     :      /ORB2/NCF,NW,PNTRIQ
*
      include 'mpif.h'
      integer  myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr 

      integer lenhost
c     character*80  host
c     CHARACTER host*(MPI_MAX_PROCESSOR_NAME)
      character*(*) host

      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond
      character msg*80

*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
*   Allocate storage for local arrays
*
cbdebug
c      print *, ' hfs 1'

      CALL ALLOC (PNTHFC,5*NVEC*NVEC,8)
*P    CALL ALLOC (PNTHFCtmp,5*NVEC*NVEC,8)
*
*   Initialise
*
      DO 2 I = 1,5
        DO 1 J = 1,NVEC*NVEC
           HFC(I,J) = 0.0D 00
*P         HFCtmp(I,J) = 0.0D 00
    1   CONTINUE
    2 CONTINUE
*
*   Calculate and save the radial integrals and angular
*   matrix elements for the two multipolarities
*
      DO 5 KT = 1,2
        DO 4 I = 1,NW
          DO 3 J = 1,NW
            IF (KT .EQ. 1) THEN
              RINTME(KT,I,J) = RINTHF (I,J,-2)
            ELSE
              RINTME(KT,I,J) = RINT (I,J,-3)
            ENDIF
            CALL MATELT (I,KT,J,APART)
            AMELT(KT,I,J) = APART
    3     CONTINUE
    4   CONTINUE
    5 CONTINUE
*
*   Set the parity of the one-body operators
*
      IPT = 1
*
*   Sweep through the Hamiltonian matrix to determine the
*   diagonal and off-diagonal hyperfine constants
*
cb MPI
c      DO 11 IC = 1,NCF
      DO 11 IC = myid+1, NCF, nprocs
*
*   Output IC on the screen to show how far the calculation has preceede
*
         IF (IC .LE. 3*nprocs .OR. IC .GE. NCF-3*nprocs
     &       .OR. MOD (IC-1,10*nprocs) .EQ. myid) THEN
cb      if (mod(IC-1,10*nprocs).eq.myid) then
        CALL CONVRT (IC,CNUM,LCNUM)
        PRINT *, 'myid=',myid, ' Col. '//CNUM(1:LCNUM)//' '
       endif
*
        DO 10 IR = 1,NCF
*
*   If LFORDR is .TRUE., a `first order' calculation is indicated;
*   only the CSFs with serial numbers exceeding IC are treated specially
*   only diagonal elements are evaluated for the `first order' CSFs
*
          IF (  LFORDR .AND.
     :         (IC .GT. ICCUT) .AND.
     :         (IC .NE. IR) ) GOTO 10
*
          ISPARC = ISPAR (IC)
          ITJPOC = ITJPO (IC)
          ITJPOR = ITJPO (IR)
*
*   Loop over the multipolarities
*
          DO 9 KT = 1,2
*
*   Initialise the accumulator
*
            ELEMNT = 0.0D 00
*
*   Consider 3 cases
*               (k)
*   (1) < J || T   || J >    , k = 1,2  and IR >= IC
*
*               (k)
*   (2) < J || T   || J-1 >  , k = 1,2
*
*               (k)
*   (3) < J || T   || J-2 >  , k = 2
*
            IF (((ITJPOC .EQ. ITJPOR) .AND. (IR. GE. IC)) .OR.
     :          ((ITJPOC - ITJPOR) .EQ. 2) .OR.
     :         (((ITJPOC - ITJPOR) .EQ. 4) .AND. (KT .EQ. 2))) THEN
*
              CALL ONEPARTICLEJJ(KT,IPT,IC,IR,IA,IB,TSHELL)
CGG              CALL TNSRJJ (KT,IPT,IC,IR,IA,IB,TSHELL)
*
*   Accumulate the contribution from the one-body operators;
*
              IF (IA .NE. 0) THEN
                IF (IA .EQ. IB) THEN
                  DO 6 IA = 1,NW
                    IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                      ELEMNT =   ELEMNT
     :                         + AMELT(KT,IA,IA)
     :                         * RINTME(KT,IA,IA)
     :                         * TSHELL(IA)
                    ENDIF
    6             CONTINUE
                ELSE
                  IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                    ELEMNT =   ELEMNT
     :                       + AMELT(KT,IA,IB)
     :                       * RINTME(KT,IA,IB)
     :                       * TSHELL(1)
                  ENDIF
                ENDIF
              ENDIF
*
*   Multiply with the configuration expansion coefficients and add the
*   contributions from the matrix elements to obtain total contributions
*
              DO 8 K = 1,NVEC
                DO 7 KK = 1,NVEC
                  LOC1 = (K-1)*NCF
                  LOC2 = (KK-1)*NCF
                  IF ((ITJPOC .EQ. ITJPOR) .AND. (IR .NE. IC)) THEN
                     CONTR =   ELEMNT
     :                       * (  EVEC(IC+LOC1) * EVEC(IR+LOC2)
     :                          + EVEC(IR+LOC1) * EVEC(IC+LOC2) )
                  ELSE
                    CONTR =    ELEMNT
     :                       *    EVEC(IC+LOC1) * EVEC(IR+LOC2)
                  ENDIF
*
*   Magnetic dipole
*
*P                IF (KT .EQ. 1) THEN
*P                  IF ((ITJPOC-ITJPOR) .EQ. 0) THEN
*P                   HFCtmp(1,NVEC*(K-1)+KK) =   HFCtmp(1,NVEC*(K-1)+KK)
*P   :                                       + CONTR
*P                  ELSEIF ((ITJPOC-ITJPOR) .EQ. 2) THEN
*P                   HFCtmp(2,NVEC*(K-1)+KK) =   HFCtmp(2,NVEC*(K-1)+KK)
*P   :                                       + CONTR
*P                  ENDIF
                  IF (KT .EQ. 1) THEN
                    IF ((ITJPOC-ITJPOR) .EQ. 0) THEN
                     HFC(1,NVEC*(K-1)+KK) =   HFC(1,NVEC*(K-1)+KK)
     :                                       + CONTR
                    ELSEIF ((ITJPOC-ITJPOR) .EQ. 2) THEN
                     HFC(2,NVEC*(K-1)+KK) =   HFC(2,NVEC*(K-1)+KK)
     :                                       + CONTR
                    ENDIF
*
*   Electric quadrupole
*
                  ELSEIF (KT .EQ. 2) THEN
                    IF ((ITJPOC-ITJPOR) .EQ. 0) THEN
*P                   HFCtmp(3,NVEC*(K-1)+KK) =   HFCtmp(3,NVEC*(K-1)+KK)
                     HFC(3,NVEC*(K-1)+KK) =   HFC(3,NVEC*(K-1)+KK)
     :                                       + CONTR
                    ELSEIF ((ITJPOC-ITJPOR) .EQ. 2) THEN
*P                   HFCtmp(4,NVEC*(K-1)+KK) =   HFCtmp(4,NVEC*(K-1)+KK)
                     HFC(4,NVEC*(K-1)+KK) =   HFC(4,NVEC*(K-1)+KK)
     :                                       + CONTR
                    ELSEIF ((ITJPOC-ITJPOR) .EQ. 4) THEN
*P                   HFCtmp(5,NVEC*(K-1)+KK) =   HFCtmp(5,NVEC*(K-1)+KK)
                     HFC(5,NVEC*(K-1)+KK) =   HFC(5,NVEC*(K-1)+KK)
     :                                       + CONTR
                    ENDIF
                  ENDIF
    7           CONTINUE
    8         CONTINUE
*
            ENDIF
    9     CONTINUE
*
   10   CONTINUE
   11 CONTINUE

*P            CALL MPI_ALLREDUCE (HFCtmp,HFC,5*NVEC*NVEC,
*P   &             MPI_DOUBLE_PRECISION,
*P   &             MPI_SUM, MPI_COMM_WORLD, ierr)
              CALL MPI_ALLREDUCE (MPI_IN_PLACE,HFC,5*NVEC*NVEC,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_SUM, MPI_COMM_WORLD, ierr)

cb 
cb from here to END only node0 works
      if (myid .eq. 0) then !node0
*
c     call MPI_Get_processor_name (host, lenhost, ierr)
      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      msg = '  start: ' // 
     &      '  Date: ' // chdate //
     &      '  Time: ' // chtime 
      msgLength = len_trim (msg)
      print *, 'Only node0 works:'
      print *, msg(1:msgLength)      

*   Printouts
*
      WRITE (24,301) CUTOFF
*
*   These are the conversion factors to obtain the hyperfine
*   constants in MHz
*
      AUMHZ = AUCM*CCMS*1.0D-06
      BARNAU = 1.0D-24/RBCM**2
      DNMAU = B1/(2.0D 00*CVAC*EMPAM)
*
      GFAC = AUMHZ*DNMAU*HFSD/HFSI
      HFAC = AUMHZ*2.0D 00*HFSQ*BARNAU
*
*   Output the hyperfine interaction constants
*
      WRITE (24,302)
*
      DO 13 I = 1,NVEC
        DO 12 II = 1,NVEC
*
          JJ = IATJPO(I)
          JJII = IATJPO(II)
*
          IF ((JJ.EQ.JJII .AND. JJII.GT.1) .OR. (JJ.GT.JJII)) THEN
*
            FJ = 0.5D 00*DBLE (JJ-1)
*
            AFA1 =   GFAC
     :             * SQRT (
     :                       1.0D 00
     :                     / (FJ*(FJ+1.0D 00))
     :                    )
            IF (JJ.EQ.2) THEN
               AFA2 = 0.D0 00
            ELSE 
               AFA2 =   GFAC
     :             * SQRT (
     :                       1.0D 00
     :                     / (FJ*(2.0D 00*FJ-1.0D 00))
     :                    )
            ENDIF
            BFA1 =   HFAC
     :             * SQRT (
     :                       (FJ*(2.0D 00*FJ-1.0D 00))
     :                     / (   (FJ+1.0D 00)
     :                         * (2.0D 00*FJ+3.0D 00) )
     :                    )
            IF (JJ.EQ.2) THEN 
               BFA2 = 0.0D 00
            ELSE
               BFA2 =   0.25D 00
     :             * HFAC
     :             * SQRT (
     :                       (FJ*(FJ-1.0D 00))
     :                     / (   (FJ+1.0D 00)
     :                         * (2.0D 00*FJ-1.0D 00) )
     :                    )
            ENDIF
            IF (JJ.EQ.4) THEN
               BFA3 = 0.0D 00
            ELSE
               BFA3 =   0.125D 00
     :             * HFAC
     :             * SQRT (
     :                       (FJ*(FJ-1.0D0 00)
     :                     * (2.0D 00*FJ-1.0D 00) )
     :                     / (2.0D 00*FJ-3.0D 00)
     :                     )
            ENDIF
*
*   Diagonal (J,J) A and B factors
*
            IF ((JJ-JJII) .EQ. 0) THEN
              IF (I .LE. II)
     :          WRITE (24,303) IVEC(I),LABJ(JJ),
     :                         LABP((IASPAR(I)+3)/2),
     :                         IVEC(II),LABJ(JJII),
     :                         LABP((IASPAR(II)+3)/2),
     :                         AFA1*HFC(1,NVEC*(I-1)+II),
     :                         BFA1*HFC(3,NVEC*(I-1)+II)
      WRITE (24,309) '                QED corrected value',
     :AFA1*HFC(1,NVEC*(I-1)+II)*2.0023193/2.0D 00
*
*   Off diagonal (J,J-1) A and B factors
*
            ELSEIF ((JJ-JJII) .EQ. 2) THEN
              WRITE (24,303) IVEC(I),LABJ(JJ),
     :                       LABP((IASPAR(I)+3)/2),
     :                       IVEC(II),LABJ(JJII),
     :                       LABP((IASPAR(II)+3)/2),
     :                       AFA2*HFC(2,NVEC*(I-1)+II),
     :                       BFA2*HFC(4,NVEC*(I-1)+II)
      WRITE (24,309) '                QED corrected value',
     :AFA2*HFC(2,NVEC*(I-1)+II)*2.0023193/2.0D 00
*
*   Off diagonal (J,J-2) B factor
*
            ELSEIF ((JJ-JJII) .EQ. 4) THEN
              WRITE (24,303) IVEC(I),LABJ(JJ),
     :                       LABP((IASPAR(I)+3)/2),
     :                       IVEC(II),LABJ(JJII),
     :                       LABP((IASPAR(II)+3)/2),0.0D 00,
     :                       BFA3*HFC(5,NVEC*(I-1)+II)
            ENDIF
*
          ENDIF
   12   CONTINUE
   13 CONTINUE
*
*   These are the factors needed to obtain the F-dependent hyperfine
*   matrix elements in Hartrees
*
      TILDE1 =   SQRT ( (HFSI+1.0D 00)
     :                  /HFSI )
     :         * HFSD
     :         * 0.5D 00
     :         * B1/(EMPAM*CVAC)
      IF (HFSI .GT. 0.6D 00) THEN
        TILDE2 = SQRT ( (HFSI+1.0D 00)
     :                 *(2.0D 00*HFSI+3.D0)
     :                 /((2.0D 00*HFSI-1.0D 00)*HFSI) )
        TILDE2 = TILDE2*HFSQ*0.5D 00*BARNAU
      ELSE
        TILDE2 = 0.0D 00
      ENDIF
*
*     II = 2*I
*
      II = NINT (2.0D 00*HFSI)
*
*   Calculate the F-dependent matrix elemnts
*
*   Loop over the states
*
      DO 16 JB = 1,NVEC
        DO 15 JA = 1,NVEC
          JJB = IATJPO(JB)-1
          JJA = IATJPO(JA)-1
          IF ((JJA .EQ. JJB .AND. JJA .GT. 0) .OR.
     :        (JJB .GT. JJA)) THEN
*
*   Determine the possible F quantum numbers for the matrix element
*
            FFMIN = MAX (ABS (JJA-II),ABS (JJB-II))
            FFMAX = MIN (JJA+II,JJB+II)
*
*   Loop over the possible F quantum numbers
*
            IFLAG = 0
            DO 14 FF = FFMIN,FFMAX,2
*
*   Phase factor
*
              IF (MOD ((II+JJA-FF)/2,2) .EQ. 1) THEN
                FACTOR1 = -1.0D 00
              ELSE
                FACTOR1 =  1.0D 00
              ENDIF
              FACTOR2 =   FACTOR1
     :                  * SQRT (
     :                            ( DBLE (JJB)
     :                             +1.0D 00)
     :                          * ( DBLE (II)
     :                             +1.0D 00)
     :                         )
*
*   Determine the Racah W coefficients.
*
CGG              CALL DRACAH (II,JJA,II,JJB,FF,2,RAC1)
              CALL GRACAH1 (II,JJA,II,JJB,FF,2,RAC1)
CGG              CALL DRACAH (II,JJA,II,JJB,FF,4,RAC2)
              CALL GRACAH1 (II,JJA,II,JJB,FF,4,RAC2)
*
*   Obtain and output matrix elements for J,J
*
              IF ((JJA-JJB) .EQ. 0) THEN
                HFSELT1 = FACTOR2*RAC1*HFC(1,NVEC*(JB-1)+JA)*TILDE1
                HFSELT2 = FACTOR2*RAC2*HFC(3,NVEC*(JB-1)+JA)*TILDE2
                IF (ABS (HFSELT1+HFSELT2) .GT.
     :              CUTOFF*10.0D-05) THEN
                  IF (IFLAG .EQ. 0) WRITE (24,304)
                  IFLAG = 1
                  WRITE (24,305) IVEC(JB),LABJ(JJB+1),
     :                           LABP((IASPAR(JB)+3)/2),
     :                           IVEC(JA),LABJ(JJA+1),
     :                           LABP((IASPAR(JA)+3)/2),
     :                           LABJ(FF+1),
     :                           HFSELT1+HFSELT2
                ENDIF
*
*   Obtain and output matrix elements for J,J-1
*
              ELSEIF (ABS (JJA-JJB) .EQ. 2) THEN
                HFSELT1 = FACTOR2*RAC1*HFC(2,NVEC*(JB-1)+JA)*TILDE1
                HFSELT2 = FACTOR2*RAC2*HFC(4,NVEC*(JB-1)+JA)*TILDE2
                IF (ABS (HFSELT1+HFSELT2) .GT.
     :              CUTOFF*10.0D-05) THEN
                  IF (IFLAG .EQ. 0) WRITE (24,304)
                  IFLAG = 1
                  WRITE (24,305) IVEC(JB),LABJ(JJB+1),
     :                           LABP((IASPAR(JB)+3)/2),
     :                           IVEC(JA),LABJ(JJA+1),
     :                           LABP((IASPAR(JA)+3)/2),
     :                           LABJ(FF+1),
     :                           HFSELT1+HFSELT2
                ENDIF
*
*   Obtain and output matrix elements for J,J-2
*
              ELSEIF (ABS (JJA-JJB) .EQ. 4) THEN
                HFSELT1 = 0.0D 00
                HFSELT2 = FACTOR2*RAC2*HFC(5,NVEC*(JB-1)+JA)*TILDE2
                IF (ABS (HFSELT1+HFSELT2) .GT.
     :              CUTOFF*10.0D-05) THEN
                  IF (IFLAG .EQ. 0) WRITE (24,304)
                  IFLAG = 1
                  WRITE (24,305) IVEC(JB),LABJ(JJB+1),
     :                           LABP((IASPAR(JB)+3)/2),
     :                           IVEC(JA),LABJ(JJA+1),
     :                           LABP((IASPAR(JA)+3)/2),
     :                           LABJ(FF+1),
     :                           HFSELT1+ HFSELT2
                ENDIF
              ELSE
                HFSELT1 = 0.0D 00
                HFSELT2 = 0.0D 00
              ENDIF
  14        CONTINUE
          ENDIF
  15    CONTINUE
  16  CONTINUE
*
      CALL DALLOC (PNTHFC)
*P    CALL DALLOC (PNTHFCtmp)
cb
      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      msg = '  stop : ' // 
     &      '  Date: ' // chdate //
     &      '  Time: ' // chtime 
      msgLength = len_trim (msg)
      print *, msg(1:msgLength)      

      endif ! node0

      RETURN
*
  301 FORMAT (//' CUTOFF set to ',1PD22.15)
  302 FORMAT (//' Interaction constants:'//
     :' Level1  J Parity  Level2  J Parity',8X,'A (MHz)',13X,'B (MHz)'/)
  303 FORMAT (1X,1I3,5X,2A4,2X,1I3,5X,2A4,1P,2D20.10)
  304 FORMAT (//' Matrix elements:'//
     :' Level1  J Parity  Level2  J Parity    F ',
     :  4X,'Matrix element (a.u.)'/)
  305 FORMAT (1X,1I3,5X,2A4,2X,1I3,5X,2A4,1X,A4,4X,1P,1D20.10)
  309 FORMAT (A35,1P,1D20.10)
*
      call MPI_Get_processor_name (host, lenhost, ierr)
      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      END
