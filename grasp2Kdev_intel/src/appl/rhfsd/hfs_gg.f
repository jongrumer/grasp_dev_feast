************************************************************************
*                                                                      *
      SUBROUTINE HFS(NAME)
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
*   Modified by Per Jonsson to evaluate g_j factors                    * 
*                                                                      *
*   Modified by Per Jonsson to read angular data from file             *
*                                                     AUGGUST 2011     *
*                                                                      *
*   Modified by G. Gaigalas by includeing the new spin-angular         *
*   libraries.                                    Vilnius May 2012     *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      PARAMETER (KEY = KEYORB)
      POINTER (PNTRIQ,RIQDUMMY)
      INTEGER FFMIN,FFMAX,FF
      CHARACTER*24 NAME
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,AVAIL
*
      DIMENSION TSHELL(NNNW),RINTME(2,NNNW,NNNW),AMELT(2,NNNW,NNNW)
      DIMENSION TSHELL_S(NNNW),IA_S(NNNW)
      DIMENSION RINTGJ(NNNW,NNNW),RINTDGJ(NNNW,NNNW)
      DIMENSION GJMELT(NNNW,NNNW),DGJMELT(NNNW,NNNW)
*
      POINTER (PNTHFC,HFC(5,*))
      POINTER (PNTGJC,GJC(*))
      POINTER (PNTDGJC,DGJC(*))
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
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
*   Allocate storage for local arrays
*
      CALL ALLOC (PNTHFC,5*NVEC*NVEC,8)
      CALL ALLOC (PNTGJC,NVEC*NVEC,8)
      CALL ALLOC (PNTDGJC,NVEC*NVEC,8)
*
*   Initialise
*
      DO 2 I = 1,5
        DO 1 J = 1,NVEC*NVEC
           HFC(I,J) = 0.0D 00
    1   CONTINUE
    2 CONTINUE
*
      DO J = 1,NVEC*NVEC
        GJC(J) = 0.0D 00
        DGJC(J) = 0.0D 00
      ENDDO

*
*   Calculate and save the radial integrals and angular
*   matrix elements for the two multipolarities
*
      DO 5 KT = 1,2
        DO 4 I = 1,NW
          DO 3 J = 1,NW
            IF (KT .EQ. 1) THEN
              RINTME(KT,I,J) = RINTHF (I,J,-2)
              RINTGJ(I,J) = RINTHF(I,J,1)
            ELSE
              RINTME(KT,I,J) = RINT (I,J,-3)
              RINTDGJ(I,J) = RINT(I,J,0)
            ENDIF
            CALL MATELT (I,KT,J,APART,GJPART,DGJPART)
            AMELT(KT,I,J) = APART
            IF (KT.EQ.1) THEN
              GJMELT(I,J) = GJPART
              DGJMELT(I,J) = DGJPART
            ENDIF
    3     CONTINUE
    4   CONTINUE
    5 CONTINUE
*
*   Set the parity of the one-body operators
*
      IPT = 1
*
*   Check if angular data is available and appropriate
*
      CALL ANGDATA(NAME,AVAIL)
     
      IF (.NOT.AVAIL) THEN
*
*   Sweep through the Hamiltonian matrix to determine the
*   diagonal and off-diagonal hyperfine constants
*
        DO 11 IC = 1,NCF
*
*   Output IC on the screen to show how far the calculation has preceede
*
          IF (MOD(IC,100) .EQ. 0) THEN
            CALL CONVRT (IC,CNUM,LCNUM)
            PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
          ENDIF
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
            IDIFF = ITJPOC - ITJPOR
*
*   Loop over the multipolarities
*
            DO 9 KT = 1,2
*
*   Initialise the accumulator
*
              ELEMNT = 0.0D 00
              ELEMNTGJ = 0.0D 00
              ELEMNTDGJ = 0.0D 00
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
!             IF (((IDIFF .EQ. 0) .AND. (IR. GE. IC)) .OR.
!    :          (IDIFF .EQ. 2) .OR.
!    :         ((IDIFF .EQ. 4) .AND. (KT .EQ. 2))) THEN
              IF ((IDIFF .EQ. 0) .AND. (IR. GE. IC)) THEN
*
                CALL ONEPARTICLEJJ(KT,IPT,IC,IR,IA,IB,TSHELL)
CGG                CALL TNSRJJ (KT,IPT,IC,IR,IA,IB,TSHELL)
*
*   Accumulate the contribution from the one-body operators;
*
                IF (IA .NE. 0) THEN
                  IF (IA .EQ. IB) THEN
                    NCOUNT = 0
                    DO 6 IA = 1,NW
                      IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                        NCOUNT = NCOUNT + 1
                        TSHELL_S(NCOUNT) = TSHELL(IA)
                        IA_S(NCOUNT) = IA
                        ELEMNT =   ELEMNT
     :                           + AMELT(KT,IA,IA)
     :                           * RINTME(KT,IA,IA)
     :                           * TSHELL(IA)
                        IF (KT.EQ.1.AND.(IDIFF .EQ. 0)) THEN
                          ELEMNTGJ = ELEMNTGJ 
     :                           + GJMELT(IA,IA)
     :                           * RINTGJ(IA,IA)
     :                           * TSHELL(IA)
                          ELEMNTDGJ = ELEMNTDGJ
     :                           + DGJMELT(IA,IA)
     :                           * RINTDGJ(IA,IA)
     :                           * TSHELL(IA)
                        ENDIF
                      ENDIF
    6               CONTINUE
                    IF (KT.EQ.1) THEN
                      WRITE(129) IC,IR,-NCOUNT
                    ELSE
                      WRITE(129) IC,IR,NCOUNT
                    END IF
                    DO I = 1,NCOUNT
                      LAB = IA_S(I)*KEY + IA_S(I)
C                      WRITE(129) TSHELL_S(I),IA_S(I),IA_S(I)
                      WRITE(129) TSHELL_S(I),LAB
                    END DO
                  ELSE
                    IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                      IF (KT.EQ.1) THEN
                        WRITE(129) IC,IR,-1
                      ELSE
                        WRITE(129) IC,IR,1
                      END IF
                      LAB = IA*KEY + IB
C                      WRITE(129) TSHELL(1),IA,IB
                      WRITE(129) TSHELL(1),LAB
                      ELEMNT =   ELEMNT
     :                       + AMELT(KT,IA,IB)
     :                       * RINTME(KT,IA,IB)
     :                       * TSHELL(1)
                      IF (KT.EQ.1.AND.(IDIFF .EQ. 0)) THEN
                        ELEMNTGJ = ELEMNTGJ
     :                       + GJMELT(IA,IB)
     :                       * RINTGJ(IA,IB)
     :                       * TSHELL(1) 
                        ELEMNTDGJ = ELEMNTDGJ
     :                       + DGJMELT(IA,IB)
     :                       * RINTDGJ(IA,IB)
     :                       * TSHELL(1) 
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
*
*   Multiply with the configuration expansion coefficients and add the
*   contributions from the matrix elements to obtain total contributions
*
                DO 8 K = 1,NVEC
                  DO 7 KK = K,K
                    LOC1 = (K-1)*NCF
                    LOC2 = (KK-1)*NCF
                    IF ((IDIFF .EQ. 0) .AND. (IR .NE. IC)) THEN
                       CONTR =   ELEMNT
     :                       * (  EVEC(IC+LOC1) * EVEC(IR+LOC2)
     :                          + EVEC(IR+LOC1) * EVEC(IC+LOC2) )
                       CONTRGJ = ELEMNTGJ
     :                       * (  EVEC(IC+LOC1) * EVEC(IR+LOC2)
     :                          + EVEC(IR+LOC1) * EVEC(IC+LOC2) )
                       CONTRDGJ = ELEMNTDGJ
     :                       * (  EVEC(IC+LOC1) * EVEC(IR+LOC2)
     :                          + EVEC(IR+LOC1) * EVEC(IC+LOC2) )
                    ELSE
                      CONTR =    ELEMNT
     :                       *    EVEC(IC+LOC1) * EVEC(IR+LOC2)
                      CONTRGJ = ELEMNTGJ
     :                       *    EVEC(IC+LOC1) * EVEC(IR+LOC2)
                      CONTRDGJ = ELEMNTDGJ
     :                       *    EVEC(IC+LOC1) * EVEC(IR+LOC2)
                    ENDIF
*
*   Magnetic dipole and the two operators of the g_j factor
*
                    IF (KT .EQ. 1) THEN
                      IF (IDIFF .EQ. 0) THEN
                        HFC(1,NVEC*(K-1)+KK) =   HFC(1,NVEC*(K-1)+KK) 
     :                                       + CONTR
                        GJC(NVEC*(K-1)+KK) = GJC(NVEC*(K-1)+KK) 
     :                                       + CONTRGJ
                        DGJC(NVEC*(K-1)+KK) = DGJC(NVEC*(K-1)+KK) 
     :                                       + CONTRDGJ
                      ELSEIF ((ITJPOC-ITJPOR) .EQ. 2) THEN
                        HFC(2,NVEC*(K-1)+KK) =   HFC(2,NVEC*(K-1)+KK)
     :                                       + CONTR
                      ENDIF
*
*   Electric quadrupole
*
                    ELSEIF (KT .EQ. 2) THEN
                      IF (IDIFF .EQ. 0) THEN
                        HFC(3,NVEC*(K-1)+KK) =   HFC(3,NVEC*(K-1)+KK)
     :                                       + CONTR
                      ELSEIF (IDIFF .EQ. 2) THEN
                        HFC(4,NVEC*(K-1)+KK) =   HFC(4,NVEC*(K-1)+KK)
     :                                       + CONTR
                      ELSEIF (IDIFF .EQ. 4) THEN
                        HFC(5,NVEC*(K-1)+KK) =   HFC(5,NVEC*(K-1)+KK)
     :                                       + CONTR
                      ENDIF
                    ENDIF
    7             CONTINUE
    8           CONTINUE
*
              ENDIF
    9       CONTINUE
*
   10     CONTINUE
   11   CONTINUE

      ELSE

*    Read angular data from file and compute matrix elements

        ICOLD = 0
        IROLD = 0
        KTOLD = 0
        DO 
          READ(129,END=999) IC,IR,NCOUNT
          IF (NCOUNT.LT.0) THEN
            KT = 1
            NCOUNT = -NCOUNT
          ELSE
            KT = 2
          END IF
*
*   Initialise the accumulator
*
          IF ((IC.NE.ICOLD).OR.(IR.NE.IROLD).OR.(KT.NE.KTOLD)) THEN
            ISPARC = ISPAR (IC)
            ITJPOC = ITJPO (IC)
            ITJPOR = ITJPO (IR)
            IDIFF = ITJPOC - ITJPOR
            ICOLD = IC
            IROLD = IR
            KTOLD = KT
            ELEMNT = 0.0D 00
            ELEMNTGJ = 0.0D 00
            ELEMNTDGJ = 0.0D 00
          END IF
*
*   Accumulate the contribution from the one-body operators;
*
          DO I = 1,NCOUNT
C            READ(129) TSHELL_R,IA,IB
            READ(129) TSHELL_R,LAB
            IA = LAB/KEY
            IB = MOD(LAB,KEY)
            ELEMNT =   ELEMNT
     :          + AMELT(KT,IA,IB)
     :          * RINTME(KT,IA,IB)
     :          * TSHELL_R
            IF (KT.EQ.1.AND.(IDIFF .EQ. 0)) THEN
              ELEMNTGJ = ELEMNTGJ
     :             + GJMELT(IA,IB)
     :             * RINTGJ(IA,IB)
     :             * TSHELL_R 
              ELEMNTDGJ = ELEMNTDGJ
     :             + DGJMELT(IA,IB)
     :             * RINTDGJ(IA,IB)
     :             * TSHELL_R
            END IF
          END DO
*
*   Multiply with the configuration expansion coefficients and add the
*   contributions from the matrix elements to obtain total contributions
*
          DO K = 1,NVEC
            DO KK = K,K
              LOC1 = (K-1)*NCF
              LOC2 = (KK-1)*NCF
              IF ((IDIFF .EQ. 0) .AND. (IR .NE. IC)) THEN
                 CONTR =   ELEMNT
     :                 * (  EVEC(IC+LOC1) * EVEC(IR+LOC2)
     :                    + EVEC(IR+LOC1) * EVEC(IC+LOC2) )
                 CONTRGJ = ELEMNTGJ
     :                 * (  EVEC(IC+LOC1) * EVEC(IR+LOC2)
     :                    + EVEC(IR+LOC1) * EVEC(IC+LOC2) )
                 CONTRDGJ = ELEMNTDGJ
     :                 * (  EVEC(IC+LOC1) * EVEC(IR+LOC2)
     :                    + EVEC(IR+LOC1) * EVEC(IC+LOC2) )
              ELSE
                 CONTR =    ELEMNT
     :                 *    EVEC(IC+LOC1) * EVEC(IR+LOC2)
                 CONTRGJ = ELEMNTGJ
     :                 *    EVEC(IC+LOC1) * EVEC(IR+LOC2)
                 CONTRDGJ = ELEMNTDGJ
     :                 *    EVEC(IC+LOC1) * EVEC(IR+LOC2)
              ENDIF
*
*   Magnetic dipole and the two operators of the g_j factor
*
              IF (KT .EQ. 1) THEN
                IF (IDIFF .EQ. 0) THEN
                  HFC(1,NVEC*(K-1)+KK) =   HFC(1,NVEC*(K-1)+KK) 
     :                                 + CONTR
                  GJC(NVEC*(K-1)+KK) = GJC(NVEC*(K-1)+KK) 
     :                                 + CONTRGJ
                  DGJC(NVEC*(K-1)+KK) = DGJC(NVEC*(K-1)+KK) 
     :                                 + CONTRDGJ
                ELSEIF ((ITJPOC-ITJPOR) .EQ. 2) THEN
                  HFC(2,NVEC*(K-1)+KK) =   HFC(2,NVEC*(K-1)+KK)
     :                                 + CONTR
                ENDIF
*
*   Electric quadrupole
*
              ELSEIF (KT .EQ. 2) THEN
                IF (IDIFF .EQ. 0) THEN
                  HFC(3,NVEC*(K-1)+KK) =   HFC(3,NVEC*(K-1)+KK)
     :                                 + CONTR
                ELSEIF (IDIFF .EQ. 2) THEN
                  HFC(4,NVEC*(K-1)+KK) =   HFC(4,NVEC*(K-1)+KK)
     :                                 + CONTR
                ELSEIF (IDIFF .EQ. 4) THEN
                  HFC(5,NVEC*(K-1)+KK) =   HFC(5,NVEC*(K-1)+KK)
     :                                 + CONTR
                ENDIF
              ENDIF
            END DO 
          END DO 
        END DO
  999   CONTINUE
      END IF

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
!     WRITE (24,302)
      WRITE (29,402)
*
      DO 13 I = 1,NVEC
        DO 12 II = I,I
*
          JJ = IATJPO(I)
          JJII = IATJPO(II)
*
          IF ((JJ.EQ.JJII .AND. JJII.GT.1) .OR. (JJ.GT.JJII)) THEN
*
            FJ = 0.5D 00*DBLE (JJ-1)
*
            GJA1 =  SQRT (
     :                       1.0D 00
     :                     / (FJ*(FJ+1.0D 00))
     :                    )

            AFA1 =   GFAC*GJA1
            
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
              IF (I .LE. II) THEN
C                GJ = CVAC*GJA1*GJC(NVEC*(I-1)+II)
C                DGJ = 0.001160D0*GJA1*DGJC(NVEC*(I-1)+II)  
!               WRITE (24,303) IVEC(I),LABJ(JJ),
!    :                         LABP((IASPAR(I)+3)/2),
!    :                         IVEC(II),LABJ(JJII),
!    :                         LABP((IASPAR(II)+3)/2),
!    :                         AFA1*HFC(1,NVEC*(I-1)+II),
!    :                         BFA1*HFC(3,NVEC*(I-1)+II)
*  
*   Output diagonal hfs and g_j factors to file <name>.h or <name>.ch
*
                IF (I .EQ. II) THEN 
                  GJ = CVAC*GJA1*GJC(NVEC*(I-1)+II)
                  DGJ = 0.001160D0*GJA1*DGJC(NVEC*(I-1)+II)  
                  WRITE (29,403) IVEC(I),LABJ(JJ),
     :                         LABP((IASPAR(I)+3)/2),
     :                         AFA1*HFC(1,NVEC*(I-1)+II),
     :                         BFA1*HFC(3,NVEC*(I-1)+II),
     :                         GJ+DGJ
                ENDIF
              ENDIF
*
*   Off diagonal (J,J-1) A and B factors
*
            ELSEIF ((JJ-JJII) .EQ. 2) THEN
!             WRITE (24,303) IVEC(I),LABJ(JJ),
!    :                       LABP((IASPAR(I)+3)/2),
!    :                       IVEC(II),LABJ(JJII),
!    :                       LABP((IASPAR(II)+3)/2),
!    :                       AFA2*HFC(2,NVEC*(I-1)+II),
!    :                       BFA2*HFC(4,NVEC*(I-1)+II)
*
*   Off diagonal (J,J-2) B factor
*
            ELSEIF ((JJ-JJII) .EQ. 4) THEN
!             WRITE (24,303) IVEC(I),LABJ(JJ),
!    :                       LABP((IASPAR(I)+3)/2),
!    :                       IVEC(II),LABJ(JJII),
!    :                       LABP((IASPAR(II)+3)/2),0.0D 00,
!    :                       BFA3*HFC(5,NVEC*(I-1)+II)
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
!     DO 16 JB = 1,NVEC
!       DO 15 JA = 1,NVEC
!         JJB = IATJPO(JB)-1
!         JJA = IATJPO(JA)-1
!         IF ((JJA .EQ. JJB .AND. JJA .GT. 0) .OR.
!    :        (JJB .GT. JJA)) THEN
*
*   Determine the possible F quantum numbers for the matrix element
*
!           FFMIN = MAX (ABS (JJA-II),ABS (JJB-II))
!           FFMAX = MIN (JJA+II,JJB+II)
*
*   Loop over the possible F quantum numbers
*
!           IFLAG = 0
!           DO 14 FF = FFMIN,FFMAX,2
*
*   Phase factor
*
!             IF (MOD ((II+JJA-FF)/2,2) .EQ. 1) THEN
!               FACTOR1 = -1.0D 00
!             ELSE
!               FACTOR1 =  1.0D 00
!             ENDIF
!             FACTOR2 =   FACTOR1
!    :                  * SQRT (
!    :                            ( DBLE (JJB)
!    :                             +1.0D 00)
!    :                          * ( DBLE (II)
!    :                             +1.0D 00)
!    :                         )
*
*   Determine the Racah W coefficients.
*
CGG              CALL DRACAH (II,JJA,II,JJB,FF,2,RAC1)
!             CALL GRACAH1 (II,JJA,II,JJB,FF,2,RAC1)
CGG              CALL DRACAH (II,JJA,II,JJB,FF,4,RAC2)
!             CALL GRACAH1 (II,JJA,II,JJB,FF,4,RAC2)
*
*   Obtain and output matrix elements for J,J
*
!             IF ((JJA-JJB) .EQ. 0) THEN
!               HFSELT1 = FACTOR2*RAC1*HFC(1,NVEC*(JB-1)+JA)*TILDE1
!               HFSELT2 = FACTOR2*RAC2*HFC(3,NVEC*(JB-1)+JA)*TILDE2
!               IF (ABS (HFSELT1+HFSELT2) .GT.
!    :              CUTOFF*10.0D-05) THEN
!                 IF (IFLAG .EQ. 0) WRITE (24,304)
!                 IFLAG = 1
!                 WRITE (24,305) IVEC(JB),LABJ(JJB+1),
!    :                           LABP((IASPAR(JB)+3)/2),
!    :                           IVEC(JA),LABJ(JJA+1),
!    :                           LABP((IASPAR(JA)+3)/2),
!    :                           LABJ(FF+1),
!    :                           HFSELT1+HFSELT2
!               ENDIF
!
!   Obtain and output matrix elements for J,J-1
!
!             ELSEIF (ABS (JJA-JJB) .EQ. 2) THEN
!               HFSELT1 = FACTOR2*RAC1*HFC(2,NVEC*(JB-1)+JA)*TILDE1
!               HFSELT2 = FACTOR2*RAC2*HFC(4,NVEC*(JB-1)+JA)*TILDE2
!               IF (ABS (HFSELT1+HFSELT2) .GT.
!    :              CUTOFF*10.0D-05) THEN
!                 IF (IFLAG .EQ. 0) WRITE (24,304)
!                 IFLAG = 1
!                 WRITE (24,305) IVEC(JB),LABJ(JJB+1),
!    :                           LABP((IASPAR(JB)+3)/2),
!    :                           IVEC(JA),LABJ(JJA+1),
!    :                           LABP((IASPAR(JA)+3)/2),
!    :                           LABJ(FF+1),
!    :                           HFSELT1+HFSELT2
!               ENDIF
*
*   Obtain and output matrix elements for J,J-2
*
!             ELSEIF (ABS (JJA-JJB) .EQ. 4) THEN
!               HFSELT1 = 0.0D 00
!               HFSELT2 = FACTOR2*RAC2*HFC(5,NVEC*(JB-1)+JA)*TILDE2
!               IF (ABS (HFSELT1+HFSELT2) .GT.
!    :              CUTOFF*10.0D-05) THEN
!                 IF (IFLAG .EQ. 0) WRITE (24,304)
!                 IFLAG = 1
!                 WRITE (24,305) IVEC(JB),LABJ(JJB+1),
!    :                           LABP((IASPAR(JB)+3)/2),
!    :                           IVEC(JA),LABJ(JJA+1),
!    :                           LABP((IASPAR(JA)+3)/2),
!    :                           LABJ(FF+1),
!    :                           HFSELT1+ HFSELT2
!               ENDIF
!             ELSE
!               HFSELT1 = 0.0D 00
!               HFSELT2 = 0.0D 00
!             ENDIF
! 14        CONTINUE
!         ENDIF
! 15    CONTINUE
! 16  CONTINUE
*
      CALL DALLOC (PNTHFC)
      RETURN
*
  302 FORMAT (//' Interaction constants:'//
     :' Level1  J Parity  Level2  J Parity',8X,'A (MHz)',13X,'B (MHz)'/)
CPJ  402 FORMAT (//' Interaction constants:'//
CPJ     :' Level1  J Parity ',8X,'A (MHz)',13X,'B (MHz)',13X,'g_J',14X,
CPJ     : 'delta g_J',11X,'total g_J'/)
  402 FORMAT (//' Interaction constants:'//
     :' Level1  J Parity ',8X,'A (MHz)',13X,'B (MHz)',13X,'g_J'/)
  303 FORMAT (1X,1I3,5X,2A4,2X,1I3,5X,2A4,1P,2D20.10)
CPJ  403 FORMAT (1X,1I3,5X,2A4,1P,5D20.10)
  403 FORMAT (1X,1I3,5X,2A4,1P,3D20.10)
  304 FORMAT (//' Matrix elements:'//
     :' Level1  J Parity  Level2  J Parity    F ',
     :  4X,'Matrix element (a.u.)'/)
  305 FORMAT (1X,1I3,5X,2A4,2X,1I3,5X,2A4,1X,A4,4X,1P,1D20.10)
*
      END
