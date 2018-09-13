************************************************************************
*                                                                      *
      SUBROUTINE HFSZEEMAN(NOFFD)
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the  hyperfine structure parameters.               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, DRACAH, ISPAR, ITJPO,  *
*                        TNSRJJ.                                       *
*               [HFSZEEMAN05]: MATELT, RINT, RINTHF.                   *
*                                                                      *
*   Written by:  Martin Andersson and Per Jonsson                      *
*                                                                      *
*   This program is a modified version of HFS written by:              *
*                              Per Jonsson and Farid A. Parpia         *
*                                                                      *
*                                         Last revision: 27 July 2006  *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CPJ      PARAMETER (NNNW = 120)
      POINTER (PNTRIQ,RIQDUMMY)
      INTEGER FFMIN,FFMAX,FF
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION TSHELL(NNNW),RINTME(2,NNNW,NNNW),AMELT(2,NNNW,NNNW)
      DIMENSION RINTGJ(NNNW,NNNW),RINTDGJ(NNNW,NNNW)
      DIMENSION GJMELT(NNNW,NNNW),DGJMELT(NNNW,NNNW)
      DIMENSION HFCMI(NVEC,NVEC), SORTHFCMI(NVEC,NVEC)
      DIMENSION HFCMII(NVEC,NVEC), SORTHFCMII(NVEC,NVEC)
      DIMENSION GJCM(NVEC,NVEC)
      DIMENSION DGJCM(NVEC,NVEC)
      DIMENSION TOTGJ(NVEC,NVEC), SORTTOTGJ(NVEC,NVEC)
      DIMENSION NSORT(NVEC), NEWSORT(NVEC)
*
      POINTER (PNTHFC,HFC(5,*))
      POINTER (PNTGJC,GJC(2,*))
      POINTER (PNTDGJC,DGJC(2,*))
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
      CALL ALLOC (PNTGJC,2*NVEC*NVEC,8)
      CALL ALLOC (PNTDGJC,2*NVEC*NVEC,8)
*
*   Initialise
*
      DO I = 1,5
        DO J = 1,NVEC*NVEC
           HFC(I,J) = 0.0D 00
        ENDDO
      ENDDO
*
      DO I = 1,2
        DO J = 1,NVEC*NVEC
          GJC(I,J) = 0.0D 00
          DGJC(I,J) = 0.0D 00
        ENDDO
      ENDDO
*
      DO I = 1,NVEC
        DO J = 1,NVEC
          HFCMI(I,J) = 0.0D 00
          HFCMII(I,J) = 0.0D 00
          GJCM(I,J) = 0.0D 00
          DGJCM(I,J) = 0.0D 00
          TOTGJ(I,J) = 0.0D 00
        ENDDO
      ENDDO
*
*   Calculate and save the radial integrals and angular
*   matrix elements for the two multipolarities
*
      DO KT = 1,2
        DO I = 1,NW
          DO J = 1,NW
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
          ENDDO
        ENDDO
      ENDDO
*
*   Set the parity of the one-body operators
*
      IPT = 1
*
*   Sweep through the Hamiltonian matrix to determine the
*   diagonal and off-diagonal hyperfine constants
*
      DO IC = 1,NCF
*
*   Output IC on the screen to show how far the calculation has preceede
*
        IF (MOD(IC,100) .EQ. 0) THEN
          CALL CONVRT (IC,CNUM,LCNUM)
          PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
        ENDIF
*
        DO IR = 1,NCF
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
          DO KT = 1,2
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
            IF (((IDIFF .EQ. 0) .AND. (IR. GE. IC)) .OR.
     :          (IDIFF .EQ. 2) .OR.
     :         ((IDIFF .EQ. 4) .AND. (KT .EQ. 2))) THEN
*
              CALL TNSRJJ (KT,IPT,IC,IR,IA,IB,TSHELL)
*
*   Accumulate the contribution from the one-body operators;
*
              IF (IA .NE. 0) THEN
                IF (IA .EQ. IB) THEN
                  DO IA = 1,NW
                    IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                      ELEMNT =   ELEMNT
     :                         + AMELT(KT,IA,IA)
     :                         * RINTME(KT,IA,IA)
     :                         * TSHELL(IA)
                      IF (KT.EQ.1) THEN
                        ELEMNTGJ = ELEMNTGJ 
     :                         + GJMELT(IA,IA)
     :                         * RINTGJ(IA,IA)
     :                         * TSHELL(IA)
                        ELEMNTDGJ = ELEMNTDGJ
     :                         + DGJMELT(IA,IA)
     :                         * RINTDGJ(IA,IA)
     :                         * TSHELL(IA)
                      ENDIF
                    ENDIF
                  ENDDO
                ELSE
                  IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                    ELEMNT =   ELEMNT
     :                       + AMELT(KT,IA,IB)
     :                       * RINTME(KT,IA,IB)
     :                       * TSHELL(1)
                    IF (KT.EQ.1) THEN
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
              DO K = 1,NVEC
*
*   These two variables are used to separate calculations where the user 
*   is only interested in the diagonal matrix elements
*
              IF (NOFFD .EQ. 0) THEN
                 NLOW = 1
                 NUPP = NVEC
              ELSE
                 NLOW = K
                 NUPP = K
              ENDIF 
                DO KK = NLOW,NUPP
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
                      GJC(1,NVEC*(K-1)+KK) = GJC(1,NVEC*(K-1)+KK) 
     :                                       + CONTRGJ
                      DGJC(1,NVEC*(K-1)+KK) = DGJC(1,NVEC*(K-1)+KK) 
     :                                       + CONTRDGJ
                    ELSEIF ((ITJPOC-ITJPOR) .EQ. 2) THEN
                      HFC(2,NVEC*(K-1)+KK) =   HFC(2,NVEC*(K-1)+KK)
     :                                       + CONTR
                      GJC(2,NVEC*(K-1)+KK) = GJC(2,NVEC*(K-1)+KK) 
     :                                       + CONTRGJ
                      DGJC(2,NVEC*(K-1)+KK) = DGJC(2,NVEC*(K-1)+KK) 
     :                                       + CONTRDGJ
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
               ENDDO
             ENDDO
*
            ENDIF
          ENDDO
*
   10   ENDDO
      ENDDO
*
*
*
*  Produce matrices for the PRINTHFSZEEMAN matlab program
*  If user only interested in the diagonal elements this
*  does not have to be done
      IF (NOFFD .EQ. 0) THEN
        DO I = 1,NVEC
          DO II=1,NVEC
            HFCMI(I,II) = HFC(1,NVEC*(I-1)+II) + HFC(2,NVEC*(I-1)+II)
            HFCMII(I,II) = HFC(3,NVEC*(I-1)+II) + HFC(4,NVEC*(I-1)+II)
     :                                   + HFC(5,NVEC*(I-1)+II)
            GJCM(I,II) = GJC(1,NVEC*(I-1)+II) + GJC(2,NVEC*(I-1)+II)
            DGJCM(I,II) = DGJC(1,NVEC*(I-1)+II) + DGJC(2,NVEC*(I-1)+II)
*   Multiply by 0.5 since the operator is 1/2N, but the matrix elements are
*  calculated for the operator N
            TOTGJ(I,II) = 0.5*(CVAC*GJCM(I,II) + 0.001160D0*DGJCM(I,II))
          ENDDO
        ENDDO
* 
        DO I = 1,NVEC
          DO J = 1,NVEC
            HFCMI(I,J) = HFCMI(J,I)
            HFCMII(I,J) = HFCMII(J,I)
            TOTGJ(I,J) = TOTGJ(J,I)
          ENDDO
        ENDDO
*
        DO I = 1,NVEC
          NSORT(I) = I
        ENDDO

        DO I = IATJPO(1),IATJPO(NVEC),2
          NSTART = 0
          NUMBER = 0
          DO J = 1,NVEC
            IF (IATJPO(J) == I) THEN
              NUMBER = NUMBER + 1
              IF (NSTART == 0) THEN
                NSTART = J
              ENDIF
            ENDIF
          ENDDO
          NSTOP = NSTART + NUMBER -1
          K = 0
          DO L = NSTART,NSTOP
            K = K + 1
            NSORT(L) = NSTOP - K + 1
          ENDDO
        ENDDO
*
        DO I = 1,NVEC
          NEWSORT(I) = NSORT(NVEC-I+1)
        ENDDO
**
*        DO I = 1,NVEC
*          DO J = 1,NVEC
*            SORTHFCMI(I,J) = HFCMI(NEWSORT(I),NEWSORT(J))
*            SORTHFCMII(I,J) = HFCMII(NEWSORT(I),NEWSORT(J))
*            SORTTOTGJ(I,J) = TOTGJ(NEWSORT(I),NEWSORT(J))
*          ENDDO
*        ENDDO
** 
*
***     Fixed so that for the lower part 
***     <JJ||H||JI> = (-1)^(JJ-JI)sqrt((2JI+1)/(2JJ+1))<JI||H||JJ>
*
        DO I = 1,NVEC
          DO J = I,NVEC
            SORTHFCMI(I,J) = HFCMI(NEWSORT(I),NEWSORT(J))
            SORTHFCMII(I,J) = HFCMII(NEWSORT(I),NEWSORT(J))
            SORTTOTGJ(I,J) = TOTGJ(NEWSORT(I),NEWSORT(J)) 
********   IATJPO(J)-(IATJPO(I) = 2JJ + 1 - 2JI +1 = 2(JJ-JI)
            IF (I .NE. J) THEN
              IF ( MOD((IATJPO(NEWSORT(J)) - 
     :                   IATJPO(NEWSORT(I))),4) .EQ. 0) THEN
                phase = 1
              ELSE
                phase = -1
              ENDIF
                factor = DSQRT(DBLE(IATJPO(NEWSORT(I)))
     :                                / DBLE(IATJPO(NEWSORT(J))))
                SORTHFCMI(J,I) = phase * factor 
     :                             * HFCMI(NEWSORT(I),NEWSORT(J))
                SORTHFCMII(J,I) = phase * factor 
     :                             * HFCMII(NEWSORT(I),NEWSORT(J))
                SORTTOTGJ(J,I) = phase * factor 
     :                             * TOTGJ(NEWSORT(I),NEWSORT(J)) 
            ENDIF
          ENDDO
        ENDDO              
*
        WRITE(111,301)
        WRITE(111,302) NVEC
        WRITE(111,303)
        DO I = 1,NVEC
          EAU = EVAL(NEWSORT(I))+EAV
          WRITE(111,304) IVEC(NEWSORT(I)),(IATJPO(NEWSORT(I))-1)/2.0,
     :                       LABP((IASPAR(NEWSORT(I))+3)/2),EAU
        ENDDO
        WRITE(111,305)
        DO I=1,NVEC
          WRITE(111,401) (SORTTOTGJ(I,J), J = 1,NVEC)
        ENDDO
        WRITE(111,306)
        DO I=1,NVEC
          WRITE(111,401) (SORTHFCMI(I,J), J = 1,NVEC)
        ENDDO
        WRITE(111,307)
        DO I=1,NVEC
          WRITE(111,401) (SORTHFCMII(I,J), J = 1,NVEC)
        ENDDO
      ENDIF
*
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
      IF (HFSI.NE.0) THEN
        WRITE (29,402)
      ELSE
        WRITE (29,404)
      ENDIF
*
      DO  I = 1,NVEC
        DO II = 1,NVEC
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

            BFA1 =   HFAC
     :             * SQRT (
     :                       (FJ*(2.0D 00*FJ-1.0D 00))
     :                     / (   (FJ+1.0D 00)
     :                         * (2.0D 00*FJ+3.0D 00) )
     :                    )
*
            IF (HFSI.NE.0) THEN
*
*   Output diagonal hfs and g_j factors to file <name>.h or <name>.ch
*
              IF ((JJ-JJII) .EQ. 0) THEN
                IF (I .EQ. II) THEN
                  GJ = CVAC*GJA1*GJC(1,NVEC*(I-1)+II)
                  DGJ = 0.001160D0*GJA1*DGJC(1,NVEC*(I-1)+II)
                  WRITE (29,403) IVEC(I),LABJ(JJ),
     :                         LABP((IASPAR(I)+3)/2),
     :                         AFA1*HFC(1,NVEC*(I-1)+II),
     :                         BFA1*HFC(3,NVEC*(I-1)+II),
     :                         GJ,
     :                         DGJ,
     :                         GJ+DGJ
                ENDIF
              ENDIF
            ELSE
*
*   Output diagonal g_j factors to file <name>.h or <name>.ch
*
              IF ((JJ-JJII) .EQ. 0) THEN
                IF (I .EQ. II) THEN
                  GJ = CVAC*GJA1*GJC(1,NVEC*(I-1)+II)
                  DGJ = 0.001160D0*GJA1*DGJC(1,NVEC*(I-1)+II)
                  WRITE (29,405) IVEC(I),LABJ(JJ),
     :                         LABP((IASPAR(I)+3)/2),
     :                         GJ,DGJ,GJ+DGJ
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
*
      CALL DALLOC (PNTHFC)
      RETURN
*
  301 FORMAT ('  Number of relativistic eigenvalues')
  302 FORMAT (I4)
  303 FORMAT  ('  Lev     J  Parity       E')
  304 FORMAT (1X,1I3,2X,F6.1,2X,1A4,F16.9)
  305 FORMAT  ('  Zeeman interaction matrix  ')
  306 FORMAT  ('  HFI-matrix for the magnetic dipole operator  ')
  307 FORMAT  ('  HFI-matrix for the electric quadrupole operator  ')
  401 FORMAT (99E13.5)
  402 FORMAT (//' Interaction constants:'//
     :' Level1  J Parity ',8X,'A (MHz)',13X,'B (MHz)',13X,'g_J',14X,
     : 'delta g_J',11X,'total g_J'/)
  403 FORMAT (1X,1I3,5X,2A4,1P,5D20.10)
  404 FORMAT (//' Interaction constants:'//
     :' Level1  J Parity ',10X,'g_J',14X,'delta g_J',11X,'total g_J'/)
  405 FORMAT (1X,1I3,5X,2A4,1P,3D20.10)
*
      END
