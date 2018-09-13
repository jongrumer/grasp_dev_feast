************************************************************************
*                                                                      *
      SUBROUTINE DENSNEW(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6)
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the  sms parameter, the electron density at the    *
*   origin.
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, GETYN          *
*                        ITJPO, RKCO, TNSRJJ                           *
*               [SMS92]: RINTISO, RINTDENS, VINTI                      *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
*                                         Last revision: 10 Nov 1995   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ,
Cww     :        PINDTE,PVALTE
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PINDTE,INDTEDUMMY)
      POINTER (PVALTE,VALTEDUMMY)
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      CHARACTER*2 CK,NH
      LOGICAL GETYN,FIRSTT,LDBPA,VSH,NUCDE,SMSSH,YES
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION TSHELL(NNNW)
      DIMENSION DINT1(NNNW,NNNW),DINT2(NNNW,NNNW),DINT3(NNNW,NNNW)
      DIMENSION DINT4(NNNW,NNNW),DINT5(NNNW,NNNW),DINT6(NNNW,NNNW)
*
      POINTER (PNSMS,SMSC(*))
      POINTER (PNDENS1,DENS1(*))
      POINTER (PNDENS2,DENS2(*))
      POINTER (PNDENS3,DENS3(*))
      POINTER (PNDENS4,DENS4(*))
      POINTER (PNDENS5,DENS5(*))
      POINTER (PNDENS6,DENS6(*))
      POINTER (PLABEL,LABEL(6,*))
      POINTER (PCOEFF,COEFF(*))
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      EXTERNAL COR,CORD
      COMMON/DEBUGA/LDBPA(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF3/EMPAM,RBCM
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /DEF11/FMTOAU,B1
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /NPAR/PARM(2),NPARM
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /OPT6/NTC(10)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /TEILST/NDTEA,NTEI,PINDTE,PVALTE,FIRSTT
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /SMS1/PNSMS,PNDENS1,PNDENS2,PNDENS3,PNDENS4,PNDENS5,PNDENS6
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
*   Set the rank (zero) and parity (even) for the one-particle
*   coefficients
*
      KA = 0
      IOPAR = 1
      INCOR = 1
*
*   Allocate storage for the arrays in BUFFER
*
      CALL ALCBUF (1)
*
*   Sweep through the Hamiltonian matrix to determine the
*   sms parameter
*
      DO 13 IC = 1,NCF
*
*   Output IC on the screen to show how far the calculation has preceede
*
        CALL CONVRT (IC,CNUM,LCNUM)
        if (mod(IC,10).eq.0) then
          PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
        end if
*
        ITJPOC = ITJPO (IC)
        DO 12 IR = IC,NCF
*
*   Matrix elements are diagonal in J
*
          IF (ITJPO(IR) .EQ. ITJPOC) THEN
*
*   Initialise the accumulator
*
            ELEMNT1 = 0.0D 00
            ELEMNT2 = 0.0D 00
            ELEMNT3 = 0.0D 00
            ELEMNT4 = 0.0D 00
            ELEMNT5 = 0.0D 00
            ELEMNT6 = 0.0D 00
*
*   Call the MCT package to compute T coefficients
*
            CALL TNSRJJ (KA,IOPAR,IC,IR,IA,IB,TSHELL)

            IF (IA .NE. 0) THEN
              IF (IA .EQ. IB) THEN
                DO 8 IA = 1,NW
                  IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                    ELEMNT1 = ELEMNT1 + DINT1(IA,IA)*TSHELL(IA)
                    ELEMNT2 = ELEMNT2 + DINT2(IA,IA)*TSHELL(IA)
                    ELEMNT3 = ELEMNT3 + DINT3(IA,IA)*TSHELL(IA)
                    ELEMNT4 = ELEMNT4 + DINT4(IA,IA)*TSHELL(IA)
                    ELEMNT5 = ELEMNT5 + DINT5(IA,IA)*TSHELL(IA)
                    ELEMNT6 = ELEMNT6 + DINT6(IA,IA)*TSHELL(IA)
                  ENDIF
    8           CONTINUE
              ELSE
                IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                  IF (NAK(IA).EQ.NAK(IB)) THEN
                    ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL(1)
                    ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL(1)
                    ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL(1)
                    ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL(1)
                    ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL(1)
                    ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL(1)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            DO 9 J = 1,NVEC
              LOC = (J-1)*NCF
              CONTRI1 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1
              CONTRI2 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT2
              CONTRI3 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT3
              CONTRI4 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT4
              CONTRI5 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT5
              CONTRI6 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT6
              IF (IR.NE.IC) THEN
                CONTRI1 = 2.0D 00 * CONTRI1
                CONTRI2 = 2.0D 00 * CONTRI2
                CONTRI3 = 2.0D 00 * CONTRI3
                CONTRI4 = 2.0D 00 * CONTRI4
                CONTRI5 = 2.0D 00 * CONTRI5
                CONTRI6 = 2.0D 00 * CONTRI6
              ENDIF
              DENS1(J) = DENS1(J) + CONTRI1
              DENS2(J) = DENS2(J) + CONTRI2
              DENS3(J) = DENS3(J) + CONTRI3
              DENS4(J) = DENS4(J) + CONTRI4
              DENS5(J) = DENS5(J) + CONTRI5
              DENS6(J) = DENS6(J) + CONTRI6
    9       CONTINUE
          ENDIF
   12   CONTINUE
   13 CONTINUE
*
*   Deallocate storage for the arrays in BUFFER
*
      CALL ALCBUF (3)
      RETURN
      END
