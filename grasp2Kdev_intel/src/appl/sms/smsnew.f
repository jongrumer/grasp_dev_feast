************************************************************************
*                                                                      *
      SUBROUTINE SMSNEW(VINT)
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
      DIMENSION VINT(NNNW,NNNW)
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
*   Call the MCP package to generate V coefficients; ac and bd
*   are the density pairs
*
*   Initialize
*
*
*   Matrix elements are diagonal in J
*
          IF (ITJPO(IR) .EQ. ITJPOC) THEN
            NVCOEF = 0
            CALL RKCO (IC,IR,COR,CORD,INCOR)
*
            DO 11 I = 1,NVCOEF
              VCOEFF = COEFF(I)
              IF (ABS (VCOEFF) .GT. CUTOFF) THEN
                IIA = LABEL(1,I)
                IIB = LABEL(2,I)
                IIC = LABEL(3,I)
                IID = LABEL(4,I)
                K   = LABEL(5,I)
*
*   Only K = 1   LABEL(5,I) .EQ. 1
*
                IF (LABEL(5,I) .EQ. 1) THEN
                  IF (LDBPA(2)) THEN
                    WRITE (99,309) K,IC,IR,
     :              NP(IIA),NH(IIA),NP(IIB),NH(IIB),
     :              NP(IIC),NH(IIC),NP(IID),NH(IID),VCOEFF
                  ENDIF
                  DO 10 J = 1,NVEC
                    LOC = (J-1)*NCF
                    CONTRI = - EVEC(IC+LOC)*EVEC(IR+LOC)*
     :                     VCOEFF*VINT(LABEL(1,I),LABEL(3,I))*
     :                     VINT(LABEL(2,I),LABEL(4,I))
                    IF (IR.NE.IC) CONTRI = 2.0D 00 * CONTRI
                    SMSC(J) = SMSC(J) + CONTRI
   10             CONTINUE
                ENDIF
              ENDIF
   11       CONTINUE
          ENDIF
   12   CONTINUE
   13 CONTINUE
*
*   Deallocate storage for the arrays in BUFFER
*
      CALL ALCBUF (3)
      RETURN
  309 FORMAT (' V^[(',1I2,')]_[',1I3,',',1I3,']',
     :   ' (',1I2,1A2,',',1I2,1A2,';',
     :        1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
*
      END
