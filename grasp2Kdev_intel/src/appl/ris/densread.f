************************************************************************
*                                                                      *
      SUBROUTINE DENSREAD(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
*                                                                      *
*   IF angular coefficients already exist                              *
*   This routine controls combines the radial and angular parts for the*
*   calculation of the NMS parameter, the electron density at the      *
*   origin and radial expectation values.                              *
*                                                                      *
*   Written by C\'edric Naz\'e                                         *
*                                                                      *
*                                         Last revision: Feb. 2011     *
************************************************************************
* DINT1 contain the density
* DINT2 contain the uncorrected NMS parameter: K^1_NMS
* DINT3 contain the expect. value <r>
* DINT4 contain the expect. value <r2>
* DINT5 contain the expect. value <r-1>
* DINT6 contain the expect. value <r-2>
* DINT7 contain the sum of NMS parameters: K^1_NMS+K^2_NMS+K^3_NMS
*

      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
*   Key is used to store the indices of a couple of orbitals
      PARAMETER (KEY = KEYORB)
*   Matrix elements smaller than CUTOFF are not accumulated
      PARAMETER (CUTOFF = 1.0D-10)

*
      DIMENSION DINT1(NNNW,NNNW),DINT2(NNNW,NNNW),DINT3(NNNW,NNNW)
      DIMENSION DINT4(NNNW,NNNW),DINT5(NNNW,NNNW),DINT6(NNNW,NNNW)
      DIMENSION DINT7(NNNW,NNNW)
      DIMENSION TSHELL_R(NNNW)
*********,IA_S(NNNW)

* The electron density and the NMS
      POINTER (PNDENS1,DENS1(*)),(PNDENS2,DENS2(*))
* The radial axpectation values R,R2
      POINTER (PNDENS3,DENS3(*)),(PNDENS4,DENS4(*))
*  and R-1 and R-2
      POINTER (PNDENS5,DENS5(*)),(PNDENS6,DENS6(*))
      POINTER (PNDENS7,DENS7(*))

*      POINTER (PNTEMT1,ELEMNT1(*)),(PNTEMT2,ELEMNT2(*))
**      POINTER (PNTEMT3,ELEMNT3(*)),(PNTEMT4,ELEMNT4(*))
**      POINTER (PNTEMT5,ELEMNT5(*)),(PNTEMT6,ELEMNT6(*))
      POINTER (PINDEX,INDEX(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/PRNT/NVEC,PNIVEC,NVECMX
     :      /EIGVEC/PNEVEC
     :      /SMS1/PNSMS1,PNSMS2,PNDENS1,PNDENS2,PNDENS3,
     :          PNDENS4,PNDENS5,PNDENS6,PNDENS7
     :      /ORB2/NCF,NW,PNTRIQ


      ICOLD = 0
      IROLD = 0

      REWIND (50)
   16 READ (50,IOSTAT = IOS) IC,IR,NCOUNT

      IF (IOS .EQ. 0) THEN
**      print*, 'ic', IC,IR,NCOUNT
*                                                        
*   Read successful; decode the labels of I(ab) 
*
*      Initialise
        IF ((IC.NE.ICOLD).OR.(IR.NE.IROLD)) THEN
**            ISPARC = ISPAR (IC)
**            ITJPOC = ITJPO (IC)
**            ITJPOR = ITJPO (IR)
**            IDIFF = ITJPOC - ITJPOR
            ICOLD = IC
            IROLD = IR
**        DO I = 1,NELMNT
          ELEMNT1 = 0.0D 00
          ELEMNT2 = 0.0D 00
          ELEMNT3 = 0.0D 00
          ELEMNT4 = 0.0D 00
          ELEMNT5 = 0.0D 00
          ELEMNT6 = 0.0D 00
          ELEMNT7 = 0.0D 00
**        ENDDO
        ENDIF
        DO I= 1,NCOUNT
           READ(50) TSHELL_R(I),LAB
           IA = LAB/KEY
           IB = MOD(LAB,KEY)
**           LOC = INDEX(I)
           ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL_R(I)
           ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL_R(I)
           ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL_R(I)
           ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL_R(I)
           ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL_R(I)
           ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL_R(I)
           ELEMNT7 = ELEMNT7 + DINT7(IA,IB)*TSHELL_R(I)
        ENDDO
*  
*   Return to the start of the loop
*   

**      ICI = 0
**      DO 21 I = 1,NELMNT
**        IRI = IROW(I)
**         IF (I .GT. IENDC(ICI)) ICI = ICI+1
         DO 22 J = 1,NVEC
            LOC = (J-1)*NCF 
            CONTRI1 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1
            CONTRI2 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT2
            CONTRI3 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT3
            CONTRI4 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT4
            CONTRI5 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT5
            CONTRI6 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT6
            CONTRI7 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT7
            IF (IR.NE.IC) THEN
               CONTRI1 = 2.0D 00 * CONTRI1
               CONTRI2 = 2.0D 00 * CONTRI2
               CONTRI3 = 2.0D 00 * CONTRI3
               CONTRI4 = 2.0D 00 * CONTRI4
               CONTRI5 = 2.0D 00 * CONTRI5
               CONTRI6 = 2.0D 00 * CONTRI6
               CONTRI7 = 2.0D 00 * CONTRI7
            ENDIF
            DENS1(J) = DENS1(J) + CONTRI1
            DENS2(J) = DENS2(J) + CONTRI2
            DENS3(J) = DENS3(J) + CONTRI3
            DENS4(J) = DENS4(J) + CONTRI4
            DENS5(J) = DENS5(J) + CONTRI5
            DENS6(J) = DENS6(J) + CONTRI6
            DENS7(J) = DENS7(J) + CONTRI7
   22    CONTINUE
**   21 CONTINUE

      GOTO 16
      ENDIF
      END
