************************************************************************
*                                                                      *
      FUNCTION BRRA (ITYPE,IA,IC,IB,ID,K)
*                                                                      *
*   This routine evaluates the transverse interaction integrals:       *
*                                                                      *
*      ITYPE = 1: General R (k; a c | b d )                            *
*            = 2: General S (k; a c | b d )                            *
*            = 3: R (k; a , b d )                                      *
*            = 4: F (k; a , b )                                        *
*            = 5: G (k; a , b ) type integral                          *
*            = 6: H (k; a , b ) type integral                          *
*                                                                      *
*   Call(s) to: [RCI92]: RKINT, SKINT.                                 *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION RAC(NNNP),RBD(NNNP)
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      MXRAC = MIN (MF(IA),MF(IC))
      DO 1 I = 1,MXRAC
         RAC(I) = PF(I,IA)*QF(I,IC)
    1 CONTINUE
*
      MXRBD = MIN (MF(IB),MF(ID))
      DO 2 I = 1,MXRBD
         RBD(I) = PF(I,IB)*QF(I,ID)
    2 CONTINUE
*
      GOTO (21,22,23,24,25,26),ITYPE
*
*   ITYPE = 1
*
   21 CONTINUE
      IF ((IA .EQ. IB) .AND. (IC .EQ. ID)) GOTO 9
      IF ((IA .EQ. ID) .AND. (IC .EQ. IB)) GOTO 10
      BRRA = ( RKINT (RAC,IA,IC,RBD,IB,ID,K,1)
     :        +RKINT (RAC,IA,IC,RBD,IB,ID,K,2)
     :        +RKINT (RBD,IB,ID,RAC,IA,IC,K,1)
     :        +RKINT (RBD,IB,ID,RAC,IA,IC,K,2))*HALF
      RETURN
*
*   ITYPE = 2
*
   22 CONTINUE
      IF ((IA .EQ. IB) .AND. (IC .EQ. ID)) GOTO 26
      IF ((IA .EQ. ID) .AND. (IC .EQ. IB)) GOTO 26
      BRRA = ( SKINT (RAC,IA,IC,RBD,IB,ID,K,1)
     :        +SKINT (RAC,IA,IC,RBD,IB,ID,K,2))*HALF
      RETURN
*
*   ITYPE = 3
*
   23 CONTINUE
      IF (IA .NE. IC) GOTO 5
      DO 4 I = 1,MXRBD
         RBD(I) = RBD(I)+PF(I,ID)*QF(I,IB)
    4 CONTINUE
      BRRA = ( RKINT (RAC,IA,IC,RBD,IB,ID,K,0)
     :        +RKINT (RBD,IB,ID,RAC,IA,IC,K,0)
     :        +RKINT (RAC,IA,IC,RBD,IB,ID,K,2)
     :        +RKINT (RBD,IB,ID,RAC,IA,IC,K,2) )*HALF
      RETURN
    5 CONTINUE
      DO 6 I = 1,MXRAC
         RAC(I) = RAC(I)+PF(I,IC)*QF(I,IA)
    6 CONTINUE
      BRRA = ( RKINT (RAC,IA,IC,RBD,IB,ID,K,1)
     :        +RKINT (RBD,IB,ID,RAC,IA,IC,K,1)
     :        +RKINT (RAC,IA,IC,RBD,IB,ID,K,0)
     :        +RKINT (RBD,IB,ID,RAC,IA,IC,K,0) )*HALF
      RETURN
*
*   ITYPE = 4
*
   24 CONTINUE
      BRRA =  RKINT (RAC,IA,IC,RBD,IB,ID,K,0)
     :       +RKINT (RBD,IB,ID,RAC,IA,IC,K,0)
      RETURN
*
*   ITYPE = 5
*
   25 CONTINUE
      IF ((IA .EQ. ID) .AND. (IC .EQ. IB)) GOTO 10
    9 BRRA = TWO*RKINT (RAC,IA,IC,RBD,IB,ID,K,1)
      RETURN
   10 BRRA =  RKINT (RAC,IA,IC,RBD,IB,ID,K,1)
     :       +RKINT (RBD,IB,ID,RAC,IA,IC,K,1)
      RETURN
*
*   ITYPE = 6
*
   26 CONTINUE
      BRRA = SKINT (RAC,IA,IC,RBD,IB,ID,K,1)
      RETURN
*
      END
