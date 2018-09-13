************************************************************************
*                                                                      *
      FUNCTION GCO (K,IR,IA,IB)
*                                                                      *
*   This routine evaluates a coefficient                               *
*                                                                      *
*                                K                                     *
*                               g   (IA,IB)                            *
*                                IR                                    *
*                                                                      *
*                                                                      *
*   Here  K  is the multipolarity, IR  is the sequence number of the   *
*   configuration, and  IA and IB are orbital sequence  numbers. See   *
*   I P Grant,  B J McKenzie,  P H Norrington,  D F  Mayers, and N C   *
*   Pyper, Computer Phys Commun 21 (1980) 207-231, Eq (7).             *
*                                                                      *
*   Call(s) to: [LIB92]: CLRX, IQ.                                     *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL FULLA,FULLB,LDBPA
      CHARACTER*2 NH
*
      COMMON/DEBUGA/LDBPA(5)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
*
      IQA = IQ (IA,IR)
      IQB = IQ (IB,IR)
*
      FULLA = IQA .EQ. NKJ(IA)+1
      FULLB = IQB .EQ. NKJ(IB)+1
*
      IF (FULLA .OR. FULLB) THEN
         QAB = DBLE (IQA*IQB)
         FAC = CLRX (NAK(IA),K,NAK(IB))
         GCO = -QAB*FAC*FAC
      ELSE
         GCO = 0.0D 00
      ENDIF
*
*   Debug printout
*
      IF (LDBPA(3) .AND. (ABS (GCO) .GT. 0.0D 00))
     :      WRITE (99,300) K,NP(IA),NH(IA),NP(IB),NH(IB),GCO,IR
*
      RETURN
*
  300 FORMAT (/'  ',1I2
     :        /' g    (',1I2,1A2,',',1I2,1A2,') = ',1PD21.14,
     :        /'  ',1I3/)
*
      END
