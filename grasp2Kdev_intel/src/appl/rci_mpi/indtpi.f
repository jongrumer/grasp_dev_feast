************************************************************************
*                                                                      *
      FUNCTION INDTPI (ITYPE,I)
*                                                                      *
*   Written by Farid A Parpia             Last revision: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6
      POINTER (PVALT1,VALT1DUMMY)
      POINTER (PVALT2,VALT2DUMMY)
      POINTER (PVALT3,VALT3DUMMY)
      POINTER (PVALT4,VALT4DUMMY)
      POINTER (PVALT5,VALT5DUMMY)
      POINTER (PVALT6,VALT6DUMMY)
      LOGICAL FIRST
*
      POINTER (PINDT1,INDTP1(*))
      POINTER (PINDT2,INDTP2(*))
      POINTER (PINDT3,INDTP3(*))
      POINTER (PINDT4,INDTP4(*))
      POINTER (PINDT5,INDTP5(*))
      POINTER (PINDT6,INDTP6(*))
*
      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
*
      IF     (ITYPE .EQ. 1) THEN
         INDTPI = INDTP1(I)
      ELSEIF (ITYPE .EQ. 2) THEN
         INDTPI = INDTP2(I)
      ELSEIF (ITYPE .EQ. 3) THEN
         INDTPI = INDTP3(I)
      ELSEIF (ITYPE .EQ. 4) THEN
         INDTPI = INDTP4(I)
      ELSEIF (ITYPE .EQ. 5) THEN
         INDTPI = INDTP5(I)
      ELSEIF (ITYPE .EQ. 6) THEN
         INDTPI = INDTP6(I)
      ENDIF
*
      RETURN
      END
