************************************************************************
*                                                                      *
      SUBROUTINE COUPLE (NCORE,J)
*                                                                      *
*   This subroutine obtains all possible couplings of two angular mo-  *
*   menta                                                              *
*                                                                      *
*   Call(s) to: RALC2D.                                                *
*                                                                      *
*   Written by Farid A Parpia and Wasantha P Wijesundera, at Oxford    *
*                                                                      *
*                                           Last update: 27 Aug 1992   *
*                                                                      *
************************************************************************
*
*   The maximum number of states of a subshell is LJLMAX; this
*   value is set to the largest value of any element of the
*   array ITAB in BLOCK DATA TERM
*
      PARAMETER (LJLMAX = 20)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWM1
CGG      PARAMETER (NNNWM1 = 119)
CGG      INTEGER NNNWM2
CGG      PARAMETER (NNNWM2 = 118)
      POINTER (PJCLST,JCLIST(NNNWM1,*))
*
      COMMON/COUBOX/PJCLST,MNJVC,LJCL(NNNWM1),ICPTR(NNNWM2)
     :      /ORBBOX/JVLIST(NNNW,LJLMAX),JWLIST(NNNW,LJLMAX),
     :              JLIST(NNNW,LJLMAX),LJL(NNNW),IOPTR(NNNW)
*
*   Perform the coupling depending on the value of J
*
      IF (J .EQ. 1) THEN
         J1 = JLIST(NCORE+1,IOPTR(NCORE+1))
         J2 = JLIST(NCORE+2,IOPTR(NCORE+2))
         JMIN = ABS (J1-J2)
         JMAX = J1+J2
         IEL = 0
         DO 1 I = JMIN,JMAX,2
            IEL = IEL+1
            IF (IEL .GT. MNJVC) THEN
               NEWSIZ = MNJVC+MNJVC/2
               CALL RALC2D (PJCLST,NNNWM1,MNJVC,NNNWM1,NEWSIZ,4)
               MNJVC = NEWSIZ
            ENDIF
            JCLIST(1,IEL) = I
    1    CONTINUE
         LJCL(1) = IEL
      ELSE
         JP1 = J+1
         JM1 = J-1
         J1 = JLIST(NCORE+JP1,IOPTR(NCORE+JP1))
         J2 = JCLIST(JM1,ICPTR(JM1))
         JMIN = ABS (J1-J2)
         JMAX = J1+J2
         IEL = 0
         DO 2 I = JMIN,JMAX,2
            IEL = IEL+1
            IF (IEL .GT. MNJVC) THEN
               NEWSIZ = MNJVC+MNJVC/2
               CALL RALC2D (PJCLST,NNNWM1,MNJVC,NNNWM1,NEWSIZ,4)
               MNJVC = NEWSIZ
            ENDIF
            JCLIST(J,IEL) = I
    2    CONTINUE
         LJCL(J) = IEL
      ENDIF
*
      RETURN
      END
