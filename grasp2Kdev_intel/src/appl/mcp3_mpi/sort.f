************************************************************************
*                                                                      *
      SUBROUTINE SORT (NFILE,NCOEFF,NTGRAL,LPRINT,nb, fhead)
      IMPLICIT REAL*8          (A-H, O-Z)
*                                                                      *
*   This routine sorts lists                                           *
*                                                                      *
*                     (ICLMN,INDEX,LABEL,COEFF)                        *
*   into lists                                                         *
*                     (LABEL,ICLMN,INDEX,COEFF)                        *
*                                                                      *
*   using Heapsort. File NFILE is closed by this routine.  NCOEFF is   *
*   the number of triads (INDEX, ...) and is an input. NTGRAL is the   *
*   number of different values of LABEL, and is an output. If LPRINT   *
*   is  .TRUE. , the contents of the sorted file are interpreted and   *
*   printed to unit 99.                                                *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL   LPRINT
      CHARACTER fhead*(*)
      CHARACTER (LEN = LEN (fhead) + 3) fullname
      CHARACTER CNUM*20, SRTLAB*8, MCPLAB*3, CK*2, NH*2
*
      POINTER (PCOEFF,COEFF(*)), (PICLMN,ICLMN(*)), (PINDEX,INDEX(*)),
     : (PLABEL,LABEL(*)), (PNSWAP,NSWAP(*))
*
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB10/NH(NNNW)

      COMMON/iounit/istdi,istdo,istde
*
*   This is the encoding key
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
*   Message
*
      CALL CONVRT (NCOEFF,CNUM,LCNUM)
      IF (nfile .GT. 31) THEN
         CALL CONVRT (NFILE-32,CK,LCK)
         PRINT *, 'Sorting '//CNUM(1:LCNUM)
     :       //' V(k='//CK(1:LCK)//') coefficients ...',nfile
      ELSE
         PRINT *, 'Sorting '//CNUM(1:LCNUM)//' T coefficients ...',nfile
      ENDIF
*
      CALL CONVRT (nfile, ck, lck)
      IF (lck .GT. 2) THEN
         WRITE (istde,*) 'sort: nfile > 99; check fullname'
         STOP
      ENDIF

      fullname = fhead//'.'//CK(1:2)

      OPEN (29, FILE = fullname, STATUS = 'OLD',
     &      FORM = 'UNFORMATTED', IOSTAT = IERR, POSITION = 'APPEND')
      IF (IERR .NE. 0) THEN
         WRITE (istde, *) ' Error when opening the file ',fullname
         STOP
      ENDIF

*     IF (NFILE.EQ.33) REWIND (20)

      READ (NFILE) MCPLAB, mb
      IF (nb .NE. mb ) THEN
         WRITE (istde,*) 'sort: nb = ', nb, '.NE. mb (=', mb,')'
         STOP
      ENDIF
*
*   Sort the list
*
      IF (NCOEFF .GT. 0) THEN
*
*   Allocate storage for all required arrays
*
         CALL ALLOC (PCOEFF,NCOEFF,8)
         CALL ALLOC (PICLMN,NCOEFF,4)
         CALL ALLOC (PINDEX,NCOEFF,4)
         CALL ALLOC (PLABEL,NCOEFF,4)
         IF (NFILE .EQ. 33) CALL ALLOC (PNSWAP,NCOEFF,4)
*
*   Read arrays into memory from NFILE
*
         DO I = 1, NCOEFF
            READ (NFILE) ICLMN(I), INDEX(I), LABEL(I), COEFF(I)
         ENDDO

*         IF (NFILE.EQ.33) THEN
*            DO 11 I = 1,NCOEFF
*               READ (20) NSWAP(I)
*   11       CONTINUE
*         ENDIF                
*
      ENDIF
*
*   Sort LABEL into ascending order using the heapsort algorithm;
*   move the associated members of COEFF and INDEX in the same
*   manner; the code below is adapted from Press et al.
*
      IF (NFILE .EQ. 33 .AND. NCOEFF .GT. 1) THEN

         L = NCOEFF/2 + 1
         IR = NCOEFF
  234    IF (L .GT. 1) THEN
            L = L - 1
            COF = COEFF(L)
            ICL = ICLMN(L)
            IND = INDEX(L)
            LAB = LABEL(L)
            NSW = NSWAP(L)
         ELSE
            COF = COEFF(IR)
            ICL = ICLMN(IR)
            IND = INDEX(IR)
            LAB = LABEL(IR)
            NSW = NSWAP(IR)
            COEFF(IR) = COEFF(1)
            ICLMN(IR) = ICLMN(1)
            INDEX(IR) = INDEX(1)
            LABEL(IR) = LABEL(1)
            NSWAP(IR) = NSWAP(1)
            IR = IR-1
            IF (IR .EQ. 1) THEN
               COEFF(1) = COF
               ICLMN(1) = ICL
               INDEX(1) = IND
               LABEL(1) = LAB
               NSWAP(1) = NSW
               GOTO 456
            ENDIF
         ENDIF
         I = L
         J = L + L
  345    IF (J .LE. IR) THEN
            IF (J .LT. IR) THEN
               IF (LABEL(J) .LT. LABEL(J+1)) J = J + 1
            ENDIF
            IF (LAB .LT. LABEL(J)) THEN
               COEFF(I) = COEFF(J)
               ICLMN(I) = ICLMN(J)
               INDEX(I) = INDEX(J)
               LABEL(I) = LABEL(J)
               NSWAP(I) = NSWAP(J)
               I = J
               J = J + J
            ELSE
               J = IR + 1
            ENDIF
            GOTO 345
         ENDIF
         COEFF(I) = COF
         ICLMN(I) = ICL
         INDEX(I) = IND
         LABEL(I) = LAB
         NSWAP(I) = NSW
         GOTO 234

      ELSE IF (NFILE .NE. 33 .AND. NCOEFF .GT. 1) THEN
*
*   Sort LABEL into ascending order using the heapsort algorithm;
*   move the associated members of COEFF and INDEX in the same
*   manner; the code below is adapted from Press et al.
*
         L = NCOEFF/2 + 1
         IR = NCOEFF
   92    IF (L .GT. 1) THEN
            L = L - 1
            COF = COEFF(L)
            ICL = ICLMN(L)
            IND = INDEX(L)
            LAB = LABEL(L)
         ELSE
            COF = COEFF(IR)
            ICL = ICLMN(IR)
            IND = INDEX(IR)
            LAB = LABEL(IR)
            COEFF(IR) = COEFF(1)
            ICLMN(IR) = ICLMN(1)
            INDEX(IR) = INDEX(1)
            LABEL(IR) = LABEL(1)
            IR = IR - 1
            IF (IR .EQ. 1) THEN
               COEFF(1) = COF
               ICLMN(1) = ICL
               INDEX(1) = IND
               LABEL(1) = LAB
               GOTO 456
            ENDIF
         ENDIF
         I = L
         J = L + L
   93    IF (J .LE. IR) THEN
            IF (J .LT. IR) THEN
               IF (LABEL(J) .LT. LABEL(J+1)) J = J + 1
            ENDIF
            IF (LAB .LT. LABEL(J)) THEN
               COEFF(I) = COEFF(J)
               ICLMN(I) = ICLMN(J)
               INDEX(I) = INDEX(J)
               LABEL(I) = LABEL(J)
               I = J
               J = J + J
            ELSE
               J = IR + 1
            ENDIF
            GOTO 93
         ENDIF
         COEFF(I) = COF
         ICLMN(I) = ICL
         INDEX(I) = IND
         LABEL(I) = LAB
         GOTO 92

      ENDIF
*
*   Sorting complete; rewrite the file header
*
*
  456 WRITE (29) 'MCP',  nb, ncf, ncoeff
*
*   Write the sorted list to mcp.xx
*
      IF (NCOEFF .GT. 0) THEN
*
         LAST = LABEL(1)
         IBEG = 1
         IEND = 1
         NTGRAL = 1
*
         DO I = 2, NCOEFF
            IF (LABEL(I) .EQ. LAST) THEN
               IEND = IEND + 1
            ELSE
               WRITE (29) LAST, IEND - IBEG + 1
               WRITE (29)
     :            (ICLMN(J), INDEX(J), COEFF(J), J = IBEG, IEND)

*              IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,IEND)
               NTGRAL = NTGRAL + 1
               LAST = LABEL(I)
               IBEG = IEND + 1
               IEND = IBEG
            ENDIF
         ENDDO
*
         IF (IBEG .LE. NCOEFF) THEN
            WRITE (29) LAST, NCOEFF - IBEG + 1
            WRITE (29) (ICLMN(J), INDEX(J), COEFF(J), J = IBEG, NCOEFF)

*           IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,NCOEFF)
         ENDIF
*
      ELSE
*
         NTGRAL = 0
*
      ENDIF
*
*     write the terminator record for this block
* 
      WRITE (29) 0, 0
      CLOSE (29)
*
*   Completion message
*
      !PRINT *, ' ... sort complete; ', ntgral, ' integrals;'
*
*   Debug printout
*
      IF (LPRINT) THEN
         WRITE (99,300)
         WRITE (6,300)
         IF (NCOEFF .GT. 0) THEN
*
            LAST = LABEL(1)
            IBEG = 1
            IEND = 1
*
            DO I = 2, NCOEFF
               IF (LABEL(I) .EQ. LAST ) THEN
                  IEND = IEND + 1
               ENDIF
  567          IF (LABEL(I) .NE. LAST .OR. I  .EQ. NCOEFF) THEN
                  LAB = LAST
                  NCONTR = IEND - IBEG + 1
                  IF (NFILE .EQ. 31) THEN
                     IA = MOD (LAB, KEY)
                     IB = LAB/KEY
                     WRITE (99,301) NP(IA), NH(IA), NP(IB), NH(IB)
                     DO J = IBEG, IEND
                        WRITE (99,302) ICLMN(J), INDEX(J), COEFF(J)
                     ENDDO
                  ELSE
                     k=NFILE - 32
                     ID = MOD (LAB, KEY)
                     LAB = LAB/KEY
                     IB  = MOD (LAB, KEY)
                     LAB = LAB/KEY
                     IC = MOD (LAB, KEY)
                     IA = LAB/KEY
                     WRITE (99,304) K, NP(IA), NH(IA), NP(IB), NH(IB),
     :                          NP(IC), NH(IC), NP(ID), NH(ID)
                     DO J = IBEG, IEND
                        WRITE (99,305) K, ICLMN(J), INDEX(J), COEFF(J)
                     ENDDO
                  ENDIF
                  LAST = LABEL(I)
                  IBEG = IEND + 1
                  IEND = IBEG
	               IF (IEND .EQ. NCOEFF) GOTO 567
               ENDIF
            ENDDO
         ENDIF
         WRITE (99,303) NTGRAL
      ENDIF
*
*   Deallocate storage
*
      IF (NCOEFF .GT. 0) THEN
         CALL DALLOC (PCOEFF)
         CALL DALLOC (PICLMN)
         CALL DALLOC (PINDEX)
         CALL DALLOC (PLABEL)
         IF (NFILE .EQ. 33) CALL DALLOC (PNSWAP)
      ENDIF

  300 FORMAT (/'From SORT:')
  301 FORMAT (' I(',1I2,1A2,',',1I2,1A2,'):')
  302 FORMAT ('  T_[',1I2,',',1I4,'] = ',1PD19.12)
  303 FORMAT ('  Number of integrals is ',1I4)
  304 FORMAT (' R^[(',1I2,')] (',1I2,1A2,',',1I2,1A2,';'
     :                          ,1I2,1A2,',',1I2,1A2,'):')
  305 FORMAT ('  V^[(',1I2,')]_[',1I2,',',1I4,'] = ',1PD19.12)

      RETURN
      END
