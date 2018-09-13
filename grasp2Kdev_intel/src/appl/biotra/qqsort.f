************************************************************************
*                                                                      *
      SUBROUTINE QQSORT (NFILE,NUMBER,KSTART,NAME,KAMAX)
*                                                                      *
*     The list of unique integrals (j,i) is formed in the order of     *
*     increasing symmetry, i.e. with j .le. i.                         *
*     With each integral (INTGRL) there is a pointer (INTPTR)          *
*     that points to the last element in the array of coefficients     *
*     for that integral.                                               *
*                                                                      *
*     On Exit                                                          *
*                                                                      *
*     INTGRL -- array of integrals (packed form of j,i)                *
*     NINTG  -- number of integrals                                    *
*     INTPTR -- array of pointers to last element in list of coeff.    *
*     CNN    -- array of coefficients                                  *
*     NCOEFF -- number of coefficients                                 *
*     JANN   -- row of matrix                                          *
*     JBNN   -- column of matrix                                       *
*                                                                      *
*     Written by Per Jonsson                                           *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (KEY = 121)
CGG begin
      PARAMETER (KEY = KEYORB)
CGG end
      PARAMETER (NF = 200)

      CHARACTER*2 NH
      CHARACTER*24 NAME
*
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW) 
     :      /ORB2/NCF,NW,PNTRIQ
     :      /DEFAULT/NDEF,NDUMP
      COMMON/MCPDATA/PJANN,PJBNN,PPPINT,PCNN,PINTPTR,NCOEFF,NINTG

      POINTER (PJANN,JANN(*)),(PPPINT,INTGRL(*)),
     :        (PJBNN,JBNN(*)),(PCNN,CNN(*)),(PINTPTR,INTPTR(*))
*
      REWIND (NFILE+80)

      NCOEFF = NUMBER
*
*   Sort the list
*
*   Allocate storage for all required arrays these arrays are then deallocated
*   in mcp
*
      CALL ALLOC (PJANN,NCOEFF,4)
      CALL ALLOC (PJBNN,NCOEFF,4)
      CALL ALLOC (PPPINT,NCOEFF,4)
      CALL ALLOC (PCNN,NCOEFF,8)
      CALL ALLOC (PINTPTR,NCOEFF,4)
*
*   Read arrays into memory from NFILE
*
      DO I = 1,NCOEFF
        READ (NFILE+80) JANN(I),JBNN(I),INTGRL(I),CNN(I)
      ENDDO
*
*   Sort INTGRL into ascending order using the heapsort algorithm;
*   (Numerical recepies page 231.) move the associated members of
*   of JANN and JBNN in the same
*   manner;
*
      IF (NCOEFF .GT. 1) THEN
*
        L = NCOEFF/2+1
        IR = NCOEFF
    2   IF (L .GT. 1) THEN
          L = L-1
          JA = JANN(L)
          JB = JBNN(L)
          INT = INTGRL(L)
          CN = CNN(L)
        ELSE
          JA = JANN(IR)
          JB = JBNN(IR)
          INT = INTGRL(IR)
          CN = CNN(IR)
          JANN(IR) = JANN(1)
          JBNN(IR) = JBNN(1)
          INTGRL(IR) = INTGRL(1)
          CNN(IR) = CNN(1)
          IR = IR-1
          IF (IR .EQ. 1) THEN
            JANN(1) = JA
            JBNN(1) = JB
            INTGRL(1) = INT
            CNN(1) = CN
            GOTO 4
          ENDIF
        ENDIF
        I = L
        J = L+L
    3   IF (J .LE. IR) THEN
          IF (J .LT. IR) THEN
            IF (INTGRL(J) .LT. INTGRL(J+1)) J = J+1
          ENDIF
          IF (INT .LT. INTGRL(J)) THEN
            JANN(I) = JANN(J)
            JBNN(I) = JBNN(J)
            INTGRL(I) = INTGRL(J)
            CNN(I) = CNN(J)
            I = J
            J = J+J
          ELSE
            J = IR+1
          ENDIF
          GOTO 3
        ENDIF
        JANN(I) = JA
        JBNN(I) = JB
        INTGRL(I) = INT
        CNN(I) = CN
        GOTO 2
      ENDIF
*

*   Sorting complete; close the file
*
Cww    4 CLOSE (80+NFILE)
    4 CONTINUE
*
*  Find the number of unique integrals and form the list of
*  pointers to the data.
*
      NINTG = 1
      INTT = INTGRL(1)
*
      DO I = 1,NCOEFF
        IF (INTGRL(I).NE.INTT) THEN
          INTPTR(NINTG) = I - 1
          NINTG = NINTG + 1
          INTT = INTGRL(I)
        ENDIF
      ENDDO

      INTPTR(NINTG) = NCOEFF
*
      DO I = 1,NINTG
        INTGRL(I) = INTGRL(INTPTR(I))
      ENDDO
*
*  If output option is set dump the data on file
*
      IF (NDUMP.EQ.1) THEN
*
*  If first set of data open the file and print
*  some data to later be able to identify the file
*
c        IF (NFILE.EQ.1.AND.IBLK.EQ.1) THEN
c         print *, ' open sorted ang. file .TB(NF)', NF
c          J = INDEX(NAME,' ')
c          OPEN (UNIT = NF,FILE=NAME(1:J-1)//'.TB',
c     :          STATUS='UNKNOWN',FORM='UNFORMATTED')
c          REWIND (NF)
c        ENDIF
        IF (NFILE.EQ.1) WRITE(NF) NCF,NW,KAMAX
*
*  Print out angular data for this kappa
*
        WRITE(NF) NINTG,NCOEFF
        DO I = 1,NINTG
          WRITE(NF) INTGRL(I),INTPTR(I)
        ENDDO
        DO I = 1,NCOEFF
          WRITE(NF) CNN(I),JANN(I),JBNN(I)
        ENDDO
      ENDIF
*
*  Has all data been processed? If so close
*  the file
*
c     IF (NFILE.EQ.KAMAX) CLOSE (NF)

*  Debug output
*
C      IF (KSTART.EQ.1) THEN
C        NNNN = 120
C      ELSE
C        NNNN = 140
C      ENDIF
C      WRITE(NNNN+NFILE,*) ' Data for '
C      WRITE(NNNN+NFILE,*) NINTG, ' Integrals'
C      DO I = 1,NINTG
C        IA = INTGRL(I)/KEY
C        IB = MOD(INTGRL(I),KEY)
C        WRITE(NNNN+NFILE,302) NP(IA),NH(IA),NP(IB),NH(IB),INTPTR(I)
C      ENDDO
C      WRITE(NNNN+NFILE,*) NCOEFF,' Coefficients'
C      WRITE(NNNN+NFILE,'(F12.8,2I6)')
C     :     (CNN(I),JANN(I),JBNN(I),I=1,NCOEFF)
*
      RETURN
*
  301 FORMAT (' T_[',1I2,',',1I2,']',
     :   ' (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
  302 FORMAT ( ' (',1I2,1A2,',',1I2,1A2,')',I6)
*
      END
