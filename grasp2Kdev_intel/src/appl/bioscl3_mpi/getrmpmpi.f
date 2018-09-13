************************************************************************
*                                                                      *
      SUBROUTINE GETRMP
*                                                                      *
*   Interactively  determines the  list of radiation multipolarities   *
*   and parities. This is loadad into COMMON/OSC6/.                    *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL LELEC,LMAGN,GETYN,YES
      CHARACTER*256 RECORD
      CHARACTER*1 RECI
*
      DIMENSION CNUM(3)
*
      POINTER (PNTRKP,KP(*))
*
      COMMON/OFFD/NOFFD1,NOFFD2
      COMMON/OSC6/NKP,PNTRKP
      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr                                           
*
*   Initial allocation for PNTRKP
*
      NDKP = 1
      CALL ALLOC (PNTRKP,NDKP,4)

      IF (myid .EQ. 0) THEN
*
*   Entry message
*
    1     PRINT *, 'Enter the list of transition specifications'
          PRINT *, ' e.g.,  E1,M2  or  E1 M2  or  E1;M2 :'
*
*   Initialise NKP
*
    2     READ (*,'(A)') RECORD
          NKP = 0
*
*   Parse RECORD from left to right
*
          ISTART = 0
          I = 1
    3     RECI = RECORD(I:I)
          IF ((RECI .NE. ' ') .AND.
     :        (RECI .NE. ',') .AND.
     :        (RECI .NE. ';')) THEN
             IF (ISTART .EQ. 0) ISTART = I
          ELSE
             IF (ISTART .NE. 0) THEN
                IEND = I-1
                RECI = RECORD(ISTART:ISTART)
                IF (RECI .EQ. 'E') THEN
                   LELEC = .TRUE.
                   LMAGN = .FALSE.
                ELSEIF (RECI .EQ. 'M') THEN
                   LELEC = .FALSE.
                   LMAGN = .TRUE.
                ELSE
                   PRINT *, 'GETRMP: Transitions must be of type'
                   PRINT *, ' E or type M; reenter ...'
                   GOTO 2
                ENDIF
                LENTH = IEND-ISTART
                IF (LENTH .NE. 1) THEN
                   PRINT *, 'GETRMP: Transition multipolarities'
                   PRINT *, ' must be integers between 1 and 9;'
                   PRINT *, ' reenter ...'
                   GOTO 2
                ENDIF
                RECI = RECORD(IEND:IEND)
                READ (RECI,'(1I1)',IOSTAT = IOS) MULT
                IF (IOS .NE. 0) THEN
                   PRINT *, 'GETRMP: Unable to decode multipolarity'
                   PRINT *, ' '//RECI//'; reenter ...'
                   GOTO 2
                ENDIF
                NKP = NKP+1
                IF (NKP .GT. NDKP) THEN
                   CALL RALLOC (PNTRKP,NDKP,NKP,4)
                   NDKP = NKP
                ENDIF
                IF (LELEC) THEN
                   KP(NKP) = MULT*(-1)**MULT
                ELSEIF (LMAGN) THEN
                   KP(NKP) = MULT*(-1)**(MULT+1)
                ENDIF
                ISTART = 0
             ENDIF
          ENDIF
*
          IF (I .LT. 256) THEN
             I = I+1
             GOTO 3
          ENDIF
*
          IF (NKP .EQ. 0) GOTO 1
       ENDIF !myid=0
       CALL MPI_Bcast (NKP, 1, MPI_INTEGER, 0,
     & MPI_COMM_WORLD, ierr)
*
*   Trim array to the exact size
*
      IF (NDKP .NE. NKP) CALL RALLOC (PNTRKP,NDKP,NKP,4)

      CALL MPI_Bcast (KP(1), NKP, MPI_INTEGER, 0,
     &    MPI_COMM_WORLD, ierr)
      
*   
*   If M1 or E2 inquire if the transitions are between levels
*   with different J quantum numbers.
*     
      IF(myid .EQ. 0) THEN
          DO I = 1,NKP
             IF (KP(I).EQ.1) THEN
                WRITE(*,*)
     :      'M1 transitions only between levels with different J?'      
                YES = GETYN ()
                IF (YES) THEN
                   NOFFD1 = 1
                ELSE
                   NOFFD1 = 0
                ENDIF
             ENDIF
             IF (KP(I).EQ.2) THEN
                WRITE(*,*)
     :      'E2 transitions only between levels with different J?'  
                YES = GETYN ()
                IF (YES) THEN
                   NOFFD2 = 1
                ELSE
                   NOFFD2 = 0
                ENDIF
             ENDIF
          ENDDO
      ENDIF !myid=0
      CALL MPI_Bcast (NOFFD1, 1, MPI_INTEGER, 0,
     &    MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (NOFFD2, 1, MPI_INTEGER, 0,
     &    MPI_COMM_WORLD, ierr)
*
*R    print *,'dbg: myid=', myid," kp=",(KP(I),i=1,nkp)
      RETURN
      END
