************************************************************************
*                                                                      *
      SUBROUTINE CUTOUT(NAME,INPCI)
*                                                                      *
*   Condenses the .csl list by eliminating CSFs that make no contri-   *
*   bution above the threshold THRESH.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
*                        LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*   Updated by Anders Ynnerman            Last revision: 31 Jan 1994   *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

Cww      INTEGER PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(*))
      CHARACTER*6 G92MIX
      CHARACTER*256 FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
      CHARACTER*2 NH
*
      EXTERNAL IQ,ISPAR,ITJPO,JCUP,JQS
*
      POINTER (PINDEX,INDEX(*))
      POINTER (PNEVECT,EVECT(*))
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
      CHARACTER*300 LINE1, LINE2, LINE3
      CHARACTER*500 string(5)
      CHARACTER*24 NAME
      common /BLK/NBLOCK,NCFBLK(30)
      print *, NBLOCK
      print *, (NCFBLK(I),I=1,NBLOCK)
*
*   Determine the cutoff criterion
*
      PRINT*, 'Enter the cut-off value for the coefficients [0--1]'
      READ *, THRESH
*
*   Open the output  .csl file; it is FORMATTED; it must not exist
*
      OPEN (UNIT = 21,FILE=TRIM(NAME)//'_cond.c',FORM='FORMATTED',
     :     STATUS='UNKNOWN')
      OPEN (UNIT = 19,FILE=TRIM(NAME)//'.c',FORM='FORMATTED',
     :     STATUS='UNKNOWN')
*
*   Now write out all CSFs in the reference list
*
    6 DO I=1,5
      READ (19, '(A)') STRING(I)
      WRITE (21, '(A)') TRIM(STRING(I))
      ENDDO


      IF (INPCI.EQ.0) THEN
        OPEN (UNIT = 68,FILE=TRIM(NAME)//'.cm',FORM='UNFORMATTED',
     :        STATUS='OLD')
        OPEN (UNIT = 78,FILE=TRIM(NAME)//'_cond.cm',FORM='UNFORMATTED',
     :        STATUS='UNKNOWN')
      ELSE
        OPEN (UNIT = 68,FILE=TRIM(NAME)//'.m',FORM='UNFORMATTED',
     :        STATUS='OLD')
        OPEN (UNIT = 78,FILE=TRIM(NAME)//'_cond.m',FORM='UNFORMATTED',
     :        STATUS='UNKNOWN')
      ENDIF
      READ(68,IOSTAT=IOS) G92MIX
      IF ((IOS.NE.0).OR.(G92MIX.NE.'G92MIX')) THEN
        WRITE(*,*) 'Not a GRASP mixing file'
        STOP
      ENDIF
      READ (68) nelec, ncftoti, nwi, nvectoti, nvecsizi, nblocki
      WRITE(78) G92MIX
      WRITE(78) nelec, ncftoti, nwi, nvectoti, nvecsizi, nblocki
      write (6,*) '   nelec  = ', nelec
      write (6,*) '   ncftoti = ', ncftoti
      write (6,*) '   nwi     = ', nwi
      write (6,*) '   nblocki = ', nblocki
      write (6,*)

      DO 2001 IBLK = 1,NBLOCK
      CALL GETMIXBLK(IBLK)    
*
*    Establish the list of CSFs
*
      CALL ALLOC (PINDEX,NCF,4)
      NSHORT = 0
      DO 2 I = 1,NCF
         DO 1 J = 1,NVEC
            IF (ABS (EVEC(I+(J-1)*NCF)) .GT. THRESH) THEN
               NSHORT = NSHORT+1
               INDEX(NSHORT) = I
               GOTO 2
            ENDIF
    1    CONTINUE
    2 CONTINUE
      CALL ALLOC (PNEVECT,NVEC*NSHORT,8)

      NCF=0
      DO 7 I = 1,NSHORT
         ICSF = INDEX(I)
   10 READ (19, '(A)') LINE1
      READ (19, '(A)') LINE2
      READ (19, '(A)') LINE3
      NCF=NCF+1
      if(ICSF.EQ.NCF) THEN 
      WRITE(21,'(A)') TRIM(LINE1)
      WRITE(21,'(A)') TRIM(LINE2)
      WRITE(21,'(A)') TRIM(LINE3)
         DO  J = 1,NVEC
            EVECT(I+(J-1)*NSHORT) = EVEC(INDEX(I)+(J-1)*NCF)
         ENDDO   

        GO TO 7
      ELSE
      GOTO 10
      ENDIF
    7 CONTINUE

      NCF=NSHORT
      WRITE(78) iblk,ncf,NVEC,iatjpo(1),iaspar(1)
      WRITE(78) (IVEC(I),I=1,NVEC)
      WRITE(78) EAV,(EVAL(I),I=1,NVEC)
      WRITE(78) ((EVECT(I+(J-1)*NCF),I=1,NCF),J=1,NVEC)

C local
      CALL DALLOC (PINDEX)
      CALL DALLOC (PNEVECT)
C allocted in getmixblk
      CALL DALLOC(PNEVAL)
      CALL DALLOC(PNEVEC)
      CALL DALLOC(PNIVEC)
      CALL DALLOC(PIATJP)
      CALL DALLOC(PIASPA)
  500   IF(LINE1(1:2).EQ.' *') then
          WRITE(21,'(A)') ' *'
          go to 2001
        else
          READ (19, '(A)',END=2001) LINE1
          go to 500
        endif

 2001 CONTINUE
*
      CALL CONVRT (NW,FILNAM,LENTH)
      PRINT *, FILNAM(1:LENTH)//' relativistic subshells;'
      CALL CONVRT (NSHORT,FILNAM,LENTH)
      PRINT *, FILNAM(1:LENTH)//' relativistic CSFs.'
*
      CLOSE (21)
      CLOSE (19)
*
      RETURN
*
  300 FORMAT (120(1X,1I2,1A2))
*
      END
