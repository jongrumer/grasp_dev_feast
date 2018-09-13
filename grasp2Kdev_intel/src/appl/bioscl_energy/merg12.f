************************************************************************
*                                                                      *
      SUBROUTINE MERG12 (NAME,NCORER,NCORE)
*                                                                      *
*   This subroutines merges the initial and final state lists          *
*   Observ that there may doublets in this list if the initial         *
*   and final states have the same parity                              *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CFF      PARAMETER (NNNW = 120)      

Cww      INTEGER PNTIQR,PNTRIQ
      POINTER (PNTIQR,IQRDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)  
      CHARACTER*24 NAME(2)
      CHARACTER*500 LINE
      CHARACTER*2 NH, NHR, NEWNH(NNNW),NEW3
      INTEGER NEWNP(0:NNNW+1),NEWNAK(NNNW),ICON(NNNW),ICON1(NNNW)
      
*
      EXTERNAL IQ,IQR,ISPAR,ISPARR,ITJPO,ITJPOR,JCUP,JCUPR,JQS,JQSR
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF1R/NELECR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB2R/NCFR,NWR,PNTIQR
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB4R/NPR(NNNW),NAKR(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORB10R/NHR(NNNW)
* NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
* NCFF(I): the end position of the Ith block for the final states in the globle CSF list
      common /IBLK/NBLOCKI,NCFI(10)
      common /FBLK/NBLOCKF,NCFF(10)

*
      OPEN (UNIT = 21,FILE='SLASK',FORM='FORMATTED',STATUS='UNKNOWN')
*
*   The same number of electrons must appear in both lists
*
      IF (NELECR .NE. NELEC) THEN
        PRINT *, 'The number of electrons is not equal in the'
        PRINT *, ' first and second GRASP92 Configuration'
        PRINT *, ' symmetry list files.'
        STOP
      ENDIF
*
*   THe core orbitals must be the same
*
      IF (NCORE.NE.NCORER) THEN
        PRINT *, 'The number of core orbitals must be the same'
        STOP
      ENDIF

      DO I = 1,NCORE
        IF (NP(I).NE.NPR(I).OR.NAK(I).NE.NAKR(I)) THEN
          PRINT *, 'The core orbitals must be the same'
          STOP
        ENDIF
      ENDDO

      DO I = 1,NW
        ICON(I) = 0
      ENDDO
*
*   For each orbital in the initial list check if there is a corresponding
*   orbital in the final list. If so give the number of the corresponding
*   orbital
*
      DO I = 1,NW
        DO J = 1,NWR
          IF (NP(I).EQ.NPR(J).AND.NAK(I).EQ.NAKR(J)) THEN
            ICON(I) = J
          ENDIF
        ENDDO
      ENDDO
*
*   Check if the ordering of the initial and final state orbitals
*   is consistent. The condition for this is that ICON is increasing
*
      J = 0
      DO I = 1,NW
        IF (ICON(I).NE.0) THEN
          J = J +1
          ICON1(J) = ICON(I)
        ENDIF
      ENDDO

      ICOMP = ICON1(1)
      DO I = 2,J
        IF (ICON1(I).LT.ICOMP) THEN
          WRITE(*,*) ' In merg12: ordering of the initial and final'
          WRITE(*,*) ' state orbitals is inconsistent. STOP'
          STOP
        ELSE
          ICOMP = ICON1(I)
        ENDIF
      ENDDO
*
*   Determine a common orbital set for the initial and final state.
*   The common set must be such that the order of both the initial and
*   final state sets are preserved. 
*
      DO I = 1,NWR
        NEWNP(I) = NPR(I)
        NEWNAK(I) = NAKR(I)
        NEWNH(I) = NHR(I)
        write(*,*) NEWNP(I),NEWNH(I)
      ENDDO
*
*   Add the initial state orbitals at the end
*
      NEW = NWR
      DO I = 1,NW
        IF (ICON(I).EQ.0) THEN
          NEW = NEW + 1
          NEWNP(NEW) = NP(I)
          NEWNAK(NEW) = NAK(I)        
          NEWNH(NEW) = NH(I)
        write(*,*) NEWNP(NEW),NEWNH(NEW),NEW
        ENDIF
      ENDDO
*
*   Now sort in the orbitals at the end in the right position
*
      DO I = NWR+1,NEW
*
*   Position in initial state list
*
        DO J = 1,NW
          IF (NEWNP(I).EQ.NP(J).AND.NEWNAK(I).EQ.NAK(J)) IPI1 = J
        ENDDO
        write(*,*) 'i,ipi1',i,ipi1
      
        DO J = 1,NWR
          IPI2 = 0
          DO K = 1,NW
            IF (NEWNP(J).EQ.NP(K).AND.NEWNAK(J).EQ.NAK(K)) IPI2 = K
          ENDDO 
          write(*,*) 'j,ipi2',j,ipi2
          IF (IPI2.NE.0) THEN
            IF (IPI1.LT.IPI2) THEN
              NEW1 = NEWNP(I)
              NEW2 = NEWNAK(I)
              NEW3 = NEWNH(I) 
              DO K = I,J+1,-1
                NEWNP(K) = NEWNP(K-1)
                NEWNAK(K) = NEWNAK(K-1) 
                NEWNH(K) = NEWNH(K-1)
              ENDDO
              NEWNP(J) = NEW1
              NEWNAK(J) = NEW2
              NEWNH(J) = NEW3  
              DO M = 1,NEW
                write(*,*) NEWNP(M),NEWNH(M)
              ENDDO
              GOTO 193
            ENDIF
          ENDIF
        ENDDO
  193   CONTINUE
      ENDDO    

      NW = NEW
      DO I = 1,NW
        NP(I) = NEWNP(I)
        NAK(I) = NEWNAK(I)
        NH(I) = NEWNH(I)
      ENDDO
*
*   Determine NKL and NKJ
*
      DO I = 1,NW
        NKJ(I) = 2*ABS(NAK(I))-1
        IF (NAK(I).GT.0) THEN
          NKL(I) = (NKJ(I)+1)/2
        ELSE
          NKL(I) = (NKJ(I)-1)/2
        ENDIF
      ENDDO
*
*   Determine the common core subshells; write out the list;
*   determine the pell subshells; write out the list; these
*   data form the first part of the header of the .csl file;
*   one additional line forms the remainder of the header of
*   the  .csl file
*
      WRITE (21,'(A)') 'Core subshells:'
      WRITE (21,301) (NP(I),NH(I),I = 1,NCORE)
      WRITE (21,'(A)') 'Peel subshells:'
      WRITE (21,301) (NP(I),NH(I),I = NCORE+1,NW)
      WRITE (21,'(A)') 'CSF(s):'
*
*   Now write out all CSFs in the initial and final state list
*
      J = INDEX(NAME(1),' ')                                               
      OPEN (UNIT = 23,FILE=NAME(1)(1:J-1)//'.c',FORM='FORMATTED',          
     :     STATUS='OLD')

      DO I = 1,5
        READ(23,'(A)') LINE
      ENDDO

      NBLOCKI = 0
      NLINE = 0
    5 READ (23,'(A)',END=98) LINE
      if(line(1:2).eq.' *') then
        NBLOCKI = NBLOCKI + 1
        NCFI(NBLOCKI) = NLINE/3
      else
        NLINE = NLINE + 1
      endif
      K = 500           
   10 IF (LINE(K:K) .EQ. ' ') THEN
        K = K-1                    
        IF (K .GT. 1) GOTO 10
      ENDIF                 
      WRITE(21,'(A)') LINE(1:K) 
      GOTO 5
   98 NBLOCKI = NBLOCKI + 1
      NCFI(NBLOCKI) = NLINE/3
      CLOSE (23)

C zou
c     if(NAME(2).EQ.NAME(1)) return
C zou
      J = INDEX(NAME(2),' ')                                
      OPEN (UNIT = 23,FILE=NAME(2)(1:J-1)//'.c',FORM='FORMATTED',
     :     STATUS='OLD')
               
      DO I = 1,5
        READ(23,'(A)') LINE 
      ENDDO                                                      
                  
      NBLOCKF = 0
      NLINE = 0
   15 READ (23,'(A)',END=99) LINE                           
      if(line(1:2).eq.' *') then
        NBLOCKF = NBLOCKF + 1
        NCFF(NBLOCKF) = NLINE/3
      else
        NLINE = NLINE + 1
      endif
      K = 500
   20 IF (LINE(K:K) .EQ. ' ') THEN  
        K = K-1                                             
        IF (K .GT. 1) GOTO 20                              
      ENDIF                                                 
      WRITE(21,'(A)') LINE(1:K)                               
      GOTO 15     
   99 NBLOCKF = NBLOCKF + 1
      NCFF(NBLOCKF) = NLINE/3
      CLOSE (23)

      CLOSE (21)
*

      print *, nblocki
      print *, (ncfi(i),i=1,nblocki)
      print *, nblockf
      print *, (ncff(i),i=1,nblockf)
  301 FORMAT (120(1X,1I2,1A2))
*
      RETURN
      END
