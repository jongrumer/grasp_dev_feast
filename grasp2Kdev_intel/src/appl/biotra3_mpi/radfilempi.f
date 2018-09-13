************************************************************************
*                                                                      *
      SUBROUTINE RADFILE(NAME)
*                                                                      *
*  This subroutine outputs the transformed radial orbitals             *
*                                                                      *
*  Written by Per Jonsson                                              *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590) 
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      CHARACTER*2 NHII,NHFF
      CHARACTER*128 NAME(2)
*  
*  Orbitals
*
      POINTER (PNTRPFII,PFII(NNNP,*)),(PNTRQFII,QFII(NNNP,*))

      POINTER (PNTRPFFF,PFFF(NNNP,*)),(PNTRQFFF,QFFF(NNNP,*))
*
*  Initial state commons
*
      COMMON/ORB1II/EII(NNNW),GAMAII(NNNW)
     :      /ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB10II/NHII(NNNW)
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*
*  Final state commons
*
      COMMON/ORB1FF/EFF(NNNW),GAMAFF(NNNW)
     :      /ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N

      J = INDEX(NAME(1),' ')
      OPEN (UNIT = 30,FILE=NAME(1)(1:J-1)//'.bw',FORM='UNFORMATTED',
     :     STATUS='UNKNOWN')

      WRITE(30) 'G92RWF'
      write(*,*) 'NWII',NWII
      DO K = 1,NWII
        WRITE(30) NPII(K),NAKII(K),EII(K),MFII(K)
        WRITE(30) PZII(K),
     :            (PFII(I,K),I=1,MFII(K)),(QFII(I,K),I=1,MFII(K))
        WRITE(30) (R(I),I = 1,MFII(K))
      ENDDO 

      CLOSE(30)

      J = INDEX(NAME(2),' ')
      OPEN (UNIT = 30,FILE=NAME(2)(1:J-1)//'.bw',FORM='UNFORMATTED',
     :     STATUS='UNKNOWN')

      WRITE(30) 'G92RWF'
      write(*,*) 'NWFF',NWFF
      DO K = 1,NWFF
        WRITE(30) NPFF(K),NAKFF(K),EFF(K),MFFF(K)
        WRITE(30) PZFF(K),
     :            (PFFF(I,K),I=1,MFFF(K)),(QFFF(I,K),I=1,MFFF(K))
        WRITE(30) (R(I),I = 1,MFFF(K))
      ENDDO

      CLOSE(30) 
     
      RETURN
      END
