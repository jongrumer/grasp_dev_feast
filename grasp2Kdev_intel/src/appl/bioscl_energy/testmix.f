************************************************************************
*                                                                      *
      SUBROUTINE TESTMIX
*                                                                      *
*   This routine checks the mixing coefficients                        *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLEPRECISION(A-H, O-Z)

      POINTER(PNEVALII,EVALII(1))
      POINTER(PNEVECII,EVECII(1))
      POINTER(PNIVECII,IVECII(1))
      POINTER(PIATJPII,IATJPOII(1)),(PIASPAII,IASPARII(1))

      POINTER(PNEVALFF,EVALFF(1))
      POINTER(PNEVECFF,EVECFF(1))
      POINTER(PNIVECFF,IVECFF(1))
      POINTER(PIATJPFF,IATJPOFF(1)),(PIASPAFF,IASPARFF(1))

      COMMON/DEF1II/EMNII,IONCTYII,NELECII,ZII
     :      /EIGVALII/EAVII,PNEVALII
     :      /EIGVECII/PNEVECII
     :      /ORB2II/NCFII,NWII
     :      /PRNTII/NVECII,PNIVECII,NVECMXII
     :      /SYMAII/PIATJPII,PIASPAII

      COMMON/DEF1FF/EMNFF,IONCTYFF,NELECFF,ZFF
     :      /EIGVALFF/EAVFF,PNEVALFF
     :      /EIGVECFF/PNEVECFF
     :      /ORB2FF/NCFFF,NWFF
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF
     :      /SYMAFF/PIATJPFF,PIASPAFF

      WRITE(*,*) ' ****************'
      WRITE(*,*) ' Entering testmix'
      WRITE(*,*) ' ****************'
      WRITE(*,*)
      WRITE(*,*) 'Initial state'
      WRITE(*,*) 'EVALII',(EAVII+ EVALII(I),I=1,NVECII)
      WRITE(*,*) 'NELECII,NCFII,NWII,NVECMXII',NELECII,NCFII,NWII,NVECMXII
      WRITE(*,*) NVECII
      WRITE(*,*) (IVECII(I),I=1,NVECII)
      WRITE(*,*) (IATJPOII(I),IASPARII(I),I=1,NVECII)
      WRITE(*,*) ((EVECII(I+(J-1)*NCFII),I=1,NCFII),J=1,NVECII)

      WRITE(*,*) 'Final state'
      WRITE(*,*) 'EVALFF',(EAVFF+EVALFF(I),I=1,NVECFF)
      WRITE(*,*) 'NELECFF,NCFFF,NWFF,NVECMXFF',NELECFF,NCFFF,NWFF,NVECMXFF
      WRITE(*,*) NVECFF
      WRITE(*,*) (IVECFF(I),I=1,NVECFF)
      WRITE(*,*) (IATJPOFF(I),IASPARFF(I),I=1,NVECFF)
      WRITE(*,*) ((EVECFF(I+(J-1)*NCFFF),I=1,NCFFF),J=1,NVECFF)
      WRITE(*,*)
      WRITE(*,*) ' ***************'
      WRITE(*,*) ' Leaving testmix'
      WRITE(*,*) ' ***************'

      RETURN
      END
