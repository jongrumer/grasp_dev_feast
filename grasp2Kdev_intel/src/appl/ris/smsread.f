************************************************************************
*                                                                      *
      SUBROUTINE SMSREAD(VINT,VINT2)
*
*   Call(s) to: [LIB92]: ALCBUF                                        *
*                                                                      *
*   WRITTEN  by C. Naz\'e  Oct. 2011                                   *
*                                                                      *
************************************************************************
*     
*

      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
*   Key is used to store the indices of a couple of orbitals
      PARAMETER (KEY = KEYORB)
*   Matrix elements smaller than CUTOFF are not accumulated
      PARAMETER (CUTOFF = 1.0D-10)

      DIMENSION VINT(NNNW,NNNW),VINT2(NNNW,NNNW)

      POINTER (PLABEL,LABEL(6,*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNSMS1,SMSC1(*))
      POINTER (PNSMS2,SMSC2(*))
      POINTER (PNTRIQ,RIQDUMMY)
      
      COMMON/EIGVEC/PNEVEC
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /ORB2/NCF,NW,PNTRIQ
     :      /SMS1/PNSMS1,PNSMS2


      CALL ALCBUF (1) 
      
      REWIND(51)
  16  READ (51,IOSTAT = IOS) IC,IR
!C      WRITE(*,*) IC,IR
      IF (IOS .EQ. 0) THEN
!cc        DO J = 1,NVEC
          READ(51) COEFFSMS,LAB
          IID = MOD (LAB, KEY)
          LAB = LAB/KEY
          IIB  = MOD (LAB, KEY)
          LAB = LAB/KEY
          IIC = MOD (LAB, KEY)
          IIA = LAB/KEY
        DO J = 1,NVEC
          LOC = (J-1)*NCF
            CONTRIK1 = - EVEC(IC+LOC)*EVEC(IR+LOC)
     :         * COEFFSMS
     :         * VINT (IIA,IIC)*VINT(IIB,IID)
            CONTRI = - EVEC(IC+LOC)*EVEC(IR+LOC)
     :         * COEFFSMS
     :         * ( VINT2(IIA,IIC)*VINT(IIB,IID)
     :         + VINT2(IIB,IID)*VINT(IIA,IIC))/2.0D 00
          IF (IR.NE.IC) THEN 
            CONTRI = 2.0D 00 * CONTRI
            CONTRIK1 = 2.0D 00 * CONTRIK1
          ENDIF
          SMSC1(J) = SMSC1(J) + CONTRIK1
          SMSC2(J) = SMSC2(J) + CONTRI
        ENDDO

      GOTO 16
      ENDIF
      CALL ALCBUF (3) 
      END
