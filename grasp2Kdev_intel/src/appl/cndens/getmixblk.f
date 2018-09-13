************************************************************************
*                                                                      *
      SUBROUTINE GETMIXBLK(IBLK)
*                                                                      *
*   Open and read the mixing coefficent files                          *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

Cww      INTEGER PNTRIQ, PNTIQR
      POINTER(PNTRIQ,RIQDUMMY)
      POINTER(PNTIQR,IQRDUMMY)

      CHARACTER*24 NAME
      CHARACTER*6 G92MIX

      POINTER(PNEVAL,EVAL(*))
      POINTER(PNEVEC,EVEC(*))
      POINTER(PNIVEC,IVEC(*))
      POINTER(PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA


*
*
      READ(68) ib,ncf,NVEC,iatjp,iaspa
      if(ib.ne.iblk) stop 'ERROR WHEN READ MIXING COEF.'
      CALL ALLOC(PNEVAL,NVEC,8)
      CALL ALLOC(PNEVEC,NCF*NVEC,8)
      CALL ALLOC(PNIVEC,NVEC,4)
      CALL ALLOC(PIATJP,NVEC,4)
      CALL ALLOC(PIASPA,NVEC,4)
      READ(68) (IVEC(I),I=1,NVEC)
      do i=1,nvec
        IATJPO(I) = iatjp
        IASPAR(I) = iaspa
      enddo
      READ(68) EAV,(EVAL(I),I=1,NVEC)

      READ(68) ((EVEC(I+(J-1)*NCF),I=1,NCF),J=1,NVEC)
*
      RETURN
      END
