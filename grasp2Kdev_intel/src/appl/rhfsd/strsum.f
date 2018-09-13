************************************************************************
*                                                                      *
      SUBROUTINE STRSUM
*                                                                      *
*   Generates the first part of  hfs92.sum  (on stream 24 and 29).     *
*                                                                      *
*   Call(s) to: ENGOUTH                                                *
*                                                                      *
*   Written by P. Jonsson                 Last revision: 20 Oct 1999   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      COMMON/EIGVAL/EAV,PNEVAL
     :      /NSMDAT/HFSI,HFSD,HFSQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
*
*    Write nuclear data
*
      IF(DABS(HFSI) .LT. 0.01D-11) THEN
         WRITE(*,'(A)') ' Nuclear spin should be greater than 0, stop'
         STOP
      ELSE IF(DABS(HFSD) .LT. 0.01D-11) THEN
         WRITE(*,'(A)')
     :          ' Nuclear dipole moment should be greater than 0, stop'
         STOP
      END IF
!     WRITE (24,302) HFSI
!     WRITE (24,303) HFSD
!     WRITE (24,304) HFSQ
      WRITE (29,302) HFSI
      WRITE (29,303) HFSD
      WRITE (29,304) HFSQ
*
*    Write the list of eigenpair indices
*
CPJ      CALL ENGOUTH (EAV,EVAL,IATJPO,IASPAR,IVEC,NVEC,3)
*
      RETURN
*
  302 FORMAT ('Nuclear spin                        ',1PD22.15,' au')
  303 FORMAT ('Nuclear magnetic dipole moment      ',1PD22.15,' n.m.')
  304 FORMAT ('Nuclear electric quadrupole moment  ',1PD22.15,' barns')
*
      END
