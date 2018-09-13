************************************************************************
      SUBROUTINE CSFWGT (LSTDIO)
*                                                                      *
*   Print  the  weights of the largest five CSF contributors to each   *
*   ASF.                                                               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, ISPAR, ITJPO.          *
*                                                                      *
*                                          Last updated: 21 Dec 1992   *
*                                          Last updated: 24 Feb 1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL LSTDIO
      CHARACTER*256 RECORD
      CHARACTER*8 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
*
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))

      COMMON/DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVEC/PNEVEC
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /SYMA/PIATJP,PIASPA

      COMMON/iounit/istdi,istdo,istde

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pncfpast, ncfpast(1))
      POINTER (pncminpast, ncminpast(1))
      POINTER (pnevecpast, nevecpast(1))
      COMMON/pos/pncfpast,pncminpast,pnevecpast,ncftot,nvecsiz

      POINTER (pidxblk, idxblk(1))   ! idx(i= 1,ncmin) is the block where 
      COMMON/blkidx/pidxblk          ! the i_th eigenvalue comes from

      DIMENSION coeff(5),ICONF(5)

      IF(lstdio) THEN
         WRITE(istdo,300)
      ELSE
         WRITE(24,300)
      ENDIF

      DO 40 IV = 1, NCMIN
!         loop over eigenvectors

         jblock = idxblk(IV)
         ncf = ncfblk(jblock)
         ncfpat = ncfpast(jblock)
         ncminpat = ncminpast(jblock)
         nevecpat = nevecpast(jblock)
         nevecoff = nevecpat + (iv - ncminpat - 1) * ncf
         NELT = MIN(5, ncf)		!Find maximum 5 ... within block
         ICF = ICCMIN(IV)
         ivtjpo = iatjpo(iv)    ! j-value related
         ivspar = iaspar(iv)    ! parity related

!         Select 5 (at most) pairs and put into arrays coeff, iconf

!         icount = 0
!         DO I = 1, ncf
!            IF(itjpo(i+ncfpat) .NE. ivtjpo .OR. 
!     &         ispar(i+ncfpat) .NE. ivspar) CYCLE
!            icount = icount + 1
!            coeff(icount) = evec(nevecoff + i)
!            iconf(icount) = i
!            IF (icount .EQ. nelt) GOTO 10
!!               normal exit
!         ENDDO
!!               abnormal exit (should not happen in block version)
!          nelt = icount
!          STOP 'csfwgt: Something wrong here'
!
!   10    CONTINUE
!         igroup = i

         DO i = 1, nelt
            coeff(i) = evec(nevecoff + i)
            iconf(i) = i
         ENDDO

!         sort the first nelt in decreasing order
         do i = 1, nelt
         do j = i+1, nelt
            if (abs(coeff(j)) .gt. abs(coeff(i))) then
               temp = coeff(i)
               coeff(i) = coeff(j)
               coeff(j) = temp
               itemp = iconf(i)
               iconf(i) = iconf(j)
               iconf(j) = itemp
            endif
         enddo
         enddo

         do i = nelt+1 , ncf
!            IF(itjpo(i+ncfpat) .NE. ivtjpo .OR. 
!     &         ispar(i+ncfpat) .NE. ivspar) GOTO 20
            w = evec(nevecoff + i)
            if (w .ne. 0.d0 .and. abs(w) .gt. abs(coeff(nelt)) ) then
!               we have a non-zero value larger than the largest so far
               do j = 1, nelt
                  if (abs(w) .gt. abs(coeff(j))) then
                     do irem = nelt, j+1, -1
                        coeff(irem) = coeff(irem-1)
                        iconf(irem) = iconf(irem-1)
                     enddo
                     coeff(j) = w
                     iconf(j) = i
                     goto 20
                  endif
               enddo
            endif
   20       CONTINUE
         enddo

        ip = (iaspar(iv) + 3) / 2

         IF (LSTDIO) THEN
            WRITE (istdo,320) jblock,ICF,LABJ(IVTJPO),LABP(IP),
     :                     (coeff(I),i=1,nelt)
            WRITE (istdo,330) (ICONF(I),I=1,NELT)
         ELSE

            WRITE (24,320) jblock,ICF,LABJ(IVTJPO),LABP(IP),
     :                     (coeff(I),i=1,nelt)
            WRITE (24,330) (ICONF(I),I=1,NELT)
         ENDIF
   40 CONTINUE

  300 FORMAT (/'Weights of major contributors to ASF:'
     :       //'Block Level J Parity      CSF contributions'/)
  310 FORMAT (1X, A14,80A)
  320 FORMAT (I3,1X,I5,2X,2A4,5(3X,F8.4))
  330 FORMAT (19X,5(3X,I8))

      RETURN
      END

