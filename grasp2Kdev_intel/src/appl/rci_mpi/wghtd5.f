************************************************************************
*                                                                      *
      SUBROUTINE WGHTD5(iatjpo, iaspar)
*                                                                      *
*   Print  the  weights of the largest five CSF contributors to each   *
*   ASF.                                                               *
*                                                                      *
*   Call(s) to: ALLOC, DALLOC, ISPAR, ITJPO.                           *
*                                                                      *
*                                          Last updated: 02 Nov 1992   *
*                                                   CFF: 01 Sep 2014   *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*4 JLBL,LABJ,LABP
*
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
*
      DIMENSION WT(5),ICONF(5)
*
      COMMON/EIGVEC/PNEVEC
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
*
*   Allocate storage for local arrays
*
*CFF  These arrays are not needed
*     CALL ALLOC (PNTRWT,NCF,8)
*     CALL ALLOC (PNNEXT,NCF,4)
*
      WRITE (24,300)
* 
*CFF  Initialized NELT for NCF >5
      NELT = 5
      IF (NCF .LT. 5) NELT = NCF
*
*CFF  New routine for finding the five largest 
*     components of an eigenvector by an algorithm
*     that contains no goto statements and is much
*     faster. Expansion coefficients with sign
*     are displayed.  
*     Begin modified code
      DO 7 IV = 1,NVEC
*
         ICF = IVEC(IV)
*        .., initialize for a new eigenvector
         wt = 0.d0;  iconf = 0
         ibegin = (IV-1)*NCF 
         iend   = IV*NCF
*
         DO 4 I = ibegin+1, iend
            awt = abs(EVEC(i))
            if (awt > wt(5)) then
*              ... the weight needs to be inserted
               do j = 1,5
                  if ( awt > wt(j) ) then
                     if (j < 5) then
*                       ... shift right
                        do jj = 5, j+1,-1
                           wt(jj) = wt(jj-1)
                           iconf(jj) = iconf(jj-1)
                        end do
                     end if
*                    ...  insert the element
                     wt(j) = awt
                     iconf(j) = i 
                     exit
                  end if
               end do
            end if
    4    CONTINUE
*
*   Print first five elements of list.
*
         IP = (iaspar + 3 ) / 2
*        NELT = MIN (I,5)
         WRITE (24,301) ICF,LABJ(IATJPO),LABP(IP),
     :          (EVEC(iconf(j)),j = 1,NELT)
         WRITE (24,302) (ICONF(j)-ibegin,j = 1,NELT)
    7 CONTINUE
*     End modified code
*
*   Deallocate storage for local arrays
*
*      CALL DALLOC (PNTRWT)
*      CALL DALLOC (PNNEXT)
*
      RETURN
*
  300 FORMAT (/'Weights of major contributors to ASF:'
     :       //'Level J Parity      CSF contributions'/)
  301 FORMAT (I3,2X,2A4,5(3X,F8.5))
  302 FORMAT (13X,5(I11))

*
      END
