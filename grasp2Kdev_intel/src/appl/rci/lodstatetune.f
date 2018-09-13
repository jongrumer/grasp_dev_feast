************************************************************************
*
      SUBROUTINE lodstate (nblock, ncfblk, idblk, nevblk, ncmaxblk)
*
*   Print block info and ask ASF serial numbers for each block
*   Give energy shifts for selected diagonal elements
*
*   Input:
*      nblock, idblk(), ncfblk()
*
*   Output:
*      nevblk(), ncmaxblk()
*      ncmin, iccmin(1:ncmin) -- via items, in common block
*
*    Memories allocated here but not de-allocated are
*        pccmin (via items, deallocated here only if ncmin=0)
*
*  Written by Xinghong He                                  Jul 17 1997
*  Updated by Xinghong He                                  Jun 10 1998
*  Shift energies by Per J                                 March  2007
*
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*8 idblk(*) ! idblk(1:nblock)

      INTEGER ncfblk(1:nblock), nevblk(1:nblock), ncmaxblk(1:nblock)

      LOGICAL GETYN, YES, lshift

      POINTER (PCCMIN,ICCMIN(*))
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX ! NCMAX not used; --> ncmaxblk()
     :      /DEFAULT/NDEF

      COMMON/iounit/istdi,istdo,istde
      COMMON/hamshiftj/nshiftj(100),nasfshift(100,100),
     :       asfenergy(100,100),lshift

      CHARACTER*256 str

!-----------------------------------------------------------------------

!      CALL ALLOC (pncmaxblk, nblock, 4)
!      CALL ALLOC (pnevblk, nblock, 4)

      WRITE (istde,*) 'There are ', nblock, ' blocks  '
     &               ,'(block   J/Parity   NCF):'
      WRITE (istde, '( 4(I3, 1X, A5, I6, 5X) )') (
     &       j, idblk(j)(1:5), ncfblk(j), j = 1, nblock )
      DO i = 1, nblock
         nevblk(i) = 0
         ncmaxblk(i) = 0
      ENDDO

      WRITE (istde,*)
      WRITE (istde,*) 'Enter ASF serial numbers for each block'

      ncmin = 0
  123 CONTINUE
      DO jblock = 1, nblock
         ncf = ncfblk(jblock)
  234    WRITE (istde,*) 'Block ', jblock, '   ncf = ', ncf
     &                , ' id = ', idblk(jblock)(1:5)

*        ...Read and parse the list of levels

         READ (*,'(A)') str
         IF (NDEF.EQ.0) WRITE(734,'(a)') trim(str)

*        ...ICCMIN is allocated and accumulated in items
*        ...ncmin is both input and output parameters to items
         ncminold = ncmin
         CALL items (ncmin, ncf, str, ierr)
         IF (ncmin .EQ. 0) THEN
            CALL dalloc (pccmin)
         ENDIF
                    IF (ierr .LT. 0) GOTO 234
         nevblk(jblock) = ncmin - ncminold

*        ...Determine ncmaxblk
         ntmp = 0
         DO I = ncminold + 1, NCMIN
            ntmp = MAX (ntmp,ICCMIN(I))
         END DO
         ncmaxblk( jblock) = ntmp

*        ...Manually shift diagonal matrix elements

         IF (lshift) THEN
           WRITE (istde,*) 'Shift diagonal matrix elements (y/n)'
           YES = GETYN ()
           IF (YES) THEN
              WRITE (istde,*) 'How many elements should be shifted'
              read (*,*) nshifted
              nshiftj(jblock) = nshifted
              DO I = 1,nshifted
                 WRITE (istde,*) 'Give ASF serial number for element',I
                 read (*,*) nasfshift(jblock,I)
                 WRITE (istde,*) 'Give shift in cm-1 for element    ',I
                 read (*,*) asfenergy(jblock,I)
                 asfenergy(jblock,I) = asfenergy(jblock,I)/(2*109737.7)
               END DO
           ELSE
              nshiftj(jblock) = 0
           ENDIF
         ENDIF   


      ENDDO

      IF (ncmin .EQ. 0) THEN
         WRITE (istde,*) 'At least one state should be selected'
         GOTO 123
      ENDIF

      RETURN
      END
