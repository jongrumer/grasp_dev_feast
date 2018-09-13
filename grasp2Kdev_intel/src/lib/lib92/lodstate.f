************************************************************************
*
      SUBROUTINE lodstate (nblock, ncfblk, idblk, nevblk, ncmaxblk)
*
*   Print block info and ask ASF serial numbers for each block
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
*
************************************************************************
      CHARACTER*8 idblk(*) ! idblk(1:nblock)

      INTEGER ncfblk(1:nblock), nevblk(1:nblock), ncmaxblk(1:nblock)

      POINTER (PCCMIN,ICCMIN(*))
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX ! NCMAX not used; --> ncmaxblk()

      COMMON/iounit/istdi,istdo,istde

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

         WRITE(734,'(a)') trim(str)    ! write to rscf.log file see, rscf

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
      ENDDO

      IF (ncmin .EQ. 0) THEN
         WRITE (istde,*) 'At least one state should be selected'
         GOTO 123
      ENDIF

      RETURN
      END
