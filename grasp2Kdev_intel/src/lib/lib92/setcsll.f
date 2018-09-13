************************************************************************
      SUBROUTINE SETCSLL (nunit, name, nblkin, nblock, ncfblk, ncftot,
     &                   idblk)
      IMPLICIT NONE
*
*  Open, read name file to get nblock, ncfblk(), idblk(), ncftot 
*
* Xinghong He 98-06-29
* 
************************************************************************
      INTEGER   nunit, nblkin, nblock, ncftot, ncfblk(*)
      CHARACTER name*(*), idblk(*)*8
! Locals
      INTEGER   i, ncsf, ios, ierr
      LOGICAL   FOUND
      CHARACTER str*15, CH*2, line3*100
!-----------------------------------------------------------------------

! Look for  <name>

      INQUIRE (FILE = name, EXIST = FOUND)
      IF (.NOT. FOUND) THEN
         PRINT *, name(1:LEN_TRIM(name)),' does not exist'
         STOP
      ENDIF

! Open it

      CALL OPENFL (nunit, name, 'FORMATTED', 'OLD', ierr)
      IF (ierr .EQ. 1) THEN
         PRINT *, 'Error when opening ',name(1:LEN_TRIM (name))
         STOP
      ENDIF

! Check the first record of the file; if not as expected, stop

      READ (nunit,'(1A15)',IOSTAT = ios) str
      IF ((ios .NE. 0) .OR.
     :    (str .NE. 'Core subshells:')) THEN
         PRINT *, 'Not a Configuration Symmetry List File;'
         CLOSE (nunit)
         STOP
      ENDIF

! Skip next 4 records

      DO i = 1, 4
         READ (nunit,*)
      ENDDO

! Determine the number of blocks in this file

      nblock    = 0
      ncsf      = 0

      ios = 0
      DO WHILE (ios .EQ. 0)
         READ (nunit,'(1A2)', IOSTAT = ios) CH
         IF (CH .EQ. ' *' .OR. IOS .NE. 0) THEN
            !.. a new block has been found
            nblock = nblock + 1
            PRINT *, 'Block ', nblock, ',  ncf = ', ncsf
            IF (nblock .GT. nblkin) THEN
               PRINT *, 'setcsll: Too many blocks(',nblock,')'
               PRINT *, 'Maximum allowed is ', nblkin
               STOP
            ENDIF
            i = LEN_TRIM (line3)
            idblk(nblock)  = line3(i-4:i)
            ncfblk(nblock) = ncsf
            ncsf = 0
            IF (IOS .EQ. 0) CYCLE
         ELSE
            READ (nunit,*)
            READ (nunit,'(A)') line3
            ncsf = ncsf + 1
         END IF
      ENDDO

! Obtain ncftot

      ncftot = 0
      DO i = 1, nblock
         ncftot = ncftot + ncfblk(i)
      ENDDO

      RETURN
      END
