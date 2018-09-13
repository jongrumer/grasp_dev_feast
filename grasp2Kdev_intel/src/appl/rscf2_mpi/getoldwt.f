************************************************************************
*
      SUBROUTINE getoldwt (ndef, ncmin, wt)
      IMPLICIT NONE
*
*   Interactively determines the weights for (E)OL calculation.
*   It's modified to always ask the question for the weight
*
*   Call(s) to: [LIB92]: GETYN.
*
*   Written by Xinghong He                Last revision: 19 Mar 1999
*
************************************************************************
      INTEGER          ndef, ncmin,  i, n159
      REAL*8           wt(ncmin),  sumwgt

      INTEGER istdi, istdo, istde
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------

! Standard weights: ncmin=1, OL calculation

      IF (ncmin .EQ. 1) THEN
         wt(1) = -1.D0
         RETURN
      ENDIF

! Select a method to assign level weights for ncmin > 1 case

      WRITE (istde,*) 'level weights (1 equal;  5 standard;  9 user)'

      ! Let user try 10 times to get the correct input. 10 is BIG
      ! enough since the idea here is to allow user mistakes and
      ! at the same time to avoid an infinity loop.

      DO i = 1, 10
         READ (istdi,*) n159
         IF (n159 .EQ. 1 .OR. n159 .EQ. 5 .OR. n159 .EQ. 9) EXIT
         WRITE (istde, *) 'Input not correct, do it again. tried=', i
      ENDDO

      IF (i .GT. 10) STOP

!------------------------------------------------------------------

      IF (n159 .EQ. 1) THEN      ! Equal weight
         DO i = 1, ncmin
            wt(i) = -2.D0
         ENDDO
      ELSEIF (n159 .EQ. 5) THEN  ! Standard weight
         DO i = 1, ncmin
            wt(i) = -1.D0
         ENDDO
      ELSEIF (n159 .EQ. 9) THEN  ! User-input weight

  123    WRITE (istde,*) 'Enter the (relative) weights of the',
     &                ncmin,' levels :'
         READ (istdi,*) (wt(i), i = 1, ncmin)

         sumwgt = 0.D0
         DO i = 1, ncmin
            IF (wt(I) .LE. 0.D0) THEN
               WRITE (istde,*) 'Weights must exceed 0;'
               GOTO 123
            ELSE
               sumwgt = sumwgt + wt(i)
            ENDIF
         ENDDO
         sumwgt = 1.D0 / sumwgt
         DO i = 1, ncmin
            wt(i) = sumwgt * wt(i)
         ENDDO

      ELSE
         WRITE (istde,*) 'Impossible ! Because it was guarded'
         STOP
      ENDIF

      RETURN 
      END
