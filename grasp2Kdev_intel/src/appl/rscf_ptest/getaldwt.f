************************************************************************
*
      SUBROUTINE getaldwt (ncf, wt)
      IMPLICIT REAL*8          (A-H,O-Z)
*
*   Interactively determines the weights.
*
*   Call(s) to: [LIB92]: ITJPO
*
*   Written by Xinghong He                Last revision: 19 Mar 1999
*
************************************************************************
      REAL*8           wt(ncf)
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------

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

      IF (i .GT. 10) STOP 'Must be running un-attended'

!------------------------------------------------------------------

      IF (n159 .EQ. 1) THEN      ! Equal weight
         DO i = 1, ncf
            wt(i) = 1.D0
         ENDDO
         sumwgt = DBLE (ncf)
      ELSEIF (n159 .EQ. 5) THEN  ! Standard weight
         sumwgt = 0.D0
         DO i = 1, ncf
            FTJPOI = DBLE (ITJPO (I))
            WT(I) = FTJPOI
            sumwgt = sumwgt + FTJPOI
         ENDDO
      ELSEIF (n159 .EQ. 9) THEN  ! User-input weight

  123    WRITE (istde,*) 'Enter the (relative) weights of the',
     &                ncf,' levels :'
         READ (istdi,*) (wt(i), i = 1, ncmin)

         sumwgt = 0.D0
         DO i = 1, ncf
            IF (wt(I) .LE. 0.D0) THEN
               WRITE (istde,*) 'Weights must exceed 0;'
               GOTO 123
            ELSE
               sumwgt = sumwgt + wt(i)
            ENDIF
         ENDDO

      ELSE
         WRITE (istde,*) 'Impossible ! Because it was guarded'
         STOP
      ENDIF

      sumwgt = 1.D0 / sumwgt
      DO i = 1, ncf
         wt(i) = sumwgt * wt(i)
      ENDDO

      RETURN 
      END
