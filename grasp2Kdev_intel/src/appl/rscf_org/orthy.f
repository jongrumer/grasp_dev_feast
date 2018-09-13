
      SUBROUTINE orthy(nw,jp,lsort)
      IMPLICIT REAL*8           (A-H, O-Z)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! nw     - input, total number of orbitals
! jp     - input, orbital jp, usually updated from  solve. It is
!          the one after IORDER
! lsort  - input, logical, .t. then sort according in decreasing
!          self-consistency order.
!
! The following quantities are not used outside this routine
!
! kfixed - output, the number of fixed orbitals of the same kappa 
!          (specified by nak(jp)) as the orbital jp.
! ktotal - output, the number of all orbitals of the same kappa as the
!          orbital jp.
! kindx  - output, an index array containing positions of the orbitals
!          in the decreasing self-consistency order. It has the size of
!          ktotal.
!
! Orbitals are grouped as fixed, unfixed (non-correlation, correlation).
! .For i=1 to kfixed, kindx(i) has the same order as array iorder (though
!  this can be modified by requiring the smallest principal quantum number
!  come first). 
! .For i=kfixed+1 to knon, kindx(i) is ordered as SCN/E**2.
! .For i=knon+1 to ktotal, kindx(i) is ordered as SCN. 
! .When there is degeneracy (same SCN/E**2 or same SCN), iorder is used.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! All common blocks are input 

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL lfix, lsort
      CHARACTER*2 NH

      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))

      COMMON/FIXD/NFIX,LFIX(NNNW)
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORBA/IORDER(NNNW)
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV

      INTEGER kindx(NNNW)
      COMMON/iounit/istdi,istdo,istde
      LOGICAL lcorre
      COMMON/corre/lcorre(NNNW)
      SAVE /corre/

!-----------------------------------------------------------------------
! For inverse - taken from orthor.f which was not used anywhere
! Not used currently
!-----------------------------------------------------------------------
      LOGICAL CHECK,NOINVT
      COMMON/INVT/NOINVT(NNNW)
      CHECK = .NOT. NOINVT(jp)

      EPS = 0.01D0*ACCY
      kfixed = 0
      nakj = nak(jp)

!-----------------------------------------------------------------------
! Find orbitals of the same kappa as orbital jp.
!-----------------------------------------------------------------------
      ktotal = 0
      DO i = 1, nw
!         Orbital jp is supposed to be an ordered one, so should be k
         ki = iorder(i)
         IF( nak(ki) .EQ. nakj ) THEN
            ktotal = ktotal + 1
            kindx(ktotal) = ki
         ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Find fixed orbitals and separate them from the rest.
!-----------------------------------------------------------------------

!      Place the fixed orbitals at the begining of the array, no 
!      sorting is done here. And find the number of non-correlation
!      orbitals.

      j = 0
      knon = 0
      DO i = 1, ktotal
         ki = kindx(i)
         IF( lfix(ki) ) THEN

            j = j + 1

             ! Shift array elements rather than exchanging
            DO ishift = i, j+1, -1
               kindx(ishift) = kindx(ishift-1)
            ENDDO
!            kindx(i) = kindx(j)

            kindx(j) = ki

         ELSE IF( .NOT. lcorre(ki) ) THEN
            knon = knon + 1
         ENDIF
      ENDDO
      kfixed = j

      IF( lsort ) THEN

!-----------------------------------------------------------------------
! Separate correlation from non-correlation orbitals
!-----------------------------------------------------------------------

         j = kfixed
         DO i = kfixed + 1, ktotal
            ki = kindx(i)
            IF( .NOT. lcorre(ki) ) THEN

               j = j + 1

                ! Shift array elements rather than exchanging
               DO ishift = i, j+1, -1
                  kindx(ishift) = kindx(ishift-1)
               ENDDO
!               kindx(i) = kindx(j)

               kindx(j) = ki
            ENDIF
         ENDDO

!-----------------------------------------------------------------------
! Sort non-correlation and correlation orbitals separately
! This can be done indepently since knon is known
!-----------------------------------------------------------------------

!           non-correlation, using criteria scnsty/E^2 
              
         DO i = kfixed + 1, kfixed + knon
            ki = kindx(i)
            DO j = i + 1, kfixed + knon
               kj = kindx(j)
               IF( scnsty(kj)/E(kj)**2 .LT. scnsty(ki)/E(ki)**2 ) THEN

                   ! No need to do shifting, exchanging is fine
                  kindx(j) = ki
                  ki = kj
               ENDIF
            ENDDO
            kindx(i) = ki
         ENDDO

!           correlation orbitals, using criteria scnsty
              
         DO i = kfixed + knon + 1, ktotal
            ki = kindx(i)
            DO j = i + 1, ktotal
               kj = kindx(j)
               IF( scnsty(kj) .LT. scnsty(ki) ) THEN

                   ! No need to do shifting, exchanging is fine
                  kindx(j) = ki
                  ki = kj
               ENDIF
            ENDDO
            kindx(i) = ki
         ENDDO
         
      ENDIF

! Finished sorting.
! Schmidt orthogonalize all orbitals of the same kappa
! The fixed orbitals are not changed

      DO Lraw = kfixed + 1, ktotal
         L = kindx(Lraw)

         NAKL = NAK(L)

         MTP0 = MF(L)

         DO  Kraw = 1, Lraw - 1
            K = kindx(Kraw)
            OVRLAP = RINT (L,K,0)

*   Schmidt orthogonalise

            PZ(L) = PZ(L)-OVRLAP*PZ(K)
            MTP = MAX (MF(L),MF(K))
            MTP0 = MAX(MTP0, MF(K))

            DO  I = 1, MTP
               PF(I,L) = PF(I,L)-OVRLAP*PF(I,K)
               QF(I,L) = QF(I,L)-OVRLAP*QF(I,K)
            ENDDO
         ENDDO

*   Normalise

         MTP = MTP0

         MF(L) = MTP
         DNORM = RINT (L,L,0)
         FACTOR = 1.D0/SQRT (DNORM)

!   Determine if inversion is necessary
!    This job was once done in dampor.f. If the purpose here is to make
!    PZ(L) positive, then it can be replaced by a simple IF structure
!    Now do it this way
!
!            IF (CHECK) THEN
!               CALL COUNT (PF(1,L),MTP,NNCFF,SGN)
!               IF (SGN .LT. 0.0D 00) THEN
!!                  INV = INV+1
!                  FACTOR = -FACTOR
!               ENDIF
!            ENDIF
!!         print *, check, np(L), nh(L), factor

         IF( PZ(L) .LT. 0.d0 ) FACTOR = -FACTOR

         PZ(L) = FACTOR*PZ(L)
         DO I = 2, MTP
            PF(I,L) = FACTOR*PF(I,L)
            QF(I,L) = FACTOR*QF(I,L)
         ENDDO

*   Find new MF(L)

         MTP = MTP+1
   20    MTP = MTP-1
         IF (ABS (PF(MTP,L)) .LT. EPS) THEN
            PF(MTP,L) = 0.D0
            QF(MTP,L) = 0.D0
            GOTO 20
         ELSE
            MF(L) = MTP
         ENDIF

      ENDDO

      RETURN
      END
