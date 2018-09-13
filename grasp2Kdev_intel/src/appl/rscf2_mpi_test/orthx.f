
      SUBROUTINE orthx(nw,jp,lsort)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! nw     - input, total number of orbitals
! jp     - input, orbital jp, usually updated from from solve. It must
!          be the one after IORDER
! lsort  - input, logical, .t. then sort according in decreasing
!          self-consistency order.
! kfixed - output, the number of fixed orbitals of the same kappa 
!          (specified by nak(jp)) as the orbital jp.
! ktotal - output, the number of all orbitals of the same kappa as the
!          orbital jp.
! kindx  - output, an index array containing positions of the orbitals
!          in the decreasing self-consistency order. It has the size of
!          ktotal.
!
! For i=1 to kfixed, kindx(i) has the same order as array iorder (though
! this can be modified by requiring the smallest principal quantum number
! come first). For i=kfixed+1 to ktotal, kindx(i) is ordered as the 
! decreasing self-consistency. For orbitals with the same value of the
! self-consistency, typically 1.D20, the same order as iorder is used
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! All common blocks are input 

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL lfix, lsort
      CHARACTER*2 NH

      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))

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

      EPS = 0.01D 00*ACCY
      nakj = nak(jp)
      kfixed = 0

! Find orbitals of the same kappa as orbital jp.

D        write(istde,*)
D        write(istde,*) '============== ',np(jp),nh(jp),' ============'
D        write(istde,*) 'i kindx(i), np, nh, scnsty/E**2,E, pos/all'
D        write(istde,*)
D        write(istde,*) ' Original order'
D        write(istde,*)

      ktotal = 0
      DO i = 1, nw
!         The orbital jp is supposed to be an ordered one, so is k
         ki = iorder(i)
         IF(nak(ki) .EQ. nakj) THEN
            ktotal = ktotal + 1
            kindx(ktotal) = ki
D           write(istde,*) ktotal, ki, np(ki), nh(ki),
D    &                     scnsty(ki)/E(ki)**2 ,E(ki), i
         ENDIF
      ENDDO

! Sort and find the fixed orbitals.

      IF (lsort) THEN

D        write(istde,*) ' New order'
D        write(istde,*)

         DO i = 1, ktotal
            ki = kindx(i)
            DO j = i+1, ktotal
               kj = kindx(j)
               IF(scnsty(kj)/E(kj)**2 .LT. scnsty(ki)/E(ki)**2) THEN

!                  To preserve the order of the "degenerate" case
                  DO ishift = j, i+1, -1
                     kindx(ishift) = kindx(ishift-1)
                  ENDDO
                  ki = kj
               ENDIF
            ENDDO
            kindx(i) = ki
            IF(lfix(ki)) kfixed = kfixed + 1
D           write(istde,*) i, ki, np(ki), nh(ki),scnsty(ki)/E(ki)**2
         ENDDO
         
      ENDIF

! Schmidt orthogonalize all orbitals of the same kappa
! The fixed orbitals are not changed

      DO Lraw = kfixed + 1, ktotal
         L = kindx(Lraw)

         NAKL = NAK(L)

         MTP0 = MF(L)

         DO  Kraw = 1, Lraw - 1
            K = kindx(Kraw)
            OVRLAP = RINT (L,K,0)
!            write(istde,*) 'state',np(k),nh(k),OVRLAP

*   Schmidt orthogonalise

            PZ(L) = PZ(L)-OVRLAP*PZ(K);pRINT*,'orthy1',pz(L);
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
         FACTOR = 1.0D 00/SQRT (DNORM)

         PZ(L) = FACTOR*PZ(L);pRINT*,'orthy2',pz(L);
         DO I = 2, MTP
            PF(I,L) = FACTOR*PF(I,L)
            QF(I,L) = FACTOR*QF(I,L)
         ENDDO

*   Find new MF(L)

         MTP = MTP+1
   20    MTP = MTP-1
         IF (ABS (PF(MTP,L)) .LT. EPS) THEN
            PF(MTP,L) = 0.0D 00
            QF(MTP,L) = 0.0D 00
            GOTO 20
         ELSE
            MF(L) = MTP
         ENDIF

!---------------------------------------------------------------
!         DO kraw = 1, Lraw
!            k= kindx(kraw)
!             write(istde,*) np(L),nh(L),np(k),nh(k),rint(L,k,0)
!         ENDDO
!---------------------------------------------------------------

      ENDDO

      RETURN
      END

