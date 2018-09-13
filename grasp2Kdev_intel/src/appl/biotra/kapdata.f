************************************************************************
*                                                                      *
      SUBROUTINE KAPDATA(NTESTG,NCORE1,NCORE2)
*                                                                      *
*   This subroutine determines the number of kappa quantum numbers     *
*   KAMAX together with the number of orbitals of each kappa.          *
*                                                                      *
*   This subroutine also checks if the orbital ordering is normal      *
*   or reversed. If reversed then JA and JB have to be permuted        *
*   In the subroutine ti1tv                                            *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      PARAMETER (NLMAX = 20)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

      LOGICAL KLAR
      DIMENSION ISORT(2*NNNW)
*
      COMMON/ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
*
      COMMON/ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
*
*   At most NLMAX orbital symmetries
*
*   At most NLMAX orbitals of each symmetry
*    
      COMMON/SBDAT/NAKINVII(NNNW),NSHLII(NLMAX),NSHLPII(NLMAX,NLMAX),
     :             NAKINVFF(NNNW),NSHLFF(NLMAX),NSHLPFF(NLMAX,NLMAX),
     :             NSHLPPII(NLMAX,NNNW),NSHLPPFF(NLMAX,NNNW),
     :             NINII(NLMAX),NINFF(NLMAX),IKAPPA(NLMAX),KAMAX
*
      COMMON/ORBORD/NORDII,NORDFF
*
      NTESTL = 00000
      NTEST = MAX0(NTESTL,NTESTG)

      DO I = 1,2*NNNW
        ISORT(I) = 100
      ENDDO

      DO I = 1,NLMAX
        IKAPPA(I)   = 0
        NSHLII(I)   = 0
        NSHLFF(I)   = 0
        NINII(I)    = 0
        NINFF(I)    = 0
      ENDDO

      DO I = 1,NNNW
        NAKINVII(I) = 0
        NAKINVFF(I) = 0
      ENDDO

      DO J = 1,NLMAX
        DO I = 1,NLMAX
          NSHLPII(I,J) = 0
          NSHLPFF(I,J) = 0
        ENDDO
      ENDDO

      DO J = 1,NNNW
        DO I = 1,NLMAX
          NSHLPPII(I,J) = 0
          NSHLPPFF(I,J) = 0
        ENDDO
      ENDDO
*
*   Sort the kappa quantum numbers
*
      DO K = 1,NWII+NWFF
        IF (K.LE.NWII) THEN
          ITAL = NAKII(K)
        ELSE
          ITAL = NAKFF(K-NWII)
        ENDIF
        I = K - 1
        KLAR = .FALSE.
   12   IF (I.GT.0 .AND. .NOT.KLAR) THEN
          IF (ITAL.LE.ISORT(I)) THEN
            ISORT(I+1) = ISORT(I)
            I = I - 1
          ELSE
            KLAR = .TRUE.
          ENDIF
        GOTO 12
        ENDIF
        ISORT(I+1) = ITAL
      ENDDO
*
*   Determine the unique set of kappa IKAPPA
*
      KAMAX = 1
      IKAPPA(1) = ISORT(1)
      DO K = 1,2*NNNW-1
        IF (ISORT(K).NE.ISORT(K+1).AND.ISORT(K+1).LT.100) THEN
          KAMAX = KAMAX + 1
          IKAPPA(KAMAX) = ISORT(K+1)
        ENDIF
      ENDDO 
*
*  Make a connection between each kappa and a number in the
*  range [1,KAMAX] as to know on which file to dump the data
*  Determine the number of shells NSHLII for each I in the
*  range [1,KAMAX] 
*      
      IF (NTEST.GE.10) THEN
        WRITE(*,*) '******************'
        WRITE(*,*) ' Entering kapdata'
        WRITE(*,*) '******************'
        WRITE(*,*)
        WRITE(*,*) 'There are',NWII,'orbitals in the initial state'
        WRITE(*,*) 'with the following n and kappa quantum numbers'
      ENDIF

      DO J = 1,NWII
        IF (NTEST.GE.10) THEN
          WRITE(*,*) 'orbital number',J,'n and kappa',NPII(J),NAKII(J)
        ENDIF
        DO I = 1,KAMAX
          IF (IKAPPA(I).EQ.NAKII(J)) THEN
            NAKINVII(J) = I
            IF (J.LE.NCORE1) NINII(I) = NINII(I) + 1
            NSHLII(I) = NSHLII(I) + 1
            NSHLPII(I,NSHLII(I)) = J
            NSHLPPII(I,J) = NSHLII(I)
          ENDIF
        ENDDO
      ENDDO         

      IF (NTEST.GE.10) THEN
        WRITE(*,*) 'There are',NWFF,'orbitals in the final state'
        WRITE(*,*) 'with the following n and kappa quantum numbers'
      ENDIF
      DO J = 1,NWFF
        IF (NTEST.GE.10) THEN
          WRITE(*,*) 'orbital number',J,'n and kappa=',NPFF(J),NAKFF(J)
        ENDIF
        DO I = 1,KAMAX
          IF (IKAPPA(I).EQ.NAKFF(J)) THEN
            NAKINVFF(J) = I
            IF (J.LE.NCORE2) NINFF(I) = NINFF(I) + 1
            NSHLFF(I) = NSHLFF(I) + 1
            NSHLPFF(I,NSHLFF(I)) = J
            NSHLPPFF(I,J) = NSHLFF(I)
          ENDIF
        ENDDO
      ENDDO         

      IF (NTEST.GE.10) THEN
        WRITE(*,*) 'Total number of different kappa',KAMAX
        DO I = 1,KAMAX
          WRITE(*,*) 'L=',I,'corresponds to kappa=',IKAPPA(I) 
          WRITE(*,*) 'nr of init.  orb. with this kappa=',NSHLII(I)
          WRITE(*,*) 'nr of final. orb. with this kappa=',NSHLFF(I)
        ENDDO
        DO I = 1,KAMAX
          WRITE(*,*) 'Position in initial state list'
          DO J = 1,NSHLII(I)
            WRITE(*,*) 'L=',I,'orb. nr',J,',position',NSHLPII(I,J)
          ENDDO
          WRITE(*,*) 'Position in final   state list'
          DO J = 1,NSHLFF(I)
            WRITE(*,*) 'L=',I,'orb. nr',J,',position',NSHLPFF(I,J)
          ENDDO
          WRITE(*,*) 'Relative positions for initial state orbitals'
          DO J = 1,NWII
            IF (NSHLPPII(I,J).NE.0) THEN
              WRITE(*,*) 'Orbital',J,'is nr',NSHLPPII(I,J),
     :        'with kappa',IKAPPA(I)
            ENDIF
          ENDDO
          WRITE(*,*) 'Relative positions for final state orbitals'
          DO J = 1,NWFF
            IF (NSHLPPFF(I,J).NE.0) THEN
              WRITE(*,*) 'Orbital',J,'is nr',NSHLPPFF(I,J),
     :        'with kappa',IKAPPA(I)
            ENDIF
          ENDDO
        ENDDO

      ENDIF
*
*  Check if the orbital ordering is normal or reversed.
*
*  For normal ordering
*
*    NPII(NSHLPII(I,1)) < NPII(NSHLPII(I,2)) < ... < NPII(NSHLPII(I,NSHLII(I))
*
*  For reversed ordering
*
*    NPII(NSHLPII(I,1)) > NPII(NSHLPII(I,2)) > ... > NPII(NSHLPII(I,NSHLII(I))
*
      NORDII = 0
      DO I = 1,KAMAX
        NREF = 0
        DO J = 1+NINII(I),NSHLII(I)
          IF (NPII(NSHLPII(I,J)).LT.NREF) NORDII = 1
          NREF = NPII(NSHLPII(I,J))
        ENDDO
      ENDDO

      NORDFF = 0
      DO I = 1,KAMAX
        NREF = 0
        DO J = 1+NINFF(I),NSHLFF(I)
          IF (NPFF(NSHLPFF(I,J)).LT.NREF) NORDFF = 1
          NREF = NPFF(NSHLPFF(I,J))
        ENDDO
      ENDDO
          
      IF (NORDII.NE.NORDFF) THEN
        WRITE(*,*) ' Orbital order of the initial and final states'
        WRITE(*,*) ' should be the same. STOP'
        STOP
      ENDIF
*
*   If not normal order check if reversed order
*
      IF (NORDII.EQ.1) THEN
        DO I = 1,KAMAX
          NREF = 0
          DO J = NSHLII(I),1+NINII(I),-1
            IF (NPII(NSHLPII(I,J)).LT.NREF) NORDII = 2
            NREF = NPII(NSHLPII(I,J))
          ENDDO
        ENDDO

        DO I = 1,KAMAX
          NREF = 0
          DO J = NSHLFF(I),1+NINFF(I),-1
            IF (NPFF(NSHLPFF(I,J)).LT.NREF) NORDFF = 2
            NREF = NPFF(NSHLPFF(I,J))
          ENDDO
        ENDDO
      ENDIF

      IF (NORDII.EQ.2.OR.NORDFF.EQ.2) THEN
        WRITE(*,*) ' The orbital order is neither normal or reversed'
        WRITE(*,*) ' STOP'
        STOP
      ENDIF

Cw      write(*,*) 'Give nordii'
Cw      read(*,*) nordii
Cw      nordff = nordii

      IF (NORDII.EQ.0) THEN
        WRITE(*,*) ' Normal orbital ordering'
      ELSE
        WRITE(*,*) ' Reverse orbital ordering'
      ENDIF

      WRITE(*,*) '*****************'
      WRITE(*,*) ' Leaving kapdata'
      WRITE(*,*) '*****************'
      WRITE(*,*)

      RETURN
      END
