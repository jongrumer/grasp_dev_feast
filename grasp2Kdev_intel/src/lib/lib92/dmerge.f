C
C fsplited from spicmv.f
C
      SUBROUTINE dmerge (n, db, dc, idy, da, dconst, dl)
C 
C  this merge version has the advantage of loading da(i)
C  and idy(i) only once.
C
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION da(n), db(*), dc(*), idy(n)

      dsum = 0.0
      DO i = 1, n
         dsum = dsum + da(i) * db(idy(i))
         dc(idy(i)) = dc(idy(i)) + dconst * da(i)
      ENDDO
      dl = dsum

      RETURN
      END
