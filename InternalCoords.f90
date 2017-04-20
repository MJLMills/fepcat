MODULE InternalCoords

  IMPLICIT NONE

  CONTAINS

    PURE REAL(8) FUNCTION distance(x,y)

      IMPLICIT NONE
      REAL(4), INTENT(IN) :: x(3), y(3)
      INTEGER :: i
      REAL(8) :: sum, diff

      sum = 0.0d0
      DO i = 1, 3
        diff = x(i) - y(i)
        sum = sum + (diff*diff)
      ENDDO

      distance = SQRT(sum)

    END FUNCTION distance

END MODULE InternalCoords
