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

!*

!  PURE REAL(8) FUNCTION angle(x,y,z)

!    IMPLICIT NONE
!    REAL(4), INTENT(IN) :: x(3), y(3), z(3)
!    REAL(4) :: diff(3))

!    angle = 0.0d0
   
    

!  END FUNCTION angle

!*

!  PURE REAL(8) FUNCTION dot_product(a,b)

!    IMPLICIT NONE
!    REAL(8), INTENT(IN) :: a(:), b(:)

!    dot_product = SUM(a(:)*b(:))

!  END FUNCTION

END MODULE InternalCoords
