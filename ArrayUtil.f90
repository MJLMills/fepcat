MODULE ArrayUtil

  ! #DES: Utility functions for generating and manipulating arrays.

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ArrayEXP, linspace

CONTAINS

!*

  FUNCTION ArrayEXP(input)

    ! #DES: Take elementwise EXP of real(8) values and return as array
    ! Domain of EXP is R^1, Range is (0,inf)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: input(:)
    REAL(8) :: ArrayEXP(SIZE(input))
    INTEGER :: i

    ArrayEXP(:) = 0.0d0
    DO i = 1, SIZE(input)
      ArrayEXP(i) = EXP(input(i))
    ENDDO

  END FUNCTION ArrayEXP

!*

  FUNCTION linspace(min, max, nPoints)

    ! #DES: Generate a linearly spaced set of nPoints points between min and max inclusive

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: min, max
    INTEGER, INTENT(IN) :: nPoints
    REAL(8) :: linspace(nPoints), stepSize
    INTEGER :: i

    stepSize = (max - min) / (nPoints - 1)

    DO i = 1, nPoints - 1
      linspace(i) = min + ((i-1) * stepSize)
    ENDDO
    ! The expression inside DO... results in error (INT*REAL8), so that the full input range is not covered.
    ! In order to ensure the correct range is returned, make the final value exact.
    linspace(nPoints) = max

  END FUNCTION linspace

!*

END MODULE ArrayUtil
