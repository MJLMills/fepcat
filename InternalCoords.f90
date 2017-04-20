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

  SUBROUTINE ReadReactionCoordinates()

    USE FileIO, ONLY : OpenFile, CloseFile
    IMPLICIT NONE
    INTEGER, PARAMETER :: unit = 22
    LOGICAL :: success
    CHARACTER(8) :: type
    INTEGER :: a, b, ios

    CALL OpenFile(unit,"rc.dat","READ",success)
    IF (success .EQV. .TRUE.) THEN

      DO WHILE (ios == 0)
        READ(unit,IOSTAT=ios) type, a, b
      ENDDO

      CALL CloseFile(unit)
    ENDIF

  END SUBROUTINE ReadReactioncoordinates

END MODULE InternalCoords
