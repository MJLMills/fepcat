MODULE Input2D

  IMPLICIT NONE

  INTEGER, PARAMETER :: maxIndices = 4       ! torsion is largest coordinate supported
  INTEGER :: dRC                             ! dimensionality of the reaction coordinate
  CHARACTER(8), ALLOCATABLE :: coordTypes(:) ! what type is each coordinate
  INTEGER, ALLOCATABLE :: coordIndices(:,:)  ! and what atoms is it composed of

  CONTAINS

!*

  SUBROUTINE ReadReactionCoordinates()

    USE FileIO, ONLY : OpenFile, CloseFile, FileLength
    IMPLICIT NONE
    INTEGER, PARAMETER :: unit = 22
    LOGICAL :: success
    INTEGER :: ios, ios2
    INTEGER :: i, j

    dRC = FileLength("rc.dat")
    ALLOCATE(coordTypes(dRC));              coordTypes(:) = ""
    ALLOCATE(coordIndices(dRC,maxIndices)); coordIndices(:,:) = 0

    CALL OpenFile(unit,"rc.dat","read",success)
    IF (success .EQV. .TRUE.) THEN

      i = 0; ios = 0
      DO WHILE (ios == 0) ! as soon as an attempt is made to read a non-existent line, loop ends

        i = i + 1
        READ(unit,'(A)',IOSTAT=ios,ADVANCE='NO') coordTypes(i)

        IF (ios == 0) THEN ! there are indices to read
          j = 0; ios2 = 0
          DO WHILE (ios2 == 0) ! once an index read fails, go to next
            j = j + 1
            READ(unit,'(I4)',IOSTAT=ios2,ADVANCE='NO') coordIndices(i,j)
          ENDDO
        ENDIF

      ENDDO

      CALL CloseFile(unit)

    ENDIF

    CALL WriteReactionCoordinates(6)

  END SUBROUTINE ReadReactionCoordinates  

!*

  SUBROUTINE WriteReactionCoordinates(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER :: coord, index

    WRITE(unit,'(A)') "*   Reaction Coordinate   *"; WRITE(unit,*)
    WRITE(unit,'(A,I3)') "Dimensionality: ", dRC

    WRITE(unit,'(A)') "Geometric Variables:"; WRITE(unit,*)
    DO coord = 1, SIZE(coordTypes)
      WRITE(unit,'(I2,1X,A,1X)',ADVANCE='NO') coord, coordTypes(coord)
      DO index = 1, SIZE(coordIndices(coord,:))
        WRITE(unit,'(I0.2,1X)',ADVANCE='NO') coordIndices(coord,index)
      ENDDO
      WRITE(unit,*)
    ENDDO

  END SUBROUTINE WriteReactionCoordinates


END MODULE Input2D
