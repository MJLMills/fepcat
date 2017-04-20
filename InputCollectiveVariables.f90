MODULE InputCollectiveVariables

  IMPLICIT NONE

!  PRIVATE
!  PUBLIC :: DetermineCollectiveVariables

  INTEGER :: dRC                             ! dimensionality of the reaction coordinate
  CHARACTER(8), ALLOCATABLE :: coordTypes(:) ! what type is each coordinate
  INTEGER, ALLOCATABLE :: coordIndices(:,:)  ! and what atoms is it composed of

  CONTAINS

!*

  SUBROUTINE DetermineCollectiveVariables(logUnit)

    USE FileIO, ONLY : FileLength
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit

    dRC = FileLength("rc.dat")
    CALL AllocateArrays(dRC)
    CALL ReadReactionCoordinates()
    CALL WriteReactionCoordinates(logUnit)

  END SUBROUTINE DetermineCollectiveVariables

!*

  SUBROUTINE AllocateArrays(N)

    IMPLICIT NONE
    INTEGER, PARAMETER :: maxIndices = 4 ! torsion is largest coordinate supported
    INTEGER, INTENT(IN) :: N

    ALLOCATE(coordTypes(N));              coordTypes(:) = ""
    ALLOCATE(coordIndices(N,maxIndices)); coordIndices(:,:) = 0

  END SUBROUTINE AllocateArrays

!*

  SUBROUTINE DeallocateArrays()

    IMPLICIT NONE

    IF (ALLOCATED(coordTypes))   DEALLOCATE(coordTypes)
    IF (ALLOCATED(coordIndices)) DEALLOCATE(coordIndices)

  END SUBROUTINE DeallocateArrays

!*

  SUBROUTINE ReadReactionCoordinates()

    USE FileIO, ONLY : OpenFile, CloseFile
    IMPLICIT NONE
    INTEGER, PARAMETER :: unit = 22
    LOGICAL :: success
    INTEGER :: ios, ios2
    INTEGER :: i, j

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
            READ(unit,'(I6)',IOSTAT=ios2,ADVANCE='NO') coordIndices(i,j)
          ENDDO
        ENDIF

      ENDDO

      CALL CloseFile(unit)

    ENDIF

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

END MODULE InputCollectiveVariables
