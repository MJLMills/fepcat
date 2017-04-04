MODULE DCDFiles

  IMPLICIT NONE
  INTEGER, PARAMETER :: dcdUnit = 14

  PRIVATE
  PUBLIC :: ReadDCD, CountAtoms

CONTAINS

  SUBROUTINE CountAtoms(fileName,nAtoms)

    USE DCDHeader, ONLY : parseHeader

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: fileName
    INTEGER, INTENT(OUT) :: nAtoms
    LOGICAL, PARAMETER :: debug = .FALSE.
    INTEGER :: nTimesteps
   
    OPEN(UNIT=dcdUnit,FILE=fileName,STATUS='old',FORM='unformatted',ACTION='read')
    CALL parseHeader(dcdUnit,debug,nAtoms,nTimesteps)
    CLOSE(dcdUnit)

  END SUBROUTINE CountAtoms

!*

  SUBROUTINE ReadDCD(fileName,coordinates)

    USE DCDHeader, ONLY : parseHeader

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: fileName
    REAL(4), INTENT(OUT) :: coordinates(:,:,:)
    LOGICAL, PARAMETER :: debug = .FALSE.
    INTEGER :: nAtoms, nTimesteps, ios
    INTEGER :: timestep, axis
    REAL(4), ALLOCATABLE :: data(:,:,:) ! 3 x nAtoms x nTimesteps

    OPEN(UNIT=dcdUnit,FILE=fileName,STATUS='old',FORM='unformatted',ACTION='read',IOSTAT=ios)
    IF (ios /= 0) STOP "Failed to open unformatted trajectory file"

    CALL parseHeader(dcdUnit,debug,nAtoms,nTimesteps)

    ALLOCATE(data(3,nAtoms,nTimesteps))
    DO timestep = 1, nTimesteps
      DO axis = 1, 3
        READ(dcdUnit,IOSTAT=ios) data(axis,1:nAtoms,timestep)
        IF (ios /= 0) STOP "Error reading coordinate line"
      ENDDO
    ENDDO

    Coordinates = data

    DEALLOCATE(data)
    CLOSE(dcdUnit)

  END SUBROUTINE ReadDCD

END MODULE DCDFiles
