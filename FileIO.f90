MODULE FileIO

  ! #DES: Subprograms for assisting with generic I/O

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: OpenFile, CloseFile, fileLengths, fileLength

CONTAINS

  SUBROUTINE OpenFile(unit,fileName,action,success)

    ! #DES: Attempt to open the file 'fileName' on unit 'unit' with action 'action'
    ! #DES: Variable 'action' can be read, write or readwrite

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    CHARACTER(*), INTENT(IN) :: fileName, action
    LOGICAL, INTENT(OUT), OPTIONAL :: success
    INTEGER :: ios

    IF (trim(adjustl(action)) /= "read" .AND.  &
&       trim(adjustl(action)) /= "write" .AND. &
&       trim(adjustl(action)) /= "readwrite")  &
&       STOP "ILLEGAL ACTION ARGUMENT PROVIDED TO OpenFile IN FileIO"

    IF (LEN(fileName) == 0) STOP "ZERO-LENGTH FILENAME PASSED TO OpenFile IN FileIO"
    IF (unit <= 0)          STOP "-VE UNIT NUMBER PASSED TO OpenFile IN FileIO"

    IF (PRESENT(success)) THEN
      IF (action == "read") THEN
        INQUIRE(FILE=trim(adjustl(fileName)),EXIST=success)
        IF (success .EQV. .FALSE.) RETURN
      ENDIF
    ENDIF

    OPEN(UNIT=unit,FILE=trim(adjustl(fileName)),IOSTAT=ios,ACTION=trim(adjustl(action)))
    IF (ios /= 0) THEN
      WRITE(*,*) "ERROR: COULD NOT OPEN ", trim(adjustl(fileName)), " FOR ", trim(adjustl(action)), " ON UNIT ", unit
      IF (PRESENT(success)) success = .FALSE.
    ENDIF

  END SUBROUTINE OpenFile

!*

  SUBROUTINE CloseFile(unit)

    ! #DES: Attempt to close the file on unit 'unit'

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER ios

    IF (unit <= 0) STOP "-VE UNIT NUMBER PASSED TO CloseFile IN FileIO"

    CLOSE(UNIT=unit,IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) "ERROR: COULD NOT CLOSE FILE ON UNIT ", unit
    ENDIF

  ENDSUBROUTINE CloseFile

!*

  INTEGER FUNCTION fileLength(fileName)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: fileName
    LOGICAL :: openSuccess
    INTEGER, PARAMETER :: eneUnit = 100
    INTEGER :: ios

    fileLength = 0

    CALL OpenFile(eneUnit,fileName,"read",success=openSuccess)
    IF (openSuccess .EQV. .TRUE.) THEN
      
      ios = 0; fileLength = 0
      DO WHILE (ios == 0)
        READ(eneUnit,*,IOSTAT=ios)
        IF (ios == 0) fileLength = fileLength + 1
      ENDDO
      CALL CloseFile(eneUnit)

    ELSE

      STOP "FileIO:fileLength: File could not be opened"

    ENDIF

  END FUNCTION fileLength

!*

  FUNCTION fileLengths(fileNames)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: fileNames(:)
    INTEGER :: fileLengths(SIZE(fileNames))
    INTEGER :: step

    fileLengths(:) = 0
    DO step = 1, SIZE(fileNames)
      fileLengths(step) = fileLength(fileNames(step))
    ENDDO

  END FUNCTION fileLengths

END MODULE FileIO
