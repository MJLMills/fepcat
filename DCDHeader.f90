MODULE DcdHeader

   ! #DES: This module deals with the dcd file format's header, which contains a lt of information

  IMPLICIT NONE

  INTEGER :: unit

  PRIVATE
  PUBLIC :: parseHeader  

CONTAINS

  SUBROUTINE parseHeader(dcdUnit,debug,nAtoms,nFrames)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dcdUnit
    LOGICAL, INTENT(IN) :: debug
    INTEGER, INTENT(OUT) :: nAtoms, nFrames

    unit = dcdUnit 
    CALL parseHeaderLineOne(debug,nFrames)
    CALL parseHeaderLineTwo(debug)
    CALL parseHeaderLineThree(debug,nAtoms)

  END SUBROUTINE parseHeader

!*

  SUBROUTINE parseHeaderLineOne(debug,nFrames)

    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: debug
    INTEGER, INTENT(OUT) :: nFrames
    INTEGER :: ios
    CHARACTER(4) :: trajectoryType
    INTEGER(4) :: nStepsBefore, interval, nSteps, intervalVelChek, na6, na7, nDegF
    INTEGER(4) :: nFixed, na10, constP, na12, na13, na14, na15, na16, na17, na18, na19, charmmVersion

    READ(unit,IOSTAT=ios) trajectoryType, nFrames, nStepsBefore, interval, nSteps,    &
                          intervalVelChek, na6, na7, nDegF, nFixed, na10, constP,     &
                          na12, na13, na14, na15, na16, na17, na18, na19, charmmVersion

    IF (ios /= 0) STOP "Error reading first line of header"
    IF (TRIM(ADJUSTL(trajectoryType)) /= 'CORD') STOP "Not a coordinate file"

    IF (debug .EQV. .TRUE.) THEN
      WRITE(*,*) "Trajectory Type:         ", trajectoryType
      WRITE(*,*) "Num. of Frames:          ", nFrames
      WRITE(*,*) "Num. Steps Before:       ", nStepsBefore
      WRITE(*,*) "Interval                 ", interval
      WRITE(*,*) "Num. of Steps:           ", nSteps
      WRITE(*,*) "Interval Velocity Check: ", intervalVelChek
      WRITE(*,*) "Degrees of Freedom:      ", nDegF
      WRITE(*,*) "Num. Fixed:              ", nFixed
      WRITE(*,*) "Constant Pressure:       ", constP
      WRITE(*,*) "CHARMM Version:          ", charmmVersion
      WRITE(*,*)
    ENDIF

  END SUBROUTINE parseHeaderLineOne

!*

  SUBROUTINE parseHeaderLineTwo(debug)

    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: debug
    INTEGER :: ios
    INTEGER, PARAMETER :: maxMaskRows = 20
    INTEGER :: maskRows
    CHARACTER(80) :: titlerow, topology, maskRow(maxMaskRows)

    READ(unit,IOSTAT=ios) maskRows, titleRow, topology, maskRow(1:maskRows-2)
    IF (ios /= 0) STOP "Error reading second line"
    IF (titleRow(1:16) /= 'Q DCD trajectory') STOP "Not a Q trajectory file"

    IF (debug .EQV. .TRUE.) THEN
      WRITE(*,*)
      WRITE(*,*) "Mask Rows: ", maskRows
      WRITE(*,*) "Title Row: ", TRIM(ADJUSTL(titleRow))
      WRITE(*,*) "Topology:  ", TRIM(ADJUSTL(topology))
      WRITE(*,*)
    ENDIF

  END SUBROUTINE parseHeaderLineTwo

!*

  SUBROUTINE parseHeaderLineThree(debug,nAtoms)

    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: debug
    INTEGER, INTENT(OUT) :: nAtoms
    INTEGER :: ios

    READ(unit,IOSTAT=ios) nAtoms
    IF (ios /= 0) STOP "Error reading third line"

    IF (debug .EQV. .TRUE.) THEN
      WRITE(*,*) "Atoms: ", nAtoms
      WRITE(*,*)
    ENDIF

  END SUBROUTINE parseHeaderLineThree

END MODULE DcdHeader

