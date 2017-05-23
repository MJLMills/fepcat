PROGRAM AveGeom

  USE Util, ONLY : Startup, Cleanup
  IMPLICIT NONE

  ! read all the structures and average each FEP step
  CALL Startup(readCoords=.TRUE.)
  CALL MeanFepStructures()
  CALL Cleanup()

  CONTAINS

!*

  SUBROUTINE MeanFepStructures()

    USE Input, ONLY : trajectory
    USE Output, ONLY : WriteXYZ
    USE FileIO, ONLY : OpenFile, CloseFile
    IMPLICIT NONE
    INTEGER :: fepstep, timestep
    CHARACTER(50) :: title, fileName
    CHARACTER(5) :: elements(SIZE(trajectory,3))
    REAL(8) :: meanStructure(SIZE(trajectory,1),3,SIZE(trajectory,3))
    INTEGER, PARAMETER :: xyzUnit = 42

    elements(:) = "H"

    DO fepstep = 1, SIZE(trajectory,1)
    
      ! sum the atomic positions over all timesteps of this simulation
      DO timestep = 1, SIZE(trajectory,4)
        meanStructure(fepstep,:,:) = meanStructure(fepstep,:,:) + trajectory(fepstep,:,:,timestep)
      ENDDO

      meanStructure(fepstep,:,:) = meanStructure(fepstep,:,:) / SIZE(trajectory,4)

      WRITE(title,'(A,I4)') "MEAN STRUCTURE AT FEP STEP ",fepstep
      WRITE(fileName,'(A,I0.4,A)') "MeanStructure_", fepstep, ".xyz"

      CALL OpenFile(xyzUnit,TRIM(ADJUSTL(fileName)),"write")
      CALL WriteXYZ(title,elements,meanStructure(fepstep,:,:),xyzUnit)
      CALL CloseFile(xyzUnit)

    ENDDO

  END SUBROUTINE MeanFepStructures

END PROGRAM AveGeom
