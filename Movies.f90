MODULE Movies

  IMPLICIT NONE

  CONTAINS

  ! DOES NOT WORK YET
  SUBROUTINE MakeFepMovie(mappingEnergies,mask,lambda,skip)

    ! #DES: Run FEP procedure for multiple points in the total simulation and
    !       print data for movie frames

    USE FileIO, ONLY : OpenFile, CloseFile
    USE Output, ONLY : WriteCSV2D
    USE FreeEnergy, ONLY : ComputeFEPProfile

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: mappingEnergies(:,:,:), lambda(:)
    LOGICAL, INTENT(IN) :: mask(:,:)
    INTEGER, INTENT(IN) :: skip
    INTEGER, PARAMETER  :: csvUnit = 86, shUnit = 87
    LOGICAL             :: localMask(SIZE(mask,1),SIZE(mask,2))
    INTEGER             :: fepstep, timestep, step
    CHARACTER(100)      :: outFileName
    CHARACTER(8)        :: timeChar
    CHARACTER(4)        :: fepChar
    CHARACTER(6)        :: head(2)
    REAL(8)             :: G_FEP(SIZE(mappingEnergies,2))
    REAL(8)             :: output(SIZE(mappingEnergies,2),2)

    head(1) = "lambda"; head(2) = "deltaG"
    DO step = 1, SIZE(mappingEnergies,2)
      output(step,1) = lambda(step)
    ENDDO

    CALL OpenFile(shUnit,"makeMovie.sh","write")

    DO fepstep = 1, SIZE(mappingEnergies,2)

      localMask(:,:) = mask(:,:)
      IF (fepstep < SIZE(mappingEnergies,2)-1) localMask(fepstep+2:,:) = .FALSE.

      DO timestep = 1, SIZE(mappingEnergies,1), skip

        localMask(fepstep+1,timestep:) = .FALSE.
        output(:,2) = 0.0d0
        CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:),localMask(:,:),G_FEP)

        DO step = 1, SIZE(mappingEnergies,2)
          output(step,1) = lambda(step)
          output(step,2) = G_FEP(step)
        ENDDO

        WRITE(fepChar,'(I0.4)') fepstep
        WRITE(timeChar,'(I0.8)') timestep
        outFileName = "FEP-FRAME_"//fepChar//"_"//timechar//".csv"
        CALL OpenFile(csvUnit,outFileName,"write")
        CALL WriteCSV2D(head,output,csvUnit)
        CALL CloseFile(csvUnit)

        WRITE(shUnit,'(A,A)') "RScript plot.r ", outFileName

      ENDDO
    ENDDO

    CALL CloseFile(shUnit)    

  END SUBROUTINE MakeFepMovie

END MODULE Movies
