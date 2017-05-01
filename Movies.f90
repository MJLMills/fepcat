MODULE Movies

  IMPLICIT NONE

  CHARACTER(1), PARAMETER :: sep = "/" ! unix-based only

  CONTAINS

!*

  SUBROUTINE MakeFepMovie(mappingEnergies,mask,lambda,skip)

    ! #DES: Run FEP procedure for multiple points in the total simulation and
    !       print data for movie frames

    USE FileIO, ONLY       : OpenFile, CloseFile
    USE FreeEnergy, ONLY   : ComputeFEPProfile
    USE MovieOptions, ONLY : movieOutputDir, plotShellScript

    IMPLICIT NONE

    REAL(8),       INTENT(IN) :: mappingEnergies(:,:,:), lambda(:)
    LOGICAL,       INTENT(IN) :: mask(:,:)
    INTEGER,       INTENT(IN) :: skip

    INTEGER,       PARAMETER  :: shUnit = 87
    CHARACTER(3),  PARAMETER  :: extension = "csv"

    LOGICAL        :: localMask(SIZE(mask,1),SIZE(mask,2))
    REAL(8)        :: output(SIZE(mappingEnergies,2),2)
    CHARACTER(100) :: outFileName
    INTEGER        :: fepstep, timestep

    output(:,1) = lambda(:)

    CALL OpenFile(shUnit,TRIM(ADJUSTL(movieOutputDir))//sep//TRIM(ADJUSTL(plotShellScript)),"write")

    DO fepstep = 1, SIZE(mappingEnergies,2)

      localMask(:,:) = mask(:,:)
      IF (fepstep < SIZE(mappingEnergies,2)-1) localMask(fepstep+2:,:) = .FALSE.

      DO timestep = 1, SIZE(mappingEnergies,1), skip

        localMask(fepstep+1,timestep:) = .FALSE.

        CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:),localMask(:,:),output(:,2))

        outFileName = createFileName(fepstep,timestep,extension)
        CALL WriteDataFrame(output,movieOutputDir//sep//TRIM(ADJUSTL(outFileName)))
        CALL WriteShellCommands(outFileName,shUnit)

      ENDDO

    ENDDO

    CALL CloseFile(shUnit)    

  END SUBROUTINE MakeFepMovie

!*

  SUBROUTINE WriteDataFrame(data,outFilePath)

    USE FileIO, ONLY : OpenFile, CloseFile
    USE Output, ONLY : WriteCSV2D

    IMPLICIT NONE

    REAL(8), INTENT(IN)      :: data(:,:)
    CHARACTER(*), INTENT(IN) :: outFilePath

    INTEGER, PARAMETER       :: csvUnit = 86
    CHARACTER(6),  PARAMETER :: head(2) = (/"lambda","deltaG"/)

      CALL OpenFile(csvUnit,TRIM(ADJUSTL(outFilePath)),"write")
      CALL WriteCSV2D(head,data,csvUnit)
      CALL CloseFile(csvUnit)

  END SUBROUTINE WriteDataFrame

!*

  SUBROUTINE WriteShellCommands(fileName,unit)

    USE MovieOptions, ONLY : genericDataFileName, plotCommand, plotScript

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: fileName
    INTEGER, INTENT(IN)      :: unit
    
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(fileName))//" "//TRIM(ADJUSTL(genericDataFileName))
    WRITE(unit,'(A)') TRIM(ADJUSTL(plotCommand))//" "//TRIM(ADJUSTL(plotScript))//" "//TRIM(ADJUSTL(genericDataFileName))
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(genericDataFileName))//" "//TRIM(ADJUSTL(fileName))

  END SUBROUTINE WriteShellCommands

!*

  CHARACTER(100) FUNCTION createFileName(fepstep,timestep,extension)

    USE MovieOptions, ONLY : outputDataFilePrefix

    IMPLICIT NONE

    INTEGER,      INTENT(IN) :: fepstep, timestep
    CHARACTER(*), INTENT(IN) :: extension
    CHARACTER(8) :: timeChar
    CHARACTER(4) :: fepChar

    WRITE(fepChar, '(I0.4)') fepstep
    WRITE(timeChar,'(I0.8)') timestep
    createFileName = TRIM(ADJUSTL(outputDataFilePrefix))//fepChar//"_"//timeChar//"."//TRIM(ADJUSTL(extension))

  END FUNCTION createFileName

END MODULE Movies
