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
    REAL(8)        :: output(SIZE(mappingEnergies,2),2), dG(SIZE(mappingEnergies,2))
    CHARACTER(100) :: outFileName, imageFileName
    INTEGER        :: fepstep, timestep

    output(:,1) = lambda(:)

    CALL OpenFile(shUnit,TRIM(ADJUSTL(movieOutputDir))//sep//TRIM(ADJUSTL(plotShellScript)),"write")
    WRITE(shUnit,*) "rm list.dat"

    output(:,2)    = 0.0d0
    localMask(:,:) = mask(:,:)

    DO fepstep = 1, SIZE(mappingEnergies,2) ! only work on the current fepstep

      DO timestep = 1, SIZE(mappingEnergies,1), skip

        CALL ComputeFEPProfile(1,fepstep,mappingEnergies(1:timestep,:,:),mask(:,1:timestep),dG(:))
        output(fepstep,2) = dG(fepstep)

        outFileName   = createFileName(fepstep,timestep,extension)
        imageFileName = createFileName(fepstep,timestep,"png")
        CALL WriteDataFrame(output,TRIM(ADJUSTL(movieOutputDir))//sep//TRIM(ADJUSTL(outFileName)))
        CALL WriteShellCommands(outFileName,imageFileName,shUnit)

      ENDDO

    ENDDO

    CALL CloseFile(shUnit)    

  END SUBROUTINE MakeFepMovie

!*

  SUBROUTINE MakeFepusMovie()

    IMPLICIT NONE

    ! similar to above but instead of plotting FEP values plot FEP/US values over time, i.e. dE vs dg#
    ! copy initially from the analysis FEP/US routine

  END SUBROUTINE MakeFepusMovie

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

  SUBROUTINE WriteShellCommands(fileName,imageName,unit)

    USE MovieOptions, ONLY : genericDataFileName, plotCommand, plotScript

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: fileName, imageName
    INTEGER, INTENT(IN)      :: unit

     WRITE(*,*) fileNAme, imageName, unit
    
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(fileName))//" "//TRIM(ADJUSTL(genericDataFileName))
    WRITE(unit,'(A)') TRIM(ADJUSTL(plotCommand))//" "//TRIM(ADJUSTL(plotScript))
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(genericDataFileName))//" "//TRIM(ADJUSTL(fileName))
    WRITE(unit,'(A)') "mv frame.png "//TRIM(ADJUSTL(imageName))
    WRITE(unit,'(A)') "echo file "//TRIM(ADJUSTL(imageName))//" >> list.dat"
    WRITE(unit,*)

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
