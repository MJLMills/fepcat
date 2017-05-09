MODULE Movies

  IMPLICIT NONE

  CHARACTER(1),  PARAMETER :: sep = "/" ! unix-based only
  CHARACTER(10), PARAMETER :: fepPrefix   = "FEP-frame_"
  CHARACTER(12), PARAMETER :: fepusPrefix = "FEPUS-frame_"
  INTEGER, PARAMETER :: shUnit = 87

  CONTAINS

!*

  SUBROUTINE SetupMovies()

    USE MovieOptions, ONLY : ProcessNameList, movieOutputDir, plotShellScript
    USE FileIO, ONLY : OpenFile
    IMPLICIT NONE

    CALL ProcessNameList()
    CALL OpenFile(shUnit,TRIM(ADJUSTL(movieOutputDir))//"/"//TRIM(ADJUSTL(plotShellScript)),"write")
    WRITE(shUnit,'(A)') "rm list.dat"    

  END SUBROUTINE SetupMovies

!*

  SUBROUTINE TeardownMovies()

    USE FileIO, ONLY : CloseFile
    IMPLICIT NONE

    CALL CloseFile(shUnit)

  END SUBROUTINE TeardownMovies

!*

  SUBROUTINE MakeFepMovie(mappingEnergies,mask,lambda)

    ! #DES: Run FEP procedure for multiple points in the total simulation and
    !       print data for movie frames

    USE FreeEnergy, ONLY   : ComputeFEPProfile
    USE MovieOptions, ONLY : movieOutputDir, plotShellScript, skip, fepScript

    IMPLICIT NONE

    REAL(8),       INTENT(IN) :: mappingEnergies(:,:,:), lambda(:)
    LOGICAL,       INTENT(IN) :: mask(:,:)

    CHARACTER(3),  PARAMETER  :: extension = "csv"

    REAL(8)        :: output(SIZE(mappingEnergies,2),2) !, dG(SIZE(mappingEnergies,2))
    LOGICAL        :: timeMask(SIZE(mask,1),SIZE(mask,2))
    CHARACTER(100) :: outFileName, imageFileName
    INTEGER        :: fepstep, timestep, nFepSteps

    output(:,1) = lambda(:); output(:,2) = 0.0d0
    timeMask(:,:) = .FALSE.
    nFepSteps = SIZE(mappingEnergies,2)

    DO fepstep = 1, nFepSteps

      DO timestep = 1, SIZE(mappingEnergies,1), skip

        WHERE(mask(1:fepstep,1:timestep) .EQV. .TRUE.) timeMask(1:fepstep,1:timestep) = .TRUE.
        CALL ComputeFEPProfile(1,fepstep,mappingEnergies(:,:,:),timeMask(:,:),output(1:fepstep,2))

        outFileName   = createFileName(fepstep,timestep,fepPrefix,extension)
        imageFileName = createFileName(fepstep,timestep,fepPrefix,"png")
        CALL WriteFepDataFrame(output,TRIM(ADJUSTL(movieOutputDir))//sep//TRIM(ADJUSTL(outFileName)))
        CALL WriteShellCommands(outFileName,imageFileName,fepScript,shUnit)

      ENDDO

    ENDDO

  END SUBROUTINE MakeFepMovie

!*

  SUBROUTINE MakeFepusMovie(energyGap,groundStateEnergy,mappingEnergies,mask,Nbins,minPop)

    USE MovieOptions, ONLY : movieOutputDir, plotShellScript, skip, fepusScript
    USE FreeEnergy, ONLY : Histogram, ComputeFepProfile, FepUS

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: energyGap(:,:), groundStateEnergy(:,:), mappingEnergies(:,:,:)
    LOGICAL, INTENT(IN) :: mask(:,:)
    INTEGER, INTENT(IN) :: Nbins, minPop

    CHARACTER(3),  PARAMETER  :: extension = "csv"

    CHARACTER(100) :: outFileName, imageFileName
    INTEGER      :: binPopulations(Nbins,SIZE(energyGap,1)), binIndices(SIZE(energyGap,1),SIZE(energyGap,2))
    REAL(8)      :: binMidpoints(Nbins)
    REAL(8)      :: binGg(Nbins)  
    REAL(8)      :: G_FEP(SIZE(energyGap,1))
    LOGICAL      :: printBin(Nbins)
    INTEGER      :: bin
    REAL(8)      :: output(Nbins,2)
    INTEGER      :: fepstep, timestep

    ! There will be Nbins bins in the histogram. For those bins that contain no points want to write zeros to the output file.
    ! Start by making the histogram using all of the data.
    ! energyGap(nFepSteps,nTimeSteps) is the reaction coordinate data for the simulation period of interest.
    ! mask(nFepSteps,nTimesteps) is the on/off switches for each piece of data.
    ! the histogram routine looks at the data you give it and returns the populations of each bin from each fepstep, 
    ! the index of which bin each point is in (fep,time -> bin) and the x-value associated with each bin (i -> RC_i)
    ! I want to force the x-axis to be the same for every run, so need the binMidpoints to reflect the full data set.

    ! These are set to the size of the complete data set and can be sliced appropriately.
    CALL Histogram(energyGap(:,:),mask(:,:),Nbins,binPopulations(:,:),binIndices(:,:),binMidpoints(:))
    output(:,1) = binMidpoints(:)

    ! Separately, can calculate the accurate dG values associated with each value of lambda via FEP
    ! This is a choice, these could also be updated with the PMF in the following loops.
    CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:),mask(:,:),profile=G_FEP(:))

    binPopulations(:,:) = 0;
    binIndices(:,:)     = 0;
    binMidpoints(:)     = 0.0d0;

    DO fepstep = 1, SIZE(mappingEnergies,2) ! only work on the current fepstep

      DO timestep = 1, SIZE(mappingEnergies,1), skip

        ! Add the contributions from this fepstep from 1:timestep to the binPopulations and binIndices arrays, forcing the range of the x-axis
        CALL Histogram(energyGap(fepstep:fepstep,1:timestep),mask(fepstep:fepstep,1:timestep),Nbins,binPopulations(:,fepstep:fepstep),binIndices(fepstep:fepstep,1:timestep),binMidpoints(:), &
&                      MINVAL(energyGap(:,:),MASK=mask .EQV. .TRUE.), &
&                      MAXVAL(energyGap(:,:),MASK=mask .EQV. .TRUE.))

        CALL FepUS(mappingEnergies(:,:,:),groundStateEnergy(:,:),G_FEP(:),binPopulations(:,:),binIndices(:,:),PMF1D=binGg(:),minPop=minPop,useBin=printBin(:))

        output(:,2) = 0.0d0
        DO bin = 1, nBins
          IF (printBin(bin) .EQV. .TRUE.) THEN
            output(bin,2) = binGg(bin)
          ENDIF
        ENDDO

        imageFileName = createFileName(fepstep,timestep,fepusPrefix,"png")
        outFileName = createFileName(fepstep,timestep,fepusPrefix,extension)
        CALL WriteFepusDataFrame(output,TRIM(ADJUSTL(movieOutputDir))//sep//outFileName)
        CALL WriteShellCommands(outFileName,imageFileName,fepusScript,shUnit)

      END DO
    END DO

  END SUBROUTINE MakeFepusMovie

!*

  SUBROUTINE MakeMDMovie

    IMPLICIT NONE

    ! for a specific timestep, write an update plot of the energy variation, E_RS, E_PS and mapping energy

  END SUBROUTINE MakeMDMovie

!*

  SUBROUTINE WriteFepDataFrame(data,outFilePath)

    USE FileIO, ONLY : OpenFile, CloseFile
    USE Output, ONLY : WriteCSV2D

    IMPLICIT NONE

    REAL(8), INTENT(IN)      :: data(:,:)
    CHARACTER(*), INTENT(IN) :: outFilePath

    INTEGER, PARAMETER       :: csvUnit = 88
    CHARACTER(6),  PARAMETER :: head(2) = (/"lambda","deltaG"/)

      CALL OpenFile(csvUnit,TRIM(ADJUSTL(outFilePath)),"write")
      CALL WriteCSV2D(head,data,csvUnit)
      CALL CloseFile(csvUnit)

  END SUBROUTINE WriteFepDataFrame

!*

  ! only differece from previous routine is contents of 'head' - make generic?
  SUBROUTINE WriteFepusDataFrame(data,outFilePath)

    USE FileIO, ONLY : OpenFile, CloseFile
    USE Output, ONLY : WriteCSV2D

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: data(:,:)
    CHARACTER(*), INTENT(IN) :: outFilePath

    INTEGER, PARAMETER :: csvUnit = 2
    CHARACTER(6), PARAMETER :: head(2) = (/"dE","dg"/)

    CALL OpenFile(csvUnit,TRIM(ADJUSTL(outFilePath)),"write")
    CALL WriteCSV2D(head,data,csvUnit)
    CALL CloseFile(csvUnit)

  END SUBROUTINE WriteFepusDataFrame

!*

  SUBROUTINE WriteShellCommands(fileName,imageName,plotScript,unit)

    USE MovieOptions, ONLY : genericDataFileName, plotCommand

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: fileName, imageName, plotScript
    INTEGER, INTENT(IN)      :: unit

    WRITE(unit,'(A)') "echo "//TRIM(ADJUSTL(fileName))
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(fileName))//" "//TRIM(ADJUSTL(genericDataFileName))
    WRITE(unit,'(A)') TRIM(ADJUSTL(plotCommand))//" "//TRIM(ADJUSTL(plotScript))
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(genericDataFileName))//" "//TRIM(ADJUSTL(fileName))
    WRITE(unit,'(A)') "mv frame.png "//TRIM(ADJUSTL(imageName))
    WRITE(unit,'(A)') "echo file "//TRIM(ADJUSTL(imageName))//" >> list.dat"
    WRITE(unit,*)

  END SUBROUTINE WriteShellCommands

!*

  CHARACTER(100) FUNCTION createFileName(fepstep,timestep,prefix,extension)

    IMPLICIT NONE

    INTEGER,      INTENT(IN) :: fepstep, timestep
    CHARACTER(*), INTENT(IN) :: prefix, extension
    CHARACTER(8) :: timeChar
    CHARACTER(4) :: fepChar

    WRITE(fepChar, '(I0.4)') fepstep
    WRITE(timeChar,'(I0.8)') timestep
    createFileName = TRIM(ADJUSTL(prefix))//fepChar//"_"//timeChar//"."//TRIM(ADJUSTL(extension))

  END FUNCTION createFileName

END MODULE Movies
