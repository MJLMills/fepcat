MODULE Movies

  IMPLICIT NONE

  CHARACTER(1),  PARAMETER :: sep = "/" ! unix-based only
  CHARACTER(10), PARAMETER :: fepPrefix   = "FEP-frame_"
  CHARACTER(12), PARAMETER :: fepusPrefix = "FEPUS-frame_"

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

    REAL(8)        :: output(SIZE(mappingEnergies,2),2), dG(SIZE(mappingEnergies,2))
    CHARACTER(100) :: outFileName, imageFileName
    INTEGER        :: fepstep, timestep

    output(:,1) = lambda(:)

    CALL OpenFile(shUnit,TRIM(ADJUSTL(movieOutputDir))//sep//TRIM(ADJUSTL(plotShellScript)),"write")
    WRITE(shUnit,*) "rm list.dat"

    DO fepstep = 1, SIZE(mappingEnergies,2) ! only work on the current fepstep

      DO timestep = 1, SIZE(mappingEnergies,1), skip

        CALL ComputeFEPProfile(1,fepstep,mappingEnergies(1:timestep,:,:),mask(:,1:timestep),dG(:))
        output(fepstep,2) = dG(fepstep)

        outFileName   = createFileName(fepstep,timestep,fepPrefix,extension)
        imageFileName = createFileName(fepstep,timestep,fepPrefix,"png")
        CALL WriteFepDataFrame(output,TRIM(ADJUSTL(movieOutputDir))//sep//TRIM(ADJUSTL(outFileName)))
        CALL WriteShellCommands(outFileName,imageFileName,shUnit)

      ENDDO

    ENDDO

    CALL CloseFile(shUnit)    

  END SUBROUTINE MakeFepMovie

!*

  SUBROUTINE MakeFepusMovie(energyGap,groundStateEnergy,mappingEnergies,mask,Nbins,minPop,skip)

    USE MovieOptions, ONLY : movieOutputDir, plotShellScript
    USE FreeEnergy, ONLY : Histogram, ComputeFepProfile, FepUS

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: energyGap(:,:), groundStateEnergy(:,:), mappingEnergies(:,:,:)
    LOGICAL, INTENT(IN) :: mask(:,:)
    INTEGER, INTENT(IN) :: Nbins, minPop, skip

    CHARACTER(3),  PARAMETER  :: extension = "csv"

    CHARACTER(100) :: outFileName
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

        ! Then call the FepUs routine on the partial data
        ! When the number of timesteps is not complete, some data comes back as infinite
        CALL FepUS(mappingEnergies(1:timestep,1:fepstep,1:fepstep),groundStateEnergy(1:timestep,1:fepstep),G_FEP(:),binPopulations(:,1:fepstep),binIndices(1:fepstep,1:timestep),PMF1D=binGg(:),minPop=minPop,useBin=printBin(:))

        output(:,2) = 0.0d0
        DO bin = 1, nBins
          IF (printBin(bin) .EQV. .TRUE.) THEN
            output(bin,2) = binGg(bin)
          ENDIF
        ENDDO

        outFileName = createFileName(fepstep,timestep,fepusPrefix,extension)
        CALL WriteFepusDataFrame(output,TRIM(ADJUSTL(movieOutputDir))//sep//outFileName)

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

  SUBROUTINE WriteShellCommands(fileName,imageName,unit)

    USE MovieOptions, ONLY : genericDataFileName, plotCommand, plotScript

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: fileName, imageName
    INTEGER, INTENT(IN)      :: unit

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
