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
    USE MovieOptions, ONLY : movieOutputDir, plotShellScript, moviestep, fepScript

    IMPLICIT NONE

    REAL(8),       INTENT(IN) :: mappingEnergies(:,:,:), lambda(:)
    LOGICAL,       INTENT(IN) :: mask(:,:)

    CHARACTER(3),  PARAMETER  :: extension = "csv"

    REAL(8)        :: output(SIZE(mappingEnergies,2),2) !, dG(SIZE(mappingEnergies,2))
    LOGICAL        :: timeMask(SIZE(mask,1),SIZE(mask,2))
    CHARACTER(100) :: outFileName, imageFileName
    INTEGER        :: fepstep, timestep, nFepSteps, i
    REAL(8)        :: minPrinted, maxPrinted

    minPrinted = HUGE(0.d0); maxPrinted = TINY(0.0d0)
    output(:,1) = lambda(:); output(:,2) = 0.0d0
    timeMask(:,:) = .FALSE.
    nFepSteps = SIZE(mappingEnergies,2)

    DO fepstep = 1, nFepSteps

      DO timestep = 1, SIZE(mappingEnergies,1), moviestep

        WHERE(mask(1:fepstep,1:timestep) .EQV. .TRUE.) timeMask(1:fepstep,1:timestep) = .TRUE.
        CALL ComputeFEPProfile(1,fepstep,mappingEnergies(:,:,:),timeMask(:,:),output(1:fepstep,2))
        outFileName   = createFileName(fepstep,timestep,fepPrefix,extension)
        DO i = 1, fepstep
          IF (output(fepstep,2) < minPrinted) minPrinted = output(fepstep,2)
          IF (output(fepstep,2) > maxPrinted) maxPrinted = output(fepstep,2)
        ENDDO
        imageFileName = createFileName(fepstep,timestep,fepPrefix,"png")
        CALL WriteFepDataFrame(output,TRIM(ADJUSTL(movieOutputDir))//sep//TRIM(ADJUSTL(outFileName)))
        CALL WriteShellCommands(outFileName,imageFileName,fepScript,shUnit)

      ENDDO

    ENDDO

    WRITE(*,*) "FEP Movie min and max", minPrinted, maxPrinted

  END SUBROUTINE MakeFepMovie

!*

  SUBROUTINE WriteHistogram(binPopulations,binMidpoints)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: binPopulations(:,:)
    REAL(8), INTENT(IN) :: binMidpoints(:)
    INTEGER :: bin, step

    WRITE(*,*) "HISTOGRAM"; WRITE(*,*)

    DO bin = 1, SIZE(binMidpoints)
      WRITE(*,'(I5,F8.2,I5)',ADVANCE='NO') bin, binMidpoints(bin), SUM(binPopulations(bin,:))
      DO step = 1, SIZE(binPopulations,2)
        WRITE(*,'(I4)',ADVANCE='NO') binPopulations(bin,step)
      ENDDO
      WRITE(*,*)
    ENDDO

    WRITE(*,*)

  END SUBROUTINE WriteHistogram

!*

  SUBROUTINE MakeFepusMovie(energyGap,groundStateEnergy,mappingEnergies,mask,Nbins,minPop)

    USE MovieOptions, ONLY : movieOutputDir, plotShellScript, moviestep, fepusScript
    USE FreeEnergy,   ONLY : Histogram, ComputeFepProfile, FepUS

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
    REAL(8)      :: output(Nbins,2), minRC, maxRC, printMin, printMax
    INTEGER      :: fepstep, timestep

    minRC = MINVAL(energyGap(:,:),MASK=mask .EQV. .TRUE.)
    maxRC = MAXVAL(energyGap(:,:),MASK=mask .EQV. .TRUE.)
    printMin = HUGE(0.0d0); printMax = TINY(0.0d0)

    CALL Histogram(energyGap(:,:),mask(:,:),Nbins,binPopulations(:,:),binIndices(:,:),output(:,1))
    CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:),mask(:,:),profile=G_FEP(:))

    binPopulations(:,:) = 0;
    binIndices(:,:)     = 0;
    binMidpoints(:)     = 0.0d0;

    DO fepstep = 1, SIZE(mappingEnergies,2) ! only work on the current fepstep

      DO timestep = 1, SIZE(mappingEnergies,1), moviestep

        CALL Histogram(energyGap(fepstep:fepstep,1:timestep),         &
        &              mask(fepstep:fepstep,1:timestep), Nbins,       &

        &              binPopulations(1:nBins,fepstep:fepstep),       & ! These are to be computed, determine which bins the fepstep data go into
        &              binIndices(fepstep:fepstep,1:timestep),        & ! and count them. binMidpoints is returned but not required - waste of time.
        &              binMidpoints(1:nBins),                         & ! Force the histogram to range over the complete RC
        &              minRC, maxRC)

        ! At this point the mask is implicit in the binIndices - no simulation data past the current fepstep,timestep should have non-zero indices
        CALL FepUS(mappingEnergies(:,:,:),groundStateEnergy(:,:),G_FEP,binPopulations(:,:),binIndices(:,:),PMF1D=binGg(:),minPop=minPop,useBin=printBin(:))

        ! Print zeros for any non-significant histogram bins
        output(:,2) = 0.0d0
        DO bin = 1, nBins
          IF (printBin(bin) .EQV. .TRUE.) THEN
            output(bin,2) = binGg(bin)
            IF (output(bin,2) > printMax) printMax = output(bin,2)
            IF (output(bin,2) < printMin) printMin = output(bin,2)
          ENDIF
        ENDDO

        imageFileName = createFileName(fepstep,timestep,fepusPrefix,"png")
        outFileName = createFileName(fepstep,timestep,fepusPrefix,extension)
        CALL WriteFepusDataFrame(output,TRIM(ADJUSTL(movieOutputDir))//sep//outFileName)
        CALL WriteShellCommands(outFileName,imageFileName,fepusScript,shUnit)

      END DO
    END DO

    WRITE(*,*) "FEP/US min and max", printMin, printMax

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
