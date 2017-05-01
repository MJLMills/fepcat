MODULE Movies

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE MakeFepMovie(mappingEnergies,mask,lambda,skip)

    ! #DES: Run FEP procedure for multiple points in the total simulation and
    !       print data for movie frames

    USE FileIO, ONLY     : OpenFile, CloseFile
    USE Output, ONLY     : WriteCSV2D
    USE FreeEnergy, ONLY : ComputeFEPProfile

    IMPLICIT NONE

    REAL(8), INTENT(IN)    :: mappingEnergies(:,:,:), lambda(:)
    LOGICAL, INTENT(IN)    :: mask(:,:)
    INTEGER, INTENT(IN)    :: skip
    INTEGER, PARAMETER     :: csvUnit = 86, shUnit = 87
    CHARACTER(15), PARAMETER :: movieDir = "fep-movie-files"

    LOGICAL             :: localMask(SIZE(mask,1),SIZE(mask,2))
    INTEGER             :: fepstep, timestep
    CHARACTER(100)      :: outFileName
    CHARACTER(6)        :: head(2)
!    REAL(8)             :: G_FEP(SIZE(mappingEnergies,2))
    REAL(8)             :: output(SIZE(mappingEnergies,2),2)

    head(1) = "lambda"; head(2) = "deltaG"
    DO fepstep = 1, SIZE(mappingEnergies,2)
      output(fepstep,1) = lambda(fepstep)
    ENDDO

    CALL OpenFile(shUnit,TRIM(ADJUSTL(movieDir))//"/makeMovie.sh","write")

    DO fepstep = 1, SIZE(mappingEnergies,2)

      localMask(:,:) = mask(:,:)
      IF (fepstep < SIZE(mappingEnergies,2)-1) localMask(fepstep+2:,:) = .FALSE.

      DO timestep = 1, SIZE(mappingEnergies,1), skip

        localMask(fepstep+1,timestep:) = .FALSE.

        output(:,2) = 0.0d0
        CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:),localMask(:,:),output(:,2))

        outFileName = createFileName(fepstep,timestep,"csv")
        CALL OpenFile(csvUnit,TRIM(ADJUSTL(movieDir))//"/"//TRIM(ADJUSTL(outFileName)),"write")
        CALL WriteCSV2D(head,output,csvUnit)
        CALL CloseFile(csvUnit)

        WRITE(shUnit,'(A,A,A)') "mv ", TRIM(ADJUSTL(outFileName)), " data.csv"
        WRITE(shUnit,'(A)')     "RScript plot.r data.csv" 
        WRITE(shUnit,'(A,A,A)') "mv data.csv ", TRIM(ADJUSTL(outFileName))


      ENDDO
    ENDDO

    CALL CloseFile(shUnit)    

  END SUBROUTINE MakeFepMovie

!*

  SUBROUTINE WriteShellScript(fileName,unit)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: fileName
    INTEGER, INTENT(IN)      :: unit
    CHARACTER(8), PARAMETER  :: genericFileName = "data.csv"
    CHARACTER(6), PARAMETER  :: plotScript = "plot.r"
    CHARACTER(7), PARAMETER  :: command = "Rscript"
    
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(fileName))//" "//TRIM(ADJUSTL(genericFileName))
    WRITE(unit,'(A)') TRIM(ADJUSTL(command))//" "//TRIM(ADJUSTL(plotScript))//" "//TRIM(ADJUSTL(genericFileName))
    WRITE(unit,'(A)') "mv "//TRIM(ADJUSTL(genericFileName))//" "//TRIM(ADJUSTL(fileName))

  END SUBROUTINE WriteShellScript

  CHARACTER(100) FUNCTION createFileName(fepstep,timestep,extension)

    IMPLICIT NONE
    INTEGER,      INTENT(IN) :: fepstep, timestep
    CHARACTER(*), INTENT(IN) :: extension    
    CHARACTER(9), PARAMETER  :: prefix = "FEP-FRAME"
    CHARACTER(8) :: timeChar
    CHARACTER(4) :: fepChar

    WRITE(fepChar, '(I0.4)') fepstep
    WRITE(timeChar,'(I0.8)') timestep
    createFileName = TRIM(ADJUSTL(prefix))//fepChar//"_"//timeChar//"."//TRIM(ADJUSTL(extension))

  END FUNCTION createFileName

END MODULE Movies
