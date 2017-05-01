MODULE MovieOptions

  IMPLICIT NONE

  NAMELIST /movie/  movieOutputDir, outputDataFilePrefix, plotShellScript, plotCommand, plotScript, genericDataFileName
  CHARACTER(500) :: movieOutputDir, outputDataFilePrefix, plotShellScript, plotCommand, plotScript, genericDataFileName

  CONTAINS

!*

  SUBROUTINE ProcessNameList

    IMPLICIT NONE

    CALL ReadNameList()
    CALL CheckNameList()
    CALL PrintNameList()

  END SUBROUTINE ProcessNameList

!*

  SUBROUTINE ReadNameList()

    USE FileIO, ONLY : OpenFile, CloseFile

    IMPLICIT NONE
    INTEGER, PARAMETER :: nmlUnit = 69

    ! usable defaults
    movieOutputDir       = "fep-movie-files"
    outputDataFilePrefix = "FEP-FRAME"
    plotShellScript      = "makeMovie.sh"
    plotCommand          = "Rscript"
    plotScript           = "plot.r"
    genericDataFileName  = "data.csv"

    CALL OpenFile(nmlUnit,"movie.nml","read")
    READ(nmlUnit,NML=movie)
    CALL CloseFile(nmlUnit)

  END SUBROUTINE ReadNameList

!*

  SUBROUTINE CheckNameList()

    IMPLICIT NONE

    IF (movieOutputDir  == "") STOP "Blank output directory for movie data not valid"
    IF (plotShellScript == "") STOP "Blank movie generation script not valid"

  END SUBROUTINE CheckNameList

!*

  SUBROUTINE PrintNameList()

    USE Log, ONLY : logUnit
    IMPLICIT NONE

    WRITE(logUnit,'(A)')   "* Movie Generation Settings *"; WRITE(logUnit,*)
    WRITE(logUnit,'(A,A)') "Output directory for movie data: ", TRIM(ADJUSTL(movieOutputDir))
    WRITE(logUnit,'(A,A)') "Prefix for output data files:    ", TRIM(ADJUSTL(outputDataFilePrefix))
    WRITE(logUnit,'(A,A)') "Script to plot all data files:   ", TRIM(ADJUSTL(plotShellScript))
    WRITE(logUnit,'(A,A)') "Command for plotting data:       ", TRIM(ADJUSTL(plotCommand))
    WRITE(logUnit,'(A,A)') "Script for plotting data:        ", TRIM(ADJUSTL(plotScript))
    WRITE(logUnit,'(A,A)') "Name for generic data file:      ", TRIM(ADJUSTL(genericDataFileName))
    WRITE(logUnit,*)

  END SUBROUTINE PrintNameList

END MODULE MovieOptions