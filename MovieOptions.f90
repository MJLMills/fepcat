MODULE MovieOptions

  IMPLICIT NONE

  NAMELIST /movie/  movieOutputDir, plotShellScript, plotCommand, fepScript, fepusScript, genericDataFileName, skip
  CHARACTER(500) :: movieOutputDir, plotShellScript, plotCommand, fepScript, fepusScript, genericDataFileName
  INTEGER        :: skip

  CONTAINS

!*

  SUBROUTINE ProcessNameList()

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
    plotShellScript      = "makeMovie.sh"
    plotCommand          = "Rscript"
    fepScript            = "plotFep.r"
    fepusScript          = "plotFepus.r"
    genericDataFileName  = "data.csv"
    skip = 1

    CALL OpenFile(nmlUnit,"movie.nml","read")
    READ(nmlUnit,NML=movie)
    CALL CloseFile(nmlUnit)

  END SUBROUTINE ReadNameList

!*

  SUBROUTINE CheckNameList()

    IMPLICIT NONE

    IF (movieOutputDir  == "") STOP "Blank output directory for movie data not valid"
    IF (plotShellScript == "") STOP "Blank movie generation script not valid"
    IF (skip <= 0)             STOP "Skip parameter must be >= 1"

  END SUBROUTINE CheckNameList

!*

  SUBROUTINE PrintNameList()

    USE Log, ONLY : logUnit
    IMPLICIT NONE

    WRITE(logUnit,'(A)')      "* Movie Generation Settings *"; WRITE(logUnit,*)
    WRITE(logUnit,'(A,A)')    "Output directory for movie data: ", TRIM(ADJUSTL(movieOutputDir))
    WRITE(logUnit,'(A,A)')    "Script to plot all data files:   ", TRIM(ADJUSTL(plotShellScript))
    WRITE(logUnit,'(A,A)')    "Command for plotting data:       ", TRIM(ADJUSTL(plotCommand))
    WRITE(logUnit,'(A,A)')    "Script for plotting FEP data:    ", TRIM(ADJUSTL(fepScript))
    WRITE(logUnit,'(A,A)')    "Script for plotting FEP/US data: ", TRIM(ADJUSTL(fepusScript))
    WRITE(logUnit,'(A,A)')    "Name for generic data file:      ", TRIM(ADJUSTL(genericDataFileName))
    WRITE(logUnit,'(A,I0.4)') "Step between data points:        ", skip
    WRITE(logUnit,*)

  END SUBROUTINE PrintNameList

END MODULE MovieOptions
