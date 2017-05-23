! These are routines common to most potential main programs that use the modules within.
! Reads input, creates log file, then post-calculation deallocates all data and closes the log.

MODULE Util

  CONTAINS

    SUBROUTINE Startup(readCoords)

      ! #DES: Setup the calculation by creating a log file and reading the necessary input.

      USE Log, ONLY   : CreateLogFile
      USE Input, ONLY : ReadInput

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: readCoords

        CALL CreateLogFile()
        CALL ReadInput(readCoords)

    END SUBROUTINE StartUp

!*

    SUBROUTINE CleanUp()

      ! #DES: Teardown the calculation by deallocating arrays and closing the log file.

      USE Input, ONLY : DeallocateInputArrays
      USE Data, ONLY  : DeallocateDataArrays
      USE Log, ONLY   : EndLogFile

      IMPLICIT NONE

        CALL DeallocateInputArrays()
        CALL DeallocateDataArrays()
        CALL EndLogFile("Normal termination of Fepcat:")

    ENDSUBROUTINE CleanUp

END MODULE Util
