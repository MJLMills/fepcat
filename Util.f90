MODULE Util

  CONTAINS

    SUBROUTINE Startup()

      ! #DES: Setup the calculation by creating a log file and reading the necessary input.

      USE Log, ONLY   : CreateLogFile
      USE Input, ONLY : ReadInput

      IMPLICIT NONE

        CALL CreateLogFile()
        CALL ReadInput()

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
