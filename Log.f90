MODULE Log

  ! #DES: Module for interacting with the log file

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CreateLogFile, EndLogFile, logUnit

  INTEGER, PARAMETER :: logUnit = 110

  CONTAINS

    SUBROUTINE CreateLogFile()

      ! #DES: Creates a new log file 'logFile' for writing and writes its head

      USE FileIO, ONLY : OpenFile      
      IMPLICIT NONE

      CALL OpenFile(logUnit,"logFile","write")
      CALL WriteLogHead()

    END SUBROUTINE CreateLogFile

!*

    SUBROUTINE EndLogFile(message)

      ! #DES: Writes 'message' to the log file, prints the tail and attempts to close it

      USE FileIO, ONLY : CloseFile
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: message

      WRITE(logUnit,'(a)',ADVANCE='NO') TRIM(ADJUSTL(message))
      CALL WriteLogTail()
      CALL CloseFile(logUnit)

    END SUBROUTINE EndLogFile

!*

    SUBROUTINE WriteLogHead()

      ! #DES: Writes calculation-agnostic header to the log file

      IMPLICIT NONE

      WRITE(logUnit,*) ""
      WRITE(logUnit,*) "*             FEPCAT             *"
      WRITE(logUnit,*) "* A Free Energy Calculation Tool *"
      WRITE(logUnit,*) "         Matthew JL Mills         "
      WRITE(logUnit,*) "     Joint BioEnergy Institute    "
      WRITE(logUnit,*) "         LICENSE GOES HERE        "
      WRITE(logUnit,*) ""

    END SUBROUTINE WriteLogHead

!*

    SUBROUTINE WriteLogTail()

      ! #DES: Writes calculation-agnostic tail to the log file (time and date of completion)

      IMPLICIT NONE
      CHARACTER(8)  :: date
      CHARACTER(10) :: time

      CALL DATE_AND_TIME(date,time)
      WRITE(logUnit,'(1X,A,2X,A,2X)') date, time

    END SUBROUTINE WriteLogTail   

END MODULE Log
