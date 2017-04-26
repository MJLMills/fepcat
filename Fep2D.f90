PROGRAM Fep2D

  ! #DES: Performs the FEP/US procedure for an n-dimensional reaction coordinate
  !       composed of geometric collective variables.

  IMPLICIT NONE

  CALL Startup()
  CALL Run2D()
  CALL CleanUp()  

  CONTAINS

!*

    SUBROUTINE Startup()

      USE Log,   ONLY : logUnit, createLogFile
      USE Data,  ONLY : ComputeDerivedData
      USE Input, ONLY : ReadInput
      USE InputCollectiveVariables, ONLY : DetermineCollectiveVariables

      IMPLICIT NONE
      LOGICAL, PARAMETER :: readCoords = .TRUE.

      CALL CreateLogFile()
      CALL ReadInput(readCoords=readCoords)
      CALL DetermineCollectiveVariables(logUnit)
      CALL ComputeDerivedData(logUnit,doTiming=.FALSE.,readCoords=readCoords,doFEPUS=.TRUE.)

    END SUBROUTINE Startup

!*

    SUBROUTINE Cleanup()

      USE Input, ONLY : DeallocateInputArrays
      USE Data,  ONLY : DeallocateDataArrays
      USE Log,   ONLY : EndLogFile
      USE InputCollectiveVariables, ONLY : DeallocateArrays

      IMPLICIT NONE

        CALL DeallocateInputArrays()
        CALL DeallocateDataArrays()
        CALL DeallocateArrays()
        CALL EndLogFile("Normal termination of Fepcat:")

    END SUBROUTINE Cleanup

!*

    SUBROUTINE Run2D()

      USE Analysis, ONLY : Fepus2D
      USE FIleIO,   ONLY : OpenFile, CloseFile
      USE Data,     ONLY : geomRC, mappingEnergies, groundStateEnergy
      USE Input,    ONLY : mask, minPop, nBins

      IMPLICIT NONE
      INTEGER, PARAMETER :: outUnit = 19
      CHARACTER(11), PARAMETER :: outFileName = "fepus2D.csv"

      CALL OpenFile(outUnit,outFileName,"write")
      CALL Fepus2D(geomRC(:,:,:),groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),nBins,minPop,outUnit)
      CALL CloseFile(outUnit)

    END SUBROUTINE Run2D

END PROGRAM Fep2D
