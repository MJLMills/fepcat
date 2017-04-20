PROGRAM Fep2D

  USE Log, ONLY  : logUnit, createLogFile

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

      CALL CreateLogFile()
      CALL ReadInput(readCoords=.TRUE.)
      CALL DetermineCollectiveVariables(logUnit)
      CALL ComputeDerivedData(logUnit,doTiming=.FALSE.,readCoords=.TRUE.)

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

      USE Data,     ONLY : geomRC, mappingEnergies, groundStateEnergy
      USE Input,    ONLY : mask, minPop, nBins
      USE Analysis, ONLY : Fepus2D

      IMPLICIT NONE

      CALL Fepus2D(geomRC,groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),nBins,minPop,logUnit)

    END SUBROUTINE Run2D

END PROGRAM Fep2D
