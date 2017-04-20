PROGRAM Fep2D

  USE Util, ONLY : Cleanup
  USE Log, ONLY  : logUnit, createLogFile

  IMPLICIT NONE

  CALL CreateLogFile()
  CALL Run2D()
  CALL CleanUp()  

  CONTAINS

    SUBROUTINE Run2D()

      USE Data, ONLY : geomRC, mappingEnergies, groundStateEnergy, ComputeDerivedData
      USE Input, ONLY : ReadInput, mask, minPop, nBins
      USE Input2D, ONLY : ReadReactionCoordinates, WriteReactionCoordinates
      USE Analysis, ONLY : Fepus2D
      USE Log, ONLY : logUnit
      IMPLICIT NONE
      LOGICAL, PARAMETER :: readCoords = .TRUE.

      ! first gotta read all the trajectory data
      CALL ReadInput(readCoords)
      ! then we need to know what the reaction coordinate collective variables are
      CALL ReadReactionCoordinates()
      CALL WriteReactionCoordinates(logUnit)

      ! and their values need to be computed
      CALL ComputeDerivedData(logUnit,doTiming=.FALSE.,readCoords = .TRUE.)

      CALL Fepus2D(geomRC,groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),nBins,minPop,logUnit)

    END SUBROUTINE Run2D

END PROGRAM Fep2D
