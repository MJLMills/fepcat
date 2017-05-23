PROGRAM Fep

  ! #DES: Program to only do FEP analysis of data

  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log,  ONLY : logUnit

  IMPLICIT NONE

    CALL Startup(readCoords=.FALSE.)
    CALL ComputeDerivedData(logUnit,doTiming=.FALSE.,readCoords =.FALSE.,doFEPUS=.FALSE.)
    CALL FepAnalysis()
    CALL CleanUp()  

  CONTAINS

  SUBROUTINE FepAnalysis

    USE Input, ONLY    : mask, energynames, outDir, stateEnergy
    USE Data, ONLY     : lambda, mappingEnergies
    USE Analysis, ONLY : WriteMeanEnergyBreakdown, FepBreakdown
    USE FileIO, ONLY : OpenFile, CloseFile

    IMPLICIT NONE
    INTEGER, PARAMETER :: outUnit = 22

      CALL OpenFile(outUnit,TRIM(ADJUSTL(outDir))//"/mean-energy.csv",'write')
      CALL WriteMeanEnergyBreakdown(stateEnergy,mask,energyNames,lambda,outUnit)
      CALL CloseFile(outUnit)

      CALL OpenFile(outUnit,TRIM(ADJUSTL(outDir))//"/fep-breakdown.csv","write")
      CALL FepBreakdown(lambda(:),mappingEnergies(:,:,:,:),mask(:,:),energyNames(:),outUnit)
      CALL CloseFile(outUnit)

  END SUBROUTINE FepAnalysis

END PROGRAM Fep
