PROGRAM Fepcat

  ! #DES: Program for free energy perturbation (FEP) calculations on molecular dynamics (MD) trajectories.
  ! #DES: For an empirical valence bond (EVB) potential of nStates states, compute the free energy change from stateA to stateB.
 
  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log, ONLY  : logUnit

  IMPLICIT NONE

    CALL Startup()
    CALL ComputeDerivedData(logUnit,doTiming=.FALSE.)
    CALL FullAnalysis()
    CALL CleanUp()

  CONTAINS

!*

    SUBROUTINE FullAnalysis()

      USE Input,    ONLY : energyNames, mask, stateEnergy, coeffs, skip, stateA, stateB
      USE Data,     ONLY : mappingEnergies, lambda
      USE Analysis, ONLY : WriteMeanEnergyBreakdown, AnalyzeSimulationConvergence, FepBreakdown, FepUsGroundState
      USE FileIO,   ONLY : OpenFile, CloseFile

      IMPLICIT NONE
      INTEGER, PARAMETER :: meanUnit = 20, convUnit = 21

      CALL OpenFile(meanUnit,"mean-energy.csv",'write')
      CALL WriteMeanEnergyBreakdown(stateEnergy,mask,energyNames,lambda,meanUnit)
      CALL CloseFile(meanUnit)

      CALL OpenFile(convUnit,"convergence.csv",'write')
      CALL AnalyzeSimulationConvergence(mappingEnergies(:,:,:,1),mask(:,:),lambda(:),coeffs(1,:,stateA:stateB),skip(:),convUnit)
      CALL CloseFile(convUnit)

      CALL FepBreakdown(lambda(:),mappingEnergies(:,:,:,:),mask(:,:),energyNames(:))
      CALL FepUsGroundState()

    END SUBROUTINE FullAnalysis

END PROGRAM Fepcat
