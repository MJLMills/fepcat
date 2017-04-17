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

      USE Input,    ONLY : energyNames, mask, stateEnergy, coeffs, skip, stateA, stateB, Nbins, minPop
      USE Data,     ONLY : mappingEnergies, lambda, groundStateEnergy, energyGap
      USE Analysis, ONLY : WriteMeanEnergyBreakdown, AnalyzeSimulationConvergence, FepBreakdown, FepUsGroundState, FepUsFreeEnergies
      USE FileIO,   ONLY : OpenFile, CloseFile
      USE DownhillSimplex, ONLY : RunNelderMead

      IMPLICIT NONE
      INTEGER, PARAMETER :: outUnit = 20
      REAL(8) :: guess(2), scale(2)

      CALL OpenFile(outUnit,"mean-energy.csv",'write')
      CALL WriteMeanEnergyBreakdown(stateEnergy,mask,energyNames,lambda,outUnit)
      CALL CloseFile(outUnit)

      CALL OpenFile(outUnit,"convergence.csv",'write')
      CALL AnalyzeSimulationConvergence(mappingEnergies(:,:,:,1),mask(:,:),lambda(:),coeffs(1,:,stateA:stateB),skip(:),outUnit)
      CALL CloseFile(outUnit)

      CALL OpenFile(outUnit,"fep-breakdown.csv","write")
      CALL FepBreakdown(lambda(:),mappingEnergies(:,:,:,:),mask(:,:),energyNames(:),outUnit)
      CALL CloseFile(outUnit)

      CALL OpenFile(outUnit,"fepus-groundstate.csv",'write')
      CALL FepUsGroundState(energyGap(:,:),groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),Nbins,minPop,outUnit)
      CALL CloseFile(outUnit)

      CALL FepUsFreeEnergies(energyGap(:,:),groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),Nbins,minPop,.TRUE.)

      guess(1) = -52.0d0; scale(1) = 3.0d0
      guess(2) = 100.0d0; scale(2) = 50.0d0
      CALL RunNelderMead(guess,scale,6,.TRUE.)

    END SUBROUTINE FullAnalysis

END PROGRAM Fepcat
