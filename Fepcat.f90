PROGRAM Fepcat

  ! #DES: Program for free energy perturbation (FEP) calculations on molecular dynamics (MD) trajectories.
  ! #DES: For an empirical valence bond (EVB) potential of nStates states, compute the free energy change from stateA to stateB.
 
  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log, ONLY  : logUnit

  IMPLICIT NONE

    CALL Startup(readCoords=.FALSE.)
    CALL ComputeDerivedData(logUnit,doTiming=.FALSE., readCoords = .FALSE., doFEPUS = .TRUE.)
    CALL FullAnalysis()
    CALL CleanUp()

  CONTAINS

!*

    SUBROUTINE FullAnalysis()

      USE Input,    ONLY : energyNames, mask, stateEnergy, coeffs, skip, stateA, stateB, Nbins, minPop, outDir
      USE Data,     ONLY : mappingEnergies, lambda, groundStateEnergy, energyGap
      USE Analysis, ONLY : WriteMeanEnergyBreakdown, AnalyzeSimulationConvergence, FepBreakdown, FepUsGroundState, FepUsFreeEnergies
      USE FileIO,   ONLY : OpenFile, CloseFile
!      USE EVBParameters, ONLY : OptimizeEVBParameters

      IMPLICIT NONE
      INTEGER, PARAMETER :: outUnit = 20

      CALL OpenFile(outUnit,TRIM(ADJUSTL(outDir))//"/mean-energy.csv",'write')
      CALL WriteMeanEnergyBreakdown(stateEnergy,mask,energyNames,lambda,outUnit)
      CALL CloseFile(outUnit)

      CALL OpenFile(outUnit,TRIM(ADJUSTL(outDir))//"/convergence.csv",'write')
      CALL AnalyzeSimulationConvergence(mappingEnergies(:,:,:,12),mask(:,:),lambda(:),coeffs(1,:,stateA:stateB),skip(:),outUnit)
      CALL CloseFile(outUnit)

      CALL OpenFile(outUnit,TRIM(ADJUSTL(outDir))//"/fep-breakdown.csv","write")
      CALL FepBreakdown(lambda(:),mappingEnergies(:,:,:,:),mask(:,:),energyNames(:),outUnit)
      CALL CloseFile(outUnit)

      CALL OpenFile(outUnit,TRIM(ADJUSTL(outDir))//"/fepus-groundstate.csv",'write')
      CALL FepUsGroundState(energyGap(:,:),groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),Nbins,minPop,outUnit)
      CALL CloseFile(outUnit)

      CALL FepUsFreeEnergies(energyGap(:,:),groundStateEnergy(:,:),mappingEnergies(:,:,:,1),mask(:,:),Nbins,minPop,.TRUE.,6)

!      CALL OptimizeEVBParameters(6,.TRUE.,.TRUE.,"CONSTANT")

    END SUBROUTINE FullAnalysis

END PROGRAM Fepcat
