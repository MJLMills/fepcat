PROGRAM Qfep

  ! #DES: Program to produce and write the FEP/US calculation data in Qfep format. (SIC throughout)

  USE Util, ONLY : Startup, Cleanup
  USE Data, ONLY : ComputeDerivedData
  USE Log,  ONLY : logUnit

  IMPLICIT NONE
  INTEGER, PARAMETER :: outUnit = 64

    CALL Startup()
    CALL ComputeDerivedData(logUnit,doTiming=.FALSE.)
    CALL QFepAnalysis()
    CALL CleanUp()

  CONTAINS

    !Run the 3-stage analysis of the simulation data as in qfep
    SUBROUTINE QfepAnalysis()

      ! #DES: Master routine: create a Qfep log file, do the analyses and write the output.

      USE FileIO, ONLY : OpenFile, CloseFile
      USE Log, ONLY : logUnit
      IMPLICIT NONE

      WRITE(logUnit,'(A)') "Performing qfep analysis - output to QFEPstyle.log"
      WRITE(logUnit,*)

      CALL OpenFile(outUnit,"QFEPstyle.log","write")

      WRITE(outUnit,'(A)') "Qfep-style output generated by Fepcat (NB: small numerical differences with Qfep output may be observed)"
      WRITE(outUnit,*) ""
      CALL WriteParameters()
      CALL AnalyzeDynamics_Q()
      CALL AnalyzeFEP_Q()
      CALL AnalyzeFepUs_Q()

      CALL CloseFile(outUnit)

    END SUBROUTINE QfepAnalysis

!*

    SUBROUTINE WriteParameters()

      USE Input, ONLY : nFepSteps, stateA, stateB, nStates, nSkip, nBins, minPop, rcCoeffA, rcCoeffB, alpha, beta
      IMPLICIT NONE

      WRITE(outUnit,'(A27)')       "--> Number of energy files:"
      WRITE(outUNit,'(A20,I4)')    "# Number of files = ", nFepSteps
      WRITE(outUnit,'(A55)')       "--> No. of states, no. of predefined off-diag elements:"
      WRITE(outUnit,'(A21,I4)')    "# Number of states = ", nStates
      WRITE(outUnit,'(A36,I4)')    "# Number of off-diagonal elements = ", 0
      WRITE(outUnit,'(A54)')       "--> Give kT & no, of pts to skip   & calculation mode:"
      WRITE(outUnit,'(A7,F7.3)')   "# kT = ", 1.0d0 / beta
      WRITE(outUnit,'(A34,I6)')    "# Number of data points to skip = ", nSkip
      WRITE(outUnit,'(A44,L1)')    "# Only QQ interactions will be considered = ", .FALSE.
      WRITE(outUnit,'(A28)')       "--> Give number of gap-bins:"
      WRITE(outUnit,'(A23,I4)')    "# Number of gap-bins = ", nBins
      WRITE(outUnit,'(A27)')       "--> Give minimum # pts/bin:"
      WRITE(outUnit,'(A37,I4)')    "# Minimum number of points per bin = ", minPop
      WRITE(outUnit,'(A27)')       "--> Give alpha for state 2:"
      WRITE(outUnit,'(A22,F7.2)')  "# Alpha for state 2 = ", alpha(stateB) - alpha(stateA)
      WRITE(outUnit,'(A16)') "-->  Hij scaling:"
      WRITE(outUNit,'(A25,F7.2)')  "# Scale factor for Hij = ", 1.0d0
      WRITE(outUnit,'(A57)')       "--> linear combination of states defining reaction coord:"
      WRITE(outUnit,'(A37,2F7.2)') "# Linear combination co-efficients = ", rcCoeffA, rcCoeffB

      WRITE(outUnit,*); WRITE(outUnit,*)

    END SUBROUTINE WriteParameters

!*

    SUBROUTINE AnalyzeDynamics_Q()

      ! #DES: Reproduce the output simulation data table of Q (same data, slightly different format)

      USE Input, ONLY : energyNames, stateEnergy, coeffs, mask, CreateFileNames, fileBase, nFepSteps
      USE StatisticalFunctions, ONLY : mean

      IMPLICIT NONE
      INTEGER :: name, fepstep, state, type
      CHARACTER(500) :: fileNames(nFepSteps)

      CALL CreateFileNames('en',fileBase,fileNames)

      WRITE(outUnit,'(A)') "# Part 0: Average energies for all states in all files"
      WRITE(outUnit,'(A9,1X,A5,1X,A8,1X,A9)',ADVANCE='NO') "FEP Index", "State", "Points", "Coeff"

      DO name = 1, SIZE(energyNames) !nEnergyTypes
        WRITE(outUnit,'(A11)',ADVANCE='NO') energyNames(name)
      ENDDO
      WRITE(outUnit,*)

      DO fepstep = 1, SIZE(stateEnergy,2) !nFepSteps
        DO state = 1, SIZE(stateEnergy,3) !nStates
          WRITE(outUnit,'(A,1X,I5,1X,I8,1X,F9.6)',ADVANCE='NO') TRIM(ADJUSTL(fileNames(fepstep))), state, COUNT(mask(fepstep,:)), coeffs(1,fepstep,state)
          DO type = 1, SIZE(energyNames)
            WRITE(outUnit,'(F11.2)',ADVANCE='NO') mean( stateEnergy(:,fepstep,state,type),mask=mask(fepstep,:) )
          ENDDO
          WRITE(outUnit,*)
        ENDDO
        WRITE(outUnit,*)
      ENDDO

    ENDSUBROUTINE AnalyzeDynamics_Q

!*

    ! Compute and print the output of the FEP procedure in Q format
    SUBROUTINE AnalyzeFep_Q()

      USE Data, ONLY : mappingEnergies
      USE Input, ONLY : stateA, coeffs, nFepSteps, mask
      USE FreeEnergy, ONLY : ComputeFEPIncrements

      IMPLICIT NONE
      REAL(8) :: forward(nFepSteps-1), reverse(nFepSteps-1), profile(nFepSteps-1)
      REAL(8) :: deltaG_f(nFepSteps), deltaG_r(nFepSteps), deltaG(nFepSteps)
      INTEGER :: fepstep

      CALL ComputeFEPIncrements(1,nFepSteps,mappingEnergies(:,:,:,1),mask(:,:),forward,reverse,profile)

      WRITE(outUnit,'(A43)') "# Part 1: Free energy perturbation summary:"
      WRITE(outUnit,*) ""
      WRITE(outUnit,'(A)') "# Calculation for full system"
      WRITE(outUnit,'(A)') "# lambda(1)      dGf sum(dGf)      dGr sum(dGr)     <dG>"

      deltaG_f(:) = 0.0d0
      deltaG(:) = 0.0d0
      !1st value is 0 in each
      DO fepstep = 2, nFepSteps
        deltaG_f(fepstep) = forward(fepstep-1)
        deltaG(fepstep)   = profile(fepstep-1)
      ENDDO

      !final value is 0
      deltaG_r(:) = 0.0d0
      DO fepstep = 1, nFepSteps - 1
        deltaG_r(fepstep) = reverse(fepstep)
      ENDDO

      DO fepstep = 1, nFepSteps
        WRITE(outUnit,'(2X,F9.6,5F9.3)') coeffs(1,fepstep,stateA), deltaG_f(fepstep), SUM(deltaG_f(1:fepstep)) , &
                                         &                         deltaG_r(fepstep), SUM(deltaG_r(fepstep:nFepSteps)), &
                                         &                         SUM(deltaG(1:fepstep))
      ENDDO

      WRITE(outUnit,*)

    ENDSUBROUTINE AnalyzeFep_Q

!*

   ! Compute and print the output of the FEP/US procedure in Q format
   SUBROUTINE AnalyzeFepUs_Q()

     USE Data,       ONLY : energyGap, mappingEnergies, groundStateEnergy, lambda
     USE Input,      ONLY : Nbins, minPop, stateA, stateB, stateEnergy, mask
     USE FreeEnergy, ONLY : Histogram, ComputeFepIncrements, FepUS

     IMPLICIT NONE
     INTEGER :: bin, step, binPopulations(Nbins,SIZE(energyGap,1)), binIndices(SIZE(energyGap,1),SIZE(energyGap,2))
     REAL(8) :: binMidpoints(Nbins), dGg(Nbins,SIZE(energyGap,1)), dGa(Nbins,SIZE(energyGap,1)), dGb(Nbins,SIZE(energyGap,1))
     REAL(8) :: binG(Nbins)
     REAL(8) :: dG_FEP(SIZE(energyGap,1)-1), G_FEP(SIZE(energyGap,1))

     INTEGER :: fepstep

     WRITE(outUnit,'(A20,F7.2)') "# Min energy-gap is:", MINVAL(energyGap(:,:))
     WRITE(outUnit,'(A20,F7.2)') "# Max energy-gap is:", MAXVAL(energyGap(:,:))
     WRITE(outUnit,*); WRITE(outUnit,*); WRITE(outUnit,*)
     WRITE(outUnit,'(A39)') "# Part 2: Reaction free energy summary:"
     WRITE(outUnit,'(A79)') "# Lambda(1)  bin Energy gap      dGa     dGb     dGg    # pts    c1**2    c2**2"

     CALL Histogram(energyGap,mask,Nbins,binPopulations,binIndices,binMidpoints)

     CALL ComputeFEPIncrements(1,SIZE(energyGap,1),mappingEnergies(:,:,:,1),mask(:,:),profile=dG_FEP)
     G_FEP(:) = 0.0d0
     DO fepstep = 2, SIZE(energyGap,1)
       G_FEP(fepstep) = SUM(dG_FEP(1:fepstep-1))
     ENDDO

     ! energyGap is the reaction coordinate, nBins is num histogram bins, binPop is histogram values,
     ! indices is which bin each point is in, binMidpoints is x values

     CALL FepUS(mappingEnergies(:,:,:,1),stateEnergy(:,:,stateA,1),G_FEP,binPopulations,binIndices,PMF2D=dGa,PMF1D=binG)
     CALL FepUS(mappingEnergies(:,:,:,1),stateEnergy(:,:,stateB,1),G_FEP,binPopulations,binIndices,PMF2D=dGb,PMF1D=binG)
     CALL FepUS(mappingEnergies(:,:,:,1),groundStateEnergy,        G_FEP,binPopulations,binIndices,PMF2D=dGg,PMF1D=binG)

     DO step = 1, SIZE(energyGap,1)
       DO bin = 1, Nbins
         IF (binPopulations(bin,step) >= minPop) WRITE(outUnit,'(2X,F9.6,I5,2X,4F9.2,2X,I5)') lambda(step), bin, binMidpoints(bin), dGa(bin,step),dGb(bin,step), dGg(bin,step), binPopulations(bin,step)
       ENDDO
     ENDDO

     WRITE(outUnit,*); WRITE(outUnit,*)
     WRITE(outUnit,'(A30)') "# Part 3: Bin-averaged summary"
     WRITE(outUnit,'(A63)') "# bin  energy gap  <dGg> <dGg norm> pts  <c1**2> <c2**2> <r_xy>"

     DO bin = 1, nBins
       IF (SUM(binPopulations(bin,:)) > 0) WRITE(outUnit,'(I4,1X,3F9.2,2X,I5)') bin, binMidpoints(bin), binG(bin), binG(bin)-MINVAL(binG), SUM(binPopulations(bin,:))
     ENDDO

    ENDSUBROUTINE AnalyzeFepUs_Q

END PROGRAM Qfep

