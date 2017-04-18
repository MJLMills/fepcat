MODULE EVBParameters

  IMPLICIT NONE

  CONTAINS

!*
    SUBROUTINE OptimizeEVBParameters(paramMask,nStates,logUnit)

      ! #Des: Guess EVB parameters and then refine them
      ! The EVB parameters are: alpha(nstates), A(nStates,nStates), mu (nStates,nStates) and eta(nStates,nStates)
      ! The optimizer needs a sensible initial guess in order to find a physically relevant minimum.

      USE DownhillSimplex, ONLY : RunNelderMead
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit, nStates
      LOGICAL, INTENT(IN) :: paramMask(:)

      REAL(8), PARAMETER :: alphaScale = 10.0d0
      REAL(8), PARAMETER :: aScale     = 50.0d0
      REAL(8), PARAMETER :: muScale    = 0.0001d0
      REAL(8), PARAMETER :: etaScale   = 0.000001d0

      CHARACTER(8) :: paramNames(4*nStates), charState
      REAL(8), ALLOCATABLE :: guess(:), scale(:)
      INTEGER :: nParams, param, state

      DO state = 1, nStates
        WRITE(charState,'(I0.2)') state
        paramNames(state) = "alpha_"//charState
        paramNames(nStates+state) = "A_"//charState
        paramNames((2*nStates)+state) = "mu_"//charState
        paramNames((3*nStates)+state) = "eta_"//charState
      ENDDO

      nParams = 0
      DO param = 1, SIZE(paramMask)
        WRITE(logUnit,*) paramNames(param), paramMask(param)
        IF (paramMask(param) .EQV. .TRUE.) nParams = nParams + 1
      ENDDO
      WRITE(logUnit,*) "Number of Optimization Parameters: ", nParams
      ALLOCATE(guess(nParams))
      ALLOCATE(scale(nParams))
      
      ! generate the guess values
      DO state = 1, nStates
      ENDDO

      guess(1) = -52.0d0; scale(1) = 10.0d0
      guess(2) = 100.0d0; scale(2) = 50.0d0

      CALL RunNelderMead(guess,scale,6,.TRUE.)

      IF (ALLOCATED(guess)) DEALLOCATE(guess)
      IF (ALLOCATED(scale)) DEALLOCATE(scale)

    ENDSUBROUTINE OptimizeEVBParameters

!*

    SUBROUTINE GuessEVBParameters(mappingEnergies,energyGap,groundStateEnergy,mask,nBins,minPop,dGTS,dGPS,alpha,mu,A) !data and targets

      USE Data, ONLY : ComputeGroundStateEnergy
      IMPLICIT NONE

      REAL(8), INTENT(IN) :: mappingEnergies(:,:,:), energyGap(:,:), groundStateEnergy(:,:)
      REAL(8), INTENT(IN) :: dGPS, dGTS
      LOGICAL, INTENT(IN) :: mask(:,:)
      INTEGER, INTENT(IN) :: nBins, minPop
      REAL(8), INTENT(OUT) :: alpha(2), mu(2,2), A(2,2)
      ! For the histogram
      INTEGER :: binPopulations(Nbins,SIZE(energyGap,1)), binIndices(SIZE(energyGap,1),SIZE(energyGap,2))
      REAL(8) :: binMidpoints(Nbins)
      ! Data needed from subroutines
      REAL(8) :: dGg(Nbins,SIZE(energyGap,1)), binG(Nbins)
      REAL(8) :: profile(SIZE(mappingEnergies,2))
      LOGICAL :: binMask(Nbins)
      INTEGER :: bin

      CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:),mask(:,:),profile=profile)
      CALL ApproximateEVBAlphas(profile,dGPS,alpha(:))

      ! Now get the FEP/US free energy barrier with zero-coupling, but WITH THE ALPHAS
      CALL ComputeGroundStateEnergy(alpha)
      ! Currently forced to get the 2D and 1D FEP/US results, but only want the 1D (binG)

      CALL Histogram(energyGap,mask,Nbins,binPopulations,binIndices,binMidpoints)
      CALL FepUS(mappingEnergies(:,:,:),groundStateEnergy,profile,binPopulations,binIndices,PMF2D=dGg,PMF1D=binG,minPop=minPop)
      binMask(:) = .TRUE.
      DO bin = 1, Nbins
        IF (SUM(binPopulations(bin,:)) < minPop) binMask(bin) = .FALSE.
      ENDDO
      ! And pass it in to get the A/mu values
      CALL ApproximateEVBCoupling(binMidpoints,binG,binMask,dGTS,mu,A)

    END SUBROUTINE GuessEVBParameters

!*

    SUBROUTINE ApproximateEVBCoupling(binMidPoints,fepusProfile,binMask,dGTS,mu,A)

      USE Input, ONLY : nStates, stateA, stateB
      USE StatisticalFunctions, ONLY : variance
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: fepusProfile(:), binMidpoints(:)
      REAL(8), INTENT(IN) :: dGTS
      LOGICAL, INTENT(IN) :: binMask(:)
      REAL(8), INTENT(OUT) :: mu(2,2), A(2,2)
      REAL(8) :: c

      mu(:,:) = 0.0d0
      A(:,:) = 0.0d0

      c = SQRT(variance(binMidpoints,binMask)) / 2.0d0
      mu(stateA,stateB) = 1.0d0 / (2.0d0 * c * c)
      A(stateA,stateB) = MAXVAL(fepusProfile) - dGTS

    END SUBROUTINE ApproximateEVBCoupling

!*

    SUBROUTINE ApproximateEVBAlphas(fepProfile,dGPS,alpha)

      ! #DES: Determine approximate values of the EVB parameters from provided target dG values.
      ! Only works for relative parameters.

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: dGPS, fepProfile(:)
      REAL(8), INTENT(OUT) :: alpha(2)

      alpha(:) = 0.0d0

      ! alpha is determined just as a constant difference between desired and experimental values
      alpha(1) = 0.0d0
      alpha(2) = dGPS - fepProfile(SIZE(fepProfile))

      WRITE(*,*) "dGPS_EXP = ", dGPS
      WRITE(*,*) "dGPS_FEP = ", fepProfile(SIZE(fepProfile))
      WRITE(*,*) "Alpha(2) = ", alpha(2)

    END SUBROUTINE ApproximateEVBAlphas

END MODULE EVBParameters
