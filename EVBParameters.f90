MODULE EVBParameters

  IMPLICIT NONE

  ! Should be input values
  REAL(8), PARAMETER :: alphaScale = 10.0d0
  REAL(8), PARAMETER :: aScale     = 50.0d0
  REAL(8), PARAMETER :: muScale    = 0.0001d0
  REAL(8), PARAMETER :: etaScale   = 0.000001d0
  CHARACTER(8), PARAMETER :: couplingTypes(3) = (/"CONSTANT","EXPONENT","GAUSSIAN"/)

  CONTAINS

    INTEGER FUNCTION countParams(optAlpha,optCoupling,couplingType)

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: optAlpha, optCoupling
      CHARACTER, INTENT(IN) :: couplingType

      countParams = 0
      IF (optAlpha) countParams = countParams + 1
      IF (optCoupling) THEN
        IF (TRIM(ADJUSTL(couplingType)) == "CONSTANT") THEN
          countParams = countParams + 1
        ELSE
          countParams = countParams + 2
        ENDIF
      ENDIF

    END FUNCTION CountParams

!*

    LOGICAL FUNCTION isCouplingTypeSupported(type)

      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: type
      INTEGER :: i
      
      isCouplingTypeSupported = .FALSE.
      DO i = 1, 3
        IF (TRIM(ADJUSTL(type)) == couplingTypes(i)) isCouplingTypeSupported = .TRUE.
      ENDDO

    END FUNCTION isCouplingTypeSupported

!*

    SUBROUTINE OptimizeEVBParameters(logUnit,optAlpha,optCoupling,couplingType)

      ! #Des: Guess EVB parameters and then refine them
      ! The EVB parameters are: alpha(nstates), A(nStates,nStates), mu (nStates,nStates) and eta(nStates,nStates)
      ! The optimizer needs a sensible initial guess in order to find a physically relevant minimum.

      USE DownhillSimplex, ONLY : RunNelderMead

      USE Input, ONLY : alpha, couplingConstant, mask, nBins, minPop, dGTS, dGPS
      USE Data, ONLY : mappingEnergies, recomputeDependentData !, EnergyGap, groundStateEnergy 
      USE FreeEnergy, ONLY : ComputeFEPProfile

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: logUnit
      LOGICAL, INTENT(IN) :: optAlpha, optCoupling
      CHARACTER(*), INTENT(IN) :: couplingType
      REAL(8) :: localAlpha(2), localA(2,2) !, localMu(2,2), localEta(2,2)
      REAL(8) :: profile(SIZE(mappingEnergies,2))
      REAL(8), ALLOCATABLE :: guess(:), scale(:) ! to pass into the optimizer
      INTEGER :: nParams, offset

      IF (isCouplingTypeSupported(couplingType) .EQV. .FALSE.) STOP "Error: Illegal value of couplingType passed to OptimizeEVBParameters"
      nParams = countParams(optAlpha,optCoupling,couplingType)
      IF (nParams <= 0) STOP "Error: Nothing to do in OptimizeEVBParameters"

      ALLOCATE(guess(nParams)); ALLOCATE(scale(nParams))

      WRITE(logUnit,*) "Auto-Determine EVB Parameters"

      localAlpha(:) = alpha(:)
      localA(:,:)   = couplingConstant(:,:)

      offset = 0
      IF (optAlpha .EQV. .TRUE.) THEN
        offset = 1
        CALL RecomputeDependentData(localAlpha,localA)
        CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:,1),mask(:,:),profile=profile)
        localAlpha(2) = localAlpha(2) + (dGPS - profile(SIZE(profile)))
        guess(1) = localAlpha(2)
        scale(1) = alphaScale
      ENDIF

      IF (optCoupling) THEN
        SELECT CASE (TRIM(ADJUSTL(couplingType)))
          CASE ("CONSTANT")
            guess(offset+1) = 100.0d0
            scale(offset+1) = aScale
          CASE ("EXPONENT")
            guess(offset+1) = 100.0d0; guess(offset+2) = 0.00001d0
            scale(offset+1) = aScale;  scale(offset+2) = muScale
          CASE ("GAUSSIAN")
            guess(offset+1) = 100.0d0; guess(offset+2) = 0.00001d0
            scale(offset+1) = aScale;  scale(offset+2) = etaScale
        END SELECT
      ENDIF

      CALL RunNelderMead(guess,scale,6,.TRUE.)

    ENDSUBROUTINE OptimizeEVBParameters

!*

    SUBROUTINE GuessEVBParameters(mappingEnergies,energyGap,groundStateEnergy,mask,nBins,minPop,dGTS,dGPS,alpha,mu,A) !data and targets

      USE Data, ONLY : ComputeGroundStateEnergy
      USE FreeEnergy, ONLY : ComputeFEPProfile, FEPUs, histogram
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

      !
      CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:),mask(:,:),profile=profile)
      CALL ApproximateEVBAlphas(profile,dGPS,alpha(:))

      ! Now get the FEP/US free energy barrier with zero-coupling, but WITH THE ALPHAS
      CALL ComputeGroundStateEnergy(alpha)
      ! Currently forced to get the 2D and 1D FEP/US results, but only want the 1D (binG)

      CALL Histogram(energyGap,mask,Nbins,binPopulations,binIndices,binMidpoints)
      CALL FepUS(mappingEnergies(:,:,:),groundStateEnergy,profile,binPopulations,binIndices,PMF2Dout=dGg,PMF1D=binG,minPop=minPop)
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
      REAL(8), INTENT(IN)  :: dGPS, fepProfile(:)
      REAL(8), INTENT(OUT) :: alpha(2)

      alpha(:) = 0.0d0

      ! alpha is determined just as a constant difference between desired and experimental values
      alpha(1) = 0.0d0
      alpha(2) = dGPS - fepProfile(SIZE(fepProfile))

    END SUBROUTINE ApproximateEVBAlphas

END MODULE EVBParameters
