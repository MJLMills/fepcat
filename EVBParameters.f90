MODULE EVBParameters

  IMPLICIT NONE

  REAL(8) :: aScale = 50.0d0, alphaScale = 10.0d0, muScale = 0.00001d0, etaScale = 0.00001d0

  CHARACTER(8), PARAMETER :: couplingTypes(3) = (/"CONSTANT","EXPONENT","GAUSSIAN"/)
  LOGICAL :: optAlpha, optCoupling

  CONTAINS

!*

    SUBROUTINE OptimizeEVBParameters(logUnit,optAlpha,optCoupling,couplingType)

      ! #Des: Guess EVB parameters and then refine them
      ! The EVB parameters are: alpha(nstates), A(nStates,nStates), mu (nStates,nStates) and eta(nStates,nStates)
      ! The optimizer needs a sensible initial guess in order to find a physically relevant minimum.

      USE DownhillSimplex, ONLY : RunNelderMead
      USE Input, ONLY      : alpha, couplingConstant, couplingGaussExpFactor, couplingExpExpFactor, mask, nBins, minPop, dGTS, dGPS
      USE Data, ONLY       : mappingEnergies, recomputeDependentData !, EnergyGap, groundStateEnergy 
      USE FreeEnergy, ONLY : ComputeFEPProfile

      IMPLICIT NONE

      INTEGER,      INTENT(IN) :: logUnit
      LOGICAL,      INTENT(IN) :: optAlpha, optCoupling
      CHARACTER(*), INTENT(IN) :: couplingType

      REAL(8) :: localAlpha(2), localA(2,2), localMu(2,2), localEta(2,2)
      REAL(8) :: profile(SIZE(mappingEnergies,2))
      INTEGER :: nParams
      REAL(8), ALLOCATABLE :: guess(:), scale(:)

      IF (isCouplingTypeSupported(couplingType) .EQV. .FALSE.) STOP "Error: Illegal value of couplingType passed to OptimizeEVBParameters"
      nParams = countParams(optAlpha,optCoupling,couplingType); IF (nParams <= 0) STOP "Error: Nothing to do in OptimizeEVBParameters"

      ALLOCATE(guess(nParams)); ALLOCATE(scale(nParams))

      WRITE(logUnit,'(A,I0.2,A)') "Auto-Determine ",nParams ," EVB Parameters"

      localAlpha(:) = 0.0d0; localA(:,:)   = 0.0d0
      localMu(:,:)  = 0.0d0; localEta(:,:) = 0.0d0

      IF (optAlpha .EQV. .TRUE.) THEN

        CALL RecomputeDependentData(alpha,couplingConstant,couplingExpExpFactor,couplingGaussExpFactor)
        CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:,1),mask(:,:),profile=profile)
        localAlpha(1) = alpha(1)
        localAlpha(2) = localAlpha(2) + (dGPS - profile(SIZE(profile)))

      ENDIF

      IF (optCoupling) THEN

        SELECT CASE (TRIM(ADJUSTL(couplingType)))
          CASE ("CONSTANT")
            localA(:,:) = couplingConstant(:,:)
          CASE ("EXPONENT")
          CASE ("GAUSSIAN")
            localA(:,:) = couplingConstant(:,:)
            localEta(:,:) = couplingGaussExpFactor(:,:)
        END SELECT

      ENDIF

      CALL LinearizeParameters(localAlpha,localA,localMu,localEta,optAlpha,optCoupling,couplingType,params=guess,scale=scale)
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


!*

    INTEGER FUNCTION countParams(optAlpha,optCoupling,couplingType)

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: optAlpha, optCoupling
      CHARACTER(*), INTENT(IN) :: couplingType

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
      DO i = 1, SIZE(couplingTypes)
        IF (TRIM(ADJUSTL(type)) == couplingTypes(i)) isCouplingTypeSupported = .TRUE.
      ENDDO

    END FUNCTION isCouplingTypeSupported

!*

  REAL(8) FUNCTION objectFunction(variables)

    ! #DES: This is the part that changes per-application. I am massively cheating by 
    !       just importing eveything here for now.

    USE Input, ONLY : dGTS, dGPS, mask, Nbins, minPop
    USE Data, ONLY : energyGap, groundStateEnergy, mappingEnergies, RecomputeDependentData
    USE Analysis, ONLY : FepUsFreeEnergies
    USE FreeEnergy, ONLY : FepUS, ScanFepUS, Histogram, computeFEPProfile

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: variables(:)

    INTEGER :: binPopulations(Nbins,SIZE(energyGap,1)), binIndices(SIZE(energyGap,1),SIZE(energyGap,2))
    REAL(8) :: binMidpoints(Nbins)
    REAL(8) :: alpha(2), A(2,2), mu(2,2), eta(2,2)
    REAL(8) :: relative(3), binGg(Nbins)
    REAL(8) :: G_FEP(SIZE(energyGap,1))
    LOGICAL :: printBin(Nbins)

    ! For EVB-parameterisation, the variables are alpha and the parameters of the off-diagonals
    ! and the object function must depend on the FEP/US free energy changes

    CALL DelinearizeParameters(variables(:),alpha(:),A(:,:),mu(:,:),eta(:,:),optAlpha,optCoupling,couplingType)

    CALL RecomputeDependentData(alpha,A,mu,eta)
    CALL ComputeFEPProfile(1,SIZE(energyGap,1),mappingEnergies(:,:,:,1),mask(:,:),profile=G_FEP)
    CALL Histogram(energyGap(:,:),mask,Nbins,binPopulations,binIndices,binMidpoints)
    CALL FepUS(mappingEnergies(:,:,:,1),groundStateEnergy(:,:),G_FEP,binPopulations,binIndices,PMF1D=binGg,minPop=minPop,useBin=printBin)
    CALL ScanFepUs(binMidpoints(:),binGg(:),mask=printBin,stationaryPoints=relative)

    relative(:) = relative(:) - relative(1)

    objectFunction = 0.0d0
    IF (optAlpha    .EQV. .TRUE.) objectFunction = objectFunction + (relative(3) - dGPS)**2.0d0
    IF (optCoupling .EQV. .TRUE.) objectFunction = objectFunction + (relative(2) - dGTS)**2.0d0
    objectFunction = SQRT(objectFunction) ! costs time but gives object function units of energy

  END FUNCTION objectFunction

!*

    SUBROUTINE LinearizeParameters(alpha,A,mu,eta,optAlpha,optCoupling,couplingType,params,scale)

      ! Knows the rules to take the parameter arrays and convert to a vector for optimisation

      IMPLICIT NONE
      REAL(8), INTENT(IN)      :: alpha(:), A(:,:), mu(:,:), eta(:,:)
      LOGICAL, INTENT(IN)      :: optAlpha, optCoupling
      CHARACTER(*), INTENT(IN) :: couplingType

      REAL(8), INTENT(OUT)     :: params(:)
      REAL(8), INTENT(OUT), OPTIONAL :: scale(:)

      INTEGER :: offset

      offset = 0
      IF (optAlpha .EQV. .TRUE.) THEN
        offset = 1
        params(1) = alpha(2)
        IF (PRESENT(scale)) scale(1) = alphaScale
      ENDIF

      IF (optCoupling .EQV. .TRUE.) THEN

        SELECT CASE (TRIM(ADJUSTL(couplingType)))

          CASE ("CONSTANT")
            params(offset+1) = A(1,2)
            IF (PRESENT(scale)) scale(offset+1) = aScale
          CASE ("EXPONENT")
            params(offset+1) = A(1,2)
            IF (PRESENT(scale)) scale(offset+1) = aScale
            params(offset+2) = mu(1,2)
            IF (PRESENT(scale)) scale(offset+2) = muScale
          CASE ("GAUSSIAN")
            params(offset+1) = A(1,2)
            IF (PRESENT(scale)) scale(offset+1) = aScale
            params(offset+2) = eta(1,2)
            IF (PRESENT(scale)) scale(offset+2) = etaScale

        END SELECT

      ENDIF

    END SUBROUTINE LinearizeParameters

!*

   SUBROUTINE DelinearizeParameters(params,alpha,A,mu,eta,optAlpha,optCoupling,couplingType)

     ! Knows the rules to take a vector of parameters and return them to their correct positions in properly-shaped arrays

     IMPLICIT NONE
     LOGICAL, INTENT(IN) :: optAlpha, optCoupling
     CHARACTER(*), INTENT(IN) :: couplingType
     REAL(8), INTENT(IN) :: params(:)
     REAL(8), INTENT(OUT) :: alpha(2), A(2,2), mu(2,2), eta(2,2)
     INTEGER :: offset

      offset = 0
      IF (optAlpha .EQV. .TRUE.) THEN
        offset = 1
        alpha(2) = params(1)
      ENDIF

      IF (optCoupling .EQV. .TRUE.) THEN

        SELECT CASE (TRIM(ADJUSTL(couplingType)))

          CASE ("CONSTANT")
            A(1,2)   = params(offset+1)
          CASE ("EXPONENT")
            A(1,2)   = params(offset+1)
            mu(1,2)  = params(offset+2)
          CASE ("GAUSSIAN")
            A(1,2)   = params(offset+1)
            eta(1,2) = params(offset+2)

        END SELECT

      ENDIF

   END SUBROUTINE DelinearizeParameters


END MODULE EVBParameters
