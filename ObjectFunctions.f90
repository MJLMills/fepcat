
MODULE ObjectFunctions

  IMPLICIT NONE

  LOGICAL      :: optAlpha     = .TRUE.
  LOGICAL      :: optCoupling  = .TRUE.
  CHARACTER(8) :: couplingType = "CONSTANT"

  CONTAINS

  REAL(8) FUNCTION objectFunction(variables)

    ! #DES: This is the part that changes per-application. I am massively cheating by 
    !       just importing eveything here for now.

    USE Input, ONLY : dGTS, dGPS, mask, Nbins, minPop
    USE Data, ONLY : energyGap, groundStateEnergy, mappingEnergies, RecomputeDependentData
    USE Analysis, ONLY : FepUsFreeEnergies
    USE EVBParameters, ONLY : DelinearizeParameters
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

END MODULE ObjectFunctions
