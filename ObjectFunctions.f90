MODULE ObjectFunctions

  IMPLICIT NONE

  CONTAINS

  REAL(8) FUNCTION objectFunction(variables)

    ! #DES: This is the part that changes per-application. I am massively cheating by 
    !       just importing eveything here for now.

    USE Input, ONLY : dGTS, dGPS, mask, Nbins, minPop
    USE Data, ONLY : energyGap, groundStateEnergy, mappingEnergies, RecomputeDependentData
    USE Analysis, ONLY : FepUsFreeEnergies
    USE FreeEnergy, ONLY : FepUS, ScanFepUS, Histogram, computeFEPProfile

    IMPLICIT NONE

    REAL(8), INTENT(IN) :: variables(:)
    REAL(8) :: alpha(2), relative(3), A(2,2)
    INTEGER :: binPopulations(Nbins,SIZE(energyGap,1)), binIndices(SIZE(energyGap,1),SIZE(energyGap,2))
    REAL(8) :: binMidpoints(Nbins)
    REAL(8) :: binGg(Nbins)
    REAL(8) :: G_FEP(SIZE(energyGap,1))
    LOGICAL :: printBin(Nbins)

    ! For EVB-parameterisation, the variables are alpha and the parameters of the off-diagonals
    ! and the object function must depend on the FEP/US free energy changes

    ! compute the free energy values with the current parameter set
    alpha(1) = 0.0d0; alpha(2) = variables(1)
    A(:,:) = 0.0d0; ! A(1,2) = variables(2); A(2,1) = A(1,2)

    CALL RecomputeDependentData(alpha,A)
    CALL ComputeFEPProfile(1,SIZE(energyGap,1),mappingEnergies(:,:,:,1),mask(:,:),profile=G_FEP)
    CALL Histogram(energyGap(:,:),mask,Nbins,binPopulations,binIndices,binMidpoints)
    CALL FepUS(mappingEnergies(:,:,:,1),groundStateEnergy(:,:),G_FEP,binPopulations,binIndices,PMF1D=binGg,minPop=minPop,useBin=printBin)
    CALL ScanFepUs(binMidpoints(:),binGg(:),mask=printBin,stationaryPoints=relative)

    relative(:) = relative(:) - relative(1)
    WRITE(*,*) relative(2), relative(3)
    objectFunction = (relative(3) - dGPS)**2.0d0 ! + (relative(2) - dGTS)**2.0d0

  END FUNCTION objectFunction

END MODULE ObjectFunctions
