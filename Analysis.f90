MODULE Analysis

  ! #DES: Routines for applying analysis to the input/derived data and writing the results
  !       to data files for visualisation. These routines should USE no data, only procedures.
  !       They should instead accept all data and options as arguments. For the sake of consistency
  !       each should also accept an already open output unit to write to.

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: AnalyzeSimulationConvergence, &
          & WriteMeanEnergyBreakdown,     &
          & DetailedFEP,                  &
          & ReportMeans,                  &
          & Fepus2D,                      &
          & FepBreakdown,                 &
          & RunLinearResponse,            & 
          & AnalyzeFep,                   &
          & FepUsGroundState,             & 
          & FepUsFreeEnergies

  CONTAINS

!*

    SUBROUTINE FepUsGroundState(energyGap,groundStateEnergy,mappingEnergies,mask,Nbins,minPop,outUnit)

      ! #DES: Compute the ground state PMF via the FEP/US method and write statistically significant 
      !       values to the output unit.

      USE FreeEnergy, ONLY : Histogram, ComputeFepProfile, FepUS
      USE Output,     ONLY : WriteCSV2D

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: energyGap(:,:), groundStateEnergy(:,:), mappingEnergies(:,:,:)
      LOGICAL, INTENT(IN) :: mask(:,:)
      INTEGER, INTENT(IN) :: Nbins, minPop, outUnit

      INTEGER      :: binPopulations(Nbins,SIZE(energyGap,1)), binIndices(SIZE(energyGap,1),SIZE(energyGap,2))
      REAL(8)      :: binMidpoints(Nbins)
      REAL(8)      :: binGg(Nbins)  
      REAL(8)      :: G_FEP(SIZE(energyGap,1))
      LOGICAL      :: printBin(Nbins)
      INTEGER      :: bin, count
      CHARACTER(3) :: head(2)
      REAL(8)      :: output(Nbins,2)

      head(1) = "dE"; head(2) = "dG"  

      CALL ComputeFEPProfile(1,SIZE(energyGap,1),mappingEnergies(:,:,:),mask(:,:),profile=G_FEP(:))

      CALL Histogram(energyGap(:,:),mask(:,:),Nbins,binPopulations(:,:),binIndices(:,:),binMidpoints(:))
      CALL FepUS(mappingEnergies(:,:,:),groundStateEnergy(:,:),G_FEP,binPopulations,binIndices,PMF1D=binGg,minPop=minPop,useBin=printBin)

      count = 0
      DO bin = 1, nBins
        IF (printBin(bin) .EQV. .TRUE.) THEN
          count = count + 1
          output(count,1) = binMidpoints(bin)
          output(count,2) = binGg(bin)
       ENDIF
      ENDDO

      CALL WriteCsv2D(head,output(1:count,:),outUnit)

    END SUBROUTINE FepUSGroundState

!*

    SUBROUTINE FepUsFreeEnergies(energyGap,groundStateEnergy,mappingEnergies,mask,Nbins,minPop,relative,outUnit)

      ! #DES: Routine to compute and return the RS, TS and PS energies on a given PMF.
      !       Assumes a general shape for the PES (through ScanFepUs) which may not be true.

      USE FreeEnergy, ONLY : ScanFepUs, ComputeFEPProfile, Histogram, FepUS

      IMPLICIT NONE

      LOGICAL, INTENT(IN), OPTIONAL :: relative
      REAL(8), INTENT(IN) :: energyGap(:,:), groundStateEnergy(:,:), mappingEnergies(:,:,:)
      LOGICAL, INTENT(IN) :: mask(:,:)
      INTEGER, INTENT(IN) :: Nbins, minPop, outUnit

      INTEGER :: binPopulations(Nbins,SIZE(energyGap,1)), binIndices(SIZE(energyGap,1),SIZE(energyGap,2))
      REAL(8) :: binMidpoints(Nbins)
      REAL(8) :: binGg(Nbins)
      REAL(8) :: G_FEP(SIZE(energyGap,1))
      LOGICAL :: printBin(Nbins)
      REAL(8) :: dG(3)

      CALL ComputeFEPProfile(1,SIZE(energyGap,1),mappingEnergies(:,:,:),mask(:,:),profile=G_FEP)
      CALL Histogram(energyGap(:,:),mask,Nbins,binPopulations,binIndices,binMidpoints)
      CALL FepUS(mappingEnergies(:,:,:),groundStateEnergy(:,:),G_FEP,binPopulations,binIndices,PMF1D=binGg,minPop=minPop,useBin=printBin)            
      CALL ScanFepUs(binMidpoints(:),binGg(:),mask=printBin,stationaryPoints=dG(:))

      IF (PRESENT(relative) .AND. relative .EQV. .TRUE.) dG(:) = dG(:) - dG(1)

      WRITE(outUnit,'(3F15.8)') dG(:)

    END SUBROUTINE FepUsFreeEnergies

!*

    SUBROUTINE FepBreakdown(lambda,mappingEnergies,mask,energyNames,outUnit)

      ! #DES: Produce and write the type-wise breakdown of the FEP result

      USE Output, ONLY     : WriteCsv2D
      USE FreeEnergy, ONLY : ComputeFEPProfile

      IMPLICIT NONE

      REAL(8), INTENT(IN)      :: lambda(:), mappingEnergies(:,:,:,:)
      LOGICAL, INTENT(IN)      :: mask(:,:)
      INTEGER, INTENT(IN)      :: outUnit
      CHARACTER(*), INTENT(IN) :: energyNames(:)

      CHARACTER(15) :: head(1+SIZE(energyNames)) ! 1 per column
      REAL(8) :: output(SIZE(mappingEnergies,2),1+SIZE(energyNames))
      INTEGER :: step, type
      REAL(8) :: profile(SIZE(mappingEnergies,2))

      head(1) = "lambda"
      
      DO type = 1, SIZE(energyNames)

        CALL ComputeFEPProfile(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:,type),mask(:,:),profile=profile)      

        head(type+1) = energyNames(type)
        DO step = 1, SIZE(mappingEnergies,2)
          output(step,1) = lambda(step)
          output(step,type+1) = profile(step)
        ENDDO

      ENDDO

      CALL WriteCSV2D(head,output,outUnit)

    END SUBROUTINE FepBreakdown

!*

    SUBROUTINE DetailedFEP(mappingEnergies,mask)

      USE Output, ONLY : WriteCSV2D
      USE Data, ONLY : lambda
      USE FreeEnergy, ONLY : ComputeFEPIncrements

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: mappingEnergies(:,:,:,:)
      LOGICAL, INTENT(IN) :: mask(:,:)

      REAL(8)       :: output(SIZE(mappingEnergies,2)-1,5)
      CHARACTER(15) :: head(5)
      INTEGER       :: step
      REAL(8)       :: forward(SIZE(mappingEnergies,2)-1), reverse(SIZE(mappingEnergies,2)-1), profile(SIZE(mappingEnergies,2)-1)

      head(1) = "lambda"
      head(2) = "forward"
      head(3) = "reverse"
      head(4) = "profile"
      head(5) = "sum"

      ! Each FEP increment being evaluated has convergence properties which can be checked

      CALL ComputeFEPIncrements(1,SIZE(mappingEnergies,2),mappingEnergies(:,:,:,3),mask(:,:),forward,reverse,profile)
      DO step = 1, SIZE(mappingEnergies,2)-1
        output(step,1) = lambda(step)
        output(step,5) = SUM(profile(1:step))
      ENDDO
      output(:,2) = forward
      output(:,3) = reverse
      output(:,4) = profile

      CALL WriteCSV2D(head,output,8)

    END SUBROUTINE DetailedFEP

!*

    ! analyze convergence of simulations in terms of mapping energies of a given type
    SUBROUTINE AnalyzeSimulationConvergence(mappingEnergies,mask,lambda,coeffs,skip,outUnit)

      USE Output, ONLY : WriteCsv2D
      USE StatisticalFunctions, ONLY : mean, varianceOfMean, FlyvbjergPetersen

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: lambda(:), coeffs(:,:), mappingEnergies(:,:,:)
      INTEGER, INTENT(IN) :: outUnit, skip(:)
      LOGICAL, INTENT(IN) :: mask(:,:)

      CHARACTER(15) :: head(8)
      REAL(8)       :: output(SIZE(mappingEnergies,2),8)
      INTEGER       :: step
      CHARACTER(4)  :: stepString

      ! for each FEP simulation, determine convergence in terms of the total mapping potential
      ! n lambda_n Ca Cb N sparkline <En> sigma(<E_n>) sparkline tau sigma(<E_n>)

      head(1) = "n"
      head(2) = "lambda_n"
      head(3) = "Ca"
      head(4) = "Cb"
      head(5) = "Nskip"
      head(6) = "<En>"
      head(7) = "sigma(<E_n>)"
      head(8) = "sigma(<E_n>)"

      DO step = 1, SIZE(mappingEnergies,2)

        output(step,1) = DBLE(step)
        output(step,2) = lambda(step)
        output(step,3) = coeffs(step,1)
        output(step,4) = coeffs(step,2)
        output(step,5) = skip(step)
        output(step,6) = mean(mappingEnergies(:,step,step),mask(step,:))

        !generate sparkline data
        WRITE(stepString,'(I0.4)') step
        CALL WriteRunningMean(mappingEnergies(:,step,step),mask(step,:),"FEP_"//stepString)
        output(step,7) = varianceOfMean(mappingEnergies(:,step,step),mask(step,:))
        output(step,8) = FlyvbjergPetersen(mappingEnergies(skip(step):,step,step))

      ENDDO

      CALL WriteCSV2D(head,output,outUnit)

    END SUBROUTINE AnalyzeSimulationConvergence

!*

    SUBROUTINE WriteMeanEnergyBreakdown(stateEnergy,mask,energyNames,lambda,outUnit)

      ! #DES: Write the analog of the QFEP table of means in plottable form.
      ! This allows you to see the contributions of various energy types along the order parameter.
      ! Assumes that the means computed from the masked data are accurate.

      USE StatisticalFunctions, ONLY : mean
      USE Output, ONLY : WriteCsv2D

      IMPLICIT NONE

      REAL(8), INTENT(IN)      :: stateEnergy(:,:,:,:), lambda(:)
      LOGICAL, INTENT(IN)      :: mask(:,:)
      CHARACTER(*), INTENT(IN) :: energyNames(:)
      INTEGER, INTENT(IN)      :: outUnit

      REAL(8)       :: output(SIZE(stateEnergy,2),1+SIZE(stateEnergy,3)*SIZE(stateEnergy,4))
      INTEGER       :: step, state, type, index
      CHARACTER(15) :: head(1+SIZE(stateEnergy,3)*SIZE(stateEnergy,4))
      CHARACTER(4)  :: stateString

      head(1) = 'lambda_n'
      DO type = 1, SIZE(stateEnergy,4)    ! nEnergyTypes
        DO state = 1, SIZE(stateEnergy,3) ! nStates
          index = 1 + (2*type-2) + state
          WRITE(stateString,'(I0.4)') state
          head(index) = TRIM(ADJUSTL(energyNames(type)))//stateString
        ENDDO
      ENDDO

      DO step = 1, SIZE(stateEnergy,2)   ! nFepSteps
        DO type = 1, SIZE(stateEnergy,4) ! nEnergyTypes
          output(step,1) = lambda(step)
          DO state = 1, SIZE(stateEnergy,3) !nStates
            index = 1 + (2*type-2) + state
            output(step,index) = mean(stateEnergy(:,step,state,type),mask(step,:))
          ENDDO
        ENDDO
      ENDDO

      CALL WriteCSV2D(head,output,outUnit)

    END SUBROUTINE WriteMeanEnergyBreakdown

!*

    SUBROUTINE AnalyzeConvergenceBasic(property,mask,output)

      ! #DES: This routine does not conform to the rules of this module - should be a subsubroutine where it is used.

      USE StatisticalFunctions, ONLY : mean, varianceOfMean

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: property(:)
      LOGICAL, INTENT(IN) :: mask(:)
      REAL(8), INTENT(OUT) :: output(4)

      output(1) = SIZE(property) - COUNT(mask)
      output(2) = SIZE(property)
      output(3) = mean(property(:),mask(:))
      output(4) = varianceOfMean(property(:),mask(:))

    END SUBROUTINE AnalyzeConvergenceBasic

!*
    ! WriteRunningMean - Write the element-by-element mean of the array 'property' to a file named 'name'_mean.csv
    ! Monopolizes unit 101, needs to be cleaned so that it takes an already open file!

    SUBROUTINE writeRunningMean(property,mask,name)

      USE StatisticalFunctions, ONLY : runningMean
      USE FileIO, ONLY : OpenFile, CloseFile
      USE Input, ONLY : outDir
      USE Output, ONLY : WriteCSV2D

      IMPLICIT NONE

      REAL(8), INTENT(IN)      :: property(:)
      LOGICAL, INTENT(IN)      :: mask(:)
      CHARACTER(*), INTENT(IN) :: name

      INTEGER, PARAMETER :: oUnit = 101
      CHARACTER(4), PARAMETER :: head(3) = (/"i   ","x_i ","mean"/)
      INTEGER :: i
      REAL(8) :: out(SIZE(property),3)

      DO i = 1, SIZE(property)
        out(i,1) = DBLE(i)
      ENDDO

      out(:,2) = property(:)
      out(:,3) = runningMean(property,mask)
      
      CALL OpenFile(oUnit,TRIM(ADJUSTL(outDir))//"/"//trim(adjustl(name))//"_runningmean.csv","write") !should already have unit of an open file
      CALL WriteCSV2D(head,out,oUnit)
      CALL CloseFile(oUnit)

    END SUBROUTINE writeRunningMean

!*

    ! This routine is not conforming to the rules of the module.
    SUBROUTINE ReportMeans()

      ! ReportMeans - csv output of FEP data in qfep style
      ! Each row describes a single FEP step in a given state and provides:
      ! step index, step lambda, number of timesteps, state, step state coefficient, mean of each defined energy type

      USE Input, ONLY : energyNames, stateEnergy, coeffs
      USE Data, ONLY : lambda
      USE StatisticalFunctions, ONLY : mean
      USE Output, ONLY : WriteCSV2D

      IMPLICIT NONE

      INTEGER                    :: nFepSteps, nStates, nColumns, nRows
      INTEGER                    :: name, fepstep, state, type, row
      CHARACTER(20), ALLOCATABLE :: head(:)
      REAL(8),       ALLOCATABLE :: out(:,:)

      nFepSteps = SIZE(stateEnergy,2)
      nStates = SIZE(stateEnergy,3)
      nColumns = SIZE(energyNames) + 5
      nRows = nFepSteps * nStates

      ALLOCATE(head(nColumns))
      ALLOCATE(out(nRows,nColumns))

      head(1) = "n";  head(2) = "lambda_n"
      head(3) = "N";  head(4) = "x"
      head(5) = "c_x"
      DO name = 1, SIZE(energyNames)
        head(5+name) = energyNames(name)
      ENDDO

      DO fepstep = 1, nFepSteps
        DO state = 1, nStates
          row = nStates*(fepstep-1) + state
          out(row,1) = fepstep
          out(row,2) = lambda(fepstep)
          out(row,3) = SIZE(stateEnergy(state,type,fepstep,:),1)
          out(row,4) = state
          out(row,5) = coeffs(1,fepstep,state)
          DO type = 1, SIZE(energyNames)
            out(row,5+type) = mean(stateEnergy(:,fepstep,state,type))
          ENDDO
        ENDDO
      ENDDO

      CALL WriteCSV2D(head,out,6) !to stdout for now

      IF (ALLOCATED(head)) DEALLOCATE(head)
      IF (ALLOCATED(out))   DEALLOCATE(out)

    END SUBROUTINE ReportMeans

!*

    SUBROUTINE Fepus2D(geomRC,groundStateEnergy,mappingEnergies,mask,N,minPop,unit)

      ! #DES: Compute and print the 2D free energy surface

      USE FreeEnergy, ONLY : Histogram2D, ComputeFEPProfile, FepUs
      USE FileIO,     ONLY : OpenFile, CloseFile
      USE Output,     ONLY : WriteCsv2D

      IMPLICIT NONE

      REAL(8), INTENT(IN) :: geomRC(:,:,:), groundStateEnergy(:,:), mappingEnergies(:,:,:)
      LOGICAL, INTENT(IN) :: mask(:,:)
      INTEGER, INTENT(IN) :: N, minPop, unit

      INTEGER :: binIndices(SIZE(geomRC,2),SIZE(geomRC,3))
      INTEGER :: binPopulations(N**SIZE(geomRC,1),SIZE(geomRC,2))
      REAL(8) :: binCoordinates(N**SIZE(geomRC,1),SIZE(geomRC,1))

      REAL(8) :: G_FEP(SIZE(geomRC,2))
      REAL(8) :: binG(N**SIZE(geomRC,1))
      LOGICAL :: printBin(N**SIZE(geomRC,1))

      REAL(8)      :: output(N**SIZE(geomRC,1),SIZE(geomRC,1)+1) ! one row for each bin to be printed, one column for each coordinate and one for the free energy
      CHARACTER(4) :: head(SIZE(geomRC,1)+1)                     ! heading for each coord, then the free energy
      INTEGER      :: dim, bin, count
      CHARACTER(1) :: dimString

      DO dim = 1, SIZE(geomRC,1)
        WRITE(dimString,'(I1)') dim
        head(dim) = "R"//dimString
      ENDDO
      head(SIZE(geomRC,1)+1) = "dG"

      G_FEP = 0.0d0
      CALL ComputeFEPProfile(1,SIZE(geomRC,2),mappingEnergies(:,:,:),mask(:,:),profile=G_FEP)
      CALL Histogram2D(geomRC(:,:,:),mask(:,:),N,binIndices(:,:),binPopulations(:,:),binCoordinates(:,:))
      printBin(:) = .FALSE.

      CALL FepUs(mappingEnergies(:,:,:),groundStateEnergy(:,:),G_FEP(:),binPopulations(:,:),binIndices(:,:),PMF1D=binG,minPop=minPop,useBin=printBin(:))

      count = 0
      DO bin = 1, N**SIZE(geomRC,1)
        IF (printBin(bin) .EQV. .TRUE.) THEN
          count = count + 1
          DO dim = 1, SIZE(geomRC,1)
            output(count,dim) = binCoordinates(bin,dim)
          ENDDO
          output(count,SIZE(geomRC,1)+1) = binG(bin)
        ENDIF
      ENDDO

      CALL WriteCsv2D(head,output(1:count,:),unit)

    END SUBROUTINE Fepus2D

!*

    !This is the initial convergence check in terms of the energy.
    !Qfep prints <E> and the means of its components, <E_i> for each FEP simulation
    !which allows the user to confirm that each component is correct and determine where 
    !issues in the total energy may come from.
    !Fepcat can do better by giving the variance of (assumed) uncorrelated data, evaluating
    !the extent of the correlation and correcting the variance appropriately.

    !Since the total energy is a sum over components, the components can be reported on in detail
    !individually too as the variances will also sum. There should also be a comparative element
    !showing the contribution of each term to the convergence of the simulation.

    ! AnalyzeConvergence - Report on the convergence of the FEP simulations
    ! The types of data written (to .csv format file) are given in the array 'head'
    ! The variance and mean of the total energies are computed. The variance of the mean is first computed
    ! assuming no correlation, then the correlation is estimated with the Flybjerg-Petersen method.
    ! Running mean and correlation plots are also produced for each step.

    SUBROUTINE AnalyzeFEP(energyType)

      USE Input, ONLY : nFepSteps, mask
      USE Data, ONLY : mappingEnergies, lambda
      USE Output, ONLY : WriteCSV2D

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: energyType
      INTEGER :: step
      CHARACTER(3) :: stepName
      CHARACTER(50):: head(6)
      REAL(8) :: out(nFepSteps,6)

      head(1) = "n"
      head(2) = "lambda_n"
      head(3) = "Mskip"
      head(4) = "M"
      head(5) = "<V>"
      head(6) = "sigma2(<V>)"

      DO step = 1, nFepSteps

        out(step,1) = step
        out(step,2) = lambda(step)
        CALL AnalyzeConvergenceBasic(MappingEnergies(:,step,step,energyType),mask=mask(step,:),output=out(step,3:6))
        WRITE(stepName,'(I0.3)') step - 1
        CALL writeRunningMean(MappingEnergies(:,step,step,energyType),mask(step,:),"mappingStep_"//stepName)

      ENDDO

      CALL WriteCsv2D(head,out,6)

    END SUBROUTINE AnalyzeFEP

!*

    SUBROUTINE RunLinearResponse()

      USE Input, ONLY : nFepSteps, energyNames, mask
      USE FileIO, ONLY : OpenFile, CloseFile
      USE Output, ONLY : WriteCSV2D
      USE FreeEnergy, ONLY : LRA
      IMPLICIT NONE
      INTEGER, PARAMETER :: lraUnit = 56
      REAL(8) :: output(1,SIZE(energyNames))
      INTEGER :: type

      CALL OpenFile(lraUnit,"lra.csv","write")

      DO type = 1, SIZE(energyNames)
        output(1,type) = LRA(type,1,nFepSteps,mask)
      ENDDO

      CALL WriteCSV2D(energyNames,output,lraUnit)
      CALL CloseFile(lraUnit)

    END SUBROUTINE

END MODULE Analysis

! PARTS TO MAKE A FULL FEP/US ANALYSIS
!     CALL FepUS(mappingEnergies(:,:,:,1),stateEnergy(:,:,stateA,1)+alpha(stateA),G_FEP,binPopulations,binIndices,PMF2D=dGa,PMF1D=binGa,minPop=minPop)
!     CALL FepUS(mappingEnergies(:,:,:,1),stateEnergy(:,:,stateB,1)+alpha(stateB),G_FEP,binPopulations,binIndices,PMF2D=dGb,PMF1D=binGb,minPop=minPop)
! REAL(8) :: binGa(Nbins), binGb(Nbins)
! REAL(8) :: dGa(Nbins,SIZE(energyGap,1)), dGb(Nbins,SIZE(energyGap,1))
!; head(3) = "Ga"; head(4) = "Gb"
!       output(bin,3) = binGa(bin)
!       output(bin,4) = binGb(bin)


