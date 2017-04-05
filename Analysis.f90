MODULE Analysis

  ! #DES: Routines for applying analysis to the input/derived data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: AnalyzeSimulationConvergence, WriteMeanEnergyBreakdown, DetailedFEP, ReportMeans, Test2D, FepBreakdown, RunLinearResponse, AnalyzeFep

  CONTAINS

!*

    SUBROUTINE FepBreakdown(lambda,mappingEnergies,mask,energyNames)

      USE Output, ONLY : WriteCsv2D
      USE FileIO, ONLY : OpenFile, CloseFile
      USE FreeEnergy, ONLY : ComputeFEPProfile
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: lambda(:), mappingEnergies(:,:,:,:)
      LOGICAL, INTENT(IN) :: mask(:,:)
      CHARACTER(*), INTENT(IN) :: energyNames(:)

      CHARACTER(15) :: head(1+SIZE(energyNames)) ! 1 per column
      REAL(8) :: output(SIZE(mappingEnergies,2),1+SIZE(energyNames))
      INTEGER :: step, type
      INTEGER, PARAMETER :: outUnit = 14
      REAL(8) :: profile(SIZE(mappingEnergies,2))

      CALL OpenFile(outUnit,"fep-breakdown.csv","write")

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

      CALL CloseFile(outUnit)

    END SUBROUTINE FepBreakdown

!*

    SUBROUTINE DetailedFEP(mappingEnergies,mask)

      USE Output, ONLY : WriteCSV2D
      USE Data, ONLY : lambda
      USE FreeEnergy, ONLY : ComputeFEPIncrements
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: mappingEnergies(:,:,:,:)
      LOGICAL, INTENT(IN) :: mask(:,:)
      REAL(8) :: output(SIZE(mappingEnergies,2)-1,5)
      CHARACTER(15) :: head(5)
      INTEGER :: step
      REAL(8) :: forward(SIZE(mappingEnergies,2)-1), reverse(SIZE(mappingEnergies,2)-1), profile(SIZE(mappingEnergies,2)-1)

      head(1) = "lambda"; head(2) = "forward"; head(3) = "reverse"; head(4) = "profile"; head(5) = "sum"

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
      REAL(8) :: output(SIZE(mappingEnergies,2),8)
      INTEGER :: step
      CHARACTER(4) :: stepString

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
      REAL(8), INTENT(IN) :: stateEnergy(:,:,:,:), lambda(:)
      LOGICAL, INTENT(IN) :: mask(:,:)
      CHARACTER(*), INTENT(IN) :: energyNames(:)
      INTEGER, INTENT(IN) :: outUnit
      REAL(8) :: output(SIZE(stateEnergy,2),1+SIZE(stateEnergy,3)*SIZE(stateEnergy,4))
      INTEGER :: step, state, type, index
      CHARACTER(15) :: head(1+SIZE(stateEnergy,3)*SIZE(stateEnergy,4))
      CHARACTER(4) :: stateString

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
    ! Monopolizes unit 101

    SUBROUTINE writeRunningMean(property,mask,name)

      USE StatisticalFunctions, ONLY : runningMean
      USe FileIO, ONLY : OpenFile, CloseFile
      USE Output, ONLY : WriteCSV2D

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: property(:)
      LOGICAL, INTENT(IN) :: mask(:)
      CHARACTER(*), INTENT(IN) :: name

      REAL(8) :: out(SIZE(property),3)
      INTEGER, PARAMETER :: oUnit = 101
      INTEGER :: i
      CHARACTER(4), PARAMETER :: head(3) = (/"i   ","x_i ","mean"/)

      DO i = 1, SIZE(property)
        out(i,1) = DBLE(i)
      ENDDO

      out(:,2) = property(:)
      out(:,3) = runningMean(property,mask)
      
      CALL OpenFile(oUnit,trim(adjustl(name))//"_runningmean.csv","write")
      CALL WriteCSV2D(head,out,oUnit)
      CALL CloseFile(oUnit)

    END SUBROUTINE writeRunningMean

    ! ReportMeans - csv output of FEP data in qfep style
    ! Each row describes a single FEP step in a given state and provides:
    ! step index, step lambda, number of timesteps, state, step state coefficient, mean of each defined energy type

    SUBROUTINE ReportMeans()

      USE Input, ONLY : energyNames, stateEnergy, coeffs
      USE Data, ONLY : lambda
      USE StatisticalFunctions, ONLY : mean
      USE Output, ONLY : WriteCSV2D

      IMPLICIT NONE
      INTEGER :: nFepSteps, nStates, nColumns, nRows
      INTEGER :: name, fepstep, state, type, row
      CHARACTER(11), ALLOCATABLE :: head(:)
      REAL(8), ALLOCATABLE :: out(:,:)

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

      DEALLOCATE(head)
      DEALLOCATE(out)

    END SUBROUTINE ReportMeans

!*

    !geomRC is the reaction coordinate to be used. It is currently hard-coded for evaluation

    SUBROUTINE Test2D

      ! This was a testing routine for the 2D plot, using the energy gap as both r1 and r2
      ! Should be cannibalized to produce the actual FEP/US driver that calls the correct histogram routine for D.

      USE Data, ONLY : geomRC, mappingEnergies, groundStateEnergy
      USE Input, ONLY : mask, minPop
      USE FreeEnergy, ONLY : Histogram2D, ComputeFEPIncrements, FepUs
      USE Output, ONLY : WriteCsv2D

      IMPLICIT NONE

      INTEGER, PARAMETER :: N = 25

      ! It would be nice to hide these in a Histogram object
      INTEGER :: binIndices(SIZE(geomRC,2),SIZE(geomRC,3))
      INTEGER :: binPopulations(N**SIZE(geomRC,1),SIZE(geomRC,2))
      REAL(8) :: binCoordinates(N**SIZE(geomRC,1),SIZE(geomRC,1))
      ! These should be hidden also so that just the profile can be obtained via a function in the call to FepUs
      REAL(8) :: G_FEP(SIZE(geomRC,2)), dG_FEP(SIZE(geomRC,2)-1)
      ! These are to be printed to the output
      REAL(8) :: binG(N**SIZE(geomRC,1)), dGg(N**SIZE(geomRC,1),SIZE(geomRC,2))

      REAL(8) :: output(N**SIZE(geomRC,1),SIZE(geomRC,1)+1)  ! one row for each bin to be printed, one column for each coordinate and one for the free energy
      CHARACTER(2) :: head(SIZE(geomRC,1)+1) ! heading for each coord, then the free energy
      INTEGER :: bin, fepstep, dim
      CHARACTER(1) :: dimString

      DO dim = 1, SIZE(geomRC,1)
        WRITE(dimString,'(I1)') dim
        head(dim) = "r"//dimString
      ENDDO
      head(SIZE(geomRC,1)+1) = "dg"

      CALL Histogram2D(geomRC,mask,N,binIndices,binPopulations,binCoordinates)

      CALL ComputeFEPIncrements(1,SIZE(geomRC,2),mappingEnergies(:,:,:,1),mask(:,:),profile=dG_FEP)
      G_FEP(:) = 0.0d0
      DO fepstep = 2, SIZE(geomRC,2)
        G_FEP(fepstep) = SUM(dG_FEP(1:fepstep-1))
      ENDDO

      CALL FepUs(mappingEnergies(:,:,:,1),groundStateEnergy,G_FEP,binPopulations,binIndices,PMF2D=dGg,PMF1D=binG,minPop=minPop)

      DO bin = 1, N**SIZE(geomRC,1)
        DO dim = 1, SIZE(geomRC,1)
          output(bin,dim) = binCoordinates(bin,dim)
        ENDDO
        output(bin,SIZE(geomRC,1)+1) = binG(bin)
      ENDDO

      CALL WriteCsv2D(head,output,6)

    END SUBROUTINE Test2D
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

