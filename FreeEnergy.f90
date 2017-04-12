MODULE FreeEnergy

  ! #DES: Routines to perform free energy calculations using the input and derived data.
  !       These routines are forbidden to USE data from other modules; they strictly receive/return data.

  USE Input, ONLY : beta
  IMPLICIT NONE

!  PRIVATE
!  PUBLIC ::

  CONTAINS

!*

    SUBROUTINE FepUs(mappingEnergies,targetEnergy,G_FEP,binPopulations,binIndices,PMF2Dout,PMF1D,minPop,useBin)

      ! #DES: Perform a FEP/US calculation of the potential of mean force (PMF) using the supplied
      !       histogram, reference free energy values and target/mapping energies.
      ! #TODO: mappingEnergies does not need to be 3D; should contain diagonals of the full mappingEnergies array

      IMPLICIT NONE
      REAL(8), INTENT(IN)            :: mappingEnergies(:,:,:), targetEnergy(:,:)
      REAL(8), INTENT(IN)            :: G_FEP(:)
      INTEGER, INTENT(IN)            :: binPopulations(:,:), binIndices(:,:)
      INTEGER, INTENT(IN)            :: minPop
      REAL(8), INTENT(OUT), OPTIONAL :: PMF2Dout(SIZE(binPopulations,1),SIZE(G_FEP))
      REAL(8), INTENT(OUT)           :: PMF1D(SIZE(binPopulations,1))
      LOGICAL, INTENT(OUT), OPTIONAL :: useBin(SIZE(binPopulations,1))
      INTEGER                        :: Nbins, NfepSteps, NtimeSteps   ! Totals
      INTEGER                        :: bin, fepstep, timestep, popSum ! Indices and Counts
      REAL(8)                        :: PMF2D(SIZE(binPopulations,1),SIZE(G_FEP))

      useBin(:)  = .FALSE.
      Nbins      = SIZE(binPopulations,1)
      NfepSteps  = SIZE(mappingEnergies,2)
      NtimeSteps = SIZE(mappingEnergies,1)

      ! Compute the free energy for the target potential for each populated bin
      ! This generates the Q-style 2D PMF where information from each fepstep is retained separately.

      PMF2D(:,:) = 0.0d0
      DO bin = 1, Nbins               ! for each range of the reaction coordinate
        DO fepstep = 1, NfepSteps     ! check each FEP simulation
          DO timestep = 1, NtimeSteps ! if this timestep of this FEP step is in the bin
            IF (binIndices(fepstep,timestep) == bin) THEN
              PMF2D(bin,fepstep) = PMF2D(bin,fepstep) + EXP( -1.0d0 * beta * (targetEnergy(timestep,fepstep) - mappingEnergies(timestep,fepstep,fepstep)) )
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      WHERE (binPopulations(:,:) > 0) PMF2D(:,:) = PMF2D(:,:) / binPopulations(:,:)

      DO bin = 1, nBins
        DO fepstep = 1, nFepSteps
          IF (binPopulations(bin,fepstep) > 0) THEN
            PMF2D(bin,fepstep) = (-1.0d0/beta) * LOG ( EXP(-1.0d0*beta*G_FEP(fepstep)) * PMF2D(bin,fepstep) )
          ENDIF
        ENDDO 
      ENDDO

      IF (PRESENT(PMF2Dout)) PMF2Dout(:,:) = PMF2D(:,:)

      ! Finally produce the 1D PMF

      PMF1D(:) = 0.0d0
      DO bin = 1, SIZE(binPopulations,1)
        DO fepstep = 1, SIZE(binPopulations,2)
          IF (binPopulations(bin,fepstep) >= minPop) THEN
            PMF1D(bin) = PMF1D(bin) + (PMF2D(bin,fepstep) * binPopulations(bin,fepstep))
          ENDIF
        ENDDO
      ENDDO

      DO bin = 1, SIZE(binPopulations,1)
        popSum = 0
        DO fepstep = 1, SIZE(binPopulations,2)
          IF (binPopulations(bin,fepstep) >= minPop) popSum = popsum + binPopulations(bin,fepstep)
        ENDDO
        IF (popSum > 0) THEN
          PMF1D(bin) = PMF1D(bin) / popSum
          useBin(bin) = .TRUE.
        ENDIF
      ENDDO

    ENDSUBROUTINE FepUs


      ! The actual result is a 1D PMF where one value is given per bin. Q and Mapping seem to disagree on how to do this.
      ! For Q, a bin average is found by sum_i[PMF_i * (N_i/N)] where N is the total points in the bin over all fep steps.
      ! For Molaris, the approach is to use the same (most populated) mapping potential to compute comtributions to the PMF
      ! from all of the points in a bin, regardless the FEP simulation the came from. This requires maintaining the binFepIndex array
      ! and using it to compute the PMF directly. Below is the Molaris algorithm.

      ! This is only needed if the MOLARIS algorithm is used
!*

    SUBROUTINE FepUsMolaris(mappingEnergies,targetEnergy,G_FEP,binPopulations,binIndices,PMF)

      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: mappingEnergies(:,:,:), targetEnergy(:,:)
      REAL(8), INTENT(IN)  :: G_FEP(:)
      INTEGER, INTENT(IN)  :: binPopulations(:,:)
      INTEGER, INTENT(IN)  :: binIndices(:,:)

      REAL(8), INTENT(OUT) :: PMF(SIZE(binPopulations,1))

      INTEGER :: binFepIndex(SIZE(binPopulations,1))
      INTEGER :: Nbins, NfepSteps, NtimeSteps, bin, fepstep, timestep

      Nbins      = SIZE(binPopulations,1)
      NfepSteps  = SIZE(mappingEnergies,2)
      NtimeSteps = SIZE(mappingEnergies,1)

      binFepIndex(:) = MAXLOC(binPopulations(:,:),2) !for a given bin, which FEP simulation produces the most members

      PMF(:) = 0.0d0
      DO bin = 1, Nbins               ! for each range of the reaction coordinate
        DO fepstep = 1, NfepSteps     ! check each FEP simulation
          DO timestep = 1, NtimeSteps ! if this timestep of this FEP step is in the bin
            IF (binIndices(fepstep,timestep) == bin) THEN
              PMF(bin) = PMF(bin) + EXP( -1.0d0 * beta * (targetEnergy(timestep,fepstep) - mappingEnergies(timestep,binFepIndex(bin),binFepIndex(bin))) )
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      PMF(:) = PMF(:) / SUM(binPopulations(bin,:)) !need to mask this for binPopulations(i,j) = 0 to prevent NaN

      DO bin = 1, nBins
        DO fepstep = 1, nFepSteps
          PMF(bin) = (-1.0d0/beta) * LOG ( EXP(-1.0d0*beta*G_FEP(binFepIndex(bin))) * PMF(bin) )
        ENDDO
      ENDDO

    END SUBROUTINE FepUsMolaris

!*

    SUBROUTINE Histogram2D(data,mask,N,linearIndices,linearPopulations,linearCoordinates)

      IMPLICIT NONE

      INTEGER, PARAMETER  :: D = 2
      REAL(8), INTENT(IN) :: data(:,:,:) ! 2 x nFepSteps x nTimesteps
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: mask(:,:)   ! nFepSteps x nTimesteps

      INTEGER, INTENT(OUT) :: linearPopulations(N**SIZE(data,1),SIZE(data,2)) !population of a bin
      INTEGER, INTENT(OUT) :: linearIndices(SIZE(data,2),SIZE(data,3))        !indexes into the linear arrays
      REAL(8), INTENT(OUT) :: linearCoordinates(N**SIZE(data,1),SIZE(data,1)) !coordinates of a bin

      INTEGER :: index(SIZE(data,1))
      REAL(8) :: binWidth(D), min(SIZE(data,1)), max(SIZE(data,1))
      INTEGER :: i, j, k, dim, bin, fepstep, timestep

      IF (SIZE(data,1) /= D) STOP "Data of incorrect dimensionality passed to Histogram2D"

      DO dim = 1, SIZE(data,1)
        min(dim) = MINVAL(data(dim,:,:),MASK=mask)
        max(dim) = MAXVAL(data(dim,:,:),MASK=mask)
        binWidth(dim) = (max(dim) - min(dim)) / N
      ENDDO

      DO dim = 1, SIZE(data,1)
        DO i = 1, N
          DO j = 1, N
            ! This is the last remaining chunk to generalize. The problem is the loop over i and j to generate the bin indices
            ! which does not easily generalize.
            ! what if, instead, you could iterate over k and, knowing N and D, generate the values of the D indices? This is the 'inverse'
            ! of the function that produces k.
            !
            ! DO k = 1, N**D
            !   DO d = 1, D
            !     index(d) = f(k,d,D,N)
            !   ENDDO
            !   DO d = 1, D
            !     linearCoordinates(k,d) = min(d) + (DBLE(index(d))-0.5d0) * binWidth(dim)
            !   ENDDO
            ! ENDDO
            !
            ! This would, I think, work for any D provided N is constant. The problem is what is the f(k,d,D,N) that generates x(D)?
            k = N*(i-1)+j
            linearCoordinates(k,1) = min(1) + (DBLE(i)-0.5d0) * binWidth(dim)
            linearCoordinates(k,2) = min(2) + (DBLE(j)-0.5d0) * binWidth(dim)
          ENDDO
        ENDDO
      ENDDO

      linearPopulations(:,:) = 0
      linearIndices(:,:) = 0
      DO fepstep = 1, SIZE(data,2)
        DO timestep = 1, SIZE(data,3)
          
          IF (mask(fepstep,timestep) .EQV. .FALSE.) CYCLE

          index(:) = 0
          DO dim = 1, D
            bin = FLOOR((data(dim,fepstep,timestep) - min(dim)) / binWidth(dim)) + 1
            IF (bin <= N) THEN
              index(dim) = bin
            ELSE
              index(dim) = N
            ENDIF
          ENDDO

          k = N * (index(1) - 1) + index(2) ! k is the linearized index, should be a function of all D indices, D and N to be general.
          linearIndices(fepstep,timestep) = k
          linearPopulations(k,fepstep) = linearPopulations(k,fepstep) + 1

        ENDDO
      ENDDO

    END SUBROUTINE Histogram2D

!*

    SUBROUTINE Histogram(data,useData,N,binPopulations,binIndices,binMidpoints)

      ! #DES: Generate a histogram along the reaction coordinate

      USE ArrayUtil, ONLY : linspace
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: N                   ! Number of bins
      REAL(8), INTENT(IN)  :: data(:,:)           ! values to histogram: nFep x nTimesteps
      LOGICAL, INTENT(IN)  :: useData(:,:)        ! mask for data values
      INTEGER, INTENT(OUT) :: binPopulations(:,:)
      INTEGER, INTENT(OUT) :: binIndices(:,:)
      REAL(8), INTENT(OUT) :: binMidpoints(:)

      LOGICAL, PARAMETER :: DEBUG = .FALSE.
      REAL(8) :: min, max, binWidth, binEdges(N+1)
      INTEGER :: bin, fepstep, timestep

      min = MINVAL(data,MASK=useData .EQV. .TRUE.)
      max = MAXVAL(data,MASK=useData .EQV. .TRUE.)
      binEdges = linspace(min,max,N+1)
      binWidth = (max - min) / N

      DO bin = 1, N
        binMidpoints(bin) = binEdges(bin) + (binWidth / 2.0d0)
      ENDDO

      binIndices(:,:) = 0
      binPopulations(:,:) = 0
      DO fepstep = 1, SIZE(data,1)
        DO timestep = 1, SIZE(data,2)
          IF (useData(fepstep,timestep)) THEN
            bin = FLOOR((data(fepstep,timestep) - min) / binWidth) + 1
            binIndices(fepstep,timestep) = bin
            IF (bin > N) THEN !max value located
              binPopulations(N,fepstep) = binPopulations(N,fepstep) + 1
              binIndices(fepstep,timestep) = N
            ELSE
              binPopulations(bin,fepstep) = binPopulations(bin,fepstep) + 1
              binIndices(fepstep,timestep) = bin
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF (DEBUG) THEN

        WRITE(*,'(A)') "1D Histogram Details"; WRITE(*,*) ""
        WRITE(*,'(A,F10.2)') "Min. Value = ", min
        WRITE(*,'(A,F10.2)') "Max. Value = ", max
        WRITE(*,'(A,I10)') "Num. Bins  = ", N 
        WRITE(*,'(A,F10.2)') "Bin width  = ", binWidth
        WRITE(*,*) ""
        WRITE(*,*) "Bin #  x_lower <= x < x_upper"
        DO bin = 1, N
          WRITE(*,'(I3,5X,3F10.2)') bin, binEdges(bin), binEdges(bin) + (binWidth/2.0d0), binEdges(bin+1)
        ENDDO

        WRITE(*,*) ""
        WRITE(*,*) "Bin Population Data"; WRITE(*,*) ""
        WRITE(*,'(A)',ADVANCE='NO') "Fep #"
        DO fepstep = 1, SIZE(data,1)
          WRITE(*,'(I4,1X)',ADVANCE='NO') fepstep
        ENDDO
        WRITE(*,*)
        WRITE(*,'(A)') "Bin #"
        DO bin = 1, N
          WRITE(*,'(I5)',ADVANCE='NO') bin
          DO fepstep = 1, SIZE(data,1)
            WRITE(*,'(I8,1X)',ADVANCE='NO') binPopulations(bin,fepstep)
          ENDDO
          WRITE(*,*)
        ENDDO

        WRITE(*,*) ""

      ENDIF

      RETURN

    END SUBROUTINE Histogram

!*

    ! Compute the free energy change from FEP step n to m
    ! That is, dG_n->m where usually m = (n+1) 
    REAL(8) FUNCTION ComputeStepFreeEnergy(E_n,E_m,mask)

      USE StatisticalFunctions, ONLY : mean
      USE ArrayUtil, ONLY : ArrayEXP
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: E_n(:), E_m(:)
      LOGICAL, INTENT(IN) :: mask(:)
      REAL(8) :: average

      IF (SIZE(E_n) /= SIZE(E_m)) STOP "Error: ComputeStepFreeEnergy"

      average = mean(ArrayEXP(-1.0d0 * beta * (E_m(:) - E_n(:))),mask(:))
      ComputeStepFreeEnergy = (-1.0d0 / beta) * LOG(average)

    END FUNCTION ComputeStepFreeEnergy

!*

    !Compute the FEP free energy increments between states
    !There are several things you might want to get from this - forward, backward, profile
    SUBROUTINE ComputeFEPIncrements(A,B,mappingEnergies,mask,forward,reverse,profile)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: A, B
      REAL(8), INTENT(IN) :: mappingEnergies(:,:,:) !timesteps, fepsteps, fepsteps
      LOGICAL, INTENT(IN) :: mask(:,:)

      REAL(8), INTENT(OUT) :: profile(:)
      REAL(8), OPTIONAL, INTENT(OUT) :: forward(:), reverse(:)

      REAL(8) :: dG_f(B-A), dG_r(B-A)
      INTEGER :: n, m

      dG_f(:) = 0.0d0
      dG_r(:) = 0.0d0
      profile(:)   = 0.0d0

      DO n = A, B - 1
        m = n + 1
        dG_f(n) = computeStepFreeEnergy(mappingEnergies(:,n,n),mappingEnergies(:,m,n),mask(n,:))
        dG_r(n) = computeStepFreeEnergy(mappingEnergies(:,m,m),mappingEnergies(:,n,m),mask(m,:))
        profile(n)   = SIGN(1.0d0,dG_f(n)) * (ABS(dG_f(n)) + ABS(dG_r(n))) / 2.0d0
      ENDDO

      IF (PRESENT(forward)) forward(:) = dG_f(:)
      IF (PRESENT(reverse)) reverse(:) = dG_r(:)

    END SUBROUTINE ComputeFEPIncrements

!*

    SUBROUTINE ComputeFEPProfile(A,B,mappingEnergies,mask,profile)

      ! #DES: This routine sums the averaged for/rev free energy increments to generate the FEP profile

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A, B
      REAL(8), INTENT(IN) :: mappingEnergies(:,:,:) !timesteps, fepsteps, fepsteps
      LOGICAL, INTENT(IN) :: mask(:,:)
      REAL(8), INTENT(OUT) :: profile(SIZE(mappingEnergies,2))
      REAL(8) :: increments(SIZE(mappingEnergies,2)-1)
      INTEGER :: step

      CALL ComputeFEPIncrements(A,B,mappingEnergies,mask,profile=increments)
      profile(:) = 0.0d0
      DO step = 2, SIZE(profile)
        profile(step) = SUM(increments(1:step-1))
      ENDDO      

    END SUBROUTINE ComputeFEPProfile

!*

    !Get the linear response approximation for moving from potential A to potential B
    !The LRA need only be computed when an energetic breakdown is required, or to justify its use by comparison of
    !its predictions of total energy changes to those of the FEP method
    REAL(8) FUNCTION LRA(energyType,A,B,mask)

      USE StatisticalFunctions, ONLY : mean
      USE Data, ONLY : mappingEnergies
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: energyType, A, B
      LOGICAL, INTENT(IN) :: mask(:,:)

      LRA = 0.5d0 * mean(mappingEnergies(:,A,B,energyType),mask(A,:)) - mean(mappingEnergies(:,A,A,energyType),mask(A,:)) + &
&                   mean(mappingEnergies(:,B,B,energyType),mask(B,:)) - mean(mappingEnergies(:,B,A,energyType),mask(B,:))

    END FUNCTION LRA

END MODULE FreeEnergy
