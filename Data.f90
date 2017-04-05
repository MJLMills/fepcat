MODULE Data

  ! #DES: Module for computing/storing data derived from the input data

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ComputeDerivedData, DeallocateDataArrays, lambda, energyGap, mappingEnergies, GroundStateEnergy, geomRC, computeGroundStateEnergy

  REAL(8), ALLOCATABLE :: lambda(:)
  REAL(8), ALLOCATABLE :: energyGap(:,:), groundStateEnergy(:,:), geomRC(:,:,:)
  REAL(8), ALLOCATABLE :: mappingEnergies(:,:,:,:), offDiagonals(:,:,:,:)

  CONTAINS

!*

    SUBROUTINE AllocateDataArrays(time)

      ! #DES: Allocate and zero memory for all derived data
      ! mappingEnergies, indices are - type, dynamicsPotential, mappingPotential, time

      USE Input, ONLY : nStates, nFepSteps, maxTimesteps, nEnergyTypes, dRC

      IMPLICIT NONE
      REAL(8), INTENT(OUT), OPTIONAL :: time
      REAL(8) :: t1, t2

      IF (PRESENT(time)) CALL CPU_TIME(t1)

      ALLOCATE(lambda(nFepSteps));                                              lambda = 0.0d0
      ALLOCATE(energyGap(nFepSteps,maxTimesteps));                              energyGap = 0.0d0
      ALLOCATE(geomRC(dRC,nFepSteps,maxTimesteps));                             geomRC = 0.0d0
      ALLOCATE(groundStateEnergy(maxTimesteps,nFepSteps));                      groundStateEnergy = 0.0d0
      ALLOCATE(mappingEnergies(maxTimesteps,nFepSteps,nFepSteps,nEnergyTypes)); !mappingEnergies = 0.0d0
      ALLOCATE(offDiagonals(maxTimesteps,nFepSteps,nStates,nStates));           offDiagonals = 0.0d0

      IF (PRESENT(time)) THEN
        CALL CPU_TIME(t2)
        time = t2 - t1
      ENDIF

    END SUBROUTINE AllocateDataArrays

!*

    SUBROUTINE DeallocateDataArrays

      ! #DES: Deallocate all memory for derived data

      IMPLICIT NONE

      IF (ALLOCATED(lambda))            DEALLOCATE(lambda)
      IF (ALLOCATED(energyGap))         DEALLOCATE(energyGap)
      IF (ALLOCATED(geomRC))            DEALLOCATE(geomRC)
      IF (ALLOCATED(groundStateEnergy)) DEALLOCATE(groundStateEnergy)
      IF (ALLOCATED(mappingEnergies))   DEALLOCATE(mappingEnergies)
      IF (ALLOCATED(offDiagonals))      DEALLOCATE(OffDiagonals)

    END SUBROUTINE DeallocateDataArrays

!*

    ! dRC - dimensionality of the reaction coordinate
    ! RCAtoms(dRC,2) - atoms involved in each RC distance

    ! The definition of the reaction coordinate has to be done in the Input module.
    SUBROUTINE ComputeGeometricRC()

      USE Input, ONLY : trajectory
      USE InternalCoords, ONLY : distance
      USE StatisticalFunctions, ONLY : mean

      IMPLICIT NONE
      INTEGER :: x, y, z !delta coord, rC-Cl1 - rC-Cl2
      REAL(8) :: r1, r2
      INTEGER :: step, timestep

      x = 1; y = 2; z = 3

      DO step = 1, SIZE(trajectory,1)
        DO timestep = 1, SIZE(trajectory,4)
          r1 = distance(trajectory(step,:,x,timestep),trajectory(step,:,y,timestep)) !1,2
          r2 = distance(trajectory(step,:,x,timestep),trajectory(step,:,z,timestep)) !1,3
          geomRC(1,step,timestep) = r1
          geomRC(2,step,timestep) = r2
        ENDDO
      ENDDO

    ENDSUBROUTINE ComputeGeometricRC

!* OFF until needed
!    SUBROUTINE RecomputeDependentData(alpha)
!
      ! #DES: Recompute all quantities that depend on the EVB parameters
!
!      IMPLICIT NONE
!      REAL(8), INTENT(IN) :: alpha(:)

!      CALL ComputeMappingEnergies(alpha)
!      CALL ComputeEnergyGap(alpha)
!      CALL ComputeGroundStateEnergy(alpha)

!    END SUBROUTINE RecomputeDependentData

!*

    SUBROUTINE ComputeDerivedData(logUnit,doTiming)

      ! #DES: Public master subroutine for computing all derived data from the inputs
      ! As this is potentially expensive, can time each piece

      USE Input, ONLY : alpha, readTrajectory, couplingA, couplingExponent, nStates
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: logUnit
      LOGICAL, INTENT(IN) :: doTiming
      REAL(8) :: time(3), t1, t2, diff

      WRITE(logUnit,'(A)') "Computing Data Derived from Input"; WRITE(logUnit,*)

      IF (doTiming) CALL CPU_TIME(t1)

      IF (doTiming) THEN 
        CALL AllocateDataArrays(time(1))
        CALL ComputeMappingEnergies(alpha,time(2))
      ELSE 
        CALL AllocateDataArrays()
        CALL ComputeMappingEnergies(alpha)
      ENDIF

      CALL ComputeLambdas()
      CALL ComputeEnergyGap(alpha)

      IF (doTiming) THEN
        CALL ComputeGroundStateEnergy(alpha,time(3))
      ELSE
        CALL ComputeGroundStateEnergy(alpha)
      ENDIF

      CALL ComputeOffDiagonals(energyGap,couplingA,couplingExponent,nStates)

      ! This should only be called when trajectory information is available
      IF (readTrajectory) CALL ComputeGeometricRC()
      
      IF (doTiming) THEN
        CALL CPU_TIME(t2); diff = t2 - t1
        WRITE(logUnit,'(A)')              "Type                  Time(s) %Time"
        WRITE(logUnit,'(A,F7.2)')         "Total               ", diff
        diff = diff * 100
        WRITE(logUnit,'(A,F7.2,3X,F5.1)') "Array Allocation    ", time(1), time(1)/diff
        WRITE(logUnit,'(A,F7.2,3X,F5.1)') "Mapping Energy      ", time(2), time(2)/diff
        WRITE(logUnit,'(A,F7.2,3X,F5.1)') "Ground State Energy ", time(3), time(3)/diff
        WRITE(logUnit,*)
      ENDIF

      WRITE(logUnit,'(A)') "Finished Computing Data Derived from Input"; WRITE(logUnit,*)


    END SUBROUTINE ComputeDerivedData

!*

    SUBROUTINE ComputeOffDiagonals(RC,couplingA,couplingExponent,nStates)

      ! #DES: Compute the off-diagonals for the EVB Hamiltonian so the GS can be evaluated

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: RC(:,:), couplingA(:,:), couplingExponent(:,:)
      INTEGER, INTENT(IN) :: nStates
      INTEGER :: step, timestep, i, j

      OffDiagonals = 0.0d0
      DO step = 1, SIZE(energyGap,1)
        DO timestep = 1, SIZE(energyGap,2)
          DO i = 1, nStates
            DO j = i+1, nStates
              OffDiagonals(timestep,step,i,j) = couplingA(i,j) * EXP(-1.0*couplingExponent(i,j)*RC(step,timestep)*RC(step,timestep))
              OffDiagonals(timestep,step,j,i) = OffDiagonals(timestep,step,i,j) ! EVB Hamiltonian must be symmetric
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE ComputeOffDiagonals

!*

    SUBROUTINE ComputeGroundStateEnergy(alpha,time)

      ! #DES: Calculate the ground state energy for each conformation

      USE Input, ONLY : stateEnergy, nStates
      USE Matrix, ONLY : Eigenvalues2DRealSymmetric, Eigenvalues3DRealSymmetric

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: alpha(:) !, sigma(:,:)
      REAL(8), INTENT(OUT), OPTIONAL :: time
      REAL(8) :: t1, t2
      INTEGER, PARAMETER :: total = 1 !gives the index of the total energy
      INTEGER :: step, timestep, i, j
      REAL(8) :: H(nStates,nStates), eigenvalues(nStates)

      IF (PRESENT(time)) CALL CPU_TIME(t1)

      DO step = 1, SIZE(stateEnergy,2)
        DO timestep = 1, SIZE(stateEnergy,1)

          ! Form the EVB Hamiltonian at this timestep
          DO i = 1, nStates
            DO j = i, nStates
              IF (i == j) THEN
                H(i,i) = stateEnergy(timestep,step,i,total) + alpha(i)
              ELSE
                H(i,j) = OffDiagonals(timestep,step,i,j)
                H(j,i) = OffDiagonals(timestep,step,j,i)
              ENDIF
            END DO
          ENDDO

          ! Obtain the eigenvalues of the Hamiltonian
          IF (nStates == 2) THEN
            eigenvalues(:) = Eigenvalues2DRealSymmetric(H)
          ELSE IF (nStates == 3) THEN
            eigenvalues(:) = Eigenvalues3DRealSymmetric(H)
          ELSE
            !numerical algorithmic solution is needed here - lib or NumRep?
            STOP "Error: Data - Support for diagonalization of EVB Hamiltonians with >3 states is pending"
          ENDIF

          groundStateEnergy(timestep,step) = MINVAL(eigenvalues(:))

        ENDDO
      ENDDO

      IF (PRESENT(time)) THEN
        CALL CPU_TIME(t2)
        time = t2 - t1
      ENDIF

    END SUBROUTINE ComputeGroundStateEnergy

!*

    SUBROUTINE ComputeMappingEnergies(alpha,time)

      ! #DES: Compute the mapping energies for each pair of states in the system.

      USE Input, ONLY : stateEnergy, coeffs, stateA, stateB, nTimesteps

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: alpha(:)
      REAL(8), INTENT(OUT), OPTIONAL :: time
      REAL(8) :: t1, t2
      INTEGER :: dyn_potential, ene_potential, type, timestep

      IF (PRESENT(time)) CALL CPU_TIME(t1)
      !This is an optimization candidate - not all of the dyn,ene entries are needed
      DO type = 1, SIZE(stateEnergy,4)
        DO dyn_potential = 1, SIZE(stateEnergy,2)
          DO ene_potential = dyn_potential-1, dyn_potential+1
            IF (ene_potential < 1 .OR. ene_potential > SIZE(stateEnergy,2)) CYCLE
            DO timestep = 1, nTimesteps(dyn_potential)
                mappingEnergies(timestep,ene_potential,dyn_potential,type) = (coeffs(timestep,ene_potential,stateA) * (stateEnergy(timestep,dyn_potential,stateA,type) + alpha(stateA)) ) + &
  &                                                                          (coeffs(timestep,ene_potential,stateB) * (stateEnergy(timestep,dyn_potential,stateB,type) + alpha(stateB)) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF (PRESENT(time)) THEN
        CALL CPU_TIME(t2)
        time = t2 - t1
      ENDIF

    END SUBROUTINE ComputeMappingEnergies

!*

    SUBROUTINE ComputeEnergyGap(alpha)

      ! #DES: Compute the value of the energy gap reaction coordinate for each configuration using the input coefficients

      USE Input, ONLY : stateEnergy, stateA, stateB, rcCoeffA, rcCoeffB

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: alpha(:)
      INTEGER :: step
      
      DO step = 1, SIZE(stateEnergy,2)
        energyGap(step,:) = (rcCoeffA * (stateEnergy(:,step,stateA,1) + alpha(stateA))) + (rcCoeffB * (stateEnergy(:,step,stateB,1) + alpha(stateB)))
      ENDDO

    END SUBROUTINE ComputeEnergyGap

!*

    SUBROUTINE ComputeLambdas()

      ! #DES: Compute the value of the FEP order parameter for each FEP step based on the input data
      ! #DES: Can be back-computed from the state coefficients in each mapping potential, provided the
      !       functional form of the order parameter is known.

      USE Input, ONLY : nFepSteps

      IMPLICIT NONE
      INTEGER :: i

      DO i = 1, nFepSteps
        lambda(i) = DBLE(i - 1) / (nFepSteps - 1.0d0)
      ENDDO

    END SUBROUTINE ComputeLambdas

END MODULE Data
