! As the derivative of the object function is not available, the optimiser must need only
! evaluations of the function to be minimized. this makes selection of an algorithm simpler,
! although potential benefits of methods that use the derivative are not available.
! The downhill simplex method is robust but potentially slow. Given that the optimization of this
! particular object function does not involve lots of variables it is unlikely a faster or more complicated
! (i.e. direction set method or PSO) algorithm is needed. It is short and self-contained.

MODULE DownhillSimplex

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RunNelderMead

  NAMELIST /opt/ c_reflect, c_expand, c_contract, c_shrink, threshold, maxSteps
  REAL(8) :: c_reflect, c_expand, c_contract, c_shrink, threshold
  INTEGER :: maxSteps

  CONTAINS

!*

  ! needs to take an object function as an argument to pass to the nelder-mead optimizer
  ! right now is just importing a fixed object function from module ObjectFunctions
  SUBROUTINE RunNelderMead(guess,lambda,logUnit,printDetails)

    ! #DES: Driver for Nelder-Mead optimizer, checks sanity of input, creates guess and informs the log.
    !       guess is the point around which the simplex should be built
    !       lambda is the characteristic length scale in each dimension
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: guess(:), lambda(:)
    INTEGER, INTENT(IN) :: logUnit
    LOGICAL, INTENT(IN) :: printDetails
    REAL(8) :: simplex(SIZE(guess)+1,SIZE(guess))
    REAL(8) :: optimum(SIZE(guess)), optValue
    LOGICAL :: converged = .FALSE.
    INTEGER :: i

    WRITE(logUnit,'(A)') "*   Initializing Nelder-Mead Downhill Simplex   *"
    WRITE(logUnit,*)
    WRITE(logUnit,'(A24,I6)')   "Parameter Space Size  = ", SIZE(guess)
    WRITE(logUnit,'(A24,I6)')   "Points in Simplex     = ", SIZE(guess) + 1

    CALL ProcessNameList(logUnit)
    CALL GenerateSimplex(guess,lambda,logUnit,simplex)
    CALL NelderMead(simplex,logUnit,printDetails,optimum,optValue,converged)

    IF (converged .EQV. .TRUE.) THEN
      WRITE(logUnit,'(A)') "Downhill Simplex Converged"
    ELSE
      WRITE(logUnit,'(A)') "Downhill Simplex Not Converged"
    ENDIF

    WRITE(logUnit,'(A,F15.8,A)',ADVANCE="NO") "Best value located = ", optValue, " at:"
    DO i = 1, SIZE(optimum)
      WRITE(logUnit,'(F15.8,1X)',ADVANCE='NO') optimum(i)
    ENDDO
    WRITE(logUnit,*)
    WRITE(logUnit,*)

  END SUBROUTINE RunNelderMead

!*

  SUBROUTINE GenerateSimplex(guess,lambda,logUnit,simplex)

    USE ObjectFunctions, ONLY : objectFunction
    IMPLICIT NONE

    REAL(8), INTENT(IN)  :: guess(:), lambda(:)
    INTEGER, INTENT(IN)  :: logUnit
    REAL(8), INTENT(OUT) :: simplex(SIZE(guess)+1,SIZE(guess))
    REAL(8) :: functionValues(SIZE(simplex,1)+1)      ! object function evaluated at each point
    INTEGER :: m, n, i

    n = SIZE(guess)
    m = n + 1

    !generate a simplex from the input point - should be done in own routine/function
    simplex(:,:)      = 0.0d0
    functionValues(:) = 0.0d0
    simplex(1,:)      = guess(:)
    functionValues(1) = objectFunction(simplex(1,:))

    DO i = 2, m
      simplex(i,:) = guess(:)
      simplex(i,i-1) = simplex(i,i-1) + lambda(i-1)
      functionValues(i) = objectFunction(simplex(i,:))
    ENDDO

    WRITE(logUnit,'(A16)') "Initial simplex:"
    CALL WriteSimplex(simplex,functionValues,logUnit)

  END SUBROUTINE GenerateSimplex

!*

  SUBROUTINE NelderMead(initialSimplex,logUnit,debug,optimum,optValue,converged)

    USE ObjectFunctions, ONLY : objectFunction
    IMPLICIT NONE

    REAL(8), INTENT(IN)  :: initialSimplex(:,:) ! m x n, each row is an n-dimensional vector
    INTEGER, INTENT(IN)  :: logUnit
    LOGICAL, INTENT(IN)  :: debug
    REAL(8), INTENT(OUT) :: optimum(SIZE(initialSimplex,2))
    REAL(8), INTENT(OUT) :: optValue
    LOGICAL, INTENT(OUT) :: converged

    REAL(8) :: simplex(SIZE(initialSimplex,1),SIZE(initialSimplex,2))
    REAL(8) :: functionValues(SIZE(initialSimplex,1)+1)      ! object function evaluated at each point
    INTEGER :: step, n, m, i
    REAL(8) :: f_reflect, f_expand, f_contract
    REAL(8) :: reflectedPoint(SIZE(initialSimplex,2)),  expandedPoint(SIZE(initialSimplex,2))
    REAL(8) :: contractedPoint(SIZE(initialSimplex,2)), centroidPoint(SIZE(initialSimplex,2))

    converged    = .FALSE.
    optimum(:)   = 0.0d0
    optValue      = 0.0d0
    simplex(:,:) = initialSimplex(:,:)
    n            = SIZE(simplex,2)
    m            = n + 1

    DO i = 1, m
      functionValues(i) = objectFunction(simplex(i,:))
    ENDDO

    !repeat until convergence or excessive compute time
    DO step = 1, maxSteps

      CALL BubbleSort(simplex,functionValues)  ! sort ascending (f(x1) smallest, f(n+1) largest)
      IF (debug .EQV. .TRUE.) WRITE(logUnit,'(A5,I0.4,A8,F10.4)') "STEP ", step, " f(x) = ", functionValues(1)

      IF (functionValues(1) < threshold) THEN ! convergence test

        converged = .TRUE.
        optimum(:) = simplex(1,:)
        optValue = functionValues(1)
        RETURN

      ENDIF

      centroidPoint(:) = centroid(simplex(1:n,:))

      !reflect
      reflectedPoint(:) = centroidPoint(:) + (c_reflect * (centroidPoint(:)-simplex(m,:)))
      f_reflect = objectFunction(reflectedPoint)

      IF ((f_reflect >= functionValues(1)) .AND. (f_reflect < functionValues(n))) THEN

        IF (debug .EQV. .TRUE.) WRITE(logUnit,*) "REFLECTING"
        simplex(m,:) = reflectedPoint(:)
        functionValues(m) = f_reflect
        CYCLE

      ELSE IF (f_reflect < functionValues(1)) THEN

        !expand
        expandedPoint(:) = centroidPoint(:) + (c_expand * (reflectedPoint(:) - centroidPoint(:)))
        f_expand = objectFunction(expandedPoint)
        IF (f_expand < f_reflect) THEN
          IF (debug .EQV. .TRUE.) WRITE(logUnit,*) "EXPANDING"
          simplex(m,:) = expandedPoint(:)
          functionValues(m) = f_expand
        ELSE
          IF (debug .EQV. .TRUE.) WRITE(logUnit,*) "REFLECTING"
          simplex(m,:) = reflectedPoint(:)
          functionValues(m) = f_reflect
        ENDIF
        CYCLE

      ELSE
        !contract
        contractedPoint(:) = centroidPoint(:) + (c_contract * (simplex(m,:) - centroidPoint(:)))
        f_contract = objectFunction(contractedPoint)
        IF (f_contract < functionValues(m)) THEN
          IF (debug .EQV. .TRUE.) WRITE(logUnit,*) "CONTRACTING"
          simplex(m,:) = contractedPoint(:)
          functionValues(m) = f_contract
          CYCLE
        ELSE
          !shrink
          IF (debug .EQV. .TRUE.) WRITE(logUnit,*) "SHRINKING"
          DO i = 2, m
            simplex(i,:) = simplex(1,:) + c_shrink * (simplex(i,:) - simplex(1,:))
          ENDDO
        ENDIF

      ENDIF

    ENDDO

    ! If optimization fails, return best value located
    optimum(:) = simplex(1,:)
    optValue = functionValues(1)
    RETURN

END SUBROUTINE NelderMead

!*

  SUBROUTINE BubbleSort(simplex,functionValues)

    ! #DES: Basic bubble sort algorithm, can be replaced if necessary.

    IMPLICIT NONE
    REAL(8) :: simplex(:,:), functionValues(:)
    reaL(8) :: f_temp, x_temp(SIZE(simplex,2))
    INTEGER :: i, j, m

    m = SIZE(simplex,1)
    DO i = 1, m !loop over vectors
      DO j = m, i+1, -1
        IF (functionValues(j-1) > functionValues(j)) THEN

          f_temp = functionValues(j-1)
          functionValues(j-1) = functionValues(j)
          functionValues(j) = f_temp

          x_temp(:) = simplex(j-1,:)
          simplex(j-1,:) = simplex(j,:)
          simplex(j,:) = x_temp(:)

        ENDIF
      ENDDO
    ENDDO       
    
  END SUBROUTINE BubbleSort

!*

  FUNCTION centroid(input)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: input(:,:) !each row a vector
    REAL(8) :: centroid(SIZE(input,2))
    INTEGER :: i, j, m, n

    m = SIZE(input,1)
    n = SIZE(input,2)

    centroid(:) = 0.0d0
    DO i = 1, m
      DO j = 1, n
        centroid(j) = centroid(j) + input(i,j)
      ENDDO
    ENDDO

    centroid(:) = centroid(:) / SIZE(input,2)

  END FUNCTION centroid

!*

  SUBROUTINE WriteSimplex(simplex,functionValues,unit)

    ! #DES: Write the provided simplex coords/values to unit.

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    REAL(8), INTENT(IN) :: simplex(:,:), functionValues(:)
    INTEGER :: i, j

    DO i = 1, SIZE(simplex,1)
      WRITE(unit,'(I0.3,1X)',ADVANCE='NO') i
      DO j = 1, SIZE(simplex,2)
        WRITE(unit,'(F10.3,1X)',ADVANCE='NO') simplex(i,j)
      ENDDO
      WRITE(unit,'(F10.3)') functionValues(i)
    ENDDO
    WRITE(unit,*)

  END SUBROUTINE WriteSimplex

!*

  SUBROUTINE ProcessNameList(logUnit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit

    CALL ReadNameList(logUnit)
    CALL CheckNameList()
    CALL PrintNameList(logUnit)

  END SUBROUTINE ProcessNameList

!*

  SUBROUTINE ReadNameList(logUnit)

    ! #DES: Open, read and close the namelist file opt.nml

    USE FileIO, ONLY : OpenFile, CloseFile
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit
    INTEGER, PARAMETER :: nmlUnit = 20
    LOGICAL :: fileRead

    CALL OpenFile(nmlUnit,"opt.nml","read",success=fileRead)

    !standard algorithm parameter values
    c_reflect  = 1.0d0
    c_expand   = 2.0d0
    c_contract = 0.5d0
    c_shrink   = 0.5d0
    threshold  = 1.0d-10
    maxSteps   = 100

    IF (fileRead .EQV. .TRUE.) THEN
      WRITE(logUnit,'(A)') "Optimization parameters will be read from opt.nml"
      READ(nmlUnit,NML=opt)
      CALL CloseFile(nmlUnit)
    ELSE
      WRITE(logUnit,'(A)') "Namelist opt.nml could not be read - default optimization parameters will be used"
    ENDIF

  END SUBROUTINE ReadNameList

!*

  SUBROUTINE CheckNameList()

    ! There must be additional constraints on these values - expand > 1, 0 < contract < 1 CHECK 

    IMPLICIT NONE

    IF (threshold  < 0.0d0) STOP "Illegal (-ve) value of threshold in DownhillSimplex:CheckNameList"
    IF (c_reflect  < 0.0d0) STOP "Illegal (-ve) value of c_reflect in DownhillSimplex:CheckNameList"
    IF (c_expand   < 0.0d0) STOP "Illegal (-ve) value of c_expand in DownhillSimplex:CheckNameList"
    IF (c_contract < 0.0d0) STOP "Illegal (-ve) value of c_contract in DownhillSimplex:CheckNameList"
    IF (c_shrink   < 0.0d0) STOP "Illegal (-ve) value of c_shrink in DownhillSimplex:CheckNameList"
    IF (maxSteps   < 1)     STOP "Illegal (-ve) value of maxSteps in DownhillSimplex:CheckNameList"

  END SUBROUTINE CheckNameList

!*

  SUBROUTINE PrintNameList(logUnit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: logUnit

    WRITE(logUnit,'(A)') "*    Optimization Parameters    *"
    WRITE(logUnit,*) ""
    WRITE(logUnit,'(A,F6.2)') "Reflection Coefficient (alpha)  = ", c_reflect
    WRITE(logUnit,'(A,F6.2)') "Expansion Coefficient (beta)    = ", c_expand
    WRITE(logUnit,'(A,F6.2)') "Contraction Coefficient (gamma) = ", c_contract
    WRITE(logUnit,'(A,F6.2)') "Shrinking Coefficient (sigma)   = ", c_shrink
    WRITE(logUnit,'(A,F6.2)') "Convergence Threshold           = ", threshold
    WRITE(logUnit,'(A,I6)')   "Maximum Steps                   = ", maxSteps
    WRITE(logUnit,*)

  END SUBROUTINE PrintNameList

END MODULE DownhillSimplex





!  SUBROUTINE WriteVector(input,unit)
!    ! #DES: Print a vector to a single line on the provided unit.
!    IMPLICIT NONE
!    INTEGER, INTENT(IN) :: unit
!    REAL(8), INTENT(IN) :: input(:)
!    INTEGER :: i
!    DO i = 1, SIZE(input)
!      WRITE(unit,'(F6.2)',ADVANCE='NO') input(i)
!    ENDDO
!    WRITE(unit,*)
!  END SUBROUTINE WriteVector

