! The goal of the trajectory analysis procedures herein is to estimate the error in each simulation.
! Assuming the data is uncorrelated, the mean, variance, standard deviation and error on the mean are given 
! by simple expressions.
! Assuming the data is correlated, it is necessary to estimate the correlation and attempt to correct the error estimate

MODULE StatisticalFunctions

  ! #DES: Subprograms for statistical analysis of MD data sets

  IMPLICIT NONE
  REAL(8), PARAMETER :: ROOT2PI = 2.50662827463d0

CONTAINS

!*

  PURE REAL(8) FUNCTION mean(x,mask)

    ! #DES: Compute the mean of a 1D array of real(8) numbers
    !       including only those for whom mask(i) is true

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: mask(:)
    INTEGER :: n

    IF (PRESENT(mask)) THEN
      n = COUNT(mask)
    ELSE
      n = SIZE(x)
    ENDIF

    IF (n > 0) THEN
      IF (PRESENT(mask)) THEN
        mean = SUM(x,mask=mask) / n
      ELSE
        mean = SUM(x) / n
      ENDIF
    ELSE
      mean = 0.0d0
    ENDIF

  END FUNCTION mean

!*

  PURE REAL(8) FUNCTION variance(x,mask)

    ! #DES: Compute the variance of a 1D array of real(8) numbers, assuming the entire population is present

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: mask(:)
    INTEGER :: n
    REAL(8) :: difference(SIZE(x))

    IF (PRESENT(mask)) THEN
      n = COUNT(mask)
    ELSE
      n = SIZE(x)
    ENDIF

    IF (n > 0) THEN

      difference(:) = x(:) - mean(x)

      IF (PRESENT(mask)) THEN
        variance = SUM(difference(:)*difference(:),mask=mask) / n
      ELSE
        variance = SUM(difference(:)*difference(:)) / n
      ENDIF
    ELSE
      variance = 0.0d0
    ENDIF

  END FUNCTION variance

!*

  PURE REAL(8) FUNCTION varianceOfMean(x,mask)

    ! #DES: Calculate the variance of the mean assuming no correlation present in data set

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(:)
    LOGICAL, INTENT(IN),OPTIONAL :: mask(:)
    INTEGER :: n

    IF (PRESENT(mask)) THEN
      n = COUNT(mask)
    ELSE
      n = SIZE(x)    
    ENDIF

    IF (n > 0) THEN
      IF (PRESENT(mask)) THEN
        varianceOfMean = variance(x,mask=mask) / n
      ELSE
        varianceOfMean = variance(x) / n
      ENDIF
    ELSE
      varianceOfMean = 0.0d0
    ENDIF

  END FUNCTION varianceOfMean

!*

  FUNCTION histogram(binEdges,input)

    ! #DES: Produce a histogram of the data 'input' using the binEdges provided

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: input(:), binEdges(:)
    INTEGER :: i, bin, nBins
    REAL(8) :: x(SIZE(binEdges)-1)
    INTEGER :: histogram(SIZE(binEdges)-1)

    nBins = SIZE(binEdges) - 1
    x(:) = binEdges(1:nBins) + ((binEdges(1:nBins) - binEdges(2:nBins+1) / 2.0d0))

    ! Count the populations
    histogram(:) = 0
    DO i = 1, SIZE(input)
      DO bin = 1, nBins
        IF ((input(i) >= binEdges(bin)) .AND. (input(i) < binEdges(bin+1))) THEN
          histogram(bin) = histogram(bin) + 1
          EXIT
        ELSE IF (input(i) == binEdges(SIZE(binEdges))) THEN
          histogram(nBins) = histogram(nBins) + 1
          EXIT
        ENDIF
      ENDDO
    ENDDO

  END FUNCTION histogram

!*

  FUNCTION gaussianDistribution(mean, var, input)

    ! #DES: generate values of the Gaussian distribution with 'mean' and 'var' at points in 'input'
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: mean, var, input(:)
    REAL(8) :: gaussianDistribution(SIZE(input))
    REAL(8) :: coeff, var2, diff(SIZE(input))
    INTEGER :: i

    IF (var > 0) THEN

      coeff = SQRT(var) * ROOT2PI
      var2 = 2.0d0 * var

      diff(:) = input(:) - mean
      DO i = 1, SIZE(input)
        gaussianDistribution(i) = EXP(-1.0d0*((diff(i)*diff(i)) / var2))
      END DO
      gaussianDistribution(:) = gaussianDistribution(:) / coeff

    ELSE
      STOP "ERROR - gaussianDistribution - variance is -ve: "
    ENDIF

  END FUNCTION gaussianDistribution

!*

  FUNCTION runningMean(input,mask)

    ! #DES: Compute the running mean of a set of data for visualising equilibration of simulations

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: input(:)
    LOGICAL, INTENT(IN), OPTIONAL :: mask(:)
    REAL(8) :: runningMean(SIZE(input))
    INTEGER :: i

    IF (PRESENT(mask) .EQV. .FALSE.) THEN
      !each output point is the mean after that many timesteps
      DO i = 1, SIZE(input)
        runningMean(i) = mean(input(1:i))
      ENDDO
    ELSE
      DO i = 1, SIZE(input)
        runningMean(i) = mean(input(1:i),mask(1:i))
      ENDDO
    ENDIF

  END FUNCTION runningMean

!*

  ! Generate values of variance of <A> for all values of T
  SUBROUTINE EstimateC(data)

    ! #DES: Compute all autocorrelation functions of the input 'data' and use to compute estimates of var(<A>)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: data(:)
    INTEGER :: cutoff, t, n
    REAL(8) :: var(0:SIZE(data(:))-2), sumTerm(1:SIZE(data(:))-2)
    REAL(8) :: c0, ct, num, den

    n = SIZE(data)

    !Compute the correlation functions and their sum
    c0 = c_t(0,data)

    !Summed term is over correlation functions with t = 1, N-1 (has no 0th entry)
    sumTerm(:) = 0.0d0
    sumTerm(1) = (1.0d0 - (1.0d0/n)) * c_t(1,data)

    DO t = 2, n-2
      ct = c_t(t,data)
      sumTerm(t) = sumTerm(t-1) + ((1.0d0 - (t/DBLE(n))) * ct)
    ENDDO

    var(0) = c0 / (n-1)
    WRITE(21,'(I10,1X,F20.10,I10,F20.10)') 0+1, SQRT(var(0))

    !Compute the variance of the mean as a function of the cutoff: T in {0,1,...,N-1}
    DO cutoff = 1, n-2

      !num / den must always be positive since sqrt is taken
      num = c0 + (2.0d0 * sumTerm(cutoff))
      den = n - 2*cutoff - 1 + ( DBLE(cutoff*(cutoff+1)) / n )
      var(cutoff) = num / den

      IF (var(cutoff) > 0.0d0) THEN
        WRITE(21,'(I10,1X,5F20.10)') cutoff+1, SQRT(var(cutoff))
      ENDIF

    ENDDO

    STOP

  END SUBROUTINE EstimateC

!*

  !Eq. 8 in Flyvbjerg & Petersen, J. Chem. Phys., 91, 461 (1989)
  REAL(8) FUNCTION c_t(t,data)

    ! #DES: Evaluate the t'th autocorrelation function of the data

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: t
    REAL(8), INTENT(IN) :: data(:)
    REAL(8) :: total, average
    INTEGER :: n, k
   
    n = SIZE(data)
    IF ((n - t) > 0) THEN
  
      average = mean(data(:))
      total = 0.0d0
      DO k = 1, n-t
        total = total + ((data(k) - average) * (data(k+t) - average))
      ENDDO
      c_t = total / DBLE(n-t)

    ELSE
      STOP "ERROR - c_t - Infinite correlation function"
    ENDIF

  END FUNCTION c_t

!*

  !Use the FP method to determine the standard deviation on the mean of correlated data
  !Can also return sigma(<A>) as function of num. of block transforms for visual inspection
  REAL(8) FUNCTION FlyvbjergPetersen(x)

    ! #DES: Use the Flyvjerg-Petersen method to determine the stdev on the mean of the data set as a function of T

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(:)
    REAL(8), ALLOCATABLE :: xPrime(:), copy(:)
    INTEGER :: nPrime, step
    REAL(8) :: T(SIZE(x)), sigma(SIZE(x)), error(SIZE(x))
    REAL(8) :: finalValue, c0prime
    LOGICAL :: converged
    !To store data per blocking transform need to know how many can be done.

    nPrime = SIZE(x)

    converged = .FALSE.
    finalValue = -1.0d0
    ALLOCATE(xPrime(SIZE(x)))
    xPrime(:) = x(:)

    step = 1
    DO WHILE (nPrime >= 2)

      c0prime = c0(xPrime)
      T(step) = 0.5d0 * ((2.0d0**(step-1)) - 1)
      sigma(step) = StDevMean(c0prime,nPrime)
      error(step) = StDevMeanError(c0prime,nPrime)

      !WRITE(*,*) T(step) + 1, sigma(step), error(step)
      IF ((converged .EQV. .FALSE.) .AND. (ABS(sigma(step) - sigma(step-1)) < error(step-1))) THEN
        finalValue = sigma(step)
        converged = .TRUE.
      ENDIF

      ALLOCATE(copy(SIZE(xPrime))); copy(:) = xPrime(:)
      DEALLOCATE(xPrime); ALLOCATE(xPrime(nPrime/2))
      xPrime = BlockTransform(copy(1:nPrime-MOD(nPrime,2))) !discard the last point if nPrime is odd
      DEALLOCATE(copy)
      nPrime = SIZE(xPrime)
      step = step + 1
      
    ENDDO

    DEALLOCATE(xPrime)

    FlyvbjergPetersen = finalValue

  END FUNCTION FlyvbjergPetersen

!*

  FUNCTION BlockTransform(x)

    ! #DES: Perform the blocking transform on the data (Eq. 20, 21)
    
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(:)
    REAL(8) :: BlockTransform(SIZE(x)/2) !rounds down
    INTEGER :: n, i

      n = SIZE(x)
      IF (MOD(n,2) /= 0) STOP "ERROR - BlockTransform"
      DO i = 1, n/2 !rounds down
        BlockTransform(i) = 0.5d0 * (x((2*i)-1)+x(2*i))
      ENDDO

  END FUNCTION BlockTransform

!*

  PURE REAL(8) FUNCTION c0(x)

    ! #DES: Compute the 0th autocorrelation function (the variance)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x(:)
    REAL(8) :: difference(SIZE(x))
    INTEGER :: n

    n = SIZE(x)
    IF (n > 0) THEN
      difference(:) = x(:) - mean(x)
      c0 = SUM(difference*difference) / n
    ELSE
      c0 = 0.0d0
    ENDIF

  END FUNCTION c0

!*

  !Compute the standard deviation of the mean from the FP quantities
  !Domain of function is (0,+inf)
  REAL(8) FUNCTION StDevMean(c_0,n)

    ! #DES: Compute the stdev of the mean from the FP quantities

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: c_0
    INTEGER, INTENT(IN) :: n

    IF (((n-1) > 0) .AND. (c_0 >= 0)) THEN
      StDevMean = SQRT(c_0 / (n-1))
    ELSE
      STOP "ERROR - StDevMean - negative variance"
    ENDIF  

  END FUNCTION StDevMean

!*

  REAL(8) FUNCTION StDevMeanError(c_0,n)

    ! #DES: Compute error in the stdev of the mean from FP quantities

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: c_0
    INTEGER, INTENT(IN) :: n

    IF ((n-1) > 0) THEN
      StDevMeanError = StDevMean(c_0,n) * (1.0d0 / SQRT(2.0d0 * (n-1)))
    ELSE
      STOP "ERROR - StDevMeanError - sqrt of -ve number"
    ENDIF

  END FUNCTION StDevMeanError

END MODULE StatisticalFunctions
