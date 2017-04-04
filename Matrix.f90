MODULE Matrix

  ! #DES: Module containing routines for manipulating matrices, focused on explicit solutions
  !       for 2- and 3-dimensional matrices. 

  IMPLICIT NONE

CONTAINS

  FUNCTION Eigenvalues2DRealSymmetric(A)

    ! #DES: Explicit solution for eigenvalues of a 2D real, symmetric matrix
    ! #DES: Obtained by solving the appropriate secular equation, det(A - lI) = 0, (a quadratic with 2 real roots)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(2,2)
    REAL(8) :: eigenvalues2DRealSymmetric(2)
    REAL(8) :: sum, diff, sqrtD !discriminant

    IF (A(1,2) /= A(2,1)) STOP "Error: Matrix - nonsymmetric matrix passed to Eigenvalues2DRealSymmetric"

    IF (A(1,2) == 0.0d0) THEN

      !trivial solution - matrix is already diagonal
      eigenvalues2DRealSymmetric(1) = A(1,1)
      eigenvalues2DRealSymmetric(2) = A(2,2)

    ELSE

      sum  = A(1,1) + A(2,2)
      diff = A(1,1) - A(2,2)
      sqrtD = SQRT(diff*diff + 4.0d0*A(1,2)*A(1,2))

      eigenvalues2DRealSymmetric(1) = 0.5d0 * (sum + sqrtD)
      eigenvalues2DRealSymmetric(2) = 0.5d0 * (sum - sqrtD)

    ENDIF

  END FUNCTION Eigenvalues2DRealSymmetric

!*

  FUNCTION Eigenvalues3DRealSymmetric(A)

    !https://gitlab.dune-project.org/lars.lubkoll/dune-common/commit/3e91e19881ec7a8708ba2dab17d7f8694ddb7e08
    !https://en.wikipedia.org/wiki/Eigenvalue_algorithm#Direct_calculation

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(3,3)
    REAL(8) :: eigenvalues3DRealSymmetric(3)
    REAL(8), PARAMETER :: pi = 3.14159d0
    REAL(8) :: p, p1, q, p2, r, phi, diff
    REAL(8) :: I(3,3), B(3,3)
    INTEGER :: j

    p1 = A(1,2)*A(1,2) + A(1,3)*A(1,3) + A(2,3)*A(2,3)

    IF (p1 == 0.0d0) THEN

      !A is already diagonal so return its diagonal entries
      DO j = 1, 3
        eigenvalues3DRealSymmetric(j) = A(j,j)
      ENDDO

    ELSE

      I(:,:) = 0.0d0
      DO j = 1, 3
        I(j,j) = 1.0d0
      ENDDO

      q = trace(A) / 3.0d0
      p2 = 2.0d0 * p1
      DO j = 1, 3
        diff = A(j,j) - q
        p2 = p2 + (diff*diff)
      ENDDO
      p = SQRT(p2 / 6.0d0)

      B(:,:) = (1.0d0 / p) * (A(:,:) - q * I(:,:)) ! I is the identity matrix

      r = determinant(B) / 2.0d0

      ! In exact arithmetic for a symmetric matrix  -1 <= r <= 1
      ! but computation error can leave it slightly outside this range.
      IF (r <= -1.0d0) THEN
        phi = pi / 3.0d0
      ELSEIF (r >= 1.0d0) THEN
        phi = 0.0d0
      ELSE
        phi = acos(r) / 3.0d0
      END IF

      ! the eigenvalues satisfy eig3 <= eig2 <= eig1
      eigenvalues3DrealSymmetric(1) = q + (2 * p * COS(phi))
      eigenvalues3DRealSymmetric(3) = q + (2 * p * COS(phi + (2*pi/3)))
      ! since trace(A) = eig1 + eig2 + eig3
      eigenvalues3DRealSymmetric(2) = (3.0d0 * q) - eigenvalues3DRealSymmetric(1) - eigenvalues3DRealSymmetric(3)

    END IF

  END FUNCTION Eigenvalues3DRealSymmetric

!*

  REAL(8) FUNCTION determinant(A)

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(3,3)

    determinant = 0.0d0
    determinant = determinant + (A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)))
    determinant = determinant + (A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)))
    determinant = determinant + (A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1)))

  END FUNCTION determinant

!*

  REAL(8) FUNCTION trace(A)

    ! #DES: Compute trace of an n-by-n square matrix, the sum of the elements on the main diagonal

    IMPLICIT NONE
    REAL(8), INTENT(IN) :: A(:,:)
    INTEGER :: i

    IF (SIZE(A,1) /= SIZE(A,2)) STOP "Error: Matrix - non-square matrix received as input to trace"

    trace = 0.0d0
    DO i = 1, SIZE(A,1)
      trace = trace + A(i,i)
    ENDDO

  END FUNCTION trace

END MODULE Matrix
