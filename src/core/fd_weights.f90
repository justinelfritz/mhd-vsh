MODULE FD_WEIGHTS
!> Generic finite-difference weight generator (Fornberg 1988, "Calculation
!> of Weights in Finite Difference Formulas", SIAM Review 40(3):685-691).
!> Given an arbitrary (possibly nonuniform) set of stencil node locations,
!> produces weights for every derivative order up to MMAX simultaneously in
!> one O(N**2) pass. This single routine is used for interior stencils,
!> one-sided boundary closures, and (via RADIAL_OPERATORS' mirrored ghost
!> nodes) the full-sphere origin-regularity rows -- uniform and stretched
!> grids alike.
USE KINDS, ONLY: dp, i4
IMPLICIT NONE
PRIVATE
PUBLIC :: FD_WEIGHTS_GENERIC

CONTAINS

!> @param X0 Point at which the derivatives are evaluated.
!> @param X Stencil node locations, size N+1 (0-indexed 0..N).
!> @param N One less than the number of stencil nodes.
!> @param MMAX Highest derivative order to produce (0=interpolation).
!> @param C Output weights, C(i,m) is the weight on X(i) for the m-th
!>   derivative, shape (0:N,0:MMAX).
SUBROUTINE FD_WEIGHTS_GENERIC(X0, X, N, MMAX, C)
  REAL(KIND=dp),    INTENT(IN)  :: X0
  INTEGER(KIND=i4), INTENT(IN)  :: N, MMAX
  REAL(KIND=dp),    INTENT(IN)  :: X(0:N)
  REAL(KIND=dp),    INTENT(OUT) :: C(0:N,0:MMAX)
  INTEGER(KIND=i4) :: I, J, K, MN
  REAL(KIND=dp) :: C1, C2, C3, C4, C5

  C = 0.0_dp
  C1 = 1.0_dp
  C4 = X(0) - X0
  C(0,0) = 1.0_dp

  DO I = 1, N
    MN = MIN(I, MMAX)
    C2 = 1.0_dp
    C5 = C4
    C4 = X(I) - X0
    DO J = 0, I-1
      C3 = X(I) - X(J)
      C2 = C2*C3
      IF (J == I-1) THEN
        DO K = MN, 1, -1
          C(I,K) = C1*(REAL(K, KIND=dp)*C(I-1,K-1) - C5*C(I-1,K))/C2
        END DO
        C(I,0) = -C1*C5*C(I-1,0)/C2
      END IF
      DO K = MN, 1, -1
        C(J,K) = (C4*C(J,K) - REAL(K, KIND=dp)*C(J,K-1))/C3
      END DO
      C(J,0) = C4*C(J,0)/C3
    END DO
    C1 = C2
  END DO
END SUBROUTINE FD_WEIGHTS_GENERIC

END MODULE FD_WEIGHTS
