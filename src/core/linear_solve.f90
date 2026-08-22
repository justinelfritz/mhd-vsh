MODULE LINEAR_SOLVE
!> Generic dense linear-system solver: factorize once, solve for
!> arbitrarily many right-hand sides against the stored factors. Wraps
!> LAPACK's DGETRF (pivoted LU factorization) and DGETRS (triangular
!> solve using those factors).
!>
!> This module has no knowledge of spherical harmonics, radial
!> operators, curvature, boundary conditions, or any physics
!> whatsoever -- A and B are just numbers to it. Every regime is free to
!> assemble its own system matrix however its own PDE requires (see
!> RADIAL_OPERATORS/BOUNDARY_CONDITIONS for optional building blocks:
!> D1/D2, ADD_CURVATURE_TERM, APPLY_ORIGIN_REGULARITY,
!> APPLY_VACUUM_BC_POLOIDAL/TOROIDAL) and only hands the finished dense
!> matrix to this module at the very last step.
!>
!> Factorizing once and solving many right-hand sides against the same
!> factors matters here specifically because a per-l radial system is
!> shared across every m at that l (2l+1 right-hand sides, one DGETRS
!> call) and typically across many timesteps too, as long as the
!> caller's own coefficients (dt, diffusivity, ...) stay fixed.
USE KINDS, ONLY: dp, i4
IMPLICIT NONE
PRIVATE
PUBLIC :: DENSE_FACTORS_T, FACTORIZE_DENSE, SOLVE_FACTORED

INTERFACE
  SUBROUTINE DGETRF(M, N, A, LDA, IPIV, INFO)
    IMPORT :: dp, i4
    INTEGER(KIND=i4), INTENT(IN)    :: M, N, LDA
    REAL(KIND=dp),    INTENT(INOUT) :: A(LDA,*)
    INTEGER(KIND=i4), INTENT(OUT)   :: IPIV(*)
    INTEGER(KIND=i4), INTENT(OUT)   :: INFO
  END SUBROUTINE DGETRF

  SUBROUTINE DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
    IMPORT :: dp, i4
    CHARACTER,        INTENT(IN)    :: TRANS
    INTEGER(KIND=i4), INTENT(IN)    :: N, NRHS, LDA, LDB
    REAL(KIND=dp),    INTENT(IN)    :: A(LDA,*)
    INTEGER(KIND=i4), INTENT(IN)    :: IPIV(*)
    REAL(KIND=dp),    INTENT(INOUT) :: B(LDB,*)
    INTEGER(KIND=i4), INTENT(OUT)   :: INFO
  END SUBROUTINE DGETRS
END INTERFACE

!> Pivoted-LU factorization of a dense (N,N) matrix, ready for repeated
!> SOLVE_FACTORED calls. INFO mirrors DGETRF's: 0 on success, >0 means
!> the factorization completed but U is exactly singular (the caller's
!> assembled system is degenerate -- e.g. a missing or inconsistent
!> boundary row), the caller's to check if it cares.
TYPE :: DENSE_FACTORS_T
  INTEGER(KIND=i4) :: N    = 0
  INTEGER(KIND=i4) :: INFO = 0
  REAL(KIND=dp),    ALLOCATABLE :: LU(:,:)
  INTEGER(KIND=i4), ALLOCATABLE :: IPIV(:)
END TYPE DENSE_FACTORS_T

CONTAINS

!> @param FACTORS Output pivoted-LU factorization of A.
!> @param A Dense system matrix, shape (N,N); not modified (factorized
!>   into an internal copy).
SUBROUTINE FACTORIZE_DENSE(FACTORS, A)
  TYPE(DENSE_FACTORS_T), INTENT(OUT) :: FACTORS
  REAL(KIND=dp),         INTENT(IN)  :: A(:,:)
  INTEGER(KIND=i4) :: N
  N = SIZE(A,1)
  FACTORS%N = N
  ALLOCATE(FACTORS%LU(N,N))
  ALLOCATE(FACTORS%IPIV(N))
  FACTORS%LU = A
  CALL DGETRF(N, N, FACTORS%LU, N, FACTORS%IPIV, FACTORS%INFO)
END SUBROUTINE FACTORIZE_DENSE

!> @param X Output solution, shape (N,NRHS): `A*X(:,k) = B(:,k)` for every
!>   column k.
!> @param FACTORS Pivoted-LU factorization from FACTORIZE_DENSE; reused
!>   unmodified, so the same FACTORS may be solved against repeatedly.
!> @param B Right-hand side(s), shape (N,NRHS); any number of columns,
!>   solved against the same stored factorization in one call.
SUBROUTINE SOLVE_FACTORED(X, FACTORS, B)
  REAL(KIND=dp),         INTENT(OUT) :: X(:,:)
  TYPE(DENSE_FACTORS_T), INTENT(IN)  :: FACTORS
  REAL(KIND=dp),         INTENT(IN)  :: B(:,:)
  INTEGER(KIND=i4) :: NRHS, INFO
  NRHS = SIZE(B,2)
  X = B
  CALL DGETRS('N', FACTORS%N, NRHS, FACTORS%LU, FACTORS%N, FACTORS%IPIV, X, FACTORS%N, INFO)
END SUBROUTINE SOLVE_FACTORED

END MODULE LINEAR_SOLVE
