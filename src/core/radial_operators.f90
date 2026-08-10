MODULE RADIAL_OPERATORS
!> Dense, l-independent radial differentiation operators (D1=d/dr,
!> D2=d^2/dr^2) built once per grid via FD_WEIGHTS_GENERIC, plus the
!> full-sphere origin-regularity correction. The l(l+1)/r**2 curvature
!> term and boundary conditions are deliberately NOT included here -- they
!> are applied later, at system-assembly time, so this module stays a
!> reusable, regime-agnostic building block.
USE KINDS,       ONLY: dp, i4
USE GRID_RADIAL, ONLY: RADIAL_GRID_T
USE FD_WEIGHTS,  ONLY: FD_WEIGHTS_GENERIC
IMPLICIT NONE
PRIVATE
PUBLIC :: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS, APPLY_ORIGIN_REGULARITY, &
          ADD_CURVATURE_TERM

TYPE :: RADIAL_OPERATOR_T
  INTEGER(KIND=i4) :: N = 0, FD_ORDER = 4
  REAL(KIND=dp), ALLOCATABLE :: D1(:,:), D2(:,:)   ! (N,N)
END TYPE RADIAL_OPERATOR_T

CONTAINS

!> Builds the baseline D1/D2 using centered stencils in the interior and
!> widened one-sided closures near either physical boundary (rows that a
!> shell regime will typically overwrite entirely via boundary-row
!> replacement, and that a full-sphere regime will overwrite near r=0 via
!> APPLY_ORIGIN_REGULARITY).
!>
!> @param OPS Output operator set.
!> @param RGRID Radial grid to differentiate on.
SUBROUTINE BUILD_RADIAL_OPERATORS(OPS, RGRID)
  TYPE(RADIAL_OPERATOR_T), INTENT(OUT) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN)  :: RGRID
  INTEGER(KIND=i4) :: N, HALF, I, LO, HI, WIDTH, K
  REAL(KIND=dp), ALLOCATABLE :: XNODES(:), C(:,:)

  N = RGRID%N
  HALF = RGRID%FD_ORDER/2

  OPS%N = N
  OPS%FD_ORDER = RGRID%FD_ORDER
  ALLOCATE(OPS%D1(N,N), OPS%D2(N,N))
  OPS%D1 = 0.0_dp
  OPS%D2 = 0.0_dp

  DO I = 1, N
    LO = MAX(1, I-HALF)
    HI = MIN(N, I+HALF)
    IF (HI - LO + 1 < 2*HALF+1) THEN
      IF (LO == 1) THEN
        HI = MIN(N, LO + 2*HALF)
      ELSE
        LO = MAX(1, HI - 2*HALF)
      END IF
    END IF
    WIDTH = HI - LO
    ALLOCATE(XNODES(0:WIDTH), C(0:WIDTH,0:2))
    XNODES = RGRID%R(LO:HI)
    CALL FD_WEIGHTS_GENERIC(RGRID%R(I), XNODES, WIDTH, 2, C)
    DO K = 0, WIDTH
      OPS%D1(I, LO+K) = C(K,1)
      OPS%D2(I, LO+K) = C(K,2)
    END DO
    DEALLOCATE(XNODES, C)
  END DO
END SUBROUTINE BUILD_RADIAL_OPERATORS

!> Rebuilds the D1/D2 rows nearest r=0 for a full-sphere grid, using an
!> augmented stencil that includes mirrored ghost nodes at negative radius
!> and the regularity closure f(-r)=(-1)**l f(r) (exact for a smooth field
!> at fixed (l,m), from the parity of Y_l^m under (r,theta,phi) ->
!> (-r,pi-theta,phi+pi)). No-op for a shell grid.
!>
!> @param D1 In/out D1 operator, modified in rows 1..FD_ORDER/2.
!> @param D2 In/out D2 operator, modified in rows 1..FD_ORDER/2.
!> @param L Spherical-harmonic degree the regularity condition applies to.
!> @param RGRID Radial grid (must be the one D1/D2 were built on).
SUBROUTINE APPLY_ORIGIN_REGULARITY(D1, D2, L, RGRID)
  REAL(KIND=dp),        INTENT(INOUT) :: D1(:,:), D2(:,:)
  INTEGER(KIND=i4),     INTENT(IN)    :: L
  TYPE(RADIAL_GRID_T),  INTENT(IN)    :: RGRID
  INTEGER(KIND=i4) :: HALF, I, K, M, P, NVIRT, COL
  REAL(KIND=dp)    :: SIGN_L
  REAL(KIND=dp), ALLOCATABLE :: XNODES(:), C(:,:)
  INTEGER(KIND=i4), ALLOCATABLE :: COLIDX(:)

  IF (.NOT. RGRID%FULL_SPHERE) RETURN
  HALF = RGRID%FD_ORDER/2
  SIGN_L = 1.0_dp - 2.0_dp*REAL(MOD(L,2), KIND=dp)

  DO I = 1, HALF
    NVIRT = HALF - I + 1
    ALLOCATE(XNODES(0:2*HALF), C(0:2*HALF,0:2), COLIDX(0:2*HALF))
    K = 0
    DO M = I-HALF, 0
      P = 1 - M
      XNODES(K) = -RGRID%R(P)
      COLIDX(K) = P
      K = K + 1
    END DO
    DO M = 1, I+HALF
      XNODES(K) = RGRID%R(M)
      COLIDX(K) = M
      K = K + 1
    END DO

    CALL FD_WEIGHTS_GENERIC(RGRID%R(I), XNODES, 2*HALF, 2, C)

    D1(I,:) = 0.0_dp
    D2(I,:) = 0.0_dp
    DO K = 0, 2*HALF
      COL = COLIDX(K)
      IF (K < NVIRT) THEN
        D1(I,COL) = D1(I,COL) + SIGN_L*C(K,1)
        D2(I,COL) = D2(I,COL) + SIGN_L*C(K,2)
      ELSE
        D1(I,COL) = D1(I,COL) + C(K,1)
        D2(I,COL) = D2(I,COL) + C(K,2)
      END IF
    END DO
    DEALLOCATE(XNODES, C, COLIDX)
  END DO
END SUBROUTINE APPLY_ORIGIN_REGULARITY

!> Adds the spherical-harmonic curvature term -l(l+1)/r**2 to the
!> diagonal of a caller-owned radial operator -- e.g. so a regime that
!> wants D2 - l(l+1)/r**2 (the radial part of the scalar Laplacian's
!> angular eigenvalue at degree l) can build it from a per-l copy of
!> RADIAL_OPERATOR_T%D2 plus this call. Purely geometric: the same term
!> at every degree l regardless of what physics the caller multiplies it
!> by, so it stays an optional building block here rather than being
!> assumed by any solver -- a regime with a different operator entirely
!> is free to skip it.
!>
!> @param OP In/out operator (a per-l copy of whatever radial operator
!>   the caller is assembling), diagonal modified in place.
!> @param L Spherical-harmonic degree, l>=0.
!> @param RGRID Radial grid OP was built on.
SUBROUTINE ADD_CURVATURE_TERM(OP, L, RGRID)
  REAL(KIND=dp),        INTENT(INOUT) :: OP(:,:)
  INTEGER(KIND=i4),     INTENT(IN)    :: L
  TYPE(RADIAL_GRID_T),  INTENT(IN)    :: RGRID
  INTEGER(KIND=i4) :: I
  REAL(KIND=dp)    :: LAMBDA_L
  LAMBDA_L = REAL(L*(L+1), KIND=dp)
  DO I = 1, RGRID%N
    OP(I,I) = OP(I,I) - LAMBDA_L/RGRID%R(I)**2
  END DO
END SUBROUTINE ADD_CURVATURE_TERM

END MODULE RADIAL_OPERATORS
