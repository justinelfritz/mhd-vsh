!> Checks HALL_REGIME end to end: after evolving, the boundary rows
!> actually solved for must satisfy the imposed vacuum conditions --
!> inner Dirichlet on both Phi and Psi, outer Dirichlet on Psi, outer
!> Robin on Phi (same checks test_diffusion_regime.f90 makes on
!> DIFFUSION_REGIME, now on the combined regime -- this is what confirms
!> HALL_REGIME's "the trailing DIFFUSION_ADVANCE call enforces the true
!> BC regardless of the Hall substeps" design argument actually holds in
!> code, not just in reasoning; see hall_regime.f90's own header) -- plus
!> a short run with no NaN/blowup, same idiom as
!> app/mhdvsh_hall_stability_experiment.f90.
PROGRAM TEST_HALL_REGIME
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_MAGNETIC_ENERGY
USE DIFFUSION_REGIME,   ONLY: DIFFUSION_STATE_T
USE HALL_REGIME,        ONLY: HALL_INIT, HALL_ADVANCE
USE TIMESTEPPER,        ONLY: RUN
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 25
REAL(KIND=dp),    PARAMETER :: R_MIN   = 0.5_dp
REAL(KIND=dp),    PARAMETER :: R_MAX   = 1.0_dp
INTEGER(KIND=i4), PARAMETER :: LMAX    = 6
REAL(KIND=dp),    PARAMETER :: ETA     = 0.05_dp
REAL(KIND=dp),    PARAMETER :: F_HALL  = 0.01_dp
REAL(KIND=dp),    PARAMETER :: DT      = 0.01_dp
INTEGER(KIND=i4), PARAMETER :: N_SUB   = 10
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 20
REAL(KIND=dp),    PARAMETER :: TOL     = 1.0E-10_dp
REAL(KIND=dp),    PARAMETER :: PI      = 3.14159265358979_dp

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(DIFFUSION_STATE_T) :: STATE
REAL(KIND=dp) :: E0, E1
INTEGER(KIND=i4) :: IR, N_FAIL

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL HALL_INIT(RGRID, OPS, LMAX, ETA, F_HALL, DT, N_SUB)

CALL ALLOC_SPECTRAL_SCALAR(STATE%PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(STATE%PSI, N_R, LMAX)

! Single seed mode, same as app/mhdvsh_hall_stability_experiment.f90 --
! sin() profile already vanishes at both ends.
DO IR = 1, N_R
  STATE%PHI%COEF(IR, YLM_INDEX(1_i4,0_i4)) = &
    CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
END DO

E0 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
CALL RUN(HALL_ADVANCE, STATE, DT, N_STEPS)
E1 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)

N_FAIL = 0
IF (IEEE_IS_NAN(E1)) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A)') "FAIL  no_blowup  E1=NaN"
ELSE
  WRITE(*,'(A,ES12.5,A,ES12.5)') "PASS  no_blowup  E0=", E0, "  E1=", E1
END IF

CALL CHECK_INNER_DIRICHLET(STATE%PHI, "phi", N_FAIL)
CALL CHECK_INNER_DIRICHLET(STATE%PSI, "psi", N_FAIL)
CALL CHECK_OUTER_DIRICHLET(STATE%PSI, N_FAIL)
CALL CHECK_OUTER_ROBIN(STATE%PHI, N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_hall_regime"
END IF

CONTAINS

SUBROUTINE CHECK_INNER_DIRICHLET(FIELD, NAME, N_FAIL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)    :: FIELD
  CHARACTER(*),             INTENT(IN)    :: NAME
  INTEGER(KIND=i4),         INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: MAX_ERR
  MAX_ERR = MAXVAL(ABS(FIELD%COEF(1,:)))
  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  inner_dirichlet_", NAME, "  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  inner_dirichlet_", NAME, "  max_err=", MAX_ERR
  END IF
END SUBROUTINE CHECK_INNER_DIRICHLET

SUBROUTINE CHECK_OUTER_DIRICHLET(PSI, N_FAIL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)    :: PSI
  INTEGER(KIND=i4),         INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: MAX_ERR
  MAX_ERR = MAXVAL(ABS(PSI%COEF(N_R,:)))
  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  outer_dirichlet_psi  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  outer_dirichlet_psi  max_err=", MAX_ERR
  END IF
END SUBROUTINE CHECK_OUTER_DIRICHLET

!> d(Phi)/dr + (l/r_out)*Phi should be ~0 at r_out for every populated
!> mode -- same check test_diffusion_regime.f90 makes, now confirming
!> the trailing DIFFUSION_ADVANCE call inside HALL_ADVANCE re-establishes
!> the true Robin BC regardless of the preceding explicit Hall substeps.
SUBROUTINE CHECK_OUTER_ROBIN(PHI, N_FAIL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN)    :: PHI
  INTEGER(KIND=i4),         INTENT(INOUT) :: N_FAIL
  COMPLEX(KIND=dp) :: RESIDUAL
  REAL(KIND=dp)    :: MAX_ERR
  INTEGER(KIND=i4) :: L, M, IDX

  MAX_ERR = 0.0_dp
  DO L = 0, LMAX
    DO M = -L, L
      IDX = YLM_INDEX(L,M)
      RESIDUAL = SUM(OPS%D1(N_R,:)*PHI%COEF(:,IDX)) + &
                 (REAL(L,KIND=dp)/RGRID%R(N_R))*PHI%COEF(N_R,IDX)
      MAX_ERR = MAX(MAX_ERR, ABS(RESIDUAL))
    END DO
  END DO
  IF (MAX_ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,ES10.3)') "FAIL  outer_robin_phi  max_err=", MAX_ERR
  ELSE
    WRITE(*,'(A,ES10.3)') "PASS  outer_robin_phi  max_err=", MAX_ERR
  END IF
END SUBROUTINE CHECK_OUTER_ROBIN

END PROGRAM TEST_HALL_REGIME
