!> Checks FIELD_DIAGNOSTICS' MAGNETIC_ENERGY_DENSITY/TOTAL_MAGNETIC_ENERGY
!> implement the formula documented in that module's header correctly --
!> NOT a check that the underlying physics coefficients are the final,
!> correct ones (they're explicitly flagged there as provisional pending
!> re-derivation). Two hand-picked modes with known closed-form profiles
!> give an exactly computable expected density at every grid point, and
!> an independently-coded trapezoidal sum gives the expected total, so
!> nothing here has to tolerate discretization error against a continuum
!> answer -- only against another discrete evaluation of the same
!> formula.
PROGRAM TEST_FIELD_DIAGNOSTICS
USE KINDS,               ONLY: dp, i4
USE VSH,                 ONLY: YLM_INDEX
USE GRID_RADIAL,         ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,    ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,         ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,   ONLY: MAGNETIC_ENERGY_DENSITY, TOTAL_MAGNETIC_ENERGY
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R   = 30
REAL(KIND=dp),    PARAMETER :: R_MIN = 0.5_dp
REAL(KIND=dp),    PARAMETER :: R_MAX = 1.0_dp
INTEGER(KIND=i4), PARAMETER :: LMAX  = 4
REAL(KIND=dp),    PARAMETER :: TOL   = 1.0E-10_dp

! Mode 1: L=2,M=1, Phi(r)=r (Phi'=1 exactly), Psi=const complex.
INTEGER(KIND=i4), PARAMETER :: L1 = 2, M1 = 1
COMPLEX(KIND=dp), PARAMETER :: PSI1 = CMPLX(1.2_dp, -0.5_dp, KIND=dp)
! Mode 2: L=3,M=-2, Phi=const complex (Phi'=0 exactly), Psi=0.
INTEGER(KIND=i4), PARAMETER :: L2 = 3, M2 = -2
COMPLEX(KIND=dp), PARAMETER :: PHI2 = CMPLX(0.8_dp, 0.6_dp, KIND=dp)

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI
REAL(KIND=dp), ALLOCATABLE :: EXPECTED_E_OF_R(:), GOT_E_OF_R(:)
REAL(KIND=dp) :: LAMBDA1, LAMBDA2, EXPECTED_E_TOTAL, GOT_E_TOTAL
REAL(KIND=dp) :: MAX_DENSITY_ERR, TOTAL_ERR
INTEGER(KIND=i4) :: IR, N_FAIL

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(PSI, N_R, LMAX)

DO IR = 1, N_R
  PHI%COEF(IR, YLM_INDEX(L1,M1)) = CMPLX(RGRID%R(IR), 0.0_dp, KIND=dp)
END DO
PSI%COEF(:, YLM_INDEX(L1,M1)) = PSI1
PHI%COEF(:, YLM_INDEX(L2,M2)) = PHI2

LAMBDA1 = REAL(L1*(L1+1), KIND=dp)
LAMBDA2 = REAL(L2*(L2+1), KIND=dp)

ALLOCATE(EXPECTED_E_OF_R(N_R))
DO IR = 1, N_R
  EXPECTED_E_OF_R(IR) = &
    LAMBDA1*(LAMBDA1*RGRID%R(IR)**2/RGRID%R(IR)**2 + 1.0_dp + ABS(PSI1)**2) + &
    LAMBDA2*(LAMBDA2*ABS(PHI2)**2/RGRID%R(IR)**2)
END DO

EXPECTED_E_TOTAL = 0.0_dp
DO IR = 1, N_R-1
  EXPECTED_E_TOTAL = EXPECTED_E_TOTAL + &
    0.5_dp*(EXPECTED_E_OF_R(IR)+EXPECTED_E_OF_R(IR+1))*(RGRID%R(IR+1)-RGRID%R(IR))
END DO
EXPECTED_E_TOTAL = 0.5_dp * EXPECTED_E_TOTAL

GOT_E_OF_R = MAGNETIC_ENERGY_DENSITY(PHI, PSI, OPS, RGRID)
GOT_E_TOTAL = TOTAL_MAGNETIC_ENERGY(PHI, PSI, OPS, RGRID)

N_FAIL = 0
MAX_DENSITY_ERR = MAXVAL(ABS(GOT_E_OF_R-EXPECTED_E_OF_R))
IF (MAX_DENSITY_ERR > TOL) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  energy_density  max_err=", MAX_DENSITY_ERR
ELSE
  WRITE(*,'(A,ES10.3)') "PASS  energy_density  max_err=", MAX_DENSITY_ERR
END IF

TOTAL_ERR = ABS(GOT_E_TOTAL-EXPECTED_E_TOTAL)
IF (TOTAL_ERR > TOL) THEN
  N_FAIL = N_FAIL + 1
  WRITE(*,'(A,ES10.3)') "FAIL  total_energy  err=", TOTAL_ERR
ELSE
  WRITE(*,'(A,ES10.3)') "PASS  total_energy  err=", TOTAL_ERR
END IF

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_field_diagnostics"
END IF

END PROGRAM TEST_FIELD_DIAGNOSTICS
