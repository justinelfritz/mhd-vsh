!> One-off diagnostic (2026-08-29, per user request): dumps the seed IC's
!> Phi(r) profile (exactly as constructed by app/mhdvsh_hall_crust_profile.f90,
!> same TOV solve/grid/IC code, not a reimplementation) across the full
!> radial domain, plus the derived B_r(r)=l(l+1)/r**2 * Phi(r) and the
!> FD-operator (OPS%D1, the actual discretized derivative the code uses,
!> not just the continuous analytic one) evaluation of the Robin-BC
!> residual at the outer boundary -- d(Phi)/dr|_Rout + (l/Rout)*Phi(Rout).
!> Built to let the user visually inspect the t=0 field near the surface
!> directly, following up on the analytic Robin-BC-violation finding
!> from earlier in this investigation (ROADMAP.md analytical task 5).
!>
!> Usage: mhdvsh_ic_boundary_check <output_data_path>
PROGRAM MHDVSH_IC_BOUNDARY_CHECK
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE TOV_SOLVER,         ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
USE CRUST_CONDUCTIVITY, ONLY: FIND_TRUNCATION_RADIUS
IMPLICIT NONE

REAL(KIND=dp),    PARAMETER :: NBAR_CENTRAL = 0.5447307_dp
CHARACTER(LEN=*), PARAMETER :: EOS_PATH = 'data/eos/APR_EOS_Cat.dat'
INTEGER(KIND=i4), PARAMETER :: NPOINTS_TOV = 414
REAL(KIND=dp),    PARAMETER :: RHOL_CGS = 2.2E14_dp
REAL(KIND=dp),    PARAMETER :: THETA_MAX_TRUNCATE = 1.3_dp
REAL(KIND=dp),    PARAMETER :: T_KELVIN = 1.0E9_dp
INTEGER(KIND=i4), PARAMETER :: N_R = 40
INTEGER(KIND=i4), PARAMETER :: LMAX = 30
REAL(KIND=dp),    PARAMETER :: PI = 3.14159265358979_dp
INTEGER(KIND=i4), PARAMETER :: L_SEED = 1_i4, M_SEED = 0_i4

TYPE(TOV_PROFILE_T)     :: PROFILE
TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(SPECTRAL_SCALAR_T) :: PHI
REAL(KIND=dp) :: R_MIN, R_MAX, LAMBDA_L, DPHI_DR_BND, ROBIN_LHS, ROBIN_RESIDUAL
REAL(KIND=dp) :: B_R
INTEGER(KIND=i4) :: IR, IDX, UNIT
CHARACTER(LEN=1024) :: OUT_PATH

IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
  WRITE(*,'(A)') 'Usage: mhdvsh_ic_boundary_check <output_data_path>'
  STOP 1
END IF
CALL GET_COMMAND_ARGUMENT(1, OUT_PATH)

CALL SOLVE_TOV_STAR(NBAR_CENTRAL, EOS_PATH, NPOINTS_TOV, PROFILE)
R_MIN = MINVAL(PROFILE%R, MASK=PROFILE%RHOCGS <= RHOL_CGS)
CALL FIND_TRUNCATION_RADIUS(PROFILE, T_KELVIN, R_MAX, THETA_MAX=THETA_MAX_TRUNCATE)

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)

CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX)
! Exact same seed IC as app/mhdvsh_hall_crust_profile.f90's own
! construction (copied verbatim, not reimplemented).
DO IR = 1, N_R
  PHI%COEF(IR, YLM_INDEX(L_SEED,M_SEED)) = &
    CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
END DO

IDX = YLM_INDEX(L_SEED, M_SEED)
LAMBDA_L = REAL(L_SEED*(L_SEED+1), KIND=dp)

OPEN(NEWUNIT=UNIT, FILE=TRIM(OUT_PATH), STATUS='REPLACE', ACTION='WRITE')
WRITE(UNIT,'(A)') '# ir  r_km  Phi_l1m0(r)  B_r(r)=l(l+1)/r**2*Phi(r)'
DO IR = 1, N_R
  B_R = LAMBDA_L / RGRID%R(IR)**2 * REAL(PHI%COEF(IR,IDX), KIND=dp)
  WRITE(UNIT,'(I8,3ES18.10)') IR, RGRID%R(IR), REAL(PHI%COEF(IR,IDX),KIND=dp), B_R
END DO
CLOSE(UNIT)

! FD-operator (OPS%D1, the ACTUAL discretized derivative used by the
! code -- not the continuous analytic derivative) evaluation of the
! Robin residual at the outer boundary: d(Phi)/dr|_Rout + (l/Rout)*Phi(Rout).
! Should be exactly zero for a BC-consistent state.
DPHI_DR_BND = SUM(OPS%D1(N_R,:) * REAL(PHI%COEF(:,IDX), KIND=dp))
ROBIN_LHS = DPHI_DR_BND + (REAL(L_SEED,KIND=dp)/RGRID%R(N_R)) * REAL(PHI%COEF(N_R,IDX),KIND=dp)
ROBIN_RESIDUAL = ROBIN_LHS

WRITE(*,'(A)')            'IC boundary check (l=1, m=0 seed mode):'
WRITE(*,'(A,ES16.8,A)')   '  R_MIN = ', R_MIN, ' km'
WRITE(*,'(A,ES16.8,A)')   '  R_MAX = ', R_MAX, ' km'
WRITE(*,'(A,ES16.8)')     '  Phi(R_MIN)          = ', REAL(PHI%COEF(1,IDX),KIND=dp)
WRITE(*,'(A,ES16.8)')     '  Phi(R_MAX)          = ', REAL(PHI%COEF(N_R,IDX),KIND=dp)
WRITE(*,'(A,ES16.8,A)')   '  dPhi/dr|_Rout (FD, OPS%D1) = ', DPHI_DR_BND, ' km**-1'
WRITE(*,'(A,ES16.8,A)')   '  required dPhi/dr|_Rout (Robin, = -(l/Rout)*Phi(Rout)) = ', &
  -(REAL(L_SEED,KIND=dp)/RGRID%R(N_R))*REAL(PHI%COEF(N_R,IDX),KIND=dp), ' km**-1'
WRITE(*,'(A,ES16.8)')     '  Robin residual (should be 0 if BC-consistent) = ', ROBIN_RESIDUAL
WRITE(*,'(A,A)')          '  full Phi(r) profile written to ', TRIM(OUT_PATH)

END PROGRAM MHDVSH_IC_BOUNDARY_CHECK
