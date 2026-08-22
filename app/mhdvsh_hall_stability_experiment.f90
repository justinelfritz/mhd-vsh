!> Stage 2 (per .claude/plans/i-would-like-to-quizzical-comet.md):
!> standalone Hall-only stability experiment. Links mhdvsh_core only --
!> no diffusion_regime.f90, confirming this is decoupled from the
!> implicit-diffusion machinery, per the plan.
!>
!> Per user direction (2026-08-20): single seed mode Phi(k=1,l=0)
!> (degree 1, axisymmetric), Psi identically zero, RK4, all 4 term-groups
!> of HALL_INDUCTION_RHS (no simplification), evolved up to LMAX=30 so
!> self-coupling of the poloidal mode can gradually populate higher
!> degrees (an axisymmetric l-cascade -- the m=l+l' selection rule keeps
!> every populated mode at order 0 for all time, since the seed and Psi
!> are both order-0/zero).
!>
!> Boundary condition: simplified homogeneous Dirichlet (Phi=Psi=0 at
!> both ends, enforced after every RK4 step) -- per the plan, the real
!> vacuum Robin BC has no existing enforcement machinery for explicit
!> (non-matrix) data; that's deferred to the eventual combined regime
!> (Stage 4). The IC's sin() profile already vanishes at both ends, so
!> it satisfies this BC from t=0.
!>
!> DT/N_STEPS/F_HALL below are a first, deliberately conservative guess
!> (Stage 2's "start small and escalate"), not a validated stable
!> configuration -- watch the printed energy/max-active-l trace and the
!> NaN check; adjust and rerun rather than trusting one pass.
PROGRAM MHDVSH_HALL_STABILITY_EXPERIMENT
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
USE KINDS,             ONLY: dp, i4
USE VSH,               ONLY: YLM_INDEX
USE GRID_RADIAL,       ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,  ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,       ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS, ONLY: TOTAL_MAGNETIC_ENERGY, &
                              POLOIDAL_MAGNETIC_ENERGY_BY_L, TOROIDAL_MAGNETIC_ENERGY_BY_L
USE HALL_INDUCTION,    ONLY: HALL_INDUCTION_RHS
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 40
REAL(KIND=dp),    PARAMETER :: R_MIN   = 0.5_dp   ! km
REAL(KIND=dp),    PARAMETER :: R_MAX   = 1.0_dp   ! km
INTEGER(KIND=i4), PARAMETER :: LMAX    = 30
REAL(KIND=dp),    PARAMETER :: F_HALL  = 0.01_dp  ! km**2/(1e12 G)/yr, "weak"
REAL(KIND=dp),    PARAMETER :: DT      = 1.0E-3_dp ! yr -- starting guess, see header
INTEGER(KIND=i4), PARAMETER :: N_STEPS = 2000
INTEGER(KIND=i4), PARAMETER :: LOG_EVERY = 20
REAL(KIND=dp),    PARAMETER :: PI = 3.14159265358979_dp
REAL(KIND=dp),    PARAMETER :: BLOWUP_FACTOR = 10.0_dp
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 22

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI
TYPE(SPECTRAL_SCALAR_T) :: K1_PHI, K1_PSI, K2_PHI, K2_PSI, K3_PHI, K3_PSI, K4_PHI, K4_PSI
TYPE(SPECTRAL_SCALAR_T) :: TMP_PHI, TMP_PSI
REAL(KIND=dp) :: E0, E_NOW, T
REAL(KIND=dp), ALLOCATABLE :: E_POL_OF_L(:), E_TOR_OF_L(:)
INTEGER(KIND=i4) :: IR, ISTEP, MAX_ACTIVE_L, L
INTEGER(KIND=i4) :: IR_MID, IDX_L0
REAL(KIND=dp) :: MAX_IMAG_FRACTION
CHARACTER(LEN=1024) :: DATA_PATH, ENERGY_PATH
CHARACTER(LEN=16) :: COL_LABEL
INTEGER(KIND=i4), PARAMETER :: ENERGY_UNIT = 23

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(PSI, N_R, LMAX)

! Single seed mode: Phi(k=1,l=0), sin() profile (vanishes at both ends,
! matching the Dirichlet BC from t=0). Psi stays exactly zero (already
! zero-initialized by ALLOC_SPECTRAL_SCALAR) -- never set.
DO IR = 1, N_R
  PHI%COEF(IR, YLM_INDEX(1_i4,0_i4)) = &
    CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
END DO

E0 = TOTAL_MAGNETIC_ENERGY(PHI, PSI, OPS, RGRID)
WRITE(*,'(A)') 'MHD-VSH Hall stability experiment (Stage 2)'
WRITE(*,'(A,I0)')     '  N_r    = ', N_R
WRITE(*,'(A,I0)')     '  Lmax   = ', LMAX
WRITE(*,'(A,ES12.5)') '  F_Hall = ', F_HALL
WRITE(*,'(A,ES12.5,A,I0,A)') '  dt     = ', DT, '  (', N_STEPS, ' steps)'
WRITE(*,'(A,ES14.6)') '  E(t=0) = ', E0
WRITE(*,'(A)') ''
WRITE(*,'(A)') '  step        t       E/E0        max_active_l'

! Per-mode Phi/Psi(l,m=0) time series: Re(Phi_l0)/Re(Psi_l0) sampled at
! the domain midpoint (not the boundary -- that's pinned to 0 by the
! Dirichlet BC, which would make every mode read exactly zero there)
! for every l=1..LMAX, one column per l per field, plus a diagnostic
! column tracking how far either field has drifted from purely real
! (both must stay exactly real for a physical field restricted to m=0
! -- Phi_{l,-m}=conj(Phi_{l,m}) with m=0 forces Phi_{l,0}/Psi_{l,0}
! real; nonzero imaginary part is a numerical/physical-consistency
! check, not something the y=Phi/Psi plots themselves show).
IR_MID = NINT(REAL(N_R+1,KIND=dp)/2.0_dp)
DATA_PATH = 'hall_mode_amplitudes.dat'
IF (COMMAND_ARGUMENT_COUNT() >= 1) CALL GET_COMMAND_ARGUMENT(1, DATA_PATH)
OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
WRITE(DATA_UNIT,'(A)', ADVANCE='NO') '# step  t'
DO L = 1, LMAX
  WRITE(COL_LABEL,'(A,I0)') '  Phi_l', L
  WRITE(DATA_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
END DO
DO L = 1, LMAX
  WRITE(COL_LABEL,'(A,I0)') '  Psi_l', L
  WRITE(DATA_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
END DO
WRITE(DATA_UNIT,'(A)') '  max_imag_fraction'

! Per-degree magnetic energy (poloidal, toroidal, both radially
! integrated -- FIELD_DIAGNOSTICS::POLOIDAL_MAGNETIC_ENERGY_BY_L/
! TOROIDAL_MAGNETIC_ENERGY_BY_L, the mhd-vsh-relations.pdf E_B,pol(t)/
! E_B,tor(t) formulas broken out by l instead of summed over it).
ENERGY_PATH = 'hall_energy_by_l.dat'
IF (COMMAND_ARGUMENT_COUNT() >= 2) CALL GET_COMMAND_ARGUMENT(2, ENERGY_PATH)
OPEN(UNIT=ENERGY_UNIT, FILE=TRIM(ENERGY_PATH), STATUS='REPLACE', ACTION='WRITE')
WRITE(ENERGY_UNIT,'(A)', ADVANCE='NO') '# step  t'
DO L = 1, LMAX
  WRITE(COL_LABEL,'(A,I0)') '  Epol_l', L
  WRITE(ENERGY_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
END DO
DO L = 1, LMAX
  WRITE(COL_LABEL,'(A,I0)') '  Etor_l', L
  WRITE(ENERGY_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
END DO
WRITE(ENERGY_UNIT,'(A)') ''

ALLOCATE(E_POL_OF_L(0:LMAX), E_TOR_OF_L(0:LMAX))
T = 0.0_dp
E_POL_OF_L = POLOIDAL_MAGNETIC_ENERGY_BY_L(PHI, OPS, RGRID)
E_TOR_OF_L = TOROIDAL_MAGNETIC_ENERGY_BY_L(PSI, RGRID)
CALL LOG_MODE_AMPLITUDES(0_i4, T)
CALL LOG_ENERGY_BY_L(0_i4, T)
DO ISTEP = 1, N_STEPS
  CALL HALL_INDUCTION_RHS(PHI, PSI, OPS, RGRID, F_HALL, K1_PHI, K1_PSI)

  CALL ALLOC_SPECTRAL_SCALAR(TMP_PHI, N_R, LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(TMP_PSI, N_R, LMAX)
  TMP_PHI%COEF = PHI%COEF + 0.5_dp*DT*K1_PHI%COEF
  TMP_PSI%COEF = PSI%COEF + 0.5_dp*DT*K1_PSI%COEF
  CALL ZERO_BOUNDARY(TMP_PHI); CALL ZERO_BOUNDARY(TMP_PSI)
  CALL HALL_INDUCTION_RHS(TMP_PHI, TMP_PSI, OPS, RGRID, F_HALL, K2_PHI, K2_PSI)

  TMP_PHI%COEF = PHI%COEF + 0.5_dp*DT*K2_PHI%COEF
  TMP_PSI%COEF = PSI%COEF + 0.5_dp*DT*K2_PSI%COEF
  CALL ZERO_BOUNDARY(TMP_PHI); CALL ZERO_BOUNDARY(TMP_PSI)
  CALL HALL_INDUCTION_RHS(TMP_PHI, TMP_PSI, OPS, RGRID, F_HALL, K3_PHI, K3_PSI)

  TMP_PHI%COEF = PHI%COEF + DT*K3_PHI%COEF
  TMP_PSI%COEF = PSI%COEF + DT*K3_PSI%COEF
  CALL ZERO_BOUNDARY(TMP_PHI); CALL ZERO_BOUNDARY(TMP_PSI)
  CALL HALL_INDUCTION_RHS(TMP_PHI, TMP_PSI, OPS, RGRID, F_HALL, K4_PHI, K4_PSI)

  PHI%COEF = PHI%COEF + (DT/6.0_dp)*(K1_PHI%COEF + 2.0_dp*K2_PHI%COEF + 2.0_dp*K3_PHI%COEF + K4_PHI%COEF)
  PSI%COEF = PSI%COEF + (DT/6.0_dp)*(K1_PSI%COEF + 2.0_dp*K2_PSI%COEF + 2.0_dp*K3_PSI%COEF + K4_PSI%COEF)
  CALL ZERO_BOUNDARY(PHI); CALL ZERO_BOUNDARY(PSI)
  T = T + DT

  IF (MOD(ISTEP, LOG_EVERY) == 0 .OR. ISTEP == N_STEPS) THEN
    E_NOW = TOTAL_MAGNETIC_ENERGY(PHI, PSI, OPS, RGRID)
    E_POL_OF_L = POLOIDAL_MAGNETIC_ENERGY_BY_L(PHI, OPS, RGRID)
    E_TOR_OF_L = TOROIDAL_MAGNETIC_ENERGY_BY_L(PSI, RGRID)
    MAX_ACTIVE_L = 0
    DO L = 0, LMAX
      IF (E_POL_OF_L(L)+E_TOR_OF_L(L) > 1.0E-10_dp*MAXVAL(E_POL_OF_L+E_TOR_OF_L)) MAX_ACTIVE_L = L
    END DO
    WRITE(*,'(I8,ES12.4,ES12.4,I8)') ISTEP, T, E_NOW/E0, MAX_ACTIVE_L
    CALL LOG_MODE_AMPLITUDES(ISTEP, T)
    CALL LOG_ENERGY_BY_L(ISTEP, T)

    IF (IEEE_IS_NAN(E_NOW)) THEN
      WRITE(*,'(A)') '  BLOWUP: energy is NaN -- stopping.'
      STOP 1
    ELSE IF (E_NOW > BLOWUP_FACTOR*E0) THEN
      WRITE(*,'(A,F6.2,A)') '  BLOWUP: energy exceeded ', BLOWUP_FACTOR, 'x initial -- stopping.'
      STOP 1
    END IF
  END IF
END DO

CLOSE(DATA_UNIT)
CLOSE(ENERGY_UNIT)
WRITE(*,'(A)') ''
WRITE(*,'(A,A)') '  per-mode amplitudes written to ', TRIM(DATA_PATH)
WRITE(*,'(A,A)') '  per-degree energy written to    ', TRIM(ENERGY_PATH)
WRITE(*,'(A)') 'RESULT: completed without blowup at this (dt, Lmax, F_Hall).'

CONTAINS

!> Writes one row to DATA_UNIT: step, t, Re(Phi_l0) then Re(Psi_l0) at
!> the domain midpoint for every l=1..LMAX, and the largest
!> |Im(.)|/|.| seen across all l and both fields (should stay at
!> round-off for a physically real, purely-axisymmetric field -- see
!> the header comment above where DATA_UNIT is opened).
SUBROUTINE LOG_MODE_AMPLITUDES(ISTEP_ARG, T_ARG)
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP_ARG
  REAL(KIND=dp),    INTENT(IN) :: T_ARG
  COMPLEX(KIND=dp) :: VAL
  INTEGER(KIND=i4) :: LL
  WRITE(DATA_UNIT,'(I8,ES14.6)', ADVANCE='NO') ISTEP_ARG, T_ARG
  MAX_IMAG_FRACTION = 0.0_dp
  DO LL = 1, LMAX
    IDX_L0 = YLM_INDEX(LL, 0_i4)
    VAL = PHI%COEF(IR_MID, IDX_L0)
    WRITE(DATA_UNIT,'(ES14.6)', ADVANCE='NO') REAL(VAL, KIND=dp)
    IF (ABS(VAL) > 1.0E-13_dp) THEN
      MAX_IMAG_FRACTION = MAX(MAX_IMAG_FRACTION, ABS(AIMAG(VAL))/ABS(VAL))
    END IF
  END DO
  DO LL = 1, LMAX
    IDX_L0 = YLM_INDEX(LL, 0_i4)
    VAL = PSI%COEF(IR_MID, IDX_L0)
    WRITE(DATA_UNIT,'(ES14.6)', ADVANCE='NO') REAL(VAL, KIND=dp)
    IF (ABS(VAL) > 1.0E-13_dp) THEN
      MAX_IMAG_FRACTION = MAX(MAX_IMAG_FRACTION, ABS(AIMAG(VAL))/ABS(VAL))
    END IF
  END DO
  WRITE(DATA_UNIT,'(ES14.6)') MAX_IMAG_FRACTION
END SUBROUTINE LOG_MODE_AMPLITUDES

!> Zeros the outer/inner radial rows (Dirichlet BC) of every mode's
!> coefficient -- HALL_INDUCTION_RHS applies no BC itself (see its own
!> docstring), this is the caller's job, per the module's contract.
SUBROUTINE ZERO_BOUNDARY(FIELD)
  TYPE(SPECTRAL_SCALAR_T), INTENT(INOUT) :: FIELD
  FIELD%COEF(1,:)   = (0.0_dp, 0.0_dp)
  FIELD%COEF(N_R,:) = (0.0_dp, 0.0_dp)
END SUBROUTINE ZERO_BOUNDARY

!> Writes one row to ENERGY_UNIT: step, t, E_POL_OF_L(1:LMAX) then
!> E_TOR_OF_L(1:LMAX) -- both already computed this call by
!> FIELD_DIAGNOSTICS::POLOIDAL_MAGNETIC_ENERGY_BY_L/
!> TOROIDAL_MAGNETIC_ENERGY_BY_L in the main loop above.
SUBROUTINE LOG_ENERGY_BY_L(ISTEP_ARG, T_ARG)
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP_ARG
  REAL(KIND=dp),    INTENT(IN) :: T_ARG
  INTEGER(KIND=i4) :: LL
  WRITE(ENERGY_UNIT,'(I8,ES14.6)', ADVANCE='NO') ISTEP_ARG, T_ARG
  DO LL = 1, LMAX
    WRITE(ENERGY_UNIT,'(ES14.6)', ADVANCE='NO') E_POL_OF_L(LL)
  END DO
  DO LL = 1, LMAX
    WRITE(ENERGY_UNIT,'(ES14.6)', ADVANCE='NO') E_TOR_OF_L(LL)
  END DO
  WRITE(ENERGY_UNIT,'(A)') ''
END SUBROUTINE LOG_ENERGY_BY_L

END PROGRAM MHDVSH_HALL_STABILITY_EXPERIMENT
