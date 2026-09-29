!> One-off verification harness (2026-08-29, per user request "understand
!> why energy conservation is being violated"). Not a permanent feature.
!>
!> Context: results/hall_crust_profile_1000yr/'s production run shows a
!> growing energy-balance residual Edot_tot that tracks Edot_S_Hall at a
!> nearly CONSTANT ratio (~-0.94) from t~150yr onward. A parallel LMAX=30
!> vs LMAX=45 rerun (results/hall_crust_profile_1000yr_lmax45/) came back
!> BIT-IDENTICAL at every logged step, ruling out spectral truncation of
!> the Hall mode-coupling term as the cause. That leaves a live hypothesis
!> that FIELD_DIAGNOSTICS::HALL_POYNTING_FLUX_RATE itself is not the
!> correct boundary-flux formula for a realistic (multi-mode-coupled)
!> field state -- e.g. a missing J-type Gaunt-coefficient term (the Hall
!> induction equations have both I- and J-coupling; the Poynting flux
!> formula in analytic_formulas/mhd-vsh-relations.tex, eq. "Energy Budget
!> in Hall Limit", has only I-coupling).
!>
!> This harness re-runs this project's OWN original validation method for
!> that exact formula (Stage 1 of .claude/plans/i-would-like-to-quizzical-
!> comet.md: step HALL_INDUCTION_RHS alone, no diffusion, at shrinking
!> dt_hall, and check the finite-difference dE/dt converges to
!> HALL_POYNTING_FLUX_RATE's prediction) -- but starting from a REAL,
!> multi-mode-coupled state (the production run's own t=1000yr checkpoint)
!> instead of the trivial single-seed-mode state Stage 1 originally used.
!> UPDATE (2026-08-29): the first run of this harness (hard-zeroing the
!> boundary during substeps, matching the then-current HALL_SUBSTEPS)
!> found the answer directly -- pure-Hall substepping grew total energy
!> 11.7x in one outer step (converged across a 128x dt_hall sweep, so a
!> real effect, not a discretization artifact), traced to the checkpoint's
!> large, physically real outer-boundary Phi value (60% of the field's
!> own global max) being hard-zeroed every RK4 stage right where F_HALL
!> is ~1e7x its inner-boundary value. hall_regime.f90's HALL_SUBSTEPS now
!> holds the boundary fixed instead of zeroing it; this harness was
!> updated to match, and is kept as a regression check that the fix
!> actually restores convergence to HALL_POYNTING_FLUX_RATE's prediction.
!>
!> Usage: mhdvsh_hall_convergence_check <checkpoint_path>
PROGRAM MHDVSH_HALL_CONVERGENCE_CHECK
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_POLOIDAL_MAGNETIC_ENERGY, TOTAL_TOROIDAL_MAGNETIC_ENERGY, &
                               HALL_POYNTING_FLUX_RATE
USE HALL_INDUCTION,     ONLY: HALL_INDUCTION_RHS
USE IO_CHECKPOINT,      ONLY: READ_CHECKPOINT
USE TOV_SOLVER,         ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
USE CRUST_CONDUCTIVITY, ONLY: ETA_AND_F_HALL_ON_GRID, FIND_TRUNCATION_RADIUS
USE UNITS,              ONLY: ENERGY_UNIT_ERG, POWER_UNIT_ERG_PER_S
IMPLICIT NONE

REAL(KIND=dp),    PARAMETER :: NBAR_CENTRAL = 0.5447307_dp
CHARACTER(LEN=*), PARAMETER :: EOS_PATH = 'data/eos/APR_EOS_Cat.dat'
INTEGER(KIND=i4), PARAMETER :: NPOINTS_TOV = 414
REAL(KIND=dp),    PARAMETER :: RHOL_CGS = 2.2E14_dp
REAL(KIND=dp),    PARAMETER :: THETA_MAX_TRUNCATE = 1.3_dp
REAL(KIND=dp),    PARAMETER :: T_KELVIN = 1.0E9_dp
REAL(KIND=dp),    PARAMETER :: DT_TOTAL = 0.192827882_dp  ! matches production run's own dt_yr at t=1000yr

TYPE(TOV_PROFILE_T)     :: PROFILE
TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(SPECTRAL_SCALAR_T) :: PHI0, PSI0, PHI, PSI
TYPE(SPECTRAL_SCALAR_T) :: K1_PHI, K1_PSI, K2_PHI, K2_PSI, K3_PHI, K3_PSI, K4_PHI, K4_PSI
TYPE(SPECTRAL_SCALAR_T) :: TMP_PHI, TMP_PSI
COMPLEX(KIND=dp), ALLOCATABLE :: PHI_IN(:), PHI_OUT(:), PSI_IN(:), PSI_OUT(:)
REAL(KIND=dp), ALLOCATABLE :: ETA_PROFILE(:), F_HALL_PROFILE(:), N_E_PROFILE(:)
REAL(KIND=dp) :: R_MIN, R_MAX, F_HALL_MID
REAL(KIND=dp) :: E0, E1, DT_HALL, EDOT_FD, EDOT_S_HALL_PRED
INTEGER(KIND=i4) :: N_R, LMAX, ISTEP, T_ISTEP, ISUB, N_SUB_TEST
REAL(KIND=dp) :: T_CK, R_MIN_CK, R_MAX_CK
CHARACTER(LEN=1024) :: CHECKPOINT_PATH

IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
  WRITE(*,'(A)') 'Usage: mhdvsh_hall_convergence_check <checkpoint_path>'
  STOP 1
END IF
CALL GET_COMMAND_ARGUMENT(1, CHECKPOINT_PATH)

CALL READ_CHECKPOINT(TRIM(CHECKPOINT_PATH), T_ISTEP, T_CK, PHI0, PSI0, N_R, LMAX, R_MIN_CK, R_MAX_CK)
WRITE(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') '  checkpoint: N_R=', N_R, ' LMAX=', LMAX, &
  ' t=', T_CK, ' R_MIN=', R_MIN_CK, ' R_MAX=', R_MAX_CK

! Rebuild the SAME TOV profile / eta(r)/f_H(r) the production run used --
! deterministic from (NBAR_CENTRAL, EOS_PATH, T_KELVIN), reproduces the
! checkpoint's own R_MIN/R_MAX exactly (checked below).
CALL SOLVE_TOV_STAR(NBAR_CENTRAL, EOS_PATH, NPOINTS_TOV, PROFILE)
R_MIN = MINVAL(PROFILE%R, MASK=PROFILE%RHOCGS <= RHOL_CGS)
CALL FIND_TRUNCATION_RADIUS(PROFILE, T_KELVIN, R_MAX, THETA_MAX=THETA_MAX_TRUNCATE)
IF (ABS(R_MIN-R_MIN_CK) > 1.0E-6_dp .OR. ABS(R_MAX-R_MAX_CK) > 1.0E-6_dp) THEN
  WRITE(*,'(A)') 'MHDVSH_HALL_CONVERGENCE_CHECK: rebuilt grid does not match checkpoint -- aborting'
  STOP 1
END IF

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)
CALL ETA_AND_F_HALL_ON_GRID(PROFILE, RGRID%R, ETA_PROFILE, F_HALL_PROFILE, N_E_PROFILE, T_KELVIN=T_KELVIN)
F_HALL_MID = F_HALL_PROFILE(N_R/2)

E0 = TOTAL_POLOIDAL_MAGNETIC_ENERGY(PHI0, OPS, RGRID) + TOTAL_TOROIDAL_MAGNETIC_ENERGY(PSI0, RGRID)
EDOT_S_HALL_PRED = HALL_POYNTING_FLUX_RATE(PHI0, PSI0, OPS, RGRID, F_HALL_MID, F_HALL_PROFILE=F_HALL_PROFILE)

! Boundary-value / F_HALL magnitude sanity check: does zeroing the
! boundary during Hall substeps discard a physically significant Phi/Psi
! value right where F_HALL is most extreme (near-surface, N_R)?
BLOCK
  INTEGER(KIND=i4) :: LL, IDX
  DO LL = 0, MIN(3_i4, LMAX)
    IDX = YLM_INDEX(LL, 0_i4)
    WRITE(*,'(A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') '  l=', LL, &
      '  PHI(r_in)=', ABS(PHI0%COEF(1,IDX)), '  PHI(r_out)=', ABS(PHI0%COEF(N_R,IDX)), &
      '  PSI(r_in)=', ABS(PSI0%COEF(1,IDX)), '  PSI(r_out)=', ABS(PSI0%COEF(N_R,IDX))
  END DO
  WRITE(*,'(A,ES14.6,A,ES14.6)') '  F_HALL_PROFILE(1)=', F_HALL_PROFILE(1), &
    '  F_HALL_PROFILE(N_R)=', F_HALL_PROFILE(N_R)
  WRITE(*,'(A,ES14.6)') '  max |PHI| over all modes at r_out (surface) = ', MAXVAL(ABS(PHI0%COEF(N_R,:)))
  WRITE(*,'(A,ES14.6)') '  max |PHI| over all modes, all r             = ', MAXVAL(ABS(PHI0%COEF))
END BLOCK
WRITE(*,'(A,ES16.8,A)') '  E(checkpoint) = ', E0*ENERGY_UNIT_ERG, ' erg'
WRITE(*,'(A,ES16.8,A)') '  HALL_POYNTING_FLUX_RATE prediction (analytic, at checkpoint state) = ', &
  EDOT_S_HALL_PRED*POWER_UNIT_ERG_PER_S, ' erg/s'
WRITE(*,'(A)') ''
WRITE(*,'(A)') '  N_SUB_TEST   dt_hall [yr]   finite-diff dE/dt [erg/s]   ' // &
  'ratio to analytic prediction'

! Hold-fixed placeholder BC during substeps, matching the 2026-08-29 fix
! in hall_regime.f90's HALL_SUBSTEPS (was a hard zero -- see this file's
! own header and hall_regime.f90's updated docstrings for why).
ALLOCATE(PHI_IN(PHI0%NLM), PHI_OUT(PHI0%NLM), PSI_IN(PSI0%NLM), PSI_OUT(PSI0%NLM))
PHI_IN = PHI0%COEF(1,:); PHI_OUT = PHI0%COEF(N_R,:)
PSI_IN = PSI0%COEF(1,:); PSI_OUT = PSI0%COEF(N_R,:)

DO N_SUB_TEST = 1, 8
  ISTEP = 2**(N_SUB_TEST-1)
  DT_HALL = DT_TOTAL / REAL(ISTEP, KIND=dp)

  CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(PSI, N_R, LMAX)
  PHI%COEF = PHI0%COEF
  PSI%COEF = PSI0%COEF

  DO ISUB = 1, ISTEP
    CALL HALL_INDUCTION_RHS(PHI, PSI, OPS, RGRID, F_HALL_MID, K1_PHI, K1_PSI, F_HALL_PROFILE=F_HALL_PROFILE)

    CALL ALLOC_SPECTRAL_SCALAR(TMP_PHI, N_R, LMAX)
    CALL ALLOC_SPECTRAL_SCALAR(TMP_PSI, N_R, LMAX)
    TMP_PHI%COEF = PHI%COEF + 0.5_dp*DT_HALL*K1_PHI%COEF
    TMP_PSI%COEF = PSI%COEF + 0.5_dp*DT_HALL*K1_PSI%COEF
    TMP_PHI%COEF(1,:) = PHI_IN; TMP_PHI%COEF(N_R,:) = PHI_OUT
    TMP_PSI%COEF(1,:) = PSI_IN; TMP_PSI%COEF(N_R,:) = PSI_OUT
    CALL HALL_INDUCTION_RHS(TMP_PHI, TMP_PSI, OPS, RGRID, F_HALL_MID, K2_PHI, K2_PSI, F_HALL_PROFILE=F_HALL_PROFILE)

    TMP_PHI%COEF = PHI%COEF + 0.5_dp*DT_HALL*K2_PHI%COEF
    TMP_PSI%COEF = PSI%COEF + 0.5_dp*DT_HALL*K2_PSI%COEF
    TMP_PHI%COEF(1,:) = PHI_IN; TMP_PHI%COEF(N_R,:) = PHI_OUT
    TMP_PSI%COEF(1,:) = PSI_IN; TMP_PSI%COEF(N_R,:) = PSI_OUT
    CALL HALL_INDUCTION_RHS(TMP_PHI, TMP_PSI, OPS, RGRID, F_HALL_MID, K3_PHI, K3_PSI, F_HALL_PROFILE=F_HALL_PROFILE)

    TMP_PHI%COEF = PHI%COEF + DT_HALL*K3_PHI%COEF
    TMP_PSI%COEF = PSI%COEF + DT_HALL*K3_PSI%COEF
    TMP_PHI%COEF(1,:) = PHI_IN; TMP_PHI%COEF(N_R,:) = PHI_OUT
    TMP_PSI%COEF(1,:) = PSI_IN; TMP_PSI%COEF(N_R,:) = PSI_OUT
    CALL HALL_INDUCTION_RHS(TMP_PHI, TMP_PSI, OPS, RGRID, F_HALL_MID, K4_PHI, K4_PSI, F_HALL_PROFILE=F_HALL_PROFILE)

    PHI%COEF = PHI%COEF + (DT_HALL/6.0_dp)*(K1_PHI%COEF + 2.0_dp*K2_PHI%COEF + 2.0_dp*K3_PHI%COEF + K4_PHI%COEF)
    PSI%COEF = PSI%COEF + (DT_HALL/6.0_dp)*(K1_PSI%COEF + 2.0_dp*K2_PSI%COEF + 2.0_dp*K3_PSI%COEF + K4_PSI%COEF)
    PHI%COEF(1,:) = PHI_IN; PHI%COEF(N_R,:) = PHI_OUT
    PSI%COEF(1,:) = PSI_IN; PSI%COEF(N_R,:) = PSI_OUT
  END DO

  E1 = TOTAL_POLOIDAL_MAGNETIC_ENERGY(PHI, OPS, RGRID) + TOTAL_TOROIDAL_MAGNETIC_ENERGY(PSI, RGRID)
  EDOT_FD = (E1-E0)/DT_TOTAL

  WRITE(*,'(I10,ES16.8,ES28.8,ES20.6,A,ES16.8,A,ES12.4)') ISTEP, DT_HALL, EDOT_FD*POWER_UNIT_ERG_PER_S, &
    (EDOT_FD/EDOT_S_HALL_PRED), '   E1=', E1*ENERGY_UNIT_ERG, ' erg   E1/E0=', (E1*ENERGY_UNIT_ERG)/(E0*ENERGY_UNIT_ERG)
END DO

END PROGRAM MHDVSH_HALL_CONVERGENCE_CHECK
