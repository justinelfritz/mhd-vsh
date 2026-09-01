!> Combined resistive+Hall regime driver with a REAL, radially-varying
!> eta(r)/f_H(r) (not a uniform toy constant, and not one physically-
!> sourced constant like mhdvsh_hall_crust_sample.f90) -- solves the TOV
!> profile fresh (same M=1.40 reference star as app/mhdvsh_tov.f90),
!> derives the crust extent from it, samples eta(r)/f_H(r) onto the
!> simulation's own radial grid via CRUST_CONDUCTIVITY::
!> ETA_AND_F_HALL_ON_GRID (new, 2026-08-24), and wires the resulting
!> arrays into HALL_INIT's existing optional ETA_PROFILE/F_HALL_PROFILE
!> arguments. Built per user request (2026-08-24) to run a `tmax=1000`
!> yr test comparable to the existing uniform-toy-model archive
!> (results/hall_dt_stability_2026-08-21/).
!>
!> @warning DYNAMIC (Hall-CFL-limited) outer timestep, not a fixed DT --
!>   built on TIMESTEPPER::RUN_ADAPTIVE + HALL_REGIME::HALL_COMPUTE_DT/
!>   HALL_SET_DT, same machinery as app/mhdvsh_hall_adaptive.f90. This
!>   was a deliberate pivot away from this driver's own original fixed-
!>   dt design (per the plan addendum this driver's header used to
!>   reference): an empirical smoke test at dt_hall=0.001 (the value
!>   already validated for BOTH the uniform-toy-constant archive,
!>   results/hall_dt_stability_2026-08-21/, AND this session's own
!>   single-representative-constant test) produced a genuine, dt_hall-
!>   dependent numerical instability -- E_toroidal jumped ~10^13x in one
!>   step -- confirmed NOT a diagnostic artifact by rerunning at
!>   dt_hall=0.00001 (100x finer), where the SAME quantity grows
!>   smoothly instead. Root cause: F_HALL climbs by ~8 orders of
!>   magnitude across just the last few grid cells near the surface
!>   (not merely the single outermost node, which IS safely isolated --
!>   see the boundary-handling note below, still true but incomplete on
!>   its own); HALL_COURANT_TIMESTEP's own `dr/(F_HALL*|current|)`
!>   estimate is exactly the tool for a spatially steep constraint like
!>   this, so this driver now uses it instead of a hand-guessed fixed
!>   dt (which made a direct fixed-dt tmax=1000yr run either unstable
!>   at dt_hall=0.001 or, at a stable-but-tiny dt_hall, revive the
!>   "1.5-2 days of compute" wall-clock problem the 2026-08-21 archive
!>   was built specifically to avoid).
!>
!> Grid/composition risk assessment, still valid on its own terms (see
!> this session's plan addendum, .claude/plans/i-would-like-to-
!> quizzical-comet.md, "Addendum (2026-08-24)"): Rmag reaches
!> ~10^5-10^6 at the single outermost grid node (the literal stellar
!> surface) but is <3 everywhere else on this project's N_R=40 uniform
!> grid -- confirmed both by direct computation and by tracing
!> DIFFUSION_INIT/HALL_SUBSTEPS/HALL_INDUCTION_RHS's own boundary
!> handling (the extreme node's ETA_PROFILE/F_HALL_PROFILE values are
!> either fully overwritten by the vacuum-BC matrix row or immediately
!> zeroed after every explicit substep, with no cross-row leakage since
!> F_HALL is applied strictly elementwise). That reasoning correctly
!> ruled out the single boundary node as a direct instability source --
!> it just didn't cover the steep-gradient INTERIOR cells next to it,
!> which is what actually bit.
!>
!> Usage: mhdvsh_hall_crust_profile <n_steps> <n_sub> [cfl_safety]
!>   [energy_data_path] [energy_by_l_path] [checkpoint_path] [resume_path]
!> cfl_safety defaults to 0.5 (matching mhdvsh_hall_adaptive.f90's own
!> choice) if omitted; argument meanings for the path arguments
!> identical to app/mhdvsh_hall.f90's own.
!>
!> @warning n_steps is now an upper SAFETY CAP, not the real target
!>   (added 2026-08-24, see TIMESTEPPER::RUN_ADAPTIVE's own T_MAX
!>   parameter): the run stops at T_MAX_TARGET=1000 simulated years
!>   first in the expected case, clipping the final step's dt so T
!>   lands exactly on 1000 yr rather than overshooting. Pick n_steps
!>   generously (comfortably more than the actual expected step count)
!>   so the safety cap doesn't cut the run short before T_MAX_TARGET.
!>
!> Units: same convention as mhdvsh_hall.f90 (UNITS::ENERGY_UNIT_ERG/
!> POWER_UNIT_ERG_PER_S applied only at this reporting layer).
PROGRAM MHDVSH_HALL_CRUST_PROFILE
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN, IEEE_IS_NAN
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE RADIAL_OPERATORS,   ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE FIELD_DIAGNOSTICS,  ONLY: TOTAL_MAGNETIC_ENERGY, TOTAL_POLOIDAL_MAGNETIC_ENERGY, &
                               TOTAL_TOROIDAL_MAGNETIC_ENERGY, JOULE_DISSIPATION_RATE, &
                               POYNTING_FLUX_RATE, HALL_POYNTING_FLUX_RATE, &
                               POLOIDAL_MAGNETIC_ENERGY_BY_L, TOROIDAL_MAGNETIC_ENERGY_BY_L, &
                               HALL_COURANT_TIMESTEP
USE DIFFUSION_REGIME,   ONLY: DIFFUSION_STATE_T
USE HALL_REGIME,        ONLY: HALL_INIT, HALL_ADVANCE, HALL_COMPUTE_DT, HALL_SET_DT
USE TIMESTEPPER,        ONLY: RUN_ADAPTIVE
USE UNITS,              ONLY: ENERGY_UNIT_ERG, POWER_UNIT_ERG_PER_S
USE IO_CHECKPOINT,      ONLY: WRITE_CHECKPOINT, READ_CHECKPOINT
USE TOV_SOLVER,         ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
USE CRUST_CONDUCTIVITY, ONLY: ETA_AND_F_HALL_ON_GRID, FIND_TRUNCATION_RADIUS
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R     = 40
INTEGER(KIND=i4), PARAMETER :: LMAX    = 30
REAL(KIND=dp),    PARAMETER :: NBAR_CENTRAL = 0.5447307_dp   ! fm**-3, M=1.40 reference star, matches mhdvsh_tov.f90
CHARACTER(LEN=*), PARAMETER :: EOS_PATH = 'data/eos/APR_EOS_Cat.dat'
INTEGER(KIND=i4), PARAMETER :: NPOINTS_TOV = 414              ! matches mhdvsh_tov.f90
REAL(KIND=dp),    PARAMETER :: RHOL_CGS = 2.2E14_dp            ! core-crust boundary
! Outer-boundary truncation (added 2026-08-24, per user request; revised
! 2026-08-24 to use CRUST_CONDUCTIVITY::FIND_TRUNCATION_RADIUS instead
! of a hardcoded density -- see ROADMAP.md, "Analytical tasks",
! "Outer-boundary EOS/radial-grid truncation", for the open-questions
! writeup and why a fixed density doesn't generalize across T_KELVIN).
! R_MAX is set to the largest radius where the degeneracy parameter
! theta=kT/E_F is still <= THETA_MAX_TRUNCATE, computed live from the
! actual PROFILE and T_KELVIN below -- NOT a fixed density, since theta
! (not density) is the physically meaningful quantity and the crossing
! density shifts by orders of magnitude with T_KELVIN (confirmed
! directly: T=1e8K needs almost no truncation, T=3e9K needs
! substantially more, see FIND_TRUNCATION_RADIUS's own header). At
! T_KELVIN=1e9K this removes only ~1.18 m (0.094%) of the crust's
! radial extent and exactly one N_R=40 grid node (the same single node
! already confirmed harmless to the implicit/explicit boundary handling
! -- see this file's own header) but drops Rmag at the new boundary by
! ~12,300x (126 vs 1.5e6 at the true surface) -- empirically found to
! relieve the dt_hall~1.6e-7 yr floor the untruncated grid forced.
REAL(KIND=dp),    PARAMETER :: THETA_MAX_TRUNCATE = 1.3_dp
REAL(KIND=dp),    PARAMETER :: T_KELVIN = 1.0E9_dp             ! isothermal-crust default, matches every prior driver
INTEGER(KIND=i4), PARAMETER :: CHECKPOINT_EVERY = 200          ! coarser than mhdvsh_hall.f90's 50 -- long production run
INTEGER(KIND=i4), PARAMETER :: DATA_UNIT = 21
INTEGER(KIND=i4), PARAMETER :: ENERGY_L_UNIT = 24
REAL(KIND=dp),    PARAMETER :: PI = 3.14159265358979_dp
REAL(KIND=dp),    PARAMETER :: GRID_TOL = 1.0E-9_dp
REAL(KIND=dp),    PARAMETER :: DEFAULT_CFL_SAFETY = 0.5_dp   ! matches mhdvsh_hall_adaptive.f90
REAL(KIND=dp),    PARAMETER :: BLOWUP_FACTOR = 10.0_dp   ! same threshold as mhdvsh_hall_dt_sweep.f90
! Target simulated time for the production run (added 2026-08-24, see
! TIMESTEPPER::RUN_ADAPTIVE's own T_MAX parameter) -- with adaptive dt
! the step-to-time mapping can't be predicted precisely in advance, so
! this drives the actual stopping point; the CLI n_steps argument is now
! an upper SAFETY CAP only (pick it generously; the run stops at
! T_MAX_TARGET yr first in the expected case).
REAL(KIND=dp),    PARAMETER :: T_MAX_TARGET = 1000.0_dp

TYPE(TOV_PROFILE_T)     :: PROFILE
TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(RADIAL_OPERATOR_T) :: OPS
TYPE(DIFFUSION_STATE_T) :: STATE
REAL(KIND=dp), ALLOCATABLE :: ETA_PROFILE(:), F_HALL_PROFILE(:), N_E_PROFILE(:)
REAL(KIND=dp) :: R_MIN, R_MAX, ETA_MID, F_HALL_MID
REAL(KIND=dp) :: DT0, TC0, CFL_SAFETY, E0, E1
REAL(KIND=dp) :: PREV_E_POL, PREV_E_TOR, LAST_T
INTEGER(KIND=i4) :: N_STEPS, N_SUB
INTEGER(KIND=i4) :: IR, L
CHARACTER(LEN=1024) :: DATA_PATH, ENERGY_L_PATH, CHECKPOINT_PATH, RESUME_PATH, ARG
CHARACTER(LEN=16) :: COL_LABEL
LOGICAL :: WRITE_DATA, WRITE_ENERGY_L, WRITE_CHECKPOINT_FLAG, RESUMING, NEED_ON_STEP, BLOWN_UP
INTEGER(KIND=i4) :: ISTEP_START, N_STEPS_REMAINING
INTEGER(KIND=i4) :: N_R_CK, LMAX_CK
REAL(KIND=dp)    :: T_START, R_MIN_CK, R_MAX_CK

IF (COMMAND_ARGUMENT_COUNT() < 2) THEN
  WRITE(*,'(A)') 'Usage: mhdvsh_hall_crust_profile <n_steps> <n_sub> [cfl_safety] ' // &
    '[energy_data_path] [energy_by_l_path] [checkpoint_path] [resume_path]'
  STOP 1
END IF
CALL GET_COMMAND_ARGUMENT(1, ARG); READ(ARG,*) N_STEPS
CALL GET_COMMAND_ARGUMENT(2, ARG); READ(ARG,*) N_SUB
CFL_SAFETY = DEFAULT_CFL_SAFETY
IF (COMMAND_ARGUMENT_COUNT() >= 3) THEN
  CALL GET_COMMAND_ARGUMENT(3, ARG); READ(ARG,*) CFL_SAFETY
END IF

! ---- Solve the TOV profile fresh and derive the crust extent from it
! (RHOCGS<=RHOL_CGS mask, matching app/mhdvsh_tov.f90's own convention)
! rather than hardcoding R_MIN/R_MAX the way app/mhdvsh_hall.f90 does --
! stays self-consistent if the EOS/TOV solve ever changes.
CALL SOLVE_TOV_STAR(NBAR_CENTRAL, EOS_PATH, NPOINTS_TOV, PROFILE)
R_MIN = MINVAL(PROFILE%R, MASK=PROFILE%RHOCGS <= RHOL_CGS)
CALL FIND_TRUNCATION_RADIUS(PROFILE, T_KELVIN, R_MAX, THETA_MAX=THETA_MAX_TRUNCATE)

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL BUILD_RADIAL_OPERATORS(OPS, RGRID)

! ---- Sample eta(r)/f_H(r)/n_e(r) onto THIS driver's own simulation
! grid (N_R=40 points), not PROFILE's own much denser log(P)-uniform
! TOV resampling (~414 points) -- see ETA_AND_F_HALL_ON_GRID's own
! header for the interpolation scheme.
CALL ETA_AND_F_HALL_ON_GRID(PROFILE, RGRID%R, ETA_PROFILE, F_HALL_PROFILE, N_E_PROFILE, T_KELVIN=T_KELVIN)

! HALL_INIT's ETA/F_HALL arguments are non-optional even when the
! profile arrays are also given -- DIFFUSION_INIT/HALL_SUBSTEPS prefer
! the profile arrays whenever allocated (confirmed by reading both
! directly), so these scalar values are effectively unused placeholders
! for the DYNAMICS. A representative mid-crust value, not a meaningless
! one, in case anything ever inspects SAVED_ETA/SAVED_F_HALL directly;
! also used as the fallback scalar in the diagnostic calls below (the
! ETA_PROFILE/F_HALL_PROFILE arguments there take priority whenever
! given, see FIELD_DIAGNOSTICS' own header).
ETA_MID = ETA_PROFILE(N_R/2)
F_HALL_MID = F_HALL_PROFILE(N_R/2)

! ---- Seed IC or resume BEFORE HALL_INIT: the initial (Hall-CFL-
! limited) DT depends on the actual starting field, same ordering as
! app/mhdvsh_hall_adaptive.f90's own (compute TC0/DT0 from the seeded
! state, then call HALL_INIT with that DT0, not an arbitrary guess).
RESUMING = (COMMAND_ARGUMENT_COUNT() >= 7)
IF (RESUMING) THEN
  CALL GET_COMMAND_ARGUMENT(7, RESUME_PATH)
  CALL READ_CHECKPOINT(TRIM(RESUME_PATH), ISTEP_START, T_START, STATE%PHI, STATE%PSI, &
    N_R_CK, LMAX_CK, R_MIN_CK, R_MAX_CK)
  IF (N_R_CK /= N_R .OR. LMAX_CK /= LMAX .OR. &
      ABS(R_MIN_CK-R_MIN) > GRID_TOL .OR. ABS(R_MAX_CK-R_MAX) > GRID_TOL) THEN
    WRITE(*,'(A)') 'MHDVSH_HALL_CRUST_PROFILE: checkpoint grid shape does not match this run''s own parameters'
    WRITE(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5)') '  checkpoint: N_R=', N_R_CK, ' LMAX=', LMAX_CK, &
      ' R_MIN=', R_MIN_CK, ' R_MAX=', R_MAX_CK
    WRITE(*,'(A,I0,A,I0,A,ES12.5,A,ES12.5)') '  this run:   N_R=', N_R, ' LMAX=', LMAX, &
      ' R_MIN=', R_MIN, ' R_MAX=', R_MAX
    STOP 1
  END IF
ELSE
  ISTEP_START = 0_i4
  T_START = 0.0_dp
  CALL ALLOC_SPECTRAL_SCALAR(STATE%PHI, N_R, LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(STATE%PSI, N_R, LMAX)

  ! Single seed mode, identical to mhdvsh_hall.f90's own IC.
  DO IR = 1, N_R
    STATE%PHI%COEF(IR, YLM_INDEX(1_i4,0_i4)) = &
      CMPLX(SIN(PI*(RGRID%R(IR)-R_MIN)/(R_MAX-R_MIN)), 0.0_dp, KIND=dp)
  END DO
END IF
LAST_T = T_START

TC0 = HALL_COURANT_TIMESTEP(STATE%PHI, STATE%PSI, OPS, RGRID, F_HALL_MID, F_HALL_PROFILE=F_HALL_PROFILE)
DT0 = CFL_SAFETY * TC0
CALL HALL_INIT(RGRID, OPS, LMAX, ETA_MID, F_HALL_MID, DT0, N_SUB, &
  ETA_PROFILE=ETA_PROFILE, F_HALL_PROFILE=F_HALL_PROFILE)

N_STEPS_REMAINING = N_STEPS - ISTEP_START
IF (N_STEPS_REMAINING <= 0) THEN
  WRITE(*,'(A,I0,A,I0,A)') 'MHDVSH_HALL_CRUST_PROFILE: checkpoint step ', ISTEP_START, &
    ' already reaches or exceeds N_STEPS=', N_STEPS, ' -- nothing to do'
  STOP 1
END IF

WRITE_CHECKPOINT_FLAG = (COMMAND_ARGUMENT_COUNT() >= 6)
IF (WRITE_CHECKPOINT_FLAG) CALL GET_COMMAND_ARGUMENT(6, CHECKPOINT_PATH)

WRITE_DATA = (COMMAND_ARGUMENT_COUNT() >= 4)
IF (WRITE_DATA) THEN
  CALL GET_COMMAND_ARGUMENT(4, DATA_PATH)
  IF (RESUMING) THEN
    OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
  ELSE
    OPEN(UNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
    WRITE(DATA_UNIT,'(A)') '# MHDVSH_HALL_CRUST_PROFILE run inputs:'
    WRITE(DATA_UNIT,'(A,I0)')       '#   N_R          = ', N_R
    WRITE(DATA_UNIT,'(A,I0)')       '#   LMAX         = ', LMAX
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   R_MIN        = ', R_MIN, ' km (TOV crust extent, live)'
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   R_MAX        = ', R_MAX, ' km (TOV crust extent, live)'
    WRITE(DATA_UNIT,'(A)')          '#   dt           = DYNAMIC, Hall-CFL-limited (see dt_yr/tc_raw_yr columns)'
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   initial dt   = ', DT0, ' yr'
    WRITE(DATA_UNIT,'(A,ES16.8)')   '#   cfl_safety   = ', CFL_SAFETY
    WRITE(DATA_UNIT,'(A,I0)')       '#   N_SUB        = ', N_SUB
    WRITE(DATA_UNIT,'(A,I0)')       '#   dt_recompute_every = ', N_SUB
    WRITE(DATA_UNIT,'(A,I0)')       '#   N_STEPS      = ', N_STEPS
    WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   T_KELVIN     = ', T_KELVIN, ' K (isothermal crust)'
    WRITE(DATA_UNIT,'(A)') '#   eta(r)/f_H(r): full radial profile via CRUST_CONDUCTIVITY::' // &
      'ETA_AND_F_HALL_ON_GRID from a fresh SOLVE_TOV_STAR (NBAR_CENTRAL=0.5447307 fm**-3, ' // &
      'data/eos/APR_EOS_Cat.dat), NOT a uniform constant'
    WRITE(DATA_UNIT,'(A)') '#   seed IC: single mode Phi(l=1,m=0)=sin(pi*(r-R_MIN)/(R_MAX-R_MIN)), Psi=0'
    WRITE(DATA_UNIT,'(A)') '# step  t_yr  dt_yr  tc_raw_yr  E_poloidal_erg  E_toroidal_erg' // &
      '  EDOT_poloidal_erg_per_s  EDOT_toroidal_erg_per_s' // &
      '  joule_dissipation_rate_erg_per_s  poynting_flux_rate_erg_per_s' // &
      '  hall_poynting_flux_rate_erg_per_s  energy_balance_residual_erg_per_s'
  END IF
  PREV_E_POL = TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)
  PREV_E_TOR = TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)
  IF (.NOT. RESUMING) THEN
    ! EDOT_poloidal/EDOT_toroidal (dE/dt of each field component -- added
    ! 2026-08-25 per user request, the "only" quantities wanted in the
    ! energy-budget plot alongside the rate terms, not the raw E_pol/
    ! E_tor values themselves) are undefined at step 0 for the same
    ! reason the residual already is: no prior sample exists yet to
    ! backward-difference against. NaN here, not 0, for the same reason
    ! given below for the residual.
    WRITE(DATA_UNIT,'(I8,11ES16.8)') 0_i4, 0.0_dp, DT0, TC0, &
      PREV_E_POL*ENERGY_UNIT_ERG, PREV_E_TOR*ENERGY_UNIT_ERG, &
      IEEE_VALUE(1.0_dp, IEEE_QUIET_NAN), IEEE_VALUE(1.0_dp, IEEE_QUIET_NAN), &
      JOULE_DISSIPATION_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA_MID, ETA_PROFILE=ETA_PROFILE) &
        *POWER_UNIT_ERG_PER_S, &
      POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA_MID, ETA_PROFILE=ETA_PROFILE) &
        *POWER_UNIT_ERG_PER_S, &
      HALL_POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, F_HALL_MID, F_HALL_PROFILE=F_HALL_PROFILE) &
        *POWER_UNIT_ERG_PER_S, &
      IEEE_VALUE(1.0_dp, IEEE_QUIET_NAN)
  END IF

  WRITE_ENERGY_L = (COMMAND_ARGUMENT_COUNT() >= 5)
  IF (WRITE_ENERGY_L) THEN
    CALL GET_COMMAND_ARGUMENT(5, ENERGY_L_PATH)
    IF (RESUMING) THEN
      OPEN(UNIT=ENERGY_L_UNIT, FILE=TRIM(ENERGY_L_PATH), STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
    ELSE
      OPEN(UNIT=ENERGY_L_UNIT, FILE=TRIM(ENERGY_L_PATH), STATUS='REPLACE', ACTION='WRITE')
      WRITE(ENERGY_L_UNIT,'(A)', ADVANCE='NO') '# step  t'
      DO L = 1, LMAX
        WRITE(COL_LABEL,'(A,I0)') '  Epol_l', L
        WRITE(ENERGY_L_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
      END DO
      DO L = 1, LMAX
        WRITE(COL_LABEL,'(A,I0)') '  Etor_l', L
        WRITE(ENERGY_L_UNIT,'(A)', ADVANCE='NO') TRIM(COL_LABEL)
      END DO
      WRITE(ENERGY_L_UNIT,'(A)') ''
      CALL LOG_ENERGY_BY_L(0_i4, 0.0_dp)
    END IF
  END IF
END IF

NEED_ON_STEP = WRITE_DATA .OR. WRITE_CHECKPOINT_FLAG
BLOWN_UP = .FALSE.

! DT_RECOMPUTE_EVERY=1, not N_SUB (unlike mhdvsh_hall_adaptive.f90's own
! choice): empirically, this profile's Courant-limited dt can collapse
! by >100x after a SINGLE outer step near the steep near-surface F_HALL
! gradient (confirmed directly -- see this file's own header). Recomputing
! only every N_SUB steps (mhdvsh_hall_adaptive.f90's convention, tuned
! for the old uniform-F_HALL toy model's much gentler CFL drift) would
! run several steps at an already-stale, too-large dt before catching
! up. Re-factorizing DIFFUSION_INIT every step is cheap at this small
! N_R=40 -- the O(Nlm**3)-ish Hall mode-coupling RHS evaluation, not the
! diffusion solve, dominates per-step cost.
E0 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)
IF (NEED_ON_STEP) THEN
  CALL RUN_ADAPTIVE(HALL_ADVANCE, STATE, HALL_COMPUTE_DT, HALL_SET_DT, N_STEPS_REMAINING, &
    DT_RECOMPUTE_EVERY=1_i4, CFL_SAFETY=CFL_SAFETY, T_START=T_START, T_MAX=T_MAX_TARGET, ON_STEP=LOG_ENERGY)
ELSE
  CALL RUN_ADAPTIVE(HALL_ADVANCE, STATE, HALL_COMPUTE_DT, HALL_SET_DT, N_STEPS_REMAINING, &
    DT_RECOMPUTE_EVERY=1_i4, CFL_SAFETY=CFL_SAFETY, T_START=T_START, T_MAX=T_MAX_TARGET)
END IF
E1 = TOTAL_MAGNETIC_ENERGY(STATE%PHI, STATE%PSI, OPS, RGRID)

IF (WRITE_CHECKPOINT_FLAG) THEN
  CALL WRITE_CHECKPOINT(TRIM(CHECKPOINT_PATH), N_STEPS, LAST_T, STATE%PHI, STATE%PSI, RGRID)
END IF

IF (WRITE_DATA) CLOSE(DATA_UNIT)
IF (WRITE_ENERGY_L) CLOSE(ENERGY_L_UNIT)

WRITE(*,'(A)')            'MHD-VSH combined resistive+Hall regime -- real crust profile eta(r)/f_H(r), adaptive dt'
WRITE(*,'(A,ES12.5,A,ES12.5,A)') '  crust extent = ', R_MIN, ' -- ', R_MAX, ' km (live TOV solve)'
WRITE(*,'(A,I0)')         '  N_r          = ', N_R
WRITE(*,'(A,I0)')         '  Lmax         = ', LMAX
WRITE(*,'(A,ES12.5,A,ES12.5,A)') '  initial dt   = ', DT0, ' yr (cfl_safety*tc, tc=', TC0, ' yr)'
WRITE(*,'(A,I0,A,I0,A)')  '  N_steps      = ', N_STEPS, ' (dt recomputed every ', N_SUB, ' steps)'
WRITE(*,'(A,ES12.5,A)')   '  final t      = ', LAST_T, ' yr'
IF (RESUMING) WRITE(*,'(A,I0,A,A)') '  resumed from step ', ISTEP_START, ' via ', TRIM(RESUME_PATH)
WRITE(*,'(A,ES12.5,A)')   '  energy(t=0)  = ', E0*ENERGY_UNIT_ERG, ' erg'
WRITE(*,'(A,ES12.5,A)')   '  energy(t=end)= ', E1*ENERGY_UNIT_ERG, ' erg'
IF (BLOWN_UP) THEN
  WRITE(*,'(A)') '  RESULT: BLOWUP DETECTED (NaN or energy grew past BLOWUP_FACTOR)'
ELSE IF (E1 < E0) THEN
  WRITE(*,'(A)') '  RESULT: energy decreased (diffusion dissipating, as expected)'
ELSE
  WRITE(*,'(A)') '  RESULT: UNEXPECTED -- energy did not decrease'
END IF
IF (WRITE_DATA) WRITE(*,'(A,A)') '  time series written to ', TRIM(DATA_PATH)
IF (WRITE_CHECKPOINT_FLAG) WRITE(*,'(A,A)') '  final checkpoint written to ', TRIM(CHECKPOINT_PATH)

CONTAINS

!> Matches TIMESTEPPER::REGIME_ON_STEP_I; see mhdvsh_hall.f90's own
!> similar routine for the general explanation. Extended for the
!> dynamic-dt case: DT_NOW is recomputed fresh each call (matching
!> mhdvsh_hall_adaptive.f90's own LOG_STATE, since DT isn't a fixed
!> outer-scope value here), used both for the residual's backward-
!> difference denominator and logged directly (dt_yr/tc_raw_yr columns)
!> so the step-size trajectory is visible, not just inferred. LAST_T is
!> tracked here (not recomputed analytically after RUN_ADAPTIVE
!> returns, since DT varies) for the final unconditional checkpoint's
!> own T argument. Blow-up detection (NaN or energy past
!> BLOWUP_FACTOR*E0, same threshold as mhdvsh_hall_dt_sweep.f90) since
!> this driver is also meant for exploratory dt/n_sub/cfl_safety checks.
SUBROUTINE LOG_ENERGY(STATE, T, ISTEP)
  CLASS(*),         INTENT(IN) :: STATE
  REAL(KIND=dp),    INTENT(IN) :: T
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP
  REAL(KIND=dp) :: E_POL, E_TOR, EDOT_POL, EDOT_TOR, EDOT_J, EDOT_S, EDOT_S_HALL, RESIDUAL, E_NOW
  REAL(KIND=dp) :: TC_NOW, DT_NOW
  INTEGER(KIND=i4) :: GLOBAL_ISTEP
  GLOBAL_ISTEP = ISTEP_START + ISTEP
  LAST_T = T
  SELECT TYPE (STATE)
  TYPE IS (DIFFUSION_STATE_T)
    E_POL = TOTAL_POLOIDAL_MAGNETIC_ENERGY(STATE%PHI, OPS, RGRID)
    E_TOR = TOTAL_TOROIDAL_MAGNETIC_ENERGY(STATE%PSI, RGRID)
    E_NOW = E_POL + E_TOR
    IF (IEEE_IS_NAN(E_NOW) .OR. E_NOW > BLOWUP_FACTOR*E0) BLOWN_UP = .TRUE.

    IF (WRITE_DATA) THEN
      TC_NOW = HALL_COMPUTE_DT(STATE)
      DT_NOW = CFL_SAFETY * TC_NOW

      EDOT_POL    = (E_POL-PREV_E_POL)/DT_NOW
      EDOT_TOR    = (E_TOR-PREV_E_TOR)/DT_NOW
      EDOT_J      = JOULE_DISSIPATION_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA_MID, ETA_PROFILE=ETA_PROFILE)
      EDOT_S      = POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, ETA_MID, ETA_PROFILE=ETA_PROFILE)
      EDOT_S_HALL = HALL_POYNTING_FLUX_RATE(STATE%PHI, STATE%PSI, OPS, RGRID, F_HALL_MID, &
        F_HALL_PROFILE=F_HALL_PROFILE)
      RESIDUAL = (EDOT_POL+EDOT_TOR) - (EDOT_J+EDOT_S+EDOT_S_HALL)

      WRITE(DATA_UNIT,'(I8,11ES16.8)') GLOBAL_ISTEP, T, DT_NOW, TC_NOW, &
        E_POL*ENERGY_UNIT_ERG, E_TOR*ENERGY_UNIT_ERG, &
        EDOT_POL*POWER_UNIT_ERG_PER_S, EDOT_TOR*POWER_UNIT_ERG_PER_S, &
        EDOT_J*POWER_UNIT_ERG_PER_S, EDOT_S*POWER_UNIT_ERG_PER_S, &
        EDOT_S_HALL*POWER_UNIT_ERG_PER_S, RESIDUAL*POWER_UNIT_ERG_PER_S

      PREV_E_POL = E_POL
      PREV_E_TOR = E_TOR
      IF (WRITE_ENERGY_L) CALL LOG_ENERGY_BY_L(GLOBAL_ISTEP, T)
    END IF

    IF (WRITE_CHECKPOINT_FLAG) THEN
      IF (MOD(GLOBAL_ISTEP, CHECKPOINT_EVERY) == 0) THEN
        CALL WRITE_CHECKPOINT(TRIM(CHECKPOINT_PATH), GLOBAL_ISTEP, T, STATE%PHI, STATE%PSI, RGRID)
      END IF
    END IF
  END SELECT
END SUBROUTINE LOG_ENERGY

!> Identical to mhdvsh_hall.f90's own -- see that file's docstring.
SUBROUTINE LOG_ENERGY_BY_L(ISTEP_ARG, T_ARG)
  INTEGER(KIND=i4), INTENT(IN) :: ISTEP_ARG
  REAL(KIND=dp),    INTENT(IN) :: T_ARG
  REAL(KIND=dp) :: E_POL_L(0:LMAX), E_TOR_L(0:LMAX)
  INTEGER(KIND=i4) :: LL
  E_POL_L = POLOIDAL_MAGNETIC_ENERGY_BY_L(STATE%PHI, OPS, RGRID)
  E_TOR_L = TOROIDAL_MAGNETIC_ENERGY_BY_L(STATE%PSI, RGRID)
  WRITE(ENERGY_L_UNIT,'(I8,ES14.6)', ADVANCE='NO') ISTEP_ARG, T_ARG
  DO LL = 1, LMAX
    WRITE(ENERGY_L_UNIT,'(ES14.6)', ADVANCE='NO') E_POL_L(LL)*ENERGY_UNIT_ERG
  END DO
  DO LL = 1, LMAX
    WRITE(ENERGY_L_UNIT,'(ES14.6)', ADVANCE='NO') E_TOR_L(LL)*ENERGY_UNIT_ERG
  END DO
  WRITE(ENERGY_L_UNIT,'(A)') ''
END SUBROUTINE LOG_ENERGY_BY_L

END PROGRAM MHDVSH_HALL_CRUST_PROFILE
