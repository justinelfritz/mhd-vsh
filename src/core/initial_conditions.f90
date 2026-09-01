!> Initial-conditions module: a selectable-by-tag library of IC models,
!> replacing the ad hoc, duplicated seed-IC loop every Hall driver
!> currently hand-rolls (`mhdvsh_hall.f90`, `mhdvsh_hall_adaptive.f90`,
!> `mhdvsh_hall_crust_profile.f90` all copy-paste the identical single
!> `Phi(l=1,m=0)=sin(pi*(r-R_min)/(R_max-R_min))` seed inline). That
!> seed is also the confirmed root cause of the first-~10yr
!> energy-balance violation investigated earlier this session (it badly
!> violates the outer Robin BC, `boundary_conditions.f90`) -- see
!> ROADMAP.md analytical task 5 / programming task 10.
!>
!> First (and currently only) model: `'bessel_riccati'`, a user-specified
!> superposition of Riccati-Bessel modes `f_n(x) = a_n*x*j_n(x) +
!> b_n*x*y_n(x)` (built on SPHERICAL_BESSEL's stable j_n/y_n recursion),
!> one per (field, l, m, n, x_scale, a, b) entry in a caller-supplied
!> mode list. This is general superposition, NOT yet the true force-free
!> eigenmode of ROADMAP.md programming task 4 (which needs a
!> transcendental-eigenvalue root-finder this codebase doesn't have yet,
!> and a specific Phi/Psi amplitude coupling) -- this module is that
!> task's direct numerical prerequisite, not its completion.
!>
!> Mode-spectrum input is a Fortran NAMELIST file (confirmed with user,
!> 2026-08-29): zero new dependencies, matches this project's
!> demonstrated minimal-dependency policy (only LAPACK + locally-vendored
!> FORTVSH, no network fetch anywhere in CMakeLists.txt). New model tags
!> add their own namelist group without disturbing this one.
!>
!> @note Composable by design (confirmed with user, 2026-08-29): `&IC_SELECT`
!>   names a LIST of active tags, not one, and BUILD_INITIAL_CONDITION
!>   allocates PHI/PSI once and lets EVERY active tag's builder ACCUMULATE
!>   its own contribution into the same fields (BUILD_BESSEL_RICCATI
!>   already adds rather than overwrites, so it composes with itself or a
!>   future second category with no further change there). Only one
!>   category exists today, so this mechanism is exercised so far only
!>   with a single-entry tag list -- genuine cross-category composition
!>   (e.g. a force-free base plus a Bessel-Riccati perturbation) is
!>   untested until a second category exists, but the plumbing is ready.
MODULE INITIAL_CONDITIONS
USE KINDS,             ONLY: dp, i4
USE VSH,               ONLY: YLM_INDEX
USE GRID_RADIAL,       ONLY: RADIAL_GRID_T
USE FIELD_TYPES,       ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE SPHERICAL_BESSEL,  ONLY: SPHERICAL_BESSEL_JN_YN
IMPLICIT NONE
PRIVATE
PUBLIC :: IC_MODE_T, READ_IC_NAMELIST, BUILD_INITIAL_CONDITION

INTEGER(KIND=i4), PARAMETER :: MAX_MODES = 100_i4
INTEGER(KIND=i4), PARAMETER :: MAX_TAGS  = 10_i4

!> One Riccati-Bessel mode: contributes a_n*x*j_n(x)+b_n*x*y_n(x),
!> x=X_SCALE*r, to either PHI or PSI at spherical-harmonic (L,M).
TYPE :: IC_MODE_T
  CHARACTER(LEN=8) :: FIELD = 'phi'   ! 'phi' or 'psi'
  INTEGER(KIND=i4) :: L = 0, M = 0, N = 0
  REAL(KIND=dp)    :: X_SCALE = 1.0_dp
  REAL(KIND=dp)    :: A = 0.0_dp, B = 0.0_dp
END TYPE IC_MODE_T

CONTAINS

!> Reads the list of active tags from the file's &IC_SELECT group, then
!> each active tag's own model-specific group (only &BESSEL_RICCATI_MODES
!> exists so far). Allocates IC_TAGS(:) to exactly N_TAGS entries and
!> MODES(:) to exactly N_MODES entries (empty, size 0, if 'bessel_riccati'
!> isn't among the active tags).
!> @param IC_TAGS Every tag named in &IC_SELECT, in file order. Every
!>   entry must be a known tag -- fails loud (STOP) on typos or a
!>   not-yet-implemented category, rather than silently doing nothing in
!>   BUILD_INITIAL_CONDITION later.
SUBROUTINE READ_IC_NAMELIST(PATH, IC_TAGS, MODES)
  CHARACTER(LEN=*),              INTENT(IN)  :: PATH
  CHARACTER(LEN=64), ALLOCATABLE,INTENT(OUT) :: IC_TAGS(:)
  TYPE(IC_MODE_T), ALLOCATABLE,  INTENT(OUT) :: MODES(:)

  ! Fixed-size buffers: Fortran NAMELIST arrays need a known declared
  ! size before the READ (no auto-allocation). TAGS is the namelist's
  ! own user-facing variable name (file syntax: `TAGS = 'a', 'b'`),
  ! distinct from the IC_TAGS(:) output argument above.
  CHARACTER(LEN=64) :: TAGS(MAX_TAGS)
  INTEGER(KIND=i4)  :: UNIT, N_TAGS, N_MODES, I, ITAG

  CHARACTER(LEN=8) :: MODE_FIELD(MAX_MODES)
  INTEGER(KIND=i4) :: MODE_L(MAX_MODES), MODE_M(MAX_MODES), MODE_N(MAX_MODES)
  REAL(KIND=dp)    :: MODE_X_SCALE(MAX_MODES), MODE_A(MAX_MODES), MODE_B(MAX_MODES)

  NAMELIST /IC_SELECT/ N_TAGS, TAGS
  NAMELIST /BESSEL_RICCATI_MODES/ N_MODES, MODE_FIELD, MODE_L, MODE_M, MODE_N, &
    MODE_X_SCALE, MODE_A, MODE_B

  N_TAGS = 0_i4
  TAGS = ''
  N_MODES = 0_i4
  MODE_FIELD = 'phi'
  MODE_L = 0_i4; MODE_M = 0_i4; MODE_N = 0_i4
  MODE_X_SCALE = 1.0_dp; MODE_A = 0.0_dp; MODE_B = 0.0_dp

  OPEN(NEWUNIT=UNIT, FILE=TRIM(PATH), STATUS='OLD', ACTION='READ')
  READ(UNIT, NML=IC_SELECT)

  IF (N_TAGS <= 0_i4 .OR. N_TAGS > MAX_TAGS) THEN
    WRITE(*,'(A,I0,A,I0)') 'READ_IC_NAMELIST: N_TAGS=', N_TAGS, ' out of range 1..', MAX_TAGS
    STOP 1
  END IF
  ALLOCATE(IC_TAGS(N_TAGS))
  DO ITAG = 1, N_TAGS
    IC_TAGS(ITAG) = TRIM(TAGS(ITAG))
    SELECT CASE (TRIM(IC_TAGS(ITAG)))
    CASE ('bessel_riccati')
      CONTINUE
    CASE DEFAULT
      WRITE(*,'(A,A,A)') 'READ_IC_NAMELIST: unknown tag "', TRIM(IC_TAGS(ITAG)), '"'
      STOP 1
    END SELECT
  END DO

  IF (ANY(IC_TAGS == 'bessel_riccati')) THEN
    REWIND(UNIT)
    READ(UNIT, NML=BESSEL_RICCATI_MODES)
    IF (N_MODES <= 0_i4 .OR. N_MODES > MAX_MODES) THEN
      WRITE(*,'(A,I0,A,I0)') 'READ_IC_NAMELIST: N_MODES=', N_MODES, &
        ' out of range 1..', MAX_MODES
      STOP 1
    END IF
    ALLOCATE(MODES(N_MODES))
    DO I = 1, N_MODES
      MODES(I)%FIELD   = TRIM(MODE_FIELD(I))
      MODES(I)%L       = MODE_L(I)
      MODES(I)%M       = MODE_M(I)
      MODES(I)%N       = MODE_N(I)
      MODES(I)%X_SCALE = MODE_X_SCALE(I)
      MODES(I)%A       = MODE_A(I)
      MODES(I)%B       = MODE_B(I)
    END DO
  ELSE
    ALLOCATE(MODES(0))
  END IF
  CLOSE(UNIT)
END SUBROUTINE READ_IC_NAMELIST

!> Allocates PHI/PSI once (size N_R=RGRID%N, degree LMAX), then lets
!> EVERY active tag's builder ACCUMULATE its own contribution into the
!> same fields, in list order -- the composability mechanism (confirmed
!> with user, 2026-08-29): a caller can activate more than one category
!> at once (e.g. a future force-free base plus a Bessel-Riccati
!> perturbation) and get their sum, not just the last one applied.
!> New tags add a new CASE calling their own accumulating builder, not a
!> restructure -- matches this project's own preference for explicit
!> dispatch over added abstraction (e.g. REGIME_INTERFACE's "no runtime
!> dispatch" choice for regimes).
SUBROUTINE BUILD_INITIAL_CONDITION(IC_TAGS, MODES, RGRID, LMAX, PHI, PSI)
  CHARACTER(LEN=*),         INTENT(IN)  :: IC_TAGS(:)
  TYPE(IC_MODE_T),          INTENT(IN)  :: MODES(:)
  TYPE(RADIAL_GRID_T),      INTENT(IN)  :: RGRID
  INTEGER(KIND=i4),         INTENT(IN)  :: LMAX
  TYPE(SPECTRAL_SCALAR_T),  INTENT(OUT) :: PHI, PSI
  INTEGER(KIND=i4) :: ITAG

  CALL ALLOC_SPECTRAL_SCALAR(PHI, RGRID%N, LMAX)
  CALL ALLOC_SPECTRAL_SCALAR(PSI, RGRID%N, LMAX)

  DO ITAG = 1, SIZE(IC_TAGS)
    SELECT CASE (TRIM(IC_TAGS(ITAG)))
    CASE ('bessel_riccati')
      CALL BUILD_BESSEL_RICCATI(MODES, RGRID, LMAX, PHI, PSI)
    CASE DEFAULT
      WRITE(*,'(A,A,A)') 'BUILD_INITIAL_CONDITION: unknown tag "', TRIM(IC_TAGS(ITAG)), '"'
      STOP 1
    END SELECT
  END DO
END SUBROUTINE BUILD_INITIAL_CONDITION

!> 'bessel_riccati' model: accumulates each mode's
!> f_n(X_SCALE*r) = A*x*j_n(x) + B*x*y_n(x) into PHI or PSI at (L,M),
!> across every radial grid point. Modes with the same (FIELD,L,M)
!> superpose (accumulated, not overwritten) -- lets a caller build a
!> genuine multi-n spectrum at one (l,m), not just one radial function
!> per mode slot.
SUBROUTINE BUILD_BESSEL_RICCATI(MODES, RGRID, LMAX, PHI, PSI)
  TYPE(IC_MODE_T),          INTENT(IN)    :: MODES(:)
  TYPE(RADIAL_GRID_T),      INTENT(IN)    :: RGRID
  INTEGER(KIND=i4),         INTENT(IN)    :: LMAX
  TYPE(SPECTRAL_SCALAR_T),  INTENT(INOUT) :: PHI, PSI

  REAL(KIND=dp), ALLOCATABLE :: J(:), Y(:)
  REAL(KIND=dp) :: X, F_N
  INTEGER(KIND=i4) :: IMODE, IR, IDX

  DO IMODE = 1, SIZE(MODES)
    ASSOCIATE (MODE => MODES(IMODE))
      IF (MODE%L < 0_i4 .OR. MODE%L > LMAX) THEN
        WRITE(*,'(A,I0,A,I0,A,I0)') 'BUILD_BESSEL_RICCATI: mode ', IMODE, &
          ' has l=', MODE%L, ' outside 0..LMAX=', LMAX
        STOP 1
      END IF
      IF (ABS(MODE%M) > MODE%L) THEN
        WRITE(*,'(A,I0,A,I0,A,I0)') 'BUILD_BESSEL_RICCATI: mode ', IMODE, &
          ' has |m|=', ABS(MODE%M), ' > l=', MODE%L
        STOP 1
      END IF

      ALLOCATE(J(0:MODE%N), Y(0:MODE%N))
      IDX = YLM_INDEX(MODE%L, MODE%M)
      DO IR = 1, RGRID%N
        X = MODE%X_SCALE * RGRID%R(IR)
        CALL SPHERICAL_BESSEL_JN_YN(X, MODE%N, J, Y)
        F_N = MODE%A*X*J(MODE%N) + MODE%B*X*Y(MODE%N)
        SELECT CASE (TRIM(MODE%FIELD))
        CASE ('phi')
          PHI%COEF(IR,IDX) = PHI%COEF(IR,IDX) + CMPLX(F_N, 0.0_dp, KIND=dp)
        CASE ('psi')
          PSI%COEF(IR,IDX) = PSI%COEF(IR,IDX) + CMPLX(F_N, 0.0_dp, KIND=dp)
        CASE DEFAULT
          WRITE(*,'(A,I0,A,A,A)') 'BUILD_BESSEL_RICCATI: mode ', IMODE, &
            ' has unknown field "', TRIM(MODE%FIELD), '" (want phi or psi)'
          STOP 1
        END SELECT
      END DO
      DEALLOCATE(J, Y)
    END ASSOCIATE
  END DO
END SUBROUTINE BUILD_BESSEL_RICCATI

END MODULE INITIAL_CONDITIONS
