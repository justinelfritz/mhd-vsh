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
!> Two models:
!>  - `'bessel_riccati'`: a user-specified superposition of Riccati-Bessel
!>    modes `f_n(x) = a_n*x*j_n(x) + b_n*x*y_n(x)` (built on
!>    SPHERICAL_BESSEL's stable j_n/y_n recursion), one per (field, l, m,
!>    n, x_scale, a, b) entry in a caller-supplied mode list -- arbitrary,
!>    not generally BC-consistent or force-free.
!>  - `'force_free_bessel'` (2026-10-15): the true force-free Bessel-Riccati
!>    eigenmode, ROADMAP.md programming task 4's completion -- solves
!>    FORCE_FREE_BESSEL::SOLVE_FORCE_FREE_EIGENVALUE for the eigenvalue
!>    `mu_l` and Riccati-Bessel mix, then sets `Psi_l(r) = -mu_l*Phi_l(r)`
!>    (the force-free coupling, Igoshev/Elfritz/Popov 2016 Eq. A10) --
!>    unlike 'bessel_riccati', this ALWAYS contributes to both Phi and Psi
!>    together (they are not independently selectable per mode).
!>    Per-mode `OUTER_BC` (2026-10-15) selects which outer condition the
!>    eigenvalue is solved against -- this project's real outer BC
!>    (Robin(Phi)+Dirichlet(Psi)=0) can't both be satisfied exactly once
!>    the force-free coupling slaves Psi to Phi (see FORCE_FREE_BESSEL's
!>    own header for the full derivation of both options):
!>    `'double_dirichlet'` (default) forces `Phi(R_out)=0` (hence
!>    `B_r(R_out)=0` too -- a special case, not the generic multipole);
!>    `'arbitrary_value'` instead normalizes `Phi(R_out)` to an arbitrary
!>    nonzero value (the paper's own Eqs. A17-A19 scheme), giving the
!>    physically typical nonzero-`B_r`-at-the-surface case, at the cost of
!>    a generally nonzero `Psi(R_out)` (and no guarantee the exterior
!>    vacuum potential decays, in either case).
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
!>   its own contribution into the same fields. Now exercised with a real
!>   second category (2026-10-15, 'force_free_bessel') -- both builders add
!>   rather than overwrite, so e.g. a force-free base plus a Bessel-Riccati
!>   perturbation compose with no further change to this dispatch.
MODULE INITIAL_CONDITIONS
USE KINDS,             ONLY: dp, i4
USE VSH,               ONLY: YLM_INDEX
USE GRID_RADIAL,       ONLY: RADIAL_GRID_T
USE FIELD_TYPES,       ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE SPHERICAL_BESSEL,  ONLY: SPHERICAL_BESSEL_JN_YN
USE FORCE_FREE_BESSEL, ONLY: SOLVE_FORCE_FREE_EIGENVALUE, &
                              FF_OUTER_BC_DOUBLE_DIRICHLET, FF_OUTER_BC_ARBITRARY_VALUE
IMPLICIT NONE
PRIVATE
PUBLIC :: IC_MODE_T, FF_MODE_T, READ_IC_NAMELIST, BUILD_INITIAL_CONDITION

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

!> One force-free eigenmode at degree (L,M): ROOT_INDEX selects which
!> eigenvalue (1=smallest/"fundamental"), AMPLITUDE scales the Riccati-
!> Bessel a-coefficient (b is fixed by the inner BC, not independently
!> settable -- see FORCE_FREE_BESSEL's own header). OUTER_BC is
!> 'double_dirichlet' (default, Phi(R_out)=0) or 'arbitrary_value'
!> (Phi(R_out) normalized to a nonzero value instead -- see this module's
!> own header and FORCE_FREE_BESSEL's for the full trade-off).
!>
!> @note ROOT_INDEX>1 is rejected (hard STOP) unless ALLOW_HIGHER_ROOTS is
!>   set (2026-10-16, per user request): ROOT_INDEX=1 is the ONLY index
!>   guaranteed free of interior zero crossings for 'double_dirichlet'
!>   (a regular Sturm-Liouville eigenvalue problem -- Sturm's own
!>   oscillation theorem guarantees the n-th eigenfunction has EXACTLY
!>   n-1 interior zeros, confirmed numerically). For 'arbitrary_value',
!>   ROOT_INDEX=1 is empirically always zero-crossing-free too (checked
!>   l=1,2,3,5,10), but that case is NOT a standard self-adjoint BVP (the
!>   outer condition is an amplitude normalization, not a homogeneous
!>   BC), so there's no theorem backing it, and -- confirmed numerically
!>   -- the correspondence isn't 1-to-1 the way it is for
!>   'double_dirichlet' (e.g. at l=10, ROOT_INDEX=1..4 are ALL
!>   zero-crossing-free; the first interior node doesn't appear until
!>   ROOT_INDEX=5). ALLOW_HIGHER_ROOTS exists for deliberate exploration
!>   (e.g. test/test_force_free_bessel.f90's own root-ladder checks call
!>   FORCE_FREE_BESSEL::SOLVE_FORCE_FREE_EIGENVALUE directly, bypassing
!>   this module's policy layer entirely) -- it does NOT itself guarantee
!>   zero crossings appear; it only removes this module's refusal to
!>   build the mode at all.
TYPE :: FF_MODE_T
  INTEGER(KIND=i4)  :: L = 1, M = 0, ROOT_INDEX = 1
  REAL(KIND=dp)     :: AMPLITUDE = 1.0_dp
  CHARACTER(LEN=32) :: OUTER_BC = 'double_dirichlet'
  LOGICAL           :: ALLOW_HIGHER_ROOTS = .FALSE.
END TYPE FF_MODE_T

CONTAINS

!> Reads the list of active tags from the file's &IC_SELECT group, then
!> each active tag's own model-specific group (&BESSEL_RICCATI_MODES
!> and/or &FORCE_FREE_BESSEL_MODES). Allocates IC_TAGS(:) to exactly
!> N_TAGS entries, and MODES(:)/FF_MODES(:) each to their own N_MODES
!> (empty, size 0, if that tag isn't among the active ones).
!> @param IC_TAGS Every tag named in &IC_SELECT, in file order. Every
!>   entry must be a known tag -- fails loud (STOP) on typos or a
!>   not-yet-implemented category, rather than silently doing nothing in
!>   BUILD_INITIAL_CONDITION later.
SUBROUTINE READ_IC_NAMELIST(PATH, IC_TAGS, MODES, FF_MODES)
  CHARACTER(LEN=*),              INTENT(IN)  :: PATH
  CHARACTER(LEN=64), ALLOCATABLE,INTENT(OUT) :: IC_TAGS(:)
  TYPE(IC_MODE_T), ALLOCATABLE,  INTENT(OUT) :: MODES(:)
  TYPE(FF_MODE_T), ALLOCATABLE,  INTENT(OUT) :: FF_MODES(:)

  ! Fixed-size buffers: Fortran NAMELIST arrays need a known declared
  ! size before the READ (no auto-allocation). TAGS is the namelist's
  ! own user-facing variable name (file syntax: `TAGS = 'a', 'b'`),
  ! distinct from the IC_TAGS(:) output argument above.
  CHARACTER(LEN=64) :: TAGS(MAX_TAGS)
  INTEGER(KIND=i4)  :: UNIT, N_TAGS, N_MODES, I, ITAG

  CHARACTER(LEN=8) :: MODE_FIELD(MAX_MODES)
  INTEGER(KIND=i4) :: MODE_L(MAX_MODES), MODE_M(MAX_MODES), MODE_N(MAX_MODES)
  REAL(KIND=dp)    :: MODE_X_SCALE(MAX_MODES), MODE_A(MAX_MODES), MODE_B(MAX_MODES)

  INTEGER(KIND=i4)  :: N_FF_MODES
  INTEGER(KIND=i4)  :: FF_MODE_L(MAX_MODES), FF_MODE_M(MAX_MODES), FF_MODE_ROOT_INDEX(MAX_MODES)
  REAL(KIND=dp)     :: FF_MODE_AMPLITUDE(MAX_MODES)
  CHARACTER(LEN=32) :: FF_MODE_OUTER_BC(MAX_MODES)
  LOGICAL           :: FF_MODE_ALLOW_HIGHER_ROOTS(MAX_MODES)

  NAMELIST /IC_SELECT/ N_TAGS, TAGS
  NAMELIST /BESSEL_RICCATI_MODES/ N_MODES, MODE_FIELD, MODE_L, MODE_M, MODE_N, &
    MODE_X_SCALE, MODE_A, MODE_B
  NAMELIST /FORCE_FREE_BESSEL_MODES/ N_FF_MODES, FF_MODE_L, FF_MODE_M, &
    FF_MODE_ROOT_INDEX, FF_MODE_AMPLITUDE, FF_MODE_OUTER_BC, FF_MODE_ALLOW_HIGHER_ROOTS

  N_TAGS = 0_i4
  TAGS = ''
  N_MODES = 0_i4
  MODE_FIELD = 'phi'
  MODE_L = 0_i4; MODE_M = 0_i4; MODE_N = 0_i4
  MODE_X_SCALE = 1.0_dp; MODE_A = 0.0_dp; MODE_B = 0.0_dp
  N_FF_MODES = 0_i4
  FF_MODE_L = 1_i4; FF_MODE_M = 0_i4; FF_MODE_ROOT_INDEX = 1_i4
  FF_MODE_AMPLITUDE = 1.0_dp
  FF_MODE_OUTER_BC = 'double_dirichlet'
  FF_MODE_ALLOW_HIGHER_ROOTS = .FALSE.

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
    CASE ('bessel_riccati', 'force_free_bessel')
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

  IF (ANY(IC_TAGS == 'force_free_bessel')) THEN
    REWIND(UNIT)
    READ(UNIT, NML=FORCE_FREE_BESSEL_MODES)
    IF (N_FF_MODES <= 0_i4 .OR. N_FF_MODES > MAX_MODES) THEN
      WRITE(*,'(A,I0,A,I0)') 'READ_IC_NAMELIST: N_FF_MODES=', N_FF_MODES, &
        ' out of range 1..', MAX_MODES
      STOP 1
    END IF
    ALLOCATE(FF_MODES(N_FF_MODES))
    DO I = 1, N_FF_MODES
      FF_MODES(I)%L          = FF_MODE_L(I)
      FF_MODES(I)%M          = FF_MODE_M(I)
      FF_MODES(I)%ROOT_INDEX = FF_MODE_ROOT_INDEX(I)
      FF_MODES(I)%AMPLITUDE  = FF_MODE_AMPLITUDE(I)
      FF_MODES(I)%OUTER_BC   = TRIM(FF_MODE_OUTER_BC(I))
      FF_MODES(I)%ALLOW_HIGHER_ROOTS = FF_MODE_ALLOW_HIGHER_ROOTS(I)
      SELECT CASE (TRIM(FF_MODES(I)%OUTER_BC))
      CASE ('double_dirichlet', 'arbitrary_value')
        CONTINUE
      CASE DEFAULT
        WRITE(*,'(A,I0,A,A,A)') 'READ_IC_NAMELIST: force_free_bessel mode ', I, &
          ' has unknown OUTER_BC "', TRIM(FF_MODES(I)%OUTER_BC), '" (want double_dirichlet or arbitrary_value)'
        STOP 1
      END SELECT
    END DO
  ELSE
    ALLOCATE(FF_MODES(0))
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
SUBROUTINE BUILD_INITIAL_CONDITION(IC_TAGS, MODES, FF_MODES, RGRID, LMAX, PHI, PSI)
  CHARACTER(LEN=*),         INTENT(IN)  :: IC_TAGS(:)
  TYPE(IC_MODE_T),          INTENT(IN)  :: MODES(:)
  TYPE(FF_MODE_T),          INTENT(IN)  :: FF_MODES(:)
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
    CASE ('force_free_bessel')
      CALL BUILD_FORCE_FREE_BESSEL(FF_MODES, RGRID, LMAX, PHI, PSI)
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

!> 'force_free_bessel' model: the true force-free eigenmode (ROADMAP.md
!> programming task 4). For each mode, solves FORCE_FREE_BESSEL::
!> SOLVE_FORCE_FREE_EIGENVALUE for the eigenvalue MU and Riccati-Bessel
!> mix (A,B) satisfying the double-Dirichlet BC at degree L, scales both
!> by MODE%AMPLITUDE, then accumulates BOTH Phi (from A,B directly) AND
!> Psi (= -MU*Phi pointwise, the force-free coupling) at (L,M) -- unlike
!> BUILD_BESSEL_RICCATI, this always touches both fields together, never
!> just one.
SUBROUTINE BUILD_FORCE_FREE_BESSEL(FF_MODES, RGRID, LMAX, PHI, PSI)
  TYPE(FF_MODE_T),          INTENT(IN)    :: FF_MODES(:)
  TYPE(RADIAL_GRID_T),      INTENT(IN)    :: RGRID
  INTEGER(KIND=i4),         INTENT(IN)    :: LMAX
  TYPE(SPECTRAL_SCALAR_T),  INTENT(INOUT) :: PHI, PSI

  REAL(KIND=dp), ALLOCATABLE :: J(:), Y(:)
  REAL(KIND=dp) :: MU, A_COEF, B_COEF, X, PHI_VAL
  INTEGER(KIND=i4) :: IMODE, IR, IDX, OUTER_BC_FLAG

  DO IMODE = 1, SIZE(FF_MODES)
    ASSOCIATE (MODE => FF_MODES(IMODE))
      IF (MODE%L < 1_i4 .OR. MODE%L > LMAX) THEN
        WRITE(*,'(A,I0,A,I0,A,I0)') 'BUILD_FORCE_FREE_BESSEL: mode ', IMODE, &
          ' has l=', MODE%L, ' outside 1..LMAX=', LMAX
        STOP 1
      END IF
      IF (ABS(MODE%M) > MODE%L) THEN
        WRITE(*,'(A,I0,A,I0,A,I0)') 'BUILD_FORCE_FREE_BESSEL: mode ', IMODE, &
          ' has |m|=', ABS(MODE%M), ' > l=', MODE%L
        STOP 1
      END IF
      IF (MODE%ROOT_INDEX > 1_i4 .AND. .NOT. MODE%ALLOW_HIGHER_ROOTS) THEN
        WRITE(*,'(A,I0,A,I0,A)') 'BUILD_FORCE_FREE_BESSEL: mode ', IMODE, &
          ' requests ROOT_INDEX=', MODE%ROOT_INDEX, ' > 1'
        WRITE(*,'(A)') '  ROOT_INDEX=1 is the only index PROVEN free of interior zero' // &
          ' crossings for double_dirichlet (Sturm-Liouville oscillation theorem), and the' // &
          ' only index EMPIRICALLY confirmed free of them for arbitrary_value (not proven' // &
          ' there -- see this module''s own header). Set ALLOW_HIGHER_ROOTS=.TRUE. to override.'
        STOP 1
      END IF

      SELECT CASE (TRIM(MODE%OUTER_BC))
      CASE ('double_dirichlet')
        OUTER_BC_FLAG = FF_OUTER_BC_DOUBLE_DIRICHLET
      CASE ('arbitrary_value')
        OUTER_BC_FLAG = FF_OUTER_BC_ARBITRARY_VALUE
      CASE DEFAULT
        WRITE(*,'(A,I0,A,A,A)') 'BUILD_FORCE_FREE_BESSEL: mode ', IMODE, &
          ' has unknown OUTER_BC "', TRIM(MODE%OUTER_BC), '" (want double_dirichlet or arbitrary_value)'
        STOP 1
      END SELECT

      CALL SOLVE_FORCE_FREE_EIGENVALUE(MODE%L, RGRID%R_MIN, RGRID%R_MAX, &
        OUTER_BC_FLAG, MODE%ROOT_INDEX, MU, A_COEF, B_COEF)
      A_COEF = MODE%AMPLITUDE * A_COEF
      B_COEF = MODE%AMPLITUDE * B_COEF

      ALLOCATE(J(0:MODE%L), Y(0:MODE%L))
      IDX = YLM_INDEX(MODE%L, MODE%M)
      DO IR = 1, RGRID%N
        X = MU * RGRID%R(IR)
        CALL SPHERICAL_BESSEL_JN_YN(X, MODE%L, J, Y)
        PHI_VAL = A_COEF*X*J(MODE%L) + B_COEF*X*Y(MODE%L)
        PHI%COEF(IR,IDX) = PHI%COEF(IR,IDX) + CMPLX(PHI_VAL, 0.0_dp, KIND=dp)
        PSI%COEF(IR,IDX) = PSI%COEF(IR,IDX) + CMPLX(-MU*PHI_VAL, 0.0_dp, KIND=dp)
      END DO
      DEALLOCATE(J, Y)
    END ASSOCIATE
  END DO
END SUBROUTINE BUILD_FORCE_FREE_BESSEL

END MODULE INITIAL_CONDITIONS
