MODULE TOV_SOLVER
!> Tolman-Oppenheimer-Volkoff stellar-structure solver -- ported from
!> ~/Desktop/EOSNS/src/nstot.f (driver logic) and derivs.f (the TOV
!> right-hand side), built on ODE_INTEGRATOR (generic RK45) and
!> EOS_TABLE (crust/low-density EOS interpolation). Faithful port: every
!> geometrized-unit conversion constant is copied byte-for-byte from the
!> original (see LOAD_EOS_TABLE's own warning) -- these are validated
!> against real, already-solved reference output (`~/Desktop/EOSNS/
!> fort.34`/`PL.DAT`, an M=1.40 star, dated 2015), not re-derived.
!>
!> Structural (non-numerical) changes from the original: SOLVE_TOV_STAR
!> returns a TOV_PROFILE_T via INTENT(OUT) instead of writing files
!> directly (this codebase's "physics module never writes files,
!> drivers do" convention -- FIELD_DIAGNOSTICS is the precedent); the
!> EOS table lookup needed inside the DERIVS callback (which
!> ODE_INTEGRATOR's DERIVS_I interface has no room to pass explicitly)
!> is threaded through as SAVE'd module-level state, set once per
!> SOLVE_TOV_STAR call -- the same pattern HALL_REGIME already uses for
!> RGRID/OPS in its own ADVANCE callback.
USE KINDS,          ONLY: dp, i4
USE ODE_INTEGRATOR, ONLY: DERIVS_I, ODEINT
USE EOS_TABLE,      ONLY: EOS_TABLE_T, LOAD_EOS_TABLE, EOS_AT_DENSITY, EOS_AT_PRESSURE
IMPLICIT NONE
PRIVATE
PUBLIC :: TOV_PROFILE_T, SOLVE_TOV_STAR

!> One radial profile of a solved TOV star, from the core-crust-facing
!> integration limit out to the surface (pressure -> ~0). Composition
!> fields (AH/ZH/XH/YN/YE) follow nstot.f's own convention exactly:
!> AH<=0 (with ZH=YE, XH=0) signals homogeneous matter (no nuclei --
!> deep crust/core, where Potekhin's lattice-conductivity physics does
!> not apply), matching how nstot.f itself gates its own conductivity
!> call.
TYPE :: TOV_PROFILE_T
  INTEGER(KIND=i4) :: N = 0
  REAL(KIND=dp) :: MASS_MSUN = 0.0_dp   !! Total gravitational mass, solar masses.
  REAL(KIND=dp) :: RADIUS_KM = 0.0_dp   !! Surface radius, km.
  REAL(KIND=dp), ALLOCATABLE :: R(:)      !! Radius, km.
  REAL(KIND=dp), ALLOCATABLE :: RHOCGS(:) !! Mass density, `g/cm**3`.
  REAL(KIND=dp), ALLOCATABLE :: PCGS(:)   !! Pressure, `dyn/cm**2`.
  REAL(KIND=dp), ALLOCATABLE :: NBFM(:)   !! Baryon number density, `fm**-3`.
  REAL(KIND=dp), ALLOCATABLE :: NEL(:)    !! Electron number density, `fm**-3`.
  REAL(KIND=dp), ALLOCATABLE :: AH(:), ZH(:), XH(:), YE(:), YN(:)
  !> Raw EOS-table mass number (nstot.f's `a`, pre-COMPOSITION), the
  !> actual field its own conductivity gate `IF (a.gt.0.d0)` tests --
  !> unlike AH (always >0: COMPOSITION substitutes a dummy AH=1 for
  !> homogeneous rows), this is <=0 exactly where nstot.f itself skips
  !> calling potekhinc (deep crust/core, no discrete nuclei). Consumed
  !> by CRUST_CONDUCTIVITY::ETA_AND_F_HALL_AT's own gate.
  REAL(KIND=dp), ALLOCATABLE :: A_TABLE(:)
END TYPE TOV_PROFILE_T

!> Per-call context for DERIVS_TOV (ODE_INTEGRATOR's DERIVS_I interface
!> has no room for extra arguments -- see module header).
TYPE(EOS_TABLE_T), SAVE :: ACTIVE_TABLE

CONTAINS

!> Solves a single TOV star of central density RHOCGS (`g/cm**3`) using
!> the EOS table at EOS_PATH, returning its radial profile sampled at
!> NPOINTS+2 points (NPOINTS interior steps in log-pressure, plus the
!> center and surface). Ported from nstot.f's full driver body.
!>
!> @param RHOCGS Central mass density, `g/cm**3` (e.g. 9.88d14 for the
!>   M=1.40 reference star nstot.f's own comments document).
!> @param EOS_PATH Path to a table in the `lowd-eos.ja.tab` format.
!> @param NPOINTS Number of log-pressure integration steps for the
!>   radial-profile pass (nstot.f prompts for this interactively; here
!>   it's an argument).
!> @param PROFILE Output: the solved radial profile (see TOV_PROFILE_T).
SUBROUTINE SOLVE_TOV_STAR(RHOCGS, EOS_PATH, NPOINTS, PROFILE)
  REAL(KIND=dp),       INTENT(IN)  :: RHOCGS
  CHARACTER(LEN=*),    INTENT(IN)  :: EOS_PATH
  INTEGER(KIND=i4),    INTENT(IN)  :: NPOINTS
  TYPE(TOV_PROFILE_T), INTENT(OUT) :: PROFILE

  REAL(KIND=dp), PARAMETER :: G = 6.67_dp
  REAL(KIND=dp), PARAMETER :: C = 2.99792458_dp
  ! nstot.f's own comments label these "c**2/G=1.347459039d28 g/cm" and
  ! "c**4/G=1.211035789d49 g*cm/s**2" (the true physical values), but the
  ! CODED literals are 1.347459039d0/1.211035789d0 -- no exponent. Ported
  ! exactly as coded (not as commented): the reference output was
  ! produced by the coded value, and other conversion factors throughout
  ! this file (1.78d12, 1.602d33, etc.) are calibrated to match it.
  REAL(KIND=dp), PARAMETER :: C2DG = 1.347459039E0_dp
  REAL(KIND=dp), PARAMETER :: C4DG = 1.211035789E0_dp
  REAL(KIND=dp) :: PI, UL, RHOC, RHOCGS_SCALE
  REAL(KIND=dp) :: P, RHO, RHO0, Z, A, XN
  REAL(KIND=dp) :: Y(5), HH, X1, X2, H1, DX
  INTEGER(KIND=i4) :: NOK, NBAD, I
  REAL(KIND=dp) :: NU0, EPHI, ELAM
  REAL(KIND=dp) :: RAD, PCGS_I, RHOCGS_I, NBFM_I
  REAL(KIND=dp) :: AH_I, ZH_I, XH_I, YE_I, YN_I

  PI = ACOS(-1.0_dp)
  ! nstot.f's own `rhocgs` variable, at the point it feeds UL/RHOC, is
  ! the central density in units of 1e14 g/cm**3 (its own comment:
  ! "give the central density in 10**14 g/cm3", e.g. `rhocgs=9.88` for
  ! the M=1.40 star -- NOT 9.88d14). This SOLVE_TOV_STAR's own RHOCGS
  ! argument is true cgs (matching every other caller in this codebase,
  ! e.g. FILL_ROW's direct use of RHOCGS below for the center row), so
  ! the 1e14 rescaling has to happen here, once, feeding only UL/RHOC --
  ! everything downstream (EOS table conversion, profile rows) already
  ! reconstructs true cgs values through UL/C4DG on its own and must NOT
  ! be rescaled again.
  RHOCGS_SCALE = RHOCGS / 1.0E14_dp
  UL = 1.0_dp / SQRT(RHOCGS_SCALE / C2DG)   ! length unit, 1e7 cm, per rhocgs

  CALL LOAD_EOS_TABLE(EOS_PATH, UL, C4DG, ACTIVE_TABLE)
  RHOC = RHOCGS_SCALE / C2DG * UL*UL   ! central density in code units (=1 by construction)

  ! ---- First integration: get R, M, central redshift -------------------
  CALL EOS_AT_DENSITY(ACTIVE_TABLE, RHOC, P, RHO0, Z, A, XN)
  RHO = RHOC

  HH = 1.0E-9_dp * P
  X1 = LOG(P - HH)
  Y(1) = SQRT((2.0_dp*HH) / ((RHOC+P)*(RHOC/3.0_dp+P)*4.0_dp*PI))
  Y(2) = 4.0_dp*PI*RHOC*Y(1)**3 / 3.0_dp
  Y(3) = 4.0_dp*PI*RHOC*Y(1)**3 / 3.0_dp
  Y(4) = 0.0_dp
  Y(5) = Y(1)

  X2 = LOG(1.0E-22_dp * P)
  H1 = 1.0E-2_dp * (X2 - X1)
  CALL ODEINT(Y, X1, X2, H1, NOK, NBAD, DERIVS_TOV)

  NU0 = -Y(4) + LOG(1.0_dp - 2.0_dp*Y(2)/Y(1))
  PROFILE%RADIUS_KM = Y(1)*UL*100.0_dp
  PROFILE%MASS_MSUN = Y(2)*UL*100.0_dp / 1.4766_dp

  ! ---- Second integration: full radial profile --------------------------
  CALL EOS_AT_DENSITY(ACTIVE_TABLE, RHOC, P, RHO0, Z, A, XN)
  RHO = RHOC
  HH = 1.0E-9_dp * P
  X1 = LOG(P - HH)
  Y(1) = SQRT((2.0_dp*HH) / ((RHOC+P)*(RHOC/3.0_dp+P)*4.0_dp*PI))
  Y(2) = 4.0_dp*PI*RHOC*Y(1)**3 / 3.0_dp
  Y(3) = 4.0_dp*PI*RHOC*Y(1)**3 / 3.0_dp
  Y(4) = NU0
  Y(5) = Y(1)

  X2 = LOG(1.0E-22_dp * P)
  DX = (X2 - X1) / REAL(NPOINTS, KIND=dp)

  PROFILE%N = NPOINTS + 2
  ALLOCATE(PROFILE%R(PROFILE%N), PROFILE%RHOCGS(PROFILE%N), PROFILE%PCGS(PROFILE%N))
  ALLOCATE(PROFILE%NBFM(PROFILE%N), PROFILE%NEL(PROFILE%N))
  ALLOCATE(PROFILE%AH(PROFILE%N), PROFILE%ZH(PROFILE%N), PROFILE%XH(PROFILE%N))
  ALLOCATE(PROFILE%YE(PROFILE%N), PROFILE%YN(PROFILE%N), PROFILE%A_TABLE(PROFILE%N))

  ! Center point (row 1)
  CALL COMPOSITION(A, Z, XN, AH_I, ZH_I, XH_I, YE_I, YN_I)
  NBFM_I = RHO0 * C2DG / UL/UL / 1.673E1_dp
  CALL FILL_ROW(PROFILE, 1, 0.0_dp, RHOCGS, P*C4DG/UL/UL*1.0E35_dp, &
    NBFM_I, A, AH_I, ZH_I, XH_I, YE_I, YN_I)

  DO I = 1, NPOINTS + 1
    X2 = X1 + DX
    H1 = 1.0E-2_dp * (X2 - X1)
    CALL ODEINT(Y, X1, X2, H1, NOK, NBAD, DERIVS_TOV)
    X1 = X2

    P = EXP(X2)
    CALL EOS_AT_PRESSURE(ACTIVE_TABLE, P, RHO, RHO0, Z, A, XN)

    RAD      = Y(1)*UL*100.0_dp
    PCGS_I   = P*C4DG/UL/UL*1.0E35_dp
    RHOCGS_I = RHO*C4DG/UL/UL*1.0E2_dp*1.78E12_dp/1.602_dp
    NBFM_I   = RHO0*C2DG/UL/UL/1.673E1_dp

    CALL COMPOSITION(A, Z, XN, AH_I, ZH_I, XH_I, YE_I, YN_I)
    CALL FILL_ROW(PROFILE, I+1, RAD, RHOCGS_I, PCGS_I, NBFM_I, &
      A, AH_I, ZH_I, XH_I, YE_I, YN_I)
  END DO

  EPHI = EXP(Y(4)/2.0_dp)   ! computed for parity with nstot.f; not yet
  ELAM = 1.0_dp/SQRT(1.0_dp - 2.0_dp*Y(2)/Y(1))   ! surfaced on PROFILE
  ASSOCIATE (UNUSED_EPHI => EPHI, UNUSED_ELAM => ELAM); END ASSOCIATE
END SUBROUTINE SOLVE_TOV_STAR

!> Composition bookkeeping, ported from nstot.f's repeated `IF (a.gt.0)`
!> block: AH<=0 (homogeneous matter, e.g. deep crust/core where nuclei
!> have dissolved) vs. a genuine (Z,A) nucleus.
SUBROUTINE COMPOSITION(A, Z, XN, AH, ZH, XH, YE, YN)
  REAL(KIND=dp), INTENT(IN)  :: A, Z, XN
  REAL(KIND=dp), INTENT(OUT) :: AH, ZH, XH, YE, YN
  IF (A > 0.0_dp) THEN
    AH = A
    ZH = Z
    XH = 1.0_dp - XN
    YE = XH*Z/A
    YN = XN
  ELSE
    AH = 1.0_dp
    YE = 1.0_dp - XN
    XH = 0.0_dp
    ZH = YE
    YN = XN
  END IF
END SUBROUTINE COMPOSITION

!> Writes one row of PROFILE (radius, density, pressure, baryon/electron
!> number density, composition) -- `NEL = NBFM*YE`, ported from nstot.f
!> line 235 (computed unconditionally, same as the original).
SUBROUTINE FILL_ROW(PROFILE, IDX, RAD, RHOCGS_I, PCGS_I, NBFM_I, A_I, AH_I, ZH_I, XH_I, YE_I, YN_I)
  TYPE(TOV_PROFILE_T), INTENT(INOUT) :: PROFILE
  INTEGER(KIND=i4),    INTENT(IN)    :: IDX
  REAL(KIND=dp),       INTENT(IN)    :: RAD, RHOCGS_I, PCGS_I, NBFM_I, A_I
  REAL(KIND=dp),       INTENT(IN)    :: AH_I, ZH_I, XH_I, YE_I, YN_I
  PROFILE%R(IDX)      = RAD
  PROFILE%RHOCGS(IDX) = RHOCGS_I
  PROFILE%PCGS(IDX)   = PCGS_I
  PROFILE%NBFM(IDX)   = NBFM_I
  PROFILE%NEL(IDX)    = NBFM_I * YE_I
  PROFILE%A_TABLE(IDX) = A_I
  PROFILE%AH(IDX)     = AH_I
  PROFILE%ZH(IDX)     = ZH_I
  PROFILE%XH(IDX)     = XH_I
  PROFILE%YE(IDX)     = YE_I
  PROFILE%YN(IDX)      = YN_I
END SUBROUTINE FILL_ROW

!> TOV right-hand side, in log(pressure) as the independent variable X
!> (`p=exp(x)`, a standard stiff-EOS trick since P drops many orders of
!> magnitude from center to surface). Y = (r, m_grav, m_baryon, nu,
!> r_star); ported from derivs.f verbatim. Matches ODE_INTEGRATOR's
!> DERIVS_I interface; reads ACTIVE_TABLE (see module header) since the
!> interface itself has no room for it.
SUBROUTINE DERIVS_TOV(X, Y, DYDX)
  REAL(KIND=dp), INTENT(IN)  :: X
  REAL(KIND=dp), INTENT(IN)  :: Y(:)
  REAL(KIND=dp), INTENT(OUT) :: DYDX(:)
  REAL(KIND=dp) :: P, RHO, RHO0, Z, A, XN, DPDR, R, MG, PI4

  PI4 = 4.0_dp * ACOS(-1.0_dp)
  R  = Y(1)
  MG = Y(2)
  P  = EXP(X)
  CALL EOS_AT_PRESSURE(ACTIVE_TABLE, P, RHO, RHO0, Z, A, XN)
  ASSOCIATE (UNUSED_Z => Z, UNUSED_A => A, UNUSED_XN => XN); END ASSOCIATE

  DPDR = -(RHO+P)*(MG + PI4*R**3*P) / (R*R - 2.0_dp*MG*R)
  DYDX(1) = P/DPDR
  DYDX(2) = PI4*R*R*RHO*(P/DPDR)
  DYDX(3) = PI4*R*R*RHO0/SQRT(1.0_dp - 2.0_dp*MG/R)*(P/DPDR)
  DYDX(4) = -2.0_dp*P/(RHO+P)
  DYDX(5) = 1.0_dp/SQRT(1.0_dp - 2.0_dp*MG/R)*(P/DPDR)
END SUBROUTINE DERIVS_TOV

END MODULE TOV_SOLVER
