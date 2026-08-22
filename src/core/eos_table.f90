MODULE EOS_TABLE
!> Tabulated crust/low-density equation-of-state loading and
!> interpolation -- ported from ~/Desktop/EOSNS/src/geteost.f (its
!> `ENTRY INITEOSTAB` for loading, both interpolation branches for
!> lookup), replacing its `COMMON /eos/` file-scope state with an
!> explicit EOS_TABLE_T passed by the caller (this codebase's own
!> convention -- RADIAL_GRID_T/RADIAL_OPERATOR_T are the precedent, not
!> module-level SAVE'd globals, except where a regime's own init/advance
!> split requires it as DIFFUSION_REGIME does).
!>
!> Table file format (unchanged from `lowd-eos.ja.tab`/
!> `lowd-eos.ja.apr.tab`): one row per row, six whitespace-separated
!> columns -- baryon number density (fm**-3), mass density (g/cm**3),
!> pressure (dyn/cm**2 x1e-33... see LOAD_EOS_TABLE's own conversion,
!> ported unchanged from the original), Z, A, neutron fraction.
!>
!> @warning Geometrized-unit conversion constants (the `1.602d0`,
!>   `1.78d12`, `1.673d15` factors below) are copied byte-for-byte from
!>   geteost.f -- these are the units bridge to nstot.f's `G=c=1` TOV
!>   integration (TOV_SOLVER), not something to re-derive or "clean up".
USE KINDS, ONLY: dp, i4
USE ODE_INTEGRATOR, ONLY: LOCATE_TABLE
IMPLICIT NONE
PRIVATE
PUBLIC :: EOS_TABLE_T, LOAD_EOS_TABLE, EOS_AT_DENSITY, EOS_AT_PRESSURE

!> One loaded, geometrized-unit-converted EOS table (see LOAD_EOS_TABLE).
!> PT/RHOT are monotonic in the same sense (both increasing with row
!> index), which LOCATE_TABLE relies on.
TYPE :: EOS_TABLE_T
  INTEGER(KIND=i4) :: N = 0
  REAL(KIND=dp), ALLOCATABLE :: PT(:), RHOT(:), RHO0T(:), ZT(:), AT(:), XNT(:)
END TYPE EOS_TABLE_T

CONTAINS

!> Loads an EOS table file (e.g. `lowd-eos.ja.tab` or the APR variant
!> `lowd-eos.ja.apr.tab`) and converts it into the geometrized (G=c=1)
!> units TOV_SOLVER integrates in, given the same central-density
!> length scale UL (`= 1/sqrt(rhocgs/c2dg)`, TOV_SOLVER's own choice)
!> and c4dg used there. Ported from geteost.f's `ENTRY INITEOSTAB`,
!> generalized to take the file path and unit-scale arguments rather
!> than a hardcoded filename and module-global units -- this is what
!> makes the APR variant (present in EOSNS but never wired to load
!> there) trivially selectable, and what lets TOV_SOLVER own its own
!> central-density-dependent length scale rather than this module
!> needing to know about it.
!>
!> @param PATH EOS table file path.
!> @param UL Length unit (central-density-dependent, from TOV_SOLVER).
!> @param C4DG `c**4/G` in the same cgs-derived units nstot.f uses.
SUBROUTINE LOAD_EOS_TABLE(PATH, UL, C4DG, TABLE)
  CHARACTER(LEN=*),   INTENT(IN)  :: PATH
  REAL(KIND=dp),      INTENT(IN)  :: UL, C4DG
  TYPE(EOS_TABLE_T),  INTENT(OUT) :: TABLE
  INTEGER(KIND=i4), PARAMETER :: MAXROWS = 2000
  INTEGER(KIND=i4), PARAMETER :: UNIT = 77
  REAL(KIND=dp) :: X1, X2, X3, X4, X5, X6
  REAL(KIND=dp) :: PT(MAXROWS), RHOT(MAXROWS), RHO0T(MAXROWS)
  REAL(KIND=dp) :: ZT(MAXROWS), AT(MAXROWS), XNT(MAXROWS)
  INTEGER(KIND=i4) :: I, IOS

  OPEN(UNIT=UNIT, FILE=TRIM(PATH), STATUS='OLD', ACTION='READ')
  DO I = 1, MAXROWS
    READ(UNIT, *, IOSTAT=IOS) X1, X2, X3, X4, X5, X6
    IF (IOS /= 0) EXIT
    PT(I)    = X3 * 1.602_dp * UL*UL / C4DG / 1.0E2_dp / 1.602E33_dp
    RHOT(I)  = X2 * 1.602_dp * UL*UL / C4DG / 1.0E2_dp / 1.78E12_dp
    RHO0T(I) = (X1 * 1.673E15_dp) * 1.602_dp * UL*UL / C4DG / 1.0E2_dp / 1.78E12_dp
    ZT(I)    = X4
    AT(I)    = X5
    XNT(I)   = X6
  END DO
  CLOSE(UNIT)

  TABLE%N = I - 1
  ALLOCATE(TABLE%PT(TABLE%N), TABLE%RHOT(TABLE%N), TABLE%RHO0T(TABLE%N))
  ALLOCATE(TABLE%ZT(TABLE%N), TABLE%AT(TABLE%N), TABLE%XNT(TABLE%N))
  TABLE%PT    = PT(1:TABLE%N)
  TABLE%RHOT  = RHOT(1:TABLE%N)
  TABLE%RHO0T = RHO0T(1:TABLE%N)
  TABLE%ZT    = ZT(1:TABLE%N)
  TABLE%AT    = AT(1:TABLE%N)
  TABLE%XNT   = XNT(1:TABLE%N)
END SUBROUTINE LOAD_EOS_TABLE

!> Interpolates the EOS at a given (geometrized-unit) mass density RHO,
!> log-log between the two bracketing table rows. Ported from
!> geteost.f's `iflag=0` branch. STOPs if RHO is outside the table
!> range (ported behavior -- the original treats this as fatal, not
!> recoverable, since it means the caller's density grid has run past
!> what the EOS covers).
SUBROUTINE EOS_AT_DENSITY(TABLE, RHO, P, RHO0, Z, A, XN)
  TYPE(EOS_TABLE_T), INTENT(IN)  :: TABLE
  REAL(KIND=dp),     INTENT(IN)  :: RHO
  REAL(KIND=dp),     INTENT(OUT) :: P, RHO0, Z, A, XN
  INTEGER(KIND=i4) :: J
  REAL(KIND=dp) :: ALFA

  CALL LOCATE_TABLE(TABLE%RHOT, RHO, J)
  IF (J == 0 .OR. J == TABLE%N) THEN
    WRITE(*,'(A)') 'EOS_AT_DENSITY: rho is out of the table'
    STOP 1
  END IF

  ALFA = LOG10(RHO/TABLE%RHOT(J)) / LOG10(TABLE%RHOT(J+1)/TABLE%RHOT(J))
  P    = TABLE%PT(J) * 10.0_dp**(LOG10(TABLE%PT(J+1)/TABLE%PT(J)) / &
         LOG10(TABLE%RHOT(J+1)/TABLE%RHOT(J)) * LOG10(RHO/TABLE%RHOT(J)))
  RHO0 = TABLE%RHO0T(J) * 10.0_dp**(LOG10(TABLE%RHO0T(J+1)/TABLE%RHO0T(J)) / &
         LOG10(TABLE%RHOT(J+1)/TABLE%RHOT(J)) * LOG10(RHO/TABLE%RHOT(J)))
  Z  = TABLE%ZT(J)  + ALFA*(TABLE%ZT(J+1)  - TABLE%ZT(J))
  A  = TABLE%AT(J)  + ALFA*(TABLE%AT(J+1)  - TABLE%AT(J))
  XN = TABLE%XNT(J) + ALFA*(TABLE%XNT(J+1) - TABLE%XNT(J))
END SUBROUTINE EOS_AT_DENSITY

!> Interpolates the EOS at a given (geometrized-unit) pressure P,
!> log-log between the two bracketing table rows, returning the
!> corresponding density RHO plus RHO0/Z/A/XN. Ported from geteost.f's
!> `iflag/=0` branch -- this is the one TOV_SOLVER's DERIVS callback
!> uses every integration step (pressure is the independent variable).
SUBROUTINE EOS_AT_PRESSURE(TABLE, P, RHO, RHO0, Z, A, XN)
  TYPE(EOS_TABLE_T), INTENT(IN)  :: TABLE
  REAL(KIND=dp),     INTENT(IN)  :: P
  REAL(KIND=dp),     INTENT(OUT) :: RHO, RHO0, Z, A, XN
  INTEGER(KIND=i4) :: J
  REAL(KIND=dp) :: ALFA

  CALL LOCATE_TABLE(TABLE%PT, P, J)
  IF (J == 0 .OR. J == TABLE%N) THEN
    WRITE(*,'(A)') 'EOS_AT_PRESSURE: p is out of the table'
    STOP 1
  END IF

  ALFA = LOG10(P/TABLE%PT(J)) / LOG10(TABLE%PT(J+1)/TABLE%PT(J))
  RHO  = TABLE%RHOT(J) * 10.0_dp**(LOG10(TABLE%RHOT(J+1)/TABLE%RHOT(J)) / &
         LOG10(TABLE%PT(J+1)/TABLE%PT(J)) * LOG10(P/TABLE%PT(J)))
  RHO0 = TABLE%RHO0T(J) * 10.0_dp**(LOG10(TABLE%RHO0T(J+1)/TABLE%RHO0T(J)) / &
         LOG10(TABLE%PT(J+1)/TABLE%PT(J)) * LOG10(P/TABLE%PT(J)))
  Z  = TABLE%ZT(J)  + ALFA*(TABLE%ZT(J+1)  - TABLE%ZT(J))
  A  = TABLE%AT(J)  + ALFA*(TABLE%AT(J+1)  - TABLE%AT(J))
  XN = TABLE%XNT(J) + ALFA*(TABLE%XNT(J+1) - TABLE%XNT(J))
END SUBROUTINE EOS_AT_PRESSURE

END MODULE EOS_TABLE
