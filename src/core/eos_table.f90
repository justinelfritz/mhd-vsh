MODULE EOS_TABLE
!> Tabulated dense-matter equation-of-state loading and interpolation --
!> ported from Dany Page's NSCool (`astroscu.unam.mx/neutrones/NSCool`,
!> ASCL entry 1609.009), `TOV/TOV.f`'s `subroutine eos` (table loading)
!> and `function ener`/`function pres`/`function rho` (interpolation),
!> replacing NSCool's own `common/eos_dat/` file-scope arrays with an
!> explicit EOS_TABLE_T passed by the caller (this codebase's own
!> convention).
!>
!> @warning Table file format: a fixed 6 header lines (`itext=6`,
!>   coded literally in `eos()`, not read from the file's own header
!>   row even though that row also happens to say 6 -- ported as coded,
!>   matching this session's "port what's coded, not what's commented"
!>   discipline), then rows of at least 3 whitespace-separated columns.
!>   `eos()` auto-detects column order via a magnitude heuristic (see
!>   LOAD_EOS_TABLE) rather than assuming a fixed order -- both branches
!>   ported, even though the one EOS table shipped with this project
!>   (`data/eos/APR_EOS_Cat.dat`) only exercises one of them.
!>
!> @warning Units: table values are read as literal cgs (`g/cm**3`,
!>   `dyn/cm**2`, `fm**-3`) and converted to this port's own code units via
!>   CONST/C2DG (see TOV_SOLVER's own header for the full unit-system
!>   derivation -- length: km, mass: solar masses, density: solar
!>   masses per `km**3`). These are NSCool's own 'cgs' input-mode
!>   constants (`const=1/1.989d18`, `c2=9.d20`, TOV.f's `cread('cgs')`
!>   branch) -- the only mode this port supports, matching what
!>   produced the bundled reference output this port is regression-
!>   tested against (`~NSCool/TOV/Profile/Prof_APR_Cat_1.4.dat`).
USE KINDS, ONLY: dp, i4
USE ODE_INTEGRATOR, ONLY: LOCATE_TABLE
IMPLICIT NONE
PRIVATE
PUBLIC :: EOS_TABLE_T, LOAD_EOS_TABLE, EOS_AT_DENSITY, EOS_AT_PRESSURE

!> One loaded, code-unit-converted EOS table. PT/RHOT are monotonic in
!> the same sense (increasing with row index); NBART is baryon number
!> density (`fm**-3`), unconverted (NSCool's own `eos()` never scales the
!> `deos`/nbar column -- only pressure and mass density get the
!> CONST/C2DG treatment). LOCATE_TABLE relies on PT's monotonicity.
TYPE :: EOS_TABLE_T
  INTEGER(KIND=i4) :: N = 0
  REAL(KIND=dp), ALLOCATABLE :: PT(:), RHOT(:), NBART(:)
END TYPE EOS_TABLE_T

CONTAINS

!> Loads an EOS table file (e.g. `data/eos/APR_EOS_Cat.dat`) and
!> converts it into this port's code units. Ported from NSCool's
!> `eos()`: skips a fixed 6 header lines (`itext=6`, coded literally),
!> reads up to 1000 rows (`limit=1000`, also coded literally --
!> APR_EOS_Cat.dat's own header row says 241 rows, but `eos()` never
!> reads that value, matching the same "coded, not commented"
!> transcription discipline as `TOV_SOLVER`), reading only the first 3
!> whitespace-separated columns of each row (the table's own additional
!> composition/species columns, present in APR_EOS_Cat.dat, are never
!> read by `eos()` -- composition here comes from CRUST_CONDUCTIVITY's
!> `OYAFORM` instead, not this table).
!>
!> Column-order auto-detection (ported from `eos()`'s own `ilist`
!> logic, checked once on the first data row): if column 3 is small
!> (<=10, i.e. a baryon density in `fm**-3`) and column 2 is large
!> (>=1e30, i.e. a pressure in `dyn/cm**2`), the order is (rho, P, nbar)
!> [`ilist=1`, APR_EOS_Cat.dat's own order]; if column 1 is small and
!> column 3 is large, the order is (nbar, rho, P) [`ilist=2`]; anything
!> else is a fatal error (matching `eos()`'s own `pause`).
!>
!> @param PATH EOS table file path.
SUBROUTINE LOAD_EOS_TABLE(PATH, TABLE)
  CHARACTER(LEN=*),  INTENT(IN)  :: PATH
  TYPE(EOS_TABLE_T), INTENT(OUT) :: TABLE
  INTEGER(KIND=i4), PARAMETER :: ITEXT = 6, MAXROWS = 1000
  INTEGER(KIND=i4), PARAMETER :: UNIT = 78
  REAL(KIND=dp), PARAMETER :: CONST = 1.0_dp/1.989E18_dp, C2DG = 9.0E20_dp
  REAL(KIND=dp) :: X1, X2, X3
  REAL(KIND=dp) :: PT(MAXROWS), RHOT(MAXROWS), NBART(MAXROWS)
  INTEGER(KIND=i4) :: I, IOS, ILIST

  OPEN(UNIT=UNIT, FILE=TRIM(PATH), STATUS='OLD', ACTION='READ')
  DO I = 1, ITEXT
    READ(UNIT, *)
  END DO

  ILIST = 0
  DO I = 1, MAXROWS
    READ(UNIT, *, IOSTAT=IOS) X1, X2, X3
    IF (IOS /= 0) EXIT
    IF (I == 1) THEN
      IF (X3 <= 10.0_dp .AND. X2 >= 1.0E30_dp) THEN
        ILIST = 1
      ELSE IF (X1 <= 10.0_dp .AND. X3 >= 1.0E30_dp) THEN
        ILIST = 2
      ELSE
        WRITE(*,'(A)') 'LOAD_EOS_TABLE: cannot determine EOS column ordering'
        STOP 1
      END IF
    END IF
    IF (ILIST == 1) THEN
      RHOT(I) = X1;  PT(I) = X2;  NBART(I) = X3
    ELSE
      RHOT(I) = X2;  PT(I) = X3;  NBART(I) = X1
    END IF
    PT(I)   = PT(I)   * CONST / C2DG
    RHOT(I) = RHOT(I) * CONST
  END DO
  CLOSE(UNIT)

  TABLE%N = I - 1
  ALLOCATE(TABLE%PT(TABLE%N), TABLE%RHOT(TABLE%N), TABLE%NBART(TABLE%N))
  TABLE%PT    = PT(1:TABLE%N)
  TABLE%RHOT  = RHOT(1:TABLE%N)
  TABLE%NBART = NBART(1:TABLE%N)
END SUBROUTINE LOAD_EOS_TABLE

!> Interpolates the EOS at a given (code-unit) baryon number density
!> NBAR, log-log between the two bracketing table rows -- ported from
!> `function pres(d)` (`d` there is NSCool's own `deos`/nbar array).
!> STOPs if NBAR is outside the table range (`pres`'s own `i2.gt.limit`
!> guard, ported as fatal here too).
!>
!> @param NBAR Baryon number density, `fm**-3` (unconverted, see module
!>   header -- NOT a mass density, despite the name matching
!>   TOV_SOLVER's usage: this is exactly what NSCool's own `rhoc` input
!>   parameter means -- `den(0)=rhoc` in TOV.f's main loop).
SUBROUTINE EOS_AT_DENSITY(TABLE, NBAR, P)
  TYPE(EOS_TABLE_T), INTENT(IN)  :: TABLE
  REAL(KIND=dp),     INTENT(IN)  :: NBAR
  REAL(KIND=dp),     INTENT(OUT) :: P
  INTEGER(KIND=i4) :: J

  CALL LOCATE_TABLE(TABLE%NBART, NBAR, J)
  IF (J == 0 .OR. J == TABLE%N) THEN
    WRITE(*,'(A)') 'EOS_AT_DENSITY: nbar is out of the table'
    STOP 1
  END IF

  P = TABLE%PT(J) * 10.0_dp**(LOG10(TABLE%PT(J+1)/TABLE%PT(J)) / &
      LOG10(TABLE%NBART(J+1)/TABLE%NBART(J)) * LOG10(NBAR/TABLE%NBART(J)))
END SUBROUTINE EOS_AT_DENSITY

!> Interpolates the EOS at a given (code-unit) pressure P, log-log
!> between the two bracketing table rows, returning mass density RHO
!> and baryon number density NBAR -- ported from `function ener(p)`
!> (-> RHO) and `function rho(p)` (-> NBAR; confusingly named in the
!> original -- it returns baryon density, not mass density, since it
!> interpolates against `deos`). This is the pair TOV_SOLVER's RHS
!> callback uses every integration step (pressure is the independent
!> variable, matching `ener`/`rho`'s own usage in `twostep`).
SUBROUTINE EOS_AT_PRESSURE(TABLE, P, RHO, NBAR)
  TYPE(EOS_TABLE_T), INTENT(IN)  :: TABLE
  REAL(KIND=dp),     INTENT(IN)  :: P
  REAL(KIND=dp),     INTENT(OUT) :: RHO, NBAR
  INTEGER(KIND=i4) :: J

  CALL LOCATE_TABLE(TABLE%PT, P, J)
  IF (J == 0 .OR. J == TABLE%N) THEN
    WRITE(*,'(A)') 'EOS_AT_PRESSURE: p is out of the table'
    STOP 1
  END IF

  RHO  = TABLE%RHOT(J) * 10.0_dp**(LOG10(TABLE%RHOT(J+1)/TABLE%RHOT(J)) / &
         LOG10(TABLE%PT(J+1)/TABLE%PT(J)) * LOG10(P/TABLE%PT(J)))
  NBAR = TABLE%NBART(J) * 10.0_dp**(LOG10(TABLE%NBART(J+1)/TABLE%NBART(J)) / &
         LOG10(TABLE%PT(J+1)/TABLE%PT(J)) * LOG10(P/TABLE%PT(J)))
END SUBROUTINE EOS_AT_PRESSURE

END MODULE EOS_TABLE
