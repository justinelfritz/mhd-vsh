!> Regression test for CRUST_CONDUCTIVITY's new CONDUCT-based engine
!> (CONDUCT_TRANSPORT/CONDUCT_CORE and its full dependency chain) against
!> TWO independent standalone harnesses built directly from Alexander
!> Potekhin's own unmodified `conduct21.f`
!> (http://www.ioffe.ru/astro/conduct/, not this port's own code):
!>
!> (1) CASE1-8: `CONDCONV`/`CONDUCT` called with the original's OWN
!>     hardcoded `xnuc=0` (its demo-driver default) at 8 densities
!>     spanning envelope through crust (including the two exact profile
!>     rows, rho=6.17686691e7 and 5.55712691e7, that exposed the old
!>     GYP/PBHY seam) plus two extra temperatures at rho=6e7 -- verifies
!>     this port's transcription of CONDUCT/ThAv18/ThAvI18/COUL19/
!>     COULAN3/CHEMPOT/CHEMP99/FERINV/EXPINT/TAUEESY is faithful to the
!>     original's own formulas, independent of this port's own xnuc-
!>     wiring design choice.
!> (2) CASE9-11: a SECOND harness, built from a COMMON-block-patched
!>     copy of the SAME unmodified `conduct21.f` (only `CONDUCT`'s
!>     internal `xnuc=0.`/`xnuct=xnuc/1.1` lines replaced with reads
!>     from a `common /xnuc_override/` block -- everything else
!>     byte-identical to the original), fed OYAFORM's own real,
!>     continuously density-dependent xnuc/xnuct at 3 densities --
!>     verifies CONDUCT_TRANSPORT's actual production path (this port's
!>     own CON_CRUST calls it with real xnuc/xnuct, not 0) end-to-end.
!>
!> Reference values were produced by compiling and running these two
!> harnesses directly; not hand-derived, not copied from this port's own
!> output. Z/A/xnuc/xnuct inputs for both harnesses came from this
!> port's own (unchanged, previously-verified) OYAFORM.
!>
!> CASE12-13: end-to-end CON_CRUST calls at the two rows that exposed the
!> old GYP/PBHY seam -- confirms the fix directly: SIGMA changes by only
!> a few percent across the old 6e7 g/cm**3 boundary now, not ~44x.
!>
!> @warning TOL is 5e-6, not full double-precision roundoff -- same root
!>   cause as this project's other F77-physics ports: bare literals in
!>   the original lack a `d0` suffix and are parsed as single precision
!>   before widening to double, while this port's `_dp`-suffixed
!>   literals are fully double-precision throughout.
PROGRAM TEST_CRUST_CONDUCTIVITY
USE KINDS,              ONLY: dp, i4
USE CRUST_CONDUCTIVITY, ONLY: OYAFORM, CONDUCT_TRANSPORT, CON_CRUST
IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: TOL = 5.0E-6_dp
INTEGER(KIND=i4) :: N_FAIL
REAL(KIND=dp) :: SIGMA, CKAPPA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH
REAL(KIND=dp) :: N_E, LAMBDA
REAL(KIND=dp) :: SIGMA_A, SIGMA_B

N_FAIL = 0

! ---- Part 1: xnuc=0 (matches the original's own demo-driver default)
! ---- at 8 densities spanning envelope through crust, T=1e9 K.
CALL CHECK_XNUC0("case1", 1.0E3_dp,          26.00617519_dp,    56.00734215_dp,    &
  3.7615773262E19_dp, 1.1230914171E16_dp, N_FAIL)
CALL CHECK_XNUC0("case2", 1.0E7_dp,          26.37184128_dp,    56.84726056_dp,    &
  3.1125947532E20_dp, 7.1862870075E16_dp, N_FAIL)
CALL CHECK_XNUC0("case3", 6.0E7_dp,          27.19082319_dp,    59.23836629_dp,    &
  6.7578316029E20_dp, 1.5672156362E17_dp, N_FAIL)
CALL CHECK_XNUC0("case4", 6.17686691E7_dp,   27.20933148_dp,    59.29758203_dp,    &
  6.8353364099E20_dp, 1.5862433762E17_dp, N_FAIL)
CALL CHECK_XNUC0("case5", 5.55712691E7_dp,   27.14276232_dp,    59.08547370_dp,    &
  6.5564369613E20_dp, 1.5178916743E17_dp, N_FAIL)
CALL CHECK_XNUC0("case6", 1.0E11_dp,         35.93642338_dp,    103.45888165_dp,   &
  1.4570939577E22_dp, 3.6254874186E18_dp, N_FAIL)
CALL CHECK_XNUC0("case7", 1.0E13_dp,         46.29287558_dp,    1163.73123056_dp,  &
  5.2093742991E22_dp, 1.3142969312E19_dp, N_FAIL)
CALL CHECK_XNUC0("case8", 1.22E14_dp,        30.94439645_dp,    1016.81844314_dp,  &
  2.6713909665E23_dp, 6.1815767995E19_dp, N_FAIL)

! ---- Part 2: real OYAFORM xnuc/xnuct wired through CONDUCT_TRANSPORT
! ---- (this port's actual production path), T=1e9 K. Confirms the
! ---- finite-nuclear-size correction (up to ~6.6x in SIGMA at the
! ---- highest density tested) is transcribed correctly, not just the
! ---- xnuc=0 point-nucleus path checked in Part 1.
CALL CONDUCT_TRANSPORT(1.0E9_dp, 6.0E7_dp,   0.0_dp, 27.19082319_dp, 59.23836629_dp,   0.0_dp, &
  0.00642206_dp, 0.00583824_dp, SIGMA, CKAPPA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH)
CALL CHECK_CLOSE("case9_sigma",  SIGMA,  6.7588144646E20_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case9_ckappa", CKAPPA, 1.5674205795E17_dp, TOL, N_FAIL)

CALL CONDUCT_TRANSPORT(1.0E9_dp, 1.0E11_dp,  0.0_dp, 35.93642338_dp, 103.45888165_dp, 0.0_dp, &
  0.07399714_dp, 0.06727013_dp, SIGMA, CKAPPA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH)
CALL CHECK_CLOSE("case10_sigma",  SIGMA,  1.5144741169E22_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case10_ckappa", CKAPPA, 3.7573854336E18_dp, TOL, N_FAIL)

CALL CONDUCT_TRANSPORT(1.0E9_dp, 1.22E14_dp, 0.0_dp, 30.94439645_dp, 1016.81844314_dp, 0.0_dp, &
  0.54701772_dp, 0.48241142_dp, SIGMA, CKAPPA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH)
CALL CHECK_CLOSE("case11_sigma",  SIGMA,  1.7549299257E24_dp, TOL, N_FAIL)
CALL CHECK_CLOSE("case11_ckappa", CKAPPA, 3.0812895028E20_dp, TOL, N_FAIL)

! ---- Part 3: end-to-end CON_CRUST across the old GYP/PBHY seam --
! ---- confirms the fix directly (row-to-row change now ~4%, not ~44x).
CALL CON_CRUST(1.0E9_dp, 6.17686691E7_dp, SIGMA_A, LAMBDA, N_E)
CALL CON_CRUST(1.0E9_dp, 5.55712691E7_dp, SIGMA_B, LAMBDA, N_E)
CALL CHECK_TRUE("case12_seam_no_longer_discontinuous", &
  ABS(SIGMA_A/SIGMA_B - 1.0_dp) < 0.10_dp, N_FAIL)
CALL CHECK_TRUE("case12_seam_sigma_a_positive_sane", SIGMA_A > 1.0E19_dp .AND. SIGMA_A < 1.0E22_dp, N_FAIL)

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_crust_conductivity"
END IF

CONTAINS

SUBROUTINE CHECK_XNUC0(NAME, RHO, Z, A, EXPECTED_SIGMA, EXPECTED_CKAPPA, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: RHO, Z, A, EXPECTED_SIGMA, EXPECTED_CKAPPA
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: SIGMA_L, CKAPPA_L, QJ_L, SIGMAT_L, CKAPPAT_L, QJT_L, SIGMAH_L, CKAPPAH_L, QJH_L
  CALL CONDUCT_TRANSPORT(1.0E9_dp, RHO, 0.0_dp, Z, A, 0.0_dp, 0.0_dp, 0.0_dp, &
    SIGMA_L, CKAPPA_L, QJ_L, SIGMAT_L, CKAPPAT_L, QJT_L, SIGMAH_L, CKAPPAH_L, QJH_L)
  CALL CHECK_CLOSE(TRIM(NAME)//"_sigma",  SIGMA_L,  EXPECTED_SIGMA,  TOL, N_FAIL)
  CALL CHECK_CLOSE(TRIM(NAME)//"_ckappa", CKAPPA_L, EXPECTED_CKAPPA, TOL, N_FAIL)
END SUBROUTINE CHECK_XNUC0

SUBROUTINE CHECK_TRUE(NAME, COND, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  LOGICAL,          INTENT(IN)    :: COND
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  IF (.NOT. COND) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A)') "FAIL  ", NAME
  ELSE
    WRITE(*,'(A,A)') "PASS  ", NAME
  END IF
END SUBROUTINE CHECK_TRUE

SUBROUTINE CHECK_CLOSE(NAME, VAL, EXPECTED, TOL, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: VAL, EXPECTED, TOL
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: ERR
  ERR = ABS(VAL-EXPECTED)/MAX(ABS(EXPECTED), 1.0E-300_dp)
  IF (ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  relerr=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  relerr=", ERR
  END IF
END SUBROUTINE CHECK_CLOSE

END PROGRAM TEST_CRUST_CONDUCTIVITY
