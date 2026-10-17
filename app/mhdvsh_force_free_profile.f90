!> One-off diagnostic (2026-10-16, per user request): dumps the radial
!> profile of Phi, Psi, and B_r (at the magnetic pole) for a force-free
!> Bessel-Riccati IC, both INSIDE the star (R_min<=r<=R_max, the actual
!> eigenmode solution) and OUTSIDE it (r>R_max, the vacuum exterior,
!> analytically continued -- not from this project's existing Robin-BC
!> vacuum-matching machinery, which assumes the DECAYING-only special
!> case, but the fully general two-term matching derived below, since
!> 'arbitrary_value' mode's Phi'(R_max) does not in general satisfy the
!> Robin relation).
!>
!> @note B_r is reported "at the magnetic pole" (theta=0): B_r(r,theta) =
!>   l(l+1)/r**2 * Phi_lm(r) * Y_lm(theta,phi), and at theta=0 the l=1,m=0
!>   angular factor is exactly 1 (P_1(cos 0)=1) -- this sidesteps needing
!>   FORTVSH's full VSH_POL_UP/DN pointwise basis-function evaluation
!>   (ROADMAP.md programming task 2, "(Phi,Psi) -> physical B converter",
!>   still not implemented anywhere in this codebase) while still giving
!>   a genuine, unambiguous physical B_r(r) profile. Psi/B_tor are NOT
!>   converted to a physical B_theta/B_phi amplitude here -- doing so
!>   needs the same not-yet-built converter, and an unresolved apparent
!>   normalization mismatch noticed between this project's own
!>   `mhd-vsh-relations.tex` "Reference expressions" (coefficient
!>   `r/(k(k+1))`) and `field_diagnostics.f90`'s stated `R_l=r/sqrt(l(l+1))`
!>   convention (coefficient `sqrt(l(l+1))/r`) was NOT chased down here --
!>   flagged rather than guessed at. Psi itself (the potential, not a
!>   claimed B_theta/B_phi magnitude) is still reported directly: it is
!>   exactly zero outside the star by construction (vacuum carries no
!>   toroidal field), which is true regardless of that unresolved
!>   normalization question.
!>
!> Exterior vacuum continuation (general two-term matching, not just the
!> decaying-only case): outside the star, Phi satisfies the mu=0 (current
!> free) limit of the same Riccati-Bessel ODE, `Phi'' - l(l+1)/r**2*Phi=0`,
!> general solution `Phi_ext(r) = A*r**(l+1) + B*r**(-l)`. Matching
!> Phi,Phi' continuously at r=R_max (P0=Phi(R_max), P1=Phi'(R_max), both
!> known from the interior solution) gives, by direct solution of the
!> resulting 2x2 linear system:
!>   A = (l*P0 + R_max*P1) / ((2l+1)*R_max**(l+1))
!>   B = R_max**l * ((l+1)*P0 - R_max*P1) / (2l+1)
!> Cross-checked (2026-10-16): this reduces EXACTLY to the already-tested
!> Robin condition's own pure-decaying case (A=0) when P1=-(l/R_max)*P0,
!> confirmed by direct substitution.
PROGRAM MHDVSH_FORCE_FREE_PROFILE
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T
USE SPHERICAL_BESSEL,   ONLY: SPHERICAL_BESSEL_JN_YN
USE FORCE_FREE_BESSEL,  ONLY: SOLVE_FORCE_FREE_EIGENVALUE, FF_OUTER_BC_ARBITRARY_VALUE
USE INITIAL_CONDITIONS, ONLY: IC_MODE_T, FF_MODE_T, READ_IC_NAMELIST, BUILD_INITIAL_CONDITION
USE TOV_SOLVER,         ONLY: TOV_PROFILE_T, SOLVE_TOV_STAR
USE CRUST_CONDUCTIVITY, ONLY: FIND_TRUNCATION_RADIUS
IMPLICIT NONE

INTEGER(KIND=i4), PARAMETER :: N_R_INTERIOR = 400_i4
INTEGER(KIND=i4), PARAMETER :: N_R_EXTERIOR = 400_i4
INTEGER(KIND=i4), PARAMETER :: LMAX = 1_i4
REAL(KIND=dp),    PARAMETER :: R_EXT_FACTOR = 4.0_dp   ! plot out to 4*R_max
REAL(KIND=dp),    PARAMETER :: NBAR_CENTRAL = 0.5447307_dp
CHARACTER(LEN=*), PARAMETER :: EOS_PATH = 'data/eos/APR_EOS_Cat.dat'
INTEGER(KIND=i4), PARAMETER :: NPOINTS_TOV = 414
REAL(KIND=dp),    PARAMETER :: RHOL_CGS = 2.2E14_dp
REAL(KIND=dp),    PARAMETER :: THETA_MAX_TRUNCATE = 1.3_dp
REAL(KIND=dp),    PARAMETER :: T_KELVIN = 1.0E9_dp

TYPE(TOV_PROFILE_T)     :: PROFILE
TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI
CHARACTER(LEN=64), ALLOCATABLE :: IC_TAGS(:)
TYPE(IC_MODE_T),   ALLOCATABLE :: MODES(:)
TYPE(FF_MODE_T),   ALLOCATABLE :: FF_MODES(:)
CHARACTER(LEN=1024) :: NML_PATH, DATA_PATH, ARG

REAL(KIND=dp) :: R_MIN, R_MAX, MU, A_COEF, B_COEF
REAL(KIND=dp) :: P0, P1, A_EXT, B_EXT
REAL(KIND=dp), ALLOCATABLE :: J(:), Y(:), DJ(:), DY(:)
REAL(KIND=dp) :: R, X, PHI_R, PSI_R, B_R_POLE
INTEGER(KIND=i4) :: IR, IDX, DATA_UNIT

IF (COMMAND_ARGUMENT_COUNT() < 2) THEN
  WRITE(*,'(A)') 'Usage: mhdvsh_force_free_profile <ic_namelist_path> <data_path>'
  STOP 1
END IF
CALL GET_COMMAND_ARGUMENT(1, ARG); NML_PATH = TRIM(ARG)
CALL GET_COMMAND_ARGUMENT(2, ARG); DATA_PATH = TRIM(ARG)

! ---- Live TOV solve, same setup as mhdvsh_hall_crust_profile.f90 ----
CALL SOLVE_TOV_STAR(NBAR_CENTRAL, EOS_PATH, NPOINTS_TOV, PROFILE)
R_MIN = MINVAL(PROFILE%R, MASK=PROFILE%RHOCGS <= RHOL_CGS)
CALL FIND_TRUNCATION_RADIUS(PROFILE, T_KELVIN, R_MAX, THETA_MAX=THETA_MAX_TRUNCATE)

! ---- Interior: build the IC via the real INITIAL_CONDITIONS pipeline ----
CALL BUILD_RADIAL_GRID(RGRID, N_R_INTERIOR, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL READ_IC_NAMELIST(TRIM(NML_PATH), IC_TAGS, MODES, FF_MODES)
CALL BUILD_INITIAL_CONDITION(IC_TAGS, MODES, FF_MODES, RGRID, LMAX, PHI, PSI)
IDX = YLM_INDEX(1_i4, 0_i4)

! Also solve directly for (mu,A_COEF,B_COEF) -- needed for the analytic
! exterior continuation below (not recoverable from PHI%COEF alone).
CALL SOLVE_FORCE_FREE_EIGENVALUE(1_i4, R_MIN, R_MAX, FF_OUTER_BC_ARBITRARY_VALUE, &
  FF_MODES(1)%ROOT_INDEX, MU, A_COEF, B_COEF)
A_COEF = FF_MODES(1)%AMPLITUDE * A_COEF
B_COEF = FF_MODES(1)%AMPLITUDE * B_COEF

OPEN(NEWUNIT=DATA_UNIT, FILE=TRIM(DATA_PATH), STATUS='REPLACE', ACTION='WRITE')
WRITE(DATA_UNIT,'(A)') '# MHDVSH_FORCE_FREE_PROFILE: l=1 force-free Bessel-Riccati IC'
WRITE(DATA_UNIT,'(A,A)') '#   ic_namelist = ', TRIM(NML_PATH)
WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   R_min = ', R_MIN, ' km (live TOV, crust/core boundary)'
WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   R_max = ', R_MAX, ' km (live TOV, theta-truncated surface)'
WRITE(DATA_UNIT,'(A,ES16.8,A)') '#   mu    = ', MU, ' km**-1 (solved eigenvalue)'
WRITE(DATA_UNIT,'(A)') '#   B_r reported AT THE MAGNETIC POLE (theta=0, l=1 angular factor =1)'
WRITE(DATA_UNIT,'(A)') '#   units: B in 10**12 G, length in km (src/core/units.f90)'
WRITE(DATA_UNIT,'(A)') '# region  r_km  phi  psi  b_r_1e12G'

ALLOCATE(J(0:1), Y(0:1), DJ(0:1), DY(0:1))
DO IR = 1, RGRID%N
  R = RGRID%R(IR)
  PHI_R = REAL(PHI%COEF(IR,IDX), KIND=dp)
  PSI_R = REAL(PSI%COEF(IR,IDX), KIND=dp)
  B_R_POLE = 2.0_dp/R**2 * PHI_R   ! l(l+1)=2 for l=1
  WRITE(DATA_UNIT,'(A,4ES18.10)') 'interior ', R, PHI_R, PSI_R, B_R_POLE
END DO

! ---- Exterior vacuum continuation (see module header for the derivation) ----
R = R_MAX
X = MU*R
CALL SPHERICAL_BESSEL_JN_YN(X, 1_i4, J, Y, DJ, DY)
P0 = A_COEF*X*J(1) + B_COEF*X*Y(1)
! d/dr[r*j_1(mu*r)] = j_1(x) + r*mu*DJ(x), x=mu*r; same for the y_1 term.
P1 = A_COEF*(J(1) + R*MU*DJ(1)) + B_COEF*(Y(1) + R*MU*DY(1))

A_EXT = (1.0_dp*P0 + R_MAX*P1) / (3.0_dp*R_MAX**2)          ! l=1: (l*P0+R*P1)/((2l+1)*R**(l+1))
B_EXT = R_MAX * (2.0_dp*P0 - R_MAX*P1) / 3.0_dp              ! l=1: R**l*((l+1)*P0-R*P1)/(2l+1)

WRITE(DATA_UNIT,'(A,ES16.8)') '#   exterior match check: P0 (interior) = ', P0
WRITE(DATA_UNIT,'(A,ES16.8)') '#   exterior match check: A_EXT*R_max+B_EXT/R_max**1 = ', &
  A_EXT*R_MAX**2 + B_EXT/R_MAX
WRITE(DATA_UNIT,'(A,ES16.8,A,ES16.8)') '#   A_ext = ', A_EXT, '   B_ext = ', B_EXT

DO IR = 1, N_R_EXTERIOR
  R = R_MAX + (R_EXT_FACTOR*R_MAX - R_MAX) * REAL(IR,KIND=dp)/REAL(N_R_EXTERIOR,KIND=dp)
  PHI_R = A_EXT*R**2 + B_EXT/R
  B_R_POLE = 2.0_dp/R**2 * PHI_R
  WRITE(DATA_UNIT,'(A,4ES18.10)') 'exterior ', R, PHI_R, 0.0_dp, B_R_POLE
END DO

CLOSE(DATA_UNIT)
WRITE(*,'(A)') 'MHDVSH_FORCE_FREE_PROFILE'
WRITE(*,'(A,ES12.5,A,ES12.5,A)') '  crust extent = ', R_MIN, ' -- ', R_MAX, ' km'
WRITE(*,'(A,ES12.5)') '  mu (eigenvalue) = ', MU
WRITE(*,'(A,ES12.5,A)') '  Phi(R_max) = ', P0, ' (should equal AMPLITUDE)'
WRITE(*,'(A,ES12.5)') '  B_r(pole,R_max) [1e12 G] = ', 2.0_dp/R_MAX**2*P0
WRITE(*,'(A,A)') '  data written to ', TRIM(DATA_PATH)
END PROGRAM MHDVSH_FORCE_FREE_PROFILE
