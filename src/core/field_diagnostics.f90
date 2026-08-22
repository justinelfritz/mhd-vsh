MODULE FIELD_DIAGNOSTICS
!> Time-resolved diagnostics computed directly from the (Phi,Psi)
!> potential representation -- never from a reconstructed physical B.
!> Magnetic energy (poloidal/toroidal, and their sum, plus PER-DEGREE-l
!> breakdowns via POLOIDAL_MAGNETIC_ENERGY_BY_L/TOROIDAL_MAGNETIC_ENERGY_BY_L,
!> e.g. for visualizing a mode cascade) is regime-agnostic;
!> the dissipation/flux terms are per-regime Ohm's law: pure-resistive
!> (diffusion, `E=eta*j`) gives JOULE_DISSIPATION_RATE and
!> POYNTING_FLUX_RATE, while the diffusion-free Hall limit
!> (`E=f_H*j x B`) gives a Joule term that is exactly zero (see
!> JOULE_DISSIPATION_RATE's own docstring -- its existing formula covers
!> that regime's resistive part unchanged) and its own, unrelated
!> HALL_POYNTING_FLUX_RATE.
!>
!> Formulas confirmed against Justin Elfritz's own independent derivation
!> (analytic_formulas/mhd-vsh-relations.{tex,pdf}, "Energy Budget in
!> Diffusion Limit" and "...in Hall Limit" sections, 2026-08-11) -- this
!> is no longer the provisional, draft-manuscript-sourced physics an
!> earlier version of this module carried; the energy formula's
!> *structure* (derived here from VSH orthonormality) matches that
!> document's diffusion-limit Eqs.3-4 exactly, up to the unit prefactor
!> now adopted to match it (Gaussian cgs, 1/8pi). The diffusion Joule/
!> Poynting formulas (Eqs.6,8) and the Hall Poynting formula (Hall-limit
!> section) are new; the latter is the first mode-coupled (Gaunt/GWI,
!> `O(Nlm**3)`) formula in this module -- every diffusion-limit formula
!> collapses to a diagonal sum thanks to VSH orthonormality, but Hall
!> Poynting is inherently cubic in the field, so coupling is unavoidable.
!>
!> Derivation sketch: `B_pol(l,m,r) = (1/R_l**2)*[Phi_lm*VSH_POL_DN(l,m) +
!> R_l*Phi_lm'*VSH_POL_UP(l,m)], B_tor(l,m,r) = -(i/R_l)*Psi_lm*
!> VSH_TOR(l,m)`, with R_l=r/sqrt(l(l+1)) (same source as
!> BOUNDARY_CONDITIONS' outer BC). FORTVSH's standard (J-coupled) VSH
!> basis is complete and orthonormal (numerically confirmed against
!> FORTVSH directly, not just asserted -- see the VSH_POL_UP/DN cross-
!> orthogonality check that also makes the Joule formula's Psi term
!> collapse to a clean diagonal sum), so every formula below is a plain
!> sum over (l,m) with no Gaunt/mode-coupling integrals needed.
!>
!> Current currents: `j_pol(l,m,r) = (1/R_l**2)*[Psi_lm*VSH_POL_DN(l,m) +
!> R_l*Psi_lm'*VSH_POL_UP(l,m)], j_tor(l,m,r) = (i/R_l)*(Phi_lm'' -
!> Phi_lm/R_l**2)*VSH_TOR(l,m)` -- i.e. Phi'' (RADIAL_OPERATOR_T%D2) is
!> needed for the Joule/Poynting terms though not for energy alone; note
!> `Phi_lm''-Phi_lm/R_l**2 = Phi_lm'' - l(l+1)*Phi_lm/r**2` is exactly
!> RADIAL_OPERATORS::ADD_CURVATURE_TERM's operator applied to Phi -- the
!> same combination DIFFUSION_REGIME already assembles for its own
!> implicit solve.
!>
!> POYNTING_FLUX_RATE assumes Phi/Psi represent a real physical field
!> (the standard l,-m <-> conj(l,m) symmetry), since its boundary terms
!> aren't manifestly real mode-by-mode the way the `|.|**2`-only energy/
!> Joule sums are -- only the m-sum as a whole is guaranteed real for a
!> real field. The real part is returned; for coefficients that don't
!> satisfy that symmetry (e.g. hand-picked test values), the discarded
!> imaginary part is not itself a sign of an error.
!>
!> HALL_COURANT_TIMESTEP (added 2026-08-21) computes the Hall-CFL-limited
!> maximum stable step size from the current state -- see its own
!> docstring for the exact formula and the RMS-vs-pointwise-max caveat
!> in how "|current|" is estimated without a physical-space B/j
!> reconstruction.
!>
!> @warning The discrete analogue of div(B)=0 (a numerical-consistency
!>   check, since it's an identity by construction in the continuous
!>   theory) is deliberately not implemented yet -- it needs its own
!>   careful derivation in terms of Phi/Psi and their FD derivatives,
!>   not a quick add-on to this module.
!>
!> @note Every formula here is written in genuinely unit-system-agnostic
!>   form: R_MIN/R_MAX/ETA/F_HALL are free REAL(dp) inputs with no unit
!>   baked into this module, so results are dimensionally correct in ANY
!>   consistent choice of length/time/field units, not just Gaussian-cgs
!>   cm/s/Gauss. See UNITS (units.f90) for this project's adopted
!>   code-unit system (B: `10**12` G, length: km, time: yr) and the
!>   multiplicative factors to convert a code-unit result here to erg/
!>   erg-per-second -- that conversion belongs at the reporting layer
!>   (a driver's printed/logged output), not inside this module.
USE KINDS,            ONLY: dp, i4
USE GLOBALS,          ONLY: pi
USE GRID_RADIAL,      ONLY: RADIAL_GRID_T
USE RADIAL_OPERATORS, ONLY: RADIAL_OPERATOR_T, ADD_CURVATURE_TERM
USE FIELD_TYPES,      ONLY: SPECTRAL_SCALAR_T
USE VSH,              ONLY: YLM_INDEX, GWI
IMPLICIT NONE
PRIVATE
PUBLIC :: POLOIDAL_ENERGY_DENSITY, TOROIDAL_ENERGY_DENSITY, &
          TOTAL_POLOIDAL_MAGNETIC_ENERGY, TOTAL_TOROIDAL_MAGNETIC_ENERGY, &
          POLOIDAL_MAGNETIC_ENERGY_BY_L, TOROIDAL_MAGNETIC_ENERGY_BY_L, &
          TOTAL_MAGNETIC_ENERGY, JOULE_DISSIPATION_RATE, POYNTING_FLUX_RATE, &
          HALL_POYNTING_FLUX_RATE, CURRENT_DENSITY_SQUARED_BY_R, HALL_COURANT_TIMESTEP

CONTAINS

!> @param PHI Poloidal potential (see module header), SPECTRAL_SCALAR_T.
!> @param OPS Radial derivative operator (for Phi'); OPS%D1 built on RGRID.
!> @param RGRID Radial grid PHI/OPS were built on.
!> Returns: poloidal magnetic energy density (`1/8pi * B_pol.B_pol*`,
!>   angle-integrated) at every radial node, size (RGRID%N). Integrate in
!>   r for total energy -- see TOTAL_POLOIDAL_MAGNETIC_ENERGY.
FUNCTION POLOIDAL_ENERGY_DENSITY(PHI, OPS, RGRID) RESULT(E_OF_R)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_OF_R(RGRID%N)
  COMPLEX(KIND=dp), ALLOCATABLE :: DPHI_DR(:,:)
  REAL(KIND=dp) :: LAMBDA_L
  INTEGER(KIND=i4) :: L, M, IDX, IR

  ALLOCATE(DPHI_DR(PHI%N_R, PHI%NLM))
  DPHI_DR = MATMUL(OPS%D1, PHI%COEF)

  E_OF_R = 0.0_dp
  DO L = 0, PHI%LMAX
    LAMBDA_L = REAL(L*(L+1), KIND=dp)
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      DO IR = 1, RGRID%N
        E_OF_R(IR) = E_OF_R(IR) + LAMBDA_L * ( &
          LAMBDA_L*ABS(PHI%COEF(IR,IDX))**2/RGRID%R(IR)**2 + ABS(DPHI_DR(IR,IDX))**2)
      END DO
    END DO
  END DO
  E_OF_R = E_OF_R / (8.0_dp*pi)
  DEALLOCATE(DPHI_DR)
END FUNCTION POLOIDAL_ENERGY_DENSITY

!> @param PSI Toroidal potential, SPECTRAL_SCALAR_T.
!> @param RGRID Radial grid PSI was built on.
!> Returns: toroidal magnetic energy density (`1/8pi * B_tor.B_tor*`,
!>   angle-integrated) at every radial node, size (RGRID%N).
FUNCTION TOROIDAL_ENERGY_DENSITY(PSI, RGRID) RESULT(E_OF_R)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PSI
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_OF_R(RGRID%N)
  REAL(KIND=dp) :: LAMBDA_L
  INTEGER(KIND=i4) :: L, M, IDX, IR

  E_OF_R = 0.0_dp
  DO L = 0, PSI%LMAX
    LAMBDA_L = REAL(L*(L+1), KIND=dp)
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      DO IR = 1, RGRID%N
        E_OF_R(IR) = E_OF_R(IR) + LAMBDA_L*ABS(PSI%COEF(IR,IDX))**2
      END DO
    END DO
  END DO
  E_OF_R = E_OF_R / (8.0_dp*pi)
END FUNCTION TOROIDAL_ENERGY_DENSITY

!> Radial (trapezoidal) integral of POLOIDAL_ENERGY_DENSITY.
FUNCTION TOTAL_POLOIDAL_MAGNETIC_ENERGY(PHI, OPS, RGRID) RESULT(E_TOTAL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_TOTAL
  E_TOTAL = RADIAL_INTEGRAL(POLOIDAL_ENERGY_DENSITY(PHI, OPS, RGRID), RGRID)
END FUNCTION TOTAL_POLOIDAL_MAGNETIC_ENERGY

!> @param PHI, OPS, RGRID as POLOIDAL_ENERGY_DENSITY.
!> Returns: total poloidal magnetic energy (mhd-vsh-relations.tex Eq. for
!>   E_B,pol(t), same formula TOTAL_POLOIDAL_MAGNETIC_ENERGY implements)
!>   broken out PER DEGREE l (summed over order m, radially integrated),
!>   size (0:PHI%LMAX) -- SUM(POLOIDAL_MAGNETIC_ENERGY_BY_L(...)) equals
!>   TOTAL_POLOIDAL_MAGNETIC_ENERGY to radial-integration precision (only
!>   the accumulation order differs -- per-l here vs. all-l-at-once
!>   there -- so this is a genuine identity, not an approximation; tested
!>   as such in test_field_diagnostics.f90).
FUNCTION POLOIDAL_MAGNETIC_ENERGY_BY_L(PHI, OPS, RGRID) RESULT(E_OF_L)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_OF_L(0:PHI%LMAX)
  COMPLEX(KIND=dp), ALLOCATABLE :: DPHI_DR(:,:)
  REAL(KIND=dp) :: LAMBDA_L, DENSITY(RGRID%N)
  INTEGER(KIND=i4) :: L, M, IDX, IR

  ALLOCATE(DPHI_DR(PHI%N_R, PHI%NLM))
  DPHI_DR = MATMUL(OPS%D1, PHI%COEF)

  DO L = 0, PHI%LMAX
    LAMBDA_L = REAL(L*(L+1), KIND=dp)
    DENSITY = 0.0_dp
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      DO IR = 1, RGRID%N
        DENSITY(IR) = DENSITY(IR) + LAMBDA_L * ( &
          LAMBDA_L*ABS(PHI%COEF(IR,IDX))**2/RGRID%R(IR)**2 + ABS(DPHI_DR(IR,IDX))**2)
      END DO
    END DO
    DENSITY = DENSITY / (8.0_dp*pi)
    E_OF_L(L) = RADIAL_INTEGRAL(DENSITY, RGRID)
  END DO
  DEALLOCATE(DPHI_DR)
END FUNCTION POLOIDAL_MAGNETIC_ENERGY_BY_L

!> Radial (trapezoidal) integral of TOROIDAL_ENERGY_DENSITY.
FUNCTION TOTAL_TOROIDAL_MAGNETIC_ENERGY(PSI, RGRID) RESULT(E_TOTAL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PSI
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_TOTAL
  E_TOTAL = RADIAL_INTEGRAL(TOROIDAL_ENERGY_DENSITY(PSI, RGRID), RGRID)
END FUNCTION TOTAL_TOROIDAL_MAGNETIC_ENERGY

!> @param PSI, RGRID as TOROIDAL_ENERGY_DENSITY.
!> Returns: total toroidal magnetic energy (mhd-vsh-relations.tex Eq. for
!>   E_B,tor(t), same formula TOTAL_TOROIDAL_MAGNETIC_ENERGY implements)
!>   broken out PER DEGREE l, size (0:PSI%LMAX) -- see
!>   POLOIDAL_MAGNETIC_ENERGY_BY_L's docstring for the same identity
!>   (SUM over l reproduces the `TOTAL_*` function exactly).
FUNCTION TOROIDAL_MAGNETIC_ENERGY_BY_L(PSI, RGRID) RESULT(E_OF_L)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PSI
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_OF_L(0:PSI%LMAX)
  REAL(KIND=dp) :: LAMBDA_L, DENSITY(RGRID%N)
  INTEGER(KIND=i4) :: L, M, IDX, IR

  DO L = 0, PSI%LMAX
    LAMBDA_L = REAL(L*(L+1), KIND=dp)
    DENSITY = 0.0_dp
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      DO IR = 1, RGRID%N
        DENSITY(IR) = DENSITY(IR) + LAMBDA_L*ABS(PSI%COEF(IR,IDX))**2
      END DO
    END DO
    DENSITY = DENSITY / (8.0_dp*pi)
    E_OF_L(L) = RADIAL_INTEGRAL(DENSITY, RGRID)
  END DO
END FUNCTION TOROIDAL_MAGNETIC_ENERGY_BY_L

!> TOTAL_POLOIDAL_MAGNETIC_ENERGY + TOTAL_TOROIDAL_MAGNETIC_ENERGY.
FUNCTION TOTAL_MAGNETIC_ENERGY(PHI, PSI, OPS, RGRID) RESULT(E_TOTAL)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: E_TOTAL
  E_TOTAL = TOTAL_POLOIDAL_MAGNETIC_ENERGY(PHI, OPS, RGRID) + &
            TOTAL_TOROIDAL_MAGNETIC_ENERGY(PSI, RGRID)
END FUNCTION TOTAL_MAGNETIC_ENERGY

!> @param PHI, PSI as above.
!> @param OPS Radial derivative operators (D1 for Psi', D2 for Phi'').
!> @param RGRID Radial grid.
!> @param ETA Resistivity (diffusivity) coefficient.
!> @note Confirmed by the user (2026-08-11) to also be the *complete*
!>   Joule term for the weak-Hall regime, unchanged -- a Hall
!>   contribution to E is always perpendicular to j (j.(j x B)=0
!>   identically), so it does no work on the current and contributes
!>   exactly zero additional dissipation; only the resistive `eta*j` term
!>   ever heats. Does NOT extend to a magnetofrictional Ohm's law this
!>   way -- that term is dissipative by design, a genuinely different
!>   mechanism, not covered by this argument.
!> Returns: volume-integrated Joule dissipation rate, <=0 (a loss).
FUNCTION JOULE_DISSIPATION_RATE(PHI, PSI, OPS, RGRID, ETA) RESULT(E_DOT_J)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp),           INTENT(IN) :: ETA
  REAL(KIND=dp) :: E_DOT_J
  E_DOT_J = -(ETA/(4.0_dp*pi)) * &
    RADIAL_INTEGRAL(CURRENT_DENSITY_SQUARED_BY_R(PHI, PSI, OPS, RGRID), RGRID)
END FUNCTION JOULE_DISSIPATION_RATE

!> Angle-integrated (over the full 4pi sphere) `|curl(B)|**2 = |j|**2` at
!> every radial node -- the j_pol/j_tor combination from the module
!> header (j_pol uses Psi, j_tor uses curvature-of-Phi, dual to how
!> B_pol/B_tor use Phi/Psi), summed via VSH orthonormality (an exact
!> angle integral, not an approximation -- same Parseval argument as
!> every energy formula above). Extracted from what used to be
!> JOULE_DISSIPATION_RATE's own inline loop; shared with
!> HALL_COURANT_TIMESTEP, which uses this per-radial-node rather than
!> integrated over r.
!> @param PHI, PSI, OPS, RGRID as JOULE_DISSIPATION_RATE.
!> Returns: `|curl(B)|**2`, angle-integrated over 4pi sr, at every radial
!>   node, size (RGRID%N) -- no ETA, no radial integral (JOULE_DISSIPATION_RATE
!>   applies both on top of this).
FUNCTION CURRENT_DENSITY_SQUARED_BY_R(PHI, PSI, OPS, RGRID) RESULT(J2_OF_R)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp) :: J2_OF_R(RGRID%N)
  COMPLEX(KIND=dp), ALLOCATABLE :: DPSI_DR(:,:), CURVATURE_PHI(:,:)
  REAL(KIND=dp)    :: LAMBDA_L, D2_CURV(RGRID%N,RGRID%N)
  INTEGER(KIND=i4) :: L, M, IDX, IR

  ALLOCATE(DPSI_DR(PSI%N_R, PSI%NLM), CURVATURE_PHI(PHI%N_R, PHI%NLM))
  DPSI_DR = MATMUL(OPS%D1, PSI%COEF)

  ! Phi'' - l(l+1)/r**2 * Phi, per l -- ADD_CURVATURE_TERM applied to a
  ! per-l copy of D2, the same combination DIFFUSION_REGIME assembles.
  J2_OF_R = 0.0_dp
  DO L = 0, PHI%LMAX
    LAMBDA_L = REAL(L*(L+1), KIND=dp)
    D2_CURV = OPS%D2
    CALL ADD_CURVATURE_TERM(D2_CURV, L, RGRID)
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      CURVATURE_PHI(:,IDX) = MATMUL(D2_CURV, PHI%COEF(:,IDX))
    END DO
    DO M = -L, L
      IDX = YLM_INDEX(L, M)
      DO IR = 1, RGRID%N
        J2_OF_R(IR) = J2_OF_R(IR) + LAMBDA_L * ( &
          LAMBDA_L*ABS(PSI%COEF(IR,IDX))**2/RGRID%R(IR)**2 + ABS(DPSI_DR(IR,IDX))**2 + &
          ABS(CURVATURE_PHI(IR,IDX))**2 )
      END DO
    END DO
  END DO
  DEALLOCATE(DPSI_DR, CURVATURE_PHI)
END FUNCTION CURRENT_DENSITY_SQUARED_BY_R

!> Hall-regime CFL/Courant-limited maximum stable step size (Justin
!> Elfritz, 2026-08-21): tc = min over radial cells of
!> `dr_i/(f_H_i*|current|_i)`, where `|current|=|curl(B)|`. A genuine
!> pointwise |curl(B)|(r,theta,phi) would need a physical-space VSH
!> synthesis this codebase deliberately never does (module header:
!> "never from a reconstructed physical B") -- |current| here is instead
!> the RMS |curl(B)| over the sphere at each radius, an EXACT quantity
!> from VSH orthonormality/Parseval (CURRENT_DENSITY_SQUARED_BY_R is the
!> angle-INTEGRATED `|curl(B)|**2` over 4pi sr; dividing by 4pi and taking
!> sqrt gives the RMS), not an approximation of some pointwise value --
!> and the standard choice for CFL estimation in spectral codes (an
!> L2-norm-based characteristic magnitude).
!>
!> Returns the RAW Courant time, with no safety margin applied --
!> callers (e.g. TIMESTEPPER::RUN_ADAPTIVE) multiply by their own
!> CFL_SAFETY<1 before using it as an actual step size.
!>
!> @param PHI, PSI, OPS, RGRID as JOULE_DISSIPATION_RATE.
!> @param F_HALL Hall prefactor, spatially uniform unless F_HALL_PROFILE
!>   is present (same convention as HALL_INDUCTION_RHS).
!> @param F_HALL_PROFILE Optional per-radial-row override (size RGRID%N).
!> Returns: TC, the minimum Courant time over all radial cells. Cells
!>   where both the local current and F_HALL vanish impose no
!>   constraint (skipped, not treated as TC=0) -- a quiescent field or
!>   F_HALL=0 region is not itself a stability problem.
FUNCTION HALL_COURANT_TIMESTEP(PHI, PSI, OPS, RGRID, F_HALL, F_HALL_PROFILE) RESULT(TC)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp),           INTENT(IN) :: F_HALL
  REAL(KIND=dp), OPTIONAL, INTENT(IN) :: F_HALL_PROFILE(:)
  REAL(KIND=dp) :: TC
  REAL(KIND=dp) :: J_RMS(RGRID%N), F_H(RGRID%N), DR, J_CELL, F_CELL
  REAL(KIND=dp), PARAMETER :: FLOOR = 1.0E-300_dp
  INTEGER(KIND=i4) :: IR

  J_RMS = SQRT(CURRENT_DENSITY_SQUARED_BY_R(PHI, PSI, OPS, RGRID) / (4.0_dp*pi))

  IF (PRESENT(F_HALL_PROFILE)) THEN
    F_H = F_HALL_PROFILE
  ELSE
    F_H = F_HALL
  END IF

  TC = HUGE(1.0_dp)
  DO IR = 1, RGRID%N - 1
    DR = RGRID%R(IR+1) - RGRID%R(IR)
    J_CELL = 0.5_dp*(J_RMS(IR)+J_RMS(IR+1))
    F_CELL = 0.5_dp*(F_H(IR)+F_H(IR+1))
    IF (J_CELL*F_CELL <= FLOOR) CYCLE   ! no Hall-CFL constraint where current/F_HALL vanish
    TC = MIN(TC, DR/(F_CELL*J_CELL))
  END DO
END FUNCTION HALL_COURANT_TIMESTEP

!> @param PHI, PSI, OPS, RGRID, ETA as JOULE_DISSIPATION_RATE.
!> @warning Pure-resistive Ohm's law only, same as JOULE_DISSIPATION_RATE
!>   -- but confirmed by the user (2026-08-11) to need an entirely new
!>   derivation for the Hall regime specifically, not just an added term:
!>   the diffusion-limit derivation this is built from does not carry
!>   over. Do not reuse this formula for a Hall/magnetofrictional Ohm's
!>   law even provisionally; wait for that derivation.
!> Returns: net Poynting flux rate through the domain's two boundaries
!>   (outer minus inner) -- a pure boundary evaluation, no radial
!>   integral. Assumes Phi/Psi represent a real field (see module
!>   header): the real part of the boundary sum is returned.
FUNCTION POYNTING_FLUX_RATE(PHI, PSI, OPS, RGRID, ETA) RESULT(E_DOT_S)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp),           INTENT(IN) :: ETA
  REAL(KIND=dp) :: E_DOT_S
  COMPLEX(KIND=dp) :: FLUX_TOTAL, DPHI_DR_BND, DPSI_DR_BND, CURV_PHI_BND
  REAL(KIND=dp)    :: LAMBDA_L, D2_CURV(RGRID%N,RGRID%N)
  INTEGER(KIND=i4) :: L, M, IDX, N

  N = RGRID%N
  FLUX_TOTAL = (0.0_dp, 0.0_dp)
  DO L = 0, PHI%LMAX
    LAMBDA_L = REAL(L*(L+1), KIND=dp)
    D2_CURV = OPS%D2
    CALL ADD_CURVATURE_TERM(D2_CURV, L, RGRID)
    DO M = -L, L
      IDX = YLM_INDEX(L, M)

      DPSI_DR_BND  = SUM(OPS%D1(N,:)*PSI%COEF(:,IDX))
      DPHI_DR_BND  = SUM(OPS%D1(N,:)*PHI%COEF(:,IDX))
      CURV_PHI_BND = SUM(D2_CURV(N,:)*PHI%COEF(:,IDX))
      FLUX_TOTAL = FLUX_TOTAL + LAMBDA_L * ( &
        DPSI_DR_BND*CONJG(PSI%COEF(N,IDX)) - CONJG(DPHI_DR_BND)*CURV_PHI_BND )

      DPSI_DR_BND  = SUM(OPS%D1(1,:)*PSI%COEF(:,IDX))
      DPHI_DR_BND  = SUM(OPS%D1(1,:)*PHI%COEF(:,IDX))
      CURV_PHI_BND = SUM(D2_CURV(1,:)*PHI%COEF(:,IDX))
      FLUX_TOTAL = FLUX_TOTAL - LAMBDA_L * ( &
        DPSI_DR_BND*CONJG(PSI%COEF(1,IDX)) - CONJG(DPHI_DR_BND)*CURV_PHI_BND )
    END DO
  END DO

  E_DOT_S = -(ETA/(4.0_dp*pi)) * REAL(FLUX_TOTAL, KIND=dp)
END FUNCTION POYNTING_FLUX_RATE

!> @param PHI, PSI, OPS, RGRID as JOULE_DISSIPATION_RATE.
!> @param F_HALL Hall prefactor (`c/(4*pi*e*n_e)`, the coefficient of the
!>   diffusion-free Hall induction equation `curl(F_HALL*j x B)`) --
!>   analytic_formulas/mhd-vsh-relations.tex, "Energy Budget in Hall
!>   Limit". A genuinely different derivation from JOULE_DISSIPATION_RATE/
!>   POYNTING_FLUX_RATE's pure-resistive one, not an extension of it: the
!>   Hall Joule term is exactly zero (j is always perpendicular to
!>   j x B, so JOULE_DISSIPATION_RATE's existing formula already covers
!>   the complete Joule term for this regime, unaffected -- see its own
!>   docstring), but the Hall Poynting flux is a genuinely new,
!>   cubic-in-field (mode-coupled) formula, computed here.
!> @note Coupling coefficient I^{nm}_{kl,k'l'} = FORTVSH's
!>   GWI(k,l,k',l',n,m) (the Gaunt coefficient
!>   `integral(Y_k^l*Y_k'^l'*conj(Y_n^m))dOmega`) -- the first use of
!>   mode-coupling in this module; every other FIELD_DIAGNOSTICS formula
!>   collapses to a diagonal sum thanks to VSH orthonormality, but this
!>   one is inherently cubic in the field (Hall is itself quadratic in
!>   B), so coupling is unavoidable.
!> Returns: net Hall-regime Poynting flux rate through the domain's two
!>   boundaries (outer minus inner) -- a pure boundary evaluation like
!>   POYNTING_FLUX_RATE. The (k,l),(k',l') pair is restricted to modes
!>   with a nonzero boundary-sampled quantity (an exact restriction, not
!>   an approximation -- see the active-mode-list comment in the body),
!>   and the target n/m are restricted via the coupling coefficient's
!>   own triangle-band/m=l+l' selection rule rather than searched, same
!>   optimization as HALL_INDUCTION_RHS -- what makes this tractable at
!>   LMAX=30 for a sparse state, where the full `O(Nlm**3)` worst case
!>   (the original, unoptimized form of this loop) would be infeasible
!>   called every logged timestep. Assumes Phi/Psi represent a real
!>   field, same as POYNTING_FLUX_RATE: the real part of the boundary
!>   sum is returned.
FUNCTION HALL_POYNTING_FLUX_RATE(PHI, PSI, OPS, RGRID, F_HALL) RESULT(E_DOT_S)
  TYPE(SPECTRAL_SCALAR_T), INTENT(IN) :: PHI, PSI
  TYPE(RADIAL_OPERATOR_T), INTENT(IN) :: OPS
  TYPE(RADIAL_GRID_T),     INTENT(IN) :: RGRID
  REAL(KIND=dp),           INTENT(IN) :: F_HALL
  REAL(KIND=dp) :: E_DOT_S

  COMPLEX(KIND=dp), ALLOCATABLE :: PHI_BND(:,:), PSI_BND(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: DPHI_BND(:,:), DPSI_BND(:,:), CURV_PHI_BND(:,:)
  REAL(KIND=dp)    :: D2_CURV(RGRID%N,RGRID%N)
  COMPLEX(KIND=dp) :: FLUX_TOTAL, COUPLING, TERM1, TERM2
  REAL(KIND=dp)    :: LAMBDA_K, LAMBDA_K2, LAMBDA_N, SGN(2)
  INTEGER(KIND=i4) :: BND_ROW(2), NLM
  INTEGER(KIND=i4) :: K, L, K2, L2, N, M, IDX_KL, IDX_K2L2, IDX_NM, IB
  INTEGER(KIND=i4), ALLOCATABLE :: ACTIVE_K(:), ACTIVE_L(:), ACTIVE_IDX(:)
  INTEGER(KIND=i4) :: N_ACTIVE, IA, IA2, KK, LL, IDX, N_LO, N_HI
  REAL(KIND=dp), PARAMETER :: ACTIVE_TOL = 1.0E-13_dp

  NLM = PHI%NLM
  BND_ROW = (/1, RGRID%N/)
  SGN      = (/-1.0_dp, 1.0_dp/)

  ALLOCATE(PHI_BND(2,NLM), PSI_BND(2,NLM))
  ALLOCATE(DPHI_BND(2,NLM), DPSI_BND(2,NLM), CURV_PHI_BND(2,NLM))

  ! Precompute every mode's boundary value/derivative/curvature once,
  ! outside the O(Nlm**3) coupling sum below, rather than re-evaluating
  ! them (potentially) millions of times.
  DO K = 0, PHI%LMAX
    D2_CURV = OPS%D2
    CALL ADD_CURVATURE_TERM(D2_CURV, K, RGRID)
    DO L = -K, K
      IDX_KL = YLM_INDEX(K, L)
      DO IB = 1, 2
        PHI_BND(IB,IDX_KL)      = PHI%COEF(BND_ROW(IB), IDX_KL)
        PSI_BND(IB,IDX_KL)      = PSI%COEF(BND_ROW(IB), IDX_KL)
        DPHI_BND(IB,IDX_KL)     = SUM(OPS%D1(BND_ROW(IB),:)*PHI%COEF(:,IDX_KL))
        DPSI_BND(IB,IDX_KL)     = SUM(OPS%D1(BND_ROW(IB),:)*PSI%COEF(:,IDX_KL))
        CURV_PHI_BND(IB,IDX_KL) = SUM(D2_CURV(BND_ROW(IB),:)*PHI%COEF(:,IDX_KL))
      END DO
    END DO
  END DO

  ! Active-mode list: any (k,l) with a nonzero boundary-sampled quantity
  ! (all of PHI_BND/PSI_BND/DPHI_BND/DPSI_BND/CURV_PHI_BND are exactly
  ! what every term below reads -- nothing else about the mode's radial
  ! profile enters this boundary-only flux, so this is an exact
  ! restriction, not an approximation, same reasoning as
  ! HALL_INDUCTION_RHS's identical optimization). Bounds the (k,l),
  ! (k',l') pair enumeration by N_ACTIVE**2 instead of NLM**2 -- what
  ! makes this tractable at LMAX=30 for a sparse/single-mode-seeded
  ! state (the full (LMAX+1)**2=961-mode space would otherwise force an
  ! infeasible ~961**3 raw loop, called every logged step).
  ALLOCATE(ACTIVE_K(NLM), ACTIVE_L(NLM), ACTIVE_IDX(NLM))
  N_ACTIVE = 0
  DO KK = 0, PHI%LMAX
    DO LL = -KK, KK
      IDX = YLM_INDEX(KK, LL)
      IF (ANY(ABS(PHI_BND(:,IDX)) > ACTIVE_TOL) .OR. ANY(ABS(PSI_BND(:,IDX)) > ACTIVE_TOL) .OR. &
          ANY(ABS(DPHI_BND(:,IDX)) > ACTIVE_TOL) .OR. ANY(ABS(DPSI_BND(:,IDX)) > ACTIVE_TOL) .OR. &
          ANY(ABS(CURV_PHI_BND(:,IDX)) > ACTIVE_TOL)) THEN
        N_ACTIVE = N_ACTIVE + 1
        ACTIVE_K(N_ACTIVE) = KK
        ACTIVE_L(N_ACTIVE) = LL
        ACTIVE_IDX(N_ACTIVE) = IDX
      END IF
    END DO
  END DO

  FLUX_TOTAL = (0.0_dp, 0.0_dp)
  DO IA = 1, N_ACTIVE
    K = ACTIVE_K(IA); L = ACTIVE_L(IA); IDX_KL = ACTIVE_IDX(IA)
    LAMBDA_K = REAL(K*(K+1), KIND=dp)
    DO IA2 = 1, N_ACTIVE
      K2 = ACTIVE_K(IA2); L2 = ACTIVE_L(IA2); IDX_K2L2 = ACTIVE_IDX(IA2)
      LAMBDA_K2 = REAL(K2*(K2+1), KIND=dp)

      ! m forced to l+l' (GWI's own selection rule); n restricted to the
      ! triangle band |k-k'|<=n<=k+k' -- both exact, not approximate.
      M = L + L2
      IF (ABS(M) > PHI%LMAX) CYCLE
      N_LO = MAX(ABS(K-K2), ABS(M))
      N_HI = MIN(PHI%LMAX, K+K2)
      DO N = N_LO, N_HI
        LAMBDA_N = REAL(N*(N+1), KIND=dp)
        COUPLING = GWI(K, L, K2, L2, N, M)
        IF (COUPLING == (0.0_dp,0.0_dp)) CYCLE
        IDX_NM = YLM_INDEX(N, M)

        DO IB = 1, 2
          TERM1 = -LAMBDA_K*(LAMBDA_K2+LAMBDA_N-LAMBDA_K) * PSI_BND(IB,IDX_KL) * &
            ( PSI_BND(IB,IDX_K2L2)*CONJG(PSI_BND(IB,IDX_NM)) + &
              DPHI_BND(IB,IDX_K2L2)*CONJG(DPHI_BND(IB,IDX_NM)) )
          TERM2 = LAMBDA_K2*(LAMBDA_K+LAMBDA_N-LAMBDA_K2) * PHI_BND(IB,IDX_K2L2) * &
            ( DPSI_BND(IB,IDX_KL)*CONJG(DPHI_BND(IB,IDX_NM)) - &
              PSI_BND(IB,IDX_NM)*CURV_PHI_BND(IB,IDX_KL) )
          FLUX_TOTAL = FLUX_TOTAL + SGN(IB)*COUPLING*(TERM1+TERM2)
        END DO
      END DO
    END DO
  END DO

  E_DOT_S = -(F_HALL/(8.0_dp*pi)) * REAL(FLUX_TOTAL, KIND=dp)
  DEALLOCATE(PHI_BND, PSI_BND, DPHI_BND, DPSI_BND, CURV_PHI_BND)
  DEALLOCATE(ACTIVE_K, ACTIVE_L, ACTIVE_IDX)
END FUNCTION HALL_POYNTING_FLUX_RATE

!> Trapezoidal radial integral of a density array sampled on RGRID%R.
FUNCTION RADIAL_INTEGRAL(F_OF_R, RGRID) RESULT(TOTAL)
  REAL(KIND=dp),        INTENT(IN) :: F_OF_R(:)
  TYPE(RADIAL_GRID_T),  INTENT(IN) :: RGRID
  REAL(KIND=dp) :: TOTAL
  INTEGER(KIND=i4) :: IR
  TOTAL = 0.0_dp
  DO IR = 1, RGRID%N-1
    TOTAL = TOTAL + 0.5_dp*(F_OF_R(IR)+F_OF_R(IR+1))*(RGRID%R(IR+1)-RGRID%R(IR))
  END DO
END FUNCTION RADIAL_INTEGRAL

END MODULE FIELD_DIAGNOSTICS
