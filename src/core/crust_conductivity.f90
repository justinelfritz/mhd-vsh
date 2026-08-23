MODULE CRUST_CONDUCTIVITY
!> Electron transport (electrical/thermal conductivity) in the neutron
!> star crust -- ported from Dany Page's NSCool
!> (`astroscu.unam.mx/neutrones/NSCool`, ASCL entry 1609.009),
!> `Code/conductivity_crust.f`'s `con_e_phon_ion_GYP` (electron-phonon
!> umklapp + electron-ion liquid-phase scattering) and its `OYAFORM`
!> dependency (nuclear structure from Oyamatsu's own density-only
!> parametrization -- ported from D.G.Yakovlev's `conrt.pas`, per that
!> subroutine's own header comment).
!>
!> `con_e_phon_ion_GYP` is explicitly cited (its own header) as "From
!> Appendix of Gnedin et al, MNRAS 324 (2001): 725, modified from
!> Potekhin et al, A&A 346 (1999): 345" -- the same two papers this
!> project's prior crust-conductivity port (from a different, private
!> research code) was already built on. This port changes which
!> published, citable implementation supplies the physics, not the
!> underlying physics itself (per user, 2026-08-23).
!>
!> @warning Every fitted constant is copied byte-for-byte from the
!>   original -- correctness here means matching NSCool's own output
!>   (see test_crust_conductivity.f90, checked against the unmodified
!>   original `conductivity_crust.f` directly), not re-deriving the
!>   physics from the cited papers.
!>
!> @warning One deliberate exception to "byte-for-byte", same finding
!>   as the prior EOSNS-based port's own identical issue: every literal
!>   here is written with an explicit `_dp` suffix (full double-
!>   precision parsing), whereas some of the original's bare literals
!>   lack a `d0` suffix and are silently parsed as single precision
!>   before being widened to double (a well-known F77 gotcha). This
!>   makes this port's arithmetic MORE precise than the original by
!>   ~1e-6 relative, not less faithful to the underlying physics --
!>   confirmed (not assumed) by rebuilding the original with
!>   `-fdefault-real-8`, which reproduces this port's output bit-for-
!>   bit. test_crust_conductivity.f90's tolerance (5e-6) accounts for
!>   exactly this effect.
!>
!> @warning Two genuine physics capabilities this port does NOT have,
!>   disclosed rather than silently dropped: (1) no magnetic-field
!>   dependence at all -- `con_e_phon_ion_GYP` takes no B argument
!>   (weaker even than this project's prior port, which at least had a
!>   dormant B-quantization branch); (2) no impurity scattering --
!>   unlike the prior port's `potekhinc.f` (which combined electron-
!>   phonon AND electron-impurity contributions via `COULIN`+`COUL99I`),
!>   `con_e_phon_ion_GYP` only ever computes electron-phonon/electron-ion
!>   scattering, with no `Z_imp`-style argument anywhere in its
!>   signature. Both are open gaps in the underlying NSCool routine
!>   itself, not something lost in translation.
!>
!> @warning `A_in`/`A1_in`/`Z_in`/`debug` dummy arguments in the
!>   original `con_e_phon_ion_GYP` are dropped in this port: its own
!>   header comment says outright "these being values defined in all
!>   other parts of this code, BUT THEY ARE RECALCULATED IN THE
!>   SUBROUTINE OYAFORM" -- confirmed by inspection, the body never
!>   reads them (Z/A1/A are freshly computed local variables from the
!>   `OYAFORM` call, never the same variables as the `_in` arguments).
!>   `debug` only gates `print*` trace statements, no physics effect.
!>   Genuinely dead inputs, not silently dropped live ones.
!>
!> @warning `N_E` is a NEW output added to this port's own
!>   `CON_E_PHON_ION_GYP` (not present in the original's argument list)
!>   -- electron number density is already computed as an internal step
!>   (`n_e=Z*n_i`) but never returned; MHD-VSH's own `ETA_AND_F_HALL_AT`
!>   needs it for `f_H`, and computing it a second time via a separate
!>   `OYAFORM` call would risk the two call sites drifting. A disclosed,
!>   purely-additive extension (SIGMA/LAMBDA/NU_E_S/NU_E_L unaffected),
!>   not a physics change.
USE KINDS,      ONLY: dp, i4
USE TOV_SOLVER, ONLY: TOV_PROFILE_T
USE UNITS,      ONLY: B_UNIT_GAUSS, LENGTH_UNIT_CM, TIME_UNIT_S
IMPLICIT NONE
PRIVATE
PUBLIC :: OYAFORM, CON_E_PHON_ION_GYP, ETA_AND_F_HALL_AT

CONTAINS

!> Nuclear structure (Z, A, form-factor parameters) from baryon density
!> alone -- ported verbatim from `subroutine OYAFORM` ("This subroutine
!> is from Oleg et al code! ... copied from 'conrt.pas' (D.G.Yakovlev)
!> and converted into Fortran. It realizes the SMOOTH COMPOSITION
!> model", per the original's own header). `SOyam` (a Fortran
!> statement function in the original) is ported as SOYAM below.
!>
!> @param BARD Baryon number density, `fm**-3`.
!> @param INDEX_PHASE 30 for densities below neutron drip, 3 for
!>   densities above (post-drip); ported argument name `Index` renamed
!>   to avoid shadowing the Fortran intrinsic INDEX.
!> Returns: Z (protons in nucleus), ANUC (baryons within the nucleus),
!>   A (baryons within the Wigner-Seitz cell), XNUC/XNUCT (effective
!>   proton-core-radius form-factor parameters, see CON_E_PHON_ION_GYP).
SUBROUTINE OYAFORM(BARD, INDEX_PHASE, Z, ANUC, A, XNUC, XNUCT)
  REAL(KIND=dp),    INTENT(IN)  :: BARD
  INTEGER(KIND=i4), INTENT(IN)  :: INDEX_PHASE
  REAL(KIND=dp),    INTENT(OUT) :: Z, ANUC, A, XNUC, XNUCT
  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp) :: F, RP, RN, NP_IN, NN_IN, NN_OUT, TP, TN, NIN, RWS
  REAL(KIND=dp) :: G_LOCAL, DN_N, NFREE
  REAL(KIND=dp) :: RP2EFF

  IF (INDEX_PHASE == 30) THEN   ! densities lower than the neutron drip
    F = LOG(1.0_dp + BARD/5.0E-9_dp)
    RP = 5.688_dp + 0.02628_dp*F + 0.009468_dp*F*F
    RN = 5.788_dp + 0.02077_dp*F + 0.01489_dp*F*F
    NP_IN = 0.0738_dp + 1.22E-4_dp*F - 1.641E-4_dp*F*F
    NN_IN = 0.0808_dp + 1.688E-4_dp*F + 9.439E-5_dp*F*F
    NN_OUT = 0.0_dp
    TP = 6.0_dp
    TN = TP
    NIN = PI/0.75_dp*RN**3*NN_IN*SOYAM(TN,1.0_dp)
    Z = PI/0.75_dp*RP**3*NP_IN*SOYAM(TP,1.0_dp)
    ANUC = Z + NIN
    A = ANUC
    RWS = (A*0.75_dp/PI/BARD)**0.333333_dp
    IF (RWS < RN) THEN
      WRITE(*,'(A)') 'OYAFORM: too large Rn for outer envelope!'; STOP 1
    END IF
  ELSE IF (INDEX_PHASE == 3) THEN   ! spheres after drip
    G_LOCAL = BARD*100.0_dp
    F = LOG(G_LOCAL)
    RWS = 31.68_dp - 8.400_dp*F - 0.2380_dp*F*F + 0.1152_dp*F**3
    TN = 1.0_dp/(0.2027_dp + 0.004506_dp*G_LOCAL)
    RN = 9.406_dp + 1.481_dp*F + 0.4625_dp*F*F + 0.05738_dp*F**3
    DN_N = (9.761_dp - 1.322_dp*F - 0.5544_dp*F*F - 0.07624_dp*F**3)/100.0_dp
    NIN = PI/0.75_dp*RN**3*DN_N*SOYAM(TN, MIN(1.0_dp,RWS/RN))
    TP = 1.0_dp/(0.1558_dp + 2.225E-3_dp*G_LOCAL + 9.452E-4_dp*G_LOCAL*G_LOCAL)
    RP = 8.345_dp + 0.7767_dp*F + 0.1333_dp*F*F + 0.008707_dp*F**3
    NP_IN = (4.040_dp - 1.097_dp*F - 0.0723_dp*F*F + 0.0225_dp*F**3)/100.0_dp
    Z = PI/0.75_dp*RP**3*NP_IN*SOYAM(TP, MIN(1.0_dp,RWS/RP))
    NFREE = BARD*PI/0.75_dp*RWS**3 - Z - NIN
    NN_OUT = NFREE/(PI/0.75_dp*RWS**3)
    NN_IN = NN_OUT + DN_N
    A = Z + NFREE + NIN
    ANUC = Z + NIN + NFREE*(RN/RWS)**3
    IF (RN > RWS) ANUC = A
  ELSE
    WRITE(*,'(A)') 'OYAFORM: invalid Index'; STOP 1
  END IF

  RP2EFF = RP*SQRT((1.0_dp - 15.0_dp/(5.0_dp+TP) + 15.0_dp/(5.0_dp+2.0_dp*TP) - &
    5.0_dp/(5.0_dp+3.0_dp*TP)) / SOYAM(TP,1.0_dp))
  XNUC = RP2EFF/RWS
  XNUCT = XNUC*TP/(0.6_dp+TP)
  ! NN_IN/NN_OUT computed above (matching the original exactly) but
  ! never read again there either -- genuinely unused past this point
  ! in the original, ported faithfully as such.
  ASSOCIATE (UNUSED_NN_IN => NN_IN, UNUSED_NN_OUT => NN_OUT); END ASSOCIATE
END SUBROUTINE OYAFORM

!> Oyamatsu's nuclear-shape volume-correction factor -- ported from
!> `SOyam`, a Fortran statement function in the original.
FUNCTION SOYAM(T, X) RESULT(RESULT_VAL)
  REAL(KIND=dp), INTENT(IN) :: T, X
  REAL(KIND=dp) :: RESULT_VAL
  RESULT_VAL = X**3 - 9.0_dp*X**(3.0_dp+T)/(3.0_dp+T) + &
    9.0_dp*X**(3.0_dp+2.0_dp*T)/(3.0_dp+2.0_dp*T) - X**(3.0_dp+3.0_dp*T)/(1.0_dp+T)
END FUNCTION SOYAM

!> Electrical/thermal conductivity from electron-phonon (crystal
!> umklapp) + electron-ion (liquid) scattering -- ported verbatim from
!> `subroutine con_e_phon_ion_GYP` (see module header for the `_in`
!> arguments/`debug` dropped and `N_E` added).
!>
!> @param T Temperature, K.
!> @param RHO Mass density, `g/cm**3`.
!> @param IFS 1: include the finite-nuclear-size correction (XNUC/XNUCT
!>   from OYAFORM); 0: XNUC=XNUCT=0.
!> Returns: SIGMA (electrical conductivity, `s**-1`), LAMBDA (thermal
!>   conductivity, erg/(K cm s)), NU_E_S/NU_E_L (longitudinal/
!>   transverse collision frequencies, `s**-1`), N_E (electron number
!>   density, `cm**-3` -- see module header's @warning).
SUBROUTINE CON_E_PHON_ION_GYP(T, RHO, IFS, SIGMA, LAMBDA, NU_E_S, NU_E_L, N_E)
  REAL(KIND=dp),    INTENT(IN)  :: T, RHO
  INTEGER(KIND=i4), INTENT(IN)  :: IFS
  REAL(KIND=dp),    INTENT(OUT) :: SIGMA, LAMBDA, NU_E_S, NU_E_L, N_E

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp), PARAMETER :: C_LIGHT = 2.99E10_dp, KB = 1.380E-16_dp, HB = 1.054E-27_dp
  REAL(KIND=dp), PARAMETER :: A_F = 1.0_dp/137.0_dp
  REAL(KIND=dp), PARAMETER :: E_CHARGE = 4.803E-10_dp, ME = 9.109E-28_dp, MU = 1.66E-24_dp
  REAL(KIND=dp), PARAMETER :: U_1 = 2.80_dp, U_2 = 13.00_dp   ! BCC lattice
  ! Standard neutron-drip density, g/cm**3 -- matches spec_heat.f's own
  ! hardcoded `parameter (rhodrip=4.3e11)` elsewhere in NSCool (the
  ! `common/rho_limits/` value this routine's original reads is set
  ! dynamically by NSCool's own full cooling-code initialization,
  ! unavailable to a standalone call like this port's own; 4.3e11 is
  ! the standard literature value both that hardcoded fallback and
  ! this project's own prior TOV/EOS port already used for the same
  ! quantity).
  REAL(KIND=dp), PARAMETER :: RHODRIP_CGS = 4.3E11_dp

  REAL(KIND=dp) :: BARD, Z, A1, A, XNUC, XNUCT
  INTEGER(KIND=i4) :: INDEX_PHASE
  REAL(KIND=dp) :: N_I, KF, PF, M_ST, EF, VF, X
  REAL(KIND=dp) :: OMEGA_P, T_P, TP, BETAZ, AI, GAMMA_ION
  REAL(KIND=dp) :: R_D, S_D, S_I, S_E, S_VAL, W, W1
  REAL(KIND=dp) :: G_S, G_L, D_FACTOR
  REAL(KIND=dp) :: LAM1A, LAM2A, LAMA, LAM1B, LAM2B, LAMB, LAM
  REAL(KIND=dp) :: LAM_S_HT, LAM_L_HT, T_U, LAM_0_LT, LAM_S_LT, LAM_L_LT
  REAL(KIND=dp) :: WW, LAM_S, LAM_L, NU0, NU_S, NU_L

  BARD = RHO/MU * 1.0E-39_dp
  IF (RHO > RHODRIP_CGS) THEN
    INDEX_PHASE = 3
  ELSE
    INDEX_PHASE = 30
  END IF
  CALL OYAFORM(BARD, INDEX_PHASE, Z, A1, A, XNUC, XNUCT)
  IF (IFS == 0) THEN
    XNUC = 0.0_dp
    XNUCT = 0.0_dp
  ELSE IF (IFS /= 1) THEN
    WRITE(*,'(A)') 'CON_E_PHON_ION_GYP: ifs badly defined'; STOP 1
  END IF

  N_I = RHO/(A*MU)
  N_E = Z*N_I
  KF = (3.0_dp*PI**2*N_E)**(1.0_dp/3.0_dp)
  PF = HB*KF
  M_ST = SQRT(ME**2 + (PF/C_LIGHT)**2)
  EF = M_ST*C_LIGHT**2
  VF = PF/M_ST
  X = PF/ME/C_LIGHT

  OMEGA_P = SQRT(4.0_dp*PI*E_CHARGE**2 * Z**2*N_I/A1/MU)
  T_P = HB*OMEGA_P/KB
  TP = T/T_P
  BETAZ = PI*A_F*Z*VF/C_LIGHT
  AI = (3.0_dp/(4.0_dp*PI*N_I))**(1.0_dp/3.0_dp)
  GAMMA_ION = Z**2*E_CHARGE**2/(KB*T*AI)

  R_D = AI/SQRT(3.0_dp*GAMMA_ION)
  S_D = 1.0_dp/(2.0_dp*KF*R_D)**2
  S_I = S_D*(1.0_dp+0.06_dp*GAMMA_ION)*EXP(-SQRT(GAMMA_ION))
  S_E = A_F/PI*C_LIGHT/VF
  S_VAL = (S_I+S_E)*EXP(-BETAZ)
  W = (U_2/S_D)*(1.0_dp+BETAZ/3.0_dp)
  W1 = 14.73_dp * XNUC**2 * (1.0_dp+Z*SQRT(XNUC)/13.0_dp) * (1.0_dp+BETAZ/3.0_dp)

  G_S = 1.0_dp/SQRT(1.0_dp+0.0361_dp/Z**(1.0_dp/3.0_dp)/TP**2) * (1.0_dp+0.122_dp*BETAZ**2)
  G_L = G_S + 0.0105_dp*TP/(TP**2+0.0081_dp)**1.5_dp * (1.0_dp+(VF/C_LIGHT)**3*BETAZ) * &
    (1.0_dp-1.0_dp/Z) * (1.0_dp+XNUCT**2*SQRT(2.0_dp*Z))

  D_FACTOR = EXP(-0.42_dp*SQRT(X/A/Z)*U_1*EXP(-9.1_dp*TP))

  W = W + W1
  CALL GET_LAM(S_VAL, W, LAM1A, LAM2A)
  LAMA = LAM1A - (VF/C_LIGHT)**2*LAM2A
  W = W1
  CALL GET_LAM(S_VAL, W, LAM1B, LAM2B)
  LAMB = LAM1B - (VF/C_LIGHT)**2*LAM2B
  LAM = LAMA - LAMB

  LAM_S_HT = LAM*G_S*D_FACTOR
  LAM_L_HT = LAM*G_L*D_FACTOR

  T_U = T_P*Z**(1.0_dp/3.0_dp)*A_F/3.0_dp/VF*C_LIGHT
  LAM_0_LT = 50.0_dp*X**0.5_dp/A1**0.5_dp/Z
  LAM_S_LT = LAM_0_LT*(4.0_dp/3.0_dp)*A_F*C_LIGHT/VF*TP**5
  LAM_L_LT = LAM_0_LT*TP**3

  WW = EXP(-T_U/T)
  LAM_S = LAM_S_HT*WW + LAM_S_LT*(1.0_dp-WW)
  LAM_L = LAM_L_HT*WW + LAM_L_LT*(1.0_dp-WW)

  NU0 = 4.0_dp*Z*EF*A_F**2/3.0_dp/PI/HB
  NU_S = NU0 * LAM_S
  NU_L = NU0 * LAM_L
  SIGMA  = N_E * E_CHARGE**2 / (M_ST * NU_S)
  LAMBDA = PI**2*KB**2 * T * N_E / (3.0_dp * M_ST * NU_L)
  NU_E_S = NU_S
  NU_E_L = NU_L
END SUBROUTINE CON_E_PHON_ION_GYP

!> Coulomb-logarithm helper, ported verbatim from `subroutine get_lam`.
SUBROUTINE GET_LAM(S_VAL, W, LAM1, LAM2)
  REAL(KIND=dp), INTENT(IN)  :: S_VAL, W
  REAL(KIND=dp), INTENT(OUT) :: LAM1, LAM2
  REAL(KIND=dp), PARAMETER :: EPS = 0.05_dp
  REAL(KIND=dp), PARAMETER :: EULER = 0.5772156_dp

  IF (S_VAL <= EPS .AND. S_VAL*W <= EPS) THEN
    LAM1 = 0.5_dp*(EXP_INT(W)+LOG(W)+EULER)
    LAM2 = (EXP(-W)-1.0_dp+W)/(2.0_dp*W)
  ELSE IF (W <= EPS) THEN
    LAM1 = W*((2.0_dp*S_VAL+1.0_dp)/(2.0_dp*S_VAL+2.0_dp) - S_VAL*LOG((S_VAL+1.0_dp)/S_VAL))
    LAM2 = W*((1.0_dp-3.0_dp*S_VAL-6.0_dp*S_VAL**2)/(4.0_dp*S_VAL+4.0_dp) + &
      1.5_dp*LOG((S_VAL+1.0_dp)/S_VAL))
  ELSE IF (W > (1.0_dp/EPS)) THEN
    LAM1 = 0.5_dp*(LOG((S_VAL+1.0_dp)/S_VAL) - 1.0_dp/(S_VAL+1.0_dp))
    LAM2 = (2.0_dp*S_VAL+1.0_dp)/(2.0_dp*S_VAL+2.0_dp) - S_VAL*LOG((S_VAL+1.0_dp)/S_VAL)
  ELSE
    LAM1 = LOG((S_VAL+1.0_dp)/S_VAL) + S_VAL/(S_VAL+1.0_dp)*(1.0_dp-EXP(-W)) - &
      (1.0_dp+S_VAL*W)*EXP(S_VAL*W)*(EXP_INT(S_VAL*W)-EXP_INT(S_VAL*W+W))
    LAM1 = 0.5_dp*LAM1
    LAM2 = (EXP(-W)-1.0_dp+W)/W - S_VAL**2/(S_VAL+1.0_dp)*(1.0_dp-EXP(-W)) - &
      2.0_dp*S_VAL*LOG((S_VAL+1.0_dp)/S_VAL) + &
      S_VAL*(2.0_dp+S_VAL*W)*EXP(S_VAL*W)*(EXP_INT(S_VAL*W)-EXP_INT(S_VAL*W+W))
    LAM2 = 0.5_dp*LAM2
  END IF
END SUBROUTINE GET_LAM

!> exp(x)*E_1(x) rational/polynomial approximation, accurate to 1e-5,
!> ported verbatim from `function exp_int` (checked by the original
!> against Abramowitz & Stegun -- NOT the same formula/algorithm as
!> this project's prior port's own EXPINT function, a different
!> exponential-integral evaluator entirely; named distinctly here to
!> avoid implying a relationship that doesn't exist).
FUNCTION EXP_INT(X) RESULT(RESULT_VAL)
  REAL(KIND=dp), INTENT(IN) :: X
  REAL(KIND=dp) :: RESULT_VAL
  REAL(KIND=dp) :: NUM, DEN

  IF (X <= 0.0_dp) THEN
    WRITE(*,'(A)') 'EXP_INT: x must be > 0 !'; STOP 1
  END IF
  IF (X >= 1.0_dp) THEN
    NUM = X**4 + 8.5733287401_dp*X**3 + 18.0590169730_dp*X**2 + 8.6347608925_dp*X + 0.2677737343_dp
    DEN = X**4 + 9.5733223454_dp*X**3 + 25.6329561486_dp*X**2 + 21.0996530827_dp*X + 3.9584969228_dp
    RESULT_VAL = NUM/DEN / (X*EXP(X))
  ELSE
    RESULT_VAL = -0.57721566_dp + 0.99999193_dp*X - 0.24991055_dp*X**2 + 0.05519968_dp*X**3 &
      - 0.00976004_dp*X**4 + 0.00107857_dp*X**5 - LOG(X)
  END IF
END FUNCTION EXP_INT

!> Self-consistent, radially-varying resistivity and Hall pre-factor
!> from a solved TOV crust profile -- new code (not a port), combining
!> CON_E_PHON_ION_GYP's SIGMA/N_E outputs into `eta = c**2/(4*pi*sigma)`
!> (standard magnetic diffusivity from conductivity) and
!> `f_H = c/(4*pi*e*n_e)` (already established in hall_induction.f90's
!> own docstring), converted from Gaussian-cgs into this project's code
!> units via UNITS. Simpler than the prior port's version of this
!> routine: CON_E_PHON_ION_GYP returns SIGMA directly, no more combining
!> a relaxation time with a separately-tabulated effective mass/electron
!> density the way `nstot.f`'s own `sigmae=3.26*tau*nel/meff` formula
!> required.
!>
!> Only rows with RHOCGS below the core-crust boundary
!> (`RHOL_CGS=2.2e14 g/cm**3`, TOV.f's own commented threshold for its
!> `rhol` variable, matching this project's prior port's identical
!> threshold) are supported -- `CON_E_PHON_ION_GYP` is meant for crust
!> densities; core rows are out of scope here as before.
!>
!> @param PROFILE Solved TOV radial profile (TOV_SOLVER::SOLVE_TOV_STAR).
!> @param ETA_PROFILE Resistivity per radial row, code units (out, allocated here).
!> @param F_HALL_PROFILE Hall pre-factor per radial row, code units (out, allocated here).
!> @param N_E_PROFILE Electron number density per radial row, `cm**-3`
!>   (out, allocated here) -- CON_E_PHON_ION_GYP's own N_E output,
!>   threaded through since it's already computed per row and a driver
!>   reporting eta(r)/f_H(r) will typically want n_e(r) alongside them.
!> @param T_KELVIN Crust temperature, K (optional; default 1e9 K,
!>   matching this project's prior port's own default, itself matching
!>   `nstot.f`'s hardcoded isothermal-crust value).
SUBROUTINE ETA_AND_F_HALL_AT(PROFILE, ETA_PROFILE, F_HALL_PROFILE, N_E_PROFILE, T_KELVIN)
  TYPE(TOV_PROFILE_T), INTENT(IN)  :: PROFILE
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: ETA_PROFILE(:), F_HALL_PROFILE(:), N_E_PROFILE(:)
  REAL(KIND=dp), OPTIONAL, INTENT(IN) :: T_KELVIN

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp), PARAMETER :: C_LIGHT_CGS = 2.99792458E10_dp   ! cm/s
  REAL(KIND=dp), PARAMETER :: E_CHARGE_ESU = 4.80320425E-10_dp ! esu
  REAL(KIND=dp), PARAMETER :: RHOL_CGS = 2.2E14_dp             ! core-crust boundary
  INTEGER(KIND=i4), PARAMETER :: IFS = 1   ! finite-nuclear-size correction ON

  REAL(KIND=dp) :: T
  REAL(KIND=dp) :: SIGMA, LAMBDA_TH, NU_E_S, NU_E_L, N_E_CGS
  REAL(KIND=dp) :: ETA_CGS, F_HALL_CGS
  INTEGER(KIND=i4) :: I

  T = 1.0E9_dp; IF (PRESENT(T_KELVIN)) T = T_KELVIN

  ALLOCATE(ETA_PROFILE(PROFILE%N), F_HALL_PROFILE(PROFILE%N), N_E_PROFILE(PROFILE%N))

  DO I = 1, PROFILE%N
    IF (PROFILE%RHOCGS(I) > RHOL_CGS) THEN
      WRITE(*,'(A)') 'ETA_AND_F_HALL_AT: core-density row (rho > 2.2e14 g/cm**3) not supported'
      STOP 1
    END IF

    CALL CON_E_PHON_ION_GYP(T, PROFILE%RHOCGS(I), IFS, SIGMA, LAMBDA_TH, NU_E_S, NU_E_L, N_E_CGS)
    ASSOCIATE (UNUSED_LAMBDA => LAMBDA_TH, UNUSED_NU_E_S => NU_E_S, UNUSED_NU_E_L => NU_E_L)
    END ASSOCIATE

    ETA_CGS    = C_LIGHT_CGS**2 / (4.0_dp*PI*SIGMA)                        ! cm**2/s
    F_HALL_CGS = C_LIGHT_CGS / (4.0_dp*PI*E_CHARGE_ESU*N_E_CGS)            ! cm**2/(G*s)

    ETA_PROFILE(I)    = ETA_CGS * TIME_UNIT_S / LENGTH_UNIT_CM**2
    F_HALL_PROFILE(I) = F_HALL_CGS * B_UNIT_GAUSS * TIME_UNIT_S / LENGTH_UNIT_CM**2
    N_E_PROFILE(I)    = N_E_CGS
  END DO
END SUBROUTINE ETA_AND_F_HALL_AT

END MODULE CRUST_CONDUCTIVITY
