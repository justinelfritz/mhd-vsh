MODULE CRUST_CONDUCTIVITY
!> Electron transport (electrical/thermal conductivity) in the neutron
!> star crust.
!>
!> `OYAFORM` (nuclear structure from Oyamatsu's own density-only
!> parametrization -- ported from D.G.Yakovlev's `conrt.pas`) is ported
!> from Dany Page's NSCool (`astroscu.unam.mx/neutrones/NSCool`, ASCL
!> entry 1609.009), `Code/conductivity_crust.f`'s own `OYAFORM`
!> dependency, and is UNCHANGED from this project's prior port.
!>
!> The conductivity ENGINE itself (given T, rho, Z, A, xnuc, xnuct ->
!> SIGMA, CKAPPA) was previously also from NSCool: two independently-fit
!> subroutines, `con_e_phon_ion_GYP` (crust) and `con_env_e_phon_ion_PBHY`
!> (envelope), spliced at a hardcoded rho=6e7 g/cm**3 threshold. That
!> splice produced a genuine ~44x non-monotonic jump in `eta(r)` at the
!> seam (found 2026-08-23 via a standalone harness against the
!> unmodified original, at the SAME density and composition on both
!> sides) -- not a porting bug, but an artifact of NSCool's own 2001-era
!> two-formula design (confirmed: at this project's T=1e9 K default the
!> two formulas only actually agree near rho~7.9e9 g/cm**3, two orders
!> of magnitude from where NSCool switches; at T>=3e9 K they never agree
!> within either formula's own stated validity range).
!>
!> Per user direction (2026-08-23), the engine is replaced with a single,
!> un-spliced formula: `CONDUCT`/`ThAv18`/`COUL19` from Alexander
!> Potekhin's own actively-maintained conductivity code, `conduct21.f`
!> (http://www.ioffe.ru/astro/conduct/, downloaded and read in full,
!> last updated 2021 there). Unlike NSCool's GYP/PBHY split, `COUL19`
!> takes finite nuclear size (`xnuc`,`xnuct`) as a CONTINUOUS INPUT
!> PARAMETER to the SAME formula body at every density -- no branch
!> selector, no seam. Confirmed via a second independent harness (this
!> port's own regression target) that wiring `OYAFORM`'s real,
!> continuously density-dependent `xnuc`/`xnuct` through this single
!> engine removes the ~44x jump entirely (row-to-row change across the
!> old seam is ~1% with this engine, not ~4400%).
!>
!> References (all confirmed directly, not from memory):
!> - A. Potekhin, `conduct21.f`, http://www.ioffe.ru/astro/conduct/
!>   -- the direct source of this port's `CONDUCT_TRANSPORT` and its
!>   full dependency chain below.
!> - Potekhin, `Astron. Astrophys.` 351 (1999): 787 -- cited directly in
!>   the original's own `CONDCONV` header for the transport-coefficient
!>   definitions this routine returns.
!> - Gnedin, Yakovlev & Potekhin, `MNRAS` 324 (2001): 725, Appendix
!>   A1.1 -- cited directly in the original's own `COULAN3` header for
!>   the finite-nuclear-size Coulomb-logarithm formula.
!> - Potekhin, Pons & Page, `Space Sci. Rev.` 191 (2015): 239 -- the
!>   code's own stated primary citation at ioffe.ru.
!>
!> @warning Every fitted constant is copied byte-for-byte from the
!>   original -- correctness here means matching `conduct21.f`'s own
!>   output (see test_crust_conductivity.f90, checked against the
!>   unmodified original directly via two independent standalone F77
!>   harnesses), not re-deriving the physics from the cited papers.
!>
!> @warning Same single-precision-literal finding as this project's
!>   other ports of F77 physics code: bare (non-`d0`) literals in the
!>   original are parsed as single precision before widening to double;
!>   this port's `_dp`-suffixed literals are more precise by ~1e-6
!>   relative. test_crust_conductivity.f90's tolerance accounts for
!>   this.
!>
!> @warning Two latent bugs/quirks in the ORIGINAL `conduct21.f` were
!>   found while reading it line-by-line (confirmed via `grep` across
!>   the whole downloaded file, not assumed) and are preserved here
!>   faithfully rather than silently "fixed":
!>   (1) `COUL19`'s own `C13` variable (used in `XW1=XW1*(1.d0+C13*
!>   BORNCOR)*...`) is referenced but never assigned anywhere in the
!>   file -- an uninitialized `SAVE`'d F77 scalar, which reads as 0 from
!>   zeroed BSS storage on essentially every real compiler (including
!>   the gfortran build this port's own reference values were generated
!>   from). Ported as an explicit `C13=0.0_dp` parameter.
!>   (2) `TAUEESY`'s high-`Y` branch (`Y>1e8`) computes a term into a
!>   variable named `CITL`, but the formula that USES it two lines later
!>   reads a DIFFERENT, similarly-spelled variable, `CILT` -- which is
!>   only ever assigned in the OTHER (`Y<=1e8`) branch. In the `Y>1e8`
!>   branch `CILT` is therefore also an uninitialized-SAVE'd (effectively
!>   zero) scalar. Ported as: compute the (unused) `CITL` term faithfully,
!>   and use `0.0_dp` for `CILT` in that branch, matching the original's
!>   actual (buggy) behavior rather than the evidently-intended one.
!>
!> @warning Genuine, disclosed scope reduction versus the full
!>   `conduct21.f`: only routines reachable from `CONDUCT` with `B0=0`
!>   (no magnetic field -- `ETA_AND_F_HALL_AT` has never had a B
!>   argument) and `Zimp=0` (no impurity scattering -- same pre-existing
!>   gap as this project's prior port) are exercised/tested. The
!>   magnetic branches inside `COUL19`/`ThAv18`/`ThAvI18` and the
!>   B0>0 branch of `CHEMPOT` ARE ported (not stubbed out) for interface
!>   completeness, but are UNVERIFIED -- no regression coverage exists
!>   for B0>0. `COUL18I` (impurity Coulomb log), `DENRFIT`/`BLINW`/
!>   `FitFERMI` (only reachable from `CHEMPOT`'s magnetic branch),
!>   `CONDI`/`HLfit8` (ion thermal conduction, disabled by default even
!>   in the original's own demo driver), and `BLOUIN20`/`BLOUIN1`/
!>   `BLOUIN2` (H/He envelope corrections, irrelevant for the iron-group
!>   crust composition `OYAFORM` produces) are confirmed unreachable for
!>   this project's usage (dependency-traced directly, not guessed) and
!>   are not ported at all. `CHEMPOT_FIT`'s own B0>0 branch explicitly
!>   `STOP`s rather than silently calling a nonexistent `DENRFIT`.
!>
!> @warning Deliberate deviation from `conduct21.f`'s own top-level
!>   `CONDUCT` interface: the original hardcodes `xnuc=0.` internally
!>   (with a commented-out, cruder placeholder formula labelled "outer
!>   cr." right next to it -- evidently a demo-driver default, not a
!>   real production treatment). This port instead exposes `XNUC`/
!>   `XNUCT` as real arguments to `CONDUCT_CORE`, fed from `OYAFORM`'s
!>   own already-tested, continuously density-dependent values -- this
!>   is the change that actually eliminates the old GYP/PBHY seam (a
!>   constant `xnuc=0` throughout would also be seam-free, but would
!>   silently drop the finite-nuclear-size correction everywhere,
!>   confirmed to matter by up to a factor of ~6.6x in `SIGMA` deep in
!>   the crust, rho=1.22e14 g/cm**3, via a second independent harness
!>   built specifically to check this). Verified against a
!>   COMMON-block-patched copy of the unmodified original (see
!>   test_crust_conductivity.f90's own header for the exact values).
!>
!> @warning `A_in`/`A1_in`/`Z_in`/`debug`-style dummy arguments are not
!>   part of this port's interfaces at all -- the composition (`Z`,`A`,
!>   `xnuc`,`xnuct`) is always freshly computed by `OYAFORM`, matching
!>   this project's own established convention (unchanged from the
!>   prior GYP/PBHY-based port).
!>
!> @warning `OYAFORM`'s own pre-drip/post-drip switch (`INDEX_PHASE` 30
!>   vs 3, at rho=4.3e11 g/cm**3, the standard neutron-drip density)
!>   still causes a real, smaller (order tens-of-percent, not ~44x)
!>   discontinuity in `A`/`xnuc` -- confirmed by direct scan. This
!>   reflects an actual physical phase transition (neutron drip), not a
!>   modeling-artifact splice, and is NOT addressed by this port
!>   (`OYAFORM` itself is unchanged).
USE KINDS,      ONLY: dp, i4
USE TOV_SOLVER, ONLY: TOV_PROFILE_T
USE UNITS,      ONLY: B_UNIT_GAUSS, LENGTH_UNIT_CM, TIME_UNIT_S
IMPLICIT NONE
PRIVATE
PUBLIC :: OYAFORM, CONDUCT_TRANSPORT, CON_CRUST, ETA_AND_F_HALL_AT

CONTAINS

!> Nuclear structure (Z, A, form-factor parameters) from baryon density
!> alone -- ported verbatim from `subroutine OYAFORM` ("This subroutine
!> is from Oleg et al code! ... copied from 'conrt.pas' (D.G.Yakovlev)
!> and converted into Fortran. It realizes the SMOOTH COMPOSITION
!> model", per the original's own header). `SOyam` (a Fortran
!> statement function in the original) is ported as SOYAM below.
!> Unchanged from this project's prior NSCool-based port.
!>
!> @param BARD Baryon number density, `fm**-3`.
!> @param INDEX_PHASE 30 for densities below neutron drip, 3 for
!>   densities above (post-drip); ported argument name `Index` renamed
!>   to avoid shadowing the Fortran intrinsic INDEX.
!> Returns: Z (protons in nucleus), ANUC (baryons within the nucleus),
!>   A (baryons within the Wigner-Seitz cell), XNUC/XNUCT (effective
!>   proton-core-radius form-factor parameters, see CONDUCT_TRANSPORT).
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

!> Public entry point mirroring the original's own `CONDCONV` -- converts
!> (T,rho,B,Zion,CMI,Zimp) into the relativistic units `CONDUCT_CORE`
!> works in, calls it, converts the outputs back to CGS. Ported from
!> `subroutine CONDCONV` verbatim (unit-conversion constants copied
!> byte-for-byte).
!>
!> @param T_KELVIN Temperature, K.
!> @param RHO Mass density, g/cm**3.
!> @param B_GAUSS Magnetic field, Gauss (0 for the field-free case --
!>   the only path this port's own tests exercise, see module header).
!> @param ZION, CMI Ion charge and mass numbers.
!> @param ZIMP Impurity parameter (0 for the pure-composition case --
!>   the only path this port's own tests exercise).
!> @param XNUC, XNUCT Finite-nuclear-size form-factor parameters (see
!>   OYAFORM) -- deliberately exposed as real arguments here rather than
!>   hardcoded to 0 as in the original's own top-level `CONDUCT`, see
!>   module header's @warning.
!> Returns: SIGMA (electrical conductivity, s**-1), CKAPPA (thermal
!>   conductivity, erg/(K cm s)), QJ (thermopower, k_B/e) -- longitudinal
!>   components; SIGMAT/CKAPPAT/QJT (transverse), SIGMAH/CKAPPAH/QJH
!>   (Hall, off-diagonal) -- all zero when B_GAUSS=0.
SUBROUTINE CONDUCT_TRANSPORT(T_KELVIN, RHO, B_GAUSS, ZION, CMI, ZIMP, XNUC, XNUCT, &
    SIGMA, CKAPPA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH)
  REAL(KIND=dp), INTENT(IN)  :: T_KELVIN, RHO, B_GAUSS, ZION, CMI, ZIMP, XNUC, XNUCT
  REAL(KIND=dp), INTENT(OUT) :: SIGMA, CKAPPA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH

  REAL(KIND=dp), PARAMETER :: AUM = 1822.9_dp, AUD = 15819.4_dp, BOHR = 137.036_dp
  REAL(KIND=dp), PARAMETER :: UNISIG = 7.763E20_dp, UNIKAP = 2.778E15_dp
  REAL(KIND=dp), PARAMETER :: UNITEMP = 5930.0_dp, UNIB = 4.414E13_dp

  REAL(KIND=dp) :: T6, TEMR, B0, DENRI
  REAL(KIND=dp) :: RSIGMA, RLET, RLTT, RQ, RKAPPA
  REAL(KIND=dp) :: RTSIGMA, RTLET, RTLTT, RTQ, RTKAPPA
  REAL(KIND=dp) :: RHSIGMA, RHLET, RHLTT, RHQ, RHKAPPA

  T6 = T_KELVIN/1.0E6_dp
  TEMR = T6/UNITEMP
  B0 = B_GAUSS/UNIB
  DENRI = RHO/(AUD*AUM*CMI)

  CALL CONDUCT_CORE(TEMR, DENRI, B0, ZION, CMI, ZIMP, XNUC, XNUCT, &
    RSIGMA, RLET, RLTT, RQ, RKAPPA, &
    RTSIGMA, RTLET, RTLTT, RTQ, RTKAPPA, &
    RHSIGMA, RHLET, RHLTT, RHQ, RHKAPPA)

  SIGMA   = RSIGMA*UNISIG
  CKAPPA  = RKAPPA*UNIKAP
  QJ      = RQ/SQRT(BOHR)
  SIGMAT  = RTSIGMA*UNISIG
  CKAPPAT = RTKAPPA*UNIKAP
  QJT     = RTQ/SQRT(BOHR)
  SIGMAH  = RHSIGMA*UNISIG
  CKAPPAH = RHKAPPA*UNIKAP
  QJH     = RHQ/SQRT(BOHR)
END SUBROUTINE CONDUCT_TRANSPORT

!> Central transport-coefficient driver -- ported from `subroutine
!> CONDUCT`. See module header's @warning for the XNUC/XNUCT interface
!> deviation from the original (hardcoded xnuc=0 there, real arguments
!> here).
SUBROUTINE CONDUCT_CORE(TEMR, DENRI, B0, ZION, CMI, ZIMP, XNUC, XNUCT, &
    RSIGMA, RLET, RLTT, RQ, RKAPPA, &
    RTSIGMA, RTLET, RTLTT, RTQ, RTKAPPA, &
    RHSIGMA, RHLET, RHLTT, RHQ, RHKAPPA)
  REAL(KIND=dp), INTENT(IN)  :: TEMR, DENRI, B0, ZION, CMI, ZIMP, XNUC, XNUCT
  REAL(KIND=dp), INTENT(OUT) :: RSIGMA, RLET, RLTT, RQ, RKAPPA
  REAL(KIND=dp), INTENT(OUT) :: RTSIGMA, RTLET, RTLTT, RTQ, RTKAPPA
  REAL(KIND=dp), INTENT(OUT) :: RHSIGMA, RHLET, RHLTT, RHQ, RHKAPPA

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp), PARAMETER :: AUM = 1822.9_dp, BOHR = 137.036_dp

  REAL(KIND=dp) :: DENR, SPHERION, GAMMA, XSR, CNL
  INTEGER(KIND=i4) :: MAGNET
  REAL(KIND=dp) :: CMU, DDENR, Q2E
  REAL(KIND=dp) :: XNIMP
  REAL(KIND=dp) :: F0B, F1B, F2B, F0C, F1C, F2C, F0D, F1D, F2D
  REAL(KIND=dp) :: D_DET
  REAL(KIND=dp) :: CTH, TAULONGT, TAUEE, VCL, TRP, BORNCOR, T0, G0, G2, THTOEL, EECOR

  IF (TEMR <= 0.0_dp .OR. DENRI <= 0.0_dp .OR. B0 < 0.0_dp .OR. ZION <= 0.0_dp .OR. CMI <= 0.0_dp) THEN
    WRITE(*,'(A,5ES14.6)') 'CONDUCT_CORE: Non-positive input parameter: ', TEMR, DENRI, B0, ZION, CMI
    STOP 1
  END IF
  IF (ZION < 0.5_dp) THEN
    WRITE(*,'(A)') 'CONDUCT_CORE: Too small ion charge'; STOP 1
  END IF
  IF (CMI < 1.0_dp) THEN
    WRITE(*,'(A)') 'CONDUCT_CORE: Too small ion mass'; STOP 1
  END IF
  IF (DENRI > 1.0E5_dp) THEN
    WRITE(*,'(A)') 'CONDUCT_CORE: Too high density'; STOP 1
  END IF

  DENR = DENRI*ZION
  SPHERION = (0.75_dp/PI/DENRI)**0.3333333_dp
  GAMMA = ZION**2/BOHR/TEMR/SPHERION
  XSR = (3.0_dp*PI**2*DENR)**0.3333333_dp

  IF (B0 > 0.0_dp) THEN
    CNL = ((SQRT(1.0_dp+XSR**2)+3.0_dp*TEMR)**2-1.0_dp)/2.0_dp/B0
  ELSE
    CNL = 10000.0_dp
  END IF
  IF (CNL < 200.0_dp) THEN
    MAGNET = 1_i4
  ELSE
    MAGNET = 0_i4
  END IF

  CMU = CHEMPOT_FIT(B0, DENR, DDENR, TEMR)
  Q2E = 4.0_dp*PI/BOHR*DDENR

  XNIMP = XNUC
  CALL THAV18(MAGNET, CMU, B0, TEMR, DENRI, &
    XSR, GAMMA, ZION, CMI, ZIMP, Q2E, &
    XNUC, XNUCT, XNIMP, &
    F0B, F1B, F2B, F0C, F1C, F2C, F0D, F1D, F2D)

  RSIGMA = F0B/BOHR
  RLET = F1B/SQRT(BOHR)
  RLTT = F2B*TEMR
  RQ = RLET/RSIGMA
  RKAPPA = RLTT - TEMR*RLET**2/RSIGMA

  RTSIGMA = F0C/BOHR
  RTLET = F1C/SQRT(BOHR)
  RTLTT = F2C*TEMR
  RHSIGMA = F0D/BOHR
  RHLET = F1D/SQRT(BOHR)
  RHLTT = F2D*TEMR
  D_DET = RTSIGMA**2 + RHSIGMA**2
  RTQ = (RTSIGMA*RTLET + RHSIGMA*RHLET)/D_DET
  RHQ = (RTSIGMA*RHLET - RHSIGMA*RTLET)/D_DET
  RTKAPPA = RTLTT - TEMR*(RTLET*RTQ - RHLET*RHQ)
  RHKAPPA = RHLTT - TEMR*(RTLET*RHQ + RHLET*RTQ)

  IF (B0 < TEMR*CMU .OR. B0 < CMU .OR. CMU < TEMR) THEN
    CTH = TEMR/9.0_dp*XSR**3/SQRT(1.0_dp+XSR**2)
    TAULONGT = RKAPPA/CTH
    CALL TAUEESY(XSR, TEMR, TAUEE)
    VCL = XSR/SQRT(1.0_dp+XSR**2)
    TRP = ZION/GAMMA*SQRT(CMI*AUM*SPHERION/3.0_dp/BOHR)
    BORNCOR = VCL*ZION*PI/BOHR
    T0 = 0.19_dp/ZION**0.16667_dp
    G0 = TRP/SQRT(TRP**2+T0**2)*(1.0_dp+(ZION/125.0_dp)**2)
    G2 = TRP/SQRT(0.0081_dp+TRP**2)**3
    THTOEL = 1.0_dp + G2/G0*(1.0_dp+BORNCOR*VCL**3)*0.0105_dp*(1.0_dp-1.0_dp/ZION)
    EECOR = TAUEE/(TAULONGT+TAUEE*THTOEL)
    RKAPPA = RKAPPA*EECOR
    RTKAPPA = RTKAPPA*EECOR
    RHKAPPA = RHKAPPA*EECOR
  END IF
END SUBROUTINE CONDUCT_CORE

!> Thermal averaging over the electron energy distribution -- ported
!> from `subroutine ThAv18`. For B=0 (`MAGNET=0`) the Landau-level loop
!> collapses to a single one-shot Simpson integration over energy (the
!> Gauss-quadrature-over-cyclotron-levels branch, labels matching the
!> original's own 11/12, is provably unreached -- confirmed by reading
!> the MAGNET branch logic directly). Two dead labels in the original
!> (11, 15 -- nothing ever jumps to them, confirmed by an independent
!> gfortran unused-label warning on the unmodified source) are dropped;
!> the "jump to this loop's own end" pattern at those two points (guarded
!> by `if (MARK.eq.1) goto <loop end>`) is ported as CYCLE. The outer
!> radial/Landau-band loop's own GOTO-driven repeat-with-adjusted-bounds
!> structure (labels 13/20/50) is preserved literally rather than
!> restructured, to avoid silently altering its iteration semantics.
SUBROUTINE THAV18(MAGNET, CMU, B, TEMR, DENRI, &
    XSR, GAMMA, ZION, CMI, ZIMP, Q2E, &
    XNUC, XNUCT, XNIMP, &
    F0, F1, F2, F0C, F1C, F2C, F0D, F1D, F2D)
  INTEGER(KIND=i4), INTENT(IN) :: MAGNET
  REAL(KIND=dp), INTENT(IN)  :: CMU, B, TEMR, DENRI, XSR, GAMMA, ZION, CMI, ZIMP, Q2E
  REAL(KIND=dp), INTENT(IN)  :: XNUC, XNUCT, XNIMP
  REAL(KIND=dp), INTENT(OUT) :: F0, F1, F2, F0C, F1C, F2C, F0D, F1D, F2D

  REAL(KIND=dp), PARAMETER :: TAIL1 = 8.0_dp, TAIL2 = 40.0_dp
  INTEGER(KIND=i4), PARAMETER :: NM = 128, NXI = 4
  REAL(KIND=dp), PARAMETER :: AXI(4) = [0.32425342340381_dp, 0.61337143270059_dp, &
    0.83603110732664_dp, 0.96816023950763_dp]
  REAL(KIND=dp), PARAMETER :: AWI(4) = [0.312347077040_dp, 0.260610696429_dp, &
    0.180648160695_dp, 0.081274388362_dp]

  REAL(KIND=dp) :: CST, EMIN, EMAX, E1, E2
  INTEGER(KIND=i4) :: NLMIN, NLMAX, NL1, NL2, NL, IK
  REAL(KIND=dp) :: EN, ENP, EFIN, ESTART, SCALE, DE, H, X, W, WI, XI, E, T
  REAL(KIND=dp) :: F, FC, FD
  REAL(KIND=dp) :: S0, S1, S2, S0C, S1C, S2C, S0D, S1D, S2D
  REAL(KIND=dp) :: DS0, DS1, DS2, DS0C, DS1C, DS2C, DS0D, DS1D, DS2D
  REAL(KIND=dp) :: DDS, DDSC, DDSD, BLINK
  INTEGER(KIND=i4) :: IX, K, MARK

  IF (MAGNET /= 0_i4 .AND. MAGNET /= 1_i4) THEN
    WRITE(*,'(A)') 'THAV18: incorrect MAGNET'; STOP 1
  END IF
  CST = 4.0_dp*3.14159265_dp*DENRI*(ZION/137.036_dp)**2
  EMIN = MAX(1.0_dp, CMU-TAIL2*TEMR)
  EMAX = CMU+TAIL2*TEMR
  E1 = MAX(1.0_dp, CMU-TAIL1*TEMR)
  E2 = MAX(1.0_dp, CMU+TAIL1*TEMR)
  IF (MAGNET == 1_i4) THEN
    NLMIN = INT((EMIN**2-1.0_dp)/2.0_dp/B, i4)
    NLMAX = INT((EMAX**2-1.0_dp)/2.0_dp/B, i4)
    NL1   = INT((E1**2-1.0_dp)/2.0_dp/B, i4)
    NL2   = INT((E2**2-1.0_dp)/2.0_dp/B, i4)
  ELSE
    NLMIN = 0_i4; NLMAX = 0_i4; NL1 = 0_i4; NL2 = 0_i4
  END IF

  S0 = 0.0_dp; S1 = 0.0_dp; S2 = 0.0_dp
  S0C = 0.0_dp; S1C = 0.0_dp; S2C = 0.0_dp
  S0D = 0.0_dp; S1D = 0.0_dp; S2D = 0.0_dp

  DO NL = NLMIN, NLMAX
    IF (MAGNET == 1_i4) THEN
      EN = SQRT(1.0_dp+2.0_dp*B*REAL(NL,dp))
      ENP = SQRT(EN**2+2.0_dp*B)
    ELSE
      EN = EMIN
      ENP = EMAX
      GOTO 20
    END IF
    IF ((NL >= NL1-1 .AND. NL <= NL2+1) .AND. NL2 < NL1+5) GOTO 20
    EFIN = ENP
    SCALE = (EFIN-EN)/TEMR
    DS0 = 0.0_dp; DS1 = 0.0_dp; DS2 = 0.0_dp
    DS0C = 0.0_dp; DS1C = 0.0_dp; DS2C = 0.0_dp
    DS0D = 0.0_dp; DS1D = 0.0_dp; DS2D = 0.0_dp
    DO IX = 1, NXI
      XI = AXI(IX)
      WI = AWI(IX)
      X = XI**2
      W = 2.0_dp*XI**2*WI
      E = EN+(EFIN-EN)*X
      CALL THAVI18(E, CMU, TEMR, CST, &
        XSR, GAMMA, B, ZION, CMI, ZIMP, Q2E, XNUC, XNUCT, XNIMP, &
        T, F, FC, FD, MARK)
      IF (MARK == 1_i4) CYCLE
      F = F/XI
      DDS = W*F*SCALE
      DS0 = DS0+DDS; DS1 = DS1+DDS*T; DS2 = DS2+DDS*T**2
      FC = FC/XI
      DDSC = W*FC*SCALE
      DS0C = DS0C+DDSC; DS1C = DS1C+DDSC*T; DS2C = DS2C+DDSC*T**2
      FD = FD/XI
      DDSD = W*FD*SCALE
      DS0D = DS0D+DDSD; DS1D = DS1D+DDSD*T; DS2D = DS2D+DDSD*T**2
    END DO
    S0 = S0+DS0; S1 = S1+DS1; S2 = S2+DS2
    S0C = S0C+DS0C; S1C = S1C+DS1C; S2C = S2C+DS2C
    S0D = S0D+DS0D; S1D = S1D+DS1D; S2D = S2D+DS2D
    GOTO 50
    20 CONTINUE
    ESTART = MAX(EN, EMIN)
    IK = 0_i4
    13 CONTINUE
    IK = IK+1_i4
    IF (IK > 1_i4) ESTART = EFIN
    EFIN = MIN(ENP, EMAX)
    IF (IK == 1_i4) THEN
      IF (NL == NL1) EFIN = E1
      IF (NL == NL2 .AND. NL /= NL1) EFIN = E2
    END IF
    IF (IK == 2_i4 .AND. NL == NL1 .AND. NL == NL2) EFIN = E2
    DE = EFIN-ESTART
    SCALE = DE/TEMR
    H = 1.0_dp/REAL(NM,dp)
    X = 0.0_dp
    BLINK = 1.0_dp
    W = H/3.0_dp*SCALE
    E = ESTART
    CALL THAVI18(E, CMU, TEMR, CST, &
      XSR, GAMMA, B, ZION, CMI, ZIMP, Q2E, XNUC, XNUCT, XNIMP, &
      T, F, FC, FD, MARK)
    DDS = W*F
    DS0 = DDS; DS1 = DDS*T; DS2 = DDS*T**2
    DDSC = W*FC
    DS0C = DDSC; DS1C = DDSC*T; DS2C = DDSC*T**2
    DDSD = W*FD
    DS0D = DDSD; DS1D = DDSD*T; DS2D = DDSD*T**2
    DO K = 1, NM
      X = X+H
      E = ESTART+DE*X
      CALL THAVI18(E, CMU, TEMR, CST, &
        XSR, GAMMA, B, ZION, CMI, ZIMP, Q2E, XNUC, XNUCT, XNIMP, &
        T, F, FC, FD, MARK)
      IF (MARK == 1_i4) CYCLE
      DDS = W*F*(3.0_dp+BLINK)
      DS0 = DS0+DDS; DS1 = DS1+DDS*T; DS2 = DS2+DDS*T**2
      DDSC = W*FC*(3.0_dp+BLINK)
      DS0C = DS0C+DDSC; DS1C = DS1C+DDSC*T; DS2C = DS2C+DDSC*T**2
      DDSD = W*FD*(3.0_dp+BLINK)
      DS0D = DS0D+DDSD; DS1D = DS1D+DDSD*T; DS2D = DS2D+DDSD*T**2
      BLINK = -BLINK
    END DO
    S0 = S0+DS0-DDS/2.0_dp; S1 = S1+DS1-DDS*T/2.0_dp; S2 = S2+DS2-DDS*T**2/2.0_dp
    S0C = S0C+DS0C-DDSC/2.0_dp; S1C = S1C+DS1C-DDSC*T/2.0_dp; S2C = S2C+DS2C-DDSC*T**2/2.0_dp
    S0D = S0D+DS0D-DDSD/2.0_dp; S1D = S1D+DS1D-DDSD*T/2.0_dp; S2D = S2D+DS2D-DDSD*T**2/2.0_dp
    IF (IK == 1_i4 .AND. (NL == NL1 .OR. NL == NL2)) GOTO 13
    IF (NL1 == NL2 .AND. IK == 2_i4) GOTO 13
    50 CONTINUE
  END DO

  F0 = S0; F1 = S1; F2 = S2
  F0C = S0C; F1C = S1C; F2C = S2C
  F0D = S0D; F1D = S1D; F2D = S2D
END SUBROUTINE THAV18

!> Integrand for thermal averaging -- ported from `subroutine ThAvI18`.
!> `XNIMP` is a genuine input carried for interface fidelity (would feed
!> `COUL18I` if impurity scattering were ported) but is unread on the
!> `ZIMP=0` path this project's usage always takes.
SUBROUTINE THAVI18(E, CMU, TEMR, CST, &
    XSR, GAMMA, B, ZION, CMI, ZIMP, Q2E, XNUC, XNUCT, XNIMP, &
    T, F, FC, FD, MARK)
  REAL(KIND=dp), INTENT(IN)  :: E, CMU, TEMR, CST
  REAL(KIND=dp), INTENT(IN)  :: XSR, GAMMA, B, ZION, CMI, ZIMP, Q2E, XNUC, XNUCT, XNIMP
  REAL(KIND=dp), INTENT(OUT) :: T, F, FC, FD
  INTEGER(KIND=i4), INTENT(OUT) :: MARK

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp) :: PCL, EX, F1
  REAL(KIND=dp) :: CLEFF, CLLONG, CLTRAN, SN, THTOEL
  REAL(KIND=dp) :: GYROM, TAULONG, TAUT, TAUTRAN, TAUHALL, C_FACTOR

  T = (E-CMU)/TEMR
  PCL = SQRT(E**2-1.0_dp)
  IF (ABS(T) > 40.0_dp) THEN
    F = 0.0_dp; FC = 0.0_dp; FD = 0.0_dp
    MARK = 1_i4
    RETURN
  ELSE
    EX = EXP(T)
    F1 = EX/(EX+1.0_dp)**2
  END IF
  IF (E-1.0_dp < 1.0E-10_dp) THEN
    F = 0.0_dp; FC = 0.0_dp; FD = 0.0_dp
    MARK = 2_i4
    RETURN
  END IF

  CALL COUL19_FIT(PCL, XSR, GAMMA, B, ZION, CMI, Q2E, XNUC, XNUCT, &
    CLEFF, CLLONG, CLTRAN, SN, THTOEL)
  ASSOCIATE (UNUSED_CLEFF => CLEFF, UNUSED_THTOEL => THTOEL); END ASSOCIATE
  GYROM = B/E
  TAULONG = PCL**3/E/CST/CLLONG/SN
  TAUT = PCL**3/E/CST/CLTRAN*SN
  IF (ZIMP > 0.0_dp) THEN
    WRITE(*,'(A)') 'THAVI18: impurity scattering (ZIMP>0) not implemented in this port'
    STOP 1
  END IF
  ASSOCIATE (UNUSED_XNIMP => XNIMP); END ASSOCIATE
  TAUTRAN = TAUT/(1.0_dp+(TAUT*GYROM)**2)
  TAUHALL = TAUT*GYROM*TAUTRAN
  C_FACTOR = SN*PCL**3/3.0_dp/PI**2/E
  F = F1*C_FACTOR*TAULONG
  FC = F1*C_FACTOR*TAUTRAN
  FD = F1*C_FACTOR*TAUHALL
  MARK = 2_i4
END SUBROUTINE THAVI18

!> Coulomb logarithm + phonon/lattice physics -- ported from `subroutine
!> COUL19`. Non-magnetic path (`PCL**2>4d2*B`, always true at B=0)
!> returns immediately after computing CLEFF; the magnetic-fit code below
!> it (unreached for this project's B=0-only current usage) is ported
!> for interface completeness but UNVERIFIED (see module header).
SUBROUTINE COUL19_FIT(PCL, XSR, GAMMA, B, ZION, CMI, Q2E, XNUC, XNUCT, &
    CLEFF, CLLONG, CLTRAN, SN, THTOEL)
  REAL(KIND=dp), INTENT(IN)  :: PCL, XSR, GAMMA, B, ZION, CMI, Q2E, XNUC, XNUCT
  REAL(KIND=dp), INTENT(OUT) :: CLEFF, CLLONG, CLTRAN, SN, THTOEL

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265359_dp
  REAL(KIND=dp), PARAMETER :: AUM = 1822.88848_dp
  REAL(KIND=dp), PARAMETER :: BOHR = 137.035999_dp
  REAL(KIND=dp), PARAMETER :: UMINUS1 = 2.79855_dp, UMINUS2 = 12.972_dp
  REAL(KIND=dp), PARAMETER :: BIG = 1.0E6_dp
  ! C13 is referenced in the original but never assigned anywhere in
  ! conduct21.f (confirmed by grep across the whole file -- this line
  ! is its only appearance). Ported as an explicit zero -- see module
  ! header's @warning.
  REAL(KIND=dp), PARAMETER :: C13 = 0.0_dp

  REAL(KIND=dp) :: DENR, DENRI, SPHERION, Q2ICL, ECL, VCL, PM2, TRP, BORNCOR
  REAL(KIND=dp) :: C_COEF, Q2S, XS, R2W, XW, XW1, CL, DLEFF
  REAL(KIND=dp) :: A0, VIBRCOR, T0, G0, GW, G2

  REAL(KIND=dp) :: ENU
  INTEGER(KIND=i4) :: NL, N
  REAL(KIND=dp) :: TEMR_LOCAL, GAMMAN1, TAUCL, GAMMAN2, GAMMAN, PB
  REAL(KIND=dp) :: XIS, ZETA, XI_LOC, XSUM, Q2M, QTM1, QTM2, QTRANM
  REAL(KIND=dp) :: QTP1, QTP2, QTRANP, Q_FACTOR
  REAL(KIND=dp) :: DNU, XS1, PN, SQB, X_LOC
  REAL(KIND=dp) :: EXW, A1_COEF, Q1, DLT, Y1, CL0, P2, Y2, PY, DT
  REAL(KIND=dp) :: DB, CL1, P1, P3, P4

  DENR = XSR**3/3.0_dp/PI**2
  DENRI = DENR/ZION
  SPHERION = (0.75_dp/PI/DENRI)**0.3333333_dp
  Q2ICL = 3.0_dp*GAMMA/SPHERION**2
  ECL = SQRT(1.0_dp+PCL**2)
  VCL = PCL/ECL
  PM2 = (2.0_dp*PCL)**2
  TRP = ZION/GAMMA*SQRT(CMI*AUM*SPHERION/3.0_dp/BOHR)
  BORNCOR = VCL*ZION*PI/BOHR

  C_COEF = (1.0_dp+0.06_dp*GAMMA)*EXP(-SQRT(GAMMA))
  Q2S = (Q2ICL*C_COEF+Q2E)*EXP(-BORNCOR)
  XS = Q2S/PM2
  R2W = UMINUS2/Q2ICL*(1.0_dp+0.333333_dp*BORNCOR)
  XW = R2W*PM2
  XW1 = 14.7327_dp*XNUC**2
  XW1 = XW1*(1.0_dp+C13*BORNCOR)*(1.0_dp+ZION/13.0_dp*SQRT(XNUC))
  XW1 = XW1*(PCL/XSR)**2
  CALL COULAN3_FIT(XS, XW, PCL, XW1, CL, DLEFF)
  ASSOCIATE (UNUSED_DLEFF => DLEFF); END ASSOCIATE
  A0 = 1.683_dp*SQRT(PCL/CMI/ZION)
  VIBRCOR = EXP(-A0/4.0_dp*UMINUS1*EXP(-9.1_dp*TRP))
  T0 = 0.19_dp/ZION**0.16667_dp
  G0 = TRP/SQRT(TRP**2+T0**2)*(1.0_dp+(ZION/125.0_dp)**2)
  GW = G0*VIBRCOR
  CLEFF = CL*GW
  G2 = TRP/SQRT(0.0081_dp+TRP**2)**3
  THTOEL = 1.0_dp+G2/G0*(1.0_dp+BORNCOR*VCL**3)*0.0105_dp*(1.0_dp-1.0_dp/ZION)* &
    (1.0_dp+XNUCT**2*SQRT(2.0_dp*ZION))

  IF (PCL**2 > 4.0E2_dp*B) THEN
    CLLONG = CLEFF
    CLTRAN = CLEFF
    SN = 1.0_dp
    RETURN
  END IF

  ! ---- Magnetic fit: only reached when B>0. UNVERIFIED by this port's
  ! own regression tests (ETA_AND_F_HALL_AT never passes a nonzero
  ! field) -- ported for interface completeness, not stubbed out, per
  ! module header.
  ENU = PCL**2/2.0_dp/B
  NL = INT(ENU, i4)
  TEMR_LOCAL = (ZION/BOHR)**2/SPHERION/GAMMA
  IF (GAMMA < 175.0_dp) THEN
    GAMMAN1 = SQRT(TEMR_LOCAL/CMI/AUM)*2.0_dp*PCL
    TAUCL = PCL**2*VCL*BOHR**2/(4.0_dp*PI*ZION**2*DENRI)
  ELSE
    GAMMAN1 = TEMR_LOCAL/TRP
    TAUCL = VCL*BOHR/(TEMR_LOCAL*(2.0_dp-VCL**2)*UMINUS2)
  END IF
  GAMMAN2 = 1.0_dp/TAUCL
  GAMMAN = MAX(GAMMAN1, GAMMAN2)
  SN = 0.0_dp
  DO N = 0, NL
    PB = SQRT(ENU-REAL(N,dp))
    SN = SN+PB
    IF (N /= 0) SN = SN+PB
  END DO
  SN = SN*1.5_dp*B*SQRT(2.0_dp*B)/PCL**3

  IF (ENU <= 1.0_dp) THEN
    XIS = Q2S/2.0_dp/B
    ZETA = R2W*2.0_dp*B
    XI_LOC = 2.0_dp*PCL**2/B
    XSUM = XI_LOC+XIS
    Q2M = (EXPINT_FIT(XSUM,1_i4)- &
      EXP(-ZETA*XI_LOC)*EXPINT_FIT((1.0_dp+ZETA)*XSUM,1_i4))/XSUM
    CLLONG = (PCL*VCL/B)**2*Q2M/1.5_dp*GW
    IF (XSUM < BIG) THEN
      QTM1 = (1.0_dp+XSUM)*EXPINT_FIT(XSUM,0_i4)-1.0_dp
      QTM2 = (1.0_dp+(1.0_dp+ZETA)*XSUM)*EXPINT_FIT((1.0_dp+ZETA)*XSUM,0_i4)-1.0_dp
    ELSE
      QTM1 = 1.0_dp/XSUM**2
      QTM2 = 1.0_dp/((1.0_dp+ZETA)*XSUM)**2
    END IF
    QTRANM = QTM1-EXP(-ZETA*XI_LOC)*QTM2
    IF (XIS < BIG) THEN
      QTP1 = (1.0_dp+XIS)*EXPINT_FIT(XIS,0_i4)-1.0_dp
      QTP2 = (1.0_dp+(1.0_dp+ZETA)*XIS)*EXPINT_FIT((1.0_dp+ZETA)*XIS,0_i4)-1.0_dp
    ELSE
      QTP1 = 1.0_dp/XIS**2
      QTP2 = 1.0_dp/((1.0_dp+ZETA)*XIS)**2
    END IF
    QTRANP = QTP1-QTP2
    Q_FACTOR = (ECL**2*QTRANP+QTRANM)*B/PCL**2
    CLTRAN = 0.375_dp*Q_FACTOR/ECL**2*GW
  ELSE
    DNU = ENU-REAL(NL,dp)
    DNU = MAX(DNU, GAMMAN)
    XS1 = (SQRT(XS)+1.0_dp/(2.0_dp+XW/2.0_dp))**2
    PN = SQRT(2.0_dp*B*DNU)
    SQB = SQRT(B)
    X_LOC = MAX(PN/SQB, 1.0E-10_dp)
    IF (XW < 0.01_dp) THEN
      EXW = 1.0_dp
    ELSE IF (XW > 50.0_dp) THEN
      EXW = 1.0_dp/XW
    ELSE
      EXW = (1.0_dp-EXP(-XW))/XW
    END IF
    A1_COEF = (30.0_dp-15.0_dp*EXW-(15.0_dp-6.0_dp*EXW)*VCL**2)/ &
      (30.0_dp-10.0_dp*EXW-(20.0_dp-5.0_dp*EXW)*VCL**2)
    Q1 = 0.25_dp*VCL**2/(1.0_dp-0.6666667_dp*VCL**2)
    DLT = SQB/PCL*(A1_COEF/X_LOC-SQRT(X_LOC)*(1.5_dp-0.5_dp*EXW+Q1)+ &
      (1.0_dp-EXW+0.75_dp*VCL**2)/(1.0_dp+VCL**2)*(X_LOC-SQRT(X_LOC))/REAL(NL,dp))
    Y1 = 1.0_dp/(1.0_dp+DLT)
    CL0 = LOG(1.0_dp+1.0_dp/XS1)
    P2 = CL0*(0.07_dp+0.2_dp*EXW)
    Y2 = 1.5_dp*CL0*(X_LOC**3-X_LOC/3.0_dp)/(REAL(NL,dp)+0.75_dp/(1.0_dp+2.0_dp*B)**2*X_LOC**2)+P2*X_LOC
    PY = 1.0_dp+0.06_dp*CL0**2/REAL(NL,dp)**2
    DT = SQRT(PY*Y1**2+Y2**2)
    CLLONG = CLEFF/DT
    DB = 1.0_dp/(1.0_dp+0.5_dp/B)
    CL1 = XS1*CL0
    P1 = 0.8_dp*(1.0_dp+CL1)+0.2_dp*CL0
    P2 = 1.42_dp-0.1_dp*DB+SQRT(CL1)/3.0_dp
    P3 = (0.68_dp-0.13_dp*DB)*CL1**0.165_dp
    P4 = (0.52_dp-0.1_dp*DB)*SQRT(SQRT(CL1))
    ! Original: `alog(NL+0.)` -- single-precision LOG of an
    ! integer-promoted-to-single-real argument (F77 `alog` is the
    ! single-precision intrinsic; `NL+0.` promotes via a single-
    ! precision literal). Ported as full double precision, LOG(REAL(NL,dp)),
    ! same class of single-precision-literal finding disclosed in the
    ! module header for other routines in this file.
    DLT = SQB/PCL*(P1/X_LOC**2*SQB/PCL+P3*LOG(REAL(NL,dp))/X_LOC- &
      (P2+P4*LOG(REAL(NL,dp)))*SQRT(X_LOC))
    CLTRAN = CLEFF*(1.0_dp+DLT)
  END IF
END SUBROUTINE COUL19_FIT

!> Analytic Coulomb logarithm -- ported from `subroutine COULAN3`.
SUBROUTINE COULAN3_FIT(XS, XW0, PCL, XW1, CLEFF, DLEFF)
  REAL(KIND=dp), INTENT(IN)  :: XS, XW0, PCL, XW1
  REAL(KIND=dp), INTENT(OUT) :: CLEFF, DLEFF

  REAL(KIND=dp), PARAMETER :: EPS = 1.0E-2_dp, EPS1 = 1.0E-3_dp, TINY_VAL = 1.0E-9_dp
  REAL(KIND=dp), PARAMETER :: EULER = 0.5772156649_dp
  REAL(KIND=dp), PARAMETER :: BIG = 1.0_dp/EPS

  REAL(KIND=dp) :: E_VAL, V_VAL, XW, B_VAL, EA, E1_VAL, E2_VAL, XS1
  REAL(KIND=dp) :: CL0, CL1, CL2, CL10, CL20, EL
  REAL(KIND=dp) :: EXW0, DL1
  INTEGER(KIND=i4) :: KEY, MI, I

  ! V_VAL is referenced by the original's own opening validity check
  ! BEFORE it's actually computed a few lines later (confirmed by
  ! reading the source directly). Same class of finding as
  ! COUL19_FIT's C13 -- an uninitialized SAVE'd F77 scalar, reads as 0
  ! from zeroed BSS in practice. Reproduced explicitly.
  V_VAL = 0.0_dp
  IF (XS < 0.0_dp .OR. XW0 < 0.0_dp .OR. V_VAL < 0.0_dp .OR. XW1 < 0.0_dp) THEN
    WRITE(*,'(A)') 'COULAN3_FIT: invalid input'; STOP 1
  END IF
  E_VAL = SQRT(1.0_dp+PCL**2)
  V_VAL = PCL/E_VAL
  KEY = 0_i4
  MI = 1_i4
  IF (XW1 < TINY_VAL) MI = 0_i4

  CL2 = 0.0_dp
  E2_VAL = 0.0_dp
  DO I = 0, MI
    IF (I == 0) THEN
      XW = XW0+XW1
      B_VAL = XS*XW
    ELSE
      XW = XW1
      B_VAL = XS*XW
    END IF
    IF (I == 0 .OR. KEY == 2_i4) THEN
      IF (XW < EPS) THEN
        KEY = 1_i4
        GOTO 50
      END IF
      IF (XW > BIG .AND. B_VAL > BIG) THEN
        KEY = 2_i4
      ELSE IF (XS < EPS1 .AND. B_VAL < EPS1/(1.0_dp+XW)) THEN
        KEY = 3_i4
      ELSE
        KEY = 4_i4
      END IF
    END IF
    50 CONTINUE
    EA = EXP(-XW)
    E1_VAL = 1.0_dp-EA
    IF (XW < EPS1) E1_VAL = XW*(1.0_dp-XW/2.0_dp)
    IF (KEY /= 1_i4) THEN
      IF (XW > EPS1) THEN
        E2_VAL = (XW-E1_VAL)/XW
      ELSE
        E2_VAL = 0.5_dp*XW*(1.0_dp-XW/3.0_dp)
      END IF
    END IF
    XS1 = XS+1.0_dp
    IF (KEY == 1_i4) THEN
      CL0 = LOG(XS1/XS)
      CL1 = 0.5_dp*XW*(2.0_dp-1.0_dp/XS1-2.0_dp*XS*CL0)
      CL2 = 0.5_dp*XW*(1.5_dp-3.0_dp*XS-1.0_dp/XS1+3.0_dp*XS**2*CL0)
    ELSE IF (KEY == 2_i4) THEN
      CL0 = LOG(XS1/XS)
      CL1 = (CL0-1.0_dp/XS1-1.0_dp/B_VAL**2)/2.0_dp
      CL2 = (2.0_dp*XS+1.0_dp)/(2.0_dp*XS1)-XS*CL0
    ELSE IF (KEY == 3_i4) THEN
      CL1 = 0.5_dp*(EA*EXPINT_FIT(XW,0_i4)+LOG(XW)+EULER)
      CL2 = 0.5_dp*E2_VAL
    ELSE IF (KEY >= 4_i4) THEN
      CL0 = LOG(XS1/XS)
      EL = EXPINT_FIT(B_VAL,0_i4)-EXPINT_FIT(B_VAL+XW,0_i4)*EA
      CL1 = 0.5_dp*(CL0+XS/XS1*E1_VAL-(1.0_dp+B_VAL)*EL)
      CL2 = 0.5_dp*(E2_VAL-XS*XS/XS1*E1_VAL-2.0_dp*XS*CL0+XS*(2.0_dp+B_VAL)*EL)
    ELSE
      WRITE(*,'(A)') 'COULAN3_FIT: invalid KEY'; STOP 1
    END IF
    IF (I == 0) THEN
      CLEFF = CL1-V_VAL**2*CL2
      CL10 = CL1
      CL20 = CL2
    ELSE
      CLEFF = CLEFF-(CL1-V_VAL**2*CL2)
      CL1 = CL10-CL1
      CL2 = CL20-CL2
    END IF
  END DO

  IF (XW0 > EPS) THEN
    EXW0 = 1.0_dp-EXP(-XW0)
  ELSE
    EXW0 = XW0
  END IF
  DL1 = EXW0/(PCL*V_VAL)/XS1**2*EXP(-XW1)
  DLEFF = DL1/E_VAL**2+2.0_dp*V_VAL**2/E_VAL*CL2
END SUBROUTINE COULAN3_FIT

!> Chemical potential of the free electron gas -- ported from `function
!> CHEMPOT`. The B0=0 path (this project's only exercised path) skips
!> straight to CHEMP99_FIT, matching the original's own early return.
!> The B0>0 iterative branch is NOT ported: it depends on `DENRFIT`,
!> confirmed unreachable for this project's current usage and not ported
!> at all (see module header) -- rather than silently mis-handle a
!> future B0>0 caller, this STOPs explicitly.
FUNCTION CHEMPOT_FIT(B0, DENR, DDENR, TEMR) RESULT(CMU_OUT)
  REAL(KIND=dp), INTENT(IN)  :: B0, DENR, TEMR
  REAL(KIND=dp), INTENT(OUT) :: DDENR
  REAL(KIND=dp) :: CMU_OUT

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp) :: PF0, CNL
  REAL(KIND=dp) :: CHEMFIT_VAL

  PF0 = (3.0_dp*PI**2*DENR)**0.3333333_dp

  IF (B0 > 0.0_dp) THEN
    CNL = ((SQRT(1.0_dp+PF0**2)+6.0_dp*TEMR)**2-1.0_dp)/2.0_dp/B0
  ELSE
    CNL = 10000.0_dp
  END IF

  IF (CNL < 200.0_dp) THEN
    WRITE(*,'(A)') 'CHEMPOT_FIT: magnetic branch (B0>0) requires DENRFIT, not implemented in this port'
    STOP 1
  END IF

  CALL CHEMP99_FIT(DENR, TEMR, CHEMFIT_VAL, DDENR)
  CMU_OUT = CHEMFIT_VAL
END FUNCTION CHEMPOT_FIT

!> Fit to the chemical potential of the free electron gas -- ported
!> verbatim from `subroutine CHEMP99`.
SUBROUTINE CHEMP99_FIT(DENR, TEMR, CHEMFIT, DDENR)
  REAL(KIND=dp), INTENT(IN)  :: DENR, TEMR
  REAL(KIND=dp), INTENT(OUT) :: CHEMFIT, DDENR

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp), PARAMETER :: AICH = 0.25954_dp, BICH = 0.072_dp, B1ICH = 0.858_dp

  REAL(KIND=dp) :: XSR, THETA, THETA32, P01, T_VAL, Y_VAL, P02, P03
  REAL(KIND=dp) :: UP, DN, CT, F_VAL, X1, CMU1
  REAL(KIND=dp) :: THETALOG, THETACOR, TB, TB1, DENOM, DMT, DTP

  XSR = (3.0_dp*PI**2*DENR)**0.33333333_dp
  IF (XSR > 0.001_dp) THEN
    THETA = TEMR/(SQRT(1.0_dp+XSR**2)-1.0_dp)
  ELSE
    THETA = 2.0_dp*TEMR/XSR**2
  END IF
  THETA32 = THETA*SQRT(THETA)
  P01 = 12.0_dp+8.0_dp/THETA32
  T_VAL = EXP(MIN(THETA,40.0_dp))
  Y_VAL = (1.0_dp/T_VAL+1.612_dp*T_VAL)/ &
    (6.192_dp*THETA**0.0944_dp/T_VAL+5.535_dp*T_VAL*THETA**0.698_dp)
  P02 = 1.3656_dp-Y_VAL
  IF (THETA > 1.0E-5_dp) THEN
    P03 = 1.5_dp/(T_VAL-1.0_dp)
  ELSE
    P03 = 1.5_dp/THETA
  END IF
  UP = 1.0_dp+P01*TEMR*P02+P03*SQRT(TEMR)
  DN = (1.0_dp+0.5_dp*TEMR/THETA)*(1.0_dp+P01*TEMR)
  CT = 1.0_dp+UP/DN*TEMR
  F_VAL = 0.66666667_dp/THETA/SQRT(THETA)
  X1 = FERINV_FIT(F_VAL,1_i4)-1.5_dp*LOG(CT)
  CMU1 = TEMR*X1
  CHEMFIT = CMU1+1.0_dp

  THETALOG = LOG(THETA)
  THETACOR = (THETA32+0.29_dp)**0.666667_dp
  TB = EXP(-B1ICH*THETALOG)
  TB1 = TB/THETA
  DENOM = (1.0_dp+AICH*TB)
  UP = (AICH*TB1+BICH*SQRT(TB1))
  DMT = -1.5_dp- &
    (B1ICH+1.0_dp)*(AICH*TB1+BICH/2.0_dp*SQRT(TB1))/DENOM+ &
    UP*B1ICH*AICH*TB/DENOM**2
  DTP = -THETA/TEMR*XSR**2/SQRT(1.0_dp+XSR**2)
  DMT = DMT- &
    1.5_dp*THETA32*TEMR**2/ &
    ((1.0_dp+2.0_dp*TEMR)*THETACOR+TEMR**2)/(THETA32+0.27_dp)
  DDENR = 3.0_dp*DENR/(DMT*DTP)/TEMR
END SUBROUTINE CHEMP99_FIT

!> Inverse Fermi integral, H.M.Antia 1993 ApJS 84:101 fit (relative
!> error 0.01%) -- ported from `function FERINV`. The original computes
!> its fit-coefficient tables once via a `SAVE`d/`KRUN`-guarded
!> assignment block (a runtime-once-init idiom); ported here as compile-
!> time PARAMETER arrays instead -- numerically identical, no runtime
!> guard needed since F90 PARAMETERs don't have the F77 DATA-vs-
!> assignment distinction that idiom worked around.
FUNCTION FERINV_FIT(F, N) RESULT(X_OUT)
  REAL(KIND=dp),    INTENT(IN) :: F
  INTEGER(KIND=i4), INTENT(IN) :: N
  REAL(KIND=dp) :: X_OUT

  REAL(KIND=dp), PARAMETER :: A0(0:3) = [785.16685_dp, 44.593646_dp, 35.954549_dp, 213.89693_dp]
  REAL(KIND=dp), PARAMETER :: A1(0:3) = [-140.34065_dp, 11.288764_dp, 13.90891_dp, 35.399035_dp]
  REAL(KIND=dp), PARAMETER :: A2(0:3) = [13.257418_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  REAL(KIND=dp), PARAMETER :: A3(0:3) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  REAL(KIND=dp), PARAMETER :: B0(0:3) = [1391.7278_dp, 39.519346_dp, 47.795853_dp, 710.85455_dp]
  REAL(KIND=dp), PARAMETER :: B1(0:3) = [-804.63066_dp, -5.7517464_dp, 12.133628_dp, 98.73747_dp]
  REAL(KIND=dp), PARAMETER :: B2(0:3) = [158.54806_dp, 0.26594291_dp, -0.23975074_dp, 1.0677555_dp]
  REAL(KIND=dp), PARAMETER :: B3(0:3) = [-10.640712_dp, 0.0_dp, 0.0_dp, -0.011827987_dp]
  REAL(KIND=dp), PARAMETER :: C0(0:3) = [0.0089742174_dp, 34.873722_dp, -0.98934493_dp, -0.51891788_dp]
  REAL(KIND=dp), PARAMETER :: C1(0:3) = [-0.10604768_dp, -26.922515_dp, 0.090731169_dp, -0.0091723019_dp]
  REAL(KIND=dp), PARAMETER :: D0(0:3) = [0.035898124_dp, 26.612832_dp, -0.68577484_dp, -0.36278896_dp]
  REAL(KIND=dp), PARAMETER :: D1(0:3) = [-0.42520975_dp, -20.452930_dp, 0.063338994_dp, -0.0061502672_dp]
  REAL(KIND=dp), PARAMETER :: D2(0:3) = [3.6612154_dp, 11.808945_dp, -0.1163584_dp, -0.03367354_dp]

  REAL(KIND=dp) :: T_VAL, UP, DOWN

  IF (N < 0_i4 .OR. N > 3_i4) THEN
    WRITE(*,'(A)') 'FERINV_FIT: Invalid subscript'; STOP 1
  END IF

  IF (F < 4.0_dp) THEN
    T_VAL = F
    UP = A0(N)+T_VAL*(A1(N)+T_VAL*(A2(N)+T_VAL*A3(N)))
    DOWN = B0(N)+T_VAL*(B1(N)+T_VAL*(B2(N)+T_VAL*B3(N)))
    X_OUT = LOG(F*UP/DOWN)
  ELSE
    IF (N == 0_i4) THEN
      T_VAL = 1.0_dp/F**2
    ELSE
      T_VAL = EXP(-LOG(F)/(0.5_dp+REAL(N,dp)))
    END IF
    UP = C0(N)+T_VAL*(C1(N)+T_VAL)
    DOWN = D0(N)+T_VAL*(D1(N)+T_VAL*D2(N))
    X_OUT = UP/DOWN/T_VAL
  END IF
END FUNCTION FERINV_FIT

!> exp(xi)*E_{L+1}(xi) exponential integral -- ported from `function
!> EXPINT` (the `conduct21.f` one, taking an order argument L; NOT the
!> same routine as this project's other, already-existing, one-argument
!> `EXPINT`-named helper in `ode_integrator.f90`/elsewhere -- named
!> distinctly, `EXPINT_FIT`, to avoid any implied relationship). Two
!> dead labels in the original (11, 21 -- confirmed unused via an
!> independent gfortran warning on the unmodified source) are dropped in
!> favor of plain DO loops.
FUNCTION EXPINT_FIT(XI, L) RESULT(Q0_OUT)
  REAL(KIND=dp),    INTENT(IN) :: XI
  INTEGER(KIND=i4), INTENT(IN) :: L
  REAL(KIND=dp) :: Q0_OUT

  REAL(KIND=dp), PARAMETER :: EULER_G = 0.5772156649_dp
  INTEGER(KIND=i4), PARAMETER :: NREP = 21

  REAL(KIND=dp) :: CL, CI, C_VAL, PSI, CMX, CM, DQ, Q0
  INTEGER(KIND=i4) :: I, K, M

  IF (XI >= 1.0_dp) THEN
    CL = REAL(L,dp)
    CI = REAL(NREP,dp)
    C_VAL = 0.0_dp
    DO I = NREP, 1, -1
      C_VAL = CI/(XI+C_VAL)
      C_VAL = (CL+CI)/(1.0_dp+C_VAL)
      CI = CI-1.0_dp
    END DO
    Q0 = 1.0_dp/(XI+C_VAL)
  ELSE
    PSI = -EULER_G
    DO K = 1, L
      PSI = PSI+1.0_dp/REAL(K,dp)
    END DO
    Q0 = 0.0_dp
    CMX = 1.0_dp
    CL = REAL(L,dp)
    CM = -1.0_dp
    DO M = 0, NREP
      CM = CM+1.0_dp
      IF (M /= 0) CMX = -CMX*XI/CM
      IF (M /= L) THEN
        DQ = CMX/(CM-CL)
      ELSE
        DQ = CMX*(LOG(XI+1.0E-20_dp)-PSI)
      END IF
      Q0 = Q0-DQ
    END DO
    Q0 = EXP(XI)*Q0
  END IF
  Q0_OUT = Q0
END FUNCTION EXPINT_FIT

!> Electron-electron collision relaxation time, Shternin & Yakovlev
!> (2006), corrected in the nondegenerate regime to match Lampe (1968)
!> -- ported from `subroutine TAUEESY`.
!>
!> @warning The original's own `Y>1e8` branch has a variable-name
!>   mismatch bug (`CITL` computed, `CILT` used) -- see module header.
!>   Ported to match the original's ACTUAL behavior (CILT=0 in that
!>   branch), not the evidently-intended one.
SUBROUTINE TAUEESY(X_VAL, TEMR, TAUEE)
  REAL(KIND=dp), INTENT(IN)  :: X_VAL, TEMR
  REAL(KIND=dp), INTENT(OUT) :: TAUEE

  REAL(KIND=dp) :: E_VAL, V_VAL, Y_VAL, CIL, CIT, CILT, CITL
  REAL(KIND=dp) :: C1_VAL, C2_VAL, A_VAL, C_VAL, YV, FI, FREQ
  REAL(KIND=dp) :: THETA, T_VAL

  E_VAL = SQRT(1.0_dp+X_VAL**2)
  V_VAL = X_VAL/E_VAL
  Y_VAL = 0.0963913_dp/TEMR*X_VAL*SQRT(V_VAL)
  IF (Y_VAL > 1.0E8_dp) THEN
    CIL = 20.4013123_dp/(V_VAL*Y_VAL**3)
    CIT = 2.404_dp*V_VAL/Y_VAL**2
    ! Original assigns this term to `CITL`, but the FI sum below reads
    ! `CILT` -- a different variable, only ever assigned in the ELSE
    ! branch. Ported faithfully: CITL computed (matching the original,
    ! otherwise unused), CILT left at its effectively-zero
    ! uninitialized-SAVE value for this branch. See module header.
    CITL = 18.52_dp*(V_VAL/Y_VAL**8)**(1.0_dp/3.0_dp)
    ASSOCIATE (UNUSED_CITL => CITL); END ASSOCIATE
    CILT = 0.0_dp
  ELSE
    C1_VAL = 0.123636_dp+0.016234_dp*V_VAL**2
    C2_VAL = 0.0762_dp+0.05714_dp*V_VAL**4
    A_VAL = 12.2_dp+25.2_dp*V_VAL**3
    C_VAL = A_VAL*EXP(C1_VAL/C2_VAL)
    YV = Y_VAL*V_VAL
    CIL = LOG(1.0_dp+128.56_dp/(37.1_dp*Y_VAL+10.83_dp*Y_VAL**2+Y_VAL**3))* &
      (0.1587_dp-0.02538_dp/(1.0_dp+0.0435_dp*Y_VAL))/V_VAL
    CIT = V_VAL**3*LOG(1.0_dp+C_VAL/(A_VAL*YV+YV**2))* &
      (2.404_dp/C_VAL+(C2_VAL-2.404_dp/C_VAL)/(1.0_dp+0.1_dp*YV))
    CILT = V_VAL*LOG(1.0_dp+C_VAL/(A_VAL*Y_VAL+10.83_dp*YV**2+YV**(8.0_dp/3.0_dp)))* &
      (18.52_dp*V_VAL**2/C_VAL+(C2_VAL-18.52_dp*V_VAL**2/C_VAL)/ &
      (1.0_dp+0.1558_dp*Y_VAL**(1.0_dp-0.75_dp*V_VAL)))
  END IF
  FI = CIL+CIT+CILT
  FREQ = 0.00021381_dp*X_VAL*Y_VAL*SQRT(V_VAL)*FI
  IF (X_VAL > 0.001_dp) THEN
    THETA = TEMR/(E_VAL-1.0_dp)
  ELSE
    THETA = 2.0_dp*TEMR/X_VAL**2
  END IF
  T_VAL = 25.0_dp*THETA
  TAUEE = (1.0_dp+T_VAL+0.4342_dp*SQRT(THETA)*T_VAL**2)/(1.0_dp+T_VAL**2)/FREQ
END SUBROUTINE TAUEESY

!> Density-regime dispatcher -- previously (2026-08-23 to 2026-08-23)
!> split NSCool's GYP/PBHY formulas at a hardcoded 6e7 g/cm**3 threshold;
!> now a thin wrapper with NO density threshold anywhere in it: `OYAFORM`
!> for composition (unchanged), `CONDUCT_TRANSPORT` for the single,
!> un-spliced conductivity formula that removes the seam entirely (see
!> module header).
!>
!> @warning `NU_E_S`/`NU_E_L` (effective collision frequencies) are no
!>   longer part of this signature -- confirmed by reading
!>   `ETA_AND_F_HALL_AT` directly, they were ALREADY unused there (wrapped
!>   in a dead ASSOCIATE block). `CONDUCT_TRANSPORT` doesn't naturally
!>   produce a directly-analogous pair of numbers; carrying forward
!>   fabricated/derived values would misrepresent the new formula's
!>   actual outputs, so they're dropped rather than faked.
SUBROUTINE CON_CRUST(T, RHO, SIGMA, LAMBDA, N_E)
  REAL(KIND=dp), INTENT(IN)  :: T, RHO
  REAL(KIND=dp), INTENT(OUT) :: SIGMA, LAMBDA, N_E

  REAL(KIND=dp), PARAMETER :: MU = 1.66E-24_dp
  ! Standard neutron-drip density, g/cm**3 -- same value/rationale as
  ! this project's prior port's identical constant (see git history).
  REAL(KIND=dp), PARAMETER :: RHODRIP_CGS = 4.3E11_dp

  REAL(KIND=dp) :: BARD, Z, A1, A, XNUC, XNUCT
  INTEGER(KIND=i4) :: INDEX_PHASE
  REAL(KIND=dp) :: N_I
  REAL(KIND=dp) :: QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH

  BARD = RHO/MU * 1.0E-39_dp
  IF (RHO > RHODRIP_CGS) THEN
    INDEX_PHASE = 3
  ELSE
    INDEX_PHASE = 30
  END IF
  CALL OYAFORM(BARD, INDEX_PHASE, Z, A1, A, XNUC, XNUCT)
  ASSOCIATE (UNUSED_A1 => A1); END ASSOCIATE

  N_I = RHO/(A*MU)
  N_E = Z*N_I

  CALL CONDUCT_TRANSPORT(T, RHO, 0.0_dp, Z, A, 0.0_dp, XNUC, XNUCT, &
    SIGMA, LAMBDA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH)
  ASSOCIATE (UNUSED_QJ => QJ, UNUSED_SIGMAT => SIGMAT, UNUSED_CKAPPAT => CKAPPAT, &
    UNUSED_QJT => QJT, UNUSED_SIGMAH => SIGMAH, UNUSED_CKAPPAH => CKAPPAH, UNUSED_QJH => QJH)
  END ASSOCIATE
END SUBROUTINE CON_CRUST

!> Self-consistent, radially-varying resistivity and Hall pre-factor
!> from a solved TOV crust profile -- new code (not a port), combining
!> CON_CRUST's SIGMA/N_E outputs into `eta = c**2/(4*pi*sigma)` (standard
!> magnetic diffusivity from conductivity) and `f_H = c/(4*pi*e*n_e)`
!> (already established in hall_induction.f90's own docstring), converted
!> from Gaussian-cgs into this project's code units via UNITS.
!>
!> Only rows with RHOCGS below the core-crust boundary
!> (`RHOL_CGS=2.2e14 g/cm**3`, matching this project's prior ports'
!> identical threshold) are supported -- `CON_CRUST` is meant for crust
!> densities; core rows are out of scope here as before.
!>
!> @param PROFILE Solved TOV radial profile (TOV_SOLVER::SOLVE_TOV_STAR).
!> @param ETA_PROFILE Resistivity per radial row, code units (out, allocated here).
!> @param F_HALL_PROFILE Hall pre-factor per radial row, code units (out, allocated here).
!> @param N_E_PROFILE Electron number density per radial row, `cm**-3`
!>   (out, allocated here) -- CON_CRUST's own N_E output, threaded
!>   through since it's already computed per row and a driver reporting
!>   eta(r)/f_H(r) will typically want n_e(r) alongside them.
!> @param T_KELVIN Crust temperature, K (optional; default 1e9 K,
!>   matching this project's prior ports' own default).
SUBROUTINE ETA_AND_F_HALL_AT(PROFILE, ETA_PROFILE, F_HALL_PROFILE, N_E_PROFILE, T_KELVIN)
  TYPE(TOV_PROFILE_T), INTENT(IN)  :: PROFILE
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: ETA_PROFILE(:), F_HALL_PROFILE(:), N_E_PROFILE(:)
  REAL(KIND=dp), OPTIONAL, INTENT(IN) :: T_KELVIN

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp), PARAMETER :: C_LIGHT_CGS = 2.99792458E10_dp   ! cm/s
  REAL(KIND=dp), PARAMETER :: E_CHARGE_ESU = 4.80320425E-10_dp ! esu
  REAL(KIND=dp), PARAMETER :: RHOL_CGS = 2.2E14_dp             ! core-crust boundary

  REAL(KIND=dp) :: T
  REAL(KIND=dp) :: SIGMA, LAMBDA_TH, N_E_CGS
  REAL(KIND=dp) :: ETA_CGS, F_HALL_CGS
  INTEGER(KIND=i4) :: I

  T = 1.0E9_dp; IF (PRESENT(T_KELVIN)) T = T_KELVIN

  ALLOCATE(ETA_PROFILE(PROFILE%N), F_HALL_PROFILE(PROFILE%N), N_E_PROFILE(PROFILE%N))

  DO I = 1, PROFILE%N
    IF (PROFILE%RHOCGS(I) > RHOL_CGS) THEN
      WRITE(*,'(A)') 'ETA_AND_F_HALL_AT: core-density row (rho > 2.2e14 g/cm**3) not supported'
      STOP 1
    END IF

    CALL CON_CRUST(T, PROFILE%RHOCGS(I), SIGMA, LAMBDA_TH, N_E_CGS)
    ASSOCIATE (UNUSED_LAMBDA => LAMBDA_TH); END ASSOCIATE

    ETA_CGS    = C_LIGHT_CGS**2 / (4.0_dp*PI*SIGMA)                        ! cm**2/s
    F_HALL_CGS = C_LIGHT_CGS / (4.0_dp*PI*E_CHARGE_ESU*N_E_CGS)            ! cm**2/(G*s)

    ETA_PROFILE(I)    = ETA_CGS * TIME_UNIT_S / LENGTH_UNIT_CM**2
    F_HALL_PROFILE(I) = F_HALL_CGS * B_UNIT_GAUSS * TIME_UNIT_S / LENGTH_UNIT_CM**2
    N_E_PROFILE(I)    = N_E_CGS
  END DO
END SUBROUTINE ETA_AND_F_HALL_AT

END MODULE CRUST_CONDUCTIVITY
