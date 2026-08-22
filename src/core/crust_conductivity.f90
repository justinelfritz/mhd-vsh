MODULE CRUST_CONDUCTIVITY
!> Electron transport coefficients (thermal conductivity, relaxation
!> time) in strongly degenerate, magnetized crustal matter -- mechanical
!> port of ~/Desktop/EOSNS/src/potekhinc.f (Potekhin 1999, A&A 351, 787;
!> nuclear-form-factor extension per Gnedin et al. 2001, MNRAS 324, 725).
!>
!> @warning This is dense, published fitting-formula physics -- every
!>   numerical constant (Coulomb-logarithm fit coefficients, magnetic
!>   quantization corrections, exponential-integral series parameters)
!>   is copied byte-for-byte from the F77 original. Do not "simplify" or
!>   re-derive any of it; correctness here means matching the original
!>   code's own output (see test_crust_conductivity.f90, which checks
!>   against the unmodified original F77 source directly), not
!>   re-deriving the physics from the cited papers.
!>
!> @warning One deliberate exception to "byte-for-byte": every literal
!>   here is written with an explicit `_dp` suffix (full double-precision
!>   parsing), whereas many of the original's bare literals (e.g.
!>   `DATA AUM/1822.9/,AUD/15819.4/`, `DATA BOHR/137.036/` -- no `d0`)
!>   are silently parsed as single precision (~7 significant digits)
!>   before being widened to double, a well-known F77 gotcha. This makes
!>   this port's arithmetic MORE precise than the original by ~1e-6
!>   relative per constant, not less faithful to the underlying
!>   Potekhin-1999 physics -- confirmed (not assumed) by rebuilding the
!>   original with `-fdefault-real-8` (forcing double-precision literal
!>   parsing), which reproduces this port's output bit-for-bit.
!>   test_crust_conductivity.f90's tolerance (5e-6) is set to account
!>   for exactly this quantified, understood effect.
!>
!> @warning One line of the original (`potekhinc`'s SIGMA/SIGMAT/SIGMAH
!>   cgs-conversion lines, from `RSIGMA`/`RTSIGMA`/`RHSIGMA`) is *not*
!>   ported: those three source variables are local to `CONDEGINc`
!>   (never part of its argument list), so in the original F77,
!>   `potekhinc`'s read of them is uninitialized memory -- and the
!>   result (`SIGMA`/`SIGMAT`/`SIGMAH`) is itself never part of
!>   `potekhinc`'s own returned arguments, so this dead, already-buggy
!>   code has zero effect on any real output. Omitted rather than
!>   faithfully reproduced as a bug.
!>
!> Internal quantities are in relativistic units (\hbar=m_e=c=1) except
!> where the entry point `POTEKHINC` converts to cgs, exactly as the
!> original documents.
USE KINDS,      ONLY: dp, i4
USE TOV_SOLVER, ONLY: TOV_PROFILE_T
USE UNITS,      ONLY: B_UNIT_GAUSS, LENGTH_UNIT_CM, TIME_UNIT_S
IMPLICIT NONE
PRIVATE
PUBLIC :: POTEKHINC, ETA_AND_F_HALL_AT

CONTAINS

!> Self-consistent, radially-varying resistivity and Hall pre-factor from
!> a solved TOV crust profile -- new code (not a port), combining
!> POTEKHINC's TAU with the profile's own n_e(r), via nstot.f's *exact*
!> `sigmae = 3.26*tau*nel/meff` combination (nstot.f:249; the same
!> formula that produced the validated fort.34/PL.DAT reference output),
!> so this reproduces the physics nstot.f itself already uses, not a
!> re-derivation. `meff` (nstot.f:236-237, `kfe=197.33*(nel*3*pi**2)**
!> (1/3)`, `meff=sqrt(0.511**2+kfe**2)/0.511`) is likewise ported
!> byte-for-byte -- it is not one of POTEKHINC's own outputs.
!>
!> `eta = c**2/(4*pi*sigma_e)` (standard magnetic diffusivity from
!> conductivity) and `f_H = c/(4*pi*e*n_e)` (already established in
!> hall_induction.f90's own docstring) are then converted from Gaussian-
!> cgs into this project's code units via UNITS -- eta needs
!> [length]**2/[time], f_H needs [length]**2/([B]*[time]), matching
!> DIFFUSION_INIT/HALL_INIT's existing scalar ETA/F_HALL convention.
!>
!> Only rows with a genuine nucleus (PROFILE%A_TABLE>0 -- NOT AH, which
!> COMPOSITION always sets >0 via a dummy substitute for homogeneous
!> rows) are supported -- deep crust/core "homogeneous matter" rows are
!> where nstot.f itself skips calling potekhinc (see nstot.f's own
!> `IF (a.gt.0.d0)` gate around its conductivity block, tested against
!> the raw table value, not the derived ah); this routine is meant for
!> the outer-crust R_MIN/R_MAX range MHD-VSH actually simulates, where
!> that always holds.
!>
!> @param PROFILE Solved TOV radial profile (TOV_SOLVER::SOLVE_TOV_STAR).
!> @param ETA_PROFILE Resistivity per radial row, code units (out, allocated here).
!> @param F_HALL_PROFILE Hall pre-factor per radial row, code units (out, allocated here).
!> @param T_KELVIN Crust temperature, K (optional; default 1e9 K, matching
!>   nstot.f's own hardcoded isothermal-crust value).
!> @param B12_GAUSS Magnetic field magnitude fed to POTEKHINC's own
!>   quantization/scattering physics, 1e12 G (optional; default 10.0,
!>   matching nstot.f's own hardcoded value -- this is POTEKHINC's input
!>   field scale, independent of the MHD-VSH-simulated field itself).
!> @param XIMP Impurity parameter (optional; default 0.1, matching
!>   nstot.f's own hardcoded value).
SUBROUTINE ETA_AND_F_HALL_AT(PROFILE, ETA_PROFILE, F_HALL_PROFILE, T_KELVIN, B12_GAUSS, XIMP)
  TYPE(TOV_PROFILE_T), INTENT(IN)  :: PROFILE
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: ETA_PROFILE(:), F_HALL_PROFILE(:)
  REAL(KIND=dp), OPTIONAL, INTENT(IN) :: T_KELVIN, B12_GAUSS, XIMP

  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp), PARAMETER :: C_LIGHT_CGS = 2.99792458E10_dp   ! cm/s
  REAL(KIND=dp), PARAMETER :: E_CHARGE_ESU = 4.80320425E-10_dp ! esu

  REAL(KIND=dp) :: T, B12, XIMP_LOCAL
  REAL(KIND=dp) :: TCOND, TCONDT, TCONDH, TAU
  REAL(KIND=dp) :: NEL_I, KFE, MEFF_I, SIGMAE_REL, SIGMAE_CGS, NE_CGS
  REAL(KIND=dp) :: ETA_CGS, F_HALL_CGS
  INTEGER(KIND=i4) :: I

  T = 1.0E9_dp;  IF (PRESENT(T_KELVIN))  T = T_KELVIN
  B12 = 10.0_dp; IF (PRESENT(B12_GAUSS)) B12 = B12_GAUSS
  XIMP_LOCAL = 0.1_dp; IF (PRESENT(XIMP)) XIMP_LOCAL = XIMP

  ALLOCATE(ETA_PROFILE(PROFILE%N), F_HALL_PROFILE(PROFILE%N))

  DO I = 1, PROFILE%N
    IF (PROFILE%A_TABLE(I) <= 0.0_dp) THEN
      WRITE(*,'(A)') 'ETA_AND_F_HALL_AT: homogeneous-matter row (no nucleus) not supported'
      STOP 1
    END IF

    CALL POTEKHINC(T, PROFILE%RHOCGS(I), B12, PROFILE%AH(I), PROFILE%ZH(I), &
      PROFILE%XH(I), XIMP_LOCAL, TCOND, TCONDT, TCONDH, TAU)
    ASSOCIATE (UNUSED_TCOND => TCOND, UNUSED_TCONDT => TCONDT, UNUSED_TCONDH => TCONDH)
    END ASSOCIATE

    NEL_I  = PROFILE%NEL(I)   ! fm**-3
    KFE    = 197.33_dp*(NEL_I*3.0_dp*PI**2)**(1.0_dp/3.0_dp)
    MEFF_I = SQRT(0.511_dp**2 + KFE**2)/0.511_dp

    SIGMAE_REL = 3.26_dp*TAU*NEL_I/MEFF_I   ! nstot.f:249, units of 1e26 s**-1
    SIGMAE_CGS = SIGMAE_REL * 1.0E26_dp     ! s**-1
    NE_CGS     = NEL_I * 1.0E39_dp          ! fm**-3 -> cm**-3

    ETA_CGS    = C_LIGHT_CGS**2 / (4.0_dp*PI*SIGMAE_CGS)              ! cm**2/s
    F_HALL_CGS = C_LIGHT_CGS / (4.0_dp*PI*E_CHARGE_ESU*NE_CGS)        ! cm**2/(G*s)

    ETA_PROFILE(I)    = ETA_CGS * TIME_UNIT_S / LENGTH_UNIT_CM**2
    F_HALL_PROFILE(I) = F_HALL_CGS * B_UNIT_GAUSS * TIME_UNIT_S / LENGTH_UNIT_CM**2
  END DO
END SUBROUTINE ETA_AND_F_HALL_AT

!> Entry point. Ported from potekhinc.f's `subroutine potekhinc` verbatim
!> (see module header re: the omitted dead SIGMA lines).
!> @param T1 Temperature, K.
!> @param RHO Mass density, g/cm**3.
!> @param B12 Magnetic field, 1e12 G.
!> @param CMI Ion mass number.
!> @param ZION Ion charge number.
!> @param XH Hydrogen-like mass fraction (nstot.f's usage: mass fraction
!>   of the nuclear species, `CMI1=CMI/XH` is nucleons-per-nucleus).
!> @param ZIMP Impurity parameter (effective Z, `Z_imp**2 = <n_j(Z-Z_j)**2>/n`).
!> @param CKAPPA, CKAPPAT, CKAPPAH Longitudinal/transverse/Hall thermal
!>   conductivity, erg/(K cm s).
!> @param TAU Longitudinal relaxation time (relativistic units).
SUBROUTINE POTEKHINC(T1, RHO, B12, CMI, ZION, XH, ZIMP, CKAPPA, CKAPPAT, CKAPPAH, TAU)
  REAL(KIND=dp), INTENT(IN)  :: T1, RHO, B12, CMI, ZION, XH, ZIMP
  REAL(KIND=dp), INTENT(OUT) :: CKAPPA, CKAPPAT, CKAPPAH, TAU
  ! AUM: atomic mass unit / electron mass. AUD: relativistic unit of
  ! density, g/cm**3 (used here only as the product AUD*AUM). DRIP
  ! (neutron drip density, g/cm**3) is documented in the original but
  ! never actually used in its arithmetic -- not ported.
  REAL(KIND=dp), PARAMETER :: AUM = 1822.9_dp, AUD = 15819.4_dp
  REAL(KIND=dp), PARAMETER :: UNIKAP = 2.778E15_dp
  REAL(KIND=dp) :: TEMP, B, CMI1, DENSI
  REAL(KIND=dp) :: RKAPPA, RTKAPPA, RHKAPPA

  TEMP = (T1/1.0E6_dp)/5930.0_dp
  B = B12/44.14_dp
  CMI1 = CMI/XH
  DENSI = RHO/(AUD*AUM*CMI1)

  CALL CONDEGINC(TEMP, DENSI, B, ZION, CMI, CMI1, ZIMP, RKAPPA, RTKAPPA, RHKAPPA, TAU)

  CKAPPA  = RKAPPA*UNIKAP
  CKAPPAT = RTKAPPA*UNIKAP
  CKAPPAH = RHKAPPA*UNIKAP
END SUBROUTINE POTEKHINC

!> Central transport-coefficient calculation. Ported from potekhinc.f's
!> `subroutine CONDEGINc` verbatim.
SUBROUTINE CONDEGINC(TEMP, DENSI, B, ZION, CMI, CMI1, ZIMP, RKAPPA, RTKAPPA, RHKAPPA, TAULONGT)
  REAL(KIND=dp), INTENT(IN)  :: TEMP, DENSI, B, ZION, CMI, CMI1, ZIMP
  REAL(KIND=dp), INTENT(OUT) :: RKAPPA, RTKAPPA, RHKAPPA, TAULONGT
  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp, BOHR = 137.036_dp
  REAL(KIND=dp) :: DENS, SPHERION, GAMMA, XSR, EF0, Q2E0, CST, XNUC
  REAL(KIND=dp) :: PCL, Q2E, CLEFF, CLLONG, CLTRAN, SN, THTOEL
  REAL(KIND=dp) :: E, TAU, GYROM, TAUT0
  REAL(KIND=dp) :: CLEFFI, CLLONGI, CLTRANI, SNI
  REAL(KIND=dp) :: TAULONG, TAUT, TAUTRAN, TAUHALL
  REAL(KIND=dp) :: TAUTT, TAUTRANT, TAUHALLT
  REAL(KIND=dp) :: C, CTH, TAUEE, EECOR

  IF (TEMP<=0.0_dp .OR. DENSI<=0.0_dp .OR. B<0.0_dp .OR. ZION<=0.0_dp .OR. CMI<=0.0_dp) THEN
    WRITE(*,'(A)') 'CONDEGIN: Non-positive input parameter'; STOP 1
  END IF
  IF (CMI1 < CMI) THEN
    WRITE(*,'(A)') 'CONDEGIN: Incorrect CMI1'; STOP 1
  END IF
  IF (ZION < 0.5_dp) THEN
    WRITE(*,'(A)') 'CONDEGIN: Too small ion charge'; STOP 1
  END IF
  IF (CMI < 1.0_dp) THEN
    WRITE(*,'(A)') 'CONDEGIN: Too small ion mass'; STOP 1
  END IF
  IF (DENSI > 1.0E6_dp) THEN
    WRITE(*,'(A)') 'CNDEGIN: Too high density'; STOP 1
  END IF

  DENS = DENSI*ZION
  SPHERION = (0.75_dp/PI/DENSI)**0.3333333_dp
  GAMMA = ZION**2/BOHR/TEMP/SPHERION
  XSR = (3.0_dp*PI**2*DENS)**0.3333333_dp
  EF0 = SQRT(1.0_dp + XSR**2)
  Q2E0 = 4.0_dp/PI/BOHR*XSR*EF0
  CST = XSR**3/0.75_dp/PI*ZION/BOHR**2   ! = 4*pi*n_i*(Z*e**2)**2
  IF (CMI == CMI1) THEN
    XNUC = 0.00155_dp*(CMI/ZION)**0.33333_dp*XSR
  ELSE
    XNUC = 0.00247_dp*XSR
  END IF

  ! Non-quantizing case only (the original's quantizing-field branch is
  ! commented out in the source itself -- ported as inactive there too).
  PCL = XSR
  Q2E = Q2E0

  CALL COULIN(PCL, XSR, GAMMA, B, ZION, CMI, Q2E, XNUC, CLEFF, CLLONG, CLTRAN, SN, THTOEL)
  E = SQRT(1.0_dp + PCL**2)
  TAU = PCL**3/E/4.0_dp/PI/DENSI/(ZION/BOHR)**2/CLLONG/SN
  GYROM = B/E
  TAUT0 = PCL**3/E/CST/CLTRAN*SN

  CLEFFI = 0.0_dp; CLLONGI = 0.0_dp; CLTRANI = 0.0_dp; SNI = 1.0_dp
  IF (ZIMP > 0.0_dp) CALL COUL99I(PCL, XSR, GAMMA, B, Q2E, CLEFFI, CLLONGI, CLTRANI, SNI)
  ! CLEFFI and SNI (from COUL99I) are not read again -- only CLLONGI/
  ! CLTRANI feed the formulas below, matching the original exactly.

  TAULONG = TAU*CLLONG*ZION**2 / (CLLONG*ZION**2 + CLLONGI*ZIMP**2)
  TAUT    = TAUT0*CLTRAN*ZION**2 / (CLTRAN*ZION**2 + CLTRANI*ZIMP**2)
  TAUTRAN = TAUT / (1.0_dp + (TAUT*GYROM)**2)
  TAUHALL = TAUT*GYROM*TAUTRAN
  ! In the original, TAUHALL's only consumer is RHSIGMA=C*TAUHALL --
  ! part of the dead RSIGMA/RTSIGMA/RHSIGMA computation this port
  ! omits (see module header), so TAUHALL is computed here (for
  ! fidelity to the original's own derivation chain) but not read
  ! further, unlike in the original where it technically was (by dead
  ! code with no effect on any real output).

  TAULONGT = TAU*CLLONG*ZION**2 / (CLLONG*ZION**2*THTOEL + CLLONGI*ZIMP**2)
  TAUTT    = TAUT0*CLTRAN*ZION**2 / (CLTRAN*ZION**2*THTOEL + CLTRANI*ZIMP**2)
  TAUTRANT = TAUTT / (1.0_dp + (TAUTT*GYROM)**2)
  TAUHALLT = TAUTT*GYROM*TAUTRANT

  C = SN*PCL**3/3.0_dp/PI**2/E/BOHR
  ! RSIGMA/RTSIGMA/RHSIGMA (=C*TAULONG/TAUTRAN/TAUHALL) intentionally
  ! not computed here -- see module header, they're dead in the caller.
  CTH = C*PI**2*TEMP/3.0_dp*BOHR
  TAUEE = RELAXEE(XSR, DENS, TEMP)
  EECOR = TAUEE/(TAULONGT+TAUEE)
  RKAPPA  = CTH*TAULONGT*EECOR
  RTKAPPA = CTH*TAUTRANT*EECOR
  RHKAPPA = CTH*TAUHALLT*EECOR
END SUBROUTINE CONDEGINC

!> Effective Coulomb logarithm (electron-phonon/ion scattering), ported
!> from potekhinc.f's `subroutine COULIN` verbatim.
SUBROUTINE COULIN(PCL, XSR, GAMMA, B, ZION, CMI, Q2E, XNUC, CLEFF, CLLONG, CLTRAN, SN, THTOEL)
  REAL(KIND=dp), INTENT(IN)  :: PCL, XSR, GAMMA, B, ZION, CMI, Q2E, XNUC
  REAL(KIND=dp), INTENT(OUT) :: CLEFF, CLLONG, CLTRAN, SN, THTOEL
  REAL(KIND=dp), PARAMETER :: UMINUS1 = 2.78_dp, UMINUS2 = 12.973_dp
  REAL(KIND=dp), PARAMETER :: AUM = 1822.9_dp, BOHR = 137.036_dp, PI = 3.14159265_dp

  REAL(KIND=dp) :: DENS, DENSI, SPHERION, Q2ICL, ECL, VCL, PM2, TRP, BORNCOR
  REAL(KIND=dp) :: C, Q2S, XS, R2W, XW, XW1, CL, A0, VIBRCOR, T0, G0, GW
  REAL(KIND=dp) :: G2, TRU, CLHIGH, EU, CLLOWK, CLLOWS
  REAL(KIND=dp) :: ENU, PB
  REAL(KIND=dp) :: XIS, ZETA, XI, XSUM, Q2M, QTRANM, QTRANP, Q
  REAL(KIND=dp) :: DNU, XS1, PN, SQB, X, EXW, A1, Q1, DLT, Y1, CL0, P2, Y2, PY, DT
  REAL(KIND=dp) :: DB, CL1, P1, P3, P4
  INTEGER(KIND=i4) :: NL, N

  DENS = XSR**3/3.0_dp/PI**2
  DENSI = DENS/ZION
  SPHERION = (0.75_dp/PI/DENSI)**0.3333333_dp
  Q2ICL = 3.0_dp*GAMMA/SPHERION**2
  ECL = SQRT(1.0_dp+PCL**2)
  VCL = PCL/ECL
  PM2 = (2.0_dp*PCL)**2
  TRP = ZION/GAMMA*SQRT(CMI*AUM*SPHERION/3.0_dp/BOHR)
  BORNCOR = VCL*ZION*PI/BOHR

  C = (1.0_dp+0.06_dp*GAMMA)*EXP(-SQRT(GAMMA))
  Q2S = (Q2ICL*C+Q2E)*EXP(-BORNCOR)
  XS = Q2S/PM2
  R2W = UMINUS2/Q2ICL*(1.0_dp+0.3333_dp*BORNCOR)
  XW = R2W*PM2
  XW1 = 14.7327_dp*XNUC**2
  XW1 = XW1*(1.0_dp+0.3333_dp*BORNCOR)*(1.0_dp+ZION/13.0_dp*SQRT(XNUC))
  CL = COULAN2(XS, XW, VCL, XW1)
  A0 = 1.683_dp*SQRT(PCL/CMI/ZION)
  VIBRCOR = EXP(-A0/4.0_dp*UMINUS1*EXP(-9.1_dp*TRP))
  T0 = 0.19_dp/ZION**0.16667_dp
  G0 = TRP/SQRT(TRP**2+T0**2)*(1.0_dp+(ZION/125.0_dp)**2)
  GW = G0*VIBRCOR
  CLEFF = CL*GW
  G2 = TRP/SQRT(0.0081_dp+TRP**2)**3
  THTOEL = 1.0_dp + G2/G0*(1.0_dp+BORNCOR*VCL**3)*0.0105_dp*(1.0_dp-1.0_dp/ZION)* &
           (1.0_dp+XNUC**2*SQRT(2.0_dp*ZION))
  TRU = TRP*3.0_dp*VCL*BOHR/ZION**0.3333333_dp
  IF (TRU < 20.0_dp) THEN
    CLHIGH = CLEFF
    EU = EXP(-1.0_dp/TRU)
    CLLOWK = 50.0_dp*SQRT(XSR/CMI)/ZION*TRP**3
    CLLOWS = CLLOWK/VCL/BOHR/0.75_dp*TRP**2
    CLEFF = CLHIGH*EU + CLLOWS*(1.0_dp-EU)
    THTOEL = (CLHIGH*THTOEL*EU + CLLOWK*(1.0_dp-EU))/CLEFF
  END IF

  ! CRITICAL: in the original, the guard `if (PCL**2.gt.4.d2*B) then`
  ! around this non-magnetic assignment is commented out, but the body
  ! (this assignment + `goto 50`, i.e. RETURN) is NOT commented -- so
  ! the original ALWAYS takes this path and jumps straight to its
  ! `50 return`, making the entire "Magnetic fit" section below
  ! genuinely dead/unreachable code in the actual, validated 2015
  ! reference run. Ported faithfully: RETURN here (matching `goto 50`),
  ! the magnetic-fit code below is kept (matching the original keeping
  ! it, just disabled) but is unreachable, not deleted, so it's
  ! available if magnetic quantization is ever deliberately re-enabled.
  CLLONG = CLEFF
  CLTRAN = CLEFF
  SN = 1.0_dp
  RETURN

  ENU = PCL**2/2.0_dp/B
  NL = INT(ENU)
  SN = 0.0_dp
  DO N = 0, NL
    PB = SQRT(ENU-REAL(N,KIND=dp))
    SN = SN + PB
    IF (N /= 0) SN = SN + PB
  END DO
  SN = SN*1.5_dp*B*SQRT(2.0_dp*B)/PCL**3

  IF (ENU <= 1.0_dp) THEN
    XIS = Q2S/2.0_dp/B
    ZETA = R2W*2.0_dp*B
    XI = 2.0_dp*PCL**2/B
    XSUM = XI+XIS
    Q2M = (EXPINT(XSUM,1) - EXP(-ZETA*XI)*EXPINT((1.0_dp+ZETA)*XSUM,1))/XSUM
    CLLONG = (PCL*VCL/B)**2*Q2M/1.5_dp*GW
    QTRANM = (1.0_dp+XSUM)*EXPINT(XSUM,0) - 1.0_dp - EXP(-ZETA*XI)* &
             ((1.0_dp+(1.0_dp+ZETA)*XSUM)*EXPINT((1.0_dp+ZETA)*XSUM,0) - 1.0_dp)
    QTRANP = (1.0_dp+XIS)*EXPINT(XIS,0) - 1.0_dp - &
             ((1.0_dp+(1.0_dp+ZETA)*XIS)*EXPINT((1.0_dp+ZETA)*XIS,0) - 1.0_dp)
    Q = (ECL**2*QTRANP + QTRANM)*B/PCL**2
    CLTRAN = 0.375_dp*Q/ECL**2*GW
  ELSE
    DNU = ENU - REAL(NL,KIND=dp)
    XS1 = (SQRT(XS) + 1.0_dp/(2.0_dp+XW/2.0_dp))**2
    PN = SQRT(2.0_dp*B*DNU)
    SQB = SQRT(B)
    X = MAX(PN/SQB, 1.0E-10_dp)
    IF (XW < 0.01_dp) THEN
      EXW = 1.0_dp
    ELSE IF (XW > 50.0_dp) THEN
      EXW = 1.0_dp/XW
    ELSE
      EXW = (1.0_dp-EXP(-XW))/XW
    END IF
    A1 = (30.0_dp-15.0_dp*EXW-(15.0_dp-6.0_dp*EXW)*VCL**2) / &
         (30.0_dp-10.0_dp*EXW-(20.0_dp-5.0_dp*EXW)*VCL**2)
    Q1 = 0.25_dp*VCL**2/(1.0_dp-0.667_dp*VCL**2)
    DLT = SQB/PCL*(A1/X - SQRT(X)*(1.5_dp-0.5_dp*EXW+Q1) + &
          (1.0_dp-EXW+0.75_dp*VCL**2)/(1.0_dp+VCL**2)*(X-SQRT(X))/REAL(NL,KIND=dp))
    Y1 = 1.0_dp/(1.0_dp+DLT)
    CL0 = LOG(1.0_dp+1.0_dp/XS1)
    P2 = CL0*(0.07_dp+0.2_dp*EXW)
    Y2 = 1.5_dp*CL0*(X**3-X/3.0_dp)/(REAL(NL,KIND=dp)+0.75_dp/(1.0_dp+2.0_dp*B)**2*X**2) + P2*X
    PY = 1.0_dp + 0.06_dp*CL0**2/REAL(NL,KIND=dp)**2
    DT = SQRT(PY*Y1**2+Y2**2)
    CLLONG = CLEFF/DT

    DB = 1.0_dp/(1.0_dp+0.5_dp/B)
    CL1 = XS1*CL0
    P1 = 0.8_dp*(1.0_dp+CL1) + 0.2_dp*CL0
    P2 = 1.42_dp - 0.1_dp*DB + SQRT(CL1)/3.0_dp
    P3 = (0.68_dp-0.13_dp*DB)*CL1**0.165_dp
    P4 = (0.52_dp-0.1_dp*DB)*SQRT(SQRT(CL1))
    DLT = SQB/PCL*(P1/X**2*SQB/PCL + P3*LOG(REAL(NL,KIND=dp))/X - &
          (P2+P4*LOG(REAL(NL,KIND=dp)))*SQRT(X))
    CLTRAN = CLEFF*(1.0_dp+DLT)
  END IF
END SUBROUTINE COULIN

!> Electron-impurity scattering Coulomb logarithm (simplified COULIN
!> variant), ported from potekhinc.f's `subroutine COUL99I` verbatim.
SUBROUTINE COUL99I(PCL, XSR, GAMMA, B, Q2E, CLEFF, CLLONG, CLTRAN, SN)
  REAL(KIND=dp), INTENT(IN)  :: PCL, XSR, GAMMA, B, Q2E
  REAL(KIND=dp), INTENT(OUT) :: CLEFF, CLLONG, CLTRAN, SN
  REAL(KIND=dp), PARAMETER :: PI = 3.14159265_dp
  REAL(KIND=dp), PARAMETER :: XW = 1.0E99_dp

  ! GAMMA is part of the original's call signature (matching COULIN's)
  ! but, unlike there, is never actually used in COUL99I's own body --
  ! ported faithfully as an unused argument, not a port error. Likewise
  ! DENS ("number density of electrons") is computed but never read
  ! again in the original; and C=(1+.06*GAMMA)*exp(-sqrt(GAMMA)) is
  ! computed in the original but bypassed (Q2S=Q2E directly, not via
  ! C) -- both omitted here as genuinely dead, matching RSIGMA/etc.
  ! (see module header), not silently dropped.
  REAL(KIND=dp) :: ECL, VCL, PM2, Q2S, XS
  REAL(KIND=dp) :: ENU, PB, XIS, XI, XSUM, Q2M, QTRANM, QTRANP, Q
  REAL(KIND=dp) :: DNU, XS1, PN, SQB, X, EXW, A1, Q1, DLT, Y1, CL0, P2, Y2, PY, DT
  REAL(KIND=dp) :: DB, CL1, P1, P3, P4
  INTEGER(KIND=i4) :: NL, N

  ECL = SQRT(1.0_dp+PCL**2)
  VCL = PCL/ECL
  PM2 = (2.0_dp*PCL)**2

  Q2S = Q2E
  XS = Q2S/PM2
  CLEFF = COULAN2(XS, XW, VCL, 0.0_dp)

  ! Same dead-code structure as COULIN (see its own detailed comment):
  ! the original's guard around this assignment is commented out, but
  ! the assignment + `goto 50` (RETURN) is not -- the "Magnetic fit"
  ! code below is genuinely unreachable in the validated reference.
  CLLONG = CLEFF
  CLTRAN = CLEFF
  SN = 1.0_dp
  RETURN

  ENU = PCL**2/2.0_dp/B
  NL = INT(ENU)
  SN = 0.0_dp
  DO N = 0, NL
    PB = SQRT(ENU-REAL(N,KIND=dp))
    SN = SN + PB
    IF (N /= 0) SN = SN + PB
  END DO
  SN = SN*1.5_dp*B*SQRT(2.0_dp*B)/PCL**3

  IF (ENU <= 1.0_dp) THEN
    XIS = Q2S/2.0_dp/B
    XI = 2.0_dp*PCL**2/B
    XSUM = XI+XIS
    Q2M = EXPINT(XSUM,1)
    CLLONG = (PCL*VCL/B)**2*Q2M/1.5_dp
    QTRANM = (1.0_dp+XSUM)*EXPINT(XSUM,0) - 1.0_dp
    QTRANP = (1.0_dp+XIS)*EXPINT(XIS,0) - 1.0_dp
    Q = (ECL**2*QTRANP + QTRANM)*B/PCL**2
    CLTRAN = 0.375_dp*Q/ECL**2
  ELSE
    DNU = ENU - REAL(NL,KIND=dp)
    XS1 = (SQRT(XS) + 1.0_dp/(2.0_dp+XW/2.0_dp))**2
    PN = SQRT(2.0_dp*B*DNU)
    SQB = SQRT(B)
    X = MAX(PN/SQB, 1.0E-10_dp)
    IF (XW < 0.01_dp) THEN
      EXW = 1.0_dp
    ELSE IF (XW > 50.0_dp) THEN
      EXW = 1.0_dp/XW
    ELSE
      EXW = (1.0_dp-EXP(-XW))/XW
    END IF
    A1 = (30.0_dp-15.0_dp*EXW-(15.0_dp-6.0_dp*EXW)*VCL**2) / &
         (30.0_dp-10.0_dp*EXW-(20.0_dp-5.0_dp*EXW)*VCL**2)
    Q1 = 0.25_dp*VCL**2/(1.0_dp-0.667_dp*VCL**2)
    DLT = SQB/PCL*(A1/X - SQRT(X)*(1.5_dp-0.5_dp*EXW+Q1) + &
          (1.0_dp-EXW+0.75_dp*VCL**2)/(1.0_dp+VCL**2)*(X-SQRT(X))/REAL(NL,KIND=dp))
    Y1 = 1.0_dp/(1.0_dp+DLT)
    CL0 = LOG(1.0_dp+1.0_dp/XS1)
    P2 = CL0*(0.07_dp+0.2_dp*EXW)
    Y2 = 1.5_dp*CL0*(X**3-X/3.0_dp)/(REAL(NL,KIND=dp)+0.75_dp/(1.0_dp+2.0_dp*B)**2*X**2) + P2*X
    PY = 1.0_dp + 0.06_dp*CL0**2/REAL(NL,KIND=dp)**2
    DT = SQRT(PY*Y1**2+Y2**2)
    CLLONG = CLEFF/DT

    DB = 1.0_dp/(1.0_dp+0.5_dp/B)
    CL1 = XS1*CL0
    P1 = 0.8_dp*(1.0_dp+CL1) + 0.2_dp*CL0
    P2 = 1.42_dp - 0.1_dp*DB + SQRT(CL1)/3.0_dp
    P3 = (0.68_dp-0.13_dp*DB)*CL1**0.165_dp
    P4 = (0.52_dp-0.1_dp*DB)*SQRT(SQRT(CL1))
    DLT = SQB/PCL*(P1/X**2*SQB/PCL + P3*LOG(REAL(NL,KIND=dp))/X - &
          (P2+P4*LOG(REAL(NL,KIND=dp)))*SQRT(X))
    CLTRAN = CLEFF*(1.0_dp+DLT)
  END IF
END SUBROUTINE COUL99I

!> Analytic Coulomb logarithm fit, ported from potekhinc.f's
!> `function COULAN2` verbatim (including its internal KEY-based
!> asymptote selection).
FUNCTION COULAN2(XS, XW0, V, XW1) RESULT(RESULT_VAL)
  REAL(KIND=dp), INTENT(IN) :: XS, XW0, V, XW1
  REAL(KIND=dp) :: RESULT_VAL
  REAL(KIND=dp), PARAMETER :: EPS = 1.0E-2_dp, EPS1 = 1.0E-3_dp, EULER = 0.5772156649_dp
  REAL(KIND=dp) :: XW, B, EA, E1, E2, CL0, CL1, CL2, EL
  INTEGER(KIND=i4) :: I, KEY

  IF (XS<0.0_dp .OR. XW0<0.0_dp .OR. V<0.0_dp .OR. XW1<0.0_dp) THEN
    WRITE(*,'(A)') 'COULAN2: invalid input'; STOP 1
  END IF

  KEY = 0
  DO I = 0, 1
    IF (I == 0) THEN
      XW = XW0+XW1
      B = XS*XW
    ELSE
      XW = XW1
      B = XS*XW
    END IF
    IF (I == 0 .OR. KEY == 2) THEN
      IF (XW < EPS) THEN
        KEY = 1
      ELSE IF (XW > 1.0_dp/EPS .AND. B > 1.0_dp/EPS) THEN
        KEY = 2
      ELSE IF (XS < EPS1 .AND. B < EPS1/(1.0_dp+XW)) THEN
        KEY = 3
      ELSE
        KEY = 4
      END IF
    END IF

    EA = EXP(-XW)
    E1 = 1.0_dp-EA
    IF (KEY /= 1) E2 = (XW-E1)/XW

    IF (KEY == 1) THEN
      CL0 = LOG((XS+1.0_dp)/XS)
      CL1 = 0.5_dp*XW*(2.0_dp-1.0_dp/(XS+1.0_dp)-2.0_dp*XS*CL0)
      CL2 = 0.5_dp*XW*(1.5_dp-3.0_dp*XS-1.0_dp/(XS+1.0_dp)+3.0_dp*XS**2*CL0)
    ELSE IF (KEY == 2) THEN
      CL0 = LOG(1.0_dp+1.0_dp/XS)
      CL1 = (CL0-1.0_dp/(1.0_dp+XS))/2.0_dp
      CL2 = (2.0_dp*XS+1.0_dp)/(2.0_dp*XS+2.0_dp) - XS*CL0
    ELSE IF (KEY == 3) THEN
      CL1 = 0.5_dp*(EA*EXPINT(XW,0)+LOG(XW)+EULER)
      CL2 = 0.5_dp*E2
    ELSE IF (KEY == 4) THEN
      CL0 = LOG((XS+1.0_dp)/XS)
      EL = EXPINT(B,0) - EXPINT(B+XW,0)*EA
      CL1 = 0.5_dp*(CL0+XS/(XS+1.0_dp)*E1-(1.0_dp+B)*EL)
      CL2 = 0.5_dp*(E2-XS*XS/(1.0_dp+XS)*E1-2.0_dp*XS*CL0+XS*(2.0_dp+B)*EL)
    ELSE
      WRITE(*,'(A)') 'COULAN2: invalid KEY'; STOP 1
    END IF

    IF (I == 0) THEN
      RESULT_VAL = CL1 - V**2*CL2
      IF (XW1 < EPS1) RETURN
    ELSE
      RESULT_VAL = RESULT_VAL - (CL1 - V**2*CL2)
    END IF
  END DO
END FUNCTION COULAN2

!> exp(XI)*E_{L+1}(XI) (exponential integral), ported from potekhinc.f's
!> `function EXPINT` verbatim.
FUNCTION EXPINT(XI, L) RESULT(RESULT_VAL)
  REAL(KIND=dp),    INTENT(IN) :: XI
  INTEGER(KIND=i4), INTENT(IN) :: L
  REAL(KIND=dp) :: RESULT_VAL
  REAL(KIND=dp),    PARAMETER :: GAMMA_E = 0.5772156649_dp
  INTEGER(KIND=i4), PARAMETER :: NREP = 21
  REAL(KIND=dp) :: CL, CI, CV, Q0, PSI, CMX, CM, DQ
  INTEGER(KIND=i4) :: I, K, M

  IF (XI >= 1.0_dp) THEN
    CL = REAL(L,KIND=dp)
    CI = REAL(NREP,KIND=dp)
    CV = 0.0_dp
    DO I = NREP, 1, -1
      CV = CI/(XI+CV)
      CV = (CL+CI)/(1.0_dp+CV)
      CI = CI - 1.0_dp
    END DO
    Q0 = 1.0_dp/(XI+CV)
  ELSE
    PSI = -GAMMA_E
    DO K = 1, L
      PSI = PSI + 1.0_dp/REAL(K,KIND=dp)
    END DO
    Q0 = 0.0_dp
    CMX = 1.0_dp
    CL = REAL(L,KIND=dp)
    CM = -1.0_dp
    DO M = 0, NREP
      CM = CM + 1.0_dp
      IF (M /= 0) CMX = -CMX*XI/CM
      IF (M /= L) THEN
        DQ = CMX/(CM-CL)
      ELSE
        DQ = CMX*(LOG(XI+1.0E-20_dp)-PSI)
      END IF
      Q0 = Q0 - DQ
    END DO
    Q0 = EXP(XI)*Q0
  END IF
  RESULT_VAL = Q0
END FUNCTION EXPINT

!> Electron-electron relaxation rate (for the thermal-conductivity
!> correction), ported from potekhinc.f's `function RELAXEE` verbatim.
FUNCTION RELAXEE(X, DENS, TEMP) RESULT(RESULT_VAL)
  REAL(KIND=dp), INTENT(IN) :: X, DENS, TEMP
  REAL(KIND=dp) :: RESULT_VAL
  REAL(KIND=dp) :: E, Y, V2, C, FJ, FREQ

  E = SQRT(1.0_dp+X**2)
  Y = 0.5245_dp*SQRT(DENS)/E/TEMP
  V2 = (X/E)**2
  C = 2.0_dp*V2 + 2.81042_dp*(1.0_dp-V2)
  FJ = (1.0_dp+1.2_dp/X**2+0.4_dp/X**4) * &
       ((Y**3/3.0_dp)*LOG((C+Y)/Y)/(1.0_dp+0.071514_dp*Y)**3 + &
       51.0033_dp*Y**4/(17.053_dp+Y)**4)
  FREQ = 0.0230117_dp*SQRT(X**3/E**5)*TEMP**2*FJ
  RESULT_VAL = 1.0_dp/FREQ
END FUNCTION RELAXEE

END MODULE CRUST_CONDUCTIVITY
