# MHD-VSH Roadmap

Master to-do list. Update this alongside the code — it's the durable,
shared record of what's done, what's blocked on analytical work, and
what's scaffoldable next, so it doesn't only live in a chat session.

## Status

Working end-to-end: grid → radial operators → boundary conditions →
`linear_solve.f90` → `regime_interface.f90`/`timestepper.f90` →
`field_diagnostics.f90`, with two real regimes running and tested —
`diffusion_regime.f90` (pure linear Ohmic diffusion) and
`hall_regime.f90` (combined resistive+Hall, the "primary interest"
case) — each producing a full time-resolved energy-*balance* plot on
every build (`E_poloidal`, `E_toroidal`, Joule dissipation, Poynting
flux — plus Hall Poynting for the combined regime — and their balance
residual). A ported TOV+crust-EOS solver (`tov_solver.f90`,
`crust_conductivity.f90`, `app/mhdvsh_tov.f90` -- ported from Dany
Page's NSCool, ASCL 1609.009, replacing an earlier port from a private,
uncredited source) supplies real `n_e(r)`/`eta(r)`/`f_H(r)` profiles,
and a dynamic Hall-CFL-limited timestep (`TIMESTEPPER::RUN_ADAPTIVE`,
`app/mhdvsh_hall_adaptive.f90`) is available alongside the fixed-`dt`
driver, which itself now supports checkpoint/restart
(`io_checkpoint.f90`) — see "Recently resolved" below for all three.
14/14 `ctest` targets passing. See `project.md` (FORD-generated API
docs, `ford project.md` → `docs/index.html`) for the module-by-module
reference.

## Recently resolved (2026-08-23, newest)

- **TOV/EOS solver re-ported from Dany Page's NSCool, replacing the
  earlier EOSNS-based port.** Per user: switch to a properly citable
  source (NSCool has an ASCL registration and associated papers;
  EOSNS was a private research code with no formal citation trail) --
  `src/core/tov_solver.f90`, `eos_table.f90`, `crust_conductivity.f90`
  now port NSCool's `TOV/TOV.f` (RK4-in-radius-r integration, a
  genuinely different numerical method than the prior RK45-in-log(P)
  port -- NSCool's own heuristic adaptive step size, not an embedded-
  error scheme) and `Code/conductivity_crust.f`'s `con_e_phon_ion_GYP`
  + `OYAFORM` (electron-phonon/electron-ion crust conductivity +
  density-only nuclear-structure fit). Investigated by downloading and
  reading NSCool's actual source directly (`NSCool.tar.gz`, public, no
  access gate), not from secondhand summaries -- see the design plan's
  own addendum for the full writeup.
  **Physics is unchanged, not upgraded**: `con_e_phon_ion_GYP`'s own
  header cites the *same two papers* (Potekhin 1999 A&A 346:345,
  Gnedin et al. 2001 MNRAS 324:725) the prior port already used --
  this replacement changes whose citable implementation supplies the
  physics, not the underlying physics itself (per user's own explicit
  choice of the GYP path over NSCool's different default Itoh+
  Yakovlev-Urpin conductivity formulation).
  **A real transcription bug was found and fixed during regression
  testing**: the first RK4 stage of the TOV integration mistakenly
  called the RHS routine instead of the original's explicit `k1=l1=
  m1=0` (which avoids a division-by-zero at r=0) -- caught because the
  post-fix `EM(1)` value matched NSCool's own bundled reference output
  to 9 significant figures, whereas the buggy version was off by ~50%
  at that same point. Final `M`/`R` now match that reference to
  `2e-7`/`3e-4` relative (the small residual `R` difference is traced,
  not assumed, to the integration's extreme low-pressure tail -- see
  `tov_solver.f90`'s own header). Also found (independently confirmed
  via `-fdefault-real-8`, same technique as before) the identical
  single-precision-literal quirk the prior port's conductivity code
  had, this time in `conductivity_crust.f`.
  **Two real physics capabilities were lost in this switch, disclosed
  rather than silently dropped**: no magnetic-field dependence at all
  (weaker even than the prior port's own dormant B-quantization
  branch), and no impurity-scattering treatment (the prior port's
  `COULIN`+`COUL99I` combination has no GYP-path equivalent).
  **`TOV_PROFILE_T`'s composition fields (`AH`/`ZH`/`XH`/`YE`/`YN`/
  `A_TABLE`) were dropped**, a deviation from the approved plan's own
  "same shape" wording, made once implementation revealed NSCool's
  `OYAFORM` derives Z/A directly from density with no need for a
  pre-tabulated per-row composition at all (and threading it through
  `TOV_PROFILE_T` would have created a circular module dependency,
  since `CRUST_CONDUCTIVITY` already depends on `TOV_SOLVER`). Crust-
  row filtering now uses a plain density threshold
  (`RHOCGS<=2.2e14 g/cm**3`, matching a threshold already present in
  both codebases) instead of the old `A_TABLE>0` gate.
  **`ETA_AND_F_HALL_AT` simplified**: `con_e_phon_ion_GYP` returns
  electrical conductivity `sigma` directly, so `eta=c**2/(4*pi*sigma)`
  no longer needs the prior port's `sigmae=3.26*tau*nel/meff`
  combination step.
  `app/mhdvsh_hall.f90`/`mhdvsh_hall_adaptive.f90`'s `R_MIN`/`R_MAX`
  updated to the new port's own crust extent (`10.3029`-`11.5633` km,
  down from `10.8033`-`11.6982` km -- a different, also-real crust
  extent from a different EOS table, not a placeholder correction).
  `data/eos/lowd-eos.ja.tab`/`.apr.tab` removed (unreadable by the new
  loader's format), replaced by `data/eos/APR_EOS_Cat.dat`. All three
  regression tests (`test_eos_table`, `test_tov_solver`,
  `test_crust_conductivity`) rewritten against the new interfaces and
  reference values; full suite still 14/14 passing.

## Recently resolved (2026-08-21, newest)

- **Checkpoint/restart implemented** (promoted to next-up per the
  previous entry below). New, regime-agnostic
  `src/core/io_checkpoint.f90` (`WRITE_CHECKPOINT`/`READ_CHECKPOINT`) --
  saves `Phi`/`Psi`'s coefficient arrays plus step/`t`/grid-shape
  metadata (`N_R`, `LMAX`, `R_MIN`, `R_MAX`) in this project's usual
  formatted-ASCII convention (one row per `(radial index, mode index)`
  pair, full double-precision round-trip precision), deliberately
  *not* saving any regime-specific cached state (`DIFFUSION_INIT`'s LU
  factorizations, `HALL_REGIME`'s saved grid/operators) since that's
  cheap and deterministic to rebuild via a fresh `*_INIT` call on
  resume -- matches the design sketch below exactly. Tested
  (`test_io_checkpoint.f90`): a distinct value per `(r,mode)` cell round-
  trips exactly.
- **Wired into `app/mhdvsh_hall.f90`** (the driver whose real ~2hr run
  actually crashed) as two new, purely-additive optional CLI args (3:
  checkpoint path, written every `CHECKPOINT_EVERY=50` steps plus once
  unconditionally at the end; 4: resume-from path) -- existing
  invocations with 0-2 args are completely unaffected. `N_STEPS`
  stays the TOTAL target step count; resuming runs only the remaining
  steps. The loaded checkpoint's grid shape is checked against this
  driver's own compiled-in `N_R`/`LMAX`/`R_MIN`/`R_MAX` before trusting
  it (a real, tested failure mode, not just planned -- see below), since
  a checkpoint carries no other physics parameters (`ETA`/`F_HALL`/`DT`/
  `N_SUB`) at all -- resuming means re-running this same driver, not a
  different configuration. The energy/per-degree-energy log files (args
  1-2) are appended to on resume rather than overwritten, with the
  already-logged resume-step row not repeated.
- **End-to-end smoke-tested, not just unit-tested**: a full 200-step run
  with checkpointing enabled; a synthetic step-150 checkpoint (crafted
  by hand from a completed run's own checkpoint, to exercise resume
  without needing an actual multi-hour crash to reproduce) resumed
  correctly to step 200 with the energy log correctly appended (no
  duplicate/missing rows, verified via `uniq -d`); the "checkpoint
  already at/past `N_STEPS`" guard and the grid-shape-mismatch guard
  both verified to fail cleanly with a clear message (not silently wrong
  output) via two more hand-crafted bad checkpoints.
- **Not done in this pass**: `mhdvsh_hall_adaptive.f90` (the dynamic-dt
  driver) has no checkpoint support yet -- only `mhdvsh_hall.f90`.
  `io_checkpoint.f90` itself is regime-agnostic and ready for it; wiring
  it in would additionally need the current adaptive `dt`/`tc` state
  preserved across a resume (not just `Phi`/`Psi`/step/`t`), not
  attempted here.

## Recently resolved (2026-08-21, most recent)

- **Dynamic, Hall-CFL-limited timestep: machinery built and wired in.**
  Per user formula: `tc = min(dr/(F_HALL*|current|))` over the radial
  domain, `dt < tc`. New `FIELD_DIAGNOSTICS::HALL_COURANT_TIMESTEP`
  computes this from the current `(Phi,Psi)` state; `|current|` is the
  RMS `|curl(B)|` over the sphere at each radius (an exact Parseval sum
  over VSH-mode coefficients -- `sum_lm j_lm.conj(j_lm)`, divided by
  `4*pi` and square-rooted -- not a physical-space reconstruction, which
  this codebase deliberately never does; the raw mode-sum alone would
  have the wrong units for a pointwise current-density magnitude, hence
  the `1/4pi`+`sqrt`. Confirmed with the user directly, 2026-08-21, after
  they asked for the derivation.). Extracted a shared
  `CURRENT_DENSITY_SQUARED_BY_R` helper out of `JOULE_DISSIPATION_RATE`'s
  own inline loop (both now call it) rather than duplicating that
  physics. New `TIMESTEPPER::RUN_ADAPTIVE` (+`REGIME_DT_I`/
  `REGIME_SET_DT_I` callback interfaces) recomputes `dt=CFL_SAFETY*tc`
  every `DT_RECOMPUTE_EVERY` steps and re-initializes whatever cached
  DT-dependent state the regime carries (amortizing the cost of
  `DIFFUSION_INIT`'s O(N_r**3)-per-l refactorization). `HALL_REGIME`
  gained `HALL_COMPUTE_DT`/`HALL_SET_DT` implementing those interfaces
  (required caching `ETA`/`ETA_PROFILE`, not previously retained). Per
  user's own choices: `CFL_SAFETY=0.5`, recompute cadence = `N_SUB`
  outer steps. New `app/mhdvsh_hall_adaptive.f90` driver exercises this
  (kept separate from `mhdvsh_hall.f90`, whose fixed-`dt` behavior stays
  untouched/reproducible for the prior dt-sweep/1000yr archive above) --
  a smoke test showed `dt` dropping from an initial `0.233` yr to `~0.021`
  yr as the seed mode's Hall-driven cascade populates new current-
  carrying modes, then stabilizing. Tested: `test_field_diagnostics.f90`
  (hand-computed `CURRENT_DENSITY_SQUARED_BY_R`/`HALL_COURANT_TIMESTEP`
  checks, plus a scalar-vs-profile-argument bit-identity check) and
  `test_timestepper.f90` (an independently-replicated recurrence check
  for `RUN_ADAPTIVE`'s recompute cadence/DT bookkeeping/`SET_DT` call
  count, via a toy regime). Directly advances analytical task 3 below,
  though via a direct Hall-CFL condition rather than literally computing
  `Rmag` as its own diagnostic (see that task's updated note).

## Recently resolved (2026-08-21, later again)

- **TOV + crust-EOS solver ported into MHD-VSH's own F90 tree** (per
  user: an existing, working F77 solver at `~/Desktop/EOSNS/` --
  `nstot.f`/`derivs.f`/`geteost.f`/`locate.f`/`odeint.f`/`potekhinc.f`
  -- ported rather than rewritten from scratch). New modules:
  `src/core/ode_integrator.f90` (generic Cash-Karp RK45 + bisection
  table lookup), `src/core/eos_table.f90` (EOS table load/interpolate),
  `src/core/tov_solver.f90` (`SOLVE_TOV_STAR` -> `TOV_PROFILE_T`
  radial profile), `src/core/crust_conductivity.f90` (Potekhin 1999
  electron-transport fitting formulas, `POTEKHINC`, plus a new
  `ETA_AND_F_HALL_AT` combining its `TAU` output with the TOV profile's
  own `n_e(r)` via `nstot.f`'s exact `sigmae=3.26*tau*nel/meff`
  formula, giving self-consistent `eta(r)`/`f_H(r)` in this project's
  code units). New driver `app/mhdvsh_tov.f90`; EOS tables
  (`lowd-eos.ja.tab`, `.apr.tab`) copied into `data/eos/`; new
  `scripts/plot_eos_profile.py` (n_e, rho, eta, f_H vs r, one figure,
  four independently-scaled log y-axes).
  **Directly resolves programming task 8** below (EOS files) and the
  `n_e(r)` piece of analytical task 3.
- **Two real bugs found and fixed during the port** (both self-caught
  during regression verification, not user-reported): (1) `nstot.f`'s
  own `rhocgs` variable, at the point it feeds the `UL`/`RHOC` length-
  scale formulas, is the central density in units of `1e14 g/cm**3`
  (its own comment: "give the central density in 10**14 g/cm3"), NOT
  raw cgs -- the first port attempt used raw cgs there directly, making
  `EOS_AT_DENSITY` immediately fail "rho is out of the table"; fixed by
  introducing an internal `RHOCGS_SCALE=RHOCGS/1e14` used only for those
  two formulas, everything else (including this module's own public
  `RHOCGS` argument, still true cgs) unaffected. (2) `TOV_PROFILE_T`'s
  `AH` field can't distinguish a homogeneous-matter row from a genuine-
  nucleus row (`COMPOSITION` always sets `AH>0`, substituting a dummy
  `1.0` for homogeneous rows) -- but `ETA_AND_F_HALL_AT`'s gate (mirroring
  `nstot.f`'s own `IF (a.gt.0.d0)` conductivity gate) needs the RAW
  table value, so a new `A_TABLE` field was added to carry it.
- **Independently regression-verified against the real, already-solved
  2015 reference output**, not just self-consistency: rebuilt the
  original, unmodified `potekhinc.f`/`nstot.f` etc. from source (a
  separate scratch harness, not touching `~/Desktop/EOSNS/` itself) and
  confirmed `mhdvsh_tov`'s own output -- `Radius=11.6982` km,
  `Mass=1.40088 Msun`, crust extent `10.8033-11.6982` km, and a spot-
  checked `n_e`/`rho` row -- matches `~/Desktop/EOSNS/fort.34`/`PL.DAT`
  (an `M=1.40` star, `rhocgs=9.88e14`) closely. Also found (and
  confirmed harmless, not a bug to reproduce) that `potekhinc.f`'s own
  `COULIN`/`COUL99I` have a commented-out-guard/active-body dead-code
  structure making their "magnetic fit" branches genuinely unreachable
  in the validated reference behavior -- ported faithfully as an
  explicit early `RETURN`, not silently dropped. And a ~1e-6-level
  discrepancy between this port and a live rebuild of the original
  traced to and confirmed as the original's un-suffixed (non-`d0`)
  literal constants being parsed as single precision before widening to
  double (a real, if tiny, F77 gotcha in the original, not a
  transliteration error) -- confirmed via `-fdefault-real-8`, not
  assumed; `test_crust_conductivity.f90`'s tolerance documents this.
  Independent test coverage: `test_eos_table.f90`, `test_tov_solver.f90`
  (regression against the real reference numbers above), and
  `test_crust_conductivity.f90` (a standalone F77 harness built from the
  *unmodified* original `potekhinc.f`, not this port's own code).
- **Optional per-radial-row `ETA_PROFILE`/`F_HALL_PROFILE` overrides**
  added to `DIFFUSION_INIT`, `HALL_INIT`, and `HALL_INDUCTION_RHS` (per
  user: "we want the option to override these values with constants
  though, for testing and validation purposes") -- absent, behavior is
  bit-identical to today's uniform-scalar case (confirmed: full 13-test
  suite unchanged after wiring). Not yet actually threaded end-to-end
  from `ETA_AND_F_HALL_AT`'s output into a driver's `HALL_INIT` call --
  the plumbing exists, a driver using it doesn't yet.
- **`mhdvsh_hall.f90`'s `R_MIN`/`R_MAX` updated to the EOS-derived crust
  extent** (`10.8033325018`-`11.6982211606` km for the `M=1.40`
  reference star), superseding the `9.0`/`10.0` km placeholder from the
  "Recently resolved (2026-08-21)" entry below -- rerunning that
  driver's own dt-sweep/1000yr archive at the corrected radii is
  deliberately deferred, not done as part of this pass.

## Recently resolved (2026-08-21, even later)

- **`dt`-stability sweep + 1000yr late-onset-instability check: both
  clean, archived permanently.** Per user: swept `dt in {0.001, 0.002,
  0.005, 0.01, 0.02}` yr over `tmax=10` yr at the current combined-regime
  parameterization (`Lmax=30`, NS-crust `R=9-10` km, `eta=1e-6`,
  `F_Hall=0.01`) -- all five stable, no instability boundary found in
  that range. Then, specifically to test the user's own hypothesis that
  a strong Hall drift building up over "a few hundred years" might
  destabilize a `dt` that looks fine over a short window, ran the
  largest swept `dt=0.02` all the way to the original `tmax=1000` yr
  target (`50000` steps, ~2hr wall-clock once system CPU contention was
  resolved -- see below). Result: **91.81% of magnetic energy remaining
  after 1000 years, zero NaN/negative-energy/non-monotonic points across
  all 1001 logged rows, largest step-to-step change at t=0 not anywhere
  later** -- the hypothesized late-onset instability did not
  materialize. Full writeup, data, and plots permanently archived at
  `results/hall_dt_stability_2026-08-21/` (not `build/artifacts/`, which
  is regenerated/ephemeral) -- see that directory's own `README.md`.
- **Real system-stability lesson from this session, not just a code
  one:** the long dt=0.02/tmax=1000yr run's OpenMP-parallelized radial
  loops pegged 12+ CPU cores continuously, which appears to have
  starved the rest of the machine badly enough to crash/reboot it
  mid-run (confirmed via `uptime` showing a fresh boot partway through
  -- not just a slow app, an actual system-level crash). The run had to
  be relaunched from scratch (no checkpoint/restart existed) and was
  deliberately reniced to the lowest CPU/IO priority for the rest of its
  execution to avoid repeating this. Directly motivates programming task
  3 below (`io_*.f90` checkpoint/restart) being promoted to next-up, not
  just a nice-to-have -- also worth remembering for any future
  long-running `Lmax=30`+ multi-hour job on this machine: renice
  proactively, don't wait for a crash.

## Recently resolved (2026-08-21, later)

- **`mhdvsh_hall`'s `ETA` corrected to a physically realistic crustal
  value**: `1e-6` km**2/yr (per user: realistic range is `1e-8` to
  `1e-5`, picked the middle), replacing the inherited `0.05` km**2/yr
  toy value the diffusion-only demo originally used -- that value was
  ~50000x too large for crust physics and caused the near-total
  (`98.8%`/`88.5%`) decay seen in the two NS-crust-radius runs above.
  With the corrected `ETA`: energy `3.936e38 -> 3.805e38` erg over
  `T=2` yr, **`96.7%` remaining** (only `3.3%` lost) -- diffusion is now
  a slow secular process relative to the `T=2` yr window, and the
  poloidal/toroidal per-degree energy plots now show the seed mode
  (`l=1`) staying essentially flat, much closer to the pure-Hall-only
  behavior, with diffusion just gently reshaping the cascade rather than
  dominating it. Also added per-degree energy logging to `mhdvsh_hall.f90`
  itself (an optional 2nd CLI arg, reusing `POLOIDAL_MAGNETIC_ENERGY_BY_L`/
  `TOROIDAL_MAGNETIC_ENERGY_BY_L` -- previously only the pure-Hall harness
  had this) and fixed two real bugs in `scripts/plot_hall_energy_by_l.py`
  found while reusing it for this driver: a hardcoded "Hall stability
  experiment" title (now inferred from the artifacts path, matching
  `plot_energy_budget.py`'s existing fix) and a fixed absolute y-axis
  range tuned for the old pure-Hall code-unit scale (`[1e-20,1e1]`) that
  silently cut off almost all the erg-scale data from this driver --
  replaced with a floor relative to each run's own peak value.

## Recently resolved (2026-08-21)

- **`mhdvsh_hall`'s domain moved into the actual neutron star crust**
  (per user): `R_MIN=9.0`/`R_MAX=10.0` km (core-crust boundary to
  surface), replacing the earlier `0.5`/`1.0` km toy shell -- all other
  parameters unchanged (`N_R=40`, `LMAX=30`, `ETA=0.05`, `F_HALL=0.01`,
  `DT=0.01`, `N_SUB=10`, `N_STEPS=200`, same single-mode `Phi(1,0)` seed).
  Result: energy `3.936e38 -> 4.545e37` erg over `T=2` yr, `11.5%`
  remaining -- much less dissipation than the old `0.5-1.0` km shell's
  `1.18%` remaining, consistent with the expected physics (diffusion
  decay rate `~eta*[(pi/dr)**2 + l(l+1)/r**2]` shrinks as both the shell
  width `dr` and the radius `r` grow). No code changes beyond the two
  parameter values; full suite still 10/10 passing.

## Recently resolved (2026-08-20, latest)

- **Combined resistive+Hall regime: first working version.** New
  `src/regimes/hall/hall_regime.f90` (`HALL_INIT`/`HALL_ADVANCE`),
  reusing `DIFFUSION_STATE_T` directly per the design plan
  (`.claude/plans/i-would-like-to-quizzical-comet.md`, Stage 4). Per
  outer step: `N_SUB` explicit Hall RK4 substeps (via
  `HALL_INDUCTION_RHS`, homogeneous-Dirichlet placeholder BC between
  them) then one implicit `DIFFUSION_ADVANCE` call.
  **Key design finding**: `DIFFUSION_REGIME::SOLVE_FIELD` already zeros
  its RHS's boundary rows before every solve, and its cached matrix
  already encodes the true vacuum BC -- so ending every outer step with
  a diffuse call enforces the real physical BC exactly, for free, with
  **zero new BC-on-raw-vector code needed** (the plan originally assumed
  this would be a hard prerequisite; verified false by reading
  `diffusion_regime.f90` directly). Confirmed empirically, not just
  reasoned: `test/test_hall_regime.f90`'s outer-Robin-BC check on `Phi`
  passes at `2.5e-15` even though the preceding Hall substeps only used
  the cheap placeholder BC. New `app/mhdvsh_hall.f90` driver logs the
  full combined energy budget (adds `hall_poynting_flux_rate` as a new
  column) -- `scripts/plot_energy_budget.py` needed zero changes, it's
  already generic against the column set.
- **`HALL_POYNTING_FLUX_RATE` optimized to match `HALL_INDUCTION_RHS`.**
  Discovered while building the above: at `LMAX=30`, calling this
  formerly-unoptimized `O(Nlm**3)` diagnostic every logged step would
  have been infeasible (`mhdvsh_hall.f90` needs it ~200 times). Applied
  the same active-mode-list + selection-rule restriction already
  validated on `HALL_INDUCTION_RHS` (exact, not approximate -- boundary-
  sampled quantities are literally the only thing this boundary-only
  flux formula reads). Confirmed lossless: `test_field_diagnostics.f90`'s
  `hall_poynting_flux_rate` check produces the bit-identical error value
  before and after.
- **`mhdvsh_hall`'s energy-budget artifact deliberately kept out of
  `ALL`** (deviates from the design plan's original "wire into ALL like
  diffusion" -- a cost tradeoff discovered only while implementing, not
  planned for): at `LMAX=30`/`N_SUB=10`, a full 200-step run costs ~5
  minutes (the real price of resolving the Hall cascade at that
  resolution), vs. diffusion's sub-second implicit solve -- forcing that
  into every `cmake --build .` would be a surprising, unwelcome
  development-loop cost. The `mhdvsh_hall` executable itself still
  builds by default; only the artifact-generating run is opt-in
  (`cmake --build . --target hall_energy_budget_plot`).

## Recently resolved (2026-08-20, still later)

- **Per-degree magnetic energy diagnostics.** New
  `FIELD_DIAGNOSTICS::POLOIDAL_MAGNETIC_ENERGY_BY_L`/
  `TOROIDAL_MAGNETIC_ENERGY_BY_L` -- the same `mhd-vsh-relations.pdf`
  `E_B,pol(t)`/`E_B,tor(t)` volume-integrated formulas
  `TOTAL_POLOIDAL_MAGNETIC_ENERGY`/`TOTAL_TOROIDAL_MAGNETIC_ENERGY`
  already implement, broken out per degree `l` instead of summed over
  it (confirmed algebraically identical to the existing per-mode terms,
  not a new derivation -- `R_l=r/sqrt(Lambda_l)` substitution collapses
  the PDF's `E_B,pol(t)=Sum_kl...`/`E_B,tor(t)=Sum_kl...` equations
  exactly onto what was already coded). Tested
  (`test_field_diagnostics.f90`): matches hand-derived per-mode values
  at machine precision, AND `SUM(*_BY_L)` reproduces the existing
  `TOTAL_*` functions exactly -- a genuine internal identity, not an
  approximation. `app/mhdvsh_hall_stability_experiment.f90`'s private,
  duplicate `MODE_ENERGY_BY_L` helper (a combined poloidal+toroidal
  version, written before this) was removed in favor of these; the
  driver now logs both separately to a second data file
  (`hall_energy_by_l.dat`).
- **First per-mode energy plots, and a genuine (if tiny) numerical
  finding.** `scripts/plot_hall_energy_by_l.py` (styled to match
  `plot_hall_mode_amplitudes.py` -- same ordinal blue ramp/colorbar,
  plain log not symlog since energy is non-negative). Both poloidal and
  toroidal energy-by-l show a clean two-band structure: the "real"
  cascade (odd-l for poloidal, even-l for toroidal, matching the
  Phi/Psi parity pattern already found) sits ~30+ orders of magnitude
  above a second band -- the "wrong-parity" modes (even-l poloidal,
  odd-l toroidal), which are exactly zero at the domain midpoint (per
  the earlier per-mode-amplitude check) but NOT exactly zero once
  integrated over the whole radial domain. Likely floating-point-level
  leakage through formally-zero coupling paths, accumulated over 2000
  RK4 steps -- not investigated further, but a concrete example of this
  new diagnostic catching something the coarser single-point amplitude
  check missed.

## Recently resolved (2026-08-20, later still)

- **First real Stage 2 Hall stability run: success, after finding and
  fixing a genuine bug in `HALL_INDUCTION_RHS`.** Per user direction: a
  single seed mode (`Phi` at degree 1, order 0, `Psi` identically zero),
  no simplified terms (all 4 term-groups), RK4, evolved to `LMAX=30` so
  self-coupling could cascade to high degree. First attempt blew up to
  NaN by step 20 -- root cause was a real singularity, not a tuning
  issue: both `Phi_dot`/`Psi_dot` carry an explicit `1/Lambda_n` factor,
  undefined at target degree `n=0` (`Lambda_0=0`), and the self-paired
  seed mode's triangle band includes `n=0` on the very first step.
  Fixed by excluding `n=0` from the target loop unconditionally
  (physically consistent -- a degree-0 poloidal/toroidal potential
  carries no field). After the fix: the full 2000-step run (`dt=0.001`,
  `T=2` code-yr) completed with no blowup, energy held within
  `[0.9999,1.0020]` of its initial value throughout (consistent with the
  Hall term's exact energy conservation -- a strong sanity check, not
  just "didn't crash"), and `max_active_l` climbed steadily from 3 to
  the full 30 by `t~1.16`, then saturated there (a truncation ceiling,
  not a physical one). New standalone driver:
  `app/mhdvsh_hall_stability_experiment.f90` (deliberately
  `EXCLUDE_FROM_ALL`, not wired to a fixed artifact plot, per the plan --
  parameters are expected to keep changing).
- **`HALL_INDUCTION_RHS` optimized: exact, not approximate.** Two
  changes, both confirmed to leave results bit-identical to the
  pre-optimization version (`test_hall_induction` output unchanged):
  (1) the outer `(k,l),(k',l')` pair now only ranges over modes actually
  nonzero in the current state (skips branches that provably contribute
  zero -- makes a sparse/single-mode-seeded run like the one above
  tractable at `LMAX=30`, where the full `(LMAX+1)**2=961`-mode space
  would otherwise force a ~961**3 raw loop); (2) target `m` is computed
  directly as `l+l'` (forced by `GWI`/`GWJ`'s own selection rule) rather
  than searched, and `n` is restricted to the triangle band
  `|k-k'|<=n<=k+k'` -- both exact consequences of the coupling
  coefficients' own vanishing conditions.
- **User correction, verified numerically: `GWI` is exactly symmetric
  and `GWJ` exactly antisymmetric under the `(k,l)<->(k',l')` swap.**
  Checked across 1392 nonzero-coupling sextuples (residuals `~1e-16`/
  `~1e-14`). `hall_induction.f90` updated to call each coefficient once
  and derive the swapped value algebraically, halving the
  coupling-coefficient cost in the inner loop.

## Recently resolved (2026-08-20, later same day)

- **Hall regime time-advance design plan written** (see
  `.claude/plans/i-would-like-to-quizzical-comet.md` for the full
  document): staged path from the `.tex`-only Hall induction equations
  to a combined resistive+Hall regime -- Stage 0 (settle explicit-
  integrator family), Stage 1 (`HALL_INDUCTION_RHS` + correctness test),
  Stage 2 (standalone stability harness), Stage 3 (hyperresistivity),
  Stage 4 (combined regime), Stage 5 (CMake/test wiring). Two scope
  decisions made: stability harness uses simplified homogeneous-Dirichlet
  BC (the real vacuum Robin BC has no existing explicit-data
  enforcement machinery -- deferred to Stage 4), and an ad-hoc
  multi-mode axisymmetric IC rather than waiting on the still-unstarted
  force-free Bessel-Riccati IC.
- **Stage 0 (toy integrator comparison): informative, but redirected --
  not a clean FE-vs-RK4 answer.** Three iterations of a hand-simplified
  2-4 mode WKB toy model each hit a different physical/numerical
  artifact: (1) fixing one mode as a constant "background" while only a
  second mode evolves creates an unphysical infinite energy reservoir
  (caught by the user); (2) even after fixing that (closed,
  self-consistent axisymmetric multi-mode system, no fixed background),
  the 2-of-4-term-group simplification used for tractability turned out
  to NOT be energy-conserving on its own -- all three integrators agreed
  with each other at small `dt` on a real (non-numerical) ~4x energy
  growth over a short time, meaning the toy problem itself, not the
  integrator, was untrustworthy. Conclusion: hand-simplifying the
  formula in a scratch script isn't a reliable way to settle this
  question -- redirected to Stage 1 (build the tested, full formula,
  then reuse it for the integrator comparison) instead of continuing to
  patch the toy model.
- **Stage 1 complete: `src/core/hall_induction.f90` (`HALL_INDUCTION_RHS`)
  implemented and tested.** All 4 term-groups of the Hall induction
  equations (`mhd-vsh-relations.tex:392-404`) transcribed, including the
  easy-to-miss detail that two of the four terms use the SWAPPED
  `(k,l)<->(k',l')` index order through the `I`/`J` Gaunt-type
  coefficients (not a typo in the source -- `GWI`/`GWJ` aren't symmetric
  under that swap). A pure, side-effect-free RHS evaluator (no BC
  applied, matches `RADIAL_OPERATORS`' convention), living in
  `mhdvsh_core` alongside `field_diagnostics.f90` rather than in a
  regime module, so it's reusable by both a future regime's `ADVANCE`
  and any standalone experiment harness. Tested (`test/test_hall_induction.f90`,
  `ctest` target `hall_induction`) against an independently-derived
  closed-form reference (two constant-radial-profile axisymmetric modes
  chosen so 2 of the 4 term-groups vanish exactly, leaving the other 2
  fully hand-computable) -- matches to `~1e-12`, with explicit
  nontrivial-value checks confirming it's not a trivial 0=0 pass. Known
  deferred optimization (noted in the module's own header): the current
  `O(Nlm**3 * N_R)` loop structure doesn't yet exploit the `l'=m-l`/
  triangle-band selection-rule collapse the plan calls for -- correctness
  came first, per the plan's own Stage ordering; Stage 2 will reveal
  whether this becomes a real bottleneck.

## Recently resolved (2026-08-20)

- **Physical units standardized.** New `src/core/units.f90` (`UNITS`
  module) names this project's adopted code-unit system: `B`: `10**12`
  Gauss, length: km, time: (Julian) year, with `eta` in `km**2/yr` and
  `f_H` in `km**2/(10**12 G)/yr` to match. No formula in
  `field_diagnostics.f90` or any regime's induction equation needed to
  change -- confirmed on inspection that every one of them is already
  genuinely unit-system-agnostic (linear/homogeneous in `B`, with
  `eta`/`f_H` as free inputs already carrying whatever length/time units
  the caller gives them), so this was a documentation-plus-conversion-
  factors task, not a physics-formula rewrite. `UNITS` provides
  `ENERGY_UNIT_ERG` (`= B_UNIT_GAUSS**2 * LENGTH_UNIT_CM**3`, by
  dimensional analysis alone -- `[energy]=[B]**2*[length]**3` in
  Gaussian cgs regardless of a formula's internal bookkeeping) and
  `POWER_UNIT_ERG_PER_S` (`= ENERGY_UNIT_ERG / TIME_UNIT_S`) as the
  multiplicative factors to convert a diagnostic's code-unit output to
  physical cgs. Per user decision, these conversions are applied only at
  the *reporting* layer (`app/mhdvsh_diffusion.f90`'s console output and
  `.dat` logging, `scripts/plot_energy_budget.py`'s axis labels) --
  `field_diagnostics.f90` itself stays untouched/unit-agnostic, matching
  the existing "physics-blind core module" pattern
  ([[project-mhdvsh-linear-solve]]-style). Tested (`test_units.f90`,
  checks the derived factors match their own documented dimensional
  relations) and exercised end-to-end: the diffusion driver's demo run
  now reports `energy(t=0) ~ 8.7e38 erg` for its toy `R_MIN=0.5`/
  `R_MAX=1.0` km problem, rather than an arbitrary O(1) number.
- **Still open / explicitly out of scope for this pass:** the demo
  driver's `R_MIN`/`R_MAX`/`ETA`/`N_STEPS` etc. are still toy numbers
  (a 0.5-1.0 km shell), just now interpreted in real units rather than
  arbitrary ones -- nobody has picked physically-motivated magnetar
  crust parameters (NS radius ~10-12 km, realistic crustal `eta`) yet;
  that's a separate, future decision, not implied by this change. Also
  untouched: the `.dat` file's `t`/`step` columns (already directly in
  years/steps, no conversion needed) and the plot's shared y-axis, which
  still mixes energy (erg) and rate (erg/s) series on one scale, a
  pre-existing quirk unrelated to this change.

## Recently resolved (2026-08-11)

- **`Phi`/`Psi` → `B_pol`/`B_tor` energy formula: confirmed, not just
  provisional.** The user's own derivation (`analytic_formulas/
  mhd-vsh-relations.{tex,pdf}`) matches `field_diagnostics.f90`'s
  existing energy formula structurally exactly (only the unit prefactor
  changed, to the proper Gaussian-cgs `1/8pi`) — cross-validated further
  by a direct numerical check of FORTVSH's VSH basis orthonormality
  (scratchpad, not a permanent test). Still open: the raw `B_pol`/`B_tor`
  coefficient split itself (`R_l`, `VSH_POL_DN`/`VSH_POL_UP`) was taken
  as given in that document rather than re-derived from first
  principles — strong indirect confirmation, not an independent proof.
- **Joule dissipation + Poynting flux: derived and implemented.** Same
  document's Eqs. 6, 8. `field_diagnostics.f90` now has
  `JOULE_DISSIPATION_RATE`/`POYNTING_FLUX_RATE`; the full balance
  identity (`Edot_B,pol + Edot_B,tor = Edot_J + Edot_S`) is plotted
  directly as a residual line, and converges to ~0 as expected (checked
  against a real run, not just unit-tested in isolation).
- **Diffusion regime's energy-budget artifact** (programming task,
  formerly "more artifacts as regimes come online") is now complete for
  this regime: all 5 series, auto-regenerated on every build via
  `build/artifacts/diffusion/energy_budget/{data,plots}/`.
- **Local git**: initialized, one commit. Still no GitHub remote (user
  said not yet) — see Infrastructure below.
- **Hall regime energy balance: fully resolved.** The user added a
  "Energy Budget in Hall Limit" section to `mhd-vsh-relations.tex`:
  confirmed the Hall Joule term is exactly zero (`j.(j x B)=0`
  identically, so `JOULE_DISSIPATION_RATE`'s existing formula already
  covers the weak-Hall regime unchanged), and derived a new Hall
  Poynting flux formula -- a genuinely harder, mode-coupled (Gaunt/
  FORTVSH's `GWI`, `O(Nlm**3)`) expression, the first coupling-based
  formula anywhere in this module. Implemented as
  `HALL_POYNTING_FLUX_RATE` and tested against an independently
  (differently-structured) re-derived version of the same formula, not
  just unit-tested against itself.
- **Hall induction equation itself (analytical task 1): derived.** Same
  `mhd-vsh-relations.tex`, new "Hall induction equations" derivation
  following the energy-budget section -- gives `Phi-dot`/`Psi-dot` in
  the weak Hall limit as explicit mode-coupled sums over `(k,l)`,
  `(k',l')` contracted through the same `I^{nm}_{klk'l'}`/`J^{nm}_{klk'l'}`
  Gaunt-type coefficients as the Hall Poynting formula (plus an
  appendix with the underlying `j x B` and `curl(f_H j x B)` VSH
  expansions the two induction-equation components were built from).
  This was the one piece analytical task 1 was waiting on -- the
  weak-Hall regime module (programming task 1 below) no longer needs a
  placeholder Hall term, it can be implemented directly from these
  formulas. Not yet cross-checked against a second, independent
  derivation the way the Hall Poynting formula was, and not yet
  transcribed into `field_diagnostics.f90`/a regime module -- still a
  `.tex`-only result at this point.

## Analytical tasks (yours)

1. **Dynamo / ion-MHD limit analytical work.** Explicitly deferred --
   don't start the dynamo regime or resolve the kinematic-vs-full-MHD
   question until this is done.
2. *(Possibly joint)* the discrete `div(B)=0` self-consistency check --
   a numerical-discretization question (does the FD scheme preserve an
   identity that's exact in the continuous theory), deferred rather
   than guessed at.
3. **Dynamic `dt` from the magnetic Reynolds number `Rmag`** (per user,
   2026-08-20) -- **substantially superseded, 2026-08-21**: rather than
   computing `Rmag = c*B/(4*pi*n_e*eta)` as its own diagnostic and
   deriving `dt` from it, the user instead specified a direct Hall-CFL
   condition (`tc = min(dr/(F_HALL*|current|))`, `dt < tc`), now fully
   implemented and wired in (see "Recently resolved (2026-08-21, most
   recent)" above -- `FIELD_DIAGNOSTICS::HALL_COURANT_TIMESTEP`,
   `TIMESTEPPER::RUN_ADAPTIVE`, `HALL_REGIME::HALL_COMPUTE_DT`/
   `HALL_SET_DT`, `app/mhdvsh_hall_adaptive.f90`). This achieves the
   operational goal (a physically-motivated, self-adjusting `dt`) without
   `Rmag` itself ever being computed as a named quantity anywhere --
   `Rmag` as a literal logged/plotted diagnostic is still open if wanted
   separately, but is no longer a blocker for anything. `n_e(r)` (the
   piece this task originally flagged as the one open gap) is also now
   available, via programming task 8's TOV/EOS port
   (`CRUST_CONDUCTIVITY::ETA_AND_F_HALL_AT`) -- not yet threaded into
   `mhdvsh_hall_adaptive.f90` itself, which still uses a uniform scalar
   `F_HALL`, not `ETA_AND_F_HALL_AT`'s per-radius profile.

## Programming tasks (scaffoldable now, independent of the above)

1. ~~**Weak-Hall electron-MHD regime module**~~ -- **DONE** (see "Recently
   resolved" above, `src/regimes/hall/hall_regime.f90`). First working
   version: fixed `N_SUB`/`DT` (no hyperresistivity yet -- still
   deferred, see the design plan's "Deferred, not skipped" section).
   Dynamic Hall-CFL `dt` now exists (2026-08-21, see analytical task 3)
   but only in the separate `app/mhdvsh_hall_adaptive.f90` driver --
   `mhdvsh_hall.f90` itself still uses fixed `DT`, deliberately (per
   user, to keep its prior dt-sweep/1000yr archive reproducible).
2. **`(Phi,Psi) -> physical B` converter**, using `VSH_POL_UP_ALL`/
   `VSH_POL_DN_ALL` (formulas already identified, just not implemented)
   -- needed eventually for I/O/visualization and for evaluating
   nonlinear terms pseudospectrally.
3. ~~**`io_checkpoint.f90`: checkpoint/restart**~~ -- **DONE** (2026-08-21,
   see "Recently resolved" above): `src/core/io_checkpoint.f90` +
   wiring into `app/mhdvsh_hall.f90`. Field-snapshot output (a general
   I/O layer for e.g. visualization, distinct from a resume checkpoint)
   is still open -- not attempted. Also still open: wiring checkpoint
   support into `mhdvsh_hall_adaptive.f90` (see its own note above).
4. **Force-free Bessel-Riccati initial conditions** -- unlike the
   analytical tasks above, this has a solid reference already (Igoshev,
   Elfritz & Popov 2016, MNRAS 462, 3689, Appendix A; arXiv:1608.08806),
   so it's an implementation task, not an open derivation: spherical
   Bessel functions `j_l`/`n_l` (FORTVSH provides none) + a root-finder
   for the transcendental eigenvalue equation.
5. **Formal regression test for the energy-balance identity itself** --
   currently only checked informally (comparing a real run's logged
   data by hand/script), not as a `CTest`. Would need a documented,
   dt-dependent tolerance (backward-Euler is first-order, so the
   residual shrinks as `dt` shrinks, not exactly zero at any finite dt).
6. **More artifacts for other diagnostics/regimes** -- e.g. a radial
   energy-density profile plot, boundary-residual plots, and the
   energy-budget artifact for every regime after diffusion -- following
   the `build/artifacts/<regime>/<diagnostic>/{data,plots}/` convention
   in `CMakeLists.txt`'s "Artifacts" section.
7. **`induction`/`magnetofriction`/`mhd_full` regime scaffolding** as
   each is actually started.
8. ~~**Equation-of-state (EOS) files for different neutron star
   configurations**~~ -- **DONE** (2026-08-21, see "Recently resolved"
   above): a full TOV+crust-EOS solver was ported from
   `~/Desktop/EOSNS/` into `src/core/{ode_integrator,eos_table,
   tov_solver,crust_conductivity}.f90` + `app/mhdvsh_tov.f90`, rather
   than just adding static data files -- gives `n_e(r)`, `eta(r)`,
   `f_H(r)` for any central density the EOS table covers, not just one
   fixed configuration. Only the `lowd-eos.ja.tab`/`.apr.tab` low-density
   (crust) tables from that source are wired in so far; other NS
   configurations (different EOS models entirely, not just different
   central densities against the same crust table) would need their own
   table file in the same 6-column format -- not attempted.

## Infrastructure

- **Documentation**: `project.md` (FORD config) is set up and verified
  (`ford project.md` runs clean, zero warnings). Doc-comment style
  throughout matches FORTVSH's own `!>` convention.
- **Git**: initialized locally (one commit so far). No GitHub remote yet
  -- user explicitly deferred that decision, don't create one or push
  without asking again.
- **CI / GitHub**: `.github/workflows/ci.yml` (build FORTVSH from
  source, build MHD-VSH against it, run `ctest`, upload
  `build/artifacts/**`) and `.github/workflows/docs.yml` (FORD build +
  GitHub Pages deploy, mirroring FORTVSH's own working setup almost
  exactly) are written and ready, but **inert until there's a GitHub
  remote** to run them against.
