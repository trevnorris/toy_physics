---
unit_id: 146
batch: IV.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 146 (post-SymPy-eps-fallback rework)

## Per-finding outcomes

Note on directive bookkeeping: `redteam/directives/stage_146.md` still does
not contain `## Applied: Fn` blocks. Verification proceeds by inspecting the
current script state and exec logs directly, plus the orchestrator notes
supplied with this task (SymPy F1/F2 eps-sample fallbacks installed,
Mathematica F1 numeric-sample fallback retained, Cluster A F3 fix retained).

### F1 — tautological_check

**Classification:** resolved

**What changed:**

- SymPy `scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:96-121`
  builds the physical integral form via
  `Sigma_eps = (1-eps)*Sigma.subs(Pi, Pi_star) + eps*varsigma_test`,
  `gbar_phys = sp.integrate(Sigma_eps*sp.cos(sp.pi*x/2), (x, 0, 1))`, and
  `Sbar_phys = sp.integrate(Sigma_eps*Kq, (x, 0, 1))` (with
  `varsigma_test = 6*x*(1-x)`). The residual expressions
  `g_eps_residual_expr` and `S_eps_residual_expr` are then evaluated at
  `eps in {1/10, 1/2}` via `sp.N(..., 30)` with tolerance `1e-15` and a
  raised `AssertionError` on failure. The old algebraic restatement is
  gone.
- Mathematica `mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:88-121`
  mirrors the integral construction (`gBarPhys`, `sBarPhys`,
  `varsigmaTest`) and chops the numeric residuals at `eps -> 1/10` and
  `eps -> 1/2` against `10^-6` (justified in-script by low-WP complex
  near-zero noise from the numeric `pStar` substitution).

**Assessment:**

The new checks are non-tautological: `gbar_phys` and `Sbar_phys` are built
by integrating `(1 - eps) Sigma_*(x) + eps varsigma(x)` against the
physical kernels `cos(pi x/2)` and `Kq(x)`, then compared to the affine
right-hand side built from `gminus` (resp. `Sformula.subs(Pi, Pi_star)`)
plus `eps*(gbar_v - gminus)`. The identity therefore exercises the
linearity of the integral plus `gminus = int Sigma_*(x) cos(pi x/2) dx` and
`Sformula(Pi_*) = int Sigma_*(x) Kq(x) dx`, which is the substantive
content the auditor wanted on the test surface. A bug in `Sigma`, `Kq`,
`Pi_*`, or `gminus` would now perturb the residual away from zero.

- SymPy transcript lines 156-161: residuals are 3.57E-18, 5.58E-31,
  1.98E-18, 6.32E-31 for the four eps/quantity pairs — all well under the
  1e-15 tolerance; both PASS lines print.
- Mathematica transcript lines 33-38: both eps samples Chop to literal `0`
  under `10^-6`; both PASS lines print.

The orchestrator-installed eps-sample fallback is documented in the SymPy
script's in-line comments at lines 106-108.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**

- SymPy `scripts/...sympy_audit.py:30-53` inserts
  `gPi_direct = sp.integrate(Sigma*sp.cos(sp.pi*x/2), (x, 0, 1))` and
  `Sq_direct = sp.integrate(Sigma*Kq, (x, 0, 1))`, each gated on
  `.has(sp.Integral)` to fall back to a 4-point numeric check at
  `Pi in {7/10, 11/10, 17/10, 23/10}` with tolerance `1e-25`.
- Mathematica `mathematica/...mathematica_audit.wl:50-51` adds
  `sDirect = FullSimplify[Integrate[sigma*kq, {x, 0, 1}], Assumptions -> p > 0]`
  followed by `expectZero["S_q(Pi) direct-formula", sDirect - sFormula]`.

**Assessment:**

Both engines now anchor the load-bearing closed-form moments against the
direct integrals of `sigma*cos(pi x/2)` and `sigma*Kq`. The Mathematica
side closes both as full symbolic identities. The SymPy `g(Pi)` half
closes symbolically; the `S_q(Pi)` half engages the documented numeric
fallback (the kernel against an exponential leaves an unevaluated
`Integral`), with `diff=0E-165` at all four samples — non-tautological
because `Sq_direct` and `Sformula` are constructed independently from the
shared `Sigma` and `Kq`. Tolerances (1e-25 numeric, exact symbolic) are
tight and the sample set avoids the closed-form pole at `Pi = kappa = pi/2`.

### F3 — insufficient_verification (banner)

**Classification:** resolved

**What changed:**

`mathematica/...mathematica_audit.wl:32` now reads
`banner["STAGE 146 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];`
and the transcript echoes the corrected stage number at line 3.

**Assessment:**

The stage number is corrected from 129 to 146. The literal still says
"FINITE-CORRECTION" rather than the directive's suggested "FIRST-ORDER",
but per the orchestrator note this is owned by Cluster A's mass-renumber
pass. The actual defect (wrong stage number) is fixed and the SymPy and
Mathematica transcripts both reflect Stage 146.

## Exec log assessment

**SymPy:** exit=0. Notable lines:

- `g(Pi) direct-formula = 0` (line 1) — F2 symbolic identity passes.
- `S_q(Pi) numeric sample Pi={7/10, 11/10, 17/10, 23/10}: diff=0E-165`
  (lines 2-5) and `PASS: S_q(Pi) direct-formula (numeric fallback at 4
  samples)` (line 6) — F2 numeric fallback fires and passes.
- `kernel check at Pi={1, 3/2, 5/2}: ... diff=0` (lines 24-26) — pre-existing
  numeric kernel check still passes.
- `Pi_* = 1.50882951349315552747043511772`, `g_* = 0.758035078944662822951954793784`,
  `S_* = 0.658075937605429271930315313437`, `g_*' = 0.0714453558083195219958997177670`,
  `S_*' = 0.0483709542125040994477978749511` (lines 27-31) — unchanged from
  prior runs and matching the notes' boxed values.
- F1 fallback (lines 156-161): all four residuals at `eps in {1/10, 1/2}`
  for `g_eps` and `S_eps` are 1e-18 to 1e-31, both PASS lines fire.
- `Stage 146 complete.` (line 163) — clean exit.

**Mathematica:** exit=0. Notable lines:

- Line 3: `STAGE 146 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS` (F3).
- Lines 7-10: `g(Pi) direct-formula = 0 / PASS` and `S_q(Pi) direct-formula = 0 / PASS` (F2).
- Lines 11-19: kernel formula samples (1, 3/2, 5/2) all diff=0, PASS.
- Lines 20-26: `Pi_* = 1.5088295134931555830055507559542749...`, then
  `g_*`, `S_*`, `g_*'`, `S_*'` printed; `Pi_* compensation point diff = 0`
  / PASS (pre-existing non-tautological assertion still holds).
- Lines 33-38: `g_eps affine law (integral form) at eps=1/10: 0`,
  `at eps=1/2: 0`, PASS; same for `S_eps`. F1 numeric-sample fallback fires
  cleanly.
- Line 40: `Stage 146 Mathematica audit passed.`

**Output freshness:**

- SymPy script mtime `2026-05-27 20:17:02`; SymPy output mtime
  `2026-05-27 20:21:14` (output post-edit by ~4 min).
- Mathematica script mtime `2026-05-27 20:02:10`; Mathematica output mtime
  `2026-05-27 20:03:01` (output post-edit by ~50 s).

Both `.txt` outputs are newer than their corresponding scripts.

## Material-change assessment

`material_change`: false.

No printed downstream constant changed across the rework. The SymPy and
Mathematica values for `Pi_*`, `g_*`, `S_*`, `g_*'`, `S_*'`, `gMinus`, and
the numeric kernel samples are identical to the audit-time values; the
symbolic `delta Pi` and `delta S` expressions are constructed from the
same `gPi`, `Sformula`, and `gminus` as before. Downstream units that
consume Stage 146 numerics see no perturbation; the orchestrator's
`upstream_stale` flag for units > 146 can stay off as far as this
verification is concerned.

## Side observations (non-blocking)

- The directive file `redteam/directives/stage_146.md` still does not have
  `## Applied: Fn` / `## Blocked: Fn` blocks. Recording the deviations
  (Mathematica F1 chop tol `10^-6`; SymPy F2 numeric fallback at
  `Pi in {7/10, 11/10, 17/10, 23/10}` with `1e-25`; SymPy F1 numeric-eps
  fallback at `eps in {1/10, 1/2}` with `1e-15`) in those blocks would make
  the deviations durable.
- The Mathematica F1 chop tolerance is `10^-6`, much looser than the SymPy
  F1 tolerance `1e-15`. The asymmetry is justified by the WP-9
  numeric-`pStar` complex-zero noise documented in the script comments;
  the actual residuals chop to literal `0`, so the headroom is unused in
  practice.
- `varsigma_test = 6*x*(1-x)` (SymPy) and `varsigmaTest = 6*x*(1-x)`
  (Mathematica) both integrate to 1 on [0,1] and are positive on (0,1),
  satisfying the convex-deformation family preconditions.
- The SymPy `expect_zero("g(Pi) direct-formula", ...)` symbolic branch on
  line 42 fires successfully for `g(Pi)`; only the `S_q(Pi)` branch needs
  the numeric fallback (`sp.integrate(Sigma*Kq, ...)` returns
  unevaluated). This is expected — the cosh-against-exponential integral
  is the harder one for SymPy.

## Verdict justification

All three findings are resolved. F1's intent — replacing tautological
distributive identities with integral-form checks that actually exercise
linearity of integration plus the `gMinus`/`Sformula` definitions — is
realized in both engines and runtime-confirmed via the documented
numeric-sample fallbacks (eps in {1/10, 1/2} at 1e-15 in SymPy, 10^-6 chop
in Mathematica). F2 anchors both closed-form moments to direct integrals
via symbolic identity (Mathematica fully; SymPy for `g(Pi)`) with a
4-sample numeric fallback for SymPy's `S_q(Pi)`. F3's banner stage number
is corrected. Both exec logs exit 0, no downstream-relevant numeric drift,
both transcripts are post-edit. Verdict: `verified`.
