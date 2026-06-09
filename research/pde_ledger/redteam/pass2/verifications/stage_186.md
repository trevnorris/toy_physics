---
unit_id: 186
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 186

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Classification:** resolved

**What changed:**
The `.wl` was re-authored (diff wl:9-61). The old hand-coded exponent lists
(`ctrCoreExponents`, `thermalExponents`, `cntPrefactorExponents`, `cntEExponents`,
`etaExponents`) and `logDriftFromExponents` are DELETED (grep confirms zero remaining).
In their place the physical micro-monomials are written as actual rational expressions
in micro-variables — `chi0Monomial = gamma0*cEtaU0/kU0`, `deltaUMonomial = Pi^2*tU0/(ell0^2*kU0)`,
`epsilonWMonomial`, `zOverOmegaMonomial`, `epsEtaMonomial = cEtaU0^2/(kU0*kEta0)` (wl:43-47),
then `ctrMonomial`/`cntMonomial` as products (wl:49-50). `logDriftFromMonomial` (wl:52-56)
rescales each micro-var `v0 -> v0*Exp[scaleS*drift]` (wl:41), takes `D[Log[scaled/orig], scaleS]`
at `scaleS->0`. `mMatDerived` is `Coefficient` of those DERIVED drifts (wl:77-80), then
checked against the hand-coded reference `mMat` (wl:83-85).

**Assessment:**
The `M_*` rows are now genuinely DERIVED from the symbolic monomial structure, not a
hand-coded-vs-hand-coded `Coefficient` round-trip. A wrong exponent in the reference `mMat`
would now fail (`mMatDerived` no longer originates from the same literals). Independently
re-traced all three rows: eps_eta row → `2*c1 - kU - kEta` = ref row3 `{0,2,0,-1,-1,0,0,0}`;
C_nt row → `{2+2E,0,2E,F-E,-1,-(2+E),1,-F}` = ref row2; C_tr row matches ref row1. Log lines
15-20 print all three "PASS: M_* row N from <monomial> matches reference". The micro→drift
positional map (wl:34-35) is faithful. `.py` `M_*` construction left unchanged per directive
(diff touches only py:114-133). Non-tautological now.

### F2 — insufficient_verification (both engines)

**Classification:** resolved

**What changed:**
SymPy (py:117-128): the old `eps_eta_logdrift = 2*C - U - eta_scaling` literal is replaced by
deriving the drift from the physical monomial `c_etaU0**2/(KU0*KEta0)` under exponential
rescaling (`c_etaU0->...exp(s*C)`, `KU0->...exp(s*U)`, `KEta0->...exp(s*eta_scaling)`), then
`diff(log(scaled/orig), s)|_{s->0}`. Mathematica (wl:153-168): reuses `logDriftFromMonomial[epsEtaMonomial]`
with `{c1->cSym, kU->uSym, kEta->etaScaling, rest->0}` instead of the old literal. Both print
`derived eps_eta log drift = ...` (sympy.log:42, math.log:55). The misleading
"Non-tautological ground check" comment is corrected to "Ground check: derive the K_eta^eff
scaling from the physical eps_eta monomial" in BOTH engines.

**Assessment:**
The `2`/`-1` coefficients now come from the monomial's exponent structure, not a `2C-U`-shaped
literal; a wrong exponent in `eps_eta = c_etaU^2/(K_U K_eta)` would change the drift and break
the `= 2C-U` match. Genuine ground check in both engines; comment no longer inverts the signal.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`delta log eps_eta = 2*c1 - kEta - kU`; `minor det(tau1,kEta,mu1) = chi + 1`;
`derived eps_eta log drift = 2*C - U - eta_scaling`; `K_eta preserving scaling matches paper 2C-U = 0`;
all `finite orbit preserves C_tr/C_nt/eps_eta = 0`; all `linearized {tau,kEta,mu} formula = 0`.

**Mathematica:** exit=0. Notable lines:
`PASS: M_* row 1/2/3 from <monomial> matches reference`; `minor det(tau1,kEta,mu1) = 1 + chi`;
all `PASS: solved {K_eta,T_U,mu_W} scaling matches paper`; `PASS: finite orbit preserves ...`;
`derived eps_eta log drift = 2*cSym - etaScaling - uSym`; `Stage 186 Mathematica audit passed.`
Banner = "STAGE 186 — EXACT MICROSCOPIC SIMILARITY ORBIT" (correct).

**Output freshness:** confirmed. Output `.txt` mtimes (1780986351 for both) are newer than the
`.py` (1780985652) and `.wl` (1780985666) source mtimes — regenerated post-fix.

## Material-change assessment

`material_change`: false. Both edits are method-only (how the central object and the ground check
are derived). All 9 deliverable values are preserved bit-for-bit vs the pre-fix run: `M_*` rows 1-3,
minor det `1+chi`, the three compatibility formulas `(tau_1, kappa_eta, mu_1)`, the orbit scaling
laws, monomial preservation. No derived result changed, so no downstream unit is affected.

## Side observations (non-blocking)

- The `dim ker M_* = 5` deliverable remains not explicitly emitted (rank-3 via nonzero minor on
  8 columns implies it). This was logged informationally by the auditor (not a finding) and is
  unchanged here — not a verification blocker.
- The new positivity assumptions on micro-variables (`positiveMicroAssumptions`, wl:36) are
  needed for `PowerExpand[Log[...]]` to be valid and are physically appropriate (couplings > 0).

## Verdict justification

Both findings are resolved. F1's re-author makes the `M_*` row-match genuinely derive each row
from the physical monomial definitions under primitive-variable rescaling, so a wrong exponent
now fails — the old hand-coded exponent lists are gone and all three rows re-traced correctly.
F2's eta-scaling round-trip is replaced in BOTH engines by a derivation from the physical
`eps_eta = c_etaU^2/(K_U K_eta)` monomial, and the misleading comment is corrected. Both scripts
exit 0 with every in-file check PASS, outputs are fresh, the `.py` changed only for F2, and all
nine deliverable values are preserved (method-only change). Verdict: verified.
