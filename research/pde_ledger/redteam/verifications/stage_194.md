---
unit_id: 194
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T12:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 194

## Per-finding outcomes

### F1 — missing_mathematica

**Classification:** resolved

**What changed:**
Codex created the new independent-engine script
`mathematica/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_mathematica_audit.wl`
(new untracked file; SymPy `.py` confirmed unmodified — `git status` shows only the `.wl` as `??`, the `.py` is not in the modified set, and `redteam/exec_logs/stage_194_diff.patch` is empty). The `.wl` re-verifies all four sections of the SymPy reference and emits 22 PASS lines, 0 FAIL, exit 0.

**Assessment:** correct, and a genuinely independent route, not a transliteration. Evidence by section:

- **Section I (Hankel fingerprint).** SymPy hand-assembles the mode from explicit rationals, `h2 = j2 + I*y2` with `j2,y2` written out in `sin/cos` (`.py` lines 44-46), then forms `Lambda_out = z*h2'/h2`. The `.wl` instead uses the Mathematica-native primitive `SphericalHankelH1[2, x]` via `FunctionExpand` (`.wl` line 63) — yielding the closed form `(I*E^(I*x)*(-3 + 3 I x + x^2))/x^3` in the output — and forms `lambdaExact = x*D[hankel2,x]/hankel2`. Different source object, native primitive, `Series[...,{x,0,7}]`+`Expand` vs SymPy's `sp.series(...,0,8)`. Both land the identical fingerprint `-3 + x^2/3 + x^4/9 + I x^5/9 - 2 x^6/27 - I x^7/27`.

- **Section II (retarded pole + χ_Q fixing).** This is the strongest independence evidence: the decomposition is *inverted*. SymPy *defines* `Omega_Q = 3 c_s/(2a)` and `sigma_Q_can = 9/(8 Omega_Q^5)` up front (`.py` lines 88-89) and then verifies the series matches. The `.wl` does the opposite — it `Solve`s for the pole from the ω² coefficient match (`Solve[Coefficient[retardedSeriesRaw,freq,2]==Coefficient[outgoingOmegaSeries,freq,2], poleScale]`, lines 107-113) and `Solve`s for sigma from the ω⁵ coefficient at `chiQ->1` (lines 115-121), each filtered through `positiveRoot`, then checks the *solved* values equal `3 c_s/2a`, `4 a^5/27 c_s^5`, and `9/(8 Omega_Q^5)`. Output confirms `poleFromEven = 3 soundSpeed/(2 radius)`, `sigmaFromOdd = 4 radius^5/(27 soundSpeed^5)`. χ_Q is isolated via `Coefficient[...,freq,5]` + `Reduce[ComplexExpand[oddRet/I==oddOut/I], chiQ]` → output `chiQ == 1` (line 35), where SymPy used the 5th `sp.diff`. Genuinely different extraction.

- **Section III (deformation algebra).** SymPy defines `L0..L5` algebraically (`.py` lines 135-138). The `.wl` builds the full operator `scaleS*(lambdaWindow5 /. x->betaStretch*x) + sigma0 + sigma2 x^2 + sigma4 x^4 + I sigma5 x^5` and *extracts* the slots via `Coefficient[deformedOperator, x, n]` (lines 172-175), then `Solve`s the two even-matching conditions for `{sigma2,sigma4}`. Both engines reproduce identical `Sigma_2/Sigma_4` formulas, the `chi_Q = 3(S β^5 + 9 Σ_5)/(3S - Σ_0)` map, and the `chi_Q-1` law.

- **Section IV (invariant tuple).** SymPy hardcodes `Kbar0 = 54 G c_s^5/(5 a^5 c^5)` (the auditor's soft observation). The `.wl` improves on this: it `Solve`s for `baseK0` from the `Gammabar_5 = 2G/(5 c^5)` relation (lines 233-238) and derives the remaining three, so even the carry-forward corollary is derived rather than asserted.

## Exec log assessment

**SymPy:** exit=0 (`redteam/exec_logs/stage_194_sympy.log` ends `# exit_code: 0`). All `expect_zero` lines print `= 0` (committed `.txt` Sections I-IV).

**Mathematica:** exit=0 (`stage_194_mathematica.log` ends `# exit_code: 0`; 22 `PASS:`, 0 `FAIL`). Notable lines:
- `positive pole from even matching = (3*soundSpeed)/(2*radius)` and `positive sigma from odd normalization = (4*radius^5)/(27*soundSpeed^5)` — derived by Solve, then PASS against the literals.
- `chi_Q condition from odd coefficient = chiQ == 1` — Reduce genuinely returns the fixing.
- `chi_Q from canonical-even deformed branch = (3*(betaStretch^5*scaleS + 9*sigma5))/(3*scaleS - sigma0))` — matches the SymPy map; `chi_Q deformation law = 0` PASS.

**Output freshness:** `.wl` mtime 1780335724 < `.wl` output mtime 1780335818 — output regenerated post-creation. SymPy `.py` mtime 1778522332 is older than its output (1780335818), consistent with the orchestrator's independent re-run refreshing both `.txt` files without touching the `.py`.

## Material-change assessment

`material_change`: false. No SymPy edit; the `.wl` is an additive independent confirmation that reproduces — does not alter — every Stage 194 result. No downstream-visible derived value changed.

## Side observations (non-blocking)

- Cosmetic only: the SymPy `banner`/`subbanner`/ledger strings still say "STAGE 177" (carried from the original report's second cosmetic note); the `.wl` correctly says "STAGE 194". Zero effect on any assertion. Not a basis to block.
- The `.wl` checks an extra relation SymPy did not (`sigma_Q^can - 9/(8 Omega_Q^5)`), strengthening rather than weakening coverage.

## Verdict justification

The single finding (missing second engine) is resolved with a genuinely independent Mathematica route: native `SphericalHankelH1` rather than hand-built `sin/cos`, and — most tellingly — an inverted decomposition that `Solve`s for the pole, sigma, and even Kbar_0 from coefficient/relation matching where SymPy asserted-then-verified them. Every load-bearing check is substantive (derived quantity vs independent literal, or two independently-derived quantities; no X−X, no hardcoded-vs-itself). Both engines exit 0 with all checks passing, outputs are fresh, the SymPy reference is untouched, and the cross-engine conclusions agree on the fingerprint coefficients, the χ_Q=1 fixing, the deformation laws, and the invariant tuple. Verdict: verified.
