---
unit_id: 208
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 208

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the new Mathematica audit at the registered target path
`mathematica/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_mathematica_audit.wl`
(248 lines). No other file was touched — the SymPy `.py` is unchanged and the
captured diff patch is empty because the new `.wl` is an untracked file
(`git status` reports it as `??`), so a tracked-file diff captures nothing.
The `## Applied: F1` block lists exactly this one file with `deviation: none`.
The orchestrator's independent re-run (exec log dated 2026-06-02T10:07:52)
exits 0 with 27 PASS / 0 FAIL; the saved `.txt` (mtime 10:08:00) is newer than
the `.wl` (mtime 10:04:18), so the output was regenerated post-fix.

**Assessment:**
The `.wl` independently establishes each manifest claim M1–M9 with substantive,
non-tautological `expectZero`/`expectTrue` checks. The helper `expectZero`
strips `ConditionalExpression[0,_]`, `FullSimplify`s under physical assumptions
(`ki>0, kj>0, r>=0, H0>0`), and `Exit[1]`s on any nonzero residual — so a wrong
closed form would fail the run rather than silently pass. Every manifest item is
represented by at least one explicit check:

- M1 → `M1 oriented slope law`, `M1 positive slope magnitude`
- M2 → `M2 slope derivative law`
- M3 → `M3 unique stationary root count`, `M3 stationarity region equals r==kj/ki`,
  `M3 gradient-optimal ratio`, `M3 gradient-optimal slope square`
- M4 → `M4 mixed curvature decomposition`, `M4 diagonal neutrality`,
  `M4 coefficient-recovered cross weight`
- M5 → `M5 cross-weight derivative`, `... at r=1`, `... value at r=1`
- M6 → six checks (gradient-optimal direction/slope/curvature, equal-mix direction/slope/curvature)
- M7 → `M7 ratio coincidence difference`, `M7 coincident directions when ki==kj`
- M8 → lower/upper Rayleigh form + lower/upper weighted form (4 checks)
- M9 → `M9 lower bracket closure quadratic`, `M9 upper bracket closure quadratic`

None of the checks is tautological: every derivative check (M2, M5) acts on an
expression genuinely depending on `r`; the M9 closure residual is a real
root-substitution into the closure quadratic (`H0 - k·τ + ½·κ·τ²`), not a
re-statement of the bracket formula; the M3 root count and `Equivalent` region
check confirm uniqueness rather than asserting the answer.

## Independent-derivation check (load-bearing)

Confirmed **independent** — native primitives, different decomposition, not a
line-by-line port. The two most port-prone spots, which the directive named as
the anti-transliteration guards, are derived by genuinely different routes:

1. **Cross weight (M4).** SymPy re-states the weight directly:
   `.py:89  w_x = sp.simplify(2 * r / (1 + r**2))`.
   Mathematica instead *recovers* it from the expanded quadratic form:
   `.wl:141-145  curvatureNumerator = Expand[Together[mixedCurvature*(1+r^2)]];
    linearHijNumerator = Normal[Series[curvatureNumerator,{hij,0,1}]];
    weightX = Coefficient[linearHijNumerator, hij]/(1+r^2)`.
   It then *checks* the recovered weight equals `2r/(1+r^2)`
   (`.wl:156 M4 coefficient-recovered cross weight`). All three weights
   (w_i, w_j, w_x) are obtained by `Coefficient` extraction, not re-stated.

2. **Gradient maximizer r_grad (M3).** SymPy substitutes the pre-baked ratio:
   `.py:68  r_grad = sp.simplify(kj / ki)`.
   Mathematica *solves* the stationarity equation:
   `.wl:106-116  Solve[{D[positiveSlope,r]==0, r>=0}, r, Reals]` plus an
   independent `Reduce[...]` of the stationarity region and an `Equivalent`
   check against `r==kj/ki` — corroboration that does not exist in the `.py`.

3. **Bracket closure (M9).** Both engines substitute the carried root map into
   the closure quadratic and simplify the residual to 0. The directive
   explicitly permitted this route ("Simplify-ing the residual to 0 OR
   Reduce-ing the quadratic"), and the auditor classified the corresponding
   SymPy check (A13) as genuine root-substitution rather than tautology. The
   `.wl` factors the map into a reusable `rootMap[h_,k_,c_]` and a
   `closureResidual[tau_,kappa_]` (`.wl:229-241`), a different code shape from
   the inline `Delta_lo`/`T_lo` construction in `.py:174-195`.

Additional structural divergence confirming this is not transliteration: the
`.wl` builds the gradient ray two independent ways and cross-checks them
(`sRay[rGrad]` from the solved root vs. `{ki,kj}/Sqrt[ki^2+kj^2]`,
`.wl:168-183`), whereas the `.py` builds `s_grad` directly and checks it against
`s.subs(r, r_grad)`; and the M8 envelope check verifies a Rayleigh form
(`sRay.loHessian.sRay`) in addition to the weighted form, using the
coefficient-recovered weights — neither present in the `.py`. Verdict on this
gate: **confirmed independent.**

## Exec log assessment

**SymPy:** n/a — no SymPy change in this finding; the `.py` was not edited.

**Mathematica:** exit=0. Notable lines:
- `Solve stationarity roots = {{r -> kj/ki}}` and `M3 unique stationary root count` PASS — maximizer derived, not assumed.
- `coefficient weights {w_i, w_x, w_j} = {(1+r^2)^(-1), (2*r)/(1+r^2), r^2/(1+r^2)}` then `M4 coefficient-recovered cross weight` PASS — cross weight recovered via Coefficient.
- `M9 lower/upper bracket closure quadratic` PASS — root-map substitution closes the quadratic.
- Footer `STAGE 208 MATHEMATICA AUDIT COMPLETED SUCCESSFULLY`, `# exit_code: 0`.

**Output freshness:** confirmed. `.wl` mtime 2026-06-02 10:04:18; output `.txt`
mtime 2026-06-02 10:08:00 (newer). 27 PASS, 0 FAIL in the saved `.txt`, matching
the exec log.

## Material-change assessment

`material_change`: false. The only edit is the addition of a brand-new
second-engine verification script. No derived value, closed form, or constant
consumed by downstream units (Stage 209's slope/curvature/bracket packet) was
altered — the `.wl` merely re-confirms the existing SymPy results on a second
engine. No downstream unit is affected.

## Side observations (non-blocking)

- The captured diff patch `stage_208_diff.patch` is empty (0 bytes). This is a
  capture artifact, not a defect: the new `.wl` is untracked (`git status` →
  `??`), so a tracked-file diff records nothing. The file demonstrably exists,
  is the sole changed artifact, and ran to exit 0.
- The SymPy `.py` banners still print "STAGE 191" (`.py:35,219`). The auditor
  already recorded this as cosmetic label drift in print strings only (no effect
  on assertions/paths/math) and explicitly did not raise it as a finding; the
  new `.wl` correctly labels itself "STAGE 208". Not blocking, not in scope here.

## Verdict justification

The sole finding F1 (missing Mathematica engine) is fully resolved: a substantive
`.wl` now exists at the registered path, exits 0 with 27 PASS / 0 FAIL, and
exercises every manifest item M1–M9 with non-tautological checks. The
load-bearing independence gate passes — the cross weight is recovered by
`Coefficient` extraction and `r_grad` by `Solve`/`Reduce`, the two routes the
directive named as the transliteration risks, and the script adds independent
corroborations (region `Reduce`, dual-construction gradient ray, Rayleigh-form
envelopes) absent from the SymPy script. No regressions, no material change.
Verdict: verified.
