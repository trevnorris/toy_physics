---
unit_id: 240
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 240

## Per-finding outcomes

### F1 — tautological_check (Omega_Q variable-independence trap)

**Classification:** resolved

**What changed:**
`scripts/...stage240..._sympy_audit.py:113-132`. The diff removes the two
no-op assertions `assert_zero(sp.diff(c0_expr, Omega_Q))` /
`assert_zero(sp.diff(c1_expr, Omega_Q))` (which differentiated expressions that
never contained `Omega_Q`) and replaces them with an extraction FROM the
`Omega_Q`-bearing precursor `Y_support`:

- `pole_denominator = 1 - omega**2/Omega_Q**2`
- guard `Y_support.has(Omega_Q)` (raises if false)
- `c1_probe = pole_denominator * Y_support`, guarded by `c1_probe.has(Omega_Q)`,
  then `c1_static = limit(c1_probe, omega, Omega_Q)`
- `static_sum = limit(Y_support, omega, 0)` (guarded `.has(Omega_Q)`), and
  `c0_static = static_sum - c1_static`
- then `assert_zero(sp.diff(c0_static, Omega_Q))`, `...(c1_static, Omega_Q)`,
  and the two cross-checks `c0_static - c0_expr == 0`, `c1_static - c1_expr == 0`.

**Assessment:**
Correct and non-vacuous. The differentiated objects are derived from
`Y_support`, which syntactically contains `Omega_Q` (the `.has(Omega_Q)` guards
on `Y_support` and on `c1_probe` fire BEFORE the limit, so the extraction path
provably touched `Omega_Q` prior to collapse — this is the exact self-test trap
the directive demanded). The limits are mathematically right:
`(1-ω²/Ω²)·Y_support → (alpha_req-alpha_mix)/alpha_req` as ω→Ω_Q, and
`Y_support|_{ω=0} = 1`, giving `c1_static = (alpha_req-alpha_mix)/alpha_req`,
`c0_static = alpha_mix/alpha_req`. These are NOT re-typed copies of
`alpha_mix/alpha_req` — they are produced by SymPy `limit` on the
pole-bearing object, and the `c0_static - c0_expr`/`c1_static - c1_expr`
cross-checks then tie the extraction back to the compiler weights. The
`diff(..., Omega_Q) == 0` is now a genuine zero: the limit collapsed a real
`Omega_Q` dependence, rather than differentiating an `Omega_Q`-free constant.
No collateral edits beyond the cited block (diff is scoped to lines 110-132).

### F2 — missing_verification_script (new Mathematica .wl, independent route)

**Classification:** resolved

**What changed:**
New file `mathematica/...stage240..._mathematica_audit.wl` (18 in-file checks
across M1-M6, all `expectZero`/`expectTrue` exiting nonzero on failure).

**Assessment — genuinely independent:**
- **M2 decomposition differs from the .py.** The SymPy script multiplies
  `pole_denominator * Y_support` and takes a limit; the Mathematica script
  instead calls `Apart[ySupport, omega]` to produce a true partial-fraction
  split (output shows two distinct simple-pole terms at `omega ± omegaQ`, NOT a
  re-typed rational), then extracts `c1` via `Limit[poleFactor*yApart, omega->omegaQ]`
  and `c0` via `Limit[yApart - c1Extracted/poleFactor, omega->0]`. That is a
  different choreography (Apart-based residue split vs. SymPy's
  `static_sum - c1_static`). The `FreeQ[poleProbe, omegaQ]` /
  `FreeQ[contactProbe, omegaQ]` guards enforce the same non-vacuousness as
  SymPy's `.has()` guards. Native primitives `Together`/`Apart`/`Limit`/
  `FullSimplify`/`Resolve` are used as the directive required, not transliterated
  SymPy calls.
- **M4 DERIVES Pi_tr, does not self-confirm.** `rhoMinimal = rhoMinimalFromC0`
  (= `1/(3/4)` from M3), then `piTrMinimal = FullSimplify[rhoMinimal*cMix]` and
  the check `piTrMinimal - (4/3)*cMix == 0`. Pi_tr is built from the M3-extracted
  loading ratio times C_mix; it is NOT defined as `(4/3)C_mix` and then
  "confirmed." This is exactly the fix to the SymPy A9 redundancy the auditor
  flagged. (Note: the SymPy A9 path was independently brought into line too —
  line 151 now derives `Pi_tr_selected = rho_min_from_c0 * C_mix` — though that
  was outside F1's named scope.)
- M1 ratio + spectral substitution, M3 specialization (4/3, 1/3), M5 selector
  reduction (`varrho = 2(1-eps_star)/3`), and M6 strict regime interval
  (`Resolve[1<4/3<2]` and `cMix < piTrMinimal < 2cMix`) all present and exact.

## Exec log assessment

**SymPy:** exit=0. 23 `[ok]` lines, 0 failures. Notable:
`[ok] Omega_Q independence of extracted c0`, `[ok] Omega_Q independence of
extracted c1`, `[ok] extracted c0 matches compiler c0`, `[ok] extracted c1
matches compiler c1` — the four new F1 labels all appear and pass.

**Mathematica:** exit=0. 18 `PASS:`, 0 `FAIL:`. The diagnostic dump confirms a
real Apart split: `M2 yApart = alphaMix/alphaReq + ((alphaMix-alphaReq)*omegaQ)/
(2*alphaReq*(omega-omegaQ)) - (...)/(omega+omegaQ)`, with
`c0Extracted = alphaMix/alphaReq`, `c1Extracted = 1 - alphaMix/alphaReq`.
Both `M2 Omega_Q independence` checks = 0. M4/M5/M6 all PASS.

**Pass counts:** SymPy 23/23 expected (10 prior + 4 new F1 + 9 unchanged tail);
Mathematica 18/18 expected (M1:2, M2:8, M3:3, M4:2, M5:1, M6:2) — every
`expectZero`/`expectTrue` in the file fired; no silent parser-skip. FAIL/guard
paths are non-vacuous (`fail` → `Exit[1]`; `assert_zero`/`.has()` raise).

**Output freshness:** confirmed. `_sympy_audit.txt` mtime 2026-06-03 08:42:38
and `_mathematica_audit.txt` mtime 08:42:38 are both newer than the source
`.py` (08:35:41) and `.wl` (08:36:32). Outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false. No derived constant or downstream-consumed result
changed. F1 re-derives the same `c0 = alpha_mix/alpha_req`,
`c1 = (alpha_req-alpha_mix)/alpha_req` by a non-vacuous route (the values are
identical, only the verification became real). F2 adds a second engine that
reproduces the existing constants (4/3, 1/3, 2(1-eps_*)/3). No unit > 240 needs
re-audit on account of unit 240.

## Side observations (non-blocking)

- The SymPy A9 path (lines 151-158) now derives `Pi_tr_selected` from
  `rho_min_from_c0 * C_mix` rather than re-typing `(4/3)C_mix`, so the
  auditor's non-blocking A9 redundancy note is effectively answered on both
  engines. Not a finding; just noting the redundancy no longer stands.
- Mathematica `expectZero` runs `FullSimplify[Together[Expand[...]]]` twice; the
  `M6` regime uses `expectTrue`/`Resolve` as intended. No correctness concern.

## Verdict justification

Both findings are resolved. F1 is genuinely de-vacuumed: the `Omega_Q`
independence is now established by extracting the weights from the
`Omega_Q`-bearing `Y_support` (guarded `.has(Omega_Q)` before the collapsing
limit) and cross-checking against the compiler weights, so the assertion can
actually fail if the static algebra were `Omega_Q`-dependent. F2 is a genuinely
independent second engine — an `Apart`-based partial-fraction extraction (not a
SymPy transliteration) with M4 deriving `Pi_tr` from `rho_alpha*C_mix` rather
than self-confirming. Both engines exit 0 with matching, non-vacuous check
counts, outputs are fresh, and no derived result changed (material_change
false).
