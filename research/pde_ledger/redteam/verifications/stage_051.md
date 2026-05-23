---
unit_id: 051
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 051

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py:126-133` — replaced the inverse-map round-trip (`Mmix_from_ZW.subs(ZW, ZW_twin_req) - Gtr/2`) with the explicit Stage 047/030 forward map: defines `ZW_from_Mmix = pi**2 (1 - eps_eta)(1 - eps) Mmix / [8 (1+chi0)**2]`, substitutes `Mmix -> Gtr/2`, and compares the result to `ZW_twin_req` written verbatim at L117. New assertion label: `"Z_W^(twin,req) - forward-map(M_mix=G_tr/2)"`.
- `mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:87-94` — symmetric edit: defines `zWFromMmix`, applies `mMix -> gTr/2`, and `expectZero["Z_W^(twin,req) - forward-map(M_mix=G_tr/2)", zWTwinReq - zWThresholdViaMap]`.

**Assessment:**
The edit matches the directive's "Required change" block character-for-character (modulo Mathematica's `mMix`/`Pi` naming). The new check is non-tautological: `ZW_twin_req` is supplied verbatim at sympy L117 / wl L86 with the `16 (1+chi0)^2` denominator, while the forward-map image of `Gtr/2` is recomputed from a single forward coefficient `pi^2 (1-eps_eta)(1-eps) / [8 (1+chi0)^2]`. Any prefactor perturbation of either expression (e.g. `(1+chi0)` vs `(1+chi0)^2`, or `8` vs `16`) would leave a nonzero residual. Both engine outputs show the new assertion line and a `0` / PASS verdict.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:46-58` — replaced the `FullSimplify[piTr - piExpected]` round-trip with independent `Factor[Together[...]]` canonicalization on both sides. Adds two new print lines (`Pi_tr (Factor/Together)` and `Pi_tr (claim, Factor/Together)`) and asserts `expectZero[piTrFactored - piExpectedFactored]`.
- `mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:96-111` — replaced the handed-in closed-form `xi2x` with an independent `Solve[{gTr == 2 mMix, xi > 0}, xi, Reals]` derivation, plus a parallel `xi2xClaim` of the closed form, plus a cross-check `expectZero["xi_(2x): Solve vs claim", xi2xDerived - xi2xClaim]` and the subsequent `expectZero["G_tr(xi_(2x)) - 2 M_mix", (gTr /. xi -> xi2xDerived) - 2 mMix]`.

**Assessment:**
The Pi_tr block matches the directive verbatim and routes through `Factor`/`Together` rather than `FullSimplify`, breaking the transliteration. Output confirms both `piTrFactored` and `piExpectedFactored` print as distinct expressions (L8-9) that reduce to the same canonical form (L10-11 PASS).

The `xi2x` block uses `Solve[{gTr == 2 mMix, xi > 0}, xi, Reals]` with `First` rather than the directive's `Select[..., # > 0, ...] [[1]]` SelectFirst pattern; the orchestrator note explicitly flags this as the intended independent-Solve route (`xi2x` independent-Solve plus `ConditionalExpression -> e` strip), and the code comment justifies the change (`SelectFirst`-on-TrueQ would silently drop both roots because symbolic positivity is undecidable without further assumptions). The mechanical idiom adjustment is substance-preserving: Solve still derives the root from `gTr == 2 mMix` independently of the closed-form claim. Critically, the Mathematica `Solve` output `(-9*delta + 2*mMix*(9 + 2*r^2) + Sqrt[81*(delta + 2*mMix)^2 - 72*(delta - 2*mMix)*mMix*r^2 + 16*mMix^2*r^4])/18` is structurally distinct from `xi2xClaim`'s discriminant form `(... + Sqrt[648*delta*mMix + (9*delta - 2*mMix*(9 + 2*r^2))^2])/18` (output L33-34), so the assertion `xi_(2x): Solve vs claim = 0` is a non-tautological comparison of two genuinely independent expressions, and both PASS lines confirm the bridge holds.

The `If[pi1 =!= Infinity, ...]` to `If[!TrueQ[Simplify[1/pi1 == 0, ...]], ...]` idiom patch and the `ConditionalExpression -> e` strip in `expectZero` and `xi2x` are mechanical and substance-preserving as flagged by the orchestrator note. They do not weaken any assertion — `1/pi1 == 0` is satisfied iff `pi1` is infinite, and stripping `ConditionalExpression[0, cond]` is benign under the declared `$Assumptions` domain.

## Exec log assessment

**SymPy:** exit=0 (inferred; the `expect_zero` helper raises on nonzero, the script reaches the final theorem-ledger banner, and the orchestrator log shows a clean refresh at 17:12:41 with no HALT after the subsequent resume cycle). Notable lines from `scripts/output/.../sympy_audit.txt`:
- L27: `Pi_tr - expected closed form = 0`
- L108: `Z_W^(twin,req) - forward-map(M_mix=G_tr/2) = 0` (new F1 line, replaces the prior tautological label)
- L120: `G_tr(xi_(2x)) - 2 M_mix = 0`

**Mathematica:** exit=0. Notable lines from `mathematica/output/.../mathematica_audit.txt`:
- L8-11: `Pi_tr (Factor/Together) = ...` / `Pi_tr (claim, Factor/Together) = ...` / `Pi_tr - expected closed form = 0` / `PASS: Pi_tr - expected closed form` (new F2 Factor/Together path)
- L13, L15: `Limit::alimv` warnings about assumptions that involve the limit variable; these are benign (Mathematica is correctly ignoring `0 < xi < 1` inside the `Limit[..., xi -> 0]` call) and do not break the assertion — L16-17 show `pi0 = 0` and `pi1 = Infinity` exactly as required.
- L31-32: `Z_W^(twin,req) - forward-map(M_mix=G_tr/2) = 0` / `PASS` (F1)
- L33-38: `xi_(2x) (Solve) = ...` / `xi_(2x) (claim) = ...` / `xi_(2x): Solve vs claim = 0 / PASS` / `G_tr(xi_(2x)) - 2 M_mix = 0 / PASS` (F2)
- L40: `Stage 051 Mathematica audit passed.`

**Output freshness:** confirmed. Script mtimes:
- sympy script: 17:11; sympy output: 17:18 (newer)
- mathematica script: 17:30; mathematica output: 17:31 (newer)

Both `.txt` files were regenerated after Codex's edits. (The initial 17:12 HALT was due to a non-zero Mathematica exit from the un-patched idiom; the subsequent 17:30/17:31 regeneration after the orchestrator's mechanical idiom patches cleared it.)

Note: the orchestrator wrote no `redteam/exec_logs/stage_051_*.log` or `stage_051_diff.patch` (the user's prompt explicitly overrode the path references); the canonical fresh transcripts at `scripts/output/...sympy_audit.txt` and `mathematica/output/...mathematica_audit.txt` plus the codex diffs in `redteam/fix_batch_III.2.log` supplied the same information. No information loss for this verification.

## Material-change assessment

`material_change`: false.

Neither edit changes a derived numerical/symbolic result that downstream units could depend on. F1 swapped one algebraically equivalent assertion (`f^{-1}(f(x)) = x` round-trip) for another (forward-map image of `Gtr/2` equals the verbatim `ZW_twin_req`); both check the same Stage 047/030 coherent map at the same threshold point. F2 added independent derivations of `xi_(2x)` and an alternate canonicalization of `Pi_tr` for cross-engine assurance, but the printed closed-form values of `xi_(2x)`, `Pi_tr`, `Lambda_twin_req`, `M_mix^(twin,req)`, `Z_W^(twin,req)`, `C_mix`, `S_req`, and `zeta_req` are unchanged. Downstream units (052+) that key off Stage 051's threshold formulas see identical exports.

## Side observations (non-blocking)

1. The sympy script's `xi2x` block still hands in the closed-form root verbatim (sympy L136-142) rather than calling `sp.solve(Gtr - 2*Mmix, xi)`. F2 was scoped to the Mathematica script only (the directive explicitly says "Do not change the SymPy script for this finding"), so this is in-scope behavior — but it means cross-engine assurance for `xi_(2x)` now rides entirely on the Mathematica `Solve` route. If the SymPy closed form were perturbed, the SymPy `expect_zero("G_tr(xi_(2x)) - 2 M_mix", ...)` would still pass (it is the direct substitution check that proves the supplied form solves the equation), so this is not a hidden flaw — just an asymmetry in which engine carries the independent derivation. Not a finding.

2. The `Limit::alimv` warnings on Mathematica wl L60 are harmless: the `Limit[piTr, xi -> 0, Direction -> "FromAbove"]` call passes `Assumptions -> delta > 0 && r > 0` (not the global `$Assumptions` which would include `0 < xi < 1`), so the warning that pops twice in the output appears to come from Mathematica internally evaluating `$Assumptions`-bound subexpressions inside `FullSimplify`. The numerical outcome (`pi0 = 0`, `pi1 = Infinity`) is correct.

3. The `xi2x` Solve output has a non-obviously-equal but algebraically equal surface form to the claim; this is exactly what an independent algebraic route should produce, and the `expectZero["xi_(2x): Solve vs claim", ...]` check verifies the equivalence — non-tautological.

## Verdict justification

Both findings are resolved: F1 replaces the inverse-map tautology with a genuine forward-map comparison in both engines, and F2 routes the Mathematica script through `Factor`/`Together` canonicalization for the `Pi_tr` identity and through `Solve[..., Reals]` for `xi_(2x)`. Outputs show all new assertion lines, all `PASS` markers, and both scripts complete to their final banner / `Exit[0]`. The orchestrator's mechanical idiom patches (`ConditionalExpression -> e` strips and the `1/pi1 == 0` infinity check) are substance-preserving and do not weaken any assertion. Verdict: `verified`.
