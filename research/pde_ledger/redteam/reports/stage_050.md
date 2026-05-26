---
unit_id: 050
batch: III.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md
  paper_appendix: present
---

# Audit unit 050 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_050.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 78; `\input` at line 218)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.txt`

## What the paper claims

The paper card's `\stagefield{Output}` says verbatim: "Lowest-twin sufficiency \eqref{eq:app-stage050-lowest-success} and higher-harmonic exclusion/softness thresholds." Four boxed equations crystallise the deliverables: (i) `S_0 = S(1; eps) = 2` for the lowest symmetric twin lane with `zeta_0 = 1` (eq:app-stage050-S0); (ii) `zeta_req <= 1 <=> S_req <= 2` (eq:app-stage050-lowest-success); (iii) the immediate exclusion test `zeta_req > 1/(2n+1)^2` for the nth twin harmonic (eq:app-stage050-higher-bound); (iv) the softness threshold `x <= x_max(n; zeta_req) = [1/((2n+1)^2 zeta_req) - 1]/[n(n+1)]` (eq:app-stage050-xmax). The paper card lists inputs as the Stage 048 `zeta_req` and the Stage 049 D/N support tower. The notes add a fifth result in Section 5 ("Exact enhancement bounds for the higher harmonics"): the ceiling `S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)`. The paper card's body and Output line do not mention this ceiling.

## What the script claims to verify

Both scripts (SymPy and Mathematica) walk through five algebraic checks: (1) `S(1; eps) - 2 = 0`, the doubling at `zeta = 1`; (2) the rational identity `zeta_req - 1 = (1 - eps)(S_req - 2) / (1 + eps(S_req - 2))`, which encodes the `zeta_req <= 1 <=> S_req <= 2` equivalence on the physical domain; (3) `d zeta_n^(twin)/dx` matches the closed form `-n(n+1)/[(2n+1)^2 (1 + xn(n+1))^2]`, and (in SymPy) `sp.solve` produces an `x_max` matching the paper's closed form, while (in Mathematica) `Solve` is also used; (4) the admissibility residual that recasts `x_max >= 0` as `(2n+1)^2 zeta_req <= 1`; (5) the enhancement ceiling identity `S_n^(max) - S_n^(twin)(x)` factors into a manifestly nonnegative form on `n >= 1, eps < 1, x > 0`. The SymPy script imports `twin_support_ratio` from the stage 049 module; the Mathematica script re-declares `zetaN` inline.

## Paper <-> script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `S_0 = S(1; eps) = 2` (eq:app-stage050-S0) | `expect_zero("S(1;eps) - 2", ...)` (sympy:47; wl:44) | partial — substitutes `zeta=1` by hand; never evaluates `twin_support_ratio(0, x)` to anchor it |
| `zeta_req <= 1 <=> S_req <= 2` (eq:app-stage050-lowest-success) | criterion identity (sympy:51-54; wl:47-50) | match |
| Higher-harmonic exclusion `zeta_req > 1/(2n+1)^2` (eq:app-stage050-higher-bound) | admissibility residual (sympy:80-85; wl:72-80) | match (via x_max sign equivalence) |
| Softness threshold `x_max` closed form (eq:app-stage050-xmax) | `x_eq` from Solve, matched to closed form (sympy:70-75; wl:53-70) | match |
| Enhancement ceiling `S_n^(max)` (notes section 5; absent from paper card) | sympy:88-107; wl:82-95 | extra — in notes only, not advertised in paper card Output |

Dominant pattern: most deliverables match; one is partial (the n=0 anchor) and one is extra relative to the paper card while present in the notes. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 47 | `expect_zero("S(1;eps) - 2", S.subs(zeta,1) - 2)` | doubling (eq S0) | partial — direct substitution, not anchored to `twin_support_ratio(0, x)` |
| A2 | sympy | 51-54 | `criterion - (1-eps)(S_req-2)/(1+eps(S_req-2)) == 0` | lowest-twin equivalence | yes |
| A3 | sympy | 65-68 | `d zeta_n/dx + n(n+1)/[(2n+1)^2(1+xn(n+1))^2] == 0` | monotonicity (sign manifest from explicit form) | yes |
| A4 | sympy | 72-75 | `x_eq - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)] == 0` (after `sp.solve`) | softness threshold | yes |
| A5 | sympy | 82-85 | admissibility residual = 0 | higher-harmonic exclusion | yes |
| A6 | sympy | 94 | `S_n^(twin)(x=0) - S_n^(max) == 0` | ceiling (notes only) | yes |
| A7 | sympy | 104-107 | factored difference of ceilings = 0 | ceiling (notes only) | yes |
| B1 | mathematica | 44 | `(sEnhance /. zeta -> 1) - 2 == 0` | doubling | partial — same direct substitution |
| B2 | mathematica | 47-50 | criterion identity | equivalence | yes |
| B3 | mathematica | 59-64 | derivative form match | monotonicity | yes |
| B4 | mathematica | 66 | `(zetaN /. x -> xEq) - zetaReq == 0` | softness threshold | yes |
| B5 | mathematica | 68-70 | `xEq - xEqClosedForm == 0` | softness threshold | yes |
| B6 | mathematica | 72-80 | admissibility residual = 0 | higher-harmonic exclusion | yes |
| B7 | mathematica | 86 | `(sN /. x -> 0) - sNMax == 0` | ceiling (notes only) | yes |
| B8 | mathematica | 92-95 | factored ceiling difference = 0 | ceiling (notes only) | yes |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:1-103`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:1-114`

**What's wrong:**
The `.wl` script is a near line-by-line port of the SymPy script. The same five-block choreography (doubling, equivalence identity, derivative form-match, Solve for `x_eq`, admissibility residual, ceiling factored form) appears in the same order with byte-for-byte equivalent target expressions. Concrete corresponding excerpts:

- sympy:42 `S = sp.simplify(1 + zeta * (1 - eps) / (1 - eps * zeta))` vs wl:40 `sEnhance = FullSimplify[1 + zeta (1 - eps)/(1 - eps zeta), ...]`
- sympy:64 `d_zeta_n_dx_target = -n * (n + 1) / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1)) ** 2)` vs wl:60 `dZetaNdxTarget = -n (n + 1) / ((2 n + 1)^2 (1 + n (n + 1) x)^2)` — identical hand-written derivative target used as a form-match
- sympy:74 `x_eq - (1/(((2*n+1)**2)*zeta_req) - 1)/(n*(n+1))` vs wl:56 + wl:68-70 — the same closed-form `xEqClosedForm` is hand-typed and subtracted from `xEq`
- sympy:100-103 `ceiling_diff_target = (1 - eps) * (2 * n + 1) ** 2 * n * (n + 1) * x / (((2 * n + 1) ** 2 - eps) * ((2 * n + 1) ** 2 * (1 + x * n * (n + 1)) - eps))` vs wl:89-91 `ceilingDiffTarget = ((1 - eps) (2 n + 1)^2 n (n + 1) x) / (((2 n + 1)^2 - eps) ((2 n + 1)^2 (1 + n (n + 1) x) - eps))` — identical factored target

In addition, the SymPy script imports `twin_support_ratio` from stage 049 (sympy:17), while the Mathematica script re-declares `zetaN = 1/((2n+1)^2(1 + x n(n+1)))` inline (wl:52). The Mathematica side does not anchor to an independent upstream representation — it just copies the same closed form. There is no independent re-derivation pathway (no Series expansion, no Limit/Maximize, no inverse-Solve of the equivalent inequality); both engines walk the same algebra with the same hand-written intermediate target expressions.

**Why this matters:**
Second-engine policy requires Mathematica to derive the claims independently so an error in the SymPy algebra (or in a hand-written target form) cannot be silently confirmed by a transliteration of the same algebra. Here the Mathematica `expectZero` calls subtract the same hand-typed targets used by SymPy; if any target form were mistyped, both engines would `PASS` identically and the engine cross-check would catch nothing.

**Required change:**
Rewrite the Mathematica script so at least three of the five load-bearing checks derive their reference expression independently rather than reusing the SymPy target form. Concretely:

- Derivative check (wl:59-64): instead of subtracting the hand-written `dZetaNdxTarget`, verify the algebraic structure of `D[zetaN, x]` by multiplying through and showing the denominator structure matches: `expectZero["d zeta_n / dx times denominator squared + n(n+1) = 0", D[zetaN, x] (2 n + 1)^2 (1 + n (n + 1) x)^2 + n (n + 1)]`. This pins the derivative to the denominator factorization of `zetaN` itself rather than to a transliterated closed form.
- `x_eq` closed-form match (wl:68-70): drop the `xEqClosedForm - xEq` subtraction and instead verify the equivalent statement `expectZero["x_eq from Solve satisfies (2n+1)^2 zeta_req (1 + n(n+1) x_eq) - 1 = 0", (2 n + 1)^2 zetaReq (1 + n (n + 1) xEq) - 1]`, which exercises the defining equation of `x_max` directly without re-typing the closed form.
- Ceiling factored-form check (wl:92-95): replace `ceilingDiff - ceilingDiffTarget` with an independent derivation, e.g., compute `Series[sNMax - sN, {x, 0, 1}]` and verify the leading coefficient is the expected positive expression, or factor via `Together[sNMax - sN]` and assert the resulting numerator equals `(1 - eps) (2 n + 1)^2 n (n + 1) x` (Mathematica derives the numerator by `Numerator[Together[...]]` rather than from a hand-typed target).

The keep-the-existing-PASS-output constraint: after the rewrite, every `expectZero` should still simplify to 0 under the same `$Assumptions`, but the target expressions in the new assertions must not be byte-equivalent to the SymPy target expressions.

**Verification:**
After fix, the Mathematica script should still exit 0 and the saved `.txt` output should still report all checks passing. A diff against the old `.wl` should show `dZetaNdxTarget`, `xEqClosedForm`, and `ceilingDiffTarget` removed in favour of structurally different reference expressions.

### F2 — paper_misalignment

**Subtype:** paper_missing_script_claim
**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_050.tex:44` (Output line; body equations end at line 42)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md:163-187` (Section 5 — enhancement bounds)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:88-112`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:82-97`

**What's wrong:**
The scripts verify a fifth claim — the higher-harmonic enhancement ceiling `S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)` — that is present in the notes (Section 5) but absent from the paper card's body and Output:

Paper card (stage_050.tex:44) Output line quote: `"Lowest-twin sufficiency \eqref{eq:app-stage050-lowest-success} and higher-harmonic exclusion/softness thresholds."` — no mention of an enhancement ceiling.

Notes (md:172-174) quote: `S_n^(twin) < S_n^(max) := 1 + (1-eps) / [ (2n+1)^2 - eps ].`

Script (sympy:90) quote: `S_n_max = sp.simplify(1 + (1 - eps) / ((2 * n + 1) ** 2 - eps))`, with assertions at sympy:94 and sympy:104-107.

Script (wl:83) quote: `sNMax = FullSimplify[1 + (1 - eps)/((2 n + 1)^2 - eps), ...]`, with assertions at wl:86 and wl:92-95.

This is a `paper_missing_script_claim`: the script tests a deliverable the paper card does not advertise. The notes do contain it, so it is not a fabrication, but the paper card's `\stagefield{Output}` line and the boxed equations in the body do not list it.

**Why this matters:**
Either the paper card is silently relying on the notes (in which case the Output line should be extended), or the script's ceiling block is verifying material that does not belong in this stage's deliverables (in which case the block should move to whichever stage owns the enhancement ceiling). Direction of resolution is the user's call.

**Required change:**
Routed to user via `## Resolve before fix_loop`. Do not auto-edit either side.

**Verification:**
After user resolves: if the paper card is extended with a fifth boxed equation `S_n^(twin) < S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)` plus an Output line referencing it, this finding closes with no script change. If the script's ceiling block is removed/migrated, the assertion inventory loses A6, A7, B7, B8 and the cross-check table's "extra" row becomes empty.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:47`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:44`

**What's wrong:**
The doubling theorem (paper eq:app-stage050-S0) rests on the upstream fact `zeta_0^(twin) = 1` (notes:84: `zeta_0^(twin) = 1.`), which originates in stage 049's `twin_support_ratio(n, x)` formula. The SymPy script imports `twin_support_ratio` from stage 049 (sympy:17) but then ignores the import when proving the doubling: it substitutes `zeta = 1` literally (`S.subs(zeta, 1) - 2`, sympy:47) without ever evaluating `twin_support_ratio(0, x)` to confirm it equals 1. The Mathematica script does not import an analogue at all and likewise substitutes `zeta -> 1` directly (wl:44).

As a result, the chain "n=0 twin gives zeta=1, therefore S_0 = S(1;eps) = 2" is broken into two unconnected halves in the script: the assertion `S(1;eps) - 2 = 0` is trivially provable by anyone willing to type `zeta = 1`, and is not anchored to the stage 049 closed form that the paper's claim relies on.

**Why this matters:**
A typo in the imported `twin_support_ratio` (or in any future change to the stage 049 module) would not be caught here. The doubling-at-n=0 chain is exactly what the paper card's first boxed equation asserts; the script should anchor that claim to the formula it already imports.

**Required change:**
Add a one-line anchor in the SymPy script just before sympy:47:

```python
expect_zero(
    "zeta_0^(twin) - 1 (anchors doubling to stage 049 import)",
    twin_support_ratio(sp.Integer(0), x) - 1,
)
```

The literal `sp.Integer(0)` is required because `twin_support_ratio` was declared with `n` a `positive=True` symbol; substituting a Python `0` into the symbolic expression is safe but using `sp.Integer(0)` makes the intent explicit and protects against assumptions changes.

In the Mathematica script, add an analogous anchor just before wl:44 (note that wl:36 sets `n >= 1` in `$Assumptions`, so the n=0 substitution should bypass the global assumption):

```mathematica
expectZero["zeta_0^(twin) - 1 (anchors doubling)", (1/((2 n + 1)^2 (1 + x n (n + 1))) /. n -> 0) - 1];
```

The substitution `n -> 0` is a literal pattern replacement; it does not interact with the integer-positivity assumption on `n` (the replacement happens before `FullSimplify` evaluates the expression), so the resulting residual `1/(1 (1 + 0)) - 1` will simplify cleanly to 0.

**Verification:**
After fix: the SymPy output should show a new line `"zeta_0^(twin) - 1 (anchors doubling to stage 049 import) = 0"` and the Mathematica output should show the analogous `PASS` line. If `twin_support_ratio` were ever broken so that the n=0 value drifted, this check would now fail at this stage rather than propagating silently.

## Independent-derivation check (Mathematica)

The Mathematica script is a near line-by-line transliteration of the SymPy script. Corresponding expressions are identical modulo Mathematica's syntax, and — most damning — the hand-written target forms inside `expectZero` (the explicit `dZetaNdxTarget`, the explicit `xEqClosedForm`, the explicit `ceilingDiffTarget`) are byte-for-byte the same expressions used as SymPy targets. There is no independent derivation pathway, and the Mathematica script does not import an analogue of `twin_support_ratio` from a stage 049 Mathematica module — it just re-declares the closed form inline. See F1.

## Engine cross-check

Both engines produce identical final structure: SymPy reports `S(1;eps) - 2 = 0`, `zeta_n^(twin) = 1/((2n+1)^2(nx(n+1)+1))`, admissibility residual `= 0`, `S_1^(max) = 2(eps-5)/(eps-9)`, `S_2^(max) = 2(eps-13)/(eps-25)`. Mathematica reports identical content modulo presentation: `S_1^(max) = 1 + (-1+eps)/(-9+eps)` (= `2(eps-5)/(eps-9)` after combining fractions), `S_2^(max) = 1 + (-1+eps)/(-25+eps)` (= `2(eps-13)/(eps-25)`). All `expectZero` calls return 0 on both sides; no engine disagreement. Caveat: because of F1, "engines agree" here really means "the two engines walked the same algebra and got the same answer." If the hand-written target forms contain a coordinated error, this engine cross-check would not catch it. That is the operational point of F1.

## Verdict justification

The math holds up under attack on its own terms: every algebraic identity the scripts assert (doubling, equivalence-identity, derivative form, x_max closed form, admissibility, ceiling factored form) simplifies to zero, and the simplifications are not hiding behind aggressive assumptions (positivity is used only where the paper's physical setup justifies it). I tried to break the equivalence direction by pushing `S_req < 1` into a regime where `1 + eps(S_req - 2)` could flip sign; on the physical domain stated in the notes (`S_req > 1`, `0 < eps < 1`) the multiplier `(1-eps)/(1+eps(S_req-2))` stays positive and the iff holds. I tried to break the monotonicity argument by checking the sign of `d zeta_n/dx`; the target form `-n(n+1)/[(2n+1)^2(1+xn(n+1))^2]` makes the sign manifest for `n >= 1, x > 0`. I tried the `n=0` corner; the script handles it by direct `zeta = 1` substitution rather than by evaluating the imported `twin_support_ratio(0, x)`, which is F3 (low severity — not an error, just a missing anchor).

The verdict is `findings` rather than `clean` because (F1) the Mathematica script is a line-by-line transliteration of the SymPy script and therefore does not satisfy the second-engine independence policy; (F2) the script verifies an enhancement-ceiling theorem that lives in the notes but is not advertised in the paper card; (F3) the doubling-at-n=0 chain is not anchored to the imported stage 049 formula. None of these mathematically invalidate any result of stage 050, so `stop_cold: null`. F2 is a paper_misalignment that needs user resolution before any edit; F1 and F3 are script-side and can be applied by Codex.

## Self-test notes

Checked variable independence (every `diff` is in `x`, and `zeta_n` genuinely depends on `x`, so no zero-derivative trap is hiding); checked the equivalence direction across the physical `S_req > 1` domain (denominator `1 + eps(S_req-2) > 1 - eps > 0` for `0 < eps < 1`); checked the n=0 substitution and confirmed `twin_support_ratio(0, x) = 1/((1)^2 (1 + x*0*1)) = 1`, so the F3 anchor would close as expected; checked output mtimes (sympy_audit.py @ 1779490749 < sympy_audit.txt @ 1779490871; mathematica.wl @ 1779492289 < mathematica.txt @ 1779492357), so no `stale_output` finding; verified that the F1 directive's proposed re-derivations target identical mathematical content via different code paths, so they do not introduce a new paper_misalignment.
