---
unit_id: 125
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage125_positive_source_theorem.md"]
  paper_appendix: present
---

# Audit unit 125 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_125.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage125_positive_source_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 29 audit-path summary; line 86 paragraph survey; `\input` line 1284)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.txt`

## What the paper claims

The stage card quotes a single boxed Output (line 16): "Positive normalized sources satisfy `0 ≤ 𝔤[σ] ≤ 1`, ruling out the upper compensated branch." The notes make this concrete with two distinct deliverables. First, the **cosine moment representation** (notes lines 33–39, boxed):
`𝔤[σ] = ∫_0^L σ(z) cos(π z/(2L)) dz`,
i.e., the mouth-bias factor is the first cosine moment of the positive normalized axial source profile on the first D/N throat interval, with `σ(z) ≥ 0` and `∫_0^L σ dz = 1`. Second, the **range bound** (notes lines 46–55, boxed): since `0 ≤ cos(π z/(2L)) ≤ 1` for `z ∈ [0,L]`, every positive normalized source has `0 ≤ 𝔤[σ] ≤ 1`. Third, the **branch-selection corollary** (notes lines 60–81, three boxed statements): the explicit Family-1 carry-forwards `𝔤_-^F1 ≈ 0.75803507894466` and `𝔤_+^F1 ≈ 2.79795199200529` then give `𝔤_+^F1 > 1` (upper branch ruled out) and `0 < 𝔤_-^F1 < 1` (lower branch is the unique admissible compensated Family-1 branch). The card's `Checks` list (lines 21–25) names three guard-rails: positivity of the mouth source, zero-flux/boundary-layer normalizations in the GNLS/localized-Maxwell reduction, and Family-1 compensation point against the lower (not the singular equal-normalized) branch.

## What the script claims to verify

Both engines verify, in two pieces. (1) **Pointwise kernel bound**: define `kernel = cos(π x/2)` on `x ∈ [0,1]` (which is the kernel `cos(π z/(2L))` after `x = z/L`); compute the minimum and maximum of this kernel on `[0,1]` and assert `min == 0`, `max == 1`. (2) **Explicit Family-1 branch values and inequalities**: define `r = sqrt(12·(37/20)² / π² − 1) = sqrt(4107 − 100π²)/(10π)`, define the two branches `g_± = (2·sqrt(4107 − 100π²) ± 37·sqrt(3))/(20π)`, assert that the closed form for `r` is consistent (`r − sqrt(4107 − 100π²)/(10π) == 0`), assert that both branches satisfy the balance relation `1 + r² − 4(g − r)² == 0`, and finally assert the three numerical inequalities `g_- > 0`, `g_- < 1`, `g_+ > 1`. Concluding print statements narrate that the upper compensated branch is impossible and the lower branch is unique under any positive localized mouth source law.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Cosine moment representation `𝔤[σ] = ∫_0^L σ(z) cos(πz/(2L)) dz` (notes lines 33–39, boxed) | (none — the integral is never set up or computed for any σ; the kernel is defined and printed but never integrated) | missing |
| Pointwise kernel bound `0 ≤ cos(πz/(2L)) ≤ 1` for `z ∈ [0,L]` (notes lines 46–49) | `expect_zero("kernel minimum on [0,1]", ...)` + `expect_zero("kernel maximum on [0,1] - 1", ...)` (sympy lines 41–44; mathematica lines 44–47) | match |
| Integral range bound `0 ≤ 𝔤[σ] ≤ 1` for any positive normalized σ (notes lines 51–55, boxed Output) — the paper card's headline claim | (not asserted on any σ family; not asserted as an integral; the connecting convex-combination step from kernel bound to integral bound is implicit/in prose only) | missing |
| Carry-forward Family-1 values `𝔤_-^F1 ≈ 0.7580…`, `𝔤_+^F1 ≈ 2.7980…` (notes lines 60–63) | Closed forms `g_± = (2·sqrt(4107 − 100π²) ± 37·sqrt(3))/(20π)` plus numerical evaluation (sympy 51–52, 61–62; mathematica 55–56, 66–67) | match (closed forms are stated, not derived; values agree with notes to printed precision) |
| Balance relation tying the two branch values to a single quadratic | `1 + r² − 4(g − r)² == 0` for both `g_-` and `g_+` (sympy 55–56; mathematica 59–60) | extra (not in the paper card or notes, but internally consistent; ties both branch values to one identity, useful sanity) |
| `𝔤_+^F1 > 1` ⇒ upper branch ruled out (notes lines 65–72, boxed) | `expect_true("g_+ > 1", ...)` (sympy 66; mathematica 71) | match |
| `0 < 𝔤_-^F1 < 1` ⇒ lower branch admissible (notes lines 74–81, boxed) | `expect_true("g_- > 0", ...)` + `expect_true("g_- < 1", ...)` (sympy 64–65; mathematica 69–70) | match |
| `Checks` line "zero-flux and boundary-layer normalizations in the GNLS/localized-Maxwell reduction" (paper card line 23) | (none — no normalization or zero-flux check in this script) | missing (but this is a heuristic guard-rail from `Checks`, not an `Output` deliverable; carry-forward responsibility, not a stage-125 verifiable) |

Set `paper_alignment: partial` — five of seven deliverables match cleanly (kernel bound, branch values, both inequalities, balance-relation extra), but the cosine moment representation and the integral range bound (the paper's stated `Output`) are not asserted in the script. The `Checks` zero-flux line is properly upstream; not a script-side gap for this unit.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43 | `expect_zero("kernel minimum on [0,1]", kernel_min)` | Pointwise kernel lower bound | yes |
| A2 | sympy | 44 | `expect_zero("kernel maximum on [0,1] - 1", kernel_max - 1)` | Pointwise kernel upper bound | yes |
| A3 | sympy | 54 | `expect_zero("r - sqrt(4107 - 100 pi^2)/(10 pi)", r - R/(10π))` | Internal consistency of `r` closed form | partial (algebraically forced by construction; see Finding discussion) |
| A4 | sympy | 55 | `expect_zero("lower branch balance relation", 1 + r**2 - 4*(g_- - r)**2)` | Branch values satisfy single quadratic | yes (non-tautological — a sign or factor error in `g_-` would break this) |
| A5 | sympy | 56 | `expect_zero("upper branch balance relation", 1 + r**2 - 4*(g_+ - r)**2)` | Branch values satisfy single quadratic | yes |
| A6 | sympy | 64 | `expect_true("g_- > 0", ...)` | Lower branch admissibility | yes |
| A7 | sympy | 65 | `expect_true("g_- < 1", ...)` | Lower branch admissibility | yes |
| A8 | sympy | 66 | `expect_true("g_+ > 1", ...)` | Upper branch ruled out | yes |
| B1 | mathematica | 46 | `expectZero["kernel minimum on [0,1]", kernelMin]` | same as A1 | yes |
| B2 | mathematica | 47 | `expectZero["kernel maximum on [0,1] - 1", kernelMax - 1]` | same as A2 | yes |
| B3 | mathematica | 58 | `expectZero["r - sqrt(4107 - 100 Pi^2)/(10 Pi)", r - rrad/(10*Pi)]` | same as A3 | partial |
| B4 | mathematica | 59 | `expectZero["lower branch balance relation", 1 + r^2 - 4*(gminus - r)^2]` | same as A4 | yes |
| B5 | mathematica | 60 | `expectZero["upper branch balance relation", 1 + r^2 - 4*(gplus - r)^2]` | same as A5 | yes |
| B6 | mathematica | 69 | `expectTrue["g_- > 0", gminus > 0]` | same as A6 | yes |
| B7 | mathematica | 70 | `expectTrue["g_- < 1", gminus < 1]` | same as A7 | yes |
| B8 | mathematica | 71 | `expectTrue["g_+ > 1", gplus > 1]` | same as A8 | yes |

Notes on A3/B3: this check tests `sqrt(12·(37/20)²/π² − 1) == sqrt(4107 − 100π²)/(10π)`, which is forced by elementary algebra (`12·(37/20)² = 4107/100`, then both sides equal `sqrt((4107 − 100π²)/(100π²))`). It is essentially a sanity-print of the two parametrizations of `r`. It does exercise that the user-typed coefficients (`12, 37, 20, 4107, 100`) are mutually consistent — so a typo in any of them would break the assertion — and it is therefore not pure tautology (e.g., if someone wrote `4108` instead of `4107`, the assertion would fail). Marked "partial" rather than "no". Does not rise to a `tautological_check` finding.

A4/A5/B4/B5 are also non-tautological in the same partial sense: they exercise that the closed-form expressions for `g_±` are the two roots of the quadratic `4(g − r)² = 1 + r²`. The quadratic itself is *not derived* in this stage; it is implicitly carried forward from an earlier Family-1 derivation. Within stage 125's scope, the assertions test internal consistency of the carry-forward, which is the appropriate scope for a "branch-admissibility" sanity ledger.

A1/A2/B1/B2 are the only assertions that exercise an honest pointwise property of the cosine kernel. They are non-tautological — a sign flip or wrong half-period would change the min/max on `[0,1]`.

A6–A8/B6–B8 are numerical inequalities on closed-form algebraic quantities; non-tautological.

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** script_missing_paper_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_125.tex:15-17` — the boxed Output statement:
  > Positive normalized sources satisfy `0 ≤ 𝔤[σ] ≤ 1`, ruling out the upper compensated branch.
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage125_positive_source_theorem.md:33-39` — the boxed cosine moment representation:
  > `𝔤[σ] = ∫_0^L σ(z) cos(π z/(2L)) dz`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage125_positive_source_theorem.md:46-55` — the boxed integral range bound:
  > Because `0 ≤ cos(π z/(2L)) ≤ 1` for `z ∈ [0,L]`, every positive normalized source law satisfies `0 ≤ 𝔤[σ] ≤ 1`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:35-47` — the kernel-only block
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:38-51` — the kernel-only block

**What's wrong:**

The paper's headline `Output` (the only quoted-block claim in the card) is an integral inequality: `0 ≤ 𝔤[σ] ≤ 1` *for any positive normalized source* `σ(z) ≥ 0` with `∫_0^L σ dz = 1`. The notes spell this out as a two-step argument: (a) the cosine moment representation `𝔤[σ] = ∫_0^L σ(z) cos(π z/(2L)) dz`, then (b) the pointwise bound `0 ≤ cos(π z/(2L)) ≤ 1` ⇒ the integral bound by `σ ≥ 0` and `∫σ = 1`.

The script verifies only step (b) at the pointwise level: it defines `kernel = cos(π x/2)` and checks `min == 0`, `max == 1` on `x ∈ [0,1]`. It never sets up the integral `𝔤[σ] = ∫ σ·kernel dz`, never enforces `σ ≥ 0` or `∫σ = 1` as symbolic constraints, and never asserts `0 ≤ 𝔤[σ] ≤ 1` for any concrete or symbolic σ.

Concretely the script (sympy lines 35–47):

```
k = sp.pi / (2 * L)
kernel = sp.cos(sp.pi * x / 2)
...
kernel_min = sp.calculus.util.minimum(kernel, x, sp.Interval(0, 1))
kernel_max = sp.calculus.util.maximum(kernel, x, sp.Interval(0, 1))
expect_zero("kernel minimum on [0,1]", kernel_min)
expect_zero("kernel maximum on [0,1] - 1", kernel_max - 1)
```

and identically in `.wl` lines 44–47.

The convex-combination step (kernel-bound + positive-normalized σ ⇒ integral-bound) is genuinely a one-liner *mathematically*, but it is the step that takes the script's pointwise fact to the paper's stated integral claim. With it absent, the script's load-bearing assertion is "the kernel is bounded in [0,1] on [0,1]"; the paper's load-bearing claim is "every positive normalized source gives `𝔤[σ] ∈ [0,1]`". These are different identities. A reader running the script sees a kernel-bound pass; they do not see a check of the paper's Output statement.

The remainder of the script (A3–A8 / B3–B8) only exercises the carry-forward Family-1 branch values and their pairwise quadratic + inequalities. Those pieces correctly support the *corollary* (`𝔤_+^F1 > 1` ⇒ upper ruled out; `0 < 𝔤_-^F1 < 1` ⇒ lower admissible), but the corollary relies on the integral bound, which is what is missing.

**Why this matters:**

Stage 125 is the unit that admits "positive localized sources" as a model class and proves that `𝔤[σ]` is constrained to `[0,1]` for that class. Downstream (Stages 133–145 per the paper card line 27) consume this constraint to fix the lower compensated branch as canonical. A future edit to the script — e.g., changing the kernel's half-period, changing the integration interval to `[−L, L]` (where `cos(π z/(2L))` is no longer nonneg over the whole interval), or changing the normalization convention — would not be caught by the existing pointwise check. A direct integral assertion against a one-parameter family of σ (or even just the delta and uniform endpoints) would close that gap.

There is also a representational asymmetry that should be either fixed or explicitly elided: the script defines two kernels — `cos(π z/(2L))` is printed informatively (sympy line 39 / .wl line 42), but the *asserted* kernel is `cos(π x/2)` on `x ∈ [0,1]`. The variable substitution `x = z/L` is correct and the bounds on `cos(π x/2)` on `[0,1]` are the same as `cos(π z/(2L))` on `[0,L]`, but the printed `Cos[k z]` form and the asserted `kernel` are never linked symbolically (e.g., via `kernel == sp.cos(k*z).subs(z, L*x)`). This is a smaller crack in the same wall.

**Required change:**

(Routed to user via the directive's `## Resolve before fix_loop` block — Codex must not auto-edit either paper or scripts to resolve this until the user picks a direction.)

The natural fix is to add an explicit moment-bound check in each script. For sympy, after the existing kernel-bound block, define a representative positive normalized source profile and assert the moment lies in [0,1]. A robust pattern that exercises the integral while staying closed-form is a beta-like one-parameter family `σ_a(z) = (a + 1) (z/L)^a / L` (positive and normalized for `a ≥ 0`), or any concrete profile such as the uniform `σ = 1/L` and the endpoint delta `σ = δ(z)` (handled symbolically via the moment-of-delta limit). The cleanest variant is the *symbolic* statement: introduce a symbolic positive σ via assumptions, set up the integral `g_func = sp.integrate(sigma*kernel, (z, 0, L))`, and assert lower/upper bounds using `sp.refine` under `sigma >= 0` and `sp.integrate(sigma, (z, 0, L)) == 1`. SymPy will not in general prove this directly, so the working version is to test on a parametric family. The mathematica mirror uses `Integrate[...]` plus the same family.

Direction (a) below in the directive lays out two concrete options: (i) parametric family with assertion over a range of the family parameter, (ii) explicit pointwise-implies-integral lemma stated as a one-line algebraic identity (the integral of `(kernel - 1)·σ` is `≤ 0` by positivity).

**Verification:**

After Codex applies (under user direction (a)), both `.py` and `.wl` should contain at least one new `expect_zero` / `expectZero` (or `expect_true` / `expectTrue`) assertion naming the integral `𝔤[σ]` for at least one concrete positive normalized σ, and asserting that the integral lies in `[0,1]`. The new check should appear in the saved output as a `= 0` or `= True` line. The paper card line 16 then has a direct script-side counterpart.

### F2 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:33-66`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:35-71`

**What's wrong:**

The Mathematica script is a line-by-line port of the SymPy script. Identical variable choreography, identical closed-form expressions, identical assertion sequence, identical numerical inequalities. Three concrete corresponding-section comparisons (out of many possible):

Sympy lines 33–44 vs Mathematica lines 35–47 — kernel and bounds:

```
# sympy
z, L = sp.symbols("z L", positive=True, real=True)
x = sp.symbols("x", real=True)
k = sp.pi / (2 * L)
kernel = sp.cos(sp.pi * x / 2)
...
kernel_min = sp.calculus.util.minimum(kernel, x, sp.Interval(0, 1))
kernel_max = sp.calculus.util.maximum(kernel, x, sp.Interval(0, 1))
expect_zero("kernel minimum on [0,1]", kernel_min)
expect_zero("kernel maximum on [0,1] - 1", kernel_max - 1)
```

```
(* mathematica *)
Clear[x, z, L];
$Assumptions = L > 0 && Element[{x, z}, Reals];
k = Pi/(2*L);
kernel = Cos[Pi*x/2];
...
kernelMin = FullSimplify[MinValue[{kernel, 0 <= x <= 1}, x], Assumptions -> $Assumptions];
kernelMax = FullSimplify[MaxValue[{kernel, 0 <= x <= 1}, x], Assumptions -> $Assumptions];
expectZero["kernel minimum on [0,1]", kernelMin];
expectZero["kernel maximum on [0,1] - 1", kernelMax - 1];
```

Sympy lines 49–52 vs Mathematica lines 53–56 — branch closed forms:

```
# sympy
r = sp.simplify(sp.sqrt(sp.Rational(12, 1) * sp.Rational(37, 20) ** 2 / sp.pi**2 - 1))
R = sp.sqrt(4107 - 100 * sp.pi**2)
gminus = sp.simplify((2 * R - 37 * sp.sqrt(3)) / (20 * sp.pi))
gplus = sp.simplify((2 * R + 37 * sp.sqrt(3)) / (20 * sp.pi))
```

```
(* mathematica *)
rrad = Sqrt[4107 - 100*Pi^2];
r = FullSimplify[Sqrt[12*(37/20)^2/Pi^2 - 1], Assumptions -> $Assumptions];
gminus = FullSimplify[(2*rrad - 37*Sqrt[3])/(20*Pi), Assumptions -> $Assumptions];
gplus = FullSimplify[(2*rrad + 37*Sqrt[3])/(20*Pi), Assumptions -> $Assumptions];
```

Sympy lines 54–56 vs Mathematica lines 58–60 — balance assertions:

```
# sympy
expect_zero("r - sqrt(4107 - 100 pi^2)/(10 pi)", r - R / (10 * sp.pi))
expect_zero("lower branch balance relation", 1 + r**2 - 4 * (gminus - r) ** 2)
expect_zero("upper branch balance relation", 1 + r**2 - 4 * (gplus - r) ** 2)
```

```
(* mathematica *)
expectZero["r - sqrt(4107 - 100 Pi^2)/(10 Pi)", r - rrad/(10*Pi)];
expectZero["lower branch balance relation", 1 + r^2 - 4*(gminus - r)^2];
expectZero["upper branch balance relation", 1 + r^2 - 4*(gplus - r)^2];
```

Both engines plug in the same closed forms for `g_±` and check the same quadratic identity. Neither derives `g_±` independently from, e.g., a `Solve[1 + r^2 - 4*(g - r)^2 == 0, g]` (mathematica) vs. `sp.solve(1 + r**2 - 4*(g - r)**2, g)` (sympy) — which would constitute genuinely independent re-derivation paths. The mathematica script is essentially `s/sp\./Mathematica /` on the sympy script.

That said: the algebraic content of this stage is narrow (kernel min/max on a closed interval + numerical inequalities on closed-form algebraic numbers), so the space of "natural independent derivation paths" is small. The transliteration is real but its forensic value is correspondingly limited. Filing as low severity.

**Why this matters:**

The second-engine policy exists to catch typos or convention slips in one engine that the other would not commit. A copy-paste mirror cannot catch a shared coefficient typo (`4107 → 4108`, `37/20 → 37/21`) that exists in both files. The risk is low here because both files were presumably written together from the same notes, but the cross-engine guarantee is correspondingly weak.

**Required change:**

Restructure either the sympy or the mathematica script (Codex's choice, but mathematica is typically the one rewritten) so that at least the Family-1 branch derivation step uses an independent route. Concretely, in the `.wl`, replace lines 55–56 with:

```
branchSolutions = Solve[1 + r^2 - 4*(g - r)^2 == 0, g];
gminus = FullSimplify[Min[g /. branchSolutions], Assumptions -> $Assumptions];
gplus  = FullSimplify[Max[g /. branchSolutions], Assumptions -> $Assumptions];
```

then add an explicit symbolic-equality check against the closed forms that sympy uses:

```
expectZero["g_- closed-form match", gminus - (2*rrad - 37*Sqrt[3])/(20*Pi)];
expectZero["g_+ closed-form match", gplus  - (2*rrad + 37*Sqrt[3])/(20*Pi)];
```

This way the Mathematica script derives `g_±` by symbolic root-finding (independent of the sympy hand-write) and verifies that its derived values match sympy's hand-written closed forms — true cross-engine corroboration on the carry-forward.

**Verification:**

The updated `.wl` should contain `Solve[1 + r^2 - 4*(g - r)^2 == 0, g]` (or equivalent independent derivation) for `g_±`, plus the two `g_± closed-form match` assertions. Output transcript should show the two new `= 0` lines. All existing assertions should continue to pass.

## Independent-derivation check (Mathematica)

See F2 above. The `.wl` is a transliteration of the `.py` — same variables, same closed forms, same assertion sequence, same conclusion text. The only structural deltas are language-mechanical (`MinValue/MaxValue` instead of `sp.calculus.util.minimum/maximum`; `FullSimplify` decorations on every line; `$Assumptions = L > 0 && Element[{x, z}, Reals]` instead of sympy positivity assumptions). No independent re-derivation of `g_±` from the underlying quadratic. Filed as F2 (low severity).

## Engine cross-check

Both scripts pass all assertions (every `... = 0` and `... = True` line in the saved transcripts). Final printed identities match across engines:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| kernel min on [0,1] | `0` | `0` |
| kernel max on [0,1] − 1 | `0` | `0` |
| `r − sqrt(4107 − 100π²)/(10π)` | `0` | `0` |
| `1 + r² − 4(g_- − r)²` | `0` | `0` |
| `1 + r² − 4(g_+ − r)²` | `0` | `0` |
| `g_-` numeric | `0.75803507894466282692` | `0.75803507894466282691968…` |
| `g_+` numeric | `2.7979519920052934101` | `2.79795199200529341011158…` |
| `g_- > 0` | `True` | `True` |
| `g_- < 1` | `True` | `True` |
| `g_+ > 1` | `True` | `True` |

Engines agree on every value to printed precision (sympy 20 digits, mathematica 20-digit precision marker). `engines_agree: true`.

Output mtimes: sympy output 2026-05-11 > sympy script 2026-05-11 (sympy output is newer than script by ~48 min — fresh). Mathematica output 2026-05-11 13:10 > mathematica script 2026-04-22 (fresh by weeks). `outputs_fresh: true`.

## Verdict justification

The script's assertions are internally consistent and pass: the cosine kernel `cos(π x/2)` does have min `0` and max `1` on `[0,1]`; the closed forms `g_± = (2·sqrt(4107 − 100π²) ± 37·sqrt(3))/(20π)` do both satisfy `1 + r² = 4(g − r)²` (verified by hand: `g_± − r = ±37·sqrt(3)/(20π)`, square is `3·37²/(400π²) = 4107/(400π²)`, four times that is `4107/(100π²) = 1 + r²` ✓); the numerical inequalities `g_- ∈ (0,1)`, `g_+ > 1` are checked correctly. Within the scope of "verify the kernel's pointwise bound and the carry-forward Family-1 branch values," the script is sound. Two issues prevent a `clean` verdict:

1. **F1 (paper_misalignment / script_missing_paper_claim, medium):** The paper's stated `Output` is the *integral* inequality `0 ≤ 𝔤[σ] ≤ 1` for any positive normalized σ, derived from the *moment representation* `𝔤[σ] = ∫_0^L σ(z) cos(π z/(2L)) dz` (notes lines 33–55). The script verifies only the pointwise kernel bound; it never sets up the integral or asserts the bound on any σ. The convex-combination step from pointwise to integral is the heart of the theorem stated by the paper card, and it is absent.

2. **F2 (mathematica_transliteration, low):** The `.wl` is a line-by-line port of the `.py`. Both hand-write the same closed forms for `g_±` rather than deriving them by independent symbolic root-finding. The forensic value of a second engine is correspondingly reduced.

Attacks tried that found no further issues:
- (a) Verified `12·(37/20)² = 12·1369/400 = 16428/400 = 4107/100`, so `r² = 4107/(100π²) − 1 = (4107 − 100π²)/(100π²)` and `r = sqrt(4107 − 100π²)/(10π)` ✓. The script's A3/B3 sanity check is non-trivial in the typo-detection sense (a wrong coefficient anywhere would break it) but algebraically forced once the coefficients are fixed.
- (b) Verified that `4107 − 100π² ≈ 4107 − 986.96 ≈ 3120.04 > 0`, so the square root is real and the branch values are real numbers. Both `g_±` are real.
- (c) Verified numerical evaluation against the notes' quoted values: `g_- ≈ 0.758035078944663…` (notes line 62) matches sympy `0.75803507894466282692` and mathematica `0.75803507894466282691968…`. `g_+ ≈ 2.79795199200529…` (notes line 63) matches sympy `2.7979519920052934101` and mathematica `2.79795199200529341011158…`. ✓
- (d) Checked symbol assumptions in sympy (`z, L` positive real; `x` real) and in mathematica (`L > 0`, `{x, z} ∈ Reals`). Both are consistent with the physical setup. No `simplify` under aggressive assumptions hides anything (the kernel and branch values are real-analytic in their declared domains).
- (e) Checked the variable rebinding `kernel = cos(π x/2)` on `x ∈ [0,1]` vs. the printed `cos(π z/(2L))` on `z ∈ [0,L]`: these are related by `x = z/L`, and the bounds on the kernel are invariant under this rescaling. Not a math finding, but the asserted form and the printed-as-physically-meaningful form are never symbolically linked — noted in F1.
- (f) Checked the cosmetic banner "STAGE 108" (sympy line 31; mathematica line 33) against the actual stage number 125. This is documentation drift from an earlier stage-numbering, mirrored in both engines. The output transcripts reflect the same banner text. Does not fit any of the ten finding categories cleanly; not filed, but worth noting on the next renumber pass.
- (g) No CRITICAL_DOWNSTREAM concern: F1's resolution does not change any numerical or symbolic value carried forward (the branch values `g_±` are not altered, the bound `[0,1]` is not altered); it adds a missing assertion that exercises the integral statement. F2's resolution is a derivation-route swap that preserves identical outputs. Downstream stages 133–145 consume the *fact* `0 < g_- < 1 < g_+`, which both findings leave intact.

## Self-test notes

- **Variable independence:** no derivatives are taken in the script; not relevant here. The proposed F1 fix uses `sp.integrate` / `Integrate` over `z`; the integrand `σ(z)·cos(π z/(2L))` does depend on `z` for any non-degenerate σ, so the integral is non-zero in general.
- **Symmetry/parity:** the proposed F1 fix integrates a positive σ against a positive kernel over `[0, L]`; both factors are nonneg on the domain, so the integral is nonneg by construction — no parity-zero trap. Upper bound `≤ 1` follows from `kernel ≤ 1` and `∫σ = 1`.
- **Trivial-case pre-check:** for the proposed F1 fix, the simplest σ is the uniform `σ = 1/L`; the integral is then `(1/L)·∫_0^L cos(π z/(2L)) dz = (1/L)·(2L/π) = 2/π ≈ 0.6366`, which is in `(0,1)` as required. ✓ For the delta limit `σ → δ(z − 0)`, the integral is `cos(0) = 1` (boundary case). ✓ For `σ → δ(z − L)`, the integral is `cos(π/2) = 0` (boundary case). ✓ The bound is tight, so a parametric family that interpolates between these endpoints exercises the full `[0,1]` range. No risk of trivially-zero collapse.
- **Path specifications:** F1 and F2 both edit existing scripts in `scripts/` and `mathematica/`; no `missing_verification_script` finding, so no new file path to specify.
- **Paper round-trip:** F1's proposed fix uses the integral `∫_0^L σ(z) cos(π z/(2L)) dz` and the bound `0 ≤ g[σ] ≤ 1`, both of which appear verbatim in the notes (lines 33–55) and indirectly in the paper card. No new paper_misalignment introduced. F2's proposed fix uses the same balance relation `1 + r² = 4(g − r)²` and the same closed forms for `g_±`; no constants are changed; no new paper_misalignment introduced.
- **Stale-output check:** sympy output mtime 1778525132 > sympy script mtime 1778522211 (output is ~48 min newer than script — fresh). Mathematica output mtime 1778526601 > mathematica script mtime 1776812696 (output is ~19 days newer than script — fresh). The "STAGE 108" banner is in the script source itself and is faithfully reflected in the output, so the output is current. No `stale_output` finding.
