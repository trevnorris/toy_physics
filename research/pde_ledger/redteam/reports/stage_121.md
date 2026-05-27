---
unit_id: 121
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage121_geometric_r_selection.md]
  paper_appendix: present
---

# Audit unit 121 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_121.tex` (present; the stage is titled "Stage 138: Geometric Selection of the Core Hybridization Ratio")
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage121_geometric_r_selection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1276 input only; no separate row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py`
- mathematica: (missing — explicitly waived by the stage card: "Mathematica audit: none yet." and by the audit unit instructions)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.txt`
- mathematica output: (missing — n/a)

## What the paper claims

The stage card body (line 16-17) gives the boxed `Output`: "Identifying the mixed D/N tube with the throat span gives `r_geom(L/a) = sqrt(12 (L/a)^2 / pi^2 - 1)`." The notes elaborate this into four concrete deliverables: (1) the exact geometric branch law `r_geom(L/a) = sqrt(12/pi^2 (L/a)^2 - 1)` derived from setting `L_W = L` in the Stage-99 mixed-tube length formula `L_W = (pi a/2) sqrt((1 + r^2)/3)`; (2) the existence condition `L/a >= pi/(2 sqrt 3) ~ 0.9069`; (3) the family-1 numerical value `r_F1 = sqrt(4107 - 168 pi^2)/(10 pi) ~ 1.77799353547498` and `r_c^F1 = r_F1^2 ~ 3.16126101219081` at the carried preferred aspect ratio `L/a = 37/20`; (4) the immediate mixed-tube consequence `Omega_W = pi c_s/(2 L)`, the first D/N half-wave of the throat corridor. The part-04 appendix simply includes the stage card via `\input{stages/stage_121}` (line 1276) — no separate appendix row to cross-check.

## What the script claims to verify

The script (41 lines) declares positive-real symbols `L, a, L_W, r`; constructs `L_W_formula = pi a/2 sqrt((1+r^2)/3)`; uses `sp.solve` to invert for `r`, then substitutes `L_W -> L` to obtain `r_geom_simplified`. It prints this closed form, then uses `expect_zero` to assert `LW_formula.subs(r, r_geom_simplified) - L == 0` — a round-trip check that the solved `r_geom` satisfies the original tube-length equation when `L_W = L`. It then substitutes `L = (37/20) a` and prints (but does NOT assert) the resulting `r_F1` value, the squared `r_c(F1)` value, and the threshold `pi/(2 sqrt 3)`.

## Paper - script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `r_geom(L/a) = sqrt(12(L/a)^2/pi^2 - 1)` (closed-form branch law) | Tube-length round-trip assertion at L:25-28 implicitly exercises it (any wrong branch from `sp.solve` would fail the round-trip) | partial (only via inverse round-trip; no direct closed-form assertion) |
| Existence threshold `L/a >= pi/(2 sqrt 3)` | `threshold = sp.simplify(sp.pi/(2*sp.sqrt(3)))` printed at L:39-40; not asserted | missing (printed only) |
| `r_F1 = sqrt(4107 - 168 pi^2)/(10 pi) ~ 1.77799353547498` at L/a = 37/20 | `r_F1 = sp.simplify(r_geom_simplified.subs({L: R*a}))` printed at L:31-33; numeric printed; not asserted against any target | missing (printed only) AND mismatch (notes give `168 pi^2`, script computes `100 pi^2`) |
| `r_c^F1 = r_F1^2 ~ 3.16126101219081` | `rc_F1 = sp.simplify(r_F1**2)` printed at L:35-37; not asserted | missing (printed only) |
| `Omega_W = pi c_s/(2 L)` (first D/N half-wave of throat corridor) | not exercised; symbol `c_s` not defined in the script | missing |

`paper_alignment: partial` — the central closed-form law is verified (via round-trip), but four of the five boxed paper-side deliverables are not asserted (three are printed only; one is absent entirely). Additionally, the notes form `sqrt(4107 - 168 pi^2)/(10 pi)` algebraically disagrees with the script's `sqrt(4107 - 100 pi^2)/(10 pi)`, although both yield the same numerical 1.77799... (see F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 25-28 | `expect_zero("tube-length relation", LW_formula.subs(r, r_geom_simplified) - L)` | Deliverable 1 (closed-form branch law), via round-trip | partial — verifies that the solved `r_geom` satisfies the original equation, but does not assert the explicit closed form `sqrt(12 L^2 - pi^2 a^2)/(pi a)`. The round-trip catches a wrong root from `sp.solve`, but not a typo in the printed closed form. |

There is only ONE assertion. Lines 30-37 (`r_F1`, `r_c(F1)`) and line 39-40 (existence threshold) are pure prints — they cannot fail. The script's output transcript reports `PASS` based on a single non-tautological but narrow check.

## Findings

### F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage121_geometric_r_selection.md:64-71` (boxed `r_F1` expression)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:30-33` (computation of `r_F1`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.txt:15` (printed result)

**What's wrong:**

Notes (line 67-71) quote: "r_F1 = sqrt(12/pi^2 (37/20)^2 - 1) = sqrt(4107 - 168 pi^2)/(10 pi) ~ 1.77799353547498."

Script output (line 15) prints: "r_F1 = sqrt(4107 - 100*pi**2)/(10*pi)" with numeric "1.7779935354749781185".

The algebraic forms are NOT equal: `sqrt(4107 - 168 pi^2)/(10 pi)` evaluates numerically to ~1.575 (since 168 pi^2 ~ 1658.1, 4107 - 1658.1 ~ 2449, sqrt ~ 49.49, /10 pi ~ 1.575), whereas `sqrt(4107 - 100 pi^2)/(10 pi)` evaluates to ~1.7780 (since 100 pi^2 ~ 986.96, 4107 - 986.96 ~ 3120, sqrt ~ 55.86, /10 pi ~ 1.778).

Direct derivation: at L/a = 37/20, `r^2 = 12 (37/20)^2 / pi^2 - 1 = 12 * 1369 / 400 / pi^2 - 1 = 16428/(400 pi^2) - 1 = (16428 - 400 pi^2)/(400 pi^2) = (4107 - 100 pi^2)/(100 pi^2)`, so `r = sqrt(4107 - 100 pi^2)/(10 pi)`. The script is mathematically correct; the notes have a typo (the boxed numerical value 1.77799... matches the script, confirming the notes' literal expression is what's miscopied).

**Why this matters:**

A reader following the notes' algebraic form will arrive at a different (and wrong) value than the script asserts. Anyone downstream who quotes `sqrt(4107 - 168 pi^2)/(10 pi)` symbolically will propagate the wrong constant.

**Resolution required (user-routed):**

Almost certainly the notes typo `168 pi^2 -> 100 pi^2`. But the orchestrator must surface the question; Codex may not silently edit prose. See `## Resolve before fix_loop` in the directive.

**Verification:**

After user-directed fix to the notes, the boxed expression reads `sqrt(4107 - 100 pi^2)/(10 pi)`, matching the script's printed form and the boxed numerical 1.77799...

### F2 — insufficient_verification

**Severity:** high

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:30-40`

**What's wrong:**

The script has exactly ONE `expect_zero` assertion (lines 25-28). The remaining "verification" content (lines 30-40) consists entirely of `print` statements computing `r_F1`, `r_c(F1)`, and the existence threshold without ever asserting them against the paper-side target values.

Specifically:
- Line 31: `r_F1 = sp.simplify(r_geom_simplified.subs({L: R*a}))` — computed but not asserted against any target form.
- Line 35: `rc_F1 = sp.simplify(r_F1**2)` — printed but not asserted to equal `(4107 - 100 pi^2)/(100 pi^2)` or `4107/(100 pi^2) - 1`.
- Line 39: `threshold = sp.simplify(sp.pi/(2*sp.sqrt(3)))` — printed but not asserted to satisfy `r_geom_simplified.subs(L, a*threshold) == 0` (the existence boundary, where `r = 0`).
- Line 32: `print("r_F1 numeric =", sp.N(r_F1, 20))` — prints `1.7779935354749781185`, but not asserted against the notes' boxed `1.77799353547498`.

These are exactly the kind of `print`-only "verification" that flunks the audit's anchoring test: any of them could silently print a wrong value (e.g., a future typo at line 30 changing `R = sp.Rational(37,20)` to `R = sp.Rational(37,10)` would print a different number and still report `PASS`).

**Why this matters:**

Four of the five paper-side boxed deliverables (closed-form law, threshold, `r_F1`, `r_c^F1`) have no script-side assertion. The single assertion that exists only verifies the consistency of `sp.solve` with its own input, not the explicit closed form, the threshold, or any numerical target.

**Required change:**

Add at least four `expect_zero` assertions, named and anchored:

1. **Closed-form branch law (direct).** Assert `r_geom_simplified**2 - (12*L**2/(sp.pi**2 * a**2) - 1) == 0`. This pins the explicit closed form of the boxed branch law, independent of `sp.solve`'s branch choice.
2. **Family-1 numerical value (symbolic).** Define `r_F1_target = sp.sqrt(4107 - 100*sp.pi**2)/(10*sp.pi)` and assert `sp.simplify(r_F1 - r_F1_target) == 0`.
3. **Family-1 squared (symbolic).** Define `rc_F1_target = sp.Rational(4107, 100)/sp.pi**2 - 1` and assert `sp.simplify(rc_F1 - rc_F1_target) == 0`.
4. **Existence threshold.** Assert `r_geom_simplified.subs(L, a*threshold).simplify() == 0` (at the threshold, `r = 0`).

These four checks each compare an independently-computed quantity to a paper-side target; none is tautological with the existing A1.

Detailed instructions in the directive.

**Verification:**

After the fix lands, the sympy transcript must show four new `expect_zero` lines printing residual `0`, in addition to the existing tube-length relation. Script must still exit `0`.

### F3 — paper_missing_script_claim

**Subtype:** script_missing_paper_claim

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage121_geometric_r_selection.md:83-91` (boxed `Omega_W = pi c_s/(2 L)`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py` (nothing — `c_s`, `Omega_W` are not declared)

**What's wrong:**

The notes contain a boxed third deliverable: "With L_W = L, the auxiliary mixed-tube pole is simply Omega_W = pi c_s/(2 L), the first D/N half-wave of the actual throat corridor." This appears in a `\boxed{...}` block in the notes and is enumerated alongside `r_geom` and `r_F1` as one of the stage's outputs. The script does not declare `c_s` or `Omega_W` and contains no related check.

**Why this matters:**

The mixed-tube pole substitution is a load-bearing identification (it ties the mixed-tube pole to the actual throat geometry rather than a detached cavity). Leaving it unasserted means a downstream consumer who quotes `Omega_W = pi c_s/(2 L)` has no script-side anchor.

**Required change:**

Add an `expect_zero` block that constructs the Stage-99 mixed-tube pole `Omega_W_formula(L_W) = pi c_s/(2 L_W)` (or whatever Stage-99 supplies; if Stage-99's pole is stated as `Omega_W = pi c_s/(2 L_W)` then identifying `L_W = L` gives `Omega_W = pi c_s/(2 L)`), and assert the substitution at `L_W -> L`:

```python
c_s, OmegaW = sp.symbols('c_s Omega_W', positive=True, real=True)
OmegaW_formula = sp.pi * c_s / (2 * LW)
OmegaW_at_LW_eq_L = OmegaW_formula.subs(LW, L)
expect_zero(
    "Omega_W identification at L_W = L",
    OmegaW_at_LW_eq_L - sp.pi*c_s/(2*L)
)
```

Note: this is essentially a definitional check, but it makes the substitution explicit in the assertion log and matches the notes' boxed claim.

**Verification:**

After fix, the sympy transcript must show one additional `expect_zero` line "Omega_W identification at L_W = L" with residual `0`. Script must exit `0`.

## Independent-derivation check (Mathematica)

N/A — Mathematica script is intentionally absent (stage card: "Mathematica audit: none yet."). Per the audit unit instructions, this absence is not a finding for this stage.

## Engine cross-check

N/A — only one engine present.

Output mtime: sympy script `Apr 18 16:42` < sympy output `May 11 12:45`. `outputs_fresh: true`. Output content matches what the current script would produce (verified by reading both).

## Verdict justification

The single existing assertion (`tube-length relation`) is non-tautological and correctly verifies (via round-trip) that the solved `r_geom` satisfies the geometric identification `L_W = L`. However, the audit fails on two fronts:

1. **Coverage (F2, F3):** Four of the five boxed paper-side deliverables (the explicit closed form of `r_geom`, the existence threshold, the F1 numerical values `r_F1` and `r_c^F1`, and the mixed-tube pole `Omega_W`) are not anchored by any assertion. Three are printed only; the fourth is not exercised at all. The `PASS` status conveys far less verification than the paper claims.

2. **Notes typo (F1):** The notes' boxed symbolic form `sqrt(4107 - 168 pi^2)/(10 pi)` does not equal `sqrt(4107 - 100 pi^2)/(10 pi)` algebraically, even though both purport to equal ~1.77799. The script is correct; the notes are wrong (almost certainly a typo). User resolution required because the framework does not auto-edit prose.

The cosmetic banner "STAGE 104 - GEOMETRIC SELECTION OF r" (script L:16; output L:11) is a relabeling artifact and is not a math finding; the underlying physics carries the right `L/a = 37/20` constant. Not filed.

Attacks tried that did NOT yield findings: (a) checked whether `sp.solve` could return the negative branch of `r = +/- sqrt(...)`; with `r` declared `positive=True`, sympy returns the positive root, which is consistent with the physics (`r >= 0`). (b) checked symbol-assumption hazards on `L_W -> L` substitution; both are positive-real, no domain crossing. (c) checked whether the round-trip check is implicitly tautological; it is not — `sp.solve` could in principle simplify `(1 + r^2)/3 = (2 L /(pi a))^2` to a wrong branch, and the assertion would catch it.

## Self-test notes

(1) Variable independence: no derivatives in the script, n/a. (2) Symmetry/parity: no integrals, n/a. (3) Trivial-case pre-check for the four proposed new assertions: at `L = a * pi/(2 sqrt 3)` (the threshold), `r_geom_simplified` should reduce to `0` (the square root argument vanishes); confirmed by direct substitution `12/(pi^2) (pi/(2 sqrt 3))^2 - 1 = 12 pi^2 / (12 pi^2) - 1 = 0`. For the proposed F2-assertion #1, `r_geom_simplified^2 = (12 L^2 - pi^2 a^2)/(pi^2 a^2) = 12 L^2 /(pi^2 a^2) - 1`, matches target — no tautology trap because `r_geom_simplified` was obtained via `sp.solve`, not by squaring the target. For F2-assertion #2, `sp.simplify(sqrt(4107 - 100 pi^2)/(10 pi) - sqrt(4107 - 100 pi^2)/(10 pi))` is structurally `0`; the test catches a sign or coefficient change in `R = Rational(37,20)` because that would change the residual to a nonzero algebraic expression. (4) Path specs: not a missing-script finding. (5) Paper round-trip: the proposed F2 target `sqrt(4107 - 100 pi^2)/(10 pi)` deliberately does NOT match the notes' (wrong) form `sqrt(4107 - 168 pi^2)/(10 pi)`. The directive flags F1 as a paper_misalignment requiring user resolution; the F2 fix uses the correct form (which matches the notes' boxed numerical value).
