---
unit_id: 018
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T18:06:25-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: misaligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 018 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_018.tex`
- notes: (none)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 58 + `\input` row 115)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.txt`

## What the paper claims

Stage 018 is the "Parent throat action bundle master." The card's `\stagefield{Output}` is verbatim: "Stage~018 exports the parent-wall bundle bridge \eqref{eq:stage018-msigma}--\eqref{eq:stage018-ksigma}." Those two equations define the branch integrals

- `M_Sigma = ∫ dw μ_η β_2²` (eq:stage018-msigma)
- `K_Sigma = ∫ dw [T_w (β_2')² + (K_η + 6 T_Ω) β_2²]` (eq:stage018-ksigma)

and the body text says these "replace abstract wall knobs in the grouped bundle" with the caveat that "the support and mixed packets still have to be computed on the same frozen branch." The card has no `\stagefield{Inputs}` line, no other body equations, and no enumeration of additional deliverables. The part-01 appendix row 58 summarises the stage as "Bundle-level parent-action identities used by the projected electromagnetic response." No `notes/stages/moving_throat_pde_stage018_*.md` file exists; the docstring of the SymPy script references a `step_16_parent_throat_action_bundle_master_notes.md` that is not present anywhere in the repository (only a related output `.out` exists under `scripts/output/stage_004_em_projected_source/004_005/`).

## What the script claims to verify

The SymPy script (mirrored by the Mathematica script) verifies a much larger set of algebraic identities. From the assertions and the docstring "Master-note audit for step_16_parent_throat_action_bundle_master_notes.md":
1. One-pole numerator equivalence `(u4 - 4 u2²) = (D0(B4+Z4) - 3(MSigma+B2+Z2)²)/D0²` for an even quartic pole with `D0 = KSigma - B0 - Z0`, `D2 = -(MSigma+B2+Z2)`, `D4 = -(B4+Z4)`.
2. Two `KSigma` closure routes (one-pole closure `K = B0 + Z0 + 3(M+B2+Z2)²/(B4+Z4)` and normalization closure `K = B0 + Z0 + N0/Ptarget`) and a compatibility identity tying them via `N0`.
3. Linearized "even-gate" 2×2 system with determinant `1/27` whose solve produces the wall-stiffness slope `dKSigma = B01+Z01+27(B41+Z41)` and wall-inertia slope `dMSigma = -(B21+Z21)+3(B41+Z41)`.
4. Residual amplitude `Xi1 = N01/N0 - 27(B41+Z41)/(KSigma-B0-Z0)` after substituting the slopes.
5. Gaussian-profile reductions `M_Sigma = √π` and `K_Sigma = 3√π/2` for `β = exp(-w²/2)` with `μ_η = T_w = K_η+6 T_Ω = 1`.

Of these five claim families, only family 5 ("M8 Gaussian inertia/stiffness integrals") directly corresponds to the paper's `\stagefield{Output}` — and even then it tests a specific concrete profile, not the abstract bridge identity the card states. Families 1–4 are not mentioned anywhere in the stage-018 paper card or the part-01 appendix row.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `M_Sigma = ∫ dw μ_η β_2²` (eq:stage018-msigma) | A16 (sympy:94) / M8 (math:127,131) — Gaussian example `∫ β² dw = √π` with μ_η implicitly 1 | partial — only one concrete profile; abstract bridge identity is not exercised |
| `K_Sigma = ∫ dw [T_w (β_2')² + (K_η + 6 T_Ω) β_2²]` (eq:stage018-ksigma) | A17 (sympy:95) / M8 (math:128,135) — Gaussian example `∫ ((β')² + β²) dw = 3√π/2` with T_w = K_η + 6 T_Ω = 1 | partial — single profile; the abstract `T_w`, `K_η`, `T_Ω` coefficients are collapsed to 1 |
| One-pole numerator identity (script families 1–2) | A1–A5 (sympy:30–40) / M1–M3-mut (math:38–58) | extra — nothing in paper card |
| Compatibility cross-closure (script family 2 cont.) | A6 (sympy:44–45) / M4 (math:60–66) | extra — nothing in paper card |
| Even-gate determinant + wall-slope solve (script family 3) | A7–A13 (sympy:62–75) / M5–M6 (math:84–108) | extra — nothing in paper card |
| Xi1 residual (script family 4) | A14–A16 (sympy:76–88) / M7–M7-mut (math:110–122) | extra — nothing in paper card |

Dominant pattern: most of the script's verification scope is not anchored to the paper card. Two partial matches at the very end (Gaussian integrals). `paper_alignment: misaligned`.

## Assertion inventory

| #  | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|----|--------|------|------|------------------------------|--------------------|
| A1 | sympy  | 30–33 | `assert_zero (u4 - 4 u2²) - (D0(B4+Z4) - 3(M+B2+Z2)²)/D0²` | none (script family 1) | yes (to script claim) |
| A2 | sympy  | 37 | `assert_zero (u4 - 4 u2²).subs(K, K_from_one_pole)` | none (script family 2) | yes |
| A3 | sympy  | 38 | `assert_zero (N0/D0).subs(K, K_from_norm) - Ptarget` | none (script family 2) | yes (definitional) |
| A4 | sympy  | 39 | `assert_nonzero (u4 - 3 u2²).subs(K, K_from_one_pole)` | none | yes (mutation) |
| A5 | sympy  | 40 | `assert_nonzero (N0/D0).subs(K_from_norm) - 2 Ptarget` | none | yes (mutation) |
| A6 | sympy  | 41–45 | `assert_zero N0_from_compat - N0_from_equality` | none | yes (post-fix of v1 F1) |
| A7 | sympy  | 62 | `assert_zero coeff_matrix.det() - 1/27` | none (script family 3) | yes |
| A8 | sympy  | 67 | `assert_zero sol[dKSigma] - (B01+Z01+27(B41+Z41))` | none | yes |
| A9 | sympy  | 68 | `assert_zero sol[dMSigma] - (-(B21+Z21)+3(B41+Z41))` | none | yes |
| A10| sympy  | 69 | `assert_nonzero coeff_matrix.det() + 1/27` | none | yes (mutation) |
| A11| sympy  | 71–72 | `assert_zero D01_comp - 27(B41+Z41)` | none | yes (corollary of A8) |
| A12| sympy  | 73 | `assert_zero K1.subs({dK: expected_dK, dM: expected_dM})` | none | yes (post-fix of v1 F2) |
| A13| sympy  | 74 | `assert_zero H_even.subs(...)` | none | yes |
| A14| sympy  | 75 | `assert_nonzero K1.subs({dK: expected_dK+1, dM: expected_dM})` | none | yes (mutation) |
| A15| sympy  | 76–79 | `assert_zero Xi1.subs(sol) - (N01/N0 - 27(B41+Z41)/D0)` | none (script family 4) | yes |
| A16| sympy  | 80–83 | `assert_nonzero Xi1.subs(sol) - (N01/N0 + 27(B41+Z41)/D0)` | none | yes (mutation) |
| A17| sympy  | 84–88 | `assert_zero Xi1_from_expected - (N01/N0 - 27(B41+Z41)/D0)` | none | yes (independent path, post-fix of v1 F4) |
| A18| sympy  | 94 | `assert_zero MSigma_example - √π` | eq:stage018-msigma (partial, μ_η=1) | partial |
| A19| sympy  | 95 | `assert_zero KSigma_example - 3√π/2` | eq:stage018-ksigma (partial, T_w=K_η+6 T_Ω=1) | partial |
| M1 | math   | 38–41 | `FullSimplify[(u4 - 4 u2²) - onePoleNumerator] === 0` | none | yes (matches A1) |
| M2 | math   | 44–47 | `FullSimplify[(u4 - 4 u2²) /. K → KFromOnePole] === 0` | none | yes (matches A2) |
| M3 | math   | 50–53 | `FullSimplify[((N0/D0) /. K → KFromNorm) - Ptarget] === 0` | none | yes |
| M3-mut | math | 55–58 | `FullSimplify[... - 2 Ptarget] === 0` flagged FAIL | none | yes (mutation) |
| M4 | math   | 60–66 | `FullSimplify[N0FromCompatibility - N0FromEquality] === 0` | none | yes |
| M5 | math   | 84–89 | `FullSimplify[detGate - 1/27] === 0` | none | yes |
| M5-mut | math | 91–94 | mutation `+ 1/27` | none | yes |
| M6 | math   | 96–108 | `FullSimplify[K1 /. closedSlopeRules] === 0` and `HEven` | none | yes (uses substitution rules rather than `Solve` — independent path) |
| M7 | math   | 110–115 | `FullSimplify[(Xi1 /. closedSlopeRules) - expectedXi1] === 0` | none | yes |
| M7-mut | math | 117–122 | sign mutation | none | yes |
| M8a | math | 127, 129–132 | `FullSimplify[massIntegral - √π] === 0` | eq:stage018-msigma (partial) | partial |
| M8b | math | 128, 134–137 | `FullSimplify[stiffnessIntegral - 3√π/2] === 0` | eq:stage018-ksigma (partial) | partial |

The "Exercises which paper claim?" column shows the misalignment graphically: 14 of 19 SymPy assertions (and the matching Mathematica blocks) trace to no paper-stated deliverable.

## Findings

### F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Severity:** medium

**Files:**
- paper side: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_018.tex:26-28`
- script side: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:1-2, 25-88`
- script side: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl:1-15, 24-122`

**What's wrong:**

Paper card `\stagefield{Output}` (stage_018.tex:26-28):

> "Stage~018 exports the parent-wall bundle bridge \eqref{eq:stage018-msigma}--\eqref{eq:stage018-ksigma}."

Body equations (stage_018.tex:14-22): only `M_Σ` and `K_Σ` branch integrals.

Part-01 appendix row 58: "Bundle-level parent-action identities used by the projected electromagnetic response." (Plural "identities" hints at more, but the card never enumerates them.)

Script docstring (sympy:1-2): `"""Master-note audit for step_16_parent_throat_action_bundle_master_notes.md."""` — references a notes file that does not exist anywhere under `/var/projects/toy_physics/research/pde_ledger/`. The Mathematica header (lines 3–15) enumerates eight independent claim families M1–M8.

The script's load-bearing assertions (A1–A17 in SymPy; M1–M7-mut in Mathematica) verify:
1. an algebraic identity about coefficients of a quartic pole expansion (one-pole numerator),
2. two `KSigma` closure conditions and their compatibility (one-pole vs. normalization closures),
3. a 2×2 even-gate determinant equal to `1/27` and the closed-form wall-stiffness / wall-inertia slopes,
4. a residual `Xi1` after applying those slopes,

— none of which are stated in the stage-018 paper card or anchored in any existing notes file (notes/stages/moving_throat_pde_stage018_*.md is empty). Only the trailing Gaussian integrals A18/A19 / M8 directly correspond to `\stagefield{Output}`, and they do so only for a specific concrete profile (`μ_η = T_w = K_η + 6 T_Ω = 1`, `β = exp(-w²/2)`) — the paper's abstract integrals are not exercised symbolically.

Conversely, the paper's two declared deliverables (`M_Σ` and `K_Σ` integral identities, treated as exports rather than identities-to-prove) get only one concrete-profile sanity check each. The paper card never restates families 1–4 of the script, and no notes file backs them either.

**Why this matters:**

This is the v2 paper-alignment gate. The script certifies a body of algebra (one-pole closures, even-gate determinant, wall-slope solve, Xi1 residual) that the published stage card does not declare as the stage's output. Three resolution directions are possible and the user must choose:

- (a) the script's extra families are genuinely part of stage 018 and the paper card is too terse (the part-01 appendix's "identities used by the projected electromagnetic response" is too thin to anchor four extra claim families on its own); the paper card needs to be expanded and the missing notes file restored.
- (b) the script's extra families belong to a different stage (e.g., 015–017 or 021's "Reduced Maxwell/mixed one-port normal form"); the script should be trimmed to just the bridge identity and the extra checks moved.
- (c) the paper card's `\stagefield{Output}` is correct and complete, the script's extra checks are scaffolding from an earlier integration step that should be pruned.

Codex must NOT auto-resolve this; the orchestrator routes the choice to the user.

**Required change:**

See `## Resolve before fix_loop` block in the directive. Do not edit either side until the user picks (a), (b), or (c).

**Verification:**

After user resolution, a follow-up directive will name the specific paper-side or script-side edits and the verifier will re-run `redteam exec-sympy 018` and `redteam exec-mathematica 018` to confirm both sides match.

## Independent-derivation check (Mathematica)

A `.wl` is present. The two scripts share the same variable names and the same overall block structure, but the algebra is partially independent:

- For `u2, u4`: SymPy computes them directly from `D0, D2, D4` as `u2 = -D2/D0`, `u4 = (D2² - D0 D4)/D0²` (script lines 28–29). Mathematica derives them via a `Series[1/(D0 + D2 x² + D4 x⁴), {x, 0, 4}]` expansion and extracts coefficients (math lines 32–35). These are genuinely different paths to the same coefficient — good independence.
- For the wall-slope solve: SymPy calls `sp.solve([Eq(K1, 0), Eq(H_even, 0)], [dK, dM])` and compares the solution to the expected closed forms (sympy:64–68). Mathematica skips `Solve` and instead applies the closed-form `closedSlopeRules` directly to `K1` and `HEven`, asserting both vanish (math:96–108). Different path; both arrive at the same identity.
- For the gate Jacobian: SymPy assembles the 2×2 matrix manually using `sp.diff(K1, dKSigma)` etc. (sympy:58–61); Mathematica uses `D[gateVector, {{dKSigma, dMSigma}}]` (math:84). Different syntactic path.
- For `Xi1`, both scripts define `Xi1 = N01/N0 - D01/D0` and substitute the closed-form slopes, asserting equality with `N01/N0 - 27(B41+Z41)/(K - B0 - Z0)`. Structurally the same here.

Overall the Mathematica script is not a line-by-line transliteration; the series-expansion path for the one-pole numerator and the substitution-rule path for the slopes are independent re-derivations. **No `mathematica_transliteration` finding.**

## Engine cross-check

Both engines pass. SymPy output (`STATUS: PASS`) and Mathematica output (each `residual = 0` for substantive checks, `STAGE 018 MATHEMATICA AUDIT PASS`) agree:
- M1 residual 0 ↔ A1 passes
- M2 residual 0 ↔ A2 passes
- M3 residual 0 ↔ A3 passes
- M3-mut residual `-Ptarget` (non-zero) ↔ A5 nonzero
- M4 residual 0 ↔ A6 passes
- M5 residual 0 ↔ A7 passes
- M5-mut residual `2/27` ↔ A10 nonzero
- M6 residuals 0 ↔ A12/A13 pass
- M7 residual 0 ↔ A15 passes
- M7-mut residual `54(B41+Z41)/(B0-K+Z0)` = `-54(B41+Z41)/(K-B0-Z0)` ↔ A16 nonzero (sign matches)
- M8 inertia 0 ↔ A18 passes
- M8 stiffness 0 ↔ A19 passes

`engines_agree: true`. Both saved outputs are newer than their scripts (`outputs_fresh: true`). No `stale_output` finding.

## Verdict justification

The script's algebra holds up under attack: I checked the determinant by hand (`(1/9)(2/3) - (-1)(-1/27) = 2/27 - 1/27 = 1/27`), the one-pole numerator factor (`u4 - 4 u2² = D0 D4/D0² - 3 D2²/D0²` cancels the higher-order pieces correctly given the sign conventions), the Gaussian moments (`∫ exp(-w²) dw = √π`, `∫ w² exp(-w²) dw = √π/2`, so the stiffness integral is `√π/2 + √π = 3√π/2`), and the mutation-row residuals. The two engines agree at every checked claim and the SymPy script has been cleaned up since the v1 audit (the three tautologies F1/F2 and the corollary-only F4 were fixed; the missing Mathematica script F3 was added). However, what the script verifies is dominantly **not what the paper card says this stage proves**. The card's `\stagefield{Output}` is the two bridge integrals; the script verifies a one-pole closure, an even-gate determinant, two wall-slope formulas, and a residual amplitude — all of which appear to belong to an earlier integration-step note (`step_16_parent_throat_action_bundle_master_notes.md`) that no longer exists in the repository. The notes file referenced by the script docstring is gone, no replacement notes file under `notes/stages/moving_throat_pde_stage018_*.md` exists, and the paper card is silent on families 1–4. Verdict is `findings` with the single finding being `paper_misalignment` (subtype `paper_missing_script_claim`). The orchestrator must halt and let the user pick the resolution direction; Codex is not authorized to silently update either side.

## Self-test notes

(1) Variable independence: the SymPy `coeff_matrix` rows differentiate `K1 = D21 + D01/9` and `H_even = D41 - (2/3) D21 - D01/27` with respect to `dKSigma` and `dMSigma`. `K1` depends on both through `D01 = dK - B01 - Z01` and `D21 = -(dM + B21 + Z21)`, so the derivatives are non-trivial (1/9, -1) and (-1/27, 2/3); determinant is `2/27 - 1/27 = 1/27`, matching A7. No silent-zero derivative trap.
(2) Parity: A18 / M8a integrate `β² = exp(-w²)` (even) over `(-∞, ∞)` → nonzero `√π`; A19 / M8b integrate `(β')² + β² = w² exp(-w²) + exp(-w²)` (even in `w`) → nonzero `3√π/2`. No parity-induced false zero.
(3) Trivial-case pre-check: A4 `(u4 - 3 u2²).subs(K, K_one_pole)` reduces to `u2² + (D0(B4+Z4) - 3(M+B2+Z2)²)/D0² = u2² + 0 = u2²` (non-zero) — mutation correctly nonzero. A16 mutation gives `-54(B41+Z41)/(K-B0-Z0)` — non-zero. The Mathematica output `(54*(B41 + Z41))/(B0 - KSigma + Z0)` matches exactly (sign flip from denominator rearrangement). (4) Paper round-trip: the only finding is `paper_misalignment`, so no Codex-applied fix is prescribed; the directive's body is a `## Resolve before fix_loop` block routed to the user. No new misalignment introduced.
