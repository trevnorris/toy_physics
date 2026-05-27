---
unit_id: 098
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage098_family1_support_is_automatic.md]
  paper_appendix: present
---

# Audit unit 098 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 1230: `\input{stages/stage_098}`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.txt`

## What the paper claims

Per the notes (the .tex card is terse), the stage establishes that on the actual isotropic branch with `rho_alpha = 4/3`, the support-ratio demand is `zeta_req^(act)(eps_blk) = (1/3)/(1 - (2/3) eps_blk) = 1/(3 - 2 eps_blk)`. The stage then proves that **any** explicit support/source family whose support ceiling satisfies `zeta_max > 1` passes the support test throughout the admissible blocked window `0 <= eps_blk < 1/zeta_max` — i.e. the gap `zeta_max - zeta_req^(act)|_{eps=1/zmax}` factors as `3 zmax (zmax-1)/(3 zmax-2) > 0` for any `zmax > 1`. The Family-1 specialization uses the carried-forward value `zeta_max^(F1) ≈ 2.46752922945601 > 1`, yielding worst-case `zeta_req^(act) < 0.456730991107963 < 2.46752922945601`. The card's central blockquote is "With `rho_alpha = 4/3`, any branch with `zeta_max > 1` passes the support test; Family-1 does." The paper card also lists three checklist items (static-limit `c_pole=1/4`, l=0/l=2 orthogonality, hypothesis carry-forward) that are **not** exercised by either script — they appear to belong to upstream/parallel stages but are listed in the card's checklist.

## What the script claims to verify

Both scripts derive the closed-form `zeta_req(eps) = 1/(3 - 2 eps)` from the literal `(4/3 - 1)/(1 - eps*(2 - 4/3))`, compute its derivative `2/(3-2 eps)^2`, evaluate it at the admissible-window edge `eps = 1/zmax` to get `zeta_edge = zmax/(3 zmax - 2)`, and verify `gap = zmax - zeta_edge = 3 zmax (zmax-1)/(3 zmax - 2)`. The Mathematica script additionally asserts `gap > 0` symbolically under the assumption `zMax > 1`, and pins the Family-1 numeric outputs `zeta_edge ≈ 0.4567309911...` and `margin ≈ 2.0107982383...` against literal targets via `expectApprox`. The SymPy script's only Family-1 check is `assert gap_F1 > 0`, a weak positivity assertion.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `zeta_req^(act) = 1/(3 - 2 eps_blk)` from `rho_alpha = 4/3` | sympy L12; math L56 (`expectZero["zeta_req - 1/(3-2 eps)"]`) | match |
| `zeta_req` strictly increasing in `eps_blk` (notes §1 "monotone increasing") | sympy L14-15 (derivative formula); math L57 | match (implicit via positivity of `2/(3-2 eps)^2`) |
| gap factorization `3 zmax (zmax-1)/(3 zmax - 2)` | sympy L21; math L58 | match |
| gap > 0 for any `zmax > 1` (the central theorem) | math L59 (`expectTrue`); sympy L27 (numeric only) | partial — sympy lacks symbolic positivity |
| Family-1 specialization `zeta_max^(F1) ≈ 2.46752922945601` | sympy L24; math L61 — both **hardcoded literals**, not derived | partial — load-bearing constant has no paper-side derivation in this stage |
| `zeta_req^(act) < 0.456730991107963 < 2.46752922945601` (notes §3) | math L67-68 (`expectApprox` targets); sympy L27 (positivity only) | match (math), partial (sympy) |
| Checks list: static limit `c_pole=1/4` | (absent in scripts) | missing (likely belongs upstream) |
| Checks list: l=0/l=2 orthogonality | (absent in scripts) | missing (likely belongs upstream) |
| Checks list: hypothesis carry-forward | n/a — ledger statement, not a check | n/a |
| Stage number labeling | scripts/banners say "Stage 81 / STAGE 081"; card title says "Stage 115"; label is `stage:098`; card body says "Stage 098" | mismatch — see F1 |

`paper_alignment: partial` because: (a) stage numbering is inconsistent across the paper card title (`Stage~115`), label (`stage:098`), body (`Stage~098`), and script banners (`STAGE 081`); (b) the Family-1 carry-forward constant is hardcoded without an in-stage derivation or paper-side anchor; (c) the paper card's Checks list lists three items that the scripts do not exercise.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 12 | `simplify(zeta_req - 1/(3-2 eps)) == 0` | closed-form `zeta_req` | yes |
| A2 | sympy | 15 | `simplify(dz - 2/(3-2 eps)^2) == 0` | monotonicity (via derivative formula) | yes |
| A3 | sympy | 21 | `simplify(gap_factored - 3*zmax*(zmax-1)/(3*zmax-2)) == 0` | gap factorization | yes |
| A4 | sympy | 27 | `assert gap_F1 > 0` | Family-1 numeric margin | partial — numerically weak; implied by A3 + zmax_F1>1 |
| M1 | math | 56 | `expectZero["zeta_req - 1/(3-2 eps)"]` | closed-form `zeta_req` | yes |
| M2 | math | 57 | `expectZero["d zeta_req / d eps exact formula"]` | monotonicity | yes |
| M3 | math | 58 | `expectZero["gap factorization"]` | gap factorization | yes |
| M4 | math | 59 | `expectTrue["automatic-support gap is positive for zmax > 1", gap > 0]` | central theorem (gap>0 ∀ zmax>1) | yes |
| M5 | math | 67 | `expectApprox["Family-1 zeta_edge numeric check", ..., 0.456730991107963169017835980412, 1e-15]` | Family-1 numeric bound from notes §3 | yes (target matches notes literal) |
| M6 | math | 68 | `expectApprox["Family-1 margin numeric check", ..., 2.01079823834804688464927835412, 1e-15]` | Family-1 numeric margin | yes |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script (+ paper card internal inconsistency)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex:1` quote: `\section[Stage~115]{Stage~115: The Explicit Family-1 Support Test Is Automatic on the Actual Isotropic Branch}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex:2` quote: `\label{stage:098}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_098.tex:7` quote: `Stage~098 is a geometry-lane firewall ledger step.`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md:1` quote: `# Moving-Throat PDE — Stage 098: The Explicit Family-1 Support Test Is Automatic on the Actual Isotropic Branch`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md:105` quote: `After Stages 80–81, the reduced theorem split is now fully sharp:`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:3` quote: `Stage 81 SymPy audit: actual isotropic branch support demand is automatic ...`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:35` quote: `print('\nSTAGE 81 AUDIT PASSED')`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:38` quote: `banner["STAGE 081 — FAMILY-1 SUPPORT IS AUTOMATIC"];`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:71` quote: `Print["Stage 098 Mathematica audit passed."];`

**What's wrong:**
The unit's stage number is reported four different ways across paper and scripts: the paper card's `\section[...]{...}` declares **Stage 115**; the `\label` and body say **Stage 098**; the SymPy docstring/print and Mathematica banner say **Stage 81 / 081**; the Mathematica final print says **Stage 098**. The notes section 4 references "Stages 80–81" suggesting the unit was originally numbered 81 before a renumbering to 098, and the paper card title was later mis-keyed to 115 (possibly a different stage's number was pasted). The math content is the same across all four; only the labeling diverges. This is a paper-side / docstring-side disagreement that has no script-math consequence but obstructs cross-referencing.

**Why this matters:**
A reader following the appendix to `Stage 098` lands on a section titled `Stage 115` whose scripts identify themselves as `Stage 81`. Provenance and downstream citation become ambiguous; any downstream stage that says "see Stage 81 / 098 / 115" can no longer be parsed unambiguously.

**Required change:**
User must pick the canonical stage number for this unit and align the four locations. Codex must not auto-resolve this — see `## Resolve before fix_loop` in the directive.

**Verification:**
After user picks the canonical number `N`, all four locations should read `Stage N`/`STAGE N`/`stage:N`, and notes section 4's reference to "Stages 80–81" should be updated if N differs from 81.

### F2 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:24` quote: `zmax_F1 = sp.N('2.46752922945601')`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:61` quote: `zMaxF1 = SetPrecision[2.46752922945601, 20];`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:67` quote: `expectApprox["Family-1 zeta_edge numeric check", zetaEdgeF1, 0.456730991107963169017835980412, 10^-15];`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:68` quote: `expectApprox["Family-1 margin numeric check", gapF1, 2.01079823834804688464927835412, 10^-15];`

**What's wrong:**
Both scripts hardcode `zeta_max^(F1) = 2.46752922945601` as a literal float. The paper card (`stage_098.tex`) does **not** mention this value at all; only the notes file does, and the notes too treat it as a carry-forward ("For the explicit Family-1 branch, `zeta_max^(F1) ≈ 2.46752922945601 > 1`") with no in-stage derivation. The two Mathematica `expectApprox` targets `0.456730991107963169017835980412` and `2.01079823834804688464927835412` are themselves arithmetic consequences of `zmax_F1 = 2.46752922945601` — substituting `2.46752922945601` into `zmax/(3 zmax - 2)` and `zmax - zmax/(3 zmax - 2)` reproduces those literals to the indicated precision. So the `expectApprox` checks against literal targets reduce to "this float, fed through this formula, equals this float," which exercises only the formula's numerical stability, not any independent fact about Family-1.

**Why this matters:**
If the carry-forward `zeta_max^(F1)` value is wrong (typo, off-by-one in a precursor stage), the script would still pass — the targets would be wrong too, in a self-consistent way. The audit's binding claim "Family-1 has `zeta_max > 1`" is robust to small errors in this value, but the specific number `2.46752922945601` is not cross-checked here against the upstream stage that derived it (the notes do not cite a stage number).

**Required change:**
Add a brief comment at the assignment site naming the source stage (or naming the precursor result that fixes `zeta_max^(F1)`). If no upstream verified result exists, the comment should so state. Replace the SymPy line 24 `sp.N('2.46752922945601')` with `sp.Rational`-style or a comment indicating the precursor source. Mathematica line 61 likewise.

This finding is **informational unless the carry-forward source cannot be located**. If a precursor stage's scripts verify `zeta_max^(F1) = 2.46752922945601`, this is acceptable as a quoted forward. Mark the source.

**Verification:**
Each assignment site has an inline comment naming the upstream stage or noting "external input, not verified here". The expectApprox targets remain unchanged.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:27` quote: `assert gap_F1 > 0`

**What's wrong:**
The SymPy script's Family-1 specialization check is just `assert gap_F1 > 0`. By that point in the script, line 21 has already proved symbolically that `gap = 3 zmax (zmax-1)/(3 zmax - 2)`, and `zmax_F1 = 2.46... > 1` makes the F1 margin trivially positive. So this assertion adds nothing — it cannot distinguish the script working from `zmax_F1` being any value > 1. The Mathematica script does better here: lines 67–68 pin both `zetaEdgeF1` and `gapF1` to specific literals with `1e-15` tolerance, which would catch a typo in `zMaxF1`. The SymPy side lacks this numeric pin.

**Why this matters:**
A typo in `sp.N('2.46752922945601')` (e.g., transposed digits) would not be detected by the SymPy script — the resulting `gap_F1` would still be positive. The notes section 3 makes specific numeric claims: `zeta_req^(act) < 0.456730991107963 < 2.46752922945601`. SymPy verifies neither bound.

**Required change:**
In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py`, after line 26, add numeric assertions matching the Mathematica side:

```python
assert abs(zeta_edge_F1 - sp.Float('0.456730991107963169017835980412', 30)) < sp.Float('1e-15')
assert abs(gap_F1 - sp.Float('2.01079823834804688464927835412', 30)) < sp.Float('1e-15')
```

Keep the existing `assert gap_F1 > 0` line, or replace it with the two pinned checks.

**Verification:**
The sympy script's output should still print the same `Family-1 zeta_edge = 0.456...` and `Family-1 margin = 2.010...` values, but a typo in `zmax_F1` would now fail at the new assertions. The script exits 0 on a fresh run.

## Independent-derivation check (Mathematica)

The Mathematica script's structure mirrors the SymPy script closely: both define `(4/3 - 1)/(1 - eps*(2 - 4/3))`, simplify, take a derivative, substitute `eps -> 1/zmax`, factor, then verify factorization. The variable names differ (`eps`/`zmax` vs `epsBlk`/`zMax`) but the algebra is identical. For an identity this short, true independent derivation is not really possible — the formula `(rho-1)/(1 - eps*(2-rho))` has a unique closed form. The Mathematica side **does** add two new checks not in SymPy: `expectTrue["... gap > 0", gap > 0]` (line 59) under explicit assumptions, and the `expectApprox` pinning of numeric targets (lines 67-68). I judge this **not** a transliteration finding — it's the natural minimum for verifying a 3-line algebraic identity, and the Mathematica side genuinely adds checks the SymPy side lacks (which is the substance of F3).

## Engine cross-check

| Quantity | SymPy output | Mathematica output | Agree? |
|---|---|---|---|
| `zeta_req(eps)` | `-1/(2*eps - 3)` | `(3 - 2*epsBlk)^(-1)` | yes (algebraically equal) |
| `d zeta_req / d eps` | `2/(2*eps - 3)**2` | `2/(3 - 2*epsBlk)^2` | yes (squared, sign-independent) |
| `zeta_edge` | `zmax/(3*zmax - 2)` | `zMax/(-2 + 3*zMax)` | yes |
| `gap` | `3*zmax*(zmax - 1)/(3*zmax - 2)` | `(3*(-1 + zMax)*zMax)/(-2 + 3*zMax)` | yes |
| Family-1 zeta_edge | `0.456730991107963169017835980412` | `0.4567309911079631649053851860788...` | agree to 16 digits; sympy prints exactly 30 digits truncated, math gives 19-digit precision string — agree within precision |
| Family-1 margin | `2.01079823834804688464927835412` | `2.0107982383480468887617291484513...` | agree to 16 digits |

Engines agree.

## Verdict justification

The math content of both scripts is correct and the central claim — `gap = 3 zmax (zmax-1)/(3 zmax-2) > 0 for zmax > 1` — is faithfully exercised by both engines. Three issues remain: (F1) the paper card has an internal stage-number inconsistency (`Stage~115` in section title vs `stage:098` label vs `Stage~098` in body) and the scripts add a fourth label `Stage 81 / 081`; (F2) the Family-1 numeric `zeta_max^(F1) = 2.46752922945601` is hardcoded with no paper-card anchor and no in-stage derivation; (F3) the SymPy script's Family-1 check is too weak to catch a typo in the hardcoded constant. F1 is a paper-side disagreement that requires user resolution (which number is canonical?). F2 and F3 are script-side and can be applied by Codex. None propagate downstream in a math-breaking way — only label propagation, which is a documentation concern. Verdict: `findings`, `stop_cold: null`.

Attacks tried that did not produce findings:
- Tried to find a sign error in the derivative — none; `2/(3-2 eps)^2` is positive for `eps < 3/2`, which is implied by `eps < 1/zmax < 1` for `zmax > 1`.
- Tried to find a domain pathology at `zmax = 2/3` (denominator zero) — excluded by `zmax > 1` assumption.
- Tried to find an `eps_blk` range hole — the admissible window `0 <= eps_blk < 1/zmax` is honored by Mathematica's `$Assumptions` (line 43).
- Tried to find a tautology in `expectZero["gap factorization", gap - gapExpected]` (line 58) — `gap` is computed via `FullSimplify[zMax - zetaEdge]` and `gapExpected` is the claimed factored form; the difference being zero is a substantive identity, not a definition.
- Tried to find a hidden assumption that makes `expectTrue["gap > 0", gap > 0]` trivially true — under `zMax > 1`, the factored form `3 zMax (zMax-1)/(3 zMax - 2)` has positive numerator (`zMax-1 > 0`) and positive denominator (`3 zMax - 2 > 1 > 0`), so the result is genuinely a consequence of the assumption, not a tautology.

## Self-test notes

Checked: (1) no `sp.diff(EXPR, VAR)` where VAR is missing from EXPR — the only diff is `sp.diff(zeta_req, eps)` and `zeta_req` does depend on `eps`. (2) No integrals, so no parity trap. (3) Trivial-case for F3's proposed assertions: substituting `zeta_edge_F1 = 0.456730991107963169017835980412` and the literal target with diff < 1e-15 reduces to `|0| < 1e-15` — passes. (4) Path specifications: F3 edits `scripts/...sympy_audit.py`, correct directory. (5) Paper round-trip: F3 introduces no new paper claim — the proposed sympy literals are copied from the existing Mathematica `expectApprox` targets, which are already in the script and align with notes §3's `0.456730991107963` quote.
