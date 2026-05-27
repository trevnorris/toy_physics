---
unit_id: 078
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md
  paper_appendix: present
---

# Audit unit 078 red-team report (v2)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_078.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 274 includes this stage in Part III)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.txt`

## What the paper claims

Stage 078 decides whether the explicit Family-1 wall-depth datum is the bottleneck. Two boxed numeric verdicts are stated: for the natural shell-weighted (chi^2) datum, `Pe_suff^(chi) ≈ 96.5285247264386 lambda_mu^2` and `Pe_fail^(chi) ≈ 11220.5441626259 lambda_mu^2`; for the conservative Jensen floor, `Pe_suff^(J) ≈ 22.0062226330754 lambda_mu^2` and `Pe_fail^(J) ≈ 2558.01892349205 lambda_mu^2`. The `\stagefield{Output}` line states: "The first explicit Family-1 support/source verdict: wall-depth is not the leading bottleneck for moderate demand." The notes derive these as ratios `Theta_w^(branch) / Theta_X_coeff` using the Stage-75 threshold-window coefficients and the Stage-77 Theta_w extraction.

## What the script claims to verify

The SymPy script computes the four `Pe / lambda_mu^2` ratios from the Stage-75 threshold-window coefficients (`Theta_fail/Pe_req ≈ 3.62606e-4`, `Theta_suff/Pe_req ≈ 4.21495e-2`) and the Stage-77 Theta extractions (`Theta_chi/lambda_mu^2 ≈ 4.06863`, `Theta_J/lambda_mu^2 ≈ 0.927552`), prints them, asserts the in-branch ordering `Pe_suff < Pe_fail` for both branches, and asserts the cross-branch nesting `Pe_suff_J < Pe_suff_chi`, `Pe_fail_J < Pe_fail_chi`, and the window overlap `Pe_suff_chi < Pe_fail_J`. The Mathematica script independently re-derives `Theta_fail` from the symbolic `Sinh/Cosh` closed form (Stage-75 line 20 of its output), bootstraps `Theta_suff` from `Theta_fail` times the decimal ratio `(4.21495e-2 / 3.62606e-4)`, adopts the Stage-77 chi/J Theta values at extended precision, computes the same four ratios, and runs `expectApprox` and `expectTrue` mirrors of the SymPy checks.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Pe_suff^(chi) ≈ 96.5285247264386 lambda_mu^2` | SymPy line 46 print, Mathematica `expectApprox` line 76 | match (sympy 96.528524726438575954; mathematica 96.52852472643852; agree to 1e-13) |
| `Pe_fail^(chi) ≈ 11220.5441626259 lambda_mu^2` | SymPy line 47 print, Mathematica `expectApprox` line 77 | match |
| `Pe_suff^(J) ≈ 22.0062226330754 lambda_mu^2` | SymPy line 48 print, Mathematica `expectApprox` line 78 | match |
| `Pe_fail^(J) ≈ 2558.01892349205 lambda_mu^2` | SymPy line 49 print, Mathematica `expectApprox` line 79 | match |
| Verdict: "wall-depth not leading bottleneck" | SymPy `Pe_suff_chi < Pe_fail_J` (window overlap), Mathematica `expectTrue` line 86 | match (the chi-branch success threshold lies below even the Jensen-floor failure ceiling, so any modest demand below 22.0 lies in the joint admissible band) |

Paper alignment: aligned. All four boxed values are reproduced to 13+ digits in both engines, and the qualitative verdict is anchored by a non-tautological window-overlap inequality.

## Assertion inventory

| #  | Script      | Line | Form                                                              | Exercises which paper claim?                            | Anchored to claim? |
|----|-------------|------|-------------------------------------------------------------------|---------------------------------------------------------|--------------------|
| A1 | sympy       | 51-52 | `Pe_suff_chi < Pe_fail_chi` (raises AssertionError)              | structural sanity (Theta_chi cancels)                  | partial            |
| A2 | sympy       | 53-54 | `Pe_suff_J < Pe_fail_J`                                           | structural sanity (Theta_J cancels)                    | partial            |
| A3 | sympy       | 60-64 | `Pe_suff_J < Pe_suff_chi` (verdict, depends on Theta_J/Theta_chi)| chi-vs-J branch nesting                                 | yes                |
| A4 | sympy       | 65-69 | `Pe_fail_J < Pe_fail_chi`                                         | chi-vs-J branch nesting                                 | yes                |
| A5 | sympy       | 70-74 | `Pe_suff_chi < Pe_fail_J` (window overlap = verdict)             | "wall-depth not leading bottleneck"                    | yes (depends on all 4 constants) |
| A6 | mathematica | 76-79 | `expectApprox` of the four Pe ratios against independently computed targets | the four numeric Pe values                  | partial (Theta_fail genuinely independent; Theta_suff and Theta_chi/J adopted) |
| A7 | mathematica | 80-81 | `expectTrue["Pe_suff^(chi) < Pe_fail^(chi)", ...]` etc.          | structural sanity (cancellation again)                 | partial            |
| A8 | mathematica | 84-86 | `expectTrue["Pe_suff^(J) < Pe_suff^(chi)", ...]` etc.            | chi-vs-J branch nesting + window overlap               | yes                |

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:44` (`thetaSuffSym = thetaFailSym * (4.21495341569977*^-2 / 3.62605617972939*^-4)`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:48-49` (`thetaChiCoeffNum`, `thetaJCoeffNum` as typed decimals)

**What's wrong:**
The v1 directive (F3) instructed Codex to derive `Theta_suff` symbolically; Codex instead bootstrapped it from `thetaFailSym` times a literal decimal ratio. The actual Stage-75 closed form is

```
Theta_suff = -(45*Pe_req*cosh(111*sqrt(5)/5) + 27*sqrt(5)*Pe_req*sinh(111*sqrt(5)/5))
              / (2500 - 2500*cosh(111*sqrt(5)/5))
```

(Stage-75 sympy output line 21), but the Mathematica script does not use it. Consequently if the literal decimal ratio `(4.21495341569977e-2 / 3.62605617972939e-4)` were corrupted, the Mathematica `expectApprox` target for `Pe_suff^*` would shift in lockstep with the literal used by SymPy and the engine cross-check would still pass.

Similarly, `thetaChiCoeffNum` and `thetaJCoeffNum` are typed-in 40-digit decimal strings; the inline comment at lines 45-47 says "we adopt them at high precision but verify their ratio chi:J matches the Stage-77 ratio" — but no such ratio check exists anywhere in the script.

**Why this matters:**
The load-bearing `Theta_fail` IS independently derived (correctly, from the symbolic `Sinh/Cosh` form). The remaining three coefficients are not, so the second engine's independence is partial. Given that all four Stage-75/77 values come from upstream stages whose audits are themselves in this red-team campaign, the cross-check is weakened but not broken. This is a residual concern, not a blocker.

**Required change:**
Replace line 44 with a direct symbolic re-derivation of `Theta_suff`:

```mathematica
(* Independent closed form for Theta_suff from Stage-75 (sympy output line 21):
   Theta_suff/Pe_req = -(45 cosh(alpha) + 27 sqrt(5) sinh(alpha))
                        / (2500 - 2500 cosh(alpha)),  alpha = 111 sqrt(5)/5 *)
thetaSuffSym = -(45 Cosh[111 Sqrt[5]/5] + 27 Sqrt[5] Sinh[111 Sqrt[5]/5])
                / (2500 - 2500 Cosh[111 Sqrt[5]/5]);
```

Either remove the misleading comment at lines 45-47 (about a ratio check that does not exist) OR add the ratio check it promises:

```mathematica
expectApprox["Theta_chi/Theta_J ratio",
  thetaChiCoeffNum / thetaJCoeffNum,
  ToExpression["4.387185...`30"],  (* the Stage-77 ratio, computed independently *)
  10^-10];
```

(If no upstream stage gives this ratio numerically, simply delete the lines 45-47 comment claim.)

**Verification:**
After fix, line 44 contains `Cosh[111 Sqrt[5]/5]` (not just a decimal ratio), and either the comment at lines 45-47 is removed or a corresponding `expectApprox` exists below the chi/J definitions. The four `... numeric check diff = ...` output lines remain PASS.

### F2 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:3` (`"""Stage 61 SymPy audit."""`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py:23` (`banner("STAGE 61 — FAMILY-1 BRANCH VERDICT")`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl:32` (`banner["STAGE 061 — FAMILY-1 BRANCH VERDICT"]`)
- notes file mentions `Stage 58`, `Stage 60`, `Stage 61` throughout when the paper-side numbering is `075`, `077`, `078`.

**What's wrong:**
The SymPy docstring/banner and Mathematica banner identify this unit as "Stage 61" (or "061"), while the paper card, filename, and pipeline manifest say `stage 078`. The notes file similarly uses the older numbering (`Stage 58`, `60`, `61`). This is purely a labeling artifact from a previous numbering scheme, but it means an unsuspecting reader of the script output sees `STAGE 61` while the paper card cites `Stage 078`, which is a documentation hazard.

**Why this matters:**
No mathematical content is affected. The risk is cosmetic: future readers (including automated audits) may confuse "Stage 61" output with a non-existent stage. The notes-side stale numbering is the user's call — notes/ is not editable by Codex per the directive policy.

**Required change (script side only, paper_misalignment fix scope is the user's call):**

(a) In `scripts/moving_throat_pde_stage078_family1_branch_verdict_sympy_audit.py`:
- Line 3: replace `Stage 61 SymPy audit.` with `Stage 078 SymPy audit.`
- Line 23: replace `banner("STAGE 61 — FAMILY-1 BRANCH VERDICT")` with `banner("STAGE 078 — FAMILY-1 BRANCH VERDICT")`

(b) In `mathematica/moving_throat_pde_stage078_family1_branch_verdict_mathematica_audit.wl`:
- Line 32: replace `banner["STAGE 061 — FAMILY-1 BRANCH VERDICT"]` with `banner["STAGE 078 — FAMILY-1 BRANCH VERDICT"]`

The notes file labeling is not the script's responsibility and is not modified.

## Resolve before fix_loop

The notes file at `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage078_family1_branch_verdict.md` uses the legacy numbering "Stage 58 / Stage 60 / Stage 61" while the paper uses "Stage 075 / 077 / 078". This is a notes-side editorial question the user should resolve (renumber notes to match paper, or leave as historical record). Codex must not edit notes/.

**Verification:**
After Codex applies the script-side renames, the script output text shows `STAGE 078 — FAMILY-1 BRANCH VERDICT` instead of `STAGE 61` / `STAGE 061`. The mathematica output's banner matches. The notes-side question is logged for user resolution.

## Independent-derivation check (Mathematica)

The Mathematica script now derives `thetaFailSym` from the explicit `Sinh/Cosh` closed form (lines 39-43), which is the load-bearing independence the v1 directive (F3) required. However, `thetaSuffSym` (line 44) is bootstrapped from `thetaFailSym` times a decimal literal ratio rather than its own symbolic closed form, and `thetaChiCoeffNum` / `thetaJCoeffNum` (lines 48-49) are typed-in decimals rather than derived. The `expectApprox` targets at lines 76-79 are now computed by Mathematica from `thetaFailSym` and the typed-in chi/J numerics, not literal-copied from SymPy output. This is a genuine improvement over v1 (where every target was a SymPy-output decimal). The remaining partial-dependency is a low-severity residual — see F1.

## Engine cross-check

| Quantity                       | SymPy                          | Mathematica                   |
|--------------------------------|--------------------------------|-------------------------------|
| `Pe_suff^(chi) / lambda_mu^2` | `96.528524726438575954`        | `96.52852472643852`           |
| `Pe_fail^(chi) / lambda_mu^2` | `11220.544162625905301`        | `11220.54416262589764...`     |
| `Pe_suff^(J)   / lambda_mu^2` | `22.006222633075413597`        | `22.006222633075424`          |
| `Pe_fail^(J)   / lambda_mu^2` | `2558.0189234920526360`        | `2558.0189234920537145...`    |

Agreement is to ~13 digits in all four. The Mathematica `expectApprox` deltas report `0` or `0``25.6...` (i.e. zero at precision ~25), so the cross-engine check is healthy. No `engine_disagreement` finding.

The Mathematica engine also independently verifies the three branch-verdict inequalities `Pe_suff^(J) < Pe_suff^(chi)`, `Pe_fail^(J) < Pe_fail^(chi)`, and `Pe_suff^(chi) < Pe_fail^(J)` via `expectTrue`, mirroring the SymPy `if not (...)` checks. Both engines emit PASS for all three. The window-overlap inequality (which depends on all four upstream constants) is the load-bearing assertion for the paper verdict and is honored by both engines.

## Verdict justification

Verdict: `findings` (2, both low-severity). The v1 directive's core requirements have been satisfied: the branch-verdict assertions now exist and depend on all four upstream constants (so a corruption in `Theta_chi`, `Theta_J`, `Theta_fail_coeff`, or `Theta_suff_coeff` would now fail the script), the Mathematica engine derives the load-bearing `Theta_fail` from its symbolic closed form, and the four boxed paper values are reproduced to ~13 digits by both engines. The paper card's verdict ("wall-depth not leading bottleneck") is anchored by the window-overlap assertion `Pe_suff_chi < Pe_fail_J`. Remaining issues: `Theta_suff` is decimal-ratio-bootstrapped rather than independently re-derived (F1), and the script banners are mis-labeled "Stage 61"/"061" instead of "Stage 078" (F2). No stop-cold: no math error, no downstream propagation, no paper↔script value mismatch.

Outputs are fresh (script .py at epoch 1779514083, .txt at 1779514311; .wl at 1779514238, .txt at 1779514318 — outputs newer than scripts in both cases).

## Self-test notes

I verified that the new branch-verdict assertion `Pe_suff_chi < Pe_fail_J` depends on all four upstream constants by reducing to `Theta_chi × Theta_fail_coeff < Theta_J × Theta_suff_coeff` (no cancellation), confirming F4 of v1 has the intended teeth. I also computed `4.21495341569977e-2 / 3.62605617972939e-4 ≈ 116.24` and `3.62605617972939e-4 × 116.24 ≈ 4.215e-2`, confirming the Mathematica `thetaSuffSym` bootstrap reproduces the Stage-75 Theta_suff coefficient — the math is right; the issue is only that the derivation is not independent. Finally, I matched the paper's four boxed values byte-for-byte against the SymPy output digits (96.5285247264386 vs 96.528524726438575954; 11220.5441626259 vs 11220.544162625905301; 22.0062226330754 vs 22.006222633075413597; 2558.01892349205 vs 2558.0189234920526360) — agreement to all stated paper digits.
