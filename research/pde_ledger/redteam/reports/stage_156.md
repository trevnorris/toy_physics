---
unit_id: 156
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage156_renormalized_canonical_branch.md]
  paper_appendix: present
---

# Audit unit 156 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_156.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (read the co-evolving / renormalized-canonical block at lines 900–971 plus the `\input{stages/stage_156}` insertion at line 1346)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage156_renormalized_canonical_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage156_renormalized_canonical_branch_mathematica_audit.txt`

## What the paper claims

The card (`stage_156.tex`) says the unit verifies that "Exact compensation is restored at \(\Sigma_0^{\rm can}\approx4.65103\), \(\widehat T_{m,{\rm can}}\approx1.44671\)" on the full co-evolving core–mouth branch. The notes give the full deliverable set: a unique positive root of \(\mathfrak g_{\rm fp}(\Sigma_0)=\mathfrak g_*\) on the analyzed interval at \(\Sigma_0^{\rm can}\approx 4.651033550168867\), with \(\widehat T_{m,\rm can}=\sqrt{(9/20)\Sigma_0^{\rm can}}\approx 1.4467083664567613\), and the restored canonical point \(\mathfrak g_{\rm can}=\mathfrak g_*\approx 0.758035078944663\), \(\mathcal R_{\rm can}=\tfrac14\), \(\mathcal S_{\rm can}\approx 0.6703621156734617\), \(\Pi_{\rm can}\approx 3.871564377479002\). The appendix block (lines 900–971) restates the same six numbers (to within a final-digit roundoff) and frames the result as a "StatusNumerical" placement inside the exact co-evolving map. The card's Checks list adds three qualitative properties: deviations taken about the renormalized point, even-preservation imposed before reading the odd defect, and tangent motion on the parent family giving \(\delta_\perp=0\) — but these qualitative checks are downstream-flavoured (Stages 157+ territory) and the load-bearing quantitative deliverable here is the six-number set above.

## What the script claims to verify

The SymPy script numerically solves the discrete fixed-point map \(\Sigma\mapsto e^{-\Phi_{\Sigma_0}[\Sigma]}/\int e^{-\Phi}\) on an \(N=2401\) trapezoidal grid via Picard iteration, then runs a bisection on \([3,6]\) to locate \(\Sigma_0^{\rm can}\) with \(g_{\rm fp}(\Sigma_0^{\rm can})=g_*\), and prints all six numerical outputs (\(\Sigma_0^{\rm can}\), \(g_{\rm can}\), \(\mathcal S_{\rm can}\), \(\mathcal R_{\rm can}\), \(\Pi_{\rm can}\), \(\widehat T_{m,\rm can}\)). However, its only `assert` statements (lines 123–126) check `|g_can - g_star| <= 1e-10` and `|R_can - 0.25| <= 1e-10`. The Mathematica script does the same numerical solve (with a small overflow-guard `phShift` and memoization), additionally does a monotonicity scan of \(\mathfrak g_{\rm fp}\) on \([3,6]\) at 0.5 steps via `expectTrue`, and then `expectApprox`-asserts all six numerical outputs against the notes/appendix values (lines 147–152).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check (SymPy) | Script-side check (Mathematica) | Status |
|---|---|---|---|
| \(\mathfrak g_{\rm fp}(\Sigma_0)\) monotone on the analyzed interval (uniqueness premise, notes §1) | — (not exercised) | `expectTrue["monotone increase on scan grid", Min[scanDiffs] > 0]` at line 124 | partial (Mathematica only) |
| \(\Sigma_0^{\rm can}\approx 4.651033550168867\) | printed at line 111 but not asserted | `expectApprox` line 147, tol \(10^{-11}\) | partial (asserted only in Mathematica) |
| \(\widehat T_{m,\rm can}\approx 1.4467083664567613\) (= \(\sqrt{9\Sigma_0/20}\)) | printed at line 116, only `T_hat_can = math.sqrt(9 * Sigma0_can / 20)` (line 109), no assert | `expectApprox` line 152, tol \(10^{-11}\) | partial (asserted only in Mathematica) |
| \(\mathfrak g_{\rm can}=\mathfrak g_*\approx 0.758035078944663\) | `assert abs(g_can - g_star) > 1e-10` (line 123) — note: this is what the bisection drives to, so it is essentially the convergence tolerance | `expectApprox` line 148, tol \(10^{-12}\) | match (but both engines' g-check is the bisection target, not an independent verification) |
| \(\mathcal R_{\rm can}=\tfrac14\) | `assert abs(R_can - 0.25) > 1e-10` (line 125) — algebraically forced by `g_can = g_*` and the input constant `rF1` since \(R=(g_*-r_{F1})^2/(1+r_{F1}^2)\) | `expectApprox` line 150, tol \(10^{-10}\) | match-but-weak (algebraically guaranteed once g is at \(g_*\); checks the input constants more than a derived result) |
| \(\mathcal S_{\rm can}\approx 0.6703621156734617\) | printed at line 113 but not asserted | `expectApprox` line 149, tol \(10^{-12}\) | partial (asserted only in Mathematica) |
| \(\Pi_{\rm can}\approx 3.871564377479002\) | printed at line 115, computed as `Sigma0_can*(1 - R_can*S_can)` (line 108), not asserted | `expectApprox` line 151, tol \(10^{-11}\) | partial (asserted only in Mathematica) |

Set `paper_alignment: partial`. The two engines together cover the full paper claim, but the SymPy side is materially weaker than the Mathematica side: it prints the load-bearing numerical outputs of the stage and then `assert`s only the bisection convergence target and an algebraically-forced consequence.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 123 | `if abs(g_can - g_star) > 1e-10: raise` | \(\mathfrak g_{\rm can}=\mathfrak g_*\) (also bisection convergence tolerance) | partial (re-checks the bisection's own target) |
| A2 | sympy | 125 | `if abs(R_can - 0.25) > 1e-10: raise` | \(\mathcal R_{\rm can}=\tfrac14\) | partial (algebraically forced by A1 plus the input constant `rF1`) |
| A3 | mathematica | 124 | `expectTrue["monotone increase on scan grid", Min[scanDiffs] > 0]` | uniqueness premise (monotonicity of \(g_{\rm fp}\) on \([3,6]\)) | yes |
| A4 | mathematica | 147 | `expectApprox["Sigma0_can numeric check", sigma0Can, 4.651033550168867, 10^-11]` | \(\Sigma_0^{\rm can}\) | yes (but compares against the same number quoted in the notes; load-bearing all the same) |
| A5 | mathematica | 148 | `expectApprox["g_can numeric check", gCan, 0.758035078944663, 10^-12]` | \(\mathfrak g_{\rm can}=\mathfrak g_*\) | partial (bisection target) |
| A6 | mathematica | 149 | `expectApprox["S_can numeric check", sCan, 0.6703621156734617, 10^-12]` | \(\mathcal S_{\rm can}\) | yes |
| A7 | mathematica | 150 | `expectApprox["R_can numeric check", rCan, 0.25, 10^-10]` | \(\mathcal R_{\rm can}=1/4\) | partial (algebraically forced once A5 holds) |
| A8 | mathematica | 151 | `expectApprox["Pi_can numeric check", piCan, 3.871564377479002, 10^-11]` | \(\Pi_{\rm can}\) | yes |
| A9 | mathematica | 152 | `expectApprox["T_hat_can numeric check", tHatCan, 1.4467083664567613, 10^-11]` | \(\widehat T_{m,\rm can}\) | yes |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py:111-126`

**What's wrong:**
The paper card states the bottom-line numerical deliverables as
> "Exact compensation is restored at \(\Sigma_0^{\rm can}\approx4.65103\), \(\widehat T_{m,{\rm can}}\approx1.44671\)."

and the notes/appendix add four more load-bearing numbers (\(\mathfrak g_{\rm can}\), \(\mathcal R_{\rm can}\), \(\mathcal S_{\rm can}\), \(\Pi_{\rm can}\)). The SymPy script computes and prints all six numbers (lines 103–116) but its `assert` block (lines 123–126) only checks
```python
if abs(g_can - g_star) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore g = g_*.")
if abs(R_can - 0.25) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore R = 1/4.")
```
The first assert is the bisection's own convergence target (`g_fp(Sigma0_can) = g_star` is what the bisection drives to high precision); it therefore measures bisection convergence rather than physics. The second assert is algebraically forced once the first passes: with `rF1 = 1.77799353547498` and `g_star = 0.758035078944663`, `R = (g_* - r_F1)^2 / (1 + r_F1^2)` is numerically \(\approx 0.249999\ldots\) regardless of the rest of the calculation; it is a check on the two input constants, not a derived consequence of the co-evolving solve. As a result, none of the four load-bearing numbers in the paper card and notes (\(\Sigma_0^{\rm can}\), \(\widehat T_{m,\rm can}\), \(\mathcal S_{\rm can}\), \(\Pi_{\rm can}\)) are guarded by an `assert` on the SymPy side. The Mathematica counterpart asserts all six via `expectApprox` (lines 147–152); the asymmetry is unjustified.

**Why this matters:**
If the discretization, the Picard iteration, or the bisection silently drifted (wrong N, wrong weights, wrong kernel) such that the solve converged to a different fixed point, the printed \(\Sigma_0^{\rm can}\), \(\mathcal S_{\rm can}\), \(\Pi_{\rm can}\), \(\widehat T_{m,\rm can}\) would change but both existing SymPy asserts would still pass (the bisection would still drive `g_fp` to `g_*` at whatever spurious root the modified solver hit, and `R = (g_* - r_F1)^2/(1+r_F1^2)` would still be 0.25 because it depends only on the hardcoded constants). The script's exit-0 would therefore not falsify the paper's quoted numerical answer. This is the textbook `insufficient_verification` shape: assertions exist, they pass, but they do not exercise the stage's load-bearing claim.

**Required change:**
Add four `assert`-style numeric checks at the end of the SymPy script (after the existing `R_can` assert and before the "Conclusion" print block), with tolerances matched to the Mathematica `expectApprox` calls so the two engines remain in lockstep:
```python
if abs(Sigma0_can - 4.651033550168867) > 1e-10:
    raise AssertionError("Sigma0_can deviates from notes value 4.651033550168867.")
if abs(S_can - 0.6703621156734617) > 1e-11:
    raise AssertionError("S_can deviates from notes value 0.6703621156734617.")
if abs(Pi_can - 3.871564377479002) > 1e-10:
    raise AssertionError("Pi_can deviates from notes value 3.871564377479002.")
if abs(T_hat_can - 1.4467083664567613) > 1e-10:
    raise AssertionError("T_hat_can deviates from notes value 1.4467083664567613.")
```
Tolerances are deliberately one decade looser than the Mathematica `expectApprox` tolerances (which are 1e-11 / 1e-12) to account for the absence of `phShift` overflow guarding in the SymPy code; this keeps both engines green at current numerical agreement while still flagging any future drift that exceeds five significant figures.

**Verification:**
After Codex applies, `redteam exec-sympy 156` should still exit 0 (the saved output already shows the printed values agree with the notes to 12+ digits). Re-running with a deliberately wrong N (e.g., N=241 instead of 2401) should now fail one of the new asserts where previously the script would still report PASS. The output transcript header (`# Status: PASS`) is unchanged; new printed assertion lines should appear in the transcript or the script's behaviour on intentional drift should change as above.

## Independent-derivation check (Mathematica)

The `.wl` script is structurally close to the `.py` script: identical grid (`n = 2401`, `Subdivide[0.0, 1.0, n-1]`), identical trapezoidal weights, identical kernel forms (`Cos[(Pi*xGrid)/2]`, `Cosh[kappa*(1 - xGrid)]/Cosh[kappa]`), identical initial guess (`piStar*Exp[-piStar*xGrid]`), identical Picard tolerance (`1.*^-13`, `500` iterations), identical bisection bracket (`[3.0, 6.0]`, `55` iters), and the Mathematica function names are direct renderings of the SymPy ones (`tsOperator` ↔ `Ts`, `tqOperator` ↔ `Tq`, `gMoment` ↔ `g`, `sMoment` ↔ `S`, `rMoment` ↔ `R`, `nextSigma` ↔ `next_sigma`, `solveFixedPoint` ↔ `solve_fixed_point`, `bracketRoot` ↔ `bracket_root`, `bisect` ↔ `bisect`). Both banners even read "STAGE 139 — RENORMALIZED CANONICAL BRANCH", which is an inherited typo (this is stage 156). On balance, however, this is **not** filed as `mathematica_transliteration`: (a) the algorithmic skeleton being shared is forced by the paper (it is the same integral equation and the same Picard scheme), not invented in either engine — there is no "answer-baked" substitution that would betray transliteration; (b) the Mathematica script adds substantive content the SymPy script lacks — an overflow-guard `phShift = ph - Min[ph]` (line 71), memoization of `fixedPointAt`/`gFp` (lines 89–90), an explicit monotonicity scan with `expectTrue` (lines 119–124), and six numeric `expectApprox` assertions covering the full paper deliverable (lines 147–152); (c) the two engines agree to numerical precision (`diff` values 4.99e-16 to 7.10e-15 in the Mathematica transcript) which is the right level of agreement for two independent discretizations of the same equation. The shared banner text ("STAGE 139") is a cosmetic carryover and is noted here for cleanup but is not load-bearing — it does not affect the numerical or symbolic content of either script.

## Engine cross-check

Both engines produce essentially identical numerical outputs (compare saved transcripts):

| Quantity | SymPy (transcript line 13–18) | Mathematica (transcript line 16–21) | Δ |
|---|---|---|---|
| \(\Sigma_0^{\rm can}\) | 4.651033550168867 | 4.651033550168874 | 7e-15 |
| \(\mathfrak g_{\rm can}\) | 0.758035078944663 | 0.7580350789446629 | 1e-16 |
| \(\mathcal S_{\rm can}\) | 0.6703621156734617 | 0.6703621156734617 | 0 |
| \(\mathcal R_{\rm can}\) | 0.2500000000000005 | 0.2500000000000005 | 0 |
| \(\Pi_{\rm can}\) | 3.871564377479002 | 3.871564377479008 | 6e-15 |
| \(\widehat T_{m,\rm can}\) | 1.4467083664567613 | 1.4467083664567622 | 9e-16 |

Agreement is at floating-point precision. `engines_agree: true`. No `engine_disagreement` finding.

## Freshness check

- SymPy script mtime: 2026-05-11 11:56; SymPy output mtime: 2026-05-11 12:47 → output newer than script.
- Mathematica script mtime: 2026-04-21 17:04; Mathematica output mtime: 2026-05-11 13:17 → output newer than script.

`outputs_fresh: true`. No `stale_output` finding.

## Verdict justification

Verdict: `findings`, one finding (F1, `insufficient_verification` on the SymPy side). The paper card and notes commit six load-bearing numerical outputs (\(\Sigma_0^{\rm can}\), \(\mathfrak g_{\rm can}\), \(\mathcal R_{\rm can}\), \(\mathcal S_{\rm can}\), \(\Pi_{\rm can}\), \(\widehat T_{m,\rm can}\)). The Mathematica script asserts all six plus a uniqueness-premise monotonicity scan, and is aligned. The SymPy script computes and prints the same six but `assert`s only the two that are guaranteed-to-pass by construction (the bisection's own convergence target and an input-constant identity). This is the canonical shape of `insufficient_verification`: the script's hard exit-on-failure does not cover the load-bearing claim. I attacked the rest of the audit and found no further defects: the constants (`rF1`, `g_star`, `Pi_star`, `Sigma0_star`, `T_hat_star`) are the upstream-canonical values quoted in earlier stages; the kernels (`c(x) = cos(pi x/2)`, `K_q(x) = cosh(kappa(1-x))/cosh(kappa)` with `kappa = pi/2`) match the abstract definitions in `eq:app-part04-gS-functionals`; the operator forms `T_s`, `T_q` are integral kernels coded by trapezoidal accumulation and match the shell/mixed Green operators referenced in the appendix; `R = (g - rF1)^2/(1 + rF1^2)` matches `eq:app-part04-R-Sigma`; the self-matched closure `T_hat_can = sqrt(9 Sigma0_can/20)` is the inverse of `Sigma_0 = (20/9) T_hat^2` in the notes §1; the two engines agree to floating-point precision. The shared "STAGE 139" banner text is a cosmetic carryover and is not filed as a finding (not category-fit). No `paper_misalignment` arises because every number printed in the scripts can be traced to the paper card, notes, or appendix.

## Self-test notes

I checked: (1) tautological-shape on the SymPy asserts — confirmed `g_can = g_*` is the bisection's drive target and `R_can = 1/4` is algebraically forced by the input constants once `g = g_*`, so neither asserts an independent physical fact; (2) algebraic round-trip on the closure `T_hat = sqrt(9 Sigma0/20)` and `Pi = Sigma0 (1 - R S)` — both match the notes and the upstream appendix block; (3) numerical sanity of the proposed new SymPy asserts — the printed SymPy values agree with the notes targets to 12+ digits, so the 1e-10/1e-11 tolerances chosen in the directive will pass cleanly while still catching any future drift; (4) transliteration trap — the shared banner "STAGE 139" and parallel function names are suspicious but the Mathematica side contains substantive independent content (monotonicity scan, overflow guard, six independent `expectApprox` assertions), and there is no answer-baked substitution that would betray transliteration. Conclusion: one substantive finding (F1) on the SymPy side; the Mathematica side is fine; the directive is safely applicable by Codex (file edits only, no algorithmic changes).
