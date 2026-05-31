---
unit_id: 184
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage184_branch_invariant_coordinates.md]
  paper_appendix: present
---

# Audit unit 184 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_184.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage184_branch_invariant_coordinates.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 99, 881-916)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.txt`

## What the paper claims

Stage 184 (\StatusExactClosure inside the coherent local D/N branch) claims that the three Stage-251/183 branch-adapted defect coordinates \((\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)\) are exactly the first grouped weak-axisymmetric logarithmic drifts of three exact branch composites. The `\stagefield{Output}` reads: "Identifies exact branch composites \(\mathfrak T_*\), \(\mathfrak N_*\), and \(\epsilon_\eta\) whose logarithmic drifts are the three normal-form coordinates." The distinct deliverables (notes §§1-6, appendix eqs `app-part05-{Rtr,Nstar,Tstar}-def` and `app-part05-branch-coordinate-laws`) are: (D1) the exact product identity \(R_{\rm target}\mathcal T^2=\Lambda_0(1-\epsilon_\eta)\); (D2) \(\delta\ln\mathfrak T_*=\Sigma_{\rm tr}\) for \(\mathfrak T_*:=R_{\rm tr}^{-C_*}\), \(C_*=(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)/(\chi_0\delta_U)\), with \(\Sigma_{\rm tr}=-C_*\Theta_1\); (D3) \(\delta\ln\mathfrak N_*=\Sigma_{\rm nt}\) for \(\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*}\), \(B_*=2(1+\chi_0+\delta_U)/\delta_U\), with \(\Sigma_{\rm nt}=\Xi_1+B_*\Theta_1\); (D4) \(\delta\ln\epsilon_\eta=\Sigma_\eta\) and the selected-branch complement \(\delta\ln[(R_{\rm target}\mathcal T^2)/\Lambda_0]=-(\epsilon_\eta/(1-\epsilon_\eta))\Sigma_\eta\); (D5) the forward zero-map: each drift vanishes when the corresponding observable drift vanishes.

## What the script claims to verify

Per the SymPy docstring (lines 8-15), the scripts verify all five deliverables: the exact product identity, the three branch-coordinate drift laws (tracking, nontracking, dressing), the selected-branch complement, and the forward zero map. The SymPy script genuinely constructs each composite (`Ttr = Rtr**(-Cstar)`, `Ntr = T2*Rtr**Bstar`, `Ecomp = (Rtarget*T2)/Lam0`, `eps_eta_var`) with `Rtr = Rtr0*exp(small*Theta1)` etc., then series-expands `log(composite/reference)` in `small` and extracts the first-order coefficient, comparing it to the Stage-251 closed forms `SigmaTr=-Cstar*Theta1`, `SigmaNT=Xi1+Bstar*Theta1`, `SigmaEta`. The Mathematica script's `expectZero` checks carry the same names, but (see findings) its load-bearing drift quantities are hardcoded rather than derived from the composites.

## Paper ↔ script cross-check

| Paper deliverable | SymPy check | Mathematica check | Status |
|---|---|---|---|
| D1 product identity \(R_{\rm target}\mathcal T^2=\Lambda_0(1-\epsilon_\eta)\) | line 62 (algebraic, genuine) | line 48 (algebraic, genuine) | match |
| D2 \(\delta\ln\mathfrak T_*=\Sigma_{\rm tr}\) | lines 65-69 (series from composite) | lines 51-55 (`dlnTtr` hardcoded `-cStar*theta1`; `tTr` unused) | partial (Mathematica) |
| D3 \(\delta\ln\mathfrak N_*=\Sigma_{\rm nt}\) | lines 72-76 (series from composite) | lines 58-62 (`dlnNtr` hardcoded `xi1+bStar*theta1`; `nTr` unused) | partial (Mathematica) |
| D4a \(\delta\ln\epsilon_\eta=\Sigma_\eta\) | lines 79-81 (series from `eps_eta_var`) | lines 65-67 (`dlnEpsEta` hardcoded `sigmaEta`) | partial (Mathematica) |
| D4b complement | lines 83-90 (series from `Ecomp`) | lines 69-73 (`dlnEcomp` hardcoded; `eComp` unused) | partial (Mathematica) |
| D5 forward zero map | lines 94-105 (subs into derived drifts) | lines 76-85 (subs into hardcoded drifts) | partial (Mathematica) |

Paper alignment is exact on constants and identities (`C_*`, `B_*`, `Σ_tr`, `Σ_nt`, product identity, complement all match paper verbatim). No `paper_misalignment`. The defect is engine-internal: the Mathematica engine does not independently exercise D2-D5. `paper_alignment: aligned` (the SymPy engine fully and faithfully exercises every deliverable; the gap is that the second engine does not independently corroborate).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `simplify(Rtarget*T2 - Lam0*(1-eps_eta_var)) == 0` | D1 | yes |
| A2 | sympy | 69 | `dln_Ttr - SigmaTr == 0`, `dln_Ttr` from `series(log(Ttr/Ttr0))` | D2 | yes |
| A3 | sympy | 76 | `dln_Ntr - SigmaNT == 0`, `dln_Ntr` from `series(log(Ntr/Ntr0))` | D3 | yes |
| A4 | sympy | 81 | `dln_eps_eta - SigmaEta == 0`, from `series(log(eps_eta_var/eps_eta))` | D4a | yes |
| A5 | sympy | 87-90 | `dln_Ecomp + eps_eta/(1-eps_eta)*SigmaEta == 0`, from `series(log(Ecomp/Ecomp0))` | D4b | yes |
| A6 | sympy | 94-96 | `SigmaTr - dln_Ttr == 0` etc. | D2-D4 (mirror) | partial (re-checks A2-A4) |
| A7 | sympy | 103-105 | zero-map subs into derived drifts | D5 | yes |
| A8 | mathematica | 48 | `expectZero[..., rTarget*t2 - lam0*(1-epsEtaVar)]` | D1 | yes |
| A9 | mathematica | 55 | `expectZero["...T_* - Sigma_tr", dlnTtr - sigmaTr]`, `dlnTtr := -cStar*theta1` | D2 | **no (tautological)** |
| A10 | mathematica | 62 | `expectZero["...N_* - Sigma_nt", dlnNtr - sigmaNT]`, `dlnNtr := xi1+bStar*theta1` | D3 | **no (tautological)** |
| A11 | mathematica | 67 | `expectZero["...eps_eta - Sigma_eta", dlnEpsEta - sigmaEta]`, `dlnEpsEta := sigmaEta` | D4a | **no (tautological)** |
| A12 | mathematica | 73 | `expectZero["...complement", dlnEcomp + epsEta*sigmaEta/(1-epsEta)]`, `dlnEcomp := -epsEta*sigmaEta/(1-epsEta)` | D4b | **no (tautological)** |
| A13 | mathematica | 76-78 | mirror checks `sigmaTr-dlnTtr` etc. | D2-D4 | **no (tautological)** |
| A14 | mathematica | 83-85 | zero-map subs into hardcoded drifts | D5 | **no (degenerate)** |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl:51-85`

**What's wrong:**
The Mathematica script's load-bearing drift quantities are hardcoded to the very expressions they are then compared against, so the four core drift-law assertions cannot fail regardless of the physics.

- Line 53: `dlnTtr = FullSimplify[-cStar*theta1, ...]`. Line 40: `sigmaTr = FullSimplify[-cStar*theta1, ...]`. Line 55 then asserts `dlnTtr - sigmaTr == 0`, i.e. `(-cStar*theta1) - (-cStar*theta1) == 0`. The composite `tTr = rTr^(-cStar)` defined at line 51 (and `tTr0` at line 52) is never used to derive `dlnTtr`.
- Line 60: `dlnNtr = FullSimplify[xi1 + bStar*theta1, ...]` vs `sigmaNT = FullSimplify[xi1 + bStar*theta1, ...]` (line 41); line 62 asserts their difference is zero. The composite `nTr = t2*rTr^bStar` (line 58, `nTr0` line 59) is never used.
- Line 65: `dlnEpsEta = FullSimplify[sigmaEta, ...]`; line 67 asserts `dlnEpsEta - sigmaEta == 0`. `epsEtaVar` (line 46) is never log-differentiated.
- Line 71: `dlnEcomp = FullSimplify[-epsEta*sigmaEta/(1 - epsEta), ...]`; line 73 asserts `dlnEcomp + epsEta*sigmaEta/(1-epsEta) == 0`. The composite `eComp = (rTarget*t2)/lam0` (line 69, `eComp0` line 70) is never used.
- Lines 76-78 (`sigmaTr - dlnTtr` etc.) and 80-85 (zero-map substitutions into the hardcoded `dlnTtr/dlnNtr/dlnEpsEta`) inherit the same defect: substituting `theta1->0` into the hardcoded `-cStar*theta1` is not a test that the composite `R_tr^(-C_*)` has the claimed drift.

By contrast, the SymPy script (lines 65-90) genuinely constructs each composite and recovers the drift via `series(log(composite/reference), small, 0, 2)`, so its checks DO depend on the composite definition and the closed forms of `C_*`/`B_*` being mutually consistent. The Mathematica engine therefore provides no independent corroboration of D2-D5; if `C_*` or the composite exponent were wrong, the Mathematica script would still PASS.

**Why this matters:**
The second-engine policy requires both engines to derive the result independently from the physical premises. As written, the Mathematica script asserts `expr == expr` for every drift law and thus verifies nothing about the branch composites; it merely confirms that two copies of the same hardcoded literal are equal. The "PASS" lines in the saved transcript (`mathematica/output/...stage184...txt` lines 27, 33, 41) are vacuous for D2-D5. Only A8 (product identity, line 48) is a genuine Mathematica check.

**Required change:**
Make the Mathematica drifts genuinely derived from the composites, mirroring the SymPy structure but via Mathematica's own series machinery (not by copying the SymPy intermediate algebra). Replace the hardcoded `dlnTtr`, `dlnNtr`, `dlnEpsEta`, `dlnEcomp` with `Series`-based extractions of the first-order `small` coefficient of `Log[composite/reference]`, using the already-defined `tTr/tTr0`, `nTr/nTr0`, `epsEtaVar`, `eComp/eComp0`. See directive F1 for the exact replacements.

**Verification:**
After the fix, lines 53/60/65/71 should compute the drift via `SeriesCoefficient[Log[.../...], {small, 0, 1}]` (or equivalent `Normal[Series[...]]` + coefficient extraction) rather than restating the target. The transcript should still print `delta ln T_* = -(((1+chi0)*(1+deltaU)*(1+chi0+deltaU)*theta1)/(chi0*deltaU))`, `delta ln N_* = (2*(1+chi0+deltaU)*theta1)/deltaU + xi1`, `delta ln eps_eta = sigmaEta`, and the complement, but now those values must be produced by the series operation. All `expectZero` lines must still report `= 0` / PASS, and the script must exit 0.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl:51-52,58-59,69-70`

**What's wrong:**
As a direct consequence of F1, the composite definitions `tTr`/`tTr0` (51-52), `nTr`/`nTr0` (58-59), and `eComp`/`eComp0` (69-70) are dead code — defined but never referenced by any assertion. The Mathematica engine therefore never exercises its own stated claim (the docstring/banner promises the branch-invariant coordinate laws, but no assertion touches the branch composites). This is the same root cause as F1; resolving F1 by routing the drifts through these composites simultaneously fixes F2 (the dead definitions become load-bearing).

**Why this matters:**
A reader auditing the Mathematica side would see composites defined and assume they are tested, when in fact the only thing tested is an `expr - expr` tautology. If F1 is fixed without consuming these definitions (e.g., by some other route), the dead code should be removed; the intended fix is to consume them.

**Required change:**
Covered by directive F1 — wire `tTr/tTr0`, `nTr/nTr0`, `eComp/eComp0`, and `epsEtaVar` into the series-based drift derivations. No separate edit; verifying F1 verifies F2.

**Verification:**
After F1's fix, `grep` for `tTr`, `nTr`, `eComp` should show each used in a `Log[.../...]` series call, not merely assigned. No composite remains write-only.

## Independent-derivation check (Mathematica)

The `.wl` is partly a transliteration in structure (same banner text, same variable names `sigmaTr/sigmaNT`, same five-section choreography as the `.py`), but the more serious issue is that for D2-D5 it is not even a transliteration of the SymPy *derivation* — the SymPy script performs `series(log(composite/reference))` to obtain each drift, whereas the Mathematica script skips the derivation entirely and assigns the known answer to `dlnTtr/dlnNtr/dlnEpsEta/dlnEcomp`. So the Mathematica side is weaker than a transliteration would be: it neither re-derives independently nor faithfully ports the SymPy algebra; it hardcodes the conclusion. The product-identity check (line 48) is a genuine, independent algebraic check and matches the SymPy line 62. I file this as `tautological_check` (F1) rather than `mathematica_transliteration` because the defect is that the assertions are vacuous, not that they echo SymPy's intermediate steps.

## Engine cross-check

Both transcripts report identical final values and all-zero residuals:
- SymPy `delta ln T_* = Theta_1*(-chi0**2*deltaU - chi0**2 - ... - 1)/(chi0*deltaU)` (output line 24) and Mathematica `delta ln T_* = -(((1+chi0)*(1+deltaU)*(1+chi0+deltaU)*theta1)/(chi0*deltaU))` (output line 25) are the same expanded-vs-factored form.
- `delta ln N_*`: SymPy `2*Theta_1*chi0/deltaU + 2*Theta_1 + 2*Theta_1/deltaU + Xi_1` (line 30) = Mathematica `(2*(1+chi0+deltaU)*theta1)/deltaU + xi1` (line 32). Same.
- All `expectZero`/`expect_zero` residuals are 0; both scripts exit 0.

The numbers agree, so there is no `engine_disagreement`. But the agreement is hollow for D2-D5: the Mathematica value agrees because it was hardcoded to the SymPy answer, not because it was independently derived. Engine agreement is therefore `true` at the level of printed values but does not constitute independent corroboration (which is why F1/F2 are filed).

## Freshness

Outputs are fresh. SymPy: script mtime 1778522333 < output mtime 1778525299. Mathematica: script mtime 1778522213 < output mtime 1778527452. No `stale_output` finding. (Note: the F1/F2 fix will require a re-run, which the verifier triggers.)

## Verdict justification

`findings` (not clean, not stop_cold). The paper claim is exactly mirrored by the SymPy script, which genuinely derives all five deliverables from the composite definitions via series expansion — I attacked the SymPy checks for tautology (the drifts are computed from `log(composite/ref)`, not asserted), for wrong constants (`C_*`, `B_*` match the paper verbatim, lines 39-40), for symbol-domain errors (`positive=True` on `chi0, deltaU, eps_eta` is consistent with the coherent positive-coupling state space the paper assumes; `eps_eta>0` plus the `1-eps_eta` denominator is fine since no assertion divides by it in a way requiring `eps_eta<1`, and the complement check is purely algebraic), and for missing branches (the forward zero-map is the only direction the paper's notes claim to verify in-script; the reverse `⟸` is asserted in prose via `C_*>0`, not required of the script) — all SymPy attacks failed, so the SymPy engine is sound and paper-aligned. The Mathematica engine, however, hardcodes the four drift quantities equal to their targets (F1) and leaves the branch composites as dead code (F2), so it provides no independent verification of D2-D5. This is a real second-engine defect, fixable by Codex (route the Mathematica drifts through `Series[Log[composite/reference]]`), hence `findings` with a directive, `stop_cold: null`. No downstream propagation: the SymPy verification already establishes the math, and the fix only strengthens the Mathematica mirror without changing any forward-carried constant.

## Self-test notes

I checked: (1) Variable independence — the proposed Mathematica series are `Log[tTr/tTr0]`, `Log[nTr/nTr0]`, `Log[epsEtaVar/epsEta]`, `Log[eComp/eComp0]`, each genuinely depending on `small` through the `Exp[small*theta1]`/`(1+small*sigmaEta)` factors, so the first-order `small` coefficient is nonzero and matches the targets (confirmed against the already-derived SymPy values in the transcript). (2) Trivial-case — substituting `theta1->0` into the genuinely-derived `dlnTtr=-cStar*theta1` gives 0 (zero-map A14 stays valid post-fix); the complement reduces to `-epsEta*sigmaEta/(1-epsEta)` so `dlnEcomp + epsEta*sigmaEta/(1-epsEta) -> 0`. (3) Paper round-trip — the fix introduces no new constants; it routes through `cStar`/`bStar` already defined to the paper's values, so no new `paper_misalignment` is created. (Minor non-blocking observation: both scripts' banner reads "STAGE 167" though this is stage 184 — cosmetic, folded into directive F1 as an optional label fix; not a math finding.)
