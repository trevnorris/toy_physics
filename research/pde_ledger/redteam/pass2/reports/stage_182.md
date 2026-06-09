---
unit_id: 182
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage182_microscopic_coherent_slippage.md]
  paper_appendix: present
---

# Audit unit 182 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_182.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage182_microscopic_coherent_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 95, 768-817 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.txt`

## What the paper claims

Stage 182 pushes the Stage-181 grouped weak-axisymmetric coherent defect down to microscopic coherent-kernel slippages. `\stagefield{Output}`: *"Reduces the defect to \((\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta)\) plus dressing slippage \(\Sigma_\eta\), and isolates \(\Sigma_{\rm tr}\)."* The notes give the full deliverable set: (1) the five microscopic slippage definitions (\(\Sigma_Z,\Sigma_\chi,\Sigma_\eta,\Sigma_\epsilon,\Sigma_\delta\)) in terms of the eight microscopic log drifts; (2) the boxed four-slippage grouped-defect law \(\Xi_1 = \Sigma_Z + \frac{2\chi_0}{1+\chi_0}\Sigma_\chi + \frac{2\epsilon_W}{1-\epsilon}\big[\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon - \frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta\big]\); (3) the selected-branch demand law \(\mathcal R_1 = -\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta - \Xi_1\); (4) the tracking-factor factorization \(\Theta_1 = -\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\Sigma_{\rm tr}\) with \(\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi\); (5) the tracking/nontracking split of \(\Xi_1\); and (6) microscopic support-blindness (no \(\lambda_\phi,K_\phi^{(\mathrm{eff})}\) drift enters). The appendix (eq:app-part05-Xi1-slippage-law, eq:app-part05-Sigma-tr-def, eq:app-part05-Theta-Sigma-tr) restates the four-slippage law, the \(\Sigma_{\rm tr}\) definition, and the \(\Theta_1\) factorization verbatim.

## What the script claims to verify

The SymPy docstring enumerates six checks matching the six deliverables. The assertions: (1) the five slippage definitions equal the corresponding physical-drift ratios (\(\Sigma_\chi=\chi_1/\chi_0\) etc.); (2) the direct \(\Xi_1\) equals the boxed slippage form; (3) the \(\mathcal R_1\) law; (4) the \(\Theta_1\) factorization via \(\Sigma_{\rm tr}\); (5) the tracking/nontracking split; (6) no support symbol (`lamphi1`, `kphi`) appears in `Xi1_direct`/`R1_direct`/`Theta1_direct` free symbols. The Mathematica script tests the same six but, for the central \(\Xi_1\) law, takes a structurally different route: it *solves* the slippage definitions for the microscopic gauges, substitutes back, confirms the gauges are fully *eliminated* (`FreeQ`), and then extracts each \(\Sigma\)-coefficient (`Coefficient[...]`) and checks it equals the paper's stated coefficient.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Five slippage defs = drift ratios | py 74-78 / wl 54-58 (`expect_zero`/`expectZero` of def − ratio) | match |
| Four-slippage \(\Xi_1\) law | py 81-92 (construct both, assert diff 0) / wl 71-89 (gauge-eliminate + per-coeff `Coefficient` vs paper coeffs) | match |
| \(\mathcal R_1\) selected-branch law | py 97-99 / wl 93-95 | match |
| \(\Theta_1\) factorization via \(\Sigma_{\rm tr}\) | py 104-113 / wl 100-110 | match |
| Tracking/nontracking split | py 120-132 / wl 115-129 | match |
| \(\Sigma_{\rm tr}\) definition | py 112 / wl 109 | match |
| Microscopic support-blindness | py 141-146 / wl 136-146 | match |

Every paper-side deliverable maps to a non-tautological script-side check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74-78 | `simplify(def_sub − ratio) == 0` | slippage defs | yes |
| A2 | sympy | 92 | `simplify(Xi1_direct − Xi1_slip.subs) == 0` | four-slippage law | yes |
| A3 | sympy | 99 | `simplify(R1_direct − R1_slip.subs) == 0` | \(\mathcal R_1\) law | yes |
| A4 | sympy | 113 | `Theta1_direct − Theta1_fact(Sigma_tr_def).subs == 0` | \(\Theta_1\) factorization | yes |
| A5 | sympy | 129-132 | `Xi1_split(Sigma_tr_def) − Xi1_slip == 0` | tracking split | yes |
| A6 | sympy | 141-146 | support symbols ∉ free_symbols | support-blindness | yes |
| B1 | math | 54-58 | `expectZero[def − ratio]` | slippage defs | yes |
| B2 | math | 78-87 | gauge `FreeQ` + `Coefficient` vs paper coeffs + constant-term 0 | four-slippage law | yes (independent) |
| B3 | math | 95 | `expectZero[r1Direct − r1Slip.subs]` | \(\mathcal R_1\) law | yes |
| B4 | math | 110 | `expectZero[theta1Direct − theta1Fact(...)]` | \(\Theta_1\) factorization | yes |
| B5 | math | 125-129 | split per-coeff `Coefficient` vs paper + slippage-form round | tracking split | yes (independent coeff route) |
| B6 | math | 136-146 | `FreeQ[form, lamphi1] && FreeQ[form, kphi]` | support-blindness | yes |

No tautological rows. Every check traces to a specific paper deliverable; no orphaned scaffolding.

## Findings

None.

I attacked the following and each held:

- **Tautology check on A2 (the central law).** The SymPy A2 constructs `Xi1_slip` (the target form) then asserts `Xi1_direct − Xi1_slip.subs(slip_subs) == 0`. This is NOT tautological because `Xi1_direct` (line 81-83) is built from the *independent* Stage-181 expression `zetaZ − omegaW + 2*chi1/(1+chi0) + 2*eps1/(1−eps)`, where `chi1`, `eps1`, `eps` are defined from the raw microscopic logs (lines 46-56), while `Xi1_slip` is written in the abstract \(\Sigma\) symbols whose `slip_subs` are the *separately stated* combinations (lines 65-71). For the assertion to pass, the Stage-181 physical form must genuinely collapse onto the paper's coefficient structure — it can fail if any coefficient is wrong. Confirmed: the Mathematica route (B2) reaches the same coefficients by an algebraically distinct path (`Solve` + `Coefficient`), so the law is not a definitional identity.

- **Mathematica B2 coefficient values vs paper.** wl 83 checks coeff(\(\Sigma_\chi\)) = `2*chi0/(1+chi0)` — matches paper. wl 85 checks coeff(\(\Sigma_\epsilon\)) = `2*epsW*(11+9*deltaU)/(11*(1-eps)*(1+deltaU))` — equals paper's \(\frac{2\epsilon_W}{1-\epsilon}\cdot\frac{11+9\delta_U}{11(1+\delta_U)}\). wl 86 checks coeff(\(\Sigma_\delta\)) = `−4*epsW*deltaU/(11*(1-eps)*(1+deltaU)^2)`, equals paper's \(\frac{2\epsilon_W}{1-\epsilon}\cdot(-\frac{2\delta_U}{11(1+\delta_U)^2})\). coeff(\(\Sigma_Z\))=1, coeff(\(\Sigma_\eta\))=0, constant term 0 — all match. No `value_mismatch`.

- **Symbol assumptions.** Background ratios `chi0, epsW, eps_eta, deltaU` declared `positive=True` (py 43) / `chi0>0 && epsW>0 && epsEta>0 && deltaU>0` (wl 29-30). Positivity is justified by the physical branch \(\chi_0>0,\delta_U>0\) the notes state (notes §6, "on the physical branch \(\chi_0>0\), \(\delta_U>0\)"); positivity of `epsW`, `epsEta` is benign for these identities (used only so denominators `1+chi0`, `1+deltaU` and `1-eps` are nonzero). No `symbol_assumption_error`. Note `eps = epsW(1 − 2δ_U/(11(1+δ_U)))` is not asserted `<1`; the `1−eps` denominator is symbolic and never numerically inverted, so no hidden singularity assumption.

- **Support-blindness scope.** Both scripts check only that no support log was *wired into* the construction (py comment 138-140, wl comment 134-135 explicitly say the \(\zeta\)-cancellation lives upstream in Stage 249). This is honestly scoped: the notes §2.2/§8 claim support-blindness "because \(\zeta\) drops out (proven Stage 181)", and the script does not over-claim to re-prove the cancellation — it verifies the weaker, in-scope statement. Not `insufficient_verification` because the paper deliverable for *this* stage is the microscopic decomposition, not a re-derivation of the \(\zeta\)-cancellation.

- **Output freshness.** sympy.py mtime 2026-05-30 01:22:19 < its .txt 01:40:24; wl mtime 01:26:17 < its .txt 01:40:40. Both outputs are fresh. No `stale_output`.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent** verification route, not a transliteration, on the load-bearing claim. Strongest evidence: the central four-slippage law. SymPy (py 84-92) writes the answer `Xi1_slip` by hand and asserts equality:
```
Xi1_slip = simplify(SigmaZ + 2*chi0/(1+chi0)*SigmaChi + 2*epsW/(1-eps)*(...))
expect_zero("Xi_1 direct - slippage form", Xi1_direct - Xi1_slip.subs(slip_subs))
```
Mathematica (wl 69-87) instead *derives* the slippage form: it solves the slippage definitions for the microscopic gauges, substitutes into `xi1Direct`, proves the gauges are eliminated (`FreeQ[xi1DirectInSigmas, Alternatives @@ logSyms]`, wl 78-81), then reads off each coefficient with `Coefficient[Expand[xi1DirectInSigmas], sigmaZ]`, `...sigmaChi`, etc. (wl 82-87) and checks each against the paper coefficient. SymPy never calls `Solve` to invert the gauges and never uses `Coefficient`; Mathematica never constructs a hand-written target form for the gauge-eliminated \(\Xi_1\). These are genuinely different algorithms (assert-equality-to-target vs. solve-eliminate-and-extract-coefficients). The auxiliary checks (slippage defs, \(\mathcal R_1\), \(\Theta_1\)) are closer in shape between the engines (both `expect_zero` of a residual), but the engine-defining claim is reached by distinct routes, so this is not a port. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass every check and agree on the final closed forms:
- \(\Xi_1\): sympy pprint (sympy.txt 19-31) and Mathematica `xi1Slip` (math.txt 38) are the same rational function in \(\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta,\chi_0,\epsilon_W,\delta_U\). Note Mathematica's printed `Xi_1` is written with `(11+9*deltaU)*epsW − 11*(1+deltaU)` in the denominator, which is exactly `−11(1+δ_U)(1−eps)` after substituting `eps`; algebraically identical.
- \(\Sigma_{\rm tr} = (1+\chi_0)\Sigma_\delta + (1+\delta_U)\Sigma_\chi\): math.txt 52 `sigmaChi + deltaU*sigmaChi + sigmaDel + chi0*sigmaDel`; sympy.txt 57-58 same. Agree.
- \(\Theta_1\): both `−χ₀δ_U Σ_tr / ((1+χ₀)(1+δ_U)(1+χ₀+δ_U))` (sympy.txt 59-62, math.txt 53). Agree.

No `engine_disagreement`.

## Verdict justification

Verdict **clean**. I read the paper card, the notes (full derivation), and the appendix rows before opening the scripts, and the scripts' verified claims match the paper's six deliverables exactly — same slippage definitions, same four-slippage \(\Xi_1\) law with the same `11+9δ_U` / `11(1+δ_U)²` coefficients, same \(\Sigma_{\rm tr}\) combination, same \(\Theta_1\) factorization, same support-blindness scope. I attacked the central A2 assertion as a possible construct-both-and-compare tautology and it survived: `Xi1_direct` is built from the independent Stage-181 physical form, and the Mathematica engine reaches the identical coefficients by a structurally different solve-and-extract route, which simultaneously refutes both the tautology and transliteration attacks. Symbol positivity is justified by the physical branch; outputs are fresh; engines agree.

## Self-test notes

Variable-independence: no `sp.diff`/`D[]` in either script, so the zero-derivative trap does not apply. Symmetry/parity: no integrals — N/A. Trivial-case: substituting \(\Sigma_Z=\Sigma_\chi=\Sigma_\epsilon=\Sigma_\delta=0\) into the \(\Xi_1\) law gives 0 (matches notes §7.1), and a single nonzero \(\Sigma_\chi\) yields the nonzero coefficient `2χ₀/(1+χ₀)` checked at wl 83 — both behave correctly. Paper round-trip: every script coefficient I traced (Σ_ε, Σ_δ) matches the paper's eq:app-part05-Xi1-slippage-law, so no fix is needed and no new misalignment is introduced.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 8 values checked, 0 misaligned

Authoritative records: SymPy source + `scripts/output/..._sympy_audit.txt`, Mathematica source + `mathematica/output/..._mathematica_audit.txt` (both fresh; not executed, reasoned from source + transcript).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(\Sigma_Z = 2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W\) | py 66 / wl 46 / sympy.txt 94 | notes:240-243 (`\Sigma_Z:=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W`) | MATCH |
| \(\Sigma_\chi = \gamma_1+c_1-\kappa_U\) | py 67 / wl 47 / sympy.txt 95 | notes:245-248 | MATCH |
| \(\Sigma_\eta = 2c_1-\kappa_U-\kappa_\eta\) | py 68 / wl 48 / sympy.txt 96 | notes:250-253 | MATCH |
| \(\Sigma_\epsilon = 2\gamma_1+2\lambda_1-\kappa_U-\kappa_W\) | py 69 / wl 49 / sympy.txt 97 | notes:255-258 | MATCH |
| \(\Sigma_\delta = \tau_1-\kappa_U\) | py 70 / wl 50 / sympy.txt 98 | notes:260-262 | MATCH |
| Four-slippage law \(\Xi_1 = \Sigma_Z + \frac{2\chi_0}{1+\chi_0}\Sigma_\chi + \frac{2\epsilon_W}{1-\epsilon}[\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon - \frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta]\) | py 84-92 / wl 82-89 / sympy.txt 19-31, math.txt 38 | tex:786-800 (eq:app-part05-Xi1-slippage-law) + notes:284-295 | MATCH |
| \(\Sigma_{\rm tr} = (1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi\) | py 112 / wl 109 / sympy.txt 57-58, math.txt 52 | tex:802-805 (eq:app-part05-Sigma-tr-def) + notes:371-374 | MATCH |
| \(\Theta_1 = -\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\Sigma_{\rm tr}\) | py 109-117 / wl 105-112 / sympy.txt 59-62, math.txt 53 | tex:807-814 (eq:app-part05-Theta-Sigma-tr) + notes:377-382 | MATCH |
| \(\mathcal R_1 = -\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta - \Xi_1\) (selected-branch law) | py 97-101 / wl 93-96 / sympy.txt 36-51, math.txt 45 | notes:319-325, 328-340 (not in tex, but in notes — natural carrier) | MATCH |

Note on \(\mathcal R_1\): the appendix and stage card do not print the \(\mathcal R_1\) law (the card's `Output` is terse), but the notes carry it explicitly (§5), so per the augmentation guard this is a MATCH (lives correctly in the `.md`), not MISSING-DELIVERABLE. Likewise the tracking/nontracking split \(\Xi_1\) form (notes §6, eq at notes:399-412) is verified by py 120-134 / wl 115-130 and lives in the notes; the terse card legitimately omits it.

INTERNAL scaffolding (no finding): per-check residuals printed as `= 0`, `PASS:`/`FAIL:` flags, `Xi_1 microscopic gauges eliminated` FreeQ flag, support-symbol-leakage `set()`/`support-blind` flags, the intermediate `eps`/`eps1`/`chi1`/`varepsW`/`deltaU1`/`zetaZ`/`omegaW`/`eta1` background-drift expressions (exist only to drive the assertions), and the human-readable carry-forward print block (which restates the five slippage defs already reconciled above).

All 8 deliverable-level values reconcile against the `.tex` and/or `.md`. No MISMATCH, no MISSING-DELIVERABLE.

## Numbering note (informational, not fixed)

No stale "Stage NNN" self-labels in either script (banners and carry-forward text say only "STAGE 182"). One historical-reference comment: py 45 / wl (implicit via `zetaZ` etc.) says "Stage-30 coherent branch definitions" (py 45 `# Stage-30 coherent branch definitions.`) and the notes §1 attribute the microscopic couplings to "Stage 047". The notes also reference "Stage 249" for the \(\zeta\)-cancellation (py 139, wl 135). These are cross-references to other stages, not this stage's self-label, and I did not verify their correctness (out of scope). `\stagefield{Purpose}` is generic boilerplate naming "Stage~182" correctly — no +17 drift observed.
</content>
</invoke>
