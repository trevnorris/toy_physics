---
unit_id: 183
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:43:16-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage183_triangular_normal_form.md]
  paper_appendix: present
---

# Audit unit 183 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_183.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage183_triangular_normal_form.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 97 status row; block narrative lines 819–916)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.txt`

## What the paper claims

Stage 183 ("Triangular normal form"; notes title "Stage 251") compresses the coherent weak-axisymmetric defect problem from five microscopic slippages `(Σ_Z, Σ_χ, Σ_ε, Σ_δ, Σ_η)` to three branch-adapted coordinates `(Σ_tr, Σ_nt, Σ_η)`, where `Σ_nt` is defined (notes eqs. boxed lines 99–108; appendix eq. `app-part05-Sigma-nt-def`) as `Σ_Z + (2ε_W/(1-ε))·(11+9δ_U)/(11(1+δ_U))·Σ_ε − [2χ₀/(1+δ_U) + 4ε_Wδ_U/(11(1-ε)(1+δ_U)²)]·Σ_δ`. The verbatim Output field reads: "Compresses the coherent weak-axisymmetric problem to \((\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)\) with triangular map to \((\Theta_1,\Xi_1,\mathcal R_1)\)." The triangular map (notes Sec. 2; appendix eq. `app-part05-triangular-normal-form`) is `Θ₁ = −C_tr·Σ_tr`, `Ξ₁ = A_tr·Σ_tr + Σ_nt`, `R₁+Ξ₁ = −ε_η/(1−ε_η)·Σ_η`, with `C_tr = χ₀δ_U/((1+χ₀)(1+δ_U)(1+χ₀+δ_U))` and `A_tr = 2χ₀/((1+χ₀)(1+δ_U))`. Three further deliverables: the exact inverse reconstruction formulas (notes Sec. 3, also Sec. 5 item 3), and the triple-rigidity theorem `Θ₁=Ξ₁=R₁=0 ⟺ Σ_tr=Σ_nt=Σ_η=0` on the constructive branch `χ₀>0, δ_U>0, 0<ε_η<1` (notes Sec. 4.4; appendix eq. `app-part05-zero-normal-form`).

## What the script claims to verify

The docstring enumerates four checks: (1) reduction to the three branch-adapted coordinates, (2) triangularity of the observable ledger, (3) correctness of the inverse reconstruction formulas, (4) the triple-rigidity theorem collapsing to `Σ_tr=Σ_nt=Σ_η=0`. Both engines import the Stage-250 expressions for `Θ₁`, `Ξ₁`, `R₁` (these are taken as given, not re-derived from PDE physics), define `Σ_nt` per the boxed formula, then assert: `Ξ₁ − (A_tr·Σ_tr + Σ_nt) = 0`; `R₁+Ξ₁ + ε_η/(1−ε_η)·Σ_η = 0`; the inverse-map round-trips `Σ_tr_inv − Σ_tr = 0`, `A_tr/C_tr − 2(1+χ₀+δ_U)/δ_U = 0`, `Σ_nt_inv − Σ_nt = 0`, `Σ_eta_inv − Σ_η = 0`; and three "triple-rigidity" forward substitutions setting the slippages to zero.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `Σ_nt` definition (boxed) | `SigmaNT_def`, sympy:67–74 / wl:47–52 | match (literal-for-literal with notes/appendix) |
| Triangular split `Ξ₁ = A_tr·Σ_tr + Σ_nt` | `expect_zero("Xi_1 - (A_tr Sigma_tr + Sigma_nt)")`, sympy:81 / wl:57 | match |
| `R₁+Ξ₁ = −ε_η/(1−ε_η)·Σ_η` | `expect_zero(...)`, sympy:82 / wl:58 | match |
| `Θ₁ = −C_tr·Σ_tr` | built in by definition `Theta1`, sympy:50–52 / wl:33–36 | match (imported; triangularity in Σ_tr only is by construction, consistent with a carry-forward stage) |
| Inverse `Σ_tr = −((1+χ₀)(1+δ_U)(1+χ₀+δ_U)/(χ₀δ_U))·Θ₁` | sympy:93–94 / wl:66–67 | match |
| `A_tr/C_tr = 2(1+χ₀+δ_U)/δ_U` | sympy:96–98 / wl:69–71 | match |
| Inverse `Σ_nt = Ξ₁ + 2(1+χ₀+δ_U)/δ_U·Θ₁` | sympy:100–101 / wl:73–74 | match |
| Inverse `Σ_η = −((1−ε_η)/ε_η)(R₁+Ξ₁)` | sympy:103–104 / wl:76–77 | match |
| Triple-rigidity `Θ₁=Ξ₁=R₁=0 ⟺ Σ_tr=Σ_nt=Σ_η=0` | forward subs only, sympy:109–114 / wl:80–85 | partial (forward "⟸" is tautological; the substantive "⟹" rigidity rides on the inverse block, not on this block) |

Every paper-side deliverable has a corresponding script-side check, and the constants and coordinate definitions match the paper literally. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 81 | `simplify(Xi1 - (A_tr·SigmaTr + SigmaNT_def)) == 0` | triangular split of Ξ₁ | yes (catches a copy error in the tr/nt split) |
| A2 | sympy | 82 | `simplify(R1 + Xi1 + eps_eta/(1-eps_eta)·SigmaEta) == 0` | R₁+Ξ₁ collapses to Σ_η | yes |
| A3 | sympy | 94 | `simplify(SigmaTr_inv - SigmaTr) == 0` | inverse Σ_tr | yes (catches wrong reciprocal of C_tr) |
| A4 | sympy | 98 | `simplify(A_tr/C_tr - 2(1+chi0+deltaU)/deltaU) == 0` | ratio identity | yes |
| A5 | sympy | 101 | `simplify(SigmaNT_inv - SigmaNT_def) == 0` | inverse Σ_nt | yes |
| A6 | sympy | 104 | `simplify(SigmaEta_inv - SigmaEta) == 0` | inverse Σ_η | yes |
| A7 | sympy | 112 | `simplify(Theta1.subs(SigmaTr→0)) == 0` | triple-rigidity (⟸) | no — tautological (0·anything) |
| A8 | sympy | 113 | `simplify((A_tr·SigmaTr + SigmaNT).subs({SigmaTr→0,SigmaNT→0})) == 0` | triple-rigidity (⟸) | no — tautological |
| A9 | sympy | 114 | `simplify((-eps_eta/(1-eps_eta)·SigmaEta).subs(SigmaEta→0)) == 0` | triple-rigidity (⟸) | no — tautological |
| B1–B9 | mathematica | 57,58,67,71,74,77,83,84,85 | `expectZero[...]` mirroring A1–A9 | same as A1–A9 | same verdicts as A1–A9 (B7–B9 tautological) |

A1–A6 (and their B-mirrors) genuinely exercise the triangular split and the invertibility; the inverse round-trips are the load-bearing proof of the rigidity `⟹` direction. A7–A9 / B7–B9 are tautological forward substitutions.

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py:106-114`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl:79-85`

**What's wrong:**
The "Triple-rigidity theorem" block claims to verify the paper's headline equivalence (notes Sec. 4.4 / appendix eq. `app-part05-zero-normal-form`):
`Θ₁ = Ξ₁ = R₁ = 0 ⟺ Σ_tr = Σ_nt = Σ_η = 0` on the branch `χ₀>0, δ_U>0, 0<ε_η<1`.
But the three assertions only substitute the slippages to zero into linear forms and confirm the result is zero:

sympy:109–114
```python
Theta_zero = sp.simplify(Theta1.subs({SigmaTr: 0}))
Xi_zero    = sp.simplify((A_tr * SigmaTr + SigmaNT).subs({SigmaTr: 0, SigmaNT: 0}))
Rsum_zero  = sp.simplify((-eps_eta / (1 - eps_eta) * SigmaEta).subs({SigmaEta: 0}))
expect_zero("Theta_1|(Sigma_tr=0)", Theta_zero)
expect_zero("Xi_1|(Sigma_tr=Sigma_nt=0)", Xi_zero)
expect_zero("(R_1+Xi_1)|(Sigma_eta=0)", Rsum_zero)
```
(wl:80–85 is the identical mirror.) These only exercise the trivial "⟸" direction (zero slippage ⟹ zero observable), which holds for any linear map regardless of the prefactors and so cannot fail no matter what the physics is. The author's comment at sympy:107–108 ("The map is triangular, so vanishing observables imply vanishing adapted slippages. We verify the forward zero map directly.") explicitly claims the non-trivial "⟹" content (vanishing observables ⟹ vanishing slippages), but that direction requires the prefactors `C_tr`, `A_tr`-channel, and `ε_η/(1−ε_η)` to be nonzero on the branch — which this block never tests. Note: bare symbols `SigmaNT` (sympy:47/wl:28) are used here, not the constructed `SigmaNT_def`, so even structurally the block is detached from the real coordinate.

**Why this matters:**
As written this block adds no verification weight; if a prefactor degenerated on the branch (breaking invertibility and hence the `⟺`), these three checks would still pass. The genuine rigidity content is in fact already carried by the inverse-reconstruction round-trips (A3/A5/A6), which prove the map is invertible. So the block is misleading scaffolding that advertises the headline theorem while testing a triviality. The fix is to make it test the actual non-trivial content — that the three prefactors are nonzero on the constructive branch — so the equivalence is genuinely exercised and the block stops being a silent pass.

**Required change:**
Replace the tautological forward substitutions with a check that the three prefactors are nonzero on the branch (the real `⟹` content of rigidity). Add a small `expect_nonzero` helper and assert the numerators do not vanish on `χ₀>0, δ_U>0, 0<ε_η<1`. Concretely, verify `C_tr ≠ 0`, the Ξ₁→Σ_nt channel is recoverable (already covered by A5, so the rigidity block need only certify `C_tr` and the dressing prefactor), and `ε_η/(1−ε_η) ≠ 0`. See the directive for the exact edit. Cite: sympy:106–114, wl:79–85.

**Verification:**
After the fix, the sympy/mathematica transcripts' "Triple-rigidity theorem" section should print the nonzero prefactor values (e.g. `C_tr = ...`, `eps_eta/(1-eps_eta) = ...`) and a new `expect_nonzero` / nonzero-confirmation line for each, rather than three lines that all read `= 0` from substituting zero in. Both scripts must still exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is a close structural port of the `.py`: identical symbol names (`chi0, epsW, deltaU`, `sigmaZ…sigmaDel`, `sigmaTr, sigmaNT`), the identical `eps = epsW(1 − (2/11)δ_U/(1+δ_U))` definition (sympy:39 / wl:32), the same `theta1`/`xi1`/`r1` constructions in the same order, the same `sigmaNTDef`, the same `aTr`/`cTr`, and then the same eight `expectZero` checks in the same order with the same names, closing with a verbatim copy of the same carry-forward print block. This is textbook transliteration shape. I did **not** raise a `mathematica_transliteration` finding because this stage is purely an algebraic bookkeeping/consistency stage: there are no independent "physical premises" to re-derive from — `Θ₁`, `Ξ₁`, `R₁` are imported from Stage 250 as givens in both engines, so any honest second-engine check necessarily starts from the same imported expressions and runs the same simplifications. The legitimate value of the second engine here is that a different CAS independently confirms the same algebraic cancellations (e.g. that `Ξ₁ − (A_tr·Σ_tr + Σ_nt)` really reduces to 0, and that `A_tr/C_tr` really equals `2(1+χ₀+δ_U)/δ_U`), which it does. For a stage with derivable physical content the same mirroring would have been a finding; here it is acceptable. (F1 applies equally to both engines' rigidity blocks.)

## Engine cross-check

Both engines agree at the level claimed. The sympy output (lines 29–62) and the mathematica output (lines 18–51) report the same residuals: `Xi_1 - (A_tr Sigma_tr + Sigma_nt) = 0`, `R_1 + Xi_1 + eps_eta/(1-eps_eta) Sigma_eta = 0`, all four inverse residuals `= 0`, `A_tr/C_tr = 2(1+chi0+deltaU)/deltaU`, and the three (tautological) rigidity lines `= 0`, with identical printed `Σ_nt`, `Θ₁`, `Ξ₁`, `R₁+Ξ₁` forms (modulo CAS formatting). Both exit 0. `engines_agree: true`. Outputs are newer than their scripts (sympy out 2026-05-11T12:48 vs script 11:58; math out 13:23 vs script 11:56), so `outputs_fresh: true`.

## Verdict justification

The scripts faithfully and non-tautologically verify the paper's main claim — the triangular split of `Ξ₁`, the `R₁+Ξ₁ → Σ_η` collapse, all four inverse-reconstruction round-trips, and the `A_tr/C_tr` ratio identity — with constants and coordinate definitions matching the notes and appendix literally, and both engines agreeing. I tried to break the split check (A1) as a copy-of-itself but it does genuinely cross-check the SigmaTr-coefficient of `Ξ₁` against `A_tr`; I tried to find a constant mismatch against the paper (the `(11+9δ_U)/11`, `2χ₀/(1+δ_U)`, `4ε_Wδ_U/(11…)` coefficients) and found none. The one real defect is the "Triple-rigidity theorem" block (F1): it advertises the headline `⟺` equivalence but its assertions are tautological forward substitutions that cannot fail — the genuine rigidity content is carried instead by the inverse block. This is a low-severity strengthening, not a correctness error in the verified result, so the verdict is `findings` (not `stop_cold`) with `paper_alignment: aligned`. I confirm I read the paper card, the notes, and the appendix block before attacking the scripts, and the script's verified claim matches the paper.

## Self-test notes

Traps checked: (1) Variable independence / derivatives — none present; the stage is pure algebra, no `sp.diff`/`D[...]`, so the "zero-derivative" trap does not apply. (2) Parity/integrals — no integrals, N/A. (3) Trivial-case pre-check on my proposed F1 fix — substituting the branch numerators: `C_tr` numerator `χ₀δ_U > 0`, `A_tr` numerator `2χ₀ > 0`, dressing numerator `ε_η > 0` on `0<ε_η<1`, so the `expect_nonzero` checks return strictly nonzero literals and would fail only if a prefactor were identically zero — i.e. genuinely non-tautological. (4) Paths — F1 edits existing files (no missing-script path to name). (5) Paper round-trip — the prescribed fix tests the prefactors `C_tr`, `A_tr`, `ε_η/(1−ε_η)` exactly as the paper states them (notes Sec. 4 / appendix eqs. `app-part05-Ctr-Atr-defs`, `app-part05-triangular-normal-form`), introducing no new constant, so it cannot create a new paper_misalignment.
