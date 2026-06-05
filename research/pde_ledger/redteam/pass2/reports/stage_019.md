---
unit_id: 019
batch: I.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 019 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_019.tex`
- notes: `(none)` (no `notes/stages/moving_throat_pde_stage019_*.md` exists)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 60)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.txt`

## What the paper claims

Stage 019 specializes the parent throat-action bundle to the isotropic case and rewrites the one-pole and outgoing-normalization targets in parent-wall variables (`K_Σ`, `M_Σ`). It exports four boxed deliverables (Output line: "Stage~019 exports \eqref{eq:stage019-parent-bundle}--\eqref{eq:stage019-isotropic-compatibility}"):
(1) the isotropic bundle moments `D0=K_Σ-B0-Z0`, `D2=-(M_Σ+B2+Z2)`, `D4=-(B4+Z4)` (eq:stage019-parent-bundle);
(2) the one-pole surface `K_Σ = B0+Z0 + 3(M_Σ+B2+Z2)²/(B4+Z4)` (eq:stage019-one-pole-parent);
(3) the outgoing normalization `K_Σ = B0+Z0 + N0/P_{0,target}` with `P_{0,target}=54Gc_s⁵/(5a⁵c⁵m̂0²)` (eq:stage019-normalization-parent);
(4) the compatibility identity `N0/P_{0,target} = 3(M_Σ+B2+Z2)²/(B4+Z4)` (eq:stage019-isotropic-compatibility).
The card additionally states the audit "checks the response-sign criterion selecting the positive one-pole branch and the constant-prefactor conditions." There is no `notes/stages` file; the `.tex` card is the sole doc carrier. Appendix row 60 is a one-line status summary (`\StatusExactClosure{}`, "Isotropic parent-action specialization for the projected bundle").

## What the script claims to verify

The SymPy script (docstring lines 2–6) verifies "the symbolic bridge from the promoted parent wall block (KSigma, MSigma) to the isotropic grouped-P2 bundle moments and target surface." Concretely it: builds `D0/D2/D4` and the normalized pole moments `u2,u4`; proves the one-pole defect equals `(D0(B4+Z4)-3(MΣ+B2+Z2)²)/D0²` (line 39); solves `K_Σ` two ways (one-pole defect=0 and `P0=P0_target`) and matches each to its closed form (lines 49–52); checks the compatibility identity directly (line 56); derives the constant-prefactor `N2`/`N4` closures, a Jacobian-determinant identity, the `P2`/`P4` factorizations, and mutation guards that the factorizations are non-trivial (lines 58–101); factorizes the one-pole numerator in `M_Σ`, verifies Vieta sum/product and the `u2`-sign discrimination of the two roots, plus three numeric stable-pole samples (lines 103–173); and anchors `M_Σ`,`K_Σ` against concrete Gaussian wall integrals (`√π`, `3√π/2`, lines 176–185). The Mathematica script (M1–M12) re-derives `u2,u4,P0,P2,P4` independently via `Series` expansion of `1/den` and `D0(N0+N2 x+N4 x²)/den²` (lines 71–82) and checks the same identities.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| eq:stage019-parent-bundle (D0,D2,D4) | sympy 31–33, output 3–5; wl 46–48 | match |
| eq:stage019-one-pole-parent (K_Σ one-pole) | sympy 49,51 (`K_from_one_pole`), output 11; wl M2 93–98 | match |
| eq:stage019-normalization-parent (K_Σ norm + P0_target) | sympy 47,50,52, output 12–13; wl M3 101–105 + line 50 | match |
| eq:stage019-isotropic-compatibility | sympy 53–56, output 14; wl: implied by M2∧M3 (both K_Σ forms equal) | match |
| response-sign criterion (positive branch) | sympy 103–173 (M-roots, u2 signs, 3 numeric samples); wl M10–M11 | match |
| constant-prefactor conditions | sympy 58–101 (N2/N4 closures, Jacobian det, factorizations, mutation guards); wl M4–M9 | match |

All paper deliverables have faithful script-side coverage. No `mismatch`, `missing`, or `extra` rows. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `assert_zero(one_pole_defect - num/D0²)` | one-pole (eq 019-one-pole) prep | yes |
| A2 | sympy | 51 | `assert_zero(K_from_one_pole - closed)` | eq:stage019-one-pole-parent | yes |
| A3 | sympy | 52 | `assert_zero(K_from_norm - closed)` | eq:stage019-normalization-parent | yes |
| A4 | sympy | 56 | `assert_zero(compatibility - identity)` | eq:stage019-isotropic-compatibility | yes |
| A5 | sympy | 85 | `assert_zero(N2_const - closed)` | constant-prefactor cond. | yes |
| A6 | sympy | 86 | `assert_zero(N4_const - closed)` | constant-prefactor cond. | yes |
| A7 | sympy | 87 | `assert_zero(det - D0³)` | constant-prefactor cond. | yes |
| A8 | sympy | 88–89 | `assert_zero(P2/P4 factorization)` | constant-prefactor cond. | yes |
| A9 | sympy | 90–91 | `assert_zero(solve N2/N4 - closed)` | constant-prefactor cond. | yes |
| A10 | sympy | 92–99 | `assert_nonzero(mutated closures ≠ 0)` | guards A8 non-tautological | yes |
| A11 | sympy | 101 | `assert_zero(N4_const@one-pole - md form)` | constant-prefactor reported value | yes |
| A12 | sympy | 107–115 | `assert_zero(M-root factorization, Vieta sum/prod)` | response-sign criterion | yes |
| A13 | sympy | 118–120 | `assert_zero(u2 on ± roots) / assert_nonzero(roots differ)` | response-sign criterion | yes |
| A14 | sympy | 123–173 | numeric stable-pole sign guards (3 samples) | response-sign criterion | yes |
| A15 | sympy | 184–185 | `assert_zero(MΣ-√π), (KΣ-3√π/2)` | concrete anchor for K_Σ/M_Σ | yes |
| M1–M12 | mathematica | 84–198 | `expectZero[...]` / `expectNonzero[...]` mirroring A1–A15 via independent Series route | same claims | yes |

Every script-side check traces to a specific paper deliverable. No orphaned assertions.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. The SymPy script writes `u2,u4,P2,P4` as explicit hand-coded rational expressions (sympy lines 35–45). The Mathematica script instead obtains them from a Taylor `Series` expansion of the closed-form generating functions:
- wl 71–75: `den = D0 + D2*x + D4*x²; poleSeries = Normal[Series[1/den,{x,0,2}]]; normalizedPoleSeries = D0*poleSeries; u2 = Coefficient[...,x,1]; u4 = Coefficient[...,x,2]`.
- wl 77–82: `bundleSeries = Normal[Series[D0(N0+N2 x+N4 x²)/den²,{x,0,2}]]; P2 = Coefficient[...,x,1]; P4 = Coefficient[...,x,2]`.
This is a structurally different route (series-coefficient extraction vs. SymPy's algebraic ratio formulas), so the second-engine independence policy is satisfied. The Jacobian/determinant, mutation guards (M9), and Gaussian-integral anchors (M12) are likewise computed natively in Mathematica. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines reach residual 0 on every corresponding identity. SymPy output prints the closed forms (output lines 3–24) with no nonzero residual; the script raises on any nonzero `assert_zero`/`assert_nonzero` and instead reaches `STATUS: PASS` (output line 37). Mathematica output: `M1…M12 OK`, every `residual = 0`, and the two mutation guards correctly report nonzero residual `-(eps/(B0 - KSigma + Z0))` (mathematica output lines 19–20), confirming the factorizations are non-tautological. `STATUS: PASS` (line 32). The engines agree: e.g. `K_from_one_pole` (sympy output 11) and wl M2 both reduce the same one-pole closed form to 0; `N4 one-pole md form` reduces to `-5 N0 (B4+Z4)²/(9(B2+MΣ+Z2)²)` in sympy (output 20–21) and the wl M5/M11 chain is residual-0 consistent. No `engine_disagreement`.

## Verdict justification

`clean`. I read the paper card, the appendix row, confirmed no notes file exists, then attacked the scripts. Attacks tried that failed to break anything: (a) tautology hunt — the load-bearing factorization checks A8/M7–M8 are guarded by explicit `assert_nonzero` mutation guards (A10/M9) that confirm a perturbed closure `N2closed+eps` / `N4closed+eps` yields nonzero residual, so the factorizations genuinely pin the closures; (b) parity/derivative trap — the Jacobian `const_prefactor_matrix` (sympy 77–80) differentiates `D0²·P2` and `D0³·P4` w.r.t. `N2,N4`, both of which genuinely depend on `N2,N4`, and the determinant check `det - D0³` is non-trivially satisfied (the variables are present, so no identically-zero-derivative pass-through bug); (c) symbol-domain check — `KSigma,MSigma,B*,Z*,N*` are `real, nonzero`, physical scales `mhat0,G,cs,a,c` are `positive`, matching the parent-wall setup; the `u2`-sign discrimination is validated numerically over three samples with `D0>0, B4+Z4>0` enforced (sympy 165–172), so the positive-branch selection is exercised, not assumed; (d) constant check — `P0_target = 54 G cs⁵/(5 a⁵ c⁵ mhat0²)` matches the paper card exactly (sympy 47 vs tex line 32), and the Gaussian wall integrals `√π`/`3√π/2` are independent absolute anchors. The script's verified claims match the paper's four boxed deliverables plus the two auxiliary criteria the card names. Both outputs are fresh (output mtimes newer than their scripts).

## Self-test notes

Checked the four directive traps even though no directive is written: (1) variable independence — the only derivatives are `sp.diff(P2zeroEq/P4zeroEq, N2/N4)` (sympy 78–79) and both polynomials genuinely contain `N2,N4`, so no identically-zero derivative; (2) symmetry/parity — the sole symmetric-domain integrals are the Gaussian wall integrals, both even-times-even integrands giving the expected nonzero `√π`/`3√π/2`; (3) trivial-case — the numeric `baseline` sample gives `u2+=0.6325>0, u2-=-0.6325<0` as asserted (output 26–27); (4) the mutation guards confirm the `assert_zero` factorizations are not vacuously true. All consistent; no finding.

## Value Reconciliation (pass-2 augmentation)

Scope note: no `notes/stages/moving_throat_pde_stage019_*.md` exists, so the `.tex` card is the sole prose carrier; reconciliation is script→`.tex`. Both saved outputs are fresh (output mtimes > script mtimes), so the committed `.txt` is the authoritative emitted record.

Deliverable-level values (the stage's boxed/Output quantities):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `D0 = K_Σ-B0-Z0`, `D2 = -(M_Σ+B2+Z2)`, `D4 = -(B4+Z4)` | py 31–33; sympy out 3–5; wl 46–48 | `stage_019.tex:17-19` (eq:stage019-parent-bundle) | MATCH |
| `K_Σ (one-pole) = B0+Z0 + 3(M_Σ+B2+Z2)²/(B4+Z4)` | py 49,51; sympy out 11; wl M2 93–98 | `stage_019.tex:24-26` (eq:stage019-one-pole-parent) | MATCH |
| `K_Σ (norm) = B0+Z0 + N0/P_{0,target}` | py 50,52; sympy out 12; wl M3 101–105 | `stage_019.tex:30` (eq:stage019-normalization-parent) | MATCH |
| `P_{0,target} = 54 G c_s⁵/(5 a⁵ c⁵ m̂0²)` | py 47; sympy out 13; wl 50 | `stage_019.tex:32` | MATCH |
| compatibility `N0/P_{0,target} = 3(M_Σ+B2+Z2)²/(B4+Z4)` | py 53–56; sympy out 14; wl M2∧M3 | `stage_019.tex:39-40` (eq:stage019-isotropic-compatibility) | MATCH |

reconciliation: complete; 5 deliverable values checked, 0 misaligned.

INTERNAL (genuine scaffolding the terse card legitimately omits; named collectively by the card as "constant-prefactor conditions" / "response-sign criterion"; not boxed deliverables, no notes file to carry them — per augmentation guards these are MATCH-by-omission, no finding):
`N2_on_constant_prefactor_branch` (2N0(B2+MΣ+Z2)/(B0-KΣ+Z0)); `N4_on_constant_prefactor_branch`; `N4_on_one_pole_plus_constant_prefactor = -5N0(B4+Z4)²/(9(B2+MΣ+Z2)²)`; `N4_md_equivalent_on_one_pole` (same value, equivalence check); `u2`/`u4` symbolic forms; `u2_on_positive_root`/`u2_on_negative_root`; Jacobian determinant `= D0³`; mutation-guard residual `-eps/D0`; numeric stable-pole samples (u2 ±0.6325, ±1.8257, ±0.0577); concrete wall-integral anchors `M_Σ=√π`, `K_Σ=3√π/2`; pass/fail flags, residual-near-zero values, sample counts.
