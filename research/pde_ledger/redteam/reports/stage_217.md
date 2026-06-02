---
unit_id: 217
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: misaligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md]
  paper_appendix: present
---

# Audit unit 217 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_217.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows/subsections at lines 65, 1155-1207)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card (`stage_217.tex`) states the `Output` verbatim: *"Finite support-five interior candidate set, preferred budget \(324\), fallback budget \(1500\), and support-five interior improvement/non-improvement theorem."* The card's derivation-ledger line claims the lifted stationary system has degree pattern \((3,3,3,3,2)\) and "proves the preferred \(179\)-candidate bound per envelope." The Part VI appendix (`stage_appendix_part06.tex`) restates the chart, the certified root function \(\tau_\star\), the lifted polynomial system, the degree pattern \((3,3,3,3,2)\), and writes the per-envelope bound as `3^4\cdot2=179` (eq. line 1200), the cross-envelope preferred budget as `2\times162=324` (eq. line 1205), and the fallback projected-chart bound as \(750\) per envelope / \(1500\) total. The appendix row (line 65) likewise says "preferred bound \(179\) per envelope." The notes file states the same scaffolding but writes the per-envelope lifted Bézout bound as `3\cdot 3\cdot 3\cdot 3\cdot 2 = 230` (lines 5, 32, 409, 616). Distinct deliverables: (D1) exact five-face reduction / face count = 5; (D2) exact interior stationary-numerator identities for \(r,s,t,u\); (D3) lifted degree pattern \((3,3,3,3,2)\) and the per-envelope lifted Bézout candidate bound; (D4) projected square-root-free degree pattern \((5,5,5,6)\) and projected bound \(750\); (D5) gradient-optimal interior ray exact under diagonal-isotropic curvature; (D6) equal-mix barycenter exact under full fivefold symmetry; (D7) the cross-envelope budgets 324 / 1500 and the improvement/non-improvement theorem.

## What the script claims to verify

The SymPy script verifies: (1) the primitive five-coordinate simplex has exactly five codimension-one quadruple faces (line 54); (2) the four exact stationary-numerator identities — that \(2\sqrt{\Delta}\,S^{3/2}\,\partial_\bullet\Phi_\star\) equals \(2M_\bullet\sqrt{\Delta}+L_\bullet\) for \(\bullet\in\{r,s,t,u\}\) (lines 92-107); (3) the lifted polynomial degrees are \((3,3,3,3,2)\) and the lifted Bézout product equals **162** (lines 133-139); (4) the projected square-root-free degrees are \((5,5,5,6)\) with product **750** (lines 162-168); (5) under diagonal-isotropic curvature, \(L_\bullet(\mathrm{diag})=2K_{\rm lin}M_\bullet\) and both \(M_\bullet\) and \(L_\bullet\) vanish at the gradient-optimal ratios \(r=k_c/k_\lambda,\dots\) (lines 211-223); (6) under full fivefold symmetry the equal-mix barycenter \(r=s=t=u=1\) makes \(M_\bullet=L_\bullet=0\) (lines 256-263). All checks pass (output exit 0). The script's verdict block self-reports "Bezout bound 162" (line 268, output line 86).

## Paper ↔ script cross-check

| Deliverable | Script check | Status |
|---|---|---|
| D1 face count = 5 | line 54 `len(faces) != 5` | match |
| D2 stationary-numerator identities (r,s,t,u) | lines 92-107 | match |
| D3 lifted pattern (3,3,3,3,2) | lines 133-134 | match |
| D3 per-envelope lifted Bézout bound | line 138 asserts **162**; paper card/appendix say **179**; notes say **230** | **mismatch** |
| D4 projected pattern (5,5,5,6) + bound 750 | lines 162-168 | match |
| D5 gradient-optimal ray exact (diag-isotropic) | lines 211-223 | match |
| D6 equal-mix barycenter exact (fivefold sym) | lines 256-263 | match |
| D7 cross-envelope budget 324 | not checked in script (script computes per-envelope only); appendix derives it as `2×162=324`, which is consistent with the script's 162, not with 179 or 230 | partial (no script check, but arithmetic only closes for 162) |
| D7 fallback budget 1500 | not checked (750 per-envelope is checked; 2×750=1500 consistent) | partial |
| D7 improvement/non-improvement theorem | not checked (status/inequality logic, carry-forward) | missing (acceptable: pure interval-comparison statement, no new computable identity beyond the candidate-set machinery) |

`paper_alignment: misaligned` — the dominant defect is the per-envelope-bound value, which the paper states two different wrong ways (179 in card+appendix, 230 in notes) versus the script's mathematically-correct 162.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `if len(faces) != 5: raise` | D1 | yes |
| A2 | sympy | 92-95 | `expect_zero(2√Δ·S^{3/2}∂_rΦ − (2M_r√Δ+L_r))` | D2 | yes (non-tautological: M_r, L_r defined independently of Φ) |
| A3 | sympy | 96-99 | s-stationary identity | D2 | yes |
| A4 | sympy | 100-103 | t-stationary identity | D2 | yes |
| A5 | sympy | 104-107 | u-stationary identity | D2 | yes |
| A6 | sympy | 133-134 | `(deg…) != (3,3,3,3,2): raise` | D3 | yes |
| A7 | sympy | 138-139 | `bezout_lift != 162: raise` | D3 (value) | yes — but paper says 179/230 → see F1 |
| A8 | sympy | 162-163 | `(deg…) != (5,5,5,6): raise` | D4 | yes |
| A9 | sympy | 167-168 | `bezout_proj != 750: raise` | D4 | yes |
| A10 | sympy | 211-214 | `L_•(diag) − 2K_lin M_• == 0` | D5 | yes |
| A11 | sympy | 216-223 | `M_•, L_• at grad ratios == 0` | D5 | yes |
| A12 | sympy | 256-263 | `M_•, L_• at barycenter == 0` | D6 | yes |
| — | mathematica | — | (absent) | all | missing — see F2 |

All SymPy assertions are non-tautological and well-anchored within the script's own scope. The defect is not internal: A7 verifies the correct arithmetic, but it disagrees with the paper-side number.

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** value_mismatch
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_217.tex:13` — "proves the preferred \(179\)-candidate bound per envelope"
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex:1200` — `3^4\cdot2=179`
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex:65` — "preferred bound \(179\) per envelope"
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction.md:409` — `3\cdot 3\cdot 3\cdot 3\cdot 2 = 230` (also lines 5, 32, 616)
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage217..._sympy_audit.py:138` — `if bezout_lift != 162: raise AssertionError`

**What's wrong:**
The per-envelope lifted Bézout bound is the product of the degree pattern \((3,3,3,3,2)\): \(3\cdot3\cdot3\cdot3\cdot2 = 3^4\cdot2 = 81\cdot2 = 162\). The script verifies exactly this (`bezout_lift != 162` raises; output line 42/86 prints "162"). The paper card and the Part VI appendix instead assert this product equals **179** (`3^4\cdot2=179`, which is arithmetically false). The notes assert it equals **230** (also false). Three sources, three values, for one product. The script's 162 is the only correct one. Moreover the appendix is internally self-contradicting: immediately after writing `3^4\cdot2=179` it writes the cross-envelope preferred budget as `2\times162=324` (line 1205) — i.e. it uses 162 as the per-envelope value in the very next equation, and \(2\times179=358\neq324\) and \(2\times230=460\neq324\). The total budget 324 (carried into the card `Output` and into the final `1140+324=1464`) closes only if the per-envelope value is 162.

**Why this matters:**
The headline deliverable of the stage's `Output` line — preferred budget 324 — is arithmetically consistent only with the per-envelope value 162 that the script verifies, while the prose number quoted alongside it (179, and 230 in the notes) is wrong and self-inconsistent. Left alone, the paper publishes a false arithmetic identity (`3^4·2=179`) as a "preferred candidate bound," and the notes carry yet a third wrong value. This is a presentation/value error, not a math error in the script; the script is right.

**Required change:**
None for Codex. This is a paper-side value disagreement and must be routed to the user. See the directive's `## Resolve before fix_loop` block. Codex must NOT edit the paper card, the appendix, or the notes, and must NOT change the script's `162` to match the paper.

**Verification:**
After the user chooses a direction, the verifier confirms either (a) the paper/notes lines now read 162 (matching script + the already-correct `2×162=324`), with no script change; or (b) — not expected, since 179/230 are arithmetically impossible for \(3^4\cdot2\) — a re-derivation justifying a different count, in which case both the degree-pattern claim and the `324` budget would also need revisiting.

### F2 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` — no `.wl` for stage 217
- target: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`

**What's wrong:**
Stage 217 is non-status-only (`is_status_only_candidate: False`) and non-checkpoint (`is_checkpoint: False`). Per the dual-engine project contract, a second-engine `.wl` is required wherever Mathematica CAN independently verify the stage. Mathematica can independently verify every computable deliverable here: the symbolic stationary-numerator identities (via `D[]`, `Together`, `Simplify`), the lifted/projected polynomial total degrees and their Bézout products (via `Exponent`/`MonomialList`/manual total-degree extraction), the diagonal-isotropic gradient-optimal reduction, and the fivefold-symmetric barycenter reduction. None of these is impossible or even hard in Mathematica. No `.wl` exists, so the second engine is entirely absent.

**Why this matters:**
The lifted candidate-set construction and the two special reductions are the load-bearing results of the stage and currently rest on a single engine. Single-engine verification cannot catch a SymPy-specific simplification artifact (e.g. a `simplify` masking a branch issue in the \(\sqrt{\Delta}\) numerator identity). The contract requires independent corroboration where it is possible, and here it is possible.

**Required change:**
Codex must DESIGN and WRITE a new Mathematica `.wl` (path above) that independently re-derives and asserts the claim manifest M1-M6 below, using native Mathematica primitives via a different decomposition than the SymPy script. See the directive for the manifest and the anti-transliteration guard.

**Verification:**
After Codex applies, the verifier runs `redteam exec-mathematica 217` and confirms the new `.wl` exists, asserts M1-M6, and exits 0; and that its degree/Bézout results agree with the SymPy script's (notably the per-envelope value 162, not the paper's 179).

## Independent-derivation check (Mathematica)

No `.wl` exists, so no transliteration check applies yet. The directive (F2) explicitly forbids a line-by-line port: the new `.wl` must derive \(\Phi_\star\) and its partials directly with `D[]` and reduce, must extract total degrees with native polynomial primitives rather than re-typing the SymPy `Poly(...).total_degree()` choreography, and must use a different substitution route for the two special reductions.

## Engine cross-check

Only one engine present; no cross-check possible. This is itself finding F2. When the `.wl` lands, the two engines must agree on: the four stationary identities reducing to 0, the lifted degree tuple `(3,3,3,3,2)` and product **162**, the projected tuple `(5,5,5,6)` and product **750**, and the vanishing of the special-reduction residuals.

## Verdict justification

The SymPy script is internally clean and adversarially sound: every assertion is non-tautological and anchored to a real deliverable; the stationary-numerator identities genuinely exercise the \(\Phi_\star\) construction (M_r/L_r are defined independently of \(\Phi\), so the identity can fail); the degree and Bézout checks compute the polynomials honestly and check the correct arithmetic (162, 750). The two defects are external to the script's correctness: (F1) the paper card, appendix, and notes quote a per-envelope bound (179 / 230) that is arithmetically wrong and inconsistent with the `324` total they themselves carry — the script's 162 is correct, so this routes to the user as a paper-side value fix; and (F2) the required second engine is absent though Mathematica can independently verify the stage. Verdict: `findings`, no stop-cold (the math is sound; F1 is a paper-side number to be corrected, F2 is a fillable gap).

## Self-test notes

Variable-independence trap: `D[Phi,r]` is nonzero (Phi depends on r via K_lin, Δ, and S), and `D[Δ,r]` is nonzero (Δ carries B r + F r² + G rs + H rt + I ru), so the manifest's stationary identity is a genuine identity, not a 0≡0 from a vanishing derivative. Parity trap: no unbounded integrals in this stage, so N/A. Trivial-case/value trap: the manifest asserts the correct arithmetic `3^4·2 = 162` (matching the .py and the appendix's own `2×162=324`), NOT the paper's 179 or the notes' 230 — so the new .wl will not bake in the disputed value and will correctly corroborate the script while the F1 value goes to the user. Path trap: target `.wl` lives under `mathematica/`, named per convention.
