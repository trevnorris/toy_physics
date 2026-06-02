---
unit_id: 214
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md"]
  paper_appendix: present
---

# Audit unit 214 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_214.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows 57/59/61, narrative 236, detailed subsection 1058-1080, `\input` at 1351)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 214 (`\stagefield{Purpose}`: "solves the three-parameter interior optimization problem for a primitive quadruple") reduces the full interior optimizer on a four-coordinate simplex to a finite algebraic candidate set. The `\stagefield{Output}` reads: "Finite quadruple-interior candidate set with preferred bound \(54\) per envelope, projected fallback bound, and interior winner/non-improvement theorem against the Stage~213 boundary ledger." The `\stagefield{Derivation ledger}` enumerates the deliverables: (1) three stationary numerator equations, (2) one auxiliary square-root variable lift, (3) the lifted degree pattern \((3,3,3,2)\), (4) the proven \(54\)-candidate preferred bound, and (5) a square-root-free projected fallback. The notes add three more: two special reductions (diagonal-isotropic curvature → gradient-optimal ray; full-symmetry → equal-mix barycenter \((1,1,1,1)/2\) is stationary) and the interior winner/non-improvement theorems. The appendix (Eq. `eq:app-part06-four-degree-pattern`/`four-bezout`) confirms the degree pattern \((3,3,3,2)\) and the preferred bound \(3\cdot3\cdot3\cdot2=54\) per envelope; it does not quote a numeric value for the projected fallback bound. The notes (section 4.3, line 435) quote the projected one-chart bound as `\boxed{5\cdot 5\cdot 6 = 218.}`.

## What the script claims to verify

The SymPy script builds the discriminant numerator `Delta` symbolically in coefficients `A..J`, forms `Phi`/`tau`, derives the three slope numerators `M_r,M_s,M_t` and discriminant-transport numerators `L_r,L_s,L_t`, and verifies the three exact derivative laws `dPhi/d{r,s,t} == N_{r,s,t}/(2(1+r^2+s^2+t^2)^{3/2} sqrt(Delta))` (section I). It forms the lifted polynomial system `F_r,F_s,F_t,F_Delta`, checks their total degrees are `(3,3,3,2)`, and asserts the lifted Bézout bound `3*3*3*2 == 54` (section II). It defines the quintic cross-consistency polynomials `C_rs,C_rt,C_st` and the sextic square conditions `S_r,S_s,S_t`, checks six elimination identities, checks the degrees `(5,5,5)`/`(6,6,6)`, and asserts the projected one-chart bound `5*5*6 == 150` (section III). It verifies the diagonal-isotropic reduction `Delta_iso` and that `tau_iso` depends only on the slope magnitude, plus the gradient-optimal-ray normalization (section IV), the equal-mix stationarity `N_{r,s,t}(1,1,1)==0` and barycenter normalization (section V), and a combinatorial integer model of the winner/non-improvement ordering (section VI).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Three stationary numerator equations / exact derivative laws | section I `expect_zero("exact dPhi/d{r,s,t} compiler")` lines 106-108 | match |
| Lifted system + degree pattern `(3,3,3,2)` | section II degree asserts lines 139-146 | match |
| Preferred bound `54` per envelope | `bezout_lift != 54` line 150-151 | match |
| Square-root-free projected fallback (degrees: quintic cross, sextic square) | section III degree asserts lines 201-204 | match |
| Projected fallback **numeric bound** (notes: `5*5*6 = 218`) | `bezout_projected != 150` lines 206-209 (script uses 5*5*6=150) | mismatch — see F1 |
| Diagonal-isotropic → gradient-optimal ray | section IV lines 233-247 | match |
| Full-symmetry → equal-mix barycenter stationary | section V lines 274-283 | match |
| Interior winner / non-improvement theorems vs Stage 213 ledger | section VI lines 292-316 (integer ordering model) | partial — exercises strict-inequality transitivity on integer samples, not the actual `tau`/`beta` brackets |
| Mathematica independent verification (dual-engine) | none | missing — see F2 |
| Section III elimination identities (six) | lines 179-185 | extra/tautological — zero by algebraic construction; see F3 |

Dominant pattern: the substantive physics deliverables (derivative laws, lifted degree ledger, the `54` bound, the two reductions, equal-mix stationarity) are faithfully and non-tautologically exercised, but (a) the notes' projected bound `218` contradicts the script's `150`, (b) the second engine is absent, and (c) section III's six identity checks are definitional tautologies. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 106 | `simplify(dPhi/dr - N_r/(...)) == 0` | derivative law / stationary numerator | yes |
| A2 | sympy | 107 | `simplify(dPhi/ds - N_s/(...)) == 0` | derivative law | yes |
| A3 | sympy | 108 | `simplify(dPhi/dt - N_t/(...)) == 0` | derivative law | yes |
| A4 | sympy | 139-146 | `deg F_{r,s,t}==3`, `deg F_Delta==2` | degree pattern (3,3,3,2) | yes |
| A5 | sympy | 150 | `bezout_lift == 54` | preferred bound 54 | yes |
| A6 | sympy | 179-181 | `Ms*Fr - Mr*Fs - Crs == 0` (×3) | projected elimination (definition) | no — tautological by construction |
| A7 | sympy | 183-185 | `Fr*(Fr-4Mr*y)+4Mr^2*FDelta - Sr == 0` (×3) | projected elimination (definition) | no — tautological by construction |
| A8 | sympy | 201-204 | `deg C==5`, `deg S==6` | square-root-free fallback degrees | yes |
| A9 | sympy | 206-209 | `bezout_projected == 150` | projected fallback numeric bound | yes (script) / conflicts notes (218) |
| A10 | sympy | 233-240 | `Delta_iso` reduction; `tau_iso` depends only on slope | diagonal-isotropic reduction | yes |
| A11 | sympy | 246-247 | gradient-ray normalization + slope value | gradient-optimal interior ray | yes |
| A12 | sympy | 278-283 | `N_{r,s,t}(1,1,1)==0`, equal-mix norm | equal-mix barycenter stationary | yes |
| A13 | sympy | 293-316 | integer ordering of `istar<bstar` / `bstar<istar` | winner / non-improvement theorems | partial — integer transitivity model |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction.md:435` quote: `\boxed{5\cdot 5\cdot 6 = 218.}`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:206-209` quote: `bezout_projected = deg_Crs * deg_Crt * deg_Sr` … `if bezout_projected != 150: raise AssertionError("Unexpected projected Bezout bound")`

**What's wrong:**
The notes state the projected one-chart Bézout bound is `5·5·6 = 218`. The script computes the same product of degrees (`deg_Crs=5`, `deg_Crt=5`, `deg_Sr=6`) and asserts the result is `150`. The two disagree. `5·5·6 = 150` is the arithmetically correct product; `218` is not the product of any of the cited degree combinations (`5·5·6=150`, `6·6·6=216`, etc. — none yield 218). So this is almost certainly an arithmetic typo on the notes side, but it is a prose document and the red-team cannot edit it. The paper card and the Part VI appendix do not quote a numeric value for the projected fallback bound (the card says only "projected fallback bound"; the appendix quotes only the preferred `54`), so the card/appendix do not themselves carry the error — the conflict is strictly notes ↔ script.

**Why this matters:**
A reader cross-referencing the notes against the verified script will see a 218-vs-150 discrepancy in a load-bearing combinatorial bound. Left unresolved, either the notes mis-state the projected fallback ceiling or (less likely) the script's degree assignment is wrong. The script's degrees (quintic cross-consistency, sextic square) match both the notes' own prose ("Each is quintic", "Each is sextic") and the appendix, so the script's `150` is the internally consistent value.

**Required change:**
None applied by Codex. Route to user — see directive `## Resolve before fix_loop`.

**Verification:**
After the user picks a direction: if notes are corrected to `150`, no script change and the assertion at line 208 already matches; if the user asserts a different degree combination is intended, the script's degree asserts (lines 201-204) and the `150` literal (line 208) get updated accordingly and the script re-run.

### F2 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `mathematica/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl` (to be created)

**What's wrong:**
Stage 214 is non-status-only (`is_status_only_candidate: False`) and non-checkpoint, so the project dual-engine contract requires both a SymPy and a Mathematica audit. Only the SymPy script exists; `\stagefield{Verification}` itself records "Mathematica audit: none yet." Every load-bearing claim in this stage is independently verifiable in Mathematica using native primitives (`D[]`, `Series`/`Coefficient` or `Exponent`/`PolynomialReduce` for total degree, `Resultant`/`GroebnerBasis` for elimination, `Solve`/`Reduce` and substitution for the reductions). There is no mathematical obstruction, so the gap is a finding, not an exemption.

**Why this matters:**
A single-engine stage has no cross-check; a transcription or algebra error in the SymPy derivative laws or degree ledger would go undetected. The contract requires the second engine wherever it is possible, which it is here.

**Required change:**
Create the `.wl` per the claim manifest in the directive (F2). It must derive the results independently via native Mathematica primitives and a different decomposition than the `.py`; a line-by-line port is rejected.

**Verification:**
Verifier runs `math -script mathematica/...stage214..._mathematica_audit.wl`, confirms it exits 0 with all manifest claims (M1-M7) checked, and confirms it is not a transliteration of the `.py`.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:179-185`

**What's wrong:**
The six section-III "elimination identity" checks are algebraically guaranteed to be zero by the way the polynomials are constructed, so they cannot fail regardless of the physics. For the cross-elimination (lines 179-181): with `Crs := Ms*Lr - Mr*Ls` (line 158), `Fr := 2*Mr*y + Lr`, `Fs := 2*Ms*y + Ls`, the tested expression `Ms*Fr - Mr*Fs - Crs = (2*Mr*Ms*y + Ms*Lr) - (2*Mr*Ms*y + Mr*Ls) - (Ms*Lr - Mr*Ls) ≡ 0` — the `y` terms cancel and the remainder is the definition of `Crs`. For the square-elimination (lines 183-185): with `Sr := Lr^2 - 4*Mr^2*Delta` (line 162) and `FDelta := y^2 - Delta`, the tested expression `Fr*(Fr - 4*Mr*y) + 4*Mr^2*FDelta - Sr = (Lr^2 - 4*Mr^2*y^2) + (4*Mr^2*y^2 - 4*Mr^2*Delta) - (Lr^2 - 4*Mr^2*Delta) ≡ 0`. Both are pure rearrangements of the definitions. The substantive content of section III (that `C` is quintic and `S` is sextic) is separately and validly checked at lines 201-204; only these six identity asserts are tautological.

**Why this matters:**
These six asserts give a false sense of verification: they confirm bookkeeping, not that the eliminants actually capture the stationary conditions. They would still pass if `M`, `L`, `Delta` were mis-defined upstream. Low severity because the load-bearing degree checks are independent and the derivative laws (section I) anchor the upstream `M`/`L`/`Delta` definitions non-tautologically.

**Required change:**
Replace the six definitional identity asserts with a non-tautological check that the eliminants actually vanish on the stationary manifold. Concretely: pick a fixed rational, non-degenerate numeric assignment for `A..J, k_i,k_j,k_k,k_l, H0`, solve the lifted stationary system (or pick a parameter family for which a stationary point is known), substitute a genuine stationary `(r,s,t,y)` into `C_rs, C_rt, C_st, S_r, S_s, S_t` and assert each evaluates to `0` numerically — and additionally assert each eliminant is NOT identically zero as a polynomial (so the vanishing is meaningful, not vacuous). See directive F3 for the acceptance criteria.

**Verification:**
The six `expect_zero` identity lines (179-185) are replaced by checks that (i) each eliminant is a nonzero polynomial and (ii) each vanishes at a constructed stationary point; the script still exits 0.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot occur yet. The directive's claim manifest carries an explicit anti-transliteration guard for the new script.

## Engine cross-check

Only one engine is present; `engines_agree: n/a`. This is the basis for F2.

## Verdict justification

The SymPy script faithfully and non-tautologically verifies the stage's central physics deliverables: the three exact derivative laws (section I attacks the stationary-numerator factorization directly and would catch a sign or factor error in `M`/`L`/`N`), the lifted `(3,3,3,2)` degree ledger and the `54` bound (matching both the card `\stagefield{Output}` and the appendix), the diagonal-isotropic and full-symmetry reductions, and the equal-mix barycenter stationarity. I attacked the derivative laws (rebuilt the `N/(...)` denominator and confirmed it is the genuine quotient, not a contrived match), the symbol domains (`r,s,t,y` nonnegative, `k`'s positive, `A..J` unrestricted — all correct; the `I` symbol shadows `sympy.I` harmlessly since the imaginary unit is unused), and the degree assignments (consistent with the appendix). Three issues survive: F1 (notes say projected bound `218`, script asserts `150`; needs user resolution since Codex cannot edit notes), F2 (no Mathematica engine, required by the dual-engine contract and fully possible here), and F3 (the six section-III elimination identities are zero by construction). Verdict `findings`, no stop-cold: the math is sound; F1 is a notes typo pending user resolution, F2/F3 are routine script-side repairs. Read paper card + notes + appendix before the script; the script's verified claim matches the paper on every deliverable except the notes' arithmetic typo for the projected bound.

## Self-test notes

Checked the variable-independence trap: the section-I checks differentiate `Phi`/`Delta` w.r.t. `r,s,t`, and `Delta` genuinely depends on `r,s,t` through `B..J` (lines 56-67), so no identically-zero-derivative trap. Checked the tautology trap and confirmed the section-III identities (F3) collapse to `0` by definition while the degree checks at 201-204 are genuinely substantive. Checked the F3 required-change for a new trap: the proposed "vanishes at a stationary point AND eliminant not identically zero" pair avoids a vacuous pass (substituting into an identically-zero polynomial), and the directive states the requirement (vanish on the manifold; nonzero as a polynomial) without prescribing the solve route. Re-checked the paper round-trip for F2/F3: neither fix introduces a new constant; the F2 manifest reproduces the same `54`/`(3,3,3,2)` the card and appendix state, and uses `150` only as the script's own consistent projected value (not `218`), so it does not bake in the F1 dispute.
