---
unit_id: 159
batch: IV.6
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage159_hybrid_outlet_projection.md]
  paper_appendix: present
---

# Audit unit 159 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_159.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage159_hybrid_outlet_projection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1352 `\input{stages/stage_159}`; no separate row content beyond the included card)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage159_hybrid_outlet_projection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage159_hybrid_outlet_projection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage159_hybrid_outlet_projection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage159_hybrid_outlet_projection_mathematica_audit.txt`

## What the paper claims

Stage 159 projects the co-evolving Family-1 linear defect onto the compensated Robin–mixed (hybrid) outlet and reduces the surviving first-order obstruction. The card's bottom-line `\stagefield{Output}` is the boxed quote: "Even preservation forces `δg=0`, `δκ_W=0`, leaving `Δ_Q=-9σ_*δγ_W/(1-σ_*)`." The notes carry the full skeleton: (i) the off-compensation scalar `δC := δρ_R - 4δσ_W`; (ii) the exact first-order outlet formulas `δE2 = (δC - 9σ_*δκ_W)/(27(1-σ_*))`, `δE4 = (5δC - 72σ_*δκ_W)/(243(1-σ_*))`, `Δ_Q = (δC - 27σ_*δγ_W)/(3(1-σ_*))`; (iii) the mouth-gain transport `δρ_R=ΞδM_s`, `δσ_W=-ΞδM_q`, giving `δC = -4ΞΣ₀δR = -16σ_*δR` with the exact cancellation of `δΣ₀`; (iv) the canonical-even gate (`δE2=δE4=0`) collapsing to `δC=δκ_W=0` because the 2×2 determinant `-27σ_*≠0`; and (v) the final `Δ_Q = -9σ_*δγ_W/(1-σ_*)`. Distinct deliverables: the three `*_expected` outlet formulas (chi/E2/E4), the mouth-transport cancellation, the σ_* substitution, the canonical-even collapse, and the final reduced defect.

## What the script claims to verify

The SymPy docstring lists four checks matching the notes: (1) the linearized hybrid-outlet defects (δE2, δE4, Δ_Q) around the compensated canonical branch; (2) the mouth-gain → outlet-loading map with the exact `δΣ₀` cancellation in `δC = δρ_R - 4δσ_W`; (3) canonical-even preservation forcing `δC = δκ_W = 0`; (4) the collapse to pure `δγ_W` on that branch. Both scripts build `L0=-3+ρ-σ`, `L2=1/3-σκ`, `L4=1/9-σκ²`, `χ=3(1-9σγ)/(3-ρ+σ)`, `E2=-L2/L0-1/9`, `E4=L2²/L0²-L4/L0-4/81` (the notes' definitions verbatim), linearize them about the compensated branch `ρ→4σ₀+δρ, σ→σ₀+δσ, κ→1/3+δκ, γ→1/9+δγ`, and assert each linearized result equals the paper's closed form. They then verify the transport cancellation, the σ_* substitution, solve the 2×2 even-gate system, and substitute `δC=0` into the Δ_Q formula.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `Δ_Q = (δC - 27σ_*δγ_W)/(3(1-σ_*))` (notes §3 / via χ) | `expect_zero("delta chi formula", delta_chi - delta_chi_expected)` (py 78, wl 56) | match |
| `δE2 = (δC - 9σ_*δκ_W)/(27(1-σ_*))` (notes l.222) | `expect_zero("delta E2 formula", ...)` (py 79, wl 57) | match |
| `δE4 = (5δC - 72σ_*δκ_W)/(243(1-σ_*))` (notes l.230) | `expect_zero("delta E4 formula", ...)` (py 80, wl 58) | match |
| `δC = -4ΞΣ₀δR`, `δΣ₀` cancels (notes §2) | `expect_zero("deltaC mouth transport", ...)` (py 92, wl 71) | match |
| `δC = -16σ_*δR` (notes l.148) | `expect_zero("sigma_star substitution", ...)` (py 95-98, wl 75-78) | match |
| Even gate → `δC=δκ_W=0`, det `-27σ_*` (notes §5) | `sp.solve([eq1,eq2])`/`Solve[...,Reals]` + det print (py 103-111, wl 84-90) | match |
| `Δ_Q = -9σ_*δγ_W/(1-σ_*)` (card Output, notes §6) | `expect_zero("final Delta_Q + 9 sigma* dgamma /(1-sigma*)", ...)` (py 115, wl 94) | match |
| numerics `r_F1≈1.778`, `7.843`, `2.614` (notes §2/§4) | not exercised (scripts are symbolic only) | extra-on-notes (intermediate, not the card Output — see reconciliation) |

`paper_alignment: aligned` — every boxed/Output deliverable of the card and notes has a faithful, non-tautological script-side check. The unverified numerics are intermediate illustrations derived from an upstream (Stage 158) input `r_F1`, not Stage 159 bottom-line deliverables.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 78 | `simplify(delta_chi - delta_chi_expected)==0` | Δ_Q closed form (notes §3) | yes |
| A2 | sympy | 79 | `E2_lin - E2_expected == 0` | δE2 formula | yes |
| A3 | sympy | 80 | `E4_lin - E4_expected == 0` | δE4 formula | yes |
| A4 | sympy | 92 | `deltaC_expr - deltaC_expected == 0` | δΣ₀ cancellation | yes |
| A5 | sympy | 95-98 | `deltaC_expected.subs(...) + 16σ*dR == 0` | δC=-16σ_*δR | yes |
| A6 | sympy | 105-108 | `sp.solve(...)==[{deltaC:0,dk:0}]` | even-gate collapse | yes |
| A7 | sympy | 114-115 | `final_defect + 9σ0 dgamma/(1-σ0)==0` | final Δ_Q (card Output) | yes |
| B1-B7 | mathematica | 56,57,58,71,75-78,86-88,93-94 | `expectZero[...]` / `Solve` mirror of A1-A7 via independent linearizer/solver | same claims | yes |

All rows non-tautological: the linearized LHS is computed from the structural `E2/E4/χ` definitions, the RHS is the paper's independently-stated closed form; a wrong coefficient in either would make the residual nonzero.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` shares the paper-side *inputs* with the `.py` — the same `L0/L2/L4/χ/E2/E4` definitions, the same `subs` branch parametrization, and the same `*_expected` target forms. That sharing is expected and correct: those are the notes' own definitions and the paper's claimed formulas, which both engines must verify against. The *derivation machinery is genuinely independent*:

- Linearizer: `.py` (lines 53-66) uses `sp.series(expr, v, 0, 2).removeO()` per variable followed by a `Poly`-monomial total-degree≤1 filter; `.wl` (lines 40-45) uses an explicit first-order multivariate Taylor expansion `(f /. vars->0) + Σ_i (D[f, vars[[i]]] /. vars->0)·vars[[i]]`. Different algorithms — one truncates a series and filters monomials, the other evaluates the Jacobian at the base point.
- Simplification: `.py` `simplify(expand(...))` vs `.wl` `FullSimplify[Together[Expand[...]]]`.
- Even-gate solver: `.py` `sp.solve(...,dict=True)` vs `.wl` `Solve[..., Reals]`.

This is two independent confirmations of the same paper claim, not a line-by-line port. NOT a `mathematica_transliteration` finding.

## Engine cross-check

The two transcripts agree term-for-term. SymPy: `delta chi formula = 0`, `delta E2 formula = 0`, `delta E4 formula = 0`, `deltaC mouth transport = 0`, `sigma_star substitution = 0`, `solution = [{deltaC: 0, dk: 0}]`, `determinant = -27*sigma0`, `final Delta_Q ... = 0`. Mathematica: identical residuals plus `PASS:` flags, `solution = {{deltaC -> 0, dk -> 0}}`, `determinant = -27*sigma0`, and "Stage 159 Mathematica audit passed." No residual, sign, or factor disagreement. Outputs are fresh: `.txt` mtimes (2026-05-28 11:30 / 11:32) postdate both scripts (2026-05-27 23:11).

## Verdict justification

Clean. I attacked the three load-bearing assertion classes and each held: (1) the `E2/E4/χ` checks are not tautological — the LHS is a linearization of the structural ledger definitions and the RHS is the paper's separately-stated closed form, so a coefficient slip would surface as a nonzero residual; (2) the even-gate `sp.solve`/`Solve` genuinely tests non-degeneracy of the 2×2 system (the `-27σ_*` determinant is printed and the paper's `σ_*≠0` caveat is honored — `σ_*=0` would degenerate but the result correctly reflects the loaded branch); (3) the final-defect check substitutes the *derived* `δC=0` into the verified Δ_Q formula rather than re-asserting it. I checked the transliteration trap and confirmed the two linearizers use distinct algorithms. The recurring `168π²`/`100π²` and Family-1 `√(4107−100π²)/(10π)` constants do not appear in either script (confirmed by grep). Both `Inputs` (compensated branch `ρ_*=4σ_*, κ_*=1/3, γ_*=1/9`) and `Output` (`Δ_Q=-9σ_*δγ_W/(1-σ_*)`) match the card and notes exactly.

## Self-test notes

Variable-independence: each `sp.diff`/`D[]` in the `.wl` linearizer is taken w.r.t. a `δ`-variable that genuinely appears in `f=expr/.subs` (ρ,σ,κ,γ all map to base+δ), so no identically-zero derivative masks a trivial pass. Trivial-case: setting all δ→0 sends every `*_lin` and `*_expected` to its base value and the residuals to 0 with no false structure; setting e.g. `δγ≠0` makes `delta_chi` carry the `-27σ₀δγ/(3(1-σ₀))` term that must cancel against `delta_chi_expected`, so the assertion is live. Even-gate: the only solution is trivial precisely because det `=-27σ₀≠0`, matching the paper's σ_*≠0 branch — not a tautology. No parity/integral traps apply (no integrals in this unit).

## Value Reconciliation (pass-2 augmentation)

The scripts emit only symbolic residuals (all `= 0`) plus three labeled non-zero pieces that are byproducts of the solve/det, not free deliverables. Enumerating the genuine RESULT-level symbolic deliverables and locating each in the docs:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δC := δρ_R - 4δσ_W` | py 90/118, wl 69/98; out l.27/33 | notes l.132 / card-derived | MATCH |
| `δE2 = (δC - 9σ_*δκ_W)/(27(1-σ_*))` | py 75/119, wl 53/99; out l.28/35 | notes l.222 | MATCH |
| `δE4 = (5δC - 72σ_*δκ_W)/(243(1-σ_*))` | py 76/120, wl 54/100; out l.29/36 | notes l.230 | MATCH |
| `Δ_Q = (δC - 27σ_*δγ_W)/(3(1-σ_*))` | py 74/121, wl 52/101; out l.30/37 | notes l.238 | MATCH |
| `δC = -16σ_*δR` (σ_* substitution) | py 97, wl 77; out l.13/17 | notes l.148 | MATCH |
| even-gate determinant `-27σ_*` | py 110-111, wl 89-90; out l.19/24 | notes l.341 (`-27σ_*`) | MATCH |
| `δC=δκ_W=0` (even gate) | py 107, wl 88; out l.18/23 | notes l.347-349 / card l.16 | MATCH |
| `Δ_Q = -9σ_*δγ_W/(1-σ_*)` (final Output) | py 115/123, wl 94/103; out l.24/29-32 | card l.16 + notes l.394 | MATCH |

Numeric values in the notes that the scripts do NOT emit (so cannot mismatch): `r_F1≈1.77799…` (notes l.109), `√(1+r_F1²)≈2.0399…` (l.165), `δC≈7.8434…σ_*δg` (l.172), `Δ_Q≈2.6144…σ_*/(1-σ_*)δg` (l.306), `δΠ_tan` coefficients `0.8324…/1.1627…` (l.451-455). These are intermediate illustrations downstream of the Stage-158 input `r_F1`; they are NOT the card's `\stagefield{Output}` deliverable (which is the symbolic `Δ_Q=-9σ_*δγ_W/(1-σ_*)`, and IS verified). Per the augmentation guard ("a terse card legitimately omits intermediate quantities"; deliverable status keys off the boxed Output), these are INTERNAL/illustrative, not MISSING-DELIVERABLEs — no finding.

INTERNAL scaffolding (no prose expected): `L0/L2/L4/χ/E2/E4` structural definitions, the `subs` branch map, the `*_lin` linearized intermediates, pass/fail flags, all `= 0` residuals.

reconciliation: complete; 8 deliverable values checked, 0 misaligned
