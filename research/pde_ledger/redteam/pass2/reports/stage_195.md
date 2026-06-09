---
unit_id: 195
batch: V.3
auditor_model: Opus 4.8 (1M context)
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md]
  paper_appendix: present
---

# Audit unit 195 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_195.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 121, 1306-1317, 1467)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_mathematica_audit.txt`

## What the paper claims

Stage 195 performs the exact source-map reduction of the canonical outgoing branch. `\stagefield{Output}`: "Factorizes the observable odd closure as \(\widehat m_0^{2}\chi_QN_Q=1\) and collapses \(\Delta_{\rm norm}\)." The notes enumerate five deliverables: (1) the exact odd ratio `Γ̄₅/Γ̄₅^target = χ_Q N_Q` (with `N_Q := P̄₀/P0^target`, `P0^target = 54Gc_s⁵/5a⁵c⁵`, `Γ̄₅^target = 2G/5c⁵`); (2) the boxed observable odd-closure factorization `m̂₀²χ_Q N_Q = 1` obtained from the observable condition `m̂₀²Γ̄₅ = Γ̄₅^target`; (3) the collapse `Δ_norm = P0^target(1/χ_Q − 1)` (equivalently `−P0^target·Δ_Q/(1+Δ_Q)` with `Δ_Q := χ_Q−1`); (4) the natural point-particle reduction `m̂₀→1 ⟹ N_Q = 1/χ_Q`, and `N_Q−1 = −(χ_Q−1)+O(²)`; (5) the source-map-reduced DtN deformation algebra `χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)`, its `N_Q`/`Δ_norm` images, the linearization `N_Q−1 = −5ε_β − δΣ₀/3S − 9δΣ₅/S + O(2)`, and the canonical-branch closure `β=1,Σ₀=Σ₅=0 ⟹ χ_Q=1, N_Q=1, Δ_norm=0`. The appendix (eqs. app-part05-source-outgoing-factorization, app-part05-natural-source-map) restates `m̂₀²χ_QN_Q=1` and `N_Q=χ_Q^{-1}` verbatim.

## What the script claims to verify

Both scripts verify all five deliverables in five sections (I–V). Section I builds the carry-forward tuple and checks the two `Γ̄₅` forms agree and that `Γ̄₅/Γ̄₅^target = χ_Q N_Q` after `P0 → N_Q·P0^target`. Section II derives the odd-closure factorization from the *observable* condition residual `(m̂₀²Γ̄₅ − Γ̄₅^target)` and checks it equals `Γ̄₅^target·(m̂₀²χ_Q N_Q − 1)`, then collapses `Δ_norm` to `P0^target(1/χ_Q − 1)` and to the `Δ_Q` form. Section III takes the point-particle limit. Section IV inserts the Stage-194 deformation algebra and checks `N_Q`, `Δ_norm^(pt)`, and their linearizations. Section V verifies the canonical branch closes. Crucially, in the `.wl` Section IV the deformation `χ_Q` is *re-derived* from the genuine `SphericalHankelH1[2,x]` outgoing operator (not asserted from the closed form), whereas the `.py` posits the closed form directly.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `Γ̄₅/Γ̄₅^target = χ_Q N_Q` | py:62-66 / wl:74-75 | match |
| (2) `m̂₀²χ_Q N_Q = 1` (from observable `m̂₀²Γ̄₅=Γ̄₅^target`) | py:87-93 / wl:80-101 | match |
| (3) `Δ_norm = P0^target(1/χ_Q−1)` and `Δ_Q` form | py:94-101 / wl:102-109 | match |
| (4) `N_Q=1/χ_Q`, `N_Q−1=−(χ_Q−1)+O(²)` | py:108-120 / wl:114-131 | match |
| (5) deformation `χ_Q`, `N_Q`, `Δ_norm^(pt)`, linearization | py:130-171 / wl:136-206 | match |
| (5c) canonical closure `χ_Q=N_Q=1, Δ_norm=0` | py:178-187 / wl:213-219 | match |

`paper_alignment: aligned`. Every boxed paper deliverable has a faithful, non-tautological script-side check in both engines, with consistent constants (`P0^target = 54Gc_s⁵/5a⁵c⁵`, `Γ̄₅^target = 2G/5c⁵`, `χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)`).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `expect_zero(Gamma5 - Gamma5_alt)` | (1) two Γ̄₅ forms agree | yes |
| A2 | sympy | 63-66 | `expect_zero(Gamma5/Gamma5_target.subs(P0,N_Q·P0_target) - chi_Q·N_Q)` | (1) odd ratio | yes |
| A3 | sympy | 87-93 | `expect_zero(oddResidual - Gamma5_target·(mhat0²χ_Q N_Q−1))` | (2) odd closure from observable cond. | yes |
| A4 | sympy | 94 | `expect_zero(Delta_norm_from_odd - P0_target(1/chi_Q−1))` | (3) collapse | yes |
| A5 | sympy | 98-101 | `expect_zero(Delta_norm_DeltaQ + P0_target·Delta_Q/(1+Delta_Q))` | (3) Δ_Q form | yes |
| A6 | sympy | 115-120 | `expect_zero(NQ_pt-1/chi_Q ; Delta_norm_pt ; NQ_pt-1 form)` | (4) point-particle | yes |
| A7 | sympy | 143-147 | `expect_zero(NQ_from_def ; Delta_norm_from_def)` | (5) deformation laws | yes |
| A8 | sympy | 164-171 | `expect_zero(linearized N_Q-1 ; linearized Delta_norm)` | (5) linearization | yes |
| A9 | sympy | 178-187 | `expect_zero(chi/N_Q/Delta canonical)` | (5c) canonical closure | yes |
| B1 | mathematica | 74-75 | `expectZero(gammaByGeometry-gammaByPole ; gammaRatio-chiQ·nQ)` | (1) | yes |
| B2 | mathematica | 98-101 | `expectZero(oddObservableResidual - gamma5Target·(mhat0²χ_Q N_Q−1))` | (2) | yes |
| B3 | mathematica | 102-109 | `expectZero(deltaAfterOdd ; Delta_Q form)` | (3) | yes |
| B4 | mathematica | 123-131 | `expectZero(nQNatural ; deltaNatural ; Delta_Q form)` | (4) | yes |
| B5 | mathematica | 146-154 | `Solve` σ2,σ4 from even moments `1/9, 4/81`; uniqueness gate | (5) even-constraint derivation | yes |
| B6 | mathematica | 156-159 | `chiFromDtn` read off x⁵ coeff of normalized Hankel DtN | (5) χ_Q DERIVED, not posited | yes |
| B7 | mathematica | 174-175 | `expectZero(nQFromDtn-target ; deltaDeformation-target)` | (5) | yes |
| B8 | mathematica | 199-206 | `expectZero(linearized laws)` | (5) | yes |
| B9 | mathematica | 213-219 | `expectZero(canonical closure)` | (5c) | yes |

All rows anchored. No tautology survives: the F1 headline factorization is, in both engines, derived from the independently-built observable residual `(m̂₀²Γ̄₅ − Γ̄₅^target)` (py:87-89 / wl:80-83) and checked to *equal* the factorized form — it is not the trivial `(m̂₀²χ_Q N_Q − 1) − (m̂₀²χ_Q N_Q − 1)` self-echo flagged in pass-1. F2's two deleted definitional-echo checks have not re-appeared.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.txt` (mtime 2026-06-01 11:53:34)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py` (mtime 2026-06-03 15:59:11)

**What's wrong:**
The committed SymPy `.txt` predates the current `.py` by ~2 days. The staleness is also visible in content: the saved output's banner reads `STAGE 178 — EXACT SOURCE-MAP REDUCTION…` and `STAGE 178 LEDGER` (lines 3, 122), but the current `.py:35` prints `STAGE 195 — …` and `.py:189` prints `STAGE 195 LEDGER`. So the captured transcript was produced by an earlier revision (pre-renumber and pre-F1/F2 de-taut). The Mathematica `.txt` is fresh relative to its `.wl` (both 2026-06-01; output `11:53` ≥ script `11:49`) and already carries the correct `STAGE 195` banner.

**Why this matters:**
The committed SymPy transcript does not reflect the current script — anyone trusting the checked-in `.txt` sees a stale STAGE-178-labelled run, not the current STAGE-195 de-tautologized assertions. Informational only; the math itself is sound (the current `.py` assertions all reduce to 0 by inspection — see Self-test notes).

**Required change:**
Re-run `python3 scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py` and overwrite the saved `.txt`. No source edit. The refreshed transcript must show the `STAGE 195` banner and all `expect_zero(...) = 0` lines.

**Verification:**
After refresh, `scripts/output/...sympy_audit.txt` mtime ≥ `.py` mtime, line 3 reads `STAGE 195 …` (not 178), line ~189 reads `STAGE 195 LEDGER`, and every `expect_zero` line ends `= 0`.

## Independent-derivation check (Mathematica)

GENUINELY INDEPENDENT on the load-bearing Section IV; necessarily structurally-similar (but not a port) on the short algebraic Sections I–III.

The heart of the V.3 retrofit concern is whether the `.wl` re-derives rather than echoes. Evidence:

- **Section IV — DIFFERENT METHOD (the discriminating section).** `.py:130`: `chi_from_def = sp.simplify(3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0))` — the closed form is *posited* as a hardcoded carry-forward from Stage 194. `.wl:136-159` instead *derives* it: `lambdaOut = FunctionExpand[x*D[SphericalHankelH1[2,x],x]/SphericalHankelH1[2,x]]` (the genuine outgoing l=2 DtN operator), Series-expanded to O(x⁵) (`lambdaWindow`), deformed (`scaleS·λ(βx)+σ0+σ2x²+σ4x⁴+iσ5x⁵`), normalized, then σ2/σ4 are `Solve`d (wl:146-153) by matching the canonical even moments `1/9` and `4/81`, with a uniqueness gate (wl:154); χ_Q is finally read off the x⁵ coefficient `chiFromDtn = Coefficient[normalizedDtn/.evenRules, x, 5]/(I/27)` (wl:156-159). The `.py` never touches a Hankel function or a Series of the operator (confirmed: no `SphericalHankel`/operator-Series in `.py`). The two engines reach the same `3(Sβ⁵+9Σ₅)/(3S−Σ₀)` by genuinely different routes — the `.wl` proves what the `.py` asserts.

- **Sections I–III — short closed-form algebra.** Both engines carry the same definitional building blocks (`Γ̄₅ = χ_Q·a⁵·P̄₀/27c_s⁵`, `Γ̄₅^target = 2G/5c⁵`) and form the same observable residual. Because these are one-line algebraic identities (not a multi-step choreography), there is no "independent" alternative route to a definition; the similarity is structural necessity, not transliteration. The `.wl` does use its own constructs (`Solve[oddClosureEquation, nQ, Reals]` at wl:85-88 to obtain `N_Q from odd closure`, vs the `.py` `1/(mhat0**2*chi_Q)` direct write at py:76) — a mild route difference even here. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every deliverable. Side-by-side of the load-bearing outputs:

- Section IV `χ_Q`: py `.txt:80-84` → `3(S·β⁵+9·Σ₅)/(3S−Σ₀)`; wl `.txt:46` → `(3*(betaStretch^5*scaleS + 9*sigma5))/(3*scaleS - sigma0)`. Identical.
- Linearized `N_Q−1`: py `.txt:100-103` → `−5·epsᵦ − dSigma₀/3S − 9·dSigma₅/S`; wl `.txt:53` → `−5*epsBeta − dSigma0/(3*scaleS) − (9*dSigma5)/scaleS`. Identical.
- `Δ_norm` after odd closure: py `.txt:50-55` → `54Gc_s⁵(1−χ_Q)/(5a⁵c⁵χ_Q)`; wl `.txt:23` → `(-54*bigG*(-1 + chiQ)*soundSpeed^5)/(5*chiQ*lightSpeed^5*radius^5)`. Identical.

All `expect_zero`/`expectZero` checks pass (`= 0` / `PASS:`). No `engine_disagreement`. (Note: the SymPy `.txt` is stale per F1, but its emitted *values* still match the current `.py`'s math; the only diff is the STAGE-178 banner label, not the results.)

## Verdict justification

`verdict: findings` with a single low-severity `stale_output` (F1). The math holds up against the paper under attack. The pass-1 F1 concern (X−X self-echo of `m̂₀²χ_Q N_Q=1`) is resolved: the headline identity is, in BOTH engines, derived from the independently-built observable odd residual `(m̂₀²Γ̄₅ − Γ̄₅^target)` and checked to equal `Γ̄₅^target·(m̂₀²χ_Q N_Q − 1)` — a real, falsifiable factorization, not a definitional tautology. The pass-1 F2 deletions did not re-surface. The V.3 retrofit concern is satisfied: the `.wl` Section IV is a genuinely independent re-derivation of `χ_Q` from the `SphericalHankelH1[2,x]` operator, not a transliteration of the `.py`'s posited closed form; Sections I–III are short algebraic identities where structural similarity is unavoidable and not a port. Constants and targets match the paper card, notes, and appendix exactly. Attacks tried: (a) hunt for surviving X−X tautology in the odd-closure factorization — failed, the residual is built from observable inputs; (b) check the `.wl` deformation route is not a hardcoded echo — failed, it is a true Hankel-operator derivation; (c) parity/derivative traps in the Series/linearization — none (see Self-test); (d) constant mismatch vs paper — none. The only defect is the stale SymPy transcript, which is informational and fixed by a re-run. A directive is written for the re-run.

## Self-test notes

Checked: (1) Variable independence — the only `D[...]` is wl:136 `D[SphericalHankelH1[2,x],x]` w.r.t. `x`, on which the Hankel function genuinely depends (nonzero), and the `.series(t,0,2)`/`Series[...,{tau,0,1}]` linearizations are in the deformation parameter `t`/`tau` on which `NQ_from_def` genuinely depends — no identically-zero-derivative trap. (2) Symmetry/parity — no unbounded integrals; the even/odd structure is a coefficient match (`x²,x⁴` even constraints, `x⁵` odd readout), self-consistent with the operator. (3) Trivial-case — the canonical substitution `β=1,Σ₀=Σ₅=0` gives `χ_Q = 3(S+0)/(3S) = 1`, `N_Q=1`, `Δ_norm=0` by hand, matching V's checks. Conclusion: no silent-pass trap; F1 is a transcript-freshness issue only.

## Value Reconciliation (pass-2 augmentation)

Every emitted deliverable value reconciles against the card and/or notes. The deliverable-level table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `P0^target = 54Gc_s⁵/(5a⁵c⁵)` | py:46 / wl:59 / py.txt:13-18 | notes:69-73 (boxed); card omits (terse) | MATCH |
| `Γ̄₅^target = 2G/(5c⁵)` | py:47 / wl:60 / py.txt:19-23 | notes:119-123 (boxed) | MATCH |
| `Γ̄₅ = χ_Q·a⁵·P̄₀/(27c_s⁵)` | py:48 / wl:62 / py.txt:24-29 | notes:97 (boxed) | MATCH |
| `Γ̄₅/Γ̄₅^target = χ_Q N_Q` | py:62-66 / wl:75 | notes:113-117 (boxed); card "Derivation ledger" | MATCH |
| `m̂₀²χ_Q N_Q = 1` (headline) | py:75,87-93 / wl:84,98-101 | card:15 `\stagefield{Output}`; notes:138-140 (boxed); appendix:1309 | MATCH |
| `N_Q = 1/(m̂₀²χ_Q)` | py:76 / wl:85-88 / py.txt:45-49 | notes:144-148 (boxed) | MATCH |
| `Δ_norm = P0^target(1/χ_Q − 1)` | py:94 / wl:104 / py.txt:50-55 | notes:174-178 (boxed); appendix:1306 | MATCH |
| `Δ_norm = −P0^target·Δ_Q/(1+Δ_Q)`, `Δ_Q:=χ_Q−1` | py:96-101 / wl:106-109 | notes:185-193 (boxed) | MATCH |
| `N_Q = 1/χ_Q` (point-particle) | py:108,115 / wl:114,123 / py.txt:63-66 | notes:214-225 (boxed); appendix:1316 | MATCH |
| `N_Q−1 = −(χ_Q−1)+O(²)` | py:117-120 / wl:128-131 | notes:236-241 (boxed) | MATCH |
| `χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` | py:130 / wl:156-159 / py.txt:80-84 | notes:263-266 (boxed); appendix:1297-1302 | MATCH |
| `N_Q = (3S−Σ₀)/(3(Sβ⁵+9Σ₅))` | py:131,145 / wl:160-161 / py.txt:85-91 | notes:277-281 (boxed) | MATCH |
| `Δ_norm^(pt) = −P0^target[3S(β⁵−1)+Σ₀+27Σ₅]/[3(Sβ⁵+9Σ₅)]` | py:133-135 / wl:163-167 / py.txt:92-97 | notes:308-314 (boxed) | MATCH |
| `N_Q−1 = −5ε_β − δΣ₀/3S − 9δΣ₅/S + O(2)` | py:166 / wl:201 / py.txt:100-103 | notes:333-340 (boxed) | MATCH |
| canonical closure `χ_Q=1, N_Q=1, Δ_norm=0` | py:178-183 / wl:213-215 | notes:382-386 (boxed) | MATCH |

INTERNAL scaffolding (no finding): `Omega_Q = 3c_s/2a`, `Gamma5_alt`/`gammaByPole`, even-moment Solve constraints `σ2,σ4` and targets `1/9`,`4/81` (wl:146-153, intermediate to the χ_Q derivation), pass/fail flags, the `Δ_norm in terms of Δ_Q` residual checks. `Σ₂,Σ₄` closed forms appear in notes:255-260 but are intermediate Stage-194 carry-forwards, not Stage-195 deliverables.

reconciliation: complete; 15 deliverable values checked, 0 misaligned.
