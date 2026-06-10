---
unit_id: 231
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 231 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_231.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 74, 848-865, 881)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.txt`

## What the paper claims

Stage 231 pulls the Stage-230 selected-branch dynamic classifier `R_ND(xi,delta) = 72 delta^2 (1-xi) / ((9 delta + 11 xi)(9 delta^2 + 18 delta xi + 11 xi^2))` back through the exact continuum-placement map so the same-charge dynamic verdict is stated in the physical kernel ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)`. The card `\stagefield{Output}` reads verbatim: "Continuum-placement version of the static-first verdict: the branch decision can be made directly from the physical placement packet, with $R_{\rm target}$ controlling the selected-class motion monotonically." The notes enumerate five deliverables: (1) the physical classifier `R_phys(delta,R_target) := R_ND(xi_req(delta,R_target),delta)` where `xi_req` solves `F(xi,delta)=R_target`; (2) the monotonicity theorem `∂_{R_target} R_phys = (∂_xi R_ND)/(∂_xi F) < 0` (equivalently `∂_{M_mix} R_phys > 0`); (3) a generic threshold compiler `P_c(xi,delta)` with `∂_xi P_c > 0`, onset `P_c(0,delta)=9 delta^2(9 c delta - 8)`, onset law `delta_c = 8/(9c)`, and pulled-back curves `R_flip(delta)=R_{R_*}(delta)`, `R_den(delta)=R_1(delta)`; (4) the continuum placement formulas (`delta=delta_0/(1-eps_eta)`, `M_mix`, `R_target`) and the product law `R_target·M_mix = 8 Lambda (1-eps_W)/pi^2`, with equivalent kernel inequalities; (5) the static-first theorem `B_dyn > B_stat` surviving the pullback. Section-5 sample slices give `R_flip/R_den` at delta=0.25/0.50/0.75. The card's `\stagefield{Verification}` line states "SymPy audit: ... Mathematica audit: none yet." (this is now stale — see F1).

## What the script claims to verify

The SymPy script verifies, against retyped published closed forms it independently re-derives with `sp.diff`/`sp.factor`/`sp.limit`/`sp.solve`: the three monotonicity numerators (`dF/dxi`, `dG/dxi`, `dR_ND/dxi`) and their signs on a 6-point stable-branch grid; the endpoint data (`F(0)=1`, `G(0)=0`, `R_ND(0)=8/(9 delta)`, `F→∞` as `xi→1⁻`); the pullback slope `dR_phys/dR_target = dR_dxi/dF_dxi < 0` numerically on the grid; the threshold polynomial `P_c`, its derivative, onset factorization, and `delta_c=8/(9c)`; the carried `R_*=s_-^den/(-s_-^num)≈1.229255438463336` and `delta_*^dyn≈0.723111617875019`; numerically-bisected pulled-back thresholds matching the section-5 sample rows, the collapse to onset for delta≥onset, and `R_flip ≤ R_den`; the placement map and product-law identity; and the carried static-first strict inequalities plus the `R_phys ⊆ [0, R_ND(0)]` subset bound on the grid. The Mathematica script verifies the same set but proves the three monotonicity sign claims and `dP_c>0` via `Resolve[ForAll[...]]` quantifier elimination, finds threshold roots via `NSolve` over Reals with a uniqueness `Select`, and checks the soft-limit pole via `Limit[...,Direction->"FromBelow"]` and a reciprocal-is-zero idiom.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) R_phys := R_ND(xi_req,delta), xi_req solves F=R_target | sympy `xi_cap`/`R_cap` bisect F=cap & eval R_ND; .wl `thresholdXi`/`thresholdR` | match |
| (2) ∂_{R_target}R_phys = dR_ND/dF < 0 | sympy L122,150-153 `dRphys_num<0` on grid; .wl M3 sampled slope<0 | match |
| (3) P_c, ∂_xi P_c>0, onset 9δ²(9cδ-8), δ_c=8/(9c) | sympy L160-169; .wl M4 (with `Resolve[ForAll]` for >0) | match |
| (3') R_flip(δ)/R_den(δ) sample slices (sec.5) | sympy L207-223 vs retyped rows; .wl M5 vs rows | match |
| (3'') collapse to onset for δ≥δ_c | sympy L234-235; .wl M5 collapse | match |
| (4) placement map + product law 8Λ(1-eps_W)/π² | sympy L247-252; .wl M6 | match |
| (4') equiv. kernel/M_mix inequalities | sympy L261-266 (printed translators) | match (algebraic rearrangement, printed) |
| (5) static-first B_dyn>B_stat survives | sympy L272-286; .wl M7 + subset bound | match |
| R_* ≈ 1.229255438463336, δ_*^dyn ≈ 0.723111617875019 | sympy L181-188; .wl rStar L130 | match |
| card Verification: "Mathematica audit: none yet" | a passing `.wl` now EXISTS | mismatch (F1) |

`paper_alignment: partial` — every math deliverable aligns across card/notes/appendix and both engines; the lone defect is the card's stale "none yet" Mathematica-verification claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 94 | `simplify(dF_dxi - expected_dF_dxi)==0` | claim 1/2 (∂_xi F) | yes |
| A2 | sympy | 95-96 | `simplify(dG_dxi - exp)==0`, `dR_dxi - exp==0` | claim 1/2 | yes |
| A3 | sympy | 102-105 | onset `F=1,G=0,R_ND=8/(9δ)`, `lim F=∞` | claim 1, geometry | yes |
| A4 | sympy | 146-153 | grid `dF>0,dG>0,dR<0,dRphys<0` | claim 2 | yes (concrete points) |
| A5 | sympy | 163,166,169 | `dP_dxi`, onset factor, `delta_c` onset=0 | claim 3 | yes |
| A6 | sympy | 187-189 | `R_*`, `delta_*^dyn`, `delta_den` close | claim 3/4 thresholds | yes |
| A7 | sympy | 220-225 | sample R_flip/R_den vs rows, R_flip≤R_den | claim 3' | yes |
| A8 | sympy | 234-235 | collapse to onset =1 | claim 3'' | yes |
| A9 | sympy | 252 | product law residual ==0 | claim 4 | yes |
| A10 | sympy | 277-286 | B_dyn>B_stat, R_phys⊆[0,R_ND(0)] | claim 5 | yes |
| M1 | math | 88-92 | `expectZero` dG/dF-poly; `Resolve[ForAll]` dF>0,dG>0,dRnd<0 | claim 1/2 | yes (QE) |
| M2 | math | 99-102 | onset + pole reciprocal=0 | claim 1 | yes |
| M3 | math | 112 | sampled pullback slope<0 | claim 2 | yes |
| M4 | math | 125-128 | dPc, Resolve dPc>0, onset, onset law | claim 3 | yes |
| M5 | math | 154-162 | sample R_flip/R_den, ≤, collapse | claim 3'/3'' | yes |
| M6 | math | 173 | product law residual=0 | claim 4 | yes |
| M7 | math | 181-190 | B_dyn>B_stat, RND subset | claim 5 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_231.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl` (exists, passes)

**What's wrong:**
The card's Verification field claims no Mathematica engine exists:

> `\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py}.  Mathematica audit: none yet.}` (stage_231.tex:11)

But a Mathematica audit `.wl` now exists for Stage 231, is committed (Jun 2), and passes all checks — its saved output ends "All Stage 231 Mathematica audits passed." The card therefore understates verification coverage: it states the second engine is absent when it is present and green. This is the identical situation pass-2 flagged on the sibling Stage 230 (report stage_230.md F1, subtype `paper_missing_script_claim`).

**Subtype:** paper_missing_script_claim

**Why this matters:**
A reader trusting the card would believe Stage 231 is single-engine and uncorroborated, when in fact both engines independently confirm every deliverable. Direction of fix (update the card vs. justify excluding the `.wl`) is a user call; Codex must not silently edit the card.

**Required change:**
See `## Resolve before fix_loop` in the directive. No script edit — the scripts and their math are correct. The expected direction under the dual-engine-required policy is a paper-side prose update to the Verification line citing the `.wl`.

**Verification:**
After user resolution, if direction (a): `stage_231.tex:11` Verification line should cite the `.wl` and `git grep "none yet"` scoped to stage 231 should be empty.

## Independent-derivation check (Mathematica)

The `.wl` is genuinely independent, in fact stronger than the `.py`, not a transliteration. Three corresponding sections:
- Monotonicity signs: the `.py` only probes a 6-point grid (`dF_num(xi,delta) > 0`, py L146-151); the `.wl` proves them universally — `expectTrue["M1 D[F,x] > 0 on stable strip", Resolve[ForAll[{x, d}, Implies[0 <= x < 1 && d > 0, dF > 0]], Reals]]` (wl L90-92). This is quantifier elimination, an operation absent from the `.py`.
- Threshold roots: the `.py` uses a hand-written bisection `bisect_increasing` (py L13-39); the `.wl` uses `NSolve[poly==0, x, Reals, WorkingPrecision->80]` then `Select[..., 0<#<1]` with a uniqueness check (wl L52-61). Different root-finding route.
- Soft limit: the `.py` uses `sp.limit(F, xi, 1, dir="-")` against `sp.oo` (py L101,105); the `.wl` uses `Limit[f, x->1, Direction->"FromBelow"]` plus the calibrated reciprocal-is-zero pole idiom `expectInfinityPole` (wl L46-50,102). Independent pole handling.
The closed-form residual checks (M1 dG/dF-poly, M4 dPc) compare Mathematica's own `Factor[D[...]]` against retyped published polynomials — non-tautological (the published poly could be wrong; both engines re-derive). Not a port.

## Engine cross-check

Both engines agree exactly. Final sample thresholds: SymPy `R_flip=1.330868539, R_den=1.393832566` (delta=0.25), `1.139956630/1.221087062` (0.50), `1.000000000/1.071471867` (0.75) (sympy output L30-32); Mathematica `R_flip=1.33086853932996..., R_den=1.39383256578151...` etc. (mathematica output L52-72) — agree to the 9-digit row precision. `R_*`: sympy `1.22925543846333...` (L25) vs `.wl` `rStar=411024574532864/334368725711457` = 1.229255438463336 (confirmed by read-only arithmetic). Product-law residual = 0 in both (sympy L39 / mathematica L80). Monotonicity-numerator closed forms identical modulo SymPy's `(xi-1)` vs Mathematica's `(-1+x)` canonicalization. No `engine_disagreement`.

## Verdict justification

Every math deliverable holds up against the paper, on both engines. Attacks tried and failed: (a) variable-independence trap — every `sp.diff`/`D[]` differentiates `F`/`G`/`R_ND`/`P_c` w.r.t. `xi`, on which they genuinely depend; the slope checks evaluate at concrete grid points giving nonzero literals, never a vacuous `0`. (b) Tautology — the `dF_dxi - expected_dF_dxi`-style asserts compare an engine-computed derivative against a retyped published form; they would fail if the published coefficients were wrong (this is exactly the pass-1 240/189→189/121 fix, now reconciled: notes lines 84/98/162 read 189/121, both scripts emit 189/121, output L6/L3 confirm). (c) Independence — the `.wl` adds `Resolve[ForAll]` global proofs the `.py` lacks. (d) Symbol assumptions — `xi,delta,c>0` and the `.wl` `0<=x<1 && d>0 && c>0` match the physical strip; no over-strong `simplify` masking. The only defect is the card's stale "Mathematica audit: none yet" line (F1), a `paper_misalignment` routed to the user gate. Verdict: `findings` (one paper_misalignment, low severity).

## Self-test notes

Checked variable-independence (all diffs over `xi`, a real dependency — no identically-zero-derivative trap; slope asserts give nonzero grid literals), trivial-case (at delta=0.75 `xi_flip=0` → `R_flip=1` exactly, output L66, the onset-collapse limit; at xi=0 `R_ND=8/(9δ)` reproduces onset), and symmetry/parity (no unbounded-domain integrals in this stage). The single finding (F1) is paper-side prose only and prescribes no script edit, so no Codex change can introduce a new paper_misalignment. Engine agreement and the pass-1 189/121 coefficient reconciliation both reconfirmed.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| R_* = 1.229255438463336 | py L183-187, sympy out L25; wl L130 | notes L152,354; appendix (carried) | MATCH |
| delta_*^dyn = 0.723111617875019 | py L184-188, sympy out L26 | notes L160,358 | MATCH |
| delta_den = 8/9 = 0.888888... | py L186-189, sympy out L27 | notes L162,398 | MATCH |
| delta_c = 8/(9c) | py L168, sympy out L22; wl L121 | notes L323,616 | MATCH |
| dF/dxi numerator poly 81δ³+189δ²ξ+72δ²+297δξ²+121ξ³ | py L83, sympy out L6; wl L81, math out L3 | notes L98 | MATCH (pass-1 fix held) |
| dR_ND/dxi numerator 81δ³+261δ²+297δξ(2-ξ)+ξ²(363-242ξ) | py L87-92; wl L78, math out L5 | notes L144 | MATCH |
| dP_c/dxi = 3(87cδ²+198cδξ+121cξ²+24δ²) | py L162, sympy out L20; wl L119, math out L43 | notes L302 | MATCH |
| P_c(0,delta) = 9δ²(9cδ-8) | py L166, sympy out L21; wl L120 | notes L313 | MATCH |
| product law R_target·M_mix = 8Λ(1-eps_W)/π² | py L250-252, sympy out L39; wl L168, math out L80 | notes L206-208,623 | MATCH |
| R_flip/R_den samples (0.25,0.50,0.75) | py L207-210, sympy out L30-32; wl L134-138, math out L52-72 | notes L447-449 | MATCH |
| B_dyn^both inf 0.967282389363822, B_stat^both 0.367930328492646 | py L272-275, sympy out L45 | notes L531-536 | MATCH |
| B_dyn^nonempty inf 0.990581810705233, B_stat^nonempty 0.737619063660757 | py L272-275, sympy out L46 | notes L539-544 | MATCH |
| placement map delta=δ0/(1-eps_eta), M_mix, R_target | py L247-249, sympy out L36-38; wl L165-167, math out L77-79 | notes L192-203 | MATCH |

INTERNAL (no finding): pass/fail flags, residuals (=0), `assert_close` tolerances, `s_minus_den`/`s_minus_num` carried slope rationals (intermediate inputs to R_*), probe-grid coordinates, bisection iteration bounds, dxi_req/dR_target & dR_phys/dR_target printed intermediate forms.

reconciliation: complete; 13 deliverable values checked, 0 misaligned. (The only paper_misalignment, F1, is a prose verification-coverage claim, not a value mismatch.)
</content>
</invoke>
