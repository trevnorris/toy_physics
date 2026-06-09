---
unit_id: 205
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T21:25:48Z
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
  notes_stage_files: [moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md]
  paper_appendix: present
---

# Audit unit 205 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_205.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows/blocks at lines 41, 676-729)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.txt`

## What the paper claims

Stage 205 upgrades the Stage 204 scalarized log-ray closure problem `\(\Phi_{\mathbf s}(\tau)=1\)` to exact second order. `\stagefield{Output}`: "Affine and logarithmic quadratic predictors, discriminant tests, curvature-corrected local expansions, and turning-point/tangency theorem." The notes enumerate seven deliverables: (1) directional Hessian operators `\(\mathcal H_{\mathbf s}=\mathcal D_{\mathbf s}^2\)`; (2) the first/second derivative scalars `\(\Phi_0,\Phi_1,\Phi_2,L_0,L_1\)` with the central bridge identity `\(\Phi_2=\Phi_0(L_1+L_0^2)\)` and `\(L_1=\Phi_2/\Phi_0-\Phi_1^2/\Phi_0^2\)`; (3) the exact quadratic affine predictor `\(\tau_{\rm quad}=2(1-\Phi_0)/(\Phi_1+\operatorname{sgn}(\Phi_1)\sqrt{\Delta_{\rm aff}})\)` and quadratic log predictor `\(\tau_{\log,2}=-2\ln\Phi_0/(L_0+\operatorname{sgn}(L_0)\sqrt{\Delta_{\log}})\)`; (4) the discriminants `\(\Delta_{\rm aff}=\Phi_1^2-2\Phi_2(\Phi_0-1)\)`, `\(\Delta_{\log}=L_0^2-2L_1\ln\Phi_0\)`; (5) the turning-point/tangency theorem: for `\(\Phi_1=0,\Phi_0\neq1\)`, real roots iff `\((1-\Phi_0)\Phi_2>0\)`, with `\(\tau_\pm=\pm\sqrt{2(1-\Phi_0)/\Phi_2}\)`, and the `\(\Phi_0=1,\Phi_1=0\)` tangency model `\(\Delta=\tfrac12\Phi_2\tau^2\)`; (6) curvature-corrected expansions `\(\tau_{\rm quad}=\tau_{\rm aff}-\tfrac{\Phi_2}{2\Phi_1}\tau_{\rm aff}^2+O(\tau_{\rm aff}^3)\)` and the log analogue; (7) the agreement theorem `\(\tau_{\log,2}-\tau_{\rm quad}=\tfrac{(L_0^2+3L_1)}{6L_0^3}\varepsilon^3+O(\varepsilon^4)\)`. The card classifies this as `\StatusExactClosure`. The card's `\stagefield{Verification}` line states "SymPy audit: ... Mathematica audit: none yet."

## What the script claims to verify

The SymPy script verifies all seven deliverables symbolically: Section I posits an abstract gradient/Hessian and confirms the bridge identity and `\(L_1\)` form; Section II posits the affine quadratic-formula predictor (both slope signs) and checks the closure residual is zero plus the `\(\Phi_2\to0\)` limit; Section III does the same for the log predictor (`\(L_1\to0\)` limit); Section IV verifies the turning-point sign bridge, the positive/negative product reality cases, the `\(\tau_\pm\)` residual, and the tangency model; Section V series-expands the curvature corrections and confirms the agreement theorem cubic coefficient `\((L_0^2+3L_1)/(6L_0^3)\)`. The Mathematica script verifies the same eleven facts (M1-M11) but by independent methods (symbolic `D[Log[...]]`, `Solve`+branch-select, `Reduce`/`Exists` region equivalence). All in-file checks emit zero residuals / True / PASS in both saved outputs.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Hessian operator `\(\mathcal H_{\mathbf s}=\mathcal D_{\mathbf s}^2\)`, `\(\Phi_2=\mathbf s^T H\mathbf s\)` | py L64 `Phi2`; wl L143-154 `chiAlong`/`logHessianByD` | match |
| Bridge `\(\Phi_2=\Phi_0(L_1+L_0^2)\)`, `\(L_1=\Phi_2/\Phi_0-\Phi_1^2/\Phi_0^2\)` | py L77-78; wl M1/M2 L156-163 | match |
| Affine quadratic predictor (both signs) + discriminant | py L86-112; wl M3/M4 L165-202 | match |
| Log quadratic predictor (both signs) + discriminant | py L119-145; wl M5/M6 L204-241 | match |
| Turning-point criterion `\((1-\Phi_0)\Phi_2>0\)`, `\(\tau_\pm\)` | py L152-179; wl M7 L243-285 | match |
| Tangency model `\(\Delta=\tfrac12\Phi_2\tau^2\)` | py L181-189; wl M8 L287-291 | match |
| Curvature corrections (ordinary + log) | py L200-217; wl M9/M10 L293-332 | match |
| Agreement theorem cubic coeff `\((L_0^2+3L_1)/(6L_0^3)\)` | py L219-226; wl M11 L334-341 | match |
| Card `\stagefield{Verification}`: "Mathematica audit: none yet" | `.wl` audit EXISTS (added pass-1 retrofit) | mismatch (F1) |

Every mathematical deliverable maps to a faithful, non-tautological check in both engines. The only discrepancy is documentary: the card's Verification line denies the existence of the Mathematica audit that is now present. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 77 | `expect_zero(L1 - (Phi2/Phi0 - Phi1^2/Phi0^2))` | deliverable 2 | yes |
| A2 | sympy | 78 | `expect_zero(Phi2 - Phi0*(L1+L0^2))` | deliverable 2 (bridge) | yes |
| A3 | sympy | 97/108 | affine residual = 0 (both signs) | deliverable 3 | yes |
| A4 | sympy | 98/111 | `limit Phi2->0 = tau_aff` (both signs) | deliverable 6/zero-curv | yes |
| A5 | sympy | 130/141 | log residual = 0 (both signs) | deliverable 3 | yes |
| A6 | sympy | 131/144 | `limit L1->0 = tau_log` | deliverable 6 | yes |
| A7 | sympy | 154-179 | turning sign-bridge + reality cases + `\(\tau_+\)` residual | deliverable 5 | yes |
| A8 | sympy | 186-189 | tangency model identity | deliverable 5 | yes |
| A9 | sympy | 210/217 | curvature correction series = 0 | deliverable 6 | yes |
| A10 | sympy | 222-226 | agreement cubic coeff `\((L_0^2+3L_1)/(6L_0^3)\)` | deliverable 7 | yes |
| M1 | math | 156 | `expectZero(logCurve - (ordinaryCurve/phiBase - ordinarySlope^2/phiBase^2))` | deliverable 2 | yes |
| M2 | math | 160 | `expectZero(ordinaryCurve - phiBase*(logCurve+logSlope^2))` | deliverable 2 | yes |
| M3 | math | 172/191 | affine residual on Solve-selected root (both signs) | deliverable 3 | yes |
| M4 | math | 176/195 | affine zero-curvature limit (both signs) | deliverable 6 | yes |
| M5 | math | 211/230 | log residual on Solve-selected root (both signs) | deliverable 3 | yes |
| M6 | math | 215/234 | log zero-curvature limit (both signs) | deliverable 6 | yes |
| M7 | math | 271-284 | `Reduce`/`Exists` region equivalence + `\(\tau_\pm\)` residuals | deliverable 5 | yes |
| M8 | math | 291 | tangency model identity | deliverable 5 | yes |
| M9/M10 | math | 324/332 | curvature correction series = 0 | deliverable 6 | yes |
| M11 | math | 335-340 | gap series eps^0..eps^3 with cubic coeff | deliverable 7 | yes |

No assertion is tautological. Every check substitutes a candidate predictor (posited in py, Solve-derived in wl) back into a closure model that was defined separately, so a wrong predictor would produce a nonzero residual. The turning-point reality checks use independent positivity reasoning (`sp.ask` / `Reduce`).

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_205.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl` (exists, 9367 bytes, passes)

**What's wrong:**
The stage card Verification line states the Mathematica audit does not exist:

> `paper/stages/stage_205.tex:11`: "SymPy audit: \StageFile{scripts/...stage205...sympy_audit.py}.  Mathematica audit: none yet."

But a Mathematica audit script now exists and passes cleanly (M1-M11 all PASS in the saved output). This is a pass-1 retrofit artifact: the `.wl` was added but the card's Verification pointer was not updated to acknowledge it. The card therefore under-reports the verification coverage that actually exists.

**Why this matters:**
The card is the authoritative statement of how the stage is verified. Leaving "none yet" makes the dual-engine coverage invisible to a reader and is internally inconsistent with the committed artifact. It is a documentary inconsistency, not a math defect — the scripts agree and are correct.

**Resolution:** routed to user (see directive `## Resolve before fix_loop`). Codex must not edit the card autonomously. Direction (update the card's Verification line to cite the `.wl`, vs. leave a deliberate "none yet" and instead delete/quarantine the `.wl`) is the user's call. Per project memory the natural direction is to update the card + the MATHEMATICA_MIRROR_POLICY tracker during the post-batch sync; this finding records the discrepancy so it is not lost.

## Independent-derivation check (Mathematica)

Verdict: **INDEPENDENT.** The `.wl` is not a transliteration of the `.py`. For every load-bearing object the extraction method differs:

1. **Log-Hessian / bridge (M1-M2 vs Section I).** SymPy POSITS the log-Hessian closed form and checks an algebraic identity:
   - py L66-67: `Hlog = sp.simplify(Hchi / Phi0 - (gvec * gvec.T) / Phi0**2)`; `L1 = (svec.T * Hlog * svec)[0]`.
   Mathematica instead DERIVES the log curvature by symbolic second differentiation of `Log` of an explicit local Taylor model:
   - wl L143-152: `chiLocal = phiBase + gradVec . vars + 1/2*vars . hessMat . vars`; `logHessianByD = Table[D[Log[chiLocal], vars[[i]], vars[[j]]] /. Thread[vars -> ConstantArray[0, 5]], ...]`.
   Derive-by-`D[Log[...]]` vs posit-the-formula → independent.

2. **Quadratic predictors (M3-M6 vs Sections II/III).** SymPy POSITS the quadratic-formula closed form:
   - py L88: `tau_quad = sp.simplify(2 * (1 - Phi0a) / (Phi1a + sp.sqrt(Delta_aff)))`.
   Mathematica calls `Solve` to obtain BOTH roots and selects the physical branch by a zero-curvature limit:
   - wl L110-116: `roots = ... tau /. Solve[ordinaryModel[p0, p1, p2, tau] == 0, tau]`; `selectByLimit["ordinary quadratic", roots, p2, target, assumptions]`.
   Solve+branch-select vs posited radical → independent.

3. **Turning-point reality (M7 vs Section IV).** SymPy substitutes positive/negative case symbols and queries `sp.ask(Q.positive(...))`:
   - py L170-174: `if not all(sp.ask(sp.Q.positive(sp.simplify(case))) for case in positive_product_cases): ...`.
   Mathematica computes the full semialgebraic regions and proves region equivalence:
   - wl L246-273: `radicandRegion = Reduce[turnRadicand >= 0 && ...]`; `rootExistRegion = Reduce[Exists[tauRoot, ... ordinaryModel[...]==0] && ...]`; `expectTrue[..., Equivalent[radicandRegion, criterionRegion]]`.
   `ask`/case-substitution vs `Reduce`/`Exists` region equivalence → independent (and the `.wl` is the stronger check here).

The shared step in Section V/M9-M11 is the final `Series` expansion — but the object being expanded is itself extracted differently (posited radical in py vs Solve-derived, limit-selected root in wl), so this is not the "both Series-expand the same closed form" port pattern the prompt warns against. The "each CAS runs its own simplifier" defense is not invoked; the methods genuinely differ at the extraction step.

## Engine cross-check

Both engines confirm identical mathematics. SymPy prints residual `= 0` for every identity (output L71-72, 90-93, 111-114, 119-134, 141-151) and the closed forms `tau_quad`, `tau_log2`, the turning radicand cases, and the cubic gap coefficient `eps³(L0e²+3L1e)/(6 L0e³)`. Mathematica prints `= 0` / `True` / PASS for the corresponding M1-M11 (output L9-72), including `M11 gap coefficient eps^3 = 0` after subtracting `(l0s²+3 l1s)/(6 l0s³)`. The agreement-theorem cubic coefficient `(L0²+3L1)/(6L0³)` matches between engines and matches the notes (§8.3) and appendix. No `engine_disagreement`.

## Verdict justification

The math holds up under attack in both engines and matches the paper card, the notes (§2-§8), and the appendix (eqs. app-part06-free-log-coords through app-part06-turning-predictors) exactly. Attacks tried and failed: (a) checked whether the predictor residual checks are tautological — they are not, because the predictor (posited or Solve-derived) is substituted into a closure model defined independently; (b) checked the negative-slope branch sign handling (py L101-104 `Phi1a_n - sqrt`; wl via Solve+select) — correct; (c) checked the turning-point reality logic for a `missing_branch` (both `(1-Φ0)Φ2>0` and `<0` cases are exercised in both engines, py L160-175 / wl M7) — covered; (d) checked the cubic agreement coefficient against the notes — exact match; (e) checked symbol domains (all Reals with the documented positivity on `Φ0,Φ1,L0` matching the notes' `Φ0>0`, `Φ1≠0`, `L0≠0` hypotheses) — justified. The single finding is documentary: the card still says "Mathematica audit: none yet" though the `.wl` exists and passes (F1, paper_misalignment, user-resolved). Verdict `findings` with `needs_user_resolution: true`; the underlying physics verification is sound and dual-engine-independent.

## Value Reconciliation (pass-2 augmentation)

This stage emits no numeric constants — every deliverable is a symbolic closed form. The reconciliation is therefore over the boxed/closed-form results.

| value (symbolic) | source (py / wl + output) | .tex / .md location | status |
|---|---|---|---|
| `Phi2 = sᵀ H s`, `Phi1 = sᵀ g` | py L63-64 / out L9-18; wl L147-148 | notes L116-122 (boxed), appendix L682-687 | MATCH |
| `L1 = Phi2/Phi0 - Phi1²/Phi0²` | py L77 / out L71; wl M1 / out L9 | notes L166-171; appendix L691-692 | MATCH |
| bridge `Phi2 = Phi0(L1+L0²)` | py L78 / out L72; wl M2 / out L11 | notes L174-179; appendix L695-697 | MATCH |
| `Delta_aff = Phi1²-2 Phi2(Phi0-1)` | py L86 / out L77-79 | notes L205-209; appendix L702 | MATCH |
| `tau_quad = 2(1-Phi0)/(Phi1+sgn√Δ)` | py L88 / out L84-89; wl M3 (Solve) | notes L226-230; appendix L708-710 | MATCH |
| `tau_aff = (1-Phi0)/Phi1` (zero-curv limit) | py L87 / out L80-83; wl M4 | notes L242-245; appendix L670 | MATCH |
| `Delta_log = L0²-2 L1 ln Phi0` | py L119 / out L98-100 | notes L212-216; appendix L704 | MATCH |
| `tau_log2 = -2 ln Phi0/(L0+sgn√Δ)` | py L121 / out L105-110; wl M5 (Solve) | notes L253-258; appendix L712-715 | MATCH |
| `tau_log = -ln Phi0/L0` (zero-curv limit) | py L120 / out L101-104; wl M6 | notes L270-273; appendix L672 | MATCH |
| turning criterion `(1-Phi0)Phi2>0` | py L153,170-175 / out L119-127; wl M7 / out L41-46 | notes L296-298; appendix L717-720 | MATCH |
| `tau_± = ±√(2(1-Phi0)/Phi2)` | py L177 / out L128; wl M7 / out L47-50 | notes L302-307; appendix L723-725 | MATCH |
| tangency `Δ = ½ Phi2 τ²` | py L183 / out L129-134; wl M8 / out L55 | notes L320-323 | MATCH |
| ordinary corr `τ_quad=τ_aff-(Phi2/2Phi1)τ_aff²` | py L205-210 / out L139-141; wl M9 / out L61 | notes L347-351; §8.1 | MATCH |
| log corr `τ_log2=τ_log-(L1/2L0)τ_log²` | py L212-217 / out L142-144; wl M10 / out L63 | notes L363-367; §8.2 | MATCH |
| agreement cubic `(L0²+3L1)/(6L0³) ε³` | py L219-226 / out L145-151; wl M11 / out L71 | notes L380-385 (boxed); §8.3 | MATCH |

INTERNAL scaffolding (no prose expected, no finding): `expect_zero`/`expectZero`/`expectTrue`/`pass`/`fail` harness, residual-equals-0 print lines, `selectByLimit`/`nearZeroBranch` branch-selection helpers, `$Assumptions`/`baseAssumptions`, `radicand_tp`/`turnRadicand` intermediates, `apos`/`bpos`/`a*`/`*Abs` case-substitution symbols, `eps`-series truncation order, `Coefficient` extraction, and the PASS/AUDIT-PASSED banners.

reconciliation: complete; 15 deliverable values checked, 0 misaligned. (The single F1 finding is a Verification-pointer staleness on the card, not a result-value mismatch; it is recorded in `## Findings` and routes to the user, but no emitted *result value* is misaligned.)

## Self-test notes

Checked variable-independence: `D[Log[chiLocal], vars[[i]], vars[[j]]]` (wl L150) genuinely depends on the `vars` differentiated, since `chiLocal` (wl L143) is an explicit function of `vars` — the second derivatives are not identically zero (confirmed by the nonzero `logHessianByD` feeding M1). Checked parity/symmetry: n/a (no unbounded-domain integrals). Trivial-case pre-check: substituting `Phi2→0` collapses `tau_quad` to `tau_aff` (py A4 / wl M4) and `L1→0` collapses `tau_log2` to `tau_log` (A6/M6), both confirmed zero, so the predictor checks are not vacuously satisfied. The lone finding is paper_misalignment (user-resolved), so the directive carries a `## Resolve before fix_loop` block and prescribes no Codex edit.
