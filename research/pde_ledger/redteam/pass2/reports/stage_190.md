---
unit_id: 190
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage190_direct_defect_vs_dressing_split.md]
  paper_appendix: present
---

# Audit unit 190 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_190.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage190_direct_defect_vs_dressing_split.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/blocks at lines 95, 111, 770-879, 1466, 1557)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` reads verbatim: "Separates direct transfer-shape drift from selected-branch dressing residual and records scalar no-go filters." (`stage_190.tex:15`). The Part-V appendix row repeats this verbatim (`stage_appendix_part05.tex:111`). The notes expand this into five concrete deliverables: (1) an exact direct-defect / dressing split of the coherent weak-axisymmetric problem; (2) an exact support-blindness theorem (`∂T²/∂ζ = ∂R_target/∂ζ = ∂N_*/∂ζ = 0`); (3) the branch-adapted coordinates `Σ_tr, Σ_nt, Σ_η` with the triangular normal form `Θ₁=−C_tr Σ_tr`, `Ξ₁=A_tr Σ_tr + Σ_nt`, `R₁+Ξ₁=−(ε_η/(1−ε_η))Σ_η`; (4) exact inverse reconstruction formulas and the triple-rigidity theorem; and (5) the no-go filter that a pure grouped real `P2` anisotropy cannot linearly feed the scalar off-bundle slippages (first invariant is quadratic, `I[x,x]=(7/10)ε²(x⁽¹⁾)²`). The notes also pin the microscopic slippage ledger (`Σ_χ, Σ_δ, Σ_Z, Σ_ε, Σ_η`) and the explicit `Ξ₁` direct-defect law. The card's body is terse (no equations); the notes and the Part-V appendix (eqs. app-part05-Xi1-slippage-law, -Sigma-tr-def, -Sigma-nt-def, -triangular-normal-form, -Ctr-Atr-defs) are the authoritative carriers of the closed forms.

## What the script claims to verify

Both engines verify, section-by-section, every notes deliverable: (I) support-blindness — they rebuild the support-loaded mass packet `M_tr = M_mix + M_supp`, recover `T², R_target, N_*, E=1−ε_η` exactly from it, and confirm the three `d(log)/dζ` derivatives vanish, plus a negative-control "spoiled" packet whose `dζ`-derivative is nonzero; (II) the five microscopic slippage laws via a logarithmic drift operator on the defining monomials; (III) the `Ξ₁ = A_tr Σ_tr + Σ_nt` decomposition with `A_tr, C_tr, Σ_tr, Σ_nt`; (IV) the triangular compiler matrix, its determinant, the three zero-derivative block-triangularity checks, the three inverse-reconstruction formulas, and the tracking-rigid / cancellation / selected-branch rigidity theorems; and (V) the grouped-`P2` projection giving `x̄=x₀`, `b_x=3a_x`, `I[x,x]=(7/10)ε²x₁²`, and vanishing linear terms. No literal answer is hardcoded; `Λ₀` is a free symbol (the script never pins its `27π²Gc_s⁵/20a⁵c⁵` value), and all results are derived in-script.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| (2) support-blindness `∂T²/∂ζ=∂R_target/∂ζ=∂N_*/∂ζ=0` | py 72-81 / wl 94-104 (reconstruct + d log/dζ) | match |
| `E = R_target T²/Λ₀ = 1−ε_η` | py 75 / wl 101 | match |
| microscopic slippage ledger `Σ_χ,Σ_δ,Σ_Z,Σ_ε,Σ_η` | py 144-148 / wl 142-152 | match |
| direct-defect law `Ξ₁` (eq app-part05-Xi1) | py 160-168, 194 / wl 158-184 | match |
| `Σ_tr` def, `A_tr`, `C_tr` (eq app-part05-Ctr-Atr) | py 169-173 / wl 165, 170-172 | match |
| `Σ_nt` def (eq app-part05-Sigma-nt) | py 174-182 / wl 173-179 | match |
| triangular normal form `(Θ₁,Ξ₁,R₁+Ξ₁)` | py 202-220 / wl 189-205 | match |
| det = `C_tr ε_η/(1−ε_η) > 0` | py 217 / wl 202 | match |
| inverse reconstruction `Σ_tr,Σ_nt,Σ_η` | py 222-232 / wl 207-234 | match |
| triple-rigidity / cancellation theorems | py 230-232 / wl 232-234 | match |
| (5) `P2` no-go: `x̄=x₀`, `b_x=3a_x`, `I=(7/10)ε²` | py 260-267 / wl 261-273 | match |

Every paper-side deliverable has a faithful, non-tautological script-side check on BOTH engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 72-74 | `expect_zero(loaded − direct)` | support-blind reconstruction (2) | yes |
| A2 | sympy | 76-81 | `expect_zero(d log/dζ)` | support-blindness (2) | yes |
| A3 | sympy | 89-90 | spoiled `dlnR ≠ 0` (negative control) | guards (2) | yes |
| A4 | sympy | 144-148 | `expect_zero(Σ − target)` | slippage ledger | yes |
| A5 | sympy | 194 | `expect_zero(Ξ − (A_tr Σ_tr+Σ_nt))` | direct-defect law (3) | yes |
| A6 | sympy | 218-220 | `expect_zero(block derivatives)` | triangularity (3) | yes |
| A7 | sympy | 226-232 | `expect_zero(inverse − orig)` | inversion (4) | yes |
| A8 | sympy | 260-267 | `expect_zero(P2 invariants)` | no-go (5) | yes |
| B1 | math | 94-104 | `expectZero(loaded−direct, d log/dζ)` | support-blindness (2) | yes |
| B2 | math | 110-117 | `expectNonzero(spoiled)+numeric witness` | guards (2) | yes |
| B3 | math | 142-152 | `expectZero(logEuler − target)` | slippage ledger | yes |
| B4 | math | 181-184 | `expectZero(projected A_tr, Σ_nt, no-Σ_χ)` | direct-defect law (3) | yes |
| B5 | math | 202-205 | `expectZero(det formula, block derivs)` | triangularity (3) | yes |
| B6 | math | 207-234 | `expectZero(LinearSolve inverse − recovered)` | inversion (4) | yes |
| B7 | math | 261-273 | `expectZero(P2 invariants + Series)` | no-go (5) | yes |

All rows are anchored (non-tautological). Each maps to a specific paper/notes deliverable; none are orphaned.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.txt:3` and `:148`

**What's wrong:**
The committed SymPy output banner reads `STAGE 173 — …` (line 3) and `STAGE 173 LEDGER` (line 148), whereas the current script prints `STAGE 190` (`...sympy_audit.py:35`, `:269`). The `.txt` mtime is `2026-06-01 11:35:23`; the `.py` mtime is `2026-06-03 15:59:11` — the output predates the script. `git log -p` confirms the ONLY change to the `.py` since the output was captured (commit e2a4780, "numbering reconciliation Phase 1") was the two banner relabels `173 → 190`; all math lines are byte-identical. So the captured numerical/symbolic content still agrees with what the current script produces — only the banner is stale.

**Why this matters:**
Cosmetic mismatch only; the math is unaffected. Filed because the freshness signal is real and the verifier's re-run will regenerate it.

**Required change:**
None for Codex. The orchestrator's independent re-run (`exec-sympy 190`) refreshes the `.txt`; the new banner will read `STAGE 190`.

**Verification:**
After re-run, `...sympy_audit.txt:3` should read `STAGE 190 — …` and `:148` `STAGE 190 LEDGER`; all `= 0` result lines remain identical.

### F2 — paper_misalignment (subtype: paper_missing_script_claim)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_190.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl` (whole file)

**What's wrong:**
The card's `\stagefield{Verification}` states: "SymPy audit: `...sympy_audit.py`.  **Mathematica audit: none yet.**" (`stage_190.tex:11`). But a fully populated Mathematica audit `...mathematica_audit.wl` exists (10.5 KB, 283 lines, mtime 2026-06-01) and produces an all-PASS output `...mathematica_audit.txt`. The card's prose is stale relative to the pass-1 dual-engine retrofit: the second engine now exists but the card still says it does not. This is a documentation-side discrepancy, not a math defect.

**Why this matters:**
The card under-reports the verification coverage. A reader citing this stage would believe only single-engine SymPy backing exists. Direction of resolution (update the card to point at the `.wl`, vs. retire the `.wl`) is a user/paper-owner call — Codex does not edit paper/. Per the prompt, this routes to the user gate, not to a Codex script edit.

**Required change:**
See `## Resolve before fix_loop` in the directive. No script edit.

**Verification:**
After user resolution: if (a), `stage_190.tex:11` should read "Mathematica audit: `\StageFile{mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl}`."; no script change.

## Independent-derivation check (Mathematica)

This is the heart of the V.3 retrofit audit. I compared each deliverable's METHOD in `.wl` vs `.py`. The `.wl` is **GENUINELY INDEPENDENT** — it uses a different operation for the same target in four of the five sections, and adds checks the `.py` omits.

1. **Slippage ledger (II).** `.py:123-125` substitutes each primitive `v → v·exp(s·v₁)` then `sp.diff(sp.log(expr_s), s).subs(s, 0)` — an exponential one-parameter drift differentiated at `s=0`. `.wl:127-134` instead applies a **logarithmic Euler operator** `logEuler[m] = Total[MapThread[#2·#1·D[Log[m], #1] &, {vars, drifts}]]` — the direct sum `Σ drift_i · v_i · ∂log/∂v_i`. Different mathematical route (Euler operator vs. exp-substitution-then-derivative) to the same five laws. INDEPENDENT.

2. **Direct-defect decomposition (III).** `.py:160-194` posits closed forms for `Xi_direct, Sigma_tr, Atr, Sigma_nt` and checks the single combination identity `Xi_direct − (Atr·Sigma_tr + Sigma_nt) == 0`. `.wl:166-168` **derives** `A_tr` by coefficient extraction `atrProjected = Coefficient[xiDirect, sigmaChi]/Coefficient[trackingCoordinate, sigmaChi]` and `sigmaNtProjected = xiDirect − atrProjected·trackingCoordinate`, then checks `Coefficient[sigmaNtProjected, sigmaChi] == 0` (`wl:182`) — a no-`Σ_χ`-residual test absent from the `.py`. Projection-extraction vs. posited-form-check: INDEPENDENT and strictly stronger on this section.

3. **Triangular compiler inversion (IV).** `.py:222-224` hand-writes the inverse closed forms `Sigma_tr_rec = −((1+χ)(1+δ)(1+χ+δ)/(χδ))·Theta`, etc., and checks each equals the original. `.wl:207` instead inverts via native linear algebra `inverseByLinearSolve = LinearSolve[compiler, {thetaVar, xiVar, dressVar}]` and checks each component (`wl:213-219`). Solver-based inverse vs. hand-written inverse: INDEPENDENT.

4. **Support-blindness (I).** Same construction on both sides (build `M_mix+M_supp`, recover transfer, differentiate in `ζ`) — here the routes are parallel. BUT the negative controls differ: `.py:84-90` spoils with `Msupp + bad·ζ·Mmix` and asserts the symbolic derivative `≠ 0`; `.wl:106-117` spoils with `loadedMass + ζ·mixMass` and pins an **exact numeric witness** `(spoiledResidual /. {eps→1/3, zeta→1/2}) + 46/133 == 0`. Different perturbation and a numeric (vs. symbolic) falsification witness.

5. **`P2` no-go (V).** `.py:241-267` builds `x20,x21,x22` componentwise and `xbar, ax, bx` from explicit weighted sums; linear-term checks via `sp.diff(...).subs(epsax,0)`. `.wl:240-273` builds the lane as `x0·{1,1,1}+epsAx·x1·{1,1/2,−1}`, applies a **projector matrix** `projectors.xLane` to extract `(x̄,a_x,b_x)`, forms `I` via `DiagonalMatrix[{4,4/5}]` quadratic form, and checks linear terms via `Coefficient[Normal[Series[...]], epsAx, 1]`. Matrix-projector + `Series` vs. componentwise + `diff`: INDEPENDENT.

Verdict: **GENUINELY INDEPENDENT** across all five sections, with consistently distinct operations (Euler operator, `Coefficient` projection, `LinearSolve`, projector matrices, `Series`) and several checks the `.py` lacks. This is NOT a transliteration. The first-pass CLEAN disposition holds up under fresh adversarial re-audit.

## Engine cross-check

Both engines exit clean (SymPy raises no AssertionError; Mathematica reaches `Exit[0]` with all `PASS:`). The shared deliverables agree: `det(compiler)` is `−χ₀δ_U ε_η/((1+χ₀)(1+δ_U)(ε_η−1)(1+χ₀+δ_U))` in SymPy (`...sympy_audit.txt:121-123`) and `−((chi·deltaU·epsEta)/((1+chi)(1+deltaU)(1+chi+deltaU)(−1+epsEta)))` in Mathematica (`...mathematica_audit.txt:60`) — identical. `E = 1−ε_η`, `A_tr = 2χ₀/((1+χ₀)(1+δ_U))`, `C_tr = χ₀δ_U/((1+χ₀)(1+δ_U)(1+χ₀+δ_U))`, and `I[x,x]−(7/10)ε²x₁² = 0` match on both. No `engine_disagreement`.

## Verdict justification

`findings` — but both findings are documentation-side and low-severity: (F1) a cosmetic stale-banner in the SymPy `.txt` (only `STAGE 173 → 190`, math byte-identical, refreshed by the verifier's re-run), and (F2) a `paper_misalignment` where the card says "Mathematica audit: none yet" while a complete, independent `.wl` now exists (routes to the user gate). The math itself is clean: I attacked the support-blindness route (tried to find a hidden `ζ` dependence that survives — the spoiled negative controls on BOTH engines confirm the route is falsifiable, so the `dζ=0` results are substantive not trivial); I checked the slippage-derivative directions (every `diff`/Euler target genuinely depends on the differentiation variable — no identically-zero-derivative trap); I checked the `P2` parity/quadratic claim (`b_x=3a_x` and `I=(7/10)ε²` reduce correctly from the `{1,1/2,−1}` signature on both engines); and I confirmed the `.wl` is an independent re-derivation, not a port. The card+notes+appendix were read and the script's verified claims match the paper's stated deliverables exactly.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable the scripts emit is symbolic except the rational coefficients `7/10`, `b_x=3a_x`, `46/133` (witness), and the grouped metric/weights. `Λ₀` is left as a FREE symbol — the script never pins its `27π²Gc_s⁵/20a⁵c⁵` numeric, so there is no script-side `Λ₀` value to reconcile (the notes/appendix carry that constant; the script correctly treats it abstractly). The card body is terse and carries no equations by design; the NOTES (`...stage190...md`) and the Part-V APPENDIX are the natural carriers, and every deliverable lives there.

| value (deliverable) | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `E = R_target T²/Λ₀ = 1−ε_η` | py:75 / wl:101; txt:32, :17 | md:106 (`𝔈:=…=1−ε_η`) | MATCH |
| `Σ_χ=γ₁+c₁−κ_U` | py:144 / wl:142; txt:49 | md:247; appx:792 region | MATCH |
| `Σ_δ=τ₁−κ_U` | py:145 / wl:143; txt:51 | md:248 | MATCH |
| `Σ_η=2c₁−κ_U−κ_η` | py:146 / wl:144; txt:53 | md:260 | MATCH |
| `Σ_ε=2γ₁+2λ₁−κ_U−κ_W` | py:147 / wl:146; txt:55 | md:255 | MATCH |
| `Σ_Z=2λ₁+μ₁−κ_η−2κ_W` | py:148 / wl:150; txt:57 | md:253 | MATCH |
| `Ξ₁` direct-defect law | py:160-168; txt:67 | md:283-293; appx:786-800 | MATCH |
| `Σ_tr=(1+χ₀)Σ_δ+(1+δ_U)Σ_χ` | py:169 / wl:165; txt:82 | md:301; appx:803-804 | MATCH |
| `A_tr=2χ₀/((1+χ₀)(1+δ_U))` | py:170 / wl:170; txt:96 | md:335; appx:861 | MATCH |
| `C_tr=χ₀δ_U/((1+χ₀)(1+δ_U)(1+χ₀+δ_U))` | py:171-173 / wl:171-172; txt:100 | md:310-311; appx:856-858 | MATCH |
| `Σ_nt` def | py:174-182 / wl:173-179; txt:84 | md:319-328; appx:824-837 | MATCH |
| triangular normal form `(Θ₁,Ξ₁,R₁+Ξ₁)` | py:207-213 / wl:189-198; txt:108 | md:383-401; appx:845-851 | MATCH |
| `det = C_tr ε_η/(1−ε_η)` | py:217 / wl:202; txt:120 | md:404-406 | MATCH |
| `Σ_tr⁻¹ = −((1+χ₀)(1+δ_U)(1+χ₀+δ_U)/(χ₀δ_U))Θ₁` | py:222 / wl:208; txt:127 | md:434-437 | MATCH |
| `A_tr/C_tr = 2(1+χ₀+δ_U)/δ_U` | py:227 / wl:215; txt:128 | md:456-459 | MATCH |
| `Σ_η⁻¹ = −((1−ε_η)/ε_η)(R₁+Ξ₁)` | py:224 / wl:211; txt:130 | md:484-487 | MATCH |
| `b_x = 3 a_x` | py:262 / wl:263; txt:140 | md:576 | MATCH |
| `I[x,x] = (7/10)ε²x₁²` | py:265 / wl:266; txt:143 | md:576 | MATCH |
| support-blindness `∂T²/∂ζ=∂R_target/∂ζ=∂N_*/∂ζ=0` | py:76-81 / wl:102-104; txt:37-39 | md:152-180; appx:817 | MATCH |

INTERNAL (scaffolding, no prose expected): `M_mix`, `M_supp`, `M_tr`, `S(ζ;ε)` support factor, `R_tr`, `B_*` (intermediate composites); the `bad`/`badPacket` spoiled-packet symbols and the `46/133` numeric falsification witness; `T²`, `R_target`, `N_*` intermediate printed forms (these ARE in the notes as boxed eqs md:65-79, but printed here only as reconstruction targets); pass/fail flags, residual `= 0` lines.

reconciliation: complete; 19 deliverable values checked, 0 misaligned. (F2 is a verification-prose mismatch about the `.wl`'s existence, not a value mismatch; no `value_mismatch` raised.)

## Self-test notes

Checked the five traps. (1) Variable-independence: every `sp.diff`/`D[]`/Euler target genuinely depends on its differentiation variable — the support-blindness `d log/dζ` acts on the support-LOADED expressions (which DO contain `ζ` before cancellation), so the `= 0` result is the substantive cancellation, not a trivially-`ζ`-free derivative; the spoiled negative controls (py:84-90, wl:106-117) confirm a non-blind packet yields nonzero, so the route can fail. (2) Parity/quadratic: the `P2` no-go `b_x=3a_x`, `I=(7/10)ε²` reduce correctly from the `{1,1/2,−1}` signature; the linear-`ε` term genuinely vanishes (even invariant). (3) Trivial-case: substituting the slippage targets confirms each `Σ − target → 0`; the determinant and inverse reconstruct to literal zeros. (4) No missing-script directive (both engines present). (5) Paper round-trip: no script-side fix prescribed (F1 is verifier-refresh, F2 is user-gated paper prose), so no new `paper_misalignment` introduced. Numbering: card `\stagefield{Purpose}` self-labels "Stage~190"/"Stage 190" (the +17 class: 190, NOT 207 — correct); notes cross-refs to Stages 188/189/191/168 are content refs and consistent; the stale `.txt` "STAGE 173" banner is the known pre-renumber label folded into F1. No numbering fix needed.
