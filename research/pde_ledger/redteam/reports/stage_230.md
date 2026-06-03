---
unit_id: 230
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 230 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_230.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows for stage 230 read: line 72 status row; lines 806-880 dynamic-window-compiler section)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` states verbatim: "Static-first theorem: the static $\Xi_1$ budget is the first kill condition; the dynamic corridor can be checked only after the static transported ceiling is satisfied." The card's purpose line says stage 230 "translates the selected-branch classifier into a dynamic-window decision rule," with inputs the Stage 229 classifier data, the rigid-split shares, and the survival-window inequality. The notes (authoritative on intent) enumerate five deliverables: (1) the exact share weights `w_num = R_ND/(1+R_ND)`, `w_den = 1/(1+R_ND)`; (2) the affine selected-branch slope compiler `S_±(R_ND) = (R_ND s_±^num + s_±^den)/(1+R_ND)`; (3) the exact sign-flip threshold `R_* = s_-^den/(-s_-^num) ≈ 1.229255438463336`; (4) the onset threshold `delta_*^dyn = 8/(9 R_*) ≈ 0.723111617875019`; (5) the static-first theorem `B_dyn^both > B_stat^both` and `B_dyn^nonempty > B_stat^nonempty` over the whole classifier half-line. The appendix (lines 839-846) restates the share weights, the dynamic sign theorem, and "the selected-branch dynamic ceilings are everywhere weaker than the universal transported static ceilings." Status is mixed (exact-closure + numerical), checkpoint: false.

## What the script claims to verify

The SymPy script verifies, in five blocks: (1) the exact Stage-229 classifier `R_ND` plus its hand-derivative form and the `8/(9δ)` onset and `xi->1^-` zero limit; (2) the rigid-split share weights, the affine compiler `S_±(R)`, and that both `dS_±/dR < 0` (monotone decreasing); (3) the sign-flip threshold `R_* = s_-^den/(-s_-^num)` and onset threshold `delta_*^dyn = 8/(9 R_*)`, each checked numerically against the notes' stated values, plus representative sign checks of `S_-` below/above `R_*` and a `delta=1` slice staying below `R_*`; (4) the dynamic ceilings `B_both`, `B_nonempty` in `|eps Xi_1|`, their endpoint/limit values, and four sample classifier points; (5) the universal static budgets and the two global strict inequalities `inf B_dyn > B_stat`. The four per-unit rigid slopes, the three `R_Q` figures, and the two static budgets are carried as decimal literals.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script check | Status |
|---|---|---|
| (1) share weights `w_num`, `w_den` | lines 85-86 define `R/(1+R)`, `1/(1+R)` | match |
| (2) affine compiler `S_±(R)` | lines 87-93 (assert vs `expected_S_±`) + monotonicity 95-98 | match |
| (3) sign-flip `R_* = s_-^den/(-s_-^num) ≈ 1.2292554...` | lines 111-112, plus `S_-(R_*)=0` (118) and sign checks (124-127) | match |
| (4) onset `delta_*^dyn = 8/(9 R_*) ≈ 0.7231116...` | lines 114-115 (numeric) + line 119 (roundtrip, tautological — see F2) | match (value) |
| (5) static-first `B_dyn > B_stat` everywhere | lines 210-211 check the infima vs budgets; "everywhere" follows from monotonicity (95-98) | match |
| R_Q dynamic figures / margins (notes §5) | lines 142-150 carry `30.199907..., 36.171186..., 21.854566...` and check `ell_±` | match |

All carried literals trace to genuine upstream derivations: Stage 228 derives the rigid R_Q slopes (`-0.52346582, 0.71358484, -0.35245541, -0.23169484` at stage228 lines 418-421), the base R_Q `[30.1999075602499, 36.1711864832695]` (stage228 line 368), the requirement threshold `21.854566296358396` (stage228 line 432), the positive-Xi_1 unit norms `1.73611234967676` / `0.692932151812037` (stage228 lines 294/298), and the static budgets `0.367930328492646` / `0.737619063660757` (stage228 lines 441-442). Stage 229 derives `R_ND` and the `8/(9δ)` onset. The stage-230 per-unit slopes (`-0.301516..., 0.411024..., -0.508643..., -0.334368...`) equal the stage-228 rigid R_Q slopes divided by the stage-228 Xi_1 unit norms; the division is performed in the notes prose, not re-derived in stage 230, but each input is individually anchored upstream. No fabricated literal. No R_Q typo of the 222/223 kind: the three R_Q figures in the script (lines 142-144) exactly match the notes (§5, lines 342-347).

Note on stage numbering: the notes prose uses the OLD numbering ("Stage 247" in its title/body, "Stage 245" = rigid, "Stage 246" = classifier), whereas the paper card, appendix, and script all use the CURRENT numbering (230 / 228 / 229). The mathematical content is identical across all four documents; only the notes' stage-number labels lag. This is a prose-only artifact in a document Codex may not edit, with zero effect on what the script verifies, so it is recorded here as an observation, not a paper_misalignment requiring user resolution. The substantive claims align exactly. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 61 | `simplify(dR_dxi - expected_dR_dxi) == 0` | classifier monotonicity (notes §1.2) | yes (autodiff vs notes hand-form) |
| A2 | sympy | 64 | `simplify(onset - 8/(9δ)) == 0` | onset form | yes |
| A3 | sympy | 66 | `soft_limit == 0` | softening limit | yes |
| A4 | sympy | 82-83 | sign guard on the four carried slopes | input-data sanity | yes (guards carry-forward signs) |
| A5 | sympy | 92-93 | `simplify(S_± - expected_S_±) == 0` | claim (2) compiler | partial (two equivalent rearrangements; see note) |
| A6 | sympy | 97-98 | `dS_± < 0` | monotonicity → "everywhere" | yes (depends on negative constant) |
| A7 | sympy | 112 | `assert_close(R_*, 1.229255438463336, 1e-15)` | claim (3) | yes (value vs notes) |
| A8 | sympy | 115 | `assert_close(delta_*, 0.723111617875019, 1e-15)` | claim (4) | yes (value vs notes) |
| A9 | sympy | 118 | `simplify(S_-(R_*)) == 0` | sign-flip definition | yes (S_- vanishes exactly at R_*) |
| A10 | sympy | 119 | `simplify(onset(δ→δ_*) - R_*) == 0` | claim (4) roundtrip | NO — tautological (F2) |
| A11 | sympy | 122-127 | `S_+(0)<0`, `S_-(1/2)>0`, `S_-(2)<0` | sign theorem | yes |
| A12 | sympy | 130-133 | `R_ND(δ=1,·) < R_*` on probe grid | denominator-like slice | yes |
| A13 | sympy | 149-150 | `assert_close(ell_±, ...)` | dynamic margins | yes |
| A14 | sympy | 161-162 | `B_both(0) ≈ 1.671...`, `B_nonempty(0)=inf` | endpoint ceiling | yes |
| A15 | sympy | 169-170 | `B_both(inf) ≈ 0.9672...`, `B_nonempty(inf) ≈ 0.9905...` | infimum ceilings | yes |
| A16 | sympy | 189-196 | four sample points `S_±, B_both, B_nonempty` | sample table (notes §7) | yes |
| A17 | sympy | 210-211 | `B_both_inf > B_stat_both`, `B_nonempty_inf > B_stat_nonempty` | claim (5) static-first | yes (infimum vs budget; "everywhere" via A6) |

Note on A5: `S_± = w_num*s_num + w_den*s_den` vs `expected_S_± = (R*s_num + s_den)/(1+R)` are two arrangements of the same construction, so A5 alone is near-tautological; however the *meaningful* anchoring of the compiler comes from A6 (monotone sign), A11 (representative signs), A14-A16 (numeric endpoints/samples vs the notes' independently quoted figures), so the compiler claim is substantively exercised overall.

## Findings

### F1 — missing_verification_script (subtype: missing_mathematica)

**Severity:** high
**Files:**
- `(missing)` — no `.wl` exists; paper card line 11 admits "Mathematica audit: none yet."

**What's wrong:**
Stage 230 is checkpoint:false but NOT status-only (`is_status_only_candidate: False`). The dual-engine rule requires a Mathematica audit wherever Mathematica CAN independently verify the stage. Every operation here is a native Mathematica primitive: rational-function `Simplify`/`Together`, symbolic `D` for `dR_ND/dxi` and `dS_±/dR`, `Limit` for the `xi->1^-` and `R->oo` endpoints, sign determination via `Reduce`/`Resolve`, `Log` for the dynamic margins, and `N` for the numeric comparisons. There is no obstacle to a fully independent second engine. The SymPy script is the only verifier.

**Why this matters:**
The static-first theorem (the card's bottom-line Output) currently rests on a single engine. The dual-engine policy exists precisely so a transliteration or a SymPy-specific quirk (e.g. an `assume`-dependent `simplify`, the `dS<0` sign resolution, or the `Limit` at `R->oo`) cannot silently pass unchallenged.

**Required change:**
Add an independent Mathematica audit `.wl` (full claim manifest and anti-transliteration guard in the directive). It must derive the compiler and thresholds from the physical premises using native primitives and a different decomposition than the SymPy script — not echo the `.py` algebra.

**Verification:**
`redteam exec-mathematica 230` produces a fresh `.txt`, the new `.wl` exists at the mandated path below, and it exits 0 with all in-file checks passing; its `R_*`, `delta_*`, endpoint ceilings, and the two static-first inequalities must agree with the SymPy values.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/.../stage230_..._sympy_audit.py:119`

**What's wrong:**
Line 114 defines `delta_dyn_star = sp.simplify(8 / (9 * R_star))`. Line 119 then asserts `sp.simplify(onset.subs(delta, delta_dyn_star) - R_star) == 0`, where `onset = 8/(9*delta)` (line 63). Substituting `delta -> 8/(9*R_star)` into `8/(9*delta)` returns `R_star` by pure construction: `8/(9 * 8/(9 R_star)) = R_star`. The subtraction is identically zero no matter what `R_star` is. The assertion cannot fail and exercises no physics — it merely confirms `8/(9·x)` is its own inverse map. The genuine onset-threshold content (notes §4.2: "requiring onset `8/(9δ) ≤ R_*` gives `δ ≥ 8/(9 R_*)`") is already captured numerically at line 115, which anchors `delta_dyn_star ≈ 0.723111617875019` against the notes' stated value.

**Why this matters:**
A reader may believe line 119 independently verifies the onset-threshold derivation, when in fact it is a definitional round-trip. It is harmless but misleading scaffolding; the surrounding A8/A7 checks already carry the real anchoring.

**Required change:**
Replace the tautological round-trip with the actual onset-threshold statement: assert that the equation `onset(delta) = R_star` (i.e. `8/(9*delta) - R_star == 0`) is solved by `delta = delta_dyn_star`, by solving for `delta` independently and checking the solution equals `delta_dyn_star` — OR assert the inequality direction `onset.subs(delta, delta_dyn_star) <= R_star` together with strict decrease, which is the physical content. See directive F2 for the precise minimal edit.

**Verification:**
Line 119 (or its replacement) must derive `delta` from `8/(9*delta) == R_star` rather than substitute the pre-defined `delta_dyn_star`; the script still exits 0.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` exists (F1).

## Engine cross-check

Not applicable — single engine present (F1).

## Verdict justification

The SymPy script is mathematically faithful to the paper card, the appendix, and the notes: the classifier, share weights, affine slope compiler, sign-flip threshold `R_* ≈ 1.2292554`, onset threshold `delta_* ≈ 0.7231116`, the dynamic ceilings, and the two static-first strict inequalities all match, and every carried literal traces to a genuine Stage-228/229 derivation (no fabricated value, no R_Q typo of the 222/223 kind). The "everywhere weaker" claim is correctly reduced to the `R->oo` infimum via the proven monotonicity `dS_±/dR < 0`, and that reduction is sound for both the `both` and `nonempty` ceilings. Attacks that failed: I checked the `dS_± < 0` symbolic comparison (resolves to True because `R` is nonnegative real and the constant numerator is negative — non-tautological); the `S_-(R_*) = 0` check (vanishes exactly by construction of `R_*`, genuine); the `dynamic_ceilings` branch logic at and above `R_*` (correct: at `R_*`, `S_-=0` routes to `improving_present` → nonempty=inf); and the infimum reduction for `B_nonempty` over the finite region (monotone-decreasing max of two decreasing functions → infimum at `R->oo`, sound). Two findings remain: a high-severity `missing_mathematica` (dual-engine required and possible) and a low-severity `tautological_check` at line 119 (definitional round-trip of the onset inverse). Verdict: `findings`, no stop-cold. Both are script-side and Codex-applicable; no paper_misalignment, so no user resolution is needed.

## Self-test notes

Variable-independence trap: the two new/edited derivative-style checks (the existing `dR_dxi` at line 54 and `dS_±` at 95-96) differentiate w.r.t. variables (`xi`, `R`) the expressions genuinely depend on, so no identically-zero-derivative trap; for the F2 fix I prescribe `solve(8/(9*delta) - R_star, delta)` (an equation solve, not a derivative), avoiding the zero-derivative pitfall entirely. Symmetry/parity trap: no unbounded-domain integrals in this unit, so N/A. Trivial-case pre-check: the F2 replacement, substituting the known `R_star`, yields `delta = 8/(9 R_star)` which equals `delta_dyn_star` (≈0.7231116) — a nonzero literal, confirming the fixed check is non-trivially satisfied. Paper round-trip: the F2 fix introduces no new constant (it re-uses `R_star` and `delta_dyn_star` already anchored to the notes), so it cannot create a new paper_misalignment; the F1 manifest reproduces only paper/notes-stated targets.
