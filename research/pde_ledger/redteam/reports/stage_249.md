---
unit_id: 249
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 249 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_249.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row 96 + items at lines 278-285)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 249 is a deliberately diagnostic stage layered on the Stage-248 event chain. Its `\stagefield{Output}`-equivalent content (the boxed packet at notes §8 and card eq:benchmark-packet) is the six-number Session-II benchmark `(Xi_turn, lambda_th, R_pk, R_int, alpha_pk, abar_h) = (0.34437471, 0.42826825, 4.94653917, 4.10920923, 0.66366992, 0.60854999)`. The distinct deliverables the stage proves are: (1) the exact subscale-helicity transfer law `partial_t h_sub + div F_{h,sub} = -2 E'·B'` obtained by subtracting the resolved ledger from the projected ledger (eq:hsub-eq); (2) the linear aligned/anti-aligned orientation closure reducing to one asymmetry scalar, `Hdot_{sub,sigma} = Gamma_0(1 + sigma alpha_h)` with `alpha_h = Gamma_1/Gamma_0` (eq:export-rate); (3) the exact peak-ratio and integrated-ratio Möbius compilers `R = (1+alpha)/(1-alpha)` and their inverses `alpha = (R-1)/(R+1)` (eq:instant-ratio, eq:integrated-ratio); (4) cancellation of any common export scale (eta_h) from the preference ratios; (5) the Session-II benchmark numbers, including the strict ordering `0 < abar_h ≈ 0.60855 < alpha_pk ≈ 0.66367 < 1`. The card and appendix both stress this is a hidden-channel unloading preference inside a declared closure, NOT a spin/Pauli theorem.

## What the script claims to verify

The docstring enumerates exactly these five deliverables. Section 1 builds abstract `sp.Function` placeholders and asserts that the projected-minus-resolved combination equals the directly-formed subscale combination. Section 2 builds the linear closure symbolically and asserts three regroupings/factorings. Section 3 forms `Gplus/Gminus`, asserts the Möbius form, and inverts it via `solve`. Section 4 forms the integrated ratio, asserts the `(1+I1/I0)/(1-I1/I0)` form, checks `eta_h` is absent from the simplified ratio, and inverts via `solve`. Section 5 hardcodes the eight reported Session-II run outputs, recomputes `R_pk`, `R_int`, `alpha_pk`, `alpha_h`, and asserts they match the published packet to tight tolerance plus the positivity/ordering inequalities.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) subscale transfer law `∂_t h_sub + ∇·F = -2 E'·B'` | §1 lines 43-59 | **mismatch/tautological** — only verifies that subtraction distributes over abstract functions; the physical transfer-law content (`-2(E·B − Ē·B̄) = -2 E'·B'`, the source-covariance relation) is never represented or exercised |
| (2) closure `Hdot = Gamma_0(1+sigma alpha_h)` | §2 lines 72-93 | partial — checks are pure distributive/factoring identities (tautological by construction); the closure structure is asserted but cannot fail |
| (3a) peak Möbius `R_pk=(1+a)/(1-a)` | §3 line 110 | match (trivial — Gamma_0 cancellation) |
| (3b) peak inverse `a=(R-1)/(R+1)` | §3 line 111 | match — genuine inversion via `solve`, can fail |
| (3c) integrated Möbius | §4 line 132 | match (trivial regrouping) |
| (3d) integrated inverse | §4 line 134 | match — genuine inversion, can fail |
| (4) scale cancellation `eta_h` drops out | §4 line 133 | match — genuine `FreeQ`-style check on a real ratio (NOT the absent-variable trap; eta_h is genuinely present in the constructed ratio at 120-121) |
| (5) Session-II benchmark packet + ordering | §5 lines 184-189 | match — values agree exactly with card eq:benchmark-packet and notes §8; checks can fail on mis-transcription |

`paper_alignment: aligned` — every numeric constant in the script matches the published card/notes to displayed precision; no value disagreement. The defects are internal verification-quality issues (tautology, missing engine), not paper-vs-script disagreements.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(eq_sub - eq_sub_expected)==0` | claim 1 | no (tautological — linearity of `diff`/subtraction over abstract functions) |
| A2 | sympy | 91 | `simplify(Hdot_sigma - Hdot_expected)==0` | claim 2 | no (distribution identity) |
| A3 | sympy | 92 | `simplify(Hdot_factored - Gamma0*(1+sigma*alpha_h))==0` | claim 2 | no (factoring of a self-substituted expr) |
| A4 | sympy | 93 | `simplify(Hdot_scale - eta_h*G0*(1+sigma*alpha_h))==0` | claim 2 | no (factoring after self-substitution) |
| A5 | sympy | 110 | `simplify(Rpk_formula - (1+a)/(1-a))==0` | claim 3a | partial (Gamma_0 cancellation, trivial) |
| A6 | sympy | 111 | `simplify(alpha_from_Rpk - (Rpk-1)/(Rpk+1))==0` | claim 3b | yes (real inversion) |
| A7 | sympy | 132 | `simplify(Rint_formula - Rint_expected)==0` | claim 3c | partial (regrouping) |
| A8 | sympy | 133 | `eta_h not in simplify(Rint_formula).free_symbols` | claim 4 | yes (real cancellation test) |
| A9 | sympy | 134 | `simplify(abar_from_Rint - (Rint-1)/(Rint+1))==0` | claim 3d | yes (real inversion) |
| A10-A24 | sympy | 178-192 | numeric positivity / ratio / ordering asserts | claim 5 | yes (real benchmark consistency, can fail on mis-transcription) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py:43-59`
- (same file) `:72-93`

**What's wrong:**
Section 1's load-bearing assertion (line 59) is `assert simplify(eq_sub - eq_sub_expected) == 0`, where `eq_sub = eq_proj - eq_res` (line 50) and `eq_sub_expected = diff(hsub,t) + Fsub + 2*Csub` with `hsub = hbar - hres`, `Fsub = Fbar - Fres`, `Csub = Cfull - Cres` (lines 46-51). Because `hbar, hres, Fbar, Fres, Cfull, Cres` are content-free `sp.Function` placeholders and `diff(hbar - hres, t) = diff(hbar,t) - diff(hres,t)` is just the linearity of differentiation, `eq_sub` and `eq_sub_expected` are *identically the same expression* — the output transcript (lines 16-17) shows them printed character-for-character identical. The assertion is algebraically guaranteed by construction and cannot fail for any physics. The paper's actual deliverable (eq:hsub-eq / notes §1) is the physical content `-2(Ē·B̄ over the full field − Ē·B̄ of projected fields) = -2 E'·B'`, i.e., that the projected-minus-resolved *source* equals minus twice the covariance helicity. None of that is represented; the script only re-checks that subtraction distributes.

Section 2's three assertions (lines 91-93) are the same defect class: `Hdot_sigma - Hdot_expected` (line 91) is distribution of `-(Phi0+sigma Phi1) - 2(C0+sigma C1)` regrouped by powers of sigma; `Hdot_factored - Gamma0*(1+sigma*alpha_h)` (line 92) is the factoring of `Gamma0 + sigma*(alpha_h*Gamma0)` after the script itself substituted `Gamma1 -> alpha_h*Gamma0` (line 79); line 93 is the same factoring after `Gamma0 -> eta_h*G0`. All three are guaranteed-true distributive/factoring identities with no physics that could break them.

**Why this matters:**
The script's banner claims it verified "the exact projected-minus-resolved subscale-helicity transfer law" and "the exact linear aligned/anti-aligned export closure." A reader trusts those as exercised. In fact the strongest physical claim of the stage (the transfer-law identity, deliverable 1) has zero non-tautological coverage, and the closure structure (deliverable 2) is asserted only as algebra that cannot fail. If the paper's source-covariance bookkeeping were wrong, these checks would still pass.

**Required change:**
Replace the abstract-placeholder identity in §1 with a check that actually carries the source content of the transfer law. Build `S_full` (full-field source `E·B`) and `S_res` (resolved-field source `Ē·B̄`) as distinct symbols, define `S_cov = S_full - S_res` as the covariance source `E'·B'`, and assert that the projected-minus-resolved equation's right-hand side equals `-2*S_cov` — i.e., that subtracting `partial_t h_res + div F_res = -2 S_res` from `partial_t hbar + div F_bar = -2 S_full` yields `partial_t h_sub + div F_sub = -2 (S_full - S_res) = -2 S_cov`. The non-trivial step is that the source on the RHS is the *covariance* `S_full - S_res`, matched to `-2 S_cov`; assert `simplify(rhs_sub - (-2*S_cov)) == 0` where `rhs_sub` is derived from the two ledgers and `S_cov` is the independently-named covariance symbol. For §2, keep the symbolic closure but additionally assert the inequality content that distinguishes the closure from trivial algebra (the positivity/preference conditions of notes §3.1: `alpha_h>0 => Gplus>Gminus` and `|alpha_h|<1 => both branches positive`) as the substantive checks; the bare factoring asserts may remain as documentation but are not the verification.

**Verification:**
After the fix, §1 must contain an assertion whose RHS is `-2*(S_full - S_res)` matched against `-2*S_cov` with `S_cov` defined as `S_full - S_res` (a non-vacuous identity that ties the transfer law's source to the covariance), and §2 must contain at least one sign/magnitude assertion on `alpha_h` tied to the branch ordering. The script still exits 0.

### F2 — missing_verification_script (subtype: missing_mathematica)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py` (sole engine)
- no `.wl` exists under `mathematica/`

**What's wrong:**
Stage 249 is SymPy-only. The manifest marks `is_status_only_candidate: False` and the stage is `checkpoint: False` but the dual-engine rule applies on the "is it possible" test. Every claim here is elementary symbolic algebra (linearity of subtraction, Möbius factoring, `Solve`-based inversion, `FreeQ`-style scale cancellation) plus numeric benchmark arithmetic. Mathematica can independently verify all of it with native primitives (`Solve`, `FullSimplify`, `FreeQ`, machine arithmetic). There is no obstruction that makes this genuinely single-engine-impossible.

**Why this matters:**
The dual-engine policy requires a second, independently-derived engine wherever Mathematica *can* verify the stage. The SymPy script additionally has tautological coverage on its two strongest deliverables (F1), so the absence of an independent engine leaves the physical content (transfer law + closure) effectively unverified by any tool.

**Required change:**
Add a new independent Mathematica audit script at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_mathematica_audit.wl`. It must use a DIFFERENT decomposition than the .py (not a transliteration): derive the Möbius inverse from `Reduce`/`Solve` on the orientation closure rather than re-stating closed forms, verify scale cancellation by symbolic limit/`FreeQ` of the constructed ratio, derive the transfer-law source as the covariance directly, and recompute the benchmark packet with `expectApprox`. Detailed claim manifest is in the directive.

**Verification:**
A `.wl` appears at the named path, uses `expectZero`/`expectTrue`/`expectApprox` (not bare prints), independently reproduces `R_pk`, `R_int`, `alpha_pk`, `abar_h` to the published values, and exits 0.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py:184-188`

**What's wrong:**
The benchmark block hardcodes the eight reported Session-II run outputs (lines 140-147) and checks recomputed ratios against the published packet — this is legitimate carry-forward (notes §5 states these are the reported run outputs, not in-stage derivations). However, the integrated-ratio path has a redundancy that weakens it: `ratio_integrated_report = 4.10920923` is a hardcoded literal (line 144), and `alpha_int_num` (line 152) is computed FROM that literal, then line 187 checks `alpha_int_num` against `0.6085499908172678` — a number that is just `(4.10920923-1)/(4.10920923+1)`, i.e., a function of the same literal. That sub-check is close to self-confirming. The genuinely independent cross-check is line 186 (`ratio_final ≈ ratio_integrated_report`), which ties the final-load ratio `20.58070146/5.00843357` to the independently-reported integrated ratio — that one is substantive and should be the anchor.

**Why this matters:**
Two of the §5 asserts (lines 187 against the derived constant, and the `alpha_int` half of the ordering) trace back to a single hardcoded literal, so they do not add independent confidence beyond line 186. Not a correctness error, but the verification is thinner than the four printed "verified objects" imply.

**Required change:**
Keep line 186 (final-load vs reported integrated ratio) as the load-bearing benchmark check. The redundant `alpha_int_num` check (line 187) may stay as documentation but should be recognized as non-independent; preferably re-derive `alpha_int` from the *final-load* ratio (`ratio_final`) rather than from the reported-ratio literal so the asymmetry scalar is anchored to an independent measurement. No paper value changes.

**Verification:**
`alpha_int_num` (or an added `alpha_final` ordering check) is computed from `ratio_final`, and the `alpha_pk > alpha_int` ordering (line 189) holds against that independently-anchored value. Script exits 0.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet apply. The directive's new-script manifest mandates an independent decomposition (Reduce/Solve-based inversion, symbolic scale-cancellation, covariance-anchored transfer law) explicitly NOT mirroring the .py choreography.

## Engine cross-check

Only SymPy is present; `engines_agree: n/a`.

## Verdict justification

`verdict: findings`, `stop_cold: null`. Paper alignment is exact — every numeric constant in the script (the six-number packet, `v_cross`, the eight reported run outputs) matches card eq:benchmark-packet and notes §8/§5 to displayed precision, so there is no `paper_misalignment`. The benchmark consistency checks (§5) and the two Möbius inversions and the eta_h cancellation (A6, A8, A9) hold up under attack and are genuinely falsifiable. What does NOT hold up: the two strongest physical deliverables (the subscale transfer law in §1 and the closure structure in §2) are verified only by tautological distributive/factoring identities over content-free placeholders (F1), and the stage has no second engine despite Mathematica being fully capable of verifying it (F2). The benchmark has a minor self-referential redundancy in the integrated-asymmetry sub-check (F3). None of the corrections changes a forwarded value, so no `CRITICAL_DOWNSTREAM`.

## Self-test notes

Variable-independence trap: checked — the one independence-style assertion (line 133, `eta_h not in free_symbols`) is NOT the absent-variable trap; `eta_h` is genuinely present in the constructed ratio (lines 120-121) and the check is that it cancels from a real division, which is legitimate. My prescribed F1 fix uses `S_cov = S_full - S_res` as a present, non-vacuous symbol (no `diff` w.r.t. an absent variable). Symmetry/parity: no integrals over unbounded domains in this stage, so parity is moot. Trivial-case pre-check: the prescribed transfer-law identity `rhs_sub = -2(S_full - S_res)` matched to `-2*S_cov` reduces to `0` only because `S_cov := S_full - S_res`, which is the intended non-trivial source identity, not a guaranteed-by-construction equality of the F1 type. Paper round-trip: the F1/F3 fixes introduce no new constants and touch no published value, so they cannot create a new paper_misalignment.
