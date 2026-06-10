---
unit_id: 249
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 249 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_249.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows at lines 96, 278-285)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_mathematica_audit.txt`

## What the paper claims

Stage 249 is an explicitly diagnostic (non-constitutive) stage attached to the Stage-248 event chain. It proves: (1) the exact subscale-helicity transfer law `\partial_t h_sub + \nabla_3·F_{h,sub} = -2 \overline{E'·B'}` obtained by subtracting the projected and resolved helicity ledgers (eq:stage249-hsub-eq); (2) the linear aligned/anti-aligned orientation closure `\dot H_{sub,sigma} = Gamma_0 + sigma Gamma_1 = Gamma_0(1+sigma alpha_h)` (eq:export-rate); (3) the exactly-invertible peak and integrated Möbius ratio compilers `R = (1+a)/(1-a)`, `a = (R-1)/(R+1)` (eq:instant-ratio, eq:integrated-ratio); (4) cancellation of any common export scale `eta_h` from the ratios; (5) the Session-II benchmark packet. The card quotes the benchmark verbatim: `(\Xi_{turn},\lambda_{th},R_pk) = (0.34437471,0.42826825,4.94653917)` and `(R_int,\alpha_pk,\bar\alpha_h)=(4.10920923,0.66366992,0.60854999)`. The `\stagefield{Verification}` line states "Mathematica audit: none yet."

## What the script claims to verify

The SymPy script's docstring claims it verifies all five deliverables. Section 1 asserts the difference-of-ledgers transfer equation and the RHS source identity `-2 S_cov`; Section 2 asserts the closure expansion, the `Gamma_0(1+sigma alpha_h)` factorization, the `eta_h` scaling form, the branch gap `2 Gamma_0 alpha_h`, and concrete-value positivity/sign checks; Section 3 derives `R_pk` from the branch definitions and *solves* for `alpha_h`, asserting both the Möbius form and the inverse; Section 4 does the integrated analogue and asserts `eta_h not in free_symbols`; Section 5 reconstructs the benchmark numbers and asserts them against the card values to tight tolerances plus the ordering `0 < \bar\alpha < \alpha_pk < 1`. The `.wl` mirrors deliverables 1-5 with leaner, differently-structured checks (M1-M5).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) subscale transfer law / `-2 S_cov` source | py S1 (eq_sub linearity + rhs_sub); wl M1 | match |
| (2) linear closure → `Gamma_0(1+sigma alpha_h)` | py S2; wl M2 | match |
| (3) peak + integrated Möbius + inverse | py S3/S4 (solve-based inverse); wl M3/M4 | match |
| (4) `eta_h` cancellation | py S4 (`free_symbols`); wl M4 (`FreeQ`) | match |
| (5) Session-II benchmark packet | py S5; wl M5 | match |
| Verification line "Mathematica audit: none yet" | a `.wl` now exists + passes | mismatch (stale paper-card claim) |

`paper_alignment: partial` — all five math deliverables align; the only discrepancy is the stale "none yet" Verification line.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(eq_sub - eq_sub_expected) == 0` | claim 1 (linearity of ledger subtraction) | yes |
| A2 | sympy | 71 | `simplify(rhs_sub - (-2*S_cov)) == 0` | claim 1 (source) | yes |
| A3 | sympy | 103-105 | closure expand / factor / eta_h scale | claim 2 | yes |
| A4 | sympy | 115-120 | branch gap + concrete sign/positivity | claim 2 / 3.1 | yes |
| A5 | sympy | 137-138 | `R_pk` form + `alpha_from_Rpk` via solve | claim 3 | yes |
| A6 | sympy | 159-161 | `R_int` form + `eta_h` absent + inverse | claim 3,4 | yes |
| A7 | sympy | 211-216 | benchmark numbers + ordering | claim 5 | yes |
| A8 | math | 82 | `expectZero[subtractedRHS-(-2 Scov)]` | claim 1 | yes |
| A9 | math | 88 | `expectZero[Hdot[sig]-Gamma0(1+sig ah)]` | claim 2 | yes |
| A10 | math | 100 | `expectZero[ahSolved-(Rpk-1)/(Rpk+1)]` | claim 3 | yes |
| A11 | math | 111-112 | `FreeQ[Rint,etah]` + ratio form | claim 3,4 | yes |
| A12 | math | 144-152 | benchmark approx + ordering | claim 5 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (direction: paper-side update of a stale Verification line)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_249.tex:4` quote: "SymPy audit: \StageFile{...sympy_audit.py}.  Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage249_..._mathematica_audit.wl` (full file; passing output in `mathematica/output/..._mathematica_audit.txt`)

**What's wrong:**
The stage card's `\stagefield{Verification}` line asserts there is no Mathematica audit ("none yet"), but pass-1 added an independent `.wl` second engine that exists and passes (M1-M5 all PASS). The card text is stale relative to the script tree.

**Why this matters:**
The paper under-reports verification coverage; a reader/auditor checking the card would believe stage 249 is single-engine when it is now dual-engine. Pure documentation lag, no math impact.

**Required change:**
Paper-side only — routed to user. See `## Resolve before fix_loop`. Codex must not edit the card autonomously.

**Verification:**
After user resolution, card line 4 names the `.wl` file instead of "none yet"; no script change.

## Independent-derivation check (Mathematica)

The `.wl` is a genuinely independent re-derivation, not a transliteration. Three contrasts: (a) The `.py` Section 1 builds the full symbolic transfer equation using `sp.Function` operators (`hbar(t)`, `Fbar(t)`, `Derivative`) and asserts the *linearity* of the ledger subtraction (`eq_sub - eq_sub_expected`); the `.wl` M1 does NOT replicate that derivative/flux algebra at all — it verifies only the scalar RHS source identity `subtractedRHS - (-2 Scov)`. (b) The `.py` defines `Phi_sigma`/`C_sigma` and expands `Hdot_sigma`, then separately substitutes `Gamma1 -> alpha_h*Gamma0`; the `.wl` M2 instead defines `Hdot[s_] := Gamma0 + s Gamma1` directly and sets `ah = Gamma1/Gamma0`, a different choreography reaching the same `Gamma0(1+sig ah)` factor. (c) The `.wl` fuses M2/M3 and reaches the peak inverse via `Solve[Rpksym == RpkExpr, ah]` with a fresh `RpkExpr = Hdot[1]/Hdot[-1]`, structurally distinct from the `.py` `solve` route. The benchmark section (M5) is necessarily numeric-data-similar, but that is shared input data, not ported algebra. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree. SymPy: `alpha_pk = 0.663669919237628`, `R_int(final) = 4.109209231260703`, `alpha_int = 0.6085499909138585`, all symbolic identities → 0. Mathematica: `alphaPeakBench = 0.6636699192376279`, `RfinalBench = 4.109209231260703`, `abarBench = 0.6085499908172678`, M1-M4 all `= 0`, benchmark diffs all `~1e-9` < `1e-7` tol. Identical numerics to displayed precision and identical symbolic residuals (0). No `engine_disagreement`.

## Pass-1 tautology/round-trip re-check

Pass-1 flagged a tautological/round-trip finding. The current fix holds:
- Section 3 inverse (`alpha_from_Rpk`) is obtained via `sp.solve(sp.Eq(Rpk, Rpk_formula), alpha_h)[0]` — SymPy independently solves the relation; the assertion then compares the *solved* result to `(Rpk-1)/(Rpk+1)`. This is a real inversion check, not `x==x`. The output confirms `alpha(R_pk) = (Rpk-1)/(Rpk+1)` was derived, not assumed.
- Section 4 likewise solves `abar` via `sp.solve` and additionally tests `eta_h not in free_symbols` (a check that can genuinely fail). The `.wl` mirrors with independent `Solve` + `FreeQ`.
- Section 1's `eq_sub - eq_sub_expected` is a linearity identity over `sp.Function` operators, not a definitional round-trip — it would catch a non-linear ledger error.
No surviving tautology or self-referential round-trip.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| R_pk = 4.94653917 | py:211 / out:51; wl:144 / out:46 | tex:72,553 ; md:443,553 | MATCH |
| alpha_pk = 0.66366992 (0.663669919...) | py:212 / out:52; wl:145 / out:48 | tex:75,555 ; md:452,555 | MATCH |
| R_int = 4.10920923 | py:171,213 / out:56; wl:129,146 / out:50 | tex:74,561 ; md:463,483,561 | MATCH |
| abar_h = 0.60854999 (0.6085499908...) | py:214 / out:57; wl:147 / out:52 | tex:76,563 ; md:473,563 | MATCH |
| Xi_turn = 0.34437471 | py:172 / out:63; wl:130 | tex:73 ; md:406,582 | MATCH |
| lambda_th = 0.42826825 | py:173 / out:64; wl:131 | tex:73 ; md:408,582 | MATCH |
| peak_aligned = 281.79830789 | py:167 / out:49; wl:125 | md:420 | MATCH (notes) |
| peak_antialigned = 56.96878122 | py:168 / out:50; wl:126 | md:423 | MATCH (notes) |
| H_final aligned = 20.58070146 | py:169 / out:53; wl:127 | md:428,571(ratio) | MATCH (notes) |
| H_final anti = 5.00843357 | py:170 / out:54; wl:128 | md:431 | MATCH (notes) |
| v_cross = 2.59221845 | py:174 / out:65; wl:132 | md:414 | MATCH (notes) |

INTERNAL (no prose expected): ratio_peak/ratio_final intermediates, peak/final difference & sum prints, Gamma0/Gamma1 symbolic forms, Phi_sigma/C_sigma, branch_gap, generic symbolic R/abar forms, pass/diff scaffolding.

reconciliation: complete; 11 deliverable values checked, 0 misaligned. (Every emitted number lands in either the card or the notes at matching precision. The only paper-side discrepancy is the prose "Mathematica audit: none yet" line, captured as F1 — not a value mismatch.)

## Verdict justification

The math holds up under attack. Both engines are present and independent (not a transliteration), they agree to displayed precision, every symbolic identity reduces to 0, and the pass-1 round-trip concern is resolved (the Möbius inverses are genuinely `solve`-derived, and `eta_h` cancellation is a falsifiable `free_symbols`/`FreeQ` test). All 11 deliverable values reconcile against the card and notes. The sole finding is a low-severity `paper_misalignment`: the card's `\stagefield{Verification}` line still says "Mathematica audit: none yet" despite the now-present, passing `.wl`. Because that is a paper-side prose discrepancy requiring user direction, `verdict: findings`, `stop_cold: null`, and the directive carries only a `## Resolve before fix_loop` block (no Codex script edits).

## Self-test notes

Variable-independence trap: the only differentiation is `sp.diff(hbar(t),t)` / `sp.diff(hres(t),t)` and `sp.diff(hsub,t)` where `hsub = hbar(t)-hres(t)` genuinely depends on `t` — derivatives are non-vanishing, no zero-derivative trap. Symmetry/parity: no unbounded integrals (integrals are symbolic `I0,I1` constants). Trivial-case: concrete substitutions in S2 (alpha_h=1/2, 3/2) give the asserted literal signs (`+1`, `-1`, positive, negative) — confirmed against output. No directive edits prescribed (only a user-resolution block), so no paper round-trip risk introduced.
