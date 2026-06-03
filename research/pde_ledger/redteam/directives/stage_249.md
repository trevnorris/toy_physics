---
unit_id: 249
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-03T12:54:57-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 249

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

There are NO paper_misalignment findings; all three are script-side and may be applied.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py:43-59` (Section 1) and `:72-93` (Section 2)

**Issue:** Section 1's assertion (line 59) compares two expressions that are identically equal by the linearity of `diff` and subtraction over content-free `sp.Function` placeholders, so it cannot fail for any physics — it does not exercise the paper's transfer-law deliverable (notes §1: the projected-minus-resolved *source* equals `-2 E'·B'`, i.e. minus twice the covariance helicity). Section 2's three asserts (lines 91-93) are guaranteed-true distributive/factoring identities (including a factoring after the script self-substitutes `Gamma1 -> alpha_h*Gamma0` at line 79), so the closure structure is asserted but cannot break.

**Required change:**

Section 1 — make the source identity load-bearing. After the existing placeholder block (keep it as documentation if you wish), add distinct source symbols and a non-vacuous identity:
- Introduce `S_full, S_res = sp.symbols("S_full S_res", real=True)` representing the full-field source `E·B` and the resolved-field source `Ē·B̄`.
- Define the covariance source `S_cov = S_full - S_res` (this is `E'·B'`).
- Form the two ledger RHS values: `rhs_proj = -2*S_full`, `rhs_res = -2*S_res`.
- Form the subtracted RHS: `rhs_sub = sp.simplify(rhs_proj - rhs_res)`.
- Assert the transfer-law source identity: `assert sp.simplify(rhs_sub - (-2*S_cov)) == 0`.
- Print `S_cov`, `rhs_sub`, and the matched `-2*S_cov` so the transcript shows the covariance-anchored source.
This is non-trivial in the sense that it ties the subtracted source to the *covariance* `S_full - S_res = S_cov`, which is exactly the physical content of eq:hsub-eq; it reduces to 0 only because `S_cov` is defined as `S_full - S_res`.

Section 2 — add the substantive sign/magnitude checks that distinguish the closure from bare algebra (notes §3.1). Keep the existing symbolic factoring asserts (they are harmless documentation) and ADD:
- Positivity/preference: assert that under `alpha_h > 0`, the aligned branch exceeds the anti-aligned branch. Concretely, with `Gamma0` positive and `alpha_h` symbolic, form `Gplus2 = Gamma0*(1+alpha_h)` and `Gminus2 = Gamma0*(1-alpha_h)` and assert `sp.simplify((Gplus2 - Gminus2) - 2*Gamma0*alpha_h) == 0` (the gap is `2 Gamma0 alpha_h`, so positive iff `alpha_h>0`); also assert via a concrete substitution that `(Gplus2 - Gminus2).subs({Gamma0: 1, alpha_h: sp.Rational(1,2)}) == 1 > 0` and `.subs({Gamma0:1, alpha_h: -sp.Rational(1,2)}) == -1 < 0`, demonstrating the sign flip tracks `alpha_h`.
- Both-positive condition: assert that at `alpha_h = sp.Rational(1,2)` (i.e. `|alpha_h|<1`) both `Gplus2.subs(...)>0` and `Gminus2.subs(...)>0` with `Gamma0->1`, and that at `alpha_h = sp.Rational(3,2)` (`|alpha_h|>1`) the anti-aligned branch goes negative (`Gminus2.subs({Gamma0:1, alpha_h: sp.Rational(3,2)}) < 0`).

**Self-test (already done by auditor):** `rhs_sub - (-2*S_cov)` with `S_cov = S_full - S_res` simplifies to `(-2 S_full + 2 S_res) - (-2(S_full - S_res)) = 0` — non-vacuous, ties to the covariance. The §2 substitutions give literal `1`, `-1`, and the both-positive case gives `1.5>0, 0.5>0`; `alpha_h=3/2` case gives `Gminus2 = 1*(1-3/2) = -1/2 < 0`. No `diff` w.r.t. an absent variable anywhere.

**Verification command:** `redteam exec-sympy 249`; the verifier confirms §1 contains an assertion matching `rhs_sub` to `-2*S_cov` with `S_cov` defined as `S_full - S_res`, §2 contains the `2*Gamma0*alpha_h` gap assertion plus the sign-flip / both-positive substitutions, and the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py`
- summary: Added the covariance-source RHS identity and aligned/anti-aligned sign and positivity checks.
- deviation: none

## F2 — missing_verification_script (subtype: missing_mathematica)

**Target:** NEW file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_mathematica_audit.wl`

**Issue:** Stage 249 is SymPy-only. Mathematica can independently verify every claim (symbolic algebra + numeric benchmark), so the dual-engine rule requires a second, independently-derived engine. The new `.wl` must NOT transliterate the `.py`; use native Mathematica primitives and a different decomposition.

**Required change:** Write the `.wl` using the project's standard harness helpers (`expectZero[expr_]`, `expectTrue[cond_]`, `expectApprox[a_, b_, tol_]` that `Print` a label and `Exit[1]` on failure — mirror the convention used by sibling `*_mathematica_audit.wl` scripts in `mathematica/`). Independently verify the claim manifest below, using `Solve`/`Reduce`/`FullSimplify`/`FreeQ` rather than restating closed forms.

**Claim manifest** (each must be an independent `expect*` check, NOT a port of the .py choreography):

- **M1 — subscale transfer-law source (covariance anchor).** Let `Sfull`, `Sres` be the full and resolved sources. Define `Scov = Sfull - Sres`. From the two ledgers with RHS `-2 Sfull` and `-2 Sres`, derive the subtracted RHS `-2 Sfull - (-2 Sres)` and assert via `expectZero[FullSimplify[(-2 Sfull + 2 Sres) - (-2 Scov)]]` that the projected-minus-resolved source equals `-2 Scov` (i.e. `-2 E'·B'`).
   Symbolic claim: `(-2 S_full) - (-2 S_res) = -2(S_full - S_res) = -2 S_cov`.

- **M2 — closure reduction.** With orientation label `sig` and `Gamma0, Gamma1`, define `Hdot[sig_] := Gamma0 + sig Gamma1` and asymmetry `ah = Gamma1/Gamma0`. Assert `expectZero[FullSimplify[Hdot[sig] - Gamma0 (1 + sig ah)]]`.
   Symbolic claim: `Gamma_0 + sigma Gamma_1 = Gamma_0(1 + sigma alpha_h)` with `alpha_h = Gamma_1/Gamma_0`.

- **M3 — peak Möbius inverse (derive, don't restate).** Set `Rpk = Hdot[1]/Hdot[-1]` expressed in `ah` (i.e. `(1+ah)/(1-ah)`). Use `Reduce`/`Solve` to solve `Rpksym == (1+ah)/(1-ah)` for `ah` and assert `expectZero[FullSimplify[ahSolved - (Rpksym - 1)/(Rpksym + 1)]]`. Solve for `ah`, do NOT hardcode the inverse.
   Symbolic claim: `alpha_h = (R_pk - 1)/(R_pk + 1)`.

- **M4 — integrated ratio + scale cancellation.** With `etah` (overall scale), `I0`, `I1`, form `Hplus = etah (I0 + I1)`, `Hminus = etah (I0 - I1)`, `Rint = FullSimplify[Hplus/Hminus]`. Assert `expectTrue[FreeQ[Rint, etah]]` (scale cancels — eta_h genuinely present then removed, NOT a diff-by-absent-variable) and `expectZero[FullSimplify[Rint - (I0 + I1)/(I0 - I1)]]`. Then `Solve` `Rintsym == (1+abar)/(1-abar)` for `abar` and assert it equals `(Rintsym - 1)/(Rintsym + 1)`.
   Symbolic claim: `R_int = (1 + I1/I0)/(1 - I1/I0)`, independent of `eta_h`; `abar = (R_int-1)/(R_int+1)`.

- **M5 — Session-II benchmark packet.** Hardcode the SAME reported run outputs as the paper card eq:benchmark-packet / notes §5 (these are carried-forward Session-II data, not in-stage derivations): `peakAligned = 281.79830789`, `peakAnti = 56.96878122`, `hAligned = 20.58070146`, `hAnti = 5.00843357`, `RintReport = 4.10920923`, `XiTurn = 0.34437471`, `lambdaTh = 0.42826825`, `vCross = 2.59221845`. Recompute and assert:
   - `expectApprox[peakAligned/peakAnti, 4.94653917, 1*^-7]`  (R_pk)
   - `expectApprox[(peakAligned/peakAnti - 1)/(peakAligned/peakAnti + 1), 0.66366992, 1*^-7]`  (alpha_pk)
   - `expectApprox[hAligned/hAnti, RintReport, 1*^-7]`  (final-load ratio matches reported integrated ratio — the independent cross-check)
   - `expectApprox[(RintReport - 1)/(RintReport + 1), 0.60854999, 1*^-7]`  (abar_h)
   - `expectTrue[0 < (RintReport-1)/(RintReport+1) < (peakAligned/peakAnti - 1)/(peakAligned/peakAnti + 1) < 1]`  (the ordering `0 < abar_h < alpha_pk < 1` from notes §5.3)
   - `expectTrue[XiTurn > 0 && lambdaTh > 0 && vCross > 0]`
   Numeric claims must match the published packet `(0.34437471, 0.42826825, 4.94653917, 4.10920923, 0.66366992, 0.60854999)` to displayed precision.

**Self-test (already done by auditor):** M1 reduces to `0` by definition of `Scov`; M3/M4 use `Solve` to derive inverses (real, can fail if forms wrong); M4's `FreeQ[Rint, etah]` is a true cancellation (etah is present in `Hplus`/`Hminus` then divides out), NOT the absent-variable trap; M5 ratios: `281.79830789/56.96878122 = 4.9465...`, `20.58070146/5.00843357 = 4.1092...`, `(4.1092-1)/(4.1092+1) = 0.60855`, ordering `0 < 0.60855 < 0.66367 < 1` holds. No integral over an unbounded domain, so parity is moot.

**Path note:** the `.wl` MUST live in `mathematica/` with the exact stem `..._mathematica_audit.wl` named above.

**Verification command:** `redteam exec-mathematica 249`; the verifier confirms the `.wl` exists at the named path, uses `expect*` helpers (no bare prints as the only check), independently reproduces the benchmark packet, and exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering the five claim-manifest checks with failing expect helpers.
- deviation: none

## F3 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py:152,187,188`

**Issue:** `alpha_int_num` (line 152) is computed from the hardcoded literal `ratio_integrated_report = 4.10920923`, and line 187 checks it against `0.6085499908172678`, which is just `(4.10920923-1)/(4.10920923+1)` — a function of the same literal, so the check is near self-confirming. The independent anchor is the final-load ratio `ratio_final = hint_aligned/hint_antialigned` (line 151).

**Required change:** Re-anchor the integrated asymmetry to an independent measurement. Change line 152 from
`alpha_int_num = (ratio_integrated_report - 1.0) / (ratio_integrated_report + 1.0)`
to derive it from the final-load ratio:
`alpha_int_num = (ratio_final - 1.0) / (ratio_final + 1.0)`
Keep line 186 (`assert abs(ratio_final - ratio_integrated_report) < 5e-9`) as the consistency tie between final-load and reported integrated ratio. The line 187 tolerance check against `0.6085499908172678` will then exercise the final-load-derived value against the published `abar_h`, making it an independent confirmation rather than a self-restatement. Do NOT change any hardcoded paper value. Confirm `alpha_final_num` (line 153) and the `alpha_peak_num > alpha_int_num` ordering (line 189) still hold.

**Self-test (already done by auditor):** `ratio_final = 20.58070146/5.00843357 = 4.109209231...`, so `alpha_int_num = (4.10920923... - 1)/(4.10920923... + 1) = 0.608549990913...`, which is within `5e-13` of `0.6085499908172678`? The final-load ratio differs from the reported `4.10920923` at the 9th decimal, so the derived alpha differs at ~1e-10. The existing line 187 tolerance is `5e-13` and would FAIL if re-anchored to `ratio_final`. Therefore: when re-anchoring, RELAX the line-187 tolerance to `5e-9` (consistent with line 186's `5e-9` on the ratios) so the now-independent check passes at the precision the reported data supports. Apply both edits together (line 152 re-anchor AND line 187 tolerance `5e-13 -> 5e-9`).

**Verification command:** `redteam exec-sympy 249`; the verifier confirms line 152 derives `alpha_int_num` from `ratio_final`, line 187 tolerance is `5e-9`, and the script exits 0 with the ordering assertions intact.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_sympy_audit.py`
- summary: Re-anchored the integrated asymmetry to the final-load ratio and relaxed the tolerance to the reported-data precision.
- deviation: none
