---
unit_id: 067
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md
  paper_appendix: present
---

# Audit unit 067 red-team report (second pass)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_067.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage067_sech_gaussian_resonance.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (the only stage 067 reference is the `\input{stages/stage_067}` line at 252; no prose summary row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage067_sech_gaussian_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.txt`

Prior-pass directive (`redteam/directives/stage_067.md`) records that four findings (F1-F4 from the first pass: duality-implication labeling, tautological self-dual scaffolding, missing sympy norm integration, mathematica numeric-target provenance) were applied at 2026-05-22T19:54:59 and the current scripts/outputs reflect those fixes.

## What the paper claims

The stage card (`paper/stages/stage_067.tex`) declares the unit a sech-Gaussian coherence resonance benchmark, with `\stagefield{Inputs}{Source profile chi_sigma=sech(y/w_f) and support profile chi_phi=e^{-y^2/w_g^2}.}` and three boxed deliverables: (i) the exact self-dual stationary point `w_g/w_f = sqrt(pi)` (eq:app-stage067-self-dual), (ii) the resonance coherence `C_res^2 ~= 0.994418836451529` (eq:app-stage067-Cres), and (iii) the penalty `P_res = 1/C_res^2 ~= 1.005612487760576` (eq:app-stage067-Pres). `\stagefield{Output}` collects (ii) and (iii): "The near-perfect benchmark coherence (eq:app-stage067-Cres) and penalty (eq:app-stage067-Pres)." `\claimstatus{}` marks the family exact-closure and the quoted coherence value numerical. The notes also enumerate (iv) the duality identity `I(r) = (r/sqrt(pi)) I(pi/r)` and the implied symmetry `C^2(r) = C^2(pi/r)`, plus (v) the uniqueness of `r_*` as the global maximum on the constructive branch (claimed there via "a numerical monotonicity audit").

## What the script claims to verify

After the first-pass fixes, both scripts now exercise: (1) the closed-form transverse norms `N_ss = 2 w_f` and `N_pp = w_g sqrt(pi/2)` via direct symbolic integration in both engines; (2) the algebraic implication `I(r) -> (r/sqrt(pi)) I(pi/r) => C^2(r) = C^2(pi/r)`, explicitly labeled as implication-only in inline comments; (3) the self-dual stationarity of `C^2` at `r_* = sqrt(pi)` from the symmetry, with the formal expect-zero assertions labeled as tautological calculus identities and the substantive evidence shifted to the broken-tangent perturbation (sympy lines 137-142) and the numerical monotonicity scan; (4) the numerical 60-dps evaluation `C_res^2 = 0.99441883645152934870...` (sympy mpmath quad; mathematica `NIntegrate` at WorkingPrecision 80); (5) numerical verification of the sech-Gaussian duality identity at five sample r-values; (6) numerical strict monotonicity on the left and right of `r_*`. The mathematica numeric targets are now annotated as cross-engine cross-checks against the sympy run.

## Paper <-> script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| Exact self-dual point `w_g/w_f = sqrt(pi)` (eq:app-stage067-self-dual) | Sec.3 of both: differentiated symmetry forces zero slope at `r_*` (labeled tautological), broken-tangent perturbation (sympy 137-142) confirms non-vacuity, Sec.6 monotonicity localizes the maximum | match |
| Exact duality `C^2(r) = C^2(pi/r)` (notes Sec.2) | Sec.2: algebraic implication asserted (now labeled implication-only); Sec.5: numerical `\|I(r) - (r/sqrt(pi)) I(pi/r)\|` at five r-samples vanishes to <= 1e-40 (sympy) / <= 1e-35 (mathematica) | match |
| `C_res^2 = 0.994418836451529...` (eq:app-stage067-Cres) | Sec.4 sympy: mpmath quad at 60 dps; Sec.4 mathematica: NIntegrate at WP 80, cross-checked against sympy value at 35-digit tolerance | match |
| `P_res = 1.005612487760576...` (eq:app-stage067-Pres) | Sec.4: `1/C_res^2` computed and (mathematica) cross-checked at 34-digit tolerance | match |
| Norms `N_ss = 2 w_f`, `N_pp = w_g sqrt(pi/2)` (notes Sec.1) | Sec.1 of both: direct symbolic integration confirms the boxed values (sympy now uses `(...).rewrite(sp.cosh)` so the integral evaluates) | match |
| Uniqueness of `r_*` as global maximum on constructive branch (notes Sec.3) | Sec.6 of both: strict-increase on `{0.55, ..., r_*}` and strict-decrease on `{r_*, ..., 4}` | match (sampled, not exhaustive — notes label it numerical too) |

`paper_alignment: aligned`. All paper-side deliverables are exercised, and at the precision/coverage the paper itself claims (exact closure for the analytic ratio, numerical for the coherence value and monotonicity).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74 | `expect_zero("N_ss integral - 2 w_f", Nss_integral - Nss)` | norms (notes Sec.1) | yes |
| A2 | sympy | 75 | `expect_zero("N_pp integral - w_g sqrt(pi/2)", ...)` | norms | yes |
| A3 | sympy | 94 | `expect_zero("C^2(r) - C^2(pi/r) under duality", C2_dual - C2_target)` | duality (notes Sec.2) | partial (algebraic implication only, labeled) |
| A4 | sympy | 118-121 | `expect_zero("self-dual overlap-slope relation", ...)` | stationary point | no (tautological, labeled at 116-117) |
| A5 | sympy | 131-134 | `expect_zero("stationary derivative of C^2 at the self-dual point", ...)` | stationary point | no (tautological, labeled at 128-130) |
| A6 | sympy | 141-142 | `if dC2_broken == 0: raise` (broken-tangent must be nonzero) | stationary-point non-vacuity | yes |
| A7 | sympy | 183-184 | `if diff > 1e-40: raise` over five r-samples | sech-Gaussian duality identity | yes |
| A8 | sympy | 198-200 | strict-increase on left grid | monotonicity / uniqueness | yes |
| A9 | sympy | 207-209 | strict-decrease on right grid | monotonicity / uniqueness | yes |
| B1 | mathematica | 52 | `expectZero["N_ss - 2 w_f", nssDirect - nssExpected]` | norms | yes |
| B2 | mathematica | 53 | `expectZero["N_pp - w_g sqrt(pi/2)", nppDirect - nppExpected]` | norms | yes |
| B3 | mathematica | 67 | `expectZero["C^2(r) - C^2(pi/r) under duality", c2Dual - c2Target]` | duality | partial (algebraic implication only, labeled at 64-66) |
| B4 | mathematica | 93-96 | `expectZero["self-dual C^2 stationary slope from symmetry solve", ...]` | stationary point | no (tautological, labeled at 89-92) |
| B5 | mathematica | 134 | `expectApprox["C_res^2 numeric check", c2Star, c2Target, 10^-35]` | C_res^2 benchmark | yes (cross-engine) |
| B6 | mathematica | 135 | `expectApprox["P_res numeric check", pres, presTarget, 10^-34]` | P_res benchmark | yes (cross-engine) |
| B7 | mathematica | 145 | `expectTrue["duality sample ...", diff <= 10^-35]` over five r-samples | duality identity | yes |
| B8 | mathematica | 157-160 | left-grid strict-increase | monotonicity / uniqueness | yes |
| B9 | mathematica | 168-171 | right-grid strict-decrease | monotonicity / uniqueness | yes |

The tautological assertions A4/A5/B4 are kept as formal assertions but explicitly annotated as such in adjacent comments. The first-pass auditor chose to label rather than demote; the second pass does not re-litigate that choice. The substantive stationary-point evidence (A6 and A8/A9 / B7-B9) covers the paper claim.

## Findings

### F1 — paper_misalignment

**Subtype:** notes_contradicts_script (script labels are pre-renumbering "Stage 50/050" while paper card and file paths are "Stage 067")

**Severity:** low
**Files:**
- paper-side: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_067.tex:1` quote: `\section[Stage 067]{Stage 067: Exact Sech--Gaussian Coherence Resonance Benchmark}`
- script-side: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:4` quote: `moving_throat_pde_stage50_sech_gaussian_sympy_audit.py`
- script-side: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py:53` quote: `banner("STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK")`
- script-side: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl:38` quote: `banner["STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"];`

**What's wrong:**
After the renumbering that produced `stage_067.tex` and the matching script filenames (`..._stage067_...py` / `..._stage067_...wl`), the script docstring and the printed banners still carry the old `Stage 50` / `STAGE 050` label. The transcripts emit:
- sympy output line 3: `STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK`
- mathematica output line 3-4: `STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK`
- mathematica output line 78: `Stage 067 Mathematica audit passed.` (this line *was* fixed, so the inconsistency is internal even to the .wl).

There is no actual math error. The notes file's title also still reads `Stage 50` (notes line 2), but the red-team must not edit notes/, so only the script-side labels are in scope. The direction of resolution is unambiguous (paper card and filenames are authoritative; "067" is the correct label everywhere), so no `## Resolve before fix_loop` block is needed.

**Why this matters:**
Audit transcripts feed downstream review. A reader cross-referencing the sympy output banner against the paper card sees "STAGE 50" and "Stage 067" and has to verify they are the same stage. The mathematica `.wl` is even inconsistent with itself: the banner says 050 but the closing print says 067. Fixing the labels keeps the transcripts aligned with the paper card and removes a stumbling block for downstream review.

**Required change:**
Bring the script labels into agreement with the file paths and paper card. Specifically:
1. In `scripts/moving_throat_pde_stage067_sech_gaussian_sympy_audit.py`, replace the docstring filename on line 4 from `moving_throat_pde_stage50_sech_gaussian_sympy_audit.py` to `moving_throat_pde_stage067_sech_gaussian_sympy_audit.py`, and replace the banner on line 53 from `"STAGE 50 — EXACT SECH–GAUSSIAN COHERENCE BENCHMARK"` to `"STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"`.
2. In `mathematica/moving_throat_pde_stage067_sech_gaussian_resonance_mathematica_audit.wl`, replace the banner on line 38 from `"STAGE 050 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"` to `"STAGE 067 — EXACT SECH-GAUSSIAN COHERENCE BENCHMARK"`.

Do not touch the notes file or the paper card. The notes' "Stage 50" header is documentation drift that the user can address out-of-band; it does not affect script correctness.

**Verification:**
After the edit, re-running the scripts emits banner lines that begin with `STAGE 067` (matching the paper card and filenames). The sympy docstring's filename line matches the actual filename. No assertion changes; the `PASS` / numeric content is unchanged.

## Independent-derivation check (Mathematica)

The Mathematica script is independent from the SymPy script on the math side (already concluded in the first-pass report, still holds):
- Norms: Mathematica integrates `Sech[y/wf]^2` and `Exp[-2 y^2/wg^2]` directly (lines 45-46); SymPy rewrites `sech` as `cosh` then integrates (lines 70-71). Same identity, different code paths.
- Stationary-point block: Mathematica works in `C^2`-space using `D[c2Fn[r] - c2Fn[Pi/r], r]` with `Solve` (lines 75-87); SymPy works in `I`-space with a hand-built `duality_tangent` and a separate hand-written `dC2_selfdual` formula (sympy lines 109-134). Different algebraic backbones reaching the same calculus identity.
- Numeric integrator: Mathematica uses `NIntegrate` with `WorkingPrecision -> 80, AccuracyGoal -> 32` over the half-line (then `*2`); SymPy uses `mpmath.quad` over `(-inf, inf)` at `dps=60`. Independent.

The mathematica numeric targets are explicitly labeled (lines 122-125) as cross-engine targets sourced from the sympy transcript, which is honest documentation of an agreement check.

No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree to >35 digits on `C_res^2` and `P_res`:
- sympy (output line 36): `C_res^2 = 0.994418836451529348706428351608877628170873348983716948813464`
- mathematica (output line 31): `C_res^2 = 0.99441883645152934870642835160887762817087334898371694998969514187456256874872`60.`
- The first 50+ digits match.
- `expectApprox["C_res^2 numeric check", c2Star, c2Target, 10^-35]` passes (mathematica output line 35).
- Duality samples at `r = 3/4, 1, 6/5, 3/2, 2` vanish to well below the tolerance in both engines.
- Monotonicity grids agree on both engines.

`engines_agree: true`.

Output freshness:
- sympy script mtime 2026-05-22 19:54, sympy output 19:56 (fresh, post-first-pass fix).
- mathematica script 19:53, mathematica output 19:56 (fresh, post-first-pass fix).

`outputs_fresh: true`.

## Verdict justification

Every paper deliverable is verified by substantive script-side evidence: the norms by direct symbolic integration in both engines (post-first-pass fix on the sympy side); the duality identity by 60-dps numerical evaluation at five sample r-values, in agreement to <= 1e-40 (sympy) / <= 1e-35 (mathematica); the stationary point at `r_* = sqrt(pi)` by the broken-tangent counter-check plus a strict-monotonicity scan on either side of `r_*`; and the boxed numerical values `C_res^2` and `P_res` by independent 60-dps integrations in two engines that agree to >35 digits. Tautological calculus identities are present in both scripts but explicitly labeled as such in the source. After the first-pass fixes (F1-F4 from the prior directive), the only remaining issue is cosmetic: the script docstring and printed banners still carry the pre-renumbering "Stage 50" / "STAGE 050" label rather than "Stage 067". This is a paper-script misalignment in the labels but not in the math; resolving it is mechanical (the paper card and filename use "067"), so it is being routed directly to Codex rather than to the user. Verdict: `findings`, one low-severity label/banner cleanup.

Attacks tried that failed in this pass: (i) re-checked whether the new `(... ).rewrite(sp.cosh)` step in sympy line 70 changes the integrand's domain — no, `sech(y/wf)^2 = 1/cosh(y/wf)^2`, both even and positive, integral evaluates to `2 w_f` for `wf > 0`; (ii) re-checked the half-line `2 * NIntegrate[... {x, 0, Infinity}]` is correct for the even integrand — yes, `sech(x) Exp[-x^2/r^2]` is even; (iii) checked that the labeled-tautology comments do not silently mask a real check the script would otherwise need — no, the `dC2_broken` block (sympy 137-142) handles the non-vacuity, and the monotonicity scan handles uniqueness on the sampled branch; (iv) re-verified `expect_zero` triggers `raise AssertionError` rather than silently passing (sympy line 50) — confirmed, behavior is correct; (v) checked Mathematica's `Solve[c2SymmetryAtRStar == 0, c2PrimeLeft, Reals]` for hidden `ConditionalExpression` wrappers — output transcript line 23 shows `{{C2PrimeLeft -> 0}}` clean, no conditional expression appears.

## Self-test notes

I checked: (1) variable independence — F1 is label-only, no derivatives or assertions are added or modified; (2) symmetry/parity — N/A for label-only fix; (3) trivial-case substitution — N/A (no new mathematical expressions); (4) path specs — both targeted files exist at the absolute paths quoted (`scripts/...stage067...py` and `mathematica/...stage067...wl`), confirmed via the file listing; (5) paper round-trip — after relabeling the scripts to "Stage 067", the script labels match the paper card (`\section{Stage 067: ...}`), the appendix `\input{stages/stage_067}`, and the file paths. No new paper_misalignment is introduced. The notes file's "Stage 50" header is the one remaining pre-renumbering label, but it is out of scope for the red-team and does not affect script behavior.
