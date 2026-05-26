---
unit_id: 036
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 036 (v2)

This verification supersedes the prior v1 verification (which covered the original 6-finding pass). This pass covers the v2 directive's two additional findings.

## Per-finding outcomes

### F1 — tautological_check (M_mix admissible/inadmissible witness)

**Classification:** resolved

**What changed:**
- SymPy (`scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`):
  - Lines 101-103: `Mmix_good = sp.simplify(G_sample - sp.Rational(1, 10))` and `Mmix_bad = sp.simplify(G_sample + sp.Rational(1, 10))` replaced with
    `Mmix_admissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 29}))`,
    `Mmix_inadmissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 1}))`,
    `G_sample_n = sp.N(G_sample)`.
  - Lines 146-155: the two `expect_true` blocks now compare `Mmix_admissible < G_sample_n` and `Mmix_inadmissible > G_sample_n`, with the admissible-side label updated to `M_mix < G(xi_req,delta)` (matching the directive's strict-`<`).
- Mathematica (`mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`):
  - Lines 113-115: `mMixGood = FullSimplify[gSample - 1/10]` / `mMixBad = FullSimplify[gSample + 1/10]` replaced with
    `mMixAdmissible = N[mMix /. {Chi -> 1, OmegaU -> 1, Delta0 -> 1, A -> 29}]`,
    `mMixInadmissible = N[mMix /. {Chi -> 1, OmegaU -> 1, Delta0 -> 1, A -> 1}]`,
    `gSampleN = N[gSample]`.
  - Lines 145-154: the two `expectTrue` blocks updated to use the new symbols and the strict-`<` admissible label.

**Assessment:**
The edit is exactly what the directive prescribed. Witnesses are now derived from the paper definition `Mmix_expr = 8 Chi^2/(pi^2 A Omega_U^2 Delta_0)` (and its Mathematica analogue `mMix`) at two independent parameter tuples — `(Chi=1, Omega_U=1, Delta_0=1, A=29)` for admissible and `(Chi=1, Omega_U=1, Delta_0=1, A=1)` for inadmissible — not by shifting `G_sample` by ±1/10. The exec logs show the new numeric values match the auditor's self-test:

- SymPy log line 33: `M_mix=0.0279506713496104, G=0.465517241379310` (admissible — matches `8/(29 pi^2) ≈ 0.02795` and `27/58 ≈ 0.46552`).
- SymPy log line 34: `M_mix=0.810569469138702, G=0.465517241379310` (inadmissible — matches `8/pi^2 ≈ 0.81057`).
- Mathematica log lines 44 and 46: same numeric values to ~17 digits.

The assertions are no longer tautological: a sign flip in the closed form of `G`, or a corruption of `Mmix_expr`'s closed form, would now actually break these checks because the two sides are computed from independent symbolic sources before being numerically compared. The two witnesses straddle `G_sample = 27/58 ≈ 0.4655` non-trivially (one well below at 0.028, one well above at 0.811). No collateral edits beyond the named line ranges.

### F2 — tautological_check (g_B,req factorization labelling)

**Classification:** resolved

**What changed:**
- SymPy: 4-line comment inserted between `print("R_target =", R_target)` and the first `expect_zero(...)` (lines 68-71); 1-line comment inserted between `print("Parametric frontier: ...")` and the second `expect_zero(...)` (line 91). Comment text matches the directive verbatim.
- Mathematica: 4-line comment (within a single `(* ... *)` block) inserted between `expectZero["G - closed form", g - gTarget];` and the next `expectZero[...]` block (lines 59-62); 1-line `(* ... *)` comment inserted between the `$Assumptions = ...` line and the next `expectZero[...]` block (line 103). Comment text matches the directive verbatim.

**Assessment:**
Comments landed at the prescribed insertion points with the prescribed wording. They flag the four definitional self-consistency assertions and point the reader to the symbolic-kappa derivation (`symbolic kappa derivation: F(xi,delta) - R_target_sym`) as the genuine anchor. The directive explicitly accepted the comment-only ("minimal") fix; the stronger derivation rewrite was optional. The `expect_zero`/`expectZero` calls themselves are unchanged, so the printed residuals stay 0 in both exec logs (sympy log line 14: `g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix) = 0`; line 29: `final-test support inequality <-> nonnegative required support loading = 0`; mathematica log lines 18-19 and 36-37: PASS on the same). The underlying tautological character of these particular assertions is unchanged — but the directive made clear this is a labelling fix, not a substance fix, and that the substance is anchored elsewhere via the symbolic kappa derivation, which the logs confirm passes independently (sympy log line 32, mathematica log line 43).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Line 14: `g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix) = 0` (F2 assertion still passes after comment insertion).
- Line 29: `final-test support inequality <-> nonnegative required support loading = 0` (F2 echo assertion still passes).
- Line 32: `symbolic kappa derivation: F(xi,delta) - R_target_sym = 0` (the genuine anchor referenced by F2's comments — independently confirms the F/G algebraic structure).
- Line 33: `admissible sample: M_mix < G(xi_req,delta) = M_mix=0.0279506713496104, G=0.465517241379310` (F1 fix — value computed from `Mmix_expr` at `A=29`, independent of `G_sample`).
- Line 34: `inadmissible sample: support deficit blocks the branch = M_mix=0.810569469138702, G=0.465517241379310` (F1 fix — value from `Mmix_expr` at `A=1`).
- Line 41: `All Stage 19 checks passed.`

**Mathematica:** exit=0. Notable lines:
- Lines 18-19, 36-37: PASS on the definitional residuals (F2 — unchanged by the comment insertion).
- Line 43: `PASS: symbolic kappa derivation: F(xi,delta) - R_target_sym`.
- Lines 44-47: admissible/inadmissible PASS lines show the same independent-parameter numeric values (`0.02795067134961042` and `0.8105694691387022`), confirming F1's fix mirrors across engines.
- Line 61: `Stage 036 Mathematica audit passed.`

Both engines exit 0 and agree on every shared value to full numerical precision.

**Output freshness:**
- Script files: `.py` mtime 2026-05-25 23:51, `.wl` mtime 2026-05-25 23:51 (updated after F1/F2 fixes).
- Exec logs: `stage_036_sympy.log` and `stage_036_mathematica.log` mtimes 2026-05-26 00:20 — newer than both script mtimes, reflect the post-fix state, and contain the new independent-parameter M_mix numeric values.
- Saved `.txt` outputs under `scripts/output/` and `mathematica/output/` are still at 2026-05-21 17:40 (older than the script mtimes). The exec logs supersede them for verification purposes (per the prompt's "use the exec logs the orchestrator captured" rule). Flagged as a side observation; not blocking.

## Material-change assessment

`material_change`: false.

The edits change only the test-witness construction (F1) and add explanatory comments (F2). No derived closed form is altered: `G(xi,delta)`, `F(xi,delta)`, `M_mix`, `R_target`, `dG/dxi`, `G_max`, the near-onset series, and the symbolic-kappa identity `F = R_target_sym` all produce identical printed forms before and after (sympy log lines 10-13, 19, 22, 39; mathematica log lines 10-13, 20-21, 29, 52). No downstream unit consuming stage-036 outputs would see a different value. The exec logs confirm engines still agree at every shared check.

## Side observations (non-blocking)

- Saved `.txt` outputs in `scripts/output/` and `mathematica/output/` for stage 036 are stale (mtime 2026-05-21 17:40, predating the 2026-05-25 23:51 script edits). Verification used the fresh exec logs at `redteam/exec_logs/stage_036_*.log` instead, which is the prompt's prescribed source. Recommend the orchestrator refresh the saved outputs at some point so that an independent reader of `output/*.txt` sees the post-fix witnesses (`0.0280` / `0.8106`) rather than the legacy `53/145` / `82/145`. Not a verification blocker.

## Verdict justification

Both findings are mechanically resolved exactly as the directive prescribed: F1 replaces the `G_sample ± 1/10` witness pre-shift with independent-parameter `M_mix` evaluations from `Mmix_expr`/`mMix`, and the exec logs confirm the new numeric witnesses (`0.0280` admissible, `0.8106` inadmissible) straddle `G_sample = 0.4655` non-trivially in both engines; F2 inserts the four definitional-self-consistency comments at the named insertion points with verbatim wording, and the underlying assertions continue to pass with residual 0. Both scripts exit 0, engine agreement holds at every shared value, and no collateral edits or regressions are visible in the diff. `material_change: false`.
