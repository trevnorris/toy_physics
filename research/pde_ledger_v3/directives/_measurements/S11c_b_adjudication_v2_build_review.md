# S11c-b adjudication layer v2 — BUILD review, two legs, consolidated (SOUND; the engines genuinely differ)

**Artifact:** the v2 extension of `scripts/S11c_b_adjudicated_comparison.py` (Bridge D = `PROFILE_GRADE_SUBS`;
strong-exact vs weak-density IBP classifier; protected-atom gating; explicit route order). **Legs (Codex-written
⇒ fresh Claude agent + Grok):** prompt `_legs/S11c_b_adjudication_v2_build_review.md`. Raw: Grok
`~/.s11_build/S11c_b_adjud_v2_build_grok.txt`; Claude-agent scratch `/tmp/.../scratchpad/` (`prod_run.out`,
`independent_verify2.out`, `scale_capability.out`, `big_scale_test.out`, `abl_drop*.out`, `final_checks.out`).
**Both legs: SOUND, NO BLOCKING FINDINGS.** Route counts (both reproduced): `MATCH=38, FLAG=12,
REPRESENTATIONAL_DIVERGENCE=0, DIVERGENCE_INCOMPLETE=0, RESIDUAL_BULK=8, PROTECTED_UNREDUCED=32,
STRUCTURE_INCOMPLETE=57, COVERAGE=84` (231; multiset-equal; JET 231/0).

## CRUX — `REPRESENTATIONAL_DIVERGENCE=0` is REAL, not a broken classifier (both legs, independent)
- The Claude agent reimplemented the Euler/divergence classifier FROM SCRATCH (fresh total-derivative + Euler
  operator; exact random Gaussian-rational point zero-test; no shared code/canonicalizer) and certified all 8
  RESIDUAL_BULK as genuine bulk — each has a nonzero variational (Euler) derivative on MULTIPLE test fields
  (v_e_W_s11cb, v_theta_s11cb, u_T_3, psi_L_s11cb); a nonzero Euler derivative is a rigorous certificate that NO
  flux V exists.
- Both legs confirmed the classifier is CAPABLE at production scale: known coupling-family total divergences
  (Grok: `div(A·v)`, `div(W_bg·A·v)`; agent: up to a 5218-op divergence) → `REPRESENTATIONAL_DIVERGENCE` with a
  verified `V` (`R − Σ∂_iV_i == 0`); the `A_T·curl` bulk control → `RESIDUAL_BULK`. So the 0 is genuine
  non-divergence, not incapacity.

## No FALSE MATCH / no over- or under-reduction (both legs)
`--drop-bridge-a` and `--drop-bridge-d` leave MATCH=38 unchanged (no bridge manufactures a MATCH). Bridge D =
imported `PROFILE_GRADE_SUBS` object identity, 12 entries, `sigma_W ≢ W_0·eta_bg`, jets retained (JET 231/0);
legacy divergence reduction DISABLED, raw `case.value` consumed, no `DIVERGENCE_REDUCED` leak; the V-certificate
uses `BridgeD(R_preD − Σ∂_iV_i)` (the non-commutation is real: `(sigma_W − W_0·eta_bg)·w1_d1·φ ≠ 0`); strong-row
HeldDiv EXPANDED not dropped (0 unexpanded after transform). Protected 32 atom-gated; ENERGY_BASIS_COUNT exact;
accounting 231 multiset-equal; bad args exit 2; hygiene clean.

## ⭐ THE MEASURED CROSS-ENGINE RESULT (the finding — rule 6)
- **AGREE (38 MATCH):** ADMISSIBILITY operator(16) + support(20) + ENERGY_BASIS_COUNT(2), rename-level.
- **Representational (32 PROTECTED_UNREDUCED):** SLAB rows + MU_THETA differ ONLY by the protected quotient
  representatives (07/10, gamma-DivGrad) — the ENERGY_BASIS non-uniqueness, not a finding.
- **20 GENUINE cross-engine differences** (independently nonzero, persist without either bridge):
  - `ADMISSIBILITY_OPERATOR_OPERAND/BODY_FORCE/THETA` (4): `κ_θ_W·σ_W/(L_W·W_0)·(η_bg·w1 − 1)·Σ_i ∂_{ii}w1`
    (a background-Laplacian `∇²w1`; PY carries it, WL operand 0).
  - `SLAB_OPERATOR_TERM_ORIGINS/ADVECTIVE` (4): `−ρ_br·(1+η_bg·w1)·Σ_i ∂_i u_{i,t}`.
  - `SLAB_OPERATOR_TERM_ORIGINS/KINETIC` (4): `μ_W·e_W_tt·W_0²·(2η_bg·w1 + η_bg²·w1²)` = `μ_W·e_W_tt·(W_bg² − W_0²)`.
  - `COUPLING_KERNEL` (8 RESIDUAL_BULK): genuine non-IBP-removable bulk `A_T·(∂v·∂w …)` differences.

## ⚠ SYSTEMATIC SIGNATURE (orchestrator, to verify next)
The 20 genuine differences are all HIGHER-BACKGROUND-ORDER terms (`η_bg²·w1²`, `∇²w1`, `W_bg²` vs `W_0²`). WL's
`truncateScalar` (`mathematica/...audit.wl` ~L102-125) does `Series[etaBg,0,1]` then `Series[sigmaW,0,1]` — WL
TRUNCATES to first background order while PY RETAINS full background-amplitude order (which §3d says is intended:
"∇²W_bg retained at background-AMPLITUDE order, not discarded by a derivative-order cap"). ⇒ hypothesis: the
disagreement is one systematic thing — WL background-order truncation vs PY retention. VERIFY by computation
(are all 20 residuals killed by dropping η_bg²/second-background-order?), then adjudicate which engine is
correct per the spec (candidate: a systematic WL truncation gap).

## Non-blocking (carry to step record)
N1 `_normalise_exact` skips `cancel` >2000 ops — INERT here (all 8 decided at the Euler stage; none reach the
homotopy path); harden with a rational-point zero-test fallback. N2 record the "3-D compact-support IBP"
assumption (the classifier uses full 3-D — the stricter test). N3 the printed `euler_signature` shows one base
(breaks after first nonzero) but the verdict is over-determined (multiple nonzero bases).

## Verdict
Layer SOUND; the adjudication is correct. The engines AGREE on admissibility and genuinely DIFFER on the
operator/coupling in 20 characterized, higher-background-order ways. NEXT: verify the systematic-truncation
hypothesis, adjudicate which engine is right (rule 13), then the honest step record + S11c card + close.
