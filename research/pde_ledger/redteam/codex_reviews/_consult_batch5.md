# Batch-5 math consult — Claude+Codex (read-only)

Stages 137, 139, 142 de-tautologization directions. One read-only `codex-chat` consult
(`-s read-only -C <root> --ephemeral`). Prompt: `redteam/tmp_prompts/consult_batch5.md`.
**All 7 questions: CONCUR. Zero DISPUTE.** Raw transcript deleted (transcript-bloat policy);
verdicts + grounding distilled below.

This is a math-coverage / de-tautologization resolution delegated to Claude+Codex
(see `feedback_claude_codex_resolve_math`); none of the items changed the conceptual nature,
so none was escalated to the user.

## 137 — core-to-mouth gain map (matrix-Schur de-tautologization)
- **Q1 (F1 Schur independence + feasibility): CONCUR.** Deriving `rho_c, sigma_c` from
  `M_core = [[K_s,λ],[λ,−K_q·D]]`, `v=(g_s,g_q)`, `δ = vᵀM_core⁻¹v` is genuinely independent of
  the hand-assigned values (entries are physical primitives), and reuses the ALREADY-VERIFIED
  owner stage 114's identical `M,c` construction (`scripts/...stage114...py:25-30`) → symbolically
  feasible under the 600s cap.
- **Q2 (F2 static-limit envelope + maps): CONCUR.** The Stage-97 core content lives in the
  concrete-core note filed as stage114 (`notes/...stage114...md:15-33, :43-52`), which defines
  `D_W^bare` and the Schur outlet; comparing the full z-dependent reduced envelope against the
  matrix source on `D_W_bare(z)` with `κ_c=κ₀/(1+r_c)`, `γ_c=γ₀/(1+r_c)`, `r_c=λ²/(K_sK_q)` is the
  right non-tautological object. (Orchestrator also hand-verified the envelope identity.)
- **Q3 (F3 outlet sign/factor): CONCUR — sign explicitly checked.** Stage 137 derives
  `Π = (L/Θ_σ)[ρ_c U_s'(0) − σ_c U_q'(0)]` (`notes/...stage137...md:40-48`), so
  `M_q = −L·σ_c/Θ` has the correct NEGATIVE sign and `L/Θ` factor for the mixed contribution.
  The directive's `M_q*S_q == −L·σ_c_schur·S_q/Θ` + non-vacuity (`+M_q` vs `−M_q` differ) is sound.

## 139 — actual Family-1 mouth gains
- **Q4 (R1 outlet residual is X−X): CONCUR.** `Ms=Π/(1−Rq·Sq)`, `Mq=−Rq·Ms`
  (`scripts/...stage139...py:11-18`) make the residual asserts (`:51-55`) identities. Resolution:
  demote to print-only + add an independent `S_q(Π_*)` recompute from the Stage-134 kernel
  (κ=π/2) vs the imported `Sq_star`.
- **Q5 (R2 R_q^comp=1/4 is definitional): CONCUR.** "The exact compensated branch sets
  `g_c = r_F1 − ½√(1+r_F1²)` and `R_q = 1/4`" (`notes/...stage139...md:91-100`). It cannot
  "emerge" non-tautologically (and is branch-blind: `g_±` both give 1/4). Resolution: RELABEL
  `R_q^comp=1/4` as definitional-consistency; load-bearing falsifiable checks = `r_F1` (vs 121),
  `R_q^nat ≠ 1/4` (real branch selection), and the **`g_-^F1` value** `0.758035078944662826919680890414`
  (discriminates the lower branch `≈0.758` from the upper `≈2.79`, which `R_q=1/4` cannot);
  Mathematica computes `g_minus` DIRECTLY as the closed form (sanctioned mirror), NOT via
  `Solve[(gc−rF)²==(1+rF²)/4]`.
  - *Orchestrator refinement (beyond the consult):* framed the `g_-` value anchor honestly as
    **branch-discrimination** (g_- is definitionally `rF−√(1+rF²)/2`, not an independent
    derivation), since both `g_±` satisfy `R_q=1/4` — the value pins the lower branch.

## 142 — self-consistent mouth branch
- **Q6 (R1 R_q(Π_*)=1/4 tautological; which Π_* anchor): CONCUR.** `Π_*=nsolve(gPi==gminus)`
  (`scripts/...stage142...py:52-58`) makes `R_q(Π_*)=1/4` solver-consistency. Codex confirmed
  142's canonical point "is the same lower Family-1 point carried from the mouth-bias stages," so
  anchoring on **stage 131's** independently-derived `Π_*` (`1.50882951349315558300555075595`,
  cleared-denominator FindRoot, batch-4-verified) — NOT 142's own nsolve output (which diverges at
  digit ~16) — is the valid non-tautological anchor. Relabel `Rq_star=1/4` as solver-consistency;
  keep the F2-kept decimal anchors.
- **Q7 (F2 projection-integral independence + feasibility): CONCUR.** Stage 129 gives the source
  law `σ_Π ∝ Π e^{−Π z}/(1−e^{−Π})` (`notes/...stage129...md:107-119`); Stage 130 defines the
  cosine projection (`notes/...stage130...md:16-30`). Replacing the self-series check with
  `gPi := ∫₀¹ σ_Π(z)·cos(πz/2)dz` is genuinely independent (built from σ_Π + cos, not from `gPi`),
  equals the closed form, and closes in elementary form. (Orchestrator hand-verified the integral
  = `2Π(2Πe^Π+π)/((4Π²+π²)(e^Π−1))` exactly.)

## Net
137 & 134 directives final as drafted. 139 F2 rewritten (definitional relabel +
branch-discrimination `g_-` anchor + direct-closed-form Mathematica `g_minus`). 142 F1 repointed
to stage-131's independent `Π_*`. No conceptual changes; no user escalation. Proceed to the fix
loop (2-parallel Codex: 134∥137, then 139∥142; never >2 concurrent `math -script`).
