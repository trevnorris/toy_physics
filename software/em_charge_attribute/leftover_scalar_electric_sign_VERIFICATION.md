# Verification — leftover-scalar electric-sign build (2026-07-21)

**Verdict: `RESULT_CONFIRMED` — 4 independent legs, no rig / no fidelity error / no target-blind leak. Bankable as stated.**

The build itself over-ran into (disregarded) self-verification; the legs below are the **independent** verification the orchestrator arranged.

## Landing (verified)
`NO_NATIVE_CLAMP; ALGEBRAIC_SIGN={V:range(attract|null|repel), M:range(attract|null|repel), J:attract}; BOLT_ON_DEFERRED_TO_R1`, Tier-A passes, Tier-B deferred.

## The four independent legs
1. **Arbiter — SymPy re-run** (orchestrator, scripts unchanged): `PY_EXIT=0`, `ENGINE_AGREE=PASS`, 13/13 teeth PASS with ablation firing at each own assert, 384-case verdict table total. Landing reproduced.
2. **Arbiter — Mathematica re-run** (independent `WolframKernel`, scripts unchanged): reproduces `Q1a=NO_NATIVE_CLAMP`, the same landing, `ENGINE_AGREE=PASS`, all teeth. (No orphan kernel / seat leak.)
3. **Grok compute-verify** (independent SymPy): `RESULT_CONFIRMED`. Reproduced `m_uu=4/7, m_ug=−2/7, m_gg=8/7`, `det m=z_g²/D>0` (positive-definite); `h`-attract / `u_L`-repel-for-V·M / interference-zero-V·M / J-strictly-attract; `NO_NATIVE_CLAMP` able-to-fail; Tier-A clean; faithful transcription; target-blindness intact.
4. **Fresh-agent adversarial-recompute** (independent from-scratch SymPy, ran every tooth's mutation): `RESULT_CONFIRMED`. Every expression/landing/Tier-A/tooth reproduced; Q1a injection flips `NO_NATIVE_CLAMP→NATIVE_CLAMP_EXISTS` (evaluator live, not hardcoded); numerically confirmed the V/M net sign flips attract→repel at `q/g ≈ 0.3`.

## The interpretive crux — both substance legs agree
**The V/M "range" is HONEST: the committed model does NOT fix the `h`-vs-`u_L` ratio.** The committed mouth source drives only `h` (the pathA_38 Coulomb, which attracts, `−g²`); the clamp datum `u₀`/`q` is the bolt-on — **absent from G0** (the zero-ledger forbids any `u_L` source/datum coupling: `r_B·u_L=0`, no direct `u_L` source, `δE_g/δu_L=0`). So `q/g` is unconstrained, and the net V/M sign `= sign(m_uu q² − m_gg g²)` is a **free calibration knob**:
- **the `u_L`-clamp CAN give like-repel `1/R²`, but only as a CALIBRATED postulate (a clamp strong enough to dominate the committed `h`-attraction) — the electric sign is NOT earned/determined.**
- the more throat-like fixed-**source** ensemble (J) is **unconditionally attract** — so if anything the analysis leans *against* native repulsion.

## What is DECIDED vs deferred
- **DECIDED (verified):** G0 has no native nonrelaxable signed `u_L` datum; the three conjugate ensemble structures + their algebraic sign ranges + the `1/R²` falloff + Tier-A consistency are fixed as above; the sign is a free clamp-to-`g_χh` ratio, not determined.
- **DEFERRED to R1/R3:** whether a consistent bolt-on holds the clamp; which mouth ensemble/normalization a throat realizes; whether the throat structure forces the clamp strong enough (`q/g` large enough) to dominate; the dressed non-perturbative two-body response; and all Tier-B gravity/sleeve/momentum checks.

## Non-blocking notes (from the reviewers)
- `z_g,z_b` are parameterized positive Schur/escape factors (not solved from the assembled Robin problem) — correctly scoped to R1; no sign conclusion depends on the deferred solve (V/M depend only on `ζ²,z_b²`; J negative-definite regardless).
- A few Tier-A checks (`transverse_decoupled`, `u2_nonselection`, `zero_ledger_preserved`) are source-consistency bookkeeping rather than computed PDE properties — acceptable at this analytic tier; each able-to-fail via an introduced token.
- Prose in the result note could state the "range is over the free clamp-to-`g_χh` ratio" more explicitly (both reviewers; optional, not a correctness issue — captured here instead).

## Bottom line
The leftover-scalar-as-charge candidate is **NOT falsified** (Tier-A consistent; like-repel is reachable) and **NOT earned** (the sign is a calibrated knob, matching the audit's I2 "sign IMPORTED"). The build decisively **proved no native clamp exists in G0** and **mapped the electric sign to a sharp R1 question**: does a consistent bolt-on hold the clamp, and does the throat force it strong enough to beat the committed `h`-attraction? Verified target-blind, so the (modest) result is trustworthy.
