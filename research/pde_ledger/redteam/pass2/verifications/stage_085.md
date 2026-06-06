---
unit_id: 085
batch: III.5
verdict: verified
overall: verified
material_change: false
verified_at: 2026-06-05
---

# Verification — unit 085 (pass 2)

**Verdict: verified — clean (0 findings).** Audit verdict was `clean`; no directive, no Codex.
Stage 085 (`quadrupole_demand_cancellation`), not a checkpoint, both engines present.

The load-bearing checks genuinely exercise the `kappa0^2 = 8/pi^2` cancellation against the
`8/(pi^2 A)` loading factors (`Pi_tr − NQ·alpha_req/beta0`, `Pi_tr/C_mix − alpha_req/alpha_mix`,
the `zeta_req` product↔loading, `rho_alpha`, and the unblocked limit) — each residual is nonzero
if the constant is wrong, so none is identity-by-construction. The only literal, `kappa0^2 =
8/pi^2`, is anchored in notes §1 (from Stages 035–036). `eps_blk` correctly real-not-positive.
Self-labels canonical ("STAGE 085" both engines); upstream cross-refs (030/035–036/052/084)
deferred per numbering policy and correct anyway.

**Value reconciliation:** 13 deliverable values checked, 0 misaligned (product identities, ratio
cancellation, spectral forms, reduced demand law, unblocked limit, `kappa0^2=8/pi^2`).

Orchestrator independent re-run (reliability gate): SymPy exit 0, Mathematica exit 0; both
committed outputs **byte-identical to HEAD**. Engines agree exactly. No paper/notes edits.
