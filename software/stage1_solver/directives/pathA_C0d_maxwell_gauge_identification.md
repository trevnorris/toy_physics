# Directive pathA_C0d — Is the Maxwell-lane near-null subspace GAUGE? (proper ∇λ basis, not one ansatz)

**Status:** READY (Claude-authored 2026-06-20; design-review = SOUND-WITH-FIXES → all fixes applied [weighted gauge term
`ξ⁻¹ grad(Z·D_A·A)` not raw ∇·A; dropped the ξ-rescale check; added a negative control; pinned the weighted-residual
tolerance 0.1]; confirm-pass = SOUND-AS-IS 2026-06-20. User GATED "C0d: confirm Maxwell modes are gauge" → EXECUTION is the
first post-/compact action; prompt pre-written at `_scratch/pathA_C0d_execute_prompt.md`.) Follows `pathA_C0c` (verdict MIXED, fidelity-verified): the worst near-null mode IS the global U(1)
PHASE gauge mode (CONFIRMED); the OTHER 4 near-null modes are ~99% in the Maxwell A-field lanes (`ar`/`aw`), labeled
`UNEXPLAINED_STIFFNESS` ONLY because the single crude `maxwell_residual_gauge` probe (one fixed `λ`) did not span them
(overlaps 0.11–0.77). Step of option C, task #78.

**Date:** 2026-06-20
**Owner:** Codex (codes + iterates until exit 0; extends the C0c diagnostics). Claude reviews after.
**Trigger:** the fidelity agent flagged that the 4 Maxwell-lane modes are almost certainly Maxwell-sector GAUGE modes a
proper basis would span — and must NOT be called "stiffness" / "production-solver-required" until tested against a FULL
A-sector gauge basis. This directive runs that test.

## Why this is the decisive, cheap next step
If the 4 Maxwell-lane modes lie in the A-sector gauge subspace, then the ENTIRE τ≈0.029 wall is a GAUGE artifact (U(1)
phase + Maxwell A-sector) — not physics — and a CHEAP combined gauge-fix/deflation could dissolve it. The test is cheap:
build the proper A-sector gauge subspace and project the measured modes onto it; no new solve.

## Key physics (read before designing)
The coded residual's vector gauge-CONTROL term is NOT raw `ξ⁻¹∇·A` — it is the WEIGHTED form
`ξ⁻¹ grad(Z(w)·D_A·A)` applied in the `ar/aw` residual blocks, where `D_A = axisymmetric_vector_divergence(ar, aw)` and
`Z(w)` is the localization weight (confirmed in C0d design-review, `operators.py:429`). The test MUST use this exact
weighted object, not raw divergence. Consequences:
- A PURE-gradient `δA = ∇λ` is generally PENALIZED by `ξ⁻¹ grad(Z·D_A·∇λ)` — so a generic `∇λ` is NOT an exact null mode
  (this is why the one crude probe gave `‖J·g‖=4.5e-3`, failing the exact gate). That does NOT mean the modes aren't gauge.
- The RESIDUAL gauge freedom the term does NOT fix is the part with near-zero WEIGHTED divergence `Z·D_A·δA ≈ 0` ⇒ NOT
  penalized ⇒ a NEAR-null mode with small-but-nonzero σ (consistent with the measured σ≈4.7e-11..2.8e-8, far above the
  phase mode's exact-zero σ≈7e-14). So the right characterization is: are the 4 modes spanned by the ∇λ subspace AND
  near-(weighted-)divergence-free — i.e. report BOTH the raw divergence `‖D_A·A[v]‖/‖A[v]‖` AND the normalized weighted
  gauge residual `‖grad(Z·D_A·A[v])‖/‖A[v]‖`, the latter being the actual penalized object.

## Stance (carry the C0c discipline)
DIAGNOSIS ONLY. Do NOT alter the faithful PDE operators, frozen physics, or `physical_export_permitted`. Do NOT implement
the gauge-fix/deflation or a re-crawl (NEXT step, gated on this). Evaluate on the EXISTING C0b/C0c saved states + assembled
Jacobian matrices (recompute the SVD modes from `runs/pathA_C0b_wall_diagnosis/matrices/*.npz` as C0c does). Reuse the C0c
overlap/projection machinery. CPU; `timeout 600` (split if needed; timeout → NOT_MEASURED); standalone `python3`; no
commentary `python3 -c`; YAML/markdown human output, JSON only for machine artifacts; chunk-1a/1b/1c gates must still pass.

## Work items

### C0d-1 — build the PROPER A-sector gauge subspace (not one ansatz)
Determine the A-lanes and the discrete gradient/divergence operators from `coupled_branch.py`/`operators.py` (what are
`a0,ar,aw`? which is the scalar potential vs the spatial `A` components?). Construct the discrete-gradient GAUGE SUBSPACE
`G = span{ discrete_grad(λ_k) embedded in the A-lanes }` for a BASIS of scalar fields `{λ_k}` (e.g. one per grid node, or a
well-conditioned smooth basis spanning the same range), respecting the actual boundary conditions. Orthonormalize `G`.
Also form the divergence-minimizing sub-basis `G_harm ⊆ G` — the residual gauge the WEIGHTED term `ξ⁻¹ grad(Z·D_A·A)` does
NOT penalize, i.e. `∇λ` with the weighted divergence `Z·D_A·∇λ ≈ 0` (NOT raw `∇²λ`). Report `dim(G)`, `dim(G_harm)`, and how
each was built.

### C0d-2 — project the 4 Maxwell-lane modes onto the gauge subspace
For each of the 4 near-null modes (recomputed from the saved matrices), report the captured fraction `‖P_G v‖²` (full ∇λ
subspace) AND `‖P_G_harm v‖²` (the residual part), plus the unexplained residual `1 − ‖P_G v‖²`. Also report BOTH the raw
A-sector divergence `‖D_A·A[v]‖/‖A[v]‖` AND the normalized WEIGHTED gauge residual `‖grad(Z·D_A·A[v])‖/‖A[v]‖` (the actual
penalized object; small ⇒ the residual gauge the penalty misses). Restrict the projection to the A-lanes (the modes are
~99% A-lane; report the small non-A remainder).

### C0d-3 — confirm these are PENALIZED-GAUGE, not genuine stiffness
Cross-check the modes are gauge-like using the C0d-2 evidence on the EXISTING saved matrices ONLY (do NOT re-assemble the
Jacobian at a different `ξ` — that changes a frozen parameter and is OUT of scope for this diagnosis-only run): a mode is
penalized-residual-gauge if it has high `‖P_G v‖²` AND a small WEIGHTED gauge residual `‖grad(Z·D_A·A[v])‖/‖A[v]‖` (it lies
in the gradient range but the `ξ⁻¹ grad(Z·D_A·A)` penalty barely sees it). State explicitly what in the data distinguishes
"penalized residual-gauge" (high P_G + small weighted gauge residual) from "genuine Maxwell-sector stiffness" (low P_G OR
large weighted gauge residual).

### C0d-4 — THE VERDICT (falsifiable, whole-cluster)
Classify each of the 4 Maxwell-lane modes `MAXWELL_GAUGE` / `GENUINE_MAXWELL_STIFFNESS` by the EXACT gate: a mode is gauge
iff `‖P_G v‖² ≥ 0.9` AND the normalized WEIGHTED gauge residual `‖grad(Z·D_A·A[v])‖/‖A[v]‖ ≤ 0.1` (state the chosen
tolerance; 0.1 is the default — also report the raw divergence for context). Overall:
  - **WALL_IS_ALL_GAUGE** — all 4 Maxwell-lane modes are MAXWELL_GAUGE (and C0c already confirmed the phase mode): the
    entire τ≈0.029 near-null subspace is gauge ⇒ recommend the COMBINED fix for the NEXT step (pin/deflate the U(1) phase +
    deflate or strengthen-gauge-fix the A-sector `∇λ` subspace), then re-crawl.
  - **MIXED_GAUGE_PLUS_RESIDUAL** — some Maxwell modes are gauge, a residual subspace is NOT (report its per-lane split +
    captured fraction) — the residual is the genuine-stiffness candidate.
  - **GENUINE_MAXWELL_STIFFNESS** — the Maxwell modes are NOT spanned by the gauge subspace (captured fraction low, not
    divergence-free) ⇒ real Maxwell-sector ill-conditioning → the harder solver question for the A-sector.
  - **DIAGNOSTIC_INCOMPLETE** — required evidence NOT_MEASURED (couldn't build `G`, divergence operator unavailable, etc.).

## Acceptance criteria (PASS/FAIL; exit-0 NECESSARY not SUFFICIENT)
1. `G` (and `G_harm`) built as a PROPER multi-`λ` discrete-gradient subspace (NOT one ansatz), with `dim` reported and the
   A-lane / gradient / divergence operators sourced from `coupled_branch.py`/`operators.py` (shown, not guessed).
2. Per-mode captured fractions `‖P_G v‖²`, `‖P_G_harm v‖²`, unexplained residual, AND BOTH the raw divergence
   `‖D_A·A[v]‖/‖A[v]‖` and the normalized WEIGHTED gauge residual `‖grad(Z·D_A·A[v])‖/‖A[v]‖`, recomputed from the saved
   matrices.
3. C0d-3 gauge-vs-stiffness discriminator stated, using ONLY the saved-matrix evidence (high P_G + small weighted gauge
   residual ⇒ gauge); NO ξ re-assembly (out of scope).
4. Exactly one C0d-4 verdict (or DIAGNOSTIC_INCOMPLETE) with falsifiable numeric `verdict_support` (captured fractions,
   raw + weighted divergences, thresholds) + the recommended combined-fix design if WALL_IS_ALL_GAUGE.
5. Anti-hardcode: BOTH a POSITIVE control (project a held-out known `∇λ` ⇒ captured fraction ~1, non-A lanes ~0) AND a
   NEGATIVE control (project a random / transverse / curl-type A-lane vector ⇒ captured fraction LOW) as unit tests.
   Faithful operators + frozen physics + export guard untouched (diff); no gauge-fix/deflation/re-crawl implemented; chunk
   gates pass; report + machine JSON.

**Fail conditions:** one-ansatz "basis"; hardcoded captured fractions/divergences; claiming MAXWELL_GAUGE without BOTH the
subspace AND divergence gates; calling modes "stiffness" without the proper-basis test; altering operators/frozen/export;
implementing the fix/re-crawl; masking NOT_MEASURED; raising the timeout cap.

## Out of scope
The combined gauge-fix/deflation itself; the re-crawl; the full-budget crawl; the production-solver decision; `pathA_22`.

## Review (orchestrator, after Codex)
Fidelity agent: is `G` a genuine multi-`λ` discrete-gradient subspace (dim ≫ 1, correct A-lane/grad/div operators)? are the
captured fractions + divergences GENUINELY computed (planted-gauge test) — not hardcoded? Adversarial agent: is the
gauge-vs-stiffness call SOUND (both gates), or is a low-divergence being over-read as gauge; if WALL_IS_ALL_GAUGE, is the
recommended combined fix correct + minimal? Diff-check faithful operators untouched. Then gate the NEXT step (the combined
gauge-fix + re-crawl).
