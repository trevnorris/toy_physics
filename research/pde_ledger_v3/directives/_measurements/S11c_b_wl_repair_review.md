# S11c-b Wolfram engine — repair round-1 review, two legs (serialized) — BOTH CLEAN

**Artifact:** `mathematica/S11c_b_brane_operator_mathematica_audit.wl` (Codex-repaired round 1, functions
`constructFullFieldBackgroundEnergy`/`backgroundBalanceFromModel`/differential-field-substitution/weak-pairing
helpers/`extractCoupling`/`simplifyWeakKernel`/`planeWaveCoefficient`). **Legs:** Grok + fresh Claude agent,
serialized (Mathematica 2-seat), prompt `_legs/S11c_b_wl_repair_review.md`. **Raw:** Grok
`~/.s11_build/S11c_b_wl_repair_grok.txt`; Claude-agent `…/scratchpad/s11cb_wl_review/`; Grok
`/tmp/s11cb_repair_review_grok/`. Each wrote an independent derivation + FORM ablations under `timeout 600`,
one kernel at a time (peak RSS ≤1 GB; no orphans). **Verdict: BOTH legs CLEAN — no blocking findings.**

## W1 — admissibility background-order balance — FIXED (both legs, ablation)
`constructFullFieldBackgroundEnergy` full-fields ONLY the gradient slot (`gradient[fullWidth] = ∇(W_bg+δW)`,
L528–576) and keeps the scalar perturbation fields (`thetaVariation`, `localEw`) as perturbations (∝
background order → vanish at first variation). `backgroundBalanceFromModel` takes the ε⁰ first variation at
𝔅⁰. Verified: emitted `E_W` body force carries only the `κ_W·σ_W·W_bg` second-jet (∇²W_bg) content; `U`→0,
`THETA`→0; uniform limit (σ_W→0) → 0 (all leaves `PossibleZeroQ`); FORM ablation stripping the background
gradient from `fullWidth` collapses the operand to 0 (proving the content IS the full-field ∇W_bg, not the
ε→0 wave limit). **No scalar over-promotion** (the SymPy B1 defect); WL is the correct §3d reference.

## W2 — weak solenoidal/irrotational kernel extraction — FIXED (both legs)
`extractCoupling` substitutes genuinely solenoidal (`curl[A]`) / irrotational (`gradient[φ]`) trials into the
operator (`linearRestrictedOperator`), so the split acts on the gradient content. Uniform limit: forward &
reverse pairing densities → 0 (the old diagonal-thickness survivor — `HAS_KW/BRHO/KAPPAW` — is gone;
THICKNESS_ROW telescopes to 0 via div-of-curl). Non-uniform: genuine transverse↔thickness coupling ∝ first
background jet (N15 spurions `gammaWidth*`/`gammaModulus*`). Breaking the operator zeroes the kernel
(extraction is operator-dependent).

## W3 — pairing-based operator-block adjointness — FIXED (both legs)
Residual is `forwardDensity − reverseRelabeledDensity` (`PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP`), NOT the
scalar-Hessian Clairaut. Verified non-tautological: on a hand-built non-adjoint operator the residual is an
irreducible bulk term (LeafCount 31253; Claude-agent's total-divergence tool confirms it is not a boundary
term); the old `D[D[E,·],·]−D[D[E,·],·]` is ≡0 for any operator. Uniform-limit residual → 0; operator-dependent.

## No regression (both legs)
N15 basis count 26 (8 WJET + 8 MUJET spurions), operator, and operator→kernel data dependence all
ablation-clean; no VERDICT/PASS/FAIL; relationals `Inactive[Equal]`; `ADMISSIBILITY_RESIDUAL` a genuine A−B
(bending force − `forceHoldEw`), not A−A.

## Interpretation notes for the step record (NOT defects)
- **W3 adjointness is a density modulo compact-support IBP** ⇒ a nonzero density ≠ non-adjoint (an adjoint
  operator yields a total-divergence density, ∫=0). The engine correctly prints the density and defers the
  reduction (rule 2); the comparator/step record must reduce modulo IBP before reading it.
- The real operator's adjoint residual is **nonzero because the `Λ_X` face response is dissipative** —
  a genuinely non-self-adjoint operator, a correct physics emission (not a bug). A guaranteed self-adjoint
  witness is not cleanly constructible (the mass constraint eliminates θ, so even the dissipation-free
  operator is not a naive symmetric Hessian).

## Status
WL round-1 repair COMPLETE and review-clean (committed `10eaa75c`). No WL round-2 needed. Awaiting the SymPy
B1 round-2 fix + its re-legs; then both engines clean → run transcripts → T7 comparator → step record.
