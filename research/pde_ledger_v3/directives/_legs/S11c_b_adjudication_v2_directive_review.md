# Independent review — S11c-b ADJUDICATION v2 build directive (decision list, pre-build; round 2)

## Artifact
`research/pde_ledger_v3/directives/S11c_b_adjudication_v2_build_directive.md` — orchestrator-written; extends the
committed adjudication layer with Bridge D (the engine's `PROFILE_GRADE_SUBS`) + an IBP/divergence classification
restricted to weak-pairing DENSITIES. DECISION-LIST review BEFORE the builder (no code yet). Round-1 legs removed
a wrong jet map (naive chain rule → now `PROFILE_GRADE_SUBS`) and an over-broad divergence classifier (now
density-only, strong operators exact). Re-derive; do not assume clean. The worst failure is a bridge or
divergence step that MANUFACTURES false agreement.

## Context (read all)
- The directive above; committed layer `scripts/S11c_b_adjudicated_comparison.py` + its review
  `_measurements/S11c_b_adjudication_build_review.md`; round-1 record `_measurements/S11c_b_adjudication_v2_directive_review.md`.
- Engine sources: `scripts/S11c_b_brane_operator_sympy_audit.py` (`PROFILE_GRADE_SUBS` ~L648-668,1862-1868;
  `sigma_W`/`eta_bg` L156-157; strong `variationalSource`) and
  `mathematica/S11c_b_brane_operator_mathematica_audit.wl` (anchoredWidth L448,535; variationalSource
  L274-277,793-795; COUPLING density L1075-1124).
- Comparator `scripts/S11c_b_cross_engine_comparator.py` (divergence/Euler machinery ~494-582; HeldDiv drop
  ~547-552); v1 `scripts/S11c_b_handcoded_comparison.py` (protected non-map L203-210).
- Baseline FLAG residuals `/tmp/S11c_b_adjudicated_run.out` (if present).

## Required checks — a finding only if it changes what gets built or what may be claimed

1. **Bridge D = PROFILE_GRADE_SUBS, correct + complete.** Confirm the directive's map is exactly the committed
   `PROFILE_GRADE_SUBS` (W_bg + mu_R_bg zero-jet via `eta_bg`; first jets via `sigma_W`; second jets via
   `sigma_W/L_W`; mu jets via `mu_R·sigma_W/W0`; density bookkeepers). Confirm `sigma_W ≢ W_0·eta_bg` is barred.
   Confirm it RETAINS all `w1_profile`/`m1_profile` jets (no order lowered). Flag any entry that differs from the
   engine dict, any missing same-class field, or any over-reach.

2. **The strong-vs-weak split — is the whitelist right?** The directive compares STRONG families exactly
   (SLAB, SLAB_ORIGINS, MU_THETA, ADMISSIBILITY_*, kinetic/container) and runs the divergence classifier ONLY on
   the COUPLING_KERNEL density. Verify from the engine construction: (a) is COUPLING_KERNEL genuinely a formed
   weak-pairing density (WL:1075-1124) where `∇·(ψV)` boundary-equivalence is valid? (b) Is any family on the
   STRONG list actually a weak density (wrongly excluded), or any OTHER family a weak density that should also be
   divergence-eligible (missing from the whitelist)? (c) Does the directive prevent any HeldDiv drop on strong
   rows?

3. **Divergence classifier rigor — can it still manufacture false agreement?** `REPRESENTATIONAL_DIVERGENCE`
   requires a printed, verified 3-component `V` (`R − Σ∂_iV_i == 0`), not just an Euler signature; unsupported →
   `DIVERGENCE_INCOMPLETE`/error. Check the three fixtures (product-rule bulk `a·φ_d1`→BULK; `a_d1·φ+a·φ_d1`→
   DIVERGENCE with `V=(aφ,0,0)`; a divergence in a strong family → ineligible) actually exclude over-reduction,
   and that one code path serves fixtures + production (fixture-local field registry). `--drop-divergence` a real
   ablation?

4. **Protected atom-gating.** Confirm any residual carrying `gamma_s11cb_{w_bg,mu_r_bg}_{07,10}` or
   `gamma{Width,Modulus}DivGrad{Theta,Ew}` is routed `PROTECTED_UNREDUCED` (raw), and ENERGY_BASIS never enters
   Bridge D or the divergence classifier. Flag any path where a protected representative could be folded.

5. **Leak (rule 5) + accounting DoD.** No expected classification/value for any case (esp. no forecast that the
   admissibility THETA or the coupling residual is representational or a finding). Case-ID multiset accounting +
   drop protection value-free and able to FAIL. Fixtures use placeholder fields.

## Method
Read the engine sources + the comparator divergence machinery FIRST; verify Bridge D against `PROFILE_GRADE_SUBS`
line-by-line and the strong-vs-weak split against the `variationalSource` construction; then judge. Show what you
read. Numbered findings (blocking vs non-blocking), file+line, concrete correction. If sound, say so and name
what you verified. ⛔ Reading + tracing; do not spawn Mathematica.
