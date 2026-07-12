# Design-review this REVISED computation directive before execution. READ-ONLY. Do NOT write or run code.

Read:
- `software/em_charge_attribute/directive_monopole_softness_estimate.md` (the REVISED directive)
- `docs/em_charge_attribute_requirements.md` (context: §4 softness drill, §5 the gate)

Context: this directive computes whether supersolid phonon/Goldstone softness pushes a 3+1D compact-U(1) Coulomb phase across the confinement boundary — the decisive gate for the "EM = decoupled charge-attribute" reframe (a 3-AI panel judged the reframe "conditional-needs-X", not a no-go). An earlier version of this directive was reviewed and got `DIRECTIVE_NEEDS_FIXES` with these catches, now supposedly folded:
1. δM/M is only a perturbativity diagnostic, not a confinement criterion → use the renormalized monopole fugacity/worldline tension vs the KNOWN compact-U(1) phase boundary.
2. Controls were partly vacuous (2+1D can't come from mass-propagation; rigid only proves δM→0) → add a rigid 3+1D coupling sweep that reproduces BOTH the weak-coupling Coulomb and strong-coupling confined phases; keep neutral-decoupled and charged→Higgs; order-one/near-boundary → INDECISIVE.
3. g∝1/a is dimensionally suspect (3+1D gauge coupling is dimensionless); M_mon∼J_ice is order-scale only; Goldstone mixes with longitudinal strain (not an independent δa); T=0 core action vs finite-T free energy must not be mixed → require a specified microscopic magnetoelastic model with a proper strain correlator + BZ cutoff.
4. The answer swings with the undetermined exponent α=d ln M_mon/d ln a → require α DERIVED from the model; output a sensitivity surface; verdict only if robust.

Assess:
A. Are catches 1–4 genuinely addressed by the revised directive (cite the sections)? Any that are only cosmetically addressed?
B. Is the revised acceptance now genuinely able-to-fail and decisive, with a valid (non-punitive) INDECISIVE outcome?
C. Any REMAINING physics error, hidden assumption, or under-determination that would make the computed verdict unreliable?
D. Is the computation tractable as specified (semi-analytic Gaussian integrate-out + known-boundary comparison, dual-engine), or does it still hide a step that needs the full lattice Monte Carlo?
E. The single most important remaining change, if any.

Output: terse, per-item A–E, then a one-line verdict: `DIRECTIVE_READY` or `DIRECTIVE_NEEDS_FIXES` (list them). Your entire final message is the deliverable.
