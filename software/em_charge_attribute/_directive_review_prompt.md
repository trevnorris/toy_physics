# Design-review this computation directive BEFORE it is executed. READ-ONLY.

You are reviewing a DIRECTIVE (not code) for a semiclassical physics estimate. Your job: is it well-posed, decisive, and genuinely able-to-fail — or will it produce a rigged/uninformative answer? Do NOT write or run code.

Read:
- `software/em_charge_attribute/directive_monopole_softness_estimate.md` (the directive under review)
- `docs/em_charge_attribute_requirements.md` (context: §4 softness drill, §5 the gate)

The physics goal: estimate the phonon-induced monopole-mass fluctuation `δM/M` for a compact-U(1) Coulomb phase on a supersolid lattice, to decide whether superfluid softness confines the gauge (kills EM) or is tolerable (EM survives). 3+1D. A unanimous 3-AI panel already judged this "not a fundamental no-go, conditional-needs-X"; this computation is the decisive gate.

Assess and report:
1. **Is the observable `δM/M` the right decisive quantity**, or is there a sharper/more correct criterion for the compact-U(1) deconfinement-under-phonons question (e.g. the fluctuation-renormalized monopole fugacity / free energy, a Lindemann-type criterion, the monopole worldline tension)? Name the most defensible criterion.
2. **Is the acceptance genuinely able-to-fail?** Are the mandated controls (2+1D must confine; rigid → survives; charged condensate → Higgs) the right able-to-fail teeth, or is one vacuous/tautological? Any missing control that would catch a wrong method?
3. **Any physics error or hidden assumption** in the directive's setup (the `g ∝ 1/a`, `M_mon ∼ J_ice`, Jensen fugacity-enhancement, δa/a from zero-point motion) that would bias the result? Is the semiclassical/Gaussian-fluctuation approach adequate for a first gate, and where exactly does it break?
4. **Is the scope honest and the model under-determined in a way that makes the answer arbitrary** (e.g. does `dM/da` / the `J_ice(a)` scaling need pinning, or the answer swings freely)? What must the directive force to be *derived* vs allow *declared* so the verdict isn't an artifact of a chosen exponent?
5. The single most important change to the directive before Codex executes it.

Output: terse, per-item, then a one-line verdict: `DIRECTIVE_READY` / `DIRECTIVE_NEEDS_FIXES` (list them). Your entire final message is the deliverable.
