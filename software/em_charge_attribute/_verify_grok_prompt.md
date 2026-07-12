# Independent physics verification of a PASS verdict. READ-ONLY. Falsification is the goal — try to break the PASS.

A computation claims `EMERGENT_EM_WITH_DUAL_SIGN`: that a model's ±w "throat" can be embedded as an emergent electric charge in a quantum-spin-ice compact-U(1) Coulomb phase, with the correct EM dual sign, non-circularly. Your job: verify it is physically valid, or find the flaw.

Read:
- `software/em_charge_attribute/reports/emergent_em_construction.md` (the claimed result)
- `software/em_charge_attribute/directive_emergent_em_construction.md` (what was required)
- the scripts: `software/em_charge_attribute/emergent_em_sympy.py`, `emergent_em_dual.wl`

Verify, adversarially:
1. **Is the emergent-charge construction genuinely non-circular?** The claim: charge = ice-rule (div E) defect, integer spectrum {-2..2}; UV global S^z U(1) is broken by a small h_i S^x perturbation so no microscopic matter U(1) is used. Is this a legitimate anti-circularity argument, or does it still smuggle in the U(1) / Gauss law somewhere?
2. **Is the ±w embedding real, not a relabel?** The claim: the physical charge is a COMPOSITE flux-dressed string endpoint (τ_σ dressed with a Wilson string of S_ℓ), additive via the divergence sum, with Gauss-preserving hopping (current from ∂H/∂A, not j=qv). Is this genuine, or does it hide a hand-imported charge? Is the incidence identity Bp_γ = (+1,…,−1,…) correct?
3. **Is the dual sign correctly derived?** S_eff = ½∫[ρ(1/εk²)ρ − j_T(μ/k²)j_T] from ONE Maxwell field → static like-charge repulsion + transverse-current magnetostatic attraction; two moving like charges NET REPULSIVE. Correct? Are the signs read as static interaction energies (not action-sign artifacts)?
4. **Are the controls physically meaningful and non-vacuous?** scalar → FAIL_SCALAR (attractive density, no transverse channel); ring-off → no photon; defects-condensed → Higgs (m_A²=g²v²); the FAIL_CHARGE_POSTULATED guards. Any that is tautological / cannot actually fail?
5. **Is the honest scope accurate, or does it overclaim anywhere?** It concedes: requirement A (the internal DOF) is a POSTULATE; the existing continuum throat PDE does NOT derive it; deconfinement is CITED (HFB) not computed; gravity/softness/4+1D deferred. Is anything claimed that isn't actually established? Is the HFB citation used correctly?
6. Any physics error, and the single biggest caveat a reader must know.

Output: terse per-item 1–6, then a one-line verdict: `VERIFICATION_CLEAN` / `VERIFICATION_ISSUES` (list). Your entire final message is the deliverable.
