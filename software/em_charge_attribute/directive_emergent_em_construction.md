# Directive — Construction B (v2): does the model's charge-attribute realize emergent EM (compact-`U(1)` + correct dual sign), with the `±w` throat as a *genuine* emergent charge?

**Revised after Codex + Grok design-review** (both `DIRECTIVE_NEEDS_FIXES`, convergent, constructive; both confirm **tractable — no Monte Carlo**). **Goal (one line):** show the model's charge-attribute realizes an **emergent deconfined compact-`U(1)` Coulomb phase (3+1D)** whose **emergent charge is the constraint *defect***, reproduces the **correct dual-sign EM signature**, and admits an **honest `±w`-throat embedding** — the "does the reframe produce EM at all" test, prior to the softness question. Ref: `docs/em_charge_attribute_requirements.md` §2 (A/B/D).

**What is actually being tested (both reviews).** The emergent-compact-`U(1)`-from-a-local-constraint *scaffold* is textbook (quantum spin ice / quantum dimers) — **CITE it, do NOT re-prove it.** The genuine able-to-fail content is whether the **model's own `±w` throat realizes that scaffold honestly** — as a *vertex defect* with additive charge from the divergence sum, flux dressing, and a Gauss-preserving hopping — or whether making it work requires importing a `U(1)` by hand (`FAIL_CHARGE_POSTULATED`). **That embedding is the centerpiece.**

---

## Pinned model (no cherry-pick)

- Pin **ONE** concrete literature Coulomb-phase Hamiltonian and state which: **3D quantum spin ice** (Hermele–Fisher–Balents 2004, `cond-mat/0305401`) **or 3D quantum dimer model** (Moessner–Sondhi 2003, `cond-mat/0307592`).
- **3+1D.** The deconfined Coulomb phase over a **finite** region is **CITED** from that model's literature — not re-derived.
- **[POSTULATE §2A]** the constituents carry the internal constrained DOF; what must be *derived* is the emergent EM + the honest `±w` embedding.

## What you must compute (requirement, not route)

1. **Emergent-gauge scaffold (analytic, from the constrained Hilbert space):** derive the emergent Gauss law, the emergent photon (2 transverse modes, 3+1D), and the emergent charge as the **constraint defect** — via *microscopic moves*: the defect/vertex operator, pair creation, the continuity equation, flux attachment, and the gauge redundancy. **NOT** a label rename; **NOT** an imposed rotor Gauss law.
2. **The dual sign, correctly formulated (decisive EM signature).** From **one** emergent Maxwell action, derive the conserved-source response
   `S_eff ~ ½ ∫ [ ρ (ε k²)⁻¹ ρ  −  j_T (μ/k²) j_T ]`,
   i.e. static like-charge **repulsion** (density–density, `1/r²` force) **and** transverse like-**current magnetostatic attraction** (`j_T`–`j_T`), **both from the one Maxwell field.** **State it correctly:** two *moving* like charges are **net repulsive**; only the **magnetic (transverse-current)** term is attractive. This requires the **quantum Coulomb-phase gapless photon** — classical/thermal spin-ice magnetostatics alone is **BANNED** as a proof. **Interpret all `S_eff` signs through the *static interaction energy* between physical sources (not bare action / metric-signature signs), to prevent sign-convention false passes.**
3. **The `±w` / throat embedding (THE crux able-to-fail).** Show the model's `±w` throat is a genuine **vertex defect** of the constrained lattice, with: **additive integer charge from the divergence sum** (not a bare `Z₂` assertion); **flux dressing**; a **Gauss-preserving hopping** operator for defect motion (so "current" = moving defect couples via the emergent `j·A` — **not** `j=qv` inserted as in `pathA_39`); net-neutrality on a closed lattice respected; and the electric-gauge-defect vs compact-`U(1)`-magnetic-monopole duality frame fixed. If additive/Gauss-consistent `±w` charge **requires importing a `U(1)` by hand → `FAIL_CHARGE_POSTULATED`.**

## Able-to-fail acceptance (MANDATORY)

- **Scalar negative control (makes the dual-sign test non-vacuous *in the script*):** implement an explicit scalar mediator `L_φ = ½(∂φ)² − g φ ρ` coupled to the scalar charge density `ρ`. It must give an **attractive** static density–density interaction and **NO transverse-current (`j_T`) / magnetic channel at all** — so it **cannot** reproduce the Maxwell dual-sign structure (repulsive static + attractive magnetic). That is `FAIL_SCALAR_SINGLE_SIGN`. An assert FIRES if the scalar case spuriously produces a *repulsive* static channel or *any* transverse-current channel. (Do NOT model the scalar as "both channels attractive" — a scalar has no magnetic sector; the discriminator is: EM = repulsive-static **+** attractive-magnetic; scalar = attractive-static **+** no-magnetic.)
- **Anti-circularity = UV criterion (not a Gauss rename):** the microscopic Hamiltonian must have **no global conserved matter `U(1)`**; the IR `U(1)` + integer charge are *derived* from the constraint's defect structure. Starting from gauge-charged matter / an imposed local `U(1)` / an un-derived rotor Gauss law → `FAIL_CHARGE_POSTULATED`. (Drop the old tautological "removing the identification breaks Gauss law" check — a relabel cannot break an operator identity.)
- **Corrected controls (replace the false "no-constraint → no gauge"):** finite-energy constraint defects are dynamical matter and *can* coexist with the Coulomb phase, so removing the constraint is NOT the right ablation. Instead: **ring-exchange removed (constraint retained) → NO propagating photon**; **defects *condensed* → screening/Higgs (photon massive)** — note it is *condensation*, not mere gap-closure/criticality, that Higgses the photon.
- **Deconfinement diagnostics = CITED** from the pinned model (photon `~1/L` scaling, flux-sector energy scaling, Wilson/static-potential). Small-ED may show finite-size **TRENDS only** — it **cannot** establish thermodynamic deconfinement or a finite phase region; do not claim it does.
- **Single-kernel assert:** no independent scalar channel mediating the static charge–charge force (the program's `pathA_38`+`pathA_39` two-kernel graveyard).
- **Structural firewall:** integer emergent charge; quantized flux; exactly 2 transverse photon modes; `1/r²` falloff; ≥2 able-to-fail structural ablations.
- No source-grep acceptance; enforce with computed runtime guards.

## Rigor / tooling

- **Dual engine** on the algebra (emergent-gauge mapping, Maxwell IR, the `S_eff` propagator signs, `1/r²`): Mathematica `.wl` verifies **algebra, not phase existence**; `ENGINE_AGREE` on the sign structure + falloff.
- Small-ED (optional) for constraint-sector / ring-exchange / defect-energy / finite-size photon **trends only**.
- Scripts under `software/em_charge_attribute/`; runners `timeout 600`; no script > 10 min.
- **Honest scope:** deconfinement **cited** (not re-proven); mapping + Maxwell IR + dual sign **analytic**; the **`±w` embedding is the genuine test**. Does **NOT** include gravity/superfluid coexistence (softness, deferred), the 4+1D brane embedding, or mode-count reconciliation with brane shear (§6 downstream). Requirement A (internal DOF) is a stated postulate.

## Output

`software/em_charge_attribute/reports/emergent_em_construction.md`: the pinned model; the emergent Gauss law + photon (from microscopic moves); the **charge = defect** derivation; the **corrected dual-sign** `S_eff` (both channels, one Maxwell field, magnetostatic-attraction correctly identified); the **`±w`/throat vertex-defect embedding** (divergence-sum additivity, flux dressing, Gauss-preserving hopping); the control results (scalar → `FAIL_SCALAR_SINGLE_SIGN`; ring-exchange-off → no photon; defect-condensed → Higgs; no second scalar kernel; deconfinement cited); the structural firewall + dual-engine logs; honest-scope caveats; and a one-line **VERDICT** — `EMERGENT_EM_WITH_DUAL_SIGN` / `FAIL_SCALAR_SINGLE_SIGN` / `FAIL_CHARGE_POSTULATED` / `CONSTRUCTION_INCONCLUSIVE`.
