# Directive — the softness gate: does supersolid phonon/Goldstone dressing confine the compact-`U(1)`?

**Goal (one line):** Decide whether **superfluid/supersolid softness pushes a 3+1D compact-`U(1)` Coulomb phase across the confinement boundary** — the make-or-break gate for the "EM = decoupled charge-attribute" reframe (`docs/em_charge_attribute_requirements.md` §5.1).

**Revised after Codex directive design-review (`DIRECTIVE_NEEDS_FIXES`).** The decisive observable is **NOT** `δM/M` (that is only a small-fluctuation *validity* check). The decisive quantity is the **phonon/Goldstone-renormalized monopole fugacity / worldline free-energy (tension), followed relative to the KNOWN 3+1D compact-`U(1)` Coulomb↔confinement phase boundary.** This is a semi-analytic estimate that *uses* the established phase boundary; it does NOT re-derive the nonperturbative transition (that is the deferred lattice model, §5.2).

---

## Fixed physics (the rest you must derive from a *specified* microscopic model)

- 3+1D **compact `U(1)`** lattice gauge theory (Villain form), spacing `a`. Monopoles = topological defects; deconfined Coulomb phase for weak effective coupling / dilute gapped monopoles; confined at strong coupling (monopole condensation). The **rigid-lattice phase boundary is a known input** (Guth 1980; Banks–Myerson–Kogut 1977; lattice compact-QED).
- The lattice is a **supersolid**: crystalline elastic modes (transverse + longitudinal phonons, elastic tensor `c_{ijkl}`) **and** a superfluid Goldstone (stiffness `ρ_s`) that **mixes with longitudinal strain** through compressibility/entrainment (do NOT treat the Goldstone as an independent geometric fluctuation).
- **The coupling must be a specified magnetoelastic term**, not `g∝1/a` by fiat. From it, DERIVE: `M_mon(a)`, the log-sensitivity `α ≡ d ln M_mon / d ln a`, and the gauge coupling's spacing/strain dependence via microscopic matching.

## What you must compute (requirement, not route)

1. **Specify one concrete microscopic model** — a compact-`U(1)`/quantum-ice Hamiltonian + an explicit magnetoelastic coupling to the supersolid elastic + superfluid-Goldstone sector. State it fully.
2. **Integrate out the phonon + Goldstone modes** (Gaussian) → an **effective monopole action** → the **renormalized monopole fugacity `y_ren` (and/or worldline tension)**. Include the **strain spectral correlator** with correct transverse/longitudinal mode structure, Goldstone–longitudinal-strain mixing, density, elastic tensor, and a physical Brillouin-zone cutoff. Handle `T=0` (core action) and finite-`T` (free energy/density) **separately** — do not mix `e^{−M/T}` with the zero-`T` core action.
3. **Locate `y_ren` relative to the known compact-`U(1)` boundary** as a function of the medium's stiffnesses (elastic moduli, `ρ_s`) and the derived `α`. The gate: does dressing keep `y_ren` on the **deconfined (Coulomb)** side across the stable-supersolid range, or push it **confined**?
4. Report **`δM/M` only as a small-fluctuation validity flag**: if `δM/M` is not `≪1` in the regime of interest, the Gaussian estimate is **out of its domain** → `INDECISIVE`.

## Able-to-fail acceptance (MANDATORY)

- **Model-pinned, not exponent-free.** `M_mon(a)` and `α` must be **derived** from the specified model, not chosen. Elastic/material parameters and their uncertainty ranges may be **declared**. Output a **sensitivity surface** over (declared elastic params, derived `α`); a physical verdict is allowed **only if robust** across the justified range — otherwise `INDECISIVE`.
- **Decisive CONTROLS (reproduce established physics or the method is wrong):**
  - **Rigid 3+1D coupling sweep (the key new control):** with phonons off, sweeping the bare coupling must reproduce the **known compact-`U(1)` transition** — weak-coupling **Coulomb** and strong-coupling **confined**. An assert FIRES if the method fails to show *both* phases.
  - **Neutral-Goldstone-decoupled control:** if the Goldstone is decoupled from the gauge sector, `y_ren` is **unshifted** (no confinement from a decoupled neutral mode).
  - **Charged-condensate control:** making the condensate gauge-*charged* produces an **Anderson–Higgs photon mass** (distinct channel) — confirms neutrality is what protects the Coulomb phase.
  - **2+1D sanity:** the same machinery in 2+1D must **not** yield a robust deconfined Coulomb phase (Polyakov) — flag if it does. (Diagnostic, not the primary gate.)
- **Indecision is a valid, non-punitive outcome:** an order-one Gaussian correction, proximity to the boundary, or the need for dislocations/nonperturbative physics → `INDECISIVE`, **not** automatically `SOFTNESS_CONFINES`.
- **Dimensional firewall:** every physical expression units-restored; ≥2 able-to-fail dimensional ablations (a dropped `ℏ`, `a`, `k_B` must trip an assert).
- No source-grep acceptance; enforce real properties with computed runtime guards.

## Rigor / tooling

- **Dual engine** where algebra permits: independent Mathematica (`.wl`) re-derivation of `y_ren`/tension + the boundary comparison, `ENGINE_AGREE` over key quantities.
- Scripts under `software/em_charge_attribute/`; runners `timeout 600`; no script > 10 min.
- **Honest scope:** Gaussian integrate-out + *use of the known phase boundary*; NOT the full nonperturbative transition (deferred lattice model). State declared-vs-derived inputs and where the Gaussian approximation breaks.

## Output

`software/em_charge_attribute/reports/monopole_softness_estimate.md`: the specified model; the derived `M_mon(a)`, `α`, strain correlator; `y_ren` and its position vs the known boundary; the **sensitivity surface** + **critical `α`**; the control results (rigid sweep shows both phases; decoupled → no shift; charged → Higgs; 2+1D not-deconfined); the `δM/M` validity flag; the dimensional-firewall + dual-engine logs; honest-scope caveats; and a one-line **VERDICT** — `SOFTNESS_TOLERABLE_COULOMB_SURVIVES` (robustly deconfined across the stable range), `SOFTNESS_CONFINES` (robustly pushed confined), or `ESTIMATE_INDECISIVE` (not robust / Gaussian out of domain / a control failed).
