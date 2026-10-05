# What the light sector establishes alongside MacCullagh — and what remains open

**Framing / write-up doc — NOT a physics authority or a spec.** This consolidates, for the eventual
prior-art / citation pass, the proposed historical comparison and the limits of the S11 results. Uniform confinement is derived
within the supplied linear model; nonuniform confinement and the material's physical admissibility remain
open. S11c is now PARTIAL; see [the closeout, MacCullagh and prohibited claims](../research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md). The literature mapping and the "where each piece sits" table live in the auto-memory
`reference-prior-art-maccullagh` (richer, but 32 days stale on the S11 state) — read both together.

⛔ **Framing discipline (see `feedback-framing-split`, `project-analog-framework-goal`):** publicly this is a
**toy analog / mathematical bridge**, ⛔ not an ontology claim that "light really is brane shear." The
survival claims below are about *self-consistency and falsifiability of the analog*, not about the ether
being real. And per the standing rule, **NOT-FOUND ≠ original** — the novelty claim is provisional until the
citation pass actually searches.

---

## The core observation

The S11 program asks whether a model with additional material motions can keep its light modes confined.
S11b answers this within the uniform linear model; S11c leaves the nonuniform answer unresolved. This is
one useful comparison with rotational-elastic light models, not a resolution of every historical objection.
[Source: S11c closeout, uniform result and MacCullagh boundary](../research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md).

The light sector's *algebra* is not new — it is MacCullagh's rotational aether (1839): curl-only stiffness
⇒ `D−1` transverse modes at `c² = μ/ρ` + a zero-restoring-force longitudinal slot (Whittaker, *A History of
the Theories of Aether and Electricity*, Ch. V). Standard ground is *solid* ground; what differentiates us
is the **architecture around** that algebra, not any individual derivation.

## Why MacCullagh's ether died (the four causes)

1. **The longitudinal-mode problem (the fatal one).** A real elastic solid carries both transverse and
   longitudinal (compression) waves; light is purely transverse. MacCullagh postulated a medium whose energy
   depends *only* on the curl of displacement, so there is **no restoring force for compression** → no
   longitudinal light. It reproduced optics but suppressed the longitudinal mode **by fiat**.
2. **No mechanical substrate.** That energy functional is not a real material (Stokes' objection): a
   *mathematical re-description* of optics with an ad-hoc constitutive law, not a medium anyone could point to.
3. **No novel predictions.** It re-encoded Fresnel's optics and stopped — nothing falsifiable beyond what
   optics already knew.
4. **Matter coupling / preferred frame.** How matter moves through the medium (drag, aberration) was never
   consistent; Michelson–Morley + relativity then removed the need for any mechanical medium at all.

(The memory's one-line version: *"19th-c. elastic aether died on longitudinal modes, preferred frames, and
matter–aether coupling."*)

## The four differentiators — mapped to the four causes

### 1. Uniform confinement is derived within a supplied material model
S11b derives zero transverse coupling on the **uniform** brane. SymPy and independent Wolfram work support
it, with the coefficient map and review limits stated in [S11b, transverse mode and comparison](../research/pde_ledger_v3/steps/S11b_interface_coupling_law.md).
This earns a conditional confinement result; it does not derive the rotational stiffness itself or establish
its mechanical admissibility. The original `S_leak` window identity is a lead, not an imported authority
([plan, S_leak](../research/pde_ledger_v3/V3_STEP_PLAN.md)). Reference/stress/torque closure remains open
([native interpretation, §§5.7, 14.3](native_light_em_and_vortex_throat_interpretation.md)).

### 2. Additional motion is retained; nonuniform conversion is unresolved
This is the direct differentiator. The thickness/breathing (longitudinal-ish) sector genuinely *exists* in
the model (S11c-c2's self-energy sector). Instead of assuming it away, S11 computes:
- **its own fate** — does it radiate or stay bound — in the homogeneous case (S11b-A/B, subsumed into the
  closed `S11b`); and
- **whether light's confinement survives a non-uniform brane** — the S11c program. The plan is explicit:
  *"the longitudinal mode's fate needs no gradients; light's confinement needs them."* S11c-d sought the
  **weak (Born) leakage law** but closes PARTIAL; S11c-e's **strong-edge interpretation** is deferred.
  [Current d record, handoff](../research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md).

The useful distinction is retaining compression, thickness and bulk channels and testing their coupling.
Uniform decoupling is earned within the supplied model; confinement around defects is not established.
Nor is the in-plane longitudinal mode the charge-carrying normal displacement. The charge connection is
an open question, not an identified holder. [Plan, S11 distinction and charge preamble](../research/pde_ledger_v3/V3_STEP_PLAN.md).

### 3. A testable leakage target, not an obtained law or magnitude bound
The intended deliverable was a flux-normalized conversion law and, with a physical throat profile, a
magnitude test. S11c has not delivered an accepted nonuniform loss number or a physical upper bound.
The omega=3 benchmark's empirical 1 ppm floor is a numerical limit, not an optics prediction.
[Canonical d record, benchmark and handoff](../research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md).
The possibility of such a test remains useful; it cannot be reported as a completed differentiator.

### 4. A proposed single-substrate architecture, with material closure still open
The intended architecture puts light, flow/gravity and throat/charge in one medium. That ambition does
not prove a complete compatible material, throat or electric-force mechanism; those obligations remain
in S12/Q2/Q3/S22 and the material audit. [Closeout, owners](../research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md). The nearest living programmes each have *half* of this:
- **Unzicker** (*ZAMM* 10.1002/zamm.202100280, 2022) — MacCullagh + topological defects ⇒ EM incl. charge,
  but **no bulk, no brane, no codimension** (an incompressible 3D solid). ⇒ the *confinement architecture is
  not his*; only the charge-from-defect move overlaps. ⚠ His sharper warning: linear-elastic MacCullagh
  *cannot* describe charge — needs a large (twist-disclination) deformation; our charge is a throat/puncture
  (also large), so this is an oracle to test, not a premise.
- **Volovik** (*The Universe in a Helium Droplet*) — one medium ⇒ gravity AND EM, but via Fermi-point
  topology / effective gauge fields, a different mechanism.

⇒ **the novel-looking residual, if any, is the ARCHITECTURE** (confinement package + stray-longitudinal-as-
charge-question + gravity-as-drain-flow-while-light-is-MacCullagh-on-the-same-sheet), ⛔ not any single
derivation.

## The bigger arc (banked, "for later")

Leakage is not purely a failure mode. The velocity leak *"lies outside the passive region ⇒ costs a named
reservoir, ⛔ not forbidden,"* and the user's **dark-energy postulate** (bulk reordering onto the brane → the
brane expands → cosmic expansion) is a later proposed role for order conversion. It is not a measured
light-leakage result or an energy-supply proof. [Plan, dark-energy postulate and S12](../research/pde_ledger_v3/V3_STEP_PLAN.md).

## Honest status (2026-10-04)

- **CONDITIONAL:** uniform decoupling within the stated linear model, supported by independent SymPy/Wolfram
  construction; stability requires nonnegative stiffness. [S11b, transverse mode/reviews](../research/pde_ledger_v3/steps/S11b_interface_coupling_law.md).
- **UNRESOLVED:** nonuniform confinement and leakage size. The later uniform speed-match check is useful
  but does not settle a defect. [d record, uniform and benchmark](../research/pde_ledger_v3/steps/S11c_d_profile_conditioned_scattering.md).
- **OPEN:** strong-edge interpretation, physical throat response and rotational material admissibility.
  These are later questions, not successes attributed to S11c. [Closeout, ownership table](../research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md).

## Where to cite (for the eventual pass)

MacCullagh 1839 & the longitudinal critique → **Whittaker Ch. V**; interface law → **structural acoustics**
(fluid-loaded plate / radiation impedance / added mass — a *free external cross-check* on S11b-A, better than
a review leg); passivity/Onsager → **odd elasticity** (Fruchart–Scheibner–Vitelli, *Annu. Rev. CMP* 14, 471,
2023); brane modes → **branons** (Cembranos–Dobado–Maroto, PRL 90, 241301, 2003); defects ⇒ long-range
fields → **Eshelby** (*Solid State Phys.* 3, 1956); living programmes → **Unzicker** (ZAMM 2022),
**Volovik**. Analog-gravity / brane-leakage template → Unruh 1981; Barceló–Liberati–Visser; Randall–Sundrum
graviton leakage. Full table + status in `reference-prior-art-maccullagh`.

## Related

Memory: `reference-prior-art-maccullagh`, `project-analog-framework-goal`, `feedback-framing-split`,
`project-puncture-deflection-charge-mechanism`, `project-s11b-interface-law-result`,
`project-native-em-mechanisms`. Plan: `research/pde_ledger_v3/V3_STEP_PLAN.md` (S11 block),
`steps/S11c_SCOPE.md`. Current scope and limits: `research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md`.
