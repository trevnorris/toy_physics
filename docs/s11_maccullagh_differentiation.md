# Why the light sector survives where MacCullagh's ether died

**Framing / write-up doc — NOT a physics authority or a spec.** This consolidates, for the eventual
prior-art / citation pass, *why* the 19th-century elastic ether died and *which specific S11 results*
address each cause of death. It is banked at the user's request (2026-09-10) while the S11c-d build is in
flight. The literature mapping and the "where each piece sits" table live in the auto-memory
`reference-prior-art-maccullagh` (richer, but 32 days stale on the S11 state) — read both together.

⛔ **Framing discipline (see `feedback-framing-split`, `project-analog-framework-goal`):** publicly this is a
**toy analog / mathematical bridge**, ⛔ not an ontology claim that "light really is brane shear." The
survival claims below are about *self-consistency and falsifiability of the analog*, not about the ether
being real. And per the standing rule, **NOT-FOUND ≠ original** — the novelty claim is provisional until the
citation pass actually searches.

---

## The core observation

The **entire S11 program is organized around the exact question that killed MacCullagh.** `V3_STEP_PLAN.md`
states S11c's headline as *"is light's confinement unconditional?"* — and confinement-vs-leakage into the
longitudinal sector is precisely the wound the elastic ether died from. We are not sidestepping the historic
failure; we are computing it head-on and letting the chips fall.

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

### 1. Confinement is *derived*, not postulated — and it is an identity. → answers cause (1) & (2)
MacCullagh *assumed* transverse-only. **S11b proved the transverse↔longitudinal coupling is identically zero
in the uniform brane** (`565b3fe8`) — light is confined as a *derived* consequence of one medium's dynamics,
not an ad-hoc energy. The `S_leak` object is an **identity**, not an ansatz: *"nobody set out to derive a
leak; it is what is left because the window has edges"* (`V3_STEP_PLAN.md`, S11b `S_leak`). That derived
substrate is exactly what cause (2) says he lacked.

### 2. The longitudinal mode is kept, computed, and shown controlled — not suppressed. → answers cause (1)
This is the direct differentiator. The thickness/breathing (longitudinal-ish) sector genuinely *exists* in
the model (S11c-c2's self-energy sector). Instead of assuming it away, S11 computes:
- **its own fate** — does it radiate or stay bound — in the homogeneous case (S11b-A/B, subsumed into the
  closed `S11b`); and
- **whether light's confinement survives a non-uniform brane** — the S11c program. The plan is explicit:
  *"the longitudinal mode's fate needs no gradients; light's confinement needs them."* S11c-d computes the
  **weak (Born) leakage law**; S11c-e the **strong slit-edge limit**.

Where MacCullagh made the compression mode vanish by assumption, we keep it and *earn* the confinement — or
get falsified. This is also where the stray longitudinal acquires a **physical role** (the deliberate anchor
for the charge sector) rather than being an embarrassment to hide.

### 3. A falsifiable leakage law + magnitude bound. → answers cause (3)
The S11c end deliverable is the **flux-normalized leakage FORM** plus a **magnitude bound** (the withheld
O(1) criterion, `N7`, which needs the throat interior `R1` — a later step). MacCullagh predicted nothing new;
this predicts *light leaks by this law, bounded by this size*, and if a generic brane feature leaked at
order-unity the model is dead. Novel, testable, killable — the thing his theory never was.

### 4. Single-substrate, multi-sector provenance. → answers cause (2) at the architecture level
Light here is *one sector* of a single medium that also yields gravity (flow) and charge/magnetism (throat),
not a standalone optics re-description. The nearest living programmes each have *half* of this:
- **Unzicker** (*ZAMM* 10.1002/zamm.202100280, 2022) — MacCullagh + topological defects ⇒ EM incl. charge,
  but **no bulk, no brane, no codimension** (an incompressible 3D solid). ⇒ the *confinement architecture is
  not his*; only the charge-from-defect move overlaps. ⚠ His sharper warning: linear-elastic MacCullagh
  *cannot* describe charge — needs a large (twist-disclination) deformation; our charge is a throat/puncture
  (also large), so this is an oracle to test, not a premise.
- **Volovik** (*The Universe in a Helium Droplet*) — one medium ⇒ gravity AND EM, but via Fermi-point
  topology / effective gauge fields, a different mechanism.

⇒ **the novel-looking residual, if any, is the ARCHITECTURE** (confinement package + stray-longitudinal-as-
charge-anchor + gravity-as-drain-flow-while-light-is-MacCullagh-on-the-same-sheet), ⛔ not any single
derivation.

## The bigger arc (banked, "for later")

Leakage is not purely a failure mode. The velocity leak *"lies outside the passive region ⇒ costs a named
reservoir, ⛔ not forbidden,"* and the user's **dark-energy postulate** (bulk reordering onto the brane → the
brane expands → cosmic expansion) is where a *controlled* leak becomes a feature. So the S11 confinement/
leakage result is also the hinge the cosmology story later hangs on.

## Honest status (as of 2026-09-10)

- **In hand:** uniform confinement = identically zero (S11b, `565b3fe8`); the longitudinal mode's homogeneous
  fate; the `S_leak` identity.
- **Being computed now:** the non-uniform *weak* leakage law (S11c-d — current build; the mixing/S-matrix/
  poles/survival/leakage-FORM engine).
- **Still ahead:** the strong slit-edge limit (S11c-e); the magnitude *bound* (needs the throat interior
  `R1`). So the differentiators are partly proven, partly the target we are building toward.

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
`steps/S11c_SCOPE.md`. Current build state: `research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`.
