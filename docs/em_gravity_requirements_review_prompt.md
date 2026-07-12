# Critical review — requirements-inversion for gravity + light + EM in one medium (toy 4D superfluid analog)

You are a critical, adversarial physics reviewer with strong condensed-matter, gauge-theory, and analog-gravity knowledge. READ-ONLY. **Falsification is the goal — a proven no-go is a first-class result; do NOT rubber-stamp.**

## Context
A toy program: one compressible superfluid in 4+1D, an ordered "brane" (our 3D space) + disordered shear-free "bulk," "throats" (punctures, `±w` = charge sign) as particles. Prior rounds: forward-derivations of magnetism (`pathA_39`) and a geon/MacCullagh reconsideration were found (by a 3-AI panel, unanimous) to be a **computed no-go on the brane** — Maxwell's dual sign needs a *first-class* gauge (Gauss) constraint, but the brane's field content is *second-class*; forcing it first-class requires killing the compressional longitudinal mode, which **is** the gravity/density (`c_s`) channel. So EM and gravity appear to demand *opposite ontological status of the same mode*.

The author has now pivoted to **requirements-inversion**: assume the model is real, work backwards from GR + Maxwell to the constituent-level requirements, and ask whether one medium can satisfy them all — treating a genuine requirements-contradiction as the falsification.

## The document under review
`docs/em_gravity_requirements_inversion.md` — read it in full. Key claims:
- **§1** recasts the no-go as: gravity needs the longitudinal/compression mode PHYSICAL (`c_s` density wave); Maxwell needs it PURE GAUGE; on the brane it is the SAME mode → contradiction, *unless* they need not share it.
- **§2** per-force constituent requirements (compressibility→gravity; orientational "arrow" rotational rigidity→light; a genuine gauge `U(1)` with first-class Gauss + Ward identity→EM).
- **§3** two candidate arenas for EM to live in a DIFFERENT DOF than gravity's compression: (3a) the brane arrows; (3b) ⭐ the **throat interior / `w`-fall** — EM as emergent gauge from the phase/orientation dynamics of constituents crossing the order→disorder transition around the `±w` topological defect (charge is already `±w`; this is the one arena the brane no-go doesn't reach; analog: emergent electrodynamics from topological defects / order-disorder transitions).
- **§4** the decisive question: can one medium carry a physical-compression sector (gravity), a rotational sector (light), AND a cleanly-gauged EM `U(1)` — with the gauge and compression sectors DECOUPLED — or is that coupling fundamental (→ deep no-go)?

## Supporting context (read as useful)
- `docs/conceptual_foundation.md` §1 (constituents/substructure), §2 (brane = ordered phase), §3 (the sectors), §7 #12 (route-c / arrows).
- `software/stage1_solver/reports/pathA_36_c5_phase_potential.md` (the brane shear/light sector: 2 transverse photons + `FAIL_CAUCHY`; the `J·θ̇·δρ_B` density-phase coupling and wrong-sign `K_θ`), `pathA_38_throat_body_electric_localization.md` (charge/`h`-branon), `pathA_35_gateL_light.md` (`FAIL_COUPLE_STRESS_NOGO` on deriving the arrow sector).

## Assess critically (the §6 questions)
1. Is §1's requirements-contradiction (compression-physical vs longitudinal-gauge, same mode) **correct and complete**? Any missing requirement, or a hidden third option (e.g. a way for one mode to be *both* without contradiction)?
2. Are §2's per-force constituent requirements right — is a **first-class gauge `U(1)`** truly non-negotiable for Maxwell's dual sign, and what is the **minimal constituent structure** that supplies one?
3. **Which arena is more promising — 3a (brane arrows) or 3b (throat `w`-fall / emergent gauge from a topological defect at an order→disorder transition)?** Is there **real condensed-matter precedent** for a genuine emergent `U(1)` gauge field (first-class, not just a Berry phase) arising from an order→disorder transition or a topological defect, that could **decouple** from a coexisting compression/gravity sector? Cite mechanisms (quantum dimer/spin liquids, spin ice emergent magnetostatics, `BF`/loop models, `θ`-vacua, Volovik's `³He` emergent gauge fields, etc.) and say honestly whether they'd survive coexisting with a physical density mode.
4. §4 decisive question: is the **decoupling possible or a genuine no-go**? Give your best physical argument each way, and your lean.
5. The **sharpest, cheapest decisive test** to settle §4 — analytic/linear where possible; honest about where it needs the nonlinear throat.

## Output
Overall assessment of the requirements-inversion framing (SOUND METHOD / FLAWED), your answer to each of #1–#5 (especially #3 — the condensed-matter precedent question, and #4 — your lean on possible-vs-no-go), and the SINGLE most important thing to compute or decide next. Be specific; use real condensed-matter / gauge-theory knowledge.
