# S11c — scope for the non-uniform transverse coupling (was "S11b-C")

⭐ **This is the SCOPING doc for S11c, the next front.** It consolidates the C requirements that were
preserved as historical input across `steps/S11b_HANDOFF.md`, the decision list `directives/S11b_unified_decisions.md`
(G14), and `V3_STEP_PLAN.md#s11b-split`. ⛔ **It is NOT the decision list.** S11c's first artifact is its OWN
decision list (orchestrator-written, two legs), which **re-validates** every requirement below — ⛔ do not
treat this package as ratified (the unified S11b decision list ratified only C's *scope*, not this specific
requirement set; `steps/S11b_HANDOFF.md:15-18`).

Naming (G1): **C → `S11c`** (a forward decision; no artifacts to migrate — C was never built). The closed
unified step is `S11b` (A+B). ⛔ Do not reuse the `S11b` slug for S11c artifacts.

## The question
**Is light's confinement unconditional?** S11b answered only the **uniform** case, where the transverse
coupling `∂²U/∂u_T ∂e_W` is **identically zero** (both engines; step record `steps/S11b_interface_coupling_law.md`).
⭐ **The gradient-driven channel is exactly what S11c computes.** The seam is real: the longitudinal mode's
fate needs no gradients (S11b did it); light's confinement needs them (S11c).

## The requirements — all named by earlier reviews, ⛔ none optional (re-validate each in the decision list)
1. **Tilted faces.** With `W₀(x)` varying, an in-plane displacement across `∇W₀` produces **normal** face
   motion of the **same gradient order** as the coupling being computed. ⛔ A frozen/flat wall (the S11b
   freeze) cannot see this — S11c must un-freeze `W₀(x)` in-plane.
2. **Eulerian vs material density.** The two differ by an advective `u·∇ρ_4D` term that ⭐ *directly changes*
   the transverse–thickness coupling. ⛔ Do not silently pick one; the choice is physics.
3. **Plane waves are NOT eigenmodes.** Treat the inhomogeneity as a **perturbation about a uniform state**
   with **mode conversion as an off-diagonal coupling**, or state exactly what restriction is granted.
   ⛔⛔ **Do NOT demand a global dispersion relation** — that error killed two S11b directive revisions.
4. **A uniform-limit control is KNOWN-VACUOUS.** B *proved* the uniform coupling is identically zero, so a
   "does it reduce to the uniform result" control tests nothing. ⇒ ⛔ that control is dead; **find another**
   (e.g. a form/gradient-order control, or a one-sided corruption of the gradient channel).
5. ⭐⭐ **Falsification FIRST — the coefficient is bounded by BENCH-TOP OPTICS.** If a slit edge converted an
   `O(1)` fraction of a photon into the thickness channel, diffraction gratings would not work. ⇒ S11c is
   **falsifiable against a lab measurement with no cosmology and no gravity sector** — ⭐ state that test
   **BEFORE** computing the coefficient (rule 5 posture: the bound is the acceptance, withheld from the
   builder; the diff happens on our side). ⇒ related: `project-double-slit-sim-idea` (a slit edge is a brane
   gradient ⇒ the coupling turns on there).

## Carry-ins (⛔ do not silently drop)
- ⛔⛔ **The background-flow correction `O(v₀ |q_n| / ω)` is uncarried and unbounded** (S11b known limit;
  fails first order where `|q v₀/ω| ≳ 1`). It exceeds first order where `k c_s0 ≫ ω` — the regime this works
  in. ⚠ `v₀` is the user's dark-energy leak. Decide in the S11c decision list whether S11c carries it or
  inherits it as a standing limit.
- The frozen-wall-width freeze (`ρ_br⁰ = rho_br`) is an S11b modeling choice; ⚠ requirement 1 (tilted faces)
  **un-freezes `W₀(x)` in-plane** — reconcile the two explicitly (in-plane variation vs the frozen slab
  thickness scale).

## Scope boundaries (G14 — do not let C swallow neighbours)
- ⭐ **In S11c** (non-uniform, linear): the full **variable-coefficient** slab spectrum, actual **leakage
  rates** in the non-uniform problem, and whether light's confinement is **unconditional**.
- ⛔ **NOT S11c — the nonlinear light program:** the DC/harmonic/sideband radiation audit and nonlinear
  intensity coupling. ⛔ Do not fold these into S11c.

## Method (the v3 pattern, and the S11b lesson)
- ⚠⚠ **SPLIT S11c FINER than S11b was.** S11b took **eleven** directive revisions before a build because the
  step surface was too large (`steps/S11b_HANDOFF.md:109`). ⭐ If the S11c surface is large, split it into
  tightly-scoped sub-steps in the decision list (⛔ decide the split WITH the two legs, do not pre-commit it
  here).
- The order is the proven v3 pattern (`.claude/skills/step-run/SKILL.md`, `steps/S11b_RUN_CHECKLIST.md`):
  **decision list (2 legs) → shared spec (2 legs) → SymPy engine + blind WL engine (each 2 legs) →
  frozen-T7 comparator → step record (2 legs) → card (2 legs)**. Whatever writes does not review.
- ⚠ **OPS (this session):** the spurious background-job killer hits codex/grok **bash** jobs generally — run
  review legs **serialized** and hand each a **mechanical digest** of only the tags/lines it needs so it
  finishes inside a sweep window; the **fresh-Agent** path stayed robust. Recover a kill mechanically
  (edits complete before the kill; review legs are idempotent — re-launch).

## Start here
Read this doc + `steps/S11b_HANDOFF.md` (C section) + `steps/S11b_interface_coupling_law.md` (the closed
uniform result S11c builds on) + `V3_STEP_PLAN.md#s11b-split`, then author the **S11c decision list**
(orchestrator-written ⇒ 2 legs Codex+Grok before any builder, rule 7 TRIGGER).
