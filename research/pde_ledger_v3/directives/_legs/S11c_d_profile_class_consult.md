# S11c-d profile-class design consult (independent physics opinion — NOT a build, NOT a leg)

You are advising on a **design decision** in a toy superfluid-analog physics program (the program builds a
math bridge for a GR+EM analog; it is not an ontology claim). I want your **independent physics opinion**. I
will synthesize your view with the other engine's and my own and bring all three back to the user. ⛔ Do not
just agree — if the choice I lean toward is wrong, say so and why.

## Program context (what the object is)
We are deriving the mode structure of a thin, variable-coefficient brane "slab" whose background thickness
field `W₀(x)` is **inhomogeneous** (varies along the in-plane coordinate `x`). The slab carries two coupled
linear sectors: a **uniform transverse wave** sector and a **thickness / self-energy** sector. A varying
background chemical-potential-like coefficient `μ_R(x)` (its gradient `∇μ_R ∝ ∇W₀`) is what **couples** the
two sectors — with a uniform background they decouple identically. The object of interest is the **linear
mixing** between the transverse and thickness sectors driven by that gradient.

## Already established and available (frozen S11c-b / S11c-c2 exports — consume as-is)
- `CLOSED_SLAB_OPERATOR` — the **diagonal** variable-coefficient (divergence-form) slab operator ⇒ for a
  given profile this gives the **local** spectrum / resolvent poles.
- `CLOSED_COUPLING_KERNEL` — the **off-diagonal** transverse↔thickness kernel ⇒ the mixing / Born amplitude,
  driven by `∇μ_R ≠ 0` (structurally `∝ k·a`).
- **Power counting (fixed):** two small parameters — wave amplitude `ε` and background inhomogeneity `η`.
  The transverse↔thickness coupling is `O(εη)`; the leakage rate `O(ε²η²)`. Every result must carry its
  `(ε,η)` order.

## The decision to advise on
S11c-d must **NAME A PROFILE CLASS** for `W₀(x)`. There is no universal spectrum for an *unspecified*
coefficient function, and a generic dispersion relation `ω(k)` for generic `W₀(x)` is **forbidden**. The
canonical named options and their matching output object:
- **Localized** profile → a **Born / scattering-and-leakage kernel**.
- **Periodic** profile → **Bloch** bands.
- **Slowly-varying** profile → **WKB** adiabatic mode-following.

And an orthogonal **regime** axis: **weak-gradient** (perturbative in `η`) vs **order-unity** (non-perturbative)
gradient.

## The endgame (why the choice matters — downstream falsification)
The falsification observable (a later sub-step) is a flux-normalized **dimensionless conversion FORM** —
photons lost at a **slit edge** in a bench-optics experiment. A slit edge is an **order-unity localized
gradient**, NOT the small-`∇W₀` regime of a Born kernel — this tension is explicit. Only the FORM is
computable now; the magnitude needs throat-interior physics that is out of scope.

## Constraints any choice must respect
- ⛔ No global dispersion `ω(k)` for a generic `W₀(x)`.
- The `(ε,η)` **polynomial** order-count (coupling `O(εη)`, leakage `O(ε²η²)`) must survive — a choice that
  sends `η → O(1)` breaks the bookkeeping.
- The output must be the **linear mixing** between the two sectors (an object, not a recipe).
- The chosen class must **consume the two frozen exports as-is** (no new upstream derivation).

## What I'm asking you
1. **Which profile class + regime should S11c-d fix, and WHY?** Weigh: internal consistency with the `(ε,η)`
   power counting and the frozen exports; alignment with the slit-edge falsification endgame; tractability now.
2. **Address the weak-gradient-Born vs order-unity-slit-edge tension directly.** If we pick the perturbative
   Born route for tractability, how should the weak→order-unity bridge be handled downstream **without
   invalidating** the eventual comparison to the lab bound? If instead you'd commit to the order-unity regime
   now, how do you keep the `(ε,η)` bookkeeping honest?
3. **Is a staged / hybrid answer right** — e.g. a Born scattering-and-leakage kernel as the S11c-d
   profile-conditioned object, with the order-unity slit edge treated as a *separate* downstream reduction?
4. **Have I mis-framed any option, or is there a better profile class I omitted?**

Give a **clear recommendation** with physics reasoning. You may read repo files under
`/var/projects/toy_physics/research/pde_ledger_v3/` for context (`directives/S11c_decisions.md`,
`directives/S11c_c2_SHARED_PHYSICS.md`, `steps/S11c_c2_self_energy_fold.md`) — but ground your recommendation
in the **physics**, not in what those docs assert.

## My current leaning (tell me if it is wrong)
Localized profile, **weak-gradient Born** scattering/leakage kernel — it keeps `(ε,η)` exact, consumes both
frozen exports directly, and defers the order-unity slit-edge bridge to the downstream falsification step
(where the order-of-magnitude reductio is withheld from the builder and diffed on our side). Push back if
that is the wrong call, or if it quietly throws away something the order-unity regime would capture.
