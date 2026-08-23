# Independent physics review — the S11c-a spec (background & geometry)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`

## What it is, and why it is load-bearing
S11c is the non-uniform transverse coupling at a brane–bulk interface — *is light's confinement
unconditional?* The prior step **S11b** (uniform) is CLOSED and proved the uniform transverse coupling is
**identically zero**; S11c computes the gradient-driven channel. S11c was split (with two prior legs) into a
five-sub-step family; **S11c-a** is the first: it **un-freezes the in-plane background** `W₀→W₀(x)` and
computes the **first-order geometry** — the tilted-face kinematics and the **shape-derivative of every
interface law** S11b wrote for a flat wall. This spec is orchestrator-written and is the physics authority
two blind engines (SymPy + Wolfram) will read; an error here corrupts both. ⛔ This is a **physics review**,
not a copy-edit.

## Derive your own view FIRST, then read the spec
Before opening the spec, form your own answer to: *what must a correct "background & geometry" sub-step set
up, given S11b's flat-wall interface laws and the folded S11c decision list?* Read:
- `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` — the closed uniform spec S11c-a builds on
  (esp. §1 geometry/fields `:63-160`, §3 relative flux `:183-207`, §4 face closure/affinity `:208-254`,
  §6 balance laws + virtual constraint + sign conventions `:301-380`, B0 bulk face response `:593-637`).
- `research/pde_ledger_v3/directives/S11c_decisions.md` — the folded decision list; S11c-a must honor
  **N3** (the FULL tilted-face shape-derivative, ⛔ not only normal face motion), **N4** (Eulerian/material
  is a representation S11b already fixed; the physics is the profile's *anchoring*), **N11b** (the
  frozen-wall reconciliation), **N12** (two-parameter ε,η power counting on every object), **N14** (F9 name
  reservations), **N6** (the η→0 limit is a known-vacuous regression, not the main control).
- `research/pde_ledger_v3/steps/S11b_interface_coupling_law.md` — the closed result (the uniform coupling is
  identically zero; the transverse mode).
- `research/pde_ledger_v3/directives/_measurements/S11c_a_SHARED_PHYSICS.md` — the rule-2 twin; spot-check
  that §1's claims about the imported S11b setup resolve to what the spec says.
There is **no do-not-read list** — read anything that helps you judge the physics.

## Required method — DOCUMENT branch
Read the sources first, form your own view, then read the spec. Quote **both sides** (source and spec line)
for every finding. ⚠ Reading order is a method instruction, not a blindness control.

## What to check — hunt each of these
1. **Faithful supply of S11b (§1).** Are the imported objects quoted correctly — the relative flux
   `J_± = ρ_m(v_w−∂_tζ_±)(±1)`, the closure `J_± = Λ_A𝒜_± + Λ_V V_±`, the affinity `𝒜_± = μ_s − δp_±/ρ_m`,
   the material virtual constraint `δ_vθ+δ_ve_W+∇·δ_vu=0` (and that `δ_vΣ_E=0` is forbidden), the sign
   conventions? Any misquote silently corrupts both engines.
2. **The full shape-derivative (N3), or a hidden truncation.** §3 lists T-a…T-g (normal, face velocity,
   relative flux, traction, face-shift, projection/window, constraint pullback). Is any interface law whose
   shape-derivative changes at `O(η)` or `O(εη)` **missing**? Conversely, does any item pre-state its result
   form (violating rule 3 / "name the object, compute the form")?
3. **Representation vs physics (N4).** Is the Eulerian/material issue correctly framed as a *representation*
   (the two must agree after `Δρ = δρ_E + u·∇ρ⁰`), with the genuine physical choice being profile
   **anchoring** (material-advected vs externally-held)? Is the representation-invariance control (§4) a
   real two-route + one-sided-corruption check, or vacuous?
4. **The freeze reconciliation (N11b, 2b).** Does the spec force the engines to *state* which density is
   constant and *compute* `∇Σ_E⁰`, rather than silently importing the freeze? Is the `rho_br`-reuse hazard
   (a varying `ρ_br⁰(x)` reusing the constant key) correctly forbidden?
5. **Admissibility (2d).** Is "the background is a static solution, or name the force holding it" the right
   obligation, and is the static-balance-residual the right computed object? An inadmissible background
   silently sources spurious coupling — is that risk closed?
6. **Power counting (N12, 2e).** Is the two-parameter `(ε,η)` bookkeeping actually binding on *every*
   emitted object, and correct (coupling `O(εη)`, and S11c-a's exports `O(η)`/`O(εη)`)? Could an engine
   discard an `O(εη)` term as "second order"?
7. **Scope boundary.** Does S11c-a correctly **defer** the coupling kernel (S11c-b), the bulk DtN closure
   (S11c-c), the spectrum/scattering (S11c-d), and the leakage/falsification (S11c-e) — or does it swallow a
   neighbour? Does anything here presuppose the FORM of a later answer (e.g. a global dispersion relation)?
8. **Under-specification (the recurring failure).** A binding kinematic relation stated in prose instead of
   an equation has fallen out of this ledger's specs repeatedly (`∇·u` four times). Is any load-bearing
   relation here prose where it must be an equation? Is any term at risk of silently vanishing?
9. **Controls (§8) and no-VERDICT.** Are C-1 (form) and C-2 (independence) genuine FORM/independence checks
   (rule 14), able to fail — not coefficient rescales or `η→0`? Is the "no terminal VERDICT / print-then-
   guard" discipline intact?

## Physics filter
Report a finding only if it catches a way the physics or the sub-step's set-up could be wrong — not
stylistic preference, not "wrong on a different input." ⛔ A leg returning "nothing survives" is weak
evidence: state what you checked, what you derived, and where you looked hardest.

## Output
1. Findings, most-serious first — each with the source quote, the spec line, and the concrete way the
   physics/set-up goes wrong if uncorrected.
2. Is any interface-law shape-derivative **missing** from T-a…T-g?
3. Is anything **under-specified** (prose where an equation is needed) or **over-specified** (a result form
   pre-stated that the engines should compute)?
4. Does the scope correctly defer S11c-b…e?
