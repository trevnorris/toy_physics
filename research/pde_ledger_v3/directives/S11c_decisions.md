# Decision list — S11c ("the non-uniform transverse coupling")

**Orchestrator-written, 2026-08-23; folded ONCE after two legs (Codex + Grok), 2026-08-23.** Settles the
structure, naming, provenance and scope questions the **S11c** specs and build directives must honour.
⛔ **Not itself a spec**, and ⛔ not a build directive. Rule 7 TRIGGER: no builder launched before these two
legs; one pass, folded, and now used — ⛔ this list is not re-reviewed (its corrections are re-checked at
each sub-step's spec, which gets its own two legs).

⭐ IDs are the **N-series** ("N" = non-uniform). The export-chain rules **`F1`–`F9`**
(`S11_export_chain_decisions_v2.md`), the three script clauses (`.claude/skills/build/SKILL.md`), and the
closed **S11b** decisions (`directives/S11b_unified_decisions.md`, `G1`–`G14`) are **inherited, not
restated** — an N-item points at them.

⚠ **The two legs converged (independently) on seven load-bearing defects in the first draft** — the five
requirements the scope doc consolidated (`steps/S11c_SCOPE.md`) were **unratified historical input**
(`steps/S11b_HANDOFF.md:12-15`) and each needed correction, not transcription. The corrections are folded
below and marked ⭐**[folded]**. Both legs verified against source by the orchestrator (rule 13); the
load-bearing physics anchors are in the `_measurements` twin.

---

## A · Structure — S11c is a STAGED FAMILY, not one step ⭐[folded]

**N1 · S11c is a staged FAMILY of export-chain sub-steps, ⛔ not one monolithic step.**
⭐[folded — both legs: one spec + one engine pair over this surface cannot be independently reviewed, and
S11b's single spec took eleven directive revisions.] Each sub-step is its own build unit: **its own spec,
two blind engines, frozen comparator, and `S11c_<x>_exports.py`** — reviewed on its own. The SymPy engines
**chain** (S11c-a imports `scripts/S11b_exports.py`, 1958 rows, `steps/S11b_interface_coupling_law.md:22`;
each later sub-step imports the prior one's exports); the Wolfram engines are **blind** and re-derive. A
single roll-up ledger **card** closes the family. `BUILD_INPUT_DIGESTS` pins, per sub-step,
`{that sub-step's audit, the imported exports, that sub-step's spec}`. ⛔ **Do not reuse the `S11b` slug.**

**N2 · The split — decided WITH the legs; ⛔ not pre-committed, ⛔ not frozen. ⭐[folded]**
Both legs independently found the draft's three-seam sketch **too coarse** — its "result" seam re-merged
spectrum + leakage + falsification + conclusion (the largest surface, and the one that re-invites the
forbidden global spectrum, `N5`). The **adopted starting family** (Codex's five, into which Grok's
kinematics/matrix-element/photon three-way maps):

| sub-step | scope — named objects, ⛔ not recipes | exports |
|---|---|---|
| **S11c-a** background & geometry | which quantities vary (`W₀(x)`, `μ_R(x)`, `ρ_4D⁰`); the profile's **physical anchoring** (material-advected vs externally-held, and the stationarity/force that holds it, `N4`); the **(ε,η) power counting** (`N12`); the **FULL tilted-face shape-derivative** of every interface law — normal, relative flux, traction, bulk-field-at-shifted-face, DtN — ⛔ not only normal face motion (`N3`) | kinematics + boundary shape-derivatives |
| **S11c-b** variable-coefficient brane operator | the divergence-form slab operator + the **off-diagonal** transverse→`{θ,e_W,u_L}` kernel; inherit S11b invariants with variable coefficients and **emit** new gradient-of-background invariants as **results** (`N15`); the **representation-invariance control** (Eulerian ≡ material after field redefinition, `N4`/`N6`) | the coupling kernel |
| **S11c-c** curved-interface bulk closure | the perturbed two-face outgoing bulk **DtN/impedance** (flux, traction, permeability/memory, face parity); the `v_bulk_normal_0` carry-or-restrict decision (`N11`) | coupled nonlocal self-energy/operator |
| **S11c-d** profile-conditioned spectrum/scattering | for an **explicitly named profile class** (localized→Born/scattering kernel; periodic→Bloch; slowly-varying→WKB) — ⛔ **not** a generic "full spectrum"; the profile-conditioned **leakage rates** live here; ⛔ **no global dispersion relation** (`N5`) | scattering amplitudes / resonances / local spectrum |
| **S11c-e** leakage observable & falsification | the flux-normalized **dimensionless conversion FORM** for a supplied profile (⭐ magnitude needs the interior `R1` — out of scope, `V3_STEP_PLAN.md:1179`); **confinement** interpreted here (`N13`); the bench-optics bound **withheld**, diffed orchestrator-side (`N7`) | the conversion observable + leakage |

⚠ This is the **starting** structure; each sub-step's SPEC refines its own boundary. ⛔ Do not treat a–e as
frozen before their specs. ⚠ Sequencing/scale (how many run now) is the user's call (rule 11).

---

## B · The five requirements — re-validated & CORRECTED ⭐[folded]

**N3 · Tilted faces — the FULL first-order boundary shape-derivative, ⛔ not only normal face motion.**
⭐[folded — both legs.] With `W₀(x)` varying, an in-plane displacement across `∇W₀` gives normal face
motion `∼ u·∇W₀` at the coupling's own gradient order (correct as far as it goes). ⛔ But the same tilt
changes, **at the same order**, every interface law: the outward normal `n_s = (−½∇W₀, s)+O(|∇W₀|²)`, hence
the normal derivative `n_s·∇₄`, the **relative flux** (`S11b_SHARED_PHYSICS.md:195`, flux is measured along
the face's outward normal), the traction/pressure work, the evaluation of bulk fields at the **shifted**
face, and the perturbed **DtN/impedance** operator. Keeping only the normal-motion term can give the wrong
coefficient or the wrong **face parity**. ⇒ S11c-a demands the full level-set/shape derivative of every
interface law.

**N4 · Eulerian vs material is a REPRESENTATION, ⛔ not a free physics choice; the physics is the profile's
anchoring.** ⭐[folded — both legs, verified against `S11b_SHARED_PHYSICS.md:73,320-341`.] S11b already
fixed the representation: `θ` is **Eulerian** (:73) and the binding constraint is **material**
(`δ_v Σ_mat = 0`; `δ_vΣ_E = 0` is explicitly forbidden, :341). Eulerian and material variables are related
by the map + Jacobian (`Δρ = δρ_E + u·∇ρ⁰`) and **must agree after that field redefinition** — treating the
choice as "physics the spec picks" can double-count or omit `u·∇ρ⁰`. ⭐ The **real** physical decision is
whether the inhomogeneous background profile is **advected with the material** or **held fixed in lab/
Eulerian space** — S11c-a states it, and requires the two derivations to map to the **same** operator (the
representation-invariance control, `N6`). ⚠ The load-bearing advective term is `u·∇Σ_E⁰ = u·∇(ρ_4D⁰W₀)`,
tied to the freeze (`N11b`); the bare `u·∇ρ_4D` can drop while the coupling does not.

**N5 · Plane waves are not eigenmodes — and the OUTPUT must be a profile-conditioned object with an
order-count, ⛔ not a "full spectrum."** ⭐[folded — both legs.] The named trap (a global dispersion
relation) is avoided, but `N10`'s "full variable-coefficient slab spectrum" **re-invites it**: there is no
universal spectrum for an **unspecified** coefficient function. ⇒ S11c-d must **name a profile class** and
fix the output as the matching object (localized→Born/scattering kernel; periodic→Bloch; slowly-varying→
WKB), ⛔ never `ω(k)` for generic `W₀(x)`. ⚠ And it must fix the **expansion**: for inhomogeneity `O(η)`,
the off-diagonal coupling is `O(η)` and the leading leakage rate `O(η²)` — so an engine cannot discard the
leading coupling as "second order" nor confuse `O(η²)` leakage with the excluded nonlinear-intensity program
(`N12`). ⛔ "Perturbation + off-diagonal coupling" is a recipe (rule 3); the **object** is the linear mixing
between the uniform transverse and thickness sectors.

**N6 · The uniform-limit control stays as a REGRESSION; the real control is the independent shape/coordinate
route + one-sided corruption, named NOW.** ⭐[folded — both legs.] The uniform limit cannot validate the
gradient coefficient/sign/parity (S11b: coupling identically zero, `steps/S11b_interface_coupling_law.md:75`)
— but it is a useful **smoke test** for a forbidden gradient-**independent** term; keep it as secondary.
⛔ A "gradient-order control" phrased as "coupling `∝ ∇W₀` ⇒ vanishes at `∇W₀=0`" **is the vacuous uniform
limit renamed**, and `∇W₀→0` is ⛔ **not** an accepted corruption. The genuine control (rule 14): derive the
off-diagonal kernel by direct level-set/graph linearization, derive it **again** after flattening faces into
material coordinates, transform Eulerian↔material exactly and compare, then **corrupt one route only** —
flip one face's slope term, or omit its advective-density contribution — and require a nonzero residual.
⚠ There are ≥2 same-order channels (tilt, `N3`; advection, `N4`); the one-sided corruption is the
**independence test between them**, ⛔ so "the gradient channel" (singular) is wrong.

**N7 · Falsification: emit the dimensionless conversion FORM; ⛔ the `O(1)`/grating reductio is withheld
from the builder.** ⭐[folded — both legs; `O(1)` STRIPPED.] The lab bounds a **dimensionless conversion**
(photons lost at an edge), ⛔ not the bare, dimensionful `∂²U/∂u_T∂e_W` — a nonzero coupling still carries an
undetermined normalization until "what couples to what" is fixed (`S11b_SHARED_PHYSICS.md:824-825`). ⇒
S11c-e emits the **flux-normalized dimensionless conversion FORM** for a supplied profile (⭐ the FORM is
computable now; the **magnitude** needs the throat interior `R1`, out of scope — `V3_STEP_PLAN.md:1179`,
which also fixes the structure: `∝ k·a`, mixing driven by `∇μ_R ≠ 0`). ⛔⛔ The `O(1)`-fraction / diffraction-
grating argument is an order-of-magnitude **hint** (rule 5: the builder iterates to any target it can see) —
it belongs in the **step record**, ⛔ never in the builder-facing acceptance text; the withheld numeric bound
is diffed orchestrator-side. ⚠ A slit edge is an order-unity localized gradient, ⛔ not the small-`∇W₀`
regime of a Born kernel — the spec must state which, so we do not compare a weak-gradient coefficient to a
non-perturbative lab situation.

---

## C · Requirements the five OMITTED — added ⭐[folded]

**N12 · Two-parameter power counting, on every result.** ⭐[folded — Codex's sixth requirement, Grok's
order-bookkeeping.] Every S11c object carries an order label in **wave amplitude `ε`** and **background
inhomogeneity `η`**: the transverse↔thickness coupling is `O(εη)`; the leakage probability/rate `O(ε²η²)`.
⛔ No result may be reported without its `(ε,η)` order. Background **admissibility**: name which quantities
vary, whether each profile is material-advected or externally held, and the **stationary equations or the
named force** that holds the background (⛔ an inadmissible background silently sources spurious coupling).

**N13 · Define "confinement of light."** ⭐[folded — Grok.] It means **survival of the transverse
polarization channel**. Conversion into a **bound** breathing/thickness mode kills the photon exactly as
bulk radiation does — the two are **distinct emitted objects**. ⛔ "Energy stays in the slab" is **not**
confinement; a grating fails if the photon becomes a breathing excitation even with zero radiation.

**N14 · Name reservations against the imported keys (the F9 hazard, made concrete).** ⭐[folded — both legs;
`G3` was only listed as a trap to hunt, ⛔ not settled.] Every new **spatially-varying** field and kernel
gets a **fresh** standard name: the profiles `W₀(x)`, `μ_R(x)`, `ρ_br⁰(x)`, the off-diagonal kernel, the
leakage observable, and the drain **`v_bulk_normal_0`**. ⛔ Never reuse an imported S11b key —
`W_0`, `e_W`, `rho_br`, `v_0` (live at `scripts/S11b_exports.py:13806` (`v_0`, S10 PREMISE) and `:13813`
(`v_bulk_normal_0`)) — for a varying object: `F9`'s object comparison proves a **false equal** (`G3`), and
for `rho_br` it is a **physics bug** — it freezes `∇Σ_E⁰ = 0` and silently drops the advective channel.

**N15 · Constitutive construction — inherit S11b's invariants, EMIT new ones as results.** ⭐[folded — both
legs; verified `S11b_SHARED_PHYSICS.md:283`.] S11b's energy basis used **in-plane translation invariance**
(`u` enters only through gradients, :283); non-uniform `W₀(x)` **breaks** it. ⇒ inherit S11b's invariants
with **variable coefficients** and **emit any new gradient-of-background invariants as RESULTS** (new
constants if they appear); ⛔ neither smuggle them in by an `W₀→W₀(x)` substitution into the uniform `U`,
nor forbid them.

---

## D · Chain wiring — inherit, do not re-invent

**N8 · The chain mechanism is S11b's, inherited verbatim, PER sub-step.** Point at `F1`–`F9`; the blind-
Wolfram control (`S9_export_chain_rebuild_directive.md:16-18`); the frozen **`T7`** comparator contract
(`G8(a)`: join by name, residual paired payloads, reject a native boolean as a residual operand, three-
valued, repoint ablation); the `D3` in-run round-trip; the `_RELATIONALS` reviver; the three script clauses;
the **no terminal `VERDICT`** rule. ⛔ Do not restate — a re-worded obligation comes out weaker
(`feedback-point-at-the-obligation`). ⭐ Each sub-step gets its **own** comparator and exports (`N1`).

**N9 · The denylist stays CUT (rule 12).** Blindness is enforced by **absence** — the blind Wolfram engine
that imports nothing is the **only** cross-engine control. ⛔ No "what neither engine may read" list may
enter any S11c spec (`G9`).

---

## E · Scope boundary & carry-ins

**N10 · Scope boundary — inherit `G14`, but the in-scope object is PROFILE-CONDITIONED, ⛔ not a generic
global spectrum. ⭐[folded — both legs, see `N5`.]**
⭐ **In S11c** (non-uniform, linear): the variable-coefficient coupling kernel, **profile-conditioned**
leakage rates for a named profile class, and whether light's confinement is **unconditional**
(`directives/S11b_unified_decisions.md:134`).
⛔ **NOT S11c — the nonlinear light program:** the DC/harmonic/sideband radiation audit and nonlinear
intensity coupling (`:135`). ⇒ the family's card names each deferral's owner; ⛔ never imply a linear result
closed a nonlinear question, and ⛔ never let a mis-ordered `O(ε²η²)` leakage term be read as the nonlinear
program (`N12`).

**N11 · Carry-ins (⛔ do not silently drop; each was lost once already). ⭐[folded]**
- **(a)** The background-flow correction is `O(v_bulk_normal_0·|q_n|/ω)` — ⭐[folded: written with the
  reserved name, ⛔ never the glyph `v₀`/key `v_0`, `N14`]. **The decision is made here** (⛔ not deferred to
  the spec, both legs): S11c **inherits it as a standing rest-frame limit** (as S11b did), and **every**
  S11c spectrum/leakage/confinement result is **conditional on the derived smallness domain**
  (`|q v_bulk_normal_0/ω| ≪ 1`, `steps/S11b_interface_coupling_law.md:159-161` — large `k c_s0/|ω|` is
  necessary, ⛔ not sufficient). ⛔ Do not make the convective bulk operator an S11c task; the coupling is
  slab-side geometry. `v_bulk_normal_0` is the user's dark-energy leak.
- **(b)** The **frozen-wall-width freeze** (`ρ_br⁰ = rho_br`) is an S11b modeling choice; `N3`/`N4`
  **un-freeze the in-plane profile**. ⭐ S11c-a reconciles them explicitly — which of `ρ_4D⁰`, `W₀`,
  `ρ_br⁰` varies in-plane — ⛔ never import the freeze silently while un-freezing the same object (`N14`:
  a varying `ρ_br⁰(x)` must NOT reuse the constant `rho_br` key).

---

## F · What the legs decided, and what this fold did not change

⭐ **Both legs' verdicts on the five (folded above):** N3 incomplete→full shape-derivative; N4 needs-
correction→representation-invariance + profile anchoring; N5 incomplete→profile class + order-count; N6
needs-correction→keep-as-regression + named independent route; N7 wrong-as-written→dimensionless FORM,
`O(1)` stripped. Plus the staged family (`N1`/`N2`), and the four added requirements (`N12`–`N15`).
⚠ **Not changed, and carried to each sub-step's spec legs:** the exact retained-order at each interface law
(`N3`); the specific profile class for S11c-d (`N5` — a spec-stage choice, not pre-committed here); the
precise reconciliation of the freeze (`N11b`). ⛔ This list is folded once and is **not** re-reviewed; those
open choices are re-checked where they are made — in the sub-step specs, each with its own two legs.
