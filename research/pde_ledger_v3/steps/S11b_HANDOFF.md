# S11b — where the work stands, 2026-08-04

> ⭐⭐⭐ **S11b IS CLOSED (`565b3fe8`, 2026-08-23); this 2026-08-04 handoff is retained ONLY as the C (now
> `S11c`) requirements source.** S11b was rebuilt as ONE unified export-chain step ("the linear brane–bulk
> interface coupling law"; decision list `ddd0ae4c`, spec `1a2395a3`, build directives `9bd2f184`): SymPy
> engine `864d6f41` + blind WL engine `ec89f9df` cross-checked (frozen T7 `17fe32c8` / re-run `fba6a34c`; X-1
> over-count corrected `53fcd98d`; WL repair `bd598ae7`), step record `8ddccb74`, card `565b3fe8` — all 13
> steps done. ⛔⛔ **The entire "where the work stands" narrative below (and the table) is HISTORICAL
> PRE-UNIFICATION state — the "immediate next is the blind Wolfram engine" and "S11b-C ▶ NEXT" are DEAD.**
> ⭐ **Live NEXT = `S11c`** (the non-uniform transverse coupling) — ⭐ start at **`steps/S11c_SCOPE.md`**; live
> state in `STATUS.md` + `steps/S11b_RUN_CHECKLIST.md`.
> ⭐ **C's material below (the four ingredients + the bench-top-optics falsification) is preserved as
> historical input, ⛔ to be revalidated in C's own decision list/spec** — the unified decision list
> ratifies only C's *scope* (the non-uniform variable-coefficient spectrum, leakage rates, and
> unconditional confinement), ⛔ not this specific requirement package.

⛔⛔ **The table below is HISTORICAL / PRE-UNIFICATION and SUPERSEDED** — its "✅ CLOSED / two blind engines /
cross-verified / card reviewed" and "S11b-C ▶ NEXT" are the OLD A/B/C execution framing. Live status: the
banner above (**S11b CLOSED; NEXT = S11c at `steps/S11c_SCOPE.md`**) + `STATUS.md` + `steps/S11b_RUN_CHECKLIST.md`.

| | state (⛔ HISTORICAL, superseded) |
|---|---|
| **S11b-A** — the bulk's response to moving faces | ✅ closed; ⚠ **its headline was later overturned** — see below |
| **S11b-B** — the homogeneous assembly | ✅ **CLOSED.** Two blind engines, repaired, cross-verified; card written and reviewed by two legs |
| **S11b-C** — the non-uniform transverse coupling | ▶ **NEXT. Not started.** |

⭐ **Read `steps/S11bB_interface_assembly.md`** — it is the closed result and supersedes A's headline.
⭐ The ledger card is `paper/steps/S11b_interface_coupling_law.tex` (one card for S11b; the A/B split was
execution history, ⛔ not ledger structure).

---

## ⭐⭐⭐ What S11b established, and what it killed

```
passive region: Λ_A⁰ ≥ 0  and  [ Λ_V⁰ = Λ_X⁰ = 0   or   (Λ_X⁰ = −Λ_V⁰ , τ_V = τ_X = 0) ]
reciprocity:    Λ_X(ω) = −Λ_V(ω)   ⚠ CONDITIONAL on microscopic time-reversibility
transverse:     coupling ≡ 0 on a UNIFORM background;  ρ_br⁰ω² = μ_R k² ;  Im ω = 0
stability:      K₀ = B_ρ⁽³⁾ − 2CW₀ + k_W W₀² > 0
```

⚠⚠ **S11b-A's headline is CONDITIONAL, ⛔ not dead — and an earlier version of this file said "DEAD."**
"The velocity coupling dissipates iff `ωτ ≠ 0`" lies **outside the passive region** ⇒ it is available to an
interface that **draws on a reservoir**, and this model supplies a candidate: the background drain `v₀`.
⛔ **The engines never emitted a prohibition.** They emitted a region, tagged *"not used to remove the
root,"* under a directive stating *"admissibility is a classification, not an acceptance gate."*

⛔⛔ **TWO REASONS THE PASSIVE REGION DOES NOT BOUND THIS MODEL:**
1. **Onsager–Casimir needs microscopic time-reversibility**, which this model ⛔ **postulates nowhere.**
   The medium is a *substructure*; ⛔ thermodynamic laws emergent from particles made **of** it are not
   inherited by it.
2. **Passivity and reciprocity describe fluctuations about an EQUILIBRIUM; our reference state carries
   `v₀ ≠ 0`** ⇒ a **driven steady state.** ⚠ `v₀` is the same neglect filed as the step's top known limit.

⇒ ⭐⭐ **THE STANDING RULE:** a non-passive coupling is admissible **only with a NAMED reservoir and a
STATED power budget.** ⛔ Unbounded growth fed by nothing is still a defect. ⭐ **The real gate is
observational, ⛔ not thermodynamic.**

⚠ **The growing root at `K₀ < 0` is NOT an energy-conservation violation** — the stored energy has no
minimum in that direction; the accounting closes exactly.

⭐ **The model gained five energy terms the specification did not know were allowed:** the engines built the
closed symmetry basis and found **ten** invariants where the spec carried five — `(∇·u)²`, `|∇θ|²`,
`∇θ·∇e_W`, `θ(∇·u)`, `e_W(∇·u)`. 24 registry rows generated in `reduction/_generated/`.

## ⭐⭐ WHAT S11b-C MUST HANDLE

**C's question: is light's confinement unconditional?** ⛔ B answered only the uniform case, where the
transverse coupling vanishes identically. ⭐ **The gradient-driven channel is exactly what C computes.**

Requirements, all named by earlier reviews and ⛔ none optional:
- **Tilted faces.** With `W₀(x)` varying, an in-plane displacement across `∇W₀` produces **normal** face
  motion of the same gradient order as the coupling being computed.
- **Eulerian vs material density.** The two differ by an advective `u·∇ρ_4D` term that ⭐ *directly changes*
  the transverse–thickness coupling.
- **Plane waves are not eigenmodes.** Treat the inhomogeneity as a perturbation about a uniform state with
  mode conversion as an off-diagonal coupling, or state what restriction is granted. ⛔ Do not demand a
  global dispersion relation — that error killed two directive revisions.
- ⚠ **A uniform-limit control is EMPTY unless the coupling can have a gradient-independent part.** B proved
  it cannot: the uniform coupling is identically zero. ⇒ ⛔ that control is now known-vacuous; find another.
- ⭐⭐ **C's coefficient is bounded by BENCH-TOP OPTICS.** If a slit edge converted an O(1) fraction of a
  photon into the thickness channel, diffraction gratings would not work. ⇒ ⭐ **C is falsifiable against a
  lab measurement with no cosmology and no gravity sector** — the strongest reason to run it, and the test
  to state **before** the coefficient is computed.

⛔⛔ **Carry into C:** the background-flow correction **`O(v₀|q_n|/ω)`** is still **uncarried and
unbounded**, and it exceeds first order where `k c_s0 ≫ ω` — the regime S11b works in. ⚠ `v₀` is the
dark-energy leak. ⚠ It was recorded by A, **lost** in B's first draft, and restored only after a review
caught it.

## ⭐ How S11b-B was actually run — the parts worth repeating

- ⭐⭐ **SUPPLY what is already verified; keep blindness only for genuinely open questions.** Every genuine
  physics confirmation came from a **review leg deriving from scratch**, ⛔ never from a blind builder.
- ⭐ **Write binding kinematic relations as EQUATIONS, never prose.** `∇·u` fell out **four times**, every
  time it was described in a sentence.
- ⛔⛔ **CUT 2026-08-12 — this bullet described quarantine, sandboxing and a tripwire, ⛔ all three of which
  `CLAUDE.md` rule 12 had already cut.** ⚠ The clean build logs it cited were evidence of nothing: they
  are equally clean when a payload is hand-typed. ⇒ ⭐ **The only blindness control is
  `research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md:17`** — the Wolfram engine
  imports nothing and re-derives; ⛔ nothing else may be built pretending to be one. ⭐ What catches the
  real defect is rule 14's **FORM ablation**.
- ⭐ **The registry generator** (`reduction/generate_rows.py`) emits a row **only where both engines agree**
  and blocks on disagreement. ⚠ It caught a defect five review rounds and a full cross-engine comparison
  had missed. ⇒ **use it for C.**
- ⚠ **Eleven directive revisions before a build.** The step surface was too large. ⭐ **Split C finer.**

## ▶ Next actions, in order

1. **S11b-C directive** — the non-uniform coupling, using the four requirements above.
   ⭐ **Split it smaller than B was.**
2. Build blind (`.wl` first, then `.py` sandboxed), reviews **on sight**, cross-engine comparison.
3. **Registry rows via the generator**, then live insertion — ⚠ `relations.yaml` is still **empty**, and
   S11b's admissibility and reciprocity relations belong in it as reviewed candidates.
4. ⚠ **Two S11b card fixes owed**, deferred because C may change the card: undefined symbols
   (`Λ_A⁰`, `τ_V`, `K₀`'s constituents) and a dropped `Λ_p⁰ = 0` qualifier.
5. ⭐ **A conceptual walkthrough is owed to the user** — the physical story, in plain language.
   ⇒ the section below is the material for it.

## ⭐⭐ THE PHYSICAL READING — material for the walkthrough

**The slab IS the brane.** Not a surface — a layer of finite thickness `W` in the `w` direction, with two
**faces** where it meets the bulk. Our 3D world is that layer.

**1 · The bulk pushes back two different ways, and frequency picks between them.** Fast enough that sound
escapes ⇒ **radiation resistance**, energy genuinely carried off. Too slow for its wavelength ⇒ the bulk
sloshes locally, carries nothing away, and just makes the face feel **heavier**. The divider is the bulk's
sound cone.

**2 · The slab can BREATHE, and that is the whole coupling.** When thickness changes, **both faces move
outward at once** and push on the bulk. A true surface would have no such mode and nothing to derive.

**3 · Mass conservation ties three things together.** Thicken the slab at fixed material and it dilutes:
`δθ + δe_W + ∇·δu = 0`. ⭐ Just "the material has to go somewhere" — and the term that kept vanishing.

**4 · You cannot get these equations from an energy principle.** The interface moves material **one way**;
an action is time-symmetric and would have it flow back. ⇒ balance laws. The symptom of getting it wrong
was a mode growing like `e^{t/τ}`.

**5 · Trapped vs leaky IS the sound-cone question.** Below the cone the disturbance cannot radiate and is
bound to the slab; above it, it leaks. ⭐ The branch ambiguity that cost three revisions was that physical
distinction wearing mathematical clothing.

**6 · What drives material across the interface?** Not bulk pressure alone — a **chemical-potential
difference**, `𝒜 = μ_s − δp/ρ_m`. Under the old closure a compressed brane sitting still next to a calm
bulk had nothing driving conversion, which cannot be right: a compressed brane *wants* to shed material.

**7 · And then thermodynamics deleted a channel.** The velocity-driven version of that conversion turns out
to be an energy **source** unless it is instantaneous — so the second law forbids it. ⭐⭐ **The model
refused an unphysical process on its own terms rather than us patching it out.**
