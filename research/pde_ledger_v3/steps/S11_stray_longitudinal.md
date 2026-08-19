# S11 · The stray longitudinal — closed homogeneous census, open interface physics

**Sector 1 (light), step 3.** Walked side by side with the user, 2026-08-02.

---

## What it is

S10 left one direction of `u` with `ω² = 0`. S11 asks what happens to that direction once the brane is
allowed to resist compression, what that costs the ledger, and what the resulting mode's physical status
actually is.

## What it does

Lifts S10's zero to a propagating longitudinal mode, enters the brane's compression modulus as a
**postulated** knob with a **named retirement condition**, closes the homogeneous three-finite-mode
census, and derives the simple \(k_w=0\) grazing locus. The free mode's actual leakage remains an
interface and spectral-boundary problem; its possible finite-amplitude photon-helper role additionally
requires nonlinear coupling and a harmonic/sideband radiation audit.

⛔ **This is not a defect being repaired.** The extra mode is a characterized departure from a
transverse-only Maxwell field. Historical drum-head discussions used the compressional response as
motivation for a material electric sector, but the current canonical ontology does **not** identify
\(u_L\) with charge: charge is the oriented throat/electric-odd boundary condition, while \(u_L\) is a
separate in-plane longitudinal material mode. ⛔ `FAIL_CAUCHY_STRAY_LONGITUDINAL` is a **misnamed** token;
never read the prefix as a verdict.

---

## The walk, move by move

**Move 1 — is there a compression channel at all?** ⭐ **The identification, flagged before use and
user-confirmed:** `u` is the **material displacement of the stuff whose density is `ρ_br`**. That is
forced by S9's own kinetic term — `½ρ_br(∂_t u)²` is only kinetic energy if `u̇` is the velocity of
material carrying inertia `ρ_br`. ⇒ Linearised continuity applies, `δρ_br = −ρ_br(∇·u)`, and since the
medium has an equation of state, **compression costs energy**. ⛔ Not a modelling choice.
⚠ **What would have made it wrong:** if `u` were an orientation/director field, `∇·u` would be a splay,
not a density change, and the channel would need separate argument.

**Move 2 — what a quadratic stiffness on the brane can be.** Decompose `∂_i u_j` into trace,
symmetric-traceless and antisymmetric parts. ⭐ **S10 did not forget a term — it kept ONE invariant of
several.** Curl-only stiffness charges for twist and nothing else, which is why the longitudinal came out
free rather than forbidden: a longitudinal wave is pure trace.

⛔⛔ **CORRECTION — the orchestrator's count in this move was WRONG, and both engines caught it.** The
walk asserted *"exactly three coefficients, for every `D ≥ 2`."* ⭐ **False under proper rotations.**

```
N_SO = {D=2 → 4, D=3 → 3, D=4 → 4, D=5 → 3}        N_O = 3 for every D
```

The extras are **reflection-odd**: `(tr G)·ε^{ij}G_{ij}` at `D=2`, `ε_{ijkl}G_{ij}G_{kl}` at `D=4`.
⚠ **And the walk's METHOD was the error, not its arithmetic** — *"decompose into pieces that do not mix
and count them"* misses **cross-pairings between isomorphic summands**, which is exactly what `D=2` and
`D=4` have (at `D=2` the trace and `Λ²` are both `SO(2)` scalars; at `D=4`, `Λ²` splits into self-dual ⊕
anti-self-dual). ⭐ Five independent derivations agree on the corrected counts.

⭐⭐ **The physical content of the extras, and it sharpens the step:** at `D=2` the extra invariant is
**NOT a total derivative** — its Euler–Lagrange operator is non-zero — so an `SO(2)`-invariant Lagrangian
**can mix the longitudinal and transverse sectors**. At `D=4` it **is** a total derivative and adds no
operator structure. ⇒ ⭐ **Sector separation is a `D=3` fact, ⛔ not a structural one.**

**Move 3 — the spectrum.** The trace term enters strictly as `B_comp k(k·a)`:

```
ρ_br ω² a = μ_R[k²a − k(k·a)] + B_comp k(k·a)
```

⇒ transverse (`k·a = 0`): the new term vanishes **identically** in the selected homogeneous \(D=3\)
quadratic action. Longitudinal (`a ∥ k`): the `μ_R` bracket vanishes and the zero is not perturbed but
**replaced**. Thus the homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero
linear cross-block under the selected action. This does not establish decoupling at nonlinear order, on
a nonuniform slab, at an interface or defect, or after additional allowed fields are introduced.

**Move 4 — where `B_comp` comes from.** ⛔ **The orchestrator's first mechanism was wrong** — it had
brane material "squeezing out into the bulk," which smuggles in a phase transition (bulk *is* the
disordered phase) and the user rejected it: *"How can the brane become unordered just by compression?"*
⭐ **The correct mechanism: the brane THICKENS.** In-plane compression raises the areal density, which
can be absorbed as a higher 4D density (charged by the EOS) **or as a wider wall** (charged by the wall
potential). ⛔ **Neither disorders anything** — the order parameter never changes; the sheet bulges into
`±w`. ⇒ Two channels accommodating additive shares of one strain, i.e. **springs in series**, so
compliances add and `B_comp` is **softer than either channel alone**.

**Move 5 — the flat direction.** `S7` records that the double-well **selects no slab width**. ⇒ If that
flatness survived, `B_wall = 0`, hence `B_comp = 0`, and S10's zero would never lift. ⭐ **It is lifted by
GRADIENTS**: a wave modulates the width, tilting the interfaces and stretching them at cost
`∝ σ_wall|∇W|²`. ⇒ flat at `k=0`, stiff as `k²`, predicting a **flexural crossover** — `ω ∝ k²` at long
wavelength. ⚠ **S11's cone below is computed with the wall width FROZEN**; unfreezing it should soften
the mode. ⛔ If it does not, **move 5 was wrong** — say so rather than reconciling.
⭐ `σ_wall` is **S6-derived**, so this predicts `B_comp` may cost no knob at all → `DEFECT_REGISTER#c12`.

**Move 6 — the bulk.** ⛔⛔ **THE ORCHESTRATOR OVERSTATED THIS MOVE TWICE, and both corrections stuck.**

Phase matching gives `k_w² = k²(c_L²/c_s0² − 1)`, and the walk read that as *"`c_L > c_s0` ⇒ it
radiates; `c_L < c_s0` ⇒ bound."* ⭐ **Neither direction is established.** Phase matching is **kinematic
only** — it determines whether a propagating bulk channel **exists**, ⛔ not whether the mode uses it,
nor whether an evanescent solution forms a **bound eigenmode**. Both require an interface coupling law
that is absent. ⚠ A third case was also omitted: `k_w = 0`, grazing.

---

## ⭐⭐ The computed result — two engines, written independently

| | value |
|---|---|
| dynamical matrix | `M_ij = μ_R(k·k)δ_ij + (B_comp − μ_R)k_i k_j` |
| transverse | `ω² = (μ_R/ρ_br)k²`, nullity `D−1`, kernel ⊥ `k` ⭐ **unchanged from S9/S10** |
| longitudinal | `ω² = (B_comp/ρ_br)k²`, nullity `1`, kernel ∥ `k` |
| cross-sector | ⊥ root independent of `B_comp`; ∥ root independent of `μ_R` — **computed residuals, both zero** |
| finite \(D=3\) census | characteristic polynomial degree 3 in `ω²`, leading coefficient `−ρ_br³ ≠ 0`; exactly three finite roots counted with multiplicity and no hidden fourth mode |
| degeneracy | **exactly** `B_comp = μ_R` |
| simple longitudinal/bulk threshold | `k_w² = k²(c_L²/c_s0² − 1)`; `KW_ZERO_LOCUS` is `B_comp = ρ_br c_s0²`, equivalently `c_L = c_s0`; kinematic only |
| dimensions | `[ρ_br]=(−D,0,1)`, `[μ_R]=[B_comp]=(2−D,−2,1)`, both speeds `(1,−1,0)` |

⭐ Every value matches predictions committed **before either script existed**
(`steps/S11_PREREGISTERED_PREDICTION.md`, `67d919bd`) — except the move-2 invariant count, where the
pre-registration was **wrong and the engines were right**.

**Controls.** ⭐ **FORM** — replacing the trace invariant with the symmetric-traceless one moves **both**
roots (`ω²_⊥ = (μ_R+μ_br)k²/ρ_br`, `ω²_∥ = 4μ_br k²/(3ρ_br)`) and the transverse root **acquires**
`μ_br`: ⇒ the trace invariant is the **unique** one that lifts the longitudinal while leaving light
untouched. ⛔ **COEFFICIENT** — rescaling `B_comp` moves only the parallel root and cannot test that
shape claim. ⚠ The form control's `4/3` independently reproduces the corpus's rejected-Cauchy-branch
coefficient, from an engine blind to that document.

## What's new

| item | class | why |
|---|---|---|
| `Q.brane.B_comp` | ⭐ **postulated**, with a **named retirement condition** | user's call: postulate now, retire visibly at S6. Knob count is an **upper bound that can only improve**. ⇒ `DEFECT_REGISTER#c12` |
| `Q.brane.c_L` | **derived** | `R5`, from `B_comp` and `ρ_br` |
| the mode itself | **derived** | nullity of the dynamical matrix |

⚠ **Naming, and it is a live hazard:** `B_comp` is a **brane** compression modulus. ⛔ It is **not**
`K_br` (a rejected elastic branch's bulk modulus), **not** `B_eff` (`= ρ_B0²/χ_c`), and not the bulk
medium's modulus. Registered with `aliases: []` and no borrowed loci.

**Registry:** ambient `10 → 12`, residue `6 → 7` (`{ħ, m, K, ρ0, ρ_br, μ_R, B_comp}`). Discrete payload
`{n_eos=5, D_brane=3}` unchanged and off the continuous axis.

---

## ⭐⭐⭐ Departure — closed homogeneous census, open interface and nonlinear roles

The audit now separates four questions that earlier prose compressed into one:

| question | present owner |
|---|---|
| where the free homogeneous longitudinal root lies relative to simple bulk sound | S11 kinematic `KW_ZERO_LOCUS` |
| whether that mode is bound, a threshold state, a bound state in the continuum, or a resonance | S11b interface and full-slab spectrum |
| whether matter directly observes the second cone | derived matter/interface coupling |
| whether the mode participates in a finite-amplitude photon helper | nonlinear intensity coupling plus harmonic and sideband response |

At homogeneous linear order, the mode's bulk leakage and direct observability require the brane–bulk
interface law. Its possible finite-amplitude photon-helper role additionally requires nonlinear intensity
coupling, the complete slab spectrum, and a harmonic and sideband radiation audit.

⛔ **What S11 does NOT deliver**, stated so it cannot be over-read: bound-versus-resonant classification ·
nonzero interface overlap · an actual leakage rate · unconditional light confinement · observability or
unobservability of the longitudinal branch · a localized or nonradiative nonlinear photon helper.

## ⭐ What this settles about the light sector as a whole

Established with the user, 2026-08-02, and load-bearing for how the departure above reads:

- ⛔⛔ **`c_L` is NOT a gravitational wave.** Gravity in this model is *"the **FLOW** between draining
  defects — carried by the flow + Bernoulli pressure, **NOT** by ripples/radiation"*
  (`docs/conceptual_foundation.md:348`). `c_L` is an in-plane displacement of brane material.
  ⚠ Separately: the corpus contains **no mechanistic account of what a gravitational wave is**.
- ⭐⭐ **A second cone is only a departure if matter COUPLES to it.** The transverse wave equation is
  Lorentz-invariant with invariant speed `c_γ` automatically, so anything built from those modes
  inherits that invariance — and this model's matter **is** built from them (a throat held open by a
  trapped brane-shear standing wave). ⛔ The orchestrator's *"a second cone means observers see no single
  invariant speed"* was wrong as stated. ⭐ **Not circular:** the wave equation supplies the invariance,
  so a standing wave built from it must restructure under boost, with the cost diverging as `v → c_γ`.
- ⭐ **Uniform flow is unobservable; GRADIENTS are observable** — and the model needs exactly that split.
  If uniform flow were observable we would detect absolute motion; if gradients were not, there would be
  no gravity. ⚠ *"The bulk flow is invisible because we are made of it"* and *"density gradients bend
  light"* are the same statement about different derivatives.
- ⭐⭐ **Light must carry stress-energy and must not dissipate unacceptably, but bulk shear-freeness alone
  proves only a narrower statement.** Absence of a transverse bulk branch removes the simplest direct
  transverse-to-transverse bulk channel. It does not exclude conversion into longitudinal, bulk-sound,
  thickness, density, order, electric, or reservoir modes. The homogeneous selected action also removes
  direct quadratic transverse–longitudinal mixing in \(D=3\), but finite-amplitude harmonics, slab
  gradients, interfaces, defects, and enlarged field content remain open radiation channels that must be
  calculated.

---

## Verification

**Two engines, written independently, agreeing on every computed value.** The blind Mathematica audit
was built **first**, imports nothing, and reads no registry; the SymPy audit was built while the `.wl`
was **quarantined out of the tree** — restored **byte-identical** to its committed blob (`55620a7f`), and
the builder's git usage was checked: `status` and `diff` only, ⛔ no `show`, no `cat-file`.

⭐⭐ **The engines DISAGREED on exactly one thing, and the disagreement was the valuable output.** The
`.wl` asserted a conclusion about transverse coupling; SymPy **refused** to conclude it. A review leg had
independently found that same line to be an unconditional literal the script *"does not earn."*
⇒ The `.wl` was repaired to **stop claiming what it does not compute** — ⛔ deliberately **not** to match
SymPy, since the blindness was already spent and steering one to match the other would have destroyed
the check retroactively. ⭐ It then **diagnosed its own circularity**, blind, and named the missing
input. **Convergence by honesty, not transcription.**

**Review legs: eight in total** — two on the build directives (which caught that the SymPy directive
omitted three tasks the `.wl` had, leaving the one surprising result with no second engine), two per
engine, and two per repair.

⭐ **A control hole was found, demonstrated, and closed.** Rewriting `R5` to `c_L − √(μ_R/ρ_br)` — ⚠ **the
exact claim this step exists to settle** — previously left **all five gates green**. Three guards
interlocked and all missed it. A new assertion now substitutes each **derived** root into the registry
residual that records it, ablation-verified against four distinct corruptions; `R4` covered too.
⇒ `DEFECT_REGISTER#f-r5`.

**Gates:** acceptance `MATCH` (`12→7, 12→7, 12→6`) · dimensional homogeneity `PASS`, `R5` HOMOGENEOUS
`[1,−1,0]`, 0 undetermined · able-to-fail `PASS` · 11 tests · both engine verdicts `PASS`.

### ⚠ Known limits — recorded, ⛔ not fixed

- **Two of six dimensional checks in the `.wl`** remain definitional identities `(X−2Y)+2Y == X` and
  ⛔ **cannot fail**; only the `ρ_br` one got an independent route. ⭐ The guard as a whole **is** able to
  fail — four independent ablations all caught, and a reviewer could not construct a single-point
  corruption that changed a printed dimension and still passed — but those two entries ⛔ **must not be
  counted as coverage**.
- **`ASSERTION_12` in the `.py`** is `c ≔ √(X)` then asserting `c² − X = 0`, identically zero by
  construction. Superseded in substance by the new registry control.
- **The `.wl` disclaims scope on the transverse tag but not the parallel one**, though both come from the
  identical substitution with the same absent operator. ⛔ Fixed **here**, in the record, rather than by a
  third repair round on pinned output lines: **the scope limit applies to both channels.**
- ⛔ **The SymPy script was repaired while UNTRACKED**, so no committed pre-repair baseline exists. A
  reviewer's `/tmp` copy and a full hand re-derivation of every load-bearing value substituted for one.
  ⚠ A straight violation of *commit before anything destructive*.
- **Disclosure, volunteered by a review leg:** a `grep` for `K_br` incidentally returned two lines of
  `V3_STEP_PLAN.md`, which was on its do-not-read list. The file was not opened and every reported result
  was derived beforehand.

---

## ⭐⭐ Strata-audit closure (2026-08-19) — the conclusion is independent of the census; the census's one physical family belongs to S11b

The exhaustive rank-drop / strata audit (spec §5 Q8a/Q8b — every `ρ×ρ` minor of each `M_r`, three
solve-variable sets, the degenerate strata) was subsequently run through **certified** census instruments
(`reduction/s11_*`; campaign closed `cbc49029`, four repair rounds, both engines, two legs + orchestrator
each). Measured result (`~/.s11_build/census_build4/`): the two engines **under-decide 917 sub-cases** and
carry finding-level gaps — 171 spurious branches, 72 omitted solve **records** (multiple missing
memberships each), 104 membership witness failures — plus the 7 registered defects (`DEFECT_REGISTER`
entries 4–7 + the obligation-4 instrument). ⚠ **Not all of these are "hard CAS questions":** the register
entries are genuine engine **logic** bugs (pointwise-vs-identical testing, joint-vs-pairwise coincidence,
false solution branches, projected witnesses); "known CAS limitations" understates them.

The census findings/under-decisions fall in **two kinds** (the first closure draft wrongly called all of
it "degeneracy strata" — an independent Codex check refuted that; corrected here):

- **(a) measure-zero eigenvalue-degeneracy strata** — `RANK_DROP` / `STACKED_DROP` / `ROOT_COINCIDENCE` /
  `STRATUM`: loci where already-known roots coincide or `M_r` drops rank. Genuine completeness
  bookkeeping. A degeneracy merges the kernels of known roots; it cannot create a new one — the
  characteristic polynomial is degree 3 in `ω²` with leading coefficient `−ρ_br³ ≠ 0`, so **exactly three
  finite roots always** (verified by independent factorization, both engines). ⚠ One nuance the first
  draft overstated: **at the coincidence `B_comp=μ_R` all three roots merge**, the full eigenspace is
  degenerate, and the longitudinal/transverse *character* dissolves on that measure-zero locus — no new
  mode, but the sector labels are not globally rigid.
- **(b) the `KW_ZERO_LOCUS` phase-matching family** (592 finding/under-decided lines, **88 on the MAIN
  package**) — its solution is `{B_comp → c_s0²ρ_br}` (and `{μ_R → c_s0²ρ_br}`), i.e. **`c_L = c_s0`**:
  exactly the `k_w = 0` **grazing threshold** of Move 6, the boundary between a kinematically closed and
  kinematically open propagating bulk-sound channel for the matched free longitudinal branch. ⭐ **This is
  physics, not bookkeeping**, but it is not by itself a bound-versus-radiating verdict. Interface overlap,
  boundary conditions, kinetic norm, bound-state-in-continuum possibilities, and leakage width belong to
  the brane–bulk interface and full-slab problem, **S11b**, not to any S11 engine round. The threshold value
  itself is decided; only the exhaustive status-token completeness on the locus is not.

⭐ **S11's own conclusion is independent of both.** The **computed result** above — the generic `M_ij`
eigenvalues (transverse `μ_R/ρ_br` nullity `D−1`, longitudinal `B_comp/ρ_br` nullity 1, no cross-modulus,
degeneracy exactly `B_comp=μ_R`, FORM-control uniqueness of the trace invariant) — is a generic-`k` fact,
its **root VALUES cross-engine agreed** (`ROOT{1,2}_N7_RESIDUAL=0`, both engines). ⚠ Scope caveats, none
touching the decoupling: (i) this is the **in-plane, frozen-wall-width** spectrum (`WALL_WIDTH_FIELDS={}`,
`hBranon` excluded, `INTERFACE_EQUATIONS_SUPPLIED={}`) — it does not decide inhomogeneous mode conversion,
leakage, or binding; (ii) root **positivity** follows trivially from `B_comp,μ_R,ρ_br>0`, `k²>0`, but the
SymPy sign **probe** leaves it `UNDECIDED`, so "both engines decide every generic verdict" is an
overstatement — the values are decided, one sign token is not; (iii) sector separation is a **`D=3`
fact** (at `D=2` the reflection-odd invariant mixes the sectors — Move 2), not structural; (iv) the
`PASS` gates live in this step record's acceptance harness, **not** in the two engine `.out` payloads.

### Ownership boundary after the grazing-threshold audit

**S11 owns:**

- the homogeneous \(D=3\) dynamical matrix;
- the two transverse roots and one longitudinal root;
- their selected quadratic linear cross-block being zero;
- the proof that the cubic characteristic polynomial has no hidden finite root for \(\rho_{\rm br}\ne0\);
- the simple **KW_ZERO_LOCUS** phase-matching threshold \(c_L=c_{s0}\).

**S11b and the nonlinear light program own:**

- the complete brane–bulk interface law and whether its overlap is nonzero;
- classification as a bound state, threshold state, bound state in the continuum, or resonance;
- actual leakage rates and behavior at and near \(k_w=0\);
- the full variable-coefficient slab spectrum;
- nonlinear intensity coupling and longitudinal participation in the photon helper;
- DC, harmonic, and sideband radiation into every receiving branch.

⇒ ⭐ **S11's stated conclusion stands closed.** Kind (a) is documented completeness bookkeeping; kind (b)
is handed to **S11b** where it already belongs. ⛔ **The S11 engine round was NOT run** — it would buy
completeness of the degeneracy census, not a change to any S11 conclusion (user's call, rule 11), and it
would not touch the KW/interface physics either. The certified instruments and all four census runs are
preserved, so the audit resumes exactly where it stands if a downstream step (or S11b's grazing case)
turns out to need a specific locus decided.
