# S9 · Light's requirement on the medium — shear on the brane, none in the bulk

**Sector 1 (light), step 1.** Walked side by side with the user, 2026-08-01.
⚠ First step banked under the **requirements-first** direction: light states what it needs; the medium's
structure is not assumed in advance.

---

## What it is

Light is a transverse wave. A GNLS superfluid cannot carry one. So light's first act is to state what
the medium must have that the GNLS does not.

## What it does

Introduces the two constants the entire light sector rests on, and fixes the light cone as their ratio.
Everything downstream in this sector is a consequence of these plus the curl-only form (S10, S11).

## The argument, in three moves

**1 — What Maxwell demands.** Two transverse polarisations · one speed · no longitudinal mode.

**2 — What a GNLS has.** Linearising Madelung about uniform `ρ₀` gives exactly **one** mode:

```
ω² = c_s²k² + (ħ²/4m²)k⁴        c_s² = (1/m)·dP/dρ = 5Kρ₀⁴/m
```

and it is **longitudinal**. The obstruction is not incidental — it is representation-theoretic:

> *"a single-component **SCALAR** superfluid cannot carry transverse light, period … `Ψ = √ρ e^{iθ}` is
> one complex scalar. Its excitations are `δρ` and `δθ`, both spin-0/longitudinal; one broken `U(1)` →
> one Goldstone phonon. **A scalar's spectrum cannot contain a spin-1 photon.** … This is
> **dimension-independent** for a pure scalar"*
> — `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:29-40`

⚠ **Status of that argument:** cited from an external review; ⛔ **no script in this repo executes it.**
A weaker in-repo form exists — `v = (ħ/m)∇θ ⇒ ∇×v ≡ 0`, so superfluid flow is curl-free — but it is
never wired into a no-shear proof (`docs/model_map.md:57` states the flow form only).

**3 — Therefore.** ⭐ **Light requires structure the GNLS does not contain.**

> *"**Brane elasticity is NOT a consequence of the GP/NLS mean-field (C1).** A GP/NLS superfluid is a
> *fluid* — zero shear modulus. The brane's shear rigidity therefore requires the **substructure**
> (constituents/cohesion *beneath* the mean-field), which the GP/NLS equation does not contain. Honest
> framing: **"GP/NLS as the effective/coarse-grained medium + a deeper substructure that supplies the
> brane elasticity."***
> — `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:230-233`

⚠⚠ **ONE SOURCE, NOT TWO — state it or the argument reads as stronger than it is.** Moves 2 and 3 both
come from `decisions/15`. ⛔ A draft of this record cited move 3 to `research/em_fields/paper/em_fields.tex:230`,
which is `\label{eq:B-vorticity}` — the phrase occurs **zero** times in that file (F4 class: right line
number, wrong file). ⇒ The whole *"a GNLS cannot carry light"* argument rests on **a single external-review
document with no executing script**, ⛔ not on two independent sources and ⛔ not on the model's own paper.

## ⭐⭐ The requirement, stated in full — and it is TWO-SIDED

The medium has **substructure** whose coarse-grained description is the GNLS (charter-level standing
postulate, ⛔ not light's to introduce). Light requires of that substructure:

| | requirement | what it buys |
|---|---|---|
| **ordered phase** (`χ_B = 1`) | carries **MacCullagh curl-only** shear | photons **exist** |
| **disordered phase** (`χ_B = 0`) | carries **no** shear | photons stay **confined** |

⭐ **Both halves are light's own bill.** ⛔ The second is not magnetism's, as an earlier draft had it:
without it a photon has somewhere else to go, light leaks off our 3D space, and that is an energy sink
with no observational room. ⚠ **Stressed hardest at the throat**, where the brane is bent into `±w`.

⭐ **A third consumer of the second half:** the recorded throat-support mechanism is an *outward trapped
standing-wave pressure* whose trapped mode is a **brane-shear standing wave**. If the bulk carried
shear that mode radiates away, the outward pressure vanishes, and the throat closes. ⇒ **Bulk
shear-freeness is what makes the geon possible at all.**

## The equations

```
L = ½ ρ_br (∂_t u)²  −  ½ μ_R (curl u)²                     (the brane's transverse sector)
c_γ² = μ_R / ρ_br                                            (R4 — the light cone)
```

## Dimensions — derived here, then checked against the register

`L` is an energy density on the 3D brane, `[L] = M L⁻¹T⁻²`. With `[u] = L`:

```
[∂_t u] = L T⁻¹                     [curl u] = L/L = 1  (dimensionless)

ρ_br (∂_t u)²  :  [ρ_br]·L²T⁻² = M L⁻¹T⁻²   ⇒   [ρ_br] = M L⁻³      = [-3, 0, 1]
μ_R  (curl u)² :  [μ_R]·1      = M L⁻¹T⁻²   ⇒   [μ_R] = M L⁻¹T⁻²    = [-1,-2, 1]

R4 :  μ_R − ρ_br = [-1-(-3), -2-0, 1-1] = [2,-2,0] = 2×[1,-1,0]  ⇒  c_γ = [1,-1,0]  ✓
```

⭐ **Independent agreement:** these were derived from the Lagrangian *before* consulting the register,
which records `μ_R` as `M L⁻¹ T⁻²` and `ρ_br` as `M L⁻³` (`research/pde_ledger_v2/notes/parameter_register.md:137`, `:138` — ⚠ **v2's** register, cited by locus; v3 does not have a `notes/` tree).

### ⭐⭐ Rebuilt 2026-08-04: the dimensions are now solved SYMBOLICALLY IN `D`, from the action

Both engines now extract the **derivative multi-order of every field factor from the Lagrangian's
expression tree** and solve the resulting linear system. ⛔ Nothing is read off a table and ⛔ no operator
is special-cased by name.

```
[ρ_br] = (−D,   0, 1)          [μ_R] = (2−D, −2, 1)
[μ_R] − [ρ_br] = (2, −2, 0)    ⭐ INDEPENDENT OF D
at D = 3:  (−3, 0, 1)  and  (−1, −2, 1)   ✓ registry
```

⭐⭐ **This retires a recorded S9 "blind spot" by explaining it.** The old audit's dimension check was
insensitive to the assumed brane dimension, filed as a weakness. It is not a weakness but an **identity**:
the speed's dimension **cannot** see `D`, because the difference is `(2,−2,0)` for every `D`.

⭐ **And the block is now able-to-fail, which it was not.** A control with a **different derivative count**
(flexural, `(∇²u)²`) moves the stiffness dimension to `(4−D,−2,1)`; a control with an **undifferentiated
field** (`−½ μ_G u·u`) yields `(−D,−2,1)`. ⚠ Before the rebuild the block emitted **byte-identical output
under a change of the action's form**, because the derivative counts were hand-encoded.

### ⭐⭐ FIVE independently-built engines derive these two dimensions, and three of them read no register

⚠ Measured 2026-08-08, and recorded here because the script that measured it was deleted the same day —
it belonged to an abandoned track, and its one physics result is this paragraph.

Five engines across three steps emit a symbolic dimension for `ρ_br` and `μ_R`, and every one emits the
**same** symbolic vector:

| engine | `[ρ_br]` | `[μ_R]` | registry references in the source |
|---|---|---|---|
| S9-py  | `(−D, 0, 1)` | `(2−D, −2, 1)` | **0** |
| S9-wl  | `(−D, 0, 1)` | `(2−D, −2, 1)` | **0** |
| S10-wl | `(−D, 0, 1)` | `(2−D, −2, 1)` | **0** |
| S10-py | `(−D, 0, 1)` | `(2−D, −2, 1)` | 14 |
| S11-py | `(−D, 0, 1)` | `(2−D, −2, 1)` | 18 |

⇒ ⭐ **three of the five could not have copied the answer**, because the register is not reachable from
their source at all. The reference count is `grep -c "registry\|reduction"` on each engine file.

⛔ **What this does NOT establish.** ⚠ The three register-free engines are not three independent
derivations of a *different* kind — S9-py and S9-wl solve the same linear system from the same posited
action, and S10-wl re-derives it in general `D`. ⇒ this corroborates the **solve**, ⛔ not the action.
A stiffness form authored identically into all five is caught by none of them ⇒ the common-mode limit
recorded under **VERIFICATION** below still stands at full width.

⚠ **`B_comp` is NOT part of this.** Only **one** engine (S11-py) emits a comparable dimension for it, so
it has no cross-engine corroboration whatever. ⛔ Do not report it alongside the two above. Its declared
`D = 3` specialisation is recorded in `REBUILD_HANDOFF.md`.

## What's new — the introduction inventory

| item | class | why |
|---|---|---|
| `ρ_br` | **postulated** | a property of the substructure; ⛔ not derivable from the GNLS |
| `μ_R` | **postulated** | ditto. ⚠ Deriving it from a polar substructure `P` returned `FAIL_COUPLE_STRESS_NOGO` (**B2**) — that route is closed, the postulate stands |
| `c_γ` | **derived** | `R4`, from the two above |
| MacCullagh curl-only **form** | **postulated** (structural) | forced by Maxwell's *no longitudinal mode* — ⭐ and **this is now COMPUTED, not merely asserted**; see the verification section |

⚠ **The far field consumes only the ratio `c_γ`.** ⛔ Do not let that become "the two are not separately
owed" — a **simulation** needs `μ_R` and `ρ_br` absolutely (→ the packet step).

## Inputs, by locus

| input | from |
|---|---|
| the standard GNLS (`ψ`, `ρ`, `θ`, Madelung, `v=(ħ/m)∇θ`) | charter §2 — comes standard |
| Maxwell's three demands | the imported theory |
| ⛔ nothing from any earlier v3 step | S9 is sector 1, step 1 |

## Regime

Linear response about an unstrained brane; small oscillations. ⛔ Everything here is quadratic — see the
packet step for what nonlinearity would require.

## Departure

⛔ **None yet at this step.** The departure enters at S11, when the medium's compressibility overrides
requirement three.

## Falsifier, and its status (charter §1.4(d))

**Requirement:** the substructure carries shear in its ordered phase and none in its disordered phase.
**Falsifier:** a substructure that supplies one and not the other — in particular, one that carries shear
in the *disordered* phase, which would unconfine light **and** kill the trapped-mode geon.
**Status: ⭐ LIVE, and it is a knit question.** Three independent consumers (photon propagation, photon
confinement, geon support), and ⛔ nothing yet shows one substructure can do both. ⚠ Bulk
shear-freeness is currently **POSTULATED** — `two_throat_simulation_handoff_spec.md:127` tokens it
`COMMITTED`; `ledger_stage006` files it under *"all POSTULATED"*; the only machine check on it is
**dimensional**.

⇒ ⭐ **This is light's first open requirement conflict**, and it is internal to light — ⛔ not a
light-versus-magnetism conflict. It is a question about the **ordering transition**: can one substructure
be ordered-and-shear-bearing in one phase and unstructured-and-shear-free in the other?

## ⭐⭐ VERIFICATION — rebuilt 2026-08-04 under the script-rebuild programme

⚠ **The physics did not change. The evidence did.** Every value the pre-rebuild engine computed is
identical today; ⛔ no result was revised, no sign flipped, no number moved.

### Two engines, independently constructed

| | engine | tags | typed conclusions |
|---|---|---|---|
| Mathematica | `mathematica/S9_light_requires_shear_mathematica_audit.wl` | 1559 | **0** |
| SymPy | `scripts/S9_light_requires_shear_sympy_audit.py` | 635 | **0** |

⛔⛔ **The SymPy engine's existence reverses a recorded exemption, and the reversal is itself a finding.**
`LAUNCH_PROMPT.md` OWED 2b claimed S9 needed no second engine because *"the cone is two lines of
algebra"*, under the **conditional** rule *"a second engine earns its place where the algebra is long
enough that it could genuinely DISAGREE."* ⚠ **The condition had expired** — 26 tags had become 1500+ —
and the conclusion was inherited anyway. ⇒ ⭐ **A conditional rule needs its CONDITION re-checked, ⛔ not
its conclusion quoted.** ⇒ [[feedback-attributions-are-my-paraphrase]].

### ⭐⭐⭐ What each verification mechanism caught — ⚠ and the two rows are DIFFERENT KINDS of claim

| defect | found by | missed by | kind of limit |
|---|---|---|---|
| ⛔ a **wrong dimension** — the implied speed dimension read only coefficient exponents, dropping an explicit wavevector factor | ⭐ **the second engine** | two review legs and a full ablation suite: the `.wl` computed it *consistently*, ablated *correctly*, and produced a number that **looked like a meaningful signal** | ⚠ **EMPIRICAL.** A leg that happened to dimension that expression by hand would have caught it. ⛔ Not a structural blindness — an **observed miss** |
| ⛔ a **wrong homogeneity test** — a `q`-substitution that silently no-ops, giving a false *positive* on a gapped root and a false **negative** on a genuinely dispersive one | ⭐ **reviewers deriving from scratch** (both legs, independently) | **cross-engine comparison** | ⭐⭐ **STRUCTURAL.** The defect came from the **shared directive**, so both engines computed it the same wrong way and **agreed**. ⛔ Comparison *cannot* see this class |

⇒ ⭐⭐ **This is the argument for running both mechanisms, demonstrated rather than asserted** — ⚠ but state
it precisely: **cross-engine blindness to a shared-directive defect is provable; review's miss is
measured, ⛔ not a theorem.** ⭐ A mechanism that missed a live defect is still not one to lean on alone.

### The checks that are now able-to-fail

- **Mode count.** `E2` = corank of `M` stacked on `kᵀ` — ⛔ **not** the `M·T = 0` test, which returns a
  **false negative** under anisotropic inertia (a transverse mode whose partner sits at a different
  frequency). `E2 = 2` at the propagating root, `0` at `ω²=0`; it drops to **1** when the inertia is made
  anisotropic. ⚠ Both the pre-rebuild engine and the first rebuild used the broken test.
- **Non-dispersiveness** by **scaling** `k → λk`: `ω²(λk) − λ²ω²(k)`. Vanishes for the main action; is
  **nonzero** for flexural and for a gapped root. ⛔ The superseded `q`-substitution family is **deleted**,
  not flagged — a flagged wrong value gets quoted without its flag.
- **The dynamical matrix has two independent routes** (position-space Euler–Lagrange vs. the quadratic
  form of the ansatz-substituted Lagrangian), residual zero in all nine actions. ⭐ **Independence proven
  by one-sided corruption:** breaking either route leaves the other byte-identical. ⚠ A reviewer also
  showed the two disagree on a **gyroscopic** term — they are structurally different enough to detect a
  real physical coupling.
- **Sign of `ω²` emitted per root**, so the sign-flip control can fail through the polarisation tests.
  ⚠ Without it, two **exponentially growing** modes carried the same E1/E2/E4 signature as two waves.
- **Nine actions**: main + 2 coefficient controls + 5 form controls + 1 sign control.

### ⭐⭐ A PHYSICS SHARPENING, and it comes from a control

⛔⛔ **Transverse propagation does NOT require the curl-only form.** The gradient-elastic control
`−½ μ_R Σ(∂_i u_j)²` — an ordinary elastic solid — carries **two transverse modes at the same
`c² = μ_R/ρ_br`**. It simply *also* propagates the longitudinal.

⇒ ⭐ **What curl-only buys, and the gradient-elastic form does not, is the ABSENCE of a propagating
longitudinal mode** — Maxwell's third
demand — ⛔ not the presence of the transverse ones. The record's `what's new` table already framed it
that way; ⭐ it is now **computed** rather than asserted, and the divergence-only control shows the
converse: the roles **swap**, the transverse pair drops to `ω² = 0` and the longitudinal propagates.

⚠ ⇒ **The defensible statement of S9 is conditional on the stiffness form**, and any prose reading
*"light requires shear"* as established by this computation is **stronger than the computation supports.**

### Automated consumption — ⚠ re-measured 2026-08-06, and the earlier figure was ~4× too strong

`reduction/engine_output_checks.py --config reduction/checks_S9.yaml` over both engines:

```
CROSS_ENGINE: agree=12 disagree=0            CROSS_ENGINE_COVERAGE: 12/12
NAMING_EFFECT: legacy_before_agree=8  declared_after_agree=12
CONTROL_RESPONSE[wl]: compared=3 responsive=3 invariant=0    CONTROL_RESPONSE[py]: 3/3 responsive
DIMENSIONS[wl]: total_tags=1559  compared=312  homogeneous=312  no_comparison=1228  unassessable=18
DIMENSIONS[py]: total_tags=635   compared=329  homogeneous=329  no_comparison=269   unassessable=37
REGISTRY_RESIDUAL: nonzero=1   (R4 — it compares a squared speed against a speed)
```

⛔⛔ **This record previously cited "dimensions 1219/1219 homogeneous".** Under the rebuilt harness that
figure is **312** actual comparisons on the Mathematica side and **329** on the SymPy side. ⇒ the cited
number counted **tags reached**, ⛔ not **comparisons made**, and overstated the check by roughly a factor
of four. ⭐ The check still passes everything it examines; ⛔ it examines about a quarter of what this
record claimed. ⚠ `DEFECT_REGISTER.md#f7` already recorded the same dent from the other direction.

⚠ **And four of the twelve agreements rest on declared spellings**, not on the engines independently
producing the same text: `legacy_before_agree=8` against `declared_after_agree=12`. Removing one
declaration at a time (`reduction/measurements/declaration_load_ablation.py`) shows exactly three are
load-bearing here, each on one row — `mu_F` (`AGREE → DISAGREE`), `lambda_scale`
(`AGREE → NAMING_MISMATCH`), and the `omega2 = omega**2` **algebraic identity** on `factored_determinant`
(`AGREE → DISAGREE`). The other four naming exceptions move nothing on this step. ⭐ Each declaration is
printed on its verdict line with a reason; ⚠ the `omega2` one is an **identity, not a spelling**, and is
declared separately for that reason.

⚠ **The control figure changed population, not value.** This record previously cited *"150 of 170 tags
respond to some control"*. The rebuilt harness scores **declared control rows** — `3` per engine, all
`RESPONSIVE` — rather than counting responsive tags. ⛔ The two numbers are not comparable and the smaller
one is **not** a regression.

⛔ Config maps **tag names only**; every comparison target is generated at compare time. ⚠ Its
`INVARIANT` and `PARITY` outputs are **triage lists, ⛔ not failures**, and a `DISAGREE` needs adjudication
too — one was a symbol-spelling false positive.

## ⛔ WHAT THIS STEP STILL DOES NOT ESTABLISH

- ⛔ **P2 — that a scalar superfluid carries no transverse mode.** Cited from one external review; ⛔ **no
  script executes it**, and the rebuild did not add one. It remains a **supplied premise**.
- ⛔ **The absence of a propagating longitudinal wave is ASSUMED, not derived.** The curl-only action sets
  the longitudinal restoring stiffness to **zero by construction**, so `ω² = 0` there is the postulate
  restated. ⛔ It removes the restoring *force*, ⛔ **not the degree of freedom** — which is why the
  longitudinal test is emitted at every root.
- ⛔ **Bulk shear-freeness** is postulated; the bulk is absent from the action, so nothing here tests it.
- ⚠ **The assumption set is computationally INERT.** Removing it entirely leaves every computed value
  byte-identical. ⭐ The physics is nonetheless correctly bounded — by generic rank **plus explicitly
  emitted exceptional loci**, and a reviewer verified that enumeration is **complete** for the anisotropic
  control (roots collide exactly where `ρ_z = ρ_br` or `k_x = k_y = 0`, and both are emitted).
  ⇒ ⭐ **Read the polarisation values as generic-locus values with the exceptions listed.**
- ⚠ **The `.py` was written AFTER the `.wl`**, so this is ⛔ not the blind-first ordering. Construction was
  independent (different naming, granularity, idioms; a hand-rolled Euler–Lagrange where SymPy has a
  built-in) and the builder was barred from reading the `.wl` — ⭐ but agreement here is **weaker evidence**
  than at a step where both engines predate any comparison.
- ⚠ One tag remains **unparsed** by the consumer: the anisotropic third root's sign is a `Piecewise`,
  because the sign is conditional on the parameter domain and the assumptions are inert.
- ⚠ Limits taken: sharp zero-width sheet · `v₀ → 0` · no dissipation · frequency-independent moduli ·
  continuum limit · amplitude → 0.

## Registry additions (executed at this step)

`Q.brane.rho_br` `[-3,0,1]` · `Q.brane.mu_R` `[-1,-2,1]` · `Q.brane.c_gamma` `[1,-1,0]` · relation
**`R4`** `c_γ = √(μ_R/ρ_br)`. ⭐ `R4` is the corpus's existing name for this relation
(`notes/parameter_register.md:271`). ⚠ `c_γ` re-enters here with provenance, having been deleted from
the *medium* block at S0.5 — that round trip is the mechanism working.
