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

> *"**Brane elasticity is NOT a consequence of the GP/NLS mean-field.** A GP/NLS superfluid is a *fluid*
> — zero shear modulus. The brane's shear rigidity therefore requires the **substructure** … which the
> GP/NLS equation does not contain."* — `research/em_fields/paper/em_fields.tex:230`

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
which records `μ_R` as `M L⁻¹ T⁻²` and `ρ_br` as `M L⁻³` (`notes/parameter_register.md:137`, `:138`).

## What's new — the introduction inventory

| item | class | why |
|---|---|---|
| `ρ_br` | **postulated** | a property of the substructure; ⛔ not derivable from the GNLS |
| `μ_R` | **postulated** | ditto. ⚠ Deriving it from a polar substructure `P` returned `FAIL_COUPLE_STRESS_NOGO` (**B2**) — that route is closed, the postulate stands |
| `c_γ` | **derived** | `R4`, from the two above |
| MacCullagh curl-only **form** | **postulated** (structural) | forced by Maxwell's *no longitudinal mode*; see S10/S11 for whether it delivers |

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

## Registry additions (executed at this step)

`Q.brane.rho_br` `[-3,0,1]` · `Q.brane.mu_R` `[-1,-2,1]` · `Q.brane.c_gamma` `[1,-1,0]` · relation
**`R4`** `c_γ = √(μ_R/ρ_br)`. ⭐ `R4` is the corpus's existing name for this relation
(`notes/parameter_register.md:271`). ⚠ `c_γ` re-enters here with provenance, having been deleted from
the *medium* block at S0.5 — that round trip is the mechanism working.
