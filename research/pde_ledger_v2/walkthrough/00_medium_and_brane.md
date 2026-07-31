# Phase 0 · Step 1 — the medium, and the brane as its ordered state

**Upstream:** nothing. This is the first step; everything here is new by construction.
**Method:** `docs/derivation_walkthrough_plan.md`. **Sources:** `docs/model_map.md:20-26`, `:49-59`.

---

## ⭐ A scoping finding, before the content

The step was asked for as "the brane". **It cannot be stated without the medium**, because the brane is
not a separate object — it is *the medium in its ordered state*. `model_map.md:22-24`: one medium, two
material states, distinguished by an order field `χ_B∈[0,1]`; brane = `χ_B=1` (ordered, shear-supporting),
bulk = `χ_B=0` (*"the SAME medium, de-structured and shear-free"*).

⇒ Recorded as one step. ⭐ **This matters for the count:** brane and bulk do **not** contribute two sets
of defining properties. They contribute one medium plus one order field. A two-object treatment would
have double-counted the substrate.

---

## What it is

**One compressible superfluid** — a Gross–Pitaevskii / nonlinear-Schrödinger (GNLS) condensate in a **4D**
bulk. Order parameter `ψ`; number density `ρ = |ψ|²`; phase `θ`; flow `v_b = (ħ/m)∇θ`. Closed by a stiff
polytropic equation of state `P = Kρ⁵`.

Our 3D space is a **domain-wall brane at `w = 0`** in that 4D bulk — the ordered, shear-supporting slice
of the same substance.

## What it does

Provides the single substrate every later sector is a mode of: a density mode (→ gravity), a
shear-supporting ordered slice (→ light and magnetism), and a puncturable wall (→ charge). ⭐ Nothing
downstream may introduce a second medium; that is the constraint this step buys.

## ⚠ What it does NOT do — the honest status

⛔ **The two-phase split is POSTULATED, not derived.** `model_map.md:26` states it plainly: the parent
potential `U(ρ)` is **single-well**, so the brane does **not emerge** from the one medium — *"it is put in
by hand"*.

⇒ The existence of an ordered phase is a **postulate**, not a consequence. Recorded here rather than
discovered at integration.

---

## What's new

⭐ Everything below must be **supplied for a simulation to run** (§1.0). Partitioned by kind, because the
kinds are not interchangeable.

### Scalars — 4 primitives
| quantity | dimension | class | why |
|---|---|---|---|
| `ħ` | `M L² T⁻¹` | **postulated** (external constant) | no defining equation in the model |
| `m_GNLS` | `M` | **postulated** | no defining equation |
| `K` | `M L¹⁸ T⁻²` | **postulated** | EOS stiffness; no defining equation |
| `ρ0` | `L⁻⁴` | **postulated** | asymptotic density datum |

`model_map.md:49` — *"Primitive constants of the medium (~4 declared universal)"*.

### Discrete / structural choices — 3
| choice | class | note |
|---|---|---|
| **bulk dimensionality `D = 4`** | **postulated** | the ambient is 4 spatial dimensions |
| **EOS polytropic exponent `= 5`** | ⭐ **calibrated** | per §1.1: chosen because the calibration was necessary for the model to work |
| **the two-phase split exists** (`U(ρ)` admits an ordered state) | **postulated** | `:26` — put in by hand, and the honest core of this step |

### Fields / profiles — 3
| field | dimension | note |
|---|---|---|
| `ψ` (hence `ρ=\|ψ\|²`, `θ`) | `ρ`: `L⁻⁴` | the order parameter itself |
| `χ_B ∈ [0,1]` | dimensionless | ⚠ **an independent scalar field, explicitly NOT `\|P_∥\|²`** (`:54`) |
| `U(ρ)` — the parent potential | — | single-well (`:26`); its **form** is an input |

⚠ A field is **not** one number. For a simulation each needs a profile or its dynamics — which is why
§1.0 partitions by kind rather than summing to an integer.

### Boundary / domain data — 1
| item | note |
|---|---|
| the wall at `w = 0`, and the asymptotic condition `ρ → ρ0` | a domain and a BC; must be supplied to run |

---

## Equations introduced

| id | relation | status |
|---|---|---|
| `EOS` | `P = K ρ⁵` | **postulated closure** (`:20`) — `K` primitive, exponent calibrated |
| `FLOW` | `v_b = (ħ/m)∇θ` | definition of the flow from the phase |
| `DENS` | `ρ = \|ψ\|²` | definition |

⛔ **No sound speed here.** `c_s² = 5Kρ⁴/m` is *derived from* the EOS (`:62`) and belongs to step 2.
Recording it here would import a derived quantity into the postulate set.

---

## Checks

| # | check | result |
|---|---|---|
| 1 | **dimensional homogeneity** | see below — `EOS` verified |
| 2 | **what's-new classification** | done above; 4 scalars, 3 discrete, 3 fields, 1 BC |
| 3 | **input provenance** | n/a — no upstream inputs; sources cited by locus |
| 4 | **step-to-step identity** | n/a — first step |
| 10–13 | fidelity · dynamical health · approximation validity · second physics leg | **not yet applicable** — no derivation performed, no dynamical block, no approximation, nothing derived to review |

**Dimensional check of the EOS, by hand, to be confirmed by the gate:**
`[K][ρ]⁵ = (M L¹⁸ T⁻²)(L⁻⁴)⁵ = M L¹⁸⁻²⁰ T⁻² = M L⁻² T⁻²`.
In 4D, pressure = energy / 4-volume = `(M L² T⁻²)/L⁴ = M L⁻² T⁻²`. ✅ consistent.

---

## Open, carried forward

1. ⚠ **`ħ`'s class.** No defining equation ⇒ not `derived`. But is it `postulated` (a property of the
   medium) or `calibrated` (an external constant we import)? `decisions/14:49` calls it *"fixed
   unit/external action constant (underived)"*. ⛔ Either way it is in the sim-input set; the label only
   moves it between tier 1 and tier 2. **User call.**
2. **`U(ρ)`'s form** is an input, and "single-well" does not fix it. How much of the later model depends
   on more than that one bit?
3. ⭐ **The two-phase postulate is the single largest item here.** If a later step derives an ordered
   phase from `U(ρ)`, this becomes `debt` rather than `postulated` — and the tier-1 core shrinks.
