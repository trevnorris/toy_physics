# DIRECTIVE — S11b-B SymPy audit and registry insertion

**Deliverables (absolute paths):**
1. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bB_interface_assembly_sympy_audit.py`
2. any additions to `/var/projects/toy_physics/research/pde_ledger_v3/reduction/quantities.yaml` and
   `relations.yaml`

Run the script and the existing gates. Iterate until everything exits 0. Then stop and exit — ⛔ do not
write a report or a summary.

## ⛔⛔ ONE OF TWO INDEPENDENT ENGINES

A blind Mathematica audit of the same physics exists. **The disagreement between the engines is the test**,
so an agreement reached through a shared source certifies nothing.
⇒ **The read-bar list is §0b of the shared specification below**, byte-identical to the other engine's.

## Registry work — this engine only, and ⛔ AFTER the physics

⚠ **Do the physics tasks FIRST, from the specification alone.** ⛔ Do not consult `reduction/` while
deriving. ⚠ **Measured hazard:** registry access lets one engine silently identify two symbols the other is
treating as independent, after which their "agreement" certifies an identification only one of them made.
⭐ §1 fixes every symbol's name and meaning explicitly, so there is nothing to infer — ⛔ do not infer
anyway.

**Then**, and only then:
1. Import the registry through `reduction/registry_read.py`. Do not re-implement its loading or validation.
2. Add the quantities and relations this step earns. Give each a `kind`, a `counting_axis`, a dimension
   with provenance, and `aliases: []` unless an alias genuinely exists. ⚠ A quantity this step must
   **postulate** because a later step owes its derivation is recorded as **postulated**, ⛔ never dressed
   as derived.
3. ⛔⛔ **EVERY NEW RELATION NEEDS ITS OWN ALGEBRA ASSERTION, because no gate checks relation algebra.**
   Substitute the expression **this script derived** into the registry residual that records it and assert
   it vanishes. Then **ablate** it — corrupt the residual, confirm the assertion fails — and report the
   literal output. ⚠ Measured: a relation was once rewritten to assert the exact claim its step existed to
   settle and **all five gates stayed green**.
4. ⛔ A check of the form `c := sqrt(X)` then asserting `c**2 - X == 0` is identically zero by
   construction. If you write one, label it definitional in its own output line.
5. **Run all gates** and report each one's literal output: acceptance, dimensional homogeneity,
   able-to-fail, pytest. ⚠ If the continuous payload changes, that is a **result to report**, ⛔ not
   something to make match by adjusting a fixture.

## Script conventions

- ⚠ **B7's dimensions must be DERIVED from the specification's equations**, ⛔ not read from
  `quantities.yaml` or the homogeneity gate. The registry is what B7 independently checks.
- One tag per line, `TAG: value`, no `WL_` prefix.
- Keep total runtime under **10 minutes**; runners get `timeout 600`.
- ⚠ **On a non-uniform background do not silently assume plane waves are eigenmodes.** Where you must
  restrict to a tractable form, ⛔ say so in the emitted value rather than in a comment.

---

# S11b-B — SHARED PHYSICS SPECIFICATION

⚠ **Inserted BYTE-IDENTICALLY into both the Mathematica and the SymPy directive.** It is the only part
they share, and it exists because on an earlier step the two engines were handed *different task lists*,
leaving the one surprising result with no second engine.

---

## 0 · What this step is

Assemble the brane's in-plane sector, the slab's thickness degree of freedom, and the bulk's response into
**one linear system on a NON-UNIFORM background**, and from it determine the physical status of the
brane's longitudinal mode and whether the transverse mode is coupled.

⛔ **Do not consult any file for the answers.** Everything needed is below. If something needed is missing,
emit the tag as `NOT_ESTABLISHED` and name what is missing.

## ⛔⛔ 0b · WHAT NEITHER ENGINE MAY READ — identical for both, and that is the point

⚠ A bar maintained separately in two headers drifts, and the engine with the shorter list becomes a
transcriber while the other derives.

⛔ **Neither engine may read, open, grep, `cat`, `git show`, or `git cat-file`:**
`research/pde_ledger_v3/scripts/` · `research/pde_ledger_v3/mathematica/` ·
`research/pde_ledger_v3/steps/` · `research/pde_ledger_v3/V3_STEP_PLAN.md` ·
`research/pde_ledger_v3/directives/S11b_*` and `S11bA_*` · `research/pde_audit/` ·
any file whose name contains `PREREGISTERED` or `PREREG` · **the other engine's deliverable for this step**.
⚠ `git status` and `git diff` are fine. ⛔ Everything you need is below; do not go looking.

## 1 · Geometry, notation, conventions

Four spatial dimensions `(x¹,x²,x³,w)`, `x = (x¹,x²,x³)`. A slab of thickness `W(x,t) = W₀ + δW` centred on
`w = 0`, faces near `w = ±W₀/2`, `D_brane = 3`. Bulk on both sides.

| symbol | meaning |
|---|---|
| `u(x,t)` | in-plane displacement, **3 components**, ⛔ no `w`-component |
| `δW(x,t)` | thickness perturbation. Faces sit at `ζ_± = ±δW/2` |
| `ρ_4D(x)` | local 4D mass density of brane material · `ρ_br ≡ ρ_4D W` (slab-integrated) |
| `μ_R` | curl-type (twist) modulus of the brane, established input |
| `B_ρ` | local 4D compression modulus of brane material · `B_ρ⁽³⁾ ≡ B_ρ W₀` |
| `μ_W`, `k_W` | thickness inertia and restoring stiffness — **symbolic**, ⛔ assume no form |
| `ρ_m`, `c_s0` | bulk mass density and sound speed |

**Conventions**, fixed so the engines produce sign-comparable results: every perturbation
`∝ exp[i(k·x − ωt)]` with `ω > 0`; face displacements measured along global `+w`; outward normals `+ŵ`
upper, `−ŵ` lower; outward face velocity `V_± ≡ (∂_tζ_±)(±1)`; response ratio `Z ≡ (pressure at a face) /
(OUTWARD normal velocity of that face)`; in each half-space retain only waves carrying energy **away** or
decaying away — ⛔ no incoming waves from infinity.

## 2 · ⛔⛔ THE BACKGROUND IS NON-UNIFORM. THIS IS NOT OPTIONAL.

`ρ_4D`, `B_ρ`, `μ_R`, `W₀` are **functions of `x`**. ⛔ Do not set them constant, and ⛔ do not linearise
about a uniform state except where a task explicitly asks for that limit as a control.

⚠ **Why it is not optional:** in a uniform brane the transverse mode's coupling vanishes by symmetry, so a
uniform calculation can return only one answer and tests nothing. One of this step's questions is whether
that decoupling is **unconditional** — which is precisely a question about non-uniformity.

## 3 · Established input from the preceding sub-step — ⛔ do NOT re-derive

The bulk's response to a face moving with outward normal velocity `V`, already verified by two independent
engines:

```
q² = ω²/c_s0² − k²            Z = ω ρ_m / q_out
q² > 0  propagating   q_out = q      Z real          ⇒ radiation resistance
q² < 0  evanescent    q_out = iα     Z = −iωρ_m/α    ⇒ inertial loading ρ_m/α per face
q² = 0  grazing       Z singular
```

**Permeable faces** carry a relative mass flux `J_± = Λ_p(ω)δp + Λ_V(ω)V_±` with
`Λ(ω) = Λ⁰/(1 − iωτ)`, `Λ_p⁰`, `Λ_V⁰`, `τ` real and free, `τ ≥ 0`, giving

```
Z_perm = (ρ_m r + Λ_V⁰)/(y r − Λ_p⁰) ,    r = 1 − iωτ ,    y = q_out/ω
```

⛔⛔ **THREE TRAPS, each measured in that sub-step:**
1. The per-face inertial loading is `+ρ_m/α` **against the outward acceleration on BOTH faces**. The
   signed pair `(ρ_m/α, −ρ_m/α)` is an artifact of the global-`w` convention. ⛔ **Do not consume a signed
   pair as an inertia.**
2. ⛔ **Propagating `Re Z` is radiation resistance and exists with impermeable faces.** It is ⛔ **not**
   evidence of leakage through the interface. Only the **evanescent** `Re Z` is created by the closure.
3. The rest-frame linearisation of the bulk discards a relative correction of order **`v₀|q_n|/ω`**, ⛔ not
   `v₀/c_s0`. ⚠ In the evanescent regime this **exceeds first order whenever `k c_s0 ≫ ω`**. `v₀` is the
   steady background normal flow speed. ⇒ ⛔ do not assume the discarded term is negligible; report where
   it is not.

## 4 · The brane's stored energy

⛔⛔ **DO NOT INTRODUCE A SINGLE IN-PLANE COMPRESSION MODULUS.** Build the compression sector from the two
physical channels below. ⚠ A step-3 calculation reported one such modulus measured with the thickness held
fixed; ⛔ writing it here alongside a dynamical thickness would count the wall channel twice, and neither
engine could detect it. ⭐ Where that quantity sits is an **output** of task **B3**, ⛔ never an input.

Twist is unchanged from earlier steps and is carried by `μ_R`. Compression is stored in **two** channels —
densification of the brane material (modulus `B_ρ`) and change of thickness (stiffness `k_W`) — and the
thickness carries inertia `μ_W`.

---

## TASKS

⛔⛔ **RULES THAT OVERRIDE EVERY TASK.**
1. Every reported value must come from a computation in the script. ⛔ A printed assertion is not a result.
2. ⛔ If a task cannot be completed, emit its tag as `NOT_ESTABLISHED` and say which input is missing.
   **A refusal is a valid and valuable output.**
3. ⛔ Never silently choose a branch, closure, convention, or expansion. If a choice is required and is not
   fixed above, introduce a **free symbol**, say so, and report the answer as a function of it.
4. ⛔⛔ **No task below tells you the FORM of its answer, and that is deliberate.** ⚠ On the preceding
   sub-step a task asked for "the order in `v₀/c_s0`", both engines returned a power of `v₀/c_s0`, and the
   true answer was a different quantity entirely — the specification had presupposed the form and the
   engines agreed on it. ⇒ **Report what you derive, ⛔ not what the phrasing suggests.**

---

**B1 · The constraint.**
Impose conservation of slab material per unit `x`-3-volume, with a source from flux across the two faces,
and **linearise about the non-uniform background of §2**. Report the resulting relation among the
thickness perturbation, the densification, the in-plane displacement, and the face flux.
⭐ State explicitly **how many independent internal degrees of freedom survive** this relation, and why.
⇒ `S11BB_CONSTRAINT`, `S11BB_INTERNAL_DOF_COUNT`

**B2 · The thickness equation.**
From §4's energy and §3's bulk response, derive the equation of motion for the thickness perturbation
including its inertia and the force exerted by the bulk on **both** faces. Report the response function
that relates the thickness perturbation to whatever drives it, and report the effective inertia in each of
§3's regimes.
⇒ `S11BB_THICKNESS_RESPONSE`, `S11BB_EFFECTIVE_INERTIA_BY_REGIME`

**B3 · The effective in-plane compression modulus.**
Eliminate the thickness degree of freedom and report the effective in-plane compressional stress response
as a function of `ω` and `k`. Report its `ω → 0` and `ω → ∞` limits.
⭐ Then report **which limit, if either, corresponds to a calculation performed with the thickness NOT a
degree of freedom.** ⛔ If neither limit yields a consistent identification, report that as the result.
⇒ `S11BB_B_EFF`, `S11BB_B_EFF_LIMITS`, `S11BB_FROZEN_THICKNESS_IDENTIFICATION`

**B4 · The longitudinal mode.**
Assemble the longitudinal in-plane equation with B3's response and report the dispersion relation. Report
whether it can be solved in closed form for `ω(k)`, and whether the root is real. If complex, report the
imaginary part and **which physical ingredient makes it nonzero**, ⛔ distinguishing the two mechanisms of
§3 trap 2.
⇒ `S11BB_LONGITUDINAL_DISPERSION`, `S11BB_ROOT_IS_REAL`, `S11BB_IMAGINARY_PART`,
  `S11BB_DISSIPATION_ORIGIN`

**B5 · ⭐⭐ The transverse channel — COMPUTE IT.**
On the non-uniform background of §2, compute the coupling between the **transverse** in-plane mode and the
thickness degree of freedom, **including the back-reaction on the transverse equation of motion**. Report
the coupling coefficient and the resulting modification to the transverse dispersion.
⭐ **Report what quantity controls its magnitude.** ⛔ Do **not** assume it is proportional to any
particular power of any particular background gradient — derive what appears.
⚠ ⛔ **Do not infer this from a divergence-free condition.** The route must be capable of returning a
nonzero result.
⇒ `S11BB_TRANSVERSE_COUPLING`, `S11BB_TRANSVERSE_DISPERSION_SHIFT`, `S11BB_COUPLING_CONTROLLED_BY`

**B6 · Does the transverse mode dissipate?**
Report whether the transverse mode acquires an imaginary frequency component through B5's coupling, on
which of `Λ_p⁰`, `Λ_V⁰`, `τ` it depends, and how it behaves across the full range of `ωτ` — small, order
unity, and large. ⛔ Report all three; ⛔ do not report one regime.
⇒ `S11BB_TRANSVERSE_DISSIPATION`, `S11BB_TRANSVERSE_DISSIPATION_VS_OMEGA_TAU`

**B7 · Dimensions.**
Derive, from the equations above and ⛔ not from any external table, the `[L,T,M]` exponents of: `B_ρ`,
`B_ρ⁽³⁾`, `μ_W`, `k_W`, B2's response function, B3's modulus, and B5's coupling coefficient. Show the route
for each and label it **independent** or **definitional**. ⚠ A route whose asserted equation *defines* the
symbol under test is definitional.
⇒ `S11BB_DIM_<name>`, `S11BB_DIM_ROUTE_KIND_<name>`

**B8 · Controls.** ⛔ FORM controls; ⛔ a coefficient rescaling tests neither.
- **A — the uniform limit.** Set every background gradient to zero and recompute B5. ⭐ **This is a
  required control**: report whether the transverse coupling vanishes **identically**, and whether B3's
  modulus is unchanged.
- **B — remove the thickness channel** (hold the thickness fixed) and recompute B3 and B4. Report what
  collapses.
- **C — impermeable faces** (`Λ_p⁰ = Λ_V⁰ = 0`) and recompute B4 and B6. Report which imaginary parts
  survive. ⚠ Per §3 trap 2, one mechanism should and one should not; ⛔ report what you find, not that.
⇒ `S11BB_CONTROL_UNIFORM`, `S11BB_CONTROL_NO_THICKNESS`, `S11BB_CONTROL_IMPERMEABLE`

**B9 · Validity.**
Report the conditions under which this step's linearisations hold, including the background-flow condition
of §3 trap 3, and **where in `(ω,k)` any of them fail**.
⇒ `S11BB_VALIDITY_CONDITIONS`, `S11BB_VALIDITY_FAILURE_REGION`

---

## OUTPUT FORMAT

One line per tag, `TAG: value`. Values must be explicit expressions wherever mathematical. End with a
single `VERDICT:` line reporting whether the script's own internal consistency checks contradicted each
other. ⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ It is not a
verdict on the physics.
