# DIRECTIVE — S11b-A SymPy audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bA_interface_response_sympy_audit.py`

Run it. Iterate until it exits 0. Then stop and exit — ⛔ do not write a report or a summary document.

## ⛔ WHAT YOU MUST NOT READ

- `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bA_*` — a blind Mathematica audit of
  this same physics exists. It has been moved out of the tree; ⛔ do not retrieve it from git either
  (`git show`, `git cat-file`, `git checkout`). **The disagreement between the two engines is the test.**
- `/var/projects/toy_physics/research/pde_audit/` — any file. ⚠ Prior work on adjacent physics whose
  current validity is **not established**. The other engine is barred from it too; an asymmetry here would
  make one engine a transcriber and the other a deriver, and their agreement would then mean nothing.
- any file whose name contains `PREREGISTERED`.

⛔ `git status` and `git diff` are fine. `git show` and `git cat-file` are not.

## ⛔⛔ NO REGISTRY WORK IN THIS SUB-STEP

⚠ Normally the SymPy engine is the only one that sees `reduction/`, and that is deliberate. ⭐ **Not here.**
This sub-step derives no brane quantity and earns no registry row; the assembly that would is deferred to
S11b-B.

⇒ ⛔ **Do not read, import, or modify `reduction/`** — not `quantities.yaml`, not `relations.yaml`, not
`registry_read.py`, not the gates. ⚠ **The reason is a measured one:** giving one engine access to
recorded relations lets it silently identify two symbols the other engine is treating as independent, and
the two then "agree" on an identification only one of them made. With the registry out of reach, both
engines face **identical** information, which is the condition their agreement is supposed to certify.

⚠ **A6's dimensions must be DERIVED from this specification's equations.** ⛔ Do not read them from any
table, registry, or gate.

## Script conventions

- Output one tag per line, `TAG: value`, with no `WL_` prefix.
- Keep total runtime under **10 minutes**; runners get `timeout 600`.
- ⛔ Do not write a test file, a fixture, or a gate — this sub-step produces one script and its printed
  tags.

---

# S11b-A — SHARED PHYSICS SPECIFICATION

⚠ **Inserted BYTE-IDENTICALLY into both the Mathematica and the SymPy directive.** It is the only part
they share, and it exists because on an earlier step the two engines were handed *different task lists*,
leaving the one surprising result with no second engine.

⚠⚠ **This is a SUB-STEP, deliberately narrow.** Two attempts to specify the whole brane–bulk interface in
one pass were rejected by review before any build ran — the second introduced new defects *caused by* its
own fixes. ⇒ ⭐ **The assembly is deferred to S11b-B and will be written against the objects THIS step
produces, not against a description of them.**

---

## 0 · Scope

**IN scope.** What the bulk does when the slab's faces move, and the identity that projects a
four-dimensional balance law onto the slab.

⛔⛔ **OUT of scope — do NOT attempt, and do NOT let any task drift into them.** The brane's in-plane
elastic sector; any modulus of the brane; the coupled brane+thickness+bulk dispersion; the transverse
mode's coupling; any statement about whether a mode radiates, is bound, is confined, or is observable.
⚠ **Those are S11b-B's.** ⭐ A tag here that answers one of them is a defect, not a bonus.

⛔ **Do not consult any file for the answers.** Everything needed is below. If something needed is
missing, emit the tag as `NOT_ESTABLISHED` and name what is missing.

## 1 · Geometry, notation, and ⭐ CONVENTIONS

Four spatial dimensions `(x¹,x²,x³,w)` and time `t`; `x = (x¹,x²,x³)`. A slab of thickness `W₀` centred on
`w = 0`, with faces at `w = ±W₀/2`. Bulk on **both** sides. `w` is normal to the slab.

⭐⭐ **Conventions, fixed here so the two engines produce comparable signs.** ⚠ Reported real and imaginary
parts are meaningless for comparison unless these are shared.

| | convention |
|---|---|
| harmonic dependence | every perturbation `∝ exp[i(k·x − ωt)]`, with `ω > 0` and `k` real |
| face displacements | `ζ₊(x,t)`, `ζ₋(x,t)` = displacement of the **upper / lower** face, **both measured along global `+w`** |
| parity combinations | **thickness** `δW ≡ ζ₊ − ζ₋` · **bodily shift** `h ≡ (ζ₊ + ζ₋)/2` |
| outward normals | upper face `n̂ = +ŵ` · lower face `n̂ = −ŵ` |
| response ratio | `Z ≡ (pressure perturbation at a face) / (OUTWARD normal velocity of that face)` |
| branch selection | energy or amplitude must not grow away from the slab; state the criterion you apply |

⛔ `x` is a position, ⛔ **never** a displacement. `Ω` (below) is a window, ⛔ never a thickness.

## 2 · The bulk sector

A scalar superfluid linearised to acoustics, with **no shear modulus**:

```
v = ∇₄ φ ,     δp = −ρ_m ∂_t φ ,     ∂_t² φ = c_s0² ∇₄² φ
```

`ρ_m` is the bulk mass density, `c_s0` the bulk sound speed. Both half-spaces `|w| > W₀/2` are bulk.

⛔ **Nothing about the slab's interior is specified and nothing about it is needed.** The slab enters this
step **only** through the motion of its faces.

## 3 · The window and the current

`Ω(w)` is a smooth **window function**, `≈1` inside the slab and `→0` outside, used to project
four-dimensional equations onto a three-dimensional description. Unless a task says otherwise take `Ω(w)`
**even** in `w`.

`j^w` is the `w`-component of the four-dimensional mass (or number) current.

⭐ **Relative flux.** Where material crosses a **moving** face, the physically meaningful flux is measured
**relative to that face and along its outward normal**:

```
J_± ≡ ρ_m ( v_w ∓ ∂_t ζ_± ) · (±1)        evaluated at the upper (+) / lower (−) face
```

⚠ Report explicitly which signed combination of `J₊` and `J₋` corresponds to **net accretion by the slab**
and which to **through-flow**. ⛔ Do not use a sum of global-`w` fluxes in place of relative ones; the two
differ, and they carry different physics.

## 4 · Background state

There is a steady background transfer of material across the interface. Treat the background as constant
over one wave period and work with the perturbation about it. ⭐ **Derive** the condition under which that
is valid, as an inequality between named timescales. ⛔ Do not assume it holds.

---

## TASKS

⛔⛔ **RULES THAT OVERRIDE EVERY TASK.**
1. **Every reported value must come from a computation in the script.** ⛔ A printed assertion is not a
   result.
2. ⛔ **If a task cannot be completed, emit its tag as `NOT_ESTABLISHED` and say which input is missing.**
   **A refusal is a valid and valuable output.**
3. ⛔ **Never silently choose a branch, a closure, or a convention.** If a choice is required and is not
   fixed above, introduce a **free symbol**, say so, and report the answer as a function of it.

---

**A1 · Projection identity.**
Integrate the four-dimensional continuity equation against `Ω(w)` over `w`, integrating by parts to
isolate the term carrying `j^w`. Report the resulting source term for a **finite** interval `[w₁, w₂]` and
for an **infinite** one.
⇒ `S11BA_PROJECTION_FINITE`, `S11BA_PROJECTION_INFINITE`

**A2 · Parity.**
With `Ω(w)` even, evaluate A1's source term for `j^w(w)` (a) **even** in `w`, (b) **odd** in `w`.
⭐ **Report the interval used and whether it is symmetric about `w = 0`**, and state for each result whether
it is **exact**, and **on what interval** — ⛔ an oddness argument does not by itself fix an integral over
an asymmetric interval.
⇒ `S11BA_PARITY_EVEN_JW`, `S11BA_PARITY_ODD_JW`, `S11BA_PARITY_INTERVAL`

**A3 · Dynamical window.**
Repeat A1 with `Ω = Ω(w; x, t)`. Enumerate **every** term present here and absent from A1.
⇒ `S11BA_DYNAMIC_WINDOW_EXTRA_TERMS`

**A4 · Bulk response to moving faces — impermeable.**
Solve §2 with both faces displaced as `ζ₊`, `ζ₋`, imposing that the bulk normal velocity at each face
equals that face's normal velocity. Treat both half-spaces. Report `Z` (as defined in §1) and its real and
imaginary parts, and any effective inertial loading per unit area, **separately for the `δW` and `h`
combinations** of §1.
⭐⭐ **There are THREE regimes** — the bulk normal wavenumber squared may be positive, negative, **or
zero**. Report all three, including the behaviour of every reported quantity **at** the zero (grazing)
case, where some may be singular. ⚠ Omitting that third case is a known prior defect of this corpus.
⇒ `S11BA_Z_IMPERMEABLE`, `S11BA_Z_BY_REGIME`, `S11BA_Z_BY_PARITY`, `S11BA_ADDED_MASS`,
  `S11BA_GRAZING_BEHAVIOUR`

**A5 · Bulk response — permeable faces.**
⭐ **Closure, fixed here so the family is finite:** take the relative flux of §3 to obey a **local,
ALGEBRAIC** (⛔ no time derivatives) linear law at each face,

```
J_±  =  Λ_p · δp|_face  +  Λ_ζ · ∂_t ζ_±
```

with `Λ_p`, `Λ_ζ` **free symbols**. Report their dimensions. Recompute A4 under this condition and report
the modified `Z` **as a function of `Λ_p` and `Λ_ζ`**. Then report, **for each of A4's three regimes and
each parity combination**, whether a dissipative part (real, in phase with velocity) is present and **on
which coefficient it depends**.
⛔ Do not select a value for either coefficient. ⛔ Do not report a single "the" permeable response.
⚠ **Record the restriction to an algebraic law as a stated scope limit** — a more general law would admit
time derivatives, and that generalisation is not attempted here.
⇒ `S11BA_PERMEABLE_COEFF_DIMS`, `S11BA_Z_PERMEABLE`, `S11BA_DISSIPATIVE_BY_REGIME_AND_PARITY`,
  `S11BA_CLOSURE_SCOPE_LIMIT`

**A6 · Dimensions.**
Derive, from the equations above and ⛔ **not** from any external table or registry, the `[L,T,M]`
exponents of: `Z`, the A4 inertial loading, `ρ_m`, `c_s0`, `Λ_p`, `Λ_ζ`, and A1's source term. Show the
route for each.
⚠ **A check reducing to an identity of the form `(X − 2Y) + 2Y == X` is worthless.** Label each route
**independent** or **definitional**.
⇒ `S11BA_DIM_<name>` per entry, plus `S11BA_DIM_ROUTE_KIND` per entry

**A7 · Controls.** ⛔ FORM controls; ⛔ do not substitute a coefficient rescaling.
- **A** — make `Ω(w)` **asymmetric** in `w`, keeping the interval symmetric; recompute A2.
- **B** — keep `Ω(w)` **even** and make the **interval** asymmetric; recompute A2.
⭐ Report which of the two moves A2's result, so that a parity selection rule is distinguished from a
domain-symmetry artifact. ⛔ Report what each control does, ⛔ not what it was expected to do.
⇒ `S11BA_CONTROL_WINDOW_PARITY`, `S11BA_CONTROL_INTERVAL_SYMMETRY`

**A8 · Validity condition.** Emit §4's derived condition.
⇒ `S11BA_BACKGROUND_VALIDITY_CONDITION`

---

## OUTPUT FORMAT

One line per tag, `TAG: value`. Values must be explicit expressions, not prose, wherever mathematical. End
with a single `VERDICT:` line reporting whether the script's own internal consistency checks contradicted
each other.
⚠ **`VERDICT: PASS` means only "my internal checks did not contradict each other."** ⛔ It is not a verdict
on the physics and must not be worded as one.
