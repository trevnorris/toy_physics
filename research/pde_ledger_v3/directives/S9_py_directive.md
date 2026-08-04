# BUILD — the S9 SymPy engine

You are building a **SymPy engine script** for step S9 of a physics ledger. Everything you need is here.

**Deliverables (absolute paths):**
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py`
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.premises`

Verify by running it yourself:
`cd /var/projects/toy_physics/research/pde_ledger_v3 && timeout 600 python3 scripts/S9_light_requires_shear_sympy_audit.py`

Exit 0, well under 600 s. Iterate until it does.

---

## 0. ⭐⭐⭐ THE RULE THAT GOVERNS THIS BUILD

> **1. The script PRINTS COMPUTED OBJECTS. It does NOT STATE CONCLUSIONS.**
> Every emitted payload is a **CAS object** — an expression, a solved root, a matrix, a list, an integer
> a computation returned, or a boolean from a **symbolic test**. ⛔ **Never a hand-authored word
> describing a result.** No `print("...: PASS")`. No `"LINEAR"`. No `f"...: {'YES' if x else 'NO'}"`.
> The branch words are hand-authored conclusions no matter how the test was computed.
>
> **2. EMIT BOTH OPERANDS AND THE RESIDUAL, then guard.** Where two **independently produced** objects
> should agree, emit both and their difference as three tags, **then** assert. ⛔ Never assert before
> emitting: a residual that is asserted always prints `0` and carries no information.
>
> **3. Interpretation belongs to the STEP RECORD.** ⛔ No verdict tag, no aggregate summary, no
> editorialising. The script prints objects; a human and an automated consumer read them.

**Two corollaries that are part of the rule:**

- ⛔⛔ **THE TAG NAME IS OUTPUT TOO.** A tag names **the object**, never its value, ratio, sign, or the
  shape of the answer. `PY_S9_TRANSVERSE_ROOT` is fine; appending the answer to the name is forbidden,
  and a genuine CAS payload does not rescue a name that gave the answer away.
- ⛔⛔ **NO TAUTOLOGICAL RESIDUAL.** If you define `q := A/B` and emit `A − q*B`, that is zero **by
  construction** and certifies nothing — it vanishes for *any* input. Before emitting a difference, ask:
  *were these two operands produced by independent routes?* If not, ⛔ do not emit the difference; emit
  the objects and say in your report that no second route exists.

⭐⭐ **THE STRUCTURAL RULE — obey it literally:**
> **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTIONS (§3.2,
> §5) and the plane-wave ansatz (§3.4). Every other expression involving them must be REACHED BY
> COMPUTATION from those.**

⇒ ⛔ Do not type a root, a dispersion, a speed, a dynamical matrix, or a dimension vector that a solver
should have returned. Each must arrive via `solve`, `det`, `rank`, `diff`, `Poly`, or the Euler–Lagrange
variation. ⭐ **Every control re-enters the chain at its ACTION**, ⛔ never at a result.

⚠ **Why:** engine scripts in this corpus were found emitting physics conclusions that no computation
produced — first as typed sentences, then, subtler, as **hand-typed algebraic expressions** that pass any
"is it a CAS object" test while having no data dependency on the derivation. Ablating the real
computation left the output byte-identical.

---

## 1. Architecture

- **Standalone.** No command-line arguments, no input files, no network.
- ⛔⛔ **DO NOT READ THE REDUCTION REGISTRY** (`reduction/quantities.yaml`, `relations.yaml`) and ⛔ **do
  not import from `reduction/`.** Comparison against the registry happens **outside** this script, in a
  consumer, afterwards. If this script read the registry it could report the registry's value instead of
  computing one.
- ⛔⛔ **DO NOT READ, RUN, OR OPEN THE MATHEMATICA ENGINE**
  `research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`, nor any `.wl` file,
  nor any other engine under `scripts/`. ⭐⭐ **A second engine is worth having ONLY if it was
  independently constructed, because the entire point is that the two can DISAGREE.** A transcription
  agrees vacuously and is worse than nothing, because it looks like confirmation.
  ⚠ If you open one anyway, ⭐ **say so plainly in your report** — that is recoverable; a silent
  transcription is not.
- **Output contract**, enforced by an automated consumer:
  - every stdout line is `TAG: value`, tag matching `[A-Za-z][A-Za-z0-9_]*`, prefixed **`PY_S9_`**;
  - ⛔ **every tag unique** — a duplicate tag breaks the consumer and is a defect;
  - ⛔ **no untagged output before the first tag**, no banners, no blank preamble, no logging;
  - render values with `sympy.srepr`-free plain `str()`/`sstr()` so they are re-parseable text.
- Also write the sidecar `...sympy_audit.premises`: one exact tag name per line, listing tags whose value
  is a **supplied premise** rather than a derived result (`#` comments and blank lines allowed). If none,
  create the file with only a comment.

---

## 2. SUPPLIED PREMISES — ⛔ NOT falsifiable by this build

| # | premise |
|---|---|
| P1 | The medium has a **substructure** whose coarse-grained description is a GNLS superfluid. Charter-level standing postulate. |
| P2 | A single-component **scalar** superfluid carries no transverse (spin-1) excitation. Cited from external review; ⛔ **no script executes this and this one does not either.** |
| P3 | The **brane** is a `D`-dimensional sheet in a higher-dimensional bulk. Treat `D` **symbolically** wherever possible; the registry records `D = 3`. |
| P4 | The **bulk carries no shear**, so it contributes nothing to the action below. Postulated. |
| P5 | The brane's transverse sector has **MacCullagh curl-only** stiffness — the stiffness depends on `u` only through `∇×u`. Postulated as a structural **form**; §5 controls it. |
| P6 | Linear response about an **unstrained brane at rest**: small oscillations, no background flow or strain, no dissipation, no relaxation time. |

---

## 3. THE PHYSICS, IN FULL

### 3.1 Field content
Three real displacement components on the brane, `u = (u1, u2, u3)`, each a function of `(t, x, y, z)`.

### 3.2 The action
```
L  =  (1/2) ρ_br (∂_t u)·(∂_t u)  −  (1/2) μ_R (∇×u)·(∇×u)
```
`ρ_br` the inertia density, `μ_R` the MacCullagh (curl-only) shear rigidity, both positive. `L` is an
**energy density on a `D`-dimensional brane**.
⚠ `μ_R` is **not** the Cauchy shear modulus `μ_br`, which is zero in this regime and absent here.

### 3.3 Equation of motion — ⛔ derive it, do not write it down
```
∂L/∂u_i  −  Σ_α ∂_α ( ∂L / ∂(∂_α u_i) )  =  0        α ∈ {t, x, y, z}
```

### 3.4 Plane-wave reduction
Substitute `u = a·exp(i(k·x − ω t))` with `a = (a1,a2,a3)`, `k = (kx,ky,kz)`, divide out the common
phase. The result is linear in `a` and defines the **dynamical matrix** `M` via `M·a = 0`.

### 3.5 Polarisation — ⭐ four DIFFERENT questions, ⛔ do not conflate them
With `q ≡ k·k`, and assuming `q ≠ 0` (impose this, do not merely comment it):
```
T = q·I − k kᵀ      (spans the TRANSVERSE subspace)
Λ = k kᵀ            (spans the LONGITUDINAL subspace)
```

| test | question | computation |
|---|---|---|
| **E1** transverse **existence** | is there `a ≠ 0` with `M a = 0` **and** `k·a = 0`? | stack `M` on the row `kᵀ`; `rank(stacked) < 3` |
| **E2** transverse **multiplicity** | how many independent transverse modes? | `3 − rank(stacked)` |
| **E3** longitudinal existence | is `k` a null vector of `M`? | `M·k == 0` |
| **E4** whole-subspace degeneracy | does `M` annihilate **all** of the transverse subspace at one `ω²`? | `M·T == 0` |

⛔⛔ **`M·T = 0` IS NOT THE EXISTENCE TEST.** It demands `M` annihilate the *entire* transverse subspace,
so a theory carrying a genuine transverse mode whose partner sits at a **different** frequency reports
**no transverse mode** — a false negative. ⭐ **E2 is the mode count**; emit it for every root of every
action. ⚠ The classes are neither exclusive nor exhaustive: a root may satisfy both, or neither, and E1
may hold while E4 fails.

### 3.6 Dimensions — ⭐ derive the derivative structure FROM the action
Convention `[L, T, M]`, ordered length, time, mass. Build from primitives:
```
acceleration = length − 2·time;  force = mass + acceleration;  energy = force + length
energy density on a D-sheet = energy − D·length
[u] = [1,0,0]
```
⭐⭐ **Read the derivative structure out of the constructed Lagrangian.** For each additive term, and each
field factor in it, extract the derivative multi-order by inspecting the expression tree (`Derivative`
nodes). The dimension of a term is
```
[coefficient] + Σ over field factors of ( [u] − dt·[T] − (dx+dy+dz)·[L] )
```
⛔ **Do not hand-write the derivative counts and ⛔ do not special-case curl or divergence by name** — the
counts must come from the tree, so that changing the action changes them. Require every term to equal the
energy density on a `D`-sheet and **solve** the linear system for the unknown coefficient dimensions.

---

## 4. THE QUESTION LIST

⛔ **This says what to compute. It does not say what anything equals** — no expected value appears
anywhere in this directive, deliberately. Your job **ends at compute-and-print**. The comparison against
the other engine and against the ledger happens afterwards, elsewhere, where a mismatch is a **finding**,
⛔ not a build failure. ⛔ **Do not tune anything to make a number come out.** If something looks wrong to
you, print it anyway and say so in your report.

Emit each as its own tag:

1. the Lagrangian as constructed; 2. the Euler–Lagrange residual vector; 3. the equation of motion;
4. the plane-wave ansatz as substituted; 5. the plane-wave residual after dividing the phase;
6. `M`; 7. `det M`, factored; 8. the complete solution set for `ω²`; 9. both generators `T` and `Λ`;
10. per root: the stacked matrix of E1, the vector `M·k`, the matrix `M·T`;
11. per root: **E1, E2, E3, E4** as computed, and the subsets of roots passing each, via a **programmed**
    selection over the full root list — ⛔ an empty subset is an empty list, never a word;
12. ⭐ the **complete root multiset with multiplicities** from the factored determinant, **and the nullity
    of `M` at each root** — a multiset stays well-posed whatever the root structure turns out to be;
13. for every root passing E1: `ω²` in terms of the scalar `q`;
14. for each such root, the candidate speed squared `ω²/q`.
    ⛔ **Do NOT emit `ω² − (that)·q`** — zero by construction, and it vanishes for a dispersive `ω²` too.
15. ⭐ instead emit the **non-circular** structural objects, per such root:
    `q*diff(ω², q) − ω²` (Euler homogeneity defect) and `diff(ω²/q, q)`, plus the degree of numerator and
    denominator of `together(ω²)` in `q`. ⭐ These vanish exactly when `ω²` is homogeneous of degree one
    in `q` and are **nonzero** for a dispersive one, so they carry real information.
16. the energy-density dimension vector, symbolic in `D`;
17. the extracted per-term derivative multi-orders and per-term dimension expressions;
18. the dimension linear system and its solution; `[ρ_br]`; `[μ_R]`; `[μ_R] − [ρ_br]`;
19. the dimension implied for the speed squared, read from the **computed** speed expression's powers of
    `ρ_br` and `μ_R`; the dimension of a squared velocity built from primitives; and their difference;
20. items 18 evaluated at `D = 3`;
21. the assumption set you imposed, as a CAS object.

### 4.1 ⭐⭐ Give `M` an independent SECOND ROUTE
`M` currently arrives one way, and nothing checks it — a sign error or dropped term there corrupts every
root, test and dimension downstream with nothing moving.

- **Route A:** Euler–Lagrange in position space → equation of motion → substitute the plane wave → read
  `M` off the amplitude coefficients.
- **Route B:** substitute the ansatz **directly into the Lagrangian**, reduce to a quadratic form in the
  amplitudes, and obtain `M` from that form by differentiating with respect to amplitude components.
  ⛔ Reuse nothing from Route A.

Emit `M_A`, `M_B`, and the residual `M_A − M_B` as three tags — ⭐ genuinely independent routes, so this
residual is **not** tautological. Emit all three, **then** guard. Do this for the main action and every
control. ⚠ If the two routes disagree by an overall factor or sign, ⛔ **do not tune** — emit the residual
and report it. A disagreement is a finding.

---

## 5. CONTROLS — ⭐ mostly FORM controls, and each re-enters at the ACTION

⚠ **A COEFFICIENT control tests arithmetic; only a FORM control tests physics**, because rescaling never
leaves the family of theories.

Structure the derivation as a **function of the action**, so a control is one call with a different
Lagrangian and the whole chain re-runs.

| id | kind | change |
|---|---|---|
| **X1** | coefficient | `ρ_br → λ_ρ ρ_br`, `λ_ρ` a positive symbol |
| **X2** | coefficient | `μ_R → λ_μ μ_R`, `λ_μ` a positive symbol |
| **X3** | ⭐ FORM | stiffness → **divergence-only**: `−(1/2) μ_R (∇·u)²` |
| **X4** | ⭐ FORM | stiffness → **isotropic gradient elasticity**: `−(1/2) μ_R Σ_{i,j} (∂_i u_j)²` |
| **X5** | sign | flip the stiffness sign: `+(1/2) μ_R (∇×u)·(∇×u)`. ⚠ Same operator, so ⛔ not a form change — it leaves the admitted family only because that family needs `μ_R > 0` |
| **X6** | ⭐ FORM | inertia → anisotropic: `(∂_t u)·diag(ρ_br, ρ_br, ρ_z)·(∂_t u)`, `ρ_z` an independent positive symbol |
| **X7** | ⭐ FORM | stiffness → **flexural**, a different derivative count: `−(1/2) μ_F (∇²u)·(∇²u)`, `μ_F` a fresh positive symbol |

⭐⭐ **Every control emits the SAME tag set as the main derivation**, under its own prefix — the upstream
chain, both routes to `M`, all four polarisation tests per root, the dispersion and homogeneity objects,
**and its own dimension block**. An automated consumer pairs each main tag with its control counterpart
and asks whether it moved; a tag with no counterpart cannot be asked.

⛔ Do not emit a `FIRED` boolean or any control verdict. ⛔ Where a control genuinely cannot produce a
quantity (a speed when no root passes E1), emit the **computed empty object**, never a placeholder.

### 5.1 Special directions
`E2` is a mode count that a later step's whole result rests on, and a symbolic `rank` returns the
**generic** value. For the **main** action and for **X6**, re-run the per-root E1–E4 with the wavevector
specialised, each as its own tag:
```
(kx,ky,kz) generic   ·   (kx,ky,0)   ·   (0,0,kz)   ·   (kx,0,0)   ·   (kx,kx,kx)
```
⛔ Emit the computed values for every case. ⛔ Do not summarise them and ⛔ do not say which are special.

---

## 6. Limits taken — ⛔ report these, do NOT emit them as tags (they are prose)

Sharp zero-width sheet · zero background flow and strain · no dissipation and no relaxation time ·
bulk absent from the action · continuum/homogenisation, so the moduli are local constants ·
frequency-independent moduli · amplitude → 0.

⛔⛔ **AND THE ONE THAT MATTERS:** a general stiffness would carry a compressional term
`−(1/2) B_L (∇·u)²`. Premise **P5** sets `B_L = 0`, sending the longitudinal restoring frequency to zero
and its timescale to infinity. ⚠ **That rate is exactly what decides whether a longitudinal wave
propagates**, which is part of what S9 claims. ⇒ ⭐ **any absence of a propagating longitudinal wave here
is ASSUMED THROUGH P5, ⛔ NOT ESTABLISHED**, and P5 removes the longitudinal *restoring force*, ⛔ not the
longitudinal *degree of freedom* — which is why §3.5 asks for E3 at every root.

---

## 7. Report back — under 40 lines

1. Deliverable paths; exit status and runtime of your own run; tag count and confirmation all are unique.
2. Confirmation you did not open the `.wl`, the registry, or any other engine. ⭐ If you did, say so.
3. For each of the three clauses, the **line number** of the mechanism that enforces it.
4. **For every emitted tag, the line number of the computation that produced its value.** If any tag's
   value is not produced by a computation, name it and say why it exists.
5. The literal `M_A − M_B` residual for the main action.
6. Anything in §3–§5 you could not compute, and what blocked it.
7. ⭐ Anything you computed that surprised you. This is wanted.

⛔ Do not commit to git. ⛔ Do not modify any file other than your two deliverables.
