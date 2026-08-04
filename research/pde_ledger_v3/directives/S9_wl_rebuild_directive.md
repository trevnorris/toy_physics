# S9 REBUILD — blind Mathematica audit of the brane transverse sector

You are building **one** self-contained Mathematica script. Everything you need is in this file.

**Deliverable (absolute path, overwrite it):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

Verify it by running it yourself:
`cd /var/projects/toy_physics && timeout 600 math -script research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

It must exit 0 and complete in well under 600 s. Iterate until it does.

---

## 0. ⭐⭐⭐ THE RULE THAT GOVERNS THIS BUILD — three clauses, non-negotiable

> **1. The script may PRINT COMPUTED OBJECTS. It may NOT STATE CONCLUSIONS.**
> Every emitted payload must be a **CAS object** — an expression, a solved root, a matrix, a list of
> roots, or a boolean returned by a **symbolic test**. ⛔ **Never emit a hand-authored sentence, word or
> label that describes a result.** `Print["TAG: PASS"]` is forbidden. `Print["TAG: NOT_APPLICABLE"]` is
> forbidden. `Print["TAG: " <> If[test, "SOME_WORD", "OTHER_WORD"]]` is forbidden — the branch words are
> hand-authored conclusions no matter how the test was computed. If you find yourself typing an answer into a
> string, stop — that is the exact defect this rebuild exists to remove.
>
> **2. EMIT BOTH OPERANDS AND THE RESIDUAL, then guard.**
> ⛔ A residual that is asserted zero always prints `0` and carries **no information**. Wherever you
> compare two things, emit **the independently computed object**, **the object it is compared against**,
> and **their difference** — three separate tags — and only *after* all three are printed may you apply
> a hard stop. **Compute → emit → then guard.** Never guard before emitting.
>
> **3. Interpretation belongs to the STEP RECORD, not to the script.**
> ⛔ The script does not editorialise, does not aggregate a verdict, and does not name what a value
> means. It prints objects; a human reads them.

**Why:** engine scripts in this corpus were found emitting physics conclusions as typed sentences that
no computation produced. Cross-engine agreement on such a tag is vacuous, and eight review legs missed
it, because *"does it say the right thing"* and *"does it depend on anything"* are different questions.

**Consequences you must honour:**
- ⛔ **No `VERDICT` tag. No `PASS`/`FAIL` tag. No aggregate check summary.** Emit each computed
  quantity on its own tag and let the reader compare them.
- ⛔ **No boolean emitted alone.** A boolean is admissible only as the value of a symbolic test, and
  only when the two operands and their residual are *also* emitted on their own tags.
- ⛔ **No tag whose value is a fixed string chosen by you.** The only strings permitted anywhere in
  output are the tag names themselves.
- ⭐ Prefer emitting a raw expression over emitting a test of it.
- ⛔⛔ **THE TAG NAME IS ALSO OUTPUT, AND IT MAY NOT STATE A CONCLUSION EITHER.** A tag name must name
  **the object**, never its value, its ratio, its sign, or the shape of the answer.
  `WL_S9_TRANSVERSE_ROOT` is fine. A name that appends the root's value, its ratio of coefficients, its
  sign, or a word for the shape of the dispersion is **forbidden** — the payload being a genuine CAS
  object does not rescue a name that already gave the answer away.
- ⛔⛔ **NO RESIDUAL MAY BE TAUTOLOGICAL.** If you define `q := A/B` and then emit `A − q·B`, that
  residual is identically zero **by construction** and carries no information. Before emitting any
  difference, ask: *were these two operands produced by independent routes?* If they were not, ⛔ do not
  emit the difference — emit the two objects and say in your report that no independent second route
  exists. ⭐ An honest "there is no second operand here" beats a zero that looks like a check.

---

## 1. Architectural constraints

- **Standalone.** No command-line arguments, no file reads, no imports, no exports, no network.
- ⛔ **This script must NOT read the reduction registry**
  (`research/pde_ledger_v3/reduction/quantities.yaml`, `relations.yaml`) — not directly, not by
  transcribing values out of it. That is a standing architectural rule: the Mathematica side is written
  blind of the registry so that the two sides can genuinely **disagree**. A disagreement is the most
  valuable output this build can produce.
- ⛔ **Do not read the predecessor version of the deliverable, and do not read any other `.wl` or `.py`
  under `research/pde_ledger_v3/`.** Build the script from this directive. If you have already read one,
  say so plainly in your final report.
- **Output discipline:** every line of stdout is exactly `TAG: <value>`, where `TAG` matches
  `[A-Za-z][A-Za-z0-9_]*`, prefixed `WL_S9_`. **Every tag must be unique** — a duplicate tag is a
  defect. Multi-line values are fine (continuation lines are attributed to the preceding tag). ⛔ No
  banner lines, no blank preamble, no untagged text before the first tag.
- Use `ToString[expr, InputForm]` so values are machine-readable.
- Set `$HistoryLength = 0` and `ClearAll["Global`*"]` at the top.

---

## 2. SUPPLIED PREMISES — true for this build, and ⛔ NOT falsifiable by it

State each of these in the script only as an input; the build cannot test any of them. They are
supplied because under-specification has cost this ledger far more than contamination ever has.

| # | premise |
|---|---|
| P1 | The medium has a **substructure** whose coarse-grained description is a GNLS superfluid. This is a charter-level standing postulate. |
| P2 | A single-component **scalar** superfluid carries no transverse (spin-1) excitation. Cited from external review; ⛔ **no script executes this and this build does not either.** |
| P3 | The **brane** is a `D`-dimensional sheet embedded in a higher-dimensional bulk. The registry records `D = 3`; this build must treat `D` symbolically wherever it can (see §4.5). |
| P4 | The **bulk carries no shear.** The bulk therefore contributes nothing to the action below. This is postulated, not derived. |
| P5 | The brane's transverse sector is governed by the **MacCullagh curl-only** stiffness, i.e. the stiffness depends on `u` only through `∇×u`. Forced by Maxwell's demand of *no longitudinal mode*; postulated as a **structural form**, not derived. §5 asks you to control it. |
| P6 | Linear response about an **unstrained brane at rest**: small oscillations, no background flow, no background strain, no dissipation, no relaxation time. |

---

## 3. THE PHYSICS, IN FULL

### 3.1 Field content

One real displacement field on the brane,

```
u(t, x) = (u₁(t,x,y,z), u₂(t,x,y,z), u₃(t,x,y,z))
```

interpreted as the brane's transverse displacement. Coordinates `{x, y, z}`; time `t`.

### 3.2 The action

The brane's transverse sector is the Lagrangian density

```
L  =  (1/2) ρ_br (∂_t u)·(∂_t u)  −  (1/2) μ_R (∇×u)·(∇×u)
```

with two positive constants:

- `ρ_br` — the brane's inertia density,
- `μ_R` — the MacCullagh (curl-only) shear rigidity.

`L` is an **energy density on the `D`-dimensional brane**.

⚠ `μ_R` is **not** the Cauchy shear modulus `μ_br`. They are different objects with near-identical
names. The Cauchy shear modulus is set to zero in this regime and does not appear here.

### 3.3 Equation of motion

Derive it. Apply the Euler–Lagrange equation in all four variables `{t, x, y, z}` to every component of
`u`:

```
∂L/∂u_i  −  Σ_α ∂_α ( ∂L / ∂(∂_α u_i) )  =  0            α ∈ {t, x, y, z}
```

⛔ Do not write the equation of motion down by hand and check it. **Vary the Lagrangian.**

### 3.4 Plane-wave reduction

Substitute

```
u(t,x) = a exp( i ( k·x − ω t ) )         a = (a₁, a₂, a₃),  k = (k_x, k_y, k_z)
```

and divide out the common phase. The result is linear in `a`, so it defines a **dynamical matrix**
`M(ω, k)` through `M · a = 0`.

### 3.5 The dispersion relation and mode selection

Solve `det M = 0` for `ω²`, treating `ω²` as a single unknown. That gives a set of roots.

**Selecting roots by polarisation — do this by computation, ⛔ never by inspection.** Define **two**
generators, and note that neither divides by `|k|²`:

```
T  =  |k|² · I  −  k kᵀ            (spans the TRANSVERSE subspace)
Λ  =  k kᵀ                          (spans the LONGITUDINAL subspace)
```

A root `ω²` is **transverse** exactly when `M(ω², k) · T = 0` as a matrix identity, and **longitudinal**
exactly when `M(ω², k) · Λ = 0`. Test every root against **both**, and emit both product matrices for
every root.

⚠ These two classes are **not** exclusive and **not** exhaustive, and you must not assume they are:
- a root may satisfy **both** (that is what happens when `M` vanishes identically at that root),
- a root may satisfy **neither** (that is what happens when no single `ω²` annihilates a whole subspace,
  as can occur once the inertia is anisotropic).

⛔ **Do not "fix" either case.** Emit the matrices and the resulting subsets exactly as computed,
including an empty subset. An empty list is a legitimate computed result; ⛔ never substitute a word for
it and never emit a hand-picked list of roots. The subsets must be produced by `Select` over the full
root list using the programmed matrix tests — ⛔ not by an `If` you wrote knowing the answer.

### 3.6 Dimensions

Use the exponent-vector convention **`[L, T, M]`**, ordered length, time, mass. So `[1,0,0]` is a
length, `[0,1,0]` a time, `[0,0,1]` a mass.

Build the derivation from first principles, not from a table:

```
acceleration    = length − 2·time
force           = mass + acceleration
energy          = force + length
energy density on a D-dimensional sheet  =  energy − D·length
```

`u` is a displacement, so `[u] = [1,0,0]`. Then `[∂_t u] = [u] − time` and `[∇×u] = [u] − length`.

Now **solve** for the two unknown dimension vectors by requiring each term of `L` to be an energy
density on the `D`-dimensional sheet:

```
[ρ_br] + 2[∂_t u]  =  energy density
[μ_R]  + 2[∇×u]    =  energy density
```

Treat `[ρ_br]` and `[μ_R]` as six unknown exponents and solve the linear system. ⛔ Do not type the
answer and verify it.

---

## 4. THE QUESTION LIST — what to compute

⛔ **This section says what to compute. It does not say what any of it equals.** No expected value
appears anywhere in this directive, deliberately. Your job **ends at compute-and-print**: the
comparison against the recorded ledger happens on the orchestrator's side afterwards, where a
disagreement is a **finding**, not a build failure. ⛔ **Do not tune anything to make a number come
out.** If a computation produces something that looks wrong to you, print it anyway and say so in your
report.

Emit each of the following as its own tag.

### 4.1 The derivation chain, printed at every stage
1. The Lagrangian density as constructed.
2. The Euler–Lagrange residual vector, before simplification of sign convention.
3. The equation of motion (the residual set equal to zero, componentwise).
4. The plane-wave ansatz as substituted.
5. The plane-wave residual after dividing out the phase.
6. The dynamical matrix `M`.
7. `det M`, factored.
8. The **complete** solution set for `ω²` — every root, not only the one you want.
9. Both generators `T` and `Λ`.
10. For **each** root: the matrix `M(ω²,k)·T`, the matrix `M(ω²,k)·Λ`, and the root itself.
    ⭐ Emit the matrices, not a verdict on them.
11. The subset of roots passing the transverse test, and the subset passing the longitudinal test —
    each produced by `Select` with the programmed test. Emit an empty subset as an empty list.
12. ⭐ **The complete root multiset with multiplicities**, taken from the factored determinant, **and the
    nullity of `M` at each root** (`Length[k] − MatrixRank[M]` evaluated at that root). ⭐ This is the
    object that stays well-posed when roots collapse or subsets come out empty, so emit it for the main
    derivation **and for every control.**

### 4.2 The transverse dispersion

Work with a **scalar** wavenumber for these tags: introduce `kMag` with `kMag² = k·k`, and express
`ω²` in terms of `kMag` so that differentiation with respect to a scalar is well defined.

13. `ω²` for the transverse root, in terms of `kMag`.
14. The candidate propagation speed squared, obtained as `ω²/kMag²`.
    ⛔ **Do NOT then emit `ω² − (that speed squared)·kMag²`.** That residual is identically zero **by
    construction** — the speed was *defined* as that quotient — so it is a tautology, not a check, and
    this rebuild exists to remove exactly that pattern. ⚠ At this step there is **no independent second
    route** to the propagation speed: the script may not read the registry (§1), and there is only one
    engine. ⭐ Say so in your report rather than manufacturing a second operand.
15. The speed-squared expression's dependence on `kMag`, emitted **as expressions**:
    `D[speedSquared, kMag]`, and `Exponent[Together[omegaSquared], kMag]` for numerator and denominator.
    ⭐ **These carry the actual structural information** about the dispersion. ⛔ Never emit a word
    describing its shape.

### 4.3 Dimensions
16. The energy-density dimension vector on a `D`-dimensional sheet, symbolic in `D`.
17. The linear system solved for the two unknown dimension vectors, printed as equations.
18. `[ρ_br]` as solved, symbolic in `D`.
19. `[μ_R]` as solved, symbolic in `D`.
20. `[μ_R] − [ρ_br]`, symbolic in `D`. ⭐ Emit the vector exactly as the solver returns it.
21. The dimension vector implied for the propagation speed squared, obtained by reading the powers of
    `ρ_br` and `μ_R` **out of the computed speed-squared expression** and combining the two solved
    dimension vectors accordingly.
22. The dimension vector of a squared velocity, built from the primitives in §3.6.
23. The **difference** between tags 21 and 22.
24. Tags 18–20 evaluated at the registry's brane dimension `D = 3`.

---

## 5. ⭐⭐ REQUIRED CONTROLS — and most of them must be FORM controls

⚠ **A COEFFICIENT control tests the arithmetic; only a FORM control tests the physics**, because
rescaling a coefficient never leaves the family of theories. The predecessor of this script had only
coefficient controls, and its own recorded limitations admit it was blind to the stiffness *form*, to a
flipped potential sign, and to the assumed brane dimension. **Fix that here.**

Structure the derivation as a **function of the action** so a control is one call with a different
Lagrangian, and the entire chain — variation, plane wave, determinant, roots, transverse selection —
re-runs. ⛔ Do not hand-write a control's expected outcome.

For **each** control below, emit exactly these five objects, and ⛔ nothing that interprets them:

1. the control's Lagrangian as constructed;
2. its factored determinant;
3. ⭐ **its complete root multiset with multiplicities, and the nullity of `M` at each root** — this is
   the primary comparison object, because it stays well defined when roots collapse into one, when a
   root disappears, and when a polarisation subset comes out empty;
4. its transverse subset and its longitudinal subset, each from the programmed matrix tests of §3.5,
   emitted as lists — ⛔ **an empty list is emitted as an empty list**;
5. the difference between its root multiset and the main derivation's root multiset, as CAS objects.

⛔⛔ **Do NOT emit "the difference between the transverse root and the non-transverse root."** That
phrasing is ill-posed and an earlier revision of this directive wrongly demanded it: under some controls
the root set **collapses to a single root**, and under others a subset is **empty**, so there is no
second operand and a builder would have to **invent** one. ⭐ The root multiset in (3) is the well-posed
replacement and it discriminates every case.

| id | kind | change |
|---|---|---|
| **X1** | coefficient | `ρ_br → λ_ρ ρ_br` with `λ_ρ` a positive symbol |
| **X2** | coefficient | `μ_R → λ_μ μ_R` with `λ_μ` a positive symbol |
| **X3** | ⭐ **FORM** | replace the curl-only stiffness by a **divergence-only** stiffness: `−(1/2) μ_R (∇·u)²` |
| **X4** | ⭐ **FORM** | replace the curl-only stiffness by **isotropic gradient elasticity**: `−(1/2) μ_R Σ_{i,j} (∂_i u_j)(∂_i u_j)` |
| **X5** | ⭐ **FORM** | flip the sign of the curl-only stiffness term: `+(1/2) μ_R (∇×u)·(∇×u)` |
| **X6** | ⭐ **FORM** | make the inertia anisotropic: replace `ρ_br (∂_t u)·(∂_t u)` by `(∂_t u)·diag(ρ_br, ρ_br, ρ_z)·(∂_t u)` with `ρ_z` an independent positive symbol |

For **X4** additionally emit the difference between its transverse root and its **non**-transverse root,
as an expression — that difference is the object that distinguishes a curl-only medium from an ordinary
gradient-elastic one.

For **X6** additionally emit the full multiset of roots so degeneracy (or its loss) is visible as data.

⛔ Do not emit a `FIRED` boolean for any control. Emit the control's computed objects and the
differences; the reader decides whether it fired.

---

## 6. ⭐ STATIC-OR-INSTANTANEOUS — state this in the script's own emitted output? NO. State it in your report.

The following limits are taken in this build. They are recorded here so a passing build is not read as
having tested them. ⛔ Do not emit these as tags — they are prose, and prose does not go in output.

- The brane's **wall width → 0** (a sharp `D`-dimensional sheet). This removes the flexural channel;
  a finite width would add a higher-order correction to the dispersion. Whether it does is a later
  step's question, not this one's.
- The background flow **v₀ → 0** and the background strain → 0. This removes all convective terms.
- **No relaxation time and no dissipation** — the action is conservative by construction.
- The bulk is **absent from the action entirely** (P4).
- **Continuum / homogenisation**: the substructure's microscopic length over the wavelength → 0, so
  `ρ_br` and `μ_R` are local constants. ⚠ This looks like an ordinary constitutive law, not like
  removing a scale.
- **The moduli are frequency-independent** — `ρ_br` and `μ_R` carry no memory kernel. This is stronger
  than "no dissipation" and is a separate assumption.
- **Amplitude → 0** (strict linearisation); every term is quadratic in `u`.

---

## 7. What you must report back

Keep it under 40 lines.

1. The absolute path of the deliverable and the exit status of your own run of it.
2. The **number of tags** emitted and confirmation that all are unique.
3. A statement, for each of the three clauses in §0, of how the script complies — with a **line number**
   for the mechanism, not a paraphrase.
4. **For every emitted tag, the line number of the computation that produced its value.** If any tag's
   value is not produced by a computation, say which one and why it exists.
5. Anything in §3 or §4 you could not compute, and what blocked it.
6. Anything you computed that surprised you. ⭐ This is wanted; it is not a failure report.

⛔ Do not commit anything to git. ⛔ Do not modify any file other than the deliverable.
