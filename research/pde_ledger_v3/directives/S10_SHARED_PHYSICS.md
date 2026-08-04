# S10 — SHARED PHYSICS

⛔⛔ **This file is read by BOTH engines. An error here defeats dual-engine agreement by construction:
both engines will compute the same wrong thing and agree.** It is reviewed as its own artifact before
either build starts.

⭐ It supplies **the setup, the field content, the governing equations, the premises, and the list of
quantities to compute.** ⛔ It supplies **no result** — no root, no count, no dimension tuple, no sign.
The engines produce those.

---

## 1 · The physical setup

A **brane** is a `D`-dimensional elastic sheet. `D` is a free integer symbol throughout; ⛔ do not fix it.

The field is `u`, the **in-plane displacement**: a `D`-component vector field on the brane,
`u_j(x_1, …, x_D, t)`, `j = 1…D`.

⚠ **SUPPLIED PREMISE — unfalsifiable within this build. State it in the output as a premise tag.**
`u` is **in-plane only**. Motion *out of* the brane is a **different field** (the `h`-branon) belonging to
a different sector, and it is **not** part of `u`. ⇒ `u` has exactly `D` components, ⛔ not `D+1`.
⭐ Nothing in this build tests that split; it is inherited.

## 2 · The action — supplied in full, as equations

```
L  =  (ρ_br / 2) · Σ_j (∂_t u_j)²   −   (μ_R / 2) · S[∂u]
```

where `S[∂u]` is the **stiffness density**, and the baseline stiffness is the **antisymmetric-derivative
(curl-only)** form, defined in every `D`:

```
S_curl[∂u]  =  (1/2) · Σ_{i=1..D} Σ_{j=1..D} ( ∂_i u_j  −  ∂_j u_i )²
```

`ρ_br > 0` and `μ_R > 0` are real constants.

⛔⛔ **DO NOT SIMPLIFY `S_curl` BY HAND, IN EITHER COORDINATE SPACE OR AFTER SUBSTITUTING THE ANSATZ.**
Whatever `S_curl` becomes on a plane wave is a **computed object** and is one of the things this build
exists to produce. ⛔ Typing its reduced form is the defect this whole rebuild removes.

## 3 · The ansatz — supplied

A plane wave with **constant** amplitude `a` and **real** wavevector `k`:

```
u_j(x, t)  =  a_j · exp( i ( Σ_m k_m x_m  −  ω t ) )
```

`a` and `k` are `D`-component symbolic vectors. Write `ω²` as a **single symbol** (`omegaSquared`), ⛔ never
as `ω` squared — the spectrum is polynomial in `ω²` and solving for `ω` introduces spurious branches.

**Assumptions in force everywhere:** `ρ_br > 0`, `μ_R > 0`, `Σ_m k_m² > 0`, and every component of `k` and
`a` real.

⚠ **SUPPLIED PREMISES — unfalsifiable within this build. Emit each as a premise tag:**
- the background brane is **unstrained and at rest**: `v₀ = 0`, so there are **no convective terms**;
- there is **no dissipation** — no relaxation time appears in `L`;
- **linear response / small oscillations**: `L` is exactly quadratic in `u`.

⭐ Why these are stated rather than assumed silently: a background drift picks out a direction, and any
statement this build makes about **degeneracy between directions** inherits `v₀ = 0`. ⛔ Do not present a
degeneracy as a property of the medium; it is a property of the medium **on a state at rest**.

## 4 · ⭐⭐⭐ THE STRUCTURAL RULE — verbatim, non-negotiable

> **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
> ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
> the chain at the ACTION, ⛔ never at a result.**

⛔⛔ **A hand-typed CAS object is still hand-typed.** `emit(FullSimplify[{0, muR k2/rhoBr}])` is a genuine
CAS object with **no data dependency on the derivation** — delete `Det` and `Solve` and the output does not
move. ⭐ **The test: if you deleted the computation, would this tag change?** If not, it is dead.

## 5 · ⭐⭐⭐ THE THREE CLAUSES

> **1. A script may PRINT computed objects. It may NOT state conclusions.** Every `emit`/`Print` payload
> must be a CAS object — an expression, a solved root, a matrix, an integer produced by a rank
> computation, a boolean returned by a symbolic test. ⛔ Never author prose describing a result.
>
> **2. EMIT BOTH OPERANDS AND THE RESIDUAL, then guard.** ⛔ A residual that is asserted zero always
> prints `0` and carries no information. Emit operand A, operand B, **and** `A − B`. Emit **first**, then
> hard-stop if you must.
>
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

### The four corollaries

1. ⛔ **A hand-typed CAS object is still hand-typed** — see §4.
2. ⛔ **The tag NAME is output too.** A name may name **the object**; ⛔ never its value, count, ratio,
   sign, or the shape of the answer. ⚠ Schematically — with **placeholder** content, ⛔ not this step's:
   `..._EQUALS_SEVEN`, `..._IS_POSITIVE`, `..._QUADRATIC_IN_Z` are forbidden **name** patterns whatever the
   payload, because the name has already answered. ⭐ `..._RANK`, `..._ROOTS`, `..._RESIDUAL` are fine —
   they name the object and leave the value to the payload.
3. ⛔ **No tautological residual.** Before emitting a difference, ask: were these two operands produced by
   **independent routes**? If `q := A/B` and you emit `A − q·B`, it is zero by construction and vanishes
   for **any** input, including a wrong one. ⭐ Where no second route exists, emit the objects and say
   there is nothing to compare against.
4. ⛔⛔ **EMISSION MUST NEVER BE CONDITIONAL ON A PAYLOAD'S VALUE.** Whether a tag appears may depend only
   on **which package** and **which quantity** it belongs to. ⚠ A value **present and identical across
   packages is a RESULT**; a value **absent** is indistinguishable from *never computed*. ⛔ Never
   de-duplicate payloads across packages. Tag **names** are unique; **payloads may legitimately repeat.**

### ⛔⛔ NO VERDICT

⛔ **There is no `VERDICT` tag, no `PASS`, no `FAIL`, and no summary judgement anywhere in either engine.**
⛔ **A physics finding must EXIT 0.** Exit non-zero **only** on operational failure — an exception, an
unsolvable system, a timeout. ⚠ Otherwise a builder iterating to exit 0 can make a genuine disagreement
disappear, and the disagreement is the most valuable output available.

---

## 6 · What to compute — the question list

⛔ **This section says what to produce. It does not say what any of it equals.** Everything below is
required for **every package** defined in §7.

### Q1 · The Lagrangian and the equations of motion

Build `L` in coordinate space from §2. Emit `L` expanded.
Take the **Euler–Lagrange** variation with respect to each `u_j`:

```
∂_t ( ∂L / ∂(∂_t u_j) )  +  Σ_i ∂_i ( ∂L / ∂(∂_i u_j) )  −  ∂L / ∂u_j  =  0
```

Emit the full system.

### Q2 · ⭐⭐ The dynamical matrix, by TWO INDEPENDENT ROUTES

`M` is the `D×D` coefficient matrix of the plane-wave amplitude equations, `M · a = 0`.

- **Route A — through the equations of motion.** Substitute the §3 ansatz into the Q1 Euler–Lagrange
  system, strip the common exponential, and extract the coefficient matrix. Call it `M_A`.
- **Route B — through the quadratic form.** Form the plane-wave quadratic density directly from the §2
  action (potential minus `ρ_br ω²/2 · Σ a_j²`) and take
  `M_B[[i,j]] = ∂²(that density)/∂a_i ∂a_j`. Call it `M_B`.

⭐ **Emit `M_A`, `M_B`, and `M_A − M_B`** — all three (clause 2). ⛔ Do not derive one from the other, and
⛔ do not define one as a relabelling of the other. They must reach `M` by genuinely different operations.

⚠ If `M_A` and `M_B` differ by an overall nonzero scalar or a sign convention, ⛔ **do not silently
normalise one to the other.** Emit both, emit the difference, and emit the ratio `M_A[[1,1]]/M_B[[1,1]]`
as an object. The orchestrator resolves conventions; the engine reports.

⭐ Use `M_B` as the matrix for everything downstream, and say so in a tag naming the route used.

### Q3 · The spectrum

- Emit `Det[M]`, **factored**.
- Solve `Det[M] = 0` for `omegaSquared`. Emit the **full solution set**, then the **distinct** roots after
  de-duplication under the assumptions. Emit **how many distinct roots** there are, as an integer produced
  by counting the list.
- ⭐ For each distinct root, emit its **SIGN** under the assumptions (`Sign[root]`, or a three-way
  positive/zero/negative symbolic test). ⚠ Without this an exponentially growing mode carries the same
  signature as a propagating wave.
- ⭐ Emit the **root-coincidence locus**: for every pair of distinct roots `r_p`, `r_q`, solve
  `r_p − r_q = 0` and emit the solution set. ⚠ This is where the generic-rank results below could fail.

### Q4 · ⭐⭐⭐ THE MODE COUNT — basis-independent tests, and this is the heart of the step

For **each distinct root** `r`, form `M_r = M` with `omegaSquared → r`, and emit:

| tag | what to compute |
|---|---|
| **N1** | `M_r` itself |
| **N2** | `rank(M_r)`, and `nu = D − rank(M_r)` |
| **N3** | ⭐⭐ `rank` of the `(D+1)×D` matrix formed by **STACKING `M_r` on top of the single row `kᵀ`**, and `nu_T = D − that rank` |
| **N4** | `nu − nu_T` |
| **N5** | the vector `M_r · k`, **in full**, componentwise |
| **N6** | a null-space basis of `M_r`; and for **each** basis vector `b`: the scalar `b · k`, and the vector residual `(Σ_m k_m²)·b − (b·k)·k` |
| **N7** | `rank(M_r) + nu`, and that value minus `D` |

⭐⭐ **WHY N3 EXISTS, and it is the reason this step is being rebuilt.** `null(M_r) ∩ k^⊥` is exactly the
null space of `M_r` stacked on `kᵀ`, so `nu_T` counts the null directions having **no component along
`k`** — and it is a **rank**, so it is **independent of which basis anything returns**.

⛔⛔ **CLASSIFYING THE VECTORS RETURNED BY `NullSpace` IS NOT A SUBSTITUTE AND MUST NOT BE THE PRIMARY
RESULT.** `NullSpace` returns an **arbitrary** basis. If a root's null space happens to contain directions
that are neither parallel nor perpendicular to `k`, then **every** returned basis vector classifies as
"neither" and the count is **unrecoverable** from the classification — while `N3` returns the right number
regardless. ⚠ Whether that situation arises in any particular package is **not stated here and is not
your concern**: compute `N3` everywhere, unconditionally.

⇒ ⭐ **N6 is retained for display only.** Emit the **residual objects**, ⛔ not a classification token as
the headline. If you also emit a word like `parallel`/`perpendicular`/`neither`, it must be **derived from
the emitted residuals in the same expression** and it must never replace them.

### Q5 · Dispersion — by SCALING, ⛔ never by substitution

For each root `r(k)`, introduce a positive scalar symbol `lambdaScale` and emit:

```
r( lambdaScale · k )  −  lambdaScale² · r( k )
```

as an expression, with **both operands** also emitted.

⛔⛔ **DO NOT test this by substituting a symbol for `Σ k_m²`.** That substitution **silently no-ops**
whenever the root is not a bare multiple of that sum, which yields a false result in **both** directions.
⭐ Scaling `k → λk` has no such failure mode.

### Q6 · Dimensions — ⛔ two rules, both learned the hard way

Work in a three-slot dimension vector `(length, time, mass)`. `D` stays **symbolic**.

- ⭐⭐ **Dimension the WHOLE expression by walking its expression tree.** ⛔ **Never** obtain a dimension by
  reading the exponents of the coefficient symbols — that silently drops every other dimensionful factor,
  a wavevector for instance, and reports a dimension short by it. ⚠ That has produced a **wrong emitted
  value that looked like a meaningful signal.**
- ⭐⭐ **Count FIELD FACTORS, not `Derivative` atoms.** ⚠ A per-term analysis that sums over derivative
  nodes contributes **nothing** for a bare undifferentiated field, so any gap or mass term is silently
  mis-dimensioned — with a clean exit and every check green.

Compute and emit:
- `[ρ_br]` and `[μ_R]` as **closed functions of `D`**, obtained by requiring `L` to be an energy density on
  a `D`-dimensional brane. Emit the equations solved and the solution.
- ⭐ For **each root** `r`, the dimension of `r / (Σ_m k_m²)`, computed **from the emitted expression for
  `r` by walking its tree** — ⛔ not by differencing coefficient dimensions.
- The dimension of **every** emitted dimensionful expression, and its **homogeneity** — a boolean per
  expression, with the per-term dimension vectors emitted alongside so a failure is readable.

⚠ Homogeneity is **blind to a wrong dimensionless coefficient.** ⇒ it is a layer, ⛔ not the answer.

### Q7 · The `D = 3` comparison against the ordinary curl

At `D = 3` **only**, form the ordinary curl vector of the amplitude,
`c = (∂_2u_3 − ∂_3u_2, ∂_3u_1 − ∂_1u_3, ∂_1u_2 − ∂_2u_1)`, and emit:
`S_curl[∂u]` at `D=3`, `c · c`, and their **difference** — all three (clause 2).

⭐ Build the gradient components as **independent symbols**, so the two sides are genuinely two
expressions and not one substituted into itself.

### Q8 · Generic rank

⚠ `MatrixRank` / `Matrix.rank()` return the **GENERIC** rank — the rank away from measure-zero loci in the
symbols. Emit, for each `M_r`:
- the locus where the generic rank **drops** — take `r` to be the rank **your own computation returned**
  for that matrix, form the `r × r` minors, and solve them **all** to zero; emit the solution set;
- and the Q3 root-coincidence locus alongside it.

⭐ The orchestrator checks the enumeration is complete; ⛔ the engine does not assert that it is.

---

## 7 · Packages — the sweep and the controls

⛔⛔ **EVERY control re-enters at the ACTION (§2), never at a result.** A control is a **different `L`**;
everything downstream is recomputed from it by the identical code path.

⭐ **Every package emits the FULL Q1–Q8 tag set.** ⛔ Emission is never conditional on a value
(corollary 4). ⚠ A quantity that is **identical** across packages **must still be emitted in each** — that
repetition is a result, and deleting it destroys the evidence.

| package | `D` | what changes, at the ACTION |
|---|---|---|
| `MAIN` | 2, 3, 4, 5, 6 | nothing — `S_curl`, isotropic `ρ_br` |
| `XFORM_FULLGRAD` | 3, 4 | `S → Σ_{i,j} (∂_i u_j)²` |
| `XFORM_DIVONLY` | 3, 4 | `S → ( Σ_i ∂_i u_i )²` |
| `XFORM_SIGNFLIP` | 3, 4 | the stiffness term enters with `+ (μ_R/2) S_curl` |
| `XFORM_ANISO` | 3, 4 | kinetic term `→ (1/2) Σ_j ρ_j (∂_t u_j)²` with **distinct** positive symbols `ρ_1 … ρ_D` |
| `XCOEF_SCALE` | 3 | `μ_R → s · μ_R`, `s` a free positive symbol |

⭐ **`MAIN` swept over `D` is the primary result of the step, ⛔ not a control.**
⭐ **`XCOEF_SCALE` is a COEFFICIENT control and tests arithmetic only** — scaling never leaves the family.
⛔ It cannot substitute for a FORM control, and ⛔ **no expectation is stated here about what it does or
does not move.** Run it and emit its full tag set like every other package.

⚠ **Runtime.** No script may exceed **10 minutes**. If `D = 6` with a control is too slow, ⭐ **narrow the
`D` list and emit a tag recording exactly which `(package, D)` pairs were run and which were skipped.**
⛔ Never silently drop one, and ⛔ never raise a timeout.

---

## 8 · Tag grammar — ⛔ both engines must produce PARALLEL tag sets

```
<ENGINE>_S10_<PACKAGE>_D<n>_<QUANTITY>
<ENGINE>_S10_<PACKAGE>_D<n>_ROOT<r>_<QUANTITY>
```

- `<ENGINE>` is `WL` or `PY`.
- `<PACKAGE>` is exactly one of the §7 names.
- `<n>` is the integer brane dimension of that run.
- `<r>` indexes the **distinct** roots, ordered by the engine's own de-duplicated list, `1`-based.
- `<QUANTITY>` names **the object**, ⛔ never its value.

⭐ Emit one line per tag: `TAG: <payload>`, payload in a fully re-parseable form
(`InputForm` in Wolfram; `sympy.srepr`-safe string or `str()` of the expression in SymPy).
⭐ Also emit, once per package-and-dimension, a tag listing the **root ordering** used, so the two engines'
`ROOT<r>` indices can be aligned by the orchestrator rather than assumed to match.

⚠ **S9 finished with 16 tag-parity gaps between engines.** ⭐ This grammar is why; follow it exactly.

---

## 9 · What is supplied vs. what this build tests

| supplied — ⛔ NOT tested here | tested here |
|---|---|
| the curl-only stiffness form (from S9) | everything in Q1–Q8, for every package |
| `u` is the in-plane `D`-vector; the `h`-branon is a separate field | |
| `v₀ = 0`, no dissipation, exactly-quadratic `L` | |
| `ρ_br > 0`, `μ_R > 0`, `Σ k_m² > 0` | |

⚠ **Emit each supplied premise as its own tag** so a reader cannot mistake a passing build for having
verified one.

---

## 10 · Report back — ⛔ under 25 lines

1. The file you wrote, its line count, and the **total number of emitted tags**.
2. The `(package, D)` pairs actually run, and any skipped, with the runtime.
3. ⭐ Anything in this specification you found **ambiguous, under-determined, ill-posed, or tautological** —
   including any requested object that **cannot be computed** from what §1–§3 supplies. ⭐ This is wanted
   and is more valuable than a clean build.
4. ⛔ **Do NOT report what any value came out to be, do not interpret any result, and do not say whether
   anything "worked".** Your job ends at compute-and-print.
