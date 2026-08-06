# S11 — SHARED PHYSICS

⛔⛔ **This file is read by BOTH engines. An error here defeats dual-engine agreement by construction:
both engines will compute the same wrong thing and agree.** It is reviewed as its own artifact before
either build starts.

⭐ It supplies **the setup, the field content, the governing equations, the premises, the list of
quantities to compute.** ⛔ It supplies **no result** — no root, no
count, no dimension tuple, no sign, no locus, no invariant, and no expectation about what any control does
or does not move.

⚠ **Question numbers are shared with `S10_SHARED_PHYSICS.md`.** `Q1`–`Q8` name the same objects there and
here. `Q9`–`Q11` are new at S11.

---

## 1 · The physical setup

A **brane** is a `D`-dimensional elastic sheet. `D` is a **free integer symbol** wherever an expression can
be written in general `D`; §7 also names explicit integer values at which the whole chain is re-run.
⛔ Do not globally specialise `D` to one value and derive everything from it.

The field is `u`, the **in-plane displacement**: a `D`-component vector field on the brane,
`u_j(x_1, …, x_D, t)`, `j = 1…D`.

⭐ **`u` is a DISPLACEMENT, so `[u] = length = (1, 0, 0)` in the `(length, time, mass)` slots.** ⛔ This is
a **supplied premise**, ⛔ not something to infer — §6's dimension solve is under-determined without it, and
two engines left to choose independently can pick differently and disagree for no physical reason.
⭐ Emit it as a premise tag.

⚠ **SUPPLIED PREMISE — unfalsifiable within this build. State it in the output as a premise tag.**
`u` is **in-plane only**. Motion *out of* the brane is a **different field** (the `h`-branon) belonging to
a different sector, and it is **not** part of `u`. ⇒ `u` has exactly `D` components, ⛔ not `D+1`.
⭐ Nothing in this build tests that split; it is inherited.

Write `G` for the `D×D` **displacement-gradient matrix**, `G_ij ≔ ∂_i u_j`. Every stiffness term below is
a quadratic form in `G`.

## 2 · The action — supplied in full, as equations

```
L_pkg  =  (ρ_br / 2) · Σ_j (∂_t u_j)²   −   W_pkg[G]
```

`W_pkg` is the package's **stiffness functional**, given in full for each package in §7. `ρ_br > 0` is a
real constant.

⚠ **`W` is a FUNCTIONAL, not a fixed coefficient times a fixed density.** ⛔ Do not assume which
coefficients a given package's `W` contains, or how many. §6 says to read them off that package's own `W`.

The stiffness densities `W` is built from, each defined in every `D`:

```
S_curl[G]   =  (1/2) · Σ_{i=1..D} Σ_{j=1..D} ( G_ij  −  G_ji )²
S_div[G]    =  ( Σ_{i=1..D} G_ii )²
S_symtl[G]  =  Σ_{i=1..D} Σ_{j=1..D} ( Stl_ij )²   ,
                 Stl  ≔  (1/2)(G + Gᵀ)  −  (1/D)·(tr G)·I      the symmetric traceless part of G
P_D[G]      =  ⭐ DEFINED BY §Q9's COMPUTED OUTPUT — see §7. ⛔ No formula for it appears in this file.
```

`μ_R > 0`, `B_comp > 0`, `μ_br > 0` and `β` are real constants. ⚠ **`β` carries no sign premise**; §7 says
why.

⛔⛔ **DO NOT SIMPLIFY ANY OF THESE BY HAND, IN EITHER COORDINATE SPACE OR AFTER SUBSTITUTING THE ANSATZ.**
Whatever a stiffness density becomes on a plane wave is a **computed object** and is one of the things this
build exists to produce. ⛔ Typing its reduced form is the defect this whole rebuild removes.

## 3 · The ansatz — supplied

⛔⛔ **USE THE REAL ANSATZ. THIS IS A BINDING EQUATION, ⛔ NOT A CONVENTION YOU MAY CHOOSE.**

```
u_j(x, t)  =  a_j · cos( Σ_m k_m x_m  −  ω t )
```

with `a` and `k` real `D`-component symbolic vectors, and **every quadratic density averaged over one
period of the PHASE**. ⛔⛔ **Average over the phase variable, ⛔ NEVER over `t` with limits containing `ω`:**

```
φ  ≡  Σ_m k_m x_m  −  ω t                 the phase
⟨ F ⟩  ≡  (1 / 2π) · ∫_0^{2π} F  dφ
```

⛔⛔ **DO NOT WRITE THIS AS `(ω/2π) ∫_0^{2π/ω} dt`.** ⚠ Those limits are a real period **only if `ω` is
real and nonzero**, assumptions not supplied here; using them would force a sign assumption on `ω²` into
the matrix that the sign test later examines. ⭐ The phase average above has **no such dependence**.

⛔ Do **not** add any assumption about the sign of `ω²` anywhere in the construction of `M`.

⛔⛔ **DO NOT use `a_j exp(i(k·x − ωt))` with real `a`.** It does **not** define a real
`(G_ij − G_ji)²` — substituting `∂ → i k` multiplies that square by `i² = −1` and **flips the sign of the
entire stiffness term**. ⇒ two engines free to choose the convention would differ by a sign and a factor,
for no physical reason.

⚠⚠ **And at S11 a convention error is no longer a common factor.** S11's `W` is a **sum of densities with
independent coefficients**; a convention applied inconsistently across the summands changes their
**relative** weight, which is exactly what §Q3 and §Q10 measure. ⛔ Apply one convention to the whole
Lagrangian.

⭐ **Emit whatever overall factor you obtain; ⛔ do not normalise it away**, and ⛔ do not assume it is
harmless — §Q1, §Q2, §Q3 and §Q7 payloads all carry it. ⚠ An overall positive constant common to the whole
quadratic form does not change the **dimension** solve of §6; it changes those payloads.

Write `ω²` as a **single symbol** (`omegaSquared`) throughout, ⛔ never as `ω` squared; solving for `ω`
introduces sign branches that are not separate physical modes.

**Assumptions in force everywhere — ⛔ BOTH ENGINES MUST CARRY THE IDENTICAL SET, and must carry it as a
JOINT assumption, not merely as per-symbol declarations:**

```
ρ_br > 0 ,  μ_R > 0 ,  B_comp > 0 ,  Σ_m k_m² > 0 ,  every k_m real ,  every a_j real ,
D a positive integer

package-domain premises, joined to the common set when that package is selected:
XFORM_TRACELESS:  μ_br > 0
XFORM_EXTRA:      β real                      ⛔ no sign premise, ⛔ no premise that β ≠ 0
XCOEF_BSCALE:     s > 0 ,  [s] = (0, 0, 0) ,  s ≠ 1
```

⚠ Declaring each symbol real/positive **individually is not equivalent** to asserting `Σ_m k_m² > 0`, and a
sign test performed under a different premise set is not the requested computation. ⇒ ⛔ Both engines
must carry the joint set.

⛔⛔ **NO PREMISE RELATING THE STIFFNESS COEFFICIENTS TO ONE ANOTHER IS SUPPLIED, AND NONE MAY BE ADDED.**
⚠ Do not assume any two stiffness coefficients are unequal, ordered, or generic with respect to each
other. ⭐ Where the coefficients may take special values is a **computed object** — §Q8 asks for it. A
genericity premise would delete that locus before it could be found.

⚠ **SUPPLIED PREMISES — unfalsifiable within this build. Emit each as a premise tag:**
- the background brane is **unstrained and at rest**: `v₀ = 0`, so there are **no convective terms** and no
  advective rate `k·v₀`;
- **no dissipation**: `L` contains **no term odd in `∂_t`** — no first time derivative, no relaxation time,
  no memory kernel;
- **linear response / small oscillations**: `L` is exactly quadratic in `u`, so no amplitude-dependent
  restoring force is retained;
- the brane's **wall width is frozen**: no field describing the sheet's thickness appears in `L`, so no
  stiffness in this file may depend on a width or its gradients.

⭐ Why these are stated rather than assumed silently: each removes something that could otherwise act, and
⛔ **any statement this build makes about the relationship between different spatial directions inherits
`v₀ = 0`** — it would then be a property of the medium **on a state at rest**, ⛔ not of the medium.

## 4 · ⭐⭐⭐ THE STRUCTURAL RULE — verbatim, non-negotiable

> **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
> ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
> the chain at the ACTION, ⛔ never at a result.**

⛔⛔ **A hand-typed CAS object is still hand-typed.** `emit[Simplify[(alpha + beta)^3]]` is a genuine CAS
object with **no data dependency on the derivation**: delete `Det` and `Solve` and the output does not
move. ⭐ **The test: if you deleted the computation, would this tag change?** If not, it is dead.

⚠ **That example uses unrelated symbols AND an unrelated SHAPE, deliberately.** ⛔ Swapping the symbols of
a real answer is **not** enough — a list of the right length, with a zero in the right slot, states the
result no matter what the other entries are called.

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
   sign, parity, or the shape of the answer. ⚠ Schematically — with **placeholder** content, ⛔ not this
   step's: `..._EQUALS_SEVEN`, `..._IS_POSITIVE`, `..._QUADRATIC_IN_Z`, `..._DEPENDS_ON_ZETA` are forbidden
   **name** patterns whatever the payload, because the name has already answered. ⭐ `..._RANK`,
   `..._ROOTS`, `..._RESIDUAL`, `..._JACOBIAN` are fine.
3. ⛔ **No tautological residual.** Before emitting a difference, ask: were these two operands produced by
   **independent routes**? ⭐⭐ **The test, and it is not the obvious one:** mutate the **source** the
   residual is supposed to police, recompute everything downstream, and see whether the residual moves. A
   residual that only moves when you corrupt an **already-computed operand** is testing your own
   arithmetic, ⛔ not the object. ⭐ Where no second route exists, emit the objects and ⛔ emit no residual.
4. ⛔⛔ **EMISSION MUST NEVER BE CONDITIONAL ON A PAYLOAD'S VALUE.** Whether a tag appears may depend only
   on **which package**, **which dimension**, and **which quantity** it belongs to. ⚠ A value **present and
   identical across packages is a RESULT**; a value **absent** is indistinguishable from *never computed*.
   ⛔ Never de-duplicate payloads across packages. Tag **names** are unique; **payloads may legitimately
   repeat.** ⇒ ⛔ **No tag for a named object, at the scope §8 gives it, may be emitted "if" anything.**

### ⛔⛔ COROLLARY 5 · A DECLARATION MUST BE WIRED TO WHAT IT DECLARES

A tag that *declares* an input — a premise, a convention, a route, a count — must be **produced from the
object the computation actually uses**. ⛔ Re-typing the same value a second time is a **hand-authored
payload with no data dependency**, and it reads to every consumer as a report of what was used.

⭐⭐ **THE TEST: perturb the thing the tag declares. The tag MUST move.**

⇒ ⭐ **Derive every declarative payload from the live object** — read it back out of the structure the
computation consumed, ⛔ never from a literal beside it.
⚠ Where a value genuinely IS a supplied constant with no in-code consumer, ⭐ **say so in the payload**, and
⛔ do not name the tag as though it reported a computation.

### ⛔⛔ NO VERDICT

⛔ **There is no `VERDICT` tag, no `PASS`, no `FAIL`, and no summary judgement anywhere in either engine.**
⛔ **A physics finding must EXIT 0.** Exit non-zero **only** on operational failure — an exception, an
unsolvable system. ⚠ Otherwise a builder iterating to exit 0 can make a genuine disagreement disappear, and
the disagreement is the most valuable output available.

⚠⚠ **This is a change from the engines being replaced.** Both as-built S11 engines end in a `VERDICT` tag,
and one renders a symbolic boolean as the typed word `TRUE`/`FALSE` with the residual discarded. ⛔ Neither
survives here. A boolean is emitted **as the CAS object the test returned**, with its operands beside it.

### ⛔⛔ THE LOCUS PROTOCOL — used by Q3, Q8 and Q11, and ⛔ never varied

⚠⚠ **"Solve this over the reals" is NOT a specification.** Measured on the CAS both engines will use:
a multivariate solve **ignores** the symbols' real declarations and returns solutions in the algebraic
closure; the single-variable real-domain solver **rejects** a tuple of variables; and the ordinary solver
returns the **same empty token** for a system that is identically satisfied and one that is inconsistent.
⇒ ⛔ An instruction to "solve over the reals" asks one engine for something it cannot do and makes
*checked-and-empty* indistinguishable from *true-everywhere*.

⭐⭐ **Therefore, wherever this file asks for a locus, emit ALL FIVE of these objects and nothing less:**

| suffix | the object |
|---|---|
| `_EQUATIONS` | the defining equation system, as a list of CAS relations, **before** any solve |
| `_SOLUTION` | the solution set exactly as your CAS returns it, with the **solve variables named** and the domain left unrestricted |
| `_IDENTICALLY_SATISFIED` | the boolean returned by testing whether **every** equation in `_EQUATIONS` simplifies to zero identically in its variables |
| `_INCONSISTENT` | the boolean returned by testing whether the system admits **no** solution at all |
| `_REAL_ADMISSIBLE` | for **each** branch in `_SOLUTION`, the symbolic test of whether that branch admits a point satisfying **every** §3 assumption, emitted **with its operands** |

⭐ `_IDENTICALLY_SATISFIED` and `_INCONSISTENT` are what separate the three degenerate cases, and they are
computed from `_EQUATIONS`, ⛔ never read off the solver's empty token.
⛔ Emit all five for every locus this file requests, unconditionally, whatever they come out to.

---

## 6 · What to compute — the question list

⛔ **This section says what to produce. It does not say what any of it equals.** Everything below is
required for **every package** defined in §7, at **every** `D` in that package's sweep.

### ⭐⭐ Ordering requirement

`§Q9` is a function of `D` alone and one package's action is built from its output. ⇒ ⭐ **Compute `Q9`
first for each `D`**, before assembling any action at that `D`.

### Q1 · The Lagrangian and the equations of motion

Build `L` in coordinate space from §2. Emit `L` expanded, and emit **each additive term of `W_pkg`
separately**, as the list `STIFFNESS_TERMS`, so a consumer can see which densities that package carries and
with which coefficients.
Take the **Euler–Lagrange** variation with respect to each `u_j`:

```
∂_t ( ∂L / ∂(∂_t u_j) )  +  Σ_i ∂_i ( ∂L / ∂(∂_i u_j) )  −  ∂L / ∂u_j  =  0
```

Emit the full system.

### Q2 · The dynamical matrix, by two routes — ⚠ and BE HONEST ABOUT WHAT THAT DOES AND DOES NOT TEST

`M` is the `D×D` coefficient matrix of the plane-wave amplitude equations, `M · a = 0`.

- **Route A — through the equations of motion.** Substitute the §3 ansatz into the Q1 Euler–Lagrange
  system, strip the common trigonometric factor, and extract the coefficient matrix. Call it `M_A`.
  ⚠ **An equation equal to zero does not fix an overall scale or row sign.** ⛔ Do not choose a
  normalisation to make `M_A` match anything — strip only the factor common to every entry, and ⭐ emit the
  factor you stripped as `M_ROUTE_A_STRIPPED_FACTOR`.
- **Route B — through the quadratic form.** Period-average the §2 quadratic density on the §3 ansatz and
  take `M_B[[i,j]] = ∂²⟨L⟩/∂a_i ∂a_j`. Call it `M_B`.

⭐ **Emit `M_A`, `M_B`, `M_A − M_B`, and the ratio `M_A[[1,1]]/M_B[[1,1]]` — all four, UNCONDITIONALLY**
(corollary 4). ⛔ Do not normalise either route to the other, and ⛔ do not make any of these four tags
conditional on the two matrices agreeing.

⛔⛔ **THESE ARE NOT TWO INDEPENDENT PHYSICS DERIVATIONS.** Both routes are built from the same action.
⇒ ⭐ **The residual tests CODING CONSISTENCY ONLY.** ⭐ Emit it, and emit a tag saying that is what it
tests. ⛔ Do not present it as verifying physics.

⭐ Use `M_B` for everything downstream, and emit `M_ROUTE_USED` naming the route, derived from the object
actually consumed (corollary 5).

⭐ Also emit `M_COEFFICIENT_JACOBIAN`: for each coefficient in the package's emitted `COEFFICIENT_ORDERING`,
the matrix `∂M/∂c`, as a list in that order. ⛔ Do not describe the pattern; emit the matrices.

### Q3 · The spectrum

- Emit `DET_M`, **factored**.
- Solve `DET_M = 0` for `omegaSquared`. Emit:
  - `ROOT_SOLUTION_SET` — the solution set **as returned, retaining multiplicity**;
  - `ROOT_DISTINCT` — the list after de-duplication under the §3 assumptions;
  - `ROOT_COUNT_ALL` and `ROOT_COUNT_DISTINCT` — each an integer produced by **counting its own list**;
  - `ROOT_ORDERING` — the ordering of `ROOT_DISTINCT` that the `ROOT<r>` index refers to.
- ⭐ For each distinct root, `ROOT<r>_VALUE` and `ROOT<r>_SIGN`. ⛔ `ROOT<r>_SIGN` is a **four-way
  symbolic test** whose payload is exactly one of the tokens `POSITIVE`, `NEGATIVE`, `ZERO`, `UNDECIDED`,
  emitted **together with** the operand expression the test was applied to. ⚠ Without a pinned token set
  the two engines emit different spellings of the same determination and the comparison reports a
  disagreement that is not one.
- ⭐ The **root-coincidence locus in the WAVEVECTOR**, `ROOT_COINCIDENCE_K_*`: for every pair of distinct
  roots, the system `r_p − r_q = 0` with solve variables `k_1 … k_D`. Emit the **full five-object locus
  protocol** from §5.
- ⭐⭐ The **root-coincidence locus in the STIFFNESS COEFFICIENTS**, `ROOT_COINCIDENCE_COEFF_*`: the same
  systems, with solve variables the coefficients in `COEFFICIENT_ORDERING`. Emit the **full five-object
  locus protocol**. ⚠ **This object does not exist at S10**: its action had one stiffness coefficient, so no
  coefficient-space coincidence was possible. ⛔ Solving only in `k` would report an empty locus for a
  coincidence that is a condition on the coefficients alone.
- ⭐ When `ROOT_COUNT_DISTINCT < 2` there are no pairs; ⛔ still emit both locus tag sets, with
  `_EQUATIONS` the empty list (corollary 4).

### Q4 · ⭐⭐⭐ THE MODE COUNT — basis-independent tests

For **each distinct root** `r`, form `M_r = M` with `omegaSquared → r`, and emit:

| tag | what to compute |
|---|---|
| `ROOT<r>_N1` | `M_r` itself |
| `ROOT<r>_N2_RANK`, `ROOT<r>_N2_NULLITY` | `rank(M_r)`, and `nu = D − rank(M_r)` |
| `ROOT<r>_N3_STACKED_RANK`, `ROOT<r>_N3_TRANSVERSE_NULLITY` | ⭐⭐ `rank` of the `(D+1)×D` matrix formed by **STACKING `M_r` on top of the single row `kᵀ`**, and `nu_T = D − that rank` |
| `ROOT<r>_N4_NULLITY_DIFFERENCE` | `nu − nu_T` |
| `ROOT<r>_N5_M_DOT_K` | the vector `M_r · k`, **in full**, componentwise |
| `ROOT<r>_N6_BASIS`, `ROOT<r>_N6_DOT_K`, `ROOT<r>_N6_RESIDUAL` | a null-space basis of `M_r`; for **each** basis vector `b`, the scalar `b · k`, and the vector residual `(Σ_m k_m²)·b − (b·k)·k` |
| `ROOT<r>_N7_BASIS_COUNT`, `ROOT<r>_N7_RESIDUAL` | ⭐ the **number of vectors the null-space routine returns**, and that count `− nu` |

⭐⭐ **WHY N3 EXISTS.** `null(M_r) ∩ k^⊥` is exactly the null space of `M_r` stacked on `kᵀ`, so `nu_T`
counts the null directions having **no component along `k`** — and it is a **rank**, so it is
**independent of which basis anything returns**.

⛔⛔ **CLASSIFYING THE VECTORS RETURNED BY `NullSpace` IS NOT A SUBSTITUTE AND MUST NOT BE THE PRIMARY
RESULT.** `NullSpace` returns an **arbitrary** basis. If a root's null space contains directions that are
neither parallel nor perpendicular to `k`, then **every** returned basis vector classifies as "neither" and
the count is **unrecoverable** from the classification — while `N3` returns the right number regardless.
⚠ Whether that arises in any particular package is **not stated here and is not your concern**: compute
`N3` everywhere, unconditionally.

⚠⚠ **§7 contains a package whose stiffness functional is not built from `S_curl` and `S_div` alone, and
nothing in this file says how its null directions lie relative to `k`.** ⛔ Any code path that assumes a
two-way parallel/perpendicular split will report a wrong count without failing.

⇒ ⭐ **`N6` is retained for display only.** Emit the **residual objects**, ⛔ not a classification token as
the headline. ⚠ `N6_BASIS` is **not** cross-engine comparable — a null-space basis is not canonical — and
it is emitted for inspection, ⛔ not as a comparison row.

⚠ **`N7` differences two DIFFERENT ALGORITHMS** — the null-space routine's basis count against the rank
routine's implied count — which is what lets it fail. ⛔ Do not "simplify" it to `rank + nu − D`: with `nu`
defined as `D − rank` that is identically zero for any rank at all.

### Q5 · Dispersion — by SCALING, ⛔ never by substitution

For each root `r(k)`, introduce a positive scalar symbol `lambdaScale` and emit `ROOT<r>_SCALED`,
`ROOT<r>_UNSCALED`, `ROOT<r>_SCALE_RATIO` and `ROOT<r>_SCALE_EXPONENT`:

```
r( lambdaScale · k )          the root at the scaled wavevector
r( k )                        the root at the original wavevector
r( lambdaScale · k ) / r( k )        ⭐ fully simplified — THE OBJECT THAT ANSWERS THIS
```

⭐⭐ **Emit the RATIO and let whatever power of `lambdaScale` appears, appear.** ⛔ **Do not difference the
scaled root against `lambdaScale` raised to a stated power** — naming a power tells you the answer before
you compute it.

`ROOT<r>_SCALE_EXPONENT` is obtained from the ratio by `Solve`/`solve` or `Exponent`/`degree` — ⛔ never
typed.

⛔⛔ **EMIT ALL FOUR TAGS FOR EVERY ROOT, UNCONDITIONALLY.** ⛔ **A missing tag is indistinguishable from a
quantity never computed**, which is the defect this whole specification exists to remove.

⭐⭐ **THE MARKER TOKEN IS PINNED, AND THIS IS NOT A STYLE CHOICE.** Where a requested object is not
defined — the ratio at a vanishing root, an exponent that is not a pure power — the payload is **exactly**
one of these tokens and nothing else:

```
UNDEFINED_RATIO          the ratio is not defined at this root
NOT_A_PURE_POWER         the ratio is defined but is not lambdaScale to any single power
```

⚠ **Why it is pinned here rather than left to each engine:** at the previous step the two engines chose
their own spellings for exactly this situation and **eleven** comparison rows reported a disagreement that
was two names for the same determination. ⛔ Do not invent a third token, and ⛔ do not emit an absence.

⛔⛔ **DO NOT test this by substituting a symbol for `Σ k_m²`.** That substitution **silently no-ops**
whenever the root is not a bare multiple of that sum, which yields a false result in **both** directions.
⭐ Scaling `k → λk` has no such failure mode.

### Q6 · Dimensions — ⛔ two rules, both learned the hard way

Work in a three-slot dimension vector `(length, time, mass)`. `D` stays **symbolic**.

```
[energy] = ( 2, −2,  1)
[L]      = [energy] · length^(−D)  = ( 2−D, −2, 1)      an energy DENSITY on a D-dimensional brane
[∂_i]    = (−1,  0,  0)
[∂_t]    = ( 0, −1,  0)
```

- ⭐⭐ **Dimension the WHOLE expression by walking its expression tree.** ⛔ **Never** obtain a dimension by
  reading the exponents of the coefficient symbols — that silently drops every other dimensionful factor,
  a wavevector for instance, and reports a dimension short by it.
- ⭐⭐ **Count FIELD FACTORS, not `Derivative` atoms.** ⚠ A per-term analysis that sums over derivative
  nodes contributes **nothing** for a bare undifferentiated field, so any gap or mass term is silently
  mis-dimensioned — with a clean exit and every check green.

⭐⭐ **THE COEFFICIENT INVENTORY, defined so both engines build the same list.** For each additive term of
`L_pkg`, take the factor multiplying the density, and take its **free symbols**. `COEFFICIENT_ORDERING` is
the union of those symbol sets over all terms, together with the inertial coefficient, sorted by the
engine's own stated rule and emitted as its own tag. ⚠ ⇒ a term whose factor is a product contributes
**every** symbol in that product, ⛔ not the product as one unknown.

Emit:
- `COEFFICIENT_ORDERING`, and `DIM_COEFFICIENTS` — the dimension of every coefficient in it, as **closed
  functions of `D`**, obtained by requiring `L` to be an energy density on a `D`-dimensional brane;
  `DIM_EQUATIONS` — the equations solved; `DIM_SOLUTION` — the solution as returned.
- ⚠ **The energy-density requirement alone does not fix every coefficient separately** — where a term
  carries a product of undetermined coefficients it fixes only their **sum of dimensions**. ⭐ Emit the
  underdetermined solution **as returned**, with its free parameters visible, ⛔ never by silently choosing
  one. Where §7 declares a control parameter **dimensionless**, use that as a stated input and emit it as a
  premise tag.
- `ROOT<r>_DIM_OVER_KSQ` — for **each root**, the dimension of `r / (Σ_m k_m²)`, computed **from the
  emitted expression for `r` by walking its tree** — ⛔ not by differencing coefficient dimensions.

⭐⭐ **THE HOMOGENEITY INVENTORY IS FINITE AND IS EXACTLY THIS LIST** — ⛔ not "every emitted expression":

```
DIM_HOMOGENEITY_ACTION     each additive term of L_pkg
DIM_HOMOGENEITY_DERIVED    each of: DET_M, every ROOT<r>_VALUE, every ROOT<r>_N5_M_DOT_K,
                           every ROOT<r>_N6_RESIDUAL, Q7_RESIDUAL, every ROOT<r>_KW_SQUARED
```

Each payload is a boolean **per expression**, with the per-term dimension vectors emitted alongside so a
failure is readable. ⛔ Emit the two classes under the two distinct tag names above, so a consumer can tell
which class it is reading.

### ⛔⛔ Q6d · MAKE THE HOMOGENEITY CHECK'S OWN VACUITY VISIBLE

⭐ **Why, and it is structural, ⛔ not a coding bug:** the coefficient dimensions are **solved for** by
requiring each action term to equal the energy density. Homogeneity is then evaluated **under that same
solution** ⇒ true by construction. ⚠ **And it cannot be repaired by adding a cleverer check built from the
coefficients:** `[u]` enters **every** coefficient with the same weight, so it **cancels identically** in
any ratio or difference of them — including in `[ω²]`, which is therefore blind to `[u]` as well.

⇒ ⭐⭐ **Emit the diagnostic instead of pretending the check is one.** All four, unconditionally:

```
DIM_EQUATION_COUNT      the number of independent dimension equations formed, as an integer
DIM_UNKNOWN_COUNT       the number of unknown coefficient dimensions solved for, as an integer
                        ⛔⛔ counted from COEFFICIENT_ORDERING, never from a fixed list of symbols
DIM_COUNT_DIFFERENCE    their difference
DIM_DETERMINACY         from the solve itself: OVER_DETERMINED | EXACTLY_DETERMINED | UNDER_DETERMINED
```

⭐ **The counting rule, and it must be stated because two engines read it oppositely otherwise:** a symbol
**declared** dimensionless in §3 or §7 is ⛔ **not** a `Q6` unknown; a symbol whose dimension the action
**determines** ⭐ **is** one, and stays one.

⛔ Emit all four whatever they come out to, and ⛔ do not add a label tag that fires only in one case
(corollary 4). ⭐ `DIM_COUNT_DIFFERENCE` is the operand; classifying it belongs to the step record.

⚠ **Homogeneity retains its value elsewhere** — on expressions the solve never touched, an inhomogeneity
**is** a real finding. That is why the two classes carry different tag names.

⚠ Homogeneity is also **blind to a wrong dimensionless coefficient.** ⇒ it is a layer, ⛔ not the answer.

⭐ **`[u]` is a PREMISE and is UNFALSIFIABLE WITHIN THIS BUILD.** ⛔ Emit it as a premise tag, and ⛔ do not
let any tag name suggest it was verified.

#### ⛔⛔ Q6r · THE REGISTRY COMPARISON — and ITS PROVENANCE IS PART OF THE PAYLOAD

⚠ **Engine-local: only the engine that loads the registry emits this**, under the `_LOCAL_` infix of §8.
The registry is the YAML quantity register under `research/pde_ledger_v3/reduction/`; read it with that
directory's own reader rather than parsing it by hand.

For **each** coefficient in `COEFFICIENT_ORDERING` that the registry also declares, emit: the **derived**
vector, the **declared** vector, their **difference**, and ⭐⭐ the declared row's **`source_locus` read out
of the registry itself**.

⛔⛔ **The `source_locus` is not decoration.** A declared dimension whose locus points at an artifact this
build replaces is **not an independent operand**, and a zero residual against it means only that the new
engine reproduces its predecessor. ⭐ Emitting the locus is what lets the step record tell the two cases
apart. ⛔ Do not filter, judge or omit any row on the basis of what its locus says — emit every row and let
the record classify.

⛔ If a derived and a declared vector disagree, **emit the disagreement and continue**. ⛔ Do not adjust the
derivation to match, and ⛔ do not exit non-zero.

### Q7 · The `D = 3` comparison against the ordinary curl

The three-vector cross product exists only at `D = 3`. Emit, at `D = 3` only, with

```
g_ij   := independent symbols standing for ∂_i u_j
c_i    := Σ_{j,k} ε_ijk g_jk
```

| tag | the object |
|---|---|
| `Q7_W_FULL` | this package's **complete** `W_pkg`, with the `g_ij` substituted in |
| `Q7_CURL_TERM` | the additive term of `W_pkg` whose density is `S_curl`, **with its coefficient**, `g_ij` substituted — ⭐ the **zero form** when the package has no such term |
| `Q7_CURL_DENSITY` | the §2 density `S_curl` itself, **unweighted**, `g_ij` substituted |
| `Q7_CURL_REFERENCE` | `c·c` |
| `Q7_RESIDUAL` | `Q7_CURL_DENSITY − Q7_CURL_REFERENCE` |

⭐ **The zero form is a legitimate payload, ⛔ not a missing one**, and it is what a package carrying no
curl term emits. ⛔ Do not invent a curl term for such a package, and ⛔ do not omit the tag.

The CAS must compute `c_i` from the Levi-Civita definition and compute `c·c` from that result. ⛔ Do not
write out its expanded polynomial. ⚠ This requirement is new relative to the engines being replaced, where
both authors typed the three components of `c` by hand.

⚠ `Q7_RESIDUAL` compares two objects that are both fixed by §2 and are the same for every package. ⇒ ⭐ it
is a check on the Levi-Civita construction and the `½` convention, ⛔ **not** evidence about this package's
action. **What ties Q7 to the action is `Q7_W_FULL` and `Q7_CURL_TERM`**, which must be the objects from
which the action was assembled: mutating the stiffness functional, with the package selector held fixed,
must move them. ⛔ It is not evidence for the mode count, which is computed in every `D`.

That mutation requirement establishes dependence on the action; it does ⛔ **not** establish that the
emitted density is the one from which the action was assembled. Nothing inside a single engine does.

### Q8 · ⭐⭐⭐ Generic rank — and RE-RUNNING Q3 AND Q4 WHERE IT DROPS

⚠ `MatrixRank` / `Matrix.rank()` return the **GENERIC** rank — the rank away from measure-zero loci in the
symbols. ⇒ ⛔⛔ **EVERY NUMBER Q4 PRODUCES IS A GENERIC NUMBER, AND SAYING SO IS NOT ENOUGH.**

⛔⛔ **AND AT S11 "THE SYMBOLS" INCLUDES THE STIFFNESS COEFFICIENTS.** S10 swept loci in the wavevector
only, because its action carried one stiffness coefficient. ⭐ Here a locus may live in the coefficients, in
the wavevector, or in both, and §3 supplies no premise excluding any of them.

**Q8a — the rank-drop loci.** For each `M_r`, take `ρ` to be the rank **your own computation returned**,
form **all** `ρ × ρ` minors, and emit `ROOT<r>_RANK_DROP_MINORS` (the minor list) plus the **full five-object locus
protocol** for each of three solve-variable sets, under the suffixes:

```
ROOT<r>_RANK_DROP_K_*         solve variables: k_1 … k_D
ROOT<r>_RANK_DROP_COEFF_*     solve variables: COEFFICIENT_ORDERING
ROOT<r>_RANK_DROP_JOINT_*     solve variables: both, jointly
```

⚠ **If the returned rank `ρ` is `0`** there are no `ρ × ρ` minors and the matrix is identically zero.
⭐ Emit `ROOT<r>_RANK_DROP_MINORS` as the **empty list** and the locus `_EQUATIONS` as the empty list; the protocol's
`_IDENTICALLY_SATISFIED` then carries the meaning. ⛔ Do not emit a `0×0` determinant, which some CAS report
as `1`.

**Q8b — ⭐⭐ STRATA, defined.** A **stratum** is a component of a rank-drop locus that **intersects the
region allowed by §3**. Emit the components your CAS returns, each with its defining equations and chosen
point, plus their ordering:

```
STRATUM_ORDERING              the list of components, in the order the STRATUM<s> index refers to
STRATUM<s>_DEFINING_EQUATIONS the equations cutting out that component
STRATUM<s>_POINT              ⭐ an explicit point on it satisfying every §3 assumption
STRATUM<s>_POINT_RESIDUAL     that point substituted back into STRATUM<s>_DEFINING_EQUATIONS
```

⭐ Then **recompute and emit the complete `Q3` and `Q4` tag sets at `STRATUM<s>_POINT`**, under that
stratum's tag scope.

⛔⛔ **THIS IS NOT BOOKKEEPING.** A generic rank does not determine the rank on an exceptional stratum.
Reporting a locus without recomputing on it leaves the spectrum and mode count there untested.

⚠ **A single point does not characterise a positive-dimensional component**, and this file does not ask it
to. ⭐ `STRATUM<s>_DEFINING_EQUATIONS` is what carries the component; the point is where `Q3`/`Q4` are
re-evaluated, and it is emitted as an operand so the orchestrator can see which point was used.

⭐ If a package has **no** stratum, emit `STRATUM_ORDERING` as the **empty list** — ⛔ do not omit the tag
(corollary 4). ⭐ That tag is what distinguishes *checked and empty* from *never checked*.

⭐ Stratum indices are aligned by the orchestrator from the emitted defining equations; they are never
assumed to correspond across engines.

### Q9 · ⭐⭐⭐ THE INVARIANT CENSUS

⛔⛔ **NAME THE OBJECT; ⛔ DO NOT FOLLOW A RECIPE.** ⚠ The engines being replaced were told to *"decompose
the matrix into pieces that do not mix under rotation and count them."* ⭐ **That instruction is
withdrawn** — it is a derivation path, ⛔ not an object, and it can miss contributions built by pairing two
summands that are isomorphic as representations. ⇒ ⭐ **Compute the object defined below by whatever route
your CAS supports, and emit what you built.**

Let `Quad(D)` be the vector space of quadratic forms in the `D²` entries of `G`, and for an orthogonal
matrix `R` let `G ↦ R · G · Rᵀ` act on it.

⭐ **The monomial ordering is pinned as follows.** Set `g_((i−1)D+j) = G_ij` for `1 ≤ i,j ≤ D`, so
`g_1, …, g_(D²)` lists the entries of `G` in row-major order. Enumerate the pairs `(p,q)` with `1 ≤ p ≤ q ≤ D²`
lexicographically, first by `p` and then by `q`, and use the corresponding list of monomials `g_p g_q`.

| tag | the object |
|---|---|
| `MONOMIAL_ORDERING` | the pinned ordering above, emitted as the coordinate list for `Quad(D)` |
| `V1_BASIS`, `V1_DIM` | a basis for the subspace invariant under that action for **every** `R` with `Rᵀ R = I` and `Det R = 1`; and its dimension, counted from the basis |
| `V2_BASIS`, `V2_DIM` | a basis for the subspace invariant for **every** `R` with `Rᵀ R = I` — the **full** orthogonal group, `Det R` unrestricted; and its dimension, likewise counted |
| `V3_DIFFERENCE` | `V1_DIM − V2_DIM` |
| `R0_MATRIX`, `R0_ORTHOGONALITY_RESIDUAL`, `R0_DETERMINANT` | an explicit real `R₀` with `R₀ᵀ R₀ = I` and `Det R₀ = −1`; `R₀ᵀ R₀ − I`; `Det R₀` |
| `V4_REFLECTED`, `V4_RESIDUAL`, `V4_SUM` | for **each** `V1` basis element `p`, in the `V1` ordering: `p(R₀ G R₀ᵀ)`; `p(R₀ G R₀ᵀ) − p(G)`; `p(R₀ G R₀ᵀ) + p(G)` |
| `V5_EULER_LAGRANGE` | for **each** `V1` basis element: its Euler–Lagrange operator, from substituting `G_ij → ∂_i u_j` and taking the §Q1 variation with respect to each `u_j` — **componentwise, in full** |
| `V6_OPERATOR`, `V6_BASIS`, `V6_DIM` | the map `p ↦ p(R₀ · G · R₀ᵀ)` **as a linear operator on the `V1` space**; a basis for its `(−1)`-eigenspace within `V1`; and that eigenspace's dimension, counted from the basis |
| `V7_RESIDUAL` | `V6_DIM − V3_DIFFERENCE` |

⛔⛔ **EVERY BASIS IS EMITTED IN A CANONICAL FORM, AND AN ORDERING IS NOT ENOUGH.** ⚠ Two engines can
return bases of the **same** subspace related by an arbitrary invertible linear map; an index ordering
aligns positions, ⛔ not vectors, so raw bases are not comparable and a comparator would report a
disagreement that is not one.
⇒ ⭐ **Emit `V1_BASIS`, `V2_BASIS` and `V6_BASIS` as the REDUCED ROW ECHELON FORM** of the coefficient
matrix of the basis vectors in the `MONOMIAL_ORDERING` coordinates. ⭐ RREF depends only on the **span**, so
two engines finding the same subspace emit the same matrix.

⛔⛔ **Emit `V4` and `V5` for EVERY `V1` basis element, unconditionally** (corollary 4) — ⚠ emitting them
only for some elements would make the tag set itself carry the answer to `V3_DIFFERENCE`.

⚠ **`V5` is what makes `V4` mean anything.** A form can differ from `V2`'s span and still contribute
nothing to the equations of motion. ⛔ Do not state which case holds — emit the operator and let it be read.

⚠ **`R₀` is a supplied kind of object, ⛔ not a supplied matrix.** Choose any explicit real `R₀` meeting the
two conditions. ⚠ Two engines choosing different `R₀` is expected.

⭐⭐ **`V7_RESIDUAL` differences two INDEPENDENTLY COMPUTED dimensions**, which is why `V6` is asked for:
`V2` imposes invariance under the whole orthogonal group; `V6` comes from eigenspaces of one reflection
acting on `V1`. ⛔ Do not obtain `V6` by subtracting `V2` from `V1`, and ⛔ do not obtain `V2` as the
`(+1)`-eigenspace of `V6_OPERATOR`; either shortcut makes the residual identically zero and destroys the
check (corollary 3).

⭐ **`Q9` is a function of `D` alone.** `W_pkg` does not enter it. ⇒ ⭐ **compute it once per `D` and emit
the same objects under each package's tag scope.** ⚠ Corollary 4 requires the **emission** in every
package; it does ⛔ not require the **computation** to be repeated.

⚠ **Cost note.** The unknowns in the `V1`/`V2` solve grow as `D²(D²+1)/2` with `D`.
⛔ That is not a reason to narrow the `D` list, and ⛔ not a reason to substitute a floating-point solve for
the exact one. ⇒ §7: emit and flush as you go.

### Q10 · The root–coefficient Jacobian

Emit `ROOT_COEFFICIENT_JACOBIAN`: for each distinct root `r_p` and each coefficient `c_q` in
`COEFFICIENT_ORDERING`, the partial derivative `∂ r_p / ∂ c_q`, as a matrix indexed by `(p, q)` in
`ROOT_ORDERING` and `COEFFICIENT_ORDERING`.

⭐ This is one object replacing a family of yes/no questions. ⛔ Do not emit a boolean asking whether a
particular root depends on a particular coefficient, and ⛔ do not name any tag after a coefficient
(corollary 2).

⚠⚠ **On a stratum the object is ambiguous unless it is split, so emit BOTH and name them apart:**

```
STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RESTRICTED   the GENERIC Jacobian above, evaluated at STRATUM<s>_POINT
STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED   the Jacobian of the roots RECOMPUTED on that stratum,
                                                  differentiated there
```

⛔ These are different objects and ⛔ neither is "the" answer. ⚠ Where root branches coalesce on a stratum
an individual branch may not be differentiable; ⭐ emit whatever the CAS returns, including a failure
object, ⛔ never an absence.

### Q11 · Brane–bulk phase matching

**Supplied.** An unbounded medium fills the space off the brane, with `w` the coordinate normal to it. The
medium supports a **scalar sound mode only**:

```
scalar field      φ(x, w, t)  =  A · cos( Σ_m k_m x_m  +  k_w w  −  ω t )
dispersion        ω²  =  c_s0² · ( Σ_m k_m²  +  k_w² )
```

`A` is a real symbolic amplitude and `c_s0 > 0` a real constant. A brane mode and a bulk disturbance in
contact share **both** `ω` and the in-plane wavevector `k`.

⚠ **That is the entire supplied content of this section.** ⛔ No interface boundary condition, no
brane–bulk interaction operator and no continuity requirement is supplied, and ⛔ **none may be invented.**

For **each** distinct root `r`, emit:

```
ROOT<r>_KW_EQUATION      the dispersion relation as a CAS relation, with omegaSquared -> r substituted
ROOT<r>_KW_SQUARED       k_w² solved from it
ROOT<r>_KW_SIGN          the four-way symbolic test on k_w² under the §3 assumptions, token set as in Q3,
                         emitted with its operand
ROOT<r>_KW_ZERO_LOCUS_*  the full five-object locus protocol for k_w² = 0, solve variables
                         COEFFICIENT_ORDERING together with c_s0
```

⛔⛔ **DO NOT emit a residual from substituting `ROOT<r>_KW_SQUARED` back into `ROOT<r>_KW_EQUATION`.**
⚠ Measured: that residual **stays zero when the source equation is mutated and re-solved**, and moves only
when the already-solved operand is corrupted — it verifies the solver against its own output. ⭐ There is no
second route to `k_w²` here, so per corollary 3 emit the two objects and **no residual**.

⭐⭐ **The closure inventory.** Emit all four, unconditionally:

| tag | the object |
|---|---|
| `C1_EQUATIONS` | the list, as CAS relations, of every equation relating `a` to `A` that is **derivable from §1–§3 and the supplied content above** — ⭐ emit the list **as built**, whatever length it has |
| `C2_UNKNOWNS`, `C2_COUNT` | the amplitude symbols appearing in the brane and bulk solutions above — that is, the components of `a` together with `A` — as a list, and its length |
| `C3_RANK` | the rank of `C1_EQUATIONS`' coefficient matrix with respect to `C2_UNKNOWNS` |
| `C4_DIFFERENCE` | `C2_COUNT − C3_RANK` |

⛔ **Do not describe `C4_DIFFERENCE` in words, do not name any tag after what it comes out to, and do not
state in advance what `C1_EQUATIONS` will contain.**

⚠⚠ **COROLLARY 3 APPLIES SHARPLY HERE and one of the engines being replaced tripped on it.** If you emit
any projection of a brane kernel vector onto a bulk polarization direction, ⭐ **emit alongside it
`<tag>_PREMISES`, the list of premise-tag names that projection was computed under**, obtained from the
premises the computation actually consumed (corollary 5). ⛔ A projection whose value follows from the
supplied premises alone computes nothing, and emitting it without them reads as a physical result.

---

## 7 · Packages — the sweep and the controls

⛔⛔ **EVERY control re-enters at the ACTION (§2), never at a result.** A control is a **different `L`**;
everything downstream is recomputed from it by the identical code path.

⭐ **Every package emits one tag per named object, at the scope §8 gives it.** ⛔ Emission is never conditional on a value
(corollary 4). ⚠ A quantity **identical** across packages **must still be emitted in each** — that
repetition is a result, and deleting it destroys the evidence.

There are seven packages. Their complete stiffness functionals are:

```
W_MAIN             =  (μ_R/2)·S_curl  +  (B_comp/2)·S_div
W_XFORM_CURLONLY   =  (μ_R/2)·S_curl
W_XFORM_DIVONLY    =  (B_comp/2)·S_div
W_XFORM_TRACELESS  =  (μ_R/2)·S_curl  +  μ_br·S_symtl
W_XFORM_EXTRA      =  (μ_R/2)·S_curl  +  (B_comp/2)·S_div  +  (β/2)·P_D
W_XCOEF_BSCALE     =  (μ_R/2)·S_curl  +  (s·B_comp/2)·S_div
W_XCOEF_BSIGN      =  (μ_R/2)·S_curl  −  (B_comp/2)·S_div
```

| package | `D` |
|---|---|
| `MAIN` | 2, 3, 4, 5 |
| `XFORM_CURLONLY` | 2, 3, 4, 5 |
| `XFORM_EXTRA` | 2, 3, 4, 5 |
| `XFORM_DIVONLY` | 3, 4 |
| `XFORM_TRACELESS` | 3, 4 |
| `XCOEF_BSCALE` | 3 |
| `XCOEF_BSIGN` | 3 |

### ⭐⭐ `P_D` — the one term this file does not write down

`P_D` is **built from §Q9's computed output at that `D`**, by this rule:

> `P_D` ≔ the sum of the basis emitted as `V6_BASIS`, taken in the emitted ordering, with
> `G_ij → ∂_i u_j` substituted.
> ⭐ If `V6_DIM` is zero its basis is empty and `P_D` is the **zero form**; the package is still built,
> still swept, and still emits its complete tag set.

⛔⛔ **DO NOT TYPE `P_D` FOR ANY `D`.** ⭐ It must be read out of the `V6` object the engine computed, so
that corrupting the census moves the package's action — corollary 5's test, applied to the one place in
this file where a package's action is not fully written out.

⚠ **The rule names `V6`, ⛔ not "the `V1` elements outside `V2`", and the difference is not cosmetic.** A
basis of `V1` need not contain a basis of `V2`, so a `V1` basis element can fail to lie in `V2`'s span while
still carrying `V2` content. ⭐ `V6`'s eigenspace is well defined however the `V1` basis was chosen.

⭐ **Emit `P_D` as computed, once per `(package, D)`, as `PD_TERM`.** ⛔ Do not normalise it, and ⛔ do not
rescale it toward any target.

⚠ **`β` carries no sign premise and no `β ≠ 0` premise**, deliberately. ⛔ A premise forcing `β ≠ 0` would
make the zero-form case an error rather than a result.

⛔⛔ **`W_XFORM_EXTRA` MAY NOT BE ASSUMED TO PRESERVE ANY SPLIT OF THE AMPLITUDE.** ⚠ §Q4's warning about a
two-way parallel/perpendicular classification is aimed at this package. ⛔ Compute `N3`; do not classify
basis vectors.

### On the other six

⚠ **`XCOEF_BSCALE`'s `s` is DIMENSIONLESS by declaration, and that is a Q6 input** — the energy-density
requirement alone fixes only the *sum* of a scale factor's dimension and its coefficient's. ⛔ That reason
does **not** hold for `β` or `μ_br`, so neither gets a dimensionless declaration and both are `Q6`
unknowns. ⭐ **This file states no value for any of them.**

For `XCOEF_BSCALE`, `s ≠ 1` excludes the unit-scale collapse into `MAIN`. Without that premise, an
implementation could satisfy every stated premise while reproducing the baseline, leaving the coefficient
control unable to police arithmetic. A premise that keeps a control distinct is ⛔ not a premise that forces
a solver to decide.

⚠ **`XCOEF_BSIGN` flips a sign in the action, ⛔ not a premise.** `B_comp > 0` stays in force; the minus
sign is part of `W`. ⛔ **No expectation is stated here about what it moves.**

⛔⛔ **THIS FILE DOES NOT CLASSIFY ANY PACKAGE AS A "FORM" OR A "COEFFICIENT" CONTROL, AND NEITHER ENGINE
MAY.** ⚠ At S10 four packages carried an `XFORM_` prefix and only two changed the stiffness functional.
⇒ ⭐ The classification is made downstream, from the emitted stiffness functional, ⛔ never from a tag
prefix. ⛔ **No expectation is stated here about what any package does or does not move.** Run each and emit
its full tag set.

⭐ **`MAIN` swept over `D` is the primary result of the step, ⛔ not a control.**

⛔⛔ **THE RUN RECORD MUST BE OBSERVED, ⛔ NOT DECLARED.** ⭐ Accumulate a `(package, D)` pair **only after
that package has finished emitting**, and emit `RUN_PAIRS` / `SKIPPED_PAIRS` **after** the sweep, with
`SKIPPED = declared ∖ completed`.

⚠⚠ **Runtime — the rule is about OBSERVABLE PROGRESS, ⛔ not total elapsed time.**

⭐⭐ **Emit every tag as it is computed, and FLUSH.** ⛔ Do not buffer output to the end of a package, of a
dimension, or of the run. ⭐ A long run that is visibly completing cells is **acceptable**; a run producing
**no output for a long stretch** is the failure, because from outside nothing distinguishes it from a solve
that will never return.

⇒ ⭐ **Prefer completing the declared sweep over narrowing it.** ⛔ Do not drop a `(package, D)` cell to
save time, and ⛔ do not reformulate a requested object into a cheaper one. If a cell genuinely cannot
complete, ⭐ record it in `SKIPPED_PAIRS` and report it in §10 — ⛔ never silently.

⛔ Still forbidden whatever the runtime: substituting a floating-point solve for an exact one, or weakening
a requested object to make it fit.

---

## 8 · Tag grammar — ⛔ both engines must produce PARALLEL tag sets

```
<ENGINE>_S11_<PACKAGE>_D<n>_<QUANTITY>
<ENGINE>_S11_<PACKAGE>_D<n>_ROOT<r>_<QUANTITY>
<ENGINE>_S11_<PACKAGE>_D<n>_STRATUM<s>_<QUANTITY>
<ENGINE>_S11_<PACKAGE>_D<n>_STRATUM<s>_ROOT<r>_<QUANTITY>
```

- `<ENGINE>` is `WL` or `PY`. ⛔ Both prefixes are mandatory.
- `<PACKAGE>` is exactly one of the §7 names.
- `<n>` is the integer brane dimension of that run.
- `<s>` indexes `STRATUM_ORDERING`, `1`-based. `<r>` indexes `ROOT_ORDERING`, `1`-based.
- `STRATUM<s>` sits immediately after `D<n>` and before `ROOT<r>` when both scopes apply.

### ⛔⛔ ONE TAG PER NAMED OBJECT

⚠⚠ **This is the defect that made the previous pair of S11 engines uncomparable.** They were built from
two directives that decomposed the work differently: one bundled a root's value, nullity and orientation
into a **single** payload where the other emitted three. Stripping the engine prefixes, the two engines
shared **one** tag suffix between them. ⇒ ⛔ A shared grammar shape is not enough; the **decomposition**
must be shared.

⭐⭐ **Emit one tag per named object, at the scope §8 gives it.** ⛔ Do not bundle two named objects into
one payload, and ⛔ do not split one named object across several payloads. Where you emit an object for
which this file supplies no tag name, choose one and list it in the §10 report.

⭐ Emit one line per tag: `TAG: <payload>`, payload in a fully re-parseable form
(`InputForm` in Wolfram; `sympy.srepr`-safe string or `str()` of the expression in SymPy).

⭐ Emit, once per package-and-dimension: `ROOT_ORDERING`, `STRATUM_ORDERING`, `COEFFICIENT_ORDERING` and
`MONOMIAL_ORDERING`, so the orchestrator can align indices rather than assume they match.

⚠ **The engines being replaced share exactly one tag suffix, and it is `VERDICT`, which §5 deletes.** ⇒ ⛔
There is no tag namespace to preserve and no reason to echo either old vocabulary.

### ⭐ Engine-local tags — declare them, so parity is meaningful

Some tags **cannot** exist in both engines: §Q6r's registry comparison has no counterpart in the engine
that does not load the registry, and each CAS emits its own solver-condition tags. ⛔ Those are **not**
disagreements, and a parity checker that reports them as gaps trains its reader to ignore it.

⭐ Give every such tag the infix `_LOCAL_` immediately after the engine prefix — `WL_S11_LOCAL_…`,
`PY_S11_LOCAL_…` — and ⭐ emit one tag per engine listing every `_LOCAL_` name it produced.

---

## 9 · What is supplied vs. what this build tests

| supplied — ⛔ NOT tested here | tested here |
|---|---|
| the curl density's form (from S9) | everything in Q1–Q11, for every package |
| that the brane resists compression at all — `S_div` is supplied, ⛔ not derived | |
| `u` is the in-plane `D`-vector; the `h`-branon is a separate field | |
| `v₀ = 0`, no dissipation, exactly-quadratic `L`, frozen wall width | |
| `ρ_br > 0`, `μ_R > 0`, `B_comp > 0`, `Σ k_m² > 0` | |
| the bulk supports a scalar sound mode **only**, with §Q11's field, ansatz and dispersion | |

⚠ **Emit each supplied premise as its own tag** so a reader cannot mistake a passing build for having
verified one.

⚠⚠ **The frozen wall width is a supplied premise with a known consequence elsewhere**, and it is listed
here so that no tag in either engine can be read as having tested it.

---

## 10 · Report back — ⛔ under 25 lines

1. The file you wrote, its line count, and the **total number of emitted tags**.
2. The `(package, D)` pairs actually run, and any skipped, with the runtime.
3. ⭐ Anything in this specification you found **ambiguous, under-determined, ill-posed, or tautological** —
   including any requested object that **cannot be computed** from what §1–§3 supplies, any `<QUANTITY>`
   name in §6/§7 you could not emit under §8's one-tag-per-object rule, and §7's `P_D` rule if it does not
   determine a unique term. ⭐ This is wanted and is more valuable than a clean build.
4. ⛔ **Do NOT report what any value came out to be, do not interpret any result, and do not say whether
   anything "worked".** Your job ends at compute-and-print.
