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

⚠ **SUPPLIED PREMISE — unfalsifiable within this build.**
`u` is **in-plane only**. Motion *out of* the brane is a **different field** (the `h`-branon) belonging to
a different sector, and it is **not** part of `u`. ⇒ `u` has exactly `D` components, ⛔ not `D+1`.
⭐ Nothing in this build tests that split; it is inherited.

Write `G` for the `D×D` **displacement-gradient matrix**, `G_ij ≔ ∂_i u_j`. Every stiffness term below is
a quadratic form in `G`.

## 2 · The action — supplied in full, as equations

Each package is an ordered pair `(T_pkg, W_pkg)`, whose two named members are its **kinetic density** and
its **stiffness functional**. Its action density is

```
L_pkg  =  T_pkg   −   W_pkg[G]
```

The kinetic density carried by the seven existing packages is named

```
T_ISO  =  (ρ_br/2)·Σ_{j=1..D} (∂_t u_j)²
```

The complete ordered pairs are given in §7. `ρ_br > 0` is a real constant.

⚠ **`W` is a FUNCTIONAL, not a fixed coefficient times a fixed density.** ⛔ Do not assume which
coefficients a given package's `W` contains, or how many. §6 builds `COEFFICIENT_ORDERING` from the
declared additive terms of `L_pkg`, together with every inertial coefficient the kinetic density carries.

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
XKIN_ANISO:       s_ρ > 0 ,  s_ρ ≠ 1
```

⛔ `s_ρ` is not declared dimensionless and no dimension is supplied for it; it is a `Q6` unknown.

⚠ Declaring each symbol real/positive **individually is not equivalent** to asserting `Σ_m k_m² > 0`, and a
sign test performed under a different premise set is not the requested computation. ⇒ ⛔ Both engines
must carry the joint set.

⛔⛔ **NO PREMISE RELATING THE STIFFNESS COEFFICIENTS TO ONE ANOTHER IS SUPPLIED, AND NONE MAY BE ADDED.**
⚠ Do not assume any two stiffness coefficients are unequal, ordered, or generic with respect to each
other. ⭐ Where the coefficients may take special values is a **computed object** — §Q8 asks for it. A
genericity premise would delete that locus before it could be found.

⚠ **SUPPLIED PREMISES — unfalsifiable within this build:**
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

### ⭐⭐ THE ONE PERMITTED INPUT — ⛔ and the block above is NOT amended for it

⚠ The quoted rule is **shared verbatim with the other steps' specifications and is not edited here.** ⭐ It
governs what an engine may **compute by hand**. §Q8c's witness pass supplies an engine something it does
**not** compute at all: an **input**, in the same sense as the action's coefficients.

⭐⭐ **The witness input is ONE FIELD — an exact point — and it is the only thing an engine may receive
from the OTHER ENGINE'S run:**

```
the exact counterpart point
```

⛔ **Nothing else may be received from the other engine** — ⛔ no equation, matrix, root, ordering, count,
result, status, certificate, scope or coverage token.
⚠ ⭐ **The receiver is NOT told which locus the point came from, and does not need to be:** ⭐ which of its
**own** loci that point lies on is something it **computes** (§Q8c). ⇒ ⛔ no selector, ⛔ no index
translation, ⛔ no donor scope.
⚠ This says nothing about §Q6r's import of the **previous step's** `LEDGER`, which is a different channel
and is unaffected.

⚠ ⭐ **Corollary 1's test does not apply to this record**, because it is an input rather than a computed
tag — ⭐ exactly as the action's coefficients are inputs. ⛔ **The exemption does not generalise:** it names
these two fields, in this pass, and nothing else. ⛔ It is **not** a §5 live-read exemption and does not
enlarge that closed list.

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
   on **which package**, **which dimension**, **which quantity** it belongs to, and — for `WITNESS<w>`
   objects alone — **how many witness inputs §Q8c supplied**. ⚠ A value **present and
   identical across packages is a RESULT**; a value **absent** is indistinguishable from *never computed*.
   ⛔ Never de-duplicate payloads across packages. Tag **names** are unique; **payloads may legitimately
   repeat.** ⇒ ⛔ **No tag for a named object, at the scope §8 gives it, may be emitted "if" anything.**
   ⚠ ⭐ **The witness axis is admissible only because `WITNESS_ORDERING` itself is unconditional**, so
   *supplied nothing* stays distinguishable from *never ran* — ⛔ and because the count of inputs is
   **not** a payload this engine computed.

### ⛔⛔ COROLLARY 5 · A DECLARATION MUST BE WIRED TO WHAT IT DECLARES

A declaration must be read out of the live object the computation consumes, ⛔ never reconstructed from a
literal beside it. ⛔ A consumer must never be manufactured merely to make a declaration appear wired.
⚠ A control that passes by construction is worth nothing.

⭐ **The live-read exemptions are closed and are exactly these four entries:**

| exempt object | exemption scope |
|---|---|
| `PREMISE_INVENTORY` | the whole tag |
| `M_ROUTE_RESIDUAL_SCOPE` | the whole tag; its payload is one fixed token |
| `STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED.FAILURE_TOKEN` | the `.FAILURE_TOKEN` field only |
| `Q6R_RESIDUAL_SCOPE` | the whole tag; its payload is one fixed token |

⛔ No holder object may be manufactured to make an exempt token appear wired, and there is no general
exemption for metadata. The `RECOMPUTED_ROOTS`, `COEFFICIENT_ORDERING`, `DEFINING_EQUATIONS` and
`EVALUATION_POINT` fields of `STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED` remain live-read
obligations.

### ⛔⛔ NO VERDICT

⛔ **There is no `VERDICT` tag, no `PASS`, no `FAIL`, and no summary judgement anywhere in either engine.**
⛔ **A physics finding must EXIT 0.** Exit non-zero **only** on operational failure — an exception, an
unsolvable system. ⚠ Otherwise a builder iterating to exit 0 can make a genuine disagreement disappear, and
the disagreement is the most valuable output available.

⚠⚠ **This is a change from the engines being replaced.** Both as-built S11 engines end in a `VERDICT` tag,
and one renders a symbolic boolean as the typed word `TRUE`/`FALSE` with the residual discarded. ⛔ Neither
survives here. A boolean-valued test is emitted **as the CAS object the test returned**, in a structured
record with its operands beside it, and ⛔ never serialised as a host-language native boolean.

### ⛔⛔ THE LOCUS PROTOCOL — used by Q3, Q8 and Q11, and ⛔ never varied

⚠⚠ **"Solve this over the reals" is NOT a specification.** Measured on the CAS both engines will use:
a multivariate solve **ignores** the symbols' real declarations and returns solutions in the algebraic
closure; the single-variable real-domain solver **rejects** a tuple of variables; and the ordinary solver
returns the **same empty token** for a system that is identically satisfied and one that is inconsistent.
⇒ ⛔ An instruction to "solve over the reals" asks one engine for something it cannot do and makes
*checked-and-empty* indistinguishable from *true-everywhere*.

⭐⭐ **Therefore, wherever this file asks for a locus, emit ALL FIVE of these core objects and nothing less:**

| suffix | the object |
|---|---|
| `_EQUATIONS` | the defining equation system, as a list of CAS relations, **before** any solve |
| `_SOLUTION` | the solution set exactly as your CAS returns it, with the **solve variables named** and the domain left unrestricted |
| `_IDENTICALLY_SATISFIED` | the typed symbolic test of whether **every** equation in `_EQUATIONS` simplifies to zero identically in its variables |
| `_INCONSISTENT` | the typed symbolic test of whether the system admits **no** solution at all |
| `_REAL_ADMISSIBLE` | an ordered entry for **every** branch in `_SOLUTION`, carrying that branch and the typed symbolic test of whether it admits a point satisfying **every** §3 assumption |

⭐ `_IDENTICALLY_SATISFIED` and `_INCONSISTENT` are what separate the three degenerate cases, and they are
computed from `_EQUATIONS`, ⛔ never read off the solver's empty token.

⭐⭐ **THE BOOLEAN-TEST PAYLOAD FORM IS PINNED.** Each `_IDENTICALLY_SATISFIED` and `_INCONSISTENT`
payload has this ordered field sequence:

```
STATUS_TOKEN: <exactly one of PROVED_TRUE · PROVED_FALSE · UNDECIDED>
TEST_OBJECT:  <the live CAS boolean-valued object, in the §8 re-parseable CAS form>
OPERANDS:     <the live operands supplied to the test>
```

Each entry of `_REAL_ADMISSIBLE` has this ordered field sequence, in `_SOLUTION` branch order:

```
BRANCH:       <the live branch>
STATUS_TOKEN: <exactly one of ADMISSIBLE · EXCLUDED · UNDECIDED>
TEST_OBJECT:  <the live CAS boolean-valued object, in the §8 re-parseable CAS form>
OPERANDS:     <the branch and every live premise supplied to the test>
```

⭐ If `_SOLUTION` is an opaque solution object rather than an exposed branch list, emit one
`_REAL_ADMISSIBLE` entry carrying that whole object; ⛔ invent no decomposition and ⛔ emit an empty
admissibility list merely because the solver did not expose branches.

⛔ **A BRANCH IS NEVER SILENTLY DROPPED.** An `EXCLUDED` branch remains in `_REAL_ADMISSIBLE` with the
test object and operands that exclude it. An undecidable test carries `UNDECIDED`; it is ⛔ never emitted
as a bare false-valued object. ⭐ `UNDECIDED` in any protocol object is an explicit **coverage finding**:
it is neither agreement nor disagreement, and it forbids the corresponding claim of real emptiness,
real non-emptiness, branch admission, branch exclusion or complete real-locus enumeration.

⭐⭐ **THE LOCUS PROTOCOL ALSO HAS THESE FOUR UNCONDITIONAL EXTENSION OBJECTS:**

| suffix | the object |
|---|---|
| `_CANONICAL_LOCUS` | when every residual `lhs − rhs` in `_EQUATIONS` is a polynomial expression in all symbols it contains, the reduced Gröbner basis of the ideal they generate in the named solve variables, for lexicographic order in the emitted solve-variable ordering and over the exact rational-function coefficient field generated by the remaining symbols; otherwise the single token `NOT_APPLICABLE` |
| `_REAL_STATUS` | exactly one of `PROVED_EMPTY` · `PROVED_NONEMPTY` · `UNDECIDED`, for the locus over the reals under the full joint premise set in force |
| `_REAL_WITNESS` | for `PROVED_NONEMPTY`, an exact point satisfying `_EQUATIONS` and the full joint premise set; for either other status, the single token `NOT_APPLICABLE` |
| `_REAL_STATUS_OPERANDS` | the live `_EQUATIONS`, named solve variables and full joint premise set on which `_REAL_STATUS` was computed |

⚠ `_CANONICAL_LOCUS` names an object, ⛔ not an admissibility algorithm. No real-domain method is
prescribed. ⛔ Do not replace it with a complex-locus claim, and ⛔ do not turn a component returned by
one CAS into a canonicalised component list.

⛔ Emit all five core objects and all four extension objects for every locus this file requests,
unconditionally, whatever they come out to.

---

## 6 · What to compute — the question list

⛔ **This section says what to produce. It does not say what any of it equals.** Everything below is
required for **every package** defined in §7, at **every** `D` in that package's sweep.

### ⭐⭐ Ordering requirement

`§Q9` is a function of `D` alone and one package's action is built from its output. ⇒ ⭐ **Compute `Q9`
first for each `D`**, before assembling any action at that `D`.

### Q1 · The Lagrangian and the equations of motion

Build `L` in coordinate space from §2. Emit `L` expanded. Also emit the package's **declared additive
terms of `T_pkg`**, in the list `KINETIC_TERMS`, and **each additive term of `W_pkg`**, in the list
`STIFFNESS_TERMS`. The two tags have the same list shape, with one entry per declared term, so a consumer
can see which densities that package carries and with which coefficients.
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
⇒ ⭐ **The residual tests CODING CONSISTENCY ONLY.** ⭐ Emit it, and emit
`M_ROUTE_RESIDUAL_SCOPE` with the single token `CODING_CONSISTENCY_ONLY`. ⛔ Do not present it as
verifying physics.

⭐ Use `M_B` for everything downstream, and emit `M_ROUTE_USED` naming the route. Its payload is exactly
the route token `M_B`, read from the same route-selection object that supplies the matrix actually
consumed downstream (corollary 5), ⛔ not from a second literal beside it.

⭐ Also emit `M_COEFFICIENT_JACOBIAN`: for each coefficient in the package's emitted `COEFFICIENT_ORDERING`,
the matrix `∂M/∂c`, as a list in that order. ⛔ Do not describe the pattern; emit the matrices.

### Q3 · The spectrum

- Emit `DET_M`, **factored**.
- Solve `DET_M = 0` for `omegaSquared`. Emit:
  - `ROOT_MULTIPLICITY_PAIRS` — an **ordered list of `(root, multiplicity)` pairs**, with each multiplicity
    a positive integer, obtained from `DET_M` as a univariate polynomial in `omegaSquared`;
  - `ROOT_SOLUTION_SET` — the solution set as the solver returns it, carrying ⛔ **no** multiplicity claim;
  - `ROOT_COUNT_ALL` — the sum of the second fields of `ROOT_MULTIPLICITY_PAIRS`;
  - `DET_M_DEGREE` — the degree of `DET_M` in `omegaSquared`, computed from the polynomial, ⛔ not from
    the root list;
  - `ROOT_DEGREE_RESIDUAL` — `DET_M_DEGREE − ROOT_COUNT_ALL`, emitted unconditionally;
  - `ROOT_DISTINCT` — the de-duplication of the first fields of `ROOT_MULTIPLICITY_PAIRS` under the §3
    assumptions;
  - `ROOT_COUNT_DISTINCT` — the integer produced by counting `ROOT_DISTINCT`;
  - `ROOT_ORDERING` — the ordering of `ROOT_DISTINCT` that the `ROOT<r>` index refers to.
- ⭐ For each distinct root, `ROOT<r>_VALUE` and `ROOT<r>_SIGN`. ⛔ `ROOT<r>_SIGN` is a **four-way
  symbolic test** whose payload is exactly one of the tokens `POSITIVE`, `NEGATIVE`, `ZERO`, `UNDECIDED`,
  emitted **together with** the operand expression the test was applied to. ⚠ Without a pinned token set
  the two engines emit different spellings of the same determination and the comparison reports a
  disagreement that is not one.

  ⭐⭐ **THE SIGN PAYLOAD FORM IS PINNED.** It is this ordered field sequence:

  ```
  SIGN_TOKEN: <one of the four tokens above>
  OPERAND:    <the live operand expression>
  ```

  Only the field names and their order are pinned here; the operand's content is ⛔ not pinned.
- ⭐ The **root-coincidence locus in the WAVEVECTOR**, `ROOT_COINCIDENCE_K_*`: for every pair of distinct
  roots, the system `r_p − r_q = 0` with solve variables `k_1 … k_D`. Emit the **full §5 locus protocol**.
- ⭐⭐ The **root-coincidence locus in the STIFFNESS COEFFICIENTS**, `ROOT_COINCIDENCE_COEFF_*`: the same
  systems, with solve variables the coefficients in `COEFFICIENT_ORDERING`. Emit the **full §5 locus
  protocol**. ⚠ **This object does not exist at S10**: its action had one stiffness coefficient, so no
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

⚠⚠ **Whether a package action preserves a two-way parallel/perpendicular split is a question about the
complete action, regardless of whether the relevant structure is carried by its kinetic density or its
stiffness functional; nothing in this file says how any package's null directions lie relative to `k`.** ⛔
Any code path that assumes a two-way parallel/perpendicular split will report a wrong count without failing.

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
the union of those symbol sets over all terms, together with every inertial coefficient the package's
kinetic density carries, sorted by the engine's own stated rule and emitted as its own tag. ⚠ ⇒ a term
whose factor is a product contributes
**every** symbol in that product, ⛔ not the product as one unknown.

⭐ Build this inventory from the package's **DECLARED additive terms** as §7 gives them, before simplifying
or evaluating any density. A declared term whose density evaluates to the zero form still contributes the
free symbols of its factor. ⛔ Reading the inventory from the expanded action instead would let
simplification change the declared unknown set; the inventory records the coefficients supplied by the
package, including every inertial coefficient its kinetic density carries, whether or not a density
survives evaluation.

Emit:
- `COEFFICIENT_ORDERING`, and `DIM_COEFFICIENTS` — the dimension of every coefficient in it, as **closed
  functions of `D`**, obtained by requiring `L` to be an energy density on a `D`-dimensional brane;
  `DIM_EQUATIONS` — the equations solved; `DIM_SOLUTION` — the solution as returned.
- ⚠ **The energy-density requirement alone does not fix every coefficient separately** — where a term
  carries a product of undetermined coefficients it fixes only their **sum of dimensions**. ⭐ Emit the
  underdetermined solution **as returned**, with its free parameters visible, ⛔ never by silently choosing
  one. Where §7 declares a control parameter **dimensionless**, use that as a stated input.
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

⭐ **`[u]` is a PREMISE and is UNFALSIFIABLE WITHIN THIS BUILD.**

#### ⛔⛔ Q6r · THE PREVIOUS-STEP EXPORT-CHAIN COMPARISON — and ITS PROVENANCE IS PART OF THE PAYLOAD

⚠ **Engine-local: only the engine that imports the previous step's `LEDGER` emits this**, under the
`_LOCAL_` infix of §8. Import `LEDGER` from `research/pde_ledger_v3/scripts/S10_exports.py`.

⭐⭐ **The coefficient name map is explicit and closed:**

```
ρ_br → rho_br        μ_R → mu_R        B_comp → B_comp       μ_br → mu_br
β    → beta          s   → s           s_ρ  → s_rho          c_s0 → c_s0
```

For each coefficient in `COEFFICIENT_ORDERING`, map its name as above and perform the live
`LEDGER[name]['dimension_key']`, then `LEDGER[that key]['value']` lookup. Once per `(package, D)`, emit
`Q6R_RESOLVED_COEFFICIENTS` with the coefficients whose lookup resolved to a dimension row, and
`Q6R_UNRESOLVED_COEFFICIENTS` with those whose lookup did not; derive both objects from those live lookups.

For each resolved coefficient, emit the **derived** vector, the **imported** vector, their **difference**,
and provenance from both rows, named apart: `class` and `step` of the coefficient row; and `class`, `step`
and `corroborated_steps` of the resolved dimension row.
⛔ Do not manufacture a placeholder vector for an unresolved coefficient.

⭐ The imported vector is the predecessor's by construction. Emit `Q6R_RESIDUAL_SCOPE` with the single
token `CROSS_STEP_CONSISTENCY_ONLY`; the residual tests cross-step consistency only and is not an
independent operand.

⛔ If a derived and an imported vector disagree, **emit the disagreement and continue**. ⛔ Do not adjust the
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

### Q8 · ⭐⭐⭐ Generic rank — and COMPONENT-SCOPED Q3 AND Q4 WHERE IT DROPS

⚠ `MatrixRank` / `Matrix.rank()` return the **GENERIC** rank — the rank away from measure-zero loci in the
symbols. ⇒ ⛔⛔ **EVERY PACKAGE-AND-DIMENSION NUMBER Q4 PRODUCES IS GENERIC; IT DOES NOT DETERMINE
THE CORRESPONDING OBJECT ON A COMPONENT.**

⛔⛔ **AND AT S11 "THE SYMBOLS" INCLUDES THE STIFFNESS COEFFICIENTS.** S10 swept loci in the wavevector
only, because its action carried one stiffness coefficient. ⭐ Here a locus may live in the coefficients, in
the wavevector, or in both, and §3 supplies no premise excluding any of them.

**Q8a — the rank-drop loci.** For each `M_r`, take `ρ` to be the rank **your own computation returned**,
form **all** `ρ × ρ` minors, and emit `ROOT<r>_RANK_DROP_MINORS` (the minor list) plus the **full §5 locus
protocol** for each of three solve-variable sets, under the suffixes:

```
ROOT<r>_RANK_DROP_K_*         solve variables: k_1 … k_D
ROOT<r>_RANK_DROP_COEFF_*     solve variables: COEFFICIENT_ORDERING
ROOT<r>_RANK_DROP_JOINT_*     solve variables: both, jointly
```

For the `(D+1)×D` matrix formed by stacking `M_r` on top of `kᵀ`, take `σ` to be the stacked rank the
engine's own `ROOT<r>_N3_STACKED_RANK` returned. Form **all** `σ × σ` minors and emit them as
`ROOT<r>_STACKED_DROP_MINORS`, plus the **full §5 locus protocol** for the same three
solve-variable sets, under the suffixes:

```
ROOT<r>_STACKED_DROP_K_*       solve variables: k_1 … k_D
ROOT<r>_STACKED_DROP_COEFF_*   solve variables: COEFFICIENT_ORDERING
ROOT<r>_STACKED_DROP_JOINT_*   solve variables: both, jointly
```

⚠ **If the returned source rank — `ρ` or `σ` — is `0`** there are no corresponding square minors and
the source matrix is identically zero. ⭐ Emit the corresponding `ROOT<r>_RANK_DROP_MINORS` or
`ROOT<r>_STACKED_DROP_MINORS` as the **empty list** and the locus `_EQUATIONS` as the empty list; the protocol's
`_IDENTICALLY_SATISFIED` then carries the meaning. ⛔ Do not emit a `0×0` determinant, which some CAS report
as `1`.

**Q8b — ⭐⭐ STRATA, defined.** All three `ROOT<r>_RANK_DROP_K_*`,
`ROOT<r>_RANK_DROP_COEFF_*`, `ROOT<r>_RANK_DROP_JOINT_*` solve-variable sets, all three
`ROOT<r>_STACKED_DROP_K_*`, `ROOT<r>_STACKED_DROP_COEFF_*`, `ROOT<r>_STACKED_DROP_JOINT_*`
solve-variable sets, and both Q3 families, `ROOT_COINCIDENCE_K_*` and `ROOT_COINCIDENCE_COEFF_*`, feed the
candidate pool. A **stratum** is a candidate component whose `_REAL_ADMISSIBLE` entry is `ADMISSIBLE` under
the region allowed by §3 and for which the engine emits the exact point below. A candidate with an
`UNDECIDED` admissibility status remains an explicit locus-protocol coverage finding and is ⛔ not
silently discarded or promoted to a stratum. Emit the admitted components your CAS returns, without
canonicalising them, each with its sources, source-locus identities, defining equations, free parameters
and chosen point, plus their ordering:

```
STRATUM_ORDERING              the list of components, in the order the STRATUM<s> index refers to
STRATUM<s>_SOURCES             the tokens for that component, drawn from exactly `RANK_DROP` ·
                               `STACKED_DROP` · `ROOT_COINCIDENCE`, in that order, each appearing at most once
STRATUM<s>_SOURCE_LOCUS_TAGS   the ordered §5 locus base-tag identities whose branches produced the component
STRATUM<s>_DEFINING_EQUATIONS the equations cutting out that component
STRATUM<s>_FREE_PARAMETERS     the free parameters in which this component's restricted objects are expressed
STRATUM<s>_POINT              ⭐ an explicit point on it satisfying every §3 assumption
```

⭐ The point is obtained from those same defining equations, so no independent route exists. Per §5
corollary 3, emit `STRATUM<s>_DEFINING_EQUATIONS` and `STRATUM<s>_POINT` and ⛔ emit no residual between
them.

⭐⭐ **First, compute the complete `Q3` and `Q4` objects ON THE COMPONENT.** Restrict the engine's own
`M` and the objects derived from it to `STRATUM<s>_DEFINING_EQUATIONS`, leave
`STRATUM<s>_FREE_PARAMETERS` free, and emit the results under the existing `STRATUM<s>` and
`STRATUM<s>_ROOT<r>` scopes. ⛔ None of these component-scoped payloads may be obtained by substituting
`STRATUM<s>_POINT`.

⚠⚠ **THE TWO ENGINES MAY DESCRIBE ONE COMPONENT IN DIFFERENT VARIABLES, AND THIS FILE DOES NOT PIN WHICH.**
⭐ Each engine emits `STRATUM<s>_FREE_PARAMETERS` as **what it actually retained**, and its component-scoped
**symbolic** payloads in those parameters.
⇒ ⭐⭐ **Those symbolic payloads are INSPECTION-ONLY: ⛔ they are not cross-engine comparison rows**, because
two faithful engines can eliminate different variables and differ for a reason that is not physics.
⭐⭐ **The COUNTS and their STATUSES ARE comparison rows** — a count is invariant under which variable was
eliminated, and a status is a property of the component, ⛔ not of its description.
⚠ **The certificate and the change locus are NOT comparison rows**: both are **expressed in** the free
parameters, so two faithful engines write them differently. ⭐ They are emitted, and the orchestrator
aligns them exactly as it aligns strata — ⛔ they are not differenced.
⚠ ⭐ **Pinning an elimination is not attempted**: no rule that names one variable is valid on every
component a CAS can return, and one that names the wrong variable **deletes a branch.**

⭐ **Every component-scoped Q3/Q4 tag whose payload is an integer count** — including a degree, rank,
nullity or basis count — has three companion object families. Write `<COUNT>` for that tag's quantity name:

```
STRATUM<s>_<COUNT>_STATUS                    exactly CONSTANT · VARIES · UNDECIDED
STRATUM<s>_<COUNT>_CONSTANCY_CERTIFICATE     an exact CAS certificate over every FREE_PARAMETER when
                                             status is CONSTANT; NOT_APPLICABLE otherwise
STRATUM<s>_<COUNT>_CHANGE_LOCUS_*            every §5 protocol tag for the sub-locus in this
                                             component where the count changes when status is VARIES;
                                             every member is NOT_APPLICABLE otherwise
```

⭐⭐ **AND `STRATUM<s>_<COUNT>` ITSELF HAS ONE PAYLOAD FORM IN ALL THREE STATUSES** — ⛔ never a bare
integer, ⛔ never a different shape per status:

```
STATUS_TOKEN: <the same token emitted under _STATUS>
VALUE:        <the count on the component where it is defined there, whatever the status;
               the single token NOT_DEFINED_ON_COMPONENT where the engine did not obtain one>
```

⚠ ⭐ Under `VARIES` the answer is the **pair** — the `VALUE` the engine obtained together with
`_CHANGE_LOCUS_*` — ⛔ neither field alone. ⚠ Where the engine obtains no single count off the change
sub-locus, `VALUE` is `NOT_DEFINED_ON_COMPONENT` and the change locus carries the object; ⛔ do not
manufacture one.

For a root-scoped count, insert the same suffixes after
`STRATUM<s>_ROOT<r>_<COUNT>`. All named tags above are emitted unconditionally. The change sub-locus is
reported only as an object; ⛔ do not promote it to another stratum, recurse on it or canonicalise its
components. ⛔ A bare generic integer is not a component-scoped answer.

⭐ Also emit `STRATUM<s>_COMPONENT_Q3_Q4_COVERAGE`, with exactly one token:
`COMPLETE_COVERAGE` when every requested component object was built and every component count has a
`CONSTANT` or `VARIES` status; `INCOMPLETE_COVERAGE` otherwise. If the engine cannot build a component-wide
Q3/Q4 object, its count status is `UNDECIDED` and its `VALUE` field is `NOT_DEFINED_ON_COMPONENT` — ⭐ the
same one record as every other status, ⛔ never a second payload shape — and ⛔ no component-level count may
be copied or inferred from a point evaluation.

⭐ **Then recompute and emit the complete `Q3` and `Q4` tag sets at `STRATUM<s>_POINT` as POINT-LOCAL
EVIDENCE.** Insert `POINT_EVIDENCE_` immediately before each unscoped quantity name and immediately after
`ROOT<r>_` for each root-scoped name; for example, `STRATUM<s>_POINT_EVIDENCE_ROOT_COUNT_ALL` and
`STRATUM<s>_ROOT<r>_POINT_EVIDENCE_N2_RANK`. These tags are evidence at the emitted point, ⛔ never the
component's answer. An engine that emits `INCOMPLETE_COVERAGE` still emits this point-evidence set, but the
set does not cure or narrow that coverage finding.

⛔⛔ **THIS IS NOT BOOKKEEPING.** A generic rank does not determine the rank on an exceptional stratum.
Reporting a locus without recomputing on it leaves the spectrum and mode count there untested.

⚠ **A single point does not characterise a positive-dimensional component.** ⭐
`STRATUM<s>_DEFINING_EQUATIONS` and `STRATUM<s>_FREE_PARAMETERS` carry the component; the point and its
point-evidence tags are emitted so their strictly local scope is visible.

⭐ If a package has **no** stratum, emit `STRATUM_ORDERING` as the **empty list** — ⛔ do not omit the tag
(corollary 4). ⭐ Read it only with the locus protocol: the typed branch and real-status objects distinguish
proved exclusion from undecided coverage and from a computation that was never emitted.

⭐ Stratum indices are aligned by the orchestrator from the emitted defining equations; they are never
assumed to correspond across engines.

**Q8c — ⭐⭐ WITNESS EXCHANGE.** ⭐⭐ **A mandatory two-pass protocol.** ⚠ It exists because the two CASes
have unequal real-domain capability (§5), so one engine can produce an admissible point where the other
cannot decide. ⇒ ⭐ **That asymmetry becomes a computation instead of an incomparability.**

1. ⭐ **Native pass.** Both engines run and emit every object this file requests. ⛔ Neither receives
   anything from the other.
2. ⭐ **Supply.** The **orchestrator** — ⛔ never either engine — collects **every exact point either engine
   emitted**, for **any** locus or stratum, and supplies each of them to the other engine as the §4 witness
   input. ⛔ **There is no pairing condition and no matching test.** ⇒ ⛔ **An empty `STRATUM_ORDERING`,
   an undecided admissibility, or a differently-presented equation cannot suppress the exchange**, because
   nothing has to match for a point to be handed over.
3. ⭐ **Witness pass.** Both engines run again, rebuild their own objects, consume only the supplied points,
   and emit the objects below. ⛔ No engine reads the other engine's output.

⭐⭐ **The receiver is not told where the point came from. It COMPUTES where the point lies among its own
objects.** ⛔ It imports no equation, matrix, root, ordering, count, result, status, scope or certificate,
and ⛔ adjusts no local object toward the counterpart.

⭐ Emit `WITNESS_ORDERING` unconditionally in the witness pass, with one entry per point supplied.
⛔ **If any point was supplied, an empty ordering is non-compliant.** For each `WITNESS<w>`, emit:

```
WITNESS<w>_RECEIVED_POINT                  the exact point supplied
WITNESS<w>_OWN_LOCUS_RESIDUALS             one entry for EVERY locus this engine emitted at this package
                                           and dimension: that locus's own base tag, and its own
                                           `lhs − rhs` for every equation, evaluated at RECEIVED_POINT
WITNESS<w>_POINT_COVERAGE                  exactly COMPLETE_POINT · INCOMPLETE_POINT, from testing whether
                                           the point assigns every symbol this engine's own `M` depends
                                           on, with omegaSquared left as the §Q3 solve indeterminate
WITNESS<w>_OWN_M_EVALUATED                 this engine's own live `M`, evaluated at RECEIVED_POINT
WITNESS<w>_POINT_EVIDENCE_<QUANTITY>       the complete Q3/Q4 point-evidence set computed from that own M
WITNESS<w>_ROOT<r>_POINT_EVIDENCE_<QUANTITY>  its root-scoped members
```

⭐⭐ **`WITNESS<w>_OWN_LOCUS_RESIDUALS` is the object that replaces a supplied selector.** ⭐ It states, by
computation, which of this engine's own loci the point satisfies — ⛔ so no engine needs to be told, and
⛔ no two engines need to present a locus the same way for the exchange to work.
⭐ The point is independent of the receiving engine's construction, so these residuals are genuine operands
and are required by §5 corollary 3.
⛔ If `WITNESS<w>_POINT_COVERAGE` is `INCOMPLETE_POINT`, the evaluated objects remain partly symbolic and
their counts are generic — ⭐ emit them, and ⛔ do not read them as point values.

⚠ ⭐ **Pairing a witness row with the counterpart's own row at the same point is the ORCHESTRATOR's job**,
from `WITNESS<w>_RECEIVED_POINT` and §8's orderings — ⭐ the same alignment it already performs for stratum
and root indices. ⛔ It is not encoded in any tag an engine receives.

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
STRATUM<s>_POINT_EVIDENCE_ROOT_COEFFICIENT_JACOBIAN_RESTRICTED
                                                  the GENERIC Jacobian above, evaluated at STRATUM<s>_POINT
                                                  ⚠ point-local EVIDENCE under §Q8b's rule, ⛔ never the
                                                  component's Jacobian
STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED   the pinned failure object below; no shared derivative of
                                                  the roots recomputed on that stratum is specified
```

⚠⚠ **This specification supplies neither tangent coordinates on a stratum nor an off-stratum extension.**
Therefore no shared construction exists: `STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED` is **always**
the explicit structured failure object below. An engine must ⛔ not substitute a construction of its own,
differentiate eliminated constants, emit a zero matrix or omit the tag. This missing construction is a
stated limitation of the specification, ⛔ not an engine defect.

⭐⭐ **THE FAILURE PAYLOAD FORM IS PINNED.** It is the following ordered field sequence, with the field
names and single condition token exactly as written:

```
FAILURE_TOKEN:        MISSING_TANGENT_COORDINATES_AND_OFF_STRATUM_EXTENSION
RECOMPUTED_ROOTS:     <the live COMPONENT-SCOPED STRATUM<s>_ROOT_DISTINCT object;
                       ⛔ never its POINT_EVIDENCE_ counterpart>
COEFFICIENT_ORDERING: <the live coefficient ordering>
DEFINING_EQUATIONS:   <the stratum's live defining equations>
EVALUATION_POINT:     <the stratum's live evaluation point>
```

Only the field names, their order and the token spelling are pinned; the four dependency contents are
⛔ not pinned.

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
ROOT<r>_KW_SIGN          the four-way symbolic test on k_w² under the governing assumption set stated
                         immediately below, using Q3's pinned token set and payload form
ROOT<r>_KW_ZERO_LOCUS_*  the full §5 locus protocol for k_w² = 0, solve variables
                         COEFFICIENT_ORDERING together with c_s0
```

⭐ The governing assumption set for `ROOT<r>_KW_SIGN` and every `_REAL_ADMISSIBLE` test in
`ROOT<r>_KW_ZERO_LOCUS_*` is the §3 joint set joined with §Q11's supplied bulk premises that `A` is real
and `c_s0` is real and positive.

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

⭐⭐ **These four objects measure the deficiency of the content supplied by §1–§3 and §Q11, ⛔ not a
property of the physical brane–bulk interface.** The step record must classify all four on that basis, so
`C4_DIFFERENCE` cannot be read as a discovered interface result.

⛔ **Do not describe `C4_DIFFERENCE` in words or name any tag after what it comes out to.**

---

## 7 · Packages — the sweep and the controls

⛔⛔ **EVERY control re-enters at the ACTION (§2), never at a result.** A control is a **different `L`**;
everything downstream is recomputed from it by the identical code path.

⭐ **Every package emits one tag per named object, at the scope §8 gives it.** ⛔ Emission is never conditional on a value
(corollary 4). ⚠ A quantity **identical** across packages **must still be emitted in each** — that
repetition is a result, and deleting it destroys the evidence.

There are eight packages. The eighth kinetic density is

```
T_ANISO  =  (ρ_br/2)·Σ_{j=2..D} (∂_t u_j)²   +   (s_ρ·ρ_br/2)·(∂_t u_1)²
```

Its distinguished axis is `j = 1`; no relation between that axis and `k` is supplied, and `k` stays
symbolic. `T_ANISO` has exactly two declared additive terms: the `j = 2..D` sum is one term and the `j = 1`
term is the other. ⛔ Neither is expanded before the §Q6 coefficient inventory is built.

The complete stiffness functionals are:

```
W_MAIN             =  (μ_R/2)·S_curl  +  (B_comp/2)·S_div
W_XFORM_CURLONLY   =  (μ_R/2)·S_curl
W_XFORM_DIVONLY    =  (B_comp/2)·S_div
W_XFORM_TRACELESS  =  (μ_R/2)·S_curl  +  μ_br·S_symtl
W_XFORM_EXTRA      =  (μ_R/2)·S_curl  +  (B_comp/2)·S_div  +  (β/2)·P_D
W_XCOEF_BSCALE     =  (μ_R/2)·S_curl  +  (s·B_comp/2)·S_div
W_XCOEF_BSIGN      =  (μ_R/2)·S_curl  −  (B_comp/2)·S_div
```

The complete package pairs, in `(kinetic density, stiffness functional)` order, are:

```
MAIN             =  (T_ISO,   W_MAIN)
XFORM_CURLONLY   =  (T_ISO,   W_XFORM_CURLONLY)
XFORM_DIVONLY    =  (T_ISO,   W_XFORM_DIVONLY)
XFORM_TRACELESS  =  (T_ISO,   W_XFORM_TRACELESS)
XFORM_EXTRA      =  (T_ISO,   W_XFORM_EXTRA)
XCOEF_BSCALE     =  (T_ISO,   W_XCOEF_BSCALE)
XCOEF_BSIGN      =  (T_ISO,   W_XCOEF_BSIGN)
XKIN_ANISO       =  (T_ANISO, W_MAIN)
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
| `XKIN_ANISO` | 2, 3, 4, 5 |

### ⭐⭐ `P_D` — the one term this file does not write down

`P_D` is **built from §Q9's computed output at that `D`**, by this rule:

> `P_D` ≔ the sum of the basis emitted as `V6_BASIS`, taken in the emitted ordering, with
> `G_ij → ∂_i u_j` substituted.
> ⭐ If `V6_DIM` is zero its basis is empty and `P_D` is the **zero form**; the package is still built,
> still swept, and still emits its complete tag set.

⛔⛔ **DO NOT TYPE `P_D` FOR ANY `D`.** ⭐ It must be read out of the `V6` object the engine computed, so
that corrupting the census moves the package's action — corollary 5's wiring requirement, applied to the one place in
this file where a package's action is not fully written out.

⚠ **The rule names `V6`, ⛔ not "the `V1` elements outside `V2`", and the difference is not cosmetic.** A
basis of `V1` need not contain a basis of `V2`, so a `V1` basis element can fail to lie in `V2`'s span while
still carrying `V2` content. ⭐ `V6`'s eigenspace is well defined however the `V1` basis was chosen.

⭐ **Emit `P_D` as computed, once per `(package, D)`, as `PD_TERM`.** ⛔ Do not normalise it, and ⛔ do not
rescale it toward any target.

⚠ **`β` carries no sign premise and no `β ≠ 0` premise**, deliberately. ⛔ A premise forcing `β ≠ 0` would
make the zero-form case an error rather than a result.

⛔⛔ **NO PACKAGE ACTION MAY BE ASSUMED TO PRESERVE ANY SPLIT OF THE AMPLITUDE.** ⚠ §Q4's warning applies
to the complete action, regardless of whether the relevant structure is carried by its kinetic density or
stiffness functional. ⛔ Compute `N3`; do not classify basis vectors.

### On the other seven

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
⇒ ⭐ The classification is made downstream, from the emitted action, ⛔ never from a tag
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
<ENGINE>_S11_<PACKAGE>_D<n>_WITNESS<w>_<QUANTITY>
<ENGINE>_S11_<PACKAGE>_D<n>_WITNESS<w>_ROOT<r>_<QUANTITY>
```

- `<ENGINE>` is `WL` or `PY`. ⛔ Both prefixes are mandatory.
- `<PACKAGE>` is exactly one of the §7 names.
- `<n>` is the integer brane dimension of that run.
- `<s>` indexes `STRATUM_ORDERING`, `1`-based. `<w>` indexes `WITNESS_ORDERING`, `1`-based.
- `<r>` is `1`-based and indexes the `ROOT_ORDERING` at its own package, component, point-evidence or witness
  scope. `POINT_EVIDENCE_` is part of `<QUANTITY>`; it does not create or move an index.
- `STRATUM<s>` or `WITNESS<w>` sits immediately after `D<n>` and before `ROOT<r>` when both scopes apply.

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

⭐ Emit, once per package-and-dimension: `ROOT_ORDERING`, `STRATUM_ORDERING`, `WITNESS_ORDERING`,
`COEFFICIENT_ORDERING` and `MONOMIAL_ORDERING`, so the orchestrator can align indices rather than assume
they match.

⚠ **The engines being replaced share exactly one tag suffix, and it is `VERDICT`, which §5 deletes.** ⇒ ⛔
There is no tag namespace to preserve and no reason to echo either old vocabulary.

### ⭐ Engine-local tags — declare them, so parity is meaningful

Some tags **cannot** exist in both engines: §Q6r's previous-step export-chain comparison has no counterpart
in the engine that does not import that chain's `LEDGER`, and each CAS emits its own solver-condition tags.
⛔ Those are **not** disagreements, and a parity checker that reports them as gaps trains its reader to
ignore it.

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
| the isotropy of the kinetic form | |
| the bulk supports a scalar sound mode **only**, with §Q11's field, ansatz and dispersion | |
| the completeness of the stratum enumeration is not established by this build | the component-scoped Q3/Q4 status, coverage and point-evidence objects; the typed locus real status and counterpart-witness evaluations |

### ⭐⭐ Premise inventory — one named object

Each engine emits **one** engine-local tag per `(package, D)`:
`<ENGINE>_S11_LOCAL_<PACKAGE>_D<n>_PREMISE_INVENTORY`. It is engine-local and is ⛔ not a cross-engine
comparison row. Each engine lists the supplied premises its run carried, in whatever form it holds them.

`PREMISE_INVENTORY` is **one named object**, so emitting it as one tag is not the bundling §8 forbids. Its
entries are declarations, ⛔ not computed objects, and the previous per-premise decomposition is what
diverged.

This lets a reader see every premise the build did not test.

⚠ **Corollary 5's closed list exempts this whole tag from the live-read requirement.** ⭐ Several supplied
premises are qualitative or assert an **absence** — there is no live CAS object to read them from, and ⛔ one
must not be manufactured. ⇒ ⭐ list them in whatever form the engine holds them.

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
