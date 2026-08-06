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

⭐ **`u` is a DISPLACEMENT, so `[u] = length = (1, 0, 0)` in the `(length, time, mass)` slots.** ⛔ This is
a **supplied premise**, ⛔ not something to infer — §6's dimension solve is under-determined without it, and
two engines left to choose independently can pick differently and disagree for no physical reason.
⭐ Emit it as a premise tag.

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
`(∂_i u_j − ∂_j u_i)²` — substituting `∂ → i k` multiplies that square by `i² = −1` and **flips the sign of
the entire stiffness term**. ⇒ two engines free to choose the convention would differ by a sign and a
factor, for no physical reason.

⭐ An overall positive constant common to the whole quadratic form changes **nothing** requested in §6.
⭐ Emit whatever factor you obtain rather than normalising it away.

Write `ω²` as a **single symbol** (`omegaSquared`) throughout, ⛔ never as `ω` squared; solving for `ω`
introduces sign branches that are not separate physical modes.

**Assumptions in force everywhere — ⛔ BOTH ENGINES MUST CARRY THE IDENTICAL SET, and must carry it as a
JOINT assumption, not merely as per-symbol declarations:**

```
ρ_br > 0 ,  μ_R > 0 ,  Σ_m k_m² > 0 ,  every k_m real ,  every a_j real ,  D a positive integer

package-domain premises, joined to the common set when that package is selected:
XFORM_ANISO:    s_ρ > 0 ,  s_ρ ≠ 1
XCOEF_SCALE:    s > 0   ,  [s]   = (0, 0, 0) ,  s ≠ 1
```

⚠ Declaring each symbol real/positive **individually is not equivalent** to asserting `Σ_m k_m² > 0`, and a
sign test performed under a different premise set is not the requested computation. ⇒ ⛔ Both engines
must carry the joint set.

⚠ **SUPPLIED PREMISES — unfalsifiable within this build. Emit each as a premise tag:**
- the background brane is **unstrained and at rest**: `v₀ = 0`, so there are **no convective terms** and no
  advective rate `k·v₀`;
- **no dissipation**: `L` contains **no term odd in `∂_t`** — no first time derivative, no relaxation time,
  no memory kernel;
- **linear response / small oscillations**: `L` is exactly quadratic in `u`, so no amplitude-dependent
  restoring force is retained.

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

### Q2 · The dynamical matrix, by two routes — ⚠ and BE HONEST ABOUT WHAT THAT DOES AND DOES NOT TEST

`M` is the `D×D` coefficient matrix of the plane-wave amplitude equations, `M · a = 0`.

- **Route A — through the equations of motion.** Substitute the §3 ansatz into the Q1 Euler–Lagrange
  system, strip the common trigonometric factor, and extract the coefficient matrix. Call it `M_A`.
- **Route B — through the quadratic form.** Period-average the §2 quadratic density on the §3 ansatz and
  take `M_B[[i,j]] = ∂²⟨L⟩/∂a_i ∂a_j`. Call it `M_B`.

⭐ **Emit `M_A`, `M_B`, and `M_A − M_B`** — all three (clause 2).

⛔⛔ **THESE ARE NOT TWO INDEPENDENT PHYSICS DERIVATIONS.** Both routes are built from the same action.
⇒ ⭐ **The residual tests CODING CONSISTENCY ONLY** — that both routes were built from that action.
⭐ Emit it, and emit a tag saying that is what it tests. ⛔ Do not present it as verifying physics.

⚠ If `M_A` and `M_B` differ by an overall scalar, ⛔ do not normalise one to the other — emit both, the
difference, and the ratio `M_A[[1,1]]/M_B[[1,1]]`.

⭐ Use `M_B` for everything downstream, and emit a tag naming the route used.

### Q3 · The spectrum

- Emit `Det[M]`, **factored**.
- Solve `Det[M] = 0` for `omegaSquared`. Emit the **full solution set**, then the **distinct** roots after
  de-duplication under the assumptions. Emit **how many distinct roots** there are, as an integer produced
  by counting the list.
- ⭐ For each distinct root, emit its **SIGN** under the assumptions (`Sign[root]`, or a three-way
  positive/zero/negative symbolic test). ⚠ Without this an exponentially growing mode carries the same
  signature as a propagating wave.
- ⭐ Emit the **root-coincidence locus**: for every pair of distinct roots `r_p`, `r_q`, solve
  `r_p − r_q = 0` **for the wavevector components `k_1 … k_D`, over the REALS**, and emit the solution set.

⛔⛔ **EVERY "SOLVE THE LOCUS" INSTRUCTION IN THIS FILE NAMES ITS SOLVE VARIABLES AND ITS DOMAIN, AND THE
DOMAIN IS ALWAYS THE REALS.** ⚠ A solve with no variable named returns solutions in whatever symbol the
CAS picks, and an unrestricted domain permits **complex** wavevectors — which §3 forbids. ⛔ Two engines
must not use different solve variables or domains for this computation.

⭐ Emit the locus **as found over the reals**, ⭐ and separately emit whether it intersects the region
allowed by the §3 assumptions.

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
| **N7** | ⭐ `nu_basis` ≔ the **number of vectors the null-space routine returns**, and `nu_basis − nu` |

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

⚠ **N7 differences two DIFFERENT ALGORITHMS** — the null-space routine's basis count against the rank
routine's implied count — which is what lets it fail. ⛔ Do not "simplify" it to `rank + nu − D`: with `nu`
defined as `D − rank` that is identically zero for any rank at all.

### Q5 · Dispersion — by SCALING, ⛔ never by substitution

For each root `r(k)`, introduce a positive scalar symbol `lambdaScale` and emit:

```
r( lambdaScale · k )          the root at the scaled wavevector
r( k )                        the root at the original wavevector
r( lambdaScale · k ) / r( k )        ⭐ fully simplified — THE OBJECT THAT ANSWERS THIS
```

⭐⭐ **Emit the RATIO and let whatever power of `lambdaScale` appears, appear.** ⛔ **Do not difference the
scaled root against `lambdaScale` raised to a stated power** — naming a power tells you the answer before
you compute it.

⭐ If the ratio is a pure power of `lambdaScale`, also emit that **exponent** as an object obtained by
`Solve`/`solve` or by `Exponent`/`degree` — ⛔ never typed.

⛔⛔ **EMIT THE ENTIRE Q5 TAG SET FOR EVERY ROOT, UNCONDITIONALLY — including when the ratio is
undefined.** ⛔ **A missing tag is indistinguishable from a quantity never computed**, which is the defect
this whole specification exists to remove.

⇒ ⭐ **The TAG SET is fixed; the PAYLOAD varies.** Where the ratio is undefined, the payload is an
**explicit marker object** saying so, ⛔ never an absence, and the remaining Q5 tags for that root are
still emitted with whatever they legitimately hold.

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

Compute and emit:
- **the dimension of every inertial coefficient and every stiffness coefficient the package's action
  actually contains**, as **closed functions of `D`**, obtained by requiring `L` to be an energy density on
  a `D`-dimensional brane. Emit the equations solved and the solution.
  ⚠ ⛔ **Which coefficients exist is PACKAGE-DEPENDENT — do not assume `ρ_br` appears.** A package whose
  kinetic term carries per-component densities has **those** as its inertial coefficients; solve for each.
- ⚠ **The energy-density requirement alone does not fix every coefficient separately** — where an action
  contains a product of two undetermined coefficients it fixes only their **sum of dimensions**. ⭐ Emit
  the underdetermined solution **as returned**, with its free parameters visible, ⛔ never by silently
  choosing one. Where §7 declares a control parameter **dimensionless**, use that as a stated input and
  emit it as a premise tag.
- ⭐ For **each root** `r`, the dimension of `r / (Σ_m k_m²)`, computed **from the emitted expression for
  `r` by walking its tree** — ⛔ not by differencing coefficient dimensions.
- The dimension of **every** emitted dimensionful expression, and its **homogeneity** — a boolean per
  expression, with the per-term dimension vectors emitted alongside so a failure is readable.

### ⛔⛔ Q6d · MAKE THE HOMOGENEITY CHECK'S OWN VACUITY VISIBLE

⭐ **Why, and it is structural, ⛔ not a coding bug:** the coefficient dimensions are **solved for** by
requiring each action term to equal the energy density. Homogeneity is then evaluated **under that same
solution** ⇒ true by construction. ⚠ **And it cannot be repaired by adding a cleverer check built from the
coefficients:** `[u]` enters **every** coefficient with the same weight, so it **cancels identically** in
any ratio or difference of them — including in `[ω²]`, which is therefore blind to `[u]` as well.

⇒ ⭐⭐ **Emit the diagnostic instead of pretending the check is one:**

- the **number of independent dimension equations** formed, as an integer;
- the **number of unknown coefficient dimensions** solved for, as an integer — ⛔⛔ **counted from the
  PACKAGE'S OWN ACTION, never from a fixed list of symbols.**
- ⭐ their **difference**;
- ⭐ and, from the solve itself, whether the system was **over-**, **exactly-**, or **under-determined**.

⭐⭐ **When the difference is `≤ 0` the homogeneity booleans for the action are VACUOUS** — there were never
more constraints than unknowns. ⭐ Emit that fact as its own tag so it appears in the output next to the
booleans it qualifies. ⛔ Do **not** suppress the booleans; ⭐ label them.

⚠ **Homogeneity retains its value elsewhere** — on expressions the solve never touched (roots, null-space
objects, minors, the Q7 comparison), an inhomogeneity **is** a real finding. ⇒ ⭐ Emit homogeneity for the
**solved** action terms and for **everything else** under **distinguishable tag names**, so a consumer can
tell which class it is reading.

⚠ Homogeneity is also **blind to a wrong dimensionless coefficient.** ⇒ it is a layer, ⛔ not the answer.

⭐ **`[u]` is a PREMISE and is UNFALSIFIABLE WITHIN THIS BUILD.** ⛔ Emit it as a premise tag, and ⛔ do not
let any tag name suggest it was verified. ⚠ The **one** place it can genuinely fail is the comparison of
**derived** against **registry-declared** coefficient dimensions — ⭐ which only engine 2 can perform.

### Q7 · The `D = 3` comparison against the ordinary curl

The stiffness is written in general `D`, because the three-vector cross product exists only at `D = 3`.
At `D = 3`, Q7 compares exactly one pair: this package's stiffness density against the ordinary
curl-squared. Emit both operands and their difference. This is a form-and-normalisation comparison.
⛔ It is not evidence for the mode count, which is computed
in every `D` and does not depend on this comparison.

Use these objects:

```
g_ij   := independent symbols standing for ∂_i u_j
S_pkg  := the stiffness-density object this package's action uses, with the g_ij substituted in
c_i    := Σ_{j,k} ε_ijk g_jk
emit:     S_pkg  ·  c·c  ·  (S_pkg − c·c)                 all three
```

The CAS must compute `c_i` from the Levi-Civita definition and compute `c·c` from that result. ⛔ Do not
write out its expanded polynomial.

The emitted `S_pkg` must be the stiffness-density object from which that package's action was assembled,
⛔ not an equal-valued rebuild. Mutating **the stiffness density in the action**, with the package selector
held **fixed**, must move the emitted `S_pkg`.

That mutation requirement establishes dependence on the action; it does ⛔ **not** establish that the
emitted density is the one from which the action was assembled. Nothing inside a single engine does.
Provenance rests on the cross-engine comparison of this emitted tag, ⛔ not on the within-engine residual.

### Q8 · ⭐⭐⭐ Generic rank — and RE-RUNNING Q4 ON EVERY STRATUM WHERE IT FAILS

⚠ `MatrixRank` / `Matrix.rank()` return the **GENERIC** rank — the rank away from measure-zero loci in the
symbols. ⇒ ⛔⛔ **EVERY NUMBER Q4 PRODUCES IS A GENERIC NUMBER, AND SAYING SO IS NOT ENOUGH.**

**Q8a — find the strata.** For each `M_r`, emit:
- the locus where the generic rank **drops** — take `r` to be the rank **your own computation returned**
  for that matrix, form the `r × r` minors, and solve them **all** to zero **for `k_1 … k_D` over the
  REALS**; emit the solution set;
- the Q3 root-coincidence locus alongside it;
- ⭐ for each locus, **whether it intersects the region allowed by §3** (`Σ k_m² > 0`, all `k_m` real, all
  control parameters in their declared ranges), emitted as a symbolic test with its operands.

**Q8b — ⭐⭐ RE-RUN Q3 AND Q4 ON EACH ALLOWED STRATUM.** For **every** stratum found in Q8a that
**intersects the allowed region**, choose an explicit point on it satisfying §3, and recompute and emit the
**complete Q3 spectrum and Q4 `N1`–`N7` set** at that point, tagged as belonging to that stratum.

⛔⛔ **THIS IS NOT BOOKKEEPING.** A generic rank does not determine the rank on an exceptional stratum.
Reporting a locus without recomputing on it leaves the spectrum and mode count there untested.

⭐ If a package has **no** allowed stratum, ⭐ **emit a tag saying the list is empty** — ⛔ do not omit the
tag (corollary 4). ⭐ The tag is what distinguishes *checked and empty* from *never checked*.

⭐ The orchestrator checks the enumeration is complete; ⛔ the engine does not assert that it is.

---

## 7 · Packages — the sweep and the controls

⛔⛔ **EVERY control re-enters at the ACTION (§2), never at a result.** A control is a **different `L`**;
everything downstream is recomputed from it by the identical code path.

⭐ **Every package emits the FULL Q1–Q8 tag set.** ⛔ Emission is never conditional on a value
(corollary 4). ⚠ A quantity that is **identical** across packages **must still be emitted in each** — that
repetition is a result, and deleting it destroys the evidence.

There are six packages: `MAIN` is the baseline, and the other five are controls. Their complete actions
are:

```
L_MAIN           = (ρ_br / 2) · Σ_j (∂_t u_j)² − (μ_R / 2) · S_curl[∂u]
L_XFORM_FULLGRAD = (ρ_br / 2) · Σ_j (∂_t u_j)² − (μ_R / 2) · Σ_{i,j} (∂_i u_j)²
L_XFORM_DIVONLY  = (ρ_br / 2) · Σ_j (∂_t u_j)² − (μ_R / 2) · (Σ_i ∂_i u_i)²
L_XFORM_SIGNFLIP = (ρ_br / 2) · Σ_j (∂_t u_j)² + (μ_R / 2) · S_curl[∂u]
L_XFORM_ANISO    = (ρ_br / 2) · [s_ρ (∂_t u_1)² + Σ_{j=2..D} (∂_t u_j)²] − (μ_R / 2) · S_curl[∂u]
L_XCOEF_SCALE    = (ρ_br / 2) · Σ_j (∂_t u_j)² − (s · μ_R / 2) · S_curl[∂u]
```

| package | `D` |
|---|---|
| `MAIN` | 2, 3, 4, 5 |
| `XFORM_FULLGRAD` | 3, 4 |
| `XFORM_DIVONLY` | 3, 4 |
| `XFORM_SIGNFLIP` | 3, 4 |
| `XFORM_ANISO` | 3, 4 |
| `XCOEF_SCALE` | 3 |

⚠ **`s` is DIMENSIONLESS by declaration, and that is a Q6 input** — for `XCOEF_SCALE` the energy-density
requirement alone fixes only the *sum* of a scale factor's dimension and its coefficient's. ⛔ **That
reason does not hold for `s_ρ`**, so `s_ρ` gets **no** dimensionless declaration.

⭐ **The counting rule, and it must be stated because two engines read it oppositely otherwise:** a symbol
**declared** dimensionless is ⛔ **not** a `Q6` unknown; a symbol whose dimension the action **determines**
⭐ **is** one, and stays one. ⇒ `s` is excluded from the unknown count for `XCOEF_SCALE`; `s_ρ` is
**included** for `XFORM_ANISO`.

⭐ **And `[s_ρ]` is then a SOLVED value like any other coefficient dimension** — emit it as `Q6` emits the
rest. ⛔ **This file states no value for it.** ⚠ If you want a second operand to difference it against,
take it from the **registry**, ⛔ never from this file — `§Q6d` already says the registry is the one place
that comparison can genuinely fail.

Their positive and distinctness conditions are package-domain premises in the joint set in §3.

For `XCOEF_SCALE` specifically, `s ≠ 1` excludes the unit-scale collapse into `MAIN`. Without that
premise, an implementation could satisfy every stated premise while reproducing the baseline, leaving the
coefficient control unable to police arithmetic. A premise that keeps a control distinct is ⛔ not a
premise that forces a solver to decide.

⚠ **`XFORM_ANISO` breaks isotropy on ONE axis only.** ⛔ Do not "improve" it by making all the densities
distinct; that would define a different control.

⚠ **`MAIN` stops at `D = 5` for ENGINE SYMMETRY.** ⛔ A `D` present in one engine only is worse than absent
from both — it silently drops that point from every cross-engine comparison.

⭐ **`MAIN` swept over `D` is the primary result of the step, ⛔ not a control.**
⭐ **`XCOEF_SCALE` is a COEFFICIENT control**: it varies a coefficient under an unchanged stiffness form.
⛔ It cannot substitute for a FORM control, and ⛔ **no expectation is stated here about what it does or
does not move.** Run it and emit its full tag set like every other package.

⛔⛔ **THE RUN RECORD MUST BE OBSERVED, ⛔ NOT DECLARED.** ⭐ Accumulate a `(package, D)` pair **only after
that package has finished emitting**, and emit `RUN_PAIRS` / `SKIPPED_PAIRS` **after** the sweep, with
`SKIPPED = declared ∖ completed`. ⭐ It is the one tag §7 requires to say what actually happened.

⚠ **Runtime.** No script may exceed **10 minutes**. If `D = 6` with a control is too slow, ⭐ **narrow the
`D` list and emit a tag recording exactly which `(package, D)` pairs were run and which were skipped.**
⛔ Never silently drop one, and ⛔ never raise a timeout.

---

## 8 · Tag grammar — ⛔ both engines must produce PARALLEL tag sets

```
<ENGINE>_S10_<PACKAGE>_D<n>_<QUANTITY>
<ENGINE>_S10_<PACKAGE>_D<n>_ROOT<r>_<QUANTITY>
<ENGINE>_S10_<PACKAGE>_D<n>_STRATUM<s>_<QUANTITY>
<ENGINE>_S10_<PACKAGE>_D<n>_STRATUM<s>_ROOT<r>_<QUANTITY>
```

- `<ENGINE>` is `WL` or `PY`.
- `<PACKAGE>` is exactly one of the §7 names.
- `<n>` is the integer brane dimension of that run.
- `<s>` indexes the allowed Q8 strata in the engine's emitted stratum ordering, `1`-based.
- `<r>` indexes the **distinct** roots, ordered by the engine's own de-duplicated list, `1`-based.
- `<QUANTITY>` names **the object**, ⛔ never its value.

`STRATUM<s>` sits immediately after `D<n>` and before `ROOT<r>` when both scopes apply. Emit, once per
package-and-dimension, a tag listing the stratum ordering used so the orchestrator can align the indices.
The as-built engines' stratum tag names are not aligned, and their comparison is manual.

⭐ Emit one line per tag: `TAG: <payload>`, payload in a fully re-parseable form
(`InputForm` in Wolfram; `sympy.srepr`-safe string or `str()` of the expression in SymPy).
⭐ Also emit, once per package-and-dimension, a tag listing the **root ordering** used, so the two engines'
`ROOT<r>` indices can be aligned by the orchestrator rather than assumed to match.

### ⭐ Engine-local tags — declare them, so parity is meaningful

Some tags **cannot** exist in both engines: engine 2's registry comparison has no engine-1 counterpart, and
each CAS emits its own solver-condition tags. ⛔ Those are **not** disagreements, and a parity checker that
reports them as gaps trains its reader to ignore it.

⭐ Give every such tag the infix `_LOCAL_` immediately after the engine prefix — `WL_S10_LOCAL_…`,
`PY_S10_LOCAL_…` — and ⭐ emit one tag per engine listing every `_LOCAL_` name it produced.

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
