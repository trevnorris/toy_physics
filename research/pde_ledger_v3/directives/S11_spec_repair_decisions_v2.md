# S11 spec repair — decision list

**Target:** `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`, **914 lines, base `fc920079` = HEAD.**
⛔ **NOT `f49a1684`** (877-line predecessor; repairing it resurrects deleted defects).

⛔ **Scope: that file only.** `scripts/S10_exports.py` is read as an input under item 3, ⛔ not edited.
⛔ **No engine is touched by this round.**

⚠ Each item states **what must become true of the spec**. ⛔ No item states, implies, or supplies a reason
from which a builder could derive what any package, root, rank, locus or dimension comes out to be.

---

## 1 · A package is a kinetic density and a stiffness functional

**1a.** `§2:41` writes one kinetic density outside the package system; `§7:737-746` defines a package as a
stiffness functional alone. `§2` must define a package as an ordered **pair**, and name both members.

**1b.** The seven existing packages keep the current kinetic density, named `T_ISO`:

```
T_ISO  =  (ρ_br/2)·Σ_{j=1..D} (∂_t u_j)²
```

**1c.** `§7` adds an eighth package, `XKIN_ANISO`:

```
kinetic     T_ANISO  =  (ρ_br/2)·Σ_{j=2..D} (∂_t u_j)²   +   (s_ρ·ρ_br/2)·(∂_t u_1)²
stiffness   W_MAIN
D sweep     2, 3, 4, 5
```

- The distinguished axis is **`j = 1`**. ⛔ No relation between it and `k` is supplied; `k` stays symbolic.
- `T_ANISO` has exactly **two declared additive terms** — the `j = 2..D` sum is one declared term, the
  `j = 1` term is the other. ⛔ Neither is expanded before the `Q6:406-412` inventory is built.

**1d.** `§3:113-116` gains the package-domain premises `s_ρ > 0`, `s_ρ ≠ 1`.
⛔ `s_ρ` is **not** declared dimensionless; it is a `Q6` unknown. ⛔ No dimension is stated for it.

**1e.** Propagation — each site currently assumes one kinetic term with one coefficient:

| locus | what must become true |
|---|---|
| `Q6:408`, `Q6:416` | *"the inertial coefficient"*, **singular**, must admit **every** inertial coefficient the package's kinetic density carries |
| `Q1:241-242` | emits `STIFFNESS_TERMS` only. It must also emit **`KINETIC_TERMS`** — the package's declared kinetic terms, same list shape, one entry per declared term |
| `§7:801-806` | the form-vs-coefficient classification is made *"from the emitted stiffness functional"*. It must be made from the emitted **action** |
| `§9:878-885` | must gain a row placing the **isotropy of the kinetic form** in the supplied column |

⭐ Once `s_ρ` is in `COEFFICIENT_ORDERING`, `Q6`, `Q8a`'s `_COEFF_`/`_JOINT_` loci, `Q10` and `Q11`'s
`_KW_ZERO_LOCUS_` reach it under rules already written. ⛔ Add no per-question clause for it.

---

## 2 · `Q8` gains the matrix that governs the transverse count

`Q4:321` defines `nu_T` from `[M_r; kᵀ]`. `Q8a:534` builds rank-drop loci from the minors of `M_r` alone.
See `DEFECT_REGISTER.md#C16` (`:517`), which allows two fixes. ⭐ **Take the stacked source.**

**2a.** `Q8a` gains a second minor family: with `σ` the **stacked rank the engine's own `Q4`
`ROOT<r>_N3_STACKED_RANK` returned**, all `σ × σ` minors of `[M_r; kᵀ]`, emitted as
`ROOT<r>_STACKED_DROP_MINORS`, with the **same three solve-variable sets** and the **same five-object
locus protocol** (`§5:204-227`). `Q8a:544-547`'s rank-`0` convention applies unchanged.

**2b.** `Q8b`'s candidate components come from **three** sources — `M_r` rank-drop, stacked degeneracy,
root coincidence. Each emitted stratum emits **`STRATUM<s>_SOURCES`**, the sources that produced it.
⛔ This file requires no deduplication rule and ⛔ forbids none; the spec says nothing either way.

**2c.** `§9` gains a row: **the completeness of the stratum enumeration is not established by this build.**

---

## 3 · `Q6r` stops requiring a deleted directory

`Q6r:476-477` names `research/pde_ledger_v3/reduction/`, deleted (63 files); `S9_REWRITE_PLAN.md:6` states
there is no registry. Running the engine gives `ModuleNotFoundError: No module named 'registry_read'`.

⭐ Repoint `Q6r` at the export chain. The engine that imports the previous step's `LEDGER` emits, per
coefficient:

**3a.** The **name map is explicit and closed**, because the spec's symbols and the `LEDGER`'s keys are
spelled differently:

```
ρ_br → rho_br        μ_R → mu_R        B_comp → B_comp       μ_br → mu_br
β    → beta          s   → s           s_ρ  → s_rho          c_s0 → c_s0
```

**3b.** Resolution path: `LEDGER[name]['dimension_key']`, then `LEDGER[that key]['value']`.
⭐ A coefficient with no entry, or an entry with no `dimension_key`, yields **no `Q6r` row** — stated
explicitly, ⛔ not left implied.

**3c.** Emitted per resolved coefficient: the **derived** vector, the **imported** vector, their
**difference**, and provenance from **both rows, named apart** — `class` and `step` of the coefficient row,
and `class`, `step` and `corroborated_steps` of the resolved dimension row.

**3d.** `Q6r:479-481` requires the declared row's `source_locus`; ⛔ no `LEDGER` entry carries that field.
Its role (`:483-487`) was to separate an independent operand from one reproducing its predecessor. Under
the chain the imported vector is the predecessor's **by construction**. ⇒ `Q6r` emits
`Q6R_RESIDUAL_SCOPE` with the single token `CROSS_STEP_CONSISTENCY_ONLY`, the construction already used at
`:268-271`, and `:483-487` is rewritten to that scope. ⭐ `Q6r` stays engine-local (`§8:865-874`); `§8:867`'s
justification, which names the registry, is restated for the chain.

---

## 4 · `Q3` names the multiplicity object, and the count is pinned

`Q3:284` asks for a solution set *"as returned, retaining multiplicity"*; a solution set does not retain
multiplicity. `Q3:286` then counts *"its own list"*, which is ambiguous on a map. Emit, with these names:

| tag | the object |
|---|---|
| `ROOT_MULTIPLICITY_PAIRS` | an **ordered list of `(root, multiplicity)` pairs**, multiplicity a positive integer, from `DET_M` as a univariate polynomial in `omegaSquared` |
| `ROOT_SOLUTION_SET` | the solution set as the solver returns it, carrying ⛔ **no** multiplicity claim |
| `ROOT_COUNT_ALL` | the **sum of the second fields** of `ROOT_MULTIPLICITY_PAIRS` |
| `DET_M_DEGREE` | the degree of `DET_M` in `omegaSquared`, computed from the polynomial, ⛔ not from the root list |
| `ROOT_DEGREE_RESIDUAL` | `DET_M_DEGREE − ROOT_COUNT_ALL` |

⭐ `ROOT_DISTINCT` is the de-duplication of `ROOT_MULTIPLICITY_PAIRS`' first fields under the `§3`
assumptions; `ROOT_COUNT_DISTINCT` counts it; `ROOT_ORDERING` is unchanged.
⚠ `ROOT_DEGREE_RESIDUAL` is emitted unconditionally, whatever it comes to (corollary 4).

---

## 5 · The live-read exemption is enumerated, closed, and field-scoped

`Corollary 5:187-191` requires a declaration to be read from the live object it declares. `:899-901`
exempts `PREMISE_INVENTORY`. The exemption must become a **closed list of exactly three entries**:

| exempt | granularity |
|---|---|
| `PREMISE_INVENTORY` | whole tag |
| `M_ROUTE_RESIDUAL_SCOPE` (`:270`) | whole tag — payload is one fixed token |
| `STRATUM<s>_ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED` `.FAILURE_TOKEN` (`:665`) | ⛔⛔ **that field only** |

⛔⛔ The other four fields of the pinned failure object at `:661-669` — `RECOMPUTED_ROOTS`,
`COEFFICIENT_ORDERING`, `DEFINING_EQUATIONS`, `EVALUATION_POINT` — **remain live-read obligations.**
⛔ No holder object may be manufactured to make any exempt token appear wired. ⛔ No general exemption for
"metadata" is created; `Q6R_RESIDUAL_SCOPE` from item 3 joins this list as a fourth whole-tag entry.

---

## What this round does not do

- ⛔ Does not touch the seven stiffness functionals, `§3`'s ansatz, `§4`, `§5`'s clauses or locus protocol,
  `Q2`, `Q5`, `Q7`, `Q9`, `Q10`, `§10`.
- ⛔ Adds no acceptance criterion, no verdict, and no statement of what any package produces.

## Registered separately — ⛔ not repaired in this round

⚠ Both are **pre-existing**; neither is created by the items above, and folding either into this round
would grow it past the defect it names.

1. **`STRATUM<s>` `Q3`/`Q4` reruns are point-dependent.** `Q8b:556` asks for *an* explicit point on the
   component; `:569` concedes a single point does not characterise a positive-dimensional one. Two allowed
   points on one component can yield different spectra, so the reruns are not a property of the component
   and are not comparable across engines by defining equations alone.
2. **`§5:213-221`'s locus protocol pins no branch normal form, admissibility algorithm, undecided policy or
   component ordering** — `_SOLUTION` is *"exactly as your CAS returns it"* (`:218`) and `_REAL_ADMISSIBLE`
   is per returned branch (`:221`) ⇒ one engine may omit an allowed component the other retains. Already visible in
   the committed outputs at `XFORM_EXTRA, D = 2` (`mathematica/out/…:561` against `scripts/out/…:524`).

⭐ Both need a `DEFECT_REGISTER` entry and their own decision. ⚠ Item 1's new package widens exposure to
the first; that is a reason to register it now, ⛔ not to repair it here.
