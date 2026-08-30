# Build directive (v2, folded from 2 decision-list legs) — S11c-b requested-truncation complete-row residual

## 0 · What this instrument is

The committed layer (`scripts/S11c_b_adjudicated_comparison.py`) and step-1 multigrade compare the two
engines' operands **per term-origin bucket**. A per-bucket difference does not decide whether the two
engines' **operators** agree, because the engines partition one operator into buckets by different
conventions. The object that decides agreement is the **complete operator per semantic row, compared at the
spec's requested truncation**. The reviewed layer already extracts and aligns those complete rows
(§2 below); this instrument applies the requested truncation to the aligned complete-row operands and
prints the per-row residual under the correct equivalence for each row's type. It performs **no** physical
classification and states **no** verdict (rule 2). Interpretation belongs to the step record.

⚠ This is v2. A v1 that instructed a hand-specified bucket sum was rejected by two decision-list legs for
sign errors, non-additive conflation, wrong coupling equivalence, an unavailable "full" residual, and leaks
(record: `directives/_measurements/S11c_b_row_residual_instrument_directive_review.md`). v2 reuses the
layer's own row extraction instead of reconstructing rows.

## 1 · Locked physics (supplied and verified — use exactly; do not re-derive)

### 1a · The requested truncation (§2a `S11c_b_SHARED_PHYSICS.md:202-207`, §3a `:248-254`)
Retain a term `ε^c η^a σ_W^b` **iff `c≤1 ∧ a≤1 ∧ b≤1`**. Grade by amplitude-factor power (bookkeeper
exponent), never by spatial-derivative order (a jet `∂^n W_bg`, any `n≥1`, carries exactly one `σ_W`).
Coefficients are Taylor-**linearized** in the background amplitude: e.g. `W_bg²→W_0²(1+2η w1)+O(η²)`,
`(1+η w1)^{-1}→1−η w1+O(η²)`. §3b's "retain full spatial dependence" means do **not** freeze a coefficient
at its constant binding before differentiation; it does **not** override this amplitude truncation. (This
reading was derived independently by both decision-list legs from §2a/§3a and is not a spec ambiguity.)

### 1b · ε handling and no untruncated residual (leg C8)
WL applies its background truncation **before emission** (`…_mathematica_audit.wl:154-160, 994-1003`), so WL
operands are already at `(a≤1 ∧ b≤1)`. Apply the requested truncation to the **PY** operands (WL is already
truncated). ⛔ Do **not** emit an "untruncated / full" cross-engine residual — it would compare a fuller PY
against a pre-truncated WL and is meaningless. The `ε` order is **not** recoverable by counting `ε` symbols
(the layer's extractor strips PY's `ε` via `coeff_epsilon`, `…comparator.py:455-466`; WL carries `ε` order as
metadata). Attach the `ε` order `c` from the object's **family metadata**: `c=1` for the wave operator rows
and the coupling kernel; `c=0` for the admissibility operand (a background-order object).

### 1c · The equivalence relation is per row type (legs G8, C6)
- **Strong slab rows** (§3b divergence-form EOM): the residual is the **exact** difference after activating
  any held in-plane divergence to strong form and applying the requested truncation. ⛔ Do **not** quotient
  a strong slab row modulo total in-plane divergence — §1d (`:159-168`) states that integrating a variable
  coefficient by parts generates first-background-jet terms that are **physics in the operator**, so
  quotienting them would manufacture false zeros. Integral linearity (splitting `∫`) is a permitted
  canonicalization; the total-divergence **quotient** is not.
- **The coupling kernel** (§3c weak variational restriction, IBP boundary term fixed to zero by compact
  support, `:292-308`): the residual is defined **modulo exact total in-plane divergence** (the full
  `∇·(cF)`, product rule included; `c∇·F` alone is not discardable per §1d). Use the layer's already-gated
  weak route (`_is_weak_scalar_density`→`classify_total_divergence`, `…adjudicated_comparison.py:779-782,
  726+`). §3c requires **both** cross-sector blocks — transverse→thickness **and** thickness→transverse —
  and their relabelled adjoint operands; compute both, do **not** manufacture the reverse block as
  "± forward".
- **The admissibility operand** (§3d, `:325-342`): an ordered pairing (bulk-DOF body force + per-face
  traction). Its components live in different slots and must **not** be summed into one row. Compare the
  full ordered association **componentwise** (each of the three momentum components, `θ`, `e_W`, and each
  per-face traction), for all four anchoring/density cases, with `c=0`.

### 1d · Faces are their own object, not folded into a strong slab row (leg C5/G3) — scope, do not fabricate an adapter
PY emits its faces as the raw S11c-a substrate bundle (`FACE_FLUX_BOUNDARY_OPERANDS`,
`…_sympy_audit.py:1434-1455,1543-1550`); WL folds faces into generalized rows before adding them to its
operator rows (`…_mathematica_audit.wl:807-814,839-847`). These are **not commensurate**, and the layer does
**not** perform the variational conversion between them. Therefore compare the strong slab rows on the
**bulk + kinetic + mass** content that both engines emit commensurately (PY `EXPANDED` and
`ADVECTIVE_MASS_OPERAND`; WL `U_MOMENTUM_ROWS`/`THICKNESS_ROW`/`MASS_EVOLUTION_ROW`). Where an aligned WL
complete row carries a folded face contribution that PY's `EXPANDED` lacks, **emit the face-attributable
piece as its own named object** (so it is visible, not silently differenced), and keep it **out** of the
bulk/kinetic/mass row residual. ⛔ Do not build a PY-generalized-face adapter here (out of scope) and do not
scalar-add PY's raw substrate bundle to a strong row.

## 2 · Reuse the layer — do not rebuild alignment or row extraction (legs C2/C4)

Import `S11c_b_adjudicated_comparison` (call it `A`) exactly as step-1 does. Its extractor already maps
**both** engines to a common semantic row schema and already uses the correctly-signed complete forms — you
do **not** sum term-origin buckets (that path caused v1's sign flip and non-additive conflation):
- `A.C.extract_slab` (`…cross_engine_comparator.py:760-800`): PY `U_BODY_BALANCE→U_MOMENTUM_ROWS` via the
  `EXPANDED` form (its own signs), `THETA_BALANCE→MU_THETA`, `E_W_BALANCE→THICKNESS_ROW`,
  `ADVECTIVE_MASS_OPERAND→MASS_EVOLUTION_ROW`, faces to a **separate** `FACE_SHAPE_SUBSTRATE`; WL's emitted
  `ROWS` (`U_MOMENTUM_ROWS`/`THICKNESS_ROW`/`MASS_EVOLUTION_ROW`/`CENTER_FACE_GENERALIZED_ROW`) and
  `DIVERGENCE_FORM_SOURCE.MU_THETA`.
- `A.transform`, `A._bridge_d`, `A.WL_TO_PY_RENAME` (rename/bridge alignment); `A._kinetic_pairs`
  (`…adjudicated_comparison.py:785+`, the thickness/u inertia pair); `A._admissibility_py_parts`
  (`…cross_engine_comparator.py:929+`); `A.classify_total_divergence` and `A._is_weak_scalar_density` (the
  weak coupling route); the layer's integral-linearity-aware zero test.

The new capability is only: **apply the requested truncation (§1a) to the aligned PY operand, then form and
print the per-row residual under the §1c equivalence for that row's type.** Reuse everything else.

## 3 · Objects to emit (per anchoring branch × density representative)

For each of the aligned semantic rows the layer exposes:
- `ROW_OPERAND_PY_TRUNC` := the aligned PY complete-row operand under the §1a requested truncation (strong
  form; any held divergence activated).
- `ROW_OPERAND_WL` := the aligned WL complete-row operand (already truncated at emission).
- `ROW_RESIDUAL` := their difference under the §1c equivalence for that row type — **exact** for strong slab
  rows (bulk+kinetic+mass), **modulo exact total in-plane divergence** for the coupling kernel (both blocks
  + relabelled adjoints), **componentwise** for the admissibility operand.
- For coupling: also emit the `classify_total_divergence` certificate object (the recovered potential `V`
  or the certified non-divergence residual) at the requested truncation.
- For any WL row carrying a folded face piece PY's `EXPANDED` lacks (§1d): emit `ROW_FACE_ATTRIBUTED` as its
  own object and exclude it from `ROW_RESIDUAL`.

Emit each as a CAS object with its `(η,σ_W)` multigrade support (the `ε` order `c` attached from metadata,
§1b). Emit an **assembly accounting** object per engine per row: the list of layer-extracted operands that
entered the row, so completeness is auditable (no operand silently dropped or double-counted).

## 4 · Guards (compute → emit → then guard; never assert a physics conclusion) (legs C9/G6)

- Validate the requested-truncation function against an **independent** second implementation (a direct
  `sp.series`/Taylor projection written a different way) on a synthetic fixture, and confirm they agree; a
  reconstruction identity `operand − polynomial − remainder ≡ 0` from a **single** route is zero by
  construction and is **not** a validation — do not hard-stop on it as a physics check.
- Exercise the truncation with a **one-sided synthetic corruption** (a fixture carrying an `η²`/`σ_W²` term
  that must be dropped, and an `η¹`/`σ_W¹` term that must survive) and emit the before/after so a reviewer
  sees the bound bites.
- An **assembly-accounting** guard (an aligned operand present in the layer's extraction but absent from a
  row's accounting list) may hard-stop — it is an instrument-integrity check, not a physics claim.
- ⛔ Do **not** assert any `ROW_RESIDUAL == 0` or `!= 0`. Its value is the measurement.

## 5 · The three clauses, delegation, and definition-of-done

**1. PRINT computed objects; never state conclusions.** Every `emit` payload is a CAS object; never a prose
result, never a family label ("REPRESENTATIONAL"/"GENUINE"), never a which-engine judgement.
**2. PRINT the residual; do not assert it.** Only the arithmetic/assembly-accounting guards (§4) may
hard-stop. `ROW_RESIDUAL` is emitted, never asserted.
**3. Interpretation belongs to the step record.**

> You construct **nothing physical**: every operand is read from the committed engine transcripts through
> the reviewed layer `A`. Every derived expression must be reached by `A`'s extraction/alignment + the
> requested-truncation function + the §1c equivalence — never hand-typed. No tag name may encode a value,
> sign, ratio, or grade of any residual. No emission may be conditional on a residual's value (only on which
> row / case / engine / block it belongs to). The per-row extraction and the coupling both-blocks handling
> are **delegated** to you: discover them from the layer's payloads and the cited adapters; the two build
> legs will form-ablate your row assembly, your truncation, your divergence activation, and your weak route.

**Definition of done (accounting):** every semantic row the layer exposes has a printed
`ROW_OPERAND_PY_TRUNC`, `ROW_OPERAND_WL`, and `ROW_RESIDUAL`; the coupling kernel has **both** cross-sector
blocks + relabelled adjoints with a divergence certificate each; the admissibility operand is compared
componentwise (all body DOFs + all per-face tractions) with `c=0`; every aligned operand is accounted into
exactly one row (assembly-accounting guard passes); the requested-truncation function passes its independent
cross-check and its one-sided corruption.

## 6 · Withheld (rule 5)
This directive states **no** expected residual value, **no** per-family verdict (representational vs
genuine), and **no** claim about which engine is spec-correct for any row. The instrument's job ends at
compute-and-print; the diff is adjudicated off-instrument.

## 7 · Deliverable and builder report
- Deliverable: `research/pde_ledger_v3/scripts/S11c_b_row_residual.py`, runnable as
  `python3 research/pde_ledger_v3/scripts/S11c_b_row_residual.py`, deterministic, exit 0 on the §4 guards,
  emitting in the step-1 grammar (a companion reporter is not required).
- No network, no committed-file writes, no engine re-run — read through the layer `A` only.
- Report (to `directives/_measurements/S11c_b_row_residual_instrument_build_directive.md`), each claim
  carrying its command: deliverable path + line count; run command + stdout size; the named
  requested-truncation function, its independent cross-check result, and its one-sided-corruption
  before/after (literal stdout); the named §1c equivalence dispatch (file:line for the strong vs weak vs
  componentwise branches); the assembly-accounting result for the mass row and the thickness/kinetic row
  (literal stdout); confirmation that both coupling blocks are emitted; and confirmation that no `assert`
  precedes any `ROW_RESIDUAL` emit (file:line of each residual emit).
