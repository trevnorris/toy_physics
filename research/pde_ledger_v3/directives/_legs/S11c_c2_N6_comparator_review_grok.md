I'll start by loading the review-legs skill and the governing N6 comparator artifacts so I can derive independently before trusting the script.I'll inventory the spec, engine tags, and reuse-base comparators first, then write independent ablations before opening the working comparator.Next I’ll pull the spec clauses and the actual engine tag inventories so the join contract is grounded in what the engines emit.SymPy streams are JSONL, not `TAG =` lines — I’ll inventory those objects and sample the inner-axis vocabularies next.I have the tag map. Next I’ll sample both engines’ inner-axis vocabularies and the column constructors so the schema bridge can be judged independently.I have the two vocabularies. Next I’ll sample WL inner keys and the inherited `residual()` machinery, then write independent scripts before opening the comparator.Independent view is in place: mechanical matches are grade/face/spelling only; block and kernel names must stay unmatched. I’ll write those derivation scripts, then open the comparator.I’ll write the independent pairing, `_NODES` reconstruction, and join-set scripts now, then open the comparator against those results.Independent join-set and schema results are filed. Next I’ll read the comparator against those mechanical pairing rules.Support values are real booleans and the comparator has no asserts or verdict tokens. Next I’ll run FORM, one-sided, PIT-seal, and schema-bridge ablations on the `/tmp` copy.The R_N6 keys still didn’t join because WL keeps a `COMPONENT` axis. I’ll collapse that (and the source `WAVE` axis) in the copy, and force PIT subtraction on a family that actually matches.# SOUND

The N6 comparator is a measurement instrument: it joins by mechanically identical typed keys, prints `operand_A` / `operand_B` / three-valued `A − B`, and decides nothing. I derived the join contract from the two engines’ committed tag streams and constructors **before** trusting the script, then ablated a `/tmp` copy. No path I probed manufactures agreement or disagreement.

Working tree was not modified. Evidence lives under `/tmp/n6_cmp_review/`.

---

## Independent view (before the script)

From the committed streams:

| Side | Shared object-name stems | Notes |
|---|---|---|
| Both | 36 families | Exact directive join set, plus `PROVENANCE` (must be surfaced, not joined) |
| SymPy-only | 79 stems | `_NODES`, `ARITHMETIC_DAG`, `REP_INVARIANCE_*`, `CONTROL_INDEPENDENCE_*`, diagnostic internals, `PRIMES` |
| WL-only | 4 | `N6_LOCAL_*` |

Prime sets **differ** (third prime): SymPy `{1000000009, 998244353, 1000000033}` vs WL `{1000000009, 998244353, 1004535809}`. A GF(p) subtraction across engines is meaningless.

Inner-axis vocabularies that are **not** mechanically identical (must stay unmatched):

- Weak/R_N6: PY `(THICKNESS_TO_TRANSVERSE, THETA)` + signature `{0,6,9,12}` vs WL `THETA_FROM_TRANSVERSE` + `LOCAL_BARE`/`FOURIER_*`, plus WL `COMPONENT`
- Reconcile `SOURCE_*`: PY `(face, grade)` (waves summed, 8 columns) vs WL `(face, wave, grade)`
- Jet vs field: PY `e_W_d1` vs WL `eW`

Mechanical identities that **should** join: grade `{1,η,σ}↔(η,σ)`, face `SUM↔0`, `delta_p_plus↔deltaPPlus`, `e_W↔eW`, axis order, `_NODES` `(srepr, "formal_jet")` round-trip, strong-row 0- vs 1-origin (`U0↔U1` as array index, **not** the string `U1`).

What this instrument **can** establish: a computed symbolic residual on keys that survive that mechanical map (carrier `C_E`/`C_M` is the live one); unmatched keys as `UndecidedResidual` + pairing-table rows; one-sided support as `NONZERO_WITNESSED` / `NO_NONZERO_FOUND`. What it **cannot**: PIT agreement, block/kernel physics identification, wave-axis collapse of the reconcile source.

Stdout: `/tmp/n6_cmp_review/00_independent_join_and_schema.stdout`

---

## Scrutiny (each with ablation stdout)

### 1. PIT seal — holds

`sealed_field` (`:623–628`) marks `primes`, `probe_numerators`/`probe_denominators`, `numerator_denominator`, sample fields. SymPy num/den is `SURFACED` with `pit_sealed`; WL SparseArrays are surfaced the same way. Support is typed from CSR **shape**, then `support_residual` (`:536–540`) — residues never enter `residual()`.

Forced PIT subtraction in a `/tmp` copy **does** emit `PIT_A_MINUS_B`; the working copy has **0** such lines.

```
CE_ORIGINAL_JOIN 1 PIT_LINES 0 HAS_FORCED False
CE_PITCOPY_JOIN 1 PIT_LINES 1 HAS_FORCED True
```

`/tmp/n6_cmp_review/04_run_patched_ablations.stdout`, `/tmp/n6_cmp_review/01_ablate_working_copy.stdout` §11.

### 2. Schema bridge — holds (this was the load-bearing call)

R_N6 PY key uses `BLOCK_DIRECTION`/`TEST_FIELD`/`KERNEL_SIGNATURE` (`:164–169`); WL uses `BLOCK_NAME`/`KERNEL_FAMILY` plus `COMPONENT` (`:196–200, :673`). Synthetic `THICKNESS_TO_TRANSVERSE,THETA` + `6` vs `THETA_FROM_TRANSVERSE` + `FOURIER_KOUT_Y`:

```
RN6_JOIN 0  UNMATCHED 2  AXIS_SET_MISMATCH 2
PAIRINGS sympy_only_column / wl_only_slot
A_minus_B UndecidedResidual('unmatched typed key; no mechanical sibling')
```

Not a count, not auto-joined.

FORM ablation of the copy (block/kernel fold **and** dropping `COMPONENT`) **manufactures** `join=1` / `matched` / `Integer(0)`. The working instrument does not have that path.

Same for reconcile `SOURCE_*`: dropping WL `WAVE` in the copy manufactures a join of summed-source vs per-wave-source; the original stays `join=0, unmatched=2`.

`e_W_d1` vs `eW` does not join. Face `0↔SUM` and grade `{1,0,1}↔(0,1)` do.

Jet encodings `thetaJet11↔theta_d1d1` and U-row 0- vs 1-origin are emitter-schema spelling (constructors: PY `U0,U1,U2` zip of `U[0:3]`; WL `U<>ToString[r]` for `r=1,2,3`). Independently: `U0_PY==U1_WL`, `U0_PY≠U2_WL`, `U1_PY≠U1_WL`.

### 3. Carrier-rep residual — computed, not sealed

`N6RC_CARRIER_EULERIAN` is in `CARRIERS`, not in `sealed_field`. Matched `C_E` prints operands and `A_minus_B`:

```
JOIN 1  A_MINUS_B Integer(0)  SEAL_SKIPPED_RESIDUAL False
disagree → Add(Mul(Integer(2), Symbol('x')), Mul(Integer(-1), Symbol('y')))
```

### 4. One-sided support — typed, not subtracted

Committed `nonzero_modular_numerator` values are JSON booleans. `NO_NONZERO_FOUND` vs `NONZERO_WITNESSED` → `UndecidedResidual`, not “structurally zero”. Native booleans are not differenced.

### 5. `_NODES` reconstruction — round-trips; missing DAG is `parse_failed`

`DAG.get` follows `ARITHMETIC_DAG.nodes[ref]` (`:489–525`). `formal_jet` uses `parse_sympy_payload` (`:479–480`). On the committed covariance DAG (6074 nodes):

```
FJ_SREPR Function('s11cc2FieldeW')(Symbol('s11cc2Y1', real=True), ...)
COMPARATOR_RECONSTRUCTED s11cc2FieldeW(...)
ROUNDTRIP_EQUAL True
OPAQUE_NEQ True
MISSING_DAG missing ARITHMETIC_DAG
```

`/tmp/n6_cmp_review/05_real_dag_formal_jet.stdout`

### 6. Join-set completeness + no leak

`SHARED` (`:49–58`) is exactly the directive join set, including `R_COV_CONTROL_DELTA`, `SOURCE_CONTROL_DELTA`, `PHI_DOMAIN_CENSUS`, `ACTUAL_CONTROL_PARAMETERS`. `_NODES` pair to WL `ARITHMETIC`, not to a WL tag. SymPy-only → `ACCOUNTING_ROW` / `no WL sibling` (`:960–965`, `:1004–1005`). `PROVENANCE` surfaced, digests sealed. `residual_target=none` (`:982`). No prime literals in the comparator.

`ACTUAL_CONTROL_PARAMETERS` field names differ (`kappa_a` vs `ADVECTION`). The script **does not** alias them (`:683–684`). Unmatched `FIELD_PATH` rows are the measurement; aliasing would be the defect.

### 7. DoD teeth

Zero extract from two nonempty shared containers increments `zero_extract_failures` (`:768–772`) and exits 2 (`:1042`). `join>0` is only a `RUN_ACCOUNTING` count (`:1036`), not the exit. All-`deferred_oversize` exits 2. Unit tests: **23 passed** (`/tmp/n6_cmp_review/06_unit_tests.stdout`).

No `assert`, no `PASS`/`AGREE`/`COVARIANT`/`VERDICT` tokens in the instrument (`:13` of ablation stdout).

---

## FORM / one-sided (mandatory)

On the unpatched `/tmp` copy, two-slot `C_E` (THETA, E_W):

- Baseline both `x` → both residuals `Integer(0)`
- Flip THETA sign **and** replace E_W with `y*z` → **both** residuals move (`Mul(Integer(2),x)` and `x-y*z`)
- Corrupt **only** SymPy THETA `x→y` → THETA residual moves, E_W stays `Integer(0)` (`ONESIDED_ONLY_THETA True`)
- Repoint WL payload under the same name → residual moves

`/tmp/n6_cmp_review/01_ablate_working_copy.stdout` §§2–4.

Axis-collapse copy vs original: `/tmp/n6_cmp_review/04_run_patched_ablations.stdout`.

---

## Ranked findings

None that survive the physics filter. Unmatched block/kernel/wave/component keys are coverage rows (`UndecidedResidual` + pairing table), not silent zeros and not counts.

**SOUND.**
