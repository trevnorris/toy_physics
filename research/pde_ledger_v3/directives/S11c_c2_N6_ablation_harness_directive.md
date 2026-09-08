# S11c-c2 N6 ablation harnesses — fixed four-engine knife list

## Role and authority

Build four committed ablation harnesses, one for each live N6 engine named below. This file fixes the complete
knife list, the construction sites, the FORM operations, the coefficient companions, and the objects to print.
The builder implements this list exactly: do not add, drop, combine, or retarget a knife.

Use the S11c-b carrier harness shape:
`scripts/S11c_b_carrier_ablation_harness_sympy.py` and
`directives/S11c_b_carrier_ablation_harness_directive.md`.

## Global harness contract

1. **Compact computed payloads only; never symbolic objects.** For the CANONICAL baseline, unablated copy, and every
   listed FORM, NO-OP/IDENTITY, and COEFFICIENT ×2 run, print tagged `{baseline, corrupted, diff}` records for every
   certified or DEAD object containing only the object's compact PIT fingerprint and SHA-256 digest. The SymPy PIT
   fingerprint is the engine-emitted `n.pit` projection keyed by object and
   `probe`/`route`: retain the emitted column keys and `nonzero_modular_numerator` bitmap/tally, but do not print the
   `numerator_denominator` sample matrix or an arithmetic DAG. The WL PIT fingerprint is a compact association keyed
   as `{case, "PROBE", cellIndex, prime}` whose value records the engine-emitted valid/rejected counts, the selected
   object's `nonzero_count`, and `circuit_leaves`; derive it only from that object's emitted `PROBE_NUMERATORS` and the
   engine's `LOCAL:PROBE` sample index. Hash the exact selected object payload before transcript projection with
   `n.sha(obj)` for SymPy (the function at `scripts/S11c_c2_N6_diagnostic_sympy.py:91`) and use the corresponding
   `emitted_object_sha256` entry from WL `RUN_PROVENANCE`. For non-PIT metadata, print the digest and any already-small
   scalar fields only; if the engine supplies no digest, compute only `n.sha(obj)` in the worker. There is no symbolic
   fallback. The `diff` member is a deterministic key-aligned arithmetic/bitmap delta or digest equality bit, not a
   symbolic difference. Small scalar residuals (`R_N6`, `R_cov`, `SPLIT_CHECK`, slot/closure guard residuals) and
   already-compact arithmetic key-aligned tables remain direct payloads. Subtract arithmetic tables on shared keys
   and record a tagged `MISSING` entry on either side for a key-set change; never coerce absence to numeric zero.
   The knife's mechanically visible payload is a moved PIT fingerprint and/or changed digest between the identical
   CANONICAL/NO-OP records and FORM/COEFFICIENT records. Payloads contain no interpretation, acceptance label, or
   physics conclusion.
2. **Wrap the live engine.** Obtain every object from the engine's production entrypoint and emitted tags. Do not
   reconstruct an engine object in the harness, and do not patch an extractor, emitter, serializer, or diff routine.
3. **One mutation, one source site.** Locate the named function or top-level assignment with `ast`/held WL source,
   count the exact old fragment as one inside that scope, replace it, and re-parse before execution. Line numbers are
   current-tree anchors; the literal fragment and named scope are authoritative.
4. **Temporary copies and compact spools only.** Put every changed engine tree under a fresh `/tmp` directory. Do not
   edit a production engine. Print a unified source diff and SHA-256 for every temporary variant. A worker may spool
   the engine's raw output only inside its fresh temporary directory long enough to compute the compact fingerprints
   and digests; never copy, echo, or embed a full symbolic object, full symbolic difference, arithmetic DAG, or raw
   engine transcript in a harness transcript, including on an error path. Discard the temporary spool with the tree.
5. **Process isolation and WL budget.** Run the canonical baseline, unablated copy, FORM variants, identity variants,
   and coefficient companions one at a time. Each SymPy run uses a fresh worker subprocess pinned as stated below.
   Every WL label, including CANONICAL and the unablated copy, means the production source plus the same declared
   build-scoped, NON-KNIFE budget patch in a fresh `/tmp` copy: retain only `{"LAB_HELD","RHOBR_CONSTANT"}` and use
   four valid PIT draws per prime/cell. Apply the exact two edits specified under Harness 1 before any knife patch;
   print the budget-patch diff and both the production and executed-source SHA-256 values separately from the knife
   diff. Use one fresh kernel under `timeout --kill-after=5 600` for each label and never overlap WL kernels. Never
   launch the unrestricted four-case engine from this harness.
6. **Imported siblings.** A SymPy worker's temporary `scripts/` directory goes first on `sys.path`. Copy the target
   engine and its locally imported N6/S11c-a/S11c-b siblings into that directory, patch only the named file, and let
   ordinary imports resolve there. Fresh workers prevent module caches and `_FACE_CACHE` from crossing variants.
7. **Canonical-copy drift print.** Run an unablated temporary copy through the same worker path and print its compact
   PIT-fingerprint/digest records, engine-source digests, emitted-object digests, and source-path provenance next to
   the canonical run. For WL, CANONICAL and the unablated copy both carry the identical NON-KNIFE budget patch and
   differ only in temporary path; compare their executed-source hashes as well as the untouched production hash.
8. **Guards follow compact payloads.** Print compact operands and residuals before schema and integrity guards. Guards
   may cover exact knife and budget-patch counts, parseability, subprocess exit, tag presence, key alignment,
   serialization, transcript byte count, and production-source immutability; they do not classify a physics payload.
   On failure, print only bounded diagnostics (exit code, stderr tail/size, spool SHA-256, and missing tags), never the
   raw object-bearing stdout or a symbolic payload.
9. **Self-tests for every retained FORM knife.** Run and print:

   - the exact FORM variant specified below;
   - a **NO-OP / IDENTITY** variant replacing the exact old fragment by itself;
   - the exact same-site **COEFFICIENT ×2** companion specified below; and
   - the knife's named one-sided **DEAD print set** from the same baseline and corrupted runs.

For every DEAD line below, print these named objects as the compact fingerprint/digest triples defined in clause 1
and add no statement about their values.

Every harness opens with a manifest recording the engine, pinned case, exact scope, old fragment, replacement,
classification, printed-object list, source hashes, Python/SymPy or Wolfram versions, and invocation.

## Deliverables

- `mathematica/S11c_c2_N6_ablation_harness.wl`
- `scripts/S11c_c2_N6_ablation_harness_wl.py`
- `_measurements/S11c_c2_N6_ablation_harness_wl.md`
- `scripts/S11c_c2_N6_covariance_ablation_harness.py`
- `_measurements/S11c_c2_N6_covariance_ablation_harness.md`
- `scripts/S11c_c2_N6_diagnostic_ablation_harness.py`
- `_measurements/S11c_c2_N6_diagnostic_ablation_harness.md`
- `scripts/S11c_c2_N6_reconcile_ablation_harness.py`
- `_measurements/S11c_c2_N6_reconcile_ablation_harness.md`

Each transcript is a KB-scale compact record of the exact invocation, manifests, source patches and hashes, PIT
fingerprint/digest triples, small scalar residuals/tables, and guards. It never contains a full symbolic object,
symbolic difference, arithmetic DAG, PIT sample matrix, or raw engine transcript.

---

## Harness 1 — blind WL engine

Engine: `mathematica/S11c_c2_N6_mathematica_audit.wl`.

Run the complete production file for the retained case; do not marker-truncate it. Static inspection finds no
environment or command-line parameter for either case selection or PIT draws, so every WL run applies these two
declared build-scoped, NON-KNIFE edits to its `/tmp` engine copy before any knife edit:

1. In the unique top-level `Do` driver at lines 1011–1014 whose body calls `buildCase[case]`, replace its exact old
   iterator

   ```wl
   {case, Tuples[{{"LAB_HELD", "MATERIAL_ADVECTED"}, {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}]}
   ```

   by

   ```wl
   {case, {{"LAB_HELD", "RHOBR_CONSTANT"}}}
   ```

   This selects the single SymPy-matched case without changing `buildCase`. The separate SETUP metadata occurrence
   of the four-case `Tuples[...]` at lines 582–583 is not the driver and is not patched.
2. In `probeCase`, replace the unique adaptive draw-count assignment at line 760

   ```wl
   drawCount = If[0 < bound < 1, Max[8, Ceiling[Log[2^-80/Max[1, unionCount]]/Log[bound]]], 8];
   ```

   by

   ```wl
   drawCount = 4;
   ```

   The emitted bounds and actual draw count remain production-computed payloads; this edit changes only the number
   of valid PIT samples, not `numericObjects`, the prime list, branch cells, or any knife site.

The held-source matcher must find exactly one driver `Do` containing `buildCase[case]` and exactly one old draw-count
assignment inside `probeCase`, then re-parse the complete held source. The case edit merely changes which argument is
passed to the unchanged `buildCase`, and the draw edit occurs after `numericObjects` has been constructed and compiled;
therefore the retained case's per-case symbolic objects are unchanged. `caseOrdinal`-derived PIT seeds may differ from
the retained case's ordinal in the former four-case schedule; every compared variant uses the same restricted schedule
and must print its actual seeds. Print the two restriction diffs, their SHA-256 values, and the manifest classification
`BUDGET_RESTRICTION / NON-KNIFE` for every run. The Python driver serializes all kernel invocations, parses the printed
`WL_S11CC2_*` assignments into compact fingerprints/digests, and runs each kernel under
`timeout --kill-after=5 600`.

Primary print set:

- `WL_S11CC2_N6RC_R_N6`
- `WL_S11CC2_N6RC_SPLIT_CHECK`
- `WL_S11CC2_N6COV_R_COV`
- `WL_S11CC2_N6_SLOT_GUARD_RESIDUAL`
- `WL_S11CC2_N6_CLOSURE_GUARD_RESIDUAL`

Also extract the operand/channel objects named in the DEAD sets below. Keep the engine's case keys and component keys
in their emitted order. Use the engine's own PIT configuration and prime list.

### K_carrier — FORM

- **Scope/site:** unique top-level assignment at line 17.
- **Old:** `materialNormalKnife = 0;`
- **FORM replacement:** `materialNormalKnife = 1;`
- **FORM:** enable the off-block material-covector contamination constructed only by
  `invT[[1, 4]] += knife jet["WBg", {2}];` in `materialGeometry` at line 289. The source-velocity construction calls
  `materialGeometry[anchor, s, 0]` separately at line 861.
- **NO-OP / IDENTITY:** `materialNormalKnife = 0;` → `materialNormalKnife = 0;`.
- **COEFFICIENT ×2, same site:** `materialNormalKnife = 0;` → `materialNormalKnife = 2;`.
- **DEAD print set, exact:** `WL_S11CC2_N6RC_CARRIER_EULERIAN`,
  `WL_S11CC2_N6RC_SOURCE_EULERIAN`, `WL_S11CC2_N6RC_SOURCE_MATERIAL`,
  `WL_S11CC2_N6RC_SOURCE_BRIDGE_RESIDUAL`, `WL_S11CC2_N6COV_SOURCE_ACTUAL`,
  `WL_S11CC2_N6COV_SOURCE_PREDICTED`, `WL_S11CC2_N6COV_R_COV`, and
  `WL_S11CC2_N6COV_FROZEN_PHI`.

The K_carrier DEAD set does not include `WL_S11CC2_N6COV_R_COV_INCREMENT`, whose construction reads `cm` at
lines 896–900, or `WL_S11CC2_N6COV_ACTUAL_CONTROL_PARAMETERS`, whose payload embeds `materialNormalKnife` at
lines 977–980.

### K_junk — FORM

- **Scope/site:** unique top-level assignment at line 20.
- **Old:** `actualJunkCoefficient = 0;`
- **FORM replacement:** `actualJunkCoefficient = 1;`
- **FORM:** enable the separate dimensioned amplitude family `jk junkMu eW` at the unique `materialAmplitude`
  construction `el[pulled] + jk junkMu eW` at line 247. `actualJunkCase` remains unchanged.
- **NO-OP / IDENTITY:** `actualJunkCoefficient = 0;` → `actualJunkCoefficient = 0;`.
- **COEFFICIENT ×2, same site:** `actualJunkCoefficient = 0;` → `actualJunkCoefficient = 2;`.
- **DEAD print set, exact:** `WL_S11CC2_N6RC_CARRIER_EULERIAN`,
  `WL_S11CC2_N6RC_CARRIER_MATERIAL`, `WL_S11CC2_N6RC_CARRIER_BRIDGE_RESIDUAL`,
  `WL_S11CC2_N6COV_SOURCE_PREDICTED`, `WL_S11CC2_N6COV_FROZEN_PHI`, and
  `WL_S11CC2_N6COV_PHI_DOMAIN_CENSUS`.

### K_split_route — FORM

- **Scope/site:** `buildCase`, unique `sourceChannel` assignment at line 892.
- **Old:**

  ```wl
  sourceChannel = joinFaceMaps[buildContraction[cm[#], es[#] - ms[#], response[#], #, False] &];
  ```

- **FORM replacement:**

  ```wl
  sourceChannel = joinFaceMaps[buildContraction[cm[#], es[#] - ms[#], response[#], #, True] &];
  ```

- **FORM:** change only the final `affine` route selector for the source channel, admitting the local affine family
  into that channel. Do not edit `carrierChannel`, `crossChannel`, `splitSum`, or `splitCheck`.
- **NO-OP / IDENTITY:** replace the old assignment by itself.
- **COEFFICIENT ×2, same site:** retain final `False` and replace only the source operand
  `es[#] - ms[#]` by `2 (es[#] - ms[#])` in this assignment.
- **DEAD print set, exact:** `WL_S11CC2_N6RC_CARRIER_EULERIAN`,
  `WL_S11CC2_N6RC_CARRIER_MATERIAL`, `WL_S11CC2_N6RC_CARRIER_BRIDGE_RESIDUAL`,
  `WL_S11CC2_N6RC_SOURCE_EULERIAN`, `WL_S11CC2_N6RC_SOURCE_MATERIAL`, and
  `WL_S11CC2_N6COV_R_COV`.

There is no diagonal-sampler mutation in this harness.

---

## Harness 2 — SymPy covariance engine

Engine: `scripts/S11c_c2_N6_covariance_sympy.py`.

Pinned invocation for every subprocess: `--anchoring LAB_HELD --density RHOBR_CONSTANT --draws 8 --seed 110603`.
Drive `run(args)` only after reproducing `main()`'s `n.emit = r.emit = emit` binding; restore both bindings in
`finally` inside the worker.

Certified print set:

- `N6COV_R_COV`
- `N6COV_SOURCE_ACTUAL`
- `N6COV_SOURCE_PREDICTED`
- `N6COV_R_COV_INCREMENT`
- `N6COV_R_COV_CONTROL_DELTA`
- `N6COV_FROZEN_PHI`
- `N6COV_PHI_DOMAIN_CENSUS`

The JSONL outer names carry the engine's `S11CC2_` prefix; select by the complete emitted name rather than a suffix
substring.

### K_junk — FORM

- **Scope/site:** unique module assignment at line 35.
- **Old:** `ACTUAL_JUNK = sp.Integer(0)`
- **FORM replacement:** `ACTUAL_JUNK = sp.Integer(1)`
- **FORM:** enable the separate `kappa_j * junk * b.e_W` amplitude appended in `actual_amplitudes` at line 145.
- **NO-OP / IDENTITY:** replace the old assignment by itself.
- **COEFFICIENT ×2, same site:** `ACTUAL_JUNK = sp.Integer(0)` → `ACTUAL_JUNK = sp.Integer(2)`.
- **DEAD print set, exact:** `N6COV_SOURCE_PREDICTED`, `N6COV_FROZEN_PHI`, and
  `N6COV_PHI_DOMAIN_CENSUS`.

### K_circular — FORM

- **Scope/site:** `run`, one contiguous block at lines 186–188. Do not patch `predicted_amplitude`.
- **Old:**

  ```python
    mu_pred = predicted_amplitude(mu_e, phi)
    mu_actual, mu_baseline = actual_amplitudes(a, b, inputs, alpha, rho,
                                               ACTUAL_A_RHO, ACTUAL_JUNK)
  ```

- **FORM replacement:**

  ```python
    mu_actual, mu_baseline = actual_amplitudes(a, b, inputs, alpha, rho,
                                               ACTUAL_A_RHO, ACTUAL_JUNK)
    mu_pred = mu_actual
  ```

- **FORM:** identify the prediction amplitude operand with the actual amplitude operand inside `run`; retain the
  existing predicted-source call at line 191 and every downstream construction.
- **NO-OP / IDENTITY:** replace the old block by itself.
- **COEFFICIENT ×2, same site:** retain the old statement order and replace only
  `mu_pred = predicted_amplitude(mu_e, phi)` by
  `mu_pred = [2 * term for term in predicted_amplitude(mu_e, phi)]`.
- **DEAD print set, exact:** `N6COV_SOURCE_ACTUAL`, `N6COV_FROZEN_PHI`, and
  `N6COV_PHI_DOMAIN_CENSUS`.

### K_rank — FORM

- **Scope/site:** enclosing top-level `FunctionDef prolonged_phi`, inside nested `image_of`; the unique
  `images[atom]` assignment at lines 92–97.
- **Old:**

  ```python
            images[atom] = b.total_derivative(image_of(parent), direction,
                                             background_depth=3)
  ```

- **FORM replacement:**

  ```python
            images[atom] = (image_of(parent) if len(paths[atom][1]) >= 2 else
                            b.total_derivative(image_of(parent), direction,
                                               background_depth=3))
  ```

- **FORM:** preserve every discovered domain key; for a path of rank at least two, identify its image with its
  parent image (applied recursively up the path, this collapses every rank≥2 image onto its rank-1 ancestor). Do not
  change `paths`, `queue`, `candidates`, `background_depth`, `uncovered`, or the coverage check.
- **NO-OP / IDENTITY:** replace the old assignment by itself.
- **COEFFICIENT ×2, same site:** replace the old assignment by

  ```python
            images[atom] = ((2 if len(paths[atom][1]) >= 2 else 1) *
                            b.total_derivative(image_of(parent), direction,
                                               background_depth=3))
  ```

- **DEAD print set, exact:** `N6COV_SOURCE_ACTUAL`.

Do not offer `background_depth` or a BFS/path-domain truncation as a K_rank alternative.

---

## Harness 3 — SymPy diagnostic engine

Engine: `scripts/S11c_c2_N6_diagnostic_sympy.py`.

Pinned invocation for every subprocess: `--anchoring LAB_HELD --density RHOBR_CONSTANT --draws 8 --seed 110602`.
Invoke the production CLI or `run(Namespace(...))` in a fresh worker and parse its JSONL.

Certified print set:

- `S11CC2_REP_INVARIANCE_EULERIAN_OPERAND`
- `S11CC2_REP_INVARIANCE_MATERIAL_OPERAND`
- `S11CC2_REP_INVARIANCE_RESIDUAL`
- `S11CC2_N6_SLOT_GUARD_RESIDUAL`, for both route labels
- `S11CC2_N6_CLOSURE_GUARD_RESIDUAL`, for both route labels

Every diagnostic run also extracts the engine's own three named control objects for each probe:

- probe `TILT`: `S11CC2_CONTROL_INDEPENDENCE_BASE`,
  `S11CC2_CONTROL_INDEPENDENCE_CORRUPTED`, `S11CC2_CONTROL_INDEPENDENCE_RESIDUAL`;
- probe `N4_ADVECTION`: `S11CC2_CONTROL_INDEPENDENCE_BASE`,
  `S11CC2_CONTROL_INDEPENDENCE_CORRUPTED`, `S11CC2_CONTROL_INDEPENDENCE_RESIDUAL`.

Print each native control group in the order `{BASE, CORRUPTED, RESIDUAL}` and retain the `probe` label. Also include
each of those six named objects in the harness's compact fingerprint/digest variant triples.

### K_EW_rowdrop — FORM

- **Scope/site:** `face_factory`, its unique row construction at line 311; `run` binds the MATERIAL invocation to
  `m_rows` at line 806.
- **Old:**

  ```python
    rows=flatten({'U':folded['U'],'E_W':folded['E_W'],'THETA':mass-correction})
  ```

- **FORM replacement:**

  ```python
    rows=flatten({'U':folded['U'],'E_W':folded['E_W'],'THETA':mass-correction})
    rows.pop('E_W')
  ```

- **FORM:** remove the complete `E_W` row family from every `face_factory` row mapping after `flatten`. Do not alter
  `flatten`, `ROWS`, `unflatten`, the Eulerian imported rows, `residual`, or the fixed output-domain completion in
  `build_increment` lines 510–512. The production path is `face_factory.rows` → `m_rows` at line 806 → `m_coeff`
  at line 809 → `M` at line 851 → `REP_INVARIANCE_RESIDUAL` at line 853.
- **NO-OP / IDENTITY:** replace the old row-construction line by itself.
- **COEFFICIENT ×2, same site:** replace only `'E_W':folded['E_W']` in the old line by
  `'E_W':2*folded['E_W']`.
- **DEAD print set, exact:** `S11CC2_REP_INVARIANCE_EULERIAN_OPERAND`,
  `S11CC2_N6_MU_RECONSTRUCTION_IMPORTED`, `S11CC2_N6_MU_RECONSTRUCTION_NATIVE`, and
  `S11CC2_N6_MU_RECONSTRUCTION_RESIDUAL`.

`K_ewsign` is **COEFFICIENT**, not FORM. Retain it only as the additional coefficient-class companion at this same
row-construction site: replace `'E_W':folded['E_W']` by `'E_W':-folded['E_W']`, and print the diagnostic certified
set. It does not replace the required literal ×2 companion above.

### K_slotdrop — FORM

- **Scope/site:** `run`, unique `slots` construction at line 788.
- **Pinned slot:** `delta_p_plus` only.
- **Old literal fragment:**

  ```python
    slots=tuple(inputs.a(prefix+label) for label in ('plus','minus') for prefix in ('delta_p_','d_w_delta_p_'))
  ```

- **FORM replacement literal fragment:**

  ```python
    slots=tuple(inputs.a(name) for name in
                ('d_w_delta_p_plus','delta_p_minus','d_w_delta_p_minus'))
  ```

- **FORM:** remove the `delta_p_plus` carrier-slot family while retaining the other three native symbols and their
  order. Do not select another slot and do not alter `build_increment`, `slot_guard`, or `closure_guard`. The same
  pinned tuple feeds the coefficient extractions at lines 797 and 809 and both increment calls at lines 850–851.
- **NO-OP / IDENTITY:** replace the old literal fragment by itself.
- **COEFFICIENT ×2, same site:** replace the old fragment by this exact local coefficient wrapper:

  ```python
    slots=tuple(inputs.a(prefix+label) for label in ('plus','minus') for prefix in ('delta_p_','d_w_delta_p_'))
    native_pressure_coefficients=globals()['pressure_coefficients']
    def pressure_coefficients(rows,active_slots):
        table=native_pressure_coefficients(rows,active_slots)
        return {key:(2*value if key[1].name=='delta_p_plus' else value)
                for key,value in table.items()}
  ```

  This companion keeps all four slots and multiplies only the already extracted `delta_p_plus` coefficients by two.
- **DEAD print set, exact:** `S11CC2_N6_MU_RECONSTRUCTION_IMPORTED`,
  `S11CC2_N6_MU_RECONSTRUCTION_NATIVE`, `S11CC2_N6_MU_RECONSTRUCTION_RESIDUAL`,
  `S11CC2_N6_MU_AMPLITUDE` for routes `EULERIAN` and `MATERIAL`, and
  `S11CC2_N6_FACE_VELOCITY` for routes `EULERIAN` and `MATERIAL`.

There is no diagnostic advection knife. Do not edit `shift`, `source_tag`, `t`, or `material_pullback`. The
production `mu_m` construction specializes the tag at `t=1` on line 800; the `t=1`/`t=0` constructions remain the
engine-owned `N4_ADVECTION` control printed above.

Do not use `N6_DIMENSIONS.failure_operand` as a DEAD probe and do not add a dimension-checker interpretation payload.

---

## Harness 4 — SymPy reconcile engine

Engine: `scripts/S11c_c2_N6_reconcile_sympy.py`.

Pinned invocation for every subprocess: `--anchoring LAB_HELD --density RHOBR_CONSTANT --draws 8 --seed 110602`.
Drive `run(args)` only after reproducing `main()`'s `n.emit = emit` binding; restore it in `finally` inside the worker.

Certified print set:

- `N6RC_R_N6`
- `N6RC_SPLIT_CHECK`
- `N6RC_CARRIER_BRIDGE_RESIDUAL`
- `N6RC_SOURCE_BRIDGE_RESIDUAL`
- `N6RC_CARRIER_CHANNEL`
- `N6RC_SOURCE_CHANNEL`
- `N6RC_CROSS_CHANNEL`
- `N6RC_CARRIER_EULERIAN`
- `N6RC_CARRIER_MATERIAL`
- `N6RC_SOURCE_EULERIAN`
- `N6RC_SOURCE_MATERIAL`
- `N6RC_EULERIAN_OPERAND`
- `N6RC_MATERIAL_OPERAND`

### K_normal — FORM, scoped S11c-a construction wrapper

- **Scope/site:** `build_material_carrier`, unique `n.face_factory` assignment at line 91. The wrapped object is the
  MATERIAL `FaceSource.normal_exact` constructed by `build_material_face_source` at S11c-a line 773.
- **Old:**

  ```python
    rows, _, provenance = n.face_factory(a, b, inputs, alpha, rho, 'MATERIAL', mu_slot)
  ```

- **FORM replacement:** replace that one line by the following scoped wrapper:

  ```python
    original_builder = a.build_material_face_source
    original_cache = dict(a._FACE_CACHE)
    def mixed_material_face_source(*args, **kwargs):
        result = original_builder(*args, **kwargs)
        components = list(result.normal_exact)
        components[0] = components[0] + a.grad_W[1]
        length = sp.sqrt(a.dot(tuple(components), tuple(components)))
        return n.replace(result, normal_exact=tuple(value / length for value in components))
    a.build_material_face_source = mixed_material_face_source
    a._FACE_CACHE.clear()
    try:
        rows, _, provenance = n.face_factory(a, b, inputs, alpha, rho, 'MATERIAL', mu_slot)
    finally:
        a.build_material_face_source = original_builder
        a._FACE_CACHE.clear()
        a._FACE_CACHE.update(original_cache)
  ```

- **FORM:** mix the second material slope channel `a.grad_W[1]` into component zero of the MATERIAL exact normal and
  renormalize. The wrapper exists only around the carrier's `n.face_factory` call; both the builder and cache are
  restored in `finally`.
- **NO-OP / IDENTITY:** replace the old `n.face_factory` assignment by itself.
- **COEFFICIENT ×2, same site:** use the same scoped wrapper and restoration block, but replace only
  `components[0] = components[0] + a.grad_W[1]` by `components[0] = 2 * components[0]`.
- **DEAD print set, exact:** `N6RC_SOURCE_EULERIAN`, `N6RC_SOURCE_MATERIAL`,
  `N6RC_SOURCE_BRIDGE_RESIDUAL`, and `N6RC_EULERIAN_OPERAND`.

Do not patch S11c-a globally. Do not patch the Eulerian `normal_exact` at S11c-a line 850. Do not patch
`material_inverse_transpose` at S11c-a line 694. `face_velocity_raw` reads `source.normal_exact` at S11c-a lines
890–891, which is why the material-normal wrapper is confined to `build_material_carrier` after
`build_material_velocity` has supplied `m_v`.

### K_source_route — FORM

- **Scope/site:** `run`, unique `source` assignment at line 260.
- **Old:**

  ```python
    source = closed_response(comp, inputs, m_coeff, ds, kernels)
  ```

- **FORM replacement:**

  ```python
    source, _ = n.build_increment(comp, inputs, m_coeff, ds, kernels, slots)
  ```

- **FORM:** replace the closed-response-only source construction by the full increment route and unpack the returned
  `(output, dimensions)` pair. Do not index the pair and do not edit line 259 or `closed_response`.
- **NO-OP / IDENTITY:** replace the old assignment by itself.
- **COEFFICIENT ×2, same site:** retain `closed_response` and replace only its `ds` argument by
  `{s: {w: 2 * value for w, value in terms.items()} for s, terms in ds.items()}`.
- **DEAD print set, exact:** `N6RC_CARRIER_EULERIAN`, `N6RC_CARRIER_MATERIAL`, and
  `N6RC_CARRIER_BRIDGE_RESIDUAL`.

### K_operand_swap — FORM operand collapse

- **Scope/site:** `run`, the unique `ds` comprehension at lines 211–212 only.
- **Old operand fragment:** `ms[s].get(w, sp.S.Zero)`
- **FORM replacement:** `es[s].get(w, sp.S.Zero)`
- **FORM:** replace only the material-source read in `es - ms` by the corresponding Eulerian-source read, producing
  an `es - es` operand identification. Do not patch the source-channel call at line 260.
- **NO-OP / IDENTITY:** replace the old operand fragment by itself.
- **COEFFICIENT ×2, same site:** replace only the old operand fragment by
  `2 * ms[s].get(w, sp.S.Zero)`. This is the RESCALE of the exact existing `ms` addend in the `ds` construction.
- **DEAD print set, exact:** `N6RC_CARRIER_EULERIAN`, `N6RC_CARRIER_MATERIAL`, and
  `N6RC_CARRIER_BRIDGE_RESIDUAL`.

The earlier whole-operand `ms`/`es` exchange is COEFFICIENT-class scalar multiplication by `-1`; it is not a FORM
knife and is not a harness variant.

---

## Builder bounds

The builder's sequence is: build the four named harnesses, run each once as specified, write the four literal
transcripts, report the produced paths, and stop. Do not invoke another model, reviewer, agent, watcher, comparator,
or orchestration step. Do not commit. Do not alter any file outside the named harnesses and transcript paths.

## CHANGE-LOG — review finding to directive fix

- Global result-language leakage and incomplete self-tests → payload-only contract plus per-FORM identity, same-site
  ×2, and named one-sided print sets.
- WL K_carrier DEAD overreach → exact DEAD list; `R_COV_INCREMENT` and `ACTUAL_CONTROL_PARAMETERS` removed from it.
- WL diagonal sampler probe and uncovered split construction → sampler mutation removed; one-site `K_split_route`
  added at `buildCase.sourceChannel`.
- Covariance circularity site → moved from `predicted_amplitude` to the ordered amplitude block in `run`.
- Covariance rank site and omitted map records → rank-image identification inside `image_of`; `FROZEN_PHI` and
  `PHI_DOMAIN_CENSUS` added to extraction.
- Diagnostic E_W classification → sign flip classified COEFFICIENT; material E_W row-family deletion is the FORM
  knife, with fixed-domain completion left intact.
- Diagnostic advection and control extraction → advection mutation removed; both native three-object control groups
  are printed for `TILT` and `N4_ADVECTION`.
- Diagnostic slot discretion and DEAD surrogate → `delta_p_plus` pinned with literal fragments; named one-sided
  operand prints replace the dimension-control surrogate.
- Reconcile normal scope → temporary `build_material_face_source` wrapper is confined to
  `build_material_carrier`, with builder and cache restoration; Eulerian/global/`material_inverse_transpose` sites
  forbidden.
- Reconcile source return type → `n.build_increment` result explicitly unpacked.
- Reconcile operand classification and vague RESCALE → one `ms` read collapses onto `es`; the coefficient companion
  doubles that exact `ms` addend at lines 211–212.
- Isolation and builder scope → fresh subprocess/kernel, temporary sibling-first imports, build→run→report→stop.
- Output-contract revision → every CANONICAL, unablated-copy, FORM, NO-OP/IDENTITY, and COEFFICIENT ×2 record now
  prints only compact PIT nonzero fingerprints, per-object SHA-256 digests, small scalar/key-aligned tables, and
  bounded guards. Full symbolic objects/differences, arithmetic DAGs, PIT sample matrices, and raw engine stdout are
  forbidden even on failures; deliverable transcripts are KB-scale.
- Coverage verification against all four engines → every certified and DEAD arithmetic object in the diagnostic and
  reconcile harnesses is inserted into `objects` and passed to `n.pit`; the five covariance arithmetic objects
  (`R_COV`, `SOURCE_ACTUAL`, `SOURCE_PREDICTED`, `R_COV_INCREMENT`, `R_COV_CONTROL_DELTA`) are likewise passed to
  `n.pit`; and every WL arithmetic object named by the primary/DEAD sets is inserted by `addNumeric` and processed by
  `probeCase`. The only named non-PIT objects are covariance/WL `FROZEN_PHI` and `PHI_DOMAIN_CENSUS`: the SymPy
  worker must compute the minimal `n.sha(obj)` digest before transcript projection, while WL uses each object's
  `emitted_object_sha256` from `RUN_PROVENANCE`. Thus every certified and DEAD object has a compact PIT fingerprint
  and/or digest, with no uncovered object and no symbolic fallback.
- WL budget revision → the engine exposes no case-selection or draw-count environment/parameter knob. Every WL
  CANONICAL/copy/knife run therefore patches only its `/tmp` source: in the unique top-level driver `Do` containing
  `buildCase[case]`, replace the four-case iterator by
  `{case, {{"LAB_HELD", "RHOBR_CONSTANT"}}}`; in the unique `probeCase` draw assignment, replace the adaptive
  eight-draw expression by `drawCount = 4;`. Static source tracing verifies that the first edit only selects the
  argument to unchanged `buildCase` and the second occurs after per-case object construction/compilation, so retained
  per-case symbolic objects are unchanged. Both edits are declared `BUDGET_RESTRICTION / NON-KNIFE`, source-diffed
  and hashed separately from every frozen knife, and every fresh serialized kernel remains under
  `timeout --kill-after=5 600`.
