CLEAR FOR THIS BOUNDED SOURCE-COMPOSITION BUILD

This is a source-only read of `worker.py`, `inputs.json`, the method and guide, the runtime-source files (launcher, supervisor, guard, hook, helpers) and the saved JSON I opened. I ran nothing. I found no scientific blocker. The items below are non-blocking.

## Blockers
None.

## Verified at source level
- **Typed factors.** `Rprod`, `E` and the whole closed-density tag are kept as three separate objects (`worker.py:493-525`).
  - The isolated-factor pairs (`direct/*-reference-factor-*`) and the full closed-density pairs (`*-actual-closed-*`) are consumed separately.
  - The inherited zero returns are asserted, not recomputed.
  - Nothing compares `density − Rprod`.
  - `Dwhole` enters once, with a unit multiplier (`:1016`).
- **Independent-D adapter** (`:527-565`). Only `reference`, `extension`, `jet_transfer`, `normal_jet` and `reference_pressure` are executed, each from its pinned `build_face` text.
  - The adapter enforces the native `kernel_apply` argument order and refuses a nonzero diagonal, second or source.
  - It gives `normal_jet == i·f·qo·P` for both signs.
  - The retained reference is `P`, and the excluded (2,1) term is isolated.
  - `T·Δ = Δ` over the full trace matrix confirms there is no [1,2] injection.
  - `reference_location` is executed from the native assignment.
  - No producer, triangular solve or response function runs.
- **Source and chemical joins.** The checks run in this order, before any success claim (`:102-121`, `:432-453`):
  1. Native chemical and density leaves are joined to the saved raw operands.
  2. The chemical amplitude, its domain record and the epsilon division are joined.
  3. `raw.subs(stage2)` then `bind` is compared with the saved combined source.
  4. The velocity normalization is joined for each face.
  - `densityRestoredAndJoined` is emitted only after those residuals pass.
- **Row census and affine reconstruction.**
  - The broad `delta_p`/`d_w_` scan runs on the constructor text with substring accounting. Attribute, dynamic or hidden forms are refused.
  - Child hashes are compared with the saved census, and all four slot ablations and the affine reconstruction are checked.
  - Grade tables use the exact quotient recurrence with a nonzero constant denominator at zero. Residuals at every retained grade are tested against zero, with the excluded remainder preserved. Epsilon is checked to appear once.
- **Placeholder substitution against the address sum** (`:821-849`).
  - Truncating to `b+c ≤ (1,1)` and the address triple set `a+b+c=g` select identical sets.
  - Tokens are keyed identically in both, so this is a real second Taylor expansion, not a restatement.
  - Normal-slot signs and momentum arguments are tested separately by `join_address_factor` (`:755-777`).
- **Coverage.** All 16 triples per row, face and slot are enumerated with explicit zero reasons (`:817`). Epsilon is recorded as homogeneity, with status `NONZERO_NOT_ASSERTED`.
- **Off-delta eligibility** (`:910-938`). The tanh-polynomial certificate establishes only that the field is nonconstant. It never asserts a Fourier value, and unknown values are excluded because the selection requires `is_zero is False` and `is_finite is True`.
- **Controls** use actual addresses with explicit supports. The depths at the control points check out as `q(3/2)=2`, `q(2)=3/2` and `q(30/13)=25/26` at `cs²=10/7`. Every movement is labelled a formal tag coefficient.
- **Saved data I opened** is consistent with the code:
  - The consumer10 normal coefficient is `∝ eta·w1_profile`.
  - The pressure consumer00 is a nonzero constant.
  - `beta = (30+9i)/109` matches the saved isolated-factor pair.
  - The speed inventory finds no `c_s` symbol in the raw operands.
- **Execution chain.**
  - `verify_gate`, `verify_invocation`, the argv tail check and the review-record route all run before `mkdir` and the SymPy import.
  - Evidence is written before each guard through `J.zero`/`J.emit`.
  - The post-hash pass covers inputs, copies and sources.
  - The launcher is hook-first, with a pipe handshake that is not a deadline.
  - `containment()` enforces the memory, swap, task and CPU limits and the absence of a CPU deadline. The guard enforces `RuntimeMaxSec=infinity` and `Restart=no`.

## Scope decisions
- The nonzero constant-height reduction is honestly marked `NOT_COMPUTED`. That is sufficient for this inventory, because the zero-profile reductions are formal and no constant-end claim is made.
- The H bound variable and its `l − hs` change of variable are joined (`:673-677`). `J` and `D` bound variables are documented but not free arguments.
- Controls cover THETA and E_W only. The U rows carry no pressure slots.

## Non-blocking items
**Validation**
1. **Excluded-(2,1) height is not joined to the native trace.** `worker.py:555-563` is arithmetic on `slot['savedHeight']` and is tautological in `height`. A sign error in the saved height, such as the lower face's `−eta·Hhat`, would not be caught. The retained [1,1] result is unaffected.
   - Correction: add `J.zero(savedHeight, trace['restoredTrace']['height'].xreplace(dict(actualHeightMap)))` per face.
2. **Address records omit dimension and normalization.** Method §4 lists them, but `addr` at `:793-811` omits them. Per-jet dimensions exist at `:465-468`.
   - Correction: copy `jetDimension`, `requiredCoefficientDimension` and the `1/(2π)` normalization into `addr`.
3. **The assembled-mixed-reference residual is circular** (`:730-736`). The tags are created and back-mapped in the same scope. The real content sits in the per-component maps and the census join at `:703`. Do not cite it as independent evidence.

**Tooling**
4. **No fallback for missing control candidates.** Per method §6, absence should be recorded and an applicable entry used. At `:960` and `:993` the worker instead fails with evidence preserved. Applicable entries appear to exist in the saved data, so I expect no impact.
5. **Pre-emit requires.** `require` at `:399`, `:529` and `:639` precedes emission of the derived data. The inputs are saved copies, so this is low impact.
6. **Exit code and status disagree after an integrity failure.** `:1073` sets `code=1` after a COMPLETED status. The `integrityFailure` flag marks it.

## Coverage limits
- I did not open every saved file. The following joins are verified in code only, and a mismatch would fail closed:
  - `units-and-scope`
  - the raw row files
  - `height-left-plus-trace-*`
  - `right-height-PV-operands`
  - the function-routes records
- I assessed no scientific content beyond the equations the worker joins.
- Source dimensions after binding are required units only. Global momentum, test-space and grazing composition remain UNRESOLVED.

## Runtime evidence still needed
- Exact-zero returns for all `J.zero` residuals, covering:
  - the 320 coverage entries
  - the actual-row versus address-sum checks per row, face and grade
  - the component, factor, source, chemical and density joins
- Structural equality of the saved Piecewise branch conditions with `expected_conditions` (`:643`).
- `sympify` of the constructor-text children, and the `exec` of the adapter assignments.
- A clean speed inventory and the `L_W=10` join.
- Control selections and `nonzero` certificates for both faces.
- Clean post-hashes, guard and cgroup logs, and an empty stderr.

None of this is result acceptance.