**NEEDS REVISION**

This is a source-only read. I read `method.md`, `implementation-guide.md`, all of `worker.py`, the `inputs.json` pins and resource block, `launcher.py`, and the `Journal` helper in `raw-helper.py`. I did not open the saved evidence JSON, `supervisor.py`, `shared-guard.py`, the completion hook, `tooling-tests.py` or `evidence-guide.md`. Every statement about saved-operand contents below is therefore a statement about what the worker assumes, not about what the files contain.

## Scientific / validation blockers

**B1. Controls do not enforce the support the method requires (`worker.py:726-735`, `760-764`, `770-783`, `794-800`).**
- `choose_source(face,(1,0),spatial=True)` accepts any jet with `spatialOrders[0]>0`. It does not require `x in item['field'].free_symbols`. Its nonzero test, `field.subs(x,0)*wave(...)`, is a pointwise value of the coefficient field, not a certificate on its transform.
- If the grade-(1,0) source field were constant, its transform is `b·δ(k-p)`. The p=3/2, k=2 control point would then be off-support, so the "p→k" sensitivity would be inert while still reported as responsive.
- The same applies to the normal consumer10 used for `q(l)→q(r)` and the lower-sign control (r=30/13 against l=2 or 3/2). Lines 762-764 require the pressure consumer to be constant, which is right. But `cn=consumer_fields[...,'d_w_delta_p_'+face][1,0]` is only checked for `cn.subs(x,0)!=0`. A constant `cn` would force r=l.
- `cn.subs(x,0)` can also reject a valid consumer whose profile derivative vanishes at 0, such as `m'`.
- **Minimum correction:** make eligibility require a nonconstant field for any control that uses off-delta support (source10 for p≠k, consumer10 for r≠l). Record the formal transform tag coefficient at the stated support. Drop the `x=0` pointwise value as the nonzero certificate.

**B2. The actual-row placeholder identity does not cover the normal factor or per-component maps (`worker.py:103-106`, `635-663`).**
- `address_sum_for` sums `consumerOriginal*responsePlaceholder*sourceOriginal*sourceAtom`. Tokens are free symbols keyed by (face, slot, b, component, c, atom).
- The comparison therefore checks only consumer×source bookkeeping and η/σ grading. It never uses `normalOriginal`, `responseOriginal` or `responseMap`.
- For the normal slot, `i·f·q(l)` is not in the substituted row or the sum. A wrong sign, a q(r) in place of q(l), or a wrong component map would pass this check.
- The method asks for a formal injection of the whole response into the actual rows.
- **Minimum correction:** substitute each token by `normalOriginal*responseOriginal`, and by the mapped `responseCoefficient*normalMultiplier`, on both sides. Then check the row coefficient against the address sum, or add an explicit residual per address of mapped against original under its map.

**B3. The broad pressure-name census is literal-AST-only (`worker.py:133-147`, `250-255`).**
- It matches `Symbol(...)`/`Function(...)` calls only when `func` is a bare `Name` and the first argument is a string constant.
- `sp.Symbol(...)`, `symbols('a b')`, `Function(var)` and f-strings are silently not hits.
- It also treats the row as one child if the top level is not an `Add(...)` call.
- The only guard is `len(entries)==old['totalChildren']` plus agreement with the saved selected list, which are both inherited.
- **Minimum correction:** also scan the raw child text for the substrings `delta_p` and `d_w_`. Refuse if the substring count exceeds the AST hits, and refuse any Symbol/Function/symbols call that has a non-literal or attribute-form constructor.

## Smaller validation gaps

- **Normal-slot ablation is missing (`worker.py:266-270`).** Ablation covers only `delta_p_plus`/`delta_p_minus`. The two `d_w_` slots are checked only through the affine reconstruction. Add ablation joins for them, since the method treats all four slots as independent.
- **The native reference assignment is not joined (`worker.py:406`, `424-425`).** The adapter sets `reference=sign*W_0/2` by hand. The original `reference` assignment text is saved but never evaluated or equated to it. Join it, or evaluate that one assignment.
- **The excluded (2,1) remainder is not required to be present (`worker.py:301-302`).** The check is `set(hg)<={(2,1)}`, which also passes an empty set. Require (2,1) present and nonzero so the old unprojected remainder is actually preserved.
- **Source dimensions are recorded but not enforced (`worker.py:336-338`).** `coefficientDimensionFromSource` is computed and never required to match anything. The "source-bound unit rules" are therefore emitted, not checked.
- **Mixed-iteration coverage is not checked.** The worker treats `rc['mixedIteration']` as one component. Nothing verifies that it contains both native first-shape assignments including the trace subtraction. This is a coverage limitation until the saved object is inspected.

## Pure tooling

I found no blockers in the gate, launcher or manifest joins I read:
- **Gate:** `verify_gate` binds worker, manifest, source pins, guard, supervisor, launcher, review record, method and authority hashes. It requires clear literal verdicts from both reviewers, and `durationLimits` must be null.
- **Launcher:** it arms the hook first. The 30 s wait is a startup handshake, not a compute deadline. It hardcodes the command through the pooled guard (4 GiB, zero swap, tasks 32). It snapshots sources.
- **Failure and posthash handling:** write-before-guard through `Journal.zero`, which emits input and raw before `require`. It uses `BaseException` handling, saves a failure record, and writes posthashes in `finally`.
- **Not read:** I could not confirm that supervisor/guard enforce the memory, no-deadline and `RuntimeMaxUSec=infinity` containment.

## Points the prompt asked about, assessed

- **Typed factors:**
  - `Rprod`, `E` and the closed-density whole tag are kept separate. `typed-factor-relation` joins `Rprod=qi*qo*E`.
  - The density-minus-factor comparison is not made.
  - The inherited zero-return pairs are required as literal zero and are not recomputed.
  - The isolated-factor and full-density joins sit on different operands.
- **Independent-D adapter:**
  - It uses only the four pinned `build_face` assignments.
  - It checks the whole trace matrix action, including the zero (1,2) injection, both face signs, and the excluded (2,1) term as `-height*i*f*qo*P`.
  - No old producer is called.
  - It is limited by the reference-assignment gap above.
- **Epsilon versus nonzero:**
  - `epsilonCount` is 0 only for exact zero and otherwise 1.
  - The status string `FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED` keeps homogeneity separate from certified nonzero.
  - The `J.zero(v, eps*cancel(v/eps))` check at line 281 is tautological, but the `require` on `eps` not in `v/eps` that follows does the real work.
- **Sheet:**
  - The positive outgoing branch and branch conditions are checked at the three control momenta.
  - H's own bound variable is recovered from `variableChange` and kept separate from J's t and D's td.
  - Unknown values are not converted to zero, because `nonzero` and `select_certified_candidate` require certified results.
- **Zero-profile and constant-height scope:**
  - The formal zero-profile source/consumer/tag reductions are checked.
  - Declaring the nonzero constant-height full-response reduction `NOT_COMPUTED` is sufficient for this bounded inventory.
  - It is correctly not claimed as global operator applicability.
  - The `constantCoefficientSource` records at lines 679-684 are emitted but never joined to anything, so they are not evidence.

## Runtime evidence still needed

- The completed run's actual row partition and broad-scan hits.
- The per-address token reconstruction residuals, including normal-slot tokens, once B2 is fixed.
- The control candidate tables showing nonconstant fields and the transform tags at the support points.
- The posthashes, `failure.json` or `checks.json`, and the guard and supervisor logs.

No execution gate exists yet. A source verdict cannot accept the runtime result.