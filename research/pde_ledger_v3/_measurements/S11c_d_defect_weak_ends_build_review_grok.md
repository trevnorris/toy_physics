I'll start with the build guide and the input packet so the assessment stays on the actual translated weak-end implementation.I’ll read the method, worker, launcher, and packet maps next, then follow those into the saved operands.The method and worker are large; I’ll keep reading the rest of the worker and the packet layout.I have the worker’s full flow; next I’ll inspect helpers, tests, and the saved schemas those joins actually consume.The helper path in the packet is still unknown, so I’ll look at the packet inventory and native sources next.I’ll keep tracing the packet names, tests, helper, and the height/normal saved operands.I’ll inspect field polynomials, local cells, factor maps, and binding context against those joins.The lower-normal control cannot bind a surviving height address on the saved operands, so this translated weak-end build is not ready.

The phase, PV/contact, native half-heights, Holder joins, T-versus-tanh field endpoints, and local P(±1) inheritance that I could inspect line up with `method.md`. The readiness failure is in the required addressed controls.

## Blockers

### 1. `omit-lower-normal-sign` has no surviving saved address (scientific validation / readiness)

`worker.py` only collects control candidates from plus-end **height** terms (`responseGrade==(1,0)` and `term!=0`), then keeps `omit-lower-normal-sign` only for `face=='minus'` and `slot=='normal'`. Absence is fatal (`require(chosen is not None, ...)`).

```451:452:input/worker.py
                    if side=='plus' and rg==(1,0) and term!=0:
                        possible_controls.append((ar,source[side]*consumer[side]*wave,rr,term))
```

```478:498:input/worker.py
    for kind in ('omit-height-contact','reverse-translation-phase','omit-lower-normal-sign'):
        ...
            if kind=='omit-lower-normal-sign' and (ar['face']!='minus' or ar['slot']!='normal'):continue
        ...
        require(chosen is not None,'responsive applicable actual-row control '+kind)
```

On the saved inventory that intersection is empty.

- Normal-slot consumers are the `d_w_delta_p_*` fields. For every row those retained grade-(0,0) pieces are `0`. The only retained nonzero piece is grade `(1,0)` (`eta*w1`).
  - `saved-inputs/inventory/THETA_BALANCE-d_w_delta_p_minus-split.json` retained `(0,0)=0`, `(1,0)=I*epsilon*w1*(6-20I)/436`
  - `saved-inputs/inventory/E_W_BALANCE-d_w_delta_p_minus-split.json` retained `(0,0)=0`, `(1,0)=-epsilon*w1/4`
  - U0/U1/U2 `d_w` splits are all zeros
- Height response is already grade `(1,0)`. Adding a grade-(1,0) consumer makes target `(2,0)`, which is outside `G` and outside the worker’s `targetGrade in G` join (`worker.py` 144–145). Those triples are not in the 16 ordered triples actually present.
- The minus-normal height route that **is** present has a zero consumer: address **10023** (`THETA_BALANCE`, `face=minus`, `slot=normal`, `component=NATIVE_HEIGHT`, `consumerTransform.coefficientId=9937a9eb…`, `status=EXACT_ZERO_CONSUMER`).
- The surviving minus-normal consumer is address **10062** (`consumerGrade=[1,0]`, field `9a4bf708…`, `FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED`), and that address is `NATIVE_FLAT`, so it never enters `possible_controls`.

`method.md` (lines 222–232) requires the lower-normal control through an actual surviving lower normal-slot address, and says to persist a coverage gap if none exists. This worker still requires a selected control. The 22 preparation tests never inspect control-address existence.

This is a changed scientific-validation claim: the build asserts three addressed row controls, and one of them cannot run on the saved schema.

### 2. Reverse-phase control does not implement the required left/right exchange (implementation)

`method.md` line 223: reverse the translation phase and verify that left/right height limits exchange.

For `reverse-translation-phase` the worker sets `mutated=sp.S.Zero` on the plus-end height term only (`worker.py` 482) and records `'height minus/plus exchange: right height becomes zero'`. It never builds the reversed-phase minus-end term `H_+=W/2` (equivalently contact plus flipped Dirichlet) or inserts that term into the minus-side cell.

Plus-end height **can** be bound: address **8034** (`THETA_BALANCE`, plus, pressure, `NATIVE_HEIGHT`, source `5/218+3I/436`, consumer `-10/109-3I/109`, `FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED`). Contact omission on that route is at least well-posed. The phase control on the same route is only half of the required exchange.

## What else was checked (no additional blockers found)

**Translation / PV / heights.** Phase increment is `(l-k)a` from `l*x-k*y` after the simultaneous shift (`worker.py` 320–321), matching `exp(i(l-k)a)`. The PV split keeps the oscillatory diagonal `E*A*fzero*chi/Q` (`328–329`). Contact `W/4` plus Dirichlet `i*pi*A0*sgn(a)` with `A0=1/(2pi)` recovers `H_-=0` and `H_+=W/2` (`330–331`). Uniform-bound prose forbids differentiating `exp(iQa)` (`345`). Both PV certificates set `normalCoefficientBoundedDirectly: true` and `qTimesTestAssumedC1: false`; on-diagonal they match `Bh` and `±I q Bh`.

**Native constant height.** Plus height `eta*w1/2` with jet `I*q_o`; minus height `-eta*w1/2` with jet `-I*q_o`. The product is `I q eta H_±` on both faces. Affine slots are `-face_affine_height*face_jet_slot+face_physical_target`. Grazing uses `cancel(candidate.subs(q,0))` on closed `q+beta` expressions. `eta^2` is kept only as `1-product**2`. Auxiliary `delta` is not a physical source law.

**Field / local joins.** Pressure `field-*-polynomial.json` stores numerator and denominator separately (`f588c943…`: numerator `T*(-80-24I)-80-24I`, denominator `6976`). Endpoints use certificate/`P0` complete quotients (`worker.py` 182, 275–284) and join `P(tanh(x/10))` to the saved reconstruction **right** (`reconstruction-input.json` for `f588c943…`). Reconstruction returns have `cancelled: 0`. Local cell U0/`u_1`/grade `(1,0)` has `P(T)=-T^2/60-9T/2-269/60` with saved `leftEndpoint=0`, `rightEndpoint=-9`, matching `T=±1`. Local identities are inherited as published zeros.

**Symbols / grades.** 400 local cells, 34 fields, 13260 addresses, 16 ordered triples, 200 end cells, pairing factor `2pi` outside `E`, source waves reduced at original `p` (`composition_p` in `new-pressure-wave-arguments.json`), output normal after diagonal is `1` or `±I q` (`factor-5` plus `I*q(l)`, `factor-12` minus `-I*q(l)`). Slope/mixed responses are tagged as ordinary-kernel Riemann–Lebesgue zeros; source/consumer cross grades still enter through flat/height target sums. Vanishing mixed kernels are not treated as pointwise absence.

**Runtime.** `main` runs `containment()` before `import sympy` (`519–521`). Helper import is AST-sliced (`definitions` + `HELPERS`). Launcher arms `codex_job_watch.py` before science, uses pooled `s11c_guarded_run.py` around `S11c_d_end_normalization_run.py`, `durationLimits: null`, `automaticRetry: false`, `Restart` implied no by standing message. `inputs.json` is `PREPARED_NO_SCIENCE_OR_GATE`; no READY gate is in the packet, which matches the stated pending assessment.

**Scope.** Real `omega=3` in `extended-binding-context.json` numeric; historical `physical-input.json` still has `omega: "1"` and `c_s0: "10"`. Control point `p=1`, `q=2`, `cs^2=180/101` satisfies `9/cs^2-(1/5)^2-(1/10)^2-1=4`.

## Inspected files / regions

- `input/build-guide.md` (full)
- `input/method.md` (full)
- `input/worker.py` (full, 543 lines)
- `input/launcher.py` (full)
- `input/inputs.json` (scope, `savedInputs` head, `sourcePins` tail, resources)
- `input/tests.py` (full)
- `input/physical-input.json`, `input/execution-authority.json`
- `input/saved-input-map.json` (alias head)
- `input/address-review-index.jsonl` (U0 0–49, 78–79, 156–157, 624–638; THETA 7956–7963, 7995–7996, 8034–8035, 8073–8074, 8112–8114, 8151, 8190, 8219, 8299; minus-normal 9945–9947, 10023–10024, 10062–10063, 10101–10102, 10140–10141; E_W 10725–10727)
- `input/address-representatives.json` (header + first full record)
- `saved-inputs/weak/`: duality, plus/minus PV certificates, analytic-conclusion, envelopes, `all-coefficient-certificates.json` (head + `f588c943`, `aa5461ab`, `6b57398f`), field `79ea957c`, `f588c943`, `18fdaed1`, `9a4bf708` polynomial/operands/reconstruction/derivative-class
- `saved-inputs/full/`: complete-weak-assembly, extended-binding-context (physical + numeric omega), all-local-cells (cells 0 and grade-(1,0)), new-pressure-factor-arguments (factors 0, 1, 2, 5, 6, 11, 12), new-pressure-wave-arguments (`e_W`, `e_W_d1`)
- `saved-inputs/reference/`: plus/minus native traces, height constants, final slots, first-shape input, restored-profile-arguments
- `saved-inputs/saved/reference/retained-response-census.json` (head)
- `saved-inputs/weak/global-profile-envelope.json`
- `saved-inputs/consumer/`: U0, THETA_BALANCE, E_W_BALANCE
- `saved-inputs/inventory/`: fields.json (head), U0-delta_p_plus-split, THETA/E_W `d_w_delta_p_minus` and THETA `delta_p_minus` splits
- `input/native/c2.py` (header through coordinate setup; not the full constructor body)

## Uninspected material relevant to the verdict

- Full 13,260 address **bodies** and the 17 MB ordered-address arrays (index + representatives + sampled records only)
- Remaining ~398 local cells and ~30 field files beyond the sampled constant, T-linear source, and T-linear consumer
- `Journal` / `decode` / `containment` text from `helperSource` (`S11c_d_defect_raw_increment.py`): not present under the packet paths tried (`inert-helpers.py`, `helper-source.py`, `S11c_d_defect_raw_increment.py`, `source/...`). Stage/decode behavior is inferred from call sites.
- Full `native/c1.py` and the rest of `native/c2.py` (worker only AST-matches `phase1` / `phase` / trace assignment strings)
- Repo originals named in the manifest (provenance only)
- Any peer report

## Nonblocking limitations

- Riemann–Lebesgue, Dirichlet, and Schwartz-seminorm arguments are recorded as assessed analysis; they are not machine measure-theory proofs (`analyticTheoremsMachineProved: false`).
- Preparation tests are standard-library JSON/AST wiring. Packet `tests.py` still points at `S11c_d_defect_weak_ends.py` and absolute `savedInputs` paths; that is packaging, not physics acceptance.
- `verify_gate` / launcher READY pins are intentionally absent until this assessment.
- Exit 0 is not scientific acceptance.
- Analytic weak-end limits do not give pencils, radiation conditions, scattering, current, loss, drain, or calibration.

The first blocker is enough to stop this build: the required lower-normal control has no applicable saved witness, and the worker will die on `require(chosen is not None, ...)` instead of closing a recorded coverage gap.

NEEDS REVISION