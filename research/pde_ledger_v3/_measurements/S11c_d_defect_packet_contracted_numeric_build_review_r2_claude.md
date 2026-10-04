**Verdict: NEEDS REVISION**

I found one blocker. The preamble cannot finish with the supplied restore routine, so no numerical work is reachable. I found no correctness fault in the numerical logic itself. I ran nothing; every statement below comes from reading source and operands.

## What I read

**Full reads:**
- `build.md`, `method.md`, `evidence-guide.md`, `review-prompt.md`
- `worker.py`, `prepare.py`, `validation.py`, `geometry.py`, `numeric.py`
- `request-index.py`, `evidence-store.py`, `launcher.py`, `restore-library.py`
- `original-source/S11c_d_defect_packet_geometry_lib.py`

**Partial reads:**
- `original-source/…preflight.py`: lines 236–280 and 296–351.
- `original-source/…weak_composition.py`: lines 150–200.
- `fourier_lib.py`: grep only.
- `containment-helper.py`: `containment()` and `verify_gate`.
- `numeric-source-contracts.json`: grep only.

**Native operands read:**
- `new-contraction-definitions.json` (full)
- `whole-tags.json` (full)
- `new-runtime-J-arithmetic-input.json` and `new-runtime-D-arithmetic-sum-input.json` (full)
- `numeric-factor-adapters.json`: lines 1–160 only
- `summand-8346-return.json`: lines 60–180
- `pressure-addresses.json`: the 8347 entry and the start of 8348, plus a grep of every `waveMultiplier`
- `physical-plan.json`: the cs and kappa entries
- The three rule files: keys, format and sizes only
- `input-manifest.json`: lines 1–110 and 15855–15906
- Test-log pass/fail lines

**Coverage limits:**
- I did not read `shared-guard.py`, `supervisor.py`, `completion-hook.py`, the test bodies, or the roughly 3,000 other JSON operands. These include the per-field proofs, per-address contraction files, `tail-plan.json` and `absolute-bounds/*`.
- Data joins for all 20 entries are therefore judged from code, with spot checks on 8346, 8347 and 8350.
- I cannot hash files, so I did not confirm that `restore-library.py` is byte-identical to the pinned `…contraction_lib.py`.

## Blocker

**B1. `restore_scalar` cannot restore operands the new joins pass to it.**

`restore-library.py:75-97` accepts only `Integer`, `Rational`, `Symbol`, `Add`, `Mul` and `Pow` calls, bare `I`, and int/str/bool constants. Anything else raises `ValueError`.

The new code applies it to operands that contain `sinh`, `pi` and `Function('…')(…)`:
- **`validation.py:46`** restores the `left` of `saved/inner/new-runtime-J-arithmetic-input.json`. Its srepr contains `Pow(sinh(Mul(Integer(5), pi, …)), Integer(-1))`. The D-arithmetic-sum `left` has the same content. This raises `ValueError: allowed scalar constructor sinh`. It is the first refusal, reached from `prepare.py:70`, before the Gaussian, geometry, tail and numeric stages.
- **`validation.py:56-57`** repeats the same two restores.
- **`validation.py:67`** restores `definition['density']` for Jwhole and Dwhole. In `whole-tags.json` (lines 88–89 and 110–111) these contain `sinh` and `pi`.
- **`validation.py:79, 82, 83`** restore `responseCoefficient`, the `permittedCalls` pairs and `actualCalls`. These contain `Function('Jwhole_plus')(…)` and `Function('common_outgoing_q')(…)` (adapters lines 15–16 and 63–64). Here `n.func` is a `Call`, not a `Name`, so it fails with "constructor call only".

The accepted contraction run used `restore_scalar` only on `a`, `mu`, `kappa`, `W`, `L`, adapter templates and CX/CY bounds (`contraction.py:138-140, 228, 272`). None of those contain these constructs. No tooling test exercises `restore_scalar` on native operands (`tooling-tests.py` has no restore or `sinh` coverage).

The result is a deterministic, fail-closed refusal that is preserved as `FAILED_PRESERVED`. The instrument cannot currently evaluate J or D.

**Fix direction:** add a restore path as strict as the original. It needs an explicit allowlist (`sinh`, `pi`, and `Function(name)(…)` with named functions and exact arity). It must be pinned and re-reviewed, and it must not weaken any `exact` zero check.

## Verified correct by reading

- **Contractions:** `inner_node` (`numeric.py:236-254`) matches every integrand in `new-contraction-definitions.json`: Y0, X1, X1T, X2T, X02T, C0–C2T, Y0CT, Y1CT, Y0C, Y1C.
- **Combination:** `combine` matches the saved outer formulas, including `Dr_wrong_root`. The preamble checks this exactly (`prepare.py:154-158`).
- **Clipping:** J/Dh/Dq clip the input k, Dr clips the output l, and the wrong-root mutant clips input with unclipped output.
- **Profile and branch:** the profile (A even, sinhc branch) and the `q` branch match the original ASTs.
- **Transform conventions:**
  - **X:** ν = k−p0 with phase e^{−iν x_u}.
  - **Y:** ν = p0−l with phase e^{−iν·5/2}, which equals e^{+i(l−p0)·5/2}. This matches `fourier_lib.py:253`.
  - **Normalization:** X carries 1/(2π) and Y carries none.
  - **Coefficients:** Xmultiplier = b·P·iⁿ and Ymultiplier = c, so alpha = b·c·P·iⁿ is applied once.
- **Derivative and moment identities:** I expanded `moment_polynomial` for n=2 by hand and it reproduces (i(ν+p0))². The B derivative mutant replaces Q with the constant (i·p0)ⁿ, so only M0 is used. Route A uses p0ⁿ.
- **Selectors and joins:** original `argumentDerivative` is required to be the integer 0. The Y interface membership, epsilon binding, unit-coefficient J/D template selection, and the explicit H=0 formal tag match the stored templates. These are `packet_H·… + packet_J` and `packet_D`.
- **Quadrature:**
  - **Route A:** the squared-map Jacobian is 2·half·u, with weight/2 for the (n+1)/2 map.
  - **Route B:** the physical map and Gauss weights are embedded correctly. The saved 15-point Kronrod and 7-point Gauss nodes carry identical tuples (spot check).
  - **Errors:** product and sum propagation are correct. The per-panel |GL24−GL48| indicators are summed and fed through `combine`.
- **Budgets:** the B pointwise target ε(1+1/|q|)/(16(3M+4)) and the weighted nested guard ε/4 are implemented as specified. The leaf ranking uses the binding guard in `outer_decision`.
- **Geometry:**
  - **Slabs:** six affinity slabs. I checked the d_g table and the clip branches by hand for all six.
  - **Coverage:** full wings, and every graph-pair intersection inside each slab.
  - **Audit:** true max/min window at three points per slab, and a no-crossing recheck.
  - **Wing mutants:** two per plan, refused with the window error.
- **Tails:** the new AST joins for `exponential_moment`, `weighted_tail` and the four contribution expressions match preflight lines 322–338. The enlarged-window J middle term includes the H overcount. I checked the Dr, Dh and Dq termwise majorants by hand.
- **Persistence:**
  - **Requests:** requests are keyed on the full canonical descriptor, with route, purpose and precision joins.
  - **Failure handling:** PENDING requests refuse reuse, the exclusive-create store refuses reopening, and partial node batches are committed on a callback failure.
  - **Disk and records:** the disk reserve and the 8 MiB cap are checked before every write.
  - **Final receipts:** posthashes and the final receipt handle failure before the leaf table exists.
- **Launcher:** the hook is armed before launch, there are no deadlines, and the containment values match (4 GiB, swap 0, pids 32, one CPU, no CPU limit).

## Non-blocking findings and shared limits

1. **Y has no independent derivative formula.** At `numeric.py:113` both routes use `amplitude=base`, so the Y A/B check is a precision check only. I verified the sign against `fourier_lib.py:253`.
2. **The Gaussian identity is not run on the real code.** The preamble check (`prepare.py:144-150`) is separate sympy code, not `moment_polynomial` or `gaussian`. Only the per-node 1e-24 A/B gates (`numeric.py:213-223`) cover the implemented recurrence.
3. **Whole definitions carry no face.** `whole-tags.json` has only H, Jwhole and Dwhole. Plus and minus bind by tag-name string (`validation.py:89`, `weak_composition.py:172-173`). This is an inherited assumption that both faces share one density.
4. **The absolute-majorant joins are partly prose.** `prepare.py:168-177` compares hard-coded polynomials to the saved certificates. The inequality from each primitive's numerator to those polynomials is prose (`prepare.py:180`). The accepted gap proofs are inherited and not rederived.
5. **Some checks are tautology-leaning.** `new-family-complete-summand` and the J `new-whole-selected-density` confirm structure and unit coefficient, not an independent factorization.
6. **The geometry audit reuses planner functions.** It shares `graphs` and `branches` with the planner (`geometry.py:109`). It is independent only for window, coverage and crossings.
7. **Storage will grow fast.** Every request descriptor embeds the full `contexts` (`numeric.py:169`). The raw descriptor is stored uncompressed in `request_index` (`request-index.py:134`), and the descriptor is re-canonicalized on every call. This is a cost issue, not a correctness fault. I cannot tell whether the 20 GiB reserve is reached before completion.
8. **K−G and inner totals are only implicit.** Outer K−G and inner totals are persisted implicitly rather than as separate named fields.

## Mandatory runtime checks (after revision)

- **Preamble:** every `exact` join returns the literal zero, and the journal shows the restored operands.
- **Gaussian gates:** the per-node full-value and stripped gates (1e-24) pass at every node.
- **Adaptive loops:** no stagnation in the inner or outer B loops.
- **Resources:** RSS stays under the 4 GiB address-space limit, and `request_index` plus record bytes stay clear of the disk reserve.
- **Failures:** any A nested-ε/4 failure, control non-response or comparison miss is preserved as is, with no retry or retuning.

## Scope limits

- This is ordinary J plus the three ADDED direct terms only, for two carriers and two windows.
- H, flat, height/PV and slope remain pending. J is not full native mixed, and J+D is not a packet value.
- The shared analytic Gaussian identity and the empirical indicators are not independent Fourier quadrature or rigorous error bounds.
- Controls are finite-window only.
- Gamma units and the accepted algebra are inherited.
- Source clearance is not runtime acceptance. Total work, storage and runtime stay unknown.