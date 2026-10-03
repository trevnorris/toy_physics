## Report: translated weak-end BUILD review

**Verdict: NEEDS REVISION.** The mathematics in `method.md` and the worker's algebra look right. The blockers are in the control construction (B1–B2) and one missing join between the phase and the saved native convention (B3).

### Blockers

**B1. The lower-normal-sign control can never be selected, so the single authorized run fails deterministically.**
- `worker.py:451-452` adds to `possible_controls` only plus-end terms with response grade (1,0) and a nonzero term.
- `worker.py:481` then restricts `omit-lower-normal-sign` to face `minus`, slot `normal` addresses from that list.
- In the saved metadata, all 585 minus-face normal-slot addresses with response grade (1,0) have status `EXACT_ZERO_CONSUMER`. Example: addressId 2067, `address-review-index.jsonl:2068`. None are `FORMAL`, and none have an `EXACT_ZERO_SOURCE_JET` status.
- Plus-face normal-slot (1,0) addresses have no `FORMAL` entries either.
- The 72 surviving minus-normal formal addresses are 48 with `NATIVE_FLAT` response and 24 with slope/mixed response. Example flat address: 10062, THETA_BALANCE, consumer grade (1,0), `address-review-index.jsonl:10063`.
- So the candidate list holds no minus-normal entry, `chosen` stays `None`, and `require` at `worker.py:498` fails after the symbols are emitted. This is a fatal wasted run with no retry authority.
- `method.md` asks for "an actual surviving lower normal-slot address". The surviving ones are flat (0,0), where the normal factor is `i f q F0`. The selector needs to draw from those. It also needs a pre-check that the end value is nonzero.
- `tests.py` has no control-availability test, which is why this was missed. This is a selection/tooling defect, not a changed physics claim.

**B2. The reverse-phase control does not exercise the phase or the PV term.**
- `worker.py:482` hard-codes `mutated = 0` for the plus-end term.
- It does not recompute the Dirichlet limit with the sign reversed, and it does not show the left and right limits exchanging as method item 4 requires.
- The result only shows that deleting the plus-end term moves the cell. The minus-end limit is not derived from the decomposition.

**B3. The translation sign is typed in, not derived from the saved native phase.**
- `worker.py:320-321` builds the phase check from a typed `l*x-k*y`.
- `E=exp(I*Q*a)` at `worker.py:328` is also typed.
- `Fourier['sourcePhase']` and `['profilePhase']` are only AST-matched against `c2.py` (`worker.py:323-324`). They are never parsed into a symbolic phase.
- If the profile convention were flipped, nothing would fail. The contact would then go to the wrong end, giving W/2 where 0 belongs.
- I checked by hand that the saved convention gives `h(Q)=W/4 δ + W/(2i) PV A/Q`. That is the transform of the increasing profile `W/4(1+tanh)` with `exp(+iQa)` for `B(v_a,u_a)`, and it matches `exp(i(l-k)a)`. The result is therefore correct, but the instrument does not enforce it.

### Checked and found sound

- **Gate and runtime pins:** helper, supervisor, guard, hook, c2, physical input, worker, launcher and authority are all in `sourcePins`.
- **Scope and authority:** the manifest scope and `execution-authority.json` scope are identical, with one run, no deadline and no retry.
- **Gate and containment:** `verify_gate` and the paired-review check are consistent. Native 4 GiB `RLIMIT_AS` is applied in `containment()` before `import sympy`.
- **Launcher:** it starts the hook first, with a 30 s startup handshake only. The command matches the guard flags (4 GiB, pool, 32 tasks, no `RuntimeMaxSec`). The supervisor stage name is free-form, so `defect_weak_ends` is fine.
- **Field endpoints:** quotient-versus-numerator storage is handled. `certs` holds the quotient and `field-*-polynomial.json` holds numerator and denominator. I traced `aa5461ab…` and the constant field `79ea957c…`. The tanh argument (`x/10`) and the T-variable are kept separate, and each endpoint is joined to the saved `reconstruction-input` and its zero return.
- **Flat and height closed forms:** `F0=mu/(q+β)` and `Bh=-i·mu·q/(q+β)` match the census, with β=(30+9i)/109. The normal Holder diagonal for the minus face is `-mu q²/(q+β)`, which equals `i·f·q·Bh` with f=-1.
- **Native trace:** the minus-face height is `-eta·w1/2`, the saved map is `w1 → 2·hhat`, and the normal jet is `-i·q_o`. The trace/inverse check preserves the excluded η² term. The affine slot gives `F0(1 − product)`, equal to `F0 + η H Bh`, on both faces. Closed grazing values are only taken from closed expressions.
- **Control point:** p=1, q=2 and cs²=180/101 give a radicand of 4, with cs in (1,4)² scope.
- **Grade triples and counts:** 16 ordered triples (4 per coordinate squared). 200 cells (2×5×5×4). Response grades (0,1) and (1,1) are zero by Riemann–Lebesgue, with source/consumer cross grades kept in the end sums. The uniform bound avoids differentiating `exp(iQa)`, and the normal slot uses the Holder certificate directly.
- **Address coverage:** controls 1 (height-contact omission) and 2 (reverse phase) have 48 plus-face pressure height addresses to draw from. Example: address 8034 has constant nonzero source and consumer fields.

### Nonblocking limitations

- Memory is untested. `raw[alias]` holds all 422 JSON files under `RLIMIT_AS` 4 GiB, including the 17 MB U0 address file.
- The response cache key `(proof, face, slot)` omits response grade. Proof 0 appears only with (0,0) in the index, and the saved response coefficient is proof-level. I did not exhaustively check the other 16 proofs.
- The launcher docstring still says "reference-response instrument".
- Final run acceptance still depends on the completion audit. The gate and build record are pending.

### Inspected

- **Fully read:** `build-guide.md`, `method.md`, `worker.py`, `launcher.py`, `runtime-source/inert-helpers.py`, `runtime-source/supervisor.py`, `execution-authority.json`.
- **Fully read from `inputs.json`:** the head, the manifest tail, and every key and pin line for the helper, c2, supervisor, guard, hook and physical input.
- **Fully read from the saved inputs:** `saved-inputs/weak/field-79ea…-{operands,reconstruction-input,derivative-class}`, `field-aa54…-{operands,reconstruction-input,polynomial}`, and `saved-inputs/reference/minus-new-native-trace`, `minus-final-native-slot-routing`, `minus-native-height-constant`.
- **Regions or greps only:**
  - `saved-inputs/weak/minus-global-height-PV-certificate.json` and `saved-inputs/weak/all-coefficient-certificates.json` (field `aa54…` entry).
  - `saved-inputs/saved/reference/retained-response-census.json` (grepped text lines).
  - `saved-inputs/full/all-local-cells.json` (first cell, lines 1-330).
  - `saved-inputs/full/new-pressure-wave-arguments.json` (grepped multiplier entries).
  - `address-review-index.jsonl` (counts via grep, lines 1-3, 2068, 8035, 10063).
  - `runtime-source/shared-guard.py` (grepped arguments and properties only).
  - `tests.py` (test names only).

### Not inspected

- Hash equality of the packet copies against the manifest pins. No shell was available.
- The remaining ~420 saved-input files and the full address arrays.
- `address-representatives.json`, `native/c1.py`, `native/c2.py` and `source/*` beyond what the worker's AST joins cite.
- The guard's main body.
- Whether the guard and supervisor are unchanged from the earlier build.
- Any peer report.

NEEDS REVISION