**NEEDS REVISION**

I found four correctable gaps in the worker and launcher. The method arithmetic and constants check out against the worker. What I could read is source only. The data behind the worker's joins is not in my directory, so those checks are unresolved.

## Coverage limit
- **Files read:** `method.md`, `build-guide.md`, `worker.py`, `launcher.py`, `tests.py`, `inputs.json` (all 1412 lines), the first 816 of 4573 lines of `packet-index.json`, and `execution-authority.json`.
- **Not read:** `prior-direct-method.md`, `prior-reference-method.md` and `review-prompt.md`. I only grepped the two prior methods for flat/depth-symbol conventions and found none.
- **Data not in my directory:** the 34 field objects, the 320-grade census, the 608 representatives, the factor proofs, the five address arrays and the helper source (`raw_increment`/`reference_grazing`). I did not read any actual address, field or helper.
- **Consequence:** every data-dependent claim below is unverified. That covers field form, jet names, constancy of source and consumer fields, flat-symbol naming, and what `Journal.zero`, `Journal.sinh_zero` and `containment()` do.

## What checks out
- **Constants:**
  - Re β_min = 3000/11101 at δ=1/10.
  - a*² = 879/400.
  - |q| ≤ |p|+4, via the gap (K+4)² − K² − 453/50 ≥ 0.
  - The weight integrals: ∫w = 10, sup w ≤ 2, per-root bound 18, Cq = 36/√a*.
  - The H bound: 5/(8π) + 25/(4π) + tail, far below 100.
  - Coefficients of the J, GD, Bc and flat/height/slope/C forms.
- **KJ and KD:** KJ = (4/5)·121·Cq/b³ and KD = 18·121·2Cq/b² follow from those inequalities. The numerator-gap polynomials have nonnegative coefficients.
- **Height-PV bounds:** the pressure and normal coefficients, and their root-difference identities, match method §5. q(l)Y is never assumed C¹, and q(l-t) is kept separate from q(k+t).
- **Controls:** the point k=3/2, l=2, t=1/2, cs²=10/7 gives qi=2, qo=3/2, qh=3/2, qs=2. The reflected-root control is therefore applicable with distinct qh and qs roots, and the H-contact point is consistent.
- **Startup order:** inert imports, then containment, then sympy, with a hook-first launcher, no deadlines and strict errors.
- **Failure handling:** exclusive-create writes, failure preservation, posthashes.
- **Manifest aliases:** every alias the worker reads is present in `savedInputs`.

## Blockers
1. **Whole-tag, face and grade joins are not implemented** (`worker.py:141-150, 183-203`).
   - `inventory/whole-tag-definitions.json` and `inventory/typed-direct-objects.json` are copied but never consumed.
   - `G` (line 26) is unused. The grade coverage is only checked as `len==320`.
   - Direct whole is checked only by `startswith('Dwhole_<face>(')` (line 197). Nothing counts that exactly one Dwhole or Hwhole appears per address, or that no resolvent or Rprod is hidden.
   - Method §6 requires whole-tag joins on both faces, independent η/σ order and "whole direct once, mixed once". The instrument would not produce that evidence, and success could not claim it.
   - **Fix:** consume both files and join the whole-tag and typed-direct signatures per face. Assert that coverage is exactly rows×slots×G×G, and that each address's (consumerGrade, sourceGrade) matches it. Assert a tag count of 1 on direct addresses, at most 1 Hwhole on mixed ones, and none elsewhere.
2. **Reused normal-q-l-to-r controls have no face or row mapping check** (`worker.py:309-323`).
   - The address comes only from `op['context']['addressId']`, and success is `len(restored)==4`.
   - A wrong face or row, or a duplicate, would still pass.
   - **Fix:** require the name set to equal `{THETA_BALANCE, E_W_BALANCE} × {plus, minus}`. For each, require that the address's face and consumer row match the name, and that the sign and `normalJet[face]` match.
3. **The Leibniz control point is unguarded, and the selectors can abort the one authorized run.**
   - Leibniz point (`worker.py:336`): `flat_point` substitutes only `qo`. Flat support makes the input depth q(l), but the code does not check that the lifted flat expression contains only `qo`. If it is written in `qi`, the movement stays symbolic and `exact_nonzero` fails or the control is wrong.
   - Selectors: they filter on grade or jet name only (lines 329, 342, 360). The constant-source and constant-consumer requirements come afterwards as `require` (lines 337, 346, 361). The first matching `min(addressId)` may therefore be nonconstant and abort after the single authorized run.
   - **Fix:** put constancy into the selector predicates. Substitute the full flat map and `require(not flat_point.free_symbols)` in the operands before the guard.
4. **Authority and containment are not bound** (tooling, but required by the identity checks).
   - `verify_gate` (lines 94-96) checks only two booleans. The manifest pins the old source-composition authority (`inputs.json:1386`), not the new `execution-authority.json`. Nothing compares `authority['scope']` to `manifest['scope']`, or checks `scienceExecutionsAuthorized==1` and `noDeadline`.
   - Line 386 saves `containment()` without asserting anything. I cannot see whether the helper raises on a violation.
   - **Fix:** pin the new authority and require those fields. Assert the 4 GiB limit, zero swap and single-thread settings in the containment record before importing sympy.

## Minor, optional
- `jets[name]=jet` (line 198) can silently merge two specs that share a name. Use `require(jets.setdefault(name, jet) == jet)`.
- `TEXT_HELPERS`, `itertools` and `SimpleNamespace` are dead.
- The review record's per-engine reports are not hash-bound.

## Runtime obligations after revision
Check every address at runtime, since I read none:
- All 34 fields reconstruct as polynomials in T with constant denominators and no leftover x.
- Real ω=3 and the actual L=10 hold.
- Factor and normal-role joins hold.
- The four prior controls restore with their actual arguments.
- The new controls persist operands and applicability before guards.
- Strict stop on any nonzero join, with posthashes and copies.

## What success would establish
- **Would establish:** a source-joined certificate that the saved coefficients, bounds and routes agree with the method's inequalities. That includes global polynomial envelopes with constants, a bounded-derivative recurrence for all 34 fields, the Schwartz weak-pairing bookkeeping, and three formal controls.
- **Would not establish:** a computed integral, a loss, or a finite solve. It also would not give a measure-theory machine proof, plane-wave or scattering results, or runtime restoration of the source inventory. The analytic arguments (induction, uniform integrability, dominated convergence) remain assessed mathematics.