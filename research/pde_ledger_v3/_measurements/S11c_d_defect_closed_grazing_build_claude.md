**NEEDS REVISION**

I found no mathematics error in the stated constants and no definite crash in the worker. The revision is needed because two of the three required controls don't test the new expressions, the per-face coefficient join is partly tautological, and the row-side map to the saved Fourier conventions isn't joined to the saved objects. I read files only and ran nothing.

## Substantive validation blockers

1. **The wrong-quadrant control is tautological** (`worker.py:300-302`).
   - It asserts `sqrt(kap²-t²)>0` and `-sqrt(...)` is not nonnegative, using hand-set `test_qs` and `test_qo`.
   - No new expression is evaluated with the flipped sign. That covers `Bc`, `limit_i`, the `qs+qo` bound, the quadrant identity and the transported sheet branch.
   - The values are not tied to the transported `Piecewise`.
   - Consequence: the control cannot expose a wrong sheet. Substitute `-qs` or `-qo` into `Bc` or `limit_i` and require a nonzero residual against `G`. Also show the flipped value fails the first-quadrant test taken from the transported sheet.

2. **The lower-jet control is tautological** (`worker.py:303-307`).
   - `wrongJet` is defined as `-actualJet`, and `actualJet` hard-codes `-I*qo*profile*limit`.
   - `normalized==1` therefore holds by construction.
   - Neither jet comes from the saved lower object (`minus-trace-domain.normalJet`, `minus-closed-before-cancel.jet`), and no `limitPoint≠0` assertion is made directly.
   - Consequence: the control would pass without any lower-face sign being correct. Build both signs from the saved plus and minus jets at the physical point, and require the exact difference to be nonzero.

3. **The face coefficient join is partly empty** (`worker.py:129-131`).
   - For the plus face, `cf=C`, so `transported == Cu` and the check compares an expression with itself.
   - The minus face reads the unnamed positional `lower-boundary-return.json[3]`.
   - Neither uses its own face's saved `plus-closed-before-cancel.json` / `minus-closed-before-cancel.json` field `mixed`. That field is `C` at ω=3, ρ=1/10 and is already in the manifest.
   - The closed-density joins at `worker.py:162-165` partly compensate, but "both faces use their actual saved coefficient" has no content for the plus face. Add `J.zero(before['mixed'].xreplace(depthmap), Cu.subs(omega,3))` per face, with a named lookup instead of `[3]`.

4. **The row-side joins and normalization are incomplete.** This is a blocker unless the claim is narrowed explicitly.
   - The D-slot is assumed to be the open raw kernel `pref·Bu`. The worker never loads `rawKernelPlus`/`rawKernelMinus`, `rawRowDensity`, `physicalRowDensity` or `before['raw']`. A single `sinh_zero` of `before['raw']` against `pref·Bu` at ω=3 would close this.
   - The Fourier join (`worker.py:252`) evaluates `source00` only at `k→kin[0]`. In `fourier-convention.json`, `kin[0]` is 0, so the join does not test the k-direction (d1) factor.
   - The u_1 channel has zero coefficient at that point, so it is also untested.
   - The multiplier `mult` and the saved `rowFactorsPerCombinedSource` / `rowAmplitudesPerWholeKernel` are never compared. The 1/60 normalization of the denominator and the row amplitudes are not tied to saved values.
   - `sourceMinus` is identical to `sourcePlus`, and the saved convention records `lowerFaceCorrectionConstructed:false`. Nothing here supports lower-face source phase or sign.
   - Without these joins, the strongest supported statement is "finite polynomial in source jets under the saved jet-name→Fourier map".

## Mathematics I checked and found sound

- **C=tB:** it follows exactly from `qh²=qi²-t(t+2k)` and `qs²=qo²+t(2l-t)`. I checked this by hand against the saved `C` and `B`. It needs no ω-specific relation, so the unrestricted-ω transport is valid.
- **Contact at t=0:** the numerator vanishes under `qh→qi`, `qs→qo`, and that identification is legitimate because `qh` is the same function as `qi` at t=0.
- **Bc and the saved reference:** `Bc = B·qi·qo/den` is correct. The saved `reference` (60+91i form) equals `qiqo/((qi+β)(qo+β))` with β=(30+9i)/109.
- **Prefactor:** the saved prefactor reproduces `WL/(4i)·A(t)·A(Q-t)`, with `5A` for the jet.
- **β constants:**
  - `D_δ`, `Re β = 0.3/D`, `Im β` and β_min=3000/11101 all check out.
  - The identity `Dmax-D=(dmax-δ)τ(2+τ(dmax+δ))` is right.
  - The identity `kn-km²=(9-δ²)(1/cs²-1/4)+(dmax²-δ²)/4` is right.
- **Envelope and tail:**
  - The envelope coefficients 18+3|t| and 43+3|t| are right, and the endpoint enclosure [-6,6] is right.
  - The tail steps all hold. They are `|Q-t|∈[|t|/2,3|t|/2]` and `|A|≤L|x|e^{-5π|x|}`. Also 61+6|t|≤12|t| and the constant `C`.
  - The antiderivative of `t³e^{-10πt}` is correct.
- **Local bound:** `2√(2m)` per endpoint is right.
- **Routes:**
  - `qs` is correctly kept distinct from `qh`.
  - For l=-k the radicands coincide (`worker.py:226`), and the envelope remains a sum, so no product of singularities arises.
  - The same-momentum and opposite-momentum endpoint lists are right.
- **Limits:** the input-only, output-only and simultaneous rational limits are right and are taken away from internal endpoints.
- **Scope:** only the saved sheet and row inputs stay at real ω=3, and only the kernel is continued.

## Concrete runtime and tooling items, none a definite crash

- **Memory and time:** `sp.cancel(sp.together())` runs on the 1.1 MB THETA/E_W increments (`worker.py:266` and `285`) with no deadline under 4 GiB. A cgroup OOM kill cannot be caught, so `failure.json` and `posthashes.json` would be lost. The incremental artifact and operation indexes would survive.
- **Piecewise handling:** the sheet code assumes `transported.args` keeps three arms in order and that `(nan, True)` is preserved (`worker.py:175-185`). A SymPy re-fold would give a preserved but spurious stop.
- **Persistence order:** in `nonnegative_polynomial`, the `require` at line 70 runs before the `-certificate` emit at line 74. The input is persisted first, so this is minor.
- **Gate:**
  - The build clearance is checked as author-writable boolean flags plus a worker/manifest/guard/supervisor hash join. The reviewers' report content is not pinned.
  - The worker does not tie `g['sharedGuard']` or `g['supervisor']` to the pinned executed paths. The launcher command uses fixed paths that are pinned in `sourcePins`, so this is mitigated.
  - The hook script and the `codex` binary are not hash-pinned.
- **Tests:** `tooling-tests.py` doesn't cover launcher gate checks, hook-first ordering, or `savedFiles ⊆ sourcePins`. That is acceptable for stdlib-only tests.
- **What passes:** containment, no-deadline enforcement, the pooled-guard environment requirement, import of scientific libraries only after containment, and unchanged-helper AST copying all pass on inspection. The launcher keeps hook-first behavior. If the hook fails to arm, the pipe closes and the coordinator refuses to launch.
- **Namespace and bookkeeping:** I traced the helper namespace (`sp`/`Str` injected before `Journal`), the `global sp` ordering, duplicate artifact names (none found), source and copy posthashes, and the failure paths. All are sound.

## Optional wording

- The nonnegativity of `dmax-δ`, `2-cs` and `dmax²-δ²` is used but only the identities are certified.
- `method.md` says "tanh factors" where the profile is sinh.
- State that the principal root is continuous on the closed upper half-plane, so a negative-real radicand gives `+i√|·|`. That is what makes the pointwise limit match the saved sheet.
- State that `|A|≤1/(2π)` on compact sets.
- The continuation of the closure leg `ρω/q` to Ω is supported only numerically at ω=3 (the saved matrix is numeric). The `ω ρ` factor in the saved lower pressure operand supports it, but it is not asserted symbolically.

## Strongest claim the future output could support

For both faces, with the stated scope, the output could support the following:
- the exact certificates, namely C=tB and the closed-density, native-law, root and sign joins;
- the β, domain and tail constants;
- the grazing limits and the L1 limit of the closed kernel;
- the per-face finite-polynomial multipliers.

Together these support an L1 limit of the analytically continued kernel, under the independently assessed uniform-integrability argument. The identification of that limit with the real-frequency grazing object, and the contact zero at every δ>0, also rest on that same assessed argument. The output would not support correctness of the controls (until blockers 1 and 2 are fixed) or the lower-face source map. It would not support an integral, a finite solve, iteration, the full operator, loss, or κ=0 or β=0 behavior.
