**NEEDS REVISION**

I found three blockers and four non-blocking concerns. This is a source-only assessment, not a runtime certificate.

## Blockers

**B1. The L_W scaling of the mapped profile field is never checked against the saved mapped field.**
- `library.py:96-97` computes `add(scale(unit(L_W),order),(-order,0,0))` from the derivative order and the registry alone. That is always zero, so `require(formal==ZERO)` can never fail.
- The saved `mapped` expression is only tested for being zero when the derivative is transverse (`library.py:95`).
- `record['L']==context['numeric']['L_W']` (`library.py:88`) only compares the stored L value.
- Nothing parses `mapped` to confirm it equals `L^r·∂ₓ^r` of the profile.
- `build.md:57-58` says "Original coefficient, every map pair, actual L_W, and saved mapped field must join", so the worker falls short of the stated join.
- For `profile-map-1325.json` (`w1_profile_d1`, L=10), the saved map `1/2 − tanh(x/10)²/2` is consistent with the rule. The worker would accept an unscaled map just as readily, which is the defect the `missing-profile-L` control is meant to catch.
- The unit conclusion "dimension zero" therefore rests on formula arithmetic, not on the saved operands.

**B2. The chemical and stage2 identity chain is joined on only one side.**
- Several `prior()` calls pass only `right=`: `worker.py:206` (raw-join), `:208` (epsilon-join), `:209` (velocity normalization) and `:218` (domain zero-grade).
- The left operand is unchecked in each.
- In the original source, the epsilon join is `bind(chemical['raw'][1]/epsilon) → chemical['amplitude']` (`S11c_d_defect_source_composition.py:117`).
- The raw-join left is the symbolic native leaf (it starts `B_rho_3*…`), while the epsilon-join left is the bound numeric `(…)/epsilon_shape`.
- No worker line ties the unit-checked symbolic raw leaf (`worker.py:202`) to the bound left of the epsilon join. At most the saved zeros imply that tie.
- The saved `*-native-chemical-amplitude-join` record (`inp['chemicalAmplitude']` against `amplitude`) and `*-native-density-raw-join` are never read by the worker. This weakens the claim that every cancellation identity keeps both operands.
- A symbol-set check on the bound left is cheap and needs no CAS: no raw native symbols should remain, and only the substituted names should appear.
- The velocity normalization left (`e_W_t/2`) happens to equal both the right side and the native velocity in the plus sample, but `worker.py:209` doesn't enforce that.

**B3. The control 1 baseline and mutant are hand-written arithmetic, not a native calculation.**
- `worker.py:305` computes `measured = L-unit·0 + (−r,0,0)` and compares it with a hard-coded `U.ZERO`.
- It never goes through the native unit engine, the saved profile map or the actual L-scaled expression.
- As written it is a tautology: any `r>0` refuses. It does not test the transport.
- It would only become responsive together with the B1 fix, by removing the L factor from the actual mapped expression and re-checking it.
- Control 2 is genuinely native. It uses the saved chemical constructor and registry, and checks for an "inhomogeneous" refusal.
- I did not confirm that `eta_bg` sits in an addition inside that constructor. The control only passes if it does, and the file I'd check was too large to read.

## Non-blocking concerns

1. **Homogeneity argument (`build.md:34-46`, `worker.py:193,246-248`).**
   - The reasoning is valid only under conditions the packet checks as flags or saved certificates rather than derives.
   - The conditions are: substitutions preserve units, the quotient is regular at eta=σ=0, and both numerator and denominator are homogeneous.
   - Not assigning a unit to the normalized numeric denominator is correct, and I agree with it.
   - The argument needs the saved denominator to be a unit-homogeneous unbound expression times a nonzero constant. That is the same as the unbound reciprocal base being homogeneous, which `reciprocal_bases` records. It also needs the quotient's Taylor coefficients to be taken in dimensionless variables, and `worker.py:230` checks this.
   - I rate this acceptable once B1 and B2 are fixed.
2. **Weak regularity test (`library.py:57`).** `any(v is True for v in row)` for `signedNonzero` accepts a row with one true entry and one false. The plus chemical-domain fraction does have this shape (`true, false`). The intended meaning is not documented. The fallback boolean state path at `library.py:53` is the exact "boolean as evidence" pattern you asked about. The plus chemical-domain case did not use it, since a fraction file exists, but I did not check the other five groups.
3. **Zero locations.** They keep only a required-unit expectation (`worker.py:269,281`), as the build states. This is not a blocker.
4. **Field ID.** `worker.py:292` treats the field ID as an extra check only and joins on `U.equal` of the actual field. That is correct.

## Coverage

**Read in full:**
- `build.md`
- `evidence-guide.md`
- `review-prompt.md`
- `worker.py`
- `library.py`

**Read in part:**
- `profile-map-1219`, `-1266` and `-1325`
- `plus-own-velocity-normalization-input`
- the head of the chemical raw-join and epsilon-join inputs
- the head of the chemical-domain fraction file
- the relevant lines of `source/S11c_d_defect_source_composition.py`

**Not examined:**
- `inputs.json` beyond grep counts
- the native unit walks and returns
- `tooling-tests.py`
- the contracts JSON
- the 544 selected addresses and 34 field records
- the consumer epsilon identities
- the other 150 or so profile maps
- the minus-face and consumer-group saved files

I make no claim about those. The three blockers above do not depend on them.

**Unknowns I'm preserving:**
- the gamma registry inference
- whether `launcher.py` and the guard behave as described
- whether the consumer slot unit returns (`worker.py:234-235`) truly match the full consumer coefficient operands

I found no ordering problem with persistence on refusal. `J.start` and `J.emit` write inputs and decision operands before each `require`, and the `finally` block writes `posthashes` and `journal-result`.

**To revise:**
- Parse each saved `mapped` value and check it against `L_W^r·dʳ/dxʳ` of the profile, or an equivalent exact operand identity.
- Join both sides of the chemical/stage2 cancellations, and use the unused amplitude-join and density-raw-join records.
- Make control 1 mutate the actual mapped expression and go through the native path.