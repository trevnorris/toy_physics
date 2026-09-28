# Literal Claude two-asymptote implementation review

Source: `/var/projects/toy_physics/_scratch/s11c/s11c-d-two-asymptote-20260927/build-review/claude.json`, field `result`; SHA-256 `d8dc6a26d5579c30d885d0f45f41dcc5865f6f693e8f7b814e35ae5baee5320e`.

**Verdict: NEEDS REVISION**

This verdict covers one limited claim: that the worker can reach `END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED` with exact first‑grade end lifts at blocks 16/17 (LEFT/RIGHT, grades 10/01), a correct formal native forcing/commutator reconstruction, and a meaningful responsive control. It does not cover any outgoing action, response, forcing domain, Green/FORM, A11/A12 or radiating claim.

The end‑lift mathematics (question 1) is sound. The revision is needed because the reconstruction guard, as written, will almost certainly fail on a correct construction, and because the control and the plane‑jet labelling don't do what the documents say. I read only packet files and executed nothing, so every "will fail" below comes from reading the source.

## Substantive findings

**S1. The reconstruction guard is wrong and will fail on a correct construction** (`S11c_d_two_asymptote_lift.py:578-586, 623-627, 642-648`)
- `lifted` substitutes the combined field `χ_L f_L + χ_R f_R` into each original integral. `linear_form` then runs `expand` inside that single integral.
- SymPy's `Add` automatically merges terms from both ends that share a monomial. For example, `C1_L/2·K·e^{ik0z'}` and `C1_R/2·K·e^{ik0z'}` become one integral with coefficient `(C1_L+C1_R)/2`.
- The side actions keep them as two separate integrals, `Integral(C1_L/2·K·e)` and `Integral(C1_R/2·K·e)`.
- These never cancel structurally, so `require(not residual.has(Integral))` at line 647 fails.
- **Smallest fix:** in `linear_form`, split each term with `term.as_independent(*bound_vars)` and emit `coeff*Integral(dependent, *limits)`. Then collect the residual by distinct Integral carrier and send each scalar coefficient to `zero_check`. This is the formal carrier‑algebra identity the documents describe.

**S2. The nonlocal mutation control is empty** (`:658-678`)
- `direct − decomposed` is zero by construction, so `actualResidual = −mutation` just restates the mutation.
- The operand used (`side − χ·plain`) is not an addend of `decomposed`, which uses `χ·reference_plane`. So nothing that is actually in the decomposition gets omitted.
- **Fix:** omit one real nonlocal addend of `decomposed` (a `side_actions[e]` term value). Pass the mutated residual through the same corrected reducer from S1, and require it to show a surviving carrier with an exact nonzero coefficient.
- Keep `integratedResponseNonzeroEstablished: False`. This shows that the formal check responds; it does not show a nonzero integrated response. It is the in‑scope test the contract asks for ("calculate the resulting residual").

**S3. Raw native plane actions are saved after the guard** (`:651-657` vs `:642-648`)
- If the guard fails (as it will today), the raw native plane‑jet operands are never persisted. That contradicts implementation.md:73-74.
- **Fix:** compute and save `raw_plane` before the reconstruction op.

**S4. The plane‑jet reduction is recorded as established when it is a premise** (`:629-632, :737`)
- The decomposition's "commutator" is `L0(χE) − χ·(P0 f − iP0′f′)e^{ik0z}`. That equals the true `[L0,χ]E` only if the native L0 acting on `z·e^{ikz}` equals `−i∂_k` of the accepted harmonic action.
- The accepted identity only covers harmonics. The extension also needs:
  - (a) L0 at zero grade is translation‑invariant and free of the regulator and profiles;
  - (b) `∂_k` commutes with the native ordered integrals and any regulator limit near k0.
- **Fix:**
  - Replace `planeJetActionUsesAcceptedFourierSourceIdentity: True` with an explicit premise record marked unverified.
  - Add the cheap check for (a): require/record `not strong_zero.has(alpha)` and no profile `AppliedUndef` in it.
- Premise (b) stays an open obligation, and that is consistent with the stop.

**S5. The "grade‑free integral" premise may not hold** (`:346-349`; contract line 36; `S11c_d_continuum_grades.py:134-168`)
- The saved grade pipeline keeps separate `factor`/`source` grade records and multiplies cell × factor × source grades together. That suggests grade dependence may sit inside the integrals.
- If it does, the guard stops the run and there is no path to `c00·∂_h I`.
- **Fix:** before gating, cite the saved grade inventory showing every `factor`/`source` record is grade (0,0,0) only. If that can't be shown, add the path that differentiates the intact original integrand in `h`, with limits unchanged and required to be free of generators.
- This is a premise to confirm or a path to add, not an error in the existing method.

## Question 1: the local coefficients check out
- **Laurent algebra:** correct for d = 0, 1, 2. With `c_{-2} = N_0/D_2` forced to zero, `B = (N_d − D_{d+1}A)/D_d`, and degree 3 is enough.
- **Branch:** the radical series recursion, positive `r0`, `r0² = radicand(k0)` and the saved residue `n·adj^{(n−1)}/det^{(n)}` all join correctly.
- **C2, C1 and the field:** `C2`, `C1` and `e^{ik0z}(C1 + izC2)` are the correct `(k−k0)^{-1}` coefficient for the full rank‑two frame.
- **Factorials and units:** the `1/j!` factor is trivial at j ≤ 1. The unit typing works out, and no channel normalization is applied.
- **Limit of the end check:** the check `P f − iP′f′ + P_h A` is a real test, but it reduces to `PA = 0` and `PB + P′A = I`. It cannot see components in the range of A, so `A P_h B` and `A P_h′ A` go untested. Those only move homogeneous plane waves between T and V.
- **Optional:** add the right‑sided identity `C2P = 0`, `C1P + C2P′ + A P_h = 0`. It is cheap and covers `A P_h B`. `A P_h′ A` rests on the derivative join at `:438`.

## Question 3: the stop itself is justified
Once S1–S4 are fixed, `END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED` is a useful and honest stop. The exact B, C1, C2 and per‑term native operands are the real deliverable. No global theorem is needed.

## Question 4: instrument
Persistence before guards, failure posthashes, `exist_ok=False`/`'x'`‑mode writes (no retry), gate hashes and the containment checks are all sound, apart from S3.

## Optional observations (no new review cycle needed)
- **Unit checker:** `UnitAudit` (`:240-294`) lacks the source's `[wm]1_profile` zero‑unit rule (source-definitions.md:16-17). An unknown Symbol raises an error, so this could stop the run spuriously.
- **Subs handling:** `canonical_origin` rewrites `Subs` with `.subs` (`:384-385`), while the source uses `Subs(*args).doit()` (source-definitions.md:314-315). The structural join at `:397` could then fail spuriously.
- **Runtime of local tail limits** (`:695-709`): 200 `Limit.doit()` calls on a local part that doesn't decay on its own, run before the later tags. This risks the 840 s budget. Drop them, or move them after all core saves.
- **Radical series cost:** reducing radical powers mod `r² = radicand` before the series would cut the `simplify` cost.
- **Codec:** `SavedCodec` allows everything in `builtins` (including `eval`). Hashes are pinned so the risk is low, but an allowlist would be safer.
- **Unused inputs:** `validationCheckpoint`, `candidateOperationIndex` and `domain` are pinned but never joined. Add `cp['productionCheckpoint']['sha256'] == actual[...]` or drop them.
- **Artifact index:** artifacts written with `j.value` are not listed in any hash index.
