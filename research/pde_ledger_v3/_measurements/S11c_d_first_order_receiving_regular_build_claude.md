NEEDS REVISION

I read `worker.py`, `build-notes.txt`, `method-plan.txt` and `guide.txt` in full, plus the `Journal.zero/nonzero` helpers and the saved chart and dual rows. I did not open most of the 2059 evidence files. The entries of the saved C3 block sit on very long lines I could not view. Items 1 and 2 below are therefore likely failures from the code's own logic, not confirmed against the data. Everything below is a source-level reading, not runtime evidence.

**Blockers**

1. **Entry denominators can be certified vacuously.** `worker.py:314` uses `sp.fraction(v)[1]` on the raw C entries. For an unjoined sum like `a/b + c/d`, `fraction` returns denominator 1. `domain_certificate` then passes trivially on `1`, and no real denominator is checked before cancellation. The guide requires every original denominator to be saved and cleared first. Fix this by using `together` (or collecting every `Pow(·,-n)` base), and persist the result for each entry.

2. **The denominator domain check is on the wrong object and the wrong criterion.**
   - It runs on C with `l` still present, not on `Cq` after `s=p²-q²`.
   - It accepts only degree-1 factors in `q` whose roots have both real and imaginary parts negative.
   - The saved chart and dual carry `K2=l²+1/20`, and the plan says so. That factor becomes `p²+1/20-q²`, with roots ±√6. Those roots are harmless: √6 > p, and the evanescent ray gives `p²+1/20+r² > 0`.
   - As written, such a factor either leaves `l` in the constant (so `is_zero is False` fails) or is an irreducible quadratic in `q` (`UNSUPPORTED_ORIGINAL_DENOMINATOR_FACTOR`). The run would stop on a legitimate denominator.
   - The same applies to `domain_certificate('determinant',D)` at `worker.py:333`.
   - Fix: apply the existing exact ray test (Bezout and Sturm, endpoints `0` and `p`, and `+∞`) to every entry denominator of `Cq` and to the raw determinant denominator `rawden`. Keep the negative-quadrant test only as an optional shortcut.

3. **Source-side `l` is never mapped to the receiving `l`.** Only `oldq→q` is substituted. At `worker.py:154` the code treats `oldl` as distinct from `l`, but the pressure rows (`:205`, `:222`, `:242`, `:245`, `:252`) never map `oldl→l`.
   - If `piece['pressure01']`, `flatResponse` or the flat-factor arguments carry `oldl`, the pressure rows keep a foreign symbol.
   - The zero checks could then pass trivially or fail spuriously.
   - `psi=aff['mu00'].free_symbols-{l}` at `:221` would also zero any foreign `l` symbol as if it were a field.
   - Add a symbol-identity join: `oldl` and `l` must be the same symbol, or `xreplace` both on every transported operand.
   - Do the same for `B`, `Bp` and `JL`, and `require` that `B.free_symbols ⊆ {l}`. Otherwise `B.subs(l,-p)` at `:295` can silently be a no-op.

**Substantive corrections (not blocking alone)**

4. **Right-end force binding.** `:181-182` checks `R3(p)*localRight == 0` but never requires `localRight == 9*U`, as the plan asks. A zero product against an unverified operand does not bind the saved right-end force. Add `zero('right=9U', right, 9*U)`.

5. **Moving-projector check is incomplete.** `:279-281` tests only the `δ` coefficient `R3(p)B' + R3'(p)U`.
   - The `δ'` coefficient `R3(p)U` and the `U·T1` term are not directly checked at `l=p`. Line 195 checks `R3*U` only symbolically, and `:276` covers only column 1.
   - Add `zero('R3(p)U', R3.subs(l,p)*U, 0)` with both columns, and include `R3(p)U·T1`.
   - The derivative-drop control compares only entry `[0,0]`. Declare it and confirm in the journal that `(R3(p)B')[0,0] ≠ 0`. Otherwise the control is vacuous.

6. **Current argument convention is unjoined.** `JL(km,kp)` is bilinear in a conjugated bra, but nothing ties the argument convention to the saved weights. `B(p)†JL(p,p)B(p)` and `B(-p)†JL(-p,-p)B(-p)` should be compared with the saved `G0` and `Gref` (a join, not a replay). Without that, `(-p,+p)` for `C_-+` rests on an assumed convention, and the wrong-bra control does not show the convention is right.

7. **Growth bookkeeping.** `sufficientReceivingPolynomialWeight` counts only the multiplier and the E3 chart power. The forcing side also grows polynomially, through the `i·q` normal factor, the flat response and `R3`. The Schwartz transforms absorb this, but the stated weight should list it explicitly or be described as a lower bound.

8. **Hard-coded values.** `cs²=3/2` is typed into `:154-155` instead of being read from the bound values. `coefficient=Lambda_A_0/(rho*(1-iωτ))` is used where the notes say `10*A_mem/rho`. Join the two forms, or state their equality from the saved operands.

**Checked and acceptable on reading**

- **Row selection and chart:** the full 5×2 eta projection against the hand table (`:178-180`), `R3`/`E3` against the literal chart, and `D5` rows 2–4 as the lower three. All five row labels are used.
- **Depth law:** `q²=6-|h|²-l²` with `p²=119/20`.
- **Both faces and normal signs:** the pressure rows join both faces and slots, and the affine transport check at `:250` fixes the ±`iq` sign.
- **Q factor:** applied to the whole step before any cusp.
- **Endpoint-zero profile test:** `localized()` divides by `1-T²`.
- **Gaussian-rational Euclid/Bezout/Sturm:** the Sturm construction, the `+∞` rule, the exact endpoint checks, and refusal on unknown signs.
- **Threshold joins:** both joins are before the adjugate.
- **Controls:** the three controls act on full operands and use `J.nonzero`, which requires an exactly nonzero, finite movement.
- **Retained scope:** non-decaying transverse components are retained in `notClaimed`.

The tooling tests and the hand table are source and metadata previews only. They are not runtime results. A pass would still leave the nonuniform power balance to be done.