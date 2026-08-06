# Stop enumerating pairs — pin the injectivity itself

**Do not commit.** One file: `reduction/test_engine_output_checks.py`.
⛔ **Do not change `reduction/engine_output_checks.py`.** The implementation is correct; this round fixes
only what the tests pin.

---

## ⭐⭐⭐ WHY THIS ROUND EXISTS — read before editing

`test_canonical_derivative_delimiters_do_not_alias_identity` has been strengthened **three times** by adding
hand-picked pairs, and a review leg has found a new escape **every time**:

| round | pairs added | weaker encodings a leg then found still passing |
|---|---|---|
| 2 | comma-in-argument | escape `,` only · escape all but `\` · drop the empty marker |
| 3 | those three | **order field raw** · no escape of `^` · escape only `\,^` · encode only `function_arguments` · minimal (`\,` + `\0`, order raw) |

⛔⛔ **A fix that bans one more spelling per round is the wrong shape.** ⭐ The next leg will find the next
escape, and the round after that another. ⚠ Note that round 3 added the `order` encoding **and its own test
does not fail an implementation that leaves `order` raw** — the pin did not even cover that round's own
change.

⇒ ⭐⭐ **Replace the enumeration with a test that pins the PROPERTY.**

## R1 · A test that fails for ANY non-injective encoding

⭐ **Name the property, ⛔ do not hand me another list of pairs.** The requirement:

> The map from `(function_name, function_arguments, differentiated)` to the rendered symbol name is
> **injective** on the declared domain — `str`, `tuple[str, ...]`, `tuple[tuple[str, int], ...]`.

⭐ **Write a test that fails for any encoding that is not.** Choose your own mechanism; two that work are a
verified left inverse (a decoder `D` with `D(N(t)) == t`) and an exhaustive sweep asserting distinct
identities never share a name. ⛔ Do not add a decoder that must be hand-maintained in lockstep with the
encoder unless a drift between them **fails the test** rather than passing silently.

⭐ **The generated domain must exercise every character the encoding treats specially** — at minimum
`\` `,` `^` `(` `)` `;` `[` `]` `0` plus an ordinary character — in **every** field, including the
`differentiated` variable names and orders, empty fields, and empty tuples.

⚠ **It is a unit test: keep it to a few seconds.** ⭐ Choose the alphabet and the length bounds to fit; a
few tens of thousands of tuples is ample. ⛔ Do not weaken the alphabet to buy speed — shorten the fields.

⭐ **Keep the existing hand-written pairs** as named regression cases for the specific historical defects.
⛔ They are no longer the pin, and ⛔ do not add more of them.

## R2 · Do the same for `OpaqueDerivative`'s equality, if it is not already property-shaped

`test_opaque_derivative_instances_keep_and_compare_their_full_identity` currently passes the weaker-variant
matrix both legs ran, so it may already be sufficient. ⭐ **Check, and change it only if a weaker
implementation passes it.** ⛔ If it is sufficient, say so and leave it alone.

---

## Acceptance — paste literal output

⭐⭐ **Build every weaker implementation below, run your new test against each, and paste the per-variant
result.** Each must **FAIL**. The applied implementation must pass.

```
no encoding at all
escape `,` only
escape `^` only
escape the separators but NOT the escape character `\`
escape everything but drop the empty-field marker
order field interpolated raw
escape only `\ , ^`  (not `( ) ; [ ]`)
encode only `function_arguments`
minimal: escape only `\` and `,` plus the empty marker, order raw
```

⚠ The last six all pass the **current** test. ⭐ If your new test fails on all nine, it is pinning the
property rather than a list.

Then:

1. The full unit suite — quote **which files you ran**; `test_engine_output_checks.py` alone is 75 tests,
   with `test_registry.py` it is 86.
2. `harness_ablation.py` — every `ACCEPTANCE` line.
3. The new test's runtime, and the number of tuples it generates.
4. `git status --short` and `git diff --stat`.

⛔ **Do not run the comparator and report an unchanged md5 as evidence.** This round touches no
implementation code, so an unchanged comparator is guaranteed and establishes nothing.

## Constraints

- ⛔ No `git add`, no `git commit`, no other git write.
- ⛔ Do not modify `engine_output_checks.py`. ⛔ Do not touch either engine or any committed `.out`.
- ⛔ Do not launch Mathematica or `wolframscript`.

## Report back — under 12 lines

1. `R1`, `R2`: done or not, with the net line change.
2. The nine-variant matrix and the acceptance output.
3. ⭐ Any weaker implementation you constructed that your new test still passes — ⛔ report it rather than
   adding a special case for it.
