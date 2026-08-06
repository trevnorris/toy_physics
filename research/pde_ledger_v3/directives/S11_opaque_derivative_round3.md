# Two items: pin the canonical test, and deliver R2's second invariant

**Do not commit.** Files: `reduction/engine_output_checks.py`, `reduction/test_engine_output_checks.py`.

⭐ **Smallest change that establishes each invariant.** ⛔ Add nothing beyond the four items below.

⚠ Two independent review legs examined round 2. ⭐ **Both confirmed the encoding itself is sound** — one
built a total decoder `D` and verified `D(N(t)) == t` on 1 280 329 tuples with zero failures, which is a
proof of injectivity, not a failed attack. ⛔ **Do not change the encoding.**

---

## R1 · ⛔ BLOCKING — the canonical test admits three weaker implementations

`test_engine_output_checks.py:286-293`. `test_canonical_derivative_delimiters_do_not_alias_identity` uses
one pair (`("x","y")` vs `("x,y",)`), so it pins only *"a comma inside an argument"*.

**Measured by both legs — reproduce before fixing.** Each of these passes that test **and** the full suite,
while leaving the false-AGREE channel open:

| weaker implementation | measured collisions |
|---|---|
| escape `,` only | 45 632 |
| escape every separator **but not the escape character `\`** | 45 632 |
| escape everything but **drop the empty-field marker** | 5 704 |

⚠ The second is the classic escaping leak, and the surviving collision is
`("f", ("\\","x"), (("x",1),))` against `("f", (",x",), (("x",1),))` — identical rendered name, `==` True,
first atom's `function_arguments` overwritten.

⭐ **Fix: add a case per weaker implementation**, so each is individually excluded — at minimum a pair that
separates only under escaping of `\`, and a pair that separates only under the empty-field marker.

⛔⛔ **AND THIS IS THE PART THAT KEEPS BEING MISSED.** ⭐ **Build each weaker implementation yourself, run
the test against it, and paste the per-variant result.** ⛔ A test that fails only on the *unfixed*
baseline is not pinning the invariant — the unfixed baseline is the easiest case to catch, and this is the
**third time in this workstream** that a half-fix passed the test, the whole suite, the full ablation
battery and produced byte-identical comparator output.

⭐ For reference the sibling test is already correct: `test_opaque_derivative_...` fails on all five
`OpaqueDerivative` weakenings. Match that standard.

## R2 · ⛔ R2 FROM THE PREVIOUS ROUND ASKED FOR TWO INVARIANTS AND DELIVERED ONE

The previous directive required `CanonicalDerivative` construction to be **non-mutating** as well as
injectively named. ⚠ Measured: it still uses the **cached** `sp.Symbol.__new__` and still rebinds attributes
on the returned object — `second is first` → `True`, with `function_arguments` rebound.

It is harmless **only** because the name is now injective, so equal names carry equal values. ⛔ That makes
the injective name a single point of failure with no second line of defence.

⭐ Give `CanonicalDerivative` the same two invariants `OpaqueDerivative` now has: construction never mutates
a previously constructed atom, and equality and hashing key on the **full identity**. ⭐ Use the same
mechanism, so the two classes do not differ for no reason.

## R3 · Remove the unenforced condition on `order`

`engine_output_checks.py:~202` interpolates `order` **raw** (`f"{encode(variable)}^{order}"`), so
injectivity is conditional on `order` being an `int`, and nothing enforces it. ⚠ Measured when violated:
`(("x",1),)` and `(("x","1"),)` render identically and `symbolic_equal` returns `True`; and `(("x",1),)`
against `(("x",1.0),)` are equal tuples rendering **differently**, a false disagree.

⚠ Both in-module construction sites coerce with `int(...)`, so it is unreachable today. ⭐ Encode the order
as well, so injectivity stops being conditional. One-line change; ⛔ do not add a type check instead.

## R4 · The import guard tests the name, not the invariant

`engine_output_checks.py:~24-25` checks that `sp.Symbol.__xnew__` **exists**. ⚠ A SymPy that keeps the name
but makes it cached passes silently, which is the exact condition the guard exists to prevent.

⭐ Make the guard test the **behaviour** — that two constructions with the same name yield distinct
objects — rather than the attribute's presence. ⛔ Keep it at import time and keep the failure loud.

---

## ⛔ Do not touch

The encoding scheme itself · the applied `OpaqueDerivative` fix and its test · either engine · any
committed `.out` · `checks_S9.yaml`, `checks_S10.yaml` · `harness_ablation.py`.

## Acceptance — paste literal output

1. `R1`: the three weaker implementations built and run against your **new** canonical test — each must
   **fail**, and the applied fix must pass. Per-variant result.
2. `R2`: a probe showing construction no longer mutates a previously constructed `CanonicalDerivative`,
   with the same probe against the pre-change module as the control.
3. `R3`: the two `order`-type collisions above, before and after.
4. `R4`: the guard firing against a stub whose `__xnew__` is cached.
5. The full unit suite.
6. `harness_ablation.py` — every `ACCEPTANCE` line.
7. Both comparators, stashed versus applied, complete stdout md5 on each.
   ⚠ ⛔ **Do not report an unchanged md5 as evidence your change is correct.** Both legs measured that no
   separator character occurs in any identity field on either config, so the encoder is the identity map on
   every real field and a knowingly-wrong variant produces the same bytes. ⭐ Report it as a regression
   check and say what it cannot establish.
8. `git status --short` and `git diff --stat`.

## Constraints

- ⛔ No `git add`, no `git commit`, no other git write except a stash you restore.
- ⛔ Do not launch Mathematica or `wolframscript`.
- ⭐ A script may run long **if it is visibly progressing**; ⛔ a long silent stretch is the failure.

## Report back — under 15 lines

1. `R1`–`R4`: done or not, with the net line change.
2. The acceptance output.
3. ⭐ Where your measurements disagreed with this directive.
4. ⛔ Do not list further defects you were not asked to fix; ⭐ say only how many you saw.
