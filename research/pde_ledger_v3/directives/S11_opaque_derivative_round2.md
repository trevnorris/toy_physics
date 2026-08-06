# Two real defects the review found, and three small corrections

**Do not commit.** Files: `reduction/engine_output_checks.py`, `reduction/test_engine_output_checks.py`,
`reduction/measurements/opaque_derivative_identity.py`.

⭐ **Prefer the smallest change that establishes the invariant.** ⛔ Do not add checks, guards or tags
beyond the five items below, and ⛔ do not improve anything you were not asked to touch.

⚠ Two independent review legs examined the previous round. ⭐ **Where this directive states a
measurement, re-run it yourself before relying on it** — the two legs contradicted each other on `R2`
and the orchestrator adjudicated by construction.

---

## R1 · ⛔ BLOCKING — the new unit test does not pin what it claims

`test_engine_output_checks.py:261-280`. The test's `conflicting` atom differs from `original` in
`orders` **and** `function_name` **and** `variables` **simultaneously**. ⇒ any implementation keying
equality on **any one** of the three passes it.

**Measured by a leg, and reproduce it before you fix it:** a variant whose `_hashable_content` adds only
`self.orders` passes this test, passes the **whole 82-test suite**, and produces **byte-identical**
`harness_ablation.py` and S10 comparator output — while a `_map_tree` function-name rename and a
variable rename both still leave `mapped == original` and `hash(mapped) == hash(original)`. Three such
weaker variants pass.

⭐ **Fix: assert on pairs that differ in EXACTLY ONE field each** — one pair differing only in `orders`,
one only in `function_name`, one only in `variables` — so each field is pinned separately.

⛔ **And drop the `hash(conflicting) != hash(original)` assertion, or restate it.** Hash inequality is
**not** implied by the invariant: unequal objects may legitimately collide, so as written a correct
implementation can fail this line. ⭐ Equality is the invariant; assert that. If you keep a hash
assertion at all, it may only assert that **equal** objects hash **equally**.

## R2 · ⛔ BLOCKING — `CanonicalDerivative` has the same defect, and it is the atom that decides verdicts

`engine_output_checks.py:174-191`. Its `sp.Symbol` name is built by joining the identity fields with
`,` and `^`, and **those separators are not escaped**, so the name is **not injective** on
`(function_name, function_arguments, differentiated)`.

⭐ **Adjudicated by the orchestrator, by construction — verify it yourself:**

```
a = CanonicalDerivative("f", ("x","y"), (("x",1),))
b = CanonicalDerivative("f", ("x,y",),  (("x",1),))
both render to  DerivativeIdentity[f(x,y);x^1]
a is b -> True ,  a == b -> True ,  a.function_arguments silently becomes ('x,y',)
symbolic_equal(a, b) -> True
```

The `^` separator collides the same way. A control with genuinely different identities and no separator
in any field stays distinct, so the probe is able to fail.

⚠⚠ **This matters more than the atom already fixed.** Comparison runs on canonicalised trees
(`~:1517`, `~:1543`), so **`CanonicalDerivative` is the actual comparison key** — a collision here makes
`symbolic_equal` return `True` for two genuinely different derivatives, which is a **false AGREE**.

⭐ **Fix: make the name injective on the full identity, and make construction non-mutating**, to the
same two invariants `OpaqueDerivative` now satisfies. ⭐ Choose your own encoding. ⛔ Do not change what
`str(atom)` is used for elsewhere without checking every consumer.

⚠ **Census, for scope:** 0 constructions on the S9 config, 2072 on S10 across 162 distinct names, and
**zero** names carrying more than one identity. ⇒ dormant on today's data. ⛔ That is a reason to fix it
cheaply, ⛔ not a reason to defer it: S11 emits Euler–Lagrange systems for seven packages.

## R3 · The `__xnew__` dependency is unpinned

`engine_output_checks.py:~167` calls `sp.Symbol.__xnew__`, an **undocumented SymPy internal**, and a leg
found **no SymPy pin anywhere in the repository**. A version that renames or removes it turns this into
an `AttributeError` on every derivative construction.

⭐ Make the failure loud and immediate rather than latent — a guard at import that fails with a clear
message if the internal is absent, **or** an implementation that does not depend on it. ⛔ Do not add a
dependency file or change packaging.

## R4 · The measurement script states a conclusion in its docstring

`reduction/measurements/opaque_derivative_identity.py` — its module docstring calls these
*"five process-global identity probes"*. ⛔ *"Process-global"* is a **conclusion about the defect**, and
post-fix it is factually wrong about what the script measures. ⭐ A script that must state no conclusion
must not state one in its docstring either. Reword to name what is printed.

## R5 · One directive claim was overstated — correct it where it is repeated

`directives/S11_opaque_derivative_fix.md:44-45` says the second channel *"needs no declared rename at
all"*. **Measured false for `OpaqueDerivative`:** `parse_derivative` builds `rendered` as an injective
serialisation (its identifiers cannot contain the separators), and `_map_tree` is the only other
construction site — a 162/162 name↔identity bijection on the real run. ⇒ within that class a
rendered-name collision can only be produced by a rename.

⭐ If any comment or docstring **you are touching** repeats that claim, correct it. ⛔ Do not edit the
committed directive itself, and ⛔ do not go looking for other places to correct.

---

## ⛔ Do not touch

The applied `OpaqueDerivative` fix itself, which both legs verified establishes its two invariants ·
either engine · any committed `.out` · `checks_S9.yaml` or `checks_S10.yaml` · `harness_ablation.py`.

## Acceptance — paste literal output

1. `R1`: the three weaker variants, built and run against your **new** test — each must now **fail**, and
   the applied fix must pass. Show the per-variant result.
2. `R2`: your own reproduction of the collision above, then the same probe after your fix, then the
   able-to-fail control.
3. `R2`: a census of `CanonicalDerivative` constructions and distinct names on **both** committed
   configs, before and after.
4. The full unit suite.
5. `harness_ablation.py` — every `ACCEPTANCE` line.
6. Both comparators, your change stashed versus applied, **complete stdout md5 on each**.
   ⚠ ⛔ **Do not report an unchanged md5 as evidence that your change is safe** — a leg measured that a
   knowingly-wrong variant also produces the identical md5, because zero rendered names carry more than
   one identity. ⭐ Report it as a regression check only.
7. `git status --short` and `git diff --stat`.

## Constraints

- ⛔ No `git add`, no `git commit`, no other git write except a stash you restore.
- ⛔ Do not launch Mathematica or `wolframscript`.
- ⭐ A script may run long **if it is visibly progressing**; ⛔ a long silent stretch is the failure.

## Report back — under 18 lines

1. `R1`–`R5`: done or not, with the net line change.
2. The acceptance output.
3. ⭐ Whether your own measurements agreed with this directive — **including where they did not**.
4. ⛔ Do not list further defects you were not asked to fix; ⭐ say only how many you saw.
