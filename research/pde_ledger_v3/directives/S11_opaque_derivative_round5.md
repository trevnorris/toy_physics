# Name the test for what it establishes, and widen its envelope where that is free

**Do not commit.** One file: `reduction/test_engine_output_checks.py`.
⛔ **Do not change `reduction/engine_output_checks.py`.** The encoding is correct and proven.

⭐ **This is the LAST round on this material.** ⛔ Do not add checks beyond `R1`–`R3`.

---

## Why

Two review legs measured round 4 independently.

⭐ **What it establishes:** the property test kills **21** distinct weaker encodings between the two legs —
every one either leg could construct by weakening the **escaping logic**, including one that was the
entirely unescaped pre-change name plus a suffix carrying `len(function_arguments)`, which had passed every
earlier round and the whole suite.

⛔ **What it does not:** a leg then built nine encodings that **pass** it. They do not weaken the escaping;
they agree with the real encoder on the swept envelope and differ outside it — truncate every field to
length ≤ 1, casefold, drop arguments beyond arity 2, mangle the digit `"1"` (which never appears as field
text anywhere in the sweep).

⚠⚠ **And that gap cannot be closed by testing.** The declared domain is
`str × tuple[str, ...] × tuple[tuple[str, int], ...]` — **infinite**. For any finite test there is an
implementation agreeing on exactly the tested set. ⛔ **Do not attempt to close it.** ⭐ Chasing it is what
consumed four rounds.

⇒ ⭐⭐ **The defect is that the test's NAME claims the thing the test cannot establish**, and this project's
governing rule is that an instrument claim is stated on an invariant, never on a sampled outcome. What
establishes domain-wide injectivity is the **alphabet argument**, which two legs have now given
independently — one with a verified left inverse over 1 280 329 tuples.

---

## R1 · Rename the test, and scope its docstring

`test_canonical_derivative_rendered_name_is_injective_on_declared_domain` — the name asserts injectivity on
the declared domain. ⛔ It does not establish that.

⭐ Rename it for what it does: verifies a **left inverse on an enumerated finite envelope**. ⭐ Give it a
docstring that states, in this order:

1. the envelope it sweeps — alphabet, field-length bounds, arities, order values, tuple count;
2. that a green run does **not** establish injectivity on the declared domain, and cannot, because that
   domain is infinite;
3. that domain-wide injectivity rests on the **alphabet argument**, ⛔ not on this test;
4. what it **is** strong against: any encoding that weakens the escaping logic. ⭐ Name the count of
   independently-constructed weaker encodings it is known to kill.

## R2 · Widen the envelope where it is free

⭐ Add to the swept domain, keeping the test **under about ten seconds**:

- **multi-character fields** — length up to 2 at least, in every field slot;
- **argument and differentiated arity 3**;
- **plain decimal order values** — `1` through `12` — as well as the existing adversarial renderings.
  ⚠ The digit `"1"` currently never appears as field text anywhere in the sweep;
- a **space**, an **uppercase letter**, and a **digit** in the field alphabet.

⛔ **Do not report a target number of surviving variants, and do not iterate until some count reaches
zero.** ⭐ Widen once, then measure.

## R3 · Put the argument where the proof is

The alphabet argument is the thing that actually establishes injectivity, and it currently lives only in
review transcripts outside the repository. ⭐ Record it as a comment or docstring **beside the encoder** in
`engine_output_checks.py` — ⚠ this is the **one** edit permitted to that file, and it must be
**comment-only**: the escape set including the escape character itself, the empty-field marker being
outside the encoder's image, and why no unescaped separator can be produced.

⛔ Change no code in that file.

---

## Acceptance — paste literal output

1. `R1`: the new test name and its full docstring, quoted.
2. `R2`: the widened envelope — alphabet, length bounds, arities, order values, total tuple count, and the
   measured runtime.
3. ⭐ **Re-run both published weaker-encoding sets against the widened test and report which pass and which
   fail.** ⛔ Report the numbers as measured; ⛔ do not tune the envelope to improve them.
   - the nine from the round-4 directive: `none`, `comma`, `caret`, `separators_not_escape`, `drop_empty`,
     `order_raw`, `backslash_comma_caret`, `arguments_only`, `minimal_order_raw`;
   - the nine a leg built that **passed** round 4: truncate fields to length ≤ 1 · escape only when the
     whole field is a single special character · drop arguments beyond arity 2 · drop differentiated beyond
     arity 2 · map characters outside the alphabet to `a` · casefold · truncate order text to its first
     character · mangle field `"1"` to `"2"` · drop spaces.
4. `R3`: the comment you added, quoted.
5. The full unit suite — quote which files you ran.
6. `git status --short` and `git diff --stat`.

⛔ Do not run the comparator or the ablation battery; this round changes no behaviour and an unchanged
result would establish nothing.

## Constraints

- ⛔ No `git add`, no `git commit`, no other git write.
- ⛔ No code change to `engine_output_checks.py` — comment only, `R3`.
- ⛔ Do not launch Mathematica or `wolframscript`.

## Report back — under 12 lines

1. `R1`–`R3`: done or not, with the net line change.
2. The two variant tables from acceptance item 3, as measured.
3. ⭐ Anything in this directive you judge wrong.
