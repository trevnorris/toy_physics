# Two edits behind a declaration that has already been made — a DECISION LIST, not a design task

**Do not commit. Do not modify any file outside the two named in D1 and D2.**

The nine `D=3` gradient-symbol naming exceptions are already written into
`research/pde_ledger_v3/reduction/checks_S10.yaml` and are not yours to change. Two things behind that
declaration are now stale or missing. Fix exactly those.

## What was declared, and where it came from

Both engines replace each first derivative of a field by a coordinate with an independent symbol before the
`§Q7` comparison, and both index that symbol **(coordinate, field)** in that order. The two spellings are
therefore related by

```
WL  g{r}x{c}   ==   PY  g{r}{c}        for r, c in 1..3        both standing for  d(u_c) / d(x_r)
```

sourced from the engine definitions — `mathematica/S10_brane_mode_spectrum_mathematica_audit.wl` lines
`148-152`, `216`, `1304-1311`, and `scripts/S10_brane_mode_spectrum_sympy_audit.py` lines `309`, `1536`,
`1538-1542`. The provenance and its limitation are written into the config above the nine entries. **Read
that comment before you start.** It is the premise both edits below rest on.

---

## D1 — `reduction/test_engine_output_checks.py`

`test_s10_declares_only_verified_dimension_and_scale_spellings` ends with an assertion that
`main_d3_q7_stiffness` reports `NAMING_MISMATCH`. That assertion recorded a declaration deliberately **not**
made. It has now been made, so the assertion is stale.

Replace it with assertions that pin what the declaration does to that row:

- its comparison status;
- the **complete set** of `naming_applied` entries it carries — written out by hand from the equation above
  and from which of the nine symbols occur in that package's payload, **not** pasted from a harness run;
- that `undeclared_spellings` is empty.

The test name now understates its scope. Rename it to cover all three declared families.

**Do not weaken any other assertion in that test**, and do not delete the
`main_d3_q2_downstream_route` `NAMING_MISMATCH` assertion — that one records a declaration still not made.

## D2 — `reduction/harness_ablation.py`

Add one numbered acceptance item, continuing the existing numbering and following the module's existing
conventions exactly (engine-free, in-memory mutation only, `require(...)` then a printed `ACCEPTANCE n PASS`
line).

**What it must establish.** Run the declared `§Q7` cross-engine rows under four in-memory variants of the
nine gradient exceptions, leaving every other exception alone:

| variant | the nine entries become |
|---|---|
| `AS_DECLARED` | unchanged |
| `TRANSPOSE` | `WL g{c}x{r} == PY g{r}{c}` |
| `DERANGED` | `WL g{r}x{(c mod 3) + 1} == PY g{r}{c}` |
| `ABSENT` | removed |

For each variant, tally the declared `§Q7` rows by comparison status and **print the tally**. Print all
four; a tally that is only asserted and never printed carries no information.

**Assert exactly this and nothing stronger:** that `DERANGED` and `ABSENT` each produce a tally that differs
from `AS_DECLARED`. That is the property that makes the declaration load-bearing rather than decorative.

**Do not assert anything about `TRANSPOSE`.** Print its tally and leave it unasserted. Whatever it comes
out to is a statement about how much the payloads can discriminate, and this battery's job there is to
publish it, not to enforce it. If you find yourself wanting to assert it, re-read the config comment.

---

## Out of scope — do not touch, even if you see something wrong

`reduction/checks_S10.yaml` · either engine · either committed `.out` · `steps/` · `paper/` ·
`REBUILD_HANDOFF.md` · any other acceptance item in the battery · any other test.

If you spot a defect in any of those, **report it and leave it alone.** That is wanted.

## Acceptance — run these and paste literal stdout

1. The full unit suite from `reduction/`. A one-line cause for any failure.
2. `harness_ablation.py` in full, literal, every acceptance line.
3. The harness on the S10 config against both committed outputs — paste the literal `CROSS_ENGINE:` and
   `CROSS_ENGINE_COVERAGE:` lines only.
4. `git status --short` and `git diff --stat`, confirming nothing outside D1 and D2 moved.

Do not iterate toward any particular tally. Whatever the battery reads is the measurement.

## Constraints

- No script over 600 seconds. A timeout means reformulate; never raise the limit.
- **Do not launch Mathematica or `wolframscript`.** Two licence seats, and the WL output is committed and
  unchanged. Read it as text.
- Do not run `git add`, `git commit`, or any other git write.

## Report back — under 20 lines

1. D1 and D2: done / partially / not, with line numbers.
2. The acceptance output.
3. Anything you measured that contradicts this directive.
4. Anything wrong that you did not fix because it was out of scope.
