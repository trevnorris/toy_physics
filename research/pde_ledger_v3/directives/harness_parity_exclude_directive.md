# Harness — add engine-local tag exclusion to the parity layer

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/engine_output_checks.py`
**Tests:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/test_engine_output_checks.py`

⛔ **Do not commit.** ⛔ Do not modify any other file. ⛔ Do not read `research/pde_ledger_v3/steps/` or
`research/pde_ledger_v3/paper/`.

## The problem

The parity layer (`check_tag_parity`) compares the tag sets of two engines and reports every name present
in one and absent from the other. ⚠ **Some tags CANNOT exist in both engines by design:** one engine reads
a registry the other is barred from, and each CAS emits its own solver-condition and route tags. The
convention already in use is an infix — a tag name containing `_LOCAL_` immediately after the engine prefix
is **engine-local by declaration**.

⛔ Reporting those as parity gaps trains a reader to ignore the parity report, which is worse than not
having one.

## What to add

1. ⭐ A config key **`parity_exclude`**: a list of substrings. Any tag whose name contains one of them is
   omitted from the parity comparison **only**. ⛔ It must still take part in every other layer —
   cross-engine, dimensions, control response.
2. ⭐ Report the exclusion **in the output**: how many tags were excluded, per engine, and by which pattern.
   ⛔ A silently smaller gap count is exactly the failure this is meant to avoid.
3. ⚠ Default when the key is absent: **exclude nothing**, and behave exactly as today.

## ⛔⛔ Constraints

- ⛔ **A physics finding must still EXIT 0.** This change must not alter any exit code.
- ⛔ Do not touch the cross-engine, dimension, or control-response layers.
- ⭐ Add tests: a tag matching an exclusion is absent from the parity gaps **and still present** in the
  other layers; the excluded count is reported; and with no `parity_exclude` key the behaviour is byte-
  identical to before.
- ⭐ Run the existing test suite and confirm it still passes.

## Report — ⛔ under 12 lines
1. What you changed, with line numbers.
2. The tests you added and that the suite passes.
3. ⭐ Anything about this design you think is wrong — in particular, whether excluding by **substring**
   is too blunt and something narrower would be safer. ⭐ This is wanted.
