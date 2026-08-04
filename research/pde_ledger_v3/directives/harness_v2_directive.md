# HARNESS v2 — parse what the engines now emit, and turn on cross-engine comparison

⛔ **Repair in place. Do not rewrite.** Extend the existing module and its tests.

**Files (absolute):**
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/engine_output_checks.py`
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/test_engine_output_checks.py`
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/checks_S9.yaml`

**Real inputs to develop against** — ⛔ do not invent fixtures for the main path:
```
cd /var/projects/toy_physics
timeout 600 math -script research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl > /tmp/s9_wl.txt
cd research/pde_ledger_v3 && timeout 600 python3 scripts/S9_light_requires_shear_sympy_audit.py > /tmp/s9_py.txt
```
⚠ Mathematica has a **two-seat licence**; run one kernel at a time.

⭐⭐ **The governing constraint is unchanged: NO EXPECTED PHYSICS VALUE MAY BE STORED** in the module, the
config, or a test. Every comparison target is generated **at compare time** — from the other engine's run,
from the same engine's control packages, or from the registry. ⛔ Config maps **names**, never values.

---

## H1 — parse the value forms the engines now emit

Three forms currently come back `UNPARSED`, and one class is mis-reported as `UNKNOWN_SYMBOL`.

- **Mathematica boolean/relational expressions**: `a > 0 && b > 0 && kx^2 + ky^2 + kz^2 > 0 &&
  Element[kx, Reals]`. Support `&&`, `||`, `!`, the relational operators, and `Element[x, Reals]`.
- **Mathematica associations**: `<|sym -> {…}, sym -> {…}|>` → a mapping.
- **SymPy containers**: lists of tuples such as `[(expr, Matrix([[…]])), …]` and
  `[(expr, True), (expr, 2)]`.

⭐⭐ **And classify by KIND, because a boolean and a mapping HAVE NO DIMENSION.** Introduce an explicit
kind for each parsed value — expression / matrix / list / mapping / boolean / relation / integer — and
run the dimensional check **only** on the dimensionful kinds. A boolean is reported as
`NON_DIMENSIONAL`, ⛔ **never** as `UNKNOWN_SYMBOL` and ⛔ never as `UNPARSED`.
⚠ ⛔ **Do NOT implement this as a list of tag names to skip.** A name-based skip list is a denylist: the
next tag that needs it will be missed silently. ⭐ Decide from the **parsed kind**.

⛔ **Everything genuinely unparseable still fails loudly** with its raw text, and an unknown symbol inside
a dimensionful expression is still named and still an operational failure.

## H2 — configuration for the symbols the engines added

Add the coefficient symbols the engines now derive (`muF`, `muG`, `rhoZ`, …) to `derived_dimensions`,
each pointing at the **tag** that carries its dimension. ⛔ No values in the config.
⚠ Some symbols exist only in one control package; a symbol absent from a package must not be an error.

## H3 — ⭐⭐ TAG-SET PARITY, as a first-class check

⚠ **This is the gap that let a real defect through undetected.** An engine suppressed tags whose payload
duplicated an earlier tag, so quantities that were *correctly invariant across controls* **disappeared
from the output entirely**. Nothing caught it: names stayed unique, no untagged output, exit 0. The
control-response layer cannot see it, because a tag absent from **both** the main and the control
packages simply drops out of the comparison.

⭐ **Add a check that, for each engine, compares the tag-suffix set of every control package against the
main package** and reports, per package, the suffixes **missing** and **extra**.

⛔⛔ **Report it; do NOT fail on it.** A missing suffix can be legitimate — a quantity may genuinely not
exist for a given control. ⭐ This is a **triage list**, like the invariant list.
⚠ **State the reason it matters in the output**: a value present-and-identical means **INVARIANT**, a
real result; a value **absent** is indistinguishable from *never computed*.

## H4 — turn on LAYER 1, cross-engine comparison

This is the point of the exercise. The two engines emit different tag names and different serialisations
for the same physics.

Extend `cross_engine` config entries to:
```yaml
cross_engine:
  - quantity: transverse_speed_squared
    wl: WL_S9_CANDIDATE_SPEED_SQUARED1
    py: PY_S9_MAIN_SPEED_SQUARED_CANDIDATES
    # optional: how to reach the comparable value inside a container
    wl_select: scalar          # scalar | list | list_of_pairs_second | ...
    py_select: list
```
⭐ **Selection describes SHAPE, never value.** Its job is to reach the comparable object inside a
container — e.g. the `.wl` emits a scalar per root while the `.py` emits a list of `(root, value)` pairs.
⛔ Never let a selector name an expected value.

Compare **symbolically** — `simplify(a - b) == 0`, elementwise for containers, with an explicit
order-insensitive mode for multisets. Report `AGREE` / `DISAGREE` / `MISSING` / `UNPARSED` /
`SHAPE_MISMATCH`.

⭐ **Populate `checks_S9.yaml` with the load-bearing quantities both engines compute**, at minimum: the
factored determinant, the full root set, the transverse multiplicity, the candidate speed squared, the
homogeneity defect, the two coefficient dimensions and their difference, the implied speed dimension and
its difference from squared velocity, the dynamical-matrix route residual, and the bare-field control's
coefficient dimension. ⛔ Do not add an entry you have not confirmed both engines actually emit.

⛔⛔ **A `DISAGREE` is a FINDING, and the harness must EXIT 0 on it.** ⭐ Physics findings never change the
exit code; only operational failures do. A build agent iterating to exit 0 must not be able to make a
disagreement disappear.

## H5 — the report stays short

Counts, then the triage lists: `DISAGREE` rows first and in full, then `INVARIANT`, `NON_HOMOGENEOUS`,
`PARITY` gaps, `UNPARSED`. ⛔ Never print hundreds of passing rows. Keep the closing line stating what the
run does **not** establish.

---

## H6 — ⭐⭐ ABLE-TO-FAIL, and it is the point of the tests

Extend `test_engine_output_checks.py`. Each check must be **proved to fire** by corrupting **real** input:

1. delete a tag from one control package ⇒ it appears in that package's `PARITY` missing list;
2. change one engine's value for a mapped quantity ⇒ `DISAGREE`, with both operands shown;
3. two engines expressing the same thing differently (`a/b` vs `a*b**-1`, `{1,2}` vs `Matrix([[1],[2]])`)
   ⇒ `AGREE` — ⭐ this proves the comparison is symbolic and shape-aware, not textual;
4. a boolean-valued tag ⇒ `NON_DIMENSIONAL`, ⛔ not `UNKNOWN_SYMBOL`, ⛔ not `UNPARSED`;
5. an unknown symbol inside a **dimensionful** expression ⇒ still named, still an operational failure;
6. a `DISAGREE` present ⇒ **exit code is 0**;
7. genuinely malformed CAS text ⇒ `UNPARSED` and a **non-zero** exit.

⭐ **Tests 3 and 6 matter most**: 3 proves the comparison is real, 6 proves a finding cannot be silenced.

⚠ The **end-to-end** test must run against the **real** S9 `.wl` and `.py` outputs, ⛔ not fixtures.

Run `python -m pytest reduction/test_engine_output_checks.py -q` from
`/var/projects/toy_physics/research/pde_ledger_v3/` and iterate until it passes.

---

## Report back — under 30 lines

1. The three deliverable paths and the pytest result.
2. The **literal CLI report** from running the harness on the real `.wl` and `.py` outputs together.
3. One line per able-to-fail case: what you corrupted, what fired.
4. **The number of `cross_engine` entries you configured, and any quantity you wanted to map but could
   not** — say why.
5. ⭐ Any `DISAGREE` the real run produced. ⛔ **Do not investigate or fix it, and ⛔ do not adjust a
   selector to make it go away** — report it verbatim. A disagreement between two independently built
   engines is the most valuable output this tool can produce.

⛔ Do not commit to git. ⛔ Do not modify any engine, any `.wl`, or any file outside the three named.
