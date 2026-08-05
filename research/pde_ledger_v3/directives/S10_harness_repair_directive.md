# S10 harness repair — make the checker able to consume a D-sweep

**Primary file:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/engine_output_checks.py`
**Config:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/checks_S10.yaml`
**One engine edit, item H1 only:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`

⛔ **Do not commit.** ⛔ **Do not edit any other file** — in particular ⛔ not
`scripts/S10_brane_mode_spectrum_sympy_audit.py`, ⛔ not `checks_S9.yaml`, ⛔ not anything under `steps/`
or `paper/`.

⚠ **Read `directives/S10_SHARED_PHYSICS.md` §8 (tag grammar) before starting.** It defines the package
and tag naming this checker has to parse.

---

## ⭐ WHAT THIS IS

The checker was written for **S9**, which ran at a single fixed spatial dimension and named its control
packages `X6`, `X7`, `X8`. **S10 sweeps `D = 2,3,4,5` in one run and names its controls in words**
(`XFORM_FULLGRAD`, `XFORM_DIVONLY`, `XFORM_SIGNFLIP`, `XFORM_ANISO`, `XCOEF_SCALE`).
⇒ Every item below is that single mismatch surfacing in a different layer.

⛔ **No physics result is in question and none of your fixes may change one.** The two engines are
finished and frozen. Your job is to make the instrument able to read them.

---

## ⛔⛔ H0 — BLOCKING AND FIRST: the ablation layer has been comparing NOTHING

```python
CONTROL_TAG_PATTERN = re.compile(r"^(?P<base>.+)_X(?P<control>[1-9][0-9]*)_(?P<suffix>.+)$")
```

⚠ **Measured: zero S10 tags match this.** `_package_layout` therefore returns empty, and both layers built
on it report:

```
CONTROL_RESPONSE: compared=0 responsive=0 invariant=0 unparsed=0 unpaired=0
TAG_PARITY: packages=0 gaps=0
```

⛔⛔ **That is a check reporting SILENCE as if it were health.** `CONTROL_RESPONSE` is the layer that
partitions main-vs-control tags into RESPONSIVE and INVARIANT — ⭐ it is the harness's **ablation** test,
the one thing that distinguishes a computed output from a typed one. It has not run.

**Fix, in two parts, and ⭐ both are required:**

**H0a — declare the packages in the config instead of guessing them from a regex.** Add to
`checks_S10.yaml`:

```yaml
main_package: MAIN
control_packages: [XFORM_FULLGRAD, XFORM_DIVONLY, XFORM_SIGNFLIP, XFORM_ANISO, XCOEF_SCALE]
```

and have `_package_layout` use them when present. ⚠ A tag is `<ENGINE>_<STEP>_<PACKAGE>_<SUFFIX>`; with the
package named explicitly the suffix is whatever follows it, so `MAIN_D3_Q1_LAGRANGIAN` and
`XFORM_ANISO_D3_Q1_LAGRANGIAN` share the suffix `D3_Q1_LAGRANGIAN` and are compared **at matching `D`**
automatically. ⭐ Keep the existing `_X<digits>_` behaviour as the fallback when the config declares
nothing, so `checks_S9.yaml` keeps working unchanged.
⚠ **Match by LONGEST declared prefix.** Package names contain underscores (`XFORM_FULLGRAD`), so a
shortest-match or greedy-base rule will split them in the wrong place and produce plausible-looking
nonsense suffixes.

**H0b — ⭐⭐ a layer that was ASKED to compare something and compared nothing must be an OPERATIONAL
FAILURE, not a clean line.**
If `control_packages` is declared and any declared package contributes **zero** tags, or if the layer ends
with `compared == 0`, fail with a message naming the packages it could not find. ⛔ The same for
`TAG_PARITY` with `packages == 0`. ⚠ This is the general defect, not a detail of H0a: **a check whose
failure mode is a quiet zero is worse than no check**, because a reader scores it as passing.

⛔⛔ **BUT SCOPE IT PRECISELY, and this is a correction to an earlier version of this directive that a
review leg was right to reject.** The rule is **declared-but-compared-nothing**, ⛔ **not** "any count of
zero". ⚠ `REGISTRY_RESIDUAL: configured=0` and `CROSS_ENGINE: configured=0` mean *nothing was asked for*,
which is a legitimate, already-visible state. ⛔ Turning those into failures converts an honest silence
into noise and would break `checks_S9.yaml`. ⇒ ⭐ **Fail only when the config names something the output
does not supply.**

**H0c — ⛔ control response and dimensions run on ONE engine only.** `run_checks` passes `default_values`
to both `check_control_response` and `check_dimensions`; only `check_tag_parity` receives every engine.
⇒ ⛔⛔ **the second engine is never ablation-checked and never dimension-checked** — half of a two-engine
corpus, unmeasured, in the two layers that matter most.
⭐ **Fix: run both layers per engine and report per engine.** ⚠ Expect this to surface findings in the
engine that has never been checked; ⛔ those are results, ⛔ not reasons to narrow the fix.

---

## ⛔ H1 — Mathematica diagnostics reach stdout (⭐ the only engine edit in this directive)

`Solve::svars` prints **10 message lines plus 10 blank lines** into engine 1's output. The harness parses
them as tags.

⭐ **Fix: `Quiet` the message at the solve site — and ⛔ only that message at that site.** The engine
already emits `..._SOLVE_SVARS_MESSAGE` tags recording that the condition occurred; ⛔ **do not remove,
weaken, or stop emitting those.** ⭐ The information stays in the tag stream; only the raw text leaves
stdout.

⛔ Do not add a blanket `Quiet` or `Off` around anything larger than the individual `Solve` call.
⛔ Do not change what is solved, or how.

⚠ **Then re-run engine 1** and overwrite `mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`.
⭐ Launch it in the **background**, never in the foreground under `timeout`. Budget **600 s**;
⛔ if it exceeds that, stop it and report, ⛔ do not wait for it. ⚠ An orphaned `WolframKernel` leaks
memory without limit — confirm with `ps -eo pid,rss,etime,comm | grep WolframKernel` that none survives
your run.

⭐ **Acceptance:** `grep -cvE '^[A-Z][A-Z0-9_]*:' <out file>` returns `0`.

---

## ⛔ H2 — `derived_dimensions` cannot take a per-coefficient list

`_dimension_table` passes the payload straight to `_as_dimension`, which requires a single
three-component vector. The coefficient-family tags carry **one vector per coefficient**, because Q6 was
amended so that *which* coefficients exist is package-dependent.

⭐ **Fix: collapse a list of dimension vectors to the one they share, and ⛔ FAIL if they do not share
one.** The config maps the family to a **single** symbol name, so a disagreement means the mapping is
ill-posed. ⛔ Do not take the first entry.

⚠ **A prototype of exactly this shape has been run and it works** — for reference, not to copy verbatim:

```python
def _collapse_dimension_list(value: object, label: str) -> object:
    sequence = _sequence_form(value)
    if sequence is None:
        return value
    entries = [item for item in sequence if _sequence_form(item) is not None]
    if not entries or len(entries) != len(sequence):
        return value          # a bare vector, or something else entirely: leave it alone
    first = _as_dimension(entries[0], f"{label}[0]")
    for index, item in enumerate(entries[1:], start=1):
        other = _as_dimension(item, f"{label}[{index}]")
        if not _dimension_equal(first, other):
            raise HarnessError(f"{label} carries disagreeing dimension vectors: ...")
    return list(first)
```

⛔⛔ **Its weakness cannot be fixed by a better shape heuristic, and ⛔ do not try.** ⚠ **Measured by a
review leg:** a 3×3 matrix whose rows happen to be equal collapses to a single row under the prototype,
and it still does under any stricter *scalar-component* or *length* test — the two inputs are genuinely
indistinguishable by shape.
⇒ ⭐ **The config must say which it is.** Give `derived_dimensions` an explicit form, e.g.

```yaml
derived_dimensions:
  rhoBr: {tag: WL_..., shape: family}   # a list of per-coefficient vectors
  muX:   {tag: WL_..., shape: vector}   # a single vector
```

⭐ Keep the plain `name: TAG` form working as `shape: vector`, so `checks_S9.yaml` is untouched.
⛔ **Collapse only when the config said `family`.** ⚠ A payload whose shape contradicts its declared form
is an error, ⛔ not something to auto-detect.

---

## ⛔ H3 — the dimension table is single-`D`, and S10 is a `D`-sweep

`checks_S10.yaml` currently points `derived_dimensions` at the **`_SPECIALIZED`** tags, whose payloads are
numeric at one `D`. ⛔ **That silently applies one package's dimensions to every other package in the
sweep.**

⭐ **Fix: point at the SYMBOLIC tags instead** — the ones whose components are expressions in the brane
dimension. `_as_dimension` already accepts `sp.Expr` components and `_dimension_equal` already compares
with `simplify(a - b) == 0`, so a symbolic table needs no new arithmetic.

⛔⛔ **DO NOT STATE THIS AS "CORRECT AT EVERY `D`". An earlier version of this directive did and a review
leg was right to reject it.** ⭐ What is true, and all that is claimed: for terms built from **this table
plus the primitives**, a homogeneity comparison is an identity in the symbol or it is not, so the verdict
does not depend on specialising it. ⚠ **A symbolic table and a numeric one CAN disagree**, and the leg
constructed the case: summands carrying `(-braneDimension, 0, 1)` against `(-3, 0, 1)` are
non-homogeneous symbolically and homogeneous once `D = 3` is substituted. ⇒ that is a reason to prefer
the symbolic table, ⛔ not a proof that it cannot be wrong.

⭐⭐ **Therefore the cross-package check is MANDATORY, not optional.** The table is built from **one
package's** tags and applied to all of them. ⭐ Emit a **report line** whenever another package's
corresponding tag disagrees with the one used, naming both packages and both vectors.
⛔ Do not make it fatal — a control that changes the action's *form* may legitimately change a
coefficient's dimension — but ⛔ **do not let it be silent either.** ⭐ Record the one-package limitation
in a comment in the config as well.

---

## ⛔ H4 — 220 payloads are UNPARSED, all from the Mathematica reader

Four distinct gaps, all the same root cause — the reader was written against a fixed 3+1 dimensional
output:

1. ⭐ **`Derivative[i,j][f][x1,x2,t]` of arbitrary arity.** The reader hardcodes
   `Derivative[i,j,k,l][f][t,x,y,z]`. **This is the bulk of the 220.**
2. ⛔ **`Element[x, Integers]`.** ⚠ **An earlier version of this directive blamed conjunctions and that
   was wrong** — `Element[k1, Reals]` inside a `&&` chain parses fine, tested four ways. The reader
   accepts **only** `Element[x, Reals]` (`:547-554`), and the engine also asserts
   `Element[braneDimension, Integers]`.
3. Association syntax `<| "key" -> value |>`; the tokenizer rejects the string key.
4. A marker head applied to a **list** argument, e.g. `suppliedInPlaneField[{u1[...], u2[...]}, ...]`.
5. ⭐ **`Piecewise[{{value, condition}}, default]`** — reported as *"argument of Piecewise is not a scalar
   expression"*. ⚠ This one comes from **S9**, whose harness run exits 2 on it alone; fixing it makes that
   run clean. ⛔ Do not change anything else about the S9 configuration.

⭐ **Fix all five.** ⛔ Where you genuinely cannot parse a construct, the payload must be reported as
UNPARSED with its reason — ⛔ never silently dropped and ⛔ never coerced into something that parses.

---

## ⛔ H5 — 794 unknown-symbol reports, in four classes that want different treatment

Run the harness and read the list yourself before deciding; the classes are:

- **Marker sentinels** the spec requires the engines to emit where a value is undefined
  (`undefinedRatio`, `notApplicable`, `decidedEmpty`, `decidedNonempty`, `exactlyDetermined`,
  `codingConsistencyOnly`, `quadraticFormRoute`, `suppliedNoDissipation`, `suppliedQuadraticResponse`, …),
  plus Mathematica's `Indeterminate`.
  ⭐ **Add a `marker_symbols:` config key.** A payload that **is** a marker is **non-dimensional** — report
  it in `NON_DIMENSIONAL`, ⛔ not as an unknown symbol.
  ⛔⛔ **A marker appearing INSIDE a larger expression is NOT the same thing and must still be reported** —
  a sentinel that has leaked into arithmetic is a defect, and ⛔ you must not silence it by declaring the
  name globally dimensionless.
  ⛔⛔ **AND THE REAL RISK, which a review leg named:** an engine can emit a **pure sentinel where a
  dimensionful answer was required** — a failed dimension arriving as `Indeterminate` — and bucketing it
  as non-dimensional makes that **drop out of the failure counts entirely.**
  ⇒ ⭐ **`Indeterminate` gets its OWN report category, separate from designed markers**, because it means
  *the engine could not determine this*, ⛔ not *this does not apply*. ⭐ Count and list both categories
  **per tag family**, so a family that is entirely sentinels is visible as such.
- **Dimension-exponent unknowns** — `dimRhoLength`, `dimRhoTime`, `dimRhoMass`, `dimMuLength`, `dimMuTime`,
  `dimMuMass`, `dimSRho*`, `dimScale*`. These are the **unknowns of the Q6 dimension solve**; they live in
  exponent position and are not quantities at all. ⭐ Give them their own config key and their own report
  line. ⛔ Do not dump them into `dimensionless`, which means something else.
- **Gradient symbols** `g1x2`, `g2x3`, … — Q7's independent displacement-gradient symbols.
  ⭐ These belong in `primitive_dimensions` as a **definitional convention**, with a comment saying which
  convention and why. ⛔ Do not derive their dimension from anything.
- ⚠ Anything left over: ⭐ **report it, do not classify it.** Bring it back in your write-up.

---

## ⚠ H6 — 9 NON_HOMOGENEOUS reports, all from Q8b's numeric witness point

Every one is a stratum tag where the engine substituted a **bare number** for a dimensionful wavevector
component (`k1 → -75`), which destroys homogeneity by construction.

⛔ **Do not change the engine, and ⛔ do not suppress the reports.** ⭐ Give them their own report
section, or an explicit declared exclusion naming the reason. ⚠ **The count must still be visible** — a
reader has to be able to see that nine tags were set aside and why.

---

## ⛔⛔ ACCEPTANCE — ⭐ ABLATE THE HARNESS. ⛔ Its own report is not evidence.

For **each** of the layers you touched — control response, tag parity, dimensions, cross-engine, and the
zero-count guard of H0b — you must demonstrate the layer is **able to fail**:

⭐ Take a **copy** of an engine output file, **corrupt one payload in it**, re-run, and show the layer
moves. ⛔ Corrupt the copy, ⛔ never the committed output. Then show the uncorrupted run does not.

⚠ For the H0b guard specifically: show that a config declaring a package **which does not exist in the
output** produces an operational failure, and that the correct config does not.

⛔ **A layer you cannot make fail is a layer you have not fixed**, whatever its report says.

### ⛔⛔ AND CORRUPT-ONE-PAYLOAD IS NOT SUFFICIENT — a review leg exhibited four ways to pass it broken

⭐ **All four of these are additionally required. ⛔ Skipping one is not reporting "partially done"; it is
reporting a layer as fixed that has not been shown to work.**

1. ⭐ **Ablate EVERY control package, one at a time — ⛔ not one of them.** A layer that only ever compares
   the first control package passes a single-corruption demo while never touching the others.
   ⚠ **I found the same hole in my own check before the leg reported it**, which is how sure you should be
   that it is real.
2. ⭐ **A SYMBOLIC-EQUIVALENCE probe.** Present the same object in two different but equivalent renderings
   — `x + y` against `y + x`, `1/2` against `Rational[1,2]` — and require the comparison to report
   **AGREE**. ⛔ A raw-string equality layer passes every corruption test and still fails this.
3. ⭐ **An EMPTY-INPUT probe.** Feed a layer an output file with none of its tags and show it complains.
   ⛔ Corrupt-and-move never establishes this, and it is the exact defect H0 exists for.
4. ⭐ **Corrupt into a WRONG VALUE, not into something UNPARSED.** A corruption that merely breaks parsing
   moves the UNPARSED count without ever exercising the comparison you are trying to test.

⭐ **Also run `python3 reduction/test_engine_output_checks.py`** (or via pytest) and report the result.
⛔ If an existing test now fails, ⛔ do not edit the test to match your change until you have told me it
failed and why.
⚠ **Expect at least one to fail:** `test_absent_parity_exclude_preserves_report_bytes` freezes the report
text byte-for-byte, including the `configured=0` lines your new report categories sit beside. ⭐ Report it
and say what changed; ⛔ do not quietly re-baseline it.

⭐ **Finally, run the S9 configuration too** and report its summary line — it is the regression check that
your changes did not break the one step whose harness config already worked:

```
timeout 900 python3 reduction/engine_output_checks.py --config reduction/checks_S9.yaml \
  --output wl=mathematica/out/S9_light_requires_shear_mathematica_audit.out \
  --output py=scripts/out/S9_light_requires_shear_sympy_audit.out
```

**Command to run, exactly:**

```
cd /var/projects/toy_physics/research/pde_ledger_v3 && \
timeout 900 python3 reduction/engine_output_checks.py \
  --config reduction/checks_S10.yaml \
  --output wl=mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out \
  --output py=scripts/out/S10_brane_mode_spectrum_sympy_audit.out
```

---

## Report back — ⛔ under 35 lines, plus the ablation table

1. One line per `H0a`,`H0b`,`H1`–`H6`: fixed / partially / not, with file and line numbers.
2. The ablation table: layer, what you corrupted, what moved. ⛔ One row per layer, no exceptions.
3. Exit code and wall-clock of the acceptance command; the unit-test result.
4. ⛔ **Do not report what any physics value came out to be.** Counts of tags, agreements and failures are
   fine; ⛔ the values themselves are not yours to summarise.
5. ⭐ Anything in this directive you believe is **wrong**, and anything you changed that I did not ask for.
   ⭐ This is wanted.
