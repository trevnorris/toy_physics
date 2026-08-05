# S10 harness repair, part B — make the reader see the ACTION

**Primary file:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/engine_output_checks.py`
**Config:** `/var/projects/toy_physics/research/pde_ledger_v3/reduction/checks_S10.yaml`

⛔ **Do not commit.** ⛔ **Do not edit** either engine, either `.out`, `checks_S9.yaml`, or anything under
`steps/` or `paper/`.
⚠ **Part A has already landed** — the control/ablation layer, the coverage guards, the per-engine
dimension and control checks, and the Mathematica `Quiet`. ⛔ Do not redo or revert any of it. ⭐ If you
believe part A got something wrong, ⭐ **say so and leave it alone.**

---

## ⭐ WHY THIS IS WORTH DOING, AND WHERE IT STOPS

⚠ **220 payloads are UNPARSED**, and **58 of the unparsed cross-engine rows include `Q1_LAGRANGIAN` and
`Q1_EULER_LAGRANGE_SYSTEM`.** ⇒ ⛔⛔ **The two engines' ACTIONS are currently never compared** — and the
action is the object every control is supposed to re-enter at. ⭐ **That is the whole reason for part B.**

⛔ **What part B is NOT:** a redesign of the parity layer, a new comparison architecture, or an attempt to
make dimensional homogeneity a strong check. ⚠ **The homogeneity test is structurally vacuous by design**
— the displacement's dimension cancels in every coefficient ratio, so no cleverer version exists — and
part A already made the vacuous portion visible. ⭐ **Label it and move on.** ⛔ Do not try to rescue it.

---

## ⛔⛔ B1 — the Mathematica reader. ⭐ This is the item that matters.

Six constructs the reader rejects, **all measured in the committed output**:

1. ⭐⭐ **`Derivative[i,j][f][x1,x2,t]` of ARBITRARY arity.** The reader hardcodes
   `Derivative[i,j,k,l][f][t,x,y,z]` — S9's fixed 3+1 shape. ⚠ S10 sweeps `D = 2,3,4,5`, so the
   derivative's arity and the argument list both vary. **This is the bulk of the 220 and it is what
   blocks the action.** ⭐ Fix it generally: any arity, any argument list, ⛔ no hardcoded coordinate
   names or ordering.
2. `Element[x, Integers]`. ⚠ The reader accepts **only** `Element[x, Reals]`, and the engine asserts
   `Element[braneDimension, Integers]`. ⛔ **Do not "fix" conjunctions** — `Element[k1, Reals]` inside a
   `&&` chain parses correctly today; I tested four forms.
3. Association syntax `<| "key" -> value |>`; the tokenizer rejects the string key.
4. A marker head applied to a **list** argument, e.g. `suppliedInPlaneField[{u1[...], u2[...]}, ...]`.
5. `Piecewise[{{value, condition}}, default]` — *"argument of Piecewise is not a scalar expression"*.
   ⚠ This one is **S9's** single UNPARSED; fixing it makes that run clean. ⛔ Change nothing else about S9.
6. ⚠ `ConditionalExpression[expr, condition]` — currently parses as a **generic function**, which
   ⛔ silently carries a condition into a dimension walk. ⭐ Handle it explicitly: take the expression,
   and emit the condition somewhere visible, ⛔ never drop it silently.

⭐ **Acceptance for B1:** report the UNPARSED count before and after, **and** confirm that
`Q1_LAGRANGIAN` and `Q1_EULER_LAGRANGE_SYSTEM` now parse **in both engines and at every `D`**.
⛔ Where a construct genuinely cannot be parsed, the payload stays UNPARSED **with its reason** —
⛔ never silently dropped, ⛔ never coerced into something that parses.

---

## ⛔ B2 — `derived_dimensions`: declare the shape, ⛔ do not guess it

**B2a — the payload is a LIST of per-coefficient dimension vectors** and `_as_dimension` requires a single
triple. ⭐ Collapse a family to the vector its members share, and ⛔ **FAIL if they do not share one** —
the config maps the family to **one** symbol name, so a disagreement means the mapping is ill-posed.
⛔ Do not take the first entry.

**B2b — ⛔⛔ the collapse must NOT be a shape heuristic.** ⚠ **Measured by a review leg:** a 3×3 matrix
whose rows happen to agree collapses *successfully* and installs a **wrong dimension with no complaint**;
a zero 3×3 becomes `(0,0,0)`. ⚠ And at `D = 3` the real inertial payload **is** exactly three
three-vectors, so this fires on live data. ⇒ ⭐ **the config declares which it is:**

```yaml
derived_dimensions:
  rhoBr: {tag: WL_..., shape: family}
  muX:   {tag: WL_..., shape: vector}
```

⭐ Keep the plain `name: TAG` form working as `shape: vector` so `checks_S9.yaml` is untouched.
⛔ **Collapse only when the config said `family`.** ⭐ Reject `sp.MatrixBase` with `rows>1 and cols>1`
**before** consulting entry shapes. ⚠ A payload whose shape contradicts its declaration is an **error**.
⚠ Note that B1 item 3 makes **Association** payloads reachable here for the first time — ⭐ handle or
reject them explicitly.

**B2c — point the config at the SYMBOLIC tags**, whose components are expressions in the brane dimension,
⛔ not the `_SPECIALIZED` ones. ⚠ The `_SPECIALIZED` payloads are numeric at one `D`, so the current
config **silently applies one `D`'s dimensions to every package in the sweep.**
⛔⛔ **Do not restate this as "correct at every `D`".** ⭐ What is true: for terms built from this table
plus the primitives, a comparison is an identity in the symbol or it is not, so the verdict does not
depend on specialising it — and a review leg verified the two tables give **byte-identical verdicts across
all 2983 tags** while also constructing a case where a numeric table passes an accidental coincidence a
symbolic one catches. ⇒ symbolic is **strictly stronger**, ⛔ not infallible.

**B2d — the table is built from ONE package and applied to all.** ⭐ **Mandatory:** emit a report line
whenever another package's corresponding tag disagrees with the one used, naming both packages and both
vectors. ⛔ Not fatal — a form control may legitimately change a coefficient's dimension — ⛔ but not
silent either. ⭐ Record the limitation in a config comment.

---

## ⚠ B3 — markers and strata: ⭐ REPORTING ONLY. ⛔ Do not build machinery.

**B3a — markers.** The engines emit explicit sentinels where a value is undefined by design
(`undefinedRatio`, `notApplicable`, `decidedEmpty`, `decidedNonempty`, `exactlyDetermined`,
`codingConsistencyOnly`, `quadraticFormRoute`, `suppliedNoDissipation`, `suppliedQuadraticResponse`, …).
⭐ Add a `marker_symbols:` config key. A payload that **is** a marker is non-dimensional.

⛔⛔ **Two things this must not become:**
- ⛔ **A marker INSIDE a larger expression is a different thing and must still be reported.** ⚠ 188 tags
  carry one. A sentinel that has leaked into arithmetic is a defect.
- ⛔⛔ **Do not bury them in an unnamed count.** ⚠ A review leg's point, and it is right: an engine that
  answered `notApplicable` to every hard question would then be indistinguishable from one that computed
  everything. ⇒ ⭐ **a named, counted report section listing the tags**, and ⭐ **`Indeterminate` gets its
  OWN category** — it means *the engine could not determine this*, ⛔ not *this does not apply*.
  ⚠ Measured: 93 whole-payload markers, `Indeterminate` among them.

**B3b — the dimension-exponent unknowns** (`dimRhoLength`, `dimMuTime`, `dimSRho*`, `dimScale*`, …) are
the **unknowns of the Q6 dimension solve**, in exponent position, ⛔ not quantities. ⭐ Their own config
key and their own report line. ⛔ Do not put them in `dimensionless`, which means something else.

**B3c — the gradient symbols** `g1x2`, `g2x3`, … are displacement gradients ⇒ dimensionless by definition.
⭐ Put them in `primitive_dimensions` with a comment saying which convention and why.

**B3d — the stratum tags.** ⚠ Part A should have made the non-homogeneous count visible. The nine records
come from **6 tags**, but there are **102 `STRATUM` tags carrying the same numeric witness substitution**
— ⇒ ⛔ **"nine" is not a count of artefacts, it is a count of nodes that survived a zero-filter.**
⭐ **Key any exclusion on the CLASS** (stratum / Q8b tags), ⭐ let the count float, ⭐ print the names.
⭐ **And add the discriminator:** re-run the same tag with the numeric witness replaced by a symbol of the
right dimension — homogeneity restored ⇒ artefact; not restored ⇒ real. ⛔ Without it the exclusion is a
declaration, not a finding.

---

## ⛔⛔ ACCEPTANCE — ⭐ ablate. ⛔ The report is not evidence.

For every layer you touch, demonstrate it is **able to fail** by corrupting a **copy** of an engine
output and showing the layer moves — ⛔ corrupt the copy, ⛔ never the committed file.

⭐ **Specifically required:**
1. ⭐ **B1:** a payload that parses only because of your fix, corrupted into a **wrong value**, must move
   the comparison. ⛔ Corrupting it into something unparsable proves nothing — that only moves the
   UNPARSED count.
2. ⭐ **B2:** one entry of a declared `family` payload changed so the members disagree **must fail
   loudly**; and a `vector`-declared payload handed a matrix must error.
3. ⭐ **B3:** a marker moved from whole-payload position into an arithmetic expression must be reported.

**The commands:**

```
cd /var/projects/toy_physics/research/pde_ledger_v3 && \
timeout 900 python3 reduction/engine_output_checks.py --config reduction/checks_S10.yaml \
  --output wl=mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out \
  --output py=scripts/out/S10_brane_mode_spectrum_sympy_audit.out

timeout 900 python3 reduction/engine_output_checks.py --config reduction/checks_S9.yaml \
  --output wl=mathematica/out/S9_light_requires_shear_mathematica_audit.out \
  --output py=scripts/out/S9_light_requires_shear_sympy_audit.out

timeout 900 python3 -m pytest reduction/test_engine_output_checks.py -q
```

⭐ Report all three. ⚠ Part A repointed the tests at the committed `out/` files; ⛔ if a test fails, ⭐ tell
me it failed and why — ⛔ do not edit a test to match your change.

## Report back — ⛔ under 30 lines, plus the ablation table

1. One line per `B1`(×6), `B2a`–`B2d`, `B3a`–`B3d`: fixed / partially / not, with line numbers.
2. **UNPARSED before and after**, and the explicit confirmation that `Q1_LAGRANGIAN` and
   `Q1_EULER_LAGRANGE_SYSTEM` parse in both engines at every `D`.
3. The ablation table: layer, what you corrupted, what moved.
4. The three command results.
5. ⛔ **Do not report what any physics value came out to be.** Counts are fine.
6. ⭐ Anything here you believe is wrong, and anything you changed that I did not ask for. ⭐ Wanted.
