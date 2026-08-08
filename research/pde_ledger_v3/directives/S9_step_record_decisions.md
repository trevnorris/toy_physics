# S9 step record — decision list

⭐ **You write the prose. I supply the decisions and the measurements.** ⛔ Do not add a measurement that is
not here; ⛔ do not soften or strengthen one. If something reads as incomplete, **say so in your report**
rather than filling it in.

## The file

`research/pde_ledger_v3/steps/S9_light_requires_shear.md` — **it already exists.** ⭐ Add a section covering
the export-chain rebuild, and ⭐ correct anything the rebuild made false. ⛔ Do not restructure the document.
⛔ **No other file.** ⛔ Do not commit.

## Context you need

S9's SymPy engine now writes `scripts/S9_exports.py`, a flat 141-entry `LEDGER` the next step imports.
⛔ **The derivation did not change** — measured: of 55 changed lines, 51 are the same code with a
`# TAG · description` comment appended; the rest are a declaration split, the emit loop, and the action line
reformatted. The Wolfram engine changed in **emitted name strings only** and re-runs byte-identical.

## What the section must say

**1 · What was built.** Annotations on every declaration; twelve objects emitted under one standard name in
both engines; the generated flat `LEDGER` (`display, value, dim?, class, step`); a knob extractor.
Class tally: **KNOB 2 · PREMISE 7 · DERIVED 132**. Three entries carry a computed `dim`.

**2 · The verification, and its real width.** Seven independent review legs each wrote their **own CAS
derivation before opening the artifact**, and all seven reproduce every S9 value. ⭐ A **FORM** ablation
(curl-only → gradient-elastic) moves **86** exported values; ⚠ **name the form**, because the count is
form-dependent — divergence-only moves 102.

**3 · ⛔⛔ THE CROSS-ENGINE TABLE IS SMALLER THAN TWELVE, and this is the most important paragraph.**
- **11 counted, 1 excluded.** The exclusion is `DYNAMICAL_MATRIX_ROUTE_RESIDUAL`, whose guard is exact
  equality against the zero matrix, so a disagreement becomes a build failure rather than a committed
  comparison.
- ⛔ **Eleven rows is eleven COMPARISONS, at most NINE independent computations.**
  `COEFFICIENT_DIMENSION_DIFFERENCE` buys only that both engines typed `mu − rho` and not `rho − mu`;
  `SPEED_DIMENSION_DIFFERENCE` buys only that both posited the same velocity-squared reference.
- ⛔ **Under the form control, 9 of the 10 exported standard rows are BYTE-IDENTICAL** — only
  `FACTORED_DETERMINANT` and `FULL_ROOT_MULTISET` move. ⇒ **the table's discriminating power against the
  ordinary-elastic alternative is TWO ROWS.**

**4 · ⛔⛔ A MEASURED COMMON-MODE BLINDNESS.** Both engines type `[u] = L` independently and **neither
emits it**. Doubling it in a scratch copy of **each** engine moves `INERTIA_`, `STIFFNESS_` and
`BARE_FIELD_COEFFICIENT_DIMENSION` to the **same wrong values**, both **exit 0**, and **they still agree**.
⇒ what those three rows buy is the **derivative-multiorder extraction and the solve**, ⛔ **not the
dimension**. ⚠ And it is an **identity, not an empirical miss**: with the field dimension left symbolic,
`[mu_R] − [rho_br]` is independent of it in all three axes. ⇒ ⛔ **no computation inside S9, in either
engine, can ever detect a wrong `[u]`.** It becomes detectable only when a consumer **imports** the entry
instead of re-declaring it — and S10 currently re-declares.

**5 · Three defects the rebuild introduced, which are open.**
- `wavevector_norm_dimension` **names the wrong object**: it denotes dim(|k|) = `[-1,0,0]`, and the value is
  `[-2,0,0]` = dim(k·k). ⭐ Right value, wrong name, in both engines and in the name-keyed file.
- **The placeholder-naming class is eight, not one**: five keys still carry the py-only `q`, and
  `dim_energy_density` / `dim_squared_velocity` have exact `.wl` counterparts under different names.
- `q_dimension` is **unpinned inside SymPy** — type it to anything and the engine exits 0 with a wrong
  value in the ledger. Two residuals detect it and are unguarded; the guard needs scoping because `X3`'s
  is legitimately nonzero under a form control. ⇒ a spec question.

**6 · What the export cannot express.** There is **no provenance field**. `wavevector_norm_dimension` shows
a consumer `PREMISE` with no indication the Wolfram engine **derives** it. ⭐ State the limitation; ⛔ do not
propose a schema.

**7 · ⛔⛔ NOTHING PERFORMS THE CROSS-ENGINE LOOKUP.** The naming convention exists so comparison is a join,
and **no artifact in the repository performs it.** ⚠ Re-pointing one line of an engine's name→object table
turns the light cone into an `omega^2` — wrong by `k²` — with the determinant, the root multiset and the
speed-dimension residual all **unmoved**, and **exit 0**. ⛔ The dimension check cannot catch it: `q_dimension`
is typed at exactly `−2L`. ⛔ The Wolfram engine has no `q` and would not move. ⇒ **the cross-engine
comparison on that row is the only instrument that would see it, and it does not exist.**

**8 · Correct what the rebuild made false**, if the existing record claims anything now untrue.

## What NOT to write

⛔ No claim that S9 derived the curl-only MacCullagh form — it is **posited**, and it is the step's live
falsifier. ⛔ No count of "independent agreements" larger than nine. ⛔ No statement that the export is
checked against the Wolfram engine. ⛔ No recommendation to build anything.

## Report — under 25 lines

What you added, what you corrected, and anything that read as incomplete.
