# S10 Mathematica engine — repair round 3 (final)

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
**Spec, ⚠ AMENDED:** `directives/S10_SHARED_PHYSICS.md` — ⭐ **re-read §5 corollary 5, Q6d, Q7, §7**.

⛔ Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file.
⛔ Do not read `steps/`, `paper/`, or the sibling SymPy engine.
⚠ Re-run after every fix: under **10 minutes**, exit **0**, and ⛔ **leave no `WolframKernel` behind.**

⭐ **A review leg ablated this engine and confirmed the whole computed chain is live** — N3, N7, both
matrix routes, per-package action re-entry, the dimension tree walk, Q6d's able-to-fail behaviour.
⛔ **Do not change any of those.** Every item below is a **declaration that is not wired to what it
declares**, or a spec amendment.

## ⛔⛔ W1 — the premise tags are second literals, not reports (corollary 5)

⚠ **Measured:** an ablation changed the analyzer's actual `[u]` and moved **281 dimension payloads**, while
`PREMISE_FIELD_DIMENSION` **kept printing its old value**. Same defect in `PREMISE_PERIOD_AVERAGE` (it
re-types the window and normalisation the code applies elsewhere) and `ANSATZ_FREQUENCY_BRANCH`.

⭐ **Fix: derive each of these payloads FROM the live object the computation consumes** — read the value
back out of the analyzer / the averaging routine / the ansatz, ⛔ never from a literal beside it.
⭐⭐ **Acceptance you must satisfy: perturb the thing each tag declares, and the tag MUST move.**

## ⛔ W2 — `Q2_DOWNSTREAM_ROUTE` is authored, not observed

It would print the same name if the other matrix were used downstream. ⭐ **Fix:** derive it from the
object actually passed to `runSpectrum`, so corrupting one route changes it.

## ⛔ W3 — `Q8_PARAMETER_SPECIALIZATIONS` can never be non-empty

Its `Cases` filter selects from a variable set **disjoint** from what `point` contains, so all instances
are `{}`. ⚠ This is a side-effect of the previous round's fix, which correctly made the parameters stay
symbolic. ⭐ **Fix:** make it report what was **actually** specialised — which, if nothing was, is an
explicit "no parameters specialised" payload. ⛔ An always-empty tag is indistinguishable from *never
computed*.

## ⛔ W4 — the run record counts loop iterations, not emissions

A package that raised and truncated would still be appended and reported as run. ⭐ **Fix per the amended
§7:** append a pair **only after that package has finished emitting**, and emit the record **after** the
sweep.

## ⛔ W5 — `lambdaScale` is never declared positive

Q5 calls it a positive scalar; it appears in no `JOINT_ASSUMPTIONS` payload. ⚠ Nothing is wrong in the
current output, but a root carrying `Sqrt`/`Abs` could leave the ratio unreduced and print a false
not-a-pure-power. ⭐ **Fix:** declare it positive and carry it in the joint assumptions.

## ⭐ W6 — Q7 now has a spec ruling; conform to it

The amended Q7 requires the comparison to use **the package's own stiffness**, which is what you already
do. ⭐ **Fix only the naming** so the tags say so, and ⚠ note that a nonzero residual for a non-curl
package is **expected** — ⛔ do not "fix" that.

## Report — ⛔ under 15 lines
1. One line per `W1`–`W6`: fixed / partially / not, with line numbers.
2. For **W1 and W2**, state the perturbation you ran and that the tag moved.
3. New tag count, wall-clock, exit code, and confirmation no kernel survived.
4. ⛔ Do not report what any value came out to be. ⭐ Anything here you believe is wrong is wanted.
