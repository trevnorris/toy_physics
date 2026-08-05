# S10 Mathematica engine — repair round 2

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
**Its specification, ⚠ AMENDED since you last saw it:**
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md` — ⭐ **re-read §3, Q5
and Q6 in full**; three sections changed.

⛔ **Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file.**
⛔ **Do not read** `research/pde_ledger_v3/steps/`, `research/pde_ledger_v3/paper/`, or the sibling SymPy
engine.

⚠ After every fix, re-run and confirm it completes under **10 minutes** and exits `0`.

---

## ⭐ SPEC AMENDMENTS — apply these first; two of them were MY errors, ⛔ not yours

**S1 · §3 — the period average is now over the PHASE.** `⟨F⟩ ≡ (1/2π)∫₀^{2π} F dφ` with `φ` the phase.
⛔ The old `t`-limits `(0, 2π/ω)` are a real period **only if `ω` is real and nonzero**, and Q3 exists
partly to find roots where `ω²` is zero or negative. ⭐ Update the averaging **and** any premise tag that
still states the old convention.

**S2 · §1 — `[u] = (1,0,0)` is now a STATED premise.** ⭐ Emit it as a premise tag. ⛔ Do not infer it.

**S3 · Q5 — emit the ENTIRE Q5 tag set for every root, unconditionally.** ⚠ The old text contradicted
corollary 4 and you followed the old text. ⭐ Where the ratio is undefined, the payload is an **explicit
marker object**; ⛔ the tag never disappears.

**S4 · Q6d — new.** Emit the number of independent dimension equations, the number of unknown coefficient
dimensions, their difference, and whether the system was over-/exactly-/under-determined. ⭐ Emit a tag
recording that the action's homogeneity booleans are **vacuous** when that difference is `≤ 0`.
⭐ Give homogeneity tags for **solved action terms** and for **everything else** distinguishable names.

---

## ⛔⛔ R1 — BLOCKING: Q7 does not re-enter at the action, so it tests nothing about the package

`runCurlComparison` (≈1125) takes only `prefix` and `dimension`. It **re-types** `S_curl` from fresh
symbols instead of substituting `∂_i u_j → g_ij` into the action's own `curlStiffness`.

⚠ **Measured:** `Q7_STIFFNESS` / `Q7_CURL_NORM` / `Q7_RESIDUAL` are **byte-identical across all six `D=3`
packages** and moved under **none of seven ablations** — including a sign flip inside the action's own curl
term. ⇒ ⭐ **Delete the entire action and those tags do not change.** That is §4's structural rule broken in
its "second source of truth" form.

⭐ **Fix:** build Q7's first operand by **substituting independent gradient symbols into the package's own
action-level stiffness**, ⛔ not by re-typing the formula. ⭐ That satisfies Q7's "independent symbols"
requirement **and** restores the data dependency.
⚠ Q7 is defined at `D=3` on the **curl** stiffness by name, so its payload may legitimately repeat across
packages that share that stiffness — ⭐ but it must **move** when the action's curl term is corrupted.

## ⛔⛔ R2 — BLOCKING: Q8b's stratum point over-pins the parameters

`pointVariables` (≈1383–1428) includes the amplitudes, `rhoBr`, `muR` and the control parameter, so
`FindInstance` fixes **those** to numbers as well as `k`.

⚠ The Q3/Q4 re-run then happens at **one arbitrary numeric parameter point**, ⛔ not at a generic point of
the stratum. ⇒ ⭐ **a rank drop caused by an accidental relation among the chosen coefficient values is
indistinguishable from a property of the stratum** — the exact generic-versus-special confusion Q8 exists
to remove.

⭐ **Fix:** solve for a point in the **wavevector components only**. ⭐ Keep `ρ_br`, `μ_R` and every control
parameter **symbolic** through the stratum re-run. ⚠ If a symbolic re-run is intractable, ⭐ specialise
**as few** parameters as possible and **emit exactly which ones were fixed and to what** as its own tag.

## ⛔ R3 — a suppressed message is hiding an incomplete solve

Line ≈56 quiets `Solve::svars` — **exactly the message that says the solve did not solve for all named
variables** — and every `_LOCUS` tag in Q3 and Q8a passes through that path.

⭐ **Fix:** ⛔ stop suppressing it. ⭐ Capture whether it fired, per solve, and **emit that as a tag beside
the locus it qualifies**. ⚠ An incomplete or parametric locus emitted with nothing recording its
incompleteness is a locus that will be read as complete.

## ⛔ R4 — emission still conditional on a payload's value, and a name that states its outcome

- `_Q8_POINT_SEARCH_UNLOCATED` versus `_Q8_SOURCE`/`_Q8_POINT` (≈1401–1412) switch on whether
  `FindInstance` succeeded. ⭐ Emit a **fixed** tag set; put the outcome in the **payload**.
- ⛔ The name `…_UNLOCATED` **states the outcome** — corollary 2 forbids that in a name. ⭐ Rename to name
  the **object** (the point search), ⛔ not its result.
- Root-pair coincidence tags vanish when a package has fewer than two roots (≈462, ≈1230). ⭐ Emit them
  with an explicit "no pairs" payload.

## ⛔ R5 — `_ANSATZ_FREQUENCY_BRANCH` needs the `_LOCAL_` infix

It is engine-specific (§8), so a parity checker currently treats it as a cross-engine quantity. ⭐ Rename
per §8. ⚠ And with S1 applied, re-check whether the branch is still taken at all.

## ⛔ R6 — stratum loci are solved over an empty variable set

`realLocusSolve[equation, {}]` / `reduceAllowed[equation, {}, …]` (≈1239–1240): at a fully numeric point
there is nothing to solve for, so those `_LOCUS` and `_ALLOWED_INTERSECTION` tags are **truth values
wearing locus names**. ⭐ With R2 applied the point is no longer fully numeric — ⭐ re-check these and give
them names matching what they actually are.

## ⛔ R7 — two payloads are self-descriptions rather than computed objects

`"Disposition" -> notReimposed / noConditionDropped` (≈431–433) asserts what the script did;
`_Q3_ROOT_LIST_STATUS` (≈447–451) restates `_Q3_ROOT_COUNT`. ⭐ Emit the **objects** — the dropped
condition itself, and the two counts — ⛔ not a sentence about them.

## ⛔ R8 — `Q8_ALLOWED_STRATA` de-duplicates on PRINTED FORM

`GatherBy` on the printed `Reduce` output (≈1358–1366) merges strata by **string equality** ⇒ two
syntactically different descriptions of the same region are counted twice, and two distinct strata that
print alike collapse to one point. ⭐ **Fix:** group by a **structural** comparison. ⚠ This makes the
stratum count a mathematical quantity rather than a syntactic one.

## ⛔ R9 — `emittedTagCount` is maintained and never emitted

(≈6, 12.) ⭐ Emit it, so the tag total is a computed output rather than something a reader must recount.

---

## ⭐ What NOT to change — all confirmed live by ablation on two independent legs

⛔ `N3` (stacked-rank transverse count), `N7` (dual-algorithm count residual), the separate construction of
`M_A` and `M_B`, per-package re-entry at the action, Q5's scaling method, the dimension **tree walk**, and
the `UnhandledHeads` machinery. ⚠ In particular `N7`'s residual held while **both its operands moved** —
⭐ that is a working two-algorithm check; ⛔ do not "simplify" it.

## Report back — ⛔ under 25 lines

1. One line per `S1`–`S4` and `R1`–`R9`: fixed / partially / not, with line numbers.
2. New tag count, wall-clock, exit code.
3. ⛔ **Do not report what any value came out to be**, and ⛔ do not say whether anything "worked".
4. ⭐ Anything in this list you believe is **wrong**, or any fix that would break something the
   "do not change" list protects. ⭐ This is wanted.
