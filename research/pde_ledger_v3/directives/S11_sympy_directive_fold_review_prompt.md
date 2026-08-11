# Independent physics review — the FOLDED S11 SymPy build directive

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_no_ledger_build_directive.md`
— **81 lines**, Codex-authored, a **fold** of two legs' findings into its own 66-line predecessor.

⚠⚠ **A builder launches from this.** ⭐ It is checked here or nowhere.
⚠⚠ **A revision that breeds defects in the material it just changed is this project's measured failure
mode** (`CLAUDE.md` rule 15). ⇒ ⭐ **the newly written text deserves more scrutiny than the text that
survived, not less.**

## Read in this order — ⛔ load-bearing

1. `/var/projects/toy_physics/CLAUDE.md`.
2. `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **the binding physics**.
   ⭐ Form your own view of what an engine must do to satisfy it. ⚠ The working-tree copy of that spec is
   **being amended in a parallel round**; ⛔ review the directive against the **committed** baseline, and
   ⛔ do not review the amendment here.
3. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_exports.py` — the only import.
4. `directives/_legs/round5_grok_build_directive.txt` and `directives/_legs/round5_claude_build_directive.md`
   — the findings this fold answers. ⚠ ⛔ **`B2` is RETRACTED by the leg that raised it.**
5. ⭐ **Only now** open the directive.

## Do not read

- The other leg's output.
- ⛔ The two S11 engines being **replaced** (`scripts/S11_stray_longitudinal_sympy_audit.py`,
  `mathematica/S11_stray_longitudinal_mathematica_audit.wl`). ⚠ They predate the spec and carry the defect
  class the rebuild removes.

## ⭐ Scope

⛔⛔ **The engine writes no ledger** — it imports `S10_exports.LEDGER` read-only for `§Q6r` and emits a
flushed stdout tag stream. ⛔ Storage, export, merge, naming-namespace and identity are **out of scope**
(`DEFECT_REGISTER.md#c20`); ⭐ a finding that the directive **strays into** them is in scope.
⛔ Four shared tag names, the decomposition rule and boolean rendering were **deliberately removed** from
this directive and are being amended into the shared spec. ⭐ Check the removal is clean; ⛔ do not ask this
file to carry them.

## What to check

### ⭐⭐ 1 · DID THE FOLD FIX WHAT IT CLAIMS, AND CAN EACH FIX FAIL?

⛔ **Construct the failure for each; do not reason about it.** For every one, report literal output.

- **Omission may never respond to a payload's content.** ⭐ The predecessor's wording put every cell in
  `SKIPPED_PAIRS` at exit 0, because `§Q11`'s `C1`–`C4` have no pinned undecided token and exist precisely
  to measure a deficiency. ⛔ Under the new wording, what happens to `C1`–`C4` when the supplied material
  does not determine them? ⭐ Walk it.
- **Scope covering every required emitted object.** ⛔ Is `PD_TERM` now reachable, **and** is it required to
  be **read out of** the computed `V6` object rather than typed? ⚠ `§7:996` exists so that corrupting the
  census **moves the package's action**. ⭐ Is any other `§7`/`§9` object still outside the scope sentence?
- **An evaluable completion condition.** ⚠ Measured on the predecessor: two implementations both satisfied
  the wording and one **certified an incomplete cell** (`24/27` tags, cell in `RUN_PAIRS`). ⛔ Build both
  against the new wording. Does the new condition separate them, and **where does its required set come
  from** — is that source independent of the emitter it checks?
- **`Q6r` resolution keyed on the lookup.** ⚠ Measured on the predecessor: stripping one metadata field
  deleted `mu_R`'s comparison — the only place this build compares its own dimension solve against the
  predecessor's. ⛔ Re-run that form control against the new text. ⭐ Are the three distinct lookup-failure
  shapes still distinguishable?
- **`M_ROUTE_USED` sourcing.** ⚠ `§Q2:344-346` requires it read from the **same route-selection object that
  supplies the matrix consumed downstream**, ⛔ not a literal beside it, and the corollary-5 exemption list
  is closed at four entries. ⛔ Can a typed string still satisfy the new wording?
- **The elapsed-time cap.** ⭐ Confirm it is gone and that **nothing replaced it**. ⚠ `§7:1041` governs by
  observable progress. ⛔ Is there any surviving sentence under which a builder could shed a cell for cost?

### ⭐⭐ 2 · WHAT DID THE FOLD BREAK?

⭐ Diff the 81 lines against the 66-line predecessor and read **every new or changed line** as if it were
new work.

- ⛔ Does any new sentence contradict another sentence in the same file?
- ⛔ Does any new sentence contradict the shared spec? ⚠ The directive says the spec wins every conflict —
  ⭐ so a contradiction is a **latent** defect, not a harmless one: it is a place two builders diverge.
- ⛔⛔ **Does any new sentence state what a computation will produce** — a value, count, membership, sign or
  outcome? ⚠ Rule 5. ⭐ Both legs scanned the predecessor and found it clean; ⛔ check the new text did not
  break that.
- ⛔ Does any new sentence **specify a recipe** where it should name an object (rule 3)?

### ⭐⭐ 3 · CAN THE ENGINE THIS DIRECTIVE DESCRIBES STATE A CONCLUSION?

⛔ The measured defect the rebuild exists to remove: a physics conclusion emitted as a **typed sentence with
no CAS object behind it**, where deleting the computation left stdout **byte-identical**. ⚠ It survived
**eight** review legs.

- ⭐ Is there any payload the directive permits to be assembled from a **literal beside the computation**
  rather than read from the live object?
- ⭐ Are both operands **and** the residual required before every guard?

### ⭐ 4 · THE REMOVAL

⭐ Four names, the decomposition rule and the rendering pin were removed.

- ⛔ Is there a **dangling reference** to any of them?
- ⛔ Does the directive still require something the removed text was the only support for?
- ⭐ Does it now leave a builder able to emit an object under **no** name at all?

### ⭐ 5 · WHAT IS MISSING

⛔ The most expensive defects here have been **absent computations**. ⭐ The spec orders residuals whose two
operands come from **genuinely different routes**. ⛔ Would an engine built from this directive emit every
one, and would a nonzero one be visible? ⚠ Name any that could be silently skipped.

## ⛔ Constraints on how you run

- ⛔ Read-only on the working tree; copy to `/tmp` to modify. ⛔ `timeout 600` on every run, ⛔ never raised.
  ⚠ ⭐ A run that exceeds it is a **failed probe** — report it and move on; ⛔ do not conclude the engine is
  infeasible from a slow route. ⚠ A previous leg made exactly that over-call and retracted it.
- ⛔ **No Mathematica.**
- ⭐ Save every script and its literal stdout to named absolute paths and report them. ⛔ A prose derivation
  will be discarded.

## Report

For each finding: **what is wrong**, **the mechanism by which a wrong result would survive**, **the literal
output that shows it**, and **what must become true instead** — ⛔ not a rewrite.
⭐ Separate blocking from non-blocking. ⛔ "Nothing found" with no ablation behind it is not a result.
