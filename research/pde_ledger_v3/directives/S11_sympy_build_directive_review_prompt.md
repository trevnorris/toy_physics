# Independent physics review — the S11 SymPy build directive

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_no_ledger_build_directive.md`
— 66 lines, **Codex-authored**.

⚠⚠ **A builder launches from this.** ⭐ It is checked here or nowhere. ⛔ Do not assume it is right because
it is short — four earlier artifacts covering adjacent material were each blocked by two legs.

## Read in this order — ⛔ load-bearing

1. `/var/projects/toy_physics/CLAUDE.md`.
2. `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **the binding physics**,
   1149 lines. ⭐ Form your own view of what an engine must do to satisfy it.
3. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_exports.py` — the only import.
4. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — the
   nearest built precedent for this pattern.
5. ⭐ **Only now** open the directive.

## Do not read

- The other leg's output.
- ⛔ `mathematica/S11_stray_longitudinal_mathematica_audit.wl` and `scripts/S11_stray_longitudinal_sympy_audit.py`
  — the engines being **replaced**. ⚠ They predate the spec and carry a defect class the rebuild removes.

## ⭐ Scope — what this build is, and is not

⛔⛔ **THE ENGINE WRITES NO LEDGER.** ⭐ It imports `S10_exports.LEDGER` read-only for `§Q6r` and emits a
flushed stdout tag stream. ⛔ There is no export, no merge, no storage key.
⛔ **Storage, naming-namespace, identity and chain questions are OUT OF SCOPE** — four designs for them were
blocked and the record is `DEFECT_REGISTER.md#c20`. ⛔ Do not reopen them here; ⭐ a finding that the
directive **strays into** them is in scope.

## What to check

### ⭐⭐ 1 · CROSS-ENGINE PARITY OF NAMES THE SPEC DOES NOT SUPPLY

`§8` requires both engines to produce **parallel** tag sets — strip `PY_`/`WL_` and the suffix is the same.
`§8:1084` also permits a builder to name an object the spec does not name, reporting it under `§10`.

⭐ **Measure this.** Which objects the directive requires does the shared spec leave **unnamed**? For each:

- Is it a **shared** tag or an engine-local one? ⭐ Only shared ones carry a parity obligation.
- ⛔ If a shared tag's name is fixed **in this directive**, what makes the Wolfram engine — which will be
  built from a **different** directive and never sees this one — produce the same suffix?
- ⭐ What does the comparator see if the two engines choose differently? Is that a **physics** disagreement
  or a **vocabulary** one, and can a reader of the comparator output tell them apart?
- ⚠ **Measured precedent:** the previous pair of S11 engines were built from two directives that decomposed
  the work differently, and **shared exactly one tag suffix between them**. ⭐ Establish whether this
  directive can reproduce that outcome, or cannot.
- ⭐ If you find a real exposure, ⭐ **say at which level it must be fixed** — the shared spec, this
  directive, or the Wolfram directive — ⛔ and do not fix it yourself.

### ⭐⭐ 2 · DOES THE DIRECTIVE COVER WHAT THE SPEC REQUIRES?

⭐ Walk `§6`'s `Q1`–`Q11` and `§7`'s packages against the directive.

- ⛔ Name any object the spec **requires** that a builder following this directive would not emit, or would
  emit at the wrong scope.
- ⛔ Name anything the directive requires that the spec **does not**, and say whether it adds physics — ⚠ the
  directive may not add physics.
- ⭐ `§5`'s locus protocol, `§Q8`'s component-scoped counts and their status tokens, and `§Q9`'s census
  feeding `P_D` are the most intricate obligations. ⭐ Is each reachable from what the directive says?
- ⚠ **`Q9` precedes action assembly** because `P_D` is built from its output. ⭐ Does the directive's ordering
  actually make that possible for every package that needs it?

### ⭐⭐ 3 · RULE 2 — ⛔ CAN THE ENGINE THIS DIRECTIVE DESCRIBES STATE A CONCLUSION?

⛔ The measured defect this whole rebuild exists to remove: an engine emitting a physics conclusion as a
**typed sentence with no CAS object behind it**, where deleting the computation left stdout
**byte-identical**. ⚠ It survived **eight** review legs.

- ⭐ Does anything in the directive permit a payload assembled from a **literal beside the computation**
  rather than read from the live object?
- ⭐ Are both operands **and** the residual required before every guard, or is there a place a bare residual
  could be emitted?
- ⭐ The directive requires symbolic tests to stay **CAS objects rather than Python booleans**. ⛔ Is that
  achievable for every test the spec orders, and what happens where it is not?
- ⭐ `RUN_PAIRS` must be **observed**, not declared. Is the directive's rule strong enough that a declared
  sweep cannot be emitted as an observed one?

### ⭐ 4 · Q6r AND THE IMPORT

- ⭐ Take the spec's coefficient name map and perform the **live** two-step lookup against the committed
  `S10_exports.py`. ⛔ Does the directive's resolved/unresolved rule handle every case the real artifact
  presents? Report what you ran.
- ⛔ Does the directive state, anywhere, **which** coefficients resolve? ⚠ That is a value a builder must
  compute, and stating it is a rule-5 leak.

### ⭐ 5 · THE RULES

- **Rule 5** — ⛔ no measured physics value, count, membership, sign, spectrum or expected outcome. ⭐ Scan
  every line.
- **Rule 3** — name the object, ⛔ not the recipe.
- **Rule 11** — ⛔ cost is never a reason to drop a control.
- ⭐⭐ **PHYSICS, ⛔ NOT PROCESS** (standing user instruction). ⛔ If an item is ceremony, that is a finding.

### ⭐ 6 · WHAT IS MISSING

⛔ The most expensive defects here have been **absent computations**. ⭐ The spec already orders residuals
whose two operands come from **genuinely different routes**. ⛔ Would an engine built from this directive
emit every one of them, and would a nonzero one be visible? ⚠ Name any that could be silently skipped.

## ⛔ Constraints on how you run

- ⛔ Read-only on the working tree; copy to `/tmp` to modify.
- ⛔ `timeout 600` on every run, ⛔ never raised. ⛔ **No Mathematica this round** — the licence has two seats
  and this artifact is a SymPy directive.
- ⭐ Save every script and its literal stdout to named absolute paths and report them. ⛔ A prose derivation
  will be discarded.

## Report

For each finding: **what is wrong**, **the mechanism by which a wrong result would survive**, **the literal
output that shows it**, and **what must become true instead** — ⛔ not a rewrite.
⭐ Separate blocking from non-blocking. ⛔ "Nothing found" with no ablation behind it is not a result.
