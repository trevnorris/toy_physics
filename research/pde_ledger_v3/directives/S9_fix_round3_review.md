# Independent review — S9 export chain, fix round 3

## Artifact

Uncommitted changes in `/var/projects/toy_physics` on top of commit `81243f1d`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
research/pde_ledger_v3/scripts/S9_exports.py                                (regenerated)
research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out  (regenerated)
```

Pre-round state: `git show 81243f1d:<path>`. `git diff` is exactly this round. The extractor and the
Wolfram engine were not changed. Directive: `directives/S9_export_chain_fix_round3.md`. Plan:
`S9_REWRITE_PLAN.md`. **Read the directives last.**

⛔ **F-2 (the standard-name→object binding) is with the user, out of scope, and must not be re-reported.**

## Context

This round was meant to be small: reclassify `root_scaling_quadratic` to `PREMISE`, make the
`dim_*` reference chain's classes consistent under one criterion, apply the entailment rule to all twelve
cross-engine rows, and enumerate the no-baseline tags mechanically. It reports an independent agreement
count of **7** and an enumeration of **22**.

## ⛔⛔ Attack this first — the orchestrator believes this round reintroduced a defect it had removed

`S9_light_requires_shear_sympy_audit.py:810`:

```python
class_tag = "PREMISE" if suffix == "ROOT_SCALING_QUADRATIC" else output_class_for(value)
```

and `:478`, inside the function `derive(` which begins at `:445`:

```python
quadratic_scaled_roots = [...]  # PREMISE · supplied quadratic scaling reference
```

**The claim to test:** the extractor scans module scope only, so that annotation is decorative and does
nothing; the class is in fact a hardcoded string comparison against a **tag name**. If so this is precisely
the defect rounds 1 and 2 removed — a class as a typed literal rather than read from a declaration — except
worse, because the annotation makes it look read.

**Establish it or refute it, with mutations:** rename the tag and see whether the class silently reverts;
change or delete the `:478` annotation and see whether anything moves; check whether the extractor reports
any finding in either case.

**Then the general question, which is the one that matters.** A leg already measured that **64 of 139
exported rows are built inside function bodies** where no annotation can be scanned. `root_scaling_quadratic`
is the first one whose correct class is known to differ from the default. So:

- Is the honest state of affairs that **the class system cannot express a class for objects constructed
  inside functions**, and a special-case merely hides that for one of them?
- Which is less dishonest in the export the knit phase will read: a **correct** class by a fragile
  name-keyed literal, or a **wrong** class (`DERIVED`) with the limitation stated plainly?
- Is there a fix that is **not machinery** — for instance carrying the class at the point where the object
  is named into the exported collection, so a rename moves both together?

⚠ **Do not implement anything.** Report the measurement and your judgement. This may be a decision for the
user alongside F-2, and saying so is a legitimate answer.

## Also attack

**① The independent count of 7.** Verify each of the five exclusions under the round's three rules —
*entailed* (algebraic consequence of counted rows), *unpublishable* (guarded, non-zero never reaches the
output), *counted*. Then attack the seven: is any of them also entailed, or unable to disagree for a reason
not yet named? Apply the rule in both directions — a row wrongly excluded is as bad as one wrongly counted.

**② The `dim_*` class chain.** The round says a supplied reference is `PREMISE` and a quantity computed on
the way to one stays `DERIVED`, and reports applying it to seven lines where the directive said five. Read
`:245`–`:251` and judge whether the criterion is applied consistently, and whether the criterion itself
separates the cases it claims to.

**③ The 22-tag enumeration.** It has been wrong three times (12, 16, 17). Recompute it mechanically from the
two tag sets and confirm both the count and the membership.

**④ Regression.** No exported value, dim, display or step moved against `81243f1d`; the only class move is
`root_scaling_quadratic`. Confirm. The Wolfram engine must still re-run byte-identical. A **form** ablation
must still move exported values — the round reports 86 on the gradient form; confirm the form and the count.

## Required method

Derive independently in CAS where a value is in question; a prose re-derivation is worth nothing. Write your
script before opening the artifact and save it and its literal stdout to named absolute paths. Ablate copies
in a scratch directory — **never the working tree, and never `git checkout`/`restore`/`stash`: the artifact
is uncommitted.** For every claim ask *which line computed this?*

## Physics filter

> Report a finding only if it catches a way the physics could be wrong, or changes what may be claimed.

**This project's failure mode is adding machinery, not omitting it.** A finding that a piece of this is
unnecessary is worth as much as one that it is insufficient. Propose a new check only if its absence lets a
specific wrong answer through — and name it.

A leg returning "nothing survives" is weak evidence alone: say what you checked, what you could not, and
what would have had to be true to find something.

## Operational

One Mathematica kernel at a time (two seats); `timeout 600` per run. Do not modify the working tree. Do not
commit. Do not restore `reduction/`. Absolute paths. Save every ablation script and its stdout; report paths.

## Output

Rank most-severe first: **claim · evidence (literal command and output, or `file:line`) · what must change.**
