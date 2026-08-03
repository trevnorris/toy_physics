# DIRECTIVE — S11b SymPy audit and registry insertion

**Deliverables (absolute paths):**
1. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_brane_bulk_interface_sympy_audit.py`
2. any additions to `/var/projects/toy_physics/research/pde_ledger_v3/reduction/quantities.yaml` and
   `/var/projects/toy_physics/research/pde_ledger_v3/reduction/relations.yaml`

Run the script and the existing gates. Iterate until everything exits 0. Then stop and exit — ⛔ do not
write a report or a summary document.

## ⛔ WHAT YOU MUST NOT READ

- `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11b_*` — a blind Mathematica audit of
  this same physics exists. It has been moved out of the tree; ⛔ do not retrieve it from git either
  (`git show`, `git cat-file`, `git checkout`). **The disagreement between the two engines is the test.**
- any file whose name contains `PREREGISTERED`.
- `/var/projects/toy_physics/research/pde_audit/` — **any file.** ⚠ That tree contains prior work on
  adjacent physics whose current validity is **not established**. The other engine is barred from it too;
  an asymmetry here would make one engine a transcriber and the other a deriver, and their agreement would
  then mean nothing.

⛔ `git status` and `git diff` are fine. `git show` and `git cat-file` are not.

## What this script does that the Mathematica one does not

The Mathematica audit is blind to the registry **by design**. This script is the **only** engine that sees
both the derivation and the recorded relations, so registry work belongs here and nowhere else.

1. **Import the registry** through
   `/var/projects/toy_physics/research/pde_ledger_v3/reduction/registry_read.py`. Do not re-implement its
   loading or validation.
2. **Add any new quantities and relations** the derivation earns. For each new quantity give a `kind`, a
   `counting_axis`, a dimension with provenance, and `aliases: []` unless an alias genuinely exists.
   ⚠ A quantity this step must **postulate** because a later step owes its derivation should be recorded
   as postulated, ⛔ never dressed as derived.
3. ⛔⛔ **EVERY NEW RELATION NEEDS ITS OWN ALGEBRA ASSERTION, because no gate checks relation algebra.**
   Substitute the root/expression **this script derived** into the registry residual that records it, and
   assert it vanishes. Then **ablate** that assertion — corrupt the registry residual and confirm the
   assertion fails — and report the ablation's literal output.
   ⚠ **This is measured, not hypothetical.** A relation on the previous step was rewritten to assert the
   exact claim that step existed to settle, and **all five gates stayed green**: dimensional homogeneity
   was blind because the two moduli share a dimension, the acceptance fixture was blind because the
   designated output stays a fresh variable, and the Mathematica engine cannot help by design.
4. ⛔ **A check of the form `c := sqrt(X)` followed by asserting `c**2 - X == 0` is identically zero by
   construction and proves nothing.** If you write one, label it as definitional in its own output line.
5. **Run all gates** and report each one's literal output: the acceptance check, the dimensional
   homogeneity gate, the able-to-fail check, and pytest. ⚠ If the continuous payload changes, that is a
   result to report, ⛔ not something to make match by adjusting a fixture.

## Script conventions

- ⛔ **There is no shared dimensions module in this ledger** — `reduction/` provides `registry_read.py`
  and the gates, nothing else. ⚠ **W9's dimensions must be DERIVED from the specification's equations**;
  ⛔ do **not** read them from `quantities.yaml`, from the homogeneity gate, or from any other table. The
  registry is what W9 is independently checking, so reading it makes the check vacuous.
- Output tag per line, `TAG: value`, with no `WL_` prefix.
- Keep the total runtime under **10 minutes**; runners get `timeout 600`.

---

