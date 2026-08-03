# DIRECTIVE — S11b-B SymPy audit and registry insertion

**Deliverables (absolute paths):**
1. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bB_interface_assembly_sympy_audit.py`
2. any additions to `/var/projects/toy_physics/research/pde_ledger_v3/reduction/quantities.yaml` and
   `relations.yaml`

Run the script and the existing gates. Iterate until everything exits 0. Then stop and exit — ⛔ do not
write a report or a summary.

## ⛔⛔ ONE OF TWO INDEPENDENT ENGINES

A blind Mathematica audit of the same physics exists. **The disagreement between the engines is the test**,
so an agreement reached through a shared source certifies nothing.
⇒ **The read-bar list is §0b of the shared specification below**, byte-identical to the other engine's.

## Registry work — this engine only, and ⛔ AFTER the physics

⚠ **Do the physics tasks FIRST, from the specification alone.** ⛔ Do not consult `reduction/` while
deriving. ⚠ **Measured hazard:** registry access lets one engine silently identify two symbols the other is
treating as independent, after which their "agreement" certifies an identification only one of them made.
⭐ §1 fixes every symbol's name and meaning explicitly, so there is nothing to infer — ⛔ do not infer
anyway.

**Then**, and only then:
1. Import the registry through `reduction/registry_read.py`. Do not re-implement its loading or validation.
2. Add the quantities and relations this step earns. Give each a `kind`, a `counting_axis`, a dimension
   with provenance, and `aliases: []` unless an alias genuinely exists. ⚠ A quantity this step must
   **postulate** because a later step owes its derivation is recorded as **postulated**, ⛔ never dressed
   as derived.
3. ⛔⛔ **EVERY NEW RELATION NEEDS ITS OWN ALGEBRA ASSERTION, because no gate checks relation algebra.**
   Substitute the expression **this script derived** into the registry residual that records it and assert
   it vanishes. Then **ablate** it — corrupt the residual, confirm the assertion fails — and report the
   literal output. ⚠ Measured: a relation was once rewritten to assert the exact claim its step existed to
   settle and **all five gates stayed green**.
4. ⛔ A check of the form `c := sqrt(X)` then asserting `c**2 - X == 0` is identically zero by
   construction. If you write one, label it definitional in its own output line.
5. **Run all gates** and report each one's literal output: acceptance, dimensional homogeneity,
   able-to-fail, pytest. ⚠ If the continuous payload changes, that is a **result to report**, ⛔ not
   something to make match by adjusting a fixture.

## Script conventions

- ⚠ **B7's dimensions must be DERIVED from the specification's equations**, ⛔ not read from
  `quantities.yaml` or the homogeneity gate. The registry is what B7 independently checks.
- One tag per line, `TAG: value`, no `WL_` prefix.
- Keep total runtime under **10 minutes**; runners get `timeout 600`.
- ⚠ **On a non-uniform background do not silently assume plane waves are eigenmodes.** Where you must
  restrict to a tractable form, ⛔ say so in the emitted value rather than in a comment.

---

