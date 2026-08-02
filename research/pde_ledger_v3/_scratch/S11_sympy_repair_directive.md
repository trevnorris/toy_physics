# DIRECTIVE — close a demonstrated hole in the S11 SymPy audit

Edit **only**:
- `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`
- `research/pde_ledger_v3/reduction/acceptance_check.py` (comment only — see R2)

⛔ Do not commit. ⛔ Do not create any other file. ⛔ Do not read
`research/pde_ledger_v3/mathematica/` — the independent second engine lives there and must stay
independent.

---

## R1 · ⭐⭐ THE HOLE — a registry relation's ALGEBRA is unchecked

**Demonstrated by an independent reviewer, not hypothesised.** Rewriting `R5` in `relations.yaml` from
`c_L − √(B_comp/ρ_br)` to `c_L − √(μ_R/ρ_br)` — i.e. silently declaring the longitudinal speed equal to
the transverse one, **the exact claim this step exists to settle** — leaves **every gate green**:

```
acceptance    exit=0   PHASE1_ACCEPTANCE: MATCH
dimgate       exit=0   DIMENSIONAL_HOMOGENEITY_GATE: PASS
able_to_fail  exit=0   ABLE_TO_FAIL_HARNESS: PASS
pytest        exit=0   11 passed
this script   exit=0   S11_VERDICT: PASS
```

Three guards interlock and all three miss it:
- **dimensional homogeneity** cannot see it — the two moduli share a dimension;
- **the acceptance fixture** cannot see it — the designated output stays a fresh variable, so the
  constraint count is unchanged;
- **the Mathematica engine** cannot see it *by design* — it must not read the registry.

⚠ **The committed `R5` is CORRECT.** This is a **missing control**, not a wrong value. ⛔ Do not "fix"
`relations.yaml`; it does not need fixing.

**Build the missing control.** Add an assertion to the audit that closes the gap between *what this
script derives* and *what the registry records*:

- take the parallel-kernel root **this script computes** in `A1`;
- load `R5` from the registry via `registry_read`;
- substitute the derived result into `R5`'s residual and **assert it vanishes symbolically**;
- print both the derived expression and the registry residual so a reader can see them side by side;
- register it in the enumerated assertion list so it participates in the verdict.

⭐ **Prove it is able to fail.** Ablate `R5`'s residual in a `/tmp` copy of the registry (or by an
in-process perturbation — your choice, but ⛔ leave the real `relations.yaml` untouched and restore
anything you perturb), confirm the new assertion **FAILs**, and report the output.

⚠ **Generalise it as far as it goes for free, and ⛔ no further.** `R4` (`c_γ`) has the identical shape
and this script has `c_γ`'s defining relation available to it. If covering `R4` costs a few lines, cover
it. ⛔ Do **not** build a framework, a schema extension, or a new gate file — the per-step audit is the
right home because it is the one place where the derivation and the registry record both exist.

## R2 · The acceptance fixture's provenance comment now UNDERSTATES

`reduction/acceptance_check.py`'s comment on `EXPECTED_MEDIUM_PAYLOAD` was rewritten to describe a
single deriving party — the builder, which also wrote the code that fixture polices. ⛔ A control whose
comment credits only the author of the thing it controls reads as weaker than it is.

**It now has three independent derivations**, all of which agreed with the committed numbers
`{12,7,5} / {12,7,5} / {12,6,6}`:
1. the orchestrator's, **pre-registered and committed before this script existed**;
2. a fresh review agent's, derived from `quantities.yaml`/`relations.yaml` alone;
3. a second independent reviewer's, by the same route.

⭐ Update the comment to record that accurately, keeping its existing standing instruction that on any
registry change the payload must be **recomputed and independently re-derived, never copied forward**.
⛔ Comment only — do not change any number.

---

## ⛔ ACCEPTANCE

- Every existing `S11_*` output line must be **unchanged**. This adds a control; it does not alter
  physics. ⛔ If you find yourself wanting to change a computed value, **stop and report it** rather than
  changing it.
- All five gates still pass on the unmutated tree.
- The new assertion appears in the printed enumerated list with its own PASS/FAIL and detail.
- The ablation you ran is reported, with the output showing the FAIL.

## Report

Under 25 lines: what the new assertion computes, the ablation and its output, whether you covered `R4`,
and confirmation that the pre-existing output lines are unchanged.
⛔ Do not summarise this directive back to me.
