# S11c-a T7 control-battery characterization (task 2) — result

The comparator's §5 control battery is VALID: the two HOLD controls hold across engines (0 genuine nonzero),
and the two BITE controls bite for every object (no dead control). No finding hides in the controls.

## HOLD controls — should be 0 across engines (route-invariance N4/N6, and S11c-a→S11b regression)
Command: `python3 S11c_a_control_hold_classify.py` (structural zero-test incl. tuples/matrices + collapse
X(args)->X + integral linearity). See `S11c_a_control_hold_classify.stdout`. Result:
- REP_INVARIANCE_RESIDUAL (88): ZERO=48, STRUCTURAL=16 (PY tuple vs WL scalar), UNJOINED=24 (one-engine-only). **NONZERO=0.**
- UNIFORM_LIMIT_RESIDUAL (412): ZERO=148, STRUCTURAL=24, UNJOINED=240 (uniform-limit pairing is partial — a
  coverage asymmetry, NOT a physics finding; RUN_ACCOUNTING families_with_unpaired=15 already flags it). **NONZERO=0.**
=> route-invariance and the uniform-limit regression hold across engines wherever the two engines join.

## BITE controls — should bite per-engine (one-sided corruption; C-1 source form ablation)
Command (per-object liveness on operand_A = PY's own base-minus-ablated / base-minus-corrupted; a real
expression = the control moved PY's operand):
```
awk '/^CASE family=CONTROL_FORM_RESIDUAL /{match($0,/OBJECT=[A-Z_0-9]+/);obj=...;g=1;next}
     g&&/^operand_A = /{ if($0 ~ /operand_A = \(?Integer\(0\)\)?(, Integer\(0\))*\)?$/) null[obj]++; else bit[obj]++; g=0 }' comparator_run.out
```
Result — EVERY object ALIVE (bit>=1):
- CONTROL_FORM (18 objects): all bit>=12; objects with `null` (FACE_SHIFT 240/480, FACE_VELOCITY 12/12,
  PROJECTION_* ~30/18, VIRTUAL_CONSTRAINT 54/18) are legitimately null where the ablated direction's W_bg jet
  does not couple to that object — each still bites elsewhere. No dead control.
- CONTROL_INDEPENDENCE (6 objects): all bit>=2 (CLOSURE 8/8, EVOLUTION 4/4, RELATIVE_FLUX 20/4, TRACTION 8/8,
  VIRTUAL_CONSTRAINT 2/6, VIRTUAL_WORK 16/0). Every object alive.

Cross-engine A_minus_B of the BITE controls is EXPECTED nonzero (each engine bites in its own representation)
and is not finding-relevant: the bite controls are form-ablations/corruptions of the SAME shape-derivative
objects whose primary cross-engine residuals were adjudicated benign (§25) — no new residual class appears.

---

## ⚠ CORRECTION (2026-08-27, post-`cccb4f9e`) — the BITE-liveness above was WRONG

The "every object bites (no dead control)" claim above came from a string awk that **counted sentinel
operands** (`TextAtom('UNDEFINED_UNJOINED')`, `Mismatch`, `<MISSING>`) as bites — a buggy liveness test.
A correct liveness test (a genuine bite = `operand_A` is a nonzero sympy expression carrying a `Symbol`;
sentinels and literal zeros excluded) run against the FINAL committed-state comparator
(`~/.s11_build/comparator_final_cccb4f9e.out`), script
`scratchpad/bite_liveness.py`, stdout `~/.s11_build/bite_liveness_cccb4f9e.out`:

- **CONTROL_FORM (∂W_bg jet ablation) bites 16 of 18 objects** — including FACE_VELOCITY (12),
  PROJECTION_{SHAPE,DYNAMIC,RESIDUAL,TERM_ORIGINS} (6 each), EVOLUTION_*, VIRTUAL_CONSTRAINT, and the
  geometry set (24–48 each). It is **DEAD (BITE=0)** for exactly **FACE_SHIFT** (0/480) and
  **PROJECTION_STATIC_OPERAND** (0/24): those operands carry the `W_bg` profile *value* (`w1_profile`),
  not its in-plane *jet* (`w1_profile_d`), so ablating a jet direction does not move them. So the earlier
  "every object bites" AND the interim narrative "only geometric objects bite" are both wrong.
- **CONTROL_INDEPENDENCE (corruption) bites all 6 objects it covers** (CLOSURE, EVOLUTION_MASS_BALANCE,
  RELATIVE_FLUX, TRACTION, VIRTUAL_CONSTRAINT, VIRTUAL_WORK).
- **HOLD controls unchanged: 0 genuine nonzero** (the HOLD result above stands).

FACE_SHIFT is therefore not exercised by the §5 form/independence battery; it is validated instead by the
cross-engine residual (joins 16/16 at 0 post-`cccb4f9e`), the density-grounding build-leg form ablation
(`_measurements/s11c_a_wl_density_grounding_fix_review_grok/ablate_density_bg_form.*`), and the
material-transport adjudication. A fully re-run CONTROL_FORM coverage table for the abstract-symbol objects
remains an owed item.
