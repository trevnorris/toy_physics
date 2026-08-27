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
