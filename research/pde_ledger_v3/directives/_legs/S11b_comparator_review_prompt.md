# Independent review — S11b cross-engine comparator (a Codex-written instrument)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`

A Python tool that joins two CAS transcripts (`PY_S11B_<QUANTITY>` and `WL_S11B_<QUANTITY>` tag streams) and
classifies each shared object AGREE / DISAGREE / UNDECIDED / UNCOMPARED, plus unpaired-tag coverage. It is
the instrument that certifies (or refutes) cross-engine agreement, so its ONE cardinal sin is to report
AGREE (or bury in a non-DISAGREE bucket) where the two payloads genuinely differ.

## What to check
Does the comparator faithfully implement its frozen contract, and — the adversarial core — **can any code
path make a genuine content/sign/dimension/structure disagreement come out as anything other than DISAGREE?**

## What you are handed
- The artifact above, and its mechanical precedent `research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py`
  (+ `S10_cross_engine_comparator_repoint_ablation.py`).
- The frozen contract it must implement: `research/pde_ledger_v3/directives/S11b_comparator_build_directive.md`
  (the T7/T5 rules and the 7-item "S11b delta"). Read this as the source of truth for what the comparator
  must do; then read the code.
- The two real engine transcripts, for running the comparator end-to-end:
  `scripts/out/S11b_interface_coupling_law_sympy_audit.out`,
  `mathematica/out/S11b_interface_coupling_law_mathematica_audit.out`. ⚠ Use these to check the tool does not
  crash and does not HIDE disagreements — ⛔ but your VERDICT on correctness must come from ablation +
  synthetic probes you construct, not from whether the real output "looks right."

## Required method — this is a SCRIPT; ablation is mandatory
Derive/construct independently. Write your own small PY/WL transcript fixtures (placeholder symbols/values)
and run the comparator on them; save every fixture and its literal stdout to named absolute paths and report
them. A prose claim about behaviour is discarded.

⛔⛔ **FORM-ABLATE each load-bearing rule and report the literal diff.** Change the STRUCTURE of a rule and
confirm a synthetic fixture catches it. At minimum:
- **Boolean rejection (delta #2):** does it reject the boolean at the LEAF only, so a differing ALGEBRAIC
  sibling in the same Association still DISAGREEs? Construct a record with `TEST_OBJECT -> True` beside a
  sibling `X -> a` (PY) vs `X -> b` (WL). If the record is reported anything but DISAGREE, that is the defect.
  Also test `False` vs `0`, `True` vs `1`, `S.true` vs `True` — none may be AGREE; and `0` vs `0` must stay
  residualable.
- **UNDECIDED vs UNCOMPARED (delta #1):** can a parse failure or an unreduced residual be buried as UNDECIDED
  (or AGREE) instead of the visible UNCOMPARED? Feed a payload the parser chokes on AND a genuinely
  different one — the difference must not vanish.
- **Structure/key (delta #3,#4):** tuple↔Association and unequal Association key-sets must be STRUCTURE/KEY
  DISAGREE. Can an intersection-only key walk drop an extra key that carries a real difference and still
  report AGREE on the rest?
- **Dimension vectors (delta #5):** are the L,T,M vectors extracted from BOTH shells and residualled? Feed
  vectors differing in one component — must be DIM DISAGREE, ⛔ never normalized/absorbed to AGREE.
- **Injective transliteration (delta #6):** feed two source symbols that collide under the map with genuinely
  different payloads — the collision must not let a real difference residual to zero.
- **`_LOCAL_` (delta #7):** a `LOCAL_` tag must never be residual-scored; a non-local one-engine tag must
  appear as unpaired coverage.
- **Repoint ablation:** substitute a different synthetic object's payload under a previously-agreeing name;
  the row must flip to DISAGREE.

Also: is the build's OWN synthetic acceptance decisive, or would a defective comparator also pass it? Name
any acceptance test a bad implementation survives.

## Physics filter
Report a finding only if it catches a way the COMPARATOR could be WRONG — most importantly, a way it could
report agreement (or a non-DISAGREE bucket) on a genuine disagreement, or crash/hide on real input. Do not
report style. A leg that finds nothing is weak evidence: show the ablations you ran and their literal diffs.

## Ablation sandbox
Copy the comparator to /tmp and ablate the COPY. ⛔ Never modify the working tree.

## Report — under ~30 lines
Numbered findings; for each: the code line, the ablation/fixture + its literal stdout path, the wrong
classification it produces, and the fix. End with a one-line verdict: is the comparator sound to FREEZE and
run, or must it be repaired first?
