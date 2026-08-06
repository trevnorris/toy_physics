# Two fixes from the second build review — G1 and G3 only

Two legs reviewed the delta. Both confirm the 26 action agreements are real — one rebuilt the comparator
independently, with a stricter key, and reproduced all 26 residuals as identically zero. Four defects were found;
this directive covers two. G2 (NAMING_MISMATCH unfalsifiable on a bare-symbol payload) and G4 (the
bijection verification is untested) are DEFERRED by scope decision and stay recorded in
`S10_harness_fix_round2.md`'s successor — do not fix them here.

Fix these; do not refactor anything else, and do not commit.

## G1 — the canonical derivative atom drops the differentiated function's argument list

`engine_output_checks.py:931-950` (`_canonicalize_derivatives.canonical`), `CanonicalDerivative` at
`:174-186`. The atom keys on function **name** plus the multiset of `(variable, order)` and never reads
what the function depends on.

Measured end-to-end against a damaged copy of committed output:

```
PY_S10_MAIN_D2_Q1_LAGRANGIAN_EXPANDED:  u1(t,x1,x2) -> u1(t,x1,x2,x3)
  CROSS_ENGINE: agree=509 disagree=32 naming_mismatch=17      (identical to baseline)
  ACTION_VERDICTS total=26 AGREE=26 DISAGREE=0
  normalized[py] and normalized[wl] print as byte-identical strings

symbolic_equal(Derivative(u1(t,x1,x2),x2), Derivative(u1(t,x1,x2,x3),x2))  = False
same pair after _canonicalize_derivatives                                   = True
```

**This is a drafting error in the previous fix list, not a misreading of it.** That list said the argument
*order* of the differentiated function is an emission convention carrying no content. That is true of the
order and false of the **membership**: which coordinates a field depends on is content.

**Required:** the canonical key includes the function's argument list as a **sorted set**, so ordering is
normalised and membership is preserved. Two derivatives of fields with different dependence sets must not
compare equal.

Also `harness_ablation.py:423,436-437` — ACCEPTANCE 16 prints `derivative_arity=observable`, but its
oracle calls `symbolic_equal` on the **raw parse**, a path the comparator no longer takes. Measured:
raw=False, canonicalised=True. The one control that names this property now asserts the opposite of the
truth. Point the oracle at the path actually used, or retire the line.

Nothing else covers this: no cross-engine row declares field content, and the parametrised tests at
`test_engine_output_checks.py:261-269` cover order, variable and function only. Add arity to them.

## G3 — a test now forbids a finding

`test_engine_output_checks.py:655-661` asserts all 26 action rows equal `AGREE`. A future genuine physics
disagreement would fail the suite rather than be reported.

The previous test was criticised for asserting only `status in VERDICTS`; the replacement over-corrected
by pinning the verdict. **Required:** assert the property that makes the verdict meaningful — that every
declared action row reaches a comparison verdict, that the canonicalisation is live, and that the
single-difference controls (order, variable, function, argument list, coefficient, form) each still reach
`DISAGREE`. A physics disagreement must remain reportable without breaking the suite.

## Constraints

- Do not change `symbolic_equal`'s equality semantics; `DEFECT_REGISTER.md#f7` stays deferred.
- A physics disagreement must not exit non-zero and must not fail the test suite.
- Report what the instrument reads. Do not iterate toward any verdict — if `G1` changes a row's verdict,
  that is the measurement.
- No script over 600 seconds. Do not launch Mathematica; `/tmp/s9_wl.txt` exists and the suite needs it.
- Do not edit `steps/`, `paper/`, `REBUILD_HANDOFF.md`, either engine, or any committed engine output.

## Out of scope — record, do not fix

A leg measured that `PY_S10_*_D3_Q7_STIFFNESS` is byte-identical across all five packages while the WL
side tracks each package's stiffness form, so `*_q7_difference` reads `py=0` against a nonzero WL
residual. **That is an engine defect, not a harness defect**, and it belongs to the separate Q7 repair.
Do not touch the engines. If it is not already in `DEFECT_REGISTER.md`, add an entry recording it with
the measurement.

## Report back — under 20 lines

1. One line per `G1` and `G3`: what changed, with file and line.
2. Literal stdout of both configs' cross-engine and action summary lines, before and after.
3. For `G1`: the arity attack re-run, showing the verdict it now produces.
5. Anything you measured that contradicts this list.

Stop there. Do not commit.
