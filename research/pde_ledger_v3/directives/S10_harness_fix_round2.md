# Eight fixes from two review legs

Both legs reviewed the build and returned NO. Every item below was measured by a leg and, where noted,
re-measured by the orchestrator. Fix these; do not refactor anything else, and do not commit.

**Do not iterate toward any verdict.** Several of these change what rows report. Report what the rebuilt
instrument reads. If a fix makes rows agree that previously disagreed, that is the measurement — but the
non-vacuity controls below are what entitle you to it, so build them first and keep them passing.

---

## F1 — No expression containing a derivative can compare equal across the engines

`engine_output_checks.py:158-172`, `:426-448`, and the `symbolic_equal` scalar branch `:670-677`.

The WL parser keeps every derivative as an opaque `sp.Symbol` named
`Derivative[0,1,0][u1][x1,x2,t]`. The SymPy side holds a real `sp.Derivative(u1(t,x1,x2), x2)`. These can
never be equal, so every action row is a guaranteed `DISAGREE` and the verdict column carries no
information.

Measured: 24 of the 26 declared action rows are `DISAGREE` with operands that are algebraically
**identical**. Independently confirmed on three rows by the orchestrator — `MAIN D3` and `MAIN D2`
Lagrangian and `XFORM_FULLGRAD D3` — residual identically zero after reconciliation.

**The fix, stated as the object rather than the mechanism:** a derivative's identity is its **function**
together with the **multiset of (variable, order) pairs**. The argument order of the differentiated
function is an emission convention and carries no content. Canonicalise both engines' derivatives to one
representation keyed on that identity.

Constraints, because this is the fix most able to manufacture agreement:

- Do not hardcode a coordinate list, coordinate names, or an argument ordering. The transform must be
  blind to which coordinate is time.
- Two derivatives that differ in function, in any variable, or in any order must remain unequal.
- **Non-vacuity, required and must be shown:** with the canonicalisation live, a coefficient control
  (`mu_R -> 2*mu_R`) and a form control (`MAIN` against `XFORM_DIVONLY`) must both still reach
  `DISAGREE`. A canonicalisation that also collapses those is erasing content, not reconciling notation.

## F2 — The action family label is unbound from the payload

`engine_output_checks.py:1145-1147`, `:1259-1276`.

`family` is accepted from the config as a free string and `action_family_reports` counts by that label
alone. Measured: a row named `main_d2_q1_lagrangian_action`, declared `family: lagrangian`, pointed at
`WL/PY_S10_MAIN_D2_Q3_DETERMINANT`, reports **AGREE** and `LAGRANGIAN_ACTION_COVERAGE: numerator=1
denominator=1 gaps=()`.

So the harness can report that the action was compared when it compared a determinant. The family label is
typed prose that no computation produced.

Bind the family to the declared tag identity: a row claiming membership in an action family whose tags are
not that family's declared tags is a declaration error, not a verdict.

## F3 — Naming mismatches are reported as physics disagreements

`engine_output_checks.py:1198-1207`. The undeclared-spelling detector fires only when
`_snake_to_lower_camel(a) == b`. Every other spelling divergence is absorbed into `DISAGREE` and printed
under `PHYSICS_DISAGREEMENTS`.

Measured on the committed S10 run: **69 of the 129** `DISAGREE` rows are pure name mismatches, and the
`NAMING_MISMATCH` count is **0** with `undeclared_spellings` raised on **0** rows. Examples: `M_B` vs
`quadraticFormRoute` (9 rows), `s` vs `coefficientScale`, `D` vs `braneDimension`, `g12` vs `g1x2`.

Two separate requirements:

1. **A name mismatch is not a physics disagreement.** Where the operands become equal under a bijection on
   leftover free-symbol names, report that as its own status, distinct from `DISAGREE`, and keep it out of
   `PHYSICS_DISAGREEMENTS`. Name the symbols involved.
2. **Declare only what is genuinely the same object, each with its reason.** `D` / `braneDimension` and
   `s` / `coefficientScale` look like the same quantity under two spellings; verify before declaring.
   ⛔ Do **not** bulk-declare every pair the bijection finds — that manufactures agreement. Anything you
   cannot justify stays reported and undeclared, which is the honest outcome.

⚠ `M_B` / `quadraticFormRoute` is a **route name**: it records which derivation route the engine took. If
the two engines genuinely took different routes that is a finding, not a spelling problem. Report it; do
not declare it away.

## F4 — The dimension layer gets quieter when it stops measuring

`engine_output_checks.py:2204-2210`. `operational_failures` inspects only `dimension.tables` and
`dimension.statuses`; it never reads `compared`, `package_disagreements`, or `proposition_sites`.

Measured on S9: deleting `dimension_sources` takes `DIMENSIONS[wl]` from `compared=312 homogeneous=312` to
`compared=0 homogeneous=0` while the operational-failure count stays at **2, identical to baseline**.
Emptying `cells` yields zero dimension findings — fewer than the healthy run. Mis-spelling the key as
`dimension_source` reproduces it with no error at all.

This is the governing failure mode still live in one layer. Give the dimension layer declaration-derived
coverage like the others: zero comparisons against a non-empty declaration is an operational failure.

## F5 — `NAMING_EFFECT` does not measure the naming rule

`engine_output_checks.py:1216-1220`, `:2374-2378`. `legacy` is computed from `selected`, which has already
had naming and identities applied, so the line reports the rule changed nothing.

Measured on S9: the rule's true effect is `AGREE 8 -> 12`, with 4 rows becoming `NAMING_MISMATCH` when it
is removed. The report prints `legacy_before_agree=12 declared_after_agree=12 changed_rows=0`.

Compute the "before" from the raw operands, before any naming or identity is applied. This line exists to
hold the naming rule accountable and currently exempts it.

## F6 — Two equation-branch mutations that no test catches

`_normalize_action` `:1035-1053`. Mutating `residuals.append(component.lhs - component.rhs)` to drop the
right-hand side, or to halve the scale, leaves the suite at `54 passed`.

The dropped-RHS variant is a semantic no-op on committed data — all 45 WL Euler-Lagrange right-hand sides
are zero — so it needs a synthetic fixture with a nonzero RHS. The scale variant is a real change and is
invisible only because those rows are `DISAGREE` either way; after `F1` it should become visible, but add
the test regardless so it does not depend on `F1` holding.

## F7 — Two dishonest acceptance items in `harness_ablation.py`

- `:101-120` ACCEPTANCE 3 claims positional residuals preserve scale and order, but its fixture is bare
  symbol lists that never reach the equation branch. With the scale mutation live it still prints `PASS`.
  Route the fixture through the real path.
- `:436-464` ACCEPTANCE 18 reads the form-versus-coefficient distinction off the config string
  `stiffness_form`, and applies the same `status changed OR operand text changed` oracle to both classes,
  which both satisfy trivially. It measures that two packages emit different text. Give it an oracle that
  can separate a form change from a coefficient change, or state plainly in its output that it cannot.

## F8 — Tests and config keys that cannot fail

- `test_engine_output_checks.py:501-507` asserts only that action rows reach `status in VERDICTS`, so it
  passes identically for a harness that returns `DISAGREE` unconditionally. Make it check something that
  distinguishes those cases.
- Mistyped config keys are silently ignored: `symbol_identities` -> `symbol_identity` and
  `dimension_sources` -> `dimension_source` both run without error, each silently disabling a layer.
  An unrecognised key at a declaration site is a declaration error.

---

## Constraints

- Do not change `symbolic_equal`'s equality semantics; `DEFECT_REGISTER.md#f7` stays deferred.
- Do not restore any hardcoded alias table in the comparator.
- A physics disagreement must not exit non-zero.
- No script over 600 seconds. Do not launch Mathematica; `/tmp/s9_wl.txt` exists and the suite needs it.
- Do not edit `steps/`, `paper/`, `REBUILD_HANDOFF.md`, either engine, or any committed engine output.

## Report back — under 30 lines

1. One line per `F1`-`F8`: what changed, with file and line.
2. Literal stdout of both configs' cross-engine and action summary lines, before and after.
3. The `F1` non-vacuity result: coefficient control and form control verdicts with the canonicalisation
   live.
4. For `F3`, the list of symbol pairs you declared and the list you left reported and undeclared, each
   with your reason.
5. Anything you measured that contradicts this list.

Stop there. Do not commit.
