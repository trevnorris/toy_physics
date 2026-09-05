# Decision-leg review — S11c-c1 comparator REPAIR directive (before the repair build)

You are an independent, adversarial reviewer of an **orchestrator-written repair directive**, BEFORE the repair
is built (rule-7 decision gate). The comparator it repairs was already re-reviewed SOUND (two legs CLEAN). This
repair applies exactly two surgical fixes. Your job: confirm each fix is correct and safe, that the repair
changes NO soundness behaviour, and that nothing is under- or over-scoped — so it is fixed in one fold before a
builder runs. Ground every judgement in a file/line or a literal command output you actually ran.

Working dir: `/var/projects/toy_physics` (branch `ledger-v3-rebuild`).

## Artifact under review
`research/pde_ledger_v3/directives/S11c_c1_comparator_repair_directive.md` (baseline comparator commit
`7141e6ad`: `research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py` + its test file).

## Context to read
- The frozen T7 contract `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:580-587` (measurement-only,
  three-valued, no pre-registered fold).
- The transliteration `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py:147` (`mechanical_lower_camel`).
- The committed comparator (its spelling map, its `Inactive[Greater]→HeldInactiveGreater` held-parse pattern,
  its 5 seals, its `raw_control_case` whitelist).
- The two `.out` streams (`datalad get` first; sample with grep — do NOT full-parse):
  `scripts/out/S11c_c1_bulk_closure_sympy_audit.out`, `mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out`.

## Check each as CLEAN / GAP / INACCURATE
1. **R1 correctness.** Verify yourself: `python3 -c "import S11b_cross_engine_comparator as m;
   print(m.mechanical_lower_camel('c_s0'))"`. Is it `cS0`? Is `Symbol('c_s0')` a real PY reserved symbol in the
   PY stream, and is `cS0` in the WL stream? Is `cS0←c_s0` therefore the same KIND of bare-symbol mechanical
   fold as the existing 11 (`rhoM←rho_m`…)? Confirm the directive's activation-gating + injectivity + no-arg-strip
   requirements are the right guards (it must fire only when PY `c_s0` is present, keep collisions=0, and not
   open any applied-head strip). Flag any way adding `cS0` could manufacture a FALSE agreement (e.g. if `cS0`
   ever names a DIFFERENT physical quantity than the sound speed `c_s0` — check the payloads).
2. **R2 correctness/safety.** The WL energy leaves `Inactive[Integrate][(integrand)]` (range-less) and
   `Inactive[Limit][Inactive[Integrate][…]…]` currently fail to parse (surfaced as `<PARSE_FAILED>` in
   UNPAIRED leaves). Is held-parsing them (a `HeldInactiveIntegrate`/`HeldInactiveLimit` AppliedUndef carrying
   raw args, like the existing `HeldInactiveGreater`) the correct, faithful fix? Confirm the directive forbids
   EVALUATING them and forbids reconciling them against the PY operand. Is "display/completeness only, no
   residual change" accurate given these leaves are join=0? Any risk the held-parse changes a residual elsewhere?
3. **Soundness preserved.** Does the "what must NOT change" list correctly protect every soundness element (the
   5 seals stay unreconciled + load-bearing, the `raw_control_case` whitelist, three-valued residual objects,
   no verdict token, per-family accounting, exit-0)? Is there any element it forgets to protect?
4. **Rule 5.** The R1 test asserts a synthetic `cS0`/`c_s0` pair residuals to `Integer(0)`. Is that a value-free
   synthetic-fixture acceptance (testing the fold MECHANISM, like a repoint ablation — legitimate per
   `S11b_comparator_build_directive.md:68-99`), or does any instruction leak an expected residual on a MEASURED
   stream? Flag any measured-stream target.
5. **Scope.** Is anything under-scoped (a needed change omitted) or over-scoped (a change that should not be in a
   surgical repair)? Would following this directive predictably require another round?

## Output
Numbered CLEAN/GAP/INACCURATE per item with the command + literal output you ran, then a one-line overall: is
this repair directive sound to hand a builder as-is, or what must be folded first. ⛔ Do NOT propose an expected
measured residual to bake in.
