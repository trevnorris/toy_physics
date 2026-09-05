# S11c-c1 WL repair DIRECTIVE — decision-leg review (run RETROACTIVELY)

The repair directive `directives/S11c_c1_wl_repair_directive.md` was built against **without** its rule-7 decision
legs. This record holds those legs, run **retroactively** (the gap that triggered the remediation;
[[feedback_directive_design_review]]). Prompt (committed): `directives/_legs/S11c_c1_wl_repair_directive_review.md`.
Orchestrator-written directive → **Codex + Grok**. A third pass (Codex **arc-verify**) checked the audit trail +
repo state.

## Commands (literal)

```bash
# Codex decision leg (retroactive)
codex exec -c model_reasoning_effort=xhigh "$(<…/_legs/S11c_c1_wl_repair_directive_review.md)" \
  > …/scratchpad/codex_wl_repairdir.log 2>&1        # exit 0; 119,266 tokens
# Grok decision leg (retroactive)
grok --prompt-file …/_legs/S11c_c1_wl_repair_directive_review.md --cwd /var/projects/toy_physics \
  --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain \
  > …/scratchpad/grok_wl_repairdir.log 2>&1         # exit 0
# Codex arc-verify (records + repo state)
codex exec … > …/scratchpad/codex_c1_arc_verify.log 2>&1
```

## Verdicts — the directive is NOT sound; the two legs caught DIFFERENT defects (one alone would have shipped it)

- **Codex** (`codex_wl_repairdir.log`): "The directive is **not fully sound**. I found two MUST-FIX design
  defects." (1) **R2 "independent data" leaves the bulk `φ` unbound from the face drive** — the routes must be
  "independent constructions from the same supplied boundary drive," not "independent data"; the branch-flip
  invariant alone would pass an arbitrarily-normalized `φ`. (2) **R4b does not bind the parity labels to the
  face→(δW,ζ_c) map** — "distinct blocks can be swapped and still satisfy every explicit R4b comparison"; needed a
  computed forward/inverse basis-binding residual. Codex: "I found no separate rule-5 leak."
- **Grok** (`grok_wl_repairdir.log`): "not sound." Its MUST-FIX = **R1's invariant leaks residual=0 and re-enters
  at a result** ("residual that is zero for the correct labelling" — rule 5). NITs: **`NONINVERTIBILITY_CONDITION`
  over-protected → the single-leg freeze persists in `fredholmFunctionSpaceOperator`** ("the one inherited weakness
  a directive-faithful re-review was structurally unable to see"); **R4b labelling** ("built code inherited the
  swap"; caught by repair-1 re-review leg-1 as a residual NIT); R2 "independent data" (NIT). Grok did NOT flag the
  R1 leak's counterpart Codex found, and Codex did NOT flag the R1 leak or the Fredholm freeze — **different legs,
  different defects**.
- **Codex arc-verify** (`codex_c1_arc_verify.log`): the rule-2 audit trail is incomplete (records point to
  ephemeral `/tmp` logs, not literal output; the "+443/−103 one file" is engine-scoped, not the 4-file commit) and
  **neither c1 engine's canonical `.out` is committed** (WL nor SymPy) — so "per-engine DONE" is real but "c1 DONE"
  overstates.

## Which reached CODE (verified against `13f0bd2c`, rule 13)

- **`NONINVERTIBILITY_CONDITION`/`fredholmFunctionSpaceOperator` single-leg freeze — PROPAGATED** (both `nZero`/
  `gZero` on `momentumOutput`, `.wl:580-597`). Emit-only; the re-review was told to treat it as byte-identical
  success, so it was structurally blind. → repair-2 R1.
- **R1 leaked-zero, R2 "independent data": directive-only** — the built code genuinely unfroze the composition
  (`operatorCompositionRecordFromDerivation`) and bound `φ` to the face velocity (`φ=(V/(iq))e^{i(k′·x+qw)}`), so
  neither reached code. Documentary erratum on the repair directive; not inherited.
- **R4b "parity swap" — the finding itself was WRONG (see correction below).**

## ⭐ CORRECTION by the repair-2 decision gate (2026-09-04) — the R4b "swap" was a naive misreading

The repair-1 re-review leg-1, both retroactive legs, and the plan's confirmed-defect #3 all reported a "parity
swap" in `PERMEABLE_PORT_HERMITIAN` (that `DELTA_W` should be the odd `(ζ₊−ζ₋)` combination). The **repair-2
decision legs** (Codex + Grok; `_measurements/S11c_c1_wl_repair2_directive_review.md`) + my own verification
(rule 13) **refute** it: `PERMEABLE_PORT_HERMITIAN` is **correct** — its diagonal blocks are the congruence
`Aᵀ·diag(P₊,P₋)·A` (`A=faceToParityMatrix=[[1/2,1],[1/2,-1]]`), both even combinations with the odd combination as
coupling, correct given the outward orientation `V_s=s∂_tζ_s` (`S11c_a:58`). The naive reading is a **category
error** (a coordinate's parity ≠ the parity of its diagonal block in a quadratic form) that would vanish the
thickness port at equal faces. **The genuine propagated parity defect is the dead parity axis in
`PERMEABLE_DISSIPATION_VS_OMEGA_TAU`** (`.wl:1713-1736`), fixed by repair-2 R2. This is the gate working twice: the
retroactive legs found a real gap (the freeze, the unbound routes), and the repair-2 gate corrected the one finding
that was itself a misreading — **before** any build acted on it.

## Disposition
⇒ the FULL REMEDIATION (`_measurements/S11c_c1_wl_remediation_plan.md`): repair-2 (R1 Fredholm two-leg; R2
DISSIPATION parity combination) with its **own** decision legs run BEFORE the build (done — this time), records +
both `.out`, DONE-overstatement corrected.
