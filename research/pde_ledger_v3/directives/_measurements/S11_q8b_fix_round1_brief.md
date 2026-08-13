# Measurements behind `S11_q8b_fix_round1_brief.md`

Commands and literal output. CLAUDE.md rule 2. Regenerate; do not transcribe.
⛔ NO S11 RESULT VALUE APPEARS HERE — rule 5. Taken at `94e14e42`, 2026-08-13.

## 1 · The mechanism (orchestrator code read)

```
$ sed -n '1751,1758p' research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
    status, certificate, change_residuals = analyzer.analyze(
        change_residuals, generic_value, recompute,
    )
    payload_value = generic_value if status in ("CONSTANT", "VARIES") else NOT_DEFINED_ON_COMPONENT
    emit_quantity(
        package, n, quantity, component_count_record(status, payload_value),
        root=root, stratum=stratum,
    )
```
The caller passes counts computed on the restricted component matrix (`:1949-1958`:
`sp.Integer(rank)`, `sp.Integer(stacked_rank)`, … from `m_r`/`stacked` of the component route) — the
same objects whose point-invariance both legs verified byte-exact, so `generic_value` is a
component-wide obtained count, not a point evaluation.

## 2 · The erasure, measured by both legs independently (leg-produced, orchestrator-read)

- Fresh agent: `/tmp/s11q8bs_leg_jK64/p06_undecided_value.py` → `.stdout` — `UNDECIDED_RECORD ...
  GENERIC_VALUE_THE_ENGINE_OBTAINED=<int>` beside emitted `VALUE: NOT_DEFINED_ON_COMPONENT`, on
  multiple strata/quantities through the real promotion path.
- Grok: `/tmp/s11q8bs_leg_R7jt/p12_undecided_value_probe.py` → `p12_wrap.stdout` —
  `CLASSIFICATION: OBTAINED_THEN_ERASED` on the same quantity families.
Both wrapped `emit_component_count` and called the original unchanged.

## 3 · The spec obligation

```
$ sed -n '737,745p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
STATUS_TOKEN: <the same token emitted under _STATUS>
VALUE:        <the count on the component where it is defined there, whatever the status;
               the single token NOT_DEFINED_ON_COMPONENT where the engine did not obtain one>
```

⚠ ⭐ Under `VARIES` the answer is the **pair** — the `VALUE` the engine obtained together with
`_CHANGE_LOCUS_*` — ⛔ neither field alone. ⚠ Where the engine obtains no single count off the change
sub-locus, `VALUE` is `NOT_DEFINED_ON_COMPONENT` and the change locus carries the object; ⛔ do not
manufacture one.

`:754-756` (the UNDECIDED+NOT_DEFINED case) begins *"If the engine cannot build a component-wide
Q3/Q4 object"* — its antecedent is cannot-build, not constancy-open. Steelman recorded by the agent
leg; adjudicated by the spec text above (`"whatever the status"`) plus `VARIES` keeping `VALUE`.

## 4 · Fold-time: both brief legs confirmed the item and broke the acceptance (verified, folded)

- Grok (`~/.s11_build/q8b_fix1_brief_leg_grok.log`, probes `/tmp/s11q8bf1_leg_H2M4/`): an undriven
  `/tmp` copy emits zero count records ⇒ item-2 fail condition never fires; VARIES keeps VALUE 6/6;
  partial fix leaves 4 erased.
- Codex (`~/.s11_build/q8b_fix1_brief_leg_codex.log`, probes `/tmp/s11q8bf1_leg_iRyz/`): a REAL driven
  promotion can be all-CONSTANT (11/11) making HEAD and the correct fix output-identical; and a wrong
  non-sentinel VALUE (0 in all 12 records) passed every written failure clause —
  `WRITTEN_EXPLICIT_ITEM2_FAILURES=0` beside `FULL_PROPERTY_FAILURES=12`.
- Orchestrator verification: both holes are visible by direct reading of the acceptance text against
  the measured facts (zero live promotions; the one-directional fail clause). Fold: item-2 now requires
  an UNDECIDED-with-obtained-count witness at baseline, the same driven component both runs, and
  bidirectional equality between `VALUE` and what the computation bound.
- Both legs independently re-adjudicated the spec (`:737-739` / `:742-745` / `:754-756` one-way
  antecedent): the erasure is not licensed; item 1's property unchanged by review.

## 5 · Everything else at `94e14e42` is leg-verified clean (not re-litigated here)

Both legs' full reports: `~/.s11_build/q8b_script_leg_agent.log`, `~/.s11_build/q8b_script_leg_grok.log`
(probes `/tmp/s11q8bs_leg_jK64/`, `/tmp/s11q8bs_leg_R7jt/`). Orchestrator spot-checks: the caller read
above; live MAIN D2 via own interceptor (5.5 s, dispositions unchanged vs pre-build).
