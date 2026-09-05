# Measurements — S11c-c1 comparator repair directive

Rule 2 (orchestrator half): the claims the repair directive makes, each with the command behind it. The repair
changes NO soundness behaviour; it applies two surgical fixes to the re-reviewed-SOUND comparator (baseline
`7141e6ad`).

## R1 — the `cS0 ← c_s0` fold is a genuine mechanical inverse (the build's directive was wrong)

```
$ python3 -c "import S11b_cross_engine_comparator as m; print(m.mechanical_lower_camel('c_s0'))"
cS0
$ grep -oE "Symbol\('c_s0'" scripts/out/S11c_c1_bulk_closure_sympy_audit.out | head -1   # real PY symbol
Symbol('c_s0'
$ grep -oc 'cS0' mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out   # WL uses cS0 (dense)
221355
```

So `cS0` IS `mechanical_lower_camel("c_s0")` for the real PY `Symbol('c_s0')` — the same kind as the 11 folds
already present. The build's directive claimed `=="cs0"` (an arithmetic error: the first piece `c` keeps case,
the second `s0`→`S0`). Corrected in the build directive's fold-1 paragraph + folded by this repair. The re-review
Claude Agent's controlled synthetic proved folding it collapses the `cS0`/`c_s0` residual to `Integer(0)` (pure
spelling artifact); it fires in zero currently-joined residuals but WL `cS0` is dense in the deferred PERMEABLE_*
families, so the fold must be in place before that heavy run.

## R2 — the WL held-integral energy leaves currently fail to parse (surfaced honestly)

Both re-review legs measured: WL `Inactive[Integrate][(integrand)]` (range-less) and
`Inactive[Limit][Inactive[Integrate][…]…]` raise `InputError: Inactive[Integrate] does not have integrand and
range` → 56 `ENERGY_RESIDUAL` + 4 `ENERGY_FACE_TRACTION` leaves shown as `<PARSE_FAILED>` in UNPAIRED (join=0)
leaves (no residual corruption). R2 held-parses them (a `HeldInactiveIntegrate`/`HeldInactiveLimit` AppliedUndef,
the existing `HeldInactiveGreater` pattern) so the operand is displayed — display/completeness only.

## Decision-leg gate (rule 7 — BEFORE the repair build)

Orchestrator-written repair directive → legs = Codex(`gpt-5.6-sol` xhigh, doc-review policy) + Grok (identical
prompt `_legs/S11c_c1_comparator_repair_directive_review.md`). Launched:

```
$ codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh "$(cat _legs/S11c_c1_comparator_repair_directive_review.md)" > <scratch>/codex_repair_dir_review.txt 2>&1   # detached; -m flag accepted (v0.153.4)
$ grok --prompt-file .../S11c_c1_comparator_repair_directive_review.md --cwd /var/projects/toy_physics --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain > <scratch>/grok_repair_dir_review.txt 2>&1   # detached
```

RESULTS (both EXIT=0; converged; `-m gpt-5.6-sol` accepted on Codex v0.153.4). Both legs: **R1 CLEAN**
(cS0←c_s0 is a genuine bare-symbol mechanical fold, injective 12/12/0-collisions, activation-gated, no
arg-strip, preserves `positive=True`; WL engine binds `"c_s0"→cS0` at `.wl:125`); **R2 INACCURATE**; soundness
envelope CLEAN; **rule 5 CLEAN** (synthetic cS0→0 is value-free); **scope GAP**.

FOLDS APPLIED (one pass): (F1) R2 form corrected — the WL integral is **4-arg TRIPLE-range**
`Inactive[Integrate][integrand,{xOne,-∞,∞},{xTwo,-∞,∞},{xThree,-∞,∞}]` + nested
`Inactive[Limit][…,outwardDistance->Infinity]`, 100 occurrences across **THREE** ENERGY families (RESIDUAL 80,
BULK_FARFIELD_FLUX 16 [was omitted], FACE_TRACTION 4); (F2) R2 is **c1-preprocessor-only** — ⛔ do not widen
S11c-a's 1-range handler, keep the inherited `Inactive[Integrate][f,{x,a,b}]→BOUND_INTEGRAL` path
(`test:223-234`) and `Inactive[Equal]→HeldEqual` intact, never a blanket hold-every-Inactive; (F3) R2 test uses
the real 4-arg form + nested limit + a BoundIntegral-preservation test; (F4) R1 test — the existing
`test_noninherited_sound_speed_spelling_is_exposed` (`test:214-219`, asserts `cS0` stays bare) is REPLACED (R1
folds it), and the positive activation fixture must carry the FULL 12-symbol vocabulary (the gate is all-or-none,
`REQUIRED_C1_PY_SYMBOLS<=names` at `:2127`) so the builder can't weaken the gate to pass a single-symbol fixture;
(F5) "keep 30/30 green" → keep all except the replaced test; the accounting SCHEMA is unchanged, only the
parse-failure COUNTS change (expected).

Leg logs (outside repo): `<scratch>/codex_repair_dir_review.txt` (sol, 177k tokens, verdict :4913-5052),
`<scratch>/grok_repair_dir_review.txt` (grok).

## Repair build + re-review (rule 9 — filled after the gate)

_pending (repair build = Codex `gpt-6-astra` high, code policy)._
