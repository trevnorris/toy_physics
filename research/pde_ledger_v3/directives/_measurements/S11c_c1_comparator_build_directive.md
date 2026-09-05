# Measurements — S11c-c1 comparator build directive

Rule 2 (orchestrator half): the factual claims the directive
`directives/S11c_c1_comparator_build_directive.md` makes about the two committed `.out` artifacts, each with the
command that produced it and its literal output. The directive states **no residual target** (rule 5); the only
claims about the artifacts are structural (tag counts, family census, payload container shapes).

## Inputs (git-annex pointers; content resolved before reading)

```
$ datalad get research/pde_ledger_v3/mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out \
              research/pde_ledger_v3/scripts/out/S11c_c1_bulk_closure_sympy_audit.out
action summary:
  get (notneeded: 2)
$ ls -lL <both>
-r--r--r-- ... 81881078 ... S11c_c1_bulk_closure_mathematica_audit.out
-r--r--r-- ... 90722857 ... S11c_c1_bulk_closure_sympy_audit.out
```

## Tag counts (directive: "PY 63 tags", "WL 51 tags")

```
$ grep -cE '^PY_S11CC1_[A-Z0-9_]+:' scripts/out/S11c_c1_bulk_closure_sympy_audit.out
63
$ grep -cE '^WL_S11CC1_[A-Z0-9_]+:' mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out
51
```

## Family census (directive: "identical 50 non-LOCAL family set"; LOCAL_ differs by design)

```
$ diff <(grep -oE '^PY_S11CC1_[A-Z0-9_]+:' <PY> | sed 's/:$//;s/^PY_S11CC1_//' | grep -v '^LOCAL_' | sort) \
       <(grep -oE '^WL_S11CC1_[A-Z0-9_]+:' <WL> | sed 's/:$//;s/^WL_S11CC1_//' | grep -v '^LOCAL_' | sort)
(no output)  => IDENTICAL non-LOCAL family sets
PY non-LOCAL: 50   WL non-LOCAL: 50
```

PY `LOCAL_` tags (13, export/fold bookkeeping): `LOCAL_CONSUMER_CLOSURE LOCAL_DELTA_MINIMAL_ASSERTION
LOCAL_EXPORT_CANDIDATE_KEY_OPERANDS LOCAL_EXPORT_F9C_WRITES LOCAL_EXPORT_ROUNDTRIP LOCAL_FOLD_AUDIT
LOCAL_LOOKUP_MANIFEST_ASSERTION LOCAL_OPERATIONAL_EXCEPTIONS LOCAL_RUN_TASKS LOCAL_SECTION8_REPORT
LOCAL_SKIPPED_TASKS LOCAL_TAGS LOCAL_TASK_TIMING_SECONDS`. WL `LOCAL_` tags (1): `LOCAL_TAG_NAMES`.

## Payload container shapes (directive: the per-family VALUE structures) — literal heads

The directive's per-family shape descriptions were read directly off the payloads; representative heads (first
~220–280 chars), which also ground the SEALED representational differences (rule 5/6):

- `PY DTN_OPERATOR: Tuple(Tuple(Tuple(Str('LAB_HELD'), Integer(1)), Tuple(Tuple(Str('VALUE'), Tuple(Tuple(Symbol('s11cc1_flat_impedance_operator', commutative=False), ...`
- `WL DTN_OPERATOR: <|"LAB_HELD" -> <|"COMPOSITION" -> <|"EXPRESSION" -> OperatorSum[FourierMultiplier[(omega*rhoM)/qOut[omega, {kOne, kTwo, kThree}]], OperatorComposition[...`
  → DtN operator ALGEBRA differs (noncommutative Symbols vs OperatorSum/…); SEALED, raw_control_case.
- `PY DTN_KERNEL: ... Pow(Symbol('s11cc1_q_out_output'), Integer(-1)) ...` vs
  `WL DTN_KERNEL: <|"LAB_HELD|FACE_1" -> <|"OUTPUT_LEG" -> {kOne,kTwo,kThree}, "INPUT_LEG" -> {kPrimeOne,...}, "KERNEL" -> <|"EXPRESSION" -> ((I/2)*omega*rhoM*(kOne*kPrime...`
  → two-momentum legs: PY opaque `s11cc1_q_out_output/_input`, WL live `qOut[omega,{k…}]`/`qOut[omega,{kPrime…}]`; SEALED.
- `PY ZERO_JET_RESIDUAL: Tuple(Tuple(Str('VALUE'), Tuple(Add(Mul(Integer(-1), Symbol('omega'), Pow(Symbol('q_out'), Integer(-1)), Symbol('rho_m', positive=True)), Mul(Symbol('omega', real=True), Pow(Symbol('q_out'), ...`
  vs `WL ZERO_JET_RESIDUAL: <|"LAB_HELD" -> 0, "MATERIAL_ADVECTED" -> 0|>`
  → PY mixes plain `Symbol('omega')`/`Symbol('q_out')` with `Symbol('omega', real=True)`; ω real-assumption difference; SEALED.
- `PY CONTROL_FORM_RESIDUAL: Tuple(Tuple(Tuple(Str('DTN'), Str('LAB_HELD'), Str('RHO4_CONSTANT'), Integer(1)), ...` vs
  `WL CONTROL_FORM_RESIDUAL: <|"DTN|LAB_HELD|RHO4_CONSTANT|DIRECTION_1" -> ...` → PY positional `Integer(1)`, WL typed `DIRECTION_1`.
- `PY PERMEABLE_PORT_HERMITIAN: (LAB_HELD, THICKNESS, PROPAGATING, PROPAGATING) -> VALUE -> FULL_BARE_DTN...` vs
  `WL PERMEABLE_PORT_HERMITIAN: <|"LAB_HELD|RHO4_CONSTANT|OUTPUT_PROPAGATING|INPUT_PROPAGATING|PARITY_DELTA_W" -> <|"PORT_MATRIX" -> ...` → parity/regime key SHAPE differs; adjudicate post-run.

(The `grep`/`awk cut -c1-280` commands that produced these heads were run against the two `.out` above.)

## Decision-leg gate (rule 7 — BEFORE the build)

Orchestrator-written directive → legs = Codex + Grok (identical prompt
`directives/_legs/S11c_c1_comparator_directive_review_prompt.md`). Launched:

```
$ codex exec -c model_reasoning_effort=xhigh "$(<_legs/S11c_c1_comparator_directive_review_prompt.md)" > <scratch>/codex_c1_comparator_directive_review.txt 2>&1   # bkzpxonv8
$ grok --prompt-file .../S11c_c1_comparator_directive_review_prompt.md --cwd /var/projects/toy_physics --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain > <scratch>/grok_c1_comparator_directive_review.txt 2>&1   # bmz9qfcq9
```

RESULTS (both EXIT=0; the gate found the directive **not sound** — do not build as-is; fold once). The two
legs converged on the majors and each added distinct catches; every finding was grounded in engine-source lines
or literal payload bytes, and I re-verified the highest-impact ones myself (rule 13):

Verified against `scripts/S11c_c1_bulk_closure_sympy_audit.py`: `FACES=(1,-1)` (:120), `DIRECTIONS=(1,2,3)`
(:122), `PARITIES=("THICKNESS","CENTRE_SHIFT")` (:123), `REGIMES=("PROPAGATING","EVANESCENT","GRAZING")` (:124);
DtN key `=(anchoring,face)` (:1070,:1185), `=(anchoring,face,density)` (:1247); CONTROL_FORM key
`=(object,anchoring,density,direction)` (:1628). WL `DIMENSIONS` field is `OBJECT_DIMENSIONS`+`DERIVATION_EQUATIONS`,
NOT `NAMED_OBJECTS` (count 0). Verified in `S11b_cross_engine_comparator.py`: three-valued residual types
`Mismatch`(:100)/`BooleanNotResidualable`(:108)/`ResidualFailure`(:114)/`UndecidedResidual`(:121), `classify_residual`(:859),
and `residual()` DOES subtract scalar Exprs (`sp.factor(sp.cancel(sp.together(py-wl)))` ~:823). Repoint precedent
`S10_cross_engine_comparator_repoint_ablation.py:38-47` + `S11b_comparator_build_directive.md:91-93` ("a symbol
rename is NOT a repoint").

FOLDS APPLIED (one pass): (F1) un-ban the three-valued residual OBJECTS — forbid only a per-case VERDICT/status
token, require the inherited `BooleanNotResidualable`/`Undecided`/`ResidualFailure`/`Mismatch` to be PRINTED as
`A_minus_B`; (F2) regime `OUTPUT_PROPAGATING` is NOT a "spelling" of `PROPAGATING` — type the literal token,
seal the regime reconciliation; (F3) fold-1 spelling list restricted to genuine `mechanical_lower_camel`
inverses of REAL PY reserved names — removed `qOut`/`kOne`/`kPrimeOne`/`BoundaryMeasure`/`momentumSwapOutput`/
`faceVelocityInput`/`w1Jet*` (applied heads stay applied; the arg-strip is the sealed defect); (F4) μ_θ registry
`mu_theta_L/M` does NOT exist here — PY is `s11cc1_mu_theta_<anchor>_<±>`+`mu_theta_drive`, WL opaque `muTheta`;
kept SEALED, no registry surgery; (F5) pin Integer axis per family (FACE {1,-1} vs DIRECTION {1,2,3}); (F6)
replace wrong hand-specified VALUE shapes with a discover-per-family mandate + the copied-sibling traps named
(`OBJECT_DIMENSIONS`, `extract_term_origins`/`extract_control` schema mismatch); (F7) `raw_control_case` scope =
EXACTLY `DTN_OPERATOR`+`NONINVERTIBILITY` OPERATOR, else descend to every reachable scalar leaf; (F8) drop
`join>0` from DoD (manufacture-joins pressure + counts outer cases not leaves) → no-silent-0-extract + leaf
coverage; (F9) repoint ablation → object-binding (substitute a neighbour's payload under a paired name), not a
map mutation; (F10) test battery adds three-valued/undecided-vs-parse-failure + nested-boolean-beside-algebraic
+ FACE≢DIRECTION + Integer(-1)-rejected-as-DIRECTION + raw_control_case-refused-on-reachable-scalar; (F11) fix
the "scalar A−B undefined" wording (residual DOES subtract scalars; undefined only for booleans/structural
mismatch); (F12) pin a closed axis vocabulary for scenario/limit/port/parity tokens (no builder axis-vs-OBJECT
choice). Rule 5 CLEAN on real-payload residual targets (both legs); ω seal generalized off ZERO_JET.

Leg logs (outside repo): `<scratch>/codex_c1_comparator_directive_review.txt` (Codex, 8 items, verdict
:8144-8421), `<scratch>/grok_c1_comparator_directive_review.txt` (Grok, 8 items).

## Build + re-review (rule 9 — filled after the gate)

**Build (Codex, detached; default model — pre-policy):** `codex exec -c model_reasoning_effort=xhigh
--sandbox danger-full-access "$(cat <directive>)"`, EXIT=0, 493,658 tokens. Deliverables:
`scripts/S11c_c1_cross_engine_comparator.py` (88 KB, ~2175 lines) + `scripts/test_S11c_c1_cross_engine_comparator.py`
(22 KB). Builder self-report: 50 families, 1,080 joins across 27 families, no family extracts 0 leaves; 11
bare-symbol folds added (rhoM…tauV — the verified list, NO qOut/kOne); raw whitelist exactly
`{DTN_OPERATOR.WHOLE_OBJECT, NONINVERTIBILITY_CONDITION.OPERATOR}`; 30/30 new + 48/48 combined tests pass; peak
RSS 436 MiB; full all-family symbolic residual NOT run (heavy — deferred, as the directive allows). Builder
flagged my directive error: `mechanical_lower_camel("c_s0")` returns `cS0`, not `cs0`.

**Re-review legs (Codex-written script → fresh Claude Agent + Grok, launched on sight):**
```
$ grok --prompt-file _legs/S11c_c1_comparator_review_prompt.md --cwd /var/projects/toy_physics --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain   # detached, EXIT=0
# + a fresh general-purpose Claude Agent (agentId a42308ce…) executing the same prompt
```
BOTH CLEAN — the comparator is SOUND (computes/prints/decides-nothing; no false agreement, no hidden
disagreement, no broken seal, no dropped family; 30/30 tests). Grok: 7/7 CLEAN. Claude Agent: checks 1–4,6,7
CLEAN, check 5 CLEAN w/ two surfaced caveats. Both ran real ablations: one-sided corruption moved only the
corrupted row; repoint moved the residual; FACE↔DIRECTION merge DROPPED joins (156→60/8, collisions surface as
`UNDEFINED_DUPLICATE_KEY`, never new zeros); every seal's name-map ablation moved a residual (seals
load-bearing); `Inactive[Greater]`→`HeldInactiveGreater` (never a native bool).

**cS0 adjudication (rule 13, self-verified):** `python3 -c "import S11b_cross_engine_comparator as m;
print(m.mechanical_lower_camel('c_s0'))"` → `cS0`. So `cS0` IS the exact mechanical inverse of the real PY
`Symbol('c_s0')` — my directive's fold-1 claim (`=="cs0"`) was WRONG. The Agent's controlled synthetic proves
folding it collapses the residual to `Integer(0)` (pure spelling artifact). It fires in ZERO currently-joined
residuals (cS0 concentrates in the deferred PERMEABLE_* families; feasible joined leaves carry no cS0) and can
only over-report DISAGREEMENT — never a false agreement. LOW/latent, but a legitimate fold to add before the
deferred heavy run.

**Two surfaced INFO caveats (not soundness defects):** (a) WL `Inactive[Integrate][(…)]` /
`Inactive[Limit][Inactive[Integrate]…]` energy leaves fail to parse (56 ENERGY_RESIDUAL + 4
ENERGY_FACE_TRACTION) — surfaced honestly as `<PARSE_FAILED>` in UNPAIRED leaves, no residual corruption; (b) 4
giant families (PERMEABLE_PORT_HERMITIAN, PERMEABLE_DISSIPATION_VS_OMEGA_TAU, UNIFORM_LIMIT_S11CC1_OPERAND,
UNIFORM_LIMIT_RESIDUAL) deferred by size (OOM risk on 30 GB) — the one place cS0 could activate (unverified).

Leg outputs (outside repo): `<scratch>/grok_c1_comparator_review.txt`, agent transcript
`<session>/tasks/a42308ce….output`; leg ablation artifacts under `/tmp/s11cc1_t7_review/` and `/tmp/s11cc1_review/`.

**Verdict:** comparator SOUND + committed as the reviewed baseline. Scoped repair to follow (gated): R1 add the
`cS0←c_s0` mechanical fold (+ fix this directive's fold-1 cS0 paragraph); R2 parse the WL held-integral energy
leaves so the operand is displayed. No soundness change.

**Correction 2026-09-05:** this directive's fold-1 `cS0` paragraph (which wrongly claimed
`mechanical_lower_camel("c_s0")=="cs0"`) is now CORRECTED to `"cS0"` — a legitimate bare-symbol mechanical fold
to add. The correction is carried by the gated repair `S11c_c1_comparator_repair_directive.md`.
