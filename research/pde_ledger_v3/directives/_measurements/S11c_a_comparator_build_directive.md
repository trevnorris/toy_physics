# Measurements backing `S11c_a_comparator_build_directive.md` (rev 2 → rev 4)

Every factual claim in the rev-2 directive was verified against the real instruments/engine outputs this
session (2026-08-26), not against a paraphrase. Commands + literal output below (rule 2 binds the
orchestrator). Paths: repo root `/var/projects/toy_physics`; refs under `~/.s11_build/`.

## 1 — 39 joinable tag families (PY ∩ WL, minus `_LOCAL_`)
```
$ grep -oE '^PY_S11CA_[A-Z0-9_]+:' ~/.s11_build/S11c_a_py_fixed_run2.out | sed 's/^PY_S11CA_//;s/:$//' | sort -u > /tmp/py_tags.txt   # 47 tags
$ grep -oE '^WL_S11CA_[A-Z0-9_]+:' research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out | sed 's/^WL_S11CA_//;s/:$//' | sort -u > /tmp/wl_tags.txt   # 40 tags
$ comm -12 /tmp/py_tags.txt /tmp/wl_tags.txt | grep -vc '^LOCAL'
39
```
PY-only = the 8 `LOCAL_EXPORT_*`/`LOCAL_RUN_TASKS`/`LOCAL_SKIPPED_TASKS`/`LOCAL_TAGS`/`LOCAL_OPERATIONAL_EXCEPTIONS`/`LOCAL_SECTION8_REPORT`;
WL-only = `LOCAL_TAG_NAMES`. Both excluded from the join by §7. The 39: physics shape-derivs (12) +
projection (4) + origins (2) + controls (15) + supplied bookkeeping (4) + computed bookkeeping (2).

## 2 — `exports.py` is a `_LEDGER`, NOT the tag stream (rev-1's wrong PY input)
```
$ sed -n '25,30p' research/pde_ledger_v3/scripts/S11c_a_exports.py
_LEDGER = {
    'A':     { 'display': 'A', 'value': _restore("Symbol('A', real=True)"),
               'value_kind': 'COMPUTED_OBJECT', 'class': 'COORDINATE', 'step': 'S11', ... },
$ grep -cE 'PY_S11CA' research/pde_ledger_v3/scripts/S11c_a_exports.py
4    # (only inside BUILD_INPUT_DIGESTS / comments — zero S11CA tag rows)
```
The tag stream is emitted to **stdout** by the engine (`Emitter.emit`: `print(f"{tag}: {rendering}")`,
engine line 354); it is captured to a transcript, never stored in `exports.py`.

## 3 — WL committed `.out` is post-fix and multi-line
```
$ git log --oneline -1 -- research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out
6fae82b8 S11c-a T7: WL background-current fix — j0=rho0*v0, built not free-premise (2 build legs PASS)
$ grep -cE 'currentWBackground|currentXBackground' <that file>
0                                  # bg-current free premises removed
$ grep -noE '^WL_S11CA_[A-Za-z0-9_]+:' <that file> | head -3
1:WL_S11CA_CONTROL_FORM_BASE_OPERAND:
962:WL_S11CA_CONTROL_FORM_ABLATED_OPERAND:      # ⇒ first tag payload spans lines 1..961
2884:WL_S11CA_BACKGROUND_DENSITY_MAP:
```
40 WL tags in 2920 lines ⇒ multi-line reassembly required (max span ~961).

## 4 — PY `.out` provenance (why a fresh committed transcript, not run2 on trust)
`scripts/out/` holds committed PY engine transcripts up to 91 MB (`S11_stray_longitudinal_sympy_audit.out`);
a 50 MB `S11c_a` transcript is within precedent. No committed `S11c_a` PY `.out` exists yet.
```
$ sha256sum research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py | cut -c1-64
0c7053aa2fdeff4149d7e25f9a4b1dda04baeeb3c0e5e72021618a3ea32c38b0
$ grep -A1 "sympy_audit.py'" research/pde_ledger_v3/scripts/S11c_a_exports.py | grep -oE '[0-9a-f]{64}'
0c7053aa2fdeff4149d7e25f9a4b1dda04baeeb3c0e5e72021618a3ea32c38b0   # committed engine == engine that wrote committed exports.py
```
Two post-fix runs (`run.out` 11:09, `run2.out` 20:40; fix commit 20:46) are byte-identical on 36 of 39
joinable tags and differ on exactly 11 = `PROJECTION_*` (5) + `CONTROL_FORM_*` (3, ablates the projection) +
`UNIFORM_LIMIT_*` (3, references the projection) — a coherent ripple of the projection fix, confirming
determinism. To remove the 20:40→20:46 doubt, the committed transcript is REGENERATED from the committed
engine on a snapshot and confirmed to reproduce run2 (per-case symbolic diff 0); committed as scripts/out/S11c_a_interface_geometry_sympy_audit.out (afdc8158).

## 5 — `canon_key` silently drops axes (the case-loss defect → full-axis keying)
`~/.s11_build/S11c_a_scratch_loader.py` `canon_key(tokens)` adds to the key ONLY tokens matching
BRANCH (`LAB_HELD`/`MATERIAL_ADVECTED`), DENSITY (`RHO4_CONSTANT`/`RHOBR_CONSTANT`), FACE (`±1`),
DOF (`DELTA_W`/`ZETA_C`/`DOF_*`), VDOF (`VIRTUAL_DOF_*`); every other token (object, direction, origin,
field, unrecognized DOF) is **dropped**, so two distinct cases collapse to one `frozenset` key and one
overwrites the other in `pck={canon_key(k):v ...}`.

## 6 — S11b reader has NO multi-line reassembly; reusable vs run/render split
```
$ grep -nE '^def (read_transcript|parse_wolfram_payload|parse_sympy_payload|_split_top_level|transliterate|residual|classify_residual|extract_dimension_vector|compare_records|render_value|main)' research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py
158:def read_transcript(...)      # per-line TAG_LINE.fullmatch — NO reassembly
380:def parse_wolfram_payload   424:def parse_sympy_payload   279:def _split_top_level
550:def transliteration_collisions  585:def transliterate  682:def residual  859:def classify_residual
964:def extract_dimension_vector    # REUSE these
999:def compare_records  1122:def render_value  1185:def render_comparison_status   # DO NOT reuse (classify-first + FINAL_OPERATIONAL_STATUS)
```

## 7 — the two smuggles + the phantom
```
$ grep -n 'MUMAP' ~/.s11_build/S11c_a_cov_all.py
12:MUMAP={sp.Symbol('mu_theta_L'):sp.Symbol('mu_theta'),sp.Symbol('mu_theta_M'):sp.Symbol('mu_theta')}
   # applied inside cov_all zero(): a many-to-one BRANCH collapse (§3a μ_θ^α) → registry, not the map
$ grep -niE 'conormal|taylor' ~/.s11_build/S11c_a_cov_all.py ~/.s11_build/S11c_a_reconcile_fixed.py
(no output)                        # fold-9 "CONORMAL Taylor" does not exist in the refs — a phantom; CONORMAL compared raw
```

## 8 — perturbation-current spellings (the one physics-bearing new fold, #11)
```
$ grep -oE "Function\('delta_j_bulk_4'\)|delta_j_bulk_[123]" ~/.s11_build/S11c_a_py_fixed_run2.out | sort | uniq -c
   32 delta_j_bulk_1     32 delta_j_bulk_2     32 delta_j_bulk_3      # in-plane: bare constant Symbols
   16 Function('delta_j_bulk_4')                                      # normal: w-field (the 49b5c525 fix)
$ awk '/PROJECTION_DYNAMIC_OPERAND/{f=1} f{print}' <WL .out> | grep -oE 'current[WX]Perturbation' | sort | uniq -c
  208 currentWPerturbation   1440 currentXPerturbation                # WL: ALL components x/w-fields
```
⇒ the name fold aligns the token only; PY's in-plane constants vs WL's in-plane x-fields is an
arg-structure mismatch the comparator surfaces as a residual (candidate finding), not a fold to build.

## 9 — DIMENSIONS + BACKGROUND_DENSITY_MAP are COMPUTED (spec §8; engine `emit_primary(..., dimensions, ...)`)
Spec §8: "Computed here: every T-object in §4, including its multigrade and **dimension**"; T-0
"construct Σ_E⁰ and its in-plane gradient" (a computed T-object). Supplied (§8): `BACKGROUND_STATE`,
`ADMISSIBILITY_PREMISE`, `FACE_MAP_{LAB_HELD,MATERIAL_ADVECTED}` (engine `emit(...)` with no dimensions,
lines 1451–1454).

## 10 — no-drift reproduction (engine state matches, pre-directive)
```
$ python3 ~/.s11_build/S11c_a_reconcile_traction_check.py
RELATIVE_FLUX  join=8  zero=8  nonzero=0
TRACTION       join=16 zero=16 nonzero=0
$ python3 ~/.s11_build/S11c_a_run_evolution.py
EVOLUTION_MASS_BALANCE  join=8 zero=8 nonz=0 err=0
```
(These 3 families are among the 36 byte-identical across run/run2, so they do not distinguish the runs; the
projection provenance is settled by the regeneration in §4.)

## 11 — per-family extractor leaf paths (rev-3; probed from the committed transcripts)
```
$ python3 (scratch_loader py_cases + wl subfields on the committed PY .out + WL .out):
KINEMATIC_BALANCE   WL {OPERAND_A_SHAPE_DERIVATIVE, OPERAND_B_SHAPE_DERIVATIVE, RESIDUAL_A_MINUS_B, BOUND_SOURCE_LAW}   PY VALUE=3-tuple (OPERAND_A, OPERAND_B, RESIDUAL)
VIRTUAL_CONSTRAINT  WL {NORMALIZED_VIRTUAL_MASS_VARIATION, EXACT_EW_BG_MAP, VIRTUAL_DISPLACEMENT_OPERAND}            PY VALUE=single Mul (eps*expr)
VIRTUAL_WORK        WL {EXACT_TRUE_AREA_SOURCE, SHAPE_DERIVATIVE, FACE_DISPLACEMENTS_FROM_MAP}                        PY VALUE=2-tuple keyed UPPER/LOWER
CONORMAL_DERIV      WL {GRAPH_EVALUATED_SOURCE, SHAPE_DERIVATIVE}                                                     PY VALUE=3-tuple (bg,deriv,total)
PROJECTION_SHAPE_DERIV WL {DYNAMIC_WINDOW, EXACT_PROJECTED_IDENTITY_ZERO_FORM, SHAPE_DERIVATIVE}                      PY VALUE=Mul (eps*Integral(Subs(Derivative(O_window...))))
EVOLUTION_TERM_ORIGINS WL {EXACT_SOURCE_TERM, SHAPE_DERIVATIVE} per origin case                                      PY VALUE=named partition map (DENSITY_TIME, VELOCITY..., ...)
PROJECTION_TERM_ORIGINS WL {DYNAMIC_SOURCE_TERM, DYNAMIC_SHAPE_DERIVATIVE, STATIC_SOURCE_TERM, STATIC_SHAPE_DERIVATIVE} PY VALUE=2-tuple (DYNAMIC map, STATIC map)
```
⇒ each family's leaf is bespoke; a single "SHAPE_DERIVATIVE.EXPRESSION" rule (rev-2) was wrong for
KINEMATIC/VIRTUAL_CONSTRAINT/VIRTUAL_WORK/projection/origins.

## 12 — keying reality (why CORRECTION 1 EXTENDS canon_key, not replaces it)
`canon_key` axis-typing joins the SIMPLE families (FACE_NORMAL 8/8, etc.) because it maps PY `DELTA_W`→(DOF,·),
WL `DOF_DELTA_W`→(DOF,·), face `±1`→(FACE,·). It FAILS where extra axes exist: it dedups the two positional
PY DOF/VDOF tokens into a frozenset (VIRTUAL_WORK), and conflates CONTROL_FORM's face `±1` with its direction
`±1`, and drops FIELD (FACE_SHIFT)/OBJECT (HOMOGENEITY). rev-2's "full tuple, face-only" REGRESSED even the
simple join to 0/8 (PY `DELTA_W` ≠ WL `DOF_DELTA_W` without the typing). rev-3 keeps the typing, extends it to
all axes, and preserves DOF/VDOF order/role.

## 13 — the perturbation-current fold hardening (two-advisor consult, 2026-08-27)
Both advisors (Codex + Grok), independently, on the projection residual: adjudicating it on the orchestrator's
unreviewed scouting tools presupposes the measurement (the inversion that hid #3); build the REVIEWED
comparator first. The current fold must be a literal AST-head rename preserving args/arity (⛔ NOT the FIELD
`AppliedUndef→bare Symbol` strip); emit the RAW residual + a separately-labeled held-context diagnostic
`J(w; held=x,t)`, ⛔ never assert either zero, ⛔ never count the reduced-zero as agreement, ⛔ never put the
held-context form in the closed map — its authorization is a from-spec adjudication (§1b defines a full 4D
current; "WL never differentiates x" is evidence, not proof). The integral canon must RETAIN bounds + binder
identity (the scouting `ITG` dropped them). The `delta_rho_4D_bulk_t` density-time term gets equal billing.
Consult prompts+outputs: `~/.s11_build/path_decision_{codex.log,grok.txt}`; SCOUT §23.
