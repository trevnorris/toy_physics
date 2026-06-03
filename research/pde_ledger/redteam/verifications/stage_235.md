---
unit_id: 235
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T22:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 235

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created a NEW second-engine script
`mathematica/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_mathematica_audit.wl`
(173 lines, untracked `??`, mtime 22:14:24). No other file was touched — `git diff HEAD`
on the `.py` is empty (mtime still 2026-05-11 11:56). The captured `stage_235_diff.patch`
is empty because the only change is a brand-new untracked file.

**Assessment:**
The `.wl` independently re-derives all six Stage-235 deliverables and exits 0 with 31
PASS lines, matching the 31 `expect*` calls in the script (M1=5, M2=7, M3=4, M4=6, M5=4,
M6=5). The fix correctly and fully addresses the sole finding (absent second engine).

- M1 (involution): `MatrixPower[mouthCompiler,2] - IdentityMatrix[2]` (line 80, native
  primitive, not SymPy's `M_rm*M_rm`), compiler map, inverse map, `Det+1`. Non-vacuous.
- M2 (projectors): builds via the directive-mandated similarity `M_rm.Diag.M_rm`
  (lines 92-93) AND — as the directive and hard-rule #1 require — independently confirms
  the defining projector PROPERTIES directly: `P_nt.P_nt==P_nt`, `P_eta.P_eta==P_eta`,
  `P_nt.P_eta==0`, `P_eta.P_nt==0`, `P_nt+P_eta==I` (lines 98-102). Closed forms checked
  against the paper targets `{{1,c},{0,0}}`, `{{0,-c},{0,1}}`. Non-vacuous.
- M3 (decomposition): `P_nt.x`, `P_eta.x`, `-Xi1` form, and recomposition. Non-vacuous.
- M4 (codim-two): native `Solve[...,Reals]` with structural `===` against `{{R1->0,E1->0}}`
  for both `{Xi1==0,E1==0}` and `{Xi1==0,R1==0}`, plus `Det!=0` boolean and a forward
  `M_rm.{0,0}` direction check, plus an EXTRA codim-one-vs-two distinction via
  `Length[Solve[{Xi1==0},{R1}]]==1` (line 132) that the `.py` does not have. Non-vacuous.
- M5 (static-blind line + norm): `Xi1` on the line = 0, `x.x == (1+c^2)q^2`, and the
  `q_eta=L/Sqrt[1+c^2] => norm^2==L^2` scaling, plus `Xi1` still 0 at size `L`. Non-vacuous.
- M6 (correction compilers): `Delta_x_static`, `Delta_x_orbit`, additive relation, and the
  full-lock identity `x_q + Delta_x_orbit == 0`. Non-vacuous.

The fail path is real: `fail` calls `Exit[1]` (lines 20-24); `expectZero`/`expectMatrixZero`
run `cleanExpr` (FullSimplify + ConditionalExpression strip) and demand strict `===0` /
`===ConstantArray[0,…]`; `expectTrue` demands `TrueQ`. A wrong projector, wrong lock
solution, or wrong norm leaves a nonzero residual and aborts with exit 1. Each LHS is
constructed independently and compared to the paper's stated closed form or zero — no
`x = expr; assert x == expr` tautology. The project idioms are honored: `stripConditional`
removes `ConditionalExpression[0,…]` before comparison (lines 26, 57); `Solve[...,Reals]`
is compared structurally rather than via a fragile `=!= Infinity` pole test.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `det(M_rm) = -1`
- `||x_eta||^2 = qeta_line**2*(c_eta**2 + 1)`
- `All Stage 235 symbolic checks passed.`
(The `.py` was not modified; this is the unchanged baseline engine, re-run for cross-check.)

**Mathematica:** exit=0. Notable lines:
- `M1 MatrixPower[M_rm,2] - I residual = {{0, 0}, {0, 0}}` → PASS
- `Solve[Xi1==0,E1==0] = {{R1 -> 0, E1 -> 0}}` → PASS (codim-two unique origin)
- `x_blind.x_blind = (1 + cEta^2)*qEta^2` → PASS
- `M6 full-lock identity residual = {0, 0}` → PASS
- `All Stage 235 Mathematica checks passed.` `# exit_code: 0`
No silent parser-skip: every `expect*` printed a residual line immediately followed by PASS;
31 PASS lines = 31 `expect*` calls. Both engines agree on `M_rm^2=I`, the projector closed
forms, the codim-two solution set `{(0,0)}`, and `‖x_η‖²=(1+c_η²)q_η²`.

**Output freshness:** confirmed. Saved outputs were regenerated post-fix:
- SymPy `.txt` mtime 2026-06-02 22:16:48 (script mtime 2026-05-11 11:56) — newer.
- Mathematica `.txt` mtime 2026-06-02 22:16:48 (script mtime 2026-06-02 22:14:24) — newer.
Both `.txt` files match their exec logs verbatim.

## Material-change assessment

`material_change`: false.

Only a second engine (`.wl`) was added; the `.py` is byte-for-byte unchanged vs HEAD and no
derived result changed. The Mathematica route confirms the exact same closed forms the
`.py` (and the paper) already stated. No downstream unit (>235) depends on a changed value,
so no `upstream_stale` propagation is warranted on substance grounds.

## Side observations (non-blocking)

- The notes file self-labels "Stage 251/252/253" — the known EM-extension renumber drift;
  canonical 235 is ground truth per the audit and the orchestrator directive. Not a finding.
- The `.wl` and its output `.txt`, like the sibling stage-231..236 second engines, are still
  untracked (`??`) — expected for a fresh batch, will be committed at batch close. Not a defect.

## Verdict justification

The sole finding (missing second engine) is fully resolved: a genuinely independent
Mathematica `.wl` now covers all six Stage-235 deliverables (M1–M6) with non-tautological,
hard-failing checks, exits 0, and its outputs are fresh. Independence is confirmed — the
script uses native primitives (`MatrixPower`, `Solve[...,Reals]` with structural `===`, an
extra codim-one-vs-two `Length` count, a forward origin-direction check) and, per the
directive, confirms the projector idempotency/complementarity/codim-two-lock PROPERTIES
directly rather than merely echoing the similarity product with renamed variables. The `.py`
is untouched (audit found it aligned), so this is a non-material change. Verdict: verified.
