---
unit_id: 089
batch: III.5
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 089

## Per-finding outcomes

### F1 — paper_misalignment

**Classification:** resolved

**What changed:**
User chose direction (a) "strengthen scripts" per `redteam/resolutions/batch_III5_paper_alignment.md` Q3.

- `scripts/.../stage089_..._sympy_audit.py:53-56` adds `Omega_at_zero = sp.simplify(sp.limit(Omega, Pe, 0))` and `zeta_F1_at_zero = sp.simplify(sp.limit(zeta_F1, Pe, 0))`, then `expect_zero("Omega(Pe -> 0) - 1", ..., tol=1e-30)` and `expect_zero("zeta_F1(Pe -> 0) - A_F1", ..., tol=1e-30)`.
- `scripts/.../stage089_..._sympy_audit.py:128-129` adds `Pe_req = sp.Integer(0)` plus `expect_zero("Pe_req (zero-bias bound from chain closure)", Pe_req)` at the end of the regime checks, closing the chain to the paper's boxed Output.
- Mathematica mirror at `wl:51-54` (`omegaAtZero = Limit[omegaPe[pe], pe -> 0]`, `zetaF1AtZero = Limit[zetaF1[pe], pe -> 0]`, two `expectApprox` calls with tol `10^-30`) and `wl:112-113` (`peReq = 0; expectApprox["Pe_req ...", peReq, 0, 10^-30]`).

**Assessment:**
Edit closes the auditor-named missing link `Omega(Pe -> 0) = 1 ⇒ zeta_F1(0) = A_F1 ⇒ Pe_req = 0`. Both engines compute the non-trivial 0/0 limit and assert it. The `Pe_req = sp.Integer(0)` / `peReq = 0` assignment is unconditional in this stage (per resolution Q3 note about pitfall #10 avoiding `sp.nsolve`), which matches the resolution doc's prescription. Banner-prose comments at `py:48-52` and `wl:46-50` reference paper eq `app-stage089-Pe-zero` correctly. The new assertions are non-tautological — sympy log shows the `zeta_F1(0) - A_F1` residual is `-2.17e-101` (not literally zero), confirming the limit is a real symbolic computation rather than a definitional substitution.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
At `wl:63-66`, the literals `peSuffChi = SetPrecision[96.5285247264386, 40]` and `peFailChi = SetPrecision[11220.5441626259, 40]` are replaced with:

```
zetaSuffTarget = SetPrecision[3.46622291347846 - 1, 40];
zetaFailTarget = SetPrecision[3.46752913273870 - 1, 40];
peSuffChi = pe /. FindRoot[zetaF1[pe] == zetaSuffTarget, {pe, 100}, WorkingPrecision -> 40];
peFailChi = pe /. FindRoot[zetaF1[pe] == zetaFailTarget, {pe, 10000}, WorkingPrecision -> 40];
```

A comment block at `wl:56-62` documents the independence rationale and notes that SymPy retains the literals (F4 path) because `sp.nsolve` is unstable near the tan(y) singularity of stage-074 closed form.

**Assessment:**
The Mathematica side now derives `peSuffChi` and `peFailChi` independently from notes-quoted `rho_target - 1` via `FindRoot`, not by transcribing the SymPy literals. The numeric agreement with SymPy (rho_suff = 3.46622291347846..., rho_fail = 3.46752913273870...) is now genuine cross-engine confirmation rather than mirrored input. Auditor's required-change spec is met: literals `96.5285247264386` and `11220.5441626259` no longer appear in the `.wl` file (confirmed via Read). Engine-independence policy met for the checkpoint bar.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
- `py:80-92` introduces a proper `expect_close(name, value, target, tol)` helper and replaces the prior `expect_zero("rho_X - (1 + zeta_X)")` lines with three `expect_close` calls comparing `rho_suff`, `rho_fail`, `rho_max` against the notes-quoted Stage-082 literals (`3.46622291347846`, `3.46752913273870`, `3.46752922945601`) at tol `1e-12`.
- `wl:79-81` mirrors with three `expectApprox` calls against the same notes-quoted targets at tol `10^-12`.

**Assessment:**
Non-tautological. The sympy log shows the residuals `4.001e-15`, `3.486e-15`, `2.233e-15` — real numerical differences caused by the 12-significant-figure carry-forward of upstream literals, exactly as the auditor's verification description predicted ("`~1e-15` rather than the trivial `0`"). Mathematica residuals are similar (`0.`, `0.`, `2.220e-15`). The `Q(zeta;0) = 1+zeta` line at `py:77` is retained as documentation, as the directive's parenthetical allowed.

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
SymPy retains the literals at `py:65-66` per the directive's "provenance-comment path" (option b), with the new comment block at `py:58-64` naming `scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt` and Stage 089 notes section 1 as upstream sources, plus an explicit reference to `notes/STAGE_VERIFICATION_COVERAGE.md pitfall #10` justifying why `sp.nsolve` was not used. The Mathematica side resolves F4 via the F2 FindRoot rederivation.

**Assessment:**
Matches the directive's allowed alternative — the resolution doc Q3 explicitly flags pitfall #10 risk for `sp.nsolve` on stage-074 closed forms, so the SymPy provenance-comment path is the prescribed route. Mathematica engine independently re-derives via `FindRoot` (F2). The combined result: both engines anchor the load-bearing Pe values to upstream notes via different mechanisms (provenance on sympy side, independent rederivation on Mathematica side). The auditor's verification criterion ("either no bare literals... or a comment block naming the upstream source") is satisfied on the sympy side.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Omega(Pe -> 0) - 1 = 0` (chain-closure check passes)
- `zeta_F1(Pe -> 0) - A_F1 = -2.17e-101` (non-trivial symbolic limit, residual at working precision floor)
- `rho_suff vs Stage-082 quote diff = 4.001e-15` (non-tautological upstream anchor)
- `Pe_req (zero-bias bound from chain closure) = 0` (boxed Output verified)

**Mathematica:** exit=0. Notable lines:
- `Omega(Pe -> 0) - 1 diff = 0` → `PASS`
- `zeta_F1(Pe -> 0) - A_F1 diff = 0\`\`49.69` → `PASS` (50-digit precision)
- `rho_suff vs Stage-082 quote diff = 0.` → `PASS` (FindRoot-derived Pe_suff_chi matches notes literal exactly at the 12-sig-fig precision of the target)
- `rho_max vs Stage-082 quote diff = 2.220e-15` → `PASS`
- `Pe_req (zero-bias bound from chain closure) diff = 0` → `PASS`
- Final: `Stage 089 Mathematica audit passed.`

**Output freshness:** confirmed. sympy script mtime `1779898875` < output mtime `1779899095`; mathematica script mtime `1779898920` < output mtime `1779899191`. Both transcripts re-generated post-fix.

## Material-change assessment

`material_change`: false.

The numeric values produced (`rho_suff`, `rho_fail`, `rho_max`, `A_F1`, `zeta_max`) are unchanged from the pre-fix run — the fixes added new assertions (chain-closure limits and non-tautological upstream anchors) but did not change any derived quantity that downstream stages would consume. `Pe_req = 0` is a new explicit assertion but matches the paper's pre-existing boxed value. No downstream re-audit required.

## Side observations (non-blocking)

1. The directive Codex was supposed to append `## Applied: F<n>` blocks under each finding (per directive header instruction at line 14). No such blocks are present in `redteam/directives/stage_089.md`; only the front-matter status changed to `applied: true`. The orchestrator's `verification_status: scripts_pass_pending_verifier` flag conveys the same information, so this is a non-blocking process observation.
2. The sympy script docstring at `py:3` still reads `moving_throat_pde_stage72_family1_minimal_isotropic_verdict_sympy_audit.py` (stale pre-renumber filename). The resolution doc at line 160 indicated docstrings would be relabelled orchestrator-direct alongside banners; the banner relabel was applied at `py:29` (`banner("STAGE 089 — ...")`) but the docstring was missed. Cosmetic only — does not affect verification. Banner sweep on `wl` side is correctly applied at `wl:32`.

Neither observation blocks verification.

## Verdict justification

All four findings are `resolved`. F1 (paper_misalignment, user-resolved) closes the link chain `Omega(Pe -> 0) = 1 ⇒ zeta_F1(0) = A_F1 ⇒ Pe_req = 0` end-to-end in both engines. F2 introduces an independent Mathematica derivation of `peSuffChi`/`peFailChi` via FindRoot from notes-quoted `zeta_target`, breaking the syntactic-port pattern. F3 replaces the tautological `rho_X - (1 + zeta_X) == 0` checks with `expect_close` cross-checks against upstream Stage-082 literals, producing non-zero `~1e-15` residuals that real numerical comparison would yield. F4 follows the directive's provenance-comment path on the SymPy side (consistent with pitfall #10) and is resolved on the Mathematica side via F2. Both engines exit 0 with output mtimes newer than script mtimes. Checkpoint bar (independent extraction, non-tautological anchors, chain closure to boxed Output) is met.
