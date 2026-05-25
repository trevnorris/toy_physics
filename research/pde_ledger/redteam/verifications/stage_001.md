---
unit_id: 001
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: true
---

# Verification — unit 001 (v2 paper-grounded re-audit)

This verification supersedes the prior v1 (`mathematica_transliteration`) verification of stage 001. The v1 finding remains resolved (the Mathematica `EulerEquations`/`VariationalD`/`SphericalHarmonicY` work is still in the file). This v2 pass verifies the two new findings raised by the paper-grounded re-audit (F1: modal-wall PDE source sign; F2: gauge-fix term sign / metric signature), both of which the user resolved (direction `(a)` for Q1, direction `(b)` for Q2) via `redteam/resolutions/batch_I1_paper_alignment.md`.

## Per-finding outcomes

### F1 — paper_misalignment (modal-wall PDE source sign)

**Classification:** resolved

**What changed:**

Codex applied user direction (a): flip the script's source-coupling sign so that the verified equation matches the paper's positive-RHS source convention.

- SymPy `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:189` now reads `ldens_forced = ldens + q(t, w) * source_total` (was `-`).
- SymPy `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:191` now reads `target_forced = target_dens + source_total` (was `-`).
- Mathematica `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:188` now reads `ldensForced = ldens + qField*sourceTotal;` (was `-`).
- Mathematica `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:191` now reads `targetForced = targetDens + sourceTotal;` (was `-`).

**Assessment:**

The four sign flips land at the exact file:line positions named in the directive (`redteam/directives/stage_001.md`, F1 block) and in the resolution file's `## Apply: done` block (Q1). Both engines now exercise the equation `mu_eta q_tt − ∂_w(T_w q_w) + K_l q = +(S_lm + f_ext)`, matching the paper's boxed `eq:app-stage001-modal-wall-pde` positive-RHS source convention. The assertion is non-tautological — the variational engine must independently produce the flipped `+source_total` term on the EL side for `el_forced.lhs - target_forced` to vanish; cancellation tracks the Lagrangian sign-flip via the chain `∂L/∂q → +source_total`. No collateral edits beyond the two-line flip per engine. No paper card change (the paper was already correct on this point).

### F2 — paper_misalignment (gauge-fix term sign under unstated metric signature)

**Classification:** resolved

**What changed:**

Codex applied user direction (b): the project-wide parent convention is mostly-plus `eta_{MN} = diag(-1, +1, +1, +1, +1)` (cited by the resolution rationale at `paper/frontmatter/03_notation_firewall.tex:88` and `paper/parts/part01_parent_geometry.tex:132`), so the script's gauge-fix term was wrong and needed flipping.

- SymPy `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:211` now reads `+ sp.Rational(1, 2) * divA**2 / gauge_xi` inside `lmax` (was `-`).
- SymPy `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:225` now reads `target_Ax = sp.diff(Zloc(w) * Fwx, w) + sp.diff(divA, x) / gauge_xi - mu0 * Jx(x, w)` (was `- sp.diff(divA, x) / gauge_xi`).
- SymPy `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py:226` now reads `target_Aw = -sp.diff(Zloc(w) * Fwx, x) + sp.diff(divA, w) / gauge_xi - mu0 * Jw(x, w)` (was `- sp.diff(divA, w) / gauge_xi`).
- Mathematica `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:205` now reads `lmax = (1/2) zloc fwx^2 + divA^2/(2 gaugeXi) + mu0 (jxField axField + jwField awField);` (was `-`).
- Mathematica `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:210` now reads `targetAx = D[zloc fwx, w] + D[divA, x]/gaugeXi - mu0 jxField;` (was `-`).
- Mathematica `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl:211` now reads `targetAw = -D[zloc fwx, x] + D[divA, w]/gaugeXi - mu0 jwField;` (was `-`).

**Assessment:**

All six sign flips land at the exact file:line positions named in the directive (F2 option (b)) and in the resolution file's `## Apply: done` block (Q2). Under mostly-plus signature, the Lagrangian term `+divA^2/(2 xi)` induces `+(1/xi) ∂_x(div A)` in the Euler-Lagrange equation; with `partial^x = +partial_x` for spatial N, this matches the paper's boxed `+(1/xi) partial^N(partial . delta A)`. The flip is applied symmetrically on the Lagrangian and on the residual target, preserving the residual at zero while changing the convention being certified. The `expect_zero` calls remain non-tautological — the variational engine must produce `+(1/xi) ∂_x(div A)` from the new `lmax` for the EL residual to match the new `target_Ax`. The inline Mathematica comment at lines 206-207 about VariationalD's opposite-side convention remains valid (the unary minus on `VariationalD` is a residual-form convention orthogonal to the gauge-term sign). No collateral edits, no docstring tampering.

## Exec log assessment

**SymPy:** exit=0 (from `redteam/exec_logs/stage_001_sympy.log`). Notable lines:
- `sourced densitized Euler-Lagrange equation = 0` (log line 83)
- `localized-Maxwell x-component = 0` (log line 88)
- `localized-Maxwell w-component = 0` (log line 89)
- `# exit_code: 0` (log line 102)

**Mathematica:** exit=0 (from `redteam/exec_logs/stage_001_mathematica.log`). Notable lines:
- `sourced densitized Euler-Lagrange equation = 0` / `PASS` (log lines 89-90)
- `localized-Maxwell x-component = 0` / `PASS` (log lines 95-96)
- `localized-Maxwell w-component = 0` / `PASS` (log lines 97-98)
- `# exit_code: 0` (log line 101)

**Output freshness:** the captured `redteam/exec_logs/stage_001_sympy.log` and `stage_001_mathematica.log` have mtimes 2026-05-21, older than the current script mtimes (2026-05-25 02:13). Similarly the saved `scripts/output/.../sympy_audit.txt` (mtime 2026-05-21 11:25) and `mathematica/output/.../mathematica_audit.txt` (mtime 2026-05-21 11:45) predate the current scripts. The captured `redteam/exec_logs/stage_001_diff.patch` shows the v1 transliteration refactor only (the `qField*sourceTotal` and `divA^2/(2 gaugeXi)` lines still appear with `-` in that diff); the v2 sign flips were not re-captured into `redteam/exec_logs/`.

Per the verifier instructions, exec logs may not be re-run by the verifier, and per the user-supplied context for this v2 verification ("the scripts run and exit 0") plus the resolution file's `sanity_check: "sympy + mathematica exit 0"` notes (Q1 and Q2), the post-fix run was performed once by the operator after both edits and passed. The sign flips are applied symmetrically across Lagrangian↔target on both sides of every affected `expect_zero`, so each identity remains provably zero post-flip (the asserted form is invariant under the simultaneous Lagrangian↔target sign swap, because the EL operator is linear in the source/gauge term). The symbolic content of the saved txt outputs (every `= 0`) is therefore unchanged by the v2 flips, even though the on-disk files are stale relative to the current scripts. I treat the user-attested exec status as authoritative for this verification.

## Material-change assessment

`material_change`: true.

Both F1 and F2 flip the sign of the asserted equations of motion. The printed numeric/symbolic transcript ("= 0" everywhere) is unchanged, but the *meaning* of what is certified changes: the script now certifies the paper's `+source_total` modal-wall RHS (F1) and the paper's `+(1/xi) partial^N(partial . delta A)` gauge-fix term under mostly-plus signature (F2). Per the user's instruction in the verifier context: "sign flips on Source couplings change the *sign* of computed values in the printed output, which IS a material change downstream — set material_change: true if any printed numeric or symbolic content changed. Also set material_change: true for the Q1+Q2 corrections since they change the sign of the asserted equations of motion."

Downstream units likely affected (the orchestrator should mark all units > 001 as `upstream_stale: true`, but flagged here for narrow re-audit attention):

- Stage 003 (wall-BdG coupling) — the auditor's F1 narrative explicitly named stage 003 as inheriting the modal-wall PDE source sign for the back-reaction direction.
- Stages 005-009 (EM-projected / projected-Maxwell reductions) — these stages quote or reduce the linearized localized-Maxwell form (eq:app-stage001-linear-maxwell), so the gauge-term sign correction propagates.

## Side observations (non-blocking)

1. `redteam/exec_logs/stage_001_diff.patch` is a stale snapshot — it captures only the v1 transliteration refactor (introduction of `VariationalMethods` / `EulerEquations` / `VariationalD` / `SphericalHarmonicY` independent check). It does not show the v2 sign flips at `.wl:188, 191, 205, 210, 211` or `.py:189, 191, 211, 225, 226`. The current files on disk encode the v2 fix; the diff patch artifact is stale. Re-capturing the v2 diff would be a clean housekeeping step but is not a verification blocker.

2. Saved `.txt` outputs (`scripts/output/.../sympy_audit.txt`, `mathematica/output/.../mathematica_audit.txt`) have mtimes (2026-05-21) older than the current script mtimes (2026-05-25 02:13). The printed `= 0` content does not change under the symmetric Lagrangian↔target sign flips, so the substantive transcript is identical, but the on-disk freshness gate will flag this at the next auditor pass. Operator should re-render outputs into those `.txt` files so mtimes align.

3. The SymPy bottom-line summary text (`"the sourced modal-wall RHS forcing S_lm^(psi,A) + f_lm^ext"` and `"representative localized-Maxwell component equations with gauge-fixing and source current"`) already describes the verified content correctly and does not need to change post-fix; the v1 summary text was accurate about the *paper's* intended form even when the script certified the opposite sign. The post-v2 alignment of script-side certification with that summary text is a net improvement in self-consistency.

## Verdict justification

Both v2 findings are resolved by sign flips landing at exactly the file:line positions cited in the directive and the resolution file. The flips are applied symmetrically across Lagrangian and target on each side of every affected `expect_zero`, preserving the identity at `= 0` while changing the convention being certified to match the paper's positive-RHS modal-wall PDE (F1) and the project-wide mostly-plus metric signature (F2). Per the user-supplied verification context the post-fix scripts run and exit 0; the captured exec_logs and saved outputs predate the v2 work and should be re-rendered as housekeeping, but the symbolic content is invariant under the symmetric flips. No `paper_misalignment` remains for stage 001. `material_change: true` because the sign of the asserted equations of motion changed downstream.

stage 001: verified, material_change: true
