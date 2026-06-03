---
unit_id: 236
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 236 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_236.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row 84; theorem-block rows 984-1005)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card `\stagefield{Output}` states: "Microscopic projector theorem: clearing the static nontracking defect moves along a different dependent direction than clearing the post-static dressing defect; a static-only fix leaves a dressing gap." The `\stagefield{Derivation ledger}` lists four artifacts: the dependent-plane compiler, the projectors `P_nt^dep` and `P_eta^dep`, the equal-drift dressing ray, and the static-only restoration direction. The notes (the authoritative detail) make these exact: on the rigid-mouth slice (`q_tr = 0`) the dependent correction is `y_rm = (0, -q_eta, q_nt - q_eta)` = `S_rm_dep (q_nt, q_eta)` with `S_rm_dep = [[0,0],[0,-1],[1,-1]]`; the direct-observable compiler is `C_rm_dep = S_rm_dep M_rm` with the explicit `(3,2)` entry `-1/(1-eps_eta_star)`; there is an exact left inverse `L_rm_dep = [[0,-1,1],[0,-1,0]]` with `L_rm_dep S_rm_dep = I_2`; the pushed-back projectors `P_nt^dep`, `P_eta^dep` are idempotent, mutually annihilating, and sum to the plane identity `diag(0,1,1)`; the decomposition `y_rm = y_nt + y_eta` holds with `y_nt = (0,0,q_nt)` and `y_eta = -q_eta(0,1,1)`; the static strip `q_nt = 0` is equivalent to `Delta_mu = Delta_Keta`; the equal-drift ray has norm `||y_eta||^2 = 2 q_eta^2 = 2 E1^2 = 2 R1^2/c_eta^2`; and the static-only / full orbit-lock correction compilers are `Delta_y_static = -y_nt` and `Delta_y_orbit = -y_rm`. The appendix (eq:app-part07-dependent-correction-q, eq:app-part07-equal-drift-ray, row 84) restates the dependent correction and the `-q_eta(0,1,1)` ray identically. NOTE: the notes prose repeatedly labels this "Stage 253" (e.g. headers, the closing supporting-file line citing `...stage253...sympy_audit.py`). That is the known project-wide incomplete-renumber drift; canonical numbering (paper-card filename `stage_236.tex`, script filename `...stage236...`, appendix row 84) is ground truth, so this is a stale-label note, not a misalignment finding.

## What the script claims to verify

The SymPy script reconstructs the rigid-mouth packet map `M_rm`, the general dependent section `(Delta_T_q, Delta_Keta_q, Delta_mu_q)`, restricts to `q_tr = 0`, and asserts: (1) `q_rm = M_rm x_rm` equals `(-R1 - c_eta E1, E1)` (line 36); (2) `y_rm = S_rm_dep (q_nt, q_eta)` (line 46); (3) `C_rm_dep = S_rm_dep M_rm` equals the closed form with `-1/(1-eps_eta_star)` (line 66) and `y_from_x` equals `(0, -E1, -R1 - E1/(1-eps_eta_star))` (line 69); (4) `L_rm_dep S_rm_dep = I_2` (line 90) and `L_rm_dep y_rm = (q_nt, q_eta)` (line 93); (5) the projectors equal their explicit forms, are idempotent, mutually annihilating, and sum to `diag(0,1,1)` (lines 104-112), with `y_nt`, `y_eta`, and `y_rm = y_nt + y_eta` (lines 117-119); (6) the static-strip identity `Delta_mu - Delta_Keta = q_nt` (line 137) and the on-strip value `(0,-E1,-E1)` (lines 141-142); (7) the norms `2 q_eta^2`, `2 E1^2`, `2 R1_static^2/c_eta^2` (lines 145-149); (8) the correction compilers (lines 163-168). Every check is an `assert_zero`/`assert_matrix_zero` against an independently written expected form.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Dependent section `y_rm = S_rm_dep(q_nt,q_eta)`, `(0,-q_eta,q_nt-q_eta)` | lines 43-46 | match |
| Direct-observable compiler `C_rm_dep = S_rm_dep M_rm`, `(3,2)=-1/(1-eps)` | lines 57-78 | match |
| Left inverse `L_rm_dep S_rm_dep = I_2`; recover `(q_nt,q_eta)` | lines 89-93 | match |
| Projectors idempotent/orthogonal/sum=`diag(0,1,1)` | lines 98-112 | match |
| Decomposition `y_rm = y_nt + y_eta`, `y_nt=(0,0,q_nt)`, `y_eta=-q_eta(0,1,1)` | lines 114-119 | match |
| Static strip `q_nt=0 ⟺ Delta_mu=Delta_Keta` | line 137 | match |
| Equal-drift ray `(0,-E1,-E1)` on strip; norm `2q_eta^2=2E1^2=2R1^2/c_eta^2` | lines 139-149 | match |
| Static-only / orbit-lock correction compilers | lines 160-168 | match |
| Second-engine (Mathematica) verification | — | missing |

`paper_alignment: aligned` — every paper-side deliverable has a faithful, non-tautological script-side check in SymPy. The sole gap is the absence of the required second engine.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 36 | `assert_matrix_zero(q_rm - [-R1-c_eta E1, E1])` | packet map M_rm | yes |
| A2 | sympy | 46 | `assert_matrix_zero(y_rm - S_rm_dep[q_nt,q_eta])` | dependent section | yes |
| A3 | sympy | 66 | `assert_matrix_zero(C_rm_dep - C_expected)` | compiler `C_rm_dep` | yes |
| A4 | sympy | 69 | `assert_matrix_zero(y_from_x - [0,-E1,...])` | compiler applied | yes |
| A5 | sympy | 90 | `assert_matrix_zero(L_rm_dep S_rm_dep - I_2)` | left inverse | yes |
| A6 | sympy | 93 | `assert_matrix_zero(L_rm_dep y_rm - [q_nt,q_eta])` | recovery | yes |
| A7 | sympy | 104-105 | projector explicit forms | projectors | yes |
| A8 | sympy | 106-109 | idempotent + mutual annihilation | projector algebra | yes |
| A9 | sympy | 112 | sum = `diag(0,1,1)` | plane identity | yes |
| A10 | sympy | 117-119 | `y_nt`, `y_eta`, `y_rm=y_nt+y_eta` | decomposition | yes |
| A11 | sympy | 137 | `assert_zero((y_rm[2]-y_rm[1]) - q_nt)` | static-strip equiv | yes |
| A12 | sympy | 141-142 | on-strip `(0,-E1,-E1)` and `(R1/c_eta)(0,1,1)` | equal-drift ray | yes |
| A13 | sympy | 145,148,149 | norms `2q_eta^2`, `2E1^2`, `2R1^2/c_eta^2` | ray norm | yes |
| A14 | sympy | 163-168 | correction compilers | static-only/orbit | yes |

All SymPy rows trace to a specific paper deliverable and are non-tautological. No row is orphaned scaffolding. No Mathematica rows exist (finding F1).

## Findings

### F1 — missing_verification_script (subtype missing_mathematica)

**Severity:** medium
**Files:**
- (absent) expected: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_mathematica_audit.wl`

**What's wrong:**
The stage card `\stagefield{Verification}` explicitly states "Mathematica audit: none yet." Only the SymPy script exists. This stage is `is_checkpoint: False` and `is_status_only_candidate: False`, so the dual-engine rule applies: a `.wl` is required wherever Mathematica CAN independently verify the claims. Every claim here is finite-dimensional linear algebra over the rational-function field `Q(eps_eta_star, chi0_star, F_star)` — matrix products, a left-inverse identity, idempotent/complementary projectors, a decomposition, a static-strip linear equivalence, and quadratic-form norm identities. Mathematica verifies all of this natively (`Dot`, `IdentityMatrix`, `Together`/`Simplify`, `Transpose`). Verification by Mathematica is clearly possible, so the test ("is it possible," not "is it necessary") is met.

**Why this matters:**
Single-engine coverage leaves no independent cross-check on the projector algebra and norm identities; an algebra slip in the SymPy choreography (e.g., a sign in `S_rm_dep` or the `-1/(1-eps)` compiler entry) would pass undetected. The second engine must derive the same results via a different decomposition.

**Required change:**
Codex writes a NEW independent-route Mathematica script (not a transliteration of the `.py`). See directive for the claim manifest and acceptance criteria.

**Verification:**
`mathematica/...stage236..._mathematica_audit.wl` exists, exits 0, and asserts the M1-M9 manifest claims via a route distinct from the SymPy decomposition.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed yet. The directive instructs an independent route and forbids line-by-line porting.

## Engine cross-check

Only one engine present; no cross-check possible. SymPy output is fresh (output mtime 1778525516 > script mtime 1778522211), exit code 0, status PASS, and its printed forms (`S_rm_dep`, `C_rm_dep`, `P_nt_dep`, `P_eta_dep`, `||y_eta||^2 = 2*q_eta**2`) match the script's expected expressions and the paper.

## Verdict justification

The SymPy script is faithful and adversarially solid: every assertion compares a constructed quantity against an independently written expected form, so none are tautological; the static-strip check (line 137) genuinely tests `Delta_mu - Delta_Keta = q_nt` rather than asserting a definition against itself; the norm checks correctly reduce because `y_eta = (0,-q_eta,-q_eta)` and `R1_static/c_eta = -E1`. Paper alignment is exact — all eight notes-§8 deliverables and the appendix equations map to checks. The single defect is the missing second engine, which the dual-engine rule requires here because Mathematica can independently verify the linear-algebra claims. Verdict: `findings` (one `missing_mathematica`); not stop-cold, since adding a `.wl` cannot invalidate downstream results.

## Self-test notes

Variable-independence: no `diff`/`D` derivatives appear in this stage (pure linear algebra), so the zero-derivative trap is N/A. Symmetry/parity: no unbounded integrals; N/A. Trivial-case pre-check: I mentally substituted `q_tr=0` and confirmed `y_rm=(0,-q_eta,q_nt-q_eta)`, the static-strip residue `(q_nt-q_eta)-(-q_eta)=q_nt`, and `y_eta·y_eta=2q_eta^2` all reduce as asserted. Path: directive targets `mathematica/` with the mandatory `_mathematica_audit.wl` suffix. Paper round-trip: the proposed `.wl` must use the paper's `S_rm_dep`, `C_rm_dep` `(3,2)=-1/(1-eps_eta_star)`, and `L_rm_dep` exactly, introducing no new constant.
