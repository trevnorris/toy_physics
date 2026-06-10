---
unit_id: 236
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 236 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_236.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row 84 + narrative lines 984-1005)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` states: "Microscopic projector theorem: clearing the static nontracking defect moves along a different dependent direction than clearing the post-static dressing defect; a static-only fix leaves a dressing gap." The derivation ledger names the deliverables: the dependent-plane compiler, projectors `P_nt^dep` and `P_eta^dep`, the equal-drift dressing ray, and the static-only restoration direction. The notes make these exact: rigid-mouth packet map `M_rm = [[-1,-c_eta],[0,1]]` with `c_eta = eps_eta_star/(1-eps_eta_star)`; dependent section `S_rm_dep = [[0,0],[0,-1],[1,-1]]` giving `y_rm = (0, -q_eta, q_nt-q_eta)`; compiler `C_rm_dep = S_rm_dep M_rm = [[0,0],[0,-1],[-1,-1/(1-eps_eta_star)]]`; left inverse `L_rm_dep = [[0,-1,1],[0,-1,0]]` with `L S = I_2`; complementary projectors `P_nt^dep = [[0,0,0],[0,0,0],[0,-1,1]]`, `P_eta^dep = [[0,0,0],[0,1,0],[0,1,0]]` (idempotent, orthogonal, summing to diag(0,1,1)); static-strip equivalence `q_nt=0 iff Delta_mu=Delta_Keta`; equal-drift ray `y_eta = -q_eta(0,1,1)^T`; norm `||y_eta||^2 = 2 q_eta^2 = 2 E1^2 = 2 R1^2/c_eta^2`; static-only correction `-y_nt = (0,0,-q_nt)` and full orbit-lock correction `-y_rm`. The appendix (eq. app-part07-dependent-correction-q, app-part07-equal-drift-ray, Thm rigid-mouth-packet) reproduces the dependent correction `(0,-q_eta,q_nt-q_eta)` and the ray `-q_eta(0,1,1)`.

## What the script claims to verify

Both scripts verify, by exact symbolic algebra over real symbols, the full chain of deliverables: the packet map action, the rigid-mouth dependent section on `q_tr=0`, the direct-observable-to-dependent compiler `C_rm_dep`, the left inverse and packet-coordinate recovery, the two complementary projectors with idempotency/orthogonality/completeness, the `y_rm = y_nt + y_eta` decomposition, the static-strip iff-equivalence, the equal-drift ray and its three equivalent norms, and the static-only vs. full orbit-lock correction compilers. There are no numeric constants, no transcendental functions, no `diff`/`D` — every check is a residual-zero of a matrix/vector identity in symbolic variables.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `M_rm` packet map → `q_nt=-R1-c_eta E1, q_eta=E1` | py L33-36 / wl M1 | match |
| `S_rm_dep`, `y_rm=(0,-q_eta,q_nt-q_eta)` | py L42-46 / wl M2 | match |
| Compiler `C_rm_dep` and its action | py L57-78 / wl M3 | match |
| Left inverse `L_rm_dep`, `L S = I_2`, recovery | py L89-93 / wl M4 | match |
| Projectors `P_nt^dep`, `P_eta^dep`, idem/orth/complete | py L98-112 / wl M5-M6 | match |
| Decomposition `y_rm = y_nt + y_eta` | py L114-119 / wl M7 | match |
| Static strip `q_nt=0 iff Delta_mu=Delta_Keta`, ray | py L137-149 / wl M8 | match |
| Norm `||y_eta||^2 = 2q_eta^2 = 2E1^2 = 2R1^2/c_eta^2` | py L144-149 / wl M8 | match |
| Static-only / orbit-lock corrections | py L160-168 / wl M9 | match |

`paper_alignment: aligned`. Every boxed deliverable in the notes/appendix has a faithful, non-tautological script check in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 36 | `assert_matrix_zero(q_rm - [...])` | M_rm map | yes |
| A2 | sympy | 46 | `y_rm - S_rm_dep·[q_nt,q_eta]` | dependent section | yes |
| A3 | sympy | 66,69 | `C_rm_dep - C_expected`, action | compiler | yes |
| A4 | sympy | 90,93 | `L·S - I_2`, `L·y_rm - [q_nt,q_eta]` | left inverse/recovery | yes |
| A5 | sympy | 104-112 | proj entries, idem, orth, sum | projectors | yes |
| A6 | sympy | 117-119 | `y_nt`,`y_eta`,decomp | decomposition | yes |
| A7 | sympy | 137 | `(y_rm[2]-y_rm[1]) - q_nt` | static-strip equiv | yes |
| A8 | sympy | 141-149 | strip vector + 3 norms | equal-drift ray/norm | yes |
| A9 | sympy | 163-168 | static/orbit corrections | correction compilers | yes |
| B1 | math | M1-M2 | `expectZero` map/section | maps | yes |
| B2 | math | M3 | basis-action compiler + product-agree | compiler (independent route) | yes |
| B3 | math | M4 | `LinearSolve` AND `Solve` recovery agree, `L·S-I` | left inverse (two independent routes) | yes |
| B4 | math | M5 | proj from recovered action, idem, orth, **eigenvalues**, **MatrixRank** | projectors (extra structural) | yes |
| B5 | math | M8 | `Equivalent[Delta_mu==Delta_Keta, q_nt==0]` iff, ray, norms | static strip/ray | yes |
| B6 | math | M9 | corrections | correction compilers | yes |

All rows Anchored=yes. No tautologies: each `assert_matrix_zero`/`expectZero` is `(computed expression) − (independently stated target)`, where the computed side flows from matrix products of inputs, not from re-substituting the target.

## Independent-derivation check (Mathematica)

The `.wl` is a **genuinely independent route**, not a transliteration. Three corresponding sections show different choreography:

1. Compiler. `.py` (L57) forms `C_rm_dep = simplify(S_rm_dep * M_rm)` directly. `.wl` (M3, L115-119) instead computes the columns by *acting on basis vectors* `sRmDep . (mRm . #) & /@ IdentityMatrix[2]`, transposes, and *separately* cross-checks `sRmDep.(mRm.obs) - cRmDep.obs == 0` — a route the `.py` never uses.

2. Left inverse. `.py` (L89) *hardcodes* `L_rm_dep` and asserts `L·S = I`. `.wl` (M4, L141-152) *derives* the recovery map two independent ways — `LinearSolve[planeBlock,{dK,dMu}]` and `Solve[sRmDep.{uNt,uEta}=={0,dK,dMu}]` — confirms they agree, then reconstructs the matrix via `Coefficient[...]` against `deltaCoords`. This is the inverse problem, not the forward assertion in `.py`.

3. Projectors. `.py` (L98) forms `P = S·Q·L` and compares to literals. `.wl` (M5) reconstructs the projector matrices from their *action on a generic vector* via `Coefficient`, then adds checks absent from `.py`: `Eigenvalues` = `{0,0,1}` and `MatrixRank` = 1. M8 also uses `Equivalent[...]` for the iff, which `.py` only checks one direction of (residual `= q_nt`).

The all-zero residual lines in the `.wl` output (e.g. M5 entries `{{0,0,0},...}`, eigenvalues `{0,0,0}`, rank `0`) are computed−expected differences, not the objects themselves; the printed matrices (`P_nt_dep = {{0,0,0},{0,0,0},{0,-1,1}}`, etc.) confirm correct nonzero structure. No `mathematica_transliteration` finding.

## Engine cross-check

Final structures agree exactly across engines: `S_rm_dep`, `C_rm_dep = {{0,0},{0,-1},{-1,(-1+epsEtaStar)^-1}}` (`.py` prints the algebraically identical `1/(eps_eta_star-1)` form), `L_rm_dep = {{0,-1,1},{0,-1,0}}`, `P_nt_dep`, `P_eta_dep`, `y_nt = (0,0,q_nt)`, `y_eta = (0,-q_eta,-q_eta)`, `Delta_y_static = (0,0,-q_nt)`, `Delta_y_orbit = (0,q_eta,q_eta-q_nt)`, `||y_eta||^2 = 2 q_eta^2`, equal-drift strip `(0,-E1,-E1)`. No `engine_disagreement`.

## Verdict justification

Clean. Attacks tried and failed: (1) variable-independence trap — no `diff`/`D` anywhere; checks are pure linear-algebra residuals, so no vacuous-derivative pass. (2) Tautology — each assertion compares a product-derived quantity to an independently written target (e.g. `C_rm_dep` from `S·M` vs. the literal `C_expected`); the `.wl` further derives `L` via Solve/LinearSolve rather than asserting it, so it cannot be circular. (3) Symbol-assumption error — `eps_eta_star` is `positive` (py) / `0 < epsEtaStar < 1` (wl), matching `c_eta = eps_eta_star/(1-eps_eta_star)` which requires `eps_eta_star != 1`; the assumption is justified and does not over-strengthen any simplify (all identities are polynomial/rational in the symbols and hold identically). (4) Paper alignment — every boxed notes/appendix result maps to a matching check in both engines; the equal-drift ray `-q_eta(0,1,1)`, dependent correction `(0,-q_eta,q_nt-q_eta)`, and norm `2q_eta^2` all match verbatim. (5) Freshness — both `.txt` (2026-06-02) are newer than both sources (py 2026-05-11, wl 2026-06-02). No stale "Stage 253" or other drifted self-label appears in either current source or output (the pass-1 heads-up residual is not present in the live files). Both engines present, independent, agree.

## Value Reconciliation (pass-2 augmentation)

All emitted deliverables are exact symbolic matrices/vectors; there are no free numeric constants. Reconciled against notes `.md` (boxed) and appendix `.tex`.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_rm = [[-1,-c_eta],[0,1]]`, `c_eta=eps_eta_star/(1-eps_eta_star)` | py L28,33 / wl L68,82 | notes L80-90 | MATCH |
| `S_rm_dep = [[0,0],[0,-1],[1,-1]]`, `y_rm=(0,-q_eta,q_nt-q_eta)` | py out L4-15 / wl out L20 | notes L126-141; appx eq L987-989 | MATCH |
| `C_rm_dep = [[0,0],[0,-1],[-1,-1/(1-eps_eta_star)]]` | py out L18-25 / wl out L31 | notes L166-173 | MATCH |
| `L_rm_dep = [[0,-1,1],[0,-1,0]]` | py out L36-39 / wl out L46 | notes L210-217 | MATCH |
| `P_nt_dep=[[0,0,0],[0,0,0],[0,-1,1]]`, `P_eta_dep=[[0,0,0],[0,1,0],[0,1,0]]` | py out L40-51 / wl out L71-72 | notes L242-256 | MATCH |
| `y_nt=(0,0,q_nt)`, `y_eta=(0,-q_eta,-q_eta)` | py out L52-63 / wl out L91-92 | notes L290-305 | MATCH |
| Static strip: `Delta_mu-Delta_Keta=q_nt`; equal-drift ray `-q_eta(0,1,1)` | py out L66-72 / wl out L97-113 | notes L322-360; appx eq L991-996 | MATCH |
| `||y_eta||^2 = 2 q_eta^2` (= `2 E1^2` = `2 R1^2/c_eta^2`) | py out L73 / wl out L114 | notes L364-366 | MATCH |
| `Delta_y_static=(0,0,-q_nt)`, `Delta_y_orbit=(0,q_eta,q_eta-q_nt)` | py out L76-87 / wl out L127-128 | notes L391-431 | MATCH |

INTERNAL scaffolding (no finding): pass/fail flags, all `expectZero`/`assert_matrix_zero` residuals (printed as `{0,...}`), eigenvalue/rank check residuals, identity matrices `I_2`/`diag(0,1,1)`, `Q_nt`/`Q_eta` selectors, `recoveredByLinearSolve`/`recoveredBySolve` intermediates.

reconciliation: complete; 9 deliverable values checked, 0 misaligned

## Self-test notes

Variable-independence trap: N/A — no `sp.diff`/`D[]` in either script; all checks are linear-algebra residuals, so no zero-derivative vacuous pass is possible. Symmetry/parity: N/A — no integrals. Trivial-case pre-check: substituting concrete `q_nt=q_eta=0` collapses every printed deliverable to the zero vector and every residual stays 0 (consistent), while generic symbolic inputs keep the targets nonzero, so the asserts are non-vacuous. Path/round-trip: no findings, no directive written.
