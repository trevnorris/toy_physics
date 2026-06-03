---
unit_id: 235
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
  notes_stage_files: [moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 235 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_235.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows/section at lines 82, 883-1005 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 235 makes the rigid-mouth packet structure exact on the direct observable plane `x_rm = (R1, E1)^T = (δln R_target, δln ε_η)^T`. `\stagefield{Output}` states: "Rigid-mouth packet projector theorem: the finite same-charge packet splits into independent nontracking and dressing axes, and orbit lock is codimension two." The distinct deliverables (card "Derivation ledger" + appendix eqs `app-part07-mrm`, `-direct-projectors`, `-rm-codim-two`, `-static-blind-line`, and notes §1–§5) are: (1) the involutive compiler `q_rm = M_rm x_rm` with `M_rm = [[-1,-c_η],[0,1]]` and `M_rm^2 = I_2`; (2) the exact direct-space projectors `P_nt = [[1,c_η],[0,0]]`, `P_η = [[0,-c_η],[0,1]]`, which are complementary idempotents (`P^2=P`, `P_nt P_η = 0`, `P_nt + P_η = I`) realizing the "independent axes" split; (3) the codimension-two orbit-lock equivalence `q_nt=0, q_η=0 ⟺ R1=0, E1=0` (also `⟺ Ξ1=0,E1=0 ⟺ Ξ1=0,R1=0`); (4) the static-blind dressing line `Ξ1=0 ⟺ x_rm = (-c_η q_η, q_η)^T`; (5) its exact norm `‖x_η‖² = (1+c_η²) q_η²` with the explicit `q_η = L/√(1+c_η²)` choice making `‖x_rm‖ = L`; and (6) the static-only / full-orbit correction compilers `Δx_static`, `Δx_orbit`. The microscopic dependent-triple carrier `(0,-q_η,q_nt-q_η)` is explicitly assigned to Stage 236 (appendix line 984), so it is out of scope here.

## What the script claims to verify

The SymPy script (docstring/print line 20) verifies the rigid-mouth packet projectors, static-blind dressing line, and codimension-two orbit-lock point. Its assertions: builds `M_rm` independently and asserts `q_rm = (Ξ1, E1)` (line 35) and `M_rm^2 = I` (line 36) and the inverse map (line 39); constructs `P_nt = M_rm Q_nt M_rm`, `P_η = M_rm Q_η M_rm` and asserts they equal the paper's closed forms plus full idempotency/complementarity (lines 60–66); asserts the direct-space decomposition `x_rm = x_nt + x_η` with explicit component forms (lines 71–74); checks `det(M_rm)=-1` and solves the two codim-two equivalence systems (lines 90–98); verifies the static-blind line vanishes `Ξ1` and has norm `(1+c_η²)q_η²` with the `L`-scaling (lines 112–120); and verifies both correction compilers and their additive relation (lines 137–149). This is exactly the set of paper deliverables (1)–(6).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `q_rm = M_rm x_rm`, `M_rm^2=I`, inverse map | lines 35, 36, 39 | match |
| (2) `P_nt`, `P_η` closed forms + idempotent/complementary | lines 60–66 | match |
| (3) codim-two `q_nt=q_η=0 ⟺ R1=E1=0` (+ Ξ1 forms) | lines 90, 92–98 | match |
| (4) static-blind line `Ξ1=0 ⟺ (-c_η q_η, q_η)` | lines 110–112 | match |
| (5) `‖x_η‖²=(1+c_η²)q_η²`, `q_η=L/√(1+c_η²)⇒‖x‖=L` | lines 114–120 | match |
| (6) `Δx_static`, `Δx_orbit`, additive split | lines 137–149 | match |
| Stage-236 microscopic triple (out of scope here) | n/a | correctly absent |
| Independent second-engine (Mathematica) verification | none | missing |

`paper_alignment: aligned` — every Stage-235 paper deliverable has a faithful, non-tautological SymPy check, and the script tests nothing beyond the paper's scope.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 35 | `assert (M_rm x_rm) - (Ξ1, E1) == 0` | (1) compiler | yes |
| A2 | sympy | 36 | `assert M_rm^2 - I == 0` | (1) involution | yes |
| A3 | sympy | 39 | `assert M_rm (q_nt,q_η) - (-q_nt-c q_η, q_η) == 0` | (1) inverse map | yes |
| A4 | sympy | 60–61 | `assert P_nt - [[1,c],[0,0]] == 0`, `P_η - [[0,-c],[0,1]] == 0` | (2) projector closed form | yes |
| A5 | sympy | 62–66 | idempotency `P^2=P`, `P_nt P_η=0`, `P_η P_nt=0`, `P_nt+P_η=I` | (2) independent axes | yes |
| A6 | sympy | 71–74 | `x_nt=(R1+c E1,0)=(-Ξ1,0)`, `x_η=(-c E1,E1)`, `x_rm=x_nt+x_η` | (2) decomposition | yes |
| A7 | sympy | 90 | `assert det(M_rm)+1==0` (det=-1) | (3) invertibility for codim-two | yes |
| A8 | sympy | 92–98 | `solve(Ξ1=0,E1=0)` and `solve(Ξ1=0,R1=0)` give `(0,0)` | (3) codim-two equivalence | yes |
| A9 | sympy | 112 | `assert Ξ1_blind == 0` on the line | (4) static-blind line | yes |
| A10 | sympy | 115 | `assert ‖x_η‖² - (1+c²)q² == 0` | (5) exact norm | yes |
| A11 | sympy | 120 | `assert ‖x_η‖²(q=L/√(1+c²)) - L² == 0` | (5) arbitrarily-far scaling | yes |
| A12 | sympy | 140–149 | `Δx_static=(q_nt,0)`, `Δx_orbit=(q_nt+c q_η,-q_η)`, split, `x_q+Δx_orbit=0` | (6) correction compilers | yes |

All twelve assertion blocks are non-tautological: in every case the left side is built from an independent construction (matrix products `M_rm Q M_rm`, solves, substitutions) and compared against the paper's stated closed form or zero. None is `x = expr; assert x == expr`. The `P_nt_expected` / `P_eta_expected` literals (lines 57–58) are the paper's *claimed* targets, not the source of the construction, so A4 is a real falsifiable check rather than a hardcoded round-trip.

## Findings

### F1 — missing_verification_script (subtype: missing_mathematica)

**Severity:** medium
**Files:**
- `(missing)` — no `mathematica/moving_throat_pde_stage235_*_mathematica_audit.wl` exists
- paper card states the gap: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_235.tex:11` — "Mathematica audit: none yet."

**What's wrong:**
The unit is `is_status_only_candidate: False` and `is_checkpoint: False`, so the dual-engine policy applies: a second-engine `.wl` is required wherever Mathematica *can* independently verify the claims. Every Stage-235 deliverable is finite linear algebra over a 2×2 matrix with one symbolic parameter `c_η` (involution, projector idempotency/complementarity, a 2-variable linear `Solve`, a quadratic-form norm, and additive correction identities). Mathematica verifies all of this natively with a genuinely different decomposition than SymPy (e.g. `MatrixPower[M,2]`, `IdentityMatrix`, `Eigen`/`Solve`, `Simplify` of inner products), so independent verification is clearly *possible*. It is currently absent.

**Why this matters:**
With only one engine, a transcription or convention error in the SymPy algebra (sign of `c_η`, projector ordering `M_rm Q M_rm` vs `S Q M`) cannot be caught by cross-check. The dual-engine requirement exists precisely to catch single-engine algebra slips.

**Required change:**
Codex writes a NEW independent-route Mathematica script (see directive F1 claim manifest). It must NOT transliterate the `.py`; it should re-derive the projectors and orbit-lock structure via native Mathematica primitives.

**Verification:**
`mathematica/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_mathematica_audit.wl` exists, exercises manifest items M1–M6, and exits 0 with all in-file checks passing; a saved output `.txt` is produced.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is not yet a risk — but the directive flags it for Codex: the new script must use a different decomposition (matrix primitives + `Solve`), not echo the SymPy variable choreography.

## Engine cross-check

Only one engine present; not applicable. After F1 lands, both engines should agree on `M_rm^2=I`, the projector closed forms, the codim-two solution set, and `‖x_η‖²=(1+c_η²)q_η²`.

## Verdict justification

The SymPy script is faithfully aligned with the Stage-235 paper card, appendix (lines 940–982), and notes §1–§5: all six deliverables are covered by non-tautological, well-anchored assertions, and nothing out-of-scope (the Stage-236 microscopic triple is correctly excluded). I attacked the projector construction (`M_rm Q M_rm` correctly equals the notes' `S_rm Q M_rm` since `M_rm` is involutive), the `P_nt_expected` literals (independent LHS construction, not a self-check), the `c_η` positivity assumption (not load-bearing for any identity — all are polynomial identities valid for any `c_η`), and the codim-two `solve` checks (genuinely unique `(0,0)` solution) — none broke. The saved output (mtime 12:51) is newer than the script (11:56), so it is fresh and matches what the current script produces. The sole finding is the absent second engine, which Mathematica can clearly provide; hence `verdict: findings` with one `missing_mathematica` item, not `clean`.

## Self-test notes

Checked trap #1 (variable independence / `diff`): no derivatives in the proposed `.wl` — pure linear algebra, so the identically-zero-derivative trap is N/A. Trap #2 (parity/unbounded integrals): no integrals; N/A. Trap #3 (trivial-case): I hand-evaluated `M_rm^2=I`, `P_nt=[[1,c],[0,0]]`, `P_nt P_η=0`, `Solve[{-R1-c E1==0,E1==0}]={0,0}`, and `(1+c²)·L²/(1+c²)=L²` — all reduce correctly. Trap #4 (paths): directive names the full `mathematica/…_mathematica_audit.wl` target. Trap #5 (paper round-trip): the new script's targets are the paper's own closed forms (eqs app-part07-mrm/-direct-projectors/-static-blind-line), so it introduces no new constant or `paper_misalignment`. Separately noted (not a finding): the notes file carries stale "Stage 251/252/253" labels and a stale `stage252_…py` supporting-file reference (notes lines 5,8,15,25,88,259,354,356,362,375,379,391,408) — the known project-wide incomplete-renumber; canonical numbering (paper-card/script/appendix = 235) is ground truth, math content matches, so this is not `paper_misalignment` and not routed to Codex.
