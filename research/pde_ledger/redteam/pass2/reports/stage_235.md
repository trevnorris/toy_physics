---
unit_id: 235
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
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
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows: line 82 status row; lines 940-982 narrative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Rigid-mouth packet projector theorem: the finite same-charge packet splits into independent nontracking and dressing axes, and orbit lock is codimension two." Concretely (card body, notes §1-§5, appendix eqs 942-981) the stage builds: the involutive compiler `q_rm = M_rm x_rm` with `M_rm = [[-1,-c_eta],[0,1]]`, `M_rm^2 = I_2`; the complementary direct-space projectors `P_nt = [[1,c_eta],[0,0]]`, `P_eta = [[0,-c_eta],[0,1]]` (idempotent, orthogonal, summing to I); the decomposition `x_rm = x_nt + x_eta`; the codimension-two orbit-lock equivalences (`q_rm=0 ⇔ R1=E1=0 ⇔ Xi1=E1=0 ⇔ Xi1=R1=0`) against the codimension-one static strip `Xi1=0`; the static-blind dressing line with exact norm `||x_eta||^2 = (1+c_eta^2) q_eta^2`; and the static-only / full-orbit correction compilers. Card `\stagefield{Verification}` states "Mathematica audit: none yet"; notes §8 lists only the SymPy file.

## What the script claims to verify

Both engines exercise the full theorem as pure 2x2 linear algebra over `c_eta>0, R1,E1,q_nt,q_eta ∈ ℝ`. SymPy (§1-§5) and Mathematica (M1-M6) each: build `M_rm`, check `M_rm x_rm` equals `(Xi1, E1)` against an independently written `Xi1 = -R1 - c_eta E1`, check `M_rm^2 = I` and `det = -1`; construct `P_nt/P_eta` via `M_rm·Q·M_rm` and check against the closed forms + idempotency/orthogonality/completeness; verify the `x_nt`/`x_eta` split; solve the four orbit-lock systems to the unique origin; verify the blind line satisfies `Xi1=0`, the exact norm, and the `q_eta = L/√(1+c_eta^2)` size-L parameterization; and verify both correction compilers and their additive relation.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `M_rm`, `q_rm=M_rm x_rm`, `M_rm^2=I` | py 31-36 / wl 78-80 | match |
| inverse map / `det=-1` | py 38-39,90 / wl 81-82 | match |
| `P_nt`,`P_eta` closed forms + projector algebra | py 54-66 / wl 92-102 | match |
| `x_rm = x_nt + x_eta` decomposition | py 68-74 / wl 109-115 | match |
| codim-two orbit-lock equivalences | py 92-98 / wl 122-132 | match |
| static-blind line + `||x_eta||^2=(1+c_eta^2)q_eta^2` + size-L | py 109-120 / wl 140-148 | match |
| correction compilers (static/orbit/additive) | py 132-149 / wl 155-166 | match |
| card: "Mathematica audit: none yet"; notes §8 SymPy-only | a verified `.wl` is present (wl + fresh output) | mismatch |

`paper_alignment: partial` — all math deliverables match; the card/notes verification-coverage statement is stale.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 35 | `q_rm - (Xi1,E1) == 0` | compiler | yes |
| A2 | sympy | 36 | `M_rm^2 - I == 0` | involution | yes |
| A3 | sympy | 60-66 | projector algebra == 0 | projectors | yes |
| A4 | sympy | 71-74 | x split == 0 | decomposition | yes |
| A5 | sympy | 92-98 | solve→{0,0} | codim-two | yes |
| A6 | sympy | 112-120 | blind line / norm / L | static-blind | yes |
| A7 | sympy | 140-149 | corrections == 0 | compilers | yes |
| M1 | math | 78-82 | matrix/Det zero | compiler/involution | yes |
| M2 | math | 96-102 | projector algebra zero | projectors | yes |
| M3 | math | 112-115 | x split zero | decomposition | yes |
| M4 | math | 127-132 | Solve→origin / strip dim | codim-two | yes |
| M5 | math | 145-148 | blind line / norm / L | static-blind | yes |
| M6 | math | 159-166 | corrections zero | compilers | yes |

All rows non-tautological: the projector closed forms (`P_nt_expected`, `{{1,cEta},{0,0}}`) are checked against the *constructed* `M_rm·Q·M_rm`, which could disagree; `q_rm` is checked against an independently authored `Xi1`. No assertion is guaranteed by construction.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (verification-coverage statement is stale)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_235.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md:395-408`

**What's wrong:**
The card states (line 11): "SymPy audit: ... Mathematica audit: none yet." The notes §8 (lines 395-408) is titled "## 8. SymPy-backed status" and lists only the SymPy file as the supporting file. But a Mathematica audit (`...stage235..._mathematica_audit.wl`) now exists, runs clean, and has a fresh saved output ("All Stage 235 Mathematica checks passed."). The paper-side coverage statement does not reflect the present second engine.

**Why this matters:**
The card under-reports verification coverage; a reader citing MTDC-T11.6 would believe the stage is single-engine when it is dual-engine. This is a prose/documentation lag, not a math error. Direction (update card+notes to reflect the `.wl`, vs. some other intent) is a user call — Codex must not silently edit the card/notes.

**Required change:**
(paper_misalignment — routed to user; see `## Resolve before fix_loop` in directive. No Codex script edit.)

**Verification:**
After user resolution, card line 11 reads "Mathematica audit: \StageFile{mathematica/...stage235..._mathematica_audit.wl}." and notes §8 references the `.wl`.

## Independent-derivation check (Mathematica)

This is a pure 2x2 linear-algebra stage; the matrix `M_rm`, the axis projectors `Q_nt/Q_eta`, and the similarity `M_rm·Q·M_rm` are the single canonical route — there is no alternative "derivation" of a fixed 2x2 involution. The `.wl` is therefore necessarily structurally parallel to the `.py`, but it is NOT a mechanical transliteration: it uses native Mathematica primitives (`MatrixPower`, `Det`, `Solve[..., Reals]` returning rule-lists, `expectTrue` on `Solve` cardinality for the strip dimension at wl:132) rather than echoing SymPy's `sp.solve(dict=True)` choreography, and adds an independent check absent from the `.py` (M4 `lockFromQ = Solve[Thread[compiledDirect=={0,0}]]` and `originFromDirect = M_rm.{0,0}` directly, wl:122-129). For a fixed-matrix algebraic identity this constitutes legitimate independent cross-engine confirmation, not `mathematica_transliteration`. Quoted correspondences: py:54 `P_nt = M_rm*Q_nt*M_rm` ↔ wl:92 `projNt = mouthCompiler.axisNt.mouthCompiler` (same construction — unavoidable); py:92 `sp.solve([Eq(Xi1,0),Eq(E1,0)], [R1,E1], dict=True)` ↔ wl:130 `solveRules[{xiPacket==0,E1==0}]` then `=== {{R1->0,E1->0}}` (different engine idiom). No finding.

## Engine cross-check

Outputs agree at full strength. SymPy: `q_rm = (-E1·c_eta - R1, E1)`, `P_nt=[[1,c_eta],[0,0]]`, `P_eta=[[0,-c_eta],[0,1]]`, `||x_eta||^2 = q_eta^2(c_eta^2+1)`, all solves → `{R1:0,E1:0}`, `Delta_x_orbit = (c_eta·q_eta + q_nt, -q_eta)`. Mathematica: `M_rm.{R1,E1} = {-(cEta*E1)-R1, E1}`, identical `P_nt/P_eta`, `x_blind.x_blind = (1+cEta^2)*qEta^2`, all `Solve` → `{{R1->0,E1->0}}`, `Delta_x_orbit = {cEta*qEta + qNt, -qEta}`. Every residual zero / every PASS. No `engine_disagreement`.

## Verdict justification

The math holds up under attack in both engines and exactly matches the appendix eqs (942-981) and notes boxes. Attacks tried and failed: (1) tautology — the projector closed forms are checked against an independently constructed similarity product, and `q_rm` against a separately written `Xi1`, so the asserts can fail if construction is wrong; (2) variable-independence/self-test trap — no `diff`/`D[]` anywhere, pure linear algebra, no vacuous-derivative path; (3) symbol-assumption — `c_eta>0` and reals are exactly the card/notes domain (`c_eta = ε*/(1-ε*) > 0`), and positivity is only used for the `√(1+c_eta^2)` norm step, which is valid; (4) hardcoded result — the closed forms are *targets verified against constructions*, not injected answers; (5) missing branch — the orbit-lock is solved as a system (unique origin), and the strip is shown codim-one via `Solve` cardinality. Outputs are fresh (both `.txt` newer than their scripts). The sole defect is documentation: the card line 11 and notes §8 still describe the stage as SymPy-only while a verified second engine exists → one `paper_misalignment` (low), routed to the user. No stale self-labels ("Stage 251/252/253" drift flagged in pass-1 history was NOT found in any file read). Verdict: findings (1, paper-side only; no script defect).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 9 deliverable values checked, 0 value-mismatched (1 coverage-statement misalignment folded into F1, not a value mismatch).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `M_rm = [[-1,-c_eta],[0,1]]` | py:31 / wl:68 / out py:5-7, wl:19 | appendix:953; notes:99-104 | MATCH |
| `M_rm^2 = I_2`, `det=-1` | py:36,90 / wl:80,82 / out py:12,33 | appendix:955; notes:108 | MATCH |
| `P_nt = [[1,c_eta],[0,0]]` | py:57 / wl:96 / out py:16-18 | appendix:960; notes:152-157 | MATCH |
| `P_eta = [[0,-c_eta],[0,1]]` | py:58 / wl:97 / out py:20-22 | appendix:962; notes:158-163 | MATCH |
| `x_nt = (R1+c_eta E1, 0)` | py:71 / wl:112 / out py:24-26 | notes:188-199 | MATCH |
| `x_eta = (-c_eta E1, E1)` | py:73 / wl:114 / out py:28-30 | notes:202-213; appendix:980 | MATCH |
| codim-two: `q_rm=0 ⇔ R1=E1=0` (+ Xi1 variants) | py:92-98 / wl:122-132 / out py:34-36 | appendix:965-973; notes:230-247 | MATCH |
| `||x_eta||^2 = (1+c_eta^2) q_eta^2` | py:115 / wl:146 / out py:44, wl:88 | notes:275-277 | MATCH |
| `Delta_x_static=(q_nt,0)`, `Delta_x_orbit=(q_nt+c_eta q_eta,-q_eta)` | py:143,145 / wl:160-161 / out py:48-55 | notes:307-335 | MATCH |
| verification coverage: "Mathematica none yet" | a verified `.wl` exists (fresh output) | card:11; notes:395-408 | MISMATCH → F1 |

INTERNAL scaffolding (no finding): `assert_zero`/`assert_matrix_zero`/`expectZero`/`expectMatrixZero`/`expectTrue`/`pass`/`fail` helpers, `cleanExpr`/`stripConditional`, `zeroTwo`, intermediate `q_vec`/`x_q`/`x_from_q`, `sol_xi_E`/`sol_xi_R`/`lockFromQ`/`originFromDirect`/`solXiE`/`solXiR` (assertion drivers), `qeta_choice`/`lengthRule`/`L` (parameterization helper).

## Self-test notes

Checked the three traps. (1) Variable-independence: no `sp.diff`/`D[]` in either script — this is fixed-matrix linear algebra, so no absent-variable-derivative vacuous-pass path exists. (2) Symmetry/parity: no integrals; the only norm is an algebraic quadratic form, sign-definite by `c_eta>0`. (3) Trivial-case: substituting `R1=E1=q_nt=q_eta=0` collapses every `assert_zero`/`expectZero` residual to a literal 0 and the `Solve` systems to the unique origin, confirming the asserts track the actual claim rather than passing vacuously; the projector idempotency/orthogonality residuals are genuine 2x2 zero matrices only because `M_rm·Q·M_rm` truly is idempotent. The single finding is paper-side (no script edit), so no directive self-test of a proposed `diff` was needed.
