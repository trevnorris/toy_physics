---
unit_id: 190
batch: V.3
created_at: 2026-06-01T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-01T17:33:02Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 190

The SymPy audit found this stage's checks substantive and correct (verdict was clean). The only required change is the dual-engine gap below. Do NOT modify the SymPy script (it is correct and is the reference engine). Do NOT touch paper.tex or notes/.

After creating the script, RUN it (`timeout 600 math -script <path>`) and iterate until it exits 0 with all checks passing. A timeout (exit 124) is a FAILURE — reformulate the math, never raise the cap. The orchestrator independently re-runs afterward.

## F1 — missing_mathematica

**Issue:** Stage 190 is dual-engine-capable (exact symbolic / rational-function algebra) but has no Mathematica `.wl`. Under the dual-engine rule, an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the SymPy script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_190.tex` and the notes file are the source of truth for the math.
- Use Mathematica-NATIVE primitives (`Series`+`Coefficient`, `D[...]`, `Integrate`, `Solve`/`Reduce`, `Inverse`/`Det`/`LinearSolve`, matrix products, etc.) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port with the same variable names and step order. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (the `expectZero` helpers, `$Assumptions` positivity declarations, `stripCE` for `ConditionalExpression`, the `math -script` entry/exit convention).
- Assert cross-engine agreement: each conclusion the `.wl` derives must match the SymPy result.

**Anti-transliteration:** a `.wl` that merely re-types the SymPy closed forms and subtracts them is a transliteration and will be REJECTED at verification. Design a genuinely independent route.

**Verification command:** the verifier runs `redteam exec-mathematica 190`, confirms exit 0 with all PASS lines, and reviews that the `.wl` is a genuinely independent route (native primitives, different decomposition) whose conclusions agree with the SymPy engine.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl`
- summary: Created the independent Mathematica audit covering support-blindness, microscopic slippage laws, direct-defect projection, triangular inverse reconstruction, and grouped P2 scalar no-go checks.
- deviation: none
