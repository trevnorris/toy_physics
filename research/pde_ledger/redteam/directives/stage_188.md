---
unit_id: 188
batch: V.3
created_at: 2026-06-01T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-01T17:19:25Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 188

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py:169-176`

**Issue:** Section VII is titled "Zero-set equivalence" and its comment claims to verify the paper's central zero-defect theorem (`Θ_1=Ξ_1=R_1=0 ⇔ δln R_tr=δln N_*=δln ε_η=0`). But its two assertions are `expect_zero(C_def_to_obs * Matrix([0,0,0]))` and `expect_zero(C_obs_to_quot * <that zero vector>)`. A matrix times the zero vector is the zero vector for ANY matrix, so these checks cannot fail and verify nothing about the equivalence. The real content of the theorem — that the maps are bijections, so the zero set is shared/unique — is already established elsewhere (invertibility checks at lines 84-85, 156-158; nonzero determinants printed at lines 79-80 and 133; round-trip at line 161). The fix replaces the vacuous trivial-direction substitution with assertions that actually exercise the bijection on the generic observable packet.

**Required change:**
Replace the block at lines 169-176, currently:
```python
subbanner("VII. Zero-set equivalence")
# Because both compilers are invertible, zero defect <-> zero observables <-> zero quotient packet.
# Verify by substituting zero defects and recovering zero observables exactly.
zero_def = sp.Matrix([0, 0, 0])
zero_obs_from_zero_def = sp.simplify(C_def_to_obs * zero_def)
zero_quot_from_zero_obs = sp.simplify(C_obs_to_quot * zero_obs_from_zero_def)
expect_zero("zero observables from zero defect", zero_obs_from_zero_def)
expect_zero("zero quotient packet from zero observables", zero_quot_from_zero_obs)
```
with:
```python
subbanner("VII. Zero-set equivalence (shared zero set via invertibility)")
# The nontrivial content of the iff is that both compilers are bijections:
# Delta_def == 0 forces Delta_obs == 0 (and conversely), and likewise for the
# quotient packet. Exercise it on the GENERIC packet, not on the zero vector
# (M * 0 == 0 holds for any M and proves nothing about uniqueness).
expect_zero(
    "C_obs->def then inverse recovers generic obs (bijection, def side)",
    sp.simplify(C_def_to_obs * (C_obs_to_def * obs)) - obs,
)
expect_zero(
    "C_obs->quot then inverse recovers generic obs (bijection, quot side)",
    sp.simplify(C_quot_to_obs * (C_obs_to_quot * obs)) - obs,
)
# Determinants are nonzero on the physical domain (also printed in II and IV):
#   det C_obs->quot = -1/C_tr,*  ,  det C_obs->def = -eps/(1-eps),  0<eps<1.
det_obs_to_def = sp.simplify(C_obs_to_def.det())
det_obs_to_quot = sp.simplify(C_obs_to_quot.det())
expect_zero(
    "1/det(C_obs->def) well-defined (nonzero det)",
    sp.simplify(det_obs_to_def * (1 / det_obs_to_def) - 1),
)
expect_zero(
    "1/det(C_obs->quot) well-defined (nonzero det)",
    sp.simplify(det_obs_to_quot * (1 / det_obs_to_quot) - 1),
)
```
Notes for Codex:
- `obs`, `C_obs_to_def`, `C_def_to_obs`, `C_obs_to_quot`, `C_quot_to_obs` are all already defined earlier in the script (lines 68, 122, 151, 71, 81). Do not redefine them.
- The two round-trip assertions are the load-bearing additions and MUST reference the generic `obs` packet (not `Matrix([0,0,0])`). They reduce to `obs - obs = 0` only because the matrices are genuine inverses; a singular compiler would fail them.
- You may keep the original two `M*0` prints as extra informational output, but they must NOT be the only checks in Section VII.

**Verification command:**
The verifier will run `redteam exec-sympy 188` and confirm (a) Section VII now contains at least two `expect_zero` calls whose first matrix-product operand references `obs` (the generic packet), and (b) the script exits 0 with all checks zero.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`
- summary: Replaced the vacuous zero-vector checks with generic observable round-trip bijection checks and determinant nonzero checks.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py:46,59`

**Issue:** Line 46 sets `Cstar = sp.simplify(1 / Ctr)`. Line 59 then asserts `expect_zero("C_* C_tr,* - 1", Cstar * Ctr - 1)`. Because `Cstar` is defined as `1/Ctr`, the expression `Cstar*Ctr - 1` is identically zero by construction and the assertion cannot fail. It does not independently exercise the paper's `C_*` closed form (appendix `eq:app-part05-Tstar-def`: `C_* = (1+χ₀,*)(1+δU,*)(1+χ₀,*+δU,*)/(χ₀,* δU,*)`). Define `C_*` from its own closed form and check it against the reciprocal of the independently built `C_tr,*`, turning the tautology into a genuine cross-check.

**Required change:**
1. Replace line 46:
```python
Cstar = sp.simplify(1 / Ctr)
```
with the independent closed form (appendix eq:app-part05-Tstar-def):
```python
Cstar = sp.simplify(
    (1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)
    / (chi0s * deltaUs)
)
```
2. Replace the line-59 assertion:
```python
expect_zero("C_* C_tr,* - 1", Cstar * Ctr - 1)
```
with a cross-check between the independently defined `C_*` and the reciprocal of the independently defined `C_tr,*`:
```python
expect_zero("C_* - 1/C_tr,*", Cstar - 1 / Ctr)
```
Notes for Codex:
- The numeric/symbolic value of `Cstar` is unchanged (the new closed form equals `1/Ctr`), so every downstream use — in particular `C_obs_to_quot = sp.diag(-Cstar, 1, 1)` at line 71 and `det(C_obs->quot)` at line 80 — is unaffected and must NOT be edited.
- Leave the adjacent line-60 check `expect_zero("A_tr,* - B_* C_tr,*", Astar - Bstar * Ctr)` exactly as is; it is the genuine load-bearing coefficient identity.

**Verification command:**
The verifier will run `redteam exec-sympy 188` and confirm (a) line 46 defines `Cstar` from the four positive symbols `chi0s, deltaUs` (not as `1 / Ctr`), (b) the Section I assertion now reads `Cstar - 1 / Ctr` (or equivalent independent comparison), and (c) the script exits 0 with `C_* - 1/C_tr,* = 0` printed.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`
- summary: Defined `Cstar` from the appendix closed form and changed the coefficient check to compare it independently against `1/Ctr`.
- deviation: none

## F3 — missing_mathematica

**Issue:** Stage 188 is dual-engine-capable (exact rational-function matrix algebra: coefficient identities, compiler inverses, determinants) but has no Mathematica `.wl`. Under the dual-engine rule, an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage188_branch_observables_completion_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the CORRECTED SymPy script (after F1/F2 above) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_188.tex` and the notes file are the source of truth. Mirror the CORRECTED checks (the de-tautologized bijection round-trips and the independent `C_*` closed form), NOT the broken originals.
- Use Mathematica-NATIVE primitives (`Inverse`/`Det`/`LinearSolve`, matrix products, `Solve`/`Reduce`, `Series`+`Coefficient`) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port with the same variable names and step order. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (`expectZero`/`expectZeroMatrix`, `$Assumptions` positivity, `stripCE`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion must match the SymPy result.

**Anti-transliteration:** a `.wl` that merely re-types the SymPy closed forms and subtracts them is a transliteration and will be REJECTED at verification. Design a genuinely independent route. RUN it (`timeout 600 math -script <path>`) and iterate to exit 0; a timeout (124) is a failure — reformulate, don't raise the cap.

**Verification command:** the verifier runs `redteam exec-mathematica 188`, confirms exit 0 with all PASS lines, and reviews that the `.wl` is a genuinely independent route whose conclusions agree with the SymPy engine.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage188_branch_observables_completion_mathematica_audit.wl`
- summary: Added a Mathematica audit that derives the observable, quotient, and defect compilers from finite log-drift equations, solves the linear maps, and checks the corrected SymPy conclusions.
- deviation: none
