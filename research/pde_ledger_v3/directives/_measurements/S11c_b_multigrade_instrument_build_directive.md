# S11c-b multigrade instrument — builder report (command-backed)

## Build + run
- Deliverable: `scripts/S11c_b_background_multigrade.py` (Codex build, log `~/.s11_build/S11c_b_multigrade_build_codex.log`, 184,077 tok). Both build legs SOUND (record: `S11c_b_multigrade_instrument_build_review.md`); the one tautological guard (`RECONSTRUCTION`) was folded out.
- Command: `cd research/pde_ledger_v3/scripts && python3 S11c_b_background_multigrade.py` (exit 0, empty stderr).
- Header (computed): `EMITTED_CASES = Integer(20)`, `GRADE_WINDOW_N = Integer(4)`.
- Consumed the committed, reviewed layer + real transcripts (printed by the run):
  - `A_MODULE_PATH = scripts/S11c_b_adjudicated_comparison.py`
  - `PY_TRANSCRIPT_PATH = scripts/out/S11c_b_brane_operator_sympy_audit.out`
  - `WL_TRANSCRIPT_PATH = mathematica/out/S11c_b_brane_operator_mathematica_audit.out`
- Guards (computed, all zero across the 20 cases): `WINDOW_CLEAN`, `GRADE_DIFFERENCE`, `REMAINDER_DIFFERENCE`.
- Raw stdout (outside repo, tree hygiene): `~/.s11_build/S11c_b_multigrade_run.out`.

## The computed `(eta_bg, sigma_W)` grade fingerprint (`a`=eta_bg power, `b`=sigma_W power)
Extracted with `scripts/S11c_b_grade_fingerprint.py` (a reporting tool over the run's emitted `MULTIGRADE_A/B/RESIDUAL`; structural zero-test on the instrument's already-normalised coefficients; full output `~/.s11_build/S11c_b_grade_fingerprint.txt`). Per case, the nonzero grade cells of the SymPy operand (A), the Wolfram operand (B), and their residual. These are COMPUTED objects; the per-case §2a/§3a/§3d adjudication (which retention is spec-correct) is the step record's, not stated here.

- **ADMISSIBILITY_OPERATOR_OPERAND (BODY_FORCE, DOF=THETA)** ×4 (both branches, both densities):
  A(PY)=`{(0,1),(1,1)}`, B(WL)=`{}`, RES=`{(0,1),(1,1)}`. WL operand identically zero.
- **SLAB_OPERATOR_TERM_ORIGINS (KINETIC)** ×4:
  U_MOMENTUM_ROWS[0..2]: A=B (RES `{}` — engines agree). THICKNESS_ROW: A(PY)=`{(0,0),(1,0),(2,0)}`, B(WL)=`{(0,0)}`, RES=`{(1,0),(2,0)}`.
- **SLAB_OPERATOR_TERM_ORIGINS (ADVECTIVE)** ×4:
  RHO4_CONSTANT: A(PY)=`{(0,1)}`, B(WL)=`{(0,0),(0,1),(1,0)}`, RES=`{(0,0),(1,0)}` (engines agree at (0,1); WL carries (0,0),(1,0) that PY does not).
  RHOBR_CONSTANT: A(PY)=`{}`, B(WL)=`{(0,0)}`, RES=`{(0,0)}`.
- **COUPLING_KERNEL (SECTOR=TRANSVERSE_TO_THICKNESS, with and without OBJECT=ADJOINTNESS_OPERAND_FORWARD)** ×8:
  A(PY)=`{(0,1),(1,1),(2,1),(3,1),(4,1)}` + exact remainder leading `{(5,1)}`; B(WL)=`{(0,1),(0,2),(1,1),(1,2)}`;
  RES=`{(0,1),(0,2),(1,1),(1,2),(2,1),(3,1),(4,1)}` + remainder `{(5,1)}`. (RHOBR variants drop the (0,2)/(1,2) or a subset — see full output.)

## What the fingerprint establishes (computed, not adjudicated)
The 20 genuine differences are NOT one systematic thing. Directions differ per family: for ADMISSIBILITY and KINETIC-THICKNESS the PY operand carries grade cells the WL operand lacks (WL zero, or WL only `(0,0)`); for ADVECTIVE the WL operand carries cells `(0,0),(1,0)` the PY operand lacks while both agree at `(0,1)`; for COUPLING the two operands each carry cells the other lacks (WL at `b=2`: `(0,2),(1,2)`; PY at `a≥2, b=1`: `(2,1),(3,1),(4,1),(5,1)…`) and differ at `(0,1),(1,1)`. First-shape-order cells (`a≤1, b≤1`) are populated in every family's residual. The original "systematic higher-background-order WL truncation" hypothesis does not uniformly explain the set; per-case adjudication (WL under-retention vs an IBP/divergence-form representational split vs a shape-order/§3a question vs a PY gap) is the step-2 task.
