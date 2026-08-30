# S11c-b step 4 — re-run + re-adjudication against the repaired WL transcript (2026-08-30)

The reviewed row-residual instrument (`scripts/S11c_b_row_residual.py`, committed `ef26084c`) re-run against
the repaired WL transcript (`mathematica/out/S11c_b_brane_operator_mathematica_audit.out`, regenerated +
committed in `b875cdde` — the admissibility θ full-field lift) and the unchanged PY transcript. Every claim
carries its command (rule 2).

## Run
```
$ python3 scripts/S11c_b_row_residual.py > ~/.s11_build/S11c_b_row_residual_rerun.out 2> …rerun.err
  # exit 0; stdout 26,758,678 B (byte-count identical to the prior fixrun); err empty; 100 cases
  #   (52 admissibility, 20 coupling, + slab/origins), no error.
```
Reproducible from the committed transcripts + instrument (ran the real consumer, not a re-read).

## Result 1 — ADMISSIBILITY θ residual → 0 (the repair, CONFIRMED by the reviewed comparator, rule 13)
All four branch×density admissibility `DOF=THETA` body-force cases now emit `ROW_RESIDUAL = Integer(0)`:
```
$ awk '/CASE.*ADMISSIBILITY.*DOF.*THETA/{...}/^ROW_RESIDUAL/{...}' rerun.out
  LAB_HELD/RHO4_CONSTANT           THETA: RESIDUAL=0
  LAB_HELD/RHOBR_CONSTANT          THETA: RESIDUAL=0
  MATERIAL_ADVECTED/RHO4_CONSTANT  THETA: RESIDUAL=0
  MATERIAL_ADVECTED/RHOBR_CONSTANT THETA: RESIDUAL=0
```
`ROW_OPERAND_WL` now equals `ROW_OPERAND_PY_TRUNC` (both `= kappa_theta_W·σ_W·(η w1 − 1)·∇²w1/(L_W·W_0)`, the
comparator's symbol-reconciled form; WL's `gradientThetaEwCoefficient` standardises to `kappa_theta_W`).
Before the repair WL's θ operand was `Integer(0)` and the residual was PY≠0 (see §1/§8 of
`S11c_b_step2_adjudication.md`). ⇒ the admissibility θ under-retention is REPAIRED; WL agrees with PY.

## Result 2 — nothing else moved (consistency vs the prior fixrun)
```
$ diff <(grep -a '^ROW_RESIDUAL = ' fixrun.out) <(grep -a '^ROW_RESIDUAL = ' rerun.out) | grep -c '^[<>]'
  8   # = 4 changed residual lines × (one '<' old + one '>' new) = exactly 4 residuals changed
```
The only four ROW_RESIDUAL lines that differ from the pre-repair fixrun are the four admissibility θ cases
(PY≠0 → 0). Every other residual — all other admissibility components (E_W/U/PER_FACE_TRACTION, still RES=0),
the entire SLAB_OPERATOR / TERM_ORIGINS, and all 20 COUPLING_KERNEL cases — is byte-identical to the prior
run. Consistent with the transcript-level scope check (only the 6 admissibility-family tags changed; §build
review) and the build legs' engine byte-identity. No unintended propagation.

## Per-family standing after the repair (rule 13)
1. **ADMISSIBILITY θ — REPAIRED.** Residual → 0 (all 4). WL now carries the §3d full-field `∇²W_bg` body force
   and agrees with PY. Mechanism (§9 + build review): the full-field lift of the MIXED `∇θ·∇e_W` invariant
   (`gradFullEw = anchoredWidth^(-1) gradient[fullWidth]`, i.e. `∇(W_bg+δW)/W_bg` — NOT a `ρ_4D` density
   weighting, NOT `/W_0`). Other admissibility components already agreed and still do.
2. **KINETIC — WL correct; NOT a gap (reversed in §9).** Residual `μ·e_W_tt·(W_bg²−W_0²)` persists (WL emits
   `μ_W W_0²`, unchanged; §1a: `e_W ≡ δW/W_0` with constant `W_0` ⇒ `δW = W_0 e_W` ⇒ the e_W-row inertia IS
   `μ_W W_0²`, ∂/∂W_bg = 0) — a representational `e_W` vs `e_W,bg=(W_0/W_bg)e_W` normalization artifact, not a
   WL under-retention. No repair. ⚠ Owed at re-adjudication of the family: confirm whether the
   sibling's `W_bg²` is a bare normalization convention (representational) or a genuine sibling bug.
3. **ADVECTIVE — representational** (continuity constraint-vs-evolution), unchanged. ⚠ owed: confirm PY
   imposes continuity as the constraint.
4. **COUPLING — genuine in-scope disagreement**, unchanged (`IN_SCOPE_WEAK_REMAINDER` nonzero ×20). §3c
   which-engine spec-adjudication owed (task #84); the ADJOINTNESS_RESIDUAL its own thread.
5. **ENERGY-BASIS (U_MOMENTUM/THICKNESS strong rows)** — deferred §1d quotient, unchanged (task #85).

⇒ STEP 3+4 CLOSE: one genuine WL under-retention (admissibility θ) repaired and re-adjudicated to agreement;
the kinetic "gap" was reversed (WL correct); advective representational; coupling a genuine surviving
disagreement (the headline — not resolved); energy-basis deferred.
