# Independent build review — S11c-b requested-truncation complete-row residual instrument

You are reviewing a **Codex-built script**. Your job is to find every way its computed objects could be
wrong, tautological, or leak/hardcode an answer. A prose re-derivation is worth nothing here — **run
ablations and report their literal stdout**, or the finding is discarded. Change nothing in the working
tree; copy to `/tmp` and ablate the copy.

Report a numbered list of findings (defect, the ablation command + its literal before/after stdout,
correction). A leg that finds nothing is weak evidence — probe hard, especially with a FORM ablation.

## Artifact
`research/pde_ledger_v3/scripts/S11c_b_row_residual.py` (run:
`python3 research/pde_ledger_v3/scripts/S11c_b_row_residual.py`). It imports the reviewed layer
`scripts/S11c_b_adjudicated_comparison.py` and reads the two committed engine transcripts
(`scripts/out/S11c_b_brane_operator_sympy_audit.out`, `mathematica/out/S11c_b_brane_operator_mathematica_audit.out`).

## What the script is supposed to compute
For each aligned semantic operator row the layer exposes (`U_MOMENTUM_ROWS`, `THICKNESS_ROW`,
`MASS_EVOLUTION_ROW`, `MU_THETA`), the admissibility operand (componentwise), and the coupling kernel (both
cross-sector blocks), per anchoring branch × density representative: the requested-truncated PY operand, the
WL operand (already truncated at emission), and their residual under the correct equivalence —
**exact** for strong slab rows, **modulo exact total in-plane divergence** for the coupling kernel,
**componentwise** for admissibility. It must PRINT operands and residuals and never assert a physics
conclusion. Full spec of intent: `directives/S11c_b_row_residual_instrument_build_directive.md`; the physics
spec it must honor: `directives/S11c_b_SHARED_PHYSICS.md` §1d, §2a, §3a, §3b, §3c, §3d.

## Derive independently FIRST (save your script + its literal stdout to a named absolute path)
Before opening the artifact, write your own small SymPy checks for these and save both script and stdout:
- The **requested truncation**: a term `ε^c η^a σ_W^b` is retained iff `c≤1 ∧ a≤1 ∧ b≤1`, grading by
  bookkeeper power (a jet `∂^n W_bg`, any `n≥1`, is `σ_W^1`), coefficients Taylor-linearized in `η`
  (`W_bg²→W_0²(1+2η w1)+O(η²)`; `1/(1+η w1)→1−η w1+O(η²)`). Confirm the artifact's truncation function agrees
  with an independent `sp.series` projection on a fixture carrying `η²`, `σ_W²`, and `η σ_W` terms.
- The product-rule identity `∇·(cF)=∇c·F + c∇·F` used to activate a held divergence to strong form.

## MANDATORY form ablations (only a FORM change tests physics; a coefficient rescale tests arithmetic)
Run each on a `/tmp` copy and report literal before/after stdout:
1. **Truncation bound.** Change the retention bound from `≤1` to `≤0` (and separately to `≤2`) in the
   requested-truncation function. Every truncated PY operand and every strong-row residual that carries an
   `η¹`/`σ_W¹` term MUST move. If any row residual is byte-identical under this change, that residual does
   not depend on the truncation — a defect (the operand was frozen/hardcoded, or the function is not wired
   in).
2. **Divergence activation.** Disable the product-rule carry-out of the held in-plane divergence (leave
   `∇·(cF)` unexpanded). The advective/mass-row residual MUST move. If it does not, the strong-form
   comparison is not actually activating divergences and the advective re-bucketing question is being
   answered on unexpanded forms.
3. **Weak coupling route.** For the coupling kernel, replace the modulo-total-in-plane-divergence residual
   with a strong-form (exact) difference. The coupling residual/certificate MUST move. Conversely, confirm
   the strong slab rows do NOT go through the divergence quotient (apply the coupling weak route to a strong
   slab row on a copy — a strong-row residual that changes under the divergence quotient would be masking
   §1d first-jet physics; report it).
4. **Row-assembly completeness.** Drop one aligned operand from a row's assembly. The assembly-accounting
   guard MUST fire (or the row residual MUST move). If nothing changes, the accounting is decorative.

## Probes (report each with the line number that computes it, or report it as uncomputed)
- **Asserts before emit.** Any `assert` that precedes the value it guards hides a disagreement (a
  perturbation that flips it kills the process, so a reviewer sees only PASS-or-crash). Report every `assert`
  that precedes a `ROW_RESIDUAL` (or operand) emit. The residual must be emitted, then optionally guarded —
  only arithmetic/assembly-accounting guards may hard-stop, never a physics residual.
- **Tautological guards (rule 2 corollary 3).** Is any emitted "check" a residual of two operands produced
  by a **single** route (zero by construction for any input)? In particular: a reconstruction identity
  `operand − polynomial − remainder`; a WL "partition check" whose target is the same local expressions the
  buckets came from. Name each and state whether an independent second route exists.
- **Hardcoded / leaked expected values.** Does any tag name encode a value/sign/grade? Is any emission
  conditional on a residual's **value** (as opposed to which row/case/engine/block)? Is any expected residual
  hardcoded anywhere? Quote it.
- **ε order.** Confirm the `ε` order `c` is attached from family metadata (`c=1` wave rows / coupling,
  `c=0` admissibility), NOT by counting `ε` symbols (WL carries `ε` as metadata; PY's `ε` is stripped by the
  layer's `coeff_epsilon`). A row misgraded to `c=0` because no `ε` symbol was found is a defect.
- **Coupling completeness.** Are BOTH cross-sector blocks (transverse→thickness and thickness→transverse)
  and their relabelled adjoints emitted, or is the reverse block manufactured as "± forward"? (§3c requires
  both.)
- **Admissibility completeness.** Is the admissibility operand compared componentwise over all body DOFs
  (3 momentum, θ, e_W) and each per-face traction, all four anchoring/density cases, with `c=0`? Or only the
  θ body force?

## Physics filter
Report a finding only if it catches a way the instrument computes the wrong object, masks a real
disagreement, manufactures a false agreement, or leaks/hardcodes the answer. Do not report "would be wrong
on a different input" style style-nits.

## Ablation sandbox
Copy the artifact (and, if you must, the layer) to `/tmp` and ablate the copy. ⛔ Never modify the working
tree. Save every ablation script AND its literal stdout to named absolute paths and report those paths.
